//--------------------------------------------------------------------------------------------------
// Aligning reads to multigenomes by extending exact matches based on read-multigenome edit distance.
// Finding exact matches betwwen reads and multigenomes using exact search with FM index.
// Exact search is perfomed with regard to a random position on reads.
// Finding inexact matches betwwen reads and multigenomes by extending FM-index based exact matches
// using edit distance between reads and multigenomes.
// Determining whether an interval on multigenomes contains a SNP position using interpolation search.
// Copyright 2014 Nam Sy Vo
//--------------------------------------------------------------------------------------------------

package isc

import (
	"github.com/vtphan/fmi" //to use FM index
	"sort"
)

//Index for SNP caller
type Index struct {
	SEQ            []byte            //store reference multigenomes
	SNP_PROF       map[int][][]byte  //hash table of SNP Profile (position, snps)
	SNP_AF         map[int][]float32 //allele frequency of SNP Profile (position, af of snps)
	SAME_LEN_SNP   map[int]int       //hash table to indicate if SNPs has same length
	SORTED_SNP_POS []int             //sorted array of SNP positions
	REV_FMI        fmi.Index         //FM-index of reverse multigenomes
}

//--------------------------------------------------------------------------------------------------
// Init function sets initial values for global variables and parameters for Index object
//--------------------------------------------------------------------------------------------------
func (I *Index) Init() {

	I.SEQ = LoadMultigenome(INPUT_INFO.Genome_file)
	PrintMemStats("memstats after loading multigenome")

	I.SNP_PROF, I.SNP_AF, I.SAME_LEN_SNP = LoadSNPLocation(INPUT_INFO.SNP_file)
	PrintMemStats("memstats after loading SNP profile")

	I.SORTED_SNP_POS = make([]int, 0, len(I.SNP_PROF))
	for k := range I.SNP_PROF {
		I.SORTED_SNP_POS = append(I.SORTED_SNP_POS, k)
	}
	sort.Sort(sort.IntSlice(I.SORTED_SNP_POS))
	PrintMemStats("memstats after loading sorted SNP postions")

	I.REV_FMI = *fmi.Load(INPUT_INFO.Rev_index_file)
	PrintMemStats("memstats after loading index of reverse multigenome")
}

//--------------------------------------------------------------------------------------------------
// Bachward Search with FM-index, start from any position on the pattern.
//--------------------------------------------------------------------------------------------------
func (I *Index) BackwardSearchFrom(index fmi.Index, pattern []byte, start_pos int) (int, int, int) {
	var sp, ep, offset uint32
	var ok bool

	c := pattern[start_pos]
	sp, ok = index.C[c]
	if !ok {
		return 0, -1, -1
	}
	ep = index.EP[c]
	var sp0, ep0 uint32
	for i := start_pos - 1; i >= 0; i-- {
		//fmt.Println("pos, # candidates: ", i, ep - sp + 1)
		c = pattern[i]
		offset, ok = index.C[c]
		if ok {
			sp0 = offset + index.OCC[c][sp-1]
			ep0 = offset + index.OCC[c][ep] - 1
			if sp0 <= ep0 {
				sp = sp0
				ep = ep0
			} else {
				return int(sp), int(ep), i + 1
			}
		} else {
			return int(sp), int(ep), i + 1
			//return 0, -1, -1
		}
	}
	return int(sp), int(ep), 0
}

//--------------------------------------------------------------------------------------------------
// FindSeeds function returns positions and distances of LCS between reads and multi-genomes.
// It uses both backward search and forward search (backward search on reverse references).
//--------------------------------------------------------------------------------------------------
func (I *Index) FindSeeds(read, rev_read []byte, p int, m_pos []int) (int, int, int, bool) {

	var rev_sp, rev_ep int = 0, INPUT_INFO.Max_snum
	var rev_s_pos, rev_e_pos, s_pos, e_pos int

	rev_s_pos = len(read) - 1 - p
	rev_sp, rev_ep, rev_e_pos = I.BackwardSearchFrom(I.REV_FMI, rev_read, rev_s_pos)
	if rev_e_pos >= 0 {
		var idx int
		//convert rev_e_pos in forward search to s_pos in backward search
		s_pos = len(read) - 1 - rev_e_pos
		e_pos = p
		if rev_ep - rev_sp + 1 <= INPUT_INFO.Max_snum && s_pos - e_pos >= INPUT_INFO.Min_slen {
			for idx = rev_sp; idx <= rev_ep; idx++ {
				m_pos[idx - rev_sp] = len(I.SEQ) - 1 - int(I.REV_FMI.SA[idx]) - (s_pos - e_pos)
			}
			return s_pos, e_pos, rev_ep - rev_sp + 1, true
		}
		return s_pos, e_pos, rev_ep - rev_sp + 1, false
	}
	return -1, -1, -1, false // will be changed later
}

//-----------------------------------------------------------------------------------------------------
// FindExtension function returns alignment (snp report) between between reads and multi-genomes.
// The alignment is built within a given threshold of distance.
//-----------------------------------------------------------------------------------------------------
func (I *Index) FindExtensions(read, qual []byte, s_pos, e_pos int, m_pos int, align_info *AlignInfo) (int, float64, 
	[]int, [][]byte, []int, []int, [][]byte, []int, int, int, bool) {

	var ref_l_flank, ref_r_flank, read_l_flank, read_r_flank, qual_l_flank, qual_r_flank []byte
	var lcs_len int = s_pos - e_pos + 1

	var isSNP, isSameLenSNP bool
	l_ext_add_len, right_ext_add_len := 0, 0
	i := 0
	for i = m_pos - e_pos; i < m_pos; i++ {
		_, isSNP = I.SNP_PROF[i]
		_, isSameLenSNP = I.SAME_LEN_SNP[i]
		if isSNP && !isSameLenSNP {
			l_ext_add_len++
		}
	}
	for i = m_pos + lcs_len; i < (m_pos + lcs_len) + (len(read) - s_pos) - 1; i++ {
		_, isSNP = I.SNP_PROF[i]
		_, isSameLenSNP = I.SAME_LEN_SNP[i]
		if isSNP && !isSameLenSNP {
			right_ext_add_len++
		}
	}
	l_most_pos := m_pos - e_pos - l_ext_add_len
	if l_most_pos >= 0 {
		ref_l_flank = I.SEQ[l_most_pos : m_pos]
	} else {
		ref_l_flank = I.SEQ[0 : m_pos]
	}
	r_most_pos := (m_pos + lcs_len) + (len(read) - s_pos) - 1 + right_ext_add_len
	if r_most_pos <= len(I.SEQ) {
		ref_r_flank = I.SEQ[m_pos+lcs_len : r_most_pos]
	} else {
		ref_r_flank = I.SEQ[m_pos+lcs_len : len(I.SEQ)]
	}

	read_l_flank = read[ : e_pos]
	qual_l_flank = qual[ : e_pos]
	left_d, left_D, l_m, l_n, l_snp_pos, l_snp_val, l_snp_idx, l_prob_h :=
		I.BackwardDistance(read_l_flank, qual_l_flank, ref_l_flank, l_most_pos, align_info.Bw_Dis, align_info.Bw_Trace)

	read_r_flank = read[s_pos + 1 : ]
	qual_r_flank = qual[s_pos + 1 : ]
	right_d, right_D, r_m, r_n, r_snp_pos, r_snp_val, r_snp_idx, r_prob_h :=
		I.ForwardDistance(read_r_flank, qual_r_flank, ref_r_flank, m_pos + lcs_len, align_info.Fw_Dis, align_info.Fw_Trace)

	dis := left_d + right_d + left_D + right_D
	prob := 0.0
	if dis <= PARA_INFO.Dist_thres {
		l_pos, l_val, l_idx, l_prob_e := I.BackwardTraceBack(read_l_flank, qual_l_flank, ref_l_flank, l_m, l_n, l_most_pos, align_info.Bw_Trace)
		r_pos, r_val, r_idx, r_prob_e := I.ForwardTraceBack(read_r_flank, qual_r_flank, ref_r_flank, r_m, r_n, m_pos + lcs_len, align_info.Fw_Trace)
		l_snp_pos = append(l_snp_pos, l_pos...)
		r_snp_pos = append(r_snp_pos, r_pos...)
		l_snp_val = append(l_snp_val, l_val...)
		r_snp_val = append(r_snp_val, r_val...)
		l_snp_idx = append(l_snp_idx, l_idx...)
		r_snp_idx = append(r_snp_idx, r_idx...)
		prob = l_prob_h * l_prob_e * r_prob_h * r_prob_e
		return dis, prob, l_snp_pos, l_snp_val, l_snp_idx, r_snp_pos, r_snp_val, r_snp_idx, l_most_pos, r_most_pos, true
	}
	return dis, prob, l_snp_pos, l_snp_val, l_snp_idx, r_snp_pos, r_snp_val, r_snp_idx, l_most_pos, r_most_pos, false
}
