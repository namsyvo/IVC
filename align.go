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
func (I *Index) FindSeeds(read, rev_read []byte, p int, match_pos []int) (int, int, int, bool) {

	var rev_sp, rev_ep int = 0, PARA_INFO.Max_match
	var rev_s_pos, rev_e_pos, s_pos, e_pos int

	rev_s_pos = len(read) - 1 - p
	rev_sp, rev_ep, rev_e_pos = I.BackwardSearchFrom(I.REV_FMI, rev_read, rev_s_pos)
	if rev_e_pos >= 0 {
		var idx int
		//convert rev_e_pos in forward search to s_pos in backward search
		s_pos = len(read) - 1 - rev_e_pos
		e_pos = p
		if rev_ep - rev_sp + 1 <= PARA_INFO.Max_match && s_pos - e_pos >= PARA_INFO.Min_seed {
			for idx = rev_sp; idx <= rev_ep; idx++ {
				match_pos[idx - rev_sp] = len(I.SEQ) - 1 - int(I.REV_FMI.SA[idx]) - (s_pos - e_pos)
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
func (I *Index) FindExtensions(read []byte, s_pos, e_pos int, match_pos int, align_info *AlignInfo) (int, []int, [][]byte, []int, []int, [][]byte, []int, int, int) {

	var ref_left_flank, ref_right_flank, read_left_flank, read_right_flank []byte
	var lcs_len int = s_pos - e_pos + 1

	var isSNP, isSameLenSNP bool
	left_ext_add_len, right_ext_add_len := 0, 0
	i := 0
	for i = match_pos - e_pos; i < match_pos; i++ {
		_, isSNP = I.SNP_PROF[i]
		_, isSameLenSNP = I.SAME_LEN_SNP[i]
		if isSNP && !isSameLenSNP {
			left_ext_add_len++
		}
	}
	for i = match_pos + lcs_len; i < (match_pos + lcs_len) + (len(read) - s_pos) - 1; i++ {
		_, isSNP = I.SNP_PROF[i]
		_, isSameLenSNP = I.SAME_LEN_SNP[i]
		if isSNP && !isSameLenSNP {
			right_ext_add_len++
		}
	}
	left_most_pos := match_pos - e_pos - left_ext_add_len
	if left_most_pos >= 0 {
		ref_left_flank = I.SEQ[left_most_pos:match_pos]
	} else {
		ref_left_flank = I.SEQ[0 : match_pos]
	}
	right_most_pos := (match_pos + lcs_len) + (len(read) - s_pos) - 1 + right_ext_add_len
	if right_most_pos <= len(I.SEQ) {
		ref_right_flank = I.SEQ[match_pos+lcs_len : right_most_pos]
	} else {
		ref_right_flank = I.SEQ[match_pos+lcs_len : len(I.SEQ)]
	}

	read_left_flank = read[ : e_pos]
	left_d, left_D, left_m, left_n, left_snp_pos, left_snp_val, left_snp_idx :=
		I.BackwardDistance(read_left_flank, ref_left_flank, left_most_pos, align_info.Bw_Dis, align_info.Bw_Trace)

	read_right_flank = read[s_pos + 1 : ]
	right_d, right_D, right_m, right_n, right_snp_pos, right_snp_val, right_snp_idx :=
		I.ForwardDistance(read_right_flank, ref_right_flank, match_pos+lcs_len, align_info.Fw_Dis, align_info.Fw_Trace)

	dis := left_d + right_d + left_D + right_D
	if dis <= PARA_INFO.Dist_thres {
		left_pos, left_val, left_idx := I.BackwardTraceBack(read_left_flank, ref_left_flank,
			left_m, left_n, left_most_pos, align_info.Bw_Trace)
		right_pos, right_val, right_idx := I.ForwardTraceBack(read_right_flank, ref_right_flank,
			right_m, right_n, match_pos + lcs_len, align_info.Fw_Trace)
		left_snp_pos = append(left_snp_pos, left_pos...)
		right_snp_pos = append(right_snp_pos, right_pos...)
		left_snp_val = append(left_snp_val, left_val...)
		right_snp_val = append(right_snp_val, right_val...)
		left_snp_idx = append(left_snp_idx, left_idx...)
		right_snp_idx = append(right_snp_idx, right_idx...)
	}
	return dis, left_snp_pos, left_snp_val, left_snp_idx, right_snp_pos, right_snp_val, right_snp_idx, left_most_pos, right_most_pos
}
