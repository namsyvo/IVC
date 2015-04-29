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
	//"fmt"
)

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
	var i int
	for i = start_pos - 1; i >= 0 && i >= start_pos - INPUT_INFO.Max_slen; i-- {
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
	return int(sp), int(ep), i + 1
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
func (S *SNP_Prof) FindExtensions(read, qual []byte, s_pos, e_pos int, m_pos int, align_info *AlignInfo) (float64, 
	[]int, [][]byte, [][]byte, []int, [][]byte, [][]byte, int, int, bool) {

	var ref_l_flank, ref_r_flank, read_l_flank, read_r_flank, qual_l_flank, qual_r_flank []byte
	var isSNP, isSameLenSNP bool

	l_ext_add_len := 0
	l_most_pos := m_pos - e_pos - l_ext_add_len
	for i := m_pos - e_pos; i < m_pos; i++ {
		_, isSNP = INDEX.SNP_PROF[i]
		_, isSameLenSNP = INDEX.SAME_LEN_SNP[i]
		if isSNP && !isSameLenSNP {
			l_ext_add_len++
		}
	}
	if l_most_pos >= 0 {
		ref_l_flank = INDEX.SEQ[l_most_pos : m_pos]
	} else {
		ref_l_flank = INDEX.SEQ[0 : m_pos]
	}
	read_l_flank, qual_l_flank = read[ : e_pos], qual[ : e_pos]
	left_d, left_D, l_bt_mat, l_m, l_n, l_snp_pos, l_snp_base, l_snp_qual :=
		S.BackwardDistance(read_l_flank, qual_l_flank, ref_l_flank, l_most_pos, align_info.Bw_Dist_D, 
			align_info.Bw_Dist_IS, align_info.Bw_Dist_IT, align_info.Bw_Trace_D, align_info.Bw_Trace_IS, align_info.Bw_Trace_IT)

	right_ext_add_len := 0
	r_most_pos := m_pos + s_pos - e_pos + 1
	for i := r_most_pos; i < r_most_pos + (len(read) - s_pos) - 1; i++ {
		_, isSNP = INDEX.SNP_PROF[i]
		_, isSameLenSNP = INDEX.SAME_LEN_SNP[i]
		if isSNP && !isSameLenSNP {
			right_ext_add_len++
		}
	}
	if r_most_pos + (len(read) - s_pos) - 1 + right_ext_add_len <= len(INDEX.SEQ) {
		ref_r_flank = INDEX.SEQ[r_most_pos : r_most_pos + (len(read) - s_pos) - 1 + right_ext_add_len]
	} else {
		ref_r_flank = INDEX.SEQ[r_most_pos : len(INDEX.SEQ)]
	}
	read_r_flank, qual_r_flank = read[s_pos + 1 : ], qual[s_pos + 1 : ]
	right_d, right_D, r_bt_mat, r_m, r_n, r_snp_pos, r_snp_base, r_snp_qual :=
		S.ForwardDistance(read_r_flank, qual_r_flank, ref_r_flank, r_most_pos, align_info.Fw_Dist_D, 
			align_info.Fw_Dist_IS, align_info.Fw_Dist_IT, align_info.Fw_Trace_D, align_info.Fw_Trace_IS, align_info.Fw_Trace_IT)

	prob := left_d + right_d + left_D + right_D
	if prob <= PARA_INFO.Prob_thres {
		if l_m < len(read_l_flank) && l_n < len(ref_l_flank) {
			l_pos, l_base, l_qual := S.BackwardTraceBack(read_l_flank, qual_l_flank, ref_l_flank, l_m, l_n, l_most_pos, l_bt_mat, 
				align_info.Bw_Trace_D, align_info.Bw_Trace_IS, align_info.Bw_Trace_IT)
			l_snp_pos = append(l_snp_pos, l_pos...)
			l_snp_base = append(l_snp_base, l_base...)
			l_snp_qual = append(l_snp_qual, l_qual...)
		}
		if r_m < len(read_r_flank) && r_n < len(ref_r_flank) {
			r_pos, r_base, r_qual := S.ForwardTraceBack(read_r_flank, qual_r_flank, ref_r_flank, r_m, r_n, r_most_pos, r_bt_mat, 
				align_info.Fw_Trace_D, align_info.Fw_Trace_IS, align_info.Fw_Trace_IT)
			r_snp_pos = append(r_snp_pos, r_pos...)
			r_snp_base = append(r_snp_base, r_base...)
			r_snp_qual = append(r_snp_qual, r_qual...)
		}
		return prob, l_snp_pos, l_snp_base, l_snp_qual, r_snp_pos, r_snp_base, r_snp_qual, l_most_pos, r_most_pos, true
	}
	return prob, l_snp_pos, l_snp_base, l_snp_qual, r_snp_pos, r_snp_base, r_snp_qual, l_most_pos, r_most_pos, false
}
