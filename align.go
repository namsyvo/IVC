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
	//"fmt"
	"sort"
	"github.com/vtphan/fmi" //to use FM index
)

//--------------------------------------------------------------------------------------------------
//Index for SNP caller
//--------------------------------------------------------------------------------------------------
type Index struct {
	SEQ            []byte            //store reference multigenomes
	SNP_PROF       map[int][][]byte  //hash table of SNP Profile (position, snps)
	SNP_AF         map[int][]float32 //allele frequency of SNP Profile (position, af of snps)
	SAME_LEN_SNP   map[int]int       //hash table to indicate if SNPs has same length
	DEL_SNP   	   map[int]int       //hash table to store length of deletions if SNPs are deletion
	SORTED_SNP_POS []int             //sorted array of SNP positions
	REV_FMI        fmi.Index         //FM-index of reverse multigenomes
}

//--------------------------------------------------------------------------------------------------
// Init function sets initial values for global variables and parameters for Index object
//--------------------------------------------------------------------------------------------------
func New_Index() *Index {

	I := new(Index)

	I.SEQ = LoadMultigenome(INPUT_INFO.Genome_file)
	PrintMemStats("memstats after loading multigenome")

	I.SNP_PROF, I.SNP_AF = LoadVarProf(INPUT_INFO.SNP_file)
	PrintMemStats("memstats after loading SNP profile")

	I.SAME_LEN_SNP = make(map[int]int)
	I.DEL_SNP = make(map[int]int)
	I.SORTED_SNP_POS = make([]int, 0, len(I.SNP_PROF))
	var same_len_flag, del_flag bool
	var snp_len int
	for snp_pos, snp_prof := range I.SNP_PROF {
		I.SORTED_SNP_POS = append(I.SORTED_SNP_POS, snp_pos)
		snp_len = len(snp_prof[0])
		same_len_flag, del_flag = true, true
		for _, val := range snp_prof[1: ] {
			if snp_len != len(val) {
				same_len_flag = false
			}
			if snp_len <= len(val) {
				del_flag = false
			}
		}
		if same_len_flag {
			I.SAME_LEN_SNP[snp_pos] = snp_len
		}
		if del_flag {
			I.DEL_SNP[snp_pos] = snp_len - 1
		}		
	}
	sort.Sort(sort.IntSlice(I.SORTED_SNP_POS))
	PrintMemStats("Memstats after creating auxiliary data structures for SNP Profile")

	I.REV_FMI = *fmi.Load(INPUT_INFO.Rev_index_file)
	PrintMemStats("memstats after loading index of reverse multigenome")

	return I
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
