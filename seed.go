//--------------------------------------------------------------------------------------------------
// IVC: seed.go - Finding seeds of alignment betwwen reads and multigenomes using FM index.
// Searching is perfomed from a random position on reads.
// Copyright 2015 Nam Sy Vo.
//--------------------------------------------------------------------------------------------------

package ivc

import (
	//"fmt"
	"github.com/vtphan/fmi" //to use FM index
)

//--------------------------------------------------------------------------------------------------
//Index represents info and index used for alignment process.
//--------------------------------------------------------------------------------------------------
type Index struct {
	Seq        []byte            //store reference multigenomes
	VarProf    map[int][][]byte  //hash table of SNP Profile (position, variants)
	VarAF      map[int][]float32 //allele frequency of SNP Profile (position, af of variants)
	SameLenVar map[int]int       //hash table to indicate if SNPs has same length
	DelVar     map[int]int       //hash table to store length of deletions if SNPs are deletion
	RevFMI     fmi.Index         //FM-index of reverse multigenomes
}

//--------------------------------------------------------------------------------------------------
// NewIndex creates an instance of Index and seta up its variables.
//--------------------------------------------------------------------------------------------------
func NewIndex() *Index {

	I := new(Index)

	I.Seq = LoadMultigenome(INPUT_INFO.Ref_file)
	PrintMemStats("Memstats after loading multigenome")

	I.VarProf, I.VarAF = LoadVarProf(INPUT_INFO.Var_prof_file)
	PrintMemStats("Memstats after loading variant profile")

	I.SameLenVar = make(map[int]int)
	I.DelVar = make(map[int]int)
	var same_len_flag, del_flag bool
	var var_len int
	for var_pos, var_prof := range I.VarProf {
		var_len = len(var_prof[0])
		same_len_flag, del_flag = true, true
		for _, val := range var_prof[1:] {
			if var_len != len(val) {
				same_len_flag = false
			}
			if var_len <= len(val) {
				del_flag = false
			}
		}
		if same_len_flag {
			I.SameLenVar[var_pos] = var_len
		}
		if del_flag {
			I.DelVar[var_pos] = var_len - 1
		}
	}
	PrintMemStats("Memstats after creating auxiliary data structures for variant profile")

	I.RevFMI = *fmi.Load(INPUT_INFO.Rev_index_file)
	PrintMemStats("Memstats after loading index of reverse multigenome")

	return I
}

//--------------------------------------------------------------------------------------------------
// BackwardSearchFrom searches for exact matches between a pattern and the reference using FM-index,
// It starts to search from any position on the pattern.
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
	for i = start_pos - 1; i >= 0 && i >= start_pos-INPUT_INFO.Max_slen; i-- {
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
// FindSeeds returns positions and distances of seeds between a read and the reference.
// It uses both backward search and forward search
// Forward search is backward search on reverse of the reference.
//--------------------------------------------------------------------------------------------------
func (I *Index) FindSeeds(read, rev_read []byte, p int, m_pos []int) (int, int, int, bool) {

	var rev_sp, rev_ep int = 0, INPUT_INFO.Max_snum
	var rev_s_pos, rev_e_pos, s_pos, e_pos int

	rev_s_pos = len(read) - 1 - p
	rev_sp, rev_ep, rev_e_pos = I.BackwardSearchFrom(I.RevFMI, rev_read, rev_s_pos)
	if rev_e_pos >= 0 {
		var idx int
		//convert rev_e_pos in forward search to s_pos in backward search
		s_pos = len(read) - 1 - rev_e_pos
		e_pos = p
		if rev_ep-rev_sp+1 <= INPUT_INFO.Max_snum && s_pos-e_pos >= INPUT_INFO.Min_slen {
			for idx = rev_sp; idx <= rev_ep; idx++ {
				m_pos[idx-rev_sp] = len(I.Seq) - 1 - int(I.RevFMI.SA[idx]) - (s_pos - e_pos)
			}
			return s_pos, e_pos, rev_ep - rev_sp + 1, true
		}
		return s_pos, e_pos, rev_ep - rev_sp + 1, false
	}
	return -1, -1, -1, false // will be changed later
}
