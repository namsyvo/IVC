//--------------------------------------------------------------------------------------------------
// IVC: seed.go - Finding seeds of alignment betwwen reads and multigenomes using FM index.
// Searching is perfomed from a random position on reads.
// Copyright 2015 Nam Sy Vo.
//--------------------------------------------------------------------------------------------------

package ivc

import (
	//"fmt"
	"github.com/vtphan/fmi" //to use FM index
	"math/rand"
	"strconv"
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
func (I *Index) BackwardSearchFrom(pattern []byte, start_pos int) (int, int, int) {
	var sp, ep, offset uint32
	var ok bool

	c := pattern[start_pos]
	sp, ok = I.RevFMI.C[c]
	if !ok {
		return 0, -1, -1
	}
	ep = I.RevFMI.EP[c]
	var sp0, ep0 uint32
	var i int
	for i = start_pos - 1; i >= 0 && i >= start_pos-INPUT_INFO.Max_slen; i-- {
		c = pattern[i]
		offset, ok = I.RevFMI.C[c]
		if ok {
			sp0 = offset + I.RevFMI.OCC[c][sp-1]
			ep0 = offset + I.RevFMI.OCC[c][ep] - 1
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
	rev_sp, rev_ep, rev_e_pos = I.BackwardSearchFrom(rev_read, rev_s_pos)
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

//---------------------------------------------------------------------------------------------------
// FindSeedsFromPairedEnds find all pairs of seeds which have proper chromosome distances.
//---------------------------------------------------------------------------------------------------
func (I *Index) FindSeedsPE(read_info *ReadInfo, seed_pos [][]int, rand_gen *rand.Rand) (*SeedInfo, *SeedInfo, bool) {

	var has_seeds_r1_or, has_seeds_r1_rc, has_seeds_r2_or, has_seeds_r2_rc bool
	var s_pos_r1_or, e_pos_r1_or, m_num_r1_or, s_pos_r1_rc, e_pos_r1_rc, m_num_r1_rc int
	var s_pos_r2_or, e_pos_r2_or, m_num_r2_or, s_pos_r2_rc, e_pos_r2_rc, m_num_r2_rc int
	var s_pos_r1, e_pos_r1, s_pos_r2, e_pos_r2, m_pos_r1, m_pos_r2 []int
	var strand_r1, strand_r2 []bool
	var i, j int

	//Take an initial position to search
	r_pos_r1_or := INPUT_INFO.Start_pos
	r_pos_r1_rc := INPUT_INFO.Start_pos
	r_pos_r2_or := INPUT_INFO.Start_pos
	r_pos_r2_rc := INPUT_INFO.Start_pos
	if INPUT_INFO.Search_mode == 1 {
		r_pos_r1_or = rand_gen.Intn(len(read_info.Read1) - 5)
		r_pos_r1_rc = rand_gen.Intn(len(read_info.Read1) - 5)
		r_pos_r2_or = rand_gen.Intn(len(read_info.Read2) - 5)
		r_pos_r2_rc = rand_gen.Intn(len(read_info.Read2) - 5)
	}

	loop_num := 1
	for loop_num <= PARA_INFO.Iter_num { //temp value, will be replaced later
		PrintLoopTraceInfo(loop_num, "FindSeedsFromPairedEnds, First:\t"+string(read_info.Read1))
		PrintLoopTraceInfo(loop_num, "FindSeedsFromPairedEnds, Second:\t"+string(read_info.Read2))

		PrintMemStats("Before FindSeeds, loop_num " + strconv.Itoa(loop_num))
		s_pos_r1_or, e_pos_r1_or, m_num_r1_or, has_seeds_r1_or =
			I.FindSeeds(read_info.Read1, read_info.Rev_read1, r_pos_r1_or, seed_pos[0])
		PrintSeedTraceInfo("r1_or", e_pos_r1_or, s_pos_r1_or, read_info.Read1)
		if has_seeds_r1_or {
			PrintExtendTraceInfo("r1_or", read_info.Read1[e_pos_r1_or:s_pos_r1_or+1],
				e_pos_r1_or, s_pos_r1_or, m_num_r1_or, seed_pos[0])
		}
		s_pos_r1_rc, e_pos_r1_rc, m_num_r1_rc, has_seeds_r1_rc =
			I.FindSeeds(read_info.Rev_comp_read1, read_info.Comp_read1, r_pos_r1_rc, seed_pos[1])
		PrintSeedTraceInfo("r1_rc", e_pos_r1_rc, s_pos_r1_rc, read_info.Rev_comp_read1)
		if has_seeds_r1_rc {
			PrintExtendTraceInfo("r1_rc", read_info.Rev_comp_read1[e_pos_r1_rc:s_pos_r1_rc+1],
				e_pos_r1_rc, s_pos_r1_rc, m_num_r1_rc, seed_pos[1])
		}
		s_pos_r2_or, e_pos_r2_or, m_num_r2_or, has_seeds_r2_or =
			I.FindSeeds(read_info.Read2, read_info.Rev_read2, r_pos_r2_or, seed_pos[2])
		PrintSeedTraceInfo("r2_or", e_pos_r2_or, s_pos_r2_or, read_info.Read2)
		if has_seeds_r2_or {
			PrintExtendTraceInfo("r2_or", read_info.Read1[e_pos_r2_or:s_pos_r2_or+1],
				e_pos_r2_or, s_pos_r2_or, m_num_r2_or, seed_pos[2])
		}
		s_pos_r2_rc, e_pos_r2_rc, m_num_r2_rc, has_seeds_r2_rc =
			I.FindSeeds(read_info.Rev_comp_read2, read_info.Comp_read2, r_pos_r2_rc, seed_pos[3])
		PrintSeedTraceInfo("r2_rc", e_pos_r2_rc, s_pos_r2_rc, read_info.Rev_comp_read2)
		if has_seeds_r2_rc {
			PrintExtendTraceInfo("r2_rc", read_info.Rev_comp_read2[e_pos_r2_rc:s_pos_r2_rc+1],
				e_pos_r2_rc, s_pos_r2_rc, m_num_r2_rc, seed_pos[3])
		}
		PrintMemStats("After FindSeeds, loop_num " + strconv.Itoa(loop_num))

		if has_seeds_r1_or && has_seeds_r2_rc {
			PrintExtendTraceInfo("r1_or(F1R2)", read_info.Read1[e_pos_r1_or:s_pos_r1_or+1],
				e_pos_r1_or, s_pos_r1_or, m_num_r1_or, seed_pos[0])
			PrintExtendTraceInfo("r2_rc(F1R2)", read_info.Read2[e_pos_r2_rc:s_pos_r2_rc+1],
				e_pos_r2_rc, s_pos_r2_rc, m_num_r2_rc, seed_pos[3])
			for i = 0; i < m_num_r1_or; i++ {
				for j = 0; j < m_num_r2_rc; j++ {
					//Check if alignments are likely pair-end alignments
					if (seed_pos[3][j]-seed_pos[0][i]) >= PARA_INFO.Read_len &&
						(seed_pos[3][j]-seed_pos[0][i]) <= PARA_INFO.Read_len+PARA_INFO.Max_ins {
						PrintPairedSeedInfo("r1_or, r2_rc, paired pos", seed_pos[0][i], seed_pos[3][j])
						s_pos_r1 = append(s_pos_r1, s_pos_r1_or)
						e_pos_r1 = append(e_pos_r1, e_pos_r1_or)
						s_pos_r2 = append(s_pos_r2, s_pos_r2_rc)
						e_pos_r2 = append(e_pos_r2, e_pos_r2_rc)
						m_pos_r1 = append(m_pos_r1, seed_pos[0][i])
						m_pos_r2 = append(m_pos_r2, seed_pos[3][j])
						strand_r1 = append(strand_r1, true)
						strand_r2 = append(strand_r2, false)
					}
				}
			}
		}
		if has_seeds_r1_rc && has_seeds_r2_or {
			PrintExtendTraceInfo("r1_rc (F2R1)", read_info.Read1[e_pos_r1_rc:s_pos_r1_rc+1],
				e_pos_r1_rc, s_pos_r1_rc, m_num_r1_rc, seed_pos[1])
			PrintExtendTraceInfo("r2_or (F2R1)", read_info.Read2[e_pos_r2_or:s_pos_r2_or+1],
				e_pos_r2_or, s_pos_r2_or, m_num_r2_or, seed_pos[2])
			for i = 0; i < m_num_r1_rc; i++ {
				for j = 0; j < m_num_r2_or; j++ {
					//Check if alignments are likely pair-end alignments
					if (seed_pos[1][i]-seed_pos[2][j]) >= PARA_INFO.Read_len &&
						(seed_pos[1][i]-seed_pos[2][j]) <= PARA_INFO.Read_len+PARA_INFO.Max_ins {
						PrintPairedSeedInfo("r1_rc, r2_or, paired pos", seed_pos[1][i], seed_pos[2][j])
						s_pos_r1 = append(s_pos_r1, s_pos_r1_rc)
						e_pos_r1 = append(e_pos_r1, e_pos_r1_rc)
						s_pos_r2 = append(s_pos_r2, s_pos_r2_or)
						e_pos_r2 = append(e_pos_r2, e_pos_r2_or)
						m_pos_r1 = append(m_pos_r1, seed_pos[1][i])
						m_pos_r2 = append(m_pos_r2, seed_pos[2][j])
						strand_r1 = append(strand_r1, false)
						strand_r2 = append(strand_r2, true)
					}
				}
			}
		}
		if len(s_pos_r1) >= 1 && len(s_pos_r1) <= INPUT_INFO.Max_psnum {
			return &SeedInfo{s_pos_r1, e_pos_r1, m_pos_r1, strand_r1}, &SeedInfo{s_pos_r2, e_pos_r2, m_pos_r2, strand_r2}, true
		}
		//Take a new position to search
		r_pos_r1_or = r_pos_r1_or + INPUT_INFO.Search_step
		r_pos_r1_rc = r_pos_r1_rc + INPUT_INFO.Search_step
		r_pos_r2_or = r_pos_r2_or + INPUT_INFO.Search_step
		r_pos_r2_rc = r_pos_r2_rc + INPUT_INFO.Search_step
		if INPUT_INFO.Search_mode == 1 {
			r_pos_r1_or = rand_gen.Intn(len(read_info.Read1) - 5)
			r_pos_r1_rc = rand_gen.Intn(len(read_info.Read1) - 5)
			r_pos_r2_or = rand_gen.Intn(len(read_info.Read2) - 5)
			r_pos_r2_rc = rand_gen.Intn(len(read_info.Read2) - 5)
		}
		loop_num++
	}
	return &SeedInfo{s_pos_r1, e_pos_r1, m_pos_r1, strand_r1}, &SeedInfo{s_pos_r2, e_pos_r2, m_pos_r2, strand_r2}, false
}
