//--------------------------------------------------------------------------------------------------
// IVC: seed.go
// Searching for seeds of alignment betwwen reads and multigenomes.
// Searching is perfomed from a random position on reads forwardly using an FM-index of reverse multigenomes.
// Copyright 2015 Nam Sy Vo.
//--------------------------------------------------------------------------------------------------

package ivc

import (
	"math/rand"
)

//--------------------------------------------------------------------------------------------------
// ForwardSearchFrom searches for exact matches between a pattern and the reference using FM-index.
// It starts to search forwardly on the pattern from any position to match backwardly on the reference.
//--------------------------------------------------------------------------------------------------
func (VC *VarCallIndex) ForwardSearchFrom(pattern []byte, s_pos int) (int, int, int) {
	var sp, ep, offset uint32
	var ok bool

	c := pattern[s_pos]
	sp, ok = VC.RevFMI.C[c]
	if !ok {
		return -1, -1, -1
	}
	ep = VC.RevFMI.EP[c]
	var sp0, ep0 uint32
	var i, L int
	for i, L = s_pos+1, len(pattern); i < L && i <= s_pos+PARA.Max_slen; i++ {
		c = pattern[i]
		offset, ok = VC.RevFMI.C[c]
		if ok {
			sp0 = offset + VC.RevFMI.OCC[c][sp-1]
			ep0 = offset + VC.RevFMI.OCC[c][ep] - 1
			if sp0 <= ep0 {
				sp = sp0
				ep = ep0
			} else {
				return int(sp), int(ep), i - 1
			}
		} else {
			return int(sp), int(ep), i - 1
			//return 0, -1, -1
		}
	}
	return int(sp), int(ep), i - 1
}

//--------------------------------------------------------------------------------------------------
// SearchSeeds returns positions and distances of seeds between a read and the reference.
// It searches forwardly on read to match backwardly on reverse of the reference.
//--------------------------------------------------------------------------------------------------
func (VC *VarCallIndex) SearchSeeds(read []byte, s_pos int, m_pos []int) (int, int, int, bool) {

	sp, ep, e_pos := VC.ForwardSearchFrom(read, s_pos)
	if e_pos >= 0 {
		if ep-sp+1 <= PARA.Max_snum && e_pos-s_pos >= PARA.Min_slen {
			for idx := sp; idx <= ep; idx++ {
				m_pos[idx-sp] = VC.SeqLen - 1 - int(VC.RevFMI.SA[idx]) - (e_pos - s_pos)
			}
			return s_pos, e_pos, ep - sp + 1, true
		}
		return s_pos, e_pos, ep - sp + 1, false
	}
	return -1, -1, -1, false // will be changed later
}

//---------------------------------------------------------------------------------------------------
// SearchSeedsPE searches for all pairs of seeds which have proper chromosome distances.
//---------------------------------------------------------------------------------------------------
func (VC *VarCallIndex) SearchSeedsPE(read_info *ReadInfo, seed_pos [][]int, rand_gen *rand.Rand) (*SeedInfo, *SeedInfo, bool) {

	var has_seeds_r1_or, has_seeds_r1_rc, has_seeds_r2_or, has_seeds_r2_rc bool
	var s_pos_r1_or, e_pos_r1_or, m_num_r1_or, s_pos_r1_rc, e_pos_r1_rc, m_num_r1_rc int
	var s_pos_r2_or, e_pos_r2_or, m_num_r2_or, s_pos_r2_rc, e_pos_r2_rc, m_num_r2_rc int
	var s_pos_r1, e_pos_r1, s_pos_r2, e_pos_r2, m_pos_r1, m_pos_r2 []int
	var strand_r1, strand_r2 []bool
	var i, j int

	var r_pos_r1_or, r_pos_r1_rc, r_pos_r2_or, r_pos_r2_rc int
	//Take an initial position to search
	if PARA.Search_mode == 1 {
		r_pos_r1_or = rand_gen.Intn(len(read_info.Read1) - PARA.Min_slen)
		r_pos_r1_rc = rand_gen.Intn(len(read_info.Read1) - PARA.Min_slen)
		r_pos_r2_or = rand_gen.Intn(len(read_info.Read2) - PARA.Min_slen)
		r_pos_r2_rc = rand_gen.Intn(len(read_info.Read2) - PARA.Min_slen)
	} else {
		r_pos_r1_or = PARA.Start_pos
		r_pos_r1_rc = PARA.Start_pos
		r_pos_r2_or = PARA.Start_pos
		r_pos_r2_rc = PARA.Start_pos
	}
	loop_num := 1
	for loop_num <= PARA.Iter_num {
		if PARA.Debug_mode {
			PrintLoopTraceInfo(loop_num, "SearchSeedsFromPairedEnds, First:\t"+string(read_info.Read1))
			PrintLoopTraceInfo(loop_num, "SearchSeedsFromPairedEnds, Second:\t"+string(read_info.Read2))
		}
		s_pos_r1_or, e_pos_r1_or, m_num_r1_or, has_seeds_r1_or =
			VC.SearchSeeds(read_info.Read1, r_pos_r1_or, seed_pos[0])
		if PARA.Debug_mode {
			PrintSeedTraceInfo("r1_or", s_pos_r1_or, e_pos_r1_or, read_info.Read1)
			if has_seeds_r1_or {
				PrintExtendTraceInfo("r1_or", read_info.Read1[s_pos_r1_or:e_pos_r1_or+1],
					s_pos_r1_or, e_pos_r1_or, m_num_r1_or, seed_pos[0])
			}
		}
		s_pos_r1_rc, e_pos_r1_rc, m_num_r1_rc, has_seeds_r1_rc =
			VC.SearchSeeds(read_info.Rev_comp_read1, r_pos_r1_rc, seed_pos[1])
		if PARA.Debug_mode {
			PrintSeedTraceInfo("r1_rc", s_pos_r1_rc, e_pos_r1_rc, read_info.Rev_comp_read1)
			if has_seeds_r1_rc {
				PrintExtendTraceInfo("r1_rc", read_info.Rev_comp_read1[s_pos_r1_rc:e_pos_r1_rc+1],
					s_pos_r1_rc, e_pos_r1_rc, m_num_r1_rc, seed_pos[1])
			}
		}
		s_pos_r2_or, e_pos_r2_or, m_num_r2_or, has_seeds_r2_or =
			VC.SearchSeeds(read_info.Read2, r_pos_r2_or, seed_pos[2])
		if PARA.Debug_mode {
			PrintSeedTraceInfo("r2_or", s_pos_r2_or, e_pos_r2_or, read_info.Read2)
			if has_seeds_r2_or {
				PrintExtendTraceInfo("r2_or", read_info.Read1[s_pos_r2_or:e_pos_r2_or+1],
					s_pos_r2_or, e_pos_r2_or, m_num_r2_or, seed_pos[2])
			}
		}
		s_pos_r2_rc, e_pos_r2_rc, m_num_r2_rc, has_seeds_r2_rc =
			VC.SearchSeeds(read_info.Rev_comp_read2, r_pos_r2_rc, seed_pos[3])
		if PARA.Debug_mode {
			PrintSeedTraceInfo("r2_rc", s_pos_r2_rc, e_pos_r2_rc, read_info.Rev_comp_read2)
			if has_seeds_r2_rc {
				PrintExtendTraceInfo("r2_rc", read_info.Rev_comp_read2[s_pos_r2_rc:e_pos_r2_rc+1],
					s_pos_r2_rc, e_pos_r2_rc, m_num_r2_rc, seed_pos[3])
			}
		}
		if has_seeds_r1_or && has_seeds_r2_rc {
			if PARA.Debug_mode {
				PrintExtendTraceInfo("r1_or(F1R2)", read_info.Read1[s_pos_r1_or:e_pos_r1_or+1],
					s_pos_r1_or, e_pos_r1_or, m_num_r1_or, seed_pos[0])
				PrintExtendTraceInfo("r2_rc(F1R2)", read_info.Read2[s_pos_r2_rc:e_pos_r2_rc+1],
					s_pos_r2_rc, e_pos_r2_rc, m_num_r2_rc, seed_pos[3])
			}
			for i = 0; i < m_num_r1_or; i++ {
				for j = 0; j < m_num_r2_rc; j++ {
					//Check if alignments are likely pair-end alignments
					if (seed_pos[3][j]-seed_pos[0][i]) >= PARA.Read_len &&
						(seed_pos[3][j]-seed_pos[0][i]) <= PARA.Read_len+PARA.Max_ins {
						if PARA.Debug_mode {
							PrintPairedSeedInfo("r1_or, r2_rc, paired pos", seed_pos[0][i], seed_pos[3][j])
						}
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
			if PARA.Debug_mode {
				PrintExtendTraceInfo("r1_rc (F2R1)", read_info.Read1[s_pos_r1_rc:e_pos_r1_rc+1],
					s_pos_r1_rc, e_pos_r1_rc, m_num_r1_rc, seed_pos[1])
				PrintExtendTraceInfo("r2_or (F2R1)", read_info.Read2[s_pos_r2_or:e_pos_r2_or+1],
					s_pos_r2_or, e_pos_r2_or, m_num_r2_or, seed_pos[2])
			}
			for i = 0; i < m_num_r1_rc; i++ {
				for j = 0; j < m_num_r2_or; j++ {
					//Check if alignments are likely pair-end alignments
					if (seed_pos[1][i]-seed_pos[2][j]) >= PARA.Read_len &&
						(seed_pos[1][i]-seed_pos[2][j]) <= PARA.Read_len+PARA.Max_ins {
						if PARA.Debug_mode {
							PrintPairedSeedInfo("r1_rc, r2_or, paired pos", seed_pos[1][i], seed_pos[2][j])
						}
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
		if len(s_pos_r1) >= 1 && len(s_pos_r1) <= PARA.Max_psnum {
			return &SeedInfo{s_pos_r1, e_pos_r1, m_pos_r1, strand_r1}, &SeedInfo{s_pos_r2, e_pos_r2, m_pos_r2, strand_r2}, true
		}
		//Take a new position to search
		if PARA.Search_mode == 1 { //random search
			r_pos_r1_or = rand_gen.Intn(len(read_info.Read1) - PARA.Min_slen)
			r_pos_r1_rc = rand_gen.Intn(len(read_info.Read1) - PARA.Min_slen)
			r_pos_r2_or = rand_gen.Intn(len(read_info.Read2) - PARA.Min_slen)
			r_pos_r2_rc = rand_gen.Intn(len(read_info.Read2) - PARA.Min_slen)
		} else {
			r_pos_r1_or = r_pos_r1_or + PARA.Search_step
			r_pos_r1_rc = r_pos_r1_rc + PARA.Search_step
			r_pos_r2_or = r_pos_r2_or + PARA.Search_step
			r_pos_r2_rc = r_pos_r2_rc + PARA.Search_step
		}
		loop_num++
	}
	return &SeedInfo{s_pos_r1, e_pos_r1, m_pos_r1, strand_r1}, &SeedInfo{s_pos_r2, e_pos_r2, m_pos_r2, strand_r2}, false
}
