//-------------------------------------------------------------------------------------------------
// Multigenome package: distance module.
// Calculating distance and determining alignment between reads and "starred" multigenomes.
// Copyright 2014 Nam Sy Vo
//-------------------------------------------------------------------------------------------------

package isc

import (
	"math"
	//"fmt"
)

//-------------------------------------------------------------------------------------------------
// Cost functions for computing distance between reads and multi-genomes.
// Input slices should have same length
//-------------------------------------------------------------------------------------------------
func AlignCostKnownLoci(read, ref, qual []byte, prob float64) float64 {
	p := 0.0
	for i := 0; i < len(read); i++ {
		if read[i] != ref[i] {
			return math.MaxFloat64
		} else {
			p = p - math.Log10(prob) - math.Log10(1.0 - math.Pow(10, -(float64(qual[i]) - 33) / 10.0))
		}
	}
	return p
}

//-------------------------------------------------------------------------------------------------
// Calculate the distance between read and ref in backward direction.
// 	read is a read.
// 	ref is part of a multi-genome.
// The reads include standard bases, the multi-genome includes standard bases and "*" characters.
//-------------------------------------------------------------------------------------------------
func (S *SNP_Prof) BackwardDistance(read, qual, ref []byte, pos int, D, IS, IT [][]float64, 
	BT_D, BT_IS, BT_IT [][][]int) (float64, float64, int, int, int, []int, [][]byte, [][]byte) {

	var snp_len int
	var snp_str string
	var is_snp, is_same_len_snp bool
	var p, min_p, snp_prob float64
	var snp_prof map[string]float64

	align_prob := 0.0
	m, n := len(read), len(ref)

	PrintEditDisInput("bw dis input: read, qual, ref", read, qual, ref)
	var snp_pos []int
	var snp_base, snp_qual [][]byte
	var sub_m float64
	for m > 0 && n > 0 {
		sub_m = NEW_SNP_RATE_LOG - math.Log10(1.0 - math.Pow(10, -(float64(qual[m - 1]) - 33) / 10.0))
		if _, is_snp = INDEX.SNP_PROF[pos + n - 1]; !is_snp {
			if read[m - 1] != ref[n - 1] {
				snp_pos = append(snp_pos, pos + n - 1)
				snp_base = append(snp_base, []byte{read[m - 1]})
				snp_qual = append(snp_qual, []byte{qual[m - 1]})
				align_prob = align_prob + sub_m
			}
			m--
			n--
		} else if snp_len, is_same_len_snp = INDEX.SAME_LEN_SNP[pos + n - 1]; is_same_len_snp {
			snp_prof, _ = S.SNP_Calls[uint32(pos + n - 1)]
			min_p = math.MaxFloat64
			for snp_str, snp_prob = range snp_prof {
				if m >= snp_len {
					p = AlignCostKnownLoci(read[m - snp_len : m], []byte(snp_str), qual[m - snp_len : m], snp_prob)
					if min_p > p {
						min_p = p
					}
				}
			}
			if min_p < math.MaxFloat64 {
				align_prob = align_prob + min_p
				snp_pos = append(snp_pos, pos + n - 1)
				snp, qlt := make([]byte, snp_len), make([]byte, snp_len)
				copy(snp, read[m - snp_len : m])
				copy(qlt, qual[m - snp_len : m])
				snp_base = append(snp_base, snp)
				snp_qual = append(snp_qual, qlt)
				m -= snp_len
				n--
			} else {
				break
			}
		} else {
			break
		}
		if align_prob > PARA_INFO.Prob_thres {
			return PARA_INFO.Prob_thres + 1, 0, 0, m, n, snp_pos, snp_base, snp_qual
		}
	}
	PrintDisInfo("bw H dis", m, n, align_prob)

	if m == len(read) || n == len(ref) {
		return align_prob, 0, 0, m, n, snp_pos, snp_base, snp_qual
	}

	/*
	Backtrace matrix, for each BT[i][j]:
	BT[i][j][0]: direction, can be 0: diagonal arrow (back to i-1,j-1), 1: up arrow (back to i-1,j),
	 	2: left arrow (back to i,j-1).
	BT[i][j][1]: matrix, can be 0: matrix for D, 1: matrix for IS, 2: matrix for IT.
	BT[i][j][2]: number of shift (equal to length of called variant) at known variant loc, 
		can be any integer number, e.g. 5 means back to i-5,j-1.
	*/
	var i, j int
	for i := 0; i <= 2 * PARA_INFO.Read_len; i++ {
		for j := 0; j <= 2 * PARA_INFO.Read_len; j++ {
			BT_D[i][j][0], BT_D[i][j][1], BT_D[i][j][2] = -1, -1, -1
			BT_IS[i][j][0], BT_IS[i][j][1], BT_IS[i][j][2] = -1, -1, -1
			BT_IT[i][j][0], BT_IT[i][j][1], BT_IT[i][j][2] = -1, -1, -1
		}
	}

	D[0][0] = 0.0
	IS[0][0] = 0.0
	IT[0][0] = 0.0
	for i = 1; i <= 2 * PARA_INFO.Read_len; i++ {
		D[i][0] = float64(math.MaxFloat32)//float64(i) * NEW_INDEL_RATE_LOG
		IS[i][0] = float64(math.MaxFloat32)//float64(i) * NEW_INDEL_RATE_LOG
		IT[i][0] = float64(math.MaxFloat32)//float64(i) * NEW_INDEL_RATE_LOG
	}
	for j = 1; j <= 2 * PARA_INFO.Read_len; j++ {
		D[0][j] = float64(math.MaxFloat32)//0.0
		IS[0][j] = float64(math.MaxFloat32)//PARA_INFO.Gap_open_cost
		IT[0][j] = 0.0//PARA_INFO.Gap_open_cost
		BT_IT[0][j][0], BT_IT[0][j][1] = 2, 2
	}

	var selected_snp_len int
	var prob_i, prob_i_snp, sub_i, mis_i, ins_i_open, ins_i_ext float64
	for i = 1; i <= m; i++ {
		prob_i = -math.Log10(1.0 - math.Pow(10, -(float64(qual[i - 1]) - 33) / 10.0))
		mis_i = PARA_INFO.Sub_cost + prob_i
		ins_i_open = PARA_INFO.Gap_open_cost + prob_i
		ins_i_ext = PARA_INFO.Gap_ext_cost + prob_i
		for j = 1; j <= n; j++ {
			if _, is_snp = INDEX.SNP_PROF[pos + j - 1]; !is_snp {
				if read[i - 1] == ref[j - 1] {
					sub_i = prob_i
				} else {
					sub_i = mis_i
				}
				D[i][j] = D[i - 1][j - 1] + sub_i
				BT_D[i][j][0], BT_D[i][j][1] = 0, 0
				if D[i][j] > IS[i - 1][j - 1] + sub_i {
					D[i][j] = IS[i - 1][j - 1] + sub_i
					BT_D[i][j][0], BT_D[i][j][1] = 0, 1
				}
				if D[i][j] > IT[i - 1][j - 1] + sub_i {
					D[i][j] = IT[i - 1][j - 1] + sub_i
					BT_D[i][j][0], BT_D[i][j][1] = 0, 2
				}

				IS[i][j] = D[i - 1][j] + ins_i_open
				BT_IS[i][j][0], BT_IS[i][j][1] = 1, 0
				if IS[i][j] > IS[i - 1][j] + ins_i_ext {
					IS[i][j] = IS[i - 1][j] + ins_i_ext
					BT_IS[i][j][0], BT_IS[i][j][1] = 1, 1
				}
				
				IT[i][j] = D[i][j - 1] + PARA_INFO.Gap_open_cost
				BT_IT[i][j][0], BT_IT[i][j][1] = 2, 0
				if IT[i][j] > IT[i][j - 1] + PARA_INFO.Gap_ext_cost {
					IT[i][j] = IT[i][j - 1] + PARA_INFO.Gap_ext_cost
					BT_IT[i][j][0], BT_IT[i][j][1] = 2, 2
				}
			} else {
				D[i][j] = float64(math.MaxFloat32)
				IS[i][j] = float64(math.MaxFloat32)
				IT[i][j] = float64(math.MaxFloat32)
				selected_snp_len = 0
				snp_prof, _ = S.SNP_Calls[uint32(pos + j - 1)]
				for snp_str, snp_prob = range snp_prof {
					snp_len = len(snp_str)
					//One possible case: i - snp_len < 0 for all k
					if i - snp_len >= 0 {
						prob_i_snp = AlignCostKnownLoci(read[i - snp_len : i], []byte(snp_str), 
							qual[i - snp_len : i], snp_prob)
						if D[i][j] > D[i - snp_len][j - 1] + prob_i_snp {
							D[i][j] = D[i - snp_len][j - 1] + prob_i_snp
							BT_D[i][j][0], BT_D[i][j][1] = 0, 0
							selected_snp_len = len(snp_str)
						}
						if D[i][j] > IS[i - snp_len][j - 1] + prob_i_snp {
							D[i][j] = IS[i - snp_len][j - 1] + prob_i_snp
							BT_D[i][j][0], BT_D[i][j][1] = 0, 1
							selected_snp_len = len(snp_str)
						}
						if D[i][j] > IT[i - snp_len][j - 1] + prob_i_snp {
							D[i][j] = IT[i - snp_len][j - 1] + prob_i_snp
							BT_D[i][j][0], BT_D[i][j][1] = 0, 2
							selected_snp_len = len(snp_str)
						}
					}
				}
				if selected_snp_len != 0 {
					BT_D[i][j][2] = selected_snp_len
				}
			}
		}
	}
	PrintDisInfo("BwEditDist, D dis", m, n, D[m][n])
	PrintDisInfo("BwEditDist, IS dis", m, n, IS[m][n])
	PrintDisInfo("BwEditDist, IT dis", m, n, IT[m][n])

	PrintEditDisMat("BwEditDist, D mat", D, m, n)
	PrintEditDisMat("BwEditDist, IS mat", IS, m, n)
	PrintEditDisMat("BwEditDist, IT mat", IT, m, n)

	PrintEditTraceMat("BwEditDist, D trace mat", BT_D, m, n)
	PrintEditTraceMat("BwEditDist, IS trace mat", BT_IS, m, n)
	PrintEditTraceMat("BwEditDist, IT trace mat", BT_IT, m, n)

	min_dist := D[m][n]
	bt_mat := 0
	if min_dist > IS[m][n] {
		min_dist = IS[m][n]
		bt_mat = 1
	}
	if min_dist > IT[m][n] {
		min_dist = IT[m][n]
		bt_mat = 2
	}

	return align_prob, min_dist, bt_mat, m, n, snp_pos, snp_base, snp_qual
}

//-------------------------------------------------------------------------------------------------
// BackwardTraceBack constructs alignment between reads and refs from BackwardDistanceMulti.
// 	read is a read.
// 	ref is part of a multi-genome.
// The reads include standard bases, the multi-genomes include standard bases and "*" characters.
//-------------------------------------------------------------------------------------------------
func (S *SNP_Prof) BackwardTraceBack(read, qual, ref []byte, m, n int, pos int, BT_Mat int, 
	BT_D, BT_IS, BT_IT [][][]int) ([]int, [][]byte, [][]byte) {

	var is_snp bool
	var snp_len int
	var snp_pos []int
	var snp_base, snp_qual [][]byte

	PrintEditDisInput("BwEditTraceBack, read, qual, ref", read[ : m], qual[ : m], ref[ : n])

	aligned_read, aligned_qual, aligned_ref := make([]byte, 0), make([]byte, 0), make([]byte, 0)
	bt_mat := BT_Mat
	i, j, k := m, n, 0
	for i > 0 || j > 0 {
		is_snp = false
		if j > 0 {
			_, is_snp = INDEX.SNP_PROF[pos + j - 1]
		}
		if !is_snp { //unknown VARIANT location
			if bt_mat == 0 {
				if read[i - 1] != ref[j - 1] {
					snp_pos = append(snp_pos, pos + j - 1)
					snp_base = append(snp_base, []byte{read[i - 1]})
					snp_qual = append(snp_qual, []byte{qual[i - 1]})
				}
				aligned_read = append(aligned_read, read[i - 1])
				aligned_qual = append(aligned_qual, qual[i - 1])
				aligned_ref = append(aligned_ref, ref[j - 1])
				GetEditTrace("0", i, j, read[i - 1], ref[j - 1])
				bt_mat = BT_D[i][j][1]
				i, j = i - 1, j - 1
			} else if bt_mat == 1 {
				aligned_read = append(aligned_read, read[i - 1])
				aligned_qual = append(aligned_qual, qual[i - 1])
				aligned_ref = append(aligned_ref, '-')
				GetEditTrace("1", i, j, read[i - 1], '-')
				bt_mat = BT_IS[i][j][1]
				i, j = i - 1, j
			} else if bt_mat == 2 {
				aligned_read = append(aligned_read, '-')
				aligned_qual = append(aligned_qual, '-')
				aligned_ref = append(aligned_ref, ref[j - 1])
				GetEditTrace("2", i, j, '-', ref[j - 1])
				bt_mat = BT_IT[i][j][1]
				i, j = i, j - 1
			}
		} else { //known VARIANT location
			if  BT_D[i][j][2] > 0 {
				snp_len = BT_D[i][j][2]
				snp_pos = append(snp_pos, pos + j - 1)
				snp, qlt := make([]byte, snp_len), make([]byte, snp_len)
				copy(snp, read[i - snp_len : i])
				copy(qlt, qual[i - snp_len : i])
				snp_base = append(snp_base, snp)
				snp_qual = append(snp_qual, qlt)
				for k = 0; k < snp_len - 1; k++ {
					aligned_read = append(aligned_read, read[i - 1 - k])
					aligned_qual = append(aligned_qual, qual[i - 1 - k])
					aligned_ref = append(aligned_ref, '+')
				}
				aligned_read = append(aligned_read, read[i - snp_len])
				aligned_qual = append(aligned_qual, qual[i - snp_len])
				aligned_ref = append(aligned_ref, ref[j - 1])
				GetEditTraceKnownLoc("3", i, j, read[i - snp_len : i], ref[j - 1])
				bt_mat = BT_D[i][j][1]
				i, j = i - snp_len, j - 1
			} else {
				aligned_read = append(aligned_read, '-')
				aligned_qual = append(aligned_qual, '-')
				aligned_ref = append(aligned_ref, ref[j - 1])
				GetEditTraceKnownLoc("4", i, j, []byte{'-'}, ref[j - 1])
				bt_mat = BT_IT[i][j][1]
				i, j = i, j - 1
			}
		}
	}

	//Put the alignment in original direction
	for i, j = 0, len(aligned_read) - 1; i < j; i, j = i + 1, j - 1 {
		aligned_read[i], aligned_read[j] = aligned_read[j], aligned_read[i]
		aligned_qual[i], aligned_qual[j] = aligned_qual[j], aligned_qual[i]
		aligned_ref[i],  aligned_ref[j]  = aligned_ref[j],  aligned_ref[i]
	}
	PrintEditAlignInfo("BwEditTraceBack, aligned read/qual/ref", aligned_read, aligned_qual, aligned_ref)

	//Get SNPs
	ref_ori_pos := -1
	for i = 0; i < len(aligned_ref); i++ {
		if aligned_read[i] != '-' {
			break
		}
		ref_ori_pos++
	}
	for i < len(aligned_ref) {
		if aligned_ref[i] != '+' && aligned_ref[i] != '-' {
			ref_ori_pos++
		}
		if aligned_ref[i] == '-' {
			if aligned_read[i - 1] != '-' {
				snp_pos = append(snp_pos, pos + ref_ori_pos)
				snp, qlt := make([]byte, 0), make([]byte, 0)
				snp = append(snp, aligned_read[i - 1])
				qlt = append(qlt, aligned_qual[i - 1])
				for j = i; j < len(aligned_ref) && aligned_ref[j] == '-'; j++ {
					snp = append(snp, aligned_read[j])
					qlt = append(qlt, aligned_qual[j])
				}
				snp_base = append(snp_base, snp)
				snp_qual = append(snp_qual, qlt)
				i = j - 1
			} else {
				for j = i; j < len(aligned_ref) && aligned_ref[j] == '-'; j++ {
				}
				i = j - 1
			}
		}
		//Need to get quality in a proper way in this case!!!
		if aligned_read[i] == '-' {
			snp_pos = append(snp_pos, pos + ref_ori_pos - 1)
			snp, qlt := make([]byte, 0), make([]byte, 0)
			snp = append(snp, aligned_ref[i - 1])
			qlt = append(qlt, aligned_qual[i - 1])
			for j = i; j < len(aligned_read) && aligned_read[j] == '-'; j++ {
				snp = append(snp, aligned_ref[j])
				ref_ori_pos++
			}
			snp_base = append(snp_base, snp)
			snp_qual = append(snp_qual, qlt)
			i = j - 1
		}
		i++
	}
	PrintVarInfo("BwEditTraceBack, variant info", snp_pos, snp_base, snp_qual)
	return snp_pos, snp_base, snp_qual
}

//-------------------------------------------------------------------------------------------------
// Calculate the distance between read and ref in forward direction.
// 	read is a read.
// 	ref is part of a multi-genome.
// The reads include standard bases, the multi-genomes include standard bases and "*" characters.
//-------------------------------------------------------------------------------------------------
func (S *SNP_Prof) ForwardDistance(read, qual, ref []byte, pos int, D [][]float64, BT [][][]int) (float64, float64, int, int, []int, [][]byte, [][]byte) {

	var snp_len int
	var snp_prof map[string]float64
	var is_snp, is_same_len_snp bool
	var snp_str string
	var p, min_p, snp_prob float64
	var snp_pos []int
	var snp_base, snp_qual [][]byte

	PrintEditDisInput("fw dis input: read, qual, ref", read, qual, ref)
	align_prob := 0.0
	M, N := len(read), len(ref)
	m, n := M, N

	for m > 0 && n > 0 {
		if _, is_snp = INDEX.SNP_PROF[pos + N - n]; !is_snp {
			if read[M - m] != ref[N - n] {
                snp_pos = append(snp_pos, pos + N - n)
                snp_base = append(snp_base, []byte{read[M - m]})
                snp_qual = append(snp_qual, []byte{qual[M - m]})
				align_prob = align_prob - math.Log10(NEW_SNP_RATE) - math.Log10(1.0 - math.Pow(10, -(float64(qual[M - m]) - 33) / 10.0))
			}
			m--
			n--
		} else if snp_len, is_same_len_snp = INDEX.SAME_LEN_SNP[pos + N - n]; is_same_len_snp {
			min_p = math.MaxFloat64
			snp_prof, is_snp = S.SNP_Calls[uint32(pos + N - n)]
			for snp_str, snp_prob = range snp_prof {
				if m >= snp_len {
					p = AlignCostKnownLoci(read[M - m : M - m + snp_len], []byte(snp_str), qual[M - m : M - m + snp_len], snp_prob)
					if min_p > p {
						min_p = p
					}
				}
			}
			if min_p < math.MaxFloat64 {
				align_prob = align_prob + min_p
				snp_pos = append(snp_pos, pos + N - n)
				snp, qlt := make([]byte, snp_len), make([]byte, snp_len)
				copy(snp, read[M - m : M - (m - snp_len)])
				copy(qlt, qual[M - m : M - (m - snp_len)])
				snp_base = append(snp_base, snp)
				snp_qual = append(snp_qual, qlt)
				m -= snp_len
				n--
			} else {
				break
			}
		} else {
			break
		}
		if align_prob > PARA_INFO.Prob_thres {
			return PARA_INFO.Prob_thres + 1, 0, m, n, snp_pos, snp_base, snp_qual
		}
	}

	PrintDisInfo("fw H dis", m, n, align_prob)

	var i, j int

	/*
	Backtrace matrix, for each BT[i][j]:
	BT[i][j][0]: direction, can be 0: d (diagonal arrow, back to i-1,j-1), 1: u (up arrow, back to i-1,j), 2: l (left arrow, back to i,j-1)
	BT[i][j][1]: number of shift (equal to length of called variant) at known variant loc, can be any integer number, for example 5 means back to i-5,j-1
	*/
	for i := 0; i <= 2 * PARA_INFO.Read_len; i++ {
		for j := 0; j <= 2 * PARA_INFO.Read_len; j++ {
			BT[i][j][0], BT[i][j][1] = 0, 0
		}
	}

	D[0][0] = 0.0
	for i = 1; i <= 2 * PARA_INFO.Read_len; i++ {
		D[i][0] = float64(i) * NEW_SNP_RATE_LOG
		BT[i][0][0] = 1
	}
	for j = 1; j <= 2 * PARA_INFO.Read_len; j++ {
		D[0][j] = 0.0
		BT[0][j][0] = 2
	}

	var temp_p, sub_i, id_i float64
	var min_snp string
	for i = 1; i <= m; i++ {
		sub_i = NEW_SNP_RATE_LOG - math.Log10(1.0 - math.Pow(10, -(float64(qual[M - i]) - 33) / 10.0))
		id_i =  NEW_INDEL_RATE_LOG - math.Log10(1.0 - math.Pow(10, -(float64(qual[M - i]) - 33) / 10.0))
		for j = 1; j <= n; j++ {
			if _, is_snp = INDEX.SNP_PROF[pos + N - j]; !is_snp {
				if read[M - i] != ref[N - j] {
					D[i][j] = D[i - 1][j - 1] + sub_i
					BT[i][j][0] = 0
					if D[i][j] > D[i - 1][j] + id_i {
						D[i][j] = D[i - 1][j] + id_i
						BT[i][j][0] = 1
					}
					if D[i][j] > D[i][j - 1] + NEW_INDEL_RATE_LOG {
						D[i][j] = D[i][j - 1] + NEW_INDEL_RATE_LOG
						BT[i][j][0] = 2
					}
				} else {
					D[i][j] = D[i - 1][j - 1]
					BT[i][j][0] = 0
				}
			} else {
				D[i][j] = float64(math.MaxFloat32)
				min_snp = ""
				snp_prof, _ = S.SNP_Calls[uint32(pos + N - j)]
				for snp_str, snp_prob = range snp_prof {
					snp_len = len(snp_str)
					//One possible case: i - snp_len < 0 for all k
					if i - snp_len >= 0 {
						temp_p = D[i - snp_len][j - 1] + AlignCostKnownLoci(read[M - i : M - i + snp_len], []byte(snp_str), qual[M - i : M - i + snp_len], snp_prob)
						if D[i][j] > temp_p {
							D[i][j] = temp_p
							min_snp = snp_str
						}
					}
				}
				if min_snp != "" {
					BT[i][j][0] = 0
					BT[i][j][1] = len(min_snp)
				}
			}
		}
	}
	PrintDisInfo("fw E dis", m, n, D[m][n])
	PrintEditDisMat("fw edit dis mat", D, m, n)
	PrintEditTraceMat("fw edit trace mat", BT, m, n)

	return align_prob, D[m][n], m, n, snp_pos, snp_base, snp_qual
}

//-------------------------------------------------------------------------------------------------
// ForwardTraceBack constructs alignment based on the results from ForwardDistanceMulti.
// 	read is a read.
// 	ref is part of a multi-genome.
// The reads include standard bases, the multi-genomes include standard bases and "*" characters.
//-------------------------------------------------------------------------------------------------
func (S *SNP_Prof) ForwardTraceBack(read, qual, ref []byte, m, n int, pos, r_pos int, BT [][][]int) ([]int, [][]byte, [][]byte) {

	var is_snp bool
	var snp_len int
	var snp_pos []int
	var snp_base, snp_qual [][]byte

	M, N := len(read), len(ref)
	PrintEditDisInput("fw E, read, qual, ref", read[M - m :], qual[M - m :], ref[N - n :])

	var aligned_read, aligned_qual, aligned_ref []byte
	i, j := m, n
	for i > 0 || j > 0 {
		_, is_snp = INDEX.SNP_PROF[pos + N - j]
		if !is_snp { //unknown VARIANT location
			if BT[i][j][0] == 0 {
				if read[M - i] != ref[N - j] {
					snp_pos = append(snp_pos, pos + N - j)
					snp_base = append(snp_base, []byte{read[M - i]})
					snp_qual = append(snp_qual, []byte{qual[M - i]})
				}
				aligned_read = append(aligned_read, read[M - i])
				aligned_qual = append(aligned_qual, qual[M - i])
				aligned_ref = append(aligned_ref, ref[N - j])
				GetEditTrace("1", i, j, read[M - i], ref[N - j])
				i, j = i - 1, j - 1
			} else if BT[i][j][0] == 1 {
				aligned_read = append(aligned_read, read[M - i])
				aligned_qual = append(aligned_qual, qual[M - i])
				aligned_ref = append(aligned_ref, '-')
				GetEditTrace("2", i, j, read[M - i], '-')
				i, j = i - 1, j
			} else {
				aligned_read = append(aligned_read, '-')
				aligned_qual = append(aligned_qual, '-')
				aligned_ref = append(aligned_ref, ref[N - j])
				GetEditTrace("3", i, j, '-', ref[N - j])
				i, j = i, j - 1
			}
		} else { //known VARIANT location
			if BT[i][j][1] > 0 {
				snp_len = BT[i][j][1]
				snp_pos = append(snp_pos, pos + N - j)
				snp, qlt := make([]byte, snp_len), make([]byte, snp_len)
				copy(snp, read[M - i : M - (i - snp_len)])
				copy(qlt, qual[M - i : M - (i - snp_len)])
				snp_base = append(snp_base, snp)
				snp_qual = append(snp_qual, qlt)
				aligned_read = append(aligned_read, read[M - i])
				aligned_qual = append(aligned_qual, qual[M - i])
				aligned_ref = append(aligned_ref, ref[N - j])
				for k := 1; k < snp_len; k++ {
					aligned_read = append(aligned_read, read[M - i + k])
					aligned_qual = append(aligned_qual, qual[M - i + k])
					aligned_ref = append(aligned_ref, '+')
				}
				GetEditTraceKnownLoc("5", i, j, read[M - i : M - i + snp_len], ref[N - j])
				i, j = i - snp_len, j - 1
			} else {
				aligned_read = append(aligned_read, '-')
				aligned_qual = append(aligned_qual, '-')
				aligned_ref = append(aligned_ref, ref[N - j])
				GetEditTrace("6", i, j, '-', ref[N - j])
				i, j = i, j - 1
			}
		}
	}
	PrintEditDisInput("fw E, read, qual, ref, after backtrace", read[M - m :], qual[M - m :], ref[N - n :])
	PrintEditAlignInfo("fw aligned read/qual/ref E", aligned_read, aligned_qual, aligned_ref)
	ref_pos := 1
	for i = 1; i < len(aligned_ref); i++ {
		if aligned_ref[i] == '-' {
			snp_pos = append(snp_pos, pos + N - n + ref_pos - 1)
			snp, qlt := make([]byte, 0), make([]byte, 0)
			snp = append(snp, aligned_read[i - 1])
			qlt = append(qlt, aligned_qual[i - 1])
			for j = i; j < len(aligned_ref) && aligned_ref[j] == '-'; j++ {
				snp = append(snp, aligned_read[j])
				qlt = append(qlt, aligned_qual[j])
			}
			snp_base = append(snp_base, snp)
			snp_qual = append(snp_qual, qlt)
			i = j - 1
		}
		ref_pos++
	}
	return snp_pos, snp_base, snp_qual
}
