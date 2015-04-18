//-------------------------------------------------------------------------------------------------
// Multigenome package: distance module.
// Calculating distance and determining alignment between reads and "starred" multigenomes.
// Copyright 2014 Nam Sy Vo
//-------------------------------------------------------------------------------------------------

package isc

import (
	"math"
)

//-------------------------------------------------------------------------------------------------
// Cost functions for computing distance between reads and multi-genomes.
// Input slices should have same length
//-------------------------------------------------------------------------------------------------
func AlignProb(read, ref, qual []byte, prob float64) float64 {
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
func (S *SNP_Prof) BackwardDistance(read, qual, ref []byte, pos int, D [][]float64, BT [][][]int) (float64, float64, int, int, []int, [][]byte, []int) {

	var snp_len int
	var snp_str string
	var is_snp, is_same_len_snp bool
	var p, min_p, snp_prob float64
	var snp_prof map[string]float64

	align_prob := 0.0
	m, n := len(read), len(ref)

	PrintEditDisInput("bw dis input: read, qual, ref", read, qual, ref)
	var snp_pos []int
 	var snp_idx []int
	var snp_val [][]byte
	var sub_m float64
	for m > 0 && n > 0 {
		sub_m = NEW_SNP_RATE_LOG - math.Log10(1.0 - math.Pow(10, -(float64(qual[m - 1]) - 33) / 10.0))
		if _, is_snp = INDEX.SNP_PROF[pos + n - 1]; !is_snp {
			if read[m - 1] != ref[n - 1] {
				snp_pos = append(snp_pos, pos + n - 1)
				snp_idx = append(snp_idx, m - 1)
				snp := read[m - 1]
				snp_val = append(snp_val, []byte{snp})
				align_prob = align_prob + sub_m
			}
			m--
			n--
		} else if snp_len, is_same_len_snp = INDEX.SAME_LEN_SNP[pos + n - 1]; is_same_len_snp {
			snp_prof, _ = S.SNP_Calls[uint32(pos + n - 1)]
			min_p = math.MaxFloat64
			for snp_str, snp_prob = range snp_prof {
				if m >= snp_len {
					p = AlignProb(read[m - snp_len : m], []byte(snp_str), qual[m - snp_len : m], snp_prob)
					if min_p > p {
						min_p = p
					}
				}
			}
			if min_p < math.MaxFloat64 {
				align_prob = align_prob + min_p
				snp_pos = append(snp_pos, pos + n - 1)
				snp_idx = append(snp_idx, m - snp_len)
				snp := make([]byte, snp_len)
				copy(snp, read[m - snp_len : m])
				snp_val = append(snp_val, snp)
				m -= snp_len
				n--
			} else {
				break
			}
		} else {
			break
		}
		if align_prob > PARA_INFO.Prob_thres {
			return PARA_INFO.Prob_thres + 1, 0, m, n, snp_pos, snp_val, snp_idx
		}
	}

	PrintEditDisInfo("bw H dis", m, n, align_prob)

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
		sub_i = NEW_SNP_RATE_LOG - math.Log10(1.0 - math.Pow(10, -(float64(qual[i - 1]) - 33) / 10.0))
		id_i =  NEW_INDEL_RATE_LOG - math.Log10(1.0 - math.Pow(10, -(float64(qual[i - 1]) - 33) / 10.0))
		for j = 1; j <= n; j++ {
			if _, is_snp = INDEX.SNP_PROF[pos + j - 1]; !is_snp {
				if read[i - 1] != ref[j - 1] {
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
				snp_prof, _ = S.SNP_Calls[uint32(pos + j - 1)]
				for snp_str, snp_prob = range snp_prof {
					snp_len = len(snp_str)
					//One possible case: i - snp_len < 0 for all k
					if i - snp_len >= 0 {
						temp_p = D[i - snp_len][j - 1] + AlignProb(read[i - snp_len : i], []byte(snp_str), qual[i - snp_len : i], snp_prob)
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
	PrintEditDisInfo("bw E dis", m, n, D[m][n])
	PrintEditDisMat("bw edit dis mat", D, m, n)
	PrintEditTraceMat("bw edit trace mat", BT, m, n)

	return align_prob, D[m][n], m, n, snp_pos, snp_val, snp_idx
}

//-------------------------------------------------------------------------------------------------
// BackwardTraceBack constructs alignment between reads and refs from BackwardDistanceMulti.
// 	read is a read.
// 	ref is part of a multi-genome.
// The reads include standard bases, the multi-genomes include standard bases and "*" characters.
//-------------------------------------------------------------------------------------------------
func (S *SNP_Prof) BackwardTraceBack(read, qual, ref []byte, m, n int, pos int, BT [][][]int) ([]int, [][]byte, []int) {

	var is_snp bool
	var snp_len int
	var snp_pos, snp_idx []int
	var snp_val [][]byte

	PrintEditDisInput("bw E, read, qual, ref", read[ : m], qual[ : m], ref[ : n])

	var aligned_read, aligned_ref []byte
	i, j := m, n
	for i > 0 || j > 0 {
		_, is_snp = INDEX.SNP_PROF[pos + j - 1]
		if !is_snp { //unknown VARIANT location
			if BT[i][j][0] == 0 {
				if read[i - 1] != ref[j - 1] {
					snp_pos = append(snp_pos, pos + j - 1)
					snp_idx = append(snp_idx, i - 1)
					snp := read[i - 1]
					snp_val = append(snp_val, []byte{snp})
				}
				aligned_read = append(aligned_read, read[i - 1])
				aligned_ref = append(aligned_ref, ref[j - 1])
				PrintEditTraceStep("1", i, j, read[i - 1], ref[j - 1])
				i, j = i - 1, j - 1
			} else if BT[i][j][0] == 1 {
				aligned_read = append(aligned_read, read[i - 1])
				aligned_ref = append(aligned_ref, '-')
				PrintEditTraceStep("2", i, j, read[i - 1], '-')
				i, j = i - 1, j
			} else {
				aligned_read = append(aligned_read, '-')
				aligned_ref = append(aligned_ref, ref[j - 1])
				PrintEditTraceStep("3", i, j, '-', ref[j - 1])
				i, j = i, j - 1
			}
		} else { //known VARIANT location
			snp_len = BT[i][j][1]
			snp_pos = append(snp_pos, pos + j - 1)
			snp_idx = append(snp_idx, i - snp_len)
			snp := make([]byte, snp_len)
			copy(snp, read[i - snp_len : i])
			snp_val = append(snp_val, snp)
			for k := 0; k < snp_len - 1; k++ {
				aligned_ref = append(aligned_ref, '+')
				aligned_read = append(aligned_read, read[i - 1 - k])
			}
			aligned_ref = append(aligned_ref, ref[j - 1])
			aligned_read = append(aligned_read, read[i - snp_len])
			PrintEditTraceStepKnownLoc("5", i, j, read[i - snp_len : i], ref[j - 1])
			i, j = i - snp_len, j - 1
		}
	}
	PrintEditDisInput("bw E, read, qual, ref, after backtrace", read[ : m], qual[ : m], ref[ : n])
    for i, j := 0, len(aligned_read) - 1; i < j; i, j = i + 1, j - 1 {
        aligned_read[i], aligned_read[j] = aligned_read[j], aligned_read[i]
        aligned_ref[i], aligned_ref[j] = aligned_ref[j], aligned_ref[i]
    }
	PrintEditAlignInfo("bw aligned read/ref E", aligned_read, aligned_ref)
	ref_pos := 1
	for i = 1; i < len(aligned_ref); i++ {
		if aligned_ref[i] == '-' {
			snp_pos = append(snp_pos, pos + ref_pos - 1)
			snp_idx = append(snp_idx, ref_pos - 1)
			snp := make([]byte, 0)
			snp = append(snp, aligned_read[ref_pos - 1])
			for j = i; j < len(aligned_ref) && aligned_ref[j] == '-'; j++ {
				snp = append(snp, aligned_read[j])
			}
			snp_val = append(snp_val, snp)
			i = j - 1
		}
		ref_pos++
	}
	return snp_pos, snp_val, snp_idx
}

//-------------------------------------------------------------------------------------------------
// Calculate the distance between read and ref in forward direction.
// 	read is a read.
// 	ref is part of a multi-genome.
// The reads include standard bases, the multi-genomes include standard bases and "*" characters.
//-------------------------------------------------------------------------------------------------
func (S *SNP_Prof) ForwardDistance(read, qual, ref []byte, pos int, D [][]float64, BT [][][]int) (float64, float64, int, int, []int, [][]byte, []int) {

	var snp_len int
	var snp_prof map[string]float64
	var is_snp, is_same_len_snp bool
	var snp_str string
	var p, min_p, snp_prob float64
	var snp_pos, snp_idx []int
	var snp_val [][]byte

	PrintEditDisInput("fw dis input: read, qual, ref", read, qual, ref)
	align_prob := 0.0
	M, N := len(read), len(ref)
	m, n := M, N

	for m > 0 && n > 0 {
		if _, is_snp = INDEX.SNP_PROF[pos + N - n]; !is_snp {
			if read[M - m] != ref[N - n] {
                snp_pos = append(snp_pos, pos + N - n)
                snp_idx = append(snp_idx, M - m)
                snp := read[M - m]
                snp_val = append(snp_val, []byte{snp})
				align_prob = align_prob - math.Log10(NEW_SNP_RATE) - math.Log10(1.0 - math.Pow(10, -(float64(qual[M - m]) - 33) / 10.0))
			}
			m--
			n--
		} else if snp_len, is_same_len_snp = INDEX.SAME_LEN_SNP[pos + N - n]; is_same_len_snp {
			min_p = math.MaxFloat64
			snp_prof, is_snp = S.SNP_Calls[uint32(pos + N - n)]
			for snp_str, snp_prob = range snp_prof {
				if m >= snp_len {
					p = AlignProb(read[M - m : M - (m - snp_len)], []byte(snp_str), qual[M - m : M - (m - snp_len)], snp_prob)
					if min_p > p {
						min_p = p
					}
				}
			}
			if min_p < math.MaxFloat64 {
				align_prob = align_prob + min_p
				snp_pos = append(snp_pos, pos + N - n)
				snp_idx = append(snp_idx, M - m)
				snp := make([]byte, snp_len)
				copy(snp, read[M - m : M - (m - snp_len)])
				snp_val = append(snp_val, snp)
				m -= snp_len
				n--
			} else {
				break
			}
		} else {
			break
		}
		if align_prob > PARA_INFO.Prob_thres {
			return PARA_INFO.Prob_thres + 1, 0, m, n, snp_pos, snp_val, snp_idx
		}
	}

	PrintEditDisInfo("fw H dis", m, n, align_prob)

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
						temp_p = D[i - snp_len][j - 1] + AlignProb(read[M - i : M - (i - snp_len)], []byte(snp_str), qual[M - i : M - (i - snp_len)], snp_prob)
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
	PrintEditDisInfo("fw E dis", m, n, D[m][n])
	PrintEditDisMat("fw edit dis mat", D, m, n)
	PrintEditTraceMat("fw edit trace mat", BT, m, n)

	return align_prob, D[m][n], m, n, snp_pos, snp_val, snp_idx
}

//-------------------------------------------------------------------------------------------------
// ForwardTraceBack constructs alignment based on the results from ForwardDistanceMulti.
// 	read is a read.
// 	ref is part of a multi-genome.
// The reads include standard bases, the multi-genomes include standard bases and "*" characters.
//-------------------------------------------------------------------------------------------------
func (S *SNP_Prof) ForwardTraceBack(read, qual, ref []byte, m, n int, pos int, BT [][][]int) ([]int, [][]byte, []int) {

	var is_snp bool
	var snp_len int
	var snp_pos, snp_idx []int
	var snp_val [][]byte

	M, N := len(read), len(ref)
	PrintEditDisInput("fw E, read, qual, ref", read[M - m :], qual[M - m :], ref[N - n :])

	var aligned_read, aligned_ref []byte
	i, j := m, n
	for i > 0 || j > 0 {
		_, is_snp = INDEX.SNP_PROF[pos + N - j]
		if !is_snp { //unknown VARIANT location
			if BT[i][j][0] == 0 {
				if read[M - i] != ref[N - j] {
					snp_pos = append(snp_pos, pos + N - j)
					snp_idx = append(snp_idx, M - i)
					snp := read[M - i]
					snp_val = append(snp_val, []byte{snp})
				}
				aligned_read = append(aligned_read, read[M - i])
				aligned_ref = append(aligned_ref, ref[N - j])
				PrintEditTraceStep("1", i, j, read[M - i], ref[N - j])
				i, j = i - 1, j - 1
			} else if BT[i][j][0] == 1 {
				aligned_read = append(aligned_read, read[M - i])
				aligned_ref = append(aligned_ref, '-')
				PrintEditTraceStep("2", i, j, read[M - i], '-')
				i, j = i - 1, j
			} else {
				aligned_read = append(aligned_read, '-')
				aligned_ref = append(aligned_ref, ref[N - j])
				PrintEditTraceStep("3", i, j, '-', ref[N - j])
				i, j = i, j - 1
			}
		} else { //known VARIANT location
			snp_len = BT[i][j][1]
			snp_pos = append(snp_pos, pos + N - j)
			snp_idx = append(snp_idx, M - i)
			snp := make([]byte, snp_len)
			copy(snp, read[M - i : M - (i - snp_len)])
			snp_val = append(snp_val, snp)
			for k := 0; k < snp_len - 1; k++ {
				aligned_ref = append(aligned_ref, '+')
				aligned_read = append(aligned_read, read[M - (i - snp_len) - 1 - k])
			}
			aligned_ref = append(aligned_ref, ref[N - j])
			aligned_read = append(aligned_read, read[M - i])
			PrintEditTraceStepKnownLoc("5", i, j, read[M - i : M - (i - snp_len)], ref[N - j])
			i, j = i - snp_len, j - 1
		}
	}
	PrintEditDisInput("fw E, read, qual, ref, after backtrace", read[M - m :], qual[M - m :], ref[N - n :])
	PrintEditAlignInfo("fw aligned read/ref E", aligned_read, aligned_ref)
	ref_pos := 1
	for i = 1; i < len(aligned_ref); i++ {
		if aligned_ref[i] == '-' {
			snp_pos = append(snp_pos, pos + N - n + ref_pos - 1)
			snp_idx = append(snp_idx, M - m + ref_pos - 1)
			snp := make([]byte, 0)
			snp = append(snp, aligned_read[ref_pos - 1])
			for j = i; j < len(aligned_ref) && aligned_ref[j] == '-'; j++ {
				snp = append(snp, aligned_read[j])
			}
			snp_val = append(snp_val, snp)
			i = j - 1
		}
		ref_pos++
	}
	return snp_pos, snp_val, snp_idx
}
