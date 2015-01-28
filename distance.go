//-------------------------------------------------------------------------------------------------
// Multigenome package: distance module.
// Calculating distance and determining alignment between reads and "starred" multigenomes.
// Copyright 2014 Nam Sy Vo
//-------------------------------------------------------------------------------------------------

package isc

import "math"

//-------------------------------------------------------------------------------------------------
// Cost functions for computing distance between reads and multi-genomes.
// Input slices should have same length
//-------------------------------------------------------------------------------------------------
func AlignProb(read, ref, qual []byte, prob float64) float64 {
	p := 0.0
	for i := 0; i < len(read); i++ {
		p = p - math.Log10(prob) - math.Log10(1.0 - math.Pow(10, -(float64(qual[i]) - 33) / 10.0))
	}
	return p
}

//-------------------------------------------------------------------------------------------------
// Calculate the distance between read and ref in backward direction.
// 	read is a read.
// 	ref is part of a multi-genome.
// The reads include standard bases, the multi-genome includes standard bases and "*" characters.
//-------------------------------------------------------------------------------------------------
func (S *SNP_Prof) BackwardDistance(read, qual, ref []byte, pos int, D [][]float64, T [][][]byte) (float64, float64, int, int, []int, [][]byte, []int) {

	var snp_len int
	var snp_str string
	var is_snp, is_same_len_snp bool
	var p, min_p, snp_prob float64
	var snp_prof map[string]float64

	align_prob := 0.0
	m, n := len(read), len(ref)
	var snp_pos []int
	var snp_idx []int
	var snp_val [][]byte
	for m > 0 && n > 0 {
		if _, is_snp = INDEX.SNP_PROF[pos + n - 1]; !is_snp {
			if read[m - 1] != ref[n - 1] {
				snp_pos = append(snp_pos, pos + n - 1)
				snp_idx = append(snp_idx, m - 1)
				snp := read[m - 1]
				snp_val = append(snp_val, []byte{snp})
				align_prob = align_prob - math.Log10(EPSILON) - math.Log10(1.0 - math.Pow(10, -(float64(qual[m - 1]) - 33) / 10.0))
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
	var i, j int
	D[0][0] = 0.0
	for i = 1; i <= 2 * PARA_INFO.Read_len; i++ {
		D[i][0] = math.MaxFloat64
	}
	for j = 1; j <= 2 * PARA_INFO.Read_len; j++ {
		D[0][j] = 0.0
	}
	var temp_p float64
	var min_snp string
	for i = 1; i <= m; i++ {
		for j = 1; j <= n; j++ {
			if _, is_snp = INDEX.SNP_PROF[pos + j - 1]; !is_snp {
				if read[i - 1] != ref[j - 1] {
					D[i][j] = D[i - 1][j - 1] - math.Log10(EPSILON) - math.Log10(1.0 - math.Pow(10, -(float64(qual[i - 1]) - 33) / 10.0))
				} else {
					D[i][j] = D[i - 1][j - 1]
				}
			} else {
				D[i][j] = math.MaxFloat64
				min_snp = ""
				snp_prof, is_snp = S.SNP_Calls[uint32(pos + j - 1)]
				for snp_str, snp_prob = range snp_prof {
					snp_len = len(snp_str)
					//One possible case: i - snp_len < 0 for all k
					if i - snp_len >= 0 {
						if snp_str != "." {
							temp_p = D[i - snp_len][j - 1] + AlignProb(read[i - snp_len : i], []byte(snp_str), qual[i - snp_len : i], snp_prob)
						} else {
							temp_p = D[i][j - 1]
						}
						if D[i][j] > temp_p {
							D[i][j] = temp_p
							min_snp = snp_str
						}
					}
				}
				T[i - 1][j - 1] = []byte(min_snp)
			}
		}
	}
	return align_prob, D[m][n], m, n, snp_pos, snp_val, snp_idx
}

//-------------------------------------------------------------------------------------------------
// BackwardTraceBack constructs alignment between reads and refs from BackwardDistanceMulti.
// 	read is a read.
// 	ref is part of a multi-genome.
// The reads include standard bases, the multi-genomes include standard bases and "*" characters.
//-------------------------------------------------------------------------------------------------
func (S *SNP_Prof) BackwardTraceBack(read, qual, ref []byte, m, n int, pos int, T [][][]byte) ([]int, [][]byte, []int) {

	var is_snp bool
	var snp_len int
	var snp_pos, snp_idx []int
	var snp_val [][]byte

	i, j := m, n
	for i > 0 || j > 0 {
		_, is_snp = INDEX.SNP_PROF[pos + j - 1]
		if i > 0 && j > 0 {
			if !is_snp {
				if read[i - 1] != ref[j - 1] {
					snp_pos = append(snp_pos, pos + j - 1)
					snp_idx = append(snp_idx, i - 1)
					snp := read[i - 1]
					snp_val = append(snp_val, []byte{snp})
				}
				i, j = i - 1, j - 1
			} else {
				if T[i - 1][j - 1][0] != '.' {
					snp_len = len(T[i - 1][j - 1])
				} else {
					snp_len = 0
				}
				snp_pos = append(snp_pos, pos + j - 1)
				snp_idx = append(snp_idx, i - snp_len)
				snp := make([]byte, snp_len)
				copy(snp, read[i - snp_len : i])
				snp_val = append(snp_val, snp)
				i, j = i - snp_len, j - 1
			}
		} else if i == 0 {
			j = j - 1
		} else if j == 0 {
			i = i - 1
		}
	}
	return snp_pos, snp_val, snp_idx
}

//-------------------------------------------------------------------------------------------------
// Calculate the distance between read and ref in forward direction.
// 	read is a read.
// 	ref is part of a multi-genome.
// The reads include standard bases, the multi-genomes include standard bases and "*" characters.
//-------------------------------------------------------------------------------------------------
func (S *SNP_Prof) ForwardDistance(read, qual, ref []byte, pos int, D [][]float64, T [][][]byte) (float64, float64, int, int, []int, [][]byte, []int) {

	var snp_len int
	var snp_prof map[string]float64
	var is_snp, is_same_len_snp bool
	var min_snp, snp_str string
	var p, min_p, snp_prob float64
	var snp_pos, snp_idx []int
	var snp_val [][]byte

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
				align_prob = align_prob - math.Log10(EPSILON) - math.Log10(1.0 - math.Pow(10, -(float64(qual[M - m]) - 33) / 10.0))
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

	var i, j int
	D[0][0] = 0.0
	for i = 1; i <= 2 * PARA_INFO.Read_len; i++ {
		D[i][0] = float64(math.MaxFloat32)
	}
	for j = 1; j <= 2 * PARA_INFO.Read_len; j++ {
		D[0][j] = 0.0
	}
	var temp_p float64
	for i = 1; i <= m; i++ {
		for j = 1; j <= n; j++ {
			if _, is_snp = INDEX.SNP_PROF[pos + N - j]; !is_snp {
				if read[M - i] != ref[N - j] {
					D[i][j] = D[i - 1][j - 1] - math.Log10(EPSILON) - math.Log10(1.0 - math.Pow(10, -(float64(qual[M - i]) - 33) / 10.0))
				} else {
					D[i][j] = D[i - 1][j - 1]
				}
			} else {
				D[i][j] = math.MaxFloat64
				min_snp = ""
				snp_prof, is_snp = S.SNP_Calls[uint32(pos + N - j)]
				for snp_str, snp_prob = range snp_prof {
					snp_len = len(snp_str)
					//One possnble case: i - snp_len < 0 for all k
					if i - snp_len >= 0 {
						if snp_str != "." {
							temp_p = D[i - snp_len][j - 1] + AlignProb(read[M - i : M - (i - snp_len)], []byte(snp_str), qual[M - i : M - (i - snp_len)], snp_prob)
						} else {
							temp_p = D[i][j - 1]
						}
						if D[i][j] > temp_p {
							D[i][j] = temp_p
							min_snp = snp_str
						}
					}
				}
				T[i - 1][j - 1] = []byte(min_snp)
			}
		}
	}
	return align_prob, D[m][n], m, n, snp_pos, snp_val, snp_idx
}

//-------------------------------------------------------------------------------------------------
// ForwardTraceBack constructs alignment based on the results from ForwardDistanceMulti.
// 	read is a read.
// 	ref is part of a multi-genome.
// The reads include standard bases, the multi-genomes include standard bases and "*" characters.
//-------------------------------------------------------------------------------------------------
func (S *SNP_Prof) ForwardTraceBack(read, qual, ref []byte, m, n int, pos int, T [][][]byte) ([]int, [][]byte, []int) {
	var is_snp bool
	var snp_len int
	var snp_pos, snp_idx []int
	var snp_val [][]byte

	M, N := len(read), len(ref)
	i, j := m, n
	for i > 0 || j > 0 {
		_, is_snp = INDEX.SNP_PROF[pos + N - j]
		if i > 0 && j > 0 {
			if !is_snp {
				if read[M - i] != ref[N - j] {
					snp_pos = append(snp_pos, pos + N - j)
					snp_idx = append(snp_idx, M - i)
					snp := read[M - i]
					snp_val = append(snp_val, []byte{snp})
				}
				i, j = i - 1, j - 1
			} else {
				if T[i - 1][j - 1][0] != '.' {
					snp_len = len(T[i - 1][j - 1])
				} else {
					snp_len = 0
				}
				snp_pos = append(snp_pos, pos + N - j)
				snp_idx = append(snp_idx, M - i)
				snp := make([]byte, snp_len)
				copy(snp, read[M - i : M - (i - snp_len)])
				snp_val = append(snp_val, snp)
				i, j = i - snp_len, j - 1
			}
		} else if i == 0 {
			j = j - 1
		} else if j == 0 {
			i = i - 1
		}
	}
	return snp_pos, snp_val, snp_idx
}
