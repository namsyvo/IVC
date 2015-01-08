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
func Cost(read, ref, qual []byte) float64 {
	cost := 0.0
	for i := 0; i < len(read); i++ {
		if read[i] != ref[i] {
			cost = cost + math.Log10(0.01) + math.Log10(1.0 - math.Pow(10, -(float64(qual[i]) - 33) / 10.0))
		} else {
			cost = cost + math.Log10(0.97) + math.Log10(1.0 - math.Pow(10, -(float64(qual[i]) - 33) / 10.0))
		}
	}
	return cost
}

//-------------------------------------------------------------------------------------------------
// Calculate the distance between read and ref in backward direction.
// 	read is a read.
// 	ref is part of a multi-genome.
// The reads include standard bases, the multi-genome includes standard bases and "*" characters.
//-------------------------------------------------------------------------------------------------
func (S *SNP_Prof) BackwardDistance(read, qual, ref []byte, pos int, D [][]float64, T [][][]byte) (float64, float64, int, int, []int, [][]byte, []int) {

	var i, j int

	var cost float64
	var m, n int
	var snp_len int
	var snp_str string
	var min_d, snp_prob float64
	var snp_prof map[string]float64
	var is_snp, is_same_len_snp bool

	p := 0.0
	m, n = len(read), len(ref)
	var snp_pos []int
	var snp_idx []int
	var snp_val [][]byte
	for m > 0 && n > 0 {
		snp_prof, is_snp = S.SNP_Calls[uint32(pos + n - 1)] //INDEX.SNP_PROF[pos + n - 1]
		snp_len, is_same_len_snp = INDEX.SAME_LEN_SNP[pos + n - 1]
		if !is_snp {
			if read[m - 1] != ref[n - 1] {
				snp_pos = append(snp_pos, pos + n - 1)
				snp_idx = append(snp_idx, m - 1)
				snp := read[m - 1]
				snp_val = append(snp_val, []byte{snp})
				p = p + math.Log10(EPSILON) + math.Log10(1.0 - math.Pow(10, -(float64(qual[m-1]) - 33) / 10.0))
			} else {
				p = p + math.Log10(1 - 3 * EPSILON) + math.Log10(1.0 - math.Pow(10, -(float64(qual[m-1]) - 33) / 10.0))
			}
			m--
			n--
		} else if is_same_len_snp {
			min_d = float64(math.MaxFloat32)
			for snp_str, snp_prob = range snp_prof {
				cost = snp_prob + Cost(read[m - snp_len : m], []byte(snp_str), qual[m - snp_len : m])
				if min_d > cost {
					min_d = cost
				}
			}
			p = p + cost
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
		//if d > PARA_INFO.Dist_thres {
		//	return PARA_INFO.Dist_thres + 1, 0, m, n, snp_pos, snp_val, snp_idx, p
		//}
	}

	D[0][0] = 0.0
	for i = 1; i <= PARA_INFO.Read_len; i++ {
		D[i][0] = 0.0
	}
	for j = 1; j <= PARA_INFO.Read_len; j++ {
		D[0][j] = -float64(math.MaxFloat32)
	}
	var temp_dis float64
	var min_snp string
	for i = 1; i <= m; i++ {
		for j = 1; j <= n; j++ {
			snp_prof, is_snp = S.SNP_Calls[uint32(pos + j - 1)]
			if !is_snp {
				if read[i - 1] != ref[j - 1] {
					D[i][j] = D[i - 1][j - 1] + math.Log10(0.01) + math.Log10(1.0 - math.Pow(10, -(float64(qual[i - 1]) - 33) / 10.0))
				} else {
					D[i][j] = D[i - 1][j - 1] + math.Log10(0.97) + math.Log10(1.0 - math.Pow(10, -(float64(qual[i - 1]) - 33) / 10.0))
				}
			} else {
				D[i][j] = 0.0
				min_snp = ""
				for snp_str, snp_prob = range snp_prof {
					snp_len = len(snp_str)
					//One possnble case: i - snp_len < 0 for all k
					if i - snp_len >= 0 {
						if snp_str[0] != '.' {
							temp_dis = D[i - snp_len][j - 1] + math.Log10(snp_prob) - Cost(read[i - snp_len : i], []byte(snp_str), qual[i - snp_len : i])
						} else {
							temp_dis = D[i][j - 1]
						}
						if D[i][j] > temp_dis {
							D[i][j] = temp_dis
							min_snp = snp_str
						}
					}
				}
				T[i - 1][j - 1] = []byte(min_snp)
			}
		}
	}
	//if D[m][n] >= INF {
	//	return d, INF, m, n, snp_pos, snp_val, snp_idx
	//}
	return p, D[m][n], m, n, snp_pos, snp_val, snp_idx
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
	var i, j int = m, n
	var snp_pos []int
	var snp_idx []int
	var snp_val [][]byte

	for i > 0 || j > 0 {
		_, is_snp = S.SNP_Calls[uint32(pos + j - 1)]
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

	var i, j int
	var M, N = len(read), len(ref)

	var cost float64
	var m, n int
	var snp_len int
	var snp_prof map[string]float64
	var is_snp, is_same_len_snp bool
	var min_snp, snp_str string
	var min_d, snp_prob float64

	m, n = M, N
	var snp_pos []int
	var snp_idx []int
	var snp_val [][]byte
	p := 0.0
	for m > 0 && n > 0 {
		snp_prof, is_snp = S.SNP_Calls[uint32(pos + N - n)]
		snp_len, is_same_len_snp = INDEX.SAME_LEN_SNP[pos + N - n]
		if !is_snp {
			if read[M - m] != ref[N - n] {
                snp_pos = append(snp_pos, pos + N - n)
                snp_idx = append(snp_idx, M - m)
                snp := read[M - m]
                snp_val = append(snp_val, []byte{snp})
				p = p + math.Log10(0.01) + math.Log10(1.0 - math.Pow(10, -(float64(qual[M - m]) - 33) / 10.0))
			} else {
				p = p + math.Log10(0.97) + math.Log10(1.0 - math.Pow(10, -(float64(qual[M - m]) - 33) / 10.0))
			}
			m--
			n--
		} else if is_same_len_snp {
			min_d = float64(math.MaxFloat32)
			for snp_str, snp_prob = range snp_prof {
				cost = math.Log10(snp_prob) + Cost(read[M - m : M - (m - snp_len)], []byte(snp_str), qual[M - m : M - (m - snp_len)])
				if min_d > cost {
					min_d = cost
				}
			}
			p = p + cost
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
		//if d > PARA_INFO.Dist_thres {
		//	return PARA_INFO.Dist_thres + 1, 0, m, n, snp_pos, snp_val, snp_idx, p
		//}
	}

	D[0][0] = 0.0
	for i = 1; i <= PARA_INFO.Read_len; i++ {
		D[i][0] = 0.0
	}
	for i = 1; i <= PARA_INFO.Read_len; i++ {
		D[0][i] = -float64(math.MaxFloat32)
	}
	var temp_dis float64
	for i = 1; i <= m; i++ {
		for j = 1; j <= n; j++ {
			snp_prof, is_snp = S.SNP_Calls[uint32(pos + N - j)]
			if !is_snp {
				if read[M - i] != ref[N - j] {
					D[i][j] = D[i - 1][j - 1] + math.Log10(0.01) - math.Log10(1.0 - math.Pow(10, -(float64(qual[M - i]) - 33) / 10.0))
				} else {
					D[i][j] = D[i - 1][j - 1] + math.Log10(0.97) - math.Log10(1.0 - math.Pow(10, -(float64(qual[M - i]) - 33) / 10.0))
				}
			} else {
				D[i][j] = 0.0
				min_snp = ""
				for snp_str, snp_prob = range snp_prof {
					snp_len = len(snp_str)
					//One possnble case: i - snp_len < 0 for all k
					if i - snp_len >= 0 {
						if snp_str[0] != '.' {
							temp_dis = D[i - snp_len][j - 1] + Cost(read[M - i : M - (i - snp_len)], []byte(snp_str), qual[M - i : M - (i - snp_len)])
						} else {
							temp_dis = D[i][j - 1]
						}
						if D[i][j] > temp_dis {
							D[i][j] = temp_dis
							min_snp = snp_str
						}
					}
				}
				T[i - 1][j - 1] = []byte(min_snp)
			}
		}
	}
	//if D[m][n] >= INF {
	//	return d, INF, m, n, snp_pos, snp_val, snp_idx
	//}
	return p, D[m][n], m, n, snp_pos, snp_val, snp_idx
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
	var i, j int = m, n
	var M, N int = len(read), len(ref)

	//Trace back from right-bottom conner of distance matrix T
	var snp_pos []int
	var snp_idx []int
	var snp_val [][]byte

	for i > 0 || j > 0 {
		_, is_snp = S.SNP_Calls[uint32(pos + j - 1)]
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
