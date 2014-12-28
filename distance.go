//-------------------------------------------------------------------------------------------------
// Multigenome package: distance module.
// Calculating distance and determining alignment between reads and "starred" multigenomes.
// Copyright 2014 Nam Sy Vo
//-------------------------------------------------------------------------------------------------

package isc

//-------------------------------------------------------------------------------------------------
// Cost functions for computing distance between reads and multi-genomes.
// Input slices should have same length
//-------------------------------------------------------------------------------------------------
func Cost(read, ref []byte) int {
	cost := 0
	for i := 0; i < len(read); i++ {
		if read[i] != ref[i] {
			cost += 1
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
func (I *Index) BackwardDistance(read, ref []byte, pos int, D [][]int, T [][][]byte) (int, int, int, int, []int, [][]byte, []int) {

	var i, j, k int

	var cost int
	var d, min_d, m, n int
	var snp_len int
	var snp_values [][]byte
	var is_snp, is_same_len_snp bool

	d = 0
	m, n = len(read), len(ref)
	var snp_pos []int
	var snp_idx []int
	var snp_val [][]byte
	for m > 0 && n > 0 {
		snp_values, is_snp = I.SNP_PROF[pos + n - 1]
		snp_len, is_same_len_snp = I.SAME_LEN_SNP[pos + n - 1]
		if !is_snp {
			if read[m - 1] != ref[n - 1] {
				snp_pos = append(snp_pos, pos + n - 1)
				snp_idx = append(snp_idx, m - 1)
				snp := read[m - 1]
				snp_val = append(snp_val, []byte{snp})
				d++
			}
			m--
			n--
		} else if is_same_len_snp {
			min_d = 1000 * INF // 1000*INF is a value for testing, will change to a better solution later
			for i = 0; i < len(snp_values); i++ {
				cost = Cost(read[m - snp_len : m], snp_values[i])
				if min_d > cost {
					min_d = cost
				}
			}
			if min_d >= INF {
				return INF, 0, m, n, snp_pos, snp_val, snp_idx
			}
			snp_pos = append(snp_pos, pos + n - 1)
			snp_idx = append(snp_idx, m - snp_len)
			snp := make([]byte, snp_len)
			copy(snp, read[m - snp_len : m])
			snp_val = append(snp_val, snp)
			d += min_d
			m -= snp_len
			n--
		} else {
			break
		}
		if d > PARA_INFO.Dist_thres {
			return PARA_INFO.Dist_thres + 1, 0, m, n, snp_pos, snp_val, snp_idx
		}
	}

	D[0][0] = 0
	for i = 1; i <= PARA_INFO.Read_len; i++ {
		D[i][0] = INF
	}
	for j = 1; j <= PARA_INFO.Read_len; j++ {
		D[0][j] = 0
	}
	var temp_dis, min_index int
	for i = 1; i <= m; i++ {
		for j = 1; j <= n; j++ {
			snp_values, is_snp = I.SNP_PROF[pos + j - 1]
			if !is_snp {
				if read[i - 1] != ref[j - 1] {
					D[i][j] = D[i - 1][j - 1] + 1
				} else {
					D[i][j] = D[i - 1][j - 1]
				}
			} else {
				D[i][j] = 1000 * INF //1000*INF is a value for testing, will change to a better solution later
				min_index = 0
				for k = 0; k < len(snp_values); k++ {
					snp_len = len(snp_values[k])
					//One possnble case: i - snp_len < 0 for all k
					if i-snp_len >= 0 {
						if snp_values[k][0] != '.' {
							temp_dis = D[i - snp_len][j - 1] + Cost(read[i - snp_len : i], snp_values[k])
						} else {
							temp_dis = D[i][j - 1]
						}
						if D[i][j] > temp_dis {
							D[i][j] = temp_dis
							min_index = k
						}
					}
				}
				T[i - 1][j - 1] = snp_values[min_index]
			}
		}
	}
	if D[m][n] >= INF {
		return d, INF, m, n, snp_pos, snp_val, snp_idx
	}
	return d, D[m][n], m, n, snp_pos, snp_val, snp_idx
}

//-------------------------------------------------------------------------------------------------
// BackwardTraceBack constructs alignment between reads and refs from BackwardDistanceMulti.
// 	read is a read.
// 	ref is part of a multi-genome.
// The reads include standard bases, the multi-genomes include standard bases and "*" characters.
//-------------------------------------------------------------------------------------------------
func (I Index) BackwardTraceBack(read, ref []byte, m, n int, pos int, T [][][]byte) ([]int, [][]byte, []int) {

	var is_snp bool
	var snp_len int
	var i, j int = m, n
	var snp_pos []int
	var snp_idx []int
	var snp_val [][]byte
	for i > 0 || j > 0 {
		_, is_snp = I.SNP_PROF[pos + j - 1]
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
func (I *Index) ForwardDistance(read, ref []byte, pos int, D [][]int, T [][][]byte) (int, int, int, int, []int, [][]byte, []int) {

	var i, j, k int
	var M, N = len(read), len(ref)

	var cost int
	var d, min_d, m, n int
	var snp_len int
	var snp_values [][]byte
	var is_snp, is_same_len_snp bool

	d = 0
	m, n = M, N
	var snp_pos []int
	var snp_idx []int
	var snp_val [][]byte
	for m > 0 && n > 0 {
		snp_values, is_snp = I.SNP_PROF[pos + N - n]
		snp_len, is_same_len_snp = I.SAME_LEN_SNP[pos + N - n]
		if !is_snp {
			if read[M - m] != ref[N - n] {
                snp_pos = append(snp_pos, pos + N - n)
                snp_idx = append(snp_idx, M - m)
                snp := read[M - m]
                snp_val = append(snp_val, []byte{snp})
				d++
			}
			m--
			n--
		} else if is_same_len_snp {
			min_d = 1000 * INF //1000*INF is a value for testing, will change to a better solution later
			for i = 0; i < len(snp_values); i++ {
				cost = Cost(read[M - m : M - (m - snp_len)], snp_values[i])
				if min_d > cost {
					min_d = cost
				}
			}
			if min_d >= INF {
				return INF, 0, m, n, snp_pos, snp_val, snp_idx
			}
			snp_pos = append(snp_pos, pos + N - n)
			snp_idx = append(snp_idx, M - m)
			snp := make([]byte, snp_len)
			copy(snp, read[M - m : M - (m - snp_len)])
			snp_val = append(snp_val, snp)
			d += min_d
			m -= snp_len
			n--
		} else {
			break
		}
		if d > PARA_INFO.Dist_thres {
			return PARA_INFO.Dist_thres + 1, 0, m, n, snp_pos, snp_val, snp_idx
		}
	}

	D[0][0] = 0
	for i = 1; i <= PARA_INFO.Read_len; i++ {
		D[i][0] = INF
	}
	for i = 1; i <= PARA_INFO.Read_len; i++ {
		D[0][i] = 0
	}
	var temp_dis, min_index int
	for i = 1; i <= m; i++ {
		for j = 1; j <= n; j++ {
			snp_values, is_snp = I.SNP_PROF[pos + N - j]
			if !is_snp {
				if read[M - i] != ref[N - j] {
					D[i][j] = D[i - 1][j - 1] + 1
				} else {
					D[i][j] = D[i - 1][j - 1]
				}
			} else {
				D[i][j] = 1000 * INF //1000*INF is a value for testing, will change to a better solution later
				min_index = 0
				for k = 0; k < len(snp_values); k++ {
					snp_len = len(snp_values[k])
					//One possnble case: i - snp_len < 0 for all k
					if i - snp_len >= 0 {
						if snp_values[k][0] != '.' {
							temp_dis = D[i - snp_len][j - 1] + Cost(read[M - i : M - (i - snp_len)], snp_values[k])
						} else {
							temp_dis = D[i][j - 1]
						}
						if D[i][j] > temp_dis {
							D[i][j] = temp_dis
							min_index = k
						}
					}
				}
				T[i - 1][j - 1] = snp_values[min_index]
			}
		}
	}
	if D[m][n] >= INF {
		return d, INF, m, n, snp_pos, snp_val, snp_idx
	}
	return d, D[m][n], m, n, snp_pos, snp_val, snp_idx
}

//-------------------------------------------------------------------------------------------------
// ForwardTraceBack constructs alignment based on the results from ForwardDistanceMulti.
// 	read is a read.
// 	ref is part of a multi-genome.
// The reads include standard bases, the multi-genomes include standard bases and "*" characters.
//-------------------------------------------------------------------------------------------------
func (I *Index) ForwardTraceBack(read, ref []byte, m, n int, pos int, T [][][]byte) ([]int, [][]byte, []int) {
	var is_snp bool
	var snp_len int
	var i, j int = m, n
	var M, N int = len(read), len(ref)

	//Trace back from right-bottom conner of distance matrix T
	var snp_pos []int
	var snp_idx []int
	var snp_val [][]byte
	for i > 0 || j > 0 {
		_, is_snp = I.SNP_PROF[pos + N - j]
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
