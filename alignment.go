//-------------------------------------------------------------------------------------------------
// ISC: alignment.go - Calculating alignment between reads and multigenomes.
// Alignment is performed in both backward and forward directions.
// Copyright 2015 Nam Sy Vo.
//-------------------------------------------------------------------------------------------------

package isc

import (
	//"fmt"
	"math"
)

//-------------------------------------------------------------------------------------------------
// AlignCostKnownLoci calculates cost of alignment between a read and the reference at known loci.
//-------------------------------------------------------------------------------------------------
func AlignCostKnownLoci(read, ref, qual []byte, prob float64) float64 {
	p := 0.0
	for i := 0; i < len(read); i++ {
		if read[i] != ref[i] {
			return math.MaxFloat64
		} else {
			p = p - math.Log10(1.0-math.Pow(10, -(float64(qual[i])-33)/10.0))
		}
	}
	return p - math.Log10(prob)
}

//-------------------------------------------------------------------------------------------------
// BackwardDistance calculates the distance between a read and a ref in backward direction.
// The read include standard bases, the ref includes standard bases and "*" characters.
//-------------------------------------------------------------------------------------------------
func (VC *VarCall) BackwardDistance(read, qual, ref []byte, pos int, D, IS, IT [][]float64,
	BT_D, BT_IS, BT_IT [][][]int, ref_pos_map []int) (float64, float64, int, int, int, []int, [][]byte, [][]byte, []int) {

	var var_len int
	var var_str string
	var is_var, is_same_len_var bool
	var p, min_p, var_prob float64
	var var_prof map[string]float64

	align_prob := 0.0
	m, n := len(read), len(ref)

	PrintEditDisInput("bw align input: read, qual, ref", read, qual, ref)
	var var_pos, var_type []int
	var var_base, var_qual [][]byte
	for m > 0 && n > 0 {
		if _, is_var = INDEX.VarProf[ref_pos_map[n-1]-PARA_INFO.Indel_backup]; is_var {
			if _, is_same_len_var = INDEX.SameLenVar[ref_pos_map[n-1]-PARA_INFO.Indel_backup]; !is_same_len_var {
				break
			}
		}
		if _, is_var = INDEX.VarProf[ref_pos_map[n-1]]; !is_var {
			if read[m-1] != ref[n-1] {
				if m+PARA_INFO.Ham_backup < len(read) && n+PARA_INFO.Ham_backup < len(ref) {
					m += PARA_INFO.Ham_backup
					n += PARA_INFO.Ham_backup
				}
				break
				/*
					var_pos = append(var_pos, ref_pos_map[n - 1])
					var_base = append(var_base, []byte{read[m - 1]})
					var_qual = append(var_qual, []byte{qual[m - 1]})
					var_type = append(var_type, 0)
				*/
				align_prob += PARA_INFO.Sub_cost - math.Log10(1.0-math.Pow(10, -(float64(qual[m-1])-33)/10.0))
			}
			m--
			n--
		} else if var_len, is_same_len_var = INDEX.SameLenVar[ref_pos_map[n-1]]; is_same_len_var {
			var_prof, _ = VC.VarProb[uint32(ref_pos_map[n-1])]
			min_p = math.MaxFloat64
			for var_str, var_prob = range var_prof {
				if m >= var_len {
					p = AlignCostKnownLoci(read[m-var_len:m], []byte(var_str), qual[m-var_len:m], var_prob)
					if min_p > p {
						min_p = p
					}
				}
			}
			if min_p < math.MaxFloat64 {
				align_prob = align_prob + min_p
				var_pos = append(var_pos, ref_pos_map[n-1])
				v, q := make([]byte, var_len), make([]byte, var_len)
				copy(v, read[m-var_len:m])
				copy(q, qual[m-var_len:m])
				var_base = append(var_base, v)
				var_qual = append(var_qual, q)
				var_type = append(var_type, 0)
				m -= var_len
				n--
			} else {
				break
			}
		} else {
			break
		}
		if align_prob > PARA_INFO.Prob_thres {
			return PARA_INFO.Prob_thres + 1, 0, -1, m, n, var_pos, var_base, var_qual, var_type
		}
	}
	PrintDisInfo("bw Ham dis", m, n, align_prob)

	if m == 0 || n == 0 {
		return align_prob, 0, -1, m, n, var_pos, var_base, var_qual, var_type
	}

	PrintEditDisInput("bw align for Edit: read, qual, ref", read[:m], qual[:m], ref[:n])

	/*
		Backtrace matrix, for each BT[i][j]:
		BT[i][j][0]: direction, can be 0: diagonal arrow (back to i-1,j-1), 1: up arrow (back to i-1,j),
		 	2: left arrow (back to i,j-1).
		BT[i][j][1]: matrix, can be 0: matrix for D, 1: matrix for IS, 2: matrix for IT.
		BT[i][j][2]: number of shift (equal to length of called variant) at known variant loc,
			can be any integer number, e.g. 5 means back to i-5,j-1.
	*/
	var i, j int
	for i := 0; i <= 2*PARA_INFO.Read_len; i++ {
		for j := 0; j <= 2*PARA_INFO.Read_len; j++ {
			BT_D[i][j][0], BT_D[i][j][1], BT_D[i][j][2] = -1, -1, -1
			BT_IS[i][j][0], BT_IS[i][j][1], BT_IS[i][j][2] = -1, -1, -1
			BT_IT[i][j][0], BT_IT[i][j][1], BT_IT[i][j][2] = -1, -1, -1
		}
	}

	D[0][0] = 0.0
	IS[0][0] = float64(math.MaxFloat32)
	IT[0][0] = float64(math.MaxFloat32)
	IS[1][0] = PARA_INFO.Gap_open_cost
	BT_IS[1][0][0], BT_IS[1][0][1] = 1, 1

	for i = 1; i <= 2*PARA_INFO.Read_len; i++ {
		D[i][0] = float64(math.MaxFloat32)
		IT[i][0] = float64(math.MaxFloat32)
	}
	for i = 2; i <= 2*PARA_INFO.Read_len; i++ {
		IS[i][0] = PARA_INFO.Gap_ext_cost
		BT_IS[i][0][0], BT_IS[i][0][1] = 1, 1
	}

	for j = 1; j <= 2*PARA_INFO.Read_len; j++ {
		D[0][j] = float64(math.MaxFloat32)
		IS[0][j] = float64(math.MaxFloat32)
		IT[0][j] = 0.0
		BT_IT[0][j][0], BT_IT[0][j][1] = 2, 2
	}

	var selected_var_len int
	var prob_i, prob_i_var, sub_i, mis_i, ins_i_open, ins_i_ext float64
	for i = 1; i <= m; i++ {
		prob_i = -math.Log10(1.0 - math.Pow(10, -(float64(qual[i-1])-33)/10.0))
		mis_i = PARA_INFO.Sub_cost + prob_i
		ins_i_open = PARA_INFO.Gap_open_cost // + prob_i
		ins_i_ext = PARA_INFO.Gap_ext_cost   // + prob_i
		for j = 1; j <= n; j++ {
			if _, is_var = INDEX.VarProf[ref_pos_map[j-1]]; !is_var {
				if read[i-1] == ref[j-1] {
					sub_i = 0 //prob_i
				} else {
					sub_i = mis_i
				}
				D[i][j] = D[i-1][j-1] + sub_i
				BT_D[i][j][0], BT_D[i][j][1] = 0, 0
				if D[i][j] > IS[i-1][j-1]+sub_i {
					D[i][j] = IS[i-1][j-1] + sub_i
					BT_D[i][j][0], BT_D[i][j][1] = 0, 1
				}
				if D[i][j] > IT[i-1][j-1]+sub_i {
					D[i][j] = IT[i-1][j-1] + sub_i
					BT_D[i][j][0], BT_D[i][j][1] = 0, 2
				}

				IS[i][j] = D[i-1][j] + ins_i_open
				BT_IS[i][j][0], BT_IS[i][j][1] = 1, 0
				if IS[i][j] > IS[i-1][j]+ins_i_ext {
					IS[i][j] = IS[i-1][j] + ins_i_ext
					BT_IS[i][j][0], BT_IS[i][j][1] = 1, 1
				}

				IT[i][j] = D[i][j-1] + PARA_INFO.Gap_open_cost
				BT_IT[i][j][0], BT_IT[i][j][1] = 2, 0
				if IT[i][j] > IT[i][j-1]+PARA_INFO.Gap_ext_cost {
					IT[i][j] = IT[i][j-1] + PARA_INFO.Gap_ext_cost
					BT_IT[i][j][0], BT_IT[i][j][1] = 2, 2
				}
			} else {
				D[i][j] = float64(math.MaxFloat32)
				IS[i][j] = float64(math.MaxFloat32)
				IT[i][j] = float64(math.MaxFloat32)
				selected_var_len = 0
				var_prof, _ = VC.VarProb[uint32(ref_pos_map[j-1])]
				for var_str, var_prob = range var_prof {
					var_len = len(var_str)
					//One possible case: i - var_len < 0 for all k
					if i-var_len >= 0 {
						prob_i_var = AlignCostKnownLoci(read[i-var_len:i], []byte(var_str),
							qual[i-var_len:i], var_prob)
						if D[i][j] > D[i-var_len][j-1]+prob_i_var {
							D[i][j] = D[i-var_len][j-1] + prob_i_var
							BT_D[i][j][0], BT_D[i][j][1] = 0, 0
							selected_var_len = len(var_str)
						}
						if D[i][j] > IS[i-var_len][j-1]+prob_i_var {
							D[i][j] = IS[i-var_len][j-1] + prob_i_var
							BT_D[i][j][0], BT_D[i][j][1] = 0, 1
							selected_var_len = len(var_str)
						}
						if D[i][j] > IT[i-var_len][j-1]+prob_i_var {
							D[i][j] = IT[i-var_len][j-1] + prob_i_var
							BT_D[i][j][0], BT_D[i][j][1] = 0, 2
							selected_var_len = len(var_str)
						}
					}
				}
				if selected_var_len != 0 {
					BT_D[i][j][2] = selected_var_len
				}
			}
		}
	}
	PrintDisInfo("BwEditDist, D dis", m, n, D[m][n])
	PrintDisInfo("BwEditDist, IS dis", m, n, IS[m][n])
	PrintDisInfo("BwEditDist, IT dis", m, n, IT[m][n])

	PrintEditDisMat("BwEditDist, D mat", D, m, n, read[:m], ref[:n])
	PrintEditDisMat("BwEditDist, IS mat", IS, m, n, read[:m], ref[:n])
	PrintEditDisMat("BwEditDist, IT mat", IT, m, n, read[:m], ref[:n])

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

	return align_prob, min_dist, bt_mat, m, n, var_pos, var_base, var_qual, var_type
}

//-------------------------------------------------------------------------------------------------
// BackwardTraceBack constructs alignment between a read and a ref from BackwardDistance.
// The read includes standard bases, the ref include standard bases and "*" characters.
//-------------------------------------------------------------------------------------------------
func (VC *VarCall) BackwardTraceBack(read, qual, ref []byte, m, n int, pos int, BT_Mat int,
	BT_D, BT_IS, BT_IT [][][]int, ref_pos_map []int) ([]int, [][]byte, [][]byte, []int) {

	var is_var, is_same_len_var, is_del bool
	var var_len int
	var var_pos, var_type []int
	var var_base, var_qual [][]byte

	PrintEditDisInput("BwEditTraceBack, read, qual, ref", read[:m], qual[:m], ref[:n])

	aligned_read, aligned_qual, aligned_ref := make([]byte, 0), make([]byte, 0), make([]byte, 0)
	bt_mat := BT_Mat
	i, j, k := m, n, 0
	for i > 0 || j > 0 {
		is_var = false
		if j > 0 {
			_, is_var = INDEX.VarProf[ref_pos_map[j-1]]
		}
		if !is_var { //unknown VARIANT location
			if bt_mat == 0 {
				if read[i-1] != ref[j-1] {
					var_pos = append(var_pos, ref_pos_map[j-1])
					var_base = append(var_base, []byte{read[i-1]})
					var_qual = append(var_qual, []byte{qual[i-1]})
					var_type = append(var_type, 0)
				}
				aligned_read = append(aligned_read, read[i-1])
				aligned_qual = append(aligned_qual, qual[i-1])
				aligned_ref = append(aligned_ref, ref[j-1])
				GetEditTrace("0", i, j, read[i-1], ref[j-1])
				bt_mat = BT_D[i][j][1]
				i, j = i-1, j-1
			} else if bt_mat == 1 {
				aligned_read = append(aligned_read, read[i-1])
				aligned_qual = append(aligned_qual, qual[i-1])
				aligned_ref = append(aligned_ref, '-')
				GetEditTrace("1", i, j, read[i-1], '-')
				bt_mat = BT_IS[i][j][1]
				i, j = i-1, j
			} else if bt_mat == 2 {
				aligned_read = append(aligned_read, '-')
				aligned_qual = append(aligned_qual, '-')
				aligned_ref = append(aligned_ref, ref[j-1])
				GetEditTrace("2", i, j, '-', ref[j-1])
				bt_mat = BT_IT[i][j][1]
				i, j = i, j-1
			}
		} else { //known VARIANT location
			if BT_D[i][j][2] > 0 {
				var_len = BT_D[i][j][2]
				var_pos = append(var_pos, ref_pos_map[j-1])
				v, q := make([]byte, var_len), make([]byte, var_len)
				copy(v, read[i-var_len:i])
				copy(q, qual[i-var_len:i])
				var_base = append(var_base, v)
				var_qual = append(var_qual, q)
				if _, is_del = INDEX.DelVar[ref_pos_map[j-1]]; is_del {
					var_type = append(var_type, 2)
				} else if _, is_same_len_var = INDEX.SameLenVar[ref_pos_map[j-1]]; is_same_len_var {var_type = append(var_type, 0)
				} else {
					var_type = append(var_type, 1)
				}
				for k = 0; k < var_len-1; k++ {
					aligned_read = append(aligned_read, read[i-1-k])
					aligned_qual = append(aligned_qual, qual[i-1-k])
					aligned_ref = append(aligned_ref, '+')
				}
				aligned_read = append(aligned_read, read[i-var_len])
				aligned_qual = append(aligned_qual, qual[i-var_len])
				aligned_ref = append(aligned_ref, ref[j-1])
				GetEditTraceKnownLoc("3", i, j, read[i-var_len:i], ref[j-1])
				bt_mat = BT_D[i][j][1]
				i, j = i-var_len, j-1
			} else {
				aligned_read = append(aligned_read, '-')
				aligned_qual = append(aligned_qual, '-')
				aligned_ref = append(aligned_ref, ref[j-1])
				GetEditTraceKnownLoc("4", i, j, []byte{'-'}, ref[j-1])
				bt_mat = BT_IT[i][j][1]
				i, j = i, j-1
			}
		}
	}

	//Put the alignment in original direction
	for i, j = 0, len(aligned_read)-1; i < j; i, j = i+1, j-1 {
		aligned_read[i], aligned_read[j] = aligned_read[j], aligned_read[i]
		aligned_qual[i], aligned_qual[j] = aligned_qual[j], aligned_qual[i]
		aligned_ref[i], aligned_ref[j] = aligned_ref[j], aligned_ref[i]
	}
	PrintEditAlignInfo("BwEditTraceBack, aligned read/qual/ref", aligned_read, aligned_qual, aligned_ref)

	//Get Vars
	ref_ori_pos := 0
	i = 0
	for i < len(aligned_ref) {
		if aligned_read[i] == '-' && aligned_ref[i] != '-' {
			ref_ori_pos++
			i++
		} else if aligned_read[i] != '-' && aligned_ref[i] == '-' {
			i++
		} else {
			break
		}
	}
	for i < len(aligned_ref) {
		if aligned_read[i] != '-' && aligned_ref[i] == '-' { //Insertions
			v, q := make([]byte, 0), make([]byte, 0)
			v = append(v, aligned_read[i-1])
			q = append(q, aligned_qual[i-1])
			for j = i; j < len(aligned_ref) && aligned_ref[j] == '-'; j++ {
				v = append(v, aligned_read[j])
				q = append(q, aligned_qual[j])
			}
			if j < len(aligned_ref) {
				var_pos = append(var_pos, ref_pos_map[ref_ori_pos-1])
				var_base = append(var_base, v)
				var_qual = append(var_qual, q)
				var_type = append(var_type, 1)
			}
			i = j
		} else if aligned_read[i] == '-' && aligned_ref[i] != '-' { //Deletions
			v, q := make([]byte, 0), make([]byte, 0)
			v = append(v, aligned_ref[i-1])
			q = append(q, aligned_qual[i-1]) //A temporary solution, need to get quality in a proper way in this case!!!
			for j = i; j < len(aligned_read) && aligned_read[j] == '-'; j++ {
				v = append(v, aligned_ref[j])
			}
			if j < len(aligned_read) {
				var_pos = append(var_pos, ref_pos_map[ref_ori_pos-1])
				var_base = append(var_base, v)
				var_qual = append(var_qual, q)
				var_type = append(var_type, 2)
			}
			ref_ori_pos += j - i
			i = j
		} else if aligned_ref[i] == '+' {
			i++
		} else {
			ref_ori_pos++
			i++
		}
	}
	PrintVarInfo("BwEditTraceBack, variant info", var_pos, var_base, var_qual)
	return var_pos, var_base, var_qual, var_type
}

//-------------------------------------------------------------------------------------------------
// ForwardDistance calculates the distance between a read and a ref in forward direction.
// The read includes standard bases, the ref includes standard bases and "*" characters.
//-------------------------------------------------------------------------------------------------
func (VC *VarCall) ForwardDistance(read, qual, ref []byte, pos int, D, IS, IT [][]float64,
	BT_D, BT_IS, BT_IT [][][]int, ref_pos_map []int) (float64, float64, int, int, int, []int, [][]byte, [][]byte, []int) {

	var var_len int
	var var_prof map[string]float64
	var is_var, is_same_len_var bool
	var var_str string
	var p, min_p, var_prob float64
	var var_pos, var_type []int
	var var_base, var_qual [][]byte

	PrintEditDisInput("fw dis input: read, qual, ref", read, qual, ref)
	align_prob := 0.0
	M, N := len(read), len(ref)
	m, n := M, N
	for m > 0 && n > 0 {
		if _, is_var = INDEX.VarProf[ref_pos_map[N-n]+PARA_INFO.Indel_backup]; is_var {
			if _, is_same_len_var = INDEX.SameLenVar[ref_pos_map[N-n]+PARA_INFO.Indel_backup]; !is_same_len_var {
				break
			}
		}
		if _, is_var = INDEX.VarProf[ref_pos_map[N-n]]; !is_var {
			if read[M-m] != ref[N-n] {
				if M-(m+PARA_INFO.Ham_backup) > 0 && N-(n+PARA_INFO.Ham_backup) > 0 {
					m += PARA_INFO.Ham_backup
					n += PARA_INFO.Ham_backup
				}
				break
				/*
				   var_pos = append(var_pos, ref_pos_map[N - n])
				   var_base = append(var_base, []byte{read[M - m]})
				   var_qual = append(var_qual, []byte{qual[M - m]})
				   var_type = append(var_type, 0)
				*/
				align_prob += PARA_INFO.Sub_cost - math.Log10(1.0-math.Pow(10, -(float64(qual[M-m])-33)/10.0))
			}
			m--
			n--
		} else if var_len, is_same_len_var = INDEX.SameLenVar[ref_pos_map[N-n]]; is_same_len_var {
			min_p = math.MaxFloat64
			var_prof, is_var = VC.VarProb[uint32(ref_pos_map[N-n])]
			for var_str, var_prob = range var_prof {
				if m >= var_len {
					p = AlignCostKnownLoci(read[M-m:M-m+var_len], []byte(var_str), qual[M-m:M-m+var_len], var_prob)
					if min_p > p {
						min_p = p
					}
				}
			}
			if min_p < math.MaxFloat64 {
				align_prob = align_prob + min_p
				var_pos = append(var_pos, ref_pos_map[N-n])
				v, q := make([]byte, var_len), make([]byte, var_len)
				copy(v, read[M-m:M-(m-var_len)])
				copy(q, qual[M-m:M-(m-var_len)])
				var_base = append(var_base, v)
				var_qual = append(var_qual, q)
				var_type = append(var_type, 0)
				m -= var_len
				n--
			} else {
				break
			}
		} else {
			break
		}
		if align_prob > PARA_INFO.Prob_thres {
			return PARA_INFO.Prob_thres + 1, 0, -1, m, n, var_pos, var_base, var_qual, var_type
		}
	}

	PrintDisInfo("fw Ham dis", m, n, align_prob)

	if m == 0 || n == 0 {
		return align_prob, 0, -1, m, n, var_pos, var_base, var_qual, var_type
	}

	PrintEditDisInput("fw align for Edit: read, qual, ref", read[M-m:M], qual[M-m:M], ref[N-n:N])

	/*
		Backtrace matrix, for each BT[i][j]:
		BT[i][j][0]: direction, can be 0: diagonal arrow (back to i-1,j-1), 1: up arrow (back to i-1,j),
		 	2: left arrow (back to i,j-1).
		BT[i][j][1]: matrix, can be 0: matrix for D, 1: matrix for IS, 2: matrix for IT.
		BT[i][j][2]: number of shift (equal to length of called variant) at known variant loc,
			can be any integer number, e.g. 5 means back to i-5,j-1.
	*/
	var i, j int
	for i := 0; i <= 2*PARA_INFO.Read_len; i++ {
		for j := 0; j <= 2*PARA_INFO.Read_len; j++ {
			BT_D[i][j][0], BT_D[i][j][1], BT_D[i][j][2] = -1, -1, -1
			BT_IS[i][j][0], BT_IS[i][j][1], BT_IS[i][j][2] = -1, -1, -1
			BT_IT[i][j][0], BT_IT[i][j][1], BT_IT[i][j][2] = -1, -1, -1
		}
	}

	D[0][0] = 0.0
	for i = 1; i <= 2*PARA_INFO.Read_len; i++ {
		D[i][0] = float64(math.MaxFloat32)
		IT[i][0] = float64(math.MaxFloat32)
	}
	IS[0][0] = float64(math.MaxFloat32)
	IS[1][0] = PARA_INFO.Gap_open_cost
	BT_IS[1][0][0], BT_IS[1][0][1] = 1, 1
	for i = 2; i <= 2*PARA_INFO.Read_len; i++ {
		IS[i][0] = PARA_INFO.Gap_ext_cost
		BT_IS[i][0][0], BT_IS[i][0][1] = 1, 1
	}

	IT[0][0] = float64(math.MaxFloat32)
	for j = 1; j <= 2*PARA_INFO.Read_len; j++ {
		D[0][j] = float64(math.MaxFloat32)
		IS[0][j] = float64(math.MaxFloat32)
		IT[0][j] = 0.0
		BT_IT[0][j][0], BT_IT[0][j][1] = 2, 2
	}

	var selected_var_len int
	var prob_i, prob_i_var, sub_i, mis_i, ins_i_open, ins_i_ext float64
	for i = 1; i <= m; i++ {
		prob_i = -math.Log10(1.0 - math.Pow(10, -(float64(qual[M-i])-33)/10.0))
		mis_i = PARA_INFO.Sub_cost + prob_i
		ins_i_open = PARA_INFO.Gap_open_cost // + prob_i
		ins_i_ext = PARA_INFO.Gap_ext_cost   // + prob_i
		for j = 1; j <= n; j++ {
			if _, is_var = INDEX.VarProf[ref_pos_map[N-j]]; !is_var {
				if read[M-i] == ref[N-j] {
					sub_i = 0 //prob_i
				} else {
					sub_i = mis_i
				}
				D[i][j] = D[i-1][j-1] + sub_i
				BT_D[i][j][0], BT_D[i][j][1] = 0, 0
				if D[i][j] > IS[i-1][j-1]+sub_i {
					D[i][j] = IS[i-1][j-1] + sub_i
					BT_D[i][j][0], BT_D[i][j][1] = 0, 1
				}
				if D[i][j] > IT[i-1][j-1]+sub_i {
					D[i][j] = IT[i-1][j-1] + sub_i
					BT_D[i][j][0], BT_D[i][j][1] = 0, 2
				}

				IS[i][j] = D[i-1][j] + ins_i_open
				BT_IS[i][j][0], BT_IS[i][j][1] = 1, 0
				if IS[i][j] > IS[i-1][j]+ins_i_ext {
					IS[i][j] = IS[i-1][j] + ins_i_ext
					BT_IS[i][j][0], BT_IS[i][j][1] = 1, 1
				}

				IT[i][j] = D[i][j-1] + PARA_INFO.Gap_open_cost
				BT_IT[i][j][0], BT_IT[i][j][1] = 2, 0
				if IT[i][j] > IT[i][j-1]+PARA_INFO.Gap_ext_cost {
					IT[i][j] = IT[i][j-1] + PARA_INFO.Gap_ext_cost
					BT_IT[i][j][0], BT_IT[i][j][1] = 2, 2
				}
			} else {
				D[i][j] = float64(math.MaxFloat32)
				IT[i][j] = float64(math.MaxFloat32)
				selected_var_len = 0
				var_prof, _ = VC.VarProb[uint32(ref_pos_map[N-j])]
				for var_str, var_prob = range var_prof {
					var_len = len(var_str)
					//One possible case: i - var_len < 0 for all k
					if i-var_len >= 0 {
						prob_i_var = AlignCostKnownLoci(read[M-i:M-i+var_len], []byte(var_str),
							qual[M-i:M-i+var_len], var_prob)
						if D[i][j] > D[i-var_len][j-1]+prob_i_var {
							D[i][j] = D[i-var_len][j-1] + prob_i_var
							BT_D[i][j][0], BT_D[i][j][1] = 0, 0
							selected_var_len = len(var_str)
						}
						/*
							if D[i][j] > IS[i - var_len][j - 1] + prob_i_var {
								D[i][j] = IS[i - var_len][j - 1] + prob_i_var
								BT_D[i][j][0], BT_D[i][j][1] = 0, 1
								selected_var_len = len(var_str)
							}
						*/
						if D[i][j] > IT[i-var_len][j-1]+prob_i_var {
							D[i][j] = IT[i-var_len][j-1] + prob_i_var
							BT_D[i][j][0], BT_D[i][j][1] = 0, 2
							selected_var_len = len(var_str)
						}
					}
				}
				if selected_var_len != 0 {
					BT_D[i][j][2] = selected_var_len
				}
				IS[i][j] = D[i-1][j] + ins_i_open
				BT_IS[i][j][0], BT_IS[i][j][1] = 1, 0
				if IS[i][j] > IS[i-1][j]+ins_i_ext {
					IS[i][j] = IS[i-1][j] + ins_i_ext
					BT_IS[i][j][0], BT_IS[i][j][1] = 1, 1
				}
			}
		}
	}
	PrintDisInfo("FwEditDist, D dis", m, n, D[m][n])
	PrintDisInfo("FwEditDist, IS dis", m, n, IS[m][n])
	PrintDisInfo("FwEditDist, IT dis", m, n, IT[m][n])

	PrintEditDisMat("FwEditDist, D mat", D, m, n, read[M-m:M], ref[N-n:N])
	PrintEditDisMat("FwEditDist, IS mat", IS, m, n, read[M-m:M], ref[N-n:N])
	PrintEditDisMat("FwEditDist, IT mat", IT, m, n, read[M-m:M], ref[N-n:N])

	PrintEditTraceMat("FwEditDist, D trace mat", BT_D, m, n)
	PrintEditTraceMat("FwEditDist, IS trace mat", BT_IS, m, n)
	PrintEditTraceMat("FwEditDist, IT trace mat", BT_IT, m, n)

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
	return align_prob, min_dist, bt_mat, m, n, var_pos, var_base, var_qual, var_type
}

//-------------------------------------------------------------------------------------------------
// ForwardTraceBack constructs alignment between a read and a ref from ForwardDistance.
// The read includes standard bases, the ref include standard bases and "*" characters.
//-------------------------------------------------------------------------------------------------
func (VC *VarCall) ForwardTraceBack(read, qual, ref []byte, m, n int, pos int, BT_Mat int,
	BT_D, BT_IS, BT_IT [][][]int, ref_pos_map []int) ([]int, [][]byte, [][]byte, []int) {

	PrintEditDisInput("FwEditTraceBack, read, qual, ref", read, qual, ref)

	var is_var, is_same_len_var, is_del bool
	var var_len int
	var var_pos, var_type []int
	var var_base, var_qual [][]byte

	aligned_read, aligned_qual, aligned_ref := make([]byte, 0), make([]byte, 0), make([]byte, 0)
	M, N := len(read), len(ref)
	bt_mat := BT_Mat
	i, j, k := m, n, 0
	for i > 0 || j > 0 {
		is_var = false
		if j > 0 {
			_, is_var = INDEX.VarProf[ref_pos_map[N-j]]
		}
		if !is_var { //unknown VARIANT location
			if bt_mat == 0 {
				if read[M-i] != ref[N-j] {
					var_pos = append(var_pos, ref_pos_map[N-j])
					var_base = append(var_base, []byte{read[M-i]})
					var_qual = append(var_qual, []byte{qual[M-i]})
					var_type = append(var_type, 0)
				}
				aligned_read = append(aligned_read, read[M-i])
				aligned_qual = append(aligned_qual, qual[M-i])
				aligned_ref = append(aligned_ref, ref[N-j])
				GetEditTrace("0", M-i, N-j, read[M-i], ref[N-j])
				bt_mat = BT_D[i][j][1]
				i, j = i-1, j-1
			} else if bt_mat == 1 {
				aligned_read = append(aligned_read, read[M-i])
				aligned_qual = append(aligned_qual, qual[M-i])
				aligned_ref = append(aligned_ref, '-')
				GetEditTrace("1", M-i, N-j, read[M-i], '-')
				bt_mat = BT_IS[i][j][1]
				i, j = i-1, j
			} else if bt_mat == 2 {
				aligned_read = append(aligned_read, '-')
				aligned_qual = append(aligned_qual, '-')
				aligned_ref = append(aligned_ref, ref[N-j])
				GetEditTrace("2", M-i, N-j, '-', ref[N-j])
				bt_mat = BT_IT[i][j][1]
				i, j = i, j-1
			}
		} else { //known VARIANT location
			if bt_mat == 0 {
				if BT_D[i][j][2] > 0 {
					var_len = BT_D[i][j][2]
					var_pos = append(var_pos, ref_pos_map[N-j])
					v, q := make([]byte, var_len), make([]byte, var_len)
					copy(v, read[M-i:M-(i-var_len)])
					copy(q, qual[M-i:M-(i-var_len)])
					var_base = append(var_base, v)
					var_qual = append(var_qual, q)
					if _, is_del = INDEX.DelVar[ref_pos_map[N-j]]; is_del {
						var_type = append(var_type, 2)
					} else if _, is_same_len_var = INDEX.SameLenVar[ref_pos_map[N-j]]; is_same_len_var {
						var_type = append(var_type, 0)
					} else {
						var_type = append(var_type, 1)
					}
					aligned_read = append(aligned_read, read[M-i])
					aligned_qual = append(aligned_qual, qual[M-i])
					aligned_ref = append(aligned_ref, ref[N-j])
					for k = 1; k < var_len; k++ {
						aligned_read = append(aligned_read, read[M-i+k])
						aligned_qual = append(aligned_qual, qual[M-i+k])
						aligned_ref = append(aligned_ref, '+')
					}
					GetEditTraceKnownLoc("3", M-i, N-j, read[M-i:M-i+var_len], ref[N-j])
					bt_mat = BT_D[i][j][1]
					i, j = i-var_len, j-1
				} else {
					aligned_read = append(aligned_read, '-')
					aligned_qual = append(aligned_qual, '-')
					aligned_ref = append(aligned_ref, ref[N-j])
					GetEditTrace("4", M-i, N-j, '-', ref[N-j])
					bt_mat = BT_IT[i][j][1]
					i, j = i, j-1
				}
			} else if bt_mat == 1 {
				aligned_read = append(aligned_read, read[M-i])
				aligned_qual = append(aligned_qual, qual[M-i])
				aligned_ref = append(aligned_ref, '-')
				GetEditTrace("1", M-i, N-j, read[M-i], '-')
				bt_mat = BT_IS[i][j][1]
				i, j = i-1, j
			} else {
				aligned_read = append(aligned_read, '-')
				aligned_qual = append(aligned_qual, '-')
				aligned_ref = append(aligned_ref, ref[N-j])
				GetEditTrace("4", M-i, N-j, '-', ref[N-j])
				bt_mat = BT_IT[i][j][1]
				i, j = i, j-1
			}
		}
	}

	PrintEditAlignInfo("FwEditTraceBack, aligned read/qual/ref", aligned_read, aligned_qual, aligned_ref)

	//Get Vars
	ref_ori_pos := N - n
	i = 0
	for i < len(aligned_ref) {
		if aligned_read[i] == '-' && aligned_ref[i] != '-' {
			ref_ori_pos++
			i++
		} else if aligned_read[i] != '-' && aligned_ref[i] == '-' {
			i++
		} else {
			break
		}
	}
	for i < len(aligned_ref) {
		if aligned_read[i] != '-' && aligned_ref[i] == '-' { //Insertions
			v, q := make([]byte, 0), make([]byte, 0)
			v = append(v, aligned_read[i-1])
			q = append(q, aligned_qual[i-1])
			for j = i; j < len(aligned_ref) && aligned_ref[j] == '-'; j++ {
				v = append(v, aligned_read[j])
				q = append(q, aligned_qual[j])
			}
			if j < len(aligned_ref) {
				var_pos = append(var_pos, ref_pos_map[ref_ori_pos-1])
				var_base = append(var_base, v)
				var_qual = append(var_qual, q)
				var_type = append(var_type, 1)
			}
			i = j
		} else if aligned_read[i] == '-' && aligned_ref[i] != '-' { //Deletions
			v, q := make([]byte, 0), make([]byte, 0)
			v = append(v, aligned_ref[i-1])
			//A temporary solution, need to get quality in a proper way in this case!!!
			q = append(q, aligned_qual[i-1])
			for j = i; j < len(aligned_read) && aligned_read[j] == '-'; j++ {
				v = append(v, aligned_ref[j])
			}
			if j < len(aligned_read) {
				var_pos = append(var_pos, ref_pos_map[ref_ori_pos-1])
				var_base = append(var_base, v)
				var_qual = append(var_qual, q)
				var_type = append(var_type, 2)
			}
			ref_ori_pos += j - i
			i = j
		} else if aligned_ref[i] == '+' {
			i++
		} else {
			ref_ori_pos++
			i++
		}
	}
	PrintVarInfo("FwEditTraceBack, variant info", var_pos, var_base, var_qual)
	return var_pos, var_base, var_qual, var_type
}