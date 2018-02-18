//-------------------------------------------------------------------------------------------------
// IVC: alignment.go
// Calculating alignment between reads and multigenomes, take into account known variants.
// Alignment is performed for left and right extensions of seeds on reads and multigenomes.
// Copyright 2015 Nam Sy Vo.
//-------------------------------------------------------------------------------------------------

package ivc

import (
	"math"
)

//-------------------------------------------------------------------------------------------------
// AlignCostVarLoci calculates cost of alignment between a read and the reference at known loci.
//-------------------------------------------------------------------------------------------------
func AlignCostVarLoci(read, ref, qual []byte, prob float64) float64 {
	//do not consider qual at this time
	if string(read) == string(ref) {
		return -0.1 * math.Log10(prob)
	} else {
		return -float64(len(ref)) * math.Log10(INDEL_ERR_RATE)
	}
}

//-------------------------------------------------------------------------------------------------
// LeftAlign calculates the distance between a read and a ref in backward direction.
// The read include standard bases, the ref includes standard bases and "*" characters.
//-------------------------------------------------------------------------------------------------
func (VC *VarCallIndex) LeftAlign(read, qual, ref []byte, pos int, D, IS, IT [][]float64,
	BT_D, BT_IS, BT_IT [][][]int, BT_K [][][]byte, ref_pos_map []int, del_ref bool) (float64, float64,
	int, int, int, []int, [][]byte, [][]byte, []int) {

	var var_len, indel_backup_pos int
	var var_val []byte
	var is_var, is_same_len_var bool
	var p, min_p, var_prob float64

	aln_dist := 0.0
	m, n := len(read), len(ref)

	if PARA.Debug_mode {
		PrintEditDisInput("LeftAlign input: read, qual, ref", pos, read, qual, ref)
	}
	var var_pos, var_type []int
	var var_base, var_qual [][]byte
	var_pos_trace := make(map[int]bool)
	var k int
	for m > 0 && n > 0 {
		indel_backup_pos = ref_pos_map[n-1] - PARA.Indel_backup
		if indel_backup_pos < 0 {
			indel_backup_pos = 0
		} else if indel_backup_pos > VC.SeqLen-1 {
			indel_backup_pos = VC.SeqLen - 1
		}
		if VC.Seq[indel_backup_pos] == '*' {
			if _, is_same_len_var = VC.SameLenVar[indel_backup_pos]; !is_same_len_var {
				break
			}
		}
		if VC.Seq[ref_pos_map[n-1]] != '*' {
			if read[m-1] != ref[n-1] {
				backup_num := PARA.Ham_backup
				if backup_num >= len(read)-m {
					backup_num = len(read) - m
				}
				for i := 0; i < backup_num; i++ {
					if _, is_var = var_pos_trace[n+i]; is_var {
						var_pos = var_pos[:len(var_pos)-1]
						var_base = var_base[:len(var_base)-1]
						var_qual = var_qual[:len(var_qual)-1]
						var_type = var_type[:len(var_type)-1]
					}
				}
				m += backup_num
				n += backup_num
				break
			}
			mapMutex.RLock()
			if _, is_var = VarCall[PARA.Proc_num*ref_pos_map[n-1]/VC.SeqLen].VarType[uint32(ref_pos_map[n-1])]; is_var {
				var_pos_trace[n-1] = true
				var_pos = append(var_pos, ref_pos_map[n-1])
				var_base = append(var_base, []byte{ref[n-1], '|', read[m-1]})
				var_qual = append(var_qual, []byte{qual[m-1]})
				var_type = append(var_type, 0)
			}
			mapMutex.RUnlock()
			m--
			n--
		} else if var_len, is_same_len_var = VC.SameLenVar[ref_pos_map[n-1]]; is_same_len_var {
			min_p = math.MaxFloat64
			for k, var_val = range VC.Variants[ref_pos_map[n-1]] {
				var_prob = float64(VC.VarAF[ref_pos_map[n-1]][k])
				if m >= var_len {
					p = AlignCostVarLoci(read[m-var_len:m], var_val, qual[m-var_len:m], var_prob)
					if min_p > p {
						min_p = p
					}
				}
			}
			if min_p < math.MaxFloat64 {
				aln_dist = aln_dist + min_p
				var_pos_trace[n-1] = true
				var_pos = append(var_pos, ref_pos_map[n-1])
				v, q := make([]byte, 2*var_len+1), make([]byte, var_len)
				copy(v[:var_len], VC.Variants[ref_pos_map[n-1]][0])
				copy(v[var_len:var_len+1], []byte{'|'})
				copy(v[var_len+1:], read[m-var_len:m])
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
		if aln_dist > PARA.Dist_thres {
			return PARA.Dist_thres + 1, 0, -1, m, n, var_pos, var_base, var_qual, var_type
		}
	}
	if PARA.Debug_mode {
		PrintDisInfo("LeftAlnHam dis", m, n, aln_dist)
	}
	if m == 0 || n == 0 {
		return aln_dist, 0, -1, m, n, var_pos, var_base, var_qual, var_type
	}
	if PARA.Debug_mode {
		PrintEditDisInput("LeftAlnEdit: read, qual, ref", pos, read[:m], qual[:m], ref[:n])
	}
	/*
		Backtrace info matrices:
		BT_K[i][j][2]: represents number of shifted bases (equal to length of called variants) at known variant locations,
			can be any integer number, e.g. 5 means back to i-5,j-1.
		For each BT_x[i][j] (x can be D, IS, or IT):
		BT_x[i][j][0]: represents direction to trace back to, can be 0: diagonal arrow (back to i-1,j-1), 1: up arrow (back to i-1,j),
		 	2: left arrow (back to i,j-1).
		BT_x[i][j][1]: represents matrix to trace back to, can be 0: trace back to matrix D, 1: trace back to matrix IS, 2: trace back to matrix IT.
	*/
	var i, j int
	for i := 0; i <= m; i++ {
		for j := 0; j <= n; j++ {
			BT_K[i][j] = nil
			BT_D[i][j][0], BT_D[i][j][1] = -1, -1
			BT_IS[i][j][0], BT_IS[i][j][1] = -1, -1
			BT_IT[i][j][0], BT_IT[i][j][1] = -1, -1
		}
	}

	D[0][0] = 0.0
	IS[0][0] = float64(math.MaxFloat32)
	IT[0][0] = float64(math.MaxFloat32)
	IS[1][0] = PARA.Gap_open
	BT_IS[1][0][0], BT_IS[1][0][1] = 1, 1

	for i = 1; i <= m; i++ {
		D[i][0] = float64(math.MaxFloat32)
		IT[i][0] = float64(math.MaxFloat32)
	}
	for i = 2; i <= m; i++ {
		IS[i][0] = PARA.Gap_ext
		BT_IS[i][0][0], BT_IS[i][0][1] = 1, 1
	}

	for j = 1; j <= n; j++ {
		D[0][j] = float64(math.MaxFloat32)
		IS[0][j] = float64(math.MaxFloat32)
		IT[0][j] = 0.0
		BT_IT[0][j][0], BT_IT[0][j][1] = 2, 2
	}

	var sel_var []byte
	var prob_i, sub_i, mis_i float64
	var is_del bool
	for i = 1; i <= m; i++ {
		mis_i = PARA.Sub_cost // + Q2C[qual[i-1]]
		for j = 1; j <= n; j++ {
			if VC.Seq[ref_pos_map[j-1]] != '*' {
				if read[i-1] == ref[j-1] {
					sub_i = 0.0
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

				IS[i][j] = D[i-1][j] + PARA.Gap_open
				BT_IS[i][j][0], BT_IS[i][j][1] = 1, 0
				if IS[i][j] > IS[i-1][j]+PARA.Gap_ext {
					IS[i][j] = IS[i-1][j] + PARA.Gap_ext
					BT_IS[i][j][0], BT_IS[i][j][1] = 1, 1
				}

				IT[i][j] = D[i][j-1] + PARA.Gap_open
				BT_IT[i][j][0], BT_IT[i][j][1] = 2, 0
				if IT[i][j] > IT[i][j-1]+PARA.Gap_ext {
					IT[i][j] = IT[i][j-1] + PARA.Gap_ext
					BT_IT[i][j][0], BT_IT[i][j][1] = 2, 2
				}
			} else {
				D[i][j] = float64(math.MaxFloat32)
				IS[i][j] = float64(math.MaxFloat32)
				IT[i][j] = float64(math.MaxFloat32)
				sel_var = nil
				for k, var_val = range VC.Variants[ref_pos_map[j-1]] {
					var_prob = float64(VC.VarAF[ref_pos_map[j-1]][k])
					var_len = len(var_val)
					if i-var_len >= 0 {
						if _, is_del = VC.DelVar[ref_pos_map[j-1]]; is_del && del_ref {
							prob_i = AlignCostVarLoci(read[i-var_len:i], var_val, qual[i-var_len:i], 1.0-var_prob)
						} else {
							prob_i = AlignCostVarLoci(read[i-var_len:i], var_val, qual[i-var_len:i], var_prob)
						}
						if D[i][j] > D[i-var_len][j-1]+prob_i {
							D[i][j] = D[i-var_len][j-1] + prob_i
							BT_D[i][j][0], BT_D[i][j][1] = 0, 0
							sel_var = var_val
						}
						if D[i][j] > IS[i-var_len][j-1]+prob_i {
							D[i][j] = IS[i-var_len][j-1] + prob_i
							BT_D[i][j][0], BT_D[i][j][1] = 0, 1
							sel_var = var_val
						}
						if D[i][j] > IT[i-var_len][j-1]+prob_i {
							D[i][j] = IT[i-var_len][j-1] + prob_i
							BT_D[i][j][0], BT_D[i][j][1] = 0, 2
							sel_var = var_val
						}
					}
				}
				if sel_var != nil {
					BT_K[i][j] = sel_var
				}
			}
		}
	}
	if PARA.Debug_mode {
		PrintDisInfo("LeftAlnEditDist, D dis", m, n, D[m][n])
		PrintDisInfo("LeftAlnEditDist, IS dis", m, n, IS[m][n])
		PrintDisInfo("LeftAlnEditDist, IT dis", m, n, IT[m][n])

		PrintEditDisMat("LeftAlnEditDist, D mat", D, m, n, read[:m], ref[:n])
		PrintEditDisMat("LeftAlnEditDist, IS mat", IS, m, n, read[:m], ref[:n])
		PrintEditDisMat("LeftAlnEditDist, IT mat", IT, m, n, read[:m], ref[:n])

		PrintEditTraceMat("LeftAlnEditDist, D trace mat", BT_D, m, n)
		PrintEditTraceMat("LeftAlnEditDist, IS trace mat", BT_IS, m, n)
		PrintEditTraceMat("LeftAlnEditDist, IT trace mat", BT_IT, m, n)
	}
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

	return aln_dist, min_dist, bt_mat, m, n, var_pos, var_base, var_qual, var_type
}

//-------------------------------------------------------------------------------------------------
// LeftAlignEditTraceBack constructs alignment between a read and a ref from LeftAlign.
// The read includes standard bases, the ref include standard bases and "*" characters.
//-------------------------------------------------------------------------------------------------
func (VC *VarCallIndex) LeftAlignEditTraceBack(read, qual, ref []byte, m, n int, pos int,
	BT_Mat int, BT_D, BT_IS, BT_IT [][][]int, BT_K [][][]byte, ref_pos_map []int, del_ref bool) ([]int, [][]byte, [][]byte, []int) {

	var var_len, ref_len int
	var var_pos, var_type []int
	var var_base, var_qual [][]byte
	var is_same_len_var, is_del bool
	if PARA.Debug_mode {
		PrintEditDisInput("LeftAlnEditTraceBack, read, qual, ref", pos, read[:m], qual[:m], ref[:n])
	}
	aln_read, aln_qual, aln_ref := make([]byte, 0), make([]byte, 0), make([]byte, 0)
	bt_mat := BT_Mat
	i, j, k := m, n, 0
	for i > 0 || j > 0 {
		if j == 0 || VC.Seq[ref_pos_map[j-1]] != '*' { //unknown VARIANT location
			if bt_mat == 0 {
				if read[i-1] != ref[j-1] {
					var_pos = append(var_pos, ref_pos_map[j-1])
					var_base = append(var_base, []byte{ref[j-1], '|', read[i-1]})
					var_qual = append(var_qual, []byte{qual[i-1]})
					var_type = append(var_type, 0)
				}
				aln_read = append(aln_read, read[i-1])
				aln_qual = append(aln_qual, qual[i-1])
				aln_ref = append(aln_ref, ref[j-1])
				//GetEditTrace("0", i, j, read[i-1], ref[j-1])
				bt_mat = BT_D[i][j][1]
				i, j = i-1, j-1
			} else if bt_mat == 1 {
				aln_read = append(aln_read, read[i-1])
				aln_qual = append(aln_qual, qual[i-1])
				aln_ref = append(aln_ref, '-')
				//GetEditTrace("1", i, j, read[i-1], '-')
				bt_mat = BT_IS[i][j][1]
				i, j = i-1, j
			} else if bt_mat == 2 {
				aln_read = append(aln_read, '-')
				aln_qual = append(aln_qual, '-')
				aln_ref = append(aln_ref, ref[j-1])
				//GetEditTrace("2", i, j, '-', ref[j-1])
				bt_mat = BT_IT[i][j][1]
				i, j = i, j-1
			}
		} else { //known VARIANT location
			if BT_K[i][j] != nil {
				var_len = len(BT_K[i][j])
				var_pos = append(var_pos, ref_pos_map[j-1])
				ref_len = len(VC.Variants[ref_pos_map[j-1]][0])
				var v []byte
				if _, is_del = VC.DelVar[ref_pos_map[j-1]]; is_del && !del_ref { //known DEL with non-reduced ref
					v = make([]byte, ref_len+ref_len+1)
					copy(v[:ref_len], VC.Variants[ref_pos_map[j-1]][0])
					copy(v[ref_len:ref_len+1], []byte{'|'})
					copy(v[ref_len+1:], VC.Variants[ref_pos_map[j-1]][0])
				} else {
					v = make([]byte, ref_len+var_len+1)
					copy(v[:ref_len], VC.Variants[ref_pos_map[j-1]][0])
					copy(v[ref_len:ref_len+1], []byte{'|'})
					copy(v[ref_len+1:], BT_K[i][j])
				}
				var_base = append(var_base, v)
				q := make([]byte, var_len)
				copy(q, qual[i-var_len:i])
				var_qual = append(var_qual, q)
				if _, is_del = VC.DelVar[ref_pos_map[j-1]]; is_del {
					var_type = append(var_type, 2)
				} else if _, is_same_len_var = VC.SameLenVar[ref_pos_map[j-1]]; is_same_len_var {
					var_type = append(var_type, 0)
				} else {
					var_type = append(var_type, 1)
				}
				for k = 0; k < var_len-1; k++ {
					aln_read = append(aln_read, read[i-1-k])
					aln_qual = append(aln_qual, qual[i-1-k])
					aln_ref = append(aln_ref, '+')
				}
				aln_read = append(aln_read, read[i-var_len])
				aln_qual = append(aln_qual, qual[i-var_len])
				aln_ref = append(aln_ref, ref[j-1])
				//GetEditTraceKnownLoc("3", i, j, read[i-var_len:i], ref[j-1])
				bt_mat = BT_D[i][j][1]
				i, j = i-var_len, j-1
			} else {
				aln_read = append(aln_read, '-')
				aln_qual = append(aln_qual, '-')
				aln_ref = append(aln_ref, ref[j-1])
				//GetEditTraceKnownLoc("4", i, j, []byte{'-'}, ref[j-1])
				bt_mat = BT_IT[i][j][1]
				i, j = i, j-1
			}
		}
	}

	//Put the alignment in original direction
	for i, j = 0, len(aln_read)-1; i < j; i, j = i+1, j-1 {
		aln_read[i], aln_read[j] = aln_read[j], aln_read[i]
		aln_qual[i], aln_qual[j] = aln_qual[j], aln_qual[i]
		aln_ref[i], aln_ref[j] = aln_ref[j], aln_ref[i]
	}
	if PARA.Debug_mode {
		PrintEditAlignInfo("LeftAlnEditTraceBack, aligned read/qual/ref", aln_read, aln_qual, aln_ref)
	}
	//Get variants
	ref_ori_pos := 0
	read_ori_pos := 0
	i = 0
	for i < len(aln_ref) {
		if aln_read[i] == '-' && aln_ref[i] != '-' {
			ref_ori_pos++
			i++
		} else if aln_read[i] != '-' && aln_ref[i] == '-' {
			read_ori_pos++
			i++
		} else {
			break
		}
	}
	for i < len(aln_ref) {
		if aln_read[i] != '-' && aln_ref[i] == '-' { //Insertions
			v, q := make([]byte, 0), make([]byte, 0)
			v = append(v, aln_ref[i-1])
			v = append(v, '|')
			v = append(v, aln_read[i-1])
			q = append(q, aln_qual[i-1])
			for j = i; j < len(aln_ref) && aln_ref[j] == '-'; j++ {
				v = append(v, aln_read[j])
				q = append(q, aln_qual[j])
			}
			if j < len(aln_ref)-1 && read_ori_pos > 1 {
				var_pos = append(var_pos, ref_pos_map[ref_ori_pos-1])
				var_base = append(var_base, v)
				var_qual = append(var_qual, q)
				var_type = append(var_type, 1)
			}
			read_ori_pos += j - i
			i = j
		} else if aln_read[i] == '-' && aln_ref[i] != '-' { //Deletions
			v, q := make([]byte, 0), make([]byte, 0)
			v = append(v, aln_ref[i-1])
			q = append(q, aln_qual[i-1]) //A temporary solution, need to get quality in a proper way in this case!!!
			for j = i; j < len(aln_read) && aln_read[j] == '-'; j++ {
				v = append(v, aln_ref[j])
			}
			if j < len(aln_read)-1 && read_ori_pos < m-1 {
				var_pos = append(var_pos, ref_pos_map[ref_ori_pos-1])
				v = append(v, '|')
				v = append(v, aln_read[i-1])
				var_base = append(var_base, v)
				var_qual = append(var_qual, q)
				var_type = append(var_type, 2)
			}
			ref_ori_pos += j - i
			i = j
		} else if aln_ref[i] == '+' {
			read_ori_pos++
			i++
		} else {
			if aln_read[i] == aln_ref[i] && i+1 < len(aln_read) && aln_read[i+1] != '-' && aln_ref[i+1] != '-' {
				if ref_pos_map != nil {
					mapMutex.RLock()
					if _, is_prof_new_var := VarCall[PARA.Proc_num*ref_pos_map[ref_ori_pos]/VC.SeqLen].VarType[uint32(ref_pos_map[ref_ori_pos])]; is_prof_new_var {
						var_pos = append(var_pos, ref_pos_map[ref_ori_pos])
						var_base = append(var_base, []byte{aln_ref[i], '|', aln_read[i]})
						var_qual = append(var_qual, []byte{aln_qual[i]})
						var_type = append(var_type, 0)
					}
					mapMutex.RUnlock()
				}
			}
			ref_ori_pos++
			read_ori_pos++
			i++
		}
	}
	return var_pos, var_base, var_qual, var_type
}

//-------------------------------------------------------------------------------------------------
// RightAlign calculates the distance between a read and a ref in forward direction.
// The read includes standard bases, the ref includes standard bases and "*" characters.
//-------------------------------------------------------------------------------------------------
func (VC *VarCallIndex) RightAlign(read, qual, ref []byte, pos int, D, IS, IT [][]float64,
	BT_D, BT_IS, BT_IT [][][]int, BT_K [][][]byte, ref_pos_map []int, del_ref bool) (float64, float64,
	int, int, int, []int, [][]byte, [][]byte, []int) {

	var var_len, indel_backup_pos int
	var is_var, is_same_len_var bool
	var var_val []byte
	var p, min_p, var_prob float64
	var var_pos, var_type []int
	var var_base, var_qual [][]byte
	var k int

	if PARA.Debug_mode {
		PrintEditDisInput("RightAlign input: read, qual, ref", pos, read, qual, ref)
	}
	aln_dist := 0.0
	M, N := len(read), len(ref)
	m, n := M, N
	var_pos_trace := make(map[int]bool)
	for m > 0 && n > 0 {
		indel_backup_pos = ref_pos_map[N-n] + PARA.Indel_backup
		if indel_backup_pos < 0 {
			indel_backup_pos = 0
		} else if indel_backup_pos > VC.SeqLen-1 {
			indel_backup_pos = VC.SeqLen - 1
		}
		if VC.Seq[indel_backup_pos] == '*' {
			if _, is_same_len_var = VC.SameLenVar[indel_backup_pos]; !is_same_len_var {
				break
			}
		}
		if VC.Seq[ref_pos_map[N-n]] != '*' {
			if read[M-m] != ref[N-n] {
				backup_num := 2 * PARA.Ham_backup
				if backup_num >= M-m {
					backup_num = M - m
				}
				for i := 0; i < backup_num; i++ {
					if _, is_var = var_pos_trace[N-(n+i+1)]; is_var {
						var_pos = var_pos[:len(var_pos)-1]
						var_base = var_base[:len(var_base)-1]
						var_qual = var_qual[:len(var_qual)-1]
						var_type = var_type[:len(var_type)-1]
					}
				}
				m += backup_num
				n += backup_num
				break
			}
			mapMutex.RLock()
			if _, is_var = VarCall[PARA.Proc_num*ref_pos_map[N-n]/VC.SeqLen].VarType[uint32(ref_pos_map[N-n])]; is_var {
				var_pos_trace[N-n] = true
				var_pos = append(var_pos, ref_pos_map[N-n])
				var_base = append(var_base, []byte{ref[N-n], '|', read[M-m]})
				var_qual = append(var_qual, []byte{qual[M-m]})
				var_type = append(var_type, 0)
			}
			mapMutex.RUnlock()
			m--
			n--
		} else if var_len, is_same_len_var = VC.SameLenVar[ref_pos_map[N-n]]; is_same_len_var {
			min_p = math.MaxFloat64
			for k, var_val = range VC.Variants[ref_pos_map[N-n]] {
				var_prob = float64(VC.VarAF[ref_pos_map[N-n]][k])
				if m >= var_len {
					p = AlignCostVarLoci(read[M-m:M-m+var_len], var_val, qual[M-m:M-m+var_len], var_prob)
					if min_p > p {
						min_p = p
					}
				}
			}
			if min_p < math.MaxFloat64 {
				aln_dist = aln_dist + min_p
				var_pos_trace[N-n] = true
				var_pos = append(var_pos, ref_pos_map[N-n])
				v, q := make([]byte, 2*var_len+1), make([]byte, var_len)
				copy(v[:var_len], VC.Variants[ref_pos_map[N-n]][0])
				copy(v[var_len:var_len+1], []byte{'|'})
				copy(v[var_len+1:], read[M-m:M-(m-var_len)])
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
		if aln_dist > PARA.Dist_thres {
			return PARA.Dist_thres + 1, 0, -1, m, n, var_pos, var_base, var_qual, var_type
		}
	}
	if PARA.Debug_mode {
		PrintDisInfo("RightAlnHam dis", m, n, aln_dist)
	}
	if m == 0 || n == 0 {
		return aln_dist, 0, -1, m, n, var_pos, var_base, var_qual, var_type
	}
	if PARA.Debug_mode {
		PrintEditDisInput("RightAlnEdit: read, qual, ref", pos, read[M-m:M], qual[M-m:M], ref[N-n:N])
	}
	//	Backtrace info matrices:
	//	BT_K[i][j]: represents number of shifted bases (equal to length of called variants) at known variant locations,
	//		can be any integer number, e.g. 5 means back to i-5,j-1.
	//  For each BT_x[i][j] (x can be D, IS, or IT):
	//	BT_x[i][j][0]: represents direction to trace back to, can be 0: diagonal arrow (back to i-1,j-1), 1: up arrow (back to i-1,j),
	//	 	2: left arrow (back to i,j-1).
	//	BT_x[i][j][1]: represents matrix to trace back to, can be 0: trace back to matrix D, 1: trace back to matrix IS, 2: trace back to matrix IT.
	var i, j int
	for i := 0; i <= m; i++ {
		for j := 0; j <= n; j++ {
			BT_K[i][j] = nil
			BT_D[i][j][0], BT_D[i][j][1] = -1, -1
			BT_IS[i][j][0], BT_IS[i][j][1] = -1, -1
			BT_IT[i][j][0], BT_IT[i][j][1] = -1, -1
		}
	}

	D[0][0] = 0.0
	for i = 1; i <= m; i++ {
		D[i][0] = float64(math.MaxFloat32)
		IT[i][0] = float64(math.MaxFloat32)
	}
	IS[0][0] = float64(math.MaxFloat32)
	IS[1][0] = PARA.Gap_open
	BT_IS[1][0][0], BT_IS[1][0][1] = 1, 1
	for i = 2; i <= m; i++ {
		IS[i][0] = PARA.Gap_ext
		BT_IS[i][0][0], BT_IS[i][0][1] = 1, 1
	}

	IT[0][0] = float64(math.MaxFloat32)
	for j = 1; j <= n; j++ {
		D[0][j] = float64(math.MaxFloat32)
		IS[0][j] = float64(math.MaxFloat32)
		IT[0][j] = 0.0
		BT_IT[0][j][0], BT_IT[0][j][1] = 2, 2
	}

	var sel_var []byte
	var prob_i, sub_i, mis_i float64
	var is_del bool
	for i = 1; i <= m; i++ {
		mis_i = PARA.Sub_cost // + Q2C[qual[M-i]]
		for j = 1; j <= n; j++ {
			if N-j < 0 || N-j >= len(ref_pos_map) {
				panic("ref_pos_map index problem")
			}
			if ref_pos_map[N-j] < 0 || ref_pos_map[N-j] > len(VC.Seq) {
				panic("VC.Seq index problem")
			}
			if VC.Seq[ref_pos_map[N-j]] != '*' {
				if read[M-i] == ref[N-j] {
					sub_i = 0.0
				} else {
					sub_i = mis_i
				}
				D[i][j] = IT[i-1][j-1] + sub_i
				BT_D[i][j][0], BT_D[i][j][1] = 0, 2
				if D[i][j] > IS[i-1][j-1]+sub_i {
					D[i][j] = IS[i-1][j-1] + sub_i
					BT_D[i][j][0], BT_D[i][j][1] = 0, 1
				}
				if D[i][j] > D[i-1][j-1]+sub_i {
					D[i][j] = D[i-1][j-1] + sub_i
					BT_D[i][j][0], BT_D[i][j][1] = 0, 0
				}
				IS[i][j] = D[i-1][j] + PARA.Gap_open
				BT_IS[i][j][0], BT_IS[i][j][1] = 1, 0
				if IS[i][j] > IS[i-1][j]+PARA.Gap_ext {
					IS[i][j] = IS[i-1][j] + PARA.Gap_ext
					BT_IS[i][j][0], BT_IS[i][j][1] = 1, 1
				}
				IT[i][j] = D[i][j-1] + PARA.Gap_open
				BT_IT[i][j][0], BT_IT[i][j][1] = 2, 0
				if IT[i][j] > IT[i][j-1]+PARA.Gap_ext {
					IT[i][j] = IT[i][j-1] + PARA.Gap_ext
					BT_IT[i][j][0], BT_IT[i][j][1] = 2, 2
				}
			} else {
				D[i][j] = float64(math.MaxFloat32)
				IT[i][j] = float64(math.MaxFloat32)
				sel_var = nil
				for k, var_val = range VC.Variants[ref_pos_map[N-j]] {
					var_prob = float64(VC.VarAF[ref_pos_map[N-j]][k])
					var_len = len(var_val)
					if i-var_len >= 0 {
						if _, is_del = VC.DelVar[ref_pos_map[N-j]]; is_del && del_ref { //convert prob with reduced-ref for known DEL
							prob_i = AlignCostVarLoci(read[M-i:M-i+var_len], var_val, qual[M-i:M-i+var_len], 1.0-var_prob)
						} else {
							prob_i = AlignCostVarLoci(read[M-i:M-i+var_len], var_val, qual[M-i:M-i+var_len], var_prob)
						}
						if D[i][j] > D[i-var_len][j-1]+prob_i {
							D[i][j] = D[i-var_len][j-1] + prob_i
							BT_D[i][j][0], BT_D[i][j][1] = 0, 0
							sel_var = var_val
						}
						/*
							if D[i][j] > IS[i - var_len][j - 1] + prob_i {
								D[i][j] = IS[i - var_len][j - 1] + prob_i
								BT_D[i][j][0], BT_D[i][j][1] = 0, 1
								sel_var = var_val
							}
						*/
						if D[i][j] > IT[i-var_len][j-1]+prob_i {
							D[i][j] = IT[i-var_len][j-1] + prob_i
							BT_D[i][j][0], BT_D[i][j][1] = 0, 2
							sel_var = var_val
						}
					}
				}
				if sel_var != nil {
					BT_K[i][j] = sel_var
				}
				IS[i][j] = D[i-1][j] + PARA.Gap_open
				BT_IS[i][j][0], BT_IS[i][j][1] = 1, 0
				if IS[i][j] > IS[i-1][j]+PARA.Gap_ext {
					IS[i][j] = IS[i-1][j] + PARA.Gap_ext
					BT_IS[i][j][0], BT_IS[i][j][1] = 1, 1
				}
			}
		}
	}
	if PARA.Debug_mode {
		PrintDisInfo("RightAlnEditDist, D dis", m, n, D[m][n])
		PrintDisInfo("RightAlnEditDist, IS dis", m, n, IS[m][n])
		PrintDisInfo("RightAlnEditDist, IT dis", m, n, IT[m][n])

		PrintEditDisMat("RightAlnEditDist, D mat", D, m, n, read[M-m:M], ref[N-n:N])
		PrintEditDisMat("RightAlnEditDist, IS mat", IS, m, n, read[M-m:M], ref[N-n:N])
		PrintEditDisMat("RightAlnEditDist, IT mat", IT, m, n, read[M-m:M], ref[N-n:N])

		PrintEditTraceMat("RightAlnEditDist, D trace mat", BT_D, m, n)
		PrintEditTraceMat("RightAlnEditDist, IS trace mat", BT_IS, m, n)
		PrintEditTraceMat("RightAlnEditDist, IT trace mat", BT_IT, m, n)
	}
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
	return aln_dist, min_dist, bt_mat, m, n, var_pos, var_base, var_qual, var_type
}

//-------------------------------------------------------------------------------------------------
// RightAlignEditTraceBack constructs alignment between a read and a ref from RightAlign.
// The read includes standard bases, the ref include standard bases and "*" characters.
//-------------------------------------------------------------------------------------------------
func (VC *VarCallIndex) RightAlignEditTraceBack(read, qual, ref []byte, m, n int, pos int,
	BT_Mat int, BT_D, BT_IS, BT_IT [][][]int, BT_K [][][]byte, ref_pos_map []int, del_ref bool) ([]int, [][]byte, [][]byte, []int) {

	if PARA.Debug_mode {
		PrintEditDisInput("RightAlnEditTraceBack, read, qual, ref", pos, read, qual, ref)
	}
	var var_len, ref_len int
	var var_pos, var_type []int
	var var_base, var_qual [][]byte
	var is_same_len_var, is_del bool

	aln_read, aln_qual, aln_ref := make([]byte, 0), make([]byte, 0), make([]byte, 0)
	M, N := len(read), len(ref)
	bt_mat := BT_Mat
	i, j, k := m, n, 0
	for i > 0 || j > 0 {
		if j == 0 || VC.Seq[ref_pos_map[N-j]] != '*' { //unknown VARIANT location
			if bt_mat == 0 {
				if read[M-i] != ref[N-j] {
					var_pos = append(var_pos, ref_pos_map[N-j])
					var_base = append(var_base, []byte{ref[N-j], '|', read[M-i]})
					var_qual = append(var_qual, []byte{qual[M-i]})
					var_type = append(var_type, 0)
				}
				aln_read = append(aln_read, read[M-i])
				aln_qual = append(aln_qual, qual[M-i])
				aln_ref = append(aln_ref, ref[N-j])
				//GetEditTrace("0", M-i, N-j, read[M-i], ref[N-j])
				bt_mat = BT_D[i][j][1]
				i, j = i-1, j-1
			} else if bt_mat == 1 {
				aln_read = append(aln_read, read[M-i])
				aln_qual = append(aln_qual, qual[M-i])
				aln_ref = append(aln_ref, '-')
				//GetEditTrace("1", M-i, N-j, read[M-i], '-')
				bt_mat = BT_IS[i][j][1]
				i, j = i-1, j
			} else if bt_mat == 2 {
				aln_read = append(aln_read, '-')
				aln_qual = append(aln_qual, '-')
				aln_ref = append(aln_ref, ref[N-j])
				//GetEditTrace("2", M-i, N-j, '-', ref[N-j])
				bt_mat = BT_IT[i][j][1]
				i, j = i, j-1
			}
		} else { //known VARIANT location
			if bt_mat == 0 {
				if BT_K[i][j] != nil {
					var_len = len(BT_K[i][j])
					var_pos = append(var_pos, ref_pos_map[N-j])
					ref_len = len(VC.Variants[ref_pos_map[N-j]][0])
					var v []byte
					if _, is_del = VC.DelVar[ref_pos_map[N-j]]; is_del && !del_ref { //known DEL with non-reduced ref
						v = make([]byte, ref_len+ref_len+1)
						copy(v[:ref_len], VC.Variants[ref_pos_map[N-j]][0])
						copy(v[ref_len:ref_len+1], []byte{'|'})
						copy(v[ref_len+1:], VC.Variants[ref_pos_map[N-j]][0])
					} else {
						v = make([]byte, ref_len+var_len+1)
						copy(v[:ref_len], VC.Variants[ref_pos_map[N-j]][0])
						copy(v[ref_len:ref_len+1], []byte{'|'})
						copy(v[ref_len+1:], BT_K[i][j])
					}
					var_base = append(var_base, v)
					q := make([]byte, var_len)
					copy(q, qual[M-i:M-(i-var_len)])
					var_qual = append(var_qual, q)
					if _, is_del = VC.DelVar[ref_pos_map[N-j]]; is_del {
						var_type = append(var_type, 2)
					} else if _, is_same_len_var = VC.SameLenVar[ref_pos_map[N-j]]; is_same_len_var {
						var_type = append(var_type, 0)
					} else {
						var_type = append(var_type, 1)
					}
					aln_read = append(aln_read, read[M-i])
					aln_qual = append(aln_qual, qual[M-i])
					aln_ref = append(aln_ref, ref[N-j])
					for k = 1; k < var_len; k++ {
						aln_read = append(aln_read, read[M-i+k])
						aln_qual = append(aln_qual, qual[M-i+k])
						aln_ref = append(aln_ref, '+')
					}
					//GetEditTraceKnownLoc("3", M-i, N-j, read[M-i:M-i+var_len], ref[N-j])
					bt_mat = BT_D[i][j][1]
					i, j = i-var_len, j-1
				} else {
					aln_read = append(aln_read, '-')
					aln_qual = append(aln_qual, '-')
					aln_ref = append(aln_ref, ref[N-j])
					//GetEditTrace("4", M-i, N-j, '-', ref[N-j])
					bt_mat = BT_IT[i][j][1]
					i, j = i, j-1
				}
			} else if bt_mat == 1 {
				aln_read = append(aln_read, read[M-i])
				aln_qual = append(aln_qual, qual[M-i])
				aln_ref = append(aln_ref, '-')
				//GetEditTrace("1", M-i, N-j, read[M-i], '-')
				bt_mat = BT_IS[i][j][1]
				i, j = i-1, j
			} else {
				aln_read = append(aln_read, '-')
				aln_qual = append(aln_qual, '-')
				aln_ref = append(aln_ref, ref[N-j])
				//GetEditTrace("4", M-i, N-j, '-', ref[N-j])
				bt_mat = BT_IT[i][j][1]
				i, j = i, j-1
			}
		}
	}
	if PARA.Debug_mode {
		PrintEditAlignInfo("RightAlnEditTraceBack, aligned read/qual/ref", aln_read, aln_qual, aln_ref)
	}
	//Get variants
	ref_ori_pos := N - n
	read_ori_pos := M - m
	i = 0
	for i < len(aln_ref) {
		if aln_read[i] == '-' && aln_ref[i] != '-' {
			ref_ori_pos++
			i++
		} else if aln_read[i] != '-' && aln_ref[i] == '-' {
			read_ori_pos++
			i++
		} else {
			break
		}
	}
	for i < len(aln_ref) {
		if aln_read[i] != '-' && aln_ref[i] == '-' { //Insertions
			v, q := make([]byte, 0), make([]byte, 0)
			v = append(v, aln_ref[i-1])
			v = append(v, '|')
			v = append(v, aln_read[i-1])
			q = append(q, aln_qual[i-1])
			for j = i; j < len(aln_ref) && aln_ref[j] == '-'; j++ {
				v = append(v, aln_read[j])
				q = append(q, aln_qual[j])
			}
			if j < len(aln_ref)-1 && read_ori_pos+j-i < M-1 && read_ori_pos > M-m+1 {
				var_pos = append(var_pos, ref_pos_map[ref_ori_pos-1])
				var_base = append(var_base, v)
				var_qual = append(var_qual, q)
				var_type = append(var_type, 1)
			}
			read_ori_pos += j - i
			i = j
		} else if aln_read[i] == '-' && aln_ref[i] != '-' { //Deletions
			v, q := make([]byte, 0), make([]byte, 0)
			v = append(v, aln_ref[i-1])
			//A temporary solution, need to get quality in a proper way in this case!!!
			q = append(q, aln_qual[i-1])
			for j = i; j < len(aln_read) && aln_read[j] == '-'; j++ {
				v = append(v, aln_ref[j])
			}
			if j < len(aln_read)-1 && read_ori_pos < M-1 && read_ori_pos > M-m+1 {
				var_pos = append(var_pos, ref_pos_map[ref_ori_pos-1])
				v = append(v, '|')
				v = append(v, aln_read[i-1])
				var_base = append(var_base, v)
				var_qual = append(var_qual, q)
				var_type = append(var_type, 2)
			}
			ref_ori_pos += j - i
			i = j
		} else if aln_ref[i] == '+' {
			read_ori_pos++
			i++
		} else {
			if aln_read[i] == aln_ref[i] && i+1 < len(aln_read) && aln_read[i+1] != '-' && aln_ref[i+1] != '-' {
				if ref_pos_map != nil {
					mapMutex.RLock()
					if _, is_prof_new_var := VarCall[PARA.Proc_num*ref_pos_map[ref_ori_pos]/VC.SeqLen].VarType[uint32(ref_pos_map[ref_ori_pos])]; is_prof_new_var {
						var_pos = append(var_pos, ref_pos_map[ref_ori_pos])
						var_base = append(var_base, []byte{aln_ref[i], '|', aln_read[i]})
						var_qual = append(var_qual, []byte{aln_qual[i]})
						var_type = append(var_type, 0)
					}
					mapMutex.RUnlock()
				}
			}
			ref_ori_pos++
			read_ori_pos++
			i++
		}
	}
	return var_pos, var_base, var_qual, var_type
}
