//---------------------------------------------------------------------------------------------------
// IVC: debug.go - Debug variables and functions.
// Copyright 2015 Nam Sy Vo.
//---------------------------------------------------------------------------------------------------

package ivc

import (
	//"bufio"
	//"bytes"
	"fmt"
	"log"
	"math"
	"os"
	"runtime"
	//"sort"
	//"strconv"
	//"strings"
	//"time"
)

//Global variable for turnning on/off info profiling
var (
	PRINT_PROCESS_MEM = true
	PRINT_MEM_STATS   = false

	PRINT_EDIT_DIST_INFO     = false
	PRINT_EDIT_DIST_MAT_INFO = false

	PRINT_VAR_CALL_INFO    = false
	PRINT_ALIGN_TRACE_INFO = false

	GET_ALIGN_READ_INFO = false
	PRINT_FN            = false
	PRINT_TPFP          = false

	GET_NO_ALIGN_READ_INFO = true
	PRINT_NA               = false
)

//Global variable for memory profiling
var (
	MEM_FILE *os.File
	MEM_STATS = new(runtime.MemStats)
)

//Printing memory information
func PrintProcessMem(mesg string) {
	if PRINT_PROCESS_MEM {
		runtime.ReadMemStats(MEM_STATS)
		log.Printf(mesg+"\t%d\t%d\t%d\t%d\t%d\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f",
			MEM_STATS.Alloc, MEM_STATS.TotalAlloc, MEM_STATS.Sys, MEM_STATS.HeapAlloc, MEM_STATS.HeapSys,
			float64(MEM_STATS.Alloc)/(math.Pow(1024, 3)), float64(MEM_STATS.TotalAlloc)/(math.Pow(1024, 3)),
			float64(MEM_STATS.Sys)/(math.Pow(1024, 3)), float64(MEM_STATS.HeapAlloc)/(math.Pow(1024, 3)),
			float64(MEM_STATS.HeapSys)/(math.Pow(1024, 3)))
	}
}

func PrintMemStats(mesg string) {
	if PRINT_MEM_STATS {
		runtime.ReadMemStats(MEM_STATS)
		log.Printf(mesg+"\t%d\t%d\t%d\t%d\t%d\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f",
			MEM_STATS.Alloc, MEM_STATS.TotalAlloc, MEM_STATS.Sys, MEM_STATS.HeapAlloc, MEM_STATS.HeapSys,
			float64(MEM_STATS.Alloc)/(math.Pow(1024, 3)), float64(MEM_STATS.TotalAlloc)/(math.Pow(1024, 3)),
			float64(MEM_STATS.Sys)/(math.Pow(1024, 3)), float64(MEM_STATS.HeapAlloc)/(math.Pow(1024, 3)),
			float64(MEM_STATS.HeapSys)/(math.Pow(1024, 3)))
	}
}

//------------------------
// Printing Alignment info
//------------------------

func PrintLoopTraceInfo(loop_num int, mess string) {
	if PRINT_ALIGN_TRACE_INFO {
		fmt.Println("Loop\t", loop_num, "\t", mess)
	}
}

func PrintSeedTraceInfo(mess string, e_pos, s_pos int, read []byte) {
	if PRINT_ALIGN_TRACE_INFO {
		fmt.Println(mess+" has seed\t", e_pos, "\t", s_pos, "\t", string(read[e_pos:s_pos+1]))
	}
}

func PrintPairedSeedInfo(mess string, match_pos_r1, match_pos_r2 int) {
	if PRINT_ALIGN_TRACE_INFO {
		fmt.Println(mess, "\t", match_pos_r1, "\t", match_pos_r2)
	}
}

func PrintExtendTraceInfo(mess string, match []byte, e_pos, s_pos, match_num int, match_pos []int) {
	if PRINT_ALIGN_TRACE_INFO {
		fmt.Println(mess, " extend\t", string(match), "\t", e_pos, "\t", s_pos, "\t", match_num)
		fmt.Print(mess, " match pos\t")
		for i := 0; i < match_num; i++ {
			fmt.Print(match_pos[i], "\t")
		}
		fmt.Println()
	}
}

func PrintMatchTraceInfo(pos, left_most_pos int, dis float64, left_var_pos []int, read []byte) {
	if PRINT_ALIGN_TRACE_INFO {
		fmt.Print("Match\t", pos, "\t", dis, "\t", left_most_pos, "\t", string(read), "\t")
		for _, pos := range left_var_pos {
			fmt.Print(pos, "\t")
		}
		fmt.Println()
	}
}

/*--------------------------
//Printing variant calling info
---------------------------*/

func PrintComparedReadRef(l_read_flank, l_ref_flank, r_read_flank, r_ref_flank []byte) {
	if PRINT_VAR_CALL_INFO {
		fmt.Println("l_read_flank", string(l_read_flank))
		fmt.Println("l_ref_flank", string(l_ref_flank))
		fmt.Println("r_read_flank", string(r_read_flank))
		fmt.Println("r_ref_flank", string(r_ref_flank))
	}
}

func PrintRefPosMap(l_ref_pos_map, r_ref_pos_map []int) {
	if PRINT_VAR_CALL_INFO {
		fmt.Println("l_ref_pos_map", l_ref_pos_map)
		fmt.Println("r_ref_pos_map", r_ref_pos_map)
	}
}

func PrintGetVariants(p_prob, m_prob1, m_prob2 float64, vars1, vars2 []VarInfo) {
	if PRINT_VAR_CALL_INFO {
		fmt.Println("dis to get vars (1st-end, 2nd-end)", p_prob, m_prob1, m_prob2)
		fmt.Println("1st-end variants")
		for _, s := range vars1 {
			fmt.Println(string(s.Bases), string(s.BQual))
		}
		fmt.Println("2nd-end variants")
		for _, s := range vars2 {
			fmt.Println(string(s.Bases), string(s.BQual))
		}
	}
}

/*-------------------------------
//Printing Dist calculation info
-------------------------------*/

func PrintEditDisInput(mess string, str_val ...[]byte) {
	if PRINT_EDIT_DIST_INFO {
		fmt.Println(mess)
		if str_val != nil {
			for _, v := range str_val {
				fmt.Println(string(v))
			}
		}
	}
}

func PrintEditDisMat(mess string, D [][]float64, m, n int, read, ref []byte) {
	if PRINT_EDIT_DIST_MAT_INFO {
		fmt.Println(mess)
		fmt.Print("\t\t")
		for j := 1; j <= n; j++ {
			fmt.Print(string(ref[j-1]), "\t")
		}
		fmt.Println()
		fmt.Print("\t")
		for j := 0; j <= n; j++ {
			fmt.Print(D[0][j], "\t")
		}
		fmt.Println()
		for i := 1; i <= m; i++ {
			fmt.Print(string(read[i-1]) + "\t")
			for j := 0; j <= n; j++ {
				fmt.Print(D[i][j], "\t")
			}
			fmt.Println()
		}
	}
}

/*
BT[i][j][0]: direction, can be 0: diagonal arrow (back to i-1,j-1), 1: up arrow (back to i-1,j), 2: left arrow (back to i,j-1)
BT[i][j][1]: matrix, can be 0: matrix for D, 1: matrix for IS, 2: matrix for IT
BT[i][j][2]: number of shift (equal to length of called variant) at known variant loc, can be any integer number, for example 5 means back to i-5,j-1
*/
func PrintEditTraceMat(mess string, BT [][][]int, m, n int) {
	if PRINT_EDIT_DIST_MAT_INFO {
		fmt.Println(mess)
		for i := 0; i <= m; i++ {
			for j := 0; j <= n; j++ {
				fmt.Print("[")
				if BT[i][j][0] == 0 {
					fmt.Print("d ")
				} else if BT[i][j][0] == 1 {
					fmt.Print("u ")
				} else if BT[i][j][0] == 2 {
					fmt.Print("l ")
				} else {
					fmt.Print("n ")
				}
				if BT[i][j][1] == 0 {
					fmt.Print("D ")
				} else if BT[i][j][1] == 1 {
					fmt.Print("S ")
				} else if BT[i][j][1] == 2 {
					fmt.Print("T ")
				} else {
					fmt.Print("N ")
				}
				fmt.Print(BT[i][j][2], "]\t")
			}
			fmt.Println()
		}
	}
}

func PrintDisInfo(mess string, i, j int, d float64) {
	if PRINT_EDIT_DIST_INFO {
		fmt.Println(mess, i, j, d)
	}
}

func GetEditTrace(mess string, i, j int, read, ref byte) {
	if PRINT_EDIT_DIST_INFO {
		fmt.Println(mess, i, j, string(read), string(ref))
	}
}

func GetEditTraceKnownLoc(mess string, i, j int, read []byte, ref byte) {
	if PRINT_EDIT_DIST_INFO {
		fmt.Println(mess, i, j, string(read), string(ref))
	}
}

func PrintEditAlignInfo(mess string, aligned_read, aligned_qual, aligned_ref []byte) {
	if PRINT_EDIT_DIST_INFO {
		fmt.Println(mess)
		fmt.Println(string(aligned_read))
		fmt.Println(string(aligned_qual))
		fmt.Println(string(aligned_ref))
	}
}

func PrintVarInfo(mess string, var_pos []int, var_val, var_qlt [][]byte) {
	if PRINT_EDIT_DIST_INFO {
		fmt.Println(mess)
		for i := 0; i < len(var_pos); i++ {
			fmt.Println(var_pos[i], string(var_val[i]), string(var_qlt[i]))
		}
	}
}

/*----------------------------------------
//Global variable for tp, fp, fn profiling
----------------------------------------*/

var (
	ALIGN_READ_INFO_CHAN    = make(chan Align_trace_info)
	ALIGN_READ_INFO_ARR     = make([]Align_trace_info, 0)
	VAR_TRACE_INFO_MAP      = make(map[uint32][]Var_trace_info)
	NO_ALIGN_READ_INFO_CHAN = make(chan Align_trace_info)
	NO_ALIGN_READ_INFO_ARR  = make([]Align_trace_info, 0)
	//TRUE_VAR_COMP           = LoadTrueVar("/data/nsvo/test-data/GRCh37_chr1/refs/ref_alt_mutant/mutant_0.3300/variant_comp.txt")
	//TRUE_VAR_PART           = LoadTrueVar("/data/nsvo/test-data/GRCh37_chr1/refs/ref_alt_mutant/mutant_0.3300/variant_part.txt")
	//TRUE_VAR_NONE           = LoadTrueVar("/data/nsvo/test-data/GRCh37_chr1/refs/ref_alt_mutant/mutant_0.3300/variant_none.txt")
	//QUAL_THRES              = 25.0
)

type Align_trace_info struct {
	read1, read2               []byte
	read_info1, read_info2     []byte
	l_align_pos1, l_align_pos2 int
	r_align_pos1, r_align_pos2 int
	align_dis1, align_dis2     int
	var_pos1, var_pos2         []uint32
	var_base1, var_base2       [][]byte
	var_baseq1, var_baseq2     [][]byte
}

type Var_trace_info struct {
	var_base, var_baseq    []byte
	align_pos1, align_pos2 int
	align_dis1, align_dis2 int
	read_info1, read_info2 []byte
	end_from               byte
}

//Reading align trace info from channel and store them
func GetAlignReadInfo() {
	if GET_ALIGN_READ_INFO {
		var i int
		var var_pos uint32
		for at := range ALIGN_READ_INFO_CHAN {
			ALIGN_READ_INFO_ARR = append(ALIGN_READ_INFO_ARR, at)
			if len(at.var_pos1) > 0 {
				for i, var_pos = range at.var_pos1 {
					var var_info Var_trace_info
					if len(at.var_base1[i]) == 0 {
						var_info = Var_trace_info{[]byte{'.'}, []byte{'I'}, at.l_align_pos1, at.l_align_pos2,
							at.align_dis1, at.align_dis2, at.read_info1, at.read_info2, '1'}
					} else {
						var_info = Var_trace_info{at.var_base1[i], at.var_baseq1[i], at.l_align_pos1, at.l_align_pos2,
							at.align_dis1, at.align_dis2, at.read_info1, at.read_info2, '1'}
					}
					VAR_TRACE_INFO_MAP[var_pos] = append(VAR_TRACE_INFO_MAP[var_pos], var_info)
				}
			}
			if len(at.var_pos2) > 0 {
				for i, var_pos = range at.var_pos2 {
					var var_info Var_trace_info
					if len(at.var_base2[i]) == 0 {
						var_info = Var_trace_info{[]byte{'.'}, []byte{'I'}, at.l_align_pos1, at.l_align_pos2,
							at.align_dis1, at.align_dis2, at.read_info1, at.read_info2, '2'}
					} else {
						var_info = Var_trace_info{at.var_base2[i], at.var_baseq2[i], at.l_align_pos1, at.l_align_pos2,
							at.align_dis1, at.align_dis2, at.read_info1, at.read_info2, '2'}
					}
					VAR_TRACE_INFO_MAP[var_pos] = append(VAR_TRACE_INFO_MAP[var_pos], var_info)
				}
			}
		}
		log.Printf("Number of aligned reads:\t%d", len(ALIGN_READ_INFO_ARR))
	}
}

//Reading noalign reads and related info from channel and store them
func GetNoAlignReadInfo() {
	if GET_NO_ALIGN_READ_INFO {
		for at := range NO_ALIGN_READ_INFO_CHAN {
			NO_ALIGN_READ_INFO_ARR = append(NO_ALIGN_READ_INFO_ARR, at)
		}
		log.Printf("Number of no-aligned reads:\t%d", len(NO_ALIGN_READ_INFO_ARR))
	}
}

//Processing noalign reads and related info
func ProcessNoAlignReadInfo(var_call map[uint32]map[string]float64) {
	if PRINT_NA {
		fmt.Println("Processing noaligned read info...")
		file, _ := os.Create(INPUT_INFO.Var_call_file + ".noalign")
		defer file.Close()
		for _, at := range NO_ALIGN_READ_INFO_ARR {
			file.WriteString(string(at.read_info1) + "\t" + string(at.read_info2) + "\n")
		}
		fmt.Println("Finish processing noaligned read info.")
	}
}

/*
 TP-FP info: Log file format:
  var_pos  true_var  var_base  var_qual end_from  align_dis1  align_dis2  align_pos_diff  l_align_pos1  l_align_pos2  true_pos_diff  true_pos1  true_pos2  read_id
  ...
 where:
  values in one line are corresponding to one variant call
  value is "None" if not exist
  var_pos is consecutive in increasing order
  end_from is the read end (1 or 2) in which Vars are called
  var_qual is probability of base being wrong (converted from encoded ASCII char in FASTQ format)

 Log file names:
  "tp_var_comp", "fp_var_comp", "tp_indel_comp", "fp_indel_comp", "tp_var_part", "fp_var_part",
  "tp_indel_part", "fp_indel_part", "tp_var_none", "fp_var_none", "tp_indel_none", "fp_indel_none",
  "fp_var_other", "fp_indel_other"
 where:
  tp: true positives vars
  fp: false positive vars
  variant: variant calls
  indel: indel calls
  comp: info at complete knowledge locations
  part: info at partial knowledge locations
  none: info at no knowledge locations
  other: info at other locations (totaly false positives)
*/
//--------------------------------------------------------------------------------------------------

/*
//Processing TP, FP variants info
func ProcessTPFPVarInfo(var_call map[uint32]map[string]float64) {
	if PRINT_TPFP {
		fmt.Println("Processing TP, FP Var info...")
		files := make([]*os.File, 14)
		file_names := []string{"tp_var_comp", "fp_var_comp", "tp_indel_comp", "fp_indel_comp",
			"tp_var_part", "fp_var_part", "tp_indel_part", "fp_indel_part",
			"tp_var_none", "fp_var_none", "tp_indel_none", "fp_indel_none",
			"fp_var_other", "fp_indel_other"}
		for i, file_name := range file_names {
			files[i], _ = os.Create(INPUT_INFO.Var_call_file + "." + file_name)
			defer files[i].Close()
		}

		var var_pos uint32
		var pos int
		var st Var_trace_info

		Var_Pos := make([]int, 0, len(VAR_TRACE_INFO_MAP))
		for var_pos, _ = range VAR_TRACE_INFO_MAP {
			Var_Pos = append(Var_Pos, int(var_pos))
		}
		sort.Ints(Var_Pos)

		prob_thres := 1 - math.Pow(10, -QUAL_THRES/10)
		var max_prob, var_prob float64
		for _, pos = range Var_Pos {
			var_pos = uint32(pos)
			for _, st = range VAR_TRACE_INFO_MAP[var_pos] {
				max_prob = 0.0
				for _, var_prob = range var_call[var_pos] {
					if max_prob < var_prob {
						max_prob = var_prob
					}
				}
				if max_prob >= prob_thres {
					OutputTPFPVarInfo(files, st, var_pos)
				}
			}
		}
		fmt.Println("Finish processing TP, FP Var info.")
	}
}

//Outputing TP, FP variants info to proper files
func OutputTPFPVarInfo(files []*os.File, st Var_trace_info, var_pos uint32) {
	var ok bool
	var val []byte
	if val, ok = TRUE_VAR_COMP[int(var_pos)]; ok {
		if len(st.var_base) == 1 {
			if st.var_base[0] == val[0] {
				WriteTPFPVarInfo(files[0], st, var_pos, val)
			} else {
				WriteTPFPVarInfo(files[1], st, var_pos, val)
			}
		} else {
			if bytes.Equal(st.var_base, val) {
				WriteTPFPVarInfo(files[2], st, var_pos, val)
			} else {
				WriteTPFPVarInfo(files[3], st, var_pos, val)
			}
		}
	} else if val, ok = TRUE_VAR_PART[int(var_pos)]; ok {
		if len(st.var_base) == 1 {
			if st.var_base[0] == val[0] {
				WriteTPFPVarInfo(files[4], st, var_pos, val)
			} else {
				WriteTPFPVarInfo(files[5], st, var_pos, val)
			}
		} else {
			if bytes.Equal(st.var_base, val) {
				WriteTPFPVarInfo(files[6], st, var_pos, val)
			} else {
				WriteTPFPVarInfo(files[7], st, var_pos, val)
			}
		}
	} else if val, ok = TRUE_VAR_NONE[int(var_pos)]; ok {
		if len(st.var_base) == 1 {
			if st.var_base[0] == val[0] {
				WriteTPFPVarInfo(files[8], st, var_pos, val)
			} else {
				WriteTPFPVarInfo(files[9], st, var_pos, val)
			}
		} else {
			if bytes.Equal(st.var_base, val) {
				WriteTPFPVarInfo(files[10], st, var_pos, val)
			} else {
				WriteTPFPVarInfo(files[11], st, var_pos, val)
			}
		}
	} else {
		if len(st.var_base) == 1 {
			WriteTPFPVarInfo(files[12], st, var_pos, val)
		} else {
			WriteTPFPVarInfo(files[13], st, var_pos, val)
		}
	}
}

//Writing TP, FP variants info to files
func WriteTPFPVarInfo(file *os.File, st Var_trace_info, var_pos uint32, true_var []byte) {
	file.WriteString(strconv.Itoa(int(var_pos)) + "\t" + string(true_var) + "\t" + string(st.var_base) + "\t")
	file.WriteString(strconv.FormatFloat(QualtoProb(st.var_baseq[0]), 'f', 5, 32) + "\t")

	file.WriteString(string(st.end_from) + "\t" + strconv.Itoa(st.align_dis1) + "\t" + strconv.Itoa(st.align_dis2) + "\t")

	if st.align_pos1 != 0 && st.align_pos2 != 0 {
		file.WriteString(strconv.Itoa(st.align_pos1-st.align_pos2) + "\t")
	} else {
		file.WriteString("None\t")
	}
	file.WriteString(strconv.Itoa(st.align_pos1) + "\t" + strconv.Itoa(st.align_pos2) + "\t")

	tokens := bytes.Split(st.read_info1, []byte{'_'})
	if len(tokens) >= 11 {
		true_pos1, err1 := strconv.ParseInt(string(tokens[2]), 10, 64)
		true_pos2, err2 := strconv.ParseInt(string(tokens[3]), 10, 64)
		if err1 == nil && err2 == nil {
			file.WriteString(strconv.FormatInt(true_pos1-true_pos2, 10) + "\t" + strconv.FormatInt(true_pos1, 10) +
				"\t" + strconv.FormatInt(true_pos2, 10) + "\t")
		} else {
			file.WriteString("None\tNone\tNone\t")
		}
		file.WriteString(string(tokens[10]) + "\n")
	} else {
		file.WriteString("None\tNone\tNone\tNone\n")
	}
	file.Sync()
}
*/
/*
 FN info: Log file format:
  fn_pos  fn_var  ref_base var_base  var_prob  align_pos_diff  true_align_pos_diff  align_pos1  align_pos2  true_pos1  true_pos2  align_dis1  align_dis2  read_id  align_base1  base_prob1  align_base2 base_prob2 ...
  ...
 where:
  values in one line are corresponding to one variant call
  value is "None" if not exist
  fn_pos is consecutive in increasing order
  var_prob, (base_prob) is probability of variant (base) being wrong (converted from encoded ASCII char in FASTQ format in case of base_prob)

 Log file names:
	"fn_snp_align_none", "fn_indel_align_none", "fn_snp_misalign_none", "fn_indel_misalign_none", "fn_snp_noalign_none", "fn_indel_noalign_none",
	"fn_snp_align_part", "fn_indel_align_part", "fn_snp_misalign_part", "fn_indel_misalign_part", "fn_snp_noalign_part", "fn_indel_noalign_part",
	"fn_snp_align_comp", "fn_indel_align_comp", "fn_snp_misalign_comp", "fn_indel_misalign_comp", "fn_snp_noalign_comp", "fn_indel_noalign_comp"
 where:
  fn: false negatives
  var: variant calls
  indel: indel calls
  align: fn pos covered by some reads and there is one base called at the pos
  misalign: fn pos covered by some reads but there is no base called at the pos
  noalign: fn pos not covered by any reads
  comp: info at complete knowledge locations
  part: info at partial knowledge locations
  none: info at no knowledge locations
  other: info at other locations (totaly false positives)
*/
//--------------------------------------------------------------------------------------------------
/*
type Var_Call struct {
	Bases string
	Prob  float64
}

var (
	FN_VAR_COMP   = make(map[int][]byte)
	FN_VAR_PART   = make(map[int][]byte)
	FN_VAR_NONE   = make(map[int][]byte)
	FN_VAR_CALL   = make(map[int]Var_Call)
	TPFP_VAR_CALL = make(map[int]bool)
)

//Processing FN variants and realted info
func ProcessFNVarInfo(var_call map[uint32]map[string]float64) {
	if PRINT_FN {
		fmt.Println("Processing FN Var info...")
		files := make([]*os.File, 18)
		file_names := []string{"fn_snp_align_none", "fn_indel_align_none", "fn_snp_misalign_none", "fn_indel_misalign_none", "fn_snp_noalign_none", "fn_indel_noalign_none",
			"fn_snp_align_part", "fn_indel_align_part", "fn_snp_misalign_part", "fn_indel_misalign_part", "fn_snp_noalign_part", "fn_indel_noalign_part",
			"fn_snp_align_comp", "fn_indel_align_comp", "fn_snp_misalign_comp", "fn_indel_misalign_comp", "fn_snp_noalign_comp", "fn_indel_noalign_comp"}
		for i, file_name := range file_names {
			files[i], _ = os.Create(INPUT_INFO.Var_call_file + "." + file_name)
			defer files[i].Close()
		}

		var var_pos int
		var var_bases []byte
		var ok bool

		prob_thres := 1 - math.Pow(10, -QUAL_THRES/10)
		var max_var_str, var_str string
		var max_prob, var_prob float64
		for pos, _ := range VAR_TRACE_INFO_MAP {
			max_prob = 0.0
			for var_str, var_prob = range var_call[pos] {
				if max_prob < var_prob {
					max_prob = var_prob
					max_var_str = var_str
				}
			}
			if max_prob >= prob_thres {
				TPFP_VAR_CALL[int(pos)] = true
			} else {
				FN_VAR_CALL[int(pos)] = Var_Call{max_var_str, max_prob}
			}
		}
		fmt.Println("# Var calls ", len(TPFP_VAR_CALL))

		for var_pos, var_bases = range TRUE_VAR_NONE {
			if _, ok = TPFP_VAR_CALL[var_pos]; !ok {
				FN_VAR_NONE[var_pos] = var_bases
			}
		}
		None_Pos := make([]int, 0)
		for var_pos, _ = range FN_VAR_NONE {
			None_Pos = append(None_Pos, var_pos)
		}
		sort.Ints(None_Pos)
		fmt.Println("# FN_VAR_NONE:", len(None_Pos))
		fmt.Println("Outputing FN_VAR_NONE:")

		OutputFNVarInfo(files[0:6], None_Pos, FN_VAR_NONE)

		for var_pos, var_bases = range TRUE_VAR_PART {
			if _, ok = TPFP_VAR_CALL[var_pos]; !ok {
				FN_VAR_PART[var_pos] = var_bases
			}
		}
		Part_Pos := make([]int, 0)
		for var_pos, _ = range FN_VAR_PART {
			Part_Pos = append(Part_Pos, var_pos)
		}
		sort.Ints(Part_Pos)
		fmt.Println("# FN_VAR_PART:", len(Part_Pos))
		fmt.Println("Outputing FN_VAR_PART:")

		OutputFNVarInfo(files[6:12], Part_Pos, FN_VAR_PART)

		for var_pos, var_bases = range TRUE_VAR_COMP {
			if _, ok = TPFP_VAR_CALL[var_pos]; !ok {
				FN_VAR_COMP[var_pos] = var_bases
			}
		}
		Comp_Pos := make([]int, 0)
		for var_pos, _ = range FN_VAR_COMP {
			Comp_Pos = append(Comp_Pos, var_pos)
		}
		sort.Ints(Comp_Pos)
		fmt.Println("# FN_VAR_COMP:", len(Comp_Pos))
		fmt.Println("Outputing FN_VAR_COMP:")
		OutputFNVarInfo(files[12:18], Comp_Pos, FN_VAR_COMP)
		fmt.Println("Finish processing FN Var info.")
	}
}

//Outputing FN variants info to proper files
func OutputFNVarInfo(files []*os.File, FN_Pos []int, FN_Var map[int][]byte) {
	var at Align_trace_info
	var true_var []byte
	var i, pos int
	var has_align_read bool
	start_time := time.Now()
	process_time := time.Since(start_time)
	for i, pos = range FN_Pos {
		start_time = time.Now()
		true_var = FN_Var[pos]
		has_align_read = false
		for _, at = range ALIGN_READ_INFO_ARR {
			if at.l_align_pos1 <= pos && at.r_align_pos1 >= pos {
				has_align_read = true
				var_call, has_align_base := FN_VAR_CALL[pos]
				if has_align_base {
					if len(true_var) == 1 {
						WriteFNVarInfo(files[0], at, pos, true_var, var_call)
					} else {
						WriteFNVarInfo(files[1], at, pos, true_var, var_call)
					}
				} else {
					if len(true_var) == 1 {
						WriteFNVarInfo(files[2], at, pos, true_var, var_call)
					} else {
						WriteFNVarInfo(files[3], at, pos, true_var, var_call)
					}
				}
			}
			if at.l_align_pos2 <= pos && at.r_align_pos2 >= pos {
				has_align_read = true
				var_call, has_align_base := FN_VAR_CALL[pos]
				if has_align_base {
					if len(true_var) == 1 {
						WriteFNVarInfo(files[0], at, pos, true_var, var_call)
					} else {
						WriteFNVarInfo(files[1], at, pos, true_var, var_call)
					}
				} else {
					if len(true_var) == 1 {
						WriteFNVarInfo(files[2], at, pos, true_var, var_call)
					} else {
						WriteFNVarInfo(files[3], at, pos, true_var, var_call)
					}
				}
			}
		}
		if !has_align_read {
			if len(true_var) == 1 {
				files[4].WriteString(strconv.Itoa(pos) + "\t" + string(true_var) + "\n")
			} else {
				files[5].WriteString(strconv.Itoa(pos) + "\t" + string(true_var) + "\n")
			}
		}
		process_time = time.Since(start_time)
		fmt.Println("time for processing each FN element:\t", i, process_time)
	}
}

//Writing to file FN variants with all reads aligned to FN locations
func WriteFNVarInfo(file *os.File, at Align_trace_info, pos int, true_var []byte, var_call Var_Call) {

	file.WriteString(strconv.Itoa(pos) + "\t" + string(true_var) + "\t")
	if len(var_call.Bases) != 0 {
		file.WriteString(var_call.Bases + "\t" + strconv.FormatFloat(var_call.Prob, 'f', 5, 32) + "\t")
	}
	if at.l_align_pos1 != 0 && at.l_align_pos2 != 0 {
		file.WriteString(strconv.Itoa(at.l_align_pos1-at.l_align_pos2) + "\t")
	} else {
		file.WriteString("None\t")
	}
	file.WriteString(strconv.Itoa(at.l_align_pos1) + "\t" + strconv.Itoa(at.l_align_pos2) + "\t")
	file.WriteString(strconv.Itoa(at.r_align_pos1) + "\t" + strconv.Itoa(at.r_align_pos2) + "\t")
	file.WriteString(strconv.Itoa(at.align_dis1) + "\t" + strconv.Itoa(at.align_dis2) + "\t")

	tokens := bytes.Split(at.read_info1, []byte{'_'})
	if len(tokens) >= 11 {
		true_pos1, err1 := strconv.ParseInt(string(tokens[2]), 10, 64)
		true_pos2, err2 := strconv.ParseInt(string(tokens[3]), 10, 64)
		if err1 == nil && err2 == nil {
			file.WriteString(strconv.FormatInt(true_pos1-true_pos2, 10) + "\t" + strconv.FormatInt(true_pos1, 10) +
				"\t" + strconv.FormatInt(true_pos2, 10) + "\t")
		} else {
			file.WriteString("None\tNone\tNone\t")
		}
		file.WriteString(string(tokens[10]) + "\t")
	} else {
		file.WriteString("None\tNone\tNone\tNone\t")
	}

	var var_info Var_trace_info
	if var_arr, ok := VAR_TRACE_INFO_MAP[uint32(pos)]; ok {
		for _, var_info = range var_arr {
			file.WriteString(strconv.Itoa(pos) + "\t" + string(var_info.var_base) + "\t" +
				strconv.FormatFloat(QualtoProb(var_info.var_baseq[0]), 'f', 5, 32) + "\t")
		}
	}
	file.WriteString("\n")
	file.Sync()
}

//Loading true variants which are used to generate the mutant genome
func LoadTrueVar(file_name string) map[int][]byte {
	barr := make(map[int][]byte)

	f, err := os.Open(file_name)
	if err != nil {
		fmt.Printf("%v\n", err)
		os.Exit(1)
	}
	defer f.Close()
	br := bufio.NewReader(f)
	for {
		line, err := br.ReadString('\n')
		if err != nil {
			break
		}
		sline := string(line[:len(line)-1])
		split := strings.Split(sline, "\t")
		k, _ := strconv.ParseInt(split[0], 10, 64)
		b := make([]byte, len(split[1]))
		copy(b, split[1])
		barr[int(k)] = b
	}
	return barr
}

/*---------------
Utility funtions
---------------*/
/*
//QualtoProb converts base qualities decoded by ASCII codes to probabilities
func QualtoProb(e byte) float64 {
	return math.Pow(10, -(float64(e)-33)/10.0)
}

//ProbtoQual converts probabilities to phred-scale quality scores
func ProbtoQual(p float64) float32 {
	return float32(-10 * math.Log10(1-p))
}
*/
