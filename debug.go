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

	PRINT_EDIT_DIST_INFO     = false
	PRINT_EDIT_DIST_MAT_INFO = false

	PRINT_VAR_CALL_INFO    = false
	PRINT_ALIGN_TRACE_INFO = false

	GET_UNALIGN_INFO   = true
	PRINT_UNALIGN_INFO = false
)

//Global variable for memory profiling
var (
	CPU_FILE  *os.File
	MEM_FILE  *os.File
	MEM_STATS *runtime.MemStats
)

//Printing memory information
func PrintProcessMem(mesg string) {
	if PRINT_PROCESS_MEM {
		PrintMemStats(mesg)
	}
}

func PrintMemStats(mesg string) {
	runtime.ReadMemStats(MEM_STATS)
	log.Printf(mesg+"\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f",
		float64(MEM_STATS.Alloc)/(math.Pow(1024, 3)), float64(MEM_STATS.TotalAlloc)/(math.Pow(1024, 3)),
		float64(MEM_STATS.Sys)/(math.Pow(1024, 3)), float64(MEM_STATS.HeapAlloc)/(math.Pow(1024, 3)),
		float64(MEM_STATS.HeapSys)/(math.Pow(1024, 3)))
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

func PrintGetVariants(mess string, paired_prob, prob1, prob2 float64, vars1, vars2 []*VarInfo) {
	if PRINT_VAR_CALL_INFO {
		fmt.Println(mess)
		fmt.Println("Dis to get vars: paired_prob, prob1 (1st-end), prob2 (2nd-end)", paired_prob, prob1, prob2)
		fmt.Println("1st-end variants")
		for _, s := range vars1 {
			fmt.Println(s.Pos, string(s.Bases), string(s.BQual))
		}
		fmt.Println("2nd-end variants")
		for _, s := range vars2 {
			fmt.Println(s.Pos, string(s.Bases), string(s.BQual))
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

//----------------------------------------
// Process unaligned-reads info
//----------------------------------------
type UnAlignInfo struct {
	read_info1, read_info2 []byte
}

var (
	UNALIGN_INFO_CHAN = make(chan UnAlignInfo)
	UNALIGN_INFO_ARR  = make([]UnAlignInfo, 0)
)

//Reading noalign reads and related info from channel and store them
func GetNoAlignReadInfo() {
	if GET_UNALIGN_INFO {
		for uai := range UNALIGN_INFO_CHAN {
			UNALIGN_INFO_ARR = append(UNALIGN_INFO_ARR, uai)
		}
		log.Printf("Number of no-aligned reads:\t%d", len(UNALIGN_INFO_ARR))
	}
}

//Processing noalign reads and related info
func ProcessNoAlignReadInfo() {
	if PRINT_UNALIGN_INFO {
		fmt.Println("Processing noaligned read info...")
		file, _ := os.Create(PARA_INFO.Var_call_file + ".unalign")
		defer file.Close()
		for _, uai := range UNALIGN_INFO_ARR {
			file.WriteString(string(uai.read_info1) + "\t" + string(uai.read_info2) + "\n")
		}
		fmt.Println("Finish processing noaligned read info.")
	}
}
