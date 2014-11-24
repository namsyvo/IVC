//---------------------------------------------------------------------------------------------------
// Calling SNPs based on read-multigenome alignment.
// Copyright 2014 Nam Sy Vo
//---------------------------------------------------------------------------------------------------

package isc

import (
	"fmt"
	"bufio"
	"math"
	"os"
	"strconv"
	"bytes"
	"strings"
	"runtime"
	"log"
	"sort"
)

//Global variable for turnning on/off info profiling
var (
	PRINT_PROCESS_MEM = true
	PRINT_MEM_STATS = false
	GET_ALIGN_TRACE_INFO = true
	PRINT_ALIGN_TRACE_INFO = false
	PRINT_TPFP = false
	PRINT_FN = false
	GET_MIS_ALIGN_TRACE_INFO = true
	PRINT_MA = false
)

//Global variable for memory profiling
var (
	MEM_STATS = new(runtime.MemStats)
)

//Printing memory information
func PrintProcessMem(mesg string) {
	if PRINT_PROCESS_MEM {
		runtime.ReadMemStats(MEM_STATS)
		log.Printf(mesg + "\t%d\t%d\t%d\t%d\t%d\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f",
			MEM_STATS.Alloc, MEM_STATS.TotalAlloc, MEM_STATS.Sys, MEM_STATS.HeapAlloc, MEM_STATS.HeapSys,
			float64(MEM_STATS.Alloc)/(math.Pow(1024, 3)), float64(MEM_STATS.TotalAlloc)/(math.Pow(1024, 3)),
			float64(MEM_STATS.Sys)/(math.Pow(1024, 3)), float64(MEM_STATS.HeapAlloc)/(math.Pow(1024, 3)),
			float64(MEM_STATS.HeapSys)/(math.Pow(1024, 3)))
	}
}

//Printing memory information
func PrintMemStats(mesg string) {
	if PRINT_MEM_STATS {
		runtime.ReadMemStats(MEM_STATS)
		log.Printf(mesg + "\t%d\t%d\t%d\t%d\t%d\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f",
			MEM_STATS.Alloc, MEM_STATS.TotalAlloc, MEM_STATS.Sys, MEM_STATS.HeapAlloc, MEM_STATS.HeapSys,
			float64(MEM_STATS.Alloc)/(math.Pow(1024, 3)), float64(MEM_STATS.TotalAlloc)/(math.Pow(1024, 3)),
			float64(MEM_STATS.Sys)/(math.Pow(1024, 3)), float64(MEM_STATS.HeapAlloc)/(math.Pow(1024, 3)),
			float64(MEM_STATS.HeapSys)/(math.Pow(1024, 3)))
	}
}

//Printing ALignment info
func PrintLoopTraceInfo(loop_num int, read []byte) {
	if PRINT_ALIGN_TRACE_INFO {
		fmt.Println("Loop\t", loop_num, "\t", string(read))
	}
}

func PrintSeedTraceInfo(mess string, e_pos, s_pos int, read []byte) {
	if PRINT_ALIGN_TRACE_INFO {
		fmt.Println(mess + " has seed\t", e_pos, "\t", s_pos, "\t", string(read))
	}
}

func PrintExtendTraceInfo(mess string, match []byte, e_pos, s_pos, match_num int, match_pos []int) {
	if PRINT_ALIGN_TRACE_INFO {
		fmt.Println(mess, " read\t", string(match), "\t", e_pos, "\t", s_pos, "\t", match_num)
		fmt.Print(mess, " match pos\t")
		for i := 0; i < match_num; i++ {
			fmt.Print(match_pos[i], "\t")
		}
		fmt.Println()
	}
}

func PrintMatchTraceInfo(i, pos, dis, left_most_pos int, left_snp_pos []int, read []byte) {
	if PRINT_ALIGN_TRACE_INFO {
		fmt.Println("Match\t", i, "\t", pos, "\t", dis, "\t", left_most_pos, "\t", string(read), "\t")
		for pos := range left_snp_pos {
			fmt.Print(pos, "\t")
		}
		fmt.Println()
	}
}

//Global variable for tp, fp, fn profiling
var (
    ALIGN_TRACE_INFO_CHAN = make(chan Align_trace_info)
	ALIGN_TRACE_INFO_ARR = make([]Align_trace_info, 0)
	SNP_TRACE_INFO_MAP = make(map[uint32][]SNP_trace_info)
    MIS_ALIGN_TRACE_INFO_CHAN = make(chan Align_trace_info)
	MIS_ALIGN_TRACE_INFO_ARR = make([]Align_trace_info, 0)
    TRUE_VAR_COMP = LoadTrueVar("/data/nsvo/test-data/GRCh37_chr1/refs/mutate-0.3300/variant_comp.txt")
    TRUE_VAR_PART = LoadTrueVar("/data/nsvo/test-data/GRCh37_chr1/refs/mutate-0.3300/variant_part.txt")
    TRUE_VAR_NONE = LoadTrueVar("/data/nsvo/test-data/GRCh37_chr1/refs/mutate-0.3300/variant_none.txt")
	QUAL_THRES = 25.0
)

type Align_trace_info struct {
	read1, read2 []byte
	read_info1, read_info2 []byte
	align_pos1, align_pos2 int
	align_right_pos1, align_right_pos2 int
	align_dis1, align_dis2 int
	snp_pos1, snp_pos2 []uint32
	snp_base1, snp_base2 [][]byte
	snp_baseq1, snp_baseq2 [][]byte
}

type SNP_trace_info struct {
	snp_base, snp_baseq []byte
	align_pos1, align_pos2 int
	align_dis1, align_dis2 int
	read_info1, read_info2 []byte
	end_from byte
}

//Read align trace info from channel and store them
func GetAlignTraceInfo() {
	if GET_ALIGN_TRACE_INFO {
		var i int
		var snp_pos uint32
		for at := range ALIGN_TRACE_INFO_CHAN {
			ALIGN_TRACE_INFO_ARR = append(ALIGN_TRACE_INFO_ARR, at)
			if len(at.snp_pos1) > 0 {
				for i, snp_pos = range at.snp_pos1 {
					var snp SNP_trace_info
					if len(at.snp_base1[i]) == 0 {
						snp = SNP_trace_info{[]byte{'.'}, []byte{'I'}, at.align_pos1, at.align_pos2, 
							at.align_dis1, at.align_dis2, at.read_info1, at.read_info2, '1'}
					} else {
						snp = SNP_trace_info{at.snp_base1[i], at.snp_baseq1[i], at.align_pos1, at.align_pos2, 
							at.align_dis1, at.align_dis2, at.read_info1, at.read_info2, '1'}
					}
					SNP_TRACE_INFO_MAP[snp_pos] = append(SNP_TRACE_INFO_MAP[snp_pos], snp)
				}
			}
			if len(at.snp_pos2) > 0 {
				for i, snp_pos = range at.snp_pos2 {
					var snp SNP_trace_info
					if len(at.snp_base2[i]) == 0 {
						snp = SNP_trace_info{[]byte{'.'}, []byte{'I'}, at.align_pos1, at.align_pos2, 
							at.align_dis1, at.align_dis2, at.read_info1, at.read_info2, '2'}
					} else {
						snp = SNP_trace_info{at.snp_base2[i], at.snp_baseq2[i], at.align_pos1, at.align_pos2, 
							at.align_dis1, at.align_dis2, at.read_info1, at.read_info2, '2'}
					}
					SNP_TRACE_INFO_MAP[snp_pos] = append(SNP_TRACE_INFO_MAP[snp_pos], snp)
				}
			}
		}
		fmt.Println("# of aligned reads:", len(ALIGN_TRACE_INFO_ARR))
	}
}

//Read misalign reads and related info from channel and store them
func GetMisAlignTraceInfo() {
	if GET_MIS_ALIGN_TRACE_INFO {
		for at := range MIS_ALIGN_TRACE_INFO_CHAN {
			MIS_ALIGN_TRACE_INFO_ARR = append(MIS_ALIGN_TRACE_INFO_ARR, at)
		}
		fmt.Println("# of misaligned reads:", len(MIS_ALIGN_TRACE_INFO_ARR))
	}
}

/*
 TP-FP info: Log file format:
  snp_pos  true_snp  snp_base  snp_qual end_from  align_dis1  align_dis2  align_pos_diff  align_pos1  align_pos2  true_pos_diff  true_pos1  true_pos2  read_id
  ...
 where:
  values in one line are corresponding to one variant call
  value is "None" if not exist
  snp_pos is consecutive in increasing order
  end_from is the read end (1 or 2) in which SNPs are called
  snp_qual is probability of base being wrong (converted from encoded ASCII char in FASTQ format)

 Log file names:
  "tp_snp_comp", "fp_snp_comp", "tp_indel_comp", "fp_indel_comp", "tp_snp_part", "fp_snp_part", 
  "tp_indel_part", "fp_indel_part", "tp_snp_none", "fp_snp_none", "tp_indel_none", "fp_indel_none", 
  "fp_snp_other", "fp_indel_other"
 where:
  tp: true positives snps
  fp: false positive snps
  snp: snp calls
  indel: indel calls
  comp: info at complete knowledge locations
  part: info at partial knowledge locations
  none: info at no knowledge locations
  other: info at other locations (totaly false positives)
*/
//--------------------------------------------------------------------------------------------------

//Process TP, FP variants info
func ProcessSNPTPFPInfo(snp_call map[uint32]map[string]float64) {
	if PRINT_TPFP {
		files := make([]*os.File, 14)
		file_names := []string{"tp_snp_comp", "fp_snp_comp", "tp_indel_comp", "fp_indel_comp", 
								"tp_snp_part", "fp_snp_part", "tp_indel_part", "fp_indel_part", 
								"tp_snp_none", "fp_snp_none", "tp_indel_none", "fp_indel_none", 
								"fp_snp_other", "fp_indel_other"}
		for i, file_name := range file_names {
			files[i], _ = os.Create(INPUT_INFO.SNP_call_file + "." + file_name)
			defer files[i].Close()
		}

		var snp_pos uint32
		var pos int
		var st SNP_trace_info
		
		SNP_Pos := make([]int, 0, len(SNP_TRACE_INFO_MAP))
		for snp_pos, _ = range SNP_TRACE_INFO_MAP {
			SNP_Pos = append(SNP_Pos, int(snp_pos))
		}
		sort.Ints(SNP_Pos)
		
		prob_thres := 1 - math.Pow(10, -QUAL_THRES/10)
		var max_prob, snp_prob float64
		for _, pos = range SNP_Pos {
			snp_pos = uint32(pos)
			for _, st = range SNP_TRACE_INFO_MAP[snp_pos] {
				max_prob = 0.0
				for _, snp_prob = range snp_call[snp_pos] {
					if max_prob < snp_prob {
						max_prob = snp_prob
					}
				}
				if max_prob >= prob_thres {
					OutputSNPTPFPInfo(files, st, snp_pos)
				}
			}
		}
	}
}

//Output TP, FP variants info to proper files
func OutputSNPTPFPInfo(files []*os.File, st SNP_trace_info, snp_pos uint32) {
	var ok bool
	var val []byte
	if val, ok = TRUE_VAR_COMP[int(snp_pos)]; ok {
		if len(st.snp_base) == 1 {
			if st.snp_base[0] == val[0] {
				WriteSNPTPFPInfo(files[0], st, snp_pos, val)
			} else {
				WriteSNPTPFPInfo(files[1], st, snp_pos, val)
			}
		} else {
			if bytes.Equal(st.snp_base, val) {
				WriteSNPTPFPInfo(files[2], st, snp_pos, val)
			} else {
				WriteSNPTPFPInfo(files[3], st, snp_pos, val)
			}
		}
	} else if val, ok = TRUE_VAR_PART[int(snp_pos)]; ok {
		if len(st.snp_base) == 1 {
			if st.snp_base[0] == val[0] {
				WriteSNPTPFPInfo(files[4], st, snp_pos, val)
			} else {
				WriteSNPTPFPInfo(files[5], st, snp_pos, val)
			}
		} else {
			if bytes.Equal(st.snp_base, val) {
				WriteSNPTPFPInfo(files[6], st, snp_pos, val)
			} else {
				WriteSNPTPFPInfo(files[7], st, snp_pos, val)
			}
		}
	} else if val, ok = TRUE_VAR_NONE[int(snp_pos)]; ok {
		if len(st.snp_base) == 1 {
			if st.snp_base[0] == val[0] {
				WriteSNPTPFPInfo(files[8], st, snp_pos, val)
			} else {
				WriteSNPTPFPInfo(files[9], st, snp_pos, val)
			}
		} else {
			if bytes.Equal(st.snp_base, val) {
				WriteSNPTPFPInfo(files[10], st, snp_pos, val)
			} else {
				WriteSNPTPFPInfo(files[11], st, snp_pos, val)
			}
		}
	} else {
		if len(st.snp_base) == 1 {
			WriteSNPTPFPInfo(files[12], st, snp_pos, val)
		} else {
			WriteSNPTPFPInfo(files[13], st, snp_pos, val)
		}
	}
}

//Write TP, FP variants info to files
func WriteSNPTPFPInfo(file *os.File, st SNP_trace_info, snp_pos uint32, true_var []byte) {
	file.WriteString(strconv.Itoa(int(snp_pos)) + "\t" + string(true_var) + "\t" + string(st.snp_base) + "\t")
	file.WriteString(strconv.FormatFloat(QualtoProb(st.snp_baseq[0]), 'f', 5, 32) + "\t")
	
	file.WriteString(string(st.end_from) + "\t" + strconv.Itoa(st.align_dis1) + "\t" + strconv.Itoa(st.align_dis2) + "\t")
	
	if st.align_pos1 != 0 && st.align_pos2 != 0 {
		file.WriteString(strconv.Itoa(st.align_pos1 - st.align_pos2) + "\t")
	} else {
		file.WriteString("None\t")
	}
	file.WriteString(strconv.Itoa(st.align_pos1) + "\t" + strconv.Itoa(st.align_pos2) + "\t")
	
	tokens := bytes.Split(st.read_info1, []byte{'_'})
	if len(tokens) >= 11 {
		true_pos1, err1 := strconv.ParseInt(string(tokens[2]), 10, 64)
		true_pos2, err2 := strconv.ParseInt(string(tokens[3]), 10, 64)
		if err1 == nil && err2 == nil {
			file.WriteString(strconv.FormatInt(true_pos1 - true_pos2, 10) + "\t" + strconv.FormatInt(true_pos1, 10) + 
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

/*
 FN info: Log file format:
  fn_pos  fn_snp  ref_base snp_base  snp_prob  align_pos_diff  true_align_pos_diff  align_pos1  align_pos2  true_pos1  true_pos2  align_dis1  align_dis2  read_id  align_base1  base_prob1  align_base2 base_prob2 ...
  ...
 where:
  values in one line are corresponding to one variant call
  value is "None" if not exist
  fn_pos is consecutive in increasing order
  snp_prob, (base_prob) is probability of snp (base) being wrong (converted from encoded ASCII char in FASTQ format in case of base_prob)

 Log file names:
	"fn_snp_align_none", "fn_indel_align_none", "fn_snp_misalign_none", "fn_indel_misalign_none", "fn_snp_noalign_none", "fn_indel_noalign_none", 
	"fn_snp_align_part", "fn_indel_align_part", "fn_snp_misalign_part", "fn_indel_misalign_part", "fn_snp_noalign_part", "fn_indel_noalign_part", 
	"fn_snp_align_comp", "fn_indel_align_comp", "fn_snp_misalign_comp", "fn_indel_misalign_comp", "fn_snp_noalign_comp", "fn_indel_noalign_comp"
 where:
  fn: false negatives
  snp: snp calls
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

type SNP_CALL struct {
	Bases string
	Prob float64
}

var (
	FN_VAR_COMP = make(map[int][]byte)
	FN_VAR_PART = make(map[int][]byte)
	FN_VAR_NONE = make(map[int][]byte)
	FN_VAR_CALL = make(map[int]SNP_CALL)
	TPFP_VAR_CALL = make(map[int]bool)
)

//Process FN variants and realted info
func ProcessSNPFNInfo(snp_call map[uint32]map[string]float64) {
	if PRINT_FN {
		files := make([]*os.File, 18)
		file_names := []string{"fn_snp_align_none", "fn_indel_align_none", "fn_snp_misalign_none", "fn_indel_misalign_none", "fn_snp_noalign_none", "fn_indel_noalign_none", 
								"fn_snp_align_part", "fn_indel_align_part", "fn_snp_misalign_part", "fn_indel_misalign_part", "fn_snp_noalign_part", "fn_indel_noalign_part", 
								"fn_snp_align_comp", "fn_indel_align_comp", "fn_snp_misalign_comp", "fn_indel_misalign_comp", "fn_snp_noalign_comp", "fn_indel_noalign_comp"}
		for i, file_name := range file_names {
			files[i], _ = os.Create(INPUT_INFO.SNP_call_file + "." + file_name)
			defer files[i].Close()
		}

		var snp_pos int
		var snp []byte
		var ok bool

		prob_thres := 1 - math.Pow(10, -QUAL_THRES/10)
		var max_snp_str, snp_str string
		var max_prob, snp_prob float64
		for pos, _ := range SNP_TRACE_INFO_MAP {
		    max_prob = 0.0
	    	for snp_str, snp_prob = range snp_call[pos] {
	        	if max_prob < snp_prob {
	            	max_prob = snp_prob
	            	max_snp_str = snp_str
	        	}
	    	}
			if max_prob >= prob_thres {
				TPFP_VAR_CALL[int(pos)] = true
			} else {
				FN_VAR_CALL[int(pos)] = SNP_CALL{max_snp_str, max_prob}
			}
		}
		fmt.Println("# SNP calls ", len(TPFP_VAR_CALL))

		for snp_pos, snp = range TRUE_VAR_NONE {
			if _, ok = TPFP_VAR_CALL[snp_pos]; !ok {
				FN_VAR_NONE[snp_pos] = snp
			}
		}
		None_Pos := make([]int, 0)
		for snp_pos, _ = range FN_VAR_NONE {
				None_Pos = append(None_Pos, snp_pos)
		}
		sort.Ints(None_Pos)
		fmt.Println("# FN_VAR_NONE:", len(None_Pos))
		fmt.Println("Outputing FN_VAR_NONE:")

		OutputSNPFNInfo(files[0:6], None_Pos, FN_VAR_NONE)

		for snp_pos, snp = range TRUE_VAR_PART {
			if _, ok = TPFP_VAR_CALL[snp_pos]; !ok {
				FN_VAR_PART[snp_pos] = snp
			}
		}
		Part_Pos := make([]int, 0)
		for snp_pos, _ = range FN_VAR_PART {
				Part_Pos = append(Part_Pos, snp_pos)
		}
		sort.Ints(Part_Pos)
		fmt.Println("# FN_VAR_PART:", len(Part_Pos))
		fmt.Println("Outputing FN_VAR_PART:")

		OutputSNPFNInfo(files[6:12], Part_Pos, FN_VAR_PART)

		for snp_pos, snp = range TRUE_VAR_COMP {
			if _, ok = TPFP_VAR_CALL[snp_pos]; !ok {
				FN_VAR_COMP[snp_pos] = snp
			}
		}
		Comp_Pos := make([]int, 0)
		for snp_pos, _ = range FN_VAR_COMP {
				Comp_Pos = append(Comp_Pos, snp_pos)
		}
		sort.Ints(Comp_Pos)
		fmt.Println("# FN_VAR_COMP:", len(Comp_Pos))
		fmt.Println("Outputing FN_VAR_COMP:")
		OutputSNPFNInfo(files[12:18], Comp_Pos, FN_VAR_COMP)
	}
}

//Output FN variants info to proper files
func OutputSNPFNInfo(files []*os.File, FN_Pos []int, FN_Var map[int][]byte) {
	var at Align_trace_info
	var true_var []byte
	var pos int
	var has_align_read bool
	for _, pos = range FN_Pos {
		true_var = FN_Var[pos]
		has_align_read = false
		for _, at = range ALIGN_TRACE_INFO_ARR {
			if at.align_pos1 <= pos && at.align_right_pos1 >= pos {
				has_align_read = true
				snp_call, has_align_base := FN_VAR_CALL[pos]
				if has_align_base {
					if len(true_var) == 1 {
						WriteSNPFNInfo(files[0], at, pos, true_var, snp_call)
					} else {
						WriteSNPFNInfo(files[1], at, pos, true_var, snp_call)
					}
				} else {
					if len(true_var) == 1 {
						WriteSNPFNInfo(files[2], at, pos, true_var, snp_call)
					} else {
						WriteSNPFNInfo(files[3], at, pos, true_var, snp_call)
					}			
				}
			}
			if at.align_pos2 <= pos && at.align_right_pos2 >= pos {
				has_align_read = true
				snp_call, has_align_base := FN_VAR_CALL[pos]
				if has_align_base {
					if len(true_var) == 1 {
						WriteSNPFNInfo(files[0], at, pos, true_var, snp_call)
					} else {
						WriteSNPFNInfo(files[1], at, pos, true_var, snp_call)
					}
				} else {
					if len(true_var) == 1 {
						WriteSNPFNInfo(files[2], at, pos, true_var, snp_call)
					} else {
						WriteSNPFNInfo(files[3], at, pos, true_var, snp_call)
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
	}
}

//Write to file FN variants with all reads aligned to FN locations
func WriteSNPFNInfo(file *os.File, at Align_trace_info, pos int, true_var []byte, snp_call SNP_CALL) {

	file.WriteString(strconv.Itoa(pos) + "\t" + string(true_var) + "\t")
	if len(snp_call.Bases) != 0 {
		file.WriteString(snp_call.Bases + "\t" + strconv.FormatFloat(snp_call.Prob, 'f', 5, 32) + "\t")
	}
	if at.align_pos1 != 0 && at.align_pos2 != 0 {
		file.WriteString(strconv.Itoa(at.align_pos1 - at.align_pos2) + "\t")
	} else {
		file.WriteString("None\t")
	}
	file.WriteString(strconv.Itoa(at.align_pos1) + "\t" + strconv.Itoa(at.align_pos2) + "\t")
	file.WriteString(strconv.Itoa(at.align_right_pos1) + "\t" + strconv.Itoa(at.align_right_pos2) + "\t")
	file.WriteString(strconv.Itoa(at.align_dis1) + "\t" + strconv.Itoa(at.align_dis2) + "\t")

	tokens := bytes.Split(at.read_info1, []byte{'_'})
	if len(tokens) >= 11 {
		true_pos1, err1 := strconv.ParseInt(string(tokens[2]), 10, 64)
		true_pos2, err2 := strconv.ParseInt(string(tokens[3]), 10, 64)
		if err1 == nil && err2 == nil {
			file.WriteString(strconv.FormatInt(true_pos1 - true_pos2, 10) + "\t" + strconv.FormatInt(true_pos1, 10) + 
				"\t" + strconv.FormatInt(true_pos2, 10) + "\t")
		} else {
			file.WriteString("None\tNone\tNone\t")
		}
		file.WriteString(string(tokens[10]) + "\t")
	} else {
		file.WriteString("None\tNone\tNone\tNone\t")
	}

	var snp SNP_trace_info
	if snp_arr, ok := SNP_TRACE_INFO_MAP[uint32(pos)]; ok {
		for _, snp = range snp_arr {
			file.WriteString(strconv.Itoa(pos) + "\t" + string(snp.snp_base) + "\t" + 
				strconv.FormatFloat(QualtoProb(snp.snp_baseq[0]), 'f', 5, 32) + "\t")
		}
	}
	file.WriteString("\n")
	file.Sync()
}

//Process misalign reads and realted info
func ProcessMisAlignInfo(snp_call map[uint32]map[string]float64) {
	if PRINT_MA {
		var pos int
		var snp []byte
		SNP_pos := make([]int, 0)
		Indel_pos := make([]int, 0)
		for pos, snp = range TRUE_VAR_COMP {
			if len(snp) == 0 || len(snp) > 1 {
				Indel_pos = append(Indel_pos, pos)
			} else {
				SNP_pos = append(SNP_pos, pos)
			}
		}
		for pos, snp = range TRUE_VAR_COMP {
			if len(snp) == 0 || len(snp) > 1 {
				Indel_pos = append(Indel_pos, pos)
			} else {
				SNP_pos = append(SNP_pos, pos)
			}
		}
		for pos, snp = range TRUE_VAR_COMP {
			if len(snp) == 0 || len(snp) > 1 {
				Indel_pos = append(Indel_pos, pos)
			} else {
				SNP_pos = append(SNP_pos, pos)
			}
		}

		fmt.Println("# SNP_pos:", len(SNP_pos))
		Sorted_snp_pos := make([]int, 0)
		for _, pos = range SNP_pos {
				Sorted_snp_pos = append(Sorted_snp_pos, pos)
		}
		sort.Ints(Sorted_snp_pos)

		fmt.Println("# Indel_pos:", len(Indel_pos))
		Sorted_indel_pos := make([]int, 0)
		for _, pos = range Indel_pos {
				Sorted_indel_pos = append(Sorted_indel_pos, pos)
		}
		sort.Ints(Sorted_indel_pos)

		file1, _ := os.Create(INPUT_INFO.SNP_call_file + ".misalign_snp")
		file2, _ := os.Create(INPUT_INFO.SNP_call_file + ".misalign_indel")
		defer file1.Close()
		defer file2.Close()

		for _, at := range MIS_ALIGN_TRACE_INFO_ARR {
			if IntervalHasSNP(Sorted_indel_pos, at.align_pos1, at.align_right_pos1) || IntervalHasSNP(Sorted_indel_pos, at.align_pos2, at.align_right_pos2){
				WriteMisALignInfo(file2, Sorted_indel_pos, at)
			} else {
				WriteMisALignInfo(file1, Sorted_snp_pos, at)
			}
		}
	}
}

//Write misalign reads and related info to file
func WriteMisALignInfo(file *os.File, var_pos []int, at Align_trace_info) {
	if at.align_pos1 != 0 && at.align_pos2 != 0 {
		file.WriteString(strconv.Itoa(at.align_pos1 - at.align_pos2) + "\t")
	} else {
		file.WriteString("None\t")
	}
	file.WriteString(strconv.Itoa(at.align_pos1) + "\t" + strconv.Itoa(at.align_pos2) + "\t")
	file.WriteString(strconv.Itoa(at.align_right_pos1) + "\t" + strconv.Itoa(at.align_right_pos2) + "\t")
	file.WriteString(strconv.Itoa(at.align_dis1) + "\t" + strconv.Itoa(at.align_dis2) + "\t")
				
	tokens := bytes.Split(at.read_info1, []byte{'_'})
	if len(tokens) >= 11 {
		true_pos1, err1 := strconv.ParseInt(string(tokens[2]), 10, 64)
		true_pos2, err2 := strconv.ParseInt(string(tokens[3]), 10, 64)
		if err1 == nil && err2 == nil {
			file.WriteString(strconv.FormatInt(true_pos1 - true_pos2, 10) + "\t" + strconv.FormatInt(true_pos1, 10) + 
				"\t" + strconv.FormatInt(true_pos2, 10) + "\t")
		} else {
			file.WriteString("None\tNone\tNone\t")
		}
		file.WriteString(string(tokens[10]) + "\t")
	} else {
		file.WriteString("None\tNone\tNone\tNone\t")
	}
	var pos int
	align_pos := make([]int, 0)
	for _, pos = range var_pos {
		if at.align_pos1 <= pos && at.align_right_pos1 > pos {
			align_pos = append(align_pos, pos)
		}
		if at.align_pos2 <= pos && at.align_right_pos2 > pos {
			align_pos = append(align_pos, pos)
		}
	}
	for _, pos = range align_pos {
		file.WriteString(strconv.Itoa(pos) + "\t")
	}
	file.WriteString("\n")
	file.Sync()
}

//Load true variants which are used to generate the mutant genome
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

//QualtoProb converts base qualities decoded by ASCII codes to probabilities
func QualtoProb(e byte) float64 {
	return math.Pow(10, -(float64(e) - 33)/10.0)
}
//ProbtoQual converts probabilities to phred-scale quality scores
func ProbtoQual(p float64) float32 {
	return float32(-10*math.Log10(1 - p))
}