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

//QualtoProb converts base qualities decoded by ASCII codes to probabilities
func QualtoProb(e byte) float64 {
	return math.Pow(10, -(float64(e) - 33)/10.0)
}
//ProbtoQual converts probabilities to phred-scale quality scores
func ProbtoQual(p float64) float32 {
	return float32(-10*math.Log10(1 - p))
}

//Global variable for memory profiling
var (
	MEM_STATS = new(runtime.MemStats)
	PRINT_MEM = false
)
//Printing memory information
func PrintMemStats(mesg string) {
	if PRINT_MEM {
		runtime.ReadMemStats(MEM_STATS)
		log.Printf(mesg + "\t%d\t%d\t%d\t%d\t%d\t%.3f%.3f%.3f%.3f%.3f",
			MEM_STATS.Alloc, MEM_STATS.TotalAlloc, MEM_STATS.Sys, MEM_STATS.HeapAlloc, MEM_STATS.HeapSys,
			float64(MEM_STATS.Alloc)/(math.Pow(1024, 3)), float64(MEM_STATS.TotalAlloc)/(math.Pow(1024, 3)),
			float64(MEM_STATS.Sys)/(math.Pow(1024, 3)), float64(MEM_STATS.HeapAlloc)/(math.Pow(1024, 3)),
			float64(MEM_STATS.HeapSys)/(math.Pow(1024, 3)))
	}
}

//Global variable for align info profiling
var (
	PRINT_ALIGN_TRACE_INFO = false
)

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
	SNP_TRACE_INFO_MAP = make(map[uint32][]SNP_trace_info)
	ALIGN_TRACE_INFO_MAP = make(map[int][]Align_trace_info)
    TRUE_VAR_COMP = LoadTrueVar("/data/nsvo/test_data/GRCh37_chr1/refs/mutate-0.3300/variant_comp.txt")
    TRUE_VAR_PART = LoadTrueVar("/data/nsvo/test_data/GRCh37_chr1/refs/mutate-0.3300/variant_part.txt")
    TRUE_VAR_NONE = LoadTrueVar("/data/nsvo/test_data/GRCh37_chr1/refs/mutate-0.3300/variant_none.txt")
	QUAL_THRES = 25.0
	GET_INFO = true
	PRINT_TPFP = false
	PRINT_FN = true
)

type Align_trace_info struct {
	read_info1, read_info2 []byte
	align_pos1, align_pos2 int
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

/*
 Log file format:
  snp_pos  true_snp  snp_base  snp_qual end_from  align_dis1  align_dis2  align_pos_diff  align_pos1  align_pos2  true_pos_diff  true_pos1  true_pos2  read_id
  ...
 where:
  values in one line are corresponding to one variant call
  value is "None" if not exist
  snp_pos is consecutive in increasing order
  end_from is the read end (1 or 2) in which SNPs are called
  snp_qual is probability of base being wrong (converted from encoded ASCII char in FASTQ format)

 Log file names:
  "tp_snp_comp", "fp_snp_comp", "tp_indel_comp", "fp_indel_comp", "tp_snp_part", "fp_snp_part", "tp_indel_part", "fp_indel_part", \
  "tp_snp_none", "fp_snp_none", "tp_indel_none", "fp_indel_none", "fp_snp_other", "fp_indel_other"
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

//Read debug info from channel and store them
func GetAlignTraceInfo() {

	if GET_INFO {
		var i int
		var snp_pos uint32
		var at Align_trace_info
		
		for at = range ALIGN_TRACE_INFO_CHAN {
			ALIGN_TRACE_INFO_MAP[at.align_pos1] = append(ALIGN_TRACE_INFO_MAP[at.align_pos1], at)
			ALIGN_TRACE_INFO_MAP[at.align_pos2] = append(ALIGN_TRACE_INFO_MAP[at.align_pos2], at)
			if len(at.snp_pos1) > 0 {
				for i, snp_pos = range at.snp_pos1 {
					var snp SNP_trace_info
					if len(at.snp_base1[i]) == 0 {
						snp = SNP_trace_info{[]byte{'.'}, []byte{'I'}, at.align_pos1, at.align_pos2, at.align_dis1, at.align_dis2, at.read_info1, at.read_info2, '1'}
					} else {
						snp = SNP_trace_info{at.snp_base1[i], at.snp_baseq1[i], at.align_pos1, at.align_pos2, at.align_dis1, at.align_dis2, at.read_info1, at.read_info2, '1'}
					}
					SNP_TRACE_INFO_MAP[snp_pos] = append(SNP_TRACE_INFO_MAP[snp_pos], snp)
				}
			}
			if len(at.snp_pos2) > 0 {
				for i, snp_pos = range at.snp_pos2 {
					var snp SNP_trace_info
					if len(at.snp_base2[i]) == 0 {
						snp = SNP_trace_info{[]byte{'.'}, []byte{'I'}, at.align_pos1, at.align_pos2, at.align_dis1, at.align_dis2, at.read_info1, at.read_info2, '2'}
					} else {
						snp = SNP_trace_info{at.snp_base2[i], at.snp_baseq2[i], at.align_pos1, at.align_pos2, at.align_dis1, at.align_dis2, at.read_info1, at.read_info2, '2'}
					}
					SNP_TRACE_INFO_MAP[snp_pos] = append(SNP_TRACE_INFO_MAP[snp_pos], snp)
				}
			}
		}
	}
}

//Process debug info
func ProcessSNPTPFPInfo(snp_call map[uint32]map[string]float64) {
	if PRINT_TPFP {
		files := make([]*os.File, 14)
		file_names := []string{"tp_snp_comp", "fp_snp_comp", "tp_indel_comp", "fp_indel_comp", "tp_snp_part", "fp_snp_part", "tp_indel_part", "fp_indel_part", "tp_snp_none", "fp_snp_none", "tp_indel_none", "fp_indel_none", "fp_snp_other", "fp_indel_other"}
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
		
		for _, pos = range SNP_Pos {
			snp_pos = uint32(pos)
			for _, st = range SNP_TRACE_INFO_MAP[snp_pos] {
				OutputSNPTPFPInfo(snp_call[snp_pos], files, st, snp_pos)
			}
		}
	}
}

//Output debug info to proper files
func OutputSNPTPFPInfo(snp_call map[string]float64, files []*os.File, st SNP_trace_info, snp_pos uint32) {
    prob := 0.0
    for _, snp_prob := range snp_call {
        if prob < snp_prob {
            prob = snp_prob
        }
    }
	qual := -10 * math.Log10(1 - prob)
	if qual >= QUAL_THRES {
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
}

//Write debug info to files
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
			file.WriteString(strconv.FormatInt(true_pos1 - true_pos2, 10) + "\t" + strconv.FormatInt(true_pos1, 10) + "\t" + strconv.FormatInt(true_pos2, 10) + "\t")
		} else {
				file.WriteString("None\tNone\tNone\t")
		}
		file.WriteString(string(tokens[10]) + "\n")
		} else {
		file.WriteString("None\tNone\tNone\tNone\n")
	}
		
	file.Sync()
}


//Process FN SNPs and realted info
func ProcessSNPFNInfo() {
	if PRINT_FN {
		file, _ := os.Create(INPUT_INFO.SNP_call_file + ".fn")
		defer file.Close()

		var snp_pos int
		var snp []byte
		SNP_FN_map := make(map[int][]byte)
		for snp_pos, snp = range TRUE_VAR_COMP {
			SNP_FN_map[snp_pos] = snp
		}
		for snp_pos, snp = range TRUE_VAR_PART {
			SNP_FN_map[snp_pos] = snp
		}
		for snp_pos, snp = range TRUE_VAR_NONE {
			SNP_FN_map[snp_pos] = snp
		}
		SNP_Pos := make([]int, 0)
		var ok bool
		for snp_pos, snp = range SNP_FN_map {
			if _, ok = SNP_TRACE_INFO_MAP[uint32(snp_pos)]; !ok {
				SNP_Pos = append(SNP_Pos, snp_pos)
			}
		}
		sort.Ints(SNP_Pos)

		var align_pos int
		var at_arr []Align_trace_info
		var at Align_trace_info
		for pos := range SNP_Pos {
			for align_pos, at_arr = range ALIGN_TRACE_INFO_MAP {
				if align_pos <= pos && align_pos + 200 > pos {
					for _, at = range at_arr {
						WriteSNPFNInfo(file, align_pos, at, pos, SNP_FN_map[pos])
					}
				}
			}
		}
	}
}

//Write to file FN SNPs with all reads aligned to those locations
func WriteSNPFNInfo(file *os.File, align_pos int, at Align_trace_info, pos int, true_snp []byte) {

	file.WriteString(strconv.Itoa(pos) + "\t" + string(true_snp) + "\t")
	if len(at.snp_pos1) > 0 {
		for i, snp_pos := range at.snp_pos1 {
			if len(at.snp_base1[i]) > 0 {
				file.WriteString(strconv.Itoa(int(snp_pos)) + "\t" + string(at.snp_base1[i]) + "\t")
				file.WriteString(strconv.FormatFloat(QualtoProb(at.snp_baseq1[i][0]), 'f', 5, 32) + "\t")
			} else {
				file.WriteString(strconv.Itoa(int(snp_pos)) + "\t" + "." + "\t")
				file.WriteString(strconv.FormatFloat(QualtoProb('I'), 'f', 5, 32) + "\t")
			}
		}
	}
	if len(at.snp_pos2) > 0 {
		for i, snp_pos := range at.snp_pos2 {
			if len(at.snp_base2[i]) > 0 {
				file.WriteString(strconv.Itoa(int(snp_pos)) + "\t" + string(at.snp_base2[i]) + "\t")
				file.WriteString(strconv.FormatFloat(QualtoProb(at.snp_baseq2[i][0]), 'f', 5, 32) + "\t")
			} else {
				file.WriteString(strconv.Itoa(int(snp_pos)) + "\t" + "." + "\t")
				file.WriteString(strconv.FormatFloat(QualtoProb('I'), 'f', 5, 32) + "\t")
			}
		}
	}
	if at.align_pos1 != 0 && at.align_pos2 != 0 {
		file.WriteString(strconv.Itoa(at.align_pos1 - at.align_pos2) + "\t")
	} else {
		file.WriteString("None\t")
	}
	file.WriteString(strconv.Itoa(at.align_pos1) + "\t" + strconv.Itoa(at.align_pos2) + "\t")
	file.WriteString(strconv.Itoa(at.align_dis1) + "\t" + strconv.Itoa(at.align_dis2) + "\t")

	tokens := bytes.Split(at.read_info1, []byte{'_'})
	if len(tokens) >= 11 {
		true_pos1, err1 := strconv.ParseInt(string(tokens[2]), 10, 64)
		true_pos2, err2 := strconv.ParseInt(string(tokens[3]), 10, 64)
		if err1 == nil && err2 == nil {
			file.WriteString(strconv.FormatInt(true_pos1 - true_pos2, 10) + "\t" + strconv.FormatInt(true_pos1, 10) + "\t" + strconv.FormatInt(true_pos2, 10) + "\t")
		} else {
			file.WriteString("None\tNone\tNone\t")
		}
		file.WriteString(string(tokens[10]))
	} else {
		file.WriteString("None\tNone\tNone\tNone")
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
