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
var MEM_STATS = new(runtime.MemStats)

//Printing memory information
func PrintMemStats(mesg string) {
    runtime.ReadMemStats(MEM_STATS)
    log.Printf(mesg + "\t%d\t%d\t%d\t%d\t%d\t%.3f%.3f%.3f%.3f%.3f",
		MEM_STATS.Alloc, MEM_STATS.TotalAlloc, MEM_STATS.Sys, MEM_STATS.HeapAlloc, MEM_STATS.HeapSys,
		float64(MEM_STATS.Alloc)/(math.Pow(1024, 3)), float64(MEM_STATS.TotalAlloc)/(math.Pow(1024, 3)),
		float64(MEM_STATS.Sys)/(math.Pow(1024, 3)), float64(MEM_STATS.HeapAlloc)/(math.Pow(1024, 3)),
		float64(MEM_STATS.HeapSys)/(math.Pow(1024, 3)))
}

//Global variable for debugging
var (
	DEBUG_INFO chan Debug_info
	TRUE_VAR_COMP, TRUE_VAR_PART, TRUE_VAR_NONE map[int][]byte
)

type Debug_info struct {
	read_info1, read_info2 []byte
	align_pos1, align_pos2 int
	snp_pos1, snp_pos2 []uint32
	snp_base1, snp_base2 [][]byte
	snp_baseq1, snp_baseq2 [][]byte
}

/*
 Log file format:
  align_pos_diff  align_pos1  align_pos2  snp_pos  snp_base  snp_base_qual  snp_qual_prob  true_pos_diff  true_pos1  true_pos2  read_id
  ...
 where:
  one line is corresponding to one variant call
  value is "None" if not exist

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
  other: info other locations
*/
//--------------------------------------------------------------------------------------------------

//Read debug info from channel and process them
func ProcessDebugInfo() {
	var i int
	var snp_pos uint32

	files := make([]*os.File, 14)
	file_names := []string{"tp_snp_comp", "fp_snp_comp", "tp_indel_comp", "fp_indel_comp", "tp_snp_part", "fp_snp_part", "tp_indel_part", "fp_indel_part", "tp_snp_none", "fp_snp_none", "tp_indel_none", "fp_indel_none", "fp_snp_other", "fp_indel_other"}
	for i, file_name := range file_names {
		files[i], _ = os.Create(INPUT_INFO.SNP_call_file + "." + file_name)
		defer files[i].Close()
	}
	for d := range DEBUG_INFO {
		if len(d.snp_pos1) > 0 {
			for i, snp_pos = range d.snp_pos1 {
				if len(d.snp_base1[i]) == 0 {
					OutputDebugInfo(files, d, snp_pos, []byte{'.'}, []byte{'I'})
				} else {
					OutputDebugInfo(files, d, snp_pos, d.snp_base1[i], d.snp_baseq1[i])
				}
			}
		}
		if len(d.snp_pos2) > 0 {
			for i, snp_pos = range d.snp_pos2 {
				if len(d.snp_base2[i]) == 0 {
					OutputDebugInfo(files, d, snp_pos, []byte{'.'}, []byte{'I'})
				} else {
					OutputDebugInfo(files, d, snp_pos, d.snp_base2[i], d.snp_baseq2[i])
				}
			}
		}
	}
}

//Output debug info to proper files
func OutputDebugInfo(files []*os.File, d Debug_info, snp_pos uint32, snp_base []byte, snp_baseq []byte) {
	var ok bool
	var val []byte
	if val, ok = TRUE_VAR_COMP[int(snp_pos)]; ok {
		if len(snp_base) == 1 {
			if snp_base[0] == val[0] {
				WriteDebugInfo(files[0], d, snp_pos, snp_base, snp_baseq)
			} else {
				WriteDebugInfo(files[1], d, snp_pos, snp_base, snp_baseq)
			}
		} else {
			if bytes.Equal(snp_base, val) {
				WriteDebugInfo(files[2], d, snp_pos, snp_base, snp_baseq)
			} else {
				WriteDebugInfo(files[3], d, snp_pos, snp_base, snp_baseq)
			}
		}
	} else if val, ok = TRUE_VAR_PART[int(snp_pos)]; ok {
		if len(snp_base) == 1 {
			if snp_base[0] == val[0] {
				WriteDebugInfo(files[4], d, snp_pos, snp_base, snp_baseq)
			} else {
				WriteDebugInfo(files[5], d, snp_pos, snp_base, snp_baseq)
			}
		} else {
			if bytes.Equal(snp_base, val) {
				WriteDebugInfo(files[6], d, snp_pos, snp_base, snp_baseq)
			} else {
				WriteDebugInfo(files[7], d, snp_pos, snp_base, snp_baseq)
			}
		}
	} else if val, ok = TRUE_VAR_NONE[int(snp_pos)]; ok {
		if len(snp_base) == 1 {
			if snp_base[0] == val[0] {
				WriteDebugInfo(files[8], d, snp_pos, snp_base, snp_baseq)
			} else {
				WriteDebugInfo(files[9], d, snp_pos, snp_base, snp_baseq)
			}
		} else {
			if bytes.Equal(snp_base, val) {
				WriteDebugInfo(files[10], d, snp_pos, snp_base, snp_baseq)
			} else {
				WriteDebugInfo(files[11], d, snp_pos, snp_base, snp_baseq)
			}
		}
	} else {
		if len(snp_base) == 1 {
			WriteDebugInfo(files[12], d, snp_pos, snp_base, snp_baseq)
		} else {
			WriteDebugInfo(files[13], d, snp_pos, snp_base, snp_baseq)
		}
	}
}

//Write debug info to files
func WriteDebugInfo(file *os.File, d Debug_info, snp_pos uint32, snp_base []byte, snp_baseq []byte) {

	if d.align_pos1 != 0 && d.align_pos2 != 0 {
		file.WriteString(strconv.Itoa(d.align_pos1 - d.align_pos2) + "\t")
	} else {
		file.WriteString("None\t")
	}
	file.WriteString(strconv.Itoa(d.align_pos1) + "\t" + strconv.Itoa(d.align_pos2) + "\t")
	file.WriteString(strconv.Itoa(int(snp_pos)) + "\t" + string(snp_base) + "\t")
	file.WriteString(strconv.FormatFloat(math.Pow(10, -(float64(snp_baseq[0]) - 33)/10.0), 'f', 5, 32) + "\t")
	//file.WriteString(string(d.read_info1[32:]) + "\n")
	
	tokens := bytes.Split(d.read_info1, []byte{'_'})
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

//Load true variants for debugging
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
