//---------------------------------------------------------------------------------------------------
// Sharing variables and functions for modules in ISC package.
// Copyright 2014 Nam Sy Vo
//---------------------------------------------------------------------------------------------------

package isc

import (
	"bufio"
	"fmt"
	"os"
	"strings"
	"strconv"
	"bytes"
	"sort"
    "github.com/vtphan/fmi"
	"math"
)

//Global constants and variables
var (
	INF int = math.MaxInt16 // Value for Infinity
	EMPTY_INT_SLICE []int = make([]int, 0)
)

//Index for SNP caller
type Index struct {
    SEQ []byte //store reference multigenomes
    SNP_PROF map[int][][]byte //hash table of SNP Profile (position, snps)
    SNP_AF map[int][]float32 //allele frequency of SNP Profile (position, af of snps)
    SAME_LEN_SNP map[int]int //hash table to indicate if SNPs has same length
    SORTED_SNP_POS []int //sorted array of SNP positions
    REV_FMI fmi.Index //FM-index of reverse multigenomes
}

//Input information
type InputInfo struct {
	//File names for:
    Genome_file string //reference multigenome
    SNP_file string //SNP profile
    Index_file string //Index of original reference genomes
    Rev_index_file string //Index of reverse reference genomes
    Read_file_1 string //first end read
    Read_file_2 string //second end read
    SNP_call_file string //store SNP call
	Search_mode int //searching mode for finding seeds
	Start_pos int //starting postion on reads for finding seeds
	Search_step int //step for searching in deterministic mode
	Proc_num int //maximum number of CPUs using by Go
	Routine_num int //number of goroutines
}

//Parameter used in alignment algorithm
type ParaInfo struct {
    Dist_thres int //threshold for distances between reads and multigenomes
    Iter_num int //number of random iterations to find proper seeds
    Max_match int //maximum number of matches
	Err_var_factor int //factor for standard variation for sequencing error
	Iter_num_factor int //factor for 
}

//"Global" variables used in alignment process (computing distance, snp call)
type AlignMem struct {
	Bw_D [][]int // Distance matrix for backward alignment
	Fw_D [][]int // Distance matrix for forward alignment
	Bw_T [][][]byte // SNP trace matrix for backward alignment
	Fw_T [][][]byte // SNP trace matrix for forward alignment
}

//Store read and information
type ReadInfo struct {
    Read1 []byte //first end read
	Read2 []byte //second end read
	Rev_read1 []byte //reverse of 1st end read
	Rev_read2 []byte //reverse of 2nd end read
	Rev_comp_read1 []byte //reverse complement of 1st end read
	Rev_comp_read2 []byte //reverse complement of 2nd end read
	Comp_read1 []byte //complement of 1st end read, ~ reverse of reverse complement of the read
	Comp_read2 []byte //complement of 2nd end read, ~ reverse of reverse complement of the read
	Qual_info_1 []byte //quality info of 1st read
	Qual_info_2 []byte //quality info of 2nd read
    Read_len int //read length
    Seq_err float32 //average sequencing error, estimated from read if real reads
}

func (read_info *ReadInfo) AllocMem() bool {
	if read_info.Read_len > 0 {
		read_info.Rev_read1 = make([]byte, read_info.Read_len)
		read_info.Rev_read2 = make([]byte, read_info.Read_len)
		read_info.Rev_comp_read1 = make([]byte, read_info.Read_len)
		read_info.Rev_comp_read2 = make([]byte, read_info.Read_len)
		read_info.Comp_read1 = make([]byte, read_info.Read_len)
		read_info.Comp_read2 = make([]byte, read_info.Read_len)
		return true
	} else {
		return false
	}
}

func (align_mem *AlignMem) AllocMem(arr_len int) {
	InitMatrix(arr_len, &align_mem.Bw_D, &align_mem.Bw_T)
	InitMatrix(arr_len, &align_mem.Fw_D, &align_mem.Fw_T)
}

//-------------------------------------------------------------------------------------------------
// Initializing variables for computing distance and alignment between reads and multi-genomes.
//-------------------------------------------------------------------------------------------------
func InitMatrix(arr_len int, D *[][]int, T *[][][]byte) {
	*D = make([][]int, arr_len + 1)
	for i:= 0; i <= arr_len; i++ {
		(*D)[i] = make([]int, arr_len + 1)
    }
	*T = make([][][]byte, arr_len)
    for i:= 0; i < arr_len; i++ {
		(*T)[i] = make([][]byte, arr_len)
    }
}


//--------------------------------------------------------------------------------------------------
// IntervalHasSNP determines whether [i, j] contains SNP positions which are stores in array A.
// This function impelements interpolation search. The array A must be sorted in increasing order.
//--------------------------------------------------------------------------------------------------
func IntervalHasSNP(A []int, i, j int) bool {
    L := 0
    R := len(A) - 1
    var m int
	for A[L] <= i && i <= A[R] && A[L] != A[R] {
		m = L + (R - L) * ((i - A[L]) / (A[R] - A[L]))  //out of range is possible here		
		if (A[m] < i) {
			L = m + 1;
		} else if A[m] > i {
			R = m - 1
		} else {
			return i <= j
		}
	}
	return i <= j && L < len(A) && i <= A[L] && j >= A[L]
}

//--------------------------------------------------------------------------------------------------
// Computing reverse, reverse complement, and complement of a read.
//--------------------------------------------------------------------------------------------------
func RevComp(read []byte, rev_read, rev_comp_read, comp_read []byte) {
	read_len := len(read)
    for i, elem := range read {
        rev_read[i] = read[read_len - 1 - i]
        if elem == 'A' {
            rev_comp_read[read_len - 1 - i] = 'T'
            comp_read[i] = 'T'
        } else if (elem == 'T') {
            rev_comp_read[read_len - 1 - i] = 'A'
            comp_read[i] = 'A'
        } else if (elem == 'C') {
            rev_comp_read[read_len - 1 - i] = 'G'
            comp_read[i] = 'G'
        } else if (elem == 'G') {
            rev_comp_read[read_len - 1 - i] = 'C'
            comp_read[i] = 'C'
        } else {
            rev_comp_read[read_len - 1 - i] = elem
            comp_read[i] = elem
        }
    }
}

//--------------------------------------------------------------------------------------------------
// ReadVCF
//--------------------------------------------------------------------------------------------------
func ReadVCF(sequence_file string) map[int]SNPProfile {
	array := make(map[int]SNPProfile)
	f,err := os.Open(sequence_file)
    if err != nil{
        fmt.Printf("%v\n", err)
        os.Exit(1)
    }

    defer f.Close()
    br := bufio.NewReader(f)
    //byte_array := bytes.Buffer{}

	for{
		line , err := br.ReadString('\n')
		if err != nil {
			//fmt.Printf("%v\n",err)
			break
		}		
		if line[0] == byte('#') {
			//fmt.Printf("%s \n",line)
		} else {
			sline := string(line)
			split := strings.Split(sline, "\t");
			//fmt.Printf("%s %s %s\n", split[1], split[3], split[4])
			pos, _ := strconv.ParseInt(split[1], 10, 64)
			pos = pos - 1
			if len(split[4]) > 1 {
				alt := strings.Split(split[4], ",")
				t := make([]string, len(alt)+1)
				t[0] = split[3]				
				for i := 0 ; i < len(alt) ; i++ {
					if alt[i] == "<DEL>" {
						t[i + 1] = "."
					} else {
						t[i + 1] = alt[i]
					}					
				}	
				//sort.Strings(t)
				//array[int(pos)] = SNP{t} // asign SNP at pos
				tmp, ok := array[int(pos)]
				if ok {
					t = append(t[ : 0], t[1 : ]...)
					tmp.Profile = append(tmp.Profile, t...)
				} else {
					tmp.Profile = append(tmp.Profile, t...)
				}
				sort.Strings(tmp.Profile)
				array[int(pos)] = tmp // append SNP at pos
				//fmt.Printf("pos=%d %q \n", pos, alt)
			} else {				
				//array[int(pos)] = SNP{[]string{split[3], split[4]}} // asign SNP at pos
				tmp, ok := array[int(pos)]
				if ok {
					if split[4] == "<DEL>" {
						tmp.Profile = append(tmp.Profile, ".")
					} else {
						tmp.Profile = append(tmp.Profile, split[4])
					}					
				} else {
					if split[4] == "<DEL>" {
						tmp.Profile = append(tmp.Profile, []string{split[3], "."}...)
					} else {
						tmp.Profile = append(tmp.Profile, []string{split[3], split[4]}...)
					}
				}
				sort.Strings(tmp.Profile)

				array[int(pos)]= tmp // append SNP at pos
				//fmt.Println(pos)
			}
		}
	}
    return array
}

//--------------------------------------------------------------------------------------------------
// ReadFASTA
//--------------------------------------------------------------------------------------------------
func ReadFASTA(sequence_file string) []byte {
    f,err := os.Open(sequence_file)
    if err != nil{
        fmt.Printf("%v\n",err)
        os.Exit(1)
    }

    defer f.Close()
    br := bufio.NewReader(f)
    byte_array := bytes.Buffer{}

    //line , err := br.ReadString('\n')
	_ , isPrefix, err := br.ReadLine()
	if err != nil || isPrefix{
		fmt.Printf("%v\n",err)
		os.Exit(1)
	}
    //fmt.Printf("%s",line)

    for {
        line, isPrefix, err := br.ReadLine()
        if err != nil || isPrefix {
            break
        } else {
            byte_array.Write([]byte(line))
        }
    }
    //byte_array.Write([]byte("$"))
    input := []byte(byte_array.String())
    return input
}
