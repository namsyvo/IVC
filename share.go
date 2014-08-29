//---------------------------------------------------------------------------------------------------
// Sharing variables and functions for modules in ISC package.
// Copyright 2014 Nam Sy Vo
//---------------------------------------------------------------------------------------------------

package isc

import (
	"github.com/vtphan/fmi"
	"math"
)

//Global constants and variables
var (
	INF             int   = math.MaxInt16 // Value for Infinity
	EMPTY_INT_SLICE []int = make([]int, 0)
)

//Index for SNP caller
type Index struct {
	SEQ            []byte            //store reference multigenomes
	SNP_PROF       map[int][][]byte  //hash table of SNP Profile (position, snps)
	SNP_AF         map[int][]float32 //allele frequency of SNP Profile (position, af of snps)
	SAME_LEN_SNP   map[int]int       //hash table to indicate if SNPs has same length
	SORTED_SNP_POS []int             //sorted array of SNP positions
	REV_FMI        fmi.Index         //FM-index of reverse multigenomes
}

//Input information
type InputInfo struct {
	//File names for:
	Genome_file    string //reference multigenome
	SNP_file       string //SNP profile
	Index_file     string //Index of original reference genomes
	Rev_index_file string //Index of reverse reference genomes
	Read_file_1    string //first end read
	Read_file_2    string //second end read
	SNP_call_file  string //store SNP call
	Search_mode    int    //searching mode for finding seeds
	Start_pos      int    //starting postion on reads for finding seeds
	Search_step    int    //step for searching in deterministic mode
	Proc_num       int    //maximum number of CPUs using by Go
	Routine_num    int    //number of goroutines
}

//Parameter used in alignment algorithm
type ParaInfo struct {
	Dist_thres      int     //threshold for distances between reads and multigenomes
	Iter_num        int     //number of random iterations to find proper seeds
	Max_match       int     //maximum number of matches
	Seq_err         float32 //average sequencing error, estmated from reads with real reads
	Err_var_factor  int     //factor for standard variation of sequencing error
	Iter_num_factor int     //factor for number of iterations 
	Read_len        int     //read length, calculated from read files
}

//Storing read and information
type ReadInfo struct {
	Read1		[]byte //first end read
	Read2		[]byte //second end read
	Qual1		[]byte //quality info of 1st read
	Qual2		[]byte //quality info of 2nd read
}

//"Global" variables used in alignment process (computing distance, snp call)
type AlignInfo struct {
	Bw_Dis   [][]int    // Distance matrix for backward alignment
	Fw_Dis   [][]int    // Distance matrix for forward alignment
	Bw_Trace [][][]byte // SNP trace matrix for backward alignment
	Fw_Trace [][][]byte // SNP trace matrix for forward alignment
}

//Assigning reads to ReadInfo.
func (read_info *ReadInfo) AssignReads(read1, read2, qual1, qual2 []byte) {
	read_info.Read1 = read1
	read_info.Read2 = read2
	read_info.Qual1 = qual1
	read_info.Qual2 = qual2
}

//Computing reverse, reverse complement, and complement of a read.
func RevComp(read []byte) ([]byte, []byte, []byte) {
	read_len := len(read)
	rev_read, rev_comp_read, comp_read :=
		make([]byte, read_len), make([]byte, read_len), make([]byte, read_len)
	for i, elem := range read {
		rev_read[i] = read[read_len-1-i]
		if elem == 'A' {
			rev_comp_read[read_len-1-i] = 'T'
			comp_read[i] = 'T'
		} else if elem == 'T' {
			rev_comp_read[read_len-1-i] = 'A'
			comp_read[i] = 'A'
		} else if elem == 'C' {
			rev_comp_read[read_len-1-i] = 'G'
			comp_read[i] = 'G'
		} else if elem == 'G' {
			rev_comp_read[read_len-1-i] = 'C'
			comp_read[i] = 'C'
		} else {
			rev_comp_read[read_len-1-i] = elem
			comp_read[i] = elem
		}
	}
	return rev_read, rev_comp_read, comp_read
}

//Allocating memory for share variables for alignment process
func (align_info *AlignInfo) InitAlignInfo(arr_len int) {
	InitMatrix(arr_len, &align_info.Bw_Dis, &align_info.Bw_Trace)
	InitMatrix(arr_len, &align_info.Fw_Dis, &align_info.Fw_Trace)
}

//Initializing variables for computing distance and alignment between reads and multi-genomes.
func InitMatrix(arr_len int, dis_mtr *[][]int, trace_mtr *[][][]byte) {
	*dis_mtr = make([][]int, arr_len+1)
	for i := 0; i <= arr_len; i++ {
		(*dis_mtr)[i] = make([]int, arr_len+1)
	}
	*trace_mtr = make([][][]byte, arr_len)
	for i := 0; i < arr_len; i++ {
		(*trace_mtr)[i] = make([][]byte, arr_len)
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
		m = L + (R-L)*((i-A[L])/(A[R]-A[L])) //out of range is possible here
		if A[m] < i {
			L = m + 1
		} else if A[m] > i {
			R = m - 1
		} else {
			return i <= j
		}
	}
	return i <= j && L < len(A) && i <= A[L] && j >= A[L]
}
