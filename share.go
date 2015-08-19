//---------------------------------------------------------------------------------------------------
// IVC: share.go - Shared variables and functions.
// Copyright 2015 Nam Sy Vo.
//---------------------------------------------------------------------------------------------------

package ivc

import (
	"bytes"
	"fmt"
	"log"
	"math"
	//"math/rand"
)

//--------------------------------------------------------------------------------------------------
// Global constants and variables
//--------------------------------------------------------------------------------------------------
var (
	STD_BASES          = []byte{'A', 'C', 'G', 'T'} //Standard bases of DNA sequences
	NEW_SNP_RATE       = 0.00001                    //Value for prior probability of new alleles
	NEW_INDEL_RATE     = 0.000001                   //Value for prior probability of new indels
	NEW_SNP_RATE_LOG   = -math.Log10(NEW_SNP_RATE)
	NEW_INDEL_RATE_LOG = -math.Log10(NEW_INDEL_RATE)
)

//--------------------------------------------------------------------------------------------------
// Global variables for read alignment and Var calling processes.
//--------------------------------------------------------------------------------------------------
var (
	INPUT_INFO *InputInfo //Input information
	PARA_INFO  *ParaInfo  //Parameters information
	INDEX      *Index     //Index for alignment
)

//--------------------------------------------------------------------------------------------------
// Input information
//--------------------------------------------------------------------------------------------------
type InputInfo struct {
	//File names for:
	Ref_file       string  //reference multigenome
	Var_prof_file  string  //variant profile
	Index_file     string  //index of original reference genomes
	Rev_index_file string  //index of reverse reference genomes
	Read_file_1    string  //first end of read
	Read_file_2    string  //second end of read
	Var_call_file  string  //store Var call
	Search_mode    int     //searching mode for finding seeds
	Start_pos      int     //starting postion on reads for finding seeds
	Search_step    int     //step for searching in deterministic mode
	Proc_num       int     //maximum number of CPUs using by Go
	Routine_num    int     //number of goroutines
	Max_snum       int     //maximum number of seeds
	Max_psnum      int     //maximum number of paired-seeds
	Min_slen       int     //minimum length of seeds
	Max_slen       int     //maximum length of seeds
	Dist_thres     int     //threshold for distances between reads and multigenomes
	Prob_thres     float64 //threshold for alignment prob between reads and multigenomes
	Iter_num       int     //number of random iterations to find proper alignments
	Cpu_prof_file  string  //file to write cpu profile
	Mem_prof_file  string  //file to write mem profile
}

//--------------------------------------------------------------------------------------------------
// Parameter used in alignment algorithm
//--------------------------------------------------------------------------------------------------
type ParaInfo struct {
	Dist_thres      int     //threshold for distances between reads and multigenomes
	Prob_thres      float64 //threshold for probabilities of correct alignment
	Iter_num        int     //number of random iterations to find proper alignments
	Max_ins         int     //maximum insert size of two aligned ends
	Err_rate        float32 //average sequencing error rate, estmated from reads with real reads
	Err_var_factor  int     //factor for standard variation of sequencing error rate
	Mut_rate        float32 //average mutation rate, estmated from reference genome
	Mut_var_factor  int     //factor for standard variation of mutation rate
	Iter_num_factor int     //factor for number of iterations
	Read_len        int     //read length, calculated from read files
	Info_len        int     //maximum size of array to store read headers
	Sub_cost        float64 //cost of substitution for Hamming and Edit distance
	Gap_open_cost   float64 //cost of gap open for Edit distance
	Gap_ext_cost    float64 //cost of gap extension for Edit distance
	Seed_backup     int     //number of backup bases from seeds
	Indel_backup    int     //number of backup bases from known indels
	Ham_backup      int     //number of backup bases from Hamming alignment
}

//--------------------------------------------------------------------------------------------------
// SetPara sets values of parameters for alignment process
//--------------------------------------------------------------------------------------------------
func SetPara(read_len, info_len int, max_ins int, err_rate, mut_rate float32, dist_thres int, prob_thres float64, iter_num int) *ParaInfo {
	para_info := new(ParaInfo)
	para_info.Seed_backup = 6
	para_info.Indel_backup = 30
	para_info.Ham_backup = 15
	para_info.Err_var_factor = 4
	para_info.Mut_var_factor = 2
	para_info.Iter_num_factor = 2
	para_info.Max_ins = max_ins   //based on simulated data, will be estimated from reads with real data
	para_info.Err_rate = err_rate //will be replaced by error rate estimated from input reads
	para_info.Mut_rate = mut_rate //will be replaced by error rate estimated from input reads
	para_info.Read_len = read_len //will be replaced by read length taken from input reads
	para_info.Info_len = info_len //will be replaced by maximum info length estimated from input reads

	err := float64(para_info.Err_rate)
	if dist_thres != 0 {
		para_info.Dist_thres = dist_thres
	} else {
		rlen := float64(para_info.Read_len)
		mut := float64(para_info.Mut_rate)
		k1 := float64(para_info.Err_var_factor)
		k2 := float64(para_info.Mut_var_factor)
		para_info.Dist_thres = int(math.Ceil(err*rlen+k1*math.Sqrt(rlen*err*(1-err)))) +
			int(math.Ceil(mut*rlen+k2*math.Sqrt(rlen*mut*(1-mut))))
	}
	if prob_thres != 0 {
		para_info.Prob_thres = prob_thres
	} else {
		para_info.Prob_thres = -float64(para_info.Dist_thres)*math.Log10(1-err) - float64(para_info.Dist_thres)*math.Log10(NEW_INDEL_RATE)
	}
	if iter_num != 0 {
		para_info.Iter_num = iter_num
	} else {
		para_info.Iter_num = para_info.Iter_num_factor * (para_info.Dist_thres + 1)
	}

	para_info.Sub_cost = 4.0      //NEW_SNP_RATE_LOG
	para_info.Gap_open_cost = 6.0 //NEW_INDEL_RATE_LOG
	para_info.Gap_ext_cost = 1.0  //NEW_SNP_RATE_LOG

	log.Printf("Prog paras:\tDist_thres=%d, Prob_thres=%.5f, Iter_num=%d, Max_ins=%d, Max_err=%.5f, Err_var_factor=%d,"+
		" Mut_rate=%.5f, Mut_var_factor=%d, Iter_num_factor=%d, Read_len=%d, Info_len=%d,"+
		" Sub_cost=%.5f, Gap_open_cost=%.5f, Gap_ext_cost=%.5f, Seed_backup=%d, Indel_backup=%d, Ham_backup=%d",
		para_info.Dist_thres, para_info.Prob_thres, para_info.Iter_num, para_info.Max_ins, para_info.Err_rate, para_info.Err_var_factor,
		para_info.Mut_rate, para_info.Mut_var_factor, para_info.Iter_num_factor, para_info.Read_len, para_info.Info_len,
		para_info.Sub_cost, para_info.Gap_open_cost, para_info.Gap_ext_cost, para_info.Seed_backup, para_info.Indel_backup, para_info.Ham_backup)

	return para_info
}

//--------------------------------------------------------------------------------------------------
// Information of input reads
//--------------------------------------------------------------------------------------------------
type ReadInfo struct {
	Read1, Read2                   []byte //first and second ends
	Qual1, Qual2                   []byte //quality info of the first read and second ends
	Rev_read1, Rev_read2           []byte //reverse of the first and second ends
	Rev_comp_read1, Rev_comp_read2 []byte //reverse complement of the first and second ends
	Comp_read1, Comp_read2         []byte //complement of the first and second ends
	Rev_qual1, Rev_qual2           []byte //quality of reverse of the first and second ends
	Info1, Info2                   []byte //info of the first and second ends
}

//--------------------------------------------------------------------------------------------------
// InitReadInfo creates a ReadInfo object and initializes its content
//--------------------------------------------------------------------------------------------------
func InitReadInfo(read_len, info_len int) *ReadInfo {
	read_info := new(ReadInfo)
	read_info.Read1, read_info.Read2 = make([]byte, read_len), make([]byte, read_len)
	read_info.Qual1, read_info.Qual2 = make([]byte, read_len), make([]byte, read_len)
	read_info.Rev_read1, read_info.Rev_read2 = make([]byte, read_len), make([]byte, read_len)
	read_info.Rev_comp_read1, read_info.Rev_comp_read2 = make([]byte, read_len), make([]byte, read_len)
	read_info.Comp_read1, read_info.Comp_read2 = make([]byte, read_len), make([]byte, read_len)
	read_info.Rev_qual1, read_info.Rev_qual2 = make([]byte, read_len), make([]byte, read_len)
	read_info.Info1, read_info.Info2 = make([]byte, info_len), make([]byte, info_len)
	return read_info
}

//--------------------------------------------------------------------------------------------------
// PrintReads prints read information
//--------------------------------------------------------------------------------------------------
func (read_info *ReadInfo) PrintReads() {
	fmt.Println("read1: ", string(read_info.Read1))
	fmt.Println("read2: ", string(read_info.Read2))
	fmt.Println("qual1: ", string(read_info.Qual1))
	fmt.Println("qual1: ", string(read_info.Qual2))
	fmt.Println("info1: ", string(read_info.Info1))
	fmt.Println("info2: ", string(read_info.Info2))
}

//--------------------------------------------------------------------------------------------------
// RevComp computes reverse, reverse complement, and complement of a read.
//--------------------------------------------------------------------------------------------------
func RevComp(read, qual []byte, rev_read, rev_comp_read, comp_read, rev_qual []byte) {
	read_len := len(read)
	for i, elem := range read {
		rev_qual[i] = qual[read_len-1-i]
		if elem == 'A' {
			rev_read[read_len-1-i] = 'A'
			rev_comp_read[read_len-1-i] = 'T'
			comp_read[i] = 'T'
		} else if elem == 'T' {
			rev_read[read_len-1-i] = 'T'
			rev_comp_read[read_len-1-i] = 'A'
			comp_read[i] = 'A'
		} else if elem == 'C' {
			rev_read[read_len-1-i] = 'C'
			rev_comp_read[read_len-1-i] = 'G'
			comp_read[i] = 'G'
		} else if elem == 'G' {
			rev_read[read_len-1-i] = 'G'
			rev_comp_read[read_len-1-i] = 'C'
			comp_read[i] = 'C'
		} else {
			rev_read[read_len-1-i] = elem
			rev_comp_read[read_len-1-i] = elem
			comp_read[i] = elem
		}
	}
}

//---------------------------------------------------------------------------------------------------
// SeedInfo represents info of seeds between reads and the multigenome.
//---------------------------------------------------------------------------------------------------
type SeedInfo struct {
	s_pos  []int  //staring position of seeds on reads.
	e_pos  []int  //ending position of seeds on reads.
	m_pos  []int  //(left-most) matching position of seeds on the reference multigenome.
	strand []bool //strand (forward or reverse) of matches on the reference multigenome.
}

//--------------------------------------------------------------------------------------------------
// Alignment information, served as shared variables between functions for alignment process
//--------------------------------------------------------------------------------------------------
type EditAlnInfo struct {
	Bw_Dist_D, Bw_Dist_IS, Bw_Dist_IT    [][]float64 // Distance matrix for backward alignment
	Bw_Trace_D, Bw_Trace_IS, Bw_Trace_IT [][][]int   // Backtrace matrix for backward alignment
	Fw_Dist_D, Fw_Dist_IS, Fw_Dist_IT    [][]float64 // Distance matrix for forward alignment
	Fw_Trace_D, Fw_Trace_IS, Fw_Trace_IT [][][]int   // Backtrace matrix for forward alignment
}

//--------------------------------------------------------------------------------------------------
// InitEditAlnInfo allocates memory for share variables for alignment process
//--------------------------------------------------------------------------------------------------
func InitEditAlnInfo(arr_len int) *EditAlnInfo {
	aln_info := new(EditAlnInfo)
	aln_info.Bw_Dist_D, aln_info.Bw_Trace_D = InitEditAlnMat(arr_len)
	aln_info.Bw_Dist_IS, aln_info.Bw_Trace_IS = InitEditAlnMat(arr_len)
	aln_info.Bw_Dist_IT, aln_info.Bw_Trace_IT = InitEditAlnMat(arr_len)
	aln_info.Fw_Dist_D, aln_info.Fw_Trace_D = InitEditAlnMat(arr_len)
	aln_info.Fw_Dist_IS, aln_info.Fw_Trace_IS = InitEditAlnMat(arr_len)
	aln_info.Fw_Dist_IT, aln_info.Fw_Trace_IT = InitEditAlnMat(arr_len)
	return aln_info
}

//--------------------------------------------------------------------------------------------------
// InitEditAlnMat initializes variables for computing distance and alignment between reads and multi-genomes.
//--------------------------------------------------------------------------------------------------
func InitEditAlnMat(arr_len int) ([][]float64, [][][]int) {
	dis_mat := make([][]float64, arr_len+1)
	for i := 0; i <= arr_len; i++ {
		dis_mat[i] = make([]float64, arr_len+1)
	}
	trace_mat := make([][][]int, arr_len+1)
	for i := 0; i <= arr_len; i++ {
		trace_mat[i] = make([][]int, arr_len+1)
		for j := 0; j <= arr_len; j++ {
			trace_mat[i][j] = make([]int, 3)
		}
	}
	return dis_mat, trace_mat
}

//--------------------------------------------------------------------------------------------------
// Utility functions
//--------------------------------------------------------------------------------------------------

//--------------------------------------------------------------------------------------------------
// IntervalHasVariants determines whether [i, j] contains variant positions which are stores in array A.
// This function implements interpolation search. The array A must be sorted in increasing order.
//--------------------------------------------------------------------------------------------------
func IntervalHasVariants(A []int, i, j int) bool {
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

//--------------------------------------------------------------------------------------------------
// SplitN splits a slice of bytes using an memory-efficient method.
//--------------------------------------------------------------------------------------------------
func SplitN(s, sep []byte, n int) ([][]byte, int) {
	first_idx, sep_idx := 0, 0
	sep_num := 0
	t := make([][]byte, 0)
	for first_idx < len(s) {
		sep_idx = bytes.Index(s[first_idx:], sep)
		if sep_idx != -1 {
			sep_num++
			tmp := make([]byte, first_idx+sep_idx-first_idx)
			copy(tmp, s[first_idx:first_idx+sep_idx])
			t = append(t, tmp)
			if sep_num == n {
				return t, sep_num
			}
			first_idx = first_idx + sep_idx + 1
		} else {
			return t, sep_num
		}
	}
	return t, sep_num
}

//--------------------------------------------------------------------------------------------------
// IndexN returns index of a pattern in a slice of bytes.
//--------------------------------------------------------------------------------------------------
func IndexN(s, sep []byte, n int) int {
	first_idx, sep_idx := 0, 0
	sep_num := 0
	for first_idx < len(s) {
		sep_idx = bytes.Index(s[first_idx:], sep)
		if sep_idx != -1 {
			sep_num++
			if sep_num == n {
				return first_idx + sep_idx
			}
			first_idx = first_idx + sep_idx + 1
		} else {
			return -1
		}
	}
	return -1
}
