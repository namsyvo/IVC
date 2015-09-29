//---------------------------------------------------------------------------------------------------
// IVC: share.go - Shared variables and functions.
// Copyright 2015 Nam Sy Vo.
//---------------------------------------------------------------------------------------------------

package ivc

import (
	"bytes"
	"log"
	"os"
	"runtime"
	"runtime/pprof"
)

//--------------------------------------------------------------------------------------------------
// Global constants
//--------------------------------------------------------------------------------------------------
var (
	NEW_SNP_RATE   = 0.00001   //probability of new alleles
	NEW_INDEL_RATE = 0.000001  //probability of new indels
	INDEL_ERR_RATE = 0.0000001 //probability of indel error
)

//--------------------------------------------------------------------------------------------------
// Global variables for whole variant calling process.
//--------------------------------------------------------------------------------------------------
var (
	PARA_INFO    *ParaInfo    //program parameters
	MULTI_GENOME *MultiGenome //multi-genome instance
)

//--------------------------------------------------------------------------------------------------
// Global variables for calculating variant quality.
//--------------------------------------------------------------------------------------------------
var (
	Q2C map[byte]float64 //pre-calculated alignment cost  based on Phred-based quality.
	Q2E map[byte]float64 //pre-calculated error probability based on Phred-based quality.
	Q2P map[byte]float64 //pre-calculated not-error probability based on Phred-based quality.
	L2E []float64        //pre-calculated indel error corresponding to length of indels.
)

//--------------------------------------------------------------------------------------------------
// Parameter information
//--------------------------------------------------------------------------------------------------
type ParaInfo struct {
	//Input file names:
	Ref_file       string //reference multigenome
	Var_prof_file  string //variant profile
	Index_file     string //index of original reference genomes
	Rev_index_file string //index of reverse reference genomes
	Read_file_1    string //first end of read
	Read_file_2    string //second end of read
	Var_call_file  string //store Var call

	//Input paras:
	Search_mode int     //searching mode for finding seeds
	Start_pos   int     //starting postion on reads for finding seeds
	Search_step int     //step for searching in deterministic mode
	Max_snum    int     //maximum number of seeds
	Max_psnum   int     //maximum number of paired-seeds
	Min_slen    int     //minimum length of seeds
	Max_slen    int     //maximum length of seeds
	Dist_thres  float64 //threshold for distances between reads and multigenomes
	Iter_num    int     //number of random iterations to find proper alignments
	Sub_cost    float64 //cost of substitution for Hamming and Edit distance
	Gap_open    float64 //cost of gap open for Edit distance
	Gap_ext     float64 //cost of gap extension for Edit distance
	Proc_num    int     //maximum number of CPUs using by Go
	Debug_mode  bool    //debug mode for output

	//Estimated paras:
	Read_len        int     //read length, calculated from read files
	Info_len        int     //maximum size of array to store read headers
	Max_ins         int     //maximum insert size of two aligned ends
	Err_rate        float32 //average sequencing error rate, estmated from reads with real reads
	Err_var_factor  int     //factor for standard variation of sequencing error rate
	Mut_rate        float32 //average mutation rate, estmated from reference genome
	Mut_var_factor  int     //factor for standard variation of mutation rate
	Iter_num_factor int     //factor for number of iterations
	Seed_backup     int     //number of backup bases from seeds
	Ham_backup      int     //number of backup bases from Hamming alignment
	Indel_backup    int     //number of backup bases from known indels
}

//--------------------------------------------------------------------------------------------------
// Read input information and set up parameters
//--------------------------------------------------------------------------------------------------
func Setup(input_para_info *ParaInfo) {

	log.Printf("----------------------------------------------------------------------------------------")
	log.Printf("Checking input information and seting up parameters...")

	//Check input files
	if _, e := os.Stat(input_para_info.Read_file_1); e != nil {
		log.Printf("Error: Read_file_1 does not exists! (err: %s)", e)
		os.Exit(1)
	}
	if _, e := os.Stat(input_para_info.Read_file_2); e != nil {
		log.Printf("Error: Read_file_2 does not exists! (err: %s)", e)
		os.Exit(1)
	}

	MEM_STATS = new(runtime.MemStats)

	PARA_INFO = SetupPara(input_para_info)
	runtime.GOMAXPROCS(PARA_INFO.Proc_num)

	if input_para_info.Debug_mode {
		var err error
		CPU_FILE, err = os.Create(input_para_info.Var_call_file + ".cprof")
		if err != nil {
			log.Fatal(err)
		}
		pprof.StartCPUProfile(CPU_FILE)
		defer pprof.StopCPUProfile()

		MEM_FILE, err = os.Create(input_para_info.Var_call_file + ".mprof")
		if err != nil {
			log.Fatal(err)
		}
		defer MEM_FILE.Close()
		log.Printf("Debug mode:\tCpu_prof_file: %s, Mem_prof_file: %s", input_para_info.Var_call_file+".cprof", input_para_info.Var_call_file+".mprof")
	}

	log.Printf("Finish checking input information and seting up parameters.")
}

//--------------------------------------------------------------------------------------------------
// SetupPara setups values of parameters for alignment process
//--------------------------------------------------------------------------------------------------
func SetupPara(input_para_info *ParaInfo) *ParaInfo {

	para_info := input_para_info

	//Estimate parameters for the variant caller
	//100 is length of testing reads, 500 is maximum length of info line of reads,
	//will be set up based on input reads
	para_info.Read_len = 100
	para_info.Info_len = 500

	//700 is maximum insert size of paired-end testing reads
	//will be estimated based on input reads
	para_info.Max_ins = 700

	//0.0015 is maximum sequencing error rate of testing reads, 0.01 is mutation rate of testing data,
	//will be replaced by based on input reads
	para_info.Err_rate = 0.0015
	para_info.Mut_rate = 0.01

	para_info.Err_var_factor = 4
	para_info.Mut_var_factor = 2
	para_info.Iter_num_factor = 2

	para_info.Seed_backup = 6
	para_info.Ham_backup = 15
	para_info.Indel_backup = 30

	//Setup input parameters if not specified
	if input_para_info.Search_mode == 0 {
		para_info.Search_mode = 1
		log.Printf("No or invalid input for searching mode, use default strategy (randomizaion).")
	} else if input_para_info.Search_mode == 1 {
		if input_para_info.Start_pos == 0 {
			para_info.Start_pos = 512
			log.Printf("Deterministic search mode: no or invalid input for start postion on reads to find seeds, use default value (%d).", para_info.Start_pos)
		}
		if input_para_info.Search_step == 0 {
			para_info.Search_step = 512
			log.Printf("Deterministic search mode: no or invalid input for searching step, use default value (%d).", para_info.Search_step)
		}
	}
	if input_para_info.Max_snum == 0 {
		para_info.Max_snum = 512
		log.Printf("No or invalid input for maximum number of seeds, use default value (%d).", para_info.Max_snum)
	}
	if input_para_info.Max_psnum == 0 {
		para_info.Max_psnum = 128
		log.Printf("No or invalid input for maximum number of paired-seeds, use default value (%d).", para_info.Max_psnum)
	}
	if input_para_info.Min_slen == 0 {
		para_info.Min_slen = 15
		log.Printf("No or invalid input for minimum length of seeds, use default value (%d).", para_info.Min_slen)
	}
	if input_para_info.Max_slen == 0 {
		para_info.Max_slen = 25
		log.Printf("No or invalid input for maximum length of seeds, use default value (%d).", para_info.Max_slen)
	}
	if input_para_info.Sub_cost == 0 {
		para_info.Sub_cost = 4
		log.Printf("No or invalid input for substitution cost of alignment, use default value (%.1f).", para_info.Sub_cost)
	}
	if input_para_info.Gap_open == 0 {
		para_info.Gap_open = 4.1
		log.Printf("No or invalid input for gap open cost of alignment, use default value (%.1f).", para_info.Gap_open)
	}
	if input_para_info.Gap_ext == 0 {
		para_info.Gap_ext = 1
		log.Printf("No or invalid input for gap extension cost of alignment, use default value (%.1f).", para_info.Gap_ext)
	}

	if input_para_info.Dist_thres == 0 {
		/*
			err := float64(para_info.Err_rate)
			rlen := float64(para_info.Read_len)
			mut := float64(para_info.Mut_rate)
			k1 := float64(para_info.Err_var_factor)
			k2 := float64(para_info.Mut_var_factor)
			var_dist = int(math.Ceil(err*rlen+k1*math.Sqrt(rlen*err*(1-err)))) + int(math.Ceil(mut*rlen+k2*math.Sqrt(rlen*mut*(1-mut))))
		    para_info.Dist_thres = -float64(var_dist)*math.Log10(1-err) - float64(var_dist)*math.Log10(NEW_INDEL_RATE)
		*/
		para_info.Dist_thres = 36
		log.Printf("No or invalid input for threshold of alignment distance, calculate based on input data (%.1f).", para_info.Dist_thres)
	}
	if input_para_info.Iter_num == 0 {
		//para_info.Iter_num = para_info.Iter_num_factor * (para_info.Dist_thres + 1)
		para_info.Iter_num = 12
		log.Printf("No or invalid input for numbers of random iterations, calculate based on input data (%d).", para_info.Iter_num)
	}

	if input_para_info.Proc_num == 0 {
		para_info.Proc_num = runtime.NumCPU()
		log.Printf("No or invalid input for number of threads, use maximum number of CPUs of the current machine (%d).", para_info.Proc_num)
	}

	log.Printf("Input files:\tGenome_file: %s, Var_file: %s, Index_file=%s, Read_file_1=%s, Read_file_2=%s, Var_call_file=%s", 
		para_info.Ref_file, para_info.Var_prof_file, para_info.Rev_index_file, para_info.Read_file_1, para_info.Read_file_2, para_info.Var_call_file)

	log.Printf("Input paras:\tSearch_mode=%d, Start_pos=%d, Search_step=%d, Max_snum=%d, Max_psnum=%d, "+
		"Min_slen=%d, Max_slen=%d, Dist_thres=%.1f, Iter_num=%d, Sub_cost=%.1f, Gap_open=%.1f, Gap_ext=%.1f, Proc_num=%d, Debug_mode=%t",
		para_info.Search_mode, para_info.Start_pos, para_info.Search_step, para_info.Max_snum, para_info.Max_psnum, para_info.Min_slen, para_info.Max_slen,
		para_info.Dist_thres, para_info.Iter_num, para_info.Sub_cost, para_info.Gap_open, para_info.Gap_ext, para_info.Proc_num, para_info.Debug_mode)

	log.Printf("Prog paras:\tMax_ins=%d, Max_err=%.5f, Mut_rate=%.5f, Err_var_factor=%d, Mut_var_factor=%d, Iter_num_factor=%d, "+
		"Read_len=%d, Info_len=%d, Seed_backup=%d, Ham_backup=%d, Indel_backup=%d", para_info.Max_ins, para_info.Err_rate, para_info.Mut_rate, 
		para_info.Err_var_factor, para_info.Mut_var_factor, para_info.Iter_num_factor, para_info.Read_len, para_info.Info_len,
		para_info.Seed_backup, para_info.Ham_backup, para_info.Indel_backup)

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
// Information of seeds between reads and the multigenome.
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
	l_Dist_D, l_Dist_IS, l_Dist_IT    [][]float64 // Distance matrix for backward alignment
	l_Trace_D, l_Trace_IS, l_Trace_IT [][][]int   // Backtrace matrix for backward alignment
	r_Dist_D, r_Dist_IS, r_Dist_IT    [][]float64 // Distance matrix for forward alignment
	r_Trace_D, r_Trace_IS, r_Trace_IT [][][]int   // Backtrace matrix for forward alignment
}

//--------------------------------------------------------------------------------------------------
// InitEditAlnInfo allocates memory for share variables for alignment process
//--------------------------------------------------------------------------------------------------
func InitEditAlnInfo(arr_len int) *EditAlnInfo {
	aln_info := new(EditAlnInfo)
	aln_info.l_Dist_D, aln_info.l_Trace_D = InitEditAlnMat(arr_len)
	aln_info.l_Dist_IS, aln_info.l_Trace_IS = InitEditAlnMat(arr_len)
	aln_info.l_Dist_IT, aln_info.l_Trace_IT = InitEditAlnMat(arr_len)
	aln_info.r_Dist_D, aln_info.r_Trace_D = InitEditAlnMat(arr_len)
	aln_info.r_Dist_IS, aln_info.r_Trace_IS = InitEditAlnMat(arr_len)
	aln_info.r_Dist_IT, aln_info.r_Trace_IT = InitEditAlnMat(arr_len)
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

//---------------------------------------------------------------------------------------------------
// Information of unaligned reads.
//---------------------------------------------------------------------------------------------------
type UnAlnInfo struct {
	read_info1, read_info2 []byte //unalgined read info.
}

//--------------------------------------------------------------------------------------------------
// SplitN splits a slice of bytes using a memory-efficient method.
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
