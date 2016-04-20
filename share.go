//---------------------------------------------------------------------------------------------------
// IVC: share.go
// Shared variables and functions for IVC modules.
// Copyright 2015 Nam Sy Vo.
//---------------------------------------------------------------------------------------------------

package ivc

import (
	"bufio"
	"bytes"
	"log"
	"os"
	"path"
	"runtime"
	"runtime/pprof"
	"sync"
)

//--------------------------------------------------------------------------------------------------
// Global constants
//--------------------------------------------------------------------------------------------------
const (
	NEW_SNP_RATE   = 0.001  // probability of new alleles
	NEW_INDEL_RATE = 0.0001 // probability of new indels
	INDEL_ERR_RATE = 0.0001 // probability of indel error
)

//--------------------------------------------------------------------------------------------------
// Global variables for calculating variant quality.
//--------------------------------------------------------------------------------------------------
var (
	PARA *ParaInfo        // all parameters of the program
	L2E  []float64        // indel error rate corresponding to lengths of indels
	Q2C  map[byte]float64 // alignment cost based on Phred-scale quality
	Q2E  map[byte]float64 // error probability based on Phred-scale quality
	Q2P  map[byte]float64 // non-error probability based on Phred-scale quality
	MUT  = &sync.Mutex{}  // mutex lock for reading/writing from/to the map of variant calls
)

//--------------------------------------------------------------------------------------------------
// Parameter information
//--------------------------------------------------------------------------------------------------
type ParaInfo struct {
	//Input file names:
	Ref_file       string // reference multigenome
	Var_prof_file  string // variant profile
	Index_file     string // index of original reference genomes
	Rev_index_file string // index of reverse reference genomes
	Read_file_1    string // first end of read
	Read_file_2    string // second end of read
	Var_call_file  string // store Var call

	// Input paras:
	Search_mode int     // searching mode for finding seeds
	Start_pos   int     // starting postion on reads for finding seeds
	Search_step int     // step for searching in deterministic mode
	Max_snum    int     // maximum number of seeds
	Max_psnum   int     // maximum number of paired-seeds
	Min_slen    int     // minimum length of seeds
	Max_slen    int     // maximum length of seeds
	Dist_thres  float64 // threshold for distances between reads and multigenomes
	Iter_num    int     // number of random iterations to find proper alignments
	Sub_cost    float64 // cost of substitution for Hamming and Edit distance
	Gap_open    float64 // cost of gap open for Edit distance
	Gap_ext     float64 // cost of gap extension for Edit distance
	Proc_num    int     // maximum number of CPUs using by Go
	Debug_mode  bool    // debug mode for output

	// Estimated paras:
	Read_len        int     // read length, calculated from read files
	Info_len        int     // maximum size of array to store read headers
	Max_ins         int     // maximum insert size of two aligned ends
	Err_rate        float32 // average sequencing error rate, estmated from reads with real reads
	Err_var_factor  int     // factor for standard variation of sequencing error rate
	Mut_rate        float32 // average mutation rate, estmated from reference genome
	Mut_var_factor  int     // factor for standard variation of mutation rate
	Iter_num_factor int     // factor for number of iterations
	Seed_backup     int     // number of backup bases from seeds
	Ham_backup      int     // number of backup bases from Hamming alignment
	Indel_backup    int     // number of backup bases from known indels
}

//--------------------------------------------------------------------------------------------------
// Read input information and set up parameters
//--------------------------------------------------------------------------------------------------
func Setup(input_para *ParaInfo) {

	log.Printf("----------------------------------------------------------------------------------------")
	log.Printf("Checking input information and seting up parameters...")

	//Check input files
	var f *os.File
	var e error
	if _, e = os.Stat(input_para.Ref_file); e != nil {
		log.Panicf("Error: %s", e)
	}
	if _, e = os.Stat(input_para.Var_prof_file); e != nil {
		log.Panicf("Error: %s", e)
	}
	if _, e = os.Stat(input_para.Rev_index_file); e != nil {
		log.Panicf("Error: %s", e)
	}
	if _, e = os.Stat(input_para.Read_file_1); e != nil {
		log.Panicf("Error: %s", e)
	}
	if _, e = os.Stat(input_para.Read_file_2); e != nil {
		log.Panicf("Error: %s", e)
	}
	PARA = SetupPara(input_para)

	if PARA.Debug_mode {
		MEM_STATS = new(runtime.MemStats)
		if CPU_FILE, e = os.Create(PARA.Var_call_file + ".cprof"); e != nil {
			log.Panicf("Error: %s", e)
		}
		pprof.StartCPUProfile(CPU_FILE)

		if MEM_FILE, e = os.Create(PARA.Var_call_file + ".mprof"); e != nil {
			log.Panicf("Error: %s", e)
		}
		log.Printf("Debug mode:\tCpu_prof_file: %s, Mem_prof_file: %s", PARA.Var_call_file+".cprof", PARA.Var_call_file+".mprof")
	}

	result_dir := path.Dir(PARA.Var_call_file)
	if _, e = os.Stat(result_dir); e != nil {
		if os.IsNotExist(e) {
			if e = os.Mkdir(result_dir, 0777); e != nil {
				log.Panicf("Error: %s", e)
			}
		} else {
			log.Panicf("Error: %s", e)
		}
	}
	if f, e = os.Create(PARA.Var_call_file); e != nil {
		log.Panicf("Error: %s", e)
	}
	w := bufio.NewWriter(f)
	if PARA.Debug_mode == false {
		w.WriteString("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" +
			"VAR_PROB\tMAP_PROB\tCOM_QUAL\tVAR_NUM\tREAD_NUM\n")
	} else {
		w.WriteString("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" +
			"VAR_PROB\tMAP_PROB\tCOM_QUAL\tVAR_NUM\tREAD_NUM\tBASE_QUAL\tCHR_DIS\tCHR_DIFF\tMAP_PROB\t" +
			"ALN_PROB\tPAIR_PROB\tS_POS1\tBRANCH1\tS_POS2\tBRANCH2\tREAD_HEADER\tALN_BASE\tBASE_NUM\n")
	}
	w.Flush()
	f.Close()

	log.Printf("Finish checking input information and seting up parameters.")
}

//--------------------------------------------------------------------------------------------------
// SetupPara setups values of parameters for alignment process
//--------------------------------------------------------------------------------------------------
func SetupPara(input_para *ParaInfo) *ParaInfo {

	para := input_para

	f, e := os.Open(para.Read_file_1)
	if e != nil {
		log.Panicf("Error: %s", e)
	}
	s := bufio.NewScanner(f)
	s.Scan()
	header := s.Bytes()
	if len(header) > 0 {
		para.Info_len = len(header) + 20 //there might be longer header, is that case, ignore the longer part
	} else {
		para.Info_len = 100
		log.Printf("Possibly missing header")
	}
	s.Scan()
	read := s.Bytes()
	if len(read) > 0 {
		para.Read_len = len(read)
	} else {
		log.Panicf("Something is wrong with input read sequence.")
	}
	f.Close()

	// 1500 is asigned based on insert size of paired-end testing reads
	// will be estimated based on input reads (= 3*avg_ins_size)
	para.Max_ins = 1500

	// 0.0015 is maximum sequencing error rate of testing reads, 0.01 is mutation rate of testing data,
	// will be set up based on input reads
	para.Err_rate = 0.0015
	para.Mut_rate = 0.01

	para.Err_var_factor = 4
	para.Mut_var_factor = 2
	para.Iter_num_factor = 2

	para.Seed_backup = 10
	para.Ham_backup = 15
	para.Indel_backup = 30

	// Setup input parameters if not specified
	if input_para.Search_mode == 0 {
		para.Search_mode = 1
		log.Printf("No or invalid input for searching mode, use default strategy (randomizaion).")
	} else if input_para.Search_mode == 1 {
		if input_para.Start_pos == 0 {
			para.Start_pos = 5
			log.Printf("Deterministic search mode: no or invalid input for start postion on reads to find seeds, use default value (%d).", para.Start_pos)
		}
		if input_para.Search_step == 0 {
			para.Search_step = 5
			log.Printf("Deterministic search mode: no or invalid input for searching step, use default value (%d).", para.Search_step)
		}
	}
	if input_para.Max_snum == 0 {
		para.Max_snum = 4096
		log.Printf("No or invalid input for maximum number of seeds, use default value (%d).", para.Max_snum)
	}
	if input_para.Max_psnum == 0 {
		para.Max_psnum = 128
		log.Printf("No or invalid input for maximum number of paired-seeds, use default value (%d).", para.Max_psnum)
	}
	if input_para.Min_slen == 0 {
		para.Min_slen = 15
		log.Printf("No or invalid input for minimum length of seeds, use default value (%d).", para.Min_slen)
	}
	if input_para.Max_slen == 0 {
		para.Max_slen = 25
		log.Printf("No or invalid input for maximum length of seeds, use default value (%d).", para.Max_slen)
	}
	if input_para.Sub_cost == 0 {
		para.Sub_cost = 4
		log.Printf("No or invalid input for substitution cost of alignment, use default value (%.1f).", para.Sub_cost)
	}
	if input_para.Gap_open == 0 {
		para.Gap_open = 4.1
		log.Printf("No or invalid input for gap open cost of alignment, use default value (%.1f).", para.Gap_open)
	}
	if input_para.Gap_ext == 0 {
		para.Gap_ext = 1
		log.Printf("No or invalid input for gap extension cost of alignment, use default value (%.1f).", para.Gap_ext)
	}

	if input_para.Dist_thres == 0 {
		/*
			err := float64(para.Err_rate)
			rlen := float64(para.Read_len)
			mut := float64(para.Mut_rate)
			k1 := float64(para.Err_var_factor)
			k2 := float64(para.Mut_var_factor)
			var_dist = int(math.Ceil(err*rlen+k1*math.Sqrt(rlen*err*(1-err)))) + int(math.Ceil(mut*rlen+k2*math.Sqrt(rlen*mut*(1-mut))))
			para.Dist_thres = -float64(var_dist)*math.Log10(1-err) - float64(var_dist)*math.Log10(NEW_INDEL_RATE)
		*/
		para.Dist_thres = 36
		log.Printf("No or invalid input for threshold of alignment distance, calculate based on input data (%.1f).", para.Dist_thres)
	}
	if input_para.Iter_num == 0 {
		//para.Iter_num = para.Iter_num_factor * (para.Dist_thres + 1)
		para.Iter_num = 12
		log.Printf("No or invalid input for numbers of random iterations, calculate based on input data (%d).", para.Iter_num)
	}

	if input_para.Proc_num == 0 {
		para.Proc_num = runtime.NumCPU()
		log.Printf("No or invalid input for number of threads, use maximum number of CPUs of the current machine (%d).", para.Proc_num)
	}

	log.Printf("Input files:\tGenome_file: %s, Var_file: %s, Index_file=%s, Read_file_1=%s, Read_file_2=%s, Var_call_file=%s",
		para.Ref_file, para.Var_prof_file, para.Rev_index_file, para.Read_file_1, para.Read_file_2, para.Var_call_file)

	log.Printf("Input paras:\tSearch_mode=%d, Start_pos=%d, Search_step=%d, Max_snum=%d, Max_psnum=%d, "+
		"Min_slen=%d, Max_slen=%d, Dist_thres=%.1f, Iter_num=%d, Sub_cost=%.1f, Gap_open=%.1f, Gap_ext=%.1f, Proc_num=%d, Debug_mode=%t",
		para.Search_mode, para.Start_pos, para.Search_step, para.Max_snum, para.Max_psnum, para.Min_slen, para.Max_slen,
		para.Dist_thres, para.Iter_num, para.Sub_cost, para.Gap_open, para.Gap_ext, para.Proc_num, para.Debug_mode)

	log.Printf("Prog paras:\tMax_ins=%d, Max_err=%.5f, Mut_rate=%.5f, Err_var_factor=%d, Mut_var_factor=%d, Iter_num_factor=%d, "+
		"Read_len=%d, Info_len=%d, Seed_backup=%d, Ham_backup=%d, Indel_backup=%d", para.Max_ins, para.Err_rate, para.Mut_rate,
		para.Err_var_factor, para.Mut_var_factor, para.Iter_num_factor, para.Read_len, para.Info_len,
		para.Seed_backup, para.Ham_backup, para.Indel_backup)

	return para
}

//--------------------------------------------------------------------------------------------------
// Information of input reads
//--------------------------------------------------------------------------------------------------
type ReadInfo struct {
	Read1, Read2                   []byte // first and second ends
	Qual1, Qual2                   []byte // quality info of the first read and second ends
	Rev_read1, Rev_read2           []byte // reverse of the first and second ends
	Rev_comp_read1, Rev_comp_read2 []byte // reverse complement of the first and second ends
	Comp_read1, Comp_read2         []byte // complement of the first and second ends
	Rev_qual1, Rev_qual2           []byte // quality of reverse of the first and second ends
	Info1, Info2                   []byte // info of the first and second ends
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
func RevComp(read, qual []byte, rev_comp_read, rev_qual []byte) {
	read_len := len(read)
	for i, elem := range read {
		rev_qual[i] = qual[read_len-1-i]
		if elem == 'A' {
			rev_comp_read[read_len-1-i] = 'T'
		} else if elem == 'T' {
			rev_comp_read[read_len-1-i] = 'A'
		} else if elem == 'C' {
			rev_comp_read[read_len-1-i] = 'G'
		} else if elem == 'G' {
			rev_comp_read[read_len-1-i] = 'C'
		} else {
			rev_comp_read[read_len-1-i] = elem
		}
	}
}

//---------------------------------------------------------------------------------------------------
// Information of seeds between reads and the multigenome.
//---------------------------------------------------------------------------------------------------
type SeedInfo struct {
	s_pos  []int  // staring position of seeds on reads
	e_pos  []int  // ending position of seeds on reads
	m_pos  []int  // (left-most) matching position of seeds on the reference multigenome
	strand []bool // strand (forward or reverse) of matches on the reference multigenome
}

//--------------------------------------------------------------------------------------------------
// Alignment information, served as shared variables between functions for alignment process
//--------------------------------------------------------------------------------------------------
type EditAlnInfo struct {
	l_Trace_K, r_Trace_K              [][][]byte  //backtrace matrix for known locations
	l_Dist_D, l_Dist_IS, l_Dist_IT    [][]float64 // distance matrix for backward alignment
	l_Trace_D, l_Trace_IS, l_Trace_IT [][][]int   // backtrace matrix for backward alignment
	r_Dist_D, r_Dist_IS, r_Dist_IT    [][]float64 // distance matrix for forward alignment
	r_Trace_D, r_Trace_IS, r_Trace_IT [][][]int   // backtrace matrix for forward alignment
}

//--------------------------------------------------------------------------------------------------
// InitEditAlnInfo allocates memory for share variables for alignment process
//--------------------------------------------------------------------------------------------------
func InitEditAlnInfo(arr_len int) *EditAlnInfo {
	aln_info := new(EditAlnInfo)
	aln_info.l_Trace_K, aln_info.r_Trace_K = InitTraceKMat(arr_len), InitTraceKMat(arr_len)
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
func InitTraceKMat(arr_len int) [][][]byte {
	trace_mat := make([][][]byte, arr_len+1)
	for i := 0; i <= arr_len; i++ {
		trace_mat[i] = make([][]byte, arr_len+1)
		for j := 0; j <= arr_len; j++ {
			trace_mat[i][j] = make([]byte, 0)
		}
	}
	return trace_mat
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
			trace_mat[i][j] = make([]int, 2)
		}
	}
	return dis_mat, trace_mat
}

//---------------------------------------------------------------------------------------------------
// Information of unaligned reads.
//---------------------------------------------------------------------------------------------------
type UnAlnInfo struct {
	read_info1, read_info2 []byte //unalgined read info
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
