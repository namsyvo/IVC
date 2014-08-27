//---------------------------------------------------------------------------------------------------
// Calling SNPs based on read-multigenome alignment results.
// Copyright 2014 Nam Sy Vo
//---------------------------------------------------------------------------------------------------

package isc

import (
	"fmt"
	"bufio"
	"log"
	"math"
	"math/rand"
	"os"
	"runtime"
	"strconv"
	"time"
	"sort"
	"sync"
)

//Global variables
var (
	INPUT_INFO	InputInfo	//Input information
	PARA_INFO	ParaInfo	//Parameters
	RAND_GEN	*rand.Rand 	//Pseudo-random number generator
	INDEX		Index      	//Index for alignment
)

//BaseInfo stores aligned bases at each position on reference multigenome
type BaseInfo struct {
	Pos uint32
	Base byte
	Qual uint16
}

//SNP stores SNP info at each position on reference multigenome
type SNP struct {
	SNP_Idx int
	SNP_Val []byte
}

//CalledSNP stores SNP calls at each position on reference multigenome
type CalledSNP struct {
	SNPVal SNP		//SNP info
	SNPQual float32 //SNP quality
}

//SNPProf stores SNP Calling information and defines SNP calling functions
type SNPProf struct {
	SNP_Prof map[int][][]byte  // to store all possible SNPs at each position
	SNP_Conf map[int][]float32 // to store quality of all possible SNPS at each position
	SNP_Call map[int][]byte    // to store SNP call at each position
	SNP_Qual map[int]float32     // to store quality score of called SNP at each position
}

//--------------------------------------------------------------------------------------------------
//InitIndex initializes indexes and parameters
//--------------------------------------------------------------------------------------------------
func (S *SNPProf) Init(input_info InputInfo) {

	INPUT_INFO = input_info
	PARA_INFO = SetPara(100, 0.01)
	RAND_GEN = rand.New(rand.NewSource(time.Now().UnixNano()))
	INDEX.Init()

	S.SNP_Prof = make(map[int][][]byte)
	S.SNP_Conf = make(map[int][]float32)
	for k, v := range INDEX.SNP_AF {
		S.SNP_Conf[k] = v
	}

	S.SNP_Call = make(map[int][]byte)
	S.SNP_Qual = make(map[int]float32)
}

//--------------------------------------------------------------------------------------------------
//Set values for parameters
//--------------------------------------------------------------------------------------------------
func SetPara(read_len int, seq_err float32) ParaInfo {
	para_info := ParaInfo{}
	para_info.Max_match = 32
	para_info.Err_var_factor = 4
	para_info.Iter_num_factor = 1
	para_info.Seq_err = seq_err //will be replaced by seq_err estimated from input reads
	para_info.Read_len = read_len //will be replaced by read length taken from input reads

	//Const for computing distance
	err := float64(para_info.Seq_err)
	rlen := float64(para_info.Read_len)
	k := float64(para_info.Err_var_factor)
	para_info.Dist_thres = int(math.Ceil(err*rlen + k*math.Sqrt(rlen*err*(1-err))))
	para_info.Iter_num = para_info.Iter_num_factor * (para_info.Dist_thres + 1)

	fmt.Println("DIST_THRES: ", para_info.Dist_thres)
	fmt.Println("ITER_NUM: ", para_info.Iter_num)

	return para_info
}

//--------------------------------------------------------------------------------------------------
//ProcessReads finds all possible SNPs from read-multigenome alignment.
//--------------------------------------------------------------------------------------------------
func (S *SNPProf) ProcessReads() uint64 {

	align_info := make([]AlignInfo, INPUT_INFO.Routine_num)
	for i := 0; i < INPUT_INFO.Routine_num; i++ {
		align_info[i].InitAlignInfo(PARA_INFO.Read_len)
	}
	match_pos := make([][]int, INPUT_INFO.Routine_num)
	for i := 0; i < INPUT_INFO.Routine_num; i++ {
		match_pos[i] = make([]int, PARA_INFO.Max_match)
	}

	data := make(chan ReadInfo, INPUT_INFO.Routine_num)
	go S.ReadReads(data)

	results := make(chan []SNP)
	var wg sync.WaitGroup
	for i := 0; i < INPUT_INFO.Routine_num; i++ {
		go S.FindSNPs(data, results, &wg, align_info[i], match_pos[i])
	}
	go func() {
		wg.Wait()
		close(results)
	}()

	//Collect SNPS from results channel and update SNPs
	var has_snp_read_num uint64 = 0
	var snp SNP
	for SNPs := range results {
		has_snp_read_num++
		for _, snp = range SNPs {
			S.SNP_Prof[snp.SNP_Idx] = append(S.SNP_Prof[snp.SNP_Idx], snp.SNP_Val)
		}
	}
	return has_snp_read_num
}

//--------------------------------------------------------------------------------------------------
//ReadReads reads all reads from input FASTQ files and put them into data channel
//--------------------------------------------------------------------------------------------------
func (S *SNPProf) ReadReads(data chan ReadInfo) {

	memstats := new(runtime.MemStats)

	fn1, fn2 := INPUT_INFO.Read_file_1, INPUT_INFO.Read_file_2
	f1, err_f1 := os.Open(fn1)
	if err_f1 != nil {
		panic("Error opening input read file " + fn1)
	}
	defer f1.Close()
	f2, err_f2 := os.Open(fn2)
	if err_f2 != nil {
		panic("Error opening input read file " + fn2)
	}
	defer f2.Close()

	read_info := ReadInfo{}
	read_num := 0
	scanner1 := bufio.NewScanner(f1)
	scanner2 := bufio.NewScanner(f2)
	var read1, read2, qual1, qual2 []byte
	for scanner1.Scan() && scanner2.Scan() { //ignore 1st lines in input FASTQ files
		scanner1.Scan()
		scanner2.Scan()
		read1 = scanner1.Bytes() //use 2nd line in input FASTQ file 1
		read2 = scanner2.Bytes() //use 2nd line in input FASTQ file 2
		scanner1.Scan() //ignore 3rd line in 1st input FASTQ file 1
		scanner2.Scan() //ignore 3rd line in 2nd input FASTQ file 2
		scanner1.Scan()
		scanner2.Scan()
		qual1 = scanner1.Bytes() //use 4th line in input FASTQ file 1
		qual2 = scanner2.Bytes() //use 4th line in input FASTQ file 2
		if len(read1) > 0 && len(read2) > 0 {
			read_num++
			read_info.AssignReads(read1, read2, qual1, qual2)
			data <- read_info
		}
		if read_num%10000 == 0 {
			runtime.ReadMemStats(memstats)
			log.Printf("isc.go: memstats after aligning each 10,000 reads:\t%d\t%d\t%d\t%d\t%d",
				memstats.Alloc,	memstats.TotalAlloc, memstats.Sys, memstats.HeapAlloc, memstats.HeapSys)
		}
	}
	close(data)
}

//--------------------------------------------------------------------------------------------------
//FindSNPs takes data from data channel, find all possible SNPs and put them into results channel
//--------------------------------------------------------------------------------------------------
func (S *SNPProf) FindSNPs(data chan ReadInfo, results chan []SNP, wg *sync.WaitGroup, align_info AlignInfo, match_pos []int) {

	wg.Add(1)
	defer wg.Done()
	var SNPs []SNP
	for read_info := range data {
		SNPs = S.FindSNPsFromReads(read_info, align_info, match_pos)
		if len(SNPs) > 0 {
			results <- SNPs
		}
	}
}

//--------------------------------------------------------------------------------------------------
// FindSNPsFromReads returns SNP profile of new genome based on SNP profile of reference multi-genomes
// and alignment between reads and multi-genomes.
// This version: find SNPs for pairend reads, treat each end separately and independently.
//--------------------------------------------------------------------------------------------------
func (S *SNPProf) FindSNPsFromReads(read_info ReadInfo, align_info AlignInfo, match_pos []int) []SNP {

	var SNPs []SNP
	//Find SNPs for the first end
	SNP1 := S.FindSNPsFromEachEnd(read_info.Read1, align_info, match_pos)

	//Find SNPs for the second end
	SNP2 := S.FindSNPsFromEachEnd(read_info.Read2, align_info, match_pos)

	//Will process constrants of two ends here
	//...

	SNPs = append(SNPs, SNP1...)
	SNPs = append(SNPs, SNP2...)

	return SNPs
}

//---------------------------------------------------------------------------------------------------
// FindSNPsFromEachEnd find SNP Call from matches between read and multigenome.
//---------------------------------------------------------------------------------------------------
func (S *SNPProf) FindSNPsFromEachEnd(read []byte, align_info AlignInfo, match_pos []int) []SNP {
	var has_seeds bool
	var p, s_pos, e_pos int
	var loop_num, match_num int
	var SNPs, snps []SNP

	//Rev_read: reverse of read
	//Rev_comp_read: reverse complement of read
	//Comp_read: complement of read, ~ reverse of reverse complement
	rev_read, rev_comp_read, comp_read := RevComp(read)

	p = INPUT_INFO.Start_pos
	loop_num = 1
	for loop_num <= PARA_INFO.Iter_num {
		//fmt.Println(loop_num, "\tread2")
		s_pos, e_pos, match_num, has_seeds = INDEX.FindSeeds(read, rev_read, p, match_pos)
		if has_seeds {
			//fmt.Println("read2, has seed\t", s_pos, "\t", e_pos, "\t", string(read_info.Read2))
			snps = S.FindSNPsFromMatch(read, s_pos, e_pos, match_pos, match_num, align_info)
			if len(snps) > 0 {
				//fmt.Println("read2, has snp\t", s_pos, "\t", e_pos, "\t", string(read_info.Read2))
				SNPs = append(SNPs, snps...)
				return SNPs
			}
		}
		//Find SNPs for the reverse complement of the second end
		s_pos, e_pos, match_num, has_seeds = INDEX.FindSeeds(rev_comp_read, comp_read, p, match_pos)
		if has_seeds {
			//fmt.Println("rc_read2, has seed\t", s_pos, "\t", e_pos, "\t", string(read_info.Rev_comp_read2))
			snps = S.FindSNPsFromMatch(rev_comp_read, s_pos, e_pos, match_pos, match_num, align_info)
			if len(snps) > 0 {
				//fmt.Println("rc_read2, has snp\t", s_pos, "\t", e_pos, "\t", string(read_info.Rev_comp_read2))
				SNPs = append(SNPs, snps...)
				return SNPs
			}
		}
		//Take a new position to search
		if INPUT_INFO.Search_mode == 1 {
			p = RAND_GEN.Intn(len(read) - 1) + 1
		} else if INPUT_INFO.Search_mode == 2 {
			p = p + INPUT_INFO.Search_step
		}
		loop_num++
	}
	return SNPs
}

//---------------------------------------------------------------------------------------------------
// FindSNPsFromMatch find SNP Call from extensions of matches between read and multigenome.
//---------------------------------------------------------------------------------------------------
func (S *SNPProf) FindSNPsFromMatch(read []byte, s_pos, e_pos int, match_pos []int, match_num int, align_info AlignInfo) []SNP {

	var pos, dis, k int
	var left_snp_idx, right_snp_idx []int
	var left_snp_val, right_snp_val [][]byte
	var snp SNP
	var snps []SNP

	for i := 0; i < match_num; i++ {
		pos = match_pos[i]
		dis, left_snp_idx, left_snp_val, right_snp_idx, right_snp_val = INDEX.FindExtensions(read, s_pos, e_pos, pos, align_info)
		if dis <= PARA_INFO.Dist_thres {
			if len(left_snp_idx) == 0 && len(right_snp_idx) == 0 {
				continue
			} else {
				for k = 0; k < len(left_snp_idx); k++ {
					snp.SNP_Idx, snp.SNP_Val = left_snp_idx[k], left_snp_val[k]
					snps = append(snps, snp)
					//fmt.Println("left snp\t", idx, "\t", string(val))
				}
				for k = 0; k < len(right_snp_idx); k++ {
					snp.SNP_Idx, snp.SNP_Val = right_snp_idx[k], right_snp_val[k]
					snps = append(snps, snp)
					//fmt.Println("right snp\t", idx, "\t", string(val))
				}
			}
		}
	}
	return snps
}

//---------------------------------------------------------------------------------------------------
// CallSNPForEachPos calls SNP and calculates SNP quality for each postion.
//---------------------------------------------------------------------------------------------------
func (S *SNPProf) CallSNPForEachPos(snp_pos int) CalledSNP {

	SNP_Qlt := make(map[string]int)
	for _, snp := range S.SNP_Prof[snp_pos] {
		SNP_Qlt[string(snp)] = SNP_Qlt[string(snp)] + 1
	}

	major_num := 0
	var major_snp string
	for snp_val, snp_num := range SNP_Qlt {
		if snp_num > major_num {
			major_num = snp_num
			major_snp = snp_val
		}
	}

	calledSNP := CalledSNP{}
	calledSNP.SNPVal = SNP{}
	calledSNP.SNPVal.SNP_Idx = snp_pos
	calledSNP.SNPVal.SNP_Val = []byte(major_snp)
	calledSNP.SNPQual = float32(major_num)/float32(len(S.SNP_Prof[snp_pos]))

	return calledSNP

	/*
	cond_prob := make([]float32, len(S.SNP_Conf[snp_pos]))
	for idx, base_conf := range S.SNP_Conf[snp_pos] {
		SNP := index.SNP_PROF[snp_pos][idx]
		for _, snp := range snp_arr {
			if bytes.Equal(snp, SNP) {
				//fmt.Println("Diff: ", string(snp), string(SNP), base_conf)
				cond_prob[idx] = base_conf * float32(1 - 0.01)
			} else {
				//fmt.Println("Same: ", string(snp), string(SNP), base_conf)
				cond_prob[idx] = base_conf * float32(0.01)
			}
		}
	}
	base_prob := float32(0)
	for _, value := range(cond_prob) {
		base_prob += value
	}
	for idx, value := range(cond_prob) {
		S.SNP_Conf[snp_pos][idx] = value / base_prob
	}
	*/

}

//-----------------------------------------------------------------------------------------------------
// CallSNPs returns called SNPs and related information.
//-----------------------------------------------------------------------------------------------------
func (S *SNPProf) CallSNPs() {

	snp_pos_chan := make(chan int, INPUT_INFO.Routine_num)
	go func() {
		for snp_pos, _ := range S.SNP_Prof {
			snp_pos_chan <- snp_pos
		}
		close(snp_pos_chan)
	}()
	/*
	//base_info_chan := make(chan BaseInfo, routine_num)
	go func() {
		base_info := []baseInfo{}
		for snp_pos, snp_val := range S.SNP_Prof {
			base_info.Pos = snp_pos
			base_info.Base = snp_val[i]
			base_info.Qual = S.SNP_Conf[snp_pos]
			base_info_chan <- base_info
		}
		close(base_info_chan)
	}()
	*/
	SNP_call_chan := make(chan CalledSNP)
	var wg sync.WaitGroup
	for i := 0; i < INPUT_INFO.Routine_num; i++ {
		go func() {
			wg.Add(1)
			defer wg.Done()
			for snp_pos := range snp_pos_chan {
				SNP_call_chan <- S.CallSNPForEachPos(snp_pos)
			}
			/*
			for base_info := range base_info_chan {
				SNP_call_chan <- CalcSNPQual(base_info)
			}
			*/
		}()
	}
	go func() {
		wg.Wait()
		close(SNP_call_chan)
	}()

	for snp_call := range SNP_call_chan {
		S.SNP_Call[snp_call.SNPVal.SNP_Idx] = snp_call.SNPVal.SNP_Val
		S.SNP_Qual[snp_call.SNPVal.SNP_Idx] = snp_call.SNPQual
	}
}

//-------------------------------------------------------------------------------------------------------
// SNPCall_tofile writes called SNPs and related information to output file in tab-delimited format
//-------------------------------------------------------------------------------------------------------
func (S *SNPProf) WriteSNPCalls() {

	file, err := os.Create(INPUT_INFO.SNP_call_file)
	if err != nil {
		return
	}
	defer file.Close()

	var snp_pos int
	var str_snp_pos, str_snp_val, str_snp_qual string

	SNP_Pos := make([]int, 0, len(S.SNP_Prof))
	for snp_pos, _ = range S.SNP_Call {
		SNP_Pos = append(SNP_Pos, snp_pos)
	}
	sort.Ints(SNP_Pos)

	for _, snp_pos = range SNP_Pos {
		str_snp_pos = strconv.Itoa(snp_pos)
		str_snp_val = string(S.SNP_Call[snp_pos])
		str_snp_qual = strconv.FormatFloat(float64(S.SNP_Qual[snp_pos]), 'f', 5, 32)

		if str_snp_val != "" {
			_, err = file.WriteString(str_snp_pos + "\t" + str_snp_val + "\t" + str_snp_qual + "\n")
		} else {
			_, err = file.WriteString(str_snp_pos + "\t" + "."		   + "\t" + str_snp_qual + "\n")
		}
		if err != nil {
			fmt.Println(err)
			break
		}
	}
}
