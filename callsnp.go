//---------------------------------------------------------------------------------------------------
// Calling SNPs based on read-multigenome alignment results.
// Copyright 2014 Nam Sy Vo
//---------------------------------------------------------------------------------------------------

package isc

import (
	"fmt"
	"bufio"
	"log"
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

//Aligned base at each position on reference multigenome
type BaseInfo struct {
	Pos uint32
	Base byte
	Qual uint16
}

//SNP info at each position on reference multigenome
type SNP struct {
	SNP_Idx int
	SNP_Val []byte
}

//SNP call at each position on reference multigenome
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

//Initialize parameters
func (S *SNPProf) InitIndex(input_info InputInfo, para_info ParaInfo) {

	INPUT_INFO = input_info
	PARA_INFO = para_info
	RAND_GEN = rand.New(rand.NewSource(time.Now().UnixNano()))
	INDEX.Init(input_info, para_info)

	S.SNP_Prof = make(map[int][][]byte)
	S.SNP_Conf = make(map[int][]float32)
	for k, v := range INDEX.SNP_AF {
		S.SNP_Conf[k] = v
	}

	S.SNP_Call = make(map[int][]byte)
	S.SNP_Qual = make(map[int]float32)
}

//--------------------------------------------------------------------------------------------------
//Read input FASTQ files and put data into data channel
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
	var line_f1, line_f2 []byte
	for scanner1.Scan() && scanner2.Scan() { //ignore 1st lines in input FASTQ files
		scanner1.Scan()
		scanner2.Scan()
		line_f1 = scanner1.Bytes() //use 2nd line in input FASTQ file 1
		line_f2 = scanner2.Bytes() //use 2nd line in input FASTQ file 2
		if len(line_f1) > 0 && len(line_f2) > 0 {
			read_num++
			read_info.AssignReads(line_f1, line_f2)
			read_info.CalcRevComp()
			data <- read_info
		}
		if read_num%10000 == 0 {
			runtime.ReadMemStats(memstats)
			log.Printf("isc.go: memstats after aligning each 10,000 reads:\t%d\t%d\t%d\t%d\t%d",
				memstats.Alloc,	memstats.TotalAlloc, memstats.Sys, memstats.HeapAlloc, memstats.HeapSys)
		}
		scanner1.Scan() //ignore 3rd line in 1st input FASTQ file 1
		scanner2.Scan() //ignore 3rd line in 2nd input FASTQ file 2
		scanner1.Scan()	//ignore 4th line in 1st input FASTQ file 1
		scanner2.Scan()	//ignore 4th line in 2nd input FASTQ file 2
	}
	close(data)
}

//--------------------------------------------------------------------------------------------------
//Take data from data channel, process them (find SNPs) and put results (SNPs) into results channel
//--------------------------------------------------------------------------------------------------
func (S *SNPProf) ProcessReads(data chan ReadInfo, results chan []SNP, wg *sync.WaitGroup, align_info AlignInfo, match_pos []int) {

	wg.Add(1)
	defer wg.Done()
	var read_info ReadInfo
	var SNPs []SNP
	for read_info = range data {
		SNPs = S.FindSNP(read_info, align_info, match_pos)
		if len(SNPs) > 0 {
			results <- SNPs
		}
	}
}

//--------------------------------------------------------------------------------------------------
//Align reads to the reference multigenome
//--------------------------------------------------------------------------------------------------
func (S *SNPProf) AlignReads() uint64 {

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
		go S.ProcessReads(data, results, &wg, align_info[i], match_pos[i])
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
// FindSNPProfile returns SNP profile of new genome based on SNP profile of reference multi-genomes
// and alignment between reads and multi-genomes.
//
//	- Find SNPs for pairend reads, treat each end separately and independently.
//--------------------------------------------------------------------------------------------------
func (S *SNPProf) FindSNP(read_info ReadInfo, align_info AlignInfo, match_pos []int) []SNP {

	var has_seeds bool
	var p, s_pos, e_pos int
	var loop_num, match_num int
	var SNPs, snps []SNP

	//Find SNPs for the first end
	p = INPUT_INFO.Start_pos
	loop_num = 1
	for loop_num <= ITER_NUM {
		//fmt.Println(loop_num, "\tread1")
		s_pos, e_pos, match_num, has_seeds = INDEX.FindSeeds(read_info.Read1, read_info.Rev_read1, p, match_pos)
		if has_seeds {
			//fmt.Println("read1, has seed\t", s_pos, "\t", e_pos, "\t", string(read_info.Read1))
			snps = S.FindSNPFromMatch(read_info.Read1, s_pos, e_pos, match_pos, match_num, align_info)
			if len(snps) > 0 {
				//fmt.Println("read1, has snp\t", s_pos, "\t", e_pos, "\t", string(read_info.Read1))
				//fmt.Println(loop_num, "\tori1\t", string(read1))
				SNPs = append(SNPs, snps...)
				break
			}
		}
		//Find SNPs for the reverse complement of the first end
		s_pos, e_pos, match_num, has_seeds = INDEX.FindSeeds(read_info.Rev_comp_read1, read_info.Comp_read1, p, match_pos)
		if has_seeds {
			//fmt.Println("rc_read1, has seed\t", s_pos, "\t", e_pos, "\t", string(read_info.Rev_comp_read1))
			snps = S.FindSNPFromMatch(read_info.Rev_comp_read1, s_pos, e_pos, match_pos, match_num, align_info)
			if len(snps) > 0 {
				//fmt.Println("rc_read1, has snp\t", s_pos, "\t", e_pos, "\t", string(read_info.Rev_comp_read1))
				//fmt.Println(loop_num, "\trev1\t", string(rev_read1))
				SNPs = append(SNPs, snps...)
				break
			}
		}
		//Take a new position to search
		if INPUT_INFO.Search_mode == 1 {
			p = RAND_GEN.Intn(read_info.Read_len-1) + 1
		} else if INPUT_INFO.Search_mode == 2 {
			p = p + INPUT_INFO.Search_step
		}
		loop_num++
	}

	//Find SNPs for the second end
	p = INPUT_INFO.Start_pos
	loop_num = 1
	for loop_num <= ITER_NUM {
		//fmt.Println(loop_num, "\tread2")
		s_pos, e_pos, match_num, has_seeds = INDEX.FindSeeds(read_info.Read2, read_info.Rev_read2, p, match_pos)
		if has_seeds {
			//fmt.Println("read2, has seed\t", s_pos, "\t", e_pos, "\t", string(read_info.Read2))
			snps = S.FindSNPFromMatch(read_info.Read2, s_pos, e_pos, match_pos, match_num, align_info)
			if len(snps) > 0 {
				//fmt.Println("read2, has snp\t", s_pos, "\t", e_pos, "\t", string(read_info.Read2))
				//fmt.Println(loop_num, "\tori2\t", string(read2))
				SNPs = append(SNPs, snps...)
				return SNPs
			}
		}
		//Find SNPs for the reverse complement of the second end
		s_pos, e_pos, match_num, has_seeds = INDEX.FindSeeds(read_info.Rev_comp_read2, read_info.Comp_read2, p, match_pos)
		if has_seeds {
			//fmt.Println("rc_read2, has seed\t", s_pos, "\t", e_pos, "\t", string(read_info.Rev_comp_read2))
			snps = S.FindSNPFromMatch(read_info.Rev_comp_read2, s_pos, e_pos, match_pos, match_num, align_info)
			if len(snps) > 0 {
				//fmt.Println("rc_read2, has snp\t", s_pos, "\t", e_pos, "\t", string(read_info.Rev_comp_read2))
				//fmt.Println(loop_num, "\trev2\t", string(rev_read2))
				SNPs = append(SNPs, snps...)
				return SNPs
			}
		}
		//Take a new position to search
		if INPUT_INFO.Search_mode == 1 {
			p = RAND_GEN.Intn(read_info.Read_len-1) + 1
		} else if INPUT_INFO.Search_mode == 2 {
			p = p + INPUT_INFO.Search_step
		}
		loop_num++
	}
	return SNPs
}

//---------------------------------------------------------------------------------------------------
// UpdateSNPCall updates SNP Call found from alignment between reads and multi-genomes.
//---------------------------------------------------------------------------------------------------
func (S *SNPProf) FindSNPFromMatch(read []byte, s_pos, e_pos int, match_pos []int, match_num int, align_info AlignInfo) []SNP {

	var pos, dis, k int
	var left_snp_idx, right_snp_idx []int
	var left_snp_val, right_snp_val [][]byte
	var snp SNP
	var snps []SNP

	for i := 0; i < match_num; i++ {
		pos = match_pos[i]
		//Call IntervalHasSNP to determine whether extension is needed
		//if index.IntervalHasSNP(A.SORTED_SNP_POS, pos - e_pos, pos - e_pos + len(read1)) {
		dis, left_snp_idx, left_snp_val, right_snp_idx, right_snp_val = INDEX.FindExtensions(read, s_pos, e_pos, pos, align_info)
		if dis <= DIST_THRES {
			//Determine SNP profile
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
		//}
	}
	return snps
}

//---------------------------------------------------------------------------------------------------
// CalcQual calculates called SNP quality.
//---------------------------------------------------------------------------------------------------
func (S *SNPProf) CalcSNPQual(snp_pos int) CalledSNP {

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
// CallSNP returns called SNPs and related information based on SNP profile constructed from
// alignment between reads and multi-genomes.
//-----------------------------------------------------------------------------------------------------
func (S *SNPProf) CallSNP() {

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
				SNP_call_chan <- S.CalcSNPQual(snp_pos)
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
// SNPCall_tofile writes called SNPs and related information to given output file in tab-delimited format
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
