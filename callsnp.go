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
	"bytes"
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
	SNP_Pos int
	SNP_Val []byte
	SNP_Qual float64
}

//CalledSNP stores SNP calls at each position on reference multigenome
type CalledSNP struct {
	SNPVal SNP		//SNP info
	SNPQual float64 //SNP quality
}

//SNPProf stores SNP Calling information and defines SNP calling functions
type SNPProf struct {
	SNP_Conf map[int][]float32 // to store prior probability of all bases at each position

	SNP_Prof map[int][][]byte  // to store all possible SNPs at each position
	SNP_BaseQ map[int][]float64// to store base quality of SNPs at each position

	SNP_Call map[int][]byte    // to store SNP call at each position
	SNP_SnpQ map[int]float64    // to store SNP call quality at each position
}

//--------------------------------------------------------------------------------------------------
//InitIndex initializes indexes and parameters
//--------------------------------------------------------------------------------------------------
func (S *SNPProf) Init(input_info InputInfo) {

	INPUT_INFO = input_info
	PARA_INFO = SetPara(100, 0.01)
	INDEX.Init()

	RAND_GEN = rand.New(rand.NewSource(time.Now().UnixNano()))
	S.SNP_Conf = make(map[int][]float32)
	for k, v := range INDEX.SNP_AF {
		S.SNP_Conf[k] = v
	}
	
	S.SNP_Prof = make(map[int][][]byte)
	S.SNP_BaseQ = make(map[int][]float64)

	S.SNP_Call = make(map[int][]byte)
	S.SNP_SnpQ = make(map[int]float64)
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

	data := make([]chan ReadInfo, INPUT_INFO.Routine_num)
	for i := range data {
		data[i] = make(chan ReadInfo)
	}
	go S.ReadReads(data)

	results := make(chan []SNP)
	var wg sync.WaitGroup
	for i := 0; i < INPUT_INFO.Routine_num; i++ {
		go S.FindSNPs(data[i], results, &wg, align_info[i], match_pos[i])
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
			if len(snp.SNP_Val) > 0 {
				S.SNP_Prof[snp.SNP_Pos] = append(S.SNP_Prof[snp.SNP_Pos], snp.SNP_Val)
				S.SNP_BaseQ[snp.SNP_Pos] = append(S.SNP_BaseQ[snp.SNP_Pos], snp.SNP_Qual)
			}
		}
	}
	return has_snp_read_num
}

//--------------------------------------------------------------------------------------------------
//ReadReads reads all reads from input FASTQ files and put them into data channel
//--------------------------------------------------------------------------------------------------
func (S *SNPProf) ReadReads(data []chan ReadInfo) {

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

	read_info := make([]ReadInfo, INPUT_INFO.Routine_num)
	r_len1, r_len2 := 100, 100
	for i := 0; i < INPUT_INFO.Routine_num; i++ {
		read_info[i].Rev_read1, read_info[i].Rev_comp_read1, read_info[i].Comp_read1 =
			make([]byte, r_len1), make([]byte, r_len1), make([]byte, r_len1)
		read_info[i].Rev_read2, read_info[i].Rev_comp_read2, read_info[i].Comp_read2 =
			make([]byte, r_len2), make([]byte, r_len2), make([]byte, r_len2)
	}
	read_num := 0
	scanner1 := bufio.NewScanner(f1)
	scanner2 := bufio.NewScanner(f2)
	var read1, read2, qual1, qual2 []byte
	Loop:
	for {
		for i := 0; i < INPUT_INFO.Routine_num; i++ {
			if !(scanner1.Scan() && scanner2.Scan()) {
				break Loop
			}
			//ignore 1st lines in input FASTQ files
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
				read_info[i].AssignReads(read1, read2, qual1, qual2)
				data[i] <- read_info[i]
			}
			if read_num%10000 == 0 {
				runtime.ReadMemStats(memstats)
				log.Printf("isc.go: memstats after distributing 10,000 more reads:\t%d\t%d\t%d\t%d\t%d",
					memstats.Alloc,	memstats.TotalAlloc, memstats.Sys, memstats.HeapAlloc, memstats.HeapSys)
			}
		}
	}
	fmt.Println("Number of reads: ", read_num)
	for i := 0; i < INPUT_INFO.Routine_num; i++ {
		close(data[i])
	}
}

//--------------------------------------------------------------------------------------------------
//FindSNPs takes data from data channel, find all possible SNPs and put them into results channel
//--------------------------------------------------------------------------------------------------
func (S *SNPProf) FindSNPs(data chan ReadInfo, results chan []SNP, wg *sync.WaitGroup, align_info AlignInfo, match_pos []int) {

	wg.Add(1)
	defer wg.Done()
	var SNPs []SNP
	for read_info := range data {
		RevComp(read_info.Read1, read_info.Rev_read1, read_info.Rev_comp_read1, read_info.Comp_read1)
		RevComp(read_info.Read2, read_info.Rev_read2, read_info.Rev_comp_read2, read_info.Comp_read2)
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
	SNP1 := S.FindSNPsFromEachEnd(read_info.Read1, read_info.Rev_read1, read_info.Rev_comp_read1, read_info.Comp_read1, align_info, match_pos)

	//Find SNPs for the second end
	SNP2 := S.FindSNPsFromEachEnd(read_info.Read2, read_info.Rev_read2, read_info.Rev_comp_read2, read_info.Comp_read2, align_info, match_pos)

	//Will process constrants of two ends here
	//...

	SNPs = append(SNPs, SNP1...)
	SNPs = append(SNPs, SNP2...)
	return SNPs
}

//---------------------------------------------------------------------------------------------------
// FindSNPsFromEachEnd find SNP Call from matches between read and multigenome.
//---------------------------------------------------------------------------------------------------
func (S *SNPProf) FindSNPsFromEachEnd(read, rev_read, rev_comp_read, comp_read []byte, align_info AlignInfo, match_pos []int) []SNP {
	var has_seeds bool
	var p, s_pos, e_pos int
	var loop_num, match_num int
	var SNPs, snps []SNP

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
	var left_snp_pos, right_snp_pos []int
	var left_snp_val, right_snp_val [][]byte
	var snp SNP
	var snps []SNP
	var q float64

	for i := 0; i < match_num; i++ {
		pos = match_pos[i]
		dis, left_snp_pos, left_snp_val, _, right_snp_pos, right_snp_val, _ =
			 INDEX.FindExtensions(read, s_pos, e_pos, pos, align_info)
		if dis <= PARA_INFO.Dist_thres {
			if len(left_snp_pos) == 0 && len(right_snp_pos) == 0 {
				continue
			} else {
				for k = 0; k < len(left_snp_pos); k++ {
					//q = rev_phred(read[left_snp_idx])
					q = 0.01
					snp.SNP_Pos, snp.SNP_Val, snp.SNP_Qual = left_snp_pos[k], left_snp_val[k], q
					snps = append(snps, snp)
				}
				for k = 0; k < len(right_snp_pos); k++ {
					//q = rev_phred(read[right_snp_idx])
					q = 0.01
					snp.SNP_Pos, snp.SNP_Val, snp.SNP_Qual = right_snp_pos[k], right_snp_val[k], q
					snps = append(snps, snp)
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

	calledSNP := CalledSNP{}
	calledSNP.SNPVal = SNP{}
	calledSNP.SNPVal.SNP_Pos = snp_pos

	has_indel := false
	for _, snp := range S.SNP_Prof[snp_pos] {
		if len(snp) > 1 {
			has_indel = true
		}
	}
	if has_indel {
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
		calledSNP.SNPVal.SNP_Val = []byte(major_snp)
		calledSNP.SNPQual = float64(major_num)/float64(len(S.SNP_Prof[snp_pos]))
		return calledSNP
	} else {
		var a []byte
		for _, a_base := range S.SNP_Prof[snp_pos] {
			//fmt.Println(snp_pos, "\t", string(a_base))
			a = append(a, a_base[0])
		}
		var b []byte
		for _, b_base := range INDEX.SNP_PROF[snp_pos] {
			b = append(b, b_base[0])
		}
		var p_b []float32
		for _, p := range S.SNP_Conf[snp_pos] {
			p_b = append(p_b, p)
		}
		for _, base := range [][]byte{[]byte{'A'}, []byte{'C'}, []byte{'G'}, []byte{'T'}} {
			if !bytes.Contains(b, base) {
				b = append(b, base[0])
				p_b = append(p_b, 0.000001)
			}
		}
		snp_call, snp_qual := CalcSNPQual(a, S.SNP_BaseQ[snp_pos], b, p_b)
		calledSNP.SNPVal.SNP_Val = []byte{snp_call}
		calledSNP.SNPQual = snp_qual
		return calledSNP
	}
}


//---------------------------------------------------------------------------------------------------
// CalcSNPQual calculates SNP quality for bases a and qualities e
//---------------------------------------------------------------------------------------------------
func CalcSNPQual(a []byte, e []float64, bases []byte, p_b []float32) (byte, float64) {
	var p_ab []float64
	for _, b := range(bases) {
		p := 1.0
		for j, w := range(a) {
			if w == b {
				p = p * (1.0 - e[j])
			} else {
				p = p * (e[j]/3)
			}
		}
		p_ab = append(p_ab, p)
	}
	p_a := 0.0
	for i, _ := range(bases) {
		p_a = p_a + float64(p_b[i])*p_ab[i]
	}
	var p_ba []float64
	for i, _ := range(bases) {
		p_ba = append(p_ba, float64(p_b[i])*p_ab[i]/p_a)
	}
	p_max := 0.0
	p_max_idx := 0
	for i, p := range p_ba {
		if p > p_max {
			p_max = p
			p_max_idx = i
		}
	}
	return bases[p_max_idx], -10*math.Log10(1 - p_max)
}

//-----------------------------------------------------------------------------------------------------
// CallSNPs returns called SNPs and related information.
//-----------------------------------------------------------------------------------------------------
func (S *SNPProf) CallSNPs() {
	/*
	for _, snp := range S.SNP_Prof[60800313] {
		println("snp", string(snp))
	}
	*/
	var snp_call CalledSNP
	for snp_pos, _ := range S.SNP_Prof {
		snp_call = S.CallSNPForEachPos(snp_pos)
		S.SNP_Call[snp_call.SNPVal.SNP_Pos] = snp_call.SNPVal.SNP_Val
		S.SNP_SnpQ[snp_call.SNPVal.SNP_Pos] = snp_call.SNPQual
	}
}

/*
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
	SNP_call_chan := make(chan CalledSNP)
	var wg sync.WaitGroup
	for i := 0; i < INPUT_INFO.Routine_num; i++ {
		go func() {
			wg.Add(1)
			defer wg.Done()
			for snp_pos := range snp_pos_chan {
				SNP_call_chan <- S.CallSNPForEachPos(snp_pos)
			}
		}()
	}
	go func() {
		wg.Wait()
		close(SNP_call_chan)
	}()

	for snp_call := range SNP_call_chan {
		S.SNP_Call[snp_call.SNPVal.SNP_Pos] = snp_call.SNPVal.SNP_Val
		S.SNP_SnpQ[snp_call.SNPVal.SNP_Pos] = snp_call.SNPQual
	}
}
*/
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
		str_snp_qual = strconv.FormatFloat(float64(S.SNP_SnpQ[snp_pos]), 'f', 5, 32)

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
