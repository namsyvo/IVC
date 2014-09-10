//---------------------------------------------------------------------------------------------------
// Calling SNPs based on read-multigenome alignment.
// Copyright 2014 Nam Sy Vo
//---------------------------------------------------------------------------------------------------

package isc

import (
	"fmt"
	"bufio"
	"math"
	"math/rand"
	"os"
	"strconv"
	"bytes"
	"time"
	"sort"
	"sync"
)

//--------------------------------------------------------------------------------------------------
// Global variables for alignment and SNP calling process.
//--------------------------------------------------------------------------------------------------
var (
	INPUT_INFO	InputInfo	//Input information
	PARA_INFO	ParaInfo	//Parameters information
	RAND_GEN	*rand.Rand 	//Pseudo-random number generator
	INDEX		Index      	//Index for alignment
)

//--------------------------------------------------------------------------------------------------
// SNP represents SNP obtained during alignment phase.
// Each variable of type SNP present SNPs at each position on reference multigenome (temporary variable).
//--------------------------------------------------------------------------------------------------
type SNP struct {
	Pos 	uint32 //SNP postion on ref
	Bases 	[]byte //bases of SNP
	BaseQ 	[]byte //quality of bases of SNP
}

//--------------------------------------------------------------------------------------------------
// SNP_Call represents SNP call obtained during SNP calling phase
// Each variable of type SNP_Call present SNPs at each position on reference multigenome (temporary variable).
//--------------------------------------------------------------------------------------------------
type SNP_Call struct {
	Pos 	uint32 //SNP postion on ref
	Bases 	[]byte //bases of SNP
	Qual 	float32 //quality of SNP
}

//--------------------------------------------------------------------------------------------------
// SNP_Prof represents SNP calls at all positions on reference multigenome.
// One variable of type SNP_Prof is created in initialization phase and stores SNP calls information
// through whole program, i.e., all phases alignment, calling SNPs, and writing SNP calls to file (permanent varialbe).
// This struct also has functions defined on it for calling SNPs.
//--------------------------------------------------------------------------------------------------
type SNP_Prof struct {
	//Prior information, obtained from ref data, without incorporating info from input reads/alignment
	pProb 			map[int][]float32 //to store prior probability of all SNPs at each position

	//Poterior information, obtained from both ref data and input reads/alignment
	SNP_Bases 		map[int][][]byte //to store all possible SNPs at each position
	SNP_BaseQ 		map[int][][]byte //to store base quality of SNPs at each position

	SNP_Call_Bases 	map[int][]byte //to store SNP call at each position
	SNP_Call_Qual 	map[int]float32 //to store SNP call quality at each position
}


////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Functions for initialization phase.
//
////////////////////////////////////////////////////////////////////////////////////////////////////


//--------------------------------------------------------------------------------------------------
// InitIndex initializes indexes and parameters.
// This function will be called from main program.
//--------------------------------------------------------------------------------------------------
func (S *SNP_Prof) Init(input_info InputInfo) {

	INPUT_INFO = input_info
	PARA_INFO = SetPara(100, 0.01)
	INDEX.Init()

	RAND_GEN = rand.New(rand.NewSource(time.Now().UnixNano()))
	S.pProb = make(map[int][]float32)
	for k, v := range INDEX.SNP_AF {
		S.pProb[k] = v
	}
	
	S.SNP_Bases = make(map[int][][]byte)
	S.SNP_BaseQ = make(map[int][][]byte)

	S.SNP_Call_Bases = make(map[int][]byte)
	S.SNP_Call_Qual = make(map[int]float32)
}

//--------------------------------------------------------------------------------------------------
//SetPara sets values for parameters
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


////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Functions for alignment phase.
//
////////////////////////////////////////////////////////////////////////////////////////////////////


//--------------------------------------------------------------------------------------------------
// ProcessReads initializes share variables, channels, reads input reads, finds all possible SNPs,
// and updates SNP information in SNP_Prof.
// This function will be called from main program.
//--------------------------------------------------------------------------------------------------
func (S *SNP_Prof) ProcessReads() (uint32, uint32) {

	r_len := PARA_INFO.Read_len
	read_info := make([]ReadInfo, INPUT_INFO.Routine_num)
	for i := 0; i < INPUT_INFO.Routine_num; i++ {
		read_info[i].Read1, read_info[i].Read2 = make([]byte, r_len), make([]byte, r_len)
		read_info[i].Qual1, read_info[i].Qual2 = make([]byte, r_len), make([]byte, r_len)
		read_info[i].Rev_read1, read_info[i].Rev_read2 = make([]byte, r_len), make([]byte, r_len)
		read_info[i].Rev_comp_read1, read_info[i].Rev_comp_read2 = make([]byte, r_len), make([]byte, r_len)
		read_info[i].Comp_read1, read_info[i].Comp_read2 = make([]byte, r_len), make([]byte, r_len)
	}
	align_info := make([]AlignInfo, INPUT_INFO.Routine_num)
	for i := 0; i < INPUT_INFO.Routine_num; i++ {
		align_info[i].InitAlignInfo(PARA_INFO.Read_len)
	}
	match_pos := make([][]int, INPUT_INFO.Routine_num)
	for i := 0; i < INPUT_INFO.Routine_num; i++ {
		match_pos[i] = make([]int, PARA_INFO.Max_match)
	}

	read_data := make(chan ReadInfo, INPUT_INFO.Routine_num)
	read_signal := make(chan bool)
	go S.ReadReads(read_data, read_signal)

	snp_results := make(chan []SNP)
	var wg sync.WaitGroup
	for i := 0; i < INPUT_INFO.Routine_num; i++ {
		go S.FindSNPs(read_data, read_signal, snp_results, &wg, &read_info[i], &align_info[i], match_pos[i])
	}
	go func() {
		wg.Wait()
		close(snp_results)
	}()

	//Collect SNPS from results channel and update SNPs
	var snp_call_num, del_num uint32 = 0, 0
	var snp SNP
	for SNPs := range snp_results {
		snp_call_num++
		for _, snp = range SNPs {
			if len(snp.Bases) > 0 {
				S.SNP_Bases[int(snp.Pos)] = append(S.SNP_Bases[int(snp.Pos)], snp.Bases)
				S.SNP_BaseQ[int(snp.Pos)] = append(S.SNP_BaseQ[int(snp.Pos)], snp.BaseQ)
			} else {
				S.SNP_Bases[int(snp.Pos)] = append(S.SNP_Bases[int(snp.Pos)], []byte{'.'})
				S.SNP_BaseQ[int(snp.Pos)] = append(S.SNP_BaseQ[int(snp.Pos)], []byte{'.'})
				del_num++
			}
		}
	}
	return snp_call_num, del_num
}

//--------------------------------------------------------------------------------------------------
// ReadReads reads all reads from input FASTQ files and put them into data channel.
//--------------------------------------------------------------------------------------------------
func (S *SNP_Prof) ReadReads(read_data chan ReadInfo, read_signal chan bool) {

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

	read_num := 0
	scanner1 := bufio.NewScanner(f1)
	scanner2 := bufio.NewScanner(f2)
	var read_info ReadInfo
	for scanner1.Scan() && scanner2.Scan() { //ignore 1st lines in input FASTQ files
		scanner1.Scan()
		scanner2.Scan()
		read_info.Read1 = scanner1.Bytes() //use 2nd line in input FASTQ file 1
		read_info.Read2 = scanner2.Bytes() //use 2nd line in input FASTQ file 2
		scanner1.Scan() //ignore 3rd line in 1st input FASTQ file 1
		scanner2.Scan() //ignore 3rd line in 2nd input FASTQ file 2
		scanner1.Scan()
		scanner2.Scan()
		read_info.Qual1 = scanner1.Bytes() //use 4th line in input FASTQ file 1
		read_info.Qual2 = scanner2.Bytes() //use 4th line in input FASTQ file 2
		if len(read_info.Read1) > 0 && len(read_info.Read2) > 0 {
			read_num++
			read_data <- read_info
			//PrintMemStats("After putting read to data " + string(read_info.Read1))
			read_signal <- true
		}
		
		if read_num == 10000 {
			PrintMemStats("Memstats after distributing 10000 reads")
		}
	}
	close(read_data)
}

//--------------------------------------------------------------------------------------------------
// FindSNPs takes data from data channel, find all possible SNPs and put them into results channel.
//--------------------------------------------------------------------------------------------------
func (S *SNP_Prof) FindSNPs(read_data chan ReadInfo, read_signal chan bool, snp_results chan []SNP, 
	wg *sync.WaitGroup, read_info *ReadInfo, align_info *AlignInfo, match_pos []int) {

	wg.Add(1)
	defer wg.Done()
	var snps []SNP
	var read ReadInfo
	for read = range read_data {
		//PrintMemStats("Before copying all info from data chan")
		copy(read_info.Read1, read.Read1)
		copy(read_info.Read2, read.Read2)
		copy(read_info.Qual1, read.Qual1)
		copy(read_info.Qual2, read.Qual2)
		<- read_signal
		//PrintMemStats("After copying all info from data chan")
		RevComp(read_info.Read1, read_info.Rev_read1, read_info.Rev_comp_read1, read_info.Comp_read1)
		//PrintMemStats("After calculating RevComp for Read1")
		RevComp(read_info.Read2, read_info.Rev_read2, read_info.Rev_comp_read2, read_info.Comp_read2)
		//PrintMemStats("After calculating RevComp for Read2")
		snps = S.FindSNPsFromReads(read_info, align_info, match_pos)
		//PrintMemStats("After finding all SNPs from reads")
		if len(snps) > 0 {
			snp_results <- snps
		}
	}
}

//--------------------------------------------------------------------------------------------------
// FindSNPsFromReads returns SNPs found from alignment between pair-end reads and the multigenome.
// This version treats each end of the reads independently.
//--------------------------------------------------------------------------------------------------
func (S *SNP_Prof) FindSNPsFromReads(read_info *ReadInfo, align_info *AlignInfo, match_pos []int) []SNP {

	var snps, snps1, snps2 []SNP

	//Find SNPs for the first end
	//PrintMemStats("Before FindSNPsFromEnd1")
	snps1 = S.FindSNPsFromEachEnd(read_info.Read1, read_info.Rev_read1, read_info.Rev_comp_read1, 
		read_info.Comp_read1, read_info.Qual1, align_info, match_pos)
	//PrintMemStats("After FindSNPsFromEnd1")

	//Find SNPs for the second end
	//PrintMemStats("Before FindSNPsFromEnd2")
	snps2 = S.FindSNPsFromEachEnd(read_info.Read2, read_info.Rev_read2, read_info.Rev_comp_read2, 
		read_info.Comp_read2, read_info.Qual2, align_info, match_pos)
	//PrintMemStats("After FindSNPsFromEnd2")

	//Will process constrants of two ends here
	//...

	snps = append(snps, snps1...)
	snps = append(snps, snps2...)
	return snps
}

//---------------------------------------------------------------------------------------------------
// FindSNPsFromEachEnd find SNPs from alignment between read (one end) and multigenome.
//---------------------------------------------------------------------------------------------------
func (S *SNP_Prof) FindSNPsFromEachEnd(read, rev_read, rev_comp_read, comp_read, qual []byte, 
	align_info *AlignInfo, match_pos []int) []SNP {
	var has_seeds bool
	var p, s_pos, e_pos int
	var loop_num, match_num int
	var snps []SNP

	p = INPUT_INFO.Start_pos
	loop_num = 1
	for loop_num <= PARA_INFO.Iter_num {
		//fmt.Println(loop_num, "\tread2")
		//PrintMemStats("Before FindSeeds, original_read, loop_num " + strconv.Itoa(loop_num))
		s_pos, e_pos, match_num, has_seeds = INDEX.FindSeeds(read, rev_read, p, match_pos)
		//PrintMemStats("After FindSeeds, original_read, loop_num " + strconv.Itoa(loop_num))
		if has_seeds {
			//fmt.Println("read2, has seed\t", s_pos, "\t", e_pos, "\t", string(read_info.Read2))
			//PrintMemStats("Before FindSNPsFromMatch, original_read, loop_num " + strconv.Itoa(loop_num))
			snps = S.FindSNPsFromMatch(read, qual, s_pos, e_pos, match_pos, match_num, align_info)
			//PrintMemStats("After FindSeeds, original_read, loop_num " + strconv.Itoa(loop_num))
			if len(snps) > 0 {
				//fmt.Println("read2, has snp\t", s_pos, "\t", e_pos, "\t", string(read_info.Read2))
				return snps
			}
		}
		//Find SNPs for the reverse complement of the second end
		//PrintMemStats("Before FindSeeds, revcomp_read, loop_num " + strconv.Itoa(loop_num))
		s_pos, e_pos, match_num, has_seeds = INDEX.FindSeeds(rev_comp_read, comp_read, p, match_pos)
		//PrintMemStats("After FindSeeds, revcomp_read, loop_num " + strconv.Itoa(loop_num))
		if has_seeds {
			//fmt.Println("rc_read2, has seed\t", s_pos, "\t", e_pos, "\t", string(read_info.Rev_comp_read2))
			//PrintMemStats("Before FindSNPsFromMatch, revcomp_read, loop_num " + strconv.Itoa(loop_num))
			snps = S.FindSNPsFromMatch(rev_comp_read, qual, s_pos, e_pos, match_pos, match_num, align_info)
			//PrintMemStats("After FindSNPsFromMatch, revcomp_read, loop_num " + strconv.Itoa(loop_num))
			if len(snps) > 0 {
				//fmt.Println("rc_read2, has snp\t", s_pos, "\t", e_pos, "\t", string(read_info.Rev_comp_read2))
				return snps
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
	return snps
}

//---------------------------------------------------------------------------------------------------
// FindSNPsFromMatch finds SNPs from extensions of matches between read (one end) and multigenome.
//---------------------------------------------------------------------------------------------------
func (S *SNP_Prof) FindSNPsFromMatch(read, qual []byte, s_pos, e_pos int, 
	match_pos []int, match_num int, align_info *AlignInfo) []SNP {

	var pos, dis, k int
	var left_snp_pos, right_snp_pos, left_snp_idx, right_snp_idx []int
	var left_snp_val, right_snp_val [][]byte
	var snps []SNP
	var snp SNP

	for i := 0; i < match_num; i++ {
		pos = match_pos[i]
		//PrintMemStats("Before FindExtensions, match_num " + strconv.Itoa(i))
		dis, left_snp_pos, left_snp_val, left_snp_idx, right_snp_pos, right_snp_val, right_snp_idx =
			 INDEX.FindExtensions(read, s_pos, e_pos, pos, align_info)
		//PrintMemStats("After FindExtensions, match_num " + strconv.Itoa(i))
		if dis <= PARA_INFO.Dist_thres {
			if len(left_snp_pos) == 0 && len(right_snp_pos) == 0 {
				continue
			} else {
				for k = 0; k < len(left_snp_pos); k++ {
					//PrintMemStats("Before GetSNP left, snp_num " + strconv.Itoa(k))
					left_snp_qual := make([]byte, len(left_snp_val[k]))
					copy(left_snp_qual, qual[left_snp_idx[k] : left_snp_idx[k] + len(left_snp_val[k])])
					snp.Pos, snp.Bases, snp.BaseQ = uint32(left_snp_pos[k]), left_snp_val[k], left_snp_qual
					snps = append(snps, snp)
					//PrintMemStats("After GetSNP left, snp_num " + strconv.Itoa(k))
				}
				for k = 0; k < len(right_snp_pos); k++ {
					//PrintMemStats("Before GetSNP right, snp_num " + strconv.Itoa(k))
					right_snp_qual := make([]byte, len(right_snp_val[k]))
					copy(right_snp_qual, qual[right_snp_idx[k] : right_snp_idx[k] + len(right_snp_val[k])])
					snp.Pos, snp.Bases, snp.BaseQ = uint32(right_snp_pos[k]), right_snp_val[k], right_snp_qual
					snps = append(snps, snp)
					//PrintMemStats("After GetSNP right, snp_num " + strconv.Itoa(k))
				}
			}
		}
	}
	return snps
}

///////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Functions for SNP calling phase.
//
///////////////////////////////////////////////////////////////////////////////////////////////////////

/*
//-----------------------------------------------------------------------------------------------------
// CallSNPs finds SNP calls and updates SNP call information in SNP_Prof (goroutines version).
// This function will be called from main program.
//-----------------------------------------------------------------------------------------------------
func (S *SNP_Prof) CallSNPs() {

	snp_pos_chan := make(chan int, INPUT_INFO.Routine_num)
	go func() {
		for snp_pos, _ := range S.SNP_Bases {
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
		S.SNP_Call_Bases[snp_call.SNPVal.SNP_Pos] = snp_call.SNPVal.SNP_Val
		S.SNP_Call_Qual[snp_call.SNPVal.SNP_Pos] = snp_call.SNPQual
	}
}
*/

//-----------------------------------------------------------------------------------------------------
// CallSNPs finds SNP calls and updates SNP call information in SNP_Prof (no-goroutines version).
// This function will be called from main program.
//-----------------------------------------------------------------------------------------------------
func (S *SNP_Prof) CallSNPs() {
	var snp_call SNP_Call
	for snp_pos, _ := range S.SNP_Bases {
		snp_call = S.CallSNPForEachPos(snp_pos)
		S.SNP_Call_Bases[int(snp_call.Pos)] = snp_call.Bases
		S.SNP_Call_Qual[int(snp_call.Pos)] = snp_call.Qual
	}
}

//---------------------------------------------------------------------------------------------------
// CallSNPForEachPos finds SNP calls and calculates SNP call quality for each postion.
//---------------------------------------------------------------------------------------------------
func (S *SNP_Prof) CallSNPForEachPos(snp_pos int) SNP_Call {

	var snp_call SNP_Call
	snp_call.Pos = uint32(snp_pos)

	has_indels := false
	for _, snp := range S.SNP_Bases[snp_pos] {
		if len(snp) > 1 || snp[0] == '.' {
			has_indels = true
			break
		}
	}
	if has_indels {
		SNP_Qlt := make(map[string]int)
		for _, snp := range S.SNP_Bases[snp_pos] {
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
		snp_call.Bases = []byte(major_snp)
		p_max := float64(major_num)/float64(len(S.SNP_Bases[snp_pos]))
		snp_call.Qual = float32(-10*math.Log10(1 - p_max))
		return snp_call
	} else {
		var a, e []byte
		for _, a_base := range S.SNP_Bases[snp_pos] {
			a = append(a, a_base[0])
		}
		for _, e_val := range S.SNP_BaseQ[snp_pos] {
			e = append(e, e_val[0])
		}
		var b []byte
		for _, b_base := range INDEX.SNP_PROF[snp_pos] {
			b = append(b, b_base[0])
		}
		var p_b []float32
		for _, p := range S.pProb[snp_pos] {
			p_b = append(p_b, p)
		}
		for _, base := range [][]byte{[]byte{'A'}, []byte{'C'}, []byte{'G'}, []byte{'T'}} {
			if !bytes.Contains(b, base) {
				b = append(b, base[0])
				p_b = append(p_b, 0.000001)
			}
		}
		snp_call_base, snp_call_qual := CalcSNPQual(a, e, b, p_b)
		snp_call.Bases = []byte{snp_call_base}
		snp_call.Qual = snp_call_qual
		return snp_call
	}
}

//---------------------------------------------------------------------------------------------------
// CalcSNPQual determines SNP call and calculates SNP call quality based on aligned bases a,
// base qualities e, and prior prob p_b.
//---------------------------------------------------------------------------------------------------
func CalcSNPQual(a, e []byte, bases []byte, p_b []float32) (byte, float32) {
	var p_ab []float64
	for _, b := range(bases) {
		p := 1.0
		for j, w := range(a) {
			q := QualtoProb(e[j])
			if w == b {
				p = p * (1.0 - q)
			} else {
				p = p * (q/3)
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
	return bases[p_max_idx], float32(-10*math.Log10(1 - p_max))
}

//-------------------------------------------------------------------------------------------------------
// WriteSNPCalls writes SNP calls and related information to output file in tab-delimited format.
//-------------------------------------------------------------------------------------------------------
func (S *SNP_Prof) WriteSNPCalls() {

	file, err := os.Create(INPUT_INFO.SNP_call_file)
	if err != nil {
		return
	}
	defer file.Close()

	var snp_pos int
	var str_snp_pos, str_snp_val, str_snp_qual string

	SNP_Pos := make([]int, 0, len(S.SNP_Bases))
	for snp_pos, _ = range S.SNP_Call_Bases {
		SNP_Pos = append(SNP_Pos, snp_pos)
	}
	sort.Ints(SNP_Pos)

	for _, snp_pos = range SNP_Pos {
		str_snp_pos = strconv.Itoa(snp_pos)
		str_snp_val = string(S.SNP_Call_Bases[snp_pos])
		str_snp_qual = strconv.FormatFloat(float64(S.SNP_Call_Qual[snp_pos]), 'f', 5, 32)
		if str_snp_qual != "+Inf" {
			_, err = file.WriteString(str_snp_pos + "\t" + str_snp_val + "\t" + str_snp_qual + "\n")
		} else {
			_, err = file.WriteString(str_snp_pos + "\t" + str_snp_val + "\t100000\n")
		}
		if err != nil {
			fmt.Println(err)
			break
		}
	}
}
