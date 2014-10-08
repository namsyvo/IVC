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
// It serves as temporary variable during SNP calling phase.
//--------------------------------------------------------------------------------------------------
type SNP struct {
	Pos 	uint32 //SNP postion on ref
	Bases 	[]byte //bases of SNP
	BaseQ 	[]byte //quality of bases of SNP
}


//--------------------------------------------------------------------------------------------------
// SNP_Prof represents all possible SNPs and their probablilties at all positions on reference multigenome.
// This struct also has functions defined on it for calling SNPs.
// SNP_Calls stores all possible variants at each position and their probablilities of being SNP calls.
// Their initial (prior) probablities will be obtained from reference genomes and SNP profiles.
// Their posterior probabilities will be updated during alignment phase based on information from read-multigenome alignment
type SNP_Prof struct {
	SNP_Calls 	map[uint32]map[string]float64
}

//--------------------------------------------------------------------------------------------------
// InitIndex initializes indexes and parameters.
// This function will be called from main program.
//--------------------------------------------------------------------------------------------------
func (S *SNP_Prof) Init(input_info InputInfo) {

	INPUT_INFO = input_info
	PARA_INFO = SetPara(100, 0.05)
	INDEX.Init()

	// Assign all possible SNPs and their prior probabilities from SNP profile.
	S.SNP_Calls = make(map[uint32]map[string]float64)
	
	RAND_GEN = rand.New(rand.NewSource(time.Now().UnixNano()))
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
	para_info.Dist_thres = int(0.02 * rlen) + int(math.Ceil(err*rlen + k*math.Sqrt(rlen*err*(1-err))))
	//factor 0.02 above is assigned based on rate of SNP and INDEL reported in SNP profile of human genome
	//it will be estimated from input info
	para_info.Iter_num = para_info.Iter_num_factor * (para_info.Dist_thres + 1)

	para_info.Iter_num = 5

	fmt.Println("DIST_THRES: ", para_info.Dist_thres)
	fmt.Println("ITER_NUM: ", para_info.Iter_num)

	return para_info
}

//--------------------------------------------------------------------------------------------------
// CallSNPs initializes share variables, channels, reads input reads, finds all possible SNPs,
// and updates SNP information in SNP_Prof.
// This function will be called from main program.
//--------------------------------------------------------------------------------------------------
func (S *SNP_Prof) CallSNPs() (int, int) {

	//Initialize inter-function share variables
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

	//The channel read_signal is used for signaling between goroutines which run ReadReads and FindSNPs,
	//when a FindSNPs goroutine finish copying a read to its own memory, it signals ReadReads goroutine to scan next reads.
	read_signal := make(chan bool)

	//Call a goroutine to read input reads
	read_data := make(chan ReadInfo, INPUT_INFO.Routine_num)
	go S.ReadReads(read_data, read_signal)

	//Call goroutines to find SNPs, pass shared variable to each goroutine
	snp_results := make(chan SNP)
	var wg sync.WaitGroup
	for i := 0; i < INPUT_INFO.Routine_num; i++ {
		go S.FindSNPs(read_data, read_signal, snp_results, &wg, &read_info[i], &align_info[i], match_pos[i])
	}
	go func() {
		wg.Wait()
		close(snp_results)
	}()

	//Collect SNPs from results channel and update SNPs and their probabilities
	var snp SNP
	for snp = range snp_results {
		if len(snp.Bases) == 1 {
			S.UpdateSNPProb(snp)
		} else {
			S.UpdateIndelProb(snp)
		}
	}
	return 0, 0
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
		if read_num%10000 == 0 {
			PrintMemStats("Memstats after distributing 10000 reads")
		}
	}
	close(read_data)
}

//--------------------------------------------------------------------------------------------------
// FindSNPs takes data from data channel, find all possible SNPs and put them into results channel.
//--------------------------------------------------------------------------------------------------
func (S *SNP_Prof) FindSNPs(read_data chan ReadInfo, read_signal chan bool, snp_results chan SNP, 
	wg *sync.WaitGroup, read_info *ReadInfo, align_info *AlignInfo, match_pos []int) {

	wg.Add(1)
	defer wg.Done()
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
		S.FindSNPsFromReads(read_info, snp_results, align_info, match_pos)
		//PrintMemStats("After finding all SNPs from reads")
	}
}

//--------------------------------------------------------------------------------------------------
// FindSNPsFromReads returns SNPs found from alignment between pair-end reads and the multigenome.
// This version treats each end of the reads independently.
//--------------------------------------------------------------------------------------------------
func (S *SNP_Prof) FindSNPsFromReads(read_info *ReadInfo, snp_results chan SNP, align_info *AlignInfo, match_pos []int) {

	var snps1, snps2 []SNP

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
	var snp SNP
	if len(snps1) > 0 {
		for _, snp = range snps1 {
			snp_results <- snp
		}
	}
	if len(snps2) > 0 {
		for _, snp = range snps2 {
			snp_results <- snp
		}
	}
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

	var pos, k, dis int
	var left_snp_pos, right_snp_pos, left_snp_idx, right_snp_idx []int
	var left_snp_val, right_snp_val [][]byte
	var snps []SNP
	var snp, snp SNP

	min_dis := INF
	for i := 0; i < match_num; i++ {
		pos = match_pos[i]
		//PrintMemStats("Before FindExtensions, match_num " + strconv.Itoa(i))
		dis, left_snp_pos, left_snp_val, left_snp_idx, right_snp_pos, right_snp_val, right_snp_idx =
			 INDEX.FindExtensions(read, s_pos, e_pos, pos, align_info)
		//PrintMemStats("After FindExtensions, match_num " + strconv.Itoa(i))
		if dis <= PARA_INFO.Dist_thres {
			if len(left_snp_pos) != 0 || len(right_snp_pos) != 0 {
				if min_dis > dis {
					min_dis = dis
					snps = make([]SNP, 0)
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
	}
	return snps
}

//---------------------------------------------------------------------------------------------------
// UpdateSNPProb updates SNP probablilities for all possible SNPs.
// Input: a snp of type SNP.
// Output: updated S.SNP_Calls[snp.Pos] based on snp.Bases and snp.BaseQ using Bayesian method.
//---------------------------------------------------------------------------------------------------
func (S *SNP_Prof) UpdateSNPProb(snp SNP) {
	pos := snp.Pos
	a := string(snp.Bases)
	q := snp.BaseQ[0]

	var p float64
	p_ab := make(map[string]float64)
	p_a := 0.0

	if _, snp_call_exist := S.SNP_Calls[pos]; !snp_call_exist {
		S.SNP_Calls[pos] = make(map[string]float64)
		if snps, snp_prof_exist := INDEX.SNP_PROF[int(pos)]; snp_prof_exist {
			snp_prof_num := len(snps)
			for idx, snp := range snps {
				S.SNP_Calls[pos][string(snp)] = float64(INDEX.SNP_AF[int(pos)][idx]) - float64(snp_prof_num) * EPSILON
			}
		} else {
			S.SNP_Calls[pos][string(INDEX.SEQ[int(pos)])] = 1 - 3 * EPSILON
		}
		for _, b := range STD_BASES {
			if _, ok := S.SNP_Calls[pos][string(b)]; !ok {
				S.SNP_Calls[pos][string(b)] = EPSILON
			}
		}
	}

	for b, p_b := range(S.SNP_Calls[pos]) {
		if a == b {
			p = 1.0 - math.Pow(10, -(float64(q) - 33) / 10.0) //Phred-encoding factor (33) need to be estimated from input data
		} else {
			p = math.Pow(10, -(float64(q) - 33) / 10.0) / 3 //need to be refined, e.g., checked with diff cases (snp vs. indel)
		}
		p_ab[b] = p
		p_a += p_b * p_ab[b]
	}
	for b, p_b := range(S.SNP_Calls[pos]) {
		S.SNP_Calls[pos][b] = p_b * (p_ab[b] / p_a)
	}
}

//---------------------------------------------------------------------------------------------------
// UpdateIndelProb updates Indel probablilities for all possible Indels.
// Input: a snp of type SNP.
// Output: updated S.SNP_Calls[snp.Pos] based on snp.Bases and snp.BaseQ using Bayesian method.
//---------------------------------------------------------------------------------------------------
func (S *SNP_Prof) UpdateIndelProb(snp SNP) {
	pos := snp.Pos
	a := string(snp.Bases)
	q := snp.BaseQ
	if len(a) == 0 {
		a = "."
		q = []byte{'I'} //need to be changed to a proper value
	}

	var p float64
	var qi byte
	p_ab := make(map[string]float64)
	p_a := 0.0

	if _, snp_call_exist := S.SNP_Calls[pos]; !snp_call_exist {
		S.SNP_Calls[pos] = make(map[string]float64)
		if snps, snp_prof_exist := INDEX.SNP_PROF[int(pos)]; snp_prof_exist {
			snp_prof_num := len(snps)
			for idx, snp := range snps {
				S.SNP_Calls[pos][string(snp)] = float64(INDEX.SNP_AF[int(pos)][idx]) - float64(snp_prof_num) * EPSILON
			}
		} else {
			S.SNP_Calls[pos][string(INDEX.SEQ[int(pos)])] = 1 - 3 * EPSILON
		}
		for _, b := range STD_BASES {
			if _, ok := S.SNP_Calls[pos][string(b)]; !ok {
				S.SNP_Calls[pos][string(b)] = EPSILON
			}
		}
	}

	if _, ok := S.SNP_Calls[pos][a]; !ok {
		S.SNP_Calls[pos][a] = EPSILON
	}

	for b, p_b := range(S.SNP_Calls[pos]) {
		p = 1
		if a == b {
			for _, qi = range q {
				p *= (1.0 - math.Pow(10, -(float64(qi) - 33) / 10.0)) //Phred-encoding factor (33) need to be estimated from input data
			}
		} else {
			for _, qi = range q {
				p *= (math.Pow(10, -(float64(qi) - 33) / 10.0) / 3) //need to be refined, e.g., checked with diff cases (snp vs. indel)
			}
		}
		p_ab[b] = p
		p_a += p_b * p_ab[b]
	}
	for b, p_b := range(S.SNP_Calls[pos]) {
		S.SNP_Calls[pos][b] = p_b * (p_ab[b] / p_a)
	}
}

//-------------------------------------------------------------------------------------------------------
// OutputSNPCalls determines SNP calls, convert their probabilities to Phred scores, and writes them to file
// in proper format (VCF-like format in this version).
//-------------------------------------------------------------------------------------------------------
func (S *SNP_Prof) OutputSNPCalls() {

	file, err := os.Create(INPUT_INFO.SNP_call_file)
	if err != nil {
		return
	}
	defer file.Close()

	var snp_pos uint32
	var str_snp_pos, snp_qual string

	SNP_Pos := make([]int, 0, len(S.SNP_Calls))
	for snp_pos, _ = range S.SNP_Calls {
		SNP_Pos = append(SNP_Pos, int(snp_pos))
	}
	sort.Ints(SNP_Pos)

	var snp_call_prob, snp_prob float64
	var snp_call, snp string
	for _, pos := range SNP_Pos {
		snp_pos = uint32(pos)
		str_snp_pos = strconv.Itoa(pos)
		snp_call_prob = 0
		for snp, snp_prob = range S.SNP_Calls[snp_pos] {
			if snp_call_prob < snp_prob {
				snp_call_prob = snp_prob
				snp_call = snp
			}
		}
		snp_qual = strconv.FormatFloat(-10 * math.Log10(1 - snp_call_prob), 'f', 5, 32)
		if snp_qual != "+Inf" {
			_, err = file.WriteString(str_snp_pos + "\t" + snp_call + "\t" + snp_qual + "\n")
		} else {
			_, err = file.WriteString(str_snp_pos + "\t" + snp_call + "\t1000\n")
		}
		if err != nil {
			fmt.Println(err)
			break
		}
	}
}