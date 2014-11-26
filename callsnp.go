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
	SNP_Bases 	map[uint32]map[string]int
}

//--------------------------------------------------------------------------------------------------
// InitIndex initializes indexes and parameters.
// This function will be called from main program.
//--------------------------------------------------------------------------------------------------
func (S *SNP_Prof) Init(input_info InputInfo) {
	INPUT_INFO = input_info
	PARA_INFO = *SetPara(100, 0.0015, 1000)
	INDEX.Init()
	S.SNP_Calls = make(map[uint32]map[string]float64)
	S.SNP_Bases = make(map[uint32]map[string]int)
	RAND_GEN = rand.New(rand.NewSource(time.Now().UnixNano()))
}

//--------------------------------------------------------------------------------------------------
// CallSNPs initializes share variables, channels, reads input reads, finds all possible SNPs,
// and updates SNP information in SNP_Prof.
// This function will be called from main program.
//--------------------------------------------------------------------------------------------------
func (S *SNP_Prof) CallSNPs() {
	//The channel read_signal is used for signaling between goroutines which run ReadReads and FindSNPs,
	//when a FindSNPs goroutine finish copying a read to its own memory, it signals ReadReads goroutine to scan next reads.
	read_signal := make(chan bool)

	//Call a goroutine to read input reads
	read_data := make(chan *ReadInfo, INPUT_INFO.Routine_num)
	go S.ReadReads(read_data, read_signal)

	//Call goroutines to find SNPs, pass shared variable to each goroutine
	snp_results := make(chan SNP)
	var wg sync.WaitGroup
	for i := 0; i < INPUT_INFO.Routine_num; i++ {
		go S.FindSNPs(read_data, read_signal, snp_results, &wg)
	}
	//------------------------
	//For debugging
	go func() {
		GetAlignTraceInfo()
	}()
	go func() {
		GetMisAlignTraceInfo()
	}()
	//------------------------
	go func() {
		wg.Wait()
		close(snp_results)
		//------------------------
		//For debugging
		close(ALIGN_TRACE_INFO_CHAN)
		close(MIS_ALIGN_TRACE_INFO_CHAN)
		//------------------------
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
	//Output SNP calls
	fmt.Println("Outputing SNP calls...")
	S.OutputSNPCalls()

	//------------------------
	//For debugging
	fmt.Println("Processing trace info...")
	ProcessSNPTPFPInfo(S.SNP_Calls)
	ProcessSNPFNInfo(S.SNP_Calls)
	ProcessMisAlignInfo(S.SNP_Calls)
	//------------------------
}

//--------------------------------------------------------------------------------------------------
// ReadReads reads all reads from input FASTQ files and put them into data channel.
//--------------------------------------------------------------------------------------------------
func (S *SNP_Prof) ReadReads(read_data chan *ReadInfo, read_signal chan bool) {

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
	read_info := InitReadInfo(PARA_INFO.Read_len, PARA_INFO.Info_len)
	for scanner1.Scan() && scanner2.Scan() {
		read_info.Info1 = read_info.Info1[ : PARA_INFO.Info_len]
		read_info.Info2 = read_info.Info2[ : PARA_INFO.Info_len]
		copy(read_info.Info1, scanner1.Bytes()) //use 1st line in input FASTQ file 1
		copy(read_info.Info2, scanner2.Bytes()) //use 1st line in input FASTQ file 2
		read_info.Info1 = read_info.Info1[ : len(scanner1.Bytes())]
		read_info.Info2 = read_info.Info2[ : len(scanner2.Bytes())]
		scanner1.Scan()
		scanner2.Scan()
		copy(read_info.Read1, scanner1.Bytes()) //use 2nd line in input FASTQ file 1
		copy(read_info.Read2, scanner2.Bytes()) //use 2nd line in input FASTQ file 2
		scanner1.Scan() //ignore 3rd line in 1st input FASTQ file 1
		scanner2.Scan() //ignore 3rd line in 2nd input FASTQ file 2
		scanner1.Scan()
		scanner2.Scan()
		copy(read_info.Qual1, scanner1.Bytes()) //use 4th line in input FASTQ file 1
		copy(read_info.Qual2, scanner2.Bytes()) //use 4th line in input FASTQ file 2
		if len(read_info.Read1) > 0 && len(read_info.Read2) > 0 {
			read_num++
			read_data <- read_info
			read_signal <- true
		}
		if read_num%10000 == 0 {
			PrintProcessMem("Memstats after distributing 10000 reads")
		}
	}
	close(read_data)
}

//--------------------------------------------------------------------------------------------------
// FindSNPs takes data from data channel, find all possible SNPs and put them into results channel.
//--------------------------------------------------------------------------------------------------
func (S *SNP_Prof) FindSNPs(read_data chan *ReadInfo, read_signal chan bool, snp_results chan SNP, 
	wg *sync.WaitGroup) {
	wg.Add(1)
	defer wg.Done()

	//Initialize inter-function share variables
	read_info := InitReadInfo(PARA_INFO.Read_len, PARA_INFO.Info_len)
	align_info := InitAlignInfo(PARA_INFO.Read_len)
	match_pos := make([]int, PARA_INFO.Max_match)

	for read := range read_data {
		PrintMemStats("Before copying all info from data chan")
		read_info.Info1 = read_info.Info1[ : PARA_INFO.Info_len]
		read_info.Info2 = read_info.Info2[ : PARA_INFO.Info_len]
		copy(read_info.Info1, read.Info1)
		copy(read_info.Info2, read.Info2)
		read_info.Info1 = read_info.Info1[ : len(read.Info1)]
		read_info.Info2 = read_info.Info2[ : len(read.Info2)]
		copy(read_info.Read1, read.Read1)
		copy(read_info.Read2, read.Read2)
		copy(read_info.Qual1, read.Qual1)
		copy(read_info.Qual2, read.Qual2)
		<- read_signal
		
		PrintMemStats("After copying all info from data chan")
		RevComp(read_info.Read1, read_info.Qual1, read_info.Rev_read1, read_info.Rev_comp_read1, read_info.Comp_read1, read_info.Rev_qual1)
		PrintMemStats("After calculating RevComp for Read1")
		RevComp(read_info.Read2, read_info.Qual2, read_info.Rev_read2, read_info.Rev_comp_read2, read_info.Comp_read2, read_info.Rev_qual2)
		PrintMemStats("After calculating RevComp for Read2")

		S.FindSNPsFromReads(read_info, snp_results, align_info, match_pos)
		PrintMemStats("After finding all SNPs from reads")
	}
}

//--------------------------------------------------------------------------------------------------
// FindSNPsFromReads returns SNPs found from alignment between pair-end reads and the multigenome.
// This version treats each end of the reads independently.
//--------------------------------------------------------------------------------------------------
func (S *SNP_Prof) FindSNPsFromReads(read_info *ReadInfo, snp_results chan SNP, align_info *AlignInfo, match_pos []int) {
	
	var snps_arr1, snps_arr2 [][]SNP
	var match_dis1, left_align_pos1, left_align_pos2, match_dis2, right_align_pos1, right_align_pos2 []int
	var strand1, strand2 bool
	loop_num := 1
	for loop_num < 5 { //an arbitrary value, will be replaced later
		//Find SNPs for the first end
		PrintMemStats("Before FindSNPsFromEnd1")
		snps_arr1, match_dis1, left_align_pos1, right_align_pos1, strand1 = S.FindSNPsFromEachEnd(read_info.Read1, read_info.Rev_read1, 
			read_info.Rev_comp_read1, read_info.Comp_read1, read_info.Qual1, read_info.Rev_qual1, align_info, match_pos)
		PrintMemStats("After FindSNPsFromEnd1")

		//Find SNPs for the second end
		PrintMemStats("Before FindSNPsFromEnd2")
		snps_arr2, match_dis2, left_align_pos2, right_align_pos2, strand2 = S.FindSNPsFromEachEnd(read_info.Read2, read_info.Rev_read2, 
			read_info.Rev_comp_read2, read_info.Comp_read2, read_info.Qual2, read_info.Rev_qual2, align_info, match_pos)
		PrintMemStats("After FindSNPsFromEnd2")

		var idx1, idx2, pos1, pos2 int
		if len(match_dis1) > 0 && len(match_dis2) > 0 && strand1 != strand2 {
			for idx1, pos1 = range left_align_pos1 {
				for idx2, pos2 = range left_align_pos2 {
					//Check if alignments are likely pair-end alignments
					if int(math.Abs(float64(pos1 - pos2))) <= PARA_INFO.Max_diff {
						var at Align_trace_info
						at.read1 = make([]byte, len(read_info.Read1))
						at.read2 = make([]byte, len(read_info.Read2))
						copy(at.read1, read_info.Read1)
						copy(at.read2, read_info.Read2)
						at.read_info1 = make([]byte, len(read_info.Info1))
						at.read_info2 = make([]byte, len(read_info.Info2))
						copy(at.read_info1, read_info.Info1)
						copy(at.read_info2, read_info.Info2)
						at.align_pos1 = pos1
						at.align_pos2 = pos2
						at.align_right_pos1 = right_align_pos1[idx1]
						at.align_right_pos2 = right_align_pos2[idx2]
						at.align_dis1 = match_dis1[idx1]
						at.align_dis2 = match_dis2[idx2]
						
						var snp SNP
						if len(snps_arr1[idx1]) > 0 {
							for _, snp = range snps_arr1[idx1] {
								snp_results <- snp
								at.snp_pos1 = append(at.snp_pos1, snp.Pos)
								at.snp_base1 = append(at.snp_base1, snp.Bases)
								at.snp_baseq1 = append(at.snp_baseq1, snp.BaseQ)
							}
						}
						if len(snps_arr2[idx2]) > 0 {
							for _, snp = range snps_arr2[idx2] {
								snp_results <- snp
								at.snp_pos2 = append(at.snp_pos2, snp.Pos)
								at.snp_base2 = append(at.snp_base2, snp.Bases)
								at.snp_baseq2 = append(at.snp_baseq2, snp.BaseQ)
							}
						}
						ALIGN_TRACE_INFO_CHAN <- at
						return
					}
				}
			}
		} else if len(match_dis1) == 0 {
			for idx2, pos2 = range left_align_pos2 {
				var at Align_trace_info
				at.read1 = make([]byte, len(read_info.Read1))
				at.read2 = make([]byte, len(read_info.Read2))
				copy(at.read1, read_info.Read1)
				copy(at.read2, read_info.Read2)
				at.read_info1 = make([]byte, len(read_info.Info1))
				at.read_info2 = make([]byte, len(read_info.Info2))
				copy(at.read_info1, read_info.Info1)
				copy(at.read_info2, read_info.Info2)
				at.align_pos1 = -1
				at.align_pos2 = pos2
				at.align_right_pos1 = -1
				at.align_right_pos2 = right_align_pos2[idx2]
				at.align_dis1 = -1
				at.align_dis2 = match_dis2[idx2]
				
				var snp SNP
				if len(snps_arr2[idx2]) > 0 {
					for _, snp = range snps_arr2[idx2] {
						snp_results <- snp
						at.snp_pos2 = append(at.snp_pos2, snp.Pos)
						at.snp_base2 = append(at.snp_base2, snp.Bases)
						at.snp_baseq2 = append(at.snp_baseq2, snp.BaseQ)
					}
				}
				ALIGN_TRACE_INFO_CHAN <- at
				return
			}
		} else if len(match_dis2) == 0 {
			for idx1, pos1 = range left_align_pos1 {
				var at Align_trace_info
				at.read1 = make([]byte, len(read_info.Read1))
				at.read2 = make([]byte, len(read_info.Read2))
				copy(at.read1, read_info.Read1)
				copy(at.read2, read_info.Read2)
				at.read_info1 = make([]byte, len(read_info.Info1))
				at.read_info2 = make([]byte, len(read_info.Info2))
				copy(at.read_info1, read_info.Info1)
				copy(at.read_info2, read_info.Info2)
				at.align_pos1 = pos1
				at.align_pos2 = -1
				at.align_right_pos1 = right_align_pos1[idx1]
				at.align_right_pos2 = -1
				at.align_dis1 = match_dis1[idx1]
				at.align_dis2 = -1
				
				var snp SNP
				if len(snps_arr1[idx1]) > 0 {
					for _, snp = range snps_arr1[idx1] {
						snp_results <- snp
						at.snp_pos1 = append(at.snp_pos1, snp.Pos)
						at.snp_base1 = append(at.snp_base1, snp.Bases)
						at.snp_baseq1 = append(at.snp_baseq1, snp.BaseQ)
					}
				}
				ALIGN_TRACE_INFO_CHAN <- at
				return
			}
		}
		
		loop_num++
	}
	/*
	var at Align_trace_info
	at.read1 = make([]byte, len(read_info.Read1))
	at.read2 = make([]byte, len(read_info.Read2))
	copy(at.read1, read_info.Read1)
	copy(at.read2, read_info.Read2)
	at.read_info1 = make([]byte, len(read_info.Info1))
	at.read_info2 = make([]byte, len(read_info.Info2))
	copy(at.read_info1, read_info.Info1)
	copy(at.read_info2, read_info.Info2)
	at.align_pos1 = left_align_pos1
	at.align_pos2 = left_align_pos2
	at.align_right_pos1 = right_align_pos1
	at.align_right_pos2 = right_align_pos2
	at.align_dis1 = min_dis1
	at.align_dis2 = min_dis2
	MIS_ALIGN_TRACE_INFO_CHAN <- at
	 */
}

//---------------------------------------------------------------------------------------------------
// FindSNPsFromEachEnd find SNPs from alignment between read (one end) and multigenome.
//---------------------------------------------------------------------------------------------------
func (S *SNP_Prof) FindSNPsFromEachEnd(read, rev_read, rev_comp_read, comp_read, qual, rev_qual []byte, 
	align_info *AlignInfo, match_pos []int) ([][]SNP, []int, []int, []int, bool) {
	var has_seeds bool
	var s_pos, e_pos, match_num int
	var snps_arr [][]SNP
	var match_dis, left_align_pos, right_align_pos []int

	p := INPUT_INFO.Start_pos
	loop_num := 1
	for loop_num <= PARA_INFO.Iter_num {
		PrintLoopTraceInfo(loop_num, read)
		PrintMemStats("Before FindSeeds, original_read, loop_num " + strconv.Itoa(loop_num))
		s_pos, e_pos, match_num, has_seeds = INDEX.FindSeeds(read, rev_read, p, match_pos)
		PrintMemStats("After FindSeeds, original_read, loop_num " + strconv.Itoa(loop_num))
		if has_seeds {
			PrintSeedTraceInfo("ori", e_pos, s_pos, read)
			PrintMemStats("Before FindSNPsFromMatch, original_read, loop_num " + strconv.Itoa(loop_num))
			snps_arr, match_dis, left_align_pos, right_align_pos = 
				S.FindSNPsFromMatch(read, qual, s_pos, e_pos, match_pos, match_num, align_info)
			PrintMemStats("After FindSeeds, original_read, loop_num " + strconv.Itoa(loop_num))
			PrintExtendTraceInfo("ori", read[e_pos : s_pos + 1], e_pos, s_pos, match_num, match_pos)
			if len(match_dis) > 0 {
				return snps_arr, match_dis, left_align_pos, right_align_pos, true				
			}
		}

		//Find SNPs for the reverse complement of the read
		PrintMemStats("Before FindSeeds, revcomp_read, loop_num " + strconv.Itoa(loop_num))
		s_pos, e_pos, match_num, has_seeds = INDEX.FindSeeds(rev_comp_read, comp_read, p, match_pos)
		PrintMemStats("After FindSeeds, revcomp_read, loop_num " + strconv.Itoa(loop_num))
		if has_seeds {
			PrintSeedTraceInfo("rev", e_pos, s_pos, rev_comp_read)
			PrintMemStats("Before FindSNPsFromMatch, revcomp_read, loop_num " + strconv.Itoa(loop_num))
			snps_arr, match_dis, left_align_pos, right_align_pos = 
				S.FindSNPsFromMatch(rev_comp_read, rev_qual, s_pos, e_pos, match_pos, match_num, align_info)
			PrintMemStats("After FindSNPsFromMatch, revcomp_read, loop_num " + strconv.Itoa(loop_num))
			PrintExtendTraceInfo("rev", read[e_pos : s_pos + 1], e_pos, s_pos, match_num, match_pos)
			if len(match_dis) > 0 {
				return snps_arr, match_dis, left_align_pos, right_align_pos, false
			}
		}
		//Take a new position to search
		if INPUT_INFO.Search_mode == 1 {
			RAND_GEN := rand.New(rand.NewSource(time.Now().UnixNano()))
			p = RAND_GEN.Intn(len(read) - 5)
		} else if INPUT_INFO.Search_mode == 2 {
			p = p + INPUT_INFO.Search_step
		}
		loop_num++
	}
	return snps_arr, match_dis, left_align_pos, right_align_pos, true
}

//---------------------------------------------------------------------------------------------------
// FindSNPsFromMatch finds SNPs from extensions of matches between read (one end) and multigenome.
//---------------------------------------------------------------------------------------------------
func (S *SNP_Prof) FindSNPsFromMatch(read, qual []byte, s_pos, e_pos int, 
	match_pos []int, match_num int, align_info *AlignInfo) ([][]SNP, []int, []int, []int) {

	var pos, k, dis, left_most_pos, right_most_pos int
	var left_snp_pos, right_snp_pos, left_snp_idx, right_snp_idx []int
	var left_snp_val, right_snp_val [][]byte
	var snps_arr [][]SNP
	//var snps []SNP
	var snp SNP

	var match_dis, left_pos, right_pos []int
	for i := 0; i < match_num; i++ {
		pos = match_pos[i]
		PrintMemStats("Before FindExtensions, match_num " + strconv.Itoa(i))
		dis, left_snp_pos, left_snp_val, left_snp_idx, right_snp_pos, right_snp_val, right_snp_idx, left_most_pos, right_most_pos =
			 INDEX.FindExtensions(read, s_pos, e_pos, pos, align_info)
		PrintMemStats("After FindExtensions, match_num " + strconv.Itoa(i))
		if dis <= PARA_INFO.Dist_thres {
			PrintMatchTraceInfo(i, pos, dis, left_most_pos, left_snp_pos, read)
			match_dis = append(match_dis, dis)
			left_pos = append(left_pos, left_most_pos)
			right_pos = append(right_pos, right_most_pos)
			snps := make([]SNP, 0)
			for k = 0; k < len(left_snp_pos); k++ {
				PrintMemStats("Before GetSNP left, snp_num " + strconv.Itoa(k))
				left_snp_qual := make([]byte, len(left_snp_val[k]))
				copy(left_snp_qual, qual[left_snp_idx[k] : left_snp_idx[k] + len(left_snp_val[k])])
				snp.Pos, snp.Bases, snp.BaseQ = uint32(left_snp_pos[k]), left_snp_val[k], left_snp_qual
				snps = append(snps, snp)
				PrintMemStats("After GetSNP left, snp_num " + strconv.Itoa(k))
			}
			for k = 0; k < len(right_snp_pos); k++ {
				PrintMemStats("Before GetSNP right, snp_num " + strconv.Itoa(k))
				right_snp_qual := make([]byte, len(right_snp_val[k]))
				copy(right_snp_qual, qual[right_snp_idx[k] : right_snp_idx[k] + len(right_snp_val[k])])
				snp.Pos, snp.Bases, snp.BaseQ = uint32(right_snp_pos[k]), right_snp_val[k], right_snp_qual
				snps = append(snps, snp)
				PrintMemStats("After GetSNP right, snp_num " + strconv.Itoa(k))
			}
			snps_arr = append(snps_arr, snps)
		}
	}
	return snps_arr, match_dis, left_pos, right_pos
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

	if _, snp_base_exist := S.SNP_Calls[pos]; !snp_base_exist {
		S.SNP_Bases[pos] = make(map[string]int)
	}
	S.SNP_Bases[pos][a] += 1

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

	if _, snp_base_exist := S.SNP_Calls[pos]; !snp_base_exist {
		S.SNP_Bases[pos] = make(map[string]int)
	}
	S.SNP_Bases[pos][a] += 1

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
	var str_snp_pos, str_snp_qual, str_base_num string

	SNP_Pos := make([]int, 0, len(S.SNP_Calls))
	for snp_pos, _ = range S.SNP_Calls {
		SNP_Pos = append(SNP_Pos, int(snp_pos))
	}
	sort.Ints(SNP_Pos)

	var snp_call_prob, snp_prob float64
	var snp_call, snp string
	for _, pos := range SNP_Pos {
		snp_pos = uint32(pos)
		str_snp_pos = strconv.Itoa(pos + 1)
		snp_call_prob = 0
		for snp, snp_prob = range S.SNP_Calls[snp_pos] {
			if snp_call_prob < snp_prob {
				snp_call_prob = snp_prob
				snp_call = snp
			}
		}
		str_snp_qual = strconv.FormatFloat(-10 * math.Log10(1 - snp_call_prob), 'f', 5, 32)
		str_base_num = strconv.Itoa(S.SNP_Bases[snp_pos][snp_call])
		if str_snp_qual != "+Inf" {
			_, err = file.WriteString(str_snp_pos + "\t" + snp_call + "\t" + str_snp_qual + "\t" + str_base_num + "\n")
		} else {
			_, err = file.WriteString(str_snp_pos + "\t" + snp_call + "\t1000\n")
		}
		if err != nil {
			fmt.Println(err)
			break
		}
	}
}