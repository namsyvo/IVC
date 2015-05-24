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
	"strings"
	"bytes"
)

/*--------------------------------------------------------------------------------------------------
SNP represents SNP obtained during alignment phase.
It serves as temporary variable during SNP calling phase.
--------------------------------------------------------------------------------------------------*/
type SNP struct {
	Pos 	uint32 //SNP postion on ref
	Bases 	[]byte //bases of SNP
	BaseQ 	[]byte //quality of bases of SNP
	Type    int    //type of SNP (sub, ins, del...)
	CDis    int    //chromosomal distance
	CDiff   int    //diff between aligned pos and true pos
	AProb   float64 //correct alignment prob
	CProb  float64 //correct paired-alignment prob
	RInfo 	[]byte //whole read
	SPos1 	int 	//starting pos on read1 of exact match (i.e. ending position from backward search with FM-index)
	SPos2 	int 	//starting pos on read2 of exact match
	Stra1 	bool 	//strand of read1 of exact match
	Stra2 	bool 	//strand of read2 of exact match
}

/*--------------------------------------------------------------------------------------------------
SNP_Prof represents info of aligned reads and bases at the SNP call positions on reference multigenome.
This struct also has functions defined on it for calling SNPs.
--------------------------------------------------------------------------------------------------*/
type SNP_Prof struct {
	/*
	SNP_Prob stores all possible SNPs at each position and their confident probablilities.
		Their initial (prior) probablities will be obtained from reference genomes and SNP profiles.
		Their posterior probabilities will be updated during alignment phase based on incomming aligned bases
	*/
	SNP_Prob 	map[uint32]map[string]float64 	//Probability of SNP calls
	SNP_BaseQ 	map[uint32]map[string][][]byte 	//Quality sequences (in FASTQ format) of aligned bases at the SNP call position
	SNP_Type    map[uint32]map[string][]int     //Type of SNPs (currently: 0:sub and known SNPs, 1:ins, 2:del)
	SNP_RNum 	map[uint32]map[string]int 		//Numer of reads (bases) aligned to the SNP call postion
	Chr_Dis     map[uint32]map[string][]int 	//Chromosomal distance between two aligned ends
	Chr_Diff    map[uint32]map[string][]int 	//Difference betwwen aligned postions and true postions (for simulated data)
	Aln_Prob    map[uint32]map[string][]float64 //Alignment probability of aligned read at the SNP call position
	Chr_Prob    map[uint32]map[string][]float64 //Paired-end mapping probability
	Start_Pos1  map[uint32]map[string][]int 	//Start position of alignment of the first end
	Start_Pos2  map[uint32]map[string][]int 	//Start position of alignment of the second end
	Strand1     map[uint32]map[string][]bool 	//Strand indicator of the first end ("true" if identicaly with ref)
	Strand2     map[uint32]map[string][]bool 	//Strand indicator of the second end
	Read_Info 	map[uint32]map[string][][]byte 	//Info of aligned reads at the SNP call position
}

/*--------------------------------------------------------------------------------------------------
InitIndex initializes indexes and parameters.
This function will be called from main program.
--------------------------------------------------------------------------------------------------*/
func New_SNP_Caller(input_info InputInfo) *SNP_Prof {

	INPUT_INFO = input_info
	/*
	SetPara: 100 is maximum length of reads, 500 is maximum length of info line of reads,
			700 is maximum insert size of paired-end simulated reads, 0.0015 is maximum sequencing error rate
			0.01 is mutation rate (currently is estimated from dbSNP of human genome)
	*/
	PARA_INFO = *SetPara(100, 500, 700, 0.0015, 0.01, input_info.Dist_thres, input_info.Iter_num)
	INDEX = *New_Index()
	RAND_GEN = rand.New(rand.NewSource(time.Now().UnixNano()))

	S := new(SNP_Prof)
	S.SNP_Prob = make(map[uint32]map[string]float64)
	S.SNP_BaseQ = make(map[uint32]map[string][][]byte)
	S.SNP_Type = make(map[uint32]map[string][]int)
	S.SNP_RNum = make(map[uint32]map[string]int)
	S.Chr_Dis = make(map[uint32]map[string][]int)
	S.Chr_Diff = make(map[uint32]map[string][]int)
	S.Aln_Prob = make(map[uint32]map[string][]float64)
	S.Chr_Prob = make(map[uint32]map[string][]float64)
	S.Read_Info = make(map[uint32]map[string][][]byte)
	S.Start_Pos1 = make(map[uint32]map[string][]int)
	S.Start_Pos2 = make(map[uint32]map[string][]int)
	S.Strand1 = make(map[uint32]map[string][]bool)
	S.Strand2 = make(map[uint32]map[string][]bool)

	var pos uint32
	var snp []byte
	var i, snp_prof_num int
	std_base_num := len(STD_BASES)
	for snp_pos, snp_value := range INDEX.SNP_PROF {
		snp_prof_num = len(snp_value)
		pos = uint32(snp_pos)
		S.SNP_Prob[pos] = make(map[string]float64)
		for i, snp = range snp_value {
			if len(snp) == 1 {
				S.SNP_Prob[pos][string(snp)] = float64(INDEX.SNP_AF[snp_pos][i]) - 
					NEW_SNP_RATE * float64(std_base_num - snp_prof_num)/float64(snp_prof_num)
				if S.SNP_Prob[pos][string(snp)] < NEW_SNP_RATE {
					S.SNP_Prob[pos][string(snp)] = NEW_SNP_RATE
				}
			} else {
				S.SNP_Prob[pos][string(snp)] = float64(INDEX.SNP_AF[snp_pos][i]) - 
					NEW_SNP_RATE * float64(snp_prof_num)
				if S.SNP_Prob[pos][string(snp)] < NEW_SNP_RATE {
					S.SNP_Prob[pos][string(snp)] = NEW_SNP_RATE
				}
			}
			for _, b := range STD_BASES {
				if _, ok := S.SNP_Prob[pos][string(b)]; !ok {
					S.SNP_Prob[pos][string(b)] = NEW_SNP_RATE
				}
			}
		}
		S.SNP_BaseQ[pos] = make(map[string][][]byte)
		S.SNP_Type[pos]  = make(map[string][]int)
		S.SNP_RNum[pos] = make(map[string]int)
		S.Chr_Dis[pos] = make(map[string][]int)
		S.Chr_Diff[pos] = make(map[string][]int)
		S.Aln_Prob[pos] = make(map[string][]float64)
		S.Chr_Prob[pos] = make(map[string][]float64)
		S.Read_Info[pos] = make(map[string][][]byte)
		S.Start_Pos1[pos] = make(map[string][]int)
		S.Start_Pos2[pos] = make(map[string][]int)
		S.Strand1[pos] = make(map[string][]bool)
		S.Strand2[pos] = make(map[string][]bool)
	}
	return S
}

/*--------------------------------------------------------------------------------------------------
CallSNPs initializes share variables, channels, reads input reads, finds all possible SNPs,
and updates SNP information in SNP_Prof.
This function will be called from main program.
--------------------------------------------------------------------------------------------------*/
func (S *SNP_Prof) CallSNPs() {
	//The channel read_signal is used for signaling between goroutines which run ReadReads and FindSNPs,
	//when a FindSNPs goroutine finish copying a read to its own memory, 
	//it signals ReadReads goroutine to scan next reads.
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
		GetNoAlignReadInfo()
	}()
	//------------------------
	go func() {
		wg.Wait()
		close(snp_results)
		//------------------------
		//For debugging
		close(NO_ALIGN_READ_INFO_CHAN)
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
	fmt.Println("Outputing SNP calls...")
	S.OutputSNPCalls()
	//------------------------
	//For debugging
	ProcessNoAlignReadInfo(S.SNP_Prob)
	//------------------------
}

/*--------------------------------------------------------------------------------------------------
ReadReads reads all reads from input FASTQ files and put them into data channel.
--------------------------------------------------------------------------------------------------*/
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

/*--------------------------------------------------------------------------------------------------
FindSNPs takes data from data channel, find all possible SNPs and put them into results channel.
--------------------------------------------------------------------------------------------------*/
func (S *SNP_Prof) FindSNPs(read_data chan *ReadInfo, read_signal chan bool, snp_results chan SNP, 
	wg *sync.WaitGroup) {
	wg.Add(1)
	defer wg.Done()

	//Initialize inter-function share variables
	read_info := InitReadInfo(PARA_INFO.Read_len, PARA_INFO.Info_len)
	align_info := InitAlignInfo(2 * PARA_INFO.Read_len)
	m_pos := make([]int, INPUT_INFO.Max_snum)

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
		RevComp(read_info.Read1, read_info.Qual1, read_info.Rev_read1, read_info.Rev_comp_read1, 
			read_info.Comp_read1, read_info.Rev_qual1)
		PrintMemStats("After calculating RevComp for Read1")
		RevComp(read_info.Read2, read_info.Qual2, read_info.Rev_read2, read_info.Rev_comp_read2, 
			read_info.Comp_read2, read_info.Rev_qual2)
		PrintMemStats("After calculating RevComp for Read2")

		S.FindSNPsFromPairedEnds(read_info, snp_results, align_info, m_pos)
		PrintMemStats("After finding all SNPs from reads")
	}
}

/*--------------------------------------------------------------------------------------------------
FindSNPsFromPairedEndReads returns SNPs found from alignment between pair-end reads and the multigenome.
This version treats each end of the reads independently.
--------------------------------------------------------------------------------------------------*/
func (S *SNP_Prof) FindSNPsFromPairedEnds(read_info *ReadInfo, snp_results chan SNP, align_info *AlignInfo, m_pos []int) {

	var snp SNP
	var snps1, snps2 []SNP
	var snps_get1, snps_get2 []SNP
	var l_align_pos1, l_align_pos2 int

	var s_pos_r1, e_pos_r1, s_pos_r2, e_pos_r2, m_pos_r1, m_pos_r2 []int
	var strand_r1, strand_r2 []bool

	//---------------------------------------------------------------------
	//get info for simulated reads, with specific format of testing dataset
	//need to be re-implemented for general data
    read_info1_tokens := bytes.Split(read_info.Info1, []byte{'_'})
	var true_pos1, true_pos2 int64
	var read_info1, read_info2 []byte
    if read_info1_tokens[0][1] != 'r' {
        true_pos1, _ = strconv.ParseInt(string(read_info1_tokens[2]), 10, 64)
		true_pos2, _ = strconv.ParseInt(string(read_info1_tokens[3]), 10, 64)
		read_info1, read_info2 = make([]byte, len(read_info.Info1)), make([]byte, len(read_info.Info2))
		copy(read_info1, read_info.Info1)
		copy(read_info2, read_info.Info2)
	} else {
        true_pos1, _ = strconv.ParseInt(string(read_info1_tokens[1]), 10, 64)
		true_pos2, _ = strconv.ParseInt(string(read_info1_tokens[2]), 10, 64)
		read_info1, read_info2 = make([]byte, len(read_info.Info1)), make([]byte, len(read_info.Info2))
		copy(read_info1, read_info.Info1)
		copy(read_info2, read_info.Info2)
	}
	//---------------------------------------------------------------------

	var p_idx, s_idx int
	var paired_prob, align_prob1, align_prob2 float64

	//Try to align both ends
	loop_num := 1
	paired_prob = math.MaxFloat64
	var has_seeds bool
	for loop_num <= PARA_INFO.Iter_num { //temp value, will be replaced later
		PrintLoopTraceInfo(loop_num, "FindSNPsFromReads")
		s_pos_r1, e_pos_r1, s_pos_r2, e_pos_r2, m_pos_r1, m_pos_r2, strand_r1, strand_r2, has_seeds = S.FindSeedsFromPairedEnds(read_info)
		if has_seeds {
			for p_idx = 0; p_idx < len(s_pos_r1); p_idx++ {
				//For conventional paired-end sequencing (i.e. Illumina) the directions should be F-R
				//For other kinds of variants (e.g inversions) or other technologies, they can be F-F or R-R
				//For mate-pair, they can be R-F (need to be confirmed)
				if strand_r1[p_idx] == strand_r2[p_idx] {
					continue
				}
				//Find SNPs for the first end
				PrintMemStats("Before FindSNPsFromEnd1")
				if strand_r1[p_idx] == true {
					snps1, l_align_pos1, _, align_prob1 = S.FindSNPsFromExtension(s_pos_r1[p_idx], e_pos_r1[p_idx], 
						m_pos_r1[p_idx], read_info.Read1, read_info.Qual1, align_info)
				} else {
					snps1, l_align_pos1, _, align_prob1 = S.FindSNPsFromExtension(s_pos_r1[p_idx], e_pos_r1[p_idx], 
						m_pos_r1[p_idx], read_info.Rev_comp_read1, read_info.Rev_qual1, align_info)
				}
				PrintMemStats("After FindSNPsFromEnd1")
				
				//Find SNPs for the second end
				PrintMemStats("Before FindSNPsFromEnd2")
				if strand_r2[p_idx] == true {
					snps2, l_align_pos2, _, align_prob2 = S.FindSNPsFromExtension(s_pos_r2[p_idx], e_pos_r2[p_idx], 
						m_pos_r2[p_idx], read_info.Read2, read_info.Qual2, align_info)
				} else {
					snps2, l_align_pos2, _, align_prob2 = S.FindSNPsFromExtension(s_pos_r2[p_idx], e_pos_r2[p_idx], 
						m_pos_r2[p_idx], read_info.Rev_comp_read2, read_info.Rev_qual2, align_info)
				}
				PrintMemStats("After FindSNPsFromEnd2")
				
				if align_prob1 != -1 && align_prob2 != -1 {
					a_prob := -math.Log10(math.Exp(-math.Pow(math.Abs(float64(l_align_pos1 - l_align_pos2)) - 400.0, 2.0) / (2*50*50)))
					if paired_prob > align_prob1 + align_prob2 {
						paired_prob = align_prob1 + align_prob2
						snps_get1 = make([]SNP, len(snps1))
						PrintGetSNP(paired_prob, align_prob1, align_prob2, snps1, snps2)
						if len(snps1) > 0 {
							for s_idx = 0; s_idx < len(snps1); s_idx++ {
								snps_get1[s_idx].Pos = snps1[s_idx].Pos
								snps_get1[s_idx].Bases = make([]byte, len(snps1[s_idx].Bases))
								snps_get1[s_idx].BaseQ = make([]byte, len(snps1[s_idx].BaseQ))
								copy(snps_get1[s_idx].Bases, snps1[s_idx].Bases)
								copy(snps_get1[s_idx].BaseQ, snps1[s_idx].BaseQ)
								snps_get1[s_idx].Type = snps1[s_idx].Type
								snps_get1[s_idx].CDis = l_align_pos1 - l_align_pos2
								snps_get1[s_idx].CDiff = l_align_pos1 - int(true_pos1)
								snps_get1[s_idx].AProb = align_prob1
								snps_get1[s_idx].CProb = a_prob
								snps_get1[s_idx].RInfo = read_info1
								snps_get1[s_idx].SPos1 = e_pos_r1[p_idx]
								snps_get1[s_idx].SPos2 = e_pos_r2[p_idx]
								snps_get1[s_idx].Stra1 = strand_r1[p_idx]
								snps_get1[s_idx].Stra2 = strand_r2[p_idx]
							}
						}
						snps_get2 = make([]SNP, len(snps2))
						if len(snps2) > 0 {
							for s_idx = 0; s_idx < len(snps2); s_idx++ {
								snps_get2[s_idx].Pos = snps2[s_idx].Pos
								snps_get2[s_idx].Bases = make([]byte, len(snps2[s_idx].Bases))
								snps_get2[s_idx].BaseQ = make([]byte, len(snps2[s_idx].BaseQ))
								copy(snps_get2[s_idx].Bases, snps2[s_idx].Bases)
								copy(snps_get2[s_idx].BaseQ, snps2[s_idx].BaseQ)
								snps_get2[s_idx].Type = snps2[s_idx].Type
								snps_get2[s_idx].CDis = l_align_pos1 - l_align_pos2
								snps_get2[s_idx].CDiff = l_align_pos2 - int(true_pos2)
								snps_get2[s_idx].AProb = align_prob2
								snps_get2[s_idx].CProb = a_prob
								snps_get2[s_idx].RInfo = read_info2
								snps_get2[s_idx].SPos1 = e_pos_r1[p_idx]
								snps_get2[s_idx].SPos2 = e_pos_r2[p_idx]
								snps_get2[s_idx].Stra1 = strand_r1[p_idx]
								snps_get2[s_idx].Stra2 = strand_r2[p_idx]
							}
						}
					}
				}
			}
		}
		if paired_prob < 1 {
			break
		}
		loop_num++
	}
	if paired_prob <= 2 * PARA_INFO.Prob_thres {
		if len(snps_get1) > 0 {
			for _, snp = range snps_get1 {
				snp_results <- snp
			}
		}
		if len(snps_get2) > 0 {
			for _, snp = range snps_get2 {
				snp_results <- snp
			}
		}
		return
	}

	//Cannot align any ends, consider as unaligned reads
	var at Align_trace_info
	at.read_info1 = read_info1
	at.read_info2 = read_info2
	NO_ALIGN_READ_INFO_CHAN <- at
}

/*--------------------------------------------------------------------------------------------------
FindSeedsFromPairedEnds find all pairs of seeds which have proper chromosome distances.
--------------------------------------------------------------------------------------------------*/
func (S *SNP_Prof) FindSeedsFromPairedEnds(read_info *ReadInfo) ([]int, []int, []int, []int, []int, 
	[]int, []bool, []bool, bool) {

	var has_seeds_r1_or, has_seeds_r1_rc, has_seeds_r2_or, has_seeds_r2_rc bool
	var s_pos_r1_or, e_pos_r1_or, m_num_r1_or, s_pos_r1_rc, e_pos_r1_rc, m_num_r1_rc int
	var s_pos_r2_or, e_pos_r2_or, m_num_r2_or, s_pos_r2_rc, e_pos_r2_rc, m_num_r2_rc int
	var s_pos_r1, e_pos_r1, s_pos_r2, e_pos_r2, m_pos_r1, m_pos_r2 []int
	var strand_r1, strand_r2 []bool
	var i, j int

	//Take an initial position to search
	r_pos_r1_or := INPUT_INFO.Start_pos
	r_pos_r1_rc := INPUT_INFO.Start_pos
	r_pos_r2_or := INPUT_INFO.Start_pos
	r_pos_r2_rc := INPUT_INFO.Start_pos
	if INPUT_INFO.Search_mode == 1 {
		RAND_GEN := rand.New(rand.NewSource(time.Now().UnixNano()))
		r_pos_r1_or = RAND_GEN.Intn(len(read_info.Read1) - 5)
		r_pos_r1_rc = RAND_GEN.Intn(len(read_info.Read1) - 5)
		r_pos_r2_or = RAND_GEN.Intn(len(read_info.Read2) - 5)
		r_pos_r2_rc = RAND_GEN.Intn(len(read_info.Read2) - 5)
	}
	m_pos_r1_or := make([]int, INPUT_INFO.Max_snum)
	m_pos_r1_rc := make([]int, INPUT_INFO.Max_snum)
	m_pos_r2_or := make([]int, INPUT_INFO.Max_snum)
	m_pos_r2_rc := make([]int, INPUT_INFO.Max_snum)

	loop_num := 1
	for loop_num <= PARA_INFO.Iter_num { //temp value, will be replaced later
		PrintLoopTraceInfo(loop_num, "FindSeedsFromPairedEnds, First:\t" + string(read_info.Read1))
		PrintLoopTraceInfo(loop_num, "FindSeedsFromPairedEnds, Second:\t" + string(read_info.Read2))

		PrintMemStats("Before FindSeeds, loop_num " + strconv.Itoa(loop_num))
		s_pos_r1_or, e_pos_r1_or, m_num_r1_or, has_seeds_r1_or = 
			INDEX.FindSeeds(read_info.Read1, read_info.Rev_read1, r_pos_r1_or, m_pos_r1_or)
		PrintSeedTraceInfo("r1_or", e_pos_r1_or, s_pos_r1_or, read_info.Read1)
		if has_seeds_r1_or {
			PrintExtendTraceInfo("r1_or", read_info.Read1[e_pos_r1_or : s_pos_r1_or + 1], 
				e_pos_r1_or, s_pos_r1_or, m_num_r1_or, m_pos_r1_or)
		}
		s_pos_r1_rc, e_pos_r1_rc, m_num_r1_rc, has_seeds_r1_rc = 
			INDEX.FindSeeds(read_info.Rev_comp_read1, read_info.Comp_read1, r_pos_r1_rc, m_pos_r1_rc)
		PrintSeedTraceInfo("r1_rc", e_pos_r1_rc, s_pos_r1_rc, read_info.Rev_comp_read1)
		if has_seeds_r1_rc {
			PrintExtendTraceInfo("r1_rc", read_info.Rev_comp_read1[e_pos_r1_rc : s_pos_r1_rc + 1], 
				e_pos_r1_rc, s_pos_r1_rc, m_num_r1_rc, m_pos_r1_rc)
		}
		s_pos_r2_or, e_pos_r2_or, m_num_r2_or, has_seeds_r2_or = 
			INDEX.FindSeeds(read_info.Read2, read_info.Rev_read2, r_pos_r2_or, m_pos_r2_or)
		PrintSeedTraceInfo("r2_or", e_pos_r2_or, s_pos_r2_or, read_info.Read2)
		if has_seeds_r2_or {
			PrintExtendTraceInfo("r2_or", read_info.Read1[e_pos_r2_or : s_pos_r2_or + 1], 
				e_pos_r2_or, s_pos_r2_or, m_num_r2_or, m_pos_r2_or)
		}
		s_pos_r2_rc, e_pos_r2_rc, m_num_r2_rc, has_seeds_r2_rc = 
			INDEX.FindSeeds(read_info.Rev_comp_read2, read_info.Comp_read2, r_pos_r2_rc, m_pos_r2_rc)
		PrintSeedTraceInfo("r2_rc", e_pos_r2_rc, s_pos_r2_rc, read_info.Rev_comp_read2)
		if has_seeds_r2_rc {
			PrintExtendTraceInfo("r2_rc", read_info.Rev_comp_read2[e_pos_r2_rc : s_pos_r2_rc + 1], 
				e_pos_r2_rc, s_pos_r2_rc, m_num_r2_rc, m_pos_r2_rc)
		}
		PrintMemStats("After FindSeeds, loop_num " + strconv.Itoa(loop_num))

		if has_seeds_r1_or && has_seeds_r2_rc {
			PrintExtendTraceInfo("r1_or(F1R2)", read_info.Read1[e_pos_r1_or : s_pos_r1_or + 1], 
				e_pos_r1_or, s_pos_r1_or, m_num_r1_or, m_pos_r1_or)
			PrintExtendTraceInfo("r2_rc(F1R2)", read_info.Read2[e_pos_r2_rc : s_pos_r2_rc + 1], 
				e_pos_r2_rc, s_pos_r2_rc, m_num_r2_rc, m_pos_r2_rc)
			for i = 0; i < m_num_r1_or; i++ {
				for j = 0; j < m_num_r2_rc; j++ {
					//Check if alignments are likely pair-end alignments
					if (m_pos_r2_rc[j] - m_pos_r1_or[i]) >= PARA_INFO.Read_len && 
						(m_pos_r2_rc[j] - m_pos_r1_or[i]) <= PARA_INFO.Read_len + PARA_INFO.Max_ins {

						PrintPairedSeedInfo("r1_or, r2_rc, paired pos", m_pos_r1_or[i], m_pos_r2_rc[j])
						s_pos_r1 = append(s_pos_r1, s_pos_r1_or)
						e_pos_r1 = append(e_pos_r1, e_pos_r1_or)
						s_pos_r2 = append(s_pos_r2, s_pos_r2_rc)
						e_pos_r2 = append(e_pos_r2, e_pos_r2_rc)
						m_pos_r1 = append(m_pos_r1, m_pos_r1_or[i])
						m_pos_r2 = append(m_pos_r2, m_pos_r2_rc[j])
						strand_r1 = append(strand_r1, true)
						strand_r2 = append(strand_r2, false)
					}
				}
			}
		}
		if has_seeds_r1_rc && has_seeds_r2_or {
			PrintExtendTraceInfo("r1_rc (F2R1)", read_info.Read1[e_pos_r1_rc : s_pos_r1_rc + 1], 
				e_pos_r1_rc, s_pos_r1_rc, m_num_r1_rc, m_pos_r1_rc)
			PrintExtendTraceInfo("r2_or (F2R1)", read_info.Read2[e_pos_r2_or : s_pos_r2_or + 1], 
				e_pos_r2_or, s_pos_r2_or, m_num_r2_or, m_pos_r2_or)
			for i = 0; i < m_num_r1_rc; i++ {
				for j = 0; j < m_num_r2_or; j++ {
					//Check if alignments are likely pair-end alignments
					if (m_pos_r1_rc[i] - m_pos_r2_or[j]) >= PARA_INFO.Read_len && 
						(m_pos_r1_rc[i] - m_pos_r2_or[j]) <= PARA_INFO.Read_len + PARA_INFO.Max_ins {

						PrintPairedSeedInfo("r1_rc, r2_or, paired pos", m_pos_r1_rc[i], m_pos_r2_or[j])
						s_pos_r1 = append(s_pos_r1, s_pos_r1_rc)
						e_pos_r1 = append(e_pos_r1, e_pos_r1_rc)
						s_pos_r2 = append(s_pos_r2, s_pos_r2_or)
						e_pos_r2 = append(e_pos_r2, e_pos_r2_or)
						m_pos_r1 = append(m_pos_r1, m_pos_r1_rc[i])
						m_pos_r2 = append(m_pos_r2, m_pos_r2_or[j])
						strand_r1 = append(strand_r1, false)
						strand_r2 = append(strand_r2, true)
					}
				}
			}
		}
		if len(s_pos_r1) >= 1 && len(s_pos_r1) <= INPUT_INFO.Max_psnum {
			return s_pos_r1, e_pos_r1, s_pos_r2, e_pos_r2, m_pos_r1, m_pos_r2, strand_r1, strand_r2, true
		}
		//Take a new position to search
		r_pos_r1_or = r_pos_r1_or + INPUT_INFO.Search_step
		r_pos_r1_rc = r_pos_r1_rc + INPUT_INFO.Search_step
		r_pos_r2_or = r_pos_r2_or + INPUT_INFO.Search_step
		r_pos_r2_rc = r_pos_r2_rc + INPUT_INFO.Search_step
		if INPUT_INFO.Search_mode == 1 {
			RAND_GEN := rand.New(rand.NewSource(time.Now().UnixNano()))
			r_pos_r1_or = RAND_GEN.Intn(len(read_info.Read1) - 5)
			r_pos_r1_rc = RAND_GEN.Intn(len(read_info.Read1) - 5)
			r_pos_r2_or = RAND_GEN.Intn(len(read_info.Read2) - 5)
			r_pos_r2_rc = RAND_GEN.Intn(len(read_info.Read2) - 5)
		}
		loop_num++
	}
	return s_pos_r1, e_pos_r1, s_pos_r2, e_pos_r2, m_pos_r1, m_pos_r2, strand_r1, strand_r2, false
}

/*--------------------------------------------------------------------------------------------------
FindSNPsFromExtension determines SNPs based on alignment between reads and multi-genomes.
	Extend read and ref from exact matches found from bachward search with FM-index
	Perform backward (for left extension of read and ref) and forward alignment (for right extension)
	between read and multigenome to determine aligned bases as candidates for SNP calls
--------------------------------------------------------------------------------------------------*/
func (S *SNP_Prof) FindSNPsFromExtension(s_pos, e_pos, m_pos int, read, qual []byte, 
	align_info *AlignInfo) ([]SNP, int, int, float64) {

	var ref_l_flank, ref_r_flank, read_l_flank, read_r_flank, qual_l_flank, qual_r_flank []byte
	var isSNP, isSameLenSNP bool

	PrintMemStats("Before FindSNPsFromExtension, m_pos " + strconv.Itoa(m_pos))
	l_ext_add_len := 0
	l_most_pos := m_pos - e_pos - l_ext_add_len
	for i := m_pos - e_pos; i < m_pos; i++ {
		_, isSNP = INDEX.SNP_PROF[i]
		_, isSameLenSNP = INDEX.SAME_LEN_SNP[i]
		if isSNP && !isSameLenSNP {
			l_ext_add_len++
		}
	}
	if l_most_pos >= 0 {
		ref_l_flank = INDEX.SEQ[l_most_pos : m_pos + PARA_INFO.Seed_backup]
	} else {
		ref_l_flank = INDEX.SEQ[0 : m_pos + PARA_INFO.Seed_backup]
	}
	read_l_flank, qual_l_flank = read[ : e_pos + PARA_INFO.Seed_backup], qual[ : e_pos + PARA_INFO.Seed_backup]
	left_d, left_D, l_bt_mat, l_m, l_n, l_snp_pos, l_snp_base, l_snp_qual, l_snp_type :=
		S.BackwardDistance(read_l_flank, qual_l_flank, ref_l_flank, l_most_pos, align_info.Bw_Dist_D, 
			align_info.Bw_Dist_IS, align_info.Bw_Dist_IT, align_info.Bw_Trace_D, align_info.Bw_Trace_IS, align_info.Bw_Trace_IT)

	right_ext_add_len := 0
	r_most_pos := m_pos + s_pos - e_pos + 1 - PARA_INFO.Seed_backup
	for i := r_most_pos; i < r_most_pos + (len(read) - s_pos) - 1; i++ {
		_, isSNP = INDEX.SNP_PROF[i]
		_, isSameLenSNP = INDEX.SAME_LEN_SNP[i]
		if isSNP && !isSameLenSNP {
			right_ext_add_len++
		}
	}
	if r_most_pos + PARA_INFO.Seed_backup + (len(read) - s_pos) - 1 + right_ext_add_len <= len(INDEX.SEQ) {
		ref_r_flank = INDEX.SEQ[r_most_pos : r_most_pos + PARA_INFO.Seed_backup + (len(read) - s_pos) - 1 + right_ext_add_len]
	} else {
		ref_r_flank = INDEX.SEQ[r_most_pos : len(INDEX.SEQ)]
	}
	read_r_flank, qual_r_flank = read[s_pos + 1 - PARA_INFO.Seed_backup : ], qual[s_pos + 1 - PARA_INFO.Seed_backup : ]
	right_d, right_D, r_bt_mat, r_m, r_n, r_snp_pos, r_snp_base, r_snp_qual, r_snp_type :=
		S.ForwardDistance(read_r_flank, qual_r_flank, ref_r_flank, r_most_pos, align_info.Fw_Dist_D, 
			align_info.Fw_Dist_IS, align_info.Fw_Dist_IT, align_info.Fw_Trace_D, align_info.Fw_Trace_IS, align_info.Fw_Trace_IT)

	var snps_arr []SNP
	prob := left_d + right_d + left_D + right_D
	if prob <= PARA_INFO.Prob_thres {
		if l_m > 0 && l_n > 0 {
			l_pos, l_base, l_qual, l_type := S.BackwardTraceBack(read_l_flank, qual_l_flank, ref_l_flank, l_m, l_n, l_most_pos, 
				l_bt_mat, align_info.Bw_Trace_D, align_info.Bw_Trace_IS, align_info.Bw_Trace_IT)
			l_snp_pos = append(l_snp_pos, l_pos...)
			l_snp_base = append(l_snp_base, l_base...)
			l_snp_qual = append(l_snp_qual, l_qual...)
			l_snp_type = append(l_snp_type, l_type...)
		}
		PrintMatchTraceInfo(m_pos, l_most_pos, prob, l_snp_pos, read)
		if r_m > 0 && r_n > 0 {
			r_pos, r_base, r_qual, r_type := S.ForwardTraceBack(read_r_flank, qual_r_flank, ref_r_flank, r_m, r_n, r_most_pos, 
				r_bt_mat, align_info.Fw_Trace_D, align_info.Fw_Trace_IS, align_info.Fw_Trace_IT)
			r_snp_pos = append(r_snp_pos, r_pos...)
			r_snp_base = append(r_snp_base, r_base...)
			r_snp_qual = append(r_snp_qual, r_qual...)
			r_snp_type = append(r_snp_type, r_type...)
		}
		PrintMatchTraceInfo(m_pos, r_most_pos, prob, r_snp_pos, read)
		var k int
		for k = 0; k < len(l_snp_pos); k++ {
			PrintMemStats("Before GetSNP left, snp_num " + strconv.Itoa(k))
			var snp SNP
			snp.Pos, snp.Bases, snp.BaseQ, snp.Type = uint32(l_snp_pos[k]), l_snp_base[k], l_snp_qual[k], l_snp_type[k]
			snps_arr = append(snps_arr, snp)
			PrintMemStats("After GetSNP left, snp_num " + strconv.Itoa(k))
		}
		for k = 0; k < len(r_snp_pos); k++ {
			PrintMemStats("Before GetSNP right, snp_num " + strconv.Itoa(k))
			var snp SNP
			snp.Pos, snp.Bases, snp.BaseQ, snp.Type = uint32(r_snp_pos[k]), r_snp_base[k], r_snp_qual[k], r_snp_type[k]
			snps_arr = append(snps_arr, snp)
			PrintMemStats("After GetSNP right, snp_num " + strconv.Itoa(k))
		}
		PrintMemStats("After FindSNPsFromExtension, m_pos " + strconv.Itoa(m_pos))
		return snps_arr, l_most_pos, r_most_pos, prob
	}
	PrintMemStats("After FindSNPsFromExtension, m_pos " + strconv.Itoa(m_pos))
	return snps_arr, -1, -1, -1
}

/*--------------------------------------------------------------------------------------------------
UpdateSNPProb updates SNP probablilities for all possible SNPs.
	Input: a snp of type SNP.
	Output: updated S.SNP_Prob[snp.Pos] based on snp.Bases and snp.BaseQ using Bayesian method.
--------------------------------------------------------------------------------------------------*/
func (S *SNP_Prof) UpdateSNPProb(snp SNP) {
	pos := snp.Pos
	a := string(snp.Bases[0])
	q := snp.BaseQ[0]

	if _, snp_exist := S.SNP_Prob[pos]; !snp_exist {
		S.SNP_Prob[pos] = make(map[string]float64)
		S.SNP_Prob[pos][string(INDEX.SEQ[int(pos)])] = 1 - 3 * NEW_SNP_RATE
		for _, b := range STD_BASES {
			if _, ok := S.SNP_Prob[pos][string(b)]; !ok {
				S.SNP_Prob[pos][string(b)] = NEW_SNP_RATE
			}
		}
		S.SNP_BaseQ[pos] = make(map[string][][]byte)
		S.SNP_Type[pos] = make(map[string][]int)
		S.SNP_RNum[pos] = make(map[string]int)
		S.Chr_Dis[pos] = make(map[string][]int)
		S.Chr_Diff[pos] = make(map[string][]int)
		S.Aln_Prob[pos] = make(map[string][]float64)
		S.Chr_Prob[pos] = make(map[string][]float64)
		S.Read_Info[pos] = make(map[string][][]byte)
		S.Start_Pos1[pos] = make(map[string][]int)
		S.Start_Pos2[pos] = make(map[string][]int)
		S.Strand1[pos] = make(map[string][]bool)
		S.Strand2[pos] = make(map[string][]bool)
	}
	S.SNP_BaseQ[pos][a] = append(S.SNP_BaseQ[pos][a], snp.BaseQ)
	S.SNP_Type[pos][a] = append(S.SNP_Type[pos][a], snp.Type)
	S.SNP_RNum[pos][a] += 1
	S.Chr_Dis[pos][a] = append(S.Chr_Dis[pos][a], snp.CDis)
	S.Chr_Diff[pos][a] = append(S.Chr_Diff[pos][a], snp.CDiff)
	S.Aln_Prob[pos][a] = append(S.Aln_Prob[pos][a], snp.AProb)
	S.Chr_Prob[pos][a] = append(S.Chr_Prob[pos][a], snp.CProb)
	S.Read_Info[pos][a] = append(S.Read_Info[pos][a], snp.RInfo)
	S.Start_Pos1[pos][a] = append(S.Start_Pos1[pos][a], snp.SPos1)
	S.Start_Pos2[pos][a] = append(S.Start_Pos2[pos][a], snp.SPos2)
	S.Strand1[pos][a] = append(S.Strand1[pos][a], snp.Stra1)
	S.Strand2[pos][a] = append(S.Strand2[pos][a], snp.Stra2)

	var p float64
	p_ab := make(map[string]float64)
	p_a := 0.0

	for b, p_b := range(S.SNP_Prob[pos]) {
		if a == b {
			p = 1.0 - math.Pow(10, -(float64(q) - 33) / 10.0) //Phred-encoding factor (33) need to be estimated from input data
		} else {
			p = math.Pow(10, -(float64(q) - 33) / 10.0) / 3 //need to be refined, e.g., checked with diff cases (snp vs. indel)
		}
		p_ab[b] = p
		p_a += p_b * p_ab[b]
	}
	for b, p_b := range(S.SNP_Prob[pos]) {
		S.SNP_Prob[pos][b] = p_b * (p_ab[b] / p_a)
	}
}

/*--------------------------------------------------------------------------------------------------
UpdateIndelProb updates Indel probablilities for all possible Indels.
	Input: a snp of type SNP.
	Output: updated S.SNP_Prob[snp.Pos] based on snp.Bases and snp.BaseQ using Bayesian method.
--------------------------------------------------------------------------------------------------*/
// Notice: Need to be corrected!
func (S *SNP_Prof) UpdateIndelProb(snp SNP) {
	pos := snp.Pos
	a := string(snp.Bases)
	q := snp.BaseQ

	if _, snp_exist := S.SNP_Prob[pos]; !snp_exist {
		S.SNP_Prob[pos] = make(map[string]float64)
		S.SNP_Prob[pos][string(INDEX.SEQ[int(pos)])] = 1 - 3 * NEW_SNP_RATE
		for _, b := range STD_BASES {
			if _, ok := S.SNP_Prob[pos][string(b)]; !ok {
				S.SNP_Prob[pos][string(b)] = NEW_SNP_RATE
			}
		}
		S.SNP_BaseQ[pos] = make(map[string][][]byte)
		S.SNP_Type[pos] = make(map[string][]int)
		S.SNP_RNum[pos] = make(map[string]int)
		S.Chr_Dis[pos] = make(map[string][]int)
		S.Chr_Diff[pos] = make(map[string][]int)
		S.Aln_Prob[pos] = make(map[string][]float64)
		S.Chr_Prob[pos] = make(map[string][]float64)
		S.Read_Info[pos] = make(map[string][][]byte)
		S.Start_Pos1[pos] = make(map[string][]int)
		S.Start_Pos2[pos] = make(map[string][]int)
		S.Strand1[pos] = make(map[string][]bool)
		S.Strand2[pos] = make(map[string][]bool)
	}
	//Notice: Using NEW_SNP_RATE, need to consider NEW_INDEL_RATE instead
	if _, ok := S.SNP_Prob[pos][a]; !ok {
		S.SNP_Prob[pos][a] = NEW_SNP_RATE
		S.SNP_Prob[pos][string(INDEX.SEQ[int(pos)])] -= NEW_SNP_RATE
	}
	S.SNP_BaseQ[pos][a] = append(S.SNP_BaseQ[pos][a], snp.BaseQ)
	S.SNP_Type[pos][a] = append(S.SNP_Type[pos][a], snp.Type)
	S.SNP_RNum[pos][a] += 1
	S.Chr_Dis[pos][a] = append(S.Chr_Dis[pos][a], snp.CDis)
	S.Chr_Diff[pos][a] = append(S.Chr_Diff[pos][a], snp.CDiff)
	S.Aln_Prob[pos][a] = append(S.Aln_Prob[pos][a], snp.AProb)
	S.Chr_Prob[pos][a] = append(S.Chr_Prob[pos][a], snp.CProb)
	S.Read_Info[pos][a] = append(S.Read_Info[pos][a], snp.RInfo)
	S.Start_Pos1[pos][a] = append(S.Start_Pos1[pos][a], snp.SPos1)
	S.Start_Pos2[pos][a] = append(S.Start_Pos2[pos][a], snp.SPos2)
	S.Strand1[pos][a] = append(S.Strand1[pos][a], snp.Stra1)
	S.Strand2[pos][a] = append(S.Strand2[pos][a], snp.Stra2)

	var p float64
	var qi byte
	p_ab := make(map[string]float64)
	p_a := 0.0

	for b, p_b := range(S.SNP_Prob[pos]) {
		p = 1
		if a == b {
			for _, qi = range q[0:1] { //temporary change
				p *= (1.0 - math.Pow(10, -(float64(qi) - 33) / 10.0)) //Phred-encoding factor (33) need to be estimated from input data
			}
		} else {
			for _, qi = range q[0:1] { //temporary change
				p *= (math.Pow(10, -(float64(qi) - 33) / 10.0) / 3) //need to be refined, e.g., checked with diff cases (snp vs. indel)
			}
		}
		p_ab[b] = p
		p_a += p_b * p_ab[b]
	}
	for b, p_b := range(S.SNP_Prob[pos]) {
		S.SNP_Prob[pos][b] = p_b * (p_ab[b] / p_a)
	}
}

/*--------------------------------------------------------------------------------------------------
OutputSNPCalls determines SNP calls, convert their probabilities to Phred scores, and writes them
to file in proper format (VCF-like format at this stage).
--------------------------------------------------------------------------------------------------*/
func (S *SNP_Prof) OutputSNPCalls() {

	file, err := os.Create(INPUT_INFO.SNP_call_file)
	if err != nil {
		return
	}
	defer file.Close()

	var snp_pos uint32
	SNP_Pos := make([]int, 0, len(S.SNP_Prob))
	for snp_pos, _ = range S.SNP_Prob {
		SNP_Pos = append(SNP_Pos, int(snp_pos))
	}
	sort.Ints(SNP_Pos)

	var snp, snp_call, str_qual string
	var str_a, str_b string
	var line_a, line_b []string
	var snp_call_prob, snp_prob float64
	var snp_num, idx int
	for _, pos := range SNP_Pos {
		snp_pos = uint32(pos)
		snp_call_prob = 0
		for snp, snp_prob = range S.SNP_Prob[snp_pos] {
			if snp_call_prob < snp_prob {
				snp_call_prob = snp_prob
				snp_call = snp
			}
		}
		if snp_call == string(INDEX.SEQ[pos]) { //do not call SNPs which are identical with Ref
			continue
		}
		if len(snp_call) == 2 && snp_call[0] == snp_call[1] { //do not call homopolymer indels with length 2
			continue
		}
		line_a = make([]string, 0)
		line_a = append(line_a, strconv.Itoa(pos + 1))
		if _, ok := S.SNP_Type[snp_pos][snp_call]; ok {
			if S.SNP_Type[snp_pos][snp_call][0] == 2 { //do not call DEL at this stage
				continue
			} else {
				line_a = append(line_a, snp_call)
			}
		} else {
			line_a = append(line_a, snp_call)
		}
		str_qual = strconv.FormatFloat(-10 * math.Log10(1 - snp_call_prob), 'f', 5, 32)
		if str_qual != "+Inf" {
			line_a = append(line_a, str_qual)
		} else {
			line_a = append(line_a, "1000")
		}
		line_a = append(line_a, strconv.FormatFloat(snp_call_prob, 'f', 5, 32))
		line_a = append(line_a, strconv.Itoa(S.SNP_RNum[snp_pos][snp_call]))
		str_a = strings.Join(line_a, "\t")
		line_b = make([]string, 0)
		for snp, snp_num = range S.SNP_RNum[snp_pos] {
			line_b = append(line_b, snp)
			line_b = append(line_b, strconv.Itoa(snp_num))
		}
		str_b = strings.Join(line_b, "\t")
		for idx, _ = range S.Chr_Dis[snp_pos][snp_call] {
			_, err = file.WriteString(str_a + "\t")
			_, err = file.WriteString(string(S.SNP_BaseQ[snp_pos][snp_call][idx]) + "\t")
			_, err = file.WriteString(strconv.Itoa(S.Chr_Dis[snp_pos][snp_call][idx]) + "\t")
			_, err = file.WriteString(strconv.Itoa(S.Chr_Diff[snp_pos][snp_call][idx]) + "\t")
			_, err = file.WriteString(strconv.FormatFloat(S.Aln_Prob[snp_pos][snp_call][idx], 'f', 20, 64) + "\t")
			_, err = file.WriteString(strconv.FormatFloat(S.Chr_Prob[snp_pos][snp_call][idx], 'f', 20, 64) + "\t")
			_, err = file.WriteString(strconv.Itoa(S.Start_Pos1[snp_pos][snp_call][idx]) + "\t")
			_, err = file.WriteString(strconv.FormatBool(S.Strand1[snp_pos][snp_call][idx]) + "\t")
			_, err = file.WriteString(strconv.Itoa(S.Start_Pos2[snp_pos][snp_call][idx]) + "\t")
			_, err = file.WriteString(strconv.FormatBool(S.Strand2[snp_pos][snp_call][idx]) + "\t")
			_, err = file.WriteString(string(S.Read_Info[snp_pos][snp_call][idx]) + "\t")
			_, err = file.WriteString(str_b + "\n")
			if err != nil {
				fmt.Println(err)
				break
			}
		}
	}
}
