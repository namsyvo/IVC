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
	"log"
	"bytes"
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
	ADis    int    //aligment distance
	CDis    int    //chromosomal distance
	CDiff   int    //diff between aligned pos and true pos
	AProb   float64 //correct alignment prob
	CProb  float64 //correct paired-alignment prob
	RID    []byte  //read ID of aligned bases at SNP pos
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
	Aln_Dis     map[uint32]map[string][]int
	Chr_Dis     map[uint32]map[string][]int
	Chr_Diff     map[uint32]map[string][]int
	Aln_Prob     map[uint32]map[string][]float64
	Chr_Prob     map[uint32]map[string][]float64
	Read_ID      map[uint32]map[string][][]byte
}

//--------------------------------------------------------------------------------------------------
// InitIndex initializes indexes and parameters.
// This function will be called from main program.
//--------------------------------------------------------------------------------------------------
func (S *SNP_Prof) Init(input_info InputInfo) {
	INPUT_INFO = input_info
	//SetPara: 100 is maximum length of reads, 500 is maximum length of info line of reads,
	//1000 is maximum insert size of paired-end simulated reads, 0.0015 is maximum sequencing error rate
	//of simulated reads, 0.01 is mutation rate (currently is estimated from dbSNP of human genome)
	PARA_INFO = *SetPara(100, 500, 700, 0.0015, 0.01)
	INDEX.Init()
	S.SNP_Calls = make(map[uint32]map[string]float64)
	S.SNP_Bases = make(map[uint32]map[string]int)
	S.Aln_Dis = make(map[uint32]map[string][]int)
	S.Chr_Dis = make(map[uint32]map[string][]int)
	S.Chr_Diff = make(map[uint32]map[string][]int)
	S.Aln_Prob = make(map[uint32]map[string][]float64)
	S.Chr_Prob = make(map[uint32]map[string][]float64)
	S.Read_ID = make(map[uint32]map[string][][]byte)
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
		GetAlignReadInfo()
	}()
	go func() {
		GetNoAlignReadInfo()
	}()
	//------------------------
	go func() {
		wg.Wait()
		close(snp_results)
		//------------------------
		//For debugging
		close(ALIGN_READ_INFO_CHAN)
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
	//Output SNP calls
	fmt.Println("Outputing SNP calls...")
	S.OutputSNPCalls()

	//------------------------
	//For debugging
	ProcessNoAlignReadInfo(S.SNP_Calls)
	ProcessFNSNPInfo(S.SNP_Calls)
	ProcessTPFPSNPInfo(S.SNP_Calls)
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
		RevComp(read_info.Read1, read_info.Qual1, read_info.Rev_read1, read_info.Rev_comp_read1, read_info.Comp_read1, read_info.Rev_qual1)
		PrintMemStats("After calculating RevComp for Read1")
		RevComp(read_info.Read2, read_info.Qual2, read_info.Rev_read2, read_info.Rev_comp_read2, read_info.Comp_read2, read_info.Rev_qual2)
		PrintMemStats("After calculating RevComp for Read2")

		S.FindSNPsFromReads(read_info, snp_results, align_info, m_pos)
		PrintMemStats("After finding all SNPs from reads")
	}
}

//--------------------------------------------------------------------------------------------------
// FindSNPsFromReads returns SNPs found from alignment between pair-end reads and the multigenome.
// This version treats each end of the reads independently.
//--------------------------------------------------------------------------------------------------
func (S *SNP_Prof) FindSNPsFromReads(read_info *ReadInfo, snp_results chan SNP, align_info *AlignInfo, m_pos []int) {

	var snp SNP
	var snps1, snps2 []SNP
	var snps_get1, snps_get2 []SNP
	var m_dis1, l_align_pos1, r_align_pos1, m_dis2, l_align_pos2, r_align_pos2 int

	var s_pos_r1, e_pos_r1, s_pos_r2, e_pos_r2, m_pos_r1, m_pos_r2 []int
	var strand_r1, strand_r2 []bool

	var at Align_trace_info
	at.read1 = make([]byte, len(read_info.Read1))
	at.read2 = make([]byte, len(read_info.Read2))
	copy(at.read1, read_info.Read1)
	copy(at.read2, read_info.Read2)
	at.read_info1 = make([]byte, len(read_info.Info1))
	at.read_info2 = make([]byte, len(read_info.Info2))
	copy(at.read_info1, read_info.Info1)
	copy(at.read_info2, read_info.Info2)

    read_info1_tokens := bytes.Split(at.read_info1, []byte{'_'})
    read_info2_tokens := bytes.Split(at.read_info2, []byte{'_'})
	var true_pos1, true_pos2 int64
	var read_id1, read_id2 []byte
    if read_info1_tokens[0][1] != 'r' {
        true_pos1, _ = strconv.ParseInt(string(read_info1_tokens[2]), 10, 64)
		true_pos2, _ = strconv.ParseInt(string(read_info1_tokens[3]), 10, 64)
		read_id1, read_id2 = read_info1_tokens[10], read_info2_tokens[10]
	}
	
	var p_idx, s_idx int
	var p_prob, m_prob1, m_prob2 float64

	//Try to align both ends
	loop_num := 1
	for loop_num <= PARA_INFO.Iter_num { //temp value, will be replaced later
		PrintLoopTraceInfo(loop_num, "FindSNPsFromReads")
		s_pos_r1, e_pos_r1, s_pos_r2, e_pos_r2, m_pos_r1, m_pos_r2, strand_r1, strand_r2, _ = S.FindSeedsFromPairedEnds(read_info)
		p_prob = 0.0
		for p_idx = 0; p_idx < len(s_pos_r1); p_idx++ {
			
			//For conventional paired-end sequencing (i.e. Illumina) the directions should be F-R
			//For other kinds of variants (e.g inversions) or other technologies, they can be F-F or R-R
			//For mate-pair, thewy can be R-F (need to be confirmed)
			//if strand_r1[p_idx] == strand_r2[p_idx] {
			//	continue
			//}
			
			//Find SNPs for the first end
			PrintMemStats("Before FindSNPsFromEnd1")
			if strand_r1[p_idx] == true {
				snps1, l_align_pos1, r_align_pos1, m_prob1 = S.FindSNPsFromExtension(s_pos_r1[p_idx], e_pos_r1[p_idx], 
					m_pos_r1[p_idx], read_info.Read1, read_info.Qual1, align_info)
			} else {
				snps1, l_align_pos1, r_align_pos1, m_prob1 = S.FindSNPsFromExtension(s_pos_r1[p_idx], e_pos_r1[p_idx], 
					m_pos_r1[p_idx], read_info.Rev_comp_read1, read_info.Rev_qual1, align_info)
			}
			PrintMemStats("After FindSNPsFromEnd1")
			
			//Find SNPs for the second end
			PrintMemStats("Before FindSNPsFromEnd2")
			if strand_r2[p_idx] == true {
				snps2, l_align_pos2, r_align_pos2, m_prob2 = S.FindSNPsFromExtension(s_pos_r2[p_idx], e_pos_r2[p_idx], 
					m_pos_r2[p_idx], read_info.Read2, read_info.Qual2, align_info)
			} else {
				snps2, l_align_pos2, r_align_pos2, m_prob2 = S.FindSNPsFromExtension(s_pos_r2[p_idx], e_pos_r2[p_idx], 
					m_pos_r2[p_idx], read_info.Rev_comp_read2, read_info.Rev_qual2, align_info)
			}
			PrintMemStats("After FindSNPsFromEnd2")
			
			if m_prob1 != -1 && m_prob2 != -1 {
				a_prob := math.Exp(-math.Pow(math.Abs(float64(l_align_pos1 - l_align_pos2)) - 400.0, 2.0) / (2*50*50) )
				m_prob := -math.Log10(m_prob1) - math.Log10(m_prob2)
				log.Printf("Out\t%1.30f\t%1.30f\t%1.30f\t%1.30f\t%1.30f", m_prob1, m_prob2, a_prob, m_prob, p_prob)
				if p_prob > m_prob {
					log.Printf("In \t%1.30f\t%1.30f\t%1.30f\t%1.30f\t%1.30f", m_prob1, m_prob2, a_prob, m_prob, p_prob)
					//fmt.Println("Min p_dis", loop_num, p_dis, m_dis1, m_dis2)
					p_prob = m_prob
					at.l_align_pos1 = l_align_pos1
					at.l_align_pos2 = l_align_pos2
					at.r_align_pos1 = r_align_pos1
					at.r_align_pos2 = r_align_pos2
					at.align_dis1 = m_dis1
					at.align_dis2 = m_dis2
					
					snps_get1 = make([]SNP, len(snps1))
					if len(snps1) > 0 {
						for s_idx = 0; s_idx < len(snps1); s_idx++ {
							snps_get1[s_idx].Pos = snps1[s_idx].Pos
							snps_get1[s_idx].Bases = make([]byte, len(snps1[s_idx].Bases))
							snps_get1[s_idx].BaseQ = make([]byte, len(snps1[s_idx].BaseQ))
							copy(snps_get1[s_idx].Bases, snps1[s_idx].Bases)
							copy(snps_get1[s_idx].BaseQ, snps1[s_idx].BaseQ)
							snps_get1[s_idx].ADis = m_dis1
							snps_get1[s_idx].CDis = l_align_pos1 - l_align_pos2
							snps_get1[s_idx].CDiff = l_align_pos1 - int(true_pos1)
							snps_get1[s_idx].AProb = m_prob1
							snps_get1[s_idx].CProb = a_prob
							snps_get1[s_idx].RID = read_id1
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
							snps_get2[s_idx].ADis = m_dis2
							snps_get2[s_idx].CDis = l_align_pos1 - l_align_pos2
							snps_get2[s_idx].CDiff = l_align_pos2 - int(true_pos2)
							snps_get2[s_idx].AProb = m_prob2
							snps_get2[s_idx].CProb = a_prob
							snps_get2[s_idx].RID = read_id2
						}
					}
				}
			}
		}
		if p_prob < 0.0 {
			//fmt.Println("Get SNP", loop_num, p_dis, len(snps_get1), len(snps_get2))
			if len(snps_get1) > 0 {
				for _, snp = range snps_get1 {
					snp_results <- snp
					at.snp_pos1 = append(at.snp_pos1, snp.Pos)
					at.snp_base1 = append(at.snp_base1, snp.Bases)
					at.snp_baseq1 = append(at.snp_baseq1, snp.BaseQ)
				}
			}
			if len(snps_get2) > 0 {
				for _, snp = range snps_get2 {
					snp_results <- snp
					at.snp_pos2 = append(at.snp_pos2, snp.Pos)
					at.snp_base2 = append(at.snp_base2, snp.Bases)
					at.snp_baseq2 = append(at.snp_baseq2, snp.BaseQ)
				}
			}
			ALIGN_READ_INFO_CHAN <- at
			return
		}
		loop_num++
	}
	/*
	//Try to align the first end
	loop_num = 1
	for loop_num <= PARA_INFO.Iter_num { //temp value, will be replaced later
		PrintLoopTraceInfo(loop_num, "FindSNPsFromReads")
		s_pos_r1, e_pos_r1, s_pos_r2, e_pos_r2, m_pos_r1, m_pos_r2, strand_r1, strand_r2 = S.FindSeedsFromPairedEnds(read_info)
		p_dis = PARA_INFO.Dist_thres + 1 //temp value
		for p_idx = 0; p_idx < len(s_pos_r1); p_idx++ {
			//Find SNPs for the first end
			PrintMemStats("Before FindSNPsFromEnd1")
			if strand_r1[p_idx] == true {
				snps1, m_dis1, l_align_pos1, r_align_pos1 = S.FindSNPsFromExtension(s_pos_r1[p_idx], e_pos_r1[p_idx], 
					m_pos_r1[p_idx], read_info.Read1, read_info.Qual1, align_info)
			} else {
				snps1, m_dis1, l_align_pos1, r_align_pos1 = S.FindSNPsFromExtension(s_pos_r1[p_idx], e_pos_r1[p_idx], 
					m_pos_r1[p_idx], read_info.Rev_comp_read1, read_info.Rev_qual1, align_info)
			}
			PrintMemStats("After FindSNPsFromEnd1")

			if m_dis1 != -1 {
				if p_dis > m_dis1 {
					//fmt.Println("Min p_dis", loop_num, p_dis, m_dis1)
					p_dis = m_dis1
					at.l_align_pos1 = l_align_pos1
					at.l_align_pos2 = -1
					at.r_align_pos1 = r_align_pos1
					at.r_align_pos2 = -1
					at.align_dis1 = m_dis1
					at.align_dis2 = -1

					snps_get1 = make([]SNP, len(snps1))
					if len(snps1) > 0 {
						for s_idx = 0; s_idx < len(snps1); s_idx++ {
							snps_get1[s_idx].Pos = snps1[s_idx].Pos
							snps_get1[s_idx].Bases = make([]byte, len(snps1[s_idx].Bases))
							snps_get1[s_idx].BaseQ = make([]byte, len(snps1[s_idx].BaseQ))
							copy(snps_get1[s_idx].Bases, snps1[s_idx].Bases)
							copy(snps_get1[s_idx].BaseQ, snps1[s_idx].BaseQ)
						}
					}
				}
			}
		}
		if p_dis <= PARA_INFO.Dist_thres { //temp value
			//fmt.Println("Get SNP", loop_num, p_dis, len(snps_get1), len(snps_get2))
			if len(snps_get1) > 0 {
				for _, snp = range snps_get1 {
					snp_results <- snp
					at.snp_pos1 = append(at.snp_pos1, snp.Pos)
					at.snp_base1 = append(at.snp_base1, snp.Bases)
					at.snp_baseq1 = append(at.snp_baseq1, snp.BaseQ)
				}
			}
			ALIGN_READ_INFO_CHAN <- at
			return
		}
		loop_num++
	}

	//Try to align the second end
	loop_num = 1
	for loop_num <= PARA_INFO.Iter_num { //temp value, will be replaced later
		PrintLoopTraceInfo(loop_num, "FindSNPsFromReads")
		s_pos_r1, e_pos_r1, s_pos_r2, e_pos_r2, m_pos_r1, m_pos_r2, strand_r1, strand_r2 = S.FindSeedsFromPairedEnds(read_info)
		p_dis = PARA_INFO.Dist_thres + 1
		for p_idx = 0; p_idx < len(s_pos_r1); p_idx++ {
			//Find SNPs for the second end
			PrintMemStats("Before FindSNPsFromEnd2")
			if strand_r2[p_idx] == true {
				snps2, m_dis2, l_align_pos2, r_align_pos2 = S.FindSNPsFromExtension(s_pos_r2[p_idx], e_pos_r2[p_idx], 
					m_pos_r2[p_idx], read_info.Read2, read_info.Qual2, align_info)
			} else {
				snps2, m_dis2, l_align_pos2, r_align_pos2 = S.FindSNPsFromExtension(s_pos_r2[p_idx], e_pos_r2[p_idx], 
					m_pos_r2[p_idx], read_info.Rev_comp_read2, read_info.Rev_qual2, align_info)
			}
			PrintMemStats("After FindSNPsFromEnd2")

			if m_dis2 != -1 {
				if p_dis > m_dis2 {
					//fmt.Println("Min p_dis", loop_num, p_dis, m_dis1)
					p_dis = m_dis1
					at.l_align_pos1 = -1
					at.l_align_pos2 = l_align_pos2
					at.r_align_pos1 = -1
					at.r_align_pos2 = r_align_pos2
					at.align_dis1 = -1
					at.align_dis2 = m_dis2

					snps_get2 = make([]SNP, len(snps2))
					if len(snps2) > 0 {
						for s_idx = 0; s_idx < len(snps2); s_idx++ {
							snps_get2[s_idx].Pos = snps2[s_idx].Pos
							snps_get2[s_idx].Bases = make([]byte, len(snps2[s_idx].Bases))
							snps_get2[s_idx].BaseQ = make([]byte, len(snps2[s_idx].BaseQ))
							copy(snps_get2[s_idx].Bases, snps2[s_idx].Bases)
							copy(snps_get2[s_idx].BaseQ, snps2[s_idx].BaseQ)
						}
					}
				}
			}
		}
		if p_dis <= PARA_INFO.Dist_thres {
			//fmt.Println("Get SNP", loop_num, p_dis, len(snps_get2), len(snps_get2))
			if len(snps_get2) > 0 {
				for _, snp = range snps_get2 {
					snp_results <- snp
					at.snp_pos2 = append(at.snp_pos2, snp.Pos)
					at.snp_base2 = append(at.snp_base2, snp.Bases)
					at.snp_baseq2 = append(at.snp_baseq2, snp.BaseQ)
				}
			}
			ALIGN_READ_INFO_CHAN <- at
			return
		}
		loop_num++
	}
	*/
	//Cannot align any ends, consider as unaligned reads
	at.l_align_pos1 = -1
	at.l_align_pos2 = -1
	at.r_align_pos1 = -1
	at.r_align_pos2 = -1
	at.align_dis1 = -1
	at.align_dis2 = -1
	NO_ALIGN_READ_INFO_CHAN <- at
}

//---------------------------------------------------------------------------------------------------
// FindSeedsFromPairedEnds find all pairs of seeds which have proper chromosome distances.
//---------------------------------------------------------------------------------------------------
func (S *SNP_Prof) FindSeedsFromPairedEnds(read_info *ReadInfo) ([]int, []int, []int, []int, []int, []int, []bool, []bool, bool) {

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
			PrintExtendTraceInfo("r1_or", read_info.Read1[e_pos_r1_or : s_pos_r1_or + 1], e_pos_r1_or, s_pos_r1_or, m_num_r1_or, m_pos_r1_or)
		}
		s_pos_r1_rc, e_pos_r1_rc, m_num_r1_rc, has_seeds_r1_rc = 
			INDEX.FindSeeds(read_info.Rev_comp_read1, read_info.Comp_read1, r_pos_r1_rc, m_pos_r1_rc)
		PrintSeedTraceInfo("r1_rc", e_pos_r1_rc, s_pos_r1_rc, read_info.Rev_comp_read1)
		if has_seeds_r1_rc {
			PrintExtendTraceInfo("r1_rc", read_info.Rev_comp_read1[e_pos_r1_rc : s_pos_r1_rc + 1], e_pos_r1_rc, s_pos_r1_rc, m_num_r1_rc, m_pos_r1_rc)
		}
		s_pos_r2_or, e_pos_r2_or, m_num_r2_or, has_seeds_r2_or = 
			INDEX.FindSeeds(read_info.Read2, read_info.Rev_read2, r_pos_r2_or, m_pos_r2_or)
		PrintSeedTraceInfo("r2_or", e_pos_r2_or, s_pos_r2_or, read_info.Read2)
		if has_seeds_r2_or {
			PrintExtendTraceInfo("r2_or", read_info.Read1[e_pos_r2_or : s_pos_r2_or + 1], e_pos_r2_or, s_pos_r2_or, m_num_r2_or, m_pos_r2_or)
		}
		s_pos_r2_rc, e_pos_r2_rc, m_num_r2_rc, has_seeds_r2_rc = 
			INDEX.FindSeeds(read_info.Rev_comp_read2, read_info.Comp_read2, r_pos_r2_rc, m_pos_r2_rc)
		PrintSeedTraceInfo("r2_rc", e_pos_r2_rc, s_pos_r2_rc, read_info.Rev_comp_read2)
		if has_seeds_r2_rc {
			PrintExtendTraceInfo("r2_rc", read_info.Rev_comp_read2[e_pos_r2_rc : s_pos_r2_rc + 1], e_pos_r2_rc, s_pos_r2_rc, m_num_r2_rc, m_pos_r2_rc)
		}
		PrintMemStats("After FindSeeds, loop_num " + strconv.Itoa(loop_num))

		if has_seeds_r1_or && has_seeds_r2_rc {
			PrintExtendTraceInfo("r1_or(F1R2)", read_info.Read1[e_pos_r1_or : s_pos_r1_or + 1], e_pos_r1_or, s_pos_r1_or, m_num_r1_or, m_pos_r1_or)
			PrintExtendTraceInfo("r2_rc(F1R2)", read_info.Read1[e_pos_r2_rc : s_pos_r2_rc + 1], e_pos_r2_rc, s_pos_r2_rc, m_num_r2_rc, m_pos_r2_rc)
			for i = 0; i < m_num_r1_or; i++ {
				for j = 0; j < m_num_r2_rc; j++ {
					//Check if alignments are likely pair-end alignments
					if (m_pos_r1_or[i] < m_pos_r2_rc[j]) && (m_pos_r2_rc[j] - m_pos_r1_or[i]) <= PARA_INFO.Max_ins {
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
			PrintExtendTraceInfo("r1_rc (F2R1)", read_info.Read1[e_pos_r1_rc : s_pos_r1_rc + 1], e_pos_r1_rc, s_pos_r1_rc, m_num_r1_rc, m_pos_r1_rc)
			PrintExtendTraceInfo("r2_or (F2R1)", read_info.Read1[e_pos_r2_or : s_pos_r2_or + 1], e_pos_r2_or, s_pos_r2_or, m_num_r2_or, m_pos_r2_or)
			for i = 0; i < m_num_r1_rc; i++ {
				for j = 0; j < m_num_r2_or; j++ {
					//Check if alignments are likely pair-end alignments
					if (m_pos_r1_rc[i] > m_pos_r2_or[j]) && (m_pos_r1_rc[i] - m_pos_r2_or[j]) <= PARA_INFO.Max_ins {
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
		/*
		if has_seeds_r1_or && has_seeds_r2_or {
			PrintExtendTraceInfo("r1_or", read_info.Read1[e_pos_r1_or : s_pos_r1_or + 1], e_pos_r1_or, s_pos_r1_or, m_num_r1_or, m_pos_r1_or)
			PrintExtendTraceInfo("r2_or", read_info.Read1[e_pos_r2_or : s_pos_r2_or + 1], e_pos_r2_or, s_pos_r2_or, m_num_r2_or, m_pos_r2_or)
			for i = 0; i < m_num_r1_or; i++ {
				for j = 0; j < m_num_r2_or; j++ {
					//Check if alignments are likely pair-end alignments
					if int(math.Abs(float64(m_pos_r1_or[i] - m_pos_r2_or[j]))) <= PARA_INFO.Max_ins {
						PrintPairedSeedInfo("r1_or, r2_or, paired pos", m_pos_r1_or[i], m_pos_r2_or[j])
						s_pos_r1 = append(s_pos_r1, s_pos_r1_or)
						e_pos_r1 = append(e_pos_r1, e_pos_r1_or)
						s_pos_r2 = append(s_pos_r2, s_pos_r2_or)
						e_pos_r2 = append(e_pos_r2, e_pos_r2_or)
						m_pos_r1 = append(m_pos_r1, m_pos_r1_or[i])
						m_pos_r2 = append(m_pos_r2, m_pos_r2_or[j])
						strand_r1 = append(strand_r1, true)
						strand_r2 = append(strand_r2, true)
					}
				}
			}
		}
		if has_seeds_r1_rc && has_seeds_r2_rc {
			PrintExtendTraceInfo("r1_rc", read_info.Read1[e_pos_r1_rc : s_pos_r1_rc + 1], e_pos_r1_rc, s_pos_r1_rc, m_num_r1_rc, m_pos_r1_rc)
			PrintExtendTraceInfo("r2_rc", read_info.Read1[e_pos_r2_rc : s_pos_r2_rc + 1], e_pos_r2_rc, s_pos_r2_rc, m_num_r2_rc, m_pos_r2_rc)
			for i = 0; i < m_num_r1_rc; i++ {
				for j = 0; j < m_num_r2_rc; j++ {
					//Check if alignments are likely pair-end alignments
					if int(math.Abs(float64(m_pos_r1_rc[i] - m_pos_r2_rc[j]))) <= PARA_INFO.Max_ins {
						PrintPairedSeedInfo("r1_rc, r2_rc, paired pos", m_pos_r1_rc[i], m_pos_r2_rc[j])
						s_pos_r1 = append(s_pos_r1, s_pos_r1_rc)
						e_pos_r1 = append(e_pos_r1, e_pos_r1_rc)
						s_pos_r2 = append(s_pos_r2, s_pos_r2_rc)
						e_pos_r2 = append(e_pos_r2, e_pos_r2_rc)
						m_pos_r1 = append(m_pos_r1, m_pos_r1_rc[i])
						m_pos_r2 = append(m_pos_r2, m_pos_r2_rc[j])
						strand_r1 = append(strand_r1, false)
						strand_r2 = append(strand_r2, false)
					}
				}
			}
		}
		*/
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

//---------------------------------------------------------------------------------------------------
// FindSNPsFromEachEnd find SNPs from alignment between read (one end) and multigenome.
//---------------------------------------------------------------------------------------------------
func (S *SNP_Prof) FindSNPsFromExtension(s_pos, e_pos, m_pos int, read, qual []byte, 
	align_info *AlignInfo) ([]SNP, int, int, float64) {

	var k int
	var prob float64

	var l_most_pos, r_most_pos int
	var l_snp_pos, r_snp_pos, l_snp_idx, r_snp_idx []int
	var l_snp_val, r_snp_val [][]byte

	var snp SNP
	var snps_arr []SNP
	var has_match bool

	PrintMemStats("Before FindExtensions, m_pos " + strconv.Itoa(m_pos))
	prob, l_snp_pos, l_snp_val, l_snp_idx, r_snp_pos, r_snp_val, r_snp_idx, l_most_pos, r_most_pos, has_match =
		S.FindExtensions(read, qual, s_pos, e_pos, m_pos, align_info)
	PrintMemStats("After FindExtensions, m_pos " + strconv.Itoa(m_pos))
	if has_match {
		PrintMatchTraceInfo(m_pos, 0, l_most_pos, l_snp_pos, read)
		for k = 0; k < len(l_snp_pos); k++ {
			PrintMemStats("Before GetSNP left, snp_num " + strconv.Itoa(k))
			l_snp_qual := make([]byte, len(l_snp_val[k]))
			copy(l_snp_qual, qual[l_snp_idx[k] : l_snp_idx[k] + len(l_snp_val[k])])
			snp.Pos, snp.Bases, snp.BaseQ = uint32(l_snp_pos[k]), l_snp_val[k], l_snp_qual
			snps_arr = append(snps_arr, snp)
			for _, q := range(snp.BaseQ) {
				if _, snp_call_exist := S.SNP_Calls[snp.Pos]; !snp_call_exist {
					prob = prob + math.Log10(1.0 - math.Pow(10, -(float64(q) - 33) / 10.0)) + math.Log10(EPSILON)
				} else {
					prob = prob + math.Log10(1.0 - math.Pow(10, -(float64(q) - 33) / 10.0)) + math.Log10(S.SNP_Calls[snp.Pos][string(snp.Bases)])
				}
			}
			PrintMemStats("After GetSNP left, snp_num " + strconv.Itoa(k))
		}
		for k = 0; k < len(r_snp_pos); k++ {
			PrintMemStats("Before GetSNP right, snp_num " + strconv.Itoa(k))
			r_snp_qual := make([]byte, len(r_snp_val[k]))
			copy(r_snp_qual, qual[r_snp_idx[k] : r_snp_idx[k] + len(r_snp_val[k])])
			snp.Pos, snp.Bases, snp.BaseQ = uint32(r_snp_pos[k]), r_snp_val[k], r_snp_qual
			snps_arr = append(snps_arr, snp)
			for _, q := range(snp.BaseQ) {
				if _, snp_call_exist := S.SNP_Calls[snp.Pos]; !snp_call_exist {
					prob = prob + math.Log10(1.0 - math.Pow(10, -(float64(q) - 33) / 10.0)) + math.Log10(EPSILON)
				} else {
					prob = prob + math.Log10(1.0 - math.Pow(10, -(float64(q) - 33) / 10.0)) + math.Log10(S.SNP_Calls[snp.Pos][string(snp.Bases)])
				}
			}
			PrintMemStats("After GetSNP right, snp_num " + strconv.Itoa(k))
		}
		return snps_arr, l_most_pos, r_most_pos, prob
	}
	return snps_arr, -1, -1, -1.0
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
		S.Aln_Dis[pos] = make(map[string][]int)
		S.Chr_Dis[pos] = make(map[string][]int)
		S.Chr_Diff[pos] = make(map[string][]int)
		S.Aln_Prob[pos] = make(map[string][]float64)
		S.Chr_Prob[pos] = make(map[string][]float64)
		S.Read_ID[pos] = make(map[string][][]byte)
	}
	S.SNP_Bases[pos][a] += 1
	S.Aln_Dis[pos][a] = append(S.Aln_Dis[pos][a], snp.ADis)
	S.Chr_Dis[pos][a] = append(S.Chr_Dis[pos][a], snp.CDis)
	S.Chr_Diff[pos][a] = append(S.Chr_Diff[pos][a], snp.CDiff)
	S.Aln_Prob[pos][a] = append(S.Aln_Prob[pos][a], snp.AProb)
	S.Chr_Prob[pos][a] = append(S.Chr_Prob[pos][a], snp.CProb)
	S.Read_ID[pos][a] = append(S.Read_ID[pos][a], snp.RID)

	var p float64
	p_ab := make(map[string]float64)
	p_a := 0.0

	if _, snp_call_exist := S.SNP_Calls[pos]; !snp_call_exist {
		S.SNP_Calls[pos] = make(map[string]float64)
		if snps, snp_prof_exist := INDEX.SNP_PROF[int(pos)]; snp_prof_exist {
			snp_prof_num := len(snps)
			std_base_num := len(STD_BASES)
			for idx, snp := range snps {
				S.SNP_Calls[pos][string(snp)] = float64(INDEX.SNP_AF[int(pos)][idx]) -  EPSILON * float64(std_base_num - snp_prof_num)/float64(snp_prof_num)
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
		S.Aln_Dis[pos] = make(map[string][]int)
		S.Chr_Dis[pos] = make(map[string][]int)
		S.Chr_Diff[pos] = make(map[string][]int)
		S.Aln_Prob[pos] = make(map[string][]float64)
		S.Chr_Prob[pos] = make(map[string][]float64)
		S.Read_ID[pos] = make(map[string][][]byte)
	}
	S.SNP_Bases[pos][a] += 1
	S.Aln_Dis[pos][a] = append(S.Aln_Dis[pos][a], snp.ADis)
	S.Chr_Dis[pos][a] = append(S.Chr_Dis[pos][a], snp.CDis)
	S.Chr_Diff[pos][a] = append(S.Chr_Diff[pos][a], snp.CDiff)
	S.Aln_Prob[pos][a] = append(S.Aln_Prob[pos][a], snp.AProb)
	S.Chr_Prob[pos][a] = append(S.Chr_Prob[pos][a], snp.CProb)
	S.Read_ID[pos][a] = append(S.Read_ID[pos][a], snp.RID)

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
	SNP_Pos := make([]int, 0, len(S.SNP_Calls))
	for snp_pos, _ = range S.SNP_Calls {
		SNP_Pos = append(SNP_Pos, int(snp_pos))
	}
	sort.Ints(SNP_Pos)

	var str_qual string
	var line_a, line_b []string
	var str_a, str_b string
	var snp_call_prob, snp_prob float64
	var snp_call, snp string
	var snp_num, aln_dis, idx int
	for _, pos := range SNP_Pos {
		snp_pos = uint32(pos)
		snp_call_prob = 0
		for snp, snp_prob = range S.SNP_Calls[snp_pos] {
			if snp_call_prob < snp_prob {
				snp_call_prob = snp_prob
				snp_call = snp
			}
		}
		line_a = make([]string, 0)
		line_a = append(line_a, strconv.Itoa(pos + 1))
		line_a = append(line_a, snp_call)
		str_qual = strconv.FormatFloat(-10 * math.Log10(1 - snp_call_prob), 'f', 5, 32)
		if str_qual != "+Inf" {
			line_a = append(line_a, str_qual)
		} else {
			line_a = append(line_a, "1000")
		}
		line_a = append(line_a, strconv.FormatFloat(snp_call_prob, 'f', 5, 32))
		line_a = append(line_a, strconv.Itoa(S.SNP_Bases[snp_pos][snp_call]))
		str_a = strings.Join(line_a, "\t")
		line_b = make([]string, 0)
		for snp, snp_num = range S.SNP_Bases[snp_pos] {
			line_b = append(line_b, snp)
			line_b = append(line_b, strconv.Itoa(snp_num))
		}
		str_b = strings.Join(line_b, "\t")
		for idx, aln_dis = range S.Aln_Dis[snp_pos][snp_call] {
			_, err = file.WriteString(str_a + "\t")
			_, err = file.WriteString(strconv.Itoa(aln_dis) + "\t")
			_, err = file.WriteString(strconv.Itoa(S.Chr_Dis[snp_pos][snp_call][idx]) + "\t")
			_, err = file.WriteString(strconv.Itoa(S.Chr_Diff[snp_pos][snp_call][idx]) + "\t")
			_, err = file.WriteString(strconv.FormatFloat(S.Aln_Prob[snp_pos][snp_call][idx], 'f', 20, 64) + "\t")
			_, err = file.WriteString(strconv.FormatFloat(S.Chr_Prob[snp_pos][snp_call][idx], 'f', 20, 64) + "\t")
			_, err = file.WriteString(string(S.Read_ID[snp_pos][snp_call][idx]) + "\t")
			_, err = file.WriteString(str_b + "\n")
			if err != nil {
				fmt.Println(err)
				break
			}
		}
	}
}