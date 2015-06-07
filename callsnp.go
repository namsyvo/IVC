//---------------------------------------------------------------------------------------------------
// IVC - Calling genomic variants based on alignment between reads and the reference multigenome.
// Variants and probability of correct variant calls is determined using Bayesian update.
// Copyright 2015 Nam Sy Vo.
//---------------------------------------------------------------------------------------------------

package isc

import (
	"bufio"
	"bytes"
	"fmt"
	"math"
	"math/rand"
	"os"
	"sort"
	"strconv"
	"strings"
	"sync"
	"time"
)

/*--------------------------------------------------------------------------------------------------
VarInfo represents variants obtained during alignment phase.
It serves as temporary variable during variant calling phase.
--------------------------------------------------------------------------------------------------*/
type VarInfo struct {
	Pos   uint32  //postion of variants on ref
	Bases []byte  //bases of variants
	BaseQ []byte  //quality of bases of variants
	Type  int     //type of variants (sub, ins, del...)
	CDis  int     //chromosomal distance of alignment of two read-ends
	CDiff int     //difference between aligned pos and true pos
	AProb float64 //prob of correct alignment
	CProb float64 //prob correct paired-mapping
	RInfo []byte  //info of reads
	SPos1 int     //starting pos on read1 of exact match (i.e. ending position from backward search with FM-index)
	SPos2 int     //starting pos on read2 of exact match
	Stra1 bool    //strand (backward/forward) of read1 of exact match
	Stra2 bool    //strand (backward/forward) of read2 of exact match
}

/*--------------------------------------------------------------------------------------------------
Var_Prof represents info of aligned reads and bases at the variant call positions on reference multigenome.
This struct also has functions defined on it for calling variants.
--------------------------------------------------------------------------------------------------*/
type Var_Prof struct {
	/*
		Var_Prob stores all possible variants at each position and their confident probablilities.
			Their prior probablities will be obtained from reference genomes and variant profiles.
			Their posterior probabilities will be updated during alignment phase based on incomming aligned bases
	*/
	Var_Prob   map[uint32]map[string]float64   //Probability of variant calls
	Var_BaseQ  map[uint32]map[string][][]byte  //Quality sequences (in FASTQ format) of aligned bases at the variant call position
	Var_Type   map[uint32]map[string][]int     //Type of variants (currently: 0:sub, 1:ins, 2:del)
	Var_RNum   map[uint32]map[string]int       //Numer of reads (bases) aligned to the variant call postion
	Chr_Dis    map[uint32]map[string][]int     //Chromosomal distance between two aligned ends
	Chr_Diff   map[uint32]map[string][]int     //Difference betwwen aligned postions and true postions (for simulated data)
	Aln_Prob   map[uint32]map[string][]float64 //Alignment probability of aligned read at the variant call position
	Chr_Prob   map[uint32]map[string][]float64 //Paired-end mapping probability
	Start_Pos1 map[uint32]map[string][]int     //Start position of alignment of the first end
	Start_Pos2 map[uint32]map[string][]int     //Start position of alignment of the second end
	Strand1    map[uint32]map[string][]bool    //Strand indicator of the first end ("true" if identicaly with ref)
	Strand2    map[uint32]map[string][]bool    //Strand indicator of the second end
	Read_Info  map[uint32]map[string][][]byte  //Info of aligned reads at the variant call position
}

/*--------------------------------------------------------------------------------------------------
InitIndex initializes indexes and parameters.
This function will be called from main program.
--------------------------------------------------------------------------------------------------*/
func NewVariantCaller(input_info InputInfo) *Var_Prof {

	INPUT_INFO = input_info
	if _, err := os.Stat(INPUT_INFO.Read_file_1); err != nil {
	    fmt.Println("Error: Read_file_1 does not exists!", err)
	    os.Exit(1)
	}
	if _, err := os.Stat(INPUT_INFO.Read_file_2); err != nil {
	    fmt.Println("Error: Read_file_2 does not exists!", err)
	    os.Exit(1)
	}
	/* SetPara: 100 is maximum length of reads, 500 is maximum length of info line of reads,
				700 is maximum insert size of paired-end simulated reads, 0.0015 is maximum sequencing error rate
				0.01 is mutation rate (currently is estimated from dbSNP of human genome)
	*/
	PARA_INFO = *SetPara(100, 500, 700, 0.0015, 0.01, input_info.Dist_thres, input_info.Iter_num)
	RAND_GEN = rand.New(rand.NewSource(time.Now().UnixNano()))

	INDEX = *NewIndex()

	S := new(Var_Prof)
	S.Var_Prob   = make(map[uint32]map[string]float64)
	S.Var_BaseQ  = make(map[uint32]map[string][][]byte)
	S.Var_Type   = make(map[uint32]map[string][]int)
	S.Var_RNum   = make(map[uint32]map[string]int)
	S.Chr_Dis    = make(map[uint32]map[string][]int)
	S.Chr_Diff   = make(map[uint32]map[string][]int)
	S.Aln_Prob   = make(map[uint32]map[string][]float64)
	S.Chr_Prob   = make(map[uint32]map[string][]float64)
	S.Read_Info  = make(map[uint32]map[string][][]byte)
	S.Start_Pos1 = make(map[uint32]map[string][]int)
	S.Start_Pos2 = make(map[uint32]map[string][]int)
	S.Strand1    = make(map[uint32]map[string][]bool)
	S.Strand2    = make(map[uint32]map[string][]bool)

	var pos uint32
	var var_bases []byte
	var i int
	for var_pos, var_prof := range INDEX.Var_Prof {
		pos = uint32(var_pos)
		S.Var_Prob[pos] = make(map[string]float64)
		//At this point, assume that all variants are biallelic
		if len(var_prof[0]) == 1 && len(var_prof[1]) == 1 {
			for i, var_bases = range var_prof {
				S.Var_Prob[pos][string(var_bases)] = float64(INDEX.Var_AF[var_pos][i]) - NEW_SNP_RATE
				if S.Var_Prob[pos][string(var_bases)] < NEW_SNP_RATE {
					S.Var_Prob[pos][string(var_bases)] = NEW_SNP_RATE
				}
			}
		} else {
			for i, var_bases = range var_prof {
				S.Var_Prob[pos][string(var_bases)] = float64(INDEX.Var_AF[var_pos][i]) - 1.5*NEW_SNP_RATE
				if S.Var_Prob[pos][string(var_bases)] < NEW_SNP_RATE {
					S.Var_Prob[pos][string(var_bases)] = NEW_SNP_RATE
				}
			}
		}
		for _, b := range STD_BASES {
			if _, ok := S.Var_Prob[pos][string(b)]; !ok {
				S.Var_Prob[pos][string(b)] = NEW_SNP_RATE
			}
		}
		S.Var_BaseQ[pos]  = make(map[string][][]byte)
		S.Var_Type[pos]   = make(map[string][]int)
		S.Var_RNum[pos]   = make(map[string]int)
		S.Chr_Dis[pos]    = make(map[string][]int)
		S.Chr_Diff[pos]   = make(map[string][]int)
		S.Aln_Prob[pos]   = make(map[string][]float64)
		S.Chr_Prob[pos]   = make(map[string][]float64)
		S.Read_Info[pos]  = make(map[string][][]byte)
		S.Start_Pos1[pos] = make(map[string][]int)
		S.Start_Pos2[pos] = make(map[string][]int)
		S.Strand1[pos]    = make(map[string][]bool)
		S.Strand2[pos]    = make(map[string][]bool)
	}
	return S
}

/*--------------------------------------------------------------------------------------------------
CallVariants initializes share variables, channels, reads input reads, finds all possible variants,
and updates variant information in Var_Prof.
This function will be called from main program.
--------------------------------------------------------------------------------------------------*/
func (S *Var_Prof) CallVariants() {
	//The channel read_signal is used for signaling between goroutines which run ReadReads and FindVariants,
	//when a FindSNPs goroutine finish copying a read to its own memory,
	//it signals ReadReads goroutine to scan next reads.
	read_signal := make(chan bool)

	//Call a goroutine to read input reads
	read_data := make(chan *ReadInfo, INPUT_INFO.Routine_num)
	go S.ReadReads(read_data, read_signal)

	//Call goroutines to find Vars, pass shared variable to each goroutine
	var_results := make(chan VarInfo)
	var wg sync.WaitGroup
	for i := 0; i < INPUT_INFO.Routine_num; i++ {
		go S.FindVariants(read_data, read_signal, var_results, &wg)
	}
	//------------------------
	//For debugging
	go func() {
		GetNoAlignReadInfo()
	}()
	//------------------------
	go func() {
		wg.Wait()
		close(var_results)
		//------------------------
		//For debugging
		close(NO_ALIGN_READ_INFO_CHAN)
		//------------------------
	}()

	//Collect Vars from results channel and update Vars and their probabilities
	var var_info VarInfo
	for var_info = range var_results {
		if var_info.Type == 0 {
			S.UpdateSNPProb(var_info)
		} else {
			S.UpdateIndelProb(var_info)
		}
	}
	//------------------------
	//For debugging
	ProcessNoAlignReadInfo(S.Var_Prob)
	//------------------------
}

/*--------------------------------------------------------------------------------------------------
ReadReads reads all reads from input FASTQ files and put them into data channel.
--------------------------------------------------------------------------------------------------*/
func (S *Var_Prof) ReadReads(read_data chan *ReadInfo, read_signal chan bool) {

	fn1, fn2 := INPUT_INFO.Read_file_1, INPUT_INFO.Read_file_2
	f1, err_f1 := os.Open(fn1)
	if err_f1 != nil {
	    fmt.Println("Error: Open read_file_1 " + fn1, err_f1)
	    os.Exit(1)
	}
	defer f1.Close()
	f2, err_f2 := os.Open(fn2)
	if err_f2 != nil {
	    fmt.Println("Error: Open read_file_2 " + fn1, err_f2)
	    os.Exit(1)
	}
	defer f2.Close()

	read_num := 0
	scanner1 := bufio.NewScanner(f1)
	scanner2 := bufio.NewScanner(f2)
	read_info := InitReadInfo(PARA_INFO.Read_len, PARA_INFO.Info_len)
	for scanner1.Scan() && scanner2.Scan() {
		read_info.Info1 = read_info.Info1[:PARA_INFO.Info_len]
		read_info.Info2 = read_info.Info2[:PARA_INFO.Info_len]
		copy(read_info.Info1, scanner1.Bytes()) //use 1st line in input FASTQ file 1
		copy(read_info.Info2, scanner2.Bytes()) //use 1st line in input FASTQ file 2
		read_info.Info1 = read_info.Info1[:len(scanner1.Bytes())]
		read_info.Info2 = read_info.Info2[:len(scanner2.Bytes())]
		scanner1.Scan()
		scanner2.Scan()
		copy(read_info.Read1, scanner1.Bytes()) //use 2nd line in input FASTQ file 1
		copy(read_info.Read2, scanner2.Bytes()) //use 2nd line in input FASTQ file 2
		scanner1.Scan()                         //ignore 3rd line in 1st input FASTQ file 1
		scanner2.Scan()                         //ignore 3rd line in 2nd input FASTQ file 2
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
FindVariants takes data from data channel, find all possible Vars and put them into results channel.
--------------------------------------------------------------------------------------------------*/
func (S *Var_Prof) FindVariants(read_data chan *ReadInfo, read_signal chan bool, var_results chan VarInfo,
	wg *sync.WaitGroup) {
	wg.Add(1)
	defer wg.Done()

	//Initialize inter-function share variables
	read_info := InitReadInfo(PARA_INFO.Read_len, PARA_INFO.Info_len)
	align_info := InitAlignInfo(2 * PARA_INFO.Read_len)
	m_pos := make([]int, INPUT_INFO.Max_snum)

	for read := range read_data {
		PrintMemStats("Before copying all info from data chan")
		read_info.Info1 = read_info.Info1[:PARA_INFO.Info_len]
		read_info.Info2 = read_info.Info2[:PARA_INFO.Info_len]
		copy(read_info.Info1, read.Info1)
		copy(read_info.Info2, read.Info2)
		read_info.Info1 = read_info.Info1[:len(read.Info1)]
		read_info.Info2 = read_info.Info2[:len(read.Info2)]
		copy(read_info.Read1, read.Read1)
		copy(read_info.Read2, read.Read2)
		copy(read_info.Qual1, read.Qual1)
		copy(read_info.Qual2, read.Qual2)
		<-read_signal

		PrintMemStats("After copying all info from data chan")
		RevComp(read_info.Read1, read_info.Qual1, read_info.Rev_read1, read_info.Rev_comp_read1,
			read_info.Comp_read1, read_info.Rev_qual1)
		PrintMemStats("After calculating RevComp for Read1")
		RevComp(read_info.Read2, read_info.Qual2, read_info.Rev_read2, read_info.Rev_comp_read2,
			read_info.Comp_read2, read_info.Rev_qual2)
		PrintMemStats("After calculating RevComp for Read2")

		S.FindVariantsFromPairedEnds(read_info, var_results, align_info, m_pos)
		PrintMemStats("After finding all Vars from reads")
	}
}

/*--------------------------------------------------------------------------------------------------
FindVariantsFromPairedEndReads returns Vars found from alignment between pair-end reads and the multigenome.
This version treats each end of the reads independently.
--------------------------------------------------------------------------------------------------*/
func (S *Var_Prof) FindVariantsFromPairedEnds(read_info *ReadInfo, var_results chan VarInfo, align_info *AlignInfo, m_pos []int) {

	var var_info VarInfo
	var vars1, vars2 []VarInfo
	var vars_get1, vars_get2 []VarInfo
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
		PrintLoopTraceInfo(loop_num, "FindVariantsFromReads")
		s_pos_r1, e_pos_r1, s_pos_r2, e_pos_r2, m_pos_r1, m_pos_r2, strand_r1, strand_r2, has_seeds = S.FindSeedsFromPairedEnds(read_info)
		if has_seeds {
			for p_idx = 0; p_idx < len(s_pos_r1); p_idx++ {
				//For conventional paired-end sequencing (i.e. Illumina) the directions should be F-R
				//For other kinds of variants (e.g inversions) or other technologies, they can be F-F or R-R
				//For mate-pair, they can be R-F (need to be confirmed)
				if strand_r1[p_idx] == strand_r2[p_idx] {
					continue
				}
				//Find Vars for the first end
				PrintMemStats("Before FindVariantsFromEnd1")
				if strand_r1[p_idx] == true {
					vars1, l_align_pos1, _, align_prob1 = S.FindVariantsFromExtension(s_pos_r1[p_idx], e_pos_r1[p_idx],
						m_pos_r1[p_idx], read_info.Read1, read_info.Qual1, align_info)
				} else {
					vars1, l_align_pos1, _, align_prob1 = S.FindVariantsFromExtension(s_pos_r1[p_idx], e_pos_r1[p_idx],
						m_pos_r1[p_idx], read_info.Rev_comp_read1, read_info.Rev_qual1, align_info)
				}
				PrintMemStats("After FindVariantsFromEnd1")

				//Find Vars for the second end
				PrintMemStats("Before FindVariantsFromEnd2")
				if strand_r2[p_idx] == true {
					vars2, l_align_pos2, _, align_prob2 = S.FindVariantsFromExtension(s_pos_r2[p_idx], e_pos_r2[p_idx],
						m_pos_r2[p_idx], read_info.Read2, read_info.Qual2, align_info)
				} else {
					vars2, l_align_pos2, _, align_prob2 = S.FindVariantsFromExtension(s_pos_r2[p_idx], e_pos_r2[p_idx],
						m_pos_r2[p_idx], read_info.Rev_comp_read2, read_info.Rev_qual2, align_info)
				}
				PrintMemStats("After FindVariantsFromEnd2")

				if align_prob1 != -1 && align_prob2 != -1 {
					a_prob := -math.Log10(math.Exp(-math.Pow(math.Abs(float64(l_align_pos1-l_align_pos2))-400.0, 2.0) / (2 * 50 * 50)))
					if paired_prob > align_prob1+align_prob2 {
						paired_prob = align_prob1 + align_prob2
						vars_get1 = make([]VarInfo, len(vars1))
						PrintGetVariants(paired_prob, align_prob1, align_prob2, vars1, vars2)
						if len(vars1) > 0 {
							for s_idx = 0; s_idx < len(vars1); s_idx++ {
								vars_get1[s_idx].Pos   = vars1[s_idx].Pos
								vars_get1[s_idx].Bases = make([]byte, len(vars1[s_idx].Bases))
								vars_get1[s_idx].BaseQ = make([]byte, len(vars1[s_idx].BaseQ))
								copy(vars_get1[s_idx].Bases, vars1[s_idx].Bases)
								copy(vars_get1[s_idx].BaseQ, vars1[s_idx].BaseQ)
								vars_get1[s_idx].Type  = vars1[s_idx].Type
								vars_get1[s_idx].CDis  = l_align_pos1 - l_align_pos2
								vars_get1[s_idx].CDiff = l_align_pos1 - int(true_pos1)
								vars_get1[s_idx].AProb = align_prob1
								vars_get1[s_idx].CProb = a_prob
								vars_get1[s_idx].RInfo = read_info1
								vars_get1[s_idx].SPos1 = e_pos_r1[p_idx]
								vars_get1[s_idx].SPos2 = e_pos_r2[p_idx]
								vars_get1[s_idx].Stra1 = strand_r1[p_idx]
								vars_get1[s_idx].Stra2 = strand_r2[p_idx]
							}
						}
						vars_get2 = make([]VarInfo, len(vars2))
						if len(vars2) > 0 {
							for s_idx = 0; s_idx < len(vars2); s_idx++ {
								vars_get2[s_idx].Pos   = vars2[s_idx].Pos
								vars_get2[s_idx].Bases = make([]byte, len(vars2[s_idx].Bases))
								vars_get2[s_idx].BaseQ = make([]byte, len(vars2[s_idx].BaseQ))
								copy(vars_get2[s_idx].Bases, vars2[s_idx].Bases)
								copy(vars_get2[s_idx].BaseQ, vars2[s_idx].BaseQ)
								vars_get2[s_idx].Type  = vars2[s_idx].Type
								vars_get2[s_idx].CDis  = l_align_pos1 - l_align_pos2
								vars_get2[s_idx].CDiff = l_align_pos2 - int(true_pos2)
								vars_get2[s_idx].AProb = align_prob2
								vars_get2[s_idx].CProb = a_prob
								vars_get2[s_idx].RInfo = read_info2
								vars_get2[s_idx].SPos1 = e_pos_r1[p_idx]
								vars_get2[s_idx].SPos2 = e_pos_r2[p_idx]
								vars_get2[s_idx].Stra1 = strand_r1[p_idx]
								vars_get2[s_idx].Stra2 = strand_r2[p_idx]
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
	if paired_prob <= 2*PARA_INFO.Prob_thres {
		if len(vars_get1) > 0 {
			for _, var_info = range vars_get1 {
				var_results <- var_info
			}
		}
		if len(vars_get2) > 0 {
			for _, var_info = range vars_get2 {
				var_results <- var_info
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
func (S *Var_Prof) FindSeedsFromPairedEnds(read_info *ReadInfo) ([]int, []int, []int, []int, []int,
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
		RAND_GEN   := rand.New(rand.NewSource(time.Now().UnixNano()))
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
		PrintLoopTraceInfo(loop_num, "FindSeedsFromPairedEnds, First:\t"+string(read_info.Read1))
		PrintLoopTraceInfo(loop_num, "FindSeedsFromPairedEnds, Second:\t"+string(read_info.Read2))

		PrintMemStats("Before FindSeeds, loop_num " + strconv.Itoa(loop_num))
		s_pos_r1_or, e_pos_r1_or, m_num_r1_or, has_seeds_r1_or =
			INDEX.FindSeeds(read_info.Read1, read_info.Rev_read1, r_pos_r1_or, m_pos_r1_or)
		PrintSeedTraceInfo("r1_or", e_pos_r1_or, s_pos_r1_or, read_info.Read1)
		if has_seeds_r1_or {
			PrintExtendTraceInfo("r1_or", read_info.Read1[e_pos_r1_or:s_pos_r1_or+1],
				e_pos_r1_or, s_pos_r1_or, m_num_r1_or, m_pos_r1_or)
		}
		s_pos_r1_rc, e_pos_r1_rc, m_num_r1_rc, has_seeds_r1_rc =
			INDEX.FindSeeds(read_info.Rev_comp_read1, read_info.Comp_read1, r_pos_r1_rc, m_pos_r1_rc)
		PrintSeedTraceInfo("r1_rc", e_pos_r1_rc, s_pos_r1_rc, read_info.Rev_comp_read1)
		if has_seeds_r1_rc {
			PrintExtendTraceInfo("r1_rc", read_info.Rev_comp_read1[e_pos_r1_rc:s_pos_r1_rc+1],
				e_pos_r1_rc, s_pos_r1_rc, m_num_r1_rc, m_pos_r1_rc)
		}
		s_pos_r2_or, e_pos_r2_or, m_num_r2_or, has_seeds_r2_or =
			INDEX.FindSeeds(read_info.Read2, read_info.Rev_read2, r_pos_r2_or, m_pos_r2_or)
		PrintSeedTraceInfo("r2_or", e_pos_r2_or, s_pos_r2_or, read_info.Read2)
		if has_seeds_r2_or {
			PrintExtendTraceInfo("r2_or", read_info.Read1[e_pos_r2_or:s_pos_r2_or+1],
				e_pos_r2_or, s_pos_r2_or, m_num_r2_or, m_pos_r2_or)
		}
		s_pos_r2_rc, e_pos_r2_rc, m_num_r2_rc, has_seeds_r2_rc =
			INDEX.FindSeeds(read_info.Rev_comp_read2, read_info.Comp_read2, r_pos_r2_rc, m_pos_r2_rc)
		PrintSeedTraceInfo("r2_rc", e_pos_r2_rc, s_pos_r2_rc, read_info.Rev_comp_read2)
		if has_seeds_r2_rc {
			PrintExtendTraceInfo("r2_rc", read_info.Rev_comp_read2[e_pos_r2_rc:s_pos_r2_rc+1],
				e_pos_r2_rc, s_pos_r2_rc, m_num_r2_rc, m_pos_r2_rc)
		}
		PrintMemStats("After FindSeeds, loop_num " + strconv.Itoa(loop_num))

		if has_seeds_r1_or && has_seeds_r2_rc {
			PrintExtendTraceInfo("r1_or(F1R2)", read_info.Read1[e_pos_r1_or:s_pos_r1_or+1],
				e_pos_r1_or, s_pos_r1_or, m_num_r1_or, m_pos_r1_or)
			PrintExtendTraceInfo("r2_rc(F1R2)", read_info.Read2[e_pos_r2_rc:s_pos_r2_rc+1],
				e_pos_r2_rc, s_pos_r2_rc, m_num_r2_rc, m_pos_r2_rc)
			for i = 0; i < m_num_r1_or; i++ {
				for j = 0; j < m_num_r2_rc; j++ {
					//Check if alignments are likely pair-end alignments
					if (m_pos_r2_rc[j]-m_pos_r1_or[i]) >= PARA_INFO.Read_len &&
						(m_pos_r2_rc[j]-m_pos_r1_or[i]) <= PARA_INFO.Read_len+PARA_INFO.Max_ins {

						PrintPairedSeedInfo("r1_or, r2_rc, paired pos", m_pos_r1_or[i], m_pos_r2_rc[j])
						s_pos_r1  = append(s_pos_r1, s_pos_r1_or)
						e_pos_r1  = append(e_pos_r1, e_pos_r1_or)
						s_pos_r2  = append(s_pos_r2, s_pos_r2_rc)
						e_pos_r2  = append(e_pos_r2, e_pos_r2_rc)
						m_pos_r1  = append(m_pos_r1, m_pos_r1_or[i])
						m_pos_r2  = append(m_pos_r2, m_pos_r2_rc[j])
						strand_r1 = append(strand_r1, true)
						strand_r2 = append(strand_r2, false)
					}
				}
			}
		}
		if has_seeds_r1_rc && has_seeds_r2_or {
			PrintExtendTraceInfo("r1_rc (F2R1)", read_info.Read1[e_pos_r1_rc:s_pos_r1_rc+1],
				e_pos_r1_rc, s_pos_r1_rc, m_num_r1_rc, m_pos_r1_rc)
			PrintExtendTraceInfo("r2_or (F2R1)", read_info.Read2[e_pos_r2_or:s_pos_r2_or+1],
				e_pos_r2_or, s_pos_r2_or, m_num_r2_or, m_pos_r2_or)
			for i = 0; i < m_num_r1_rc; i++ {
				for j = 0; j < m_num_r2_or; j++ {
					//Check if alignments are likely pair-end alignments
					if (m_pos_r1_rc[i]-m_pos_r2_or[j]) >= PARA_INFO.Read_len &&
						(m_pos_r1_rc[i]-m_pos_r2_or[j]) <= PARA_INFO.Read_len+PARA_INFO.Max_ins {

						PrintPairedSeedInfo("r1_rc, r2_or, paired pos", m_pos_r1_rc[i], m_pos_r2_or[j])
						s_pos_r1  = append(s_pos_r1, s_pos_r1_rc)
						e_pos_r1  = append(e_pos_r1, e_pos_r1_rc)
						s_pos_r2  = append(s_pos_r2, s_pos_r2_or)
						e_pos_r2  = append(e_pos_r2, e_pos_r2_or)
						m_pos_r1  = append(m_pos_r1, m_pos_r1_rc[i])
						m_pos_r2  = append(m_pos_r2, m_pos_r2_or[j])
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
			RAND_GEN   := rand.New(rand.NewSource(time.Now().UnixNano()))
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
FindVariantsFromExtension determines Vars based on alignment between reads and multi-genomes.
	Extend read and ref from exact matches found from bachward search with FM-index
	Perform backward (for left extension of read and ref) and forward alignment (for right extension)
	between read and multigenome to determine aligned bases as candidates for variant calls
--------------------------------------------------------------------------------------------------*/
func (S *Var_Prof) FindVariantsFromExtension(s_pos, e_pos, m_pos int, read, qual []byte,
	align_info *AlignInfo) ([]VarInfo, int, int, float64) {

	PrintMemStats("Before FindVariantsFromExtension, m_pos " + strconv.Itoa(m_pos))

	var i, j, del_len int
	var is_var, is_del bool

	l_read_flank_len := e_pos + PARA_INFO.Seed_backup
	l_read_flank, l_qual_flank := read[:l_read_flank_len], qual[:l_read_flank_len]

	l_ref_flank := make([]byte, 0)
	l_ref_pos_map := make([]int, 0)
	l_align_e_pos := m_pos - 1 + PARA_INFO.Seed_backup
	i = l_align_e_pos
	j = 0 //to check length of l_ref_flank
	for j < l_read_flank_len && i >= 0 {
		if _, is_var = INDEX.Var_Prof[i]; is_var {
			if del_len, is_del = INDEX.Del_Var[i]; is_del {
				if i+del_len <= l_align_e_pos && del_len < len(l_ref_flank) {
					l_ref_flank = l_ref_flank[:len(l_ref_flank)-del_len]
					l_ref_pos_map = l_ref_pos_map[:len(l_ref_pos_map)-del_len]
					j -= del_len
				}
			}
		}
		l_ref_pos_map = append(l_ref_pos_map, i)
		l_ref_flank = append(l_ref_flank, INDEX.Seq[i])
		j++
		i--
	}
	l_align_s_pos := i + 1

	//Reverse l_ref_pos_map and l_ref_flank to get them in original direction
	for i, j = 0, len(l_ref_pos_map)-1; i < j; i, j = i+1, j-1 {
		l_ref_pos_map[i], l_ref_pos_map[j] = l_ref_pos_map[j], l_ref_pos_map[i]
	}
	for i, j = 0, len(l_ref_flank)-1; i < j; i, j = i+1, j-1 {
		l_ref_flank[i], l_ref_flank[j] = l_ref_flank[j], l_ref_flank[i]
	}

	seed_len := s_pos - e_pos + 1
	r_read_flank_len := len(read) - s_pos - 1 + PARA_INFO.Seed_backup
	r_read_flank, r_qual_flank := read[len(read)-r_read_flank_len:], qual[len(read)-r_read_flank_len:]

	r_ref_flank := make([]byte, 0)
	r_ref_pos_map := make([]int, 0)
	r_align_s_pos := m_pos + seed_len - PARA_INFO.Seed_backup
	i = r_align_s_pos
	j = 0 //to check length of r_ref_flank
	for j < r_read_flank_len && i < len(INDEX.Seq) {
		r_ref_pos_map = append(r_ref_pos_map, i)
		r_ref_flank = append(r_ref_flank, INDEX.Seq[i])
		if _, is_var = INDEX.Var_Prof[i]; is_var {
			if del_len, is_del = INDEX.Del_Var[i]; is_del {
				if i+del_len < len(INDEX.Seq) {
					i += del_len
				}
			}
		}
		j++
		i++
	}

	PrintComparedReadRef(l_read_flank, l_ref_flank, r_read_flank, r_ref_flank)
	PrintRefPosMap(l_ref_pos_map, r_ref_pos_map)

	l_Ham_dist, l_Edit_dist, l_bt_mat, l_m, l_n, l_var_pos, l_var_base, l_var_qual, l_var_type :=
		S.BackwardDistance(l_read_flank, l_qual_flank, l_ref_flank, l_align_s_pos, align_info.Bw_Dist_D,
			align_info.Bw_Dist_IS, align_info.Bw_Dist_IT, align_info.Bw_Trace_D, align_info.Bw_Trace_IS, align_info.Bw_Trace_IT, l_ref_pos_map)
	r_Ham_dist, r_Edit_dist, r_bt_mat, r_m, r_n, r_var_pos, r_var_base, r_var_qual, r_var_type :=
		S.ForwardDistance(r_read_flank, r_qual_flank, r_ref_flank, r_align_s_pos, align_info.Fw_Dist_D,
			align_info.Fw_Dist_IS, align_info.Fw_Dist_IT, align_info.Fw_Trace_D, align_info.Fw_Trace_IS, align_info.Fw_Trace_IT, r_ref_pos_map)

	var vars_arr []VarInfo
	prob := l_Ham_dist + r_Ham_dist + l_Edit_dist + r_Edit_dist
	if prob <= PARA_INFO.Prob_thres {
		if l_m > 0 && l_n > 0 {
			l_pos, l_base, l_qual, l_type := S.BackwardTraceBack(l_read_flank, l_qual_flank, l_ref_flank, l_m, l_n, l_align_s_pos,
				l_bt_mat, align_info.Bw_Trace_D, align_info.Bw_Trace_IS, align_info.Bw_Trace_IT, l_ref_pos_map)
			l_var_pos  = append(l_var_pos, l_pos...)
			l_var_base = append(l_var_base, l_base...)
			l_var_qual = append(l_var_qual, l_qual...)
			l_var_type = append(l_var_type, l_type...)
		}
		PrintMatchTraceInfo(m_pos, l_align_s_pos, prob, l_var_pos, read)
		if r_m > 0 && r_n > 0 {
			r_pos, r_base, r_qual, r_type := S.ForwardTraceBack(r_read_flank, r_qual_flank, r_ref_flank, r_m, r_n, r_align_s_pos,
				r_bt_mat, align_info.Fw_Trace_D, align_info.Fw_Trace_IS, align_info.Fw_Trace_IT, r_ref_pos_map)
			r_var_pos  = append(r_var_pos, r_pos...)
			r_var_base = append(r_var_base, r_base...)
			r_var_qual = append(r_var_qual, r_qual...)
			r_var_type = append(r_var_type, r_type...)
		}
		PrintMatchTraceInfo(m_pos, r_align_s_pos, prob, r_var_pos, read)
		var k int
		for k = 0; k < len(l_var_pos); k++ {
			PrintMemStats("Before GetVar left, var_num " + strconv.Itoa(k))
			var var_info VarInfo
			var_info.Pos, var_info.Bases, var_info.BaseQ, var_info.Type = uint32(l_var_pos[k]), l_var_base[k], l_var_qual[k], l_var_type[k]
			vars_arr = append(vars_arr, var_info)
			PrintMemStats("After GetVar left, var_num " + strconv.Itoa(k))
		}
		for k = 0; k < len(r_var_pos); k++ {
			PrintMemStats("Before GetVar right, var_num " + strconv.Itoa(k))
			var var_info VarInfo
			var_info.Pos, var_info.Bases, var_info.BaseQ, var_info.Type = uint32(r_var_pos[k]), r_var_base[k], r_var_qual[k], r_var_type[k]
			vars_arr = append(vars_arr, var_info)
			PrintMemStats("After GetVar right, var_num " + strconv.Itoa(k))
		}
		PrintMemStats("After FindVariantsFromExtension, m_pos " + strconv.Itoa(m_pos))
		return vars_arr, l_align_s_pos, r_align_s_pos, prob
	}
	PrintMemStats("After FindVariantsFromExtension, m_pos " + strconv.Itoa(m_pos))
	return vars_arr, -1, -1, -1
}

/*--------------------------------------------------------------------------------------------------
UpdateSNPProb updates variant probablilities for all possible variants.
	Input: a variant of type VarInfo which is a SNP.
	Output: update of all SNPs in S.Var_Prob[var_info.Pos].
--------------------------------------------------------------------------------------------------*/
func (S *Var_Prof) UpdateSNPProb(var_info VarInfo) {
	pos := var_info.Pos
	a := string(var_info.Bases[0])
	q := var_info.BaseQ[0]

	if _, var_exist := S.Var_Prob[pos]; !var_exist {
		S.Var_Prob[pos] = make(map[string]float64)
		S.Var_Prob[pos][string(INDEX.Seq[int(pos)])] = 1 - 3*NEW_SNP_RATE
		for _, b := range STD_BASES {
			if _, ok := S.Var_Prob[pos][string(b)]; !ok {
				S.Var_Prob[pos][string(b)] = NEW_SNP_RATE
			}
		}
		S.Var_BaseQ[pos]  = make(map[string][][]byte)
		S.Var_Type[pos]   = make(map[string][]int)
		S.Var_RNum[pos]   = make(map[string]int)
		S.Chr_Dis[pos]    = make(map[string][]int)
		S.Chr_Diff[pos]   = make(map[string][]int)
		S.Aln_Prob[pos]   = make(map[string][]float64)
		S.Chr_Prob[pos]   = make(map[string][]float64)
		S.Read_Info[pos]  = make(map[string][][]byte)
		S.Start_Pos1[pos] = make(map[string][]int)
		S.Start_Pos2[pos] = make(map[string][]int)
		S.Strand1[pos]    = make(map[string][]bool)
		S.Strand2[pos]    = make(map[string][]bool)
	}
	S.Var_BaseQ[pos][a]  = append(S.Var_BaseQ[pos][a], var_info.BaseQ)
	S.Var_Type[pos][a]   = append(S.Var_Type[pos][a], var_info.Type)
	S.Var_RNum[pos][a]   += 1
	S.Chr_Dis[pos][a]    = append(S.Chr_Dis[pos][a], var_info.CDis)
	S.Chr_Diff[pos][a]   = append(S.Chr_Diff[pos][a], var_info.CDiff)
	S.Aln_Prob[pos][a]   = append(S.Aln_Prob[pos][a], var_info.AProb)
	S.Chr_Prob[pos][a]   = append(S.Chr_Prob[pos][a], var_info.CProb)
	S.Read_Info[pos][a]  = append(S.Read_Info[pos][a], var_info.RInfo)
	S.Start_Pos1[pos][a] = append(S.Start_Pos1[pos][a], var_info.SPos1)
	S.Start_Pos2[pos][a] = append(S.Start_Pos2[pos][a], var_info.SPos2)
	S.Strand1[pos][a]    = append(S.Strand1[pos][a], var_info.Stra1)
	S.Strand2[pos][a]    = append(S.Strand2[pos][a], var_info.Stra2)

	var p float64
	p_ab := make(map[string]float64)
	p_a := 0.0

	for b, p_b := range S.Var_Prob[pos] {
		if a == b {
			p = 1.0 - math.Pow(10, -(float64(q)-33)/10.0) //Phred-encoding factor (33) need to be estimated from input data
		} else {
			p = math.Pow(10, -(float64(q)-33)/10.0) / 3 //need to be refined, e.g., checked with diff cases (snp vs. indel)
		}
		p_ab[b] = p
		p_a += p_b * p_ab[b]
	}
	for b, p_b := range S.Var_Prob[pos] {
		S.Var_Prob[pos][b] = p_b * (p_ab[b] / p_a)
	}
}

/*--------------------------------------------------------------------------------------------------
UpdateIndelProb updates Indel probablilities for all possible Indels.
	Input: a variant of type VarInfo which is an Indel.
	Output: update of all Indels in S.Var_Prob[var_info.Pos].
--------------------------------------------------------------------------------------------------*/
func (S *Var_Prof) UpdateIndelProb(var_info VarInfo) {
	pos := var_info.Pos
	a := string(var_info.Bases)
	q := var_info.BaseQ

	//Notice: Need to be corrected!
	if _, var_exist := S.Var_Prob[pos]; !var_exist {
		S.Var_Prob[pos] = make(map[string]float64)
		S.Var_Prob[pos][string(INDEX.Seq[int(pos)])] = 1 - 4*NEW_SNP_RATE
		for _, b := range STD_BASES {
			if _, ok := S.Var_Prob[pos][string(b)]; !ok {
				S.Var_Prob[pos][string(b)] = NEW_SNP_RATE
			}
		}
		S.Var_BaseQ[pos]  = make(map[string][][]byte)
		S.Var_Type[pos]   = make(map[string][]int)
		S.Var_RNum[pos]   = make(map[string]int)
		S.Chr_Dis[pos]    = make(map[string][]int)
		S.Chr_Diff[pos]   = make(map[string][]int)
		S.Aln_Prob[pos]   = make(map[string][]float64)
		S.Chr_Prob[pos]   = make(map[string][]float64)
		S.Read_Info[pos]  = make(map[string][][]byte)
		S.Start_Pos1[pos] = make(map[string][]int)
		S.Start_Pos2[pos] = make(map[string][]int)
		S.Strand1[pos]    = make(map[string][]bool)
		S.Strand2[pos]    = make(map[string][]bool)
	}
	//Notice: Using NEW_SNP_RATE, need to consider NEW_INDEL_RATE instead
	if _, ok := S.Var_Prob[pos][a]; !ok {
		S.Var_Prob[pos][a] = NEW_SNP_RATE
	}
	S.Var_BaseQ[pos][a]  = append(S.Var_BaseQ[pos][a], var_info.BaseQ)
	S.Var_Type[pos][a]   = append(S.Var_Type[pos][a], var_info.Type)
	S.Var_RNum[pos][a]   += 1
	S.Chr_Dis[pos][a]    = append(S.Chr_Dis[pos][a], var_info.CDis)
	S.Chr_Diff[pos][a]   = append(S.Chr_Diff[pos][a], var_info.CDiff)
	S.Aln_Prob[pos][a]   = append(S.Aln_Prob[pos][a], var_info.AProb)
	S.Chr_Prob[pos][a]   = append(S.Chr_Prob[pos][a], var_info.CProb)
	S.Read_Info[pos][a]  = append(S.Read_Info[pos][a], var_info.RInfo)
	S.Start_Pos1[pos][a] = append(S.Start_Pos1[pos][a], var_info.SPos1)
	S.Start_Pos2[pos][a] = append(S.Start_Pos2[pos][a], var_info.SPos2)
	S.Strand1[pos][a]    = append(S.Strand1[pos][a], var_info.Stra1)
	S.Strand2[pos][a]    = append(S.Strand2[pos][a], var_info.Stra2)

	var p float64
	var qi byte
	p_ab := make(map[string]float64)
	p_a := 0.0

	for b, p_b := range S.Var_Prob[pos] {
		p = 1
		if a == b {
			for _, qi = range q[0:1] { //temporary change
				p *= (1.0 - math.Pow(10, -(float64(qi)-33)/10.0)) //Phred-encoding factor (33) need to be estimated from input data
			}
		} else {
			for _, qi = range q[0:1] { //temporary change
				p *= (math.Pow(10, -(float64(qi)-33)/10.0) / 3) //need to be refined, e.g., checked with diff cases (snp vs. indel)
			}
		}
		p_ab[b] = p
		p_a += p_b * p_ab[b]
	}
	for b, p_b := range S.Var_Prob[pos] {
		S.Var_Prob[pos][b] = p_b * (p_ab[b] / p_a)
	}
}

/*--------------------------------------------------------------------------------------------------
OutputVarCalls determines variant calls, convert their probabilities to Phred scores, and writes them
to file in proper format (VCF-like format at this stage).
--------------------------------------------------------------------------------------------------*/
func (S *Var_Prof) OutputVarCalls() {

	file, err := os.Create(INPUT_INFO.Var_call_file)
	if err != nil {
		fmt.Println("Error: Create output file", err)
		os.Exit(1)
	}
	defer file.Close()

	fmt.Println("Outputing Variant Calls...")
	file.WriteString("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tVAR_PROB\tBASE_NUM\tBASE_QUAL\tCHR_DIS\tCHR_DIFF\tALN_PROB\tPAIR_PROB\tS_POS1\tBRANCH1\tS_POS2\tBRANCH2\tREAD_HEADER\tALN_BASE\tBASE_NUM\t\n")
	var var_pos uint32
	Var_Pos := make([]int, 0, len(S.Var_Prob))
	for var_pos, _ = range S.Var_Prob {
		Var_Pos = append(Var_Pos, int(var_pos))
	}
	sort.Ints(Var_Pos)

	var var_base, var_call, str_qual string
	var str_a, str_b string
	var line_a, line_b []string
	var var_call_prob, var_prob float64
	var var_num, idx int
	var is_var bool
	for _, pos := range Var_Pos {
		var_pos = uint32(pos)
		//Get variant call by considering maximum prob
		var_call_prob = 0
		for var_base, var_prob = range S.Var_Prob[var_pos] {
			if var_call_prob < var_prob {
				var_call_prob = var_prob
				var_call = var_base
			}
		}
		//Start getting variant call info
		line_a = make([]string, 0)
		//#CHROM
		line_a = append(line_a, ".")
		//POS
		line_a = append(line_a, strconv.Itoa(pos+1))
		//ID
		line_a = append(line_a, ".")
		//REF & ALT
		if _, is_var = INDEX.Var_Prof[pos]; is_var {
			line_a = append(line_a, string(INDEX.Var_Prof[pos][0]))
			line_a = append(line_a, var_call)
		} else {
			if len(S.Var_Type[var_pos][var_call]) != 0 {
				if S.Var_Type[var_pos][var_call][0] == 2 { //DEL
					//Ignore indel calls of length 2 that are homopolymer
					if len(var_call) == 2 && var_call[0] == var_call[1] {
						continue
					}
					line_a = append(line_a, var_call)
					line_a = append(line_a, string(INDEX.Seq[pos]))
				} else if S.Var_Type[var_pos][var_call][0] == 1 { //INS
					//Ignore indel calls of length 2 that are homopolymer
					if len(var_call) == 2 && var_call[0] == var_call[1] {
						continue
					}
					line_a = append(line_a, string(INDEX.Seq[pos]))
					line_a = append(line_a, var_call)
				} else { //SUB
					//Ignore variants that are identical with ref
					if var_call == string(INDEX.Seq[pos]) {
						continue
					}
					line_a = append(line_a, string(INDEX.Seq[pos]))
					line_a = append(line_a, var_call)
				}
			} else {
				fmt.Println("Var calls not from alignment", pos, var_call, string(INDEX.Seq[pos]))
				continue
			}
		}
		//QUAL
		str_qual = strconv.FormatFloat(-10*math.Log10(1-var_call_prob), 'f', 5, 32)
		if str_qual != "+Inf" {
			line_a = append(line_a, str_qual)
		} else {
			line_a = append(line_a, "1000")
		}
		//FILTER
		line_a = append(line_a, ".")
		//INFO
		line_a = append(line_a, ".")
		//FORMAT
		line_a = append(line_a, ".")
		//ISC-INFO
		line_a = append(line_a, strconv.FormatFloat(var_call_prob, 'f', 5, 32))
		line_a = append(line_a, strconv.Itoa(S.Var_RNum[var_pos][var_call]))
		str_a = strings.Join(line_a, "\t")
		line_b = make([]string, 0)
		for var_base, var_num = range S.Var_RNum[var_pos] {
			line_b = append(line_b, var_base)
			line_b = append(line_b, strconv.Itoa(var_num))
		}
		str_b = strings.Join(line_b, "\t")
		//Write variant calls to file
		for idx, _ = range S.Var_BaseQ[var_pos][var_call] {
			_, err = file.WriteString(str_a + "\t")
			_, err = file.WriteString(string(S.Var_BaseQ[var_pos][var_call][idx]) + "\t")
			_, err = file.WriteString(strconv.Itoa(S.Chr_Dis[var_pos][var_call][idx]) + "\t")
			_, err = file.WriteString(strconv.Itoa(S.Chr_Diff[var_pos][var_call][idx]) + "\t")
			_, err = file.WriteString(strconv.FormatFloat(S.Aln_Prob[var_pos][var_call][idx], 'f', 20, 64) + "\t")
			_, err = file.WriteString(strconv.FormatFloat(S.Chr_Prob[var_pos][var_call][idx], 'f', 20, 64) + "\t")
			_, err = file.WriteString(strconv.Itoa(S.Start_Pos1[var_pos][var_call][idx]) + "\t")
			_, err = file.WriteString(strconv.FormatBool(S.Strand1[var_pos][var_call][idx]) + "\t")
			_, err = file.WriteString(strconv.Itoa(S.Start_Pos2[var_pos][var_call][idx]) + "\t")
			_, err = file.WriteString(strconv.FormatBool(S.Strand2[var_pos][var_call][idx]) + "\t")
			_, err = file.WriteString(string(S.Read_Info[var_pos][var_call][idx]) + "\t")
			_, err = file.WriteString(str_b + "\n")
			if err != nil {
				fmt.Println("Error: Write variant calls to file", err)
				break
			}
		}
	}
}
