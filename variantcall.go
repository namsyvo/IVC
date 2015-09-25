//---------------------------------------------------------------------------------------------------
// IVC: variantcall.go - Calling genomic variants based on alignment between reads and the reference multigenome.
// Variants and probability of correct variant calls is determined using Bayesian update.
// Copyright 2015 Nam Sy Vo.
//---------------------------------------------------------------------------------------------------

package ivc

import (
	"bufio"
	"bytes"
	"log"
	"math"
	"math/rand"
	"os"
	"runtime"
	"runtime/pprof"
	"sort"
	"strconv"
	"strings"
	"sync"
	"time"
)

var (
	Q2C map[byte]float64 //Pre-calculated alignment cost which is based on Phred-based quality
	Q2E map[byte]float64 //Pre-calculated error rate which is corresponding to Phred-based quality
	Q2P map[byte]float64 //Pre-calculated probability which is corresponding to Phred-based quality
	L2E []float64        //Pre-calculated indel error which is corresponding to length of indels
)

//---------------------------------------------------------------------------------------------------
// VarInfo represents variants obtained during alignment phase.
// It serves as temporary variable during variant calling phase.
//---------------------------------------------------------------------------------------------------
type VarInfo struct {
	Pos     uint32  //Postion of variant (on the reference)
	Bases   []byte  //Aligned bases to be the variant
	BQual   []byte  //Quality sequences (in FASTQ format) of bases to be the variant
	Type    int     //Type of the variant (0: sub, 1: ins, 2: del; other types will be considered in future)
	CDis    int     //Chromosomal distance between alignment positions of two read-ends
	CDiff   int     //Chromosomal distance between aligned pos and true pos
	MProb   float64 //Probability of mapping read corectly (mapping quality)
	AProb   float64 //Probability of aligning read correctly (alignment quality)
	IProb   float64 //Probability of insert size to be correct (for pair-end reads)
	SPos1   int     //Starting position on read1 of exact match (or ending position from backward search with FM-index)
	SPos2   int     //Starting position on read2 of exact match (or ending position from backward search with FM-index)
	Strand1 bool    //Strand (backward/forward) of read1 of exact match
	Strand2 bool    //Strand (backward/forward) of read2 of exact match
	RInfo   []byte  //Information sequences (in FASTQ format) of aligned reads (header of reads in FASTQ format)
}

//---------------------------------------------------------------------------------------------------
// VarCall represents info of aligned reads and bases at the variant call positions on reference multigenome.
// This struct also has functions defined on it for calling variants.
//---------------------------------------------------------------------------------------------------
type VarCall struct {
	//VarProb stores all possible variants at each position and their confident probablilities.
	//Prior probablities will be obtained from reference genomes and variant profiles.
	//Posterior probabilities will be updated during alignment phase based on incomming aligned bases
	VarProb   map[uint32]map[string]float64   //Probability of the variant call
	VarType   map[uint32]map[string]int       //Type of variants (0: sub, 1: ins, 2: del; other types will be considered in future)
	VarRNum   map[uint32]map[string]int       //Numer of reads (bases) aligned to the variant call postion
	ChrDis    map[uint32]map[string][]int     //Chromosomal distance between two aligned ends
	ChrDiff   map[uint32]map[string][]int     //Chromosomal distance betwwen the aligned postion and true postion (for simulated data)
	MapProb   map[uint32]map[string][]float64 //Probability of mapping read corectly (mapping quality)
	AlnProb   map[uint32]map[string][]float64 //Probability of aligning read correctly (alignment quality)
	ChrProb   map[uint32]map[string][]float64 //Probability of insert size to be correct (for pair-end reads)
	StartPos1 map[uint32]map[string][]int     //Start position of alignment of the first end
	StartPos2 map[uint32]map[string][]int     //Start position of alignment of the second end
	Strand1   map[uint32]map[string][]bool    //Strand indicator of the first end ("true" if read has same strand with ref, "false" otherwise)
	Strand2   map[uint32]map[string][]bool    //Strand indicator of the second end ("true" if read has same strand with ref, "false" otherwise)
	VarBQual  map[uint32]map[string][][]byte  //Quality sequences (in FASTQ format) of aligned bases at the variant call position
	ReadInfo  map[uint32]map[string][][]byte  //Information sequences (in FASTQ format) of aligned reads (header of reads in FASTQ format)
}

//---------------------------------------------------------------------------------------------------
// NewVariantCaller creates an instance of VarCall and sets up its variables.
// This function will be called from the main program.
//---------------------------------------------------------------------------------------------------
func NewVariantCaller(input_info *InputInfo) *VarCall {

	//Initialize global variables
	INPUT_INFO = input_info
	if _, e := os.Stat(INPUT_INFO.Read_file_1); e != nil {
		log.Printf("Error: Read_file_1 does not exists! (err: %s)", e)
		os.Exit(1)
	}
	if _, e := os.Stat(INPUT_INFO.Read_file_2); e != nil {
		log.Printf("Error: Read_file_2 does not exists! (err: %s)", e)
		os.Exit(1)
	}
	log.Printf("Initializing indexes and parameters...")
	start_time := time.Now()
	//SetPara: 100 is maximum length of reads, 500 is maximum length of info line of reads,
	//700 is maximum insert size of paired-end simulated reads, 0.0015 is maximum sequencing error rate
	//0.01 is mutation rate (currently is estimated from dbSNP of human genome)
	PARA_INFO = SetPara(100, 500, 700, 0.0015, 0.01, INPUT_INFO.Dist_thres, INPUT_INFO.Prob_thres, INPUT_INFO.Iter_num)
	runtime.GOMAXPROCS(INPUT_INFO.Proc_num)

	//Initialize Index object for finding seeds
	INDEX = NewIndex()

	//"Local-global" variable which shared by all downstream functions from CallVariants function.
	//May consider a better solution later.
	Q2C = make(map[byte]float64)
	Q2E = make(map[byte]float64)
	Q2P = make(map[byte]float64)
	L2E = make([]float64, 50)
	var q byte
	for i := 33; i < 74; i++ {
		q = byte(i)
		//Phred-encoding factor (33) need to be estimated from input data
		Q2C[q] = -math.Log10(1.0 - math.Pow(10, -(float64(q)-33)/10.0))
		Q2E[q] = math.Pow(10, -(float64(q)-33)/10.0) / 3.0
		Q2P[q] = 1.0 - math.Pow(10, -(float64(q)-33)/10.0)
	}
	for i := 1; i < 50; i++ {
		L2E[i] = math.Pow(INDEL_ERR, float64(i))
	}

	//Initialize VarCall object for calling variants
	VC := new(VarCall)
	VC.VarProb = make(map[uint32]map[string]float64)
	VC.VarType = make(map[uint32]map[string]int)
	VC.VarRNum = make(map[uint32]map[string]int)
	if INPUT_INFO.Debug_mode {
		VC.ChrDis = make(map[uint32]map[string][]int)
		VC.ChrDiff = make(map[uint32]map[string][]int)
		VC.MapProb = make(map[uint32]map[string][]float64)
		VC.AlnProb = make(map[uint32]map[string][]float64)
		VC.ChrProb = make(map[uint32]map[string][]float64)
		VC.StartPos1 = make(map[uint32]map[string][]int)
		VC.StartPos2 = make(map[uint32]map[string][]int)
		VC.Strand1 = make(map[uint32]map[string][]bool)
		VC.Strand2 = make(map[uint32]map[string][]bool)
		VC.VarBQual = make(map[uint32]map[string][][]byte)
		VC.ReadInfo = make(map[uint32]map[string][][]byte)
	}
	var pos uint32
	var var_bases []byte
	var i int
	for var_pos, var_prof := range INDEX.VarProf {
		pos = uint32(var_pos)
		VC.VarProb[pos] = make(map[string]float64)
		//At this point, assume that all variants are biallelic
		if len(var_prof[0]) == 1 && len(var_prof[1]) == 1 {
			for i, var_bases = range var_prof {
				VC.VarProb[pos][string(var_bases)] = float64(INDEX.VarAF[var_pos][i]) - NEW_SNP_RATE
				if VC.VarProb[pos][string(var_bases)] < NEW_SNP_RATE {
					VC.VarProb[pos][string(var_bases)] = NEW_SNP_RATE
				}
			}
		} else {
			for i, var_bases = range var_prof {
				VC.VarProb[pos][string(var_bases)] = float64(INDEX.VarAF[var_pos][i]) - 1.5*NEW_SNP_RATE
				if VC.VarProb[pos][string(var_bases)] < NEW_SNP_RATE {
					VC.VarProb[pos][string(var_bases)] = NEW_SNP_RATE
				}
			}
		}
		for _, b := range STD_BASES {
			if _, ok := VC.VarProb[pos][string(b)]; !ok {
				VC.VarProb[pos][string(b)] = NEW_SNP_RATE
			}
		}
		VC.VarType[pos] = make(map[string]int)
		VC.VarRNum[pos] = make(map[string]int)
		if INPUT_INFO.Debug_mode {
			VC.ChrDis[pos] = make(map[string][]int)
			VC.ChrDiff[pos] = make(map[string][]int)
			VC.MapProb[pos] = make(map[string][]float64)
			VC.AlnProb[pos] = make(map[string][]float64)
			VC.ChrProb[pos] = make(map[string][]float64)
			VC.StartPos1[pos] = make(map[string][]int)
			VC.StartPos2[pos] = make(map[string][]int)
			VC.Strand1[pos] = make(map[string][]bool)
			VC.Strand2[pos] = make(map[string][]bool)
			VC.VarBQual[pos] = make(map[string][][]byte)
			VC.ReadInfo[pos] = make(map[string][][]byte)
		}
	}
	PrintProcessMem("Memstats after initializing the variant caller")
	index_time := time.Since(start_time)
	log.Printf("Time for initializing the variant caller:\t%s", index_time)
	log.Printf("Finish initializing indexes and parameters.")
	return VC
}

//---------------------------------------------------------------------------------------------------
// CallVariants initializes share variables, channels, reads input reads, finds all possible variants,
// and updates variant information in VarCall.
// This function will be called from main program.
//---------------------------------------------------------------------------------------------------
func (VC *VarCall) CallVariants() {
	log.Printf("Calling variants...")
	start_time := time.Now()
	//The channel read_signal is used for signaling between goroutines which run ReadReads and FindVariants,
	//when a FindSNPs goroutine finish copying a read to its own memory,
	//it signals ReadReads goroutine to scan next reads.
	read_signal := make(chan bool)

	//Call a goroutine to read input reads
	read_data := make(chan *ReadInfo, INPUT_INFO.Proc_num)
	go VC.ReadReads(read_data, read_signal)

	//Call goroutines to find Vars, pass shared variable to each goroutine
	var_results := make(chan *VarInfo)
	var wg sync.WaitGroup
	for i := 0; i < INPUT_INFO.Proc_num; i++ {
		wg.Add(1)
		go VC.FindVariants(read_data, read_signal, var_results, &wg)
	}
	go GetNoAlignReadInfo()

	go func() {
		wg.Wait()
		close(var_results)
		close(UNALIGN_INFO_CHAN)
	}()

	//Collect variants from results channel and update variant probabilities
	var var_info *VarInfo
	for var_info = range var_results {
		VC.UpdateVariantProb(var_info)
	}
	ProcessNoAlignReadInfo()
	PrintProcessMem("Memstats after calling variants")
	call_var_time := time.Since(start_time)
	log.Printf("Time for calling variants:\t%s", call_var_time)
	log.Printf("Finish calling variants.")
}

//---------------------------------------------------------------------------------------------------
// ReadReads reads all reads from input FASTQ files and put them into data channel.
//---------------------------------------------------------------------------------------------------
func (VC *VarCall) ReadReads(read_data chan *ReadInfo, read_signal chan bool) {

	fn1, fn2 := INPUT_INFO.Read_file_1, INPUT_INFO.Read_file_2
	f1, e1 := os.Open(fn1)
	if e1 != nil {
		log.Printf("Error: Open read_file_1 %s, (err: %s)", fn1, e1)
		os.Exit(1)
	}
	defer f1.Close()
	f2, e2 := os.Open(fn2)
	if e2 != nil {
		log.Printf("Error: Open read_file_2 %s, (err: %s)", fn1, e2)
		os.Exit(1)
	}
	defer f2.Close()

	read_num := 0
	scanner1 := bufio.NewScanner(f1)
	scanner2 := bufio.NewScanner(f2)
	read_info := InitReadInfo(PARA_INFO.Read_len, PARA_INFO.Info_len)
	for scanner1.Scan() && scanner2.Scan() {
		read_info.Info1 = read_info.Info1[:len(scanner1.Bytes())]
		read_info.Info2 = read_info.Info2[:len(scanner2.Bytes())]
		copy(read_info.Info1, scanner1.Bytes()) //use 1st line in input FASTQ file 1
		copy(read_info.Info2, scanner2.Bytes()) //use 1st line in input FASTQ file 2
		scanner1.Scan()
		scanner2.Scan()
		read_info.Read1 = read_info.Read1[:len(scanner1.Bytes())]
		read_info.Read2 = read_info.Read2[:len(scanner2.Bytes())]
		copy(read_info.Read1, scanner1.Bytes()) //use 2nd line in input FASTQ file 1
		copy(read_info.Read2, scanner2.Bytes()) //use 2nd line in input FASTQ file 2
		scanner1.Scan()                         //ignore 3rd line in 1st input FASTQ file 1
		scanner2.Scan()                         //ignore 3rd line in 2nd input FASTQ file 2
		scanner1.Scan()
		scanner2.Scan()
		read_info.Qual1 = read_info.Qual1[:len(scanner1.Bytes())]
		read_info.Qual2 = read_info.Qual2[:len(scanner2.Bytes())]
		copy(read_info.Qual1, scanner1.Bytes()) //use 4th line in input FASTQ file 1
		copy(read_info.Qual2, scanner2.Bytes()) //use 4th line in input FASTQ file 2
		if len(read_info.Read1) > 0 && len(read_info.Read2) > 0 {
			read_num++
			read_data <- read_info
			read_signal <- true
		}
		if read_num%10000 == 0 {
			PrintProcessMem("Memstats after distributing " + strconv.Itoa(read_num) + " reads")
			if INPUT_INFO.Debug_mode {
				pprof.WriteHeapProfile(MEM_FILE)
			}
		}
	}
	close(read_data)
}

//---------------------------------------------------------------------------------------------------
// FindVariants takes data from data channel, find all possible Vars and put them into results channel.
//---------------------------------------------------------------------------------------------------
func (VC *VarCall) FindVariants(read_data chan *ReadInfo, read_signal chan bool, var_results chan *VarInfo,
	wg *sync.WaitGroup) {

	defer wg.Done()

	//Initialize inter-function share variables
	read_info := InitReadInfo(PARA_INFO.Read_len, PARA_INFO.Info_len)
	edit_aln_info := InitEditAlnInfo(2 * PARA_INFO.Read_len)
	seed_pos := make([][]int, 4)
	for i := 0; i < 4; i++ {
		seed_pos[i] = make([]int, INPUT_INFO.Max_snum)
	}
	rand_gen := rand.New(rand.NewSource(time.Now().UnixNano()))
	for read := range read_data {
		read_info.Info1 = read_info.Info1[:len(read.Info1)]
		read_info.Info2 = read_info.Info2[:len(read.Info2)]
		copy(read_info.Info1, read.Info1)
		copy(read_info.Info2, read.Info2)
		read_info.Read1 = read_info.Read1[:len(read.Read1)]
		read_info.Read2 = read_info.Read2[:len(read.Read2)]
		copy(read_info.Read1, read.Read1)
		copy(read_info.Read2, read.Read2)
		read_info.Qual1 = read_info.Qual1[:len(read.Qual1)]
		read_info.Qual2 = read_info.Qual2[:len(read.Qual2)]
		copy(read_info.Qual1, read.Qual1)
		copy(read_info.Qual2, read.Qual2)
		<-read_signal

		RevComp(read_info.Read1, read_info.Qual1, read_info.Rev_read1, read_info.Rev_comp_read1,
			read_info.Comp_read1, read_info.Rev_qual1)
		RevComp(read_info.Read2, read_info.Qual2, read_info.Rev_read2, read_info.Rev_comp_read2,
			read_info.Comp_read2, read_info.Rev_qual2)

		VC.FindVariantsPE(read_info, edit_aln_info, seed_pos, rand_gen, var_results)
	}
}

//---------------------------------------------------------------------------------------------------
// FindVariantsPE finds variants from alignment between pair-end reads and the multigenome.
// It uses seed-and-extend strategy and finds best alignment candidates through several iterations.
//---------------------------------------------------------------------------------------------------
func (VC *VarCall) FindVariantsPE(read_info *ReadInfo, edit_aln_info *EditAlnInfo, seed_pos [][]int,
	rand_gen *rand.Rand, var_results chan *VarInfo) {
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
	var vars1, vars2, vars_get1, vars_get2 []*VarInfo
	var l_align_pos1, l_align_pos2 int
	var seed_info1, seed_info2 *SeedInfo
	var has_seeds bool
	var align_prob1, align_prob2 float64
	var cand_num []int
	var p_idx, s_idx, c_num int

	paired_prob := math.MaxFloat64
	loop_has_cand := 0
	for loop_num := 1; loop_num <= PARA_INFO.Iter_num; loop_num++ {
		if loop_num == 1 {
			//seed_info1, seed_info2, has_seeds = INDEX.FindSeedsPE(read_info, seed_pos, rand_gen, 0) //Search from benginning
			seed_info1, seed_info2, has_seeds = INDEX.FindSeedsPE(read_info, seed_pos, rand_gen, 1) //Random search
		} else {
			seed_info1, seed_info2, has_seeds = INDEX.FindSeedsPE(read_info, seed_pos, rand_gen, 1) //Random search
		}
		if !has_seeds {
			cand_num = append(cand_num, 0)
			continue
		}
		c_num = 0
		for p_idx = 0; p_idx < len(seed_info1.s_pos); p_idx++ {
			//For conventional paired-end sequencing (i.e. Illumina) the directions should be F-R
			//For other kinds of variants (e.g inversions) or other technologies, they can be F-F or R-R
			//For mate-pair, they can be R-F (need to be confirmed)
			if seed_info1.strand[p_idx] == seed_info2.strand[p_idx] {
				continue
			}
			//Find variants for the first end
			if seed_info1.strand[p_idx] == true {
				vars1, _, _, align_prob1 = VC.ExtendSeeds(seed_info1.s_pos[p_idx], seed_info1.e_pos[p_idx],
					seed_info1.m_pos[p_idx], read_info.Read1, read_info.Qual1, edit_aln_info)
			} else {
				vars1, _, _, align_prob1 = VC.ExtendSeeds(seed_info1.s_pos[p_idx], seed_info1.e_pos[p_idx],
					seed_info1.m_pos[p_idx], read_info.Rev_comp_read1, read_info.Rev_qual1, edit_aln_info)
			}
			//Find variants for the second end
			if seed_info2.strand[p_idx] == true {
				vars2, _, _, align_prob2 = VC.ExtendSeeds(seed_info2.s_pos[p_idx], seed_info2.e_pos[p_idx],
					seed_info2.m_pos[p_idx], read_info.Read2, read_info.Qual2, edit_aln_info)
			} else {
				vars2, _, _, align_prob2 = VC.ExtendSeeds(seed_info2.s_pos[p_idx], seed_info2.e_pos[p_idx],
					seed_info2.m_pos[p_idx], read_info.Rev_comp_read2, read_info.Rev_qual2, edit_aln_info)
			}
			//Currently, variants can be called iff both read-ends can be aligned
			if align_prob1 != -1 && align_prob2 != -1 {
				c_num++
				ins_prob := -math.Log10(math.Exp(-math.Pow(math.Abs(float64(l_align_pos1-l_align_pos2))-400.0, 2.0) / (2 * 50 * 50)))
				if paired_prob > align_prob1+align_prob2 {
					paired_prob = align_prob1 + align_prob2
					PrintGetVariants("Find_min", paired_prob, align_prob1, align_prob2, vars1, vars2)
					vars_get1 = make([]*VarInfo, len(vars1)) //need to reset vars_get1 here
					vars_get2 = make([]*VarInfo, len(vars2)) //need to reset vars_get2 here
					loop_has_cand = loop_num
					for s_idx = 0; s_idx < len(vars1); s_idx++ {
						vars_get1[s_idx] = vars1[s_idx]
						if INPUT_INFO.Debug_mode {
							//Update vars_get1 with other info
							vars_get1[s_idx].CDis = l_align_pos1 - l_align_pos2
							vars_get1[s_idx].CDiff = l_align_pos1 - int(true_pos1)
							vars_get1[s_idx].AProb = align_prob1
							vars_get1[s_idx].IProb = ins_prob
							vars_get1[s_idx].SPos1 = seed_info1.e_pos[p_idx]
							vars_get1[s_idx].SPos2 = seed_info2.e_pos[p_idx]
							vars_get1[s_idx].Strand1 = seed_info1.strand[p_idx]
							vars_get1[s_idx].Strand2 = seed_info2.strand[p_idx]
							vars_get1[s_idx].RInfo = read_info1
						}
					}
					for s_idx = 0; s_idx < len(vars2); s_idx++ {
						vars_get2[s_idx] = vars2[s_idx]
						if INPUT_INFO.Debug_mode {
							//Update vars_get2 with other info
							vars_get2[s_idx].CDis = l_align_pos1 - l_align_pos2
							vars_get2[s_idx].CDiff = l_align_pos2 - int(true_pos2)
							vars_get2[s_idx].AProb = align_prob2
							vars_get2[s_idx].IProb = ins_prob
							vars_get2[s_idx].SPos1 = seed_info1.e_pos[p_idx]
							vars_get2[s_idx].SPos2 = seed_info2.e_pos[p_idx]
							vars_get2[s_idx].Strand1 = seed_info1.strand[p_idx]
							vars_get2[s_idx].Strand2 = seed_info2.strand[p_idx]
							vars_get2[s_idx].RInfo = read_info2
						}
					}
				}
			}
		}
		cand_num = append(cand_num, c_num)
		if paired_prob < 1 { //in this case, it is likely the best candidates
			break
		}
	}
	if loop_has_cand != 0 {
		map_qual := 1.0 / float64(cand_num[loop_has_cand-1])
		PrintGetVariants("Final_var", paired_prob, align_prob1, align_prob2, vars_get1, vars_get2)
		for _, var_info := range vars_get1 {
			var_info.MProb = map_qual
			var_results <- var_info
		}
		for _, var_info := range vars_get2 {
			var_info.MProb = map_qual
			var_results <- var_info
		}
		return
	}
	//Get unaligned paired-end reads
	uai := UnAlignInfo{}
	if INPUT_INFO.Debug_mode {
		uai.read_info1 = read_info1
		uai.read_info2 = read_info2
	}
	UNALIGN_INFO_CHAN <- uai
}

//---------------------------------------------------------------------------------------------------
// FindVariantsFromExtension determines variants based on alignment between reads and multi-genomes.
// Extend read and ref from exact matches (found from searching with FM-index).
// Perform backward alignment (for left extension of read and ref) and forward alignment (for right extension)
// between read and multigenome to determine aligned bases as candidates for variant calls.
//---------------------------------------------------------------------------------------------------
func (VC *VarCall) ExtendSeeds(s_pos, e_pos, m_pos int, read, qual []byte, edit_aln_info *EditAlnInfo) ([]*VarInfo, int, int, float64) {

	var i, j, del_len int
	var is_var, is_del bool

	l_read_flank_len := e_pos + PARA_INFO.Seed_backup
	l_read_flank, l_qual_flank := read[:l_read_flank_len], qual[:l_read_flank_len]

	l_ref_flank := make([]byte, 0)
	l_ref_pos_map := make([]int, 0)
	l_align_e_pos := m_pos - 1 + PARA_INFO.Seed_backup
	i = l_align_e_pos
	j = 0 //to check length of l_ref_flank
	for j < l_read_flank_len+PARA_INFO.Indel_backup && i >= 0 {
		if _, is_var = INDEX.VarProf[i]; is_var {
			if del_len, is_del = INDEX.DelVar[i]; is_del {
				if del_len < j && del_len < len(l_ref_flank) {
					l_ref_flank = l_ref_flank[:len(l_ref_flank)-del_len]
					l_ref_pos_map = l_ref_pos_map[:len(l_ref_pos_map)-del_len]
					j -= del_len
				} else {
					return nil, -1, -1, -1
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
	for j < r_read_flank_len+PARA_INFO.Indel_backup && i < len(INDEX.Seq) {
		r_ref_pos_map = append(r_ref_pos_map, i)
		r_ref_flank = append(r_ref_flank, INDEX.Seq[i])
		if _, is_var = INDEX.VarProf[i]; is_var {
			if del_len, is_del = INDEX.DelVar[i]; is_del {
				if del_len < r_read_flank_len-j && i+del_len < len(INDEX.Seq) {
					i += del_len
				} else {
					r_ref_flank = r_ref_flank[:len(r_ref_flank)-1]
					break
				}
			}
		}
		j++
		i++
	}

	PrintComparedReadRef(l_read_flank, l_ref_flank, r_read_flank, r_ref_flank)
	//PrintRefPosMap(l_ref_pos_map, r_ref_pos_map)

	l_Ham_dist, l_Edit_dist, l_bt_mat, l_m, l_n, l_var_pos, l_var_base, l_var_qual, l_var_type :=
		VC.LeftAlign(l_read_flank, l_qual_flank, l_ref_flank, l_align_s_pos, edit_aln_info.Bw_Dist_D,
			edit_aln_info.Bw_Dist_IS, edit_aln_info.Bw_Dist_IT, edit_aln_info.Bw_Trace_D, edit_aln_info.Bw_Trace_IS, edit_aln_info.Bw_Trace_IT, l_ref_pos_map)
	r_Ham_dist, r_Edit_dist, r_bt_mat, r_m, r_n, r_var_pos, r_var_base, r_var_qual, r_var_type :=
		VC.RightAlign(r_read_flank, r_qual_flank, r_ref_flank, r_align_s_pos, edit_aln_info.Fw_Dist_D,
			edit_aln_info.Fw_Dist_IS, edit_aln_info.Fw_Dist_IT, edit_aln_info.Fw_Trace_D, edit_aln_info.Fw_Trace_IS, edit_aln_info.Fw_Trace_IT, r_ref_pos_map)

	prob := l_Ham_dist + l_Edit_dist + r_Ham_dist + r_Edit_dist
	if prob <= PARA_INFO.Prob_thres {
		if l_m > 0 && l_n > 0 {
			l_pos, l_base, l_qual, l_type := VC.LeftAlignEditTraceBack(l_read_flank, l_qual_flank, l_ref_flank, l_m, l_n, l_align_s_pos,
				l_bt_mat, edit_aln_info.Bw_Trace_D, edit_aln_info.Bw_Trace_IS, edit_aln_info.Bw_Trace_IT, l_ref_pos_map)
			l_var_pos = append(l_var_pos, l_pos...)
			l_var_base = append(l_var_base, l_base...)
			l_var_qual = append(l_var_qual, l_qual...)
			l_var_type = append(l_var_type, l_type...)
		}
		PrintMatchTraceInfo(m_pos, l_align_s_pos, prob, l_var_pos, read)
		if r_m > 0 && r_n > 0 {
			r_pos, r_base, r_qual, r_type := VC.RightAlignEditTraceBack(r_read_flank, r_qual_flank, r_ref_flank, r_m, r_n, r_align_s_pos,
				r_bt_mat, edit_aln_info.Fw_Trace_D, edit_aln_info.Fw_Trace_IS, edit_aln_info.Fw_Trace_IT, r_ref_pos_map)
			r_var_pos = append(r_var_pos, r_pos...)
			r_var_base = append(r_var_base, r_base...)
			r_var_qual = append(r_var_qual, r_qual...)
			r_var_type = append(r_var_type, r_type...)
		}
		PrintMatchTraceInfo(m_pos, r_align_s_pos, prob, r_var_pos, read)
		var k int
		var vars_arr []*VarInfo
		for k = 0; k < len(l_var_pos); k++ {
			var_info := new(VarInfo)
			var_info.Pos, var_info.Bases, var_info.BQual, var_info.Type = uint32(l_var_pos[k]), l_var_base[k], l_var_qual[k], l_var_type[k]
			vars_arr = append(vars_arr, var_info)
		}
		for k = 0; k < len(r_var_pos); k++ {
			var_info := new(VarInfo)
			var_info.Pos, var_info.Bases, var_info.BQual, var_info.Type = uint32(r_var_pos[k]), r_var_base[k], r_var_qual[k], r_var_type[k]
			vars_arr = append(vars_arr, var_info)
		}
		return vars_arr, l_align_s_pos, r_align_s_pos, prob
	}
	return nil, -1, -1, -1
}

//---------------------------------------------------------------------------------------------------
// UpdateVariantProb updates probablilities of variants at the location.
//	Input: a variant of type VarInfo.
//---------------------------------------------------------------------------------------------------
func (VC *VarCall) UpdateVariantProb(var_info *VarInfo) {
	pos := var_info.Pos
	a := string(var_info.Bases)
	t := var_info.Type
	//if found new variant locations
	if _, var_prof_exist := VC.VarProb[pos]; !var_prof_exist {
		VC.VarProb[pos] = make(map[string]float64)
		if t == 0 {
			VC.VarProb[pos][string(INDEX.Seq[int(pos)])] = 1 - NEW_SNP_RATE
			VC.VarProb[pos][a] = NEW_SNP_RATE
		} else {
			VC.VarProb[pos][string(INDEX.Seq[int(pos)])] = 1 - NEW_INDEL_RATE
			VC.VarProb[pos][a] = NEW_INDEL_RATE
		}
		VC.VarType[pos] = make(map[string]int)
		VC.VarRNum[pos] = make(map[string]int)
		if INPUT_INFO.Debug_mode {
			VC.ChrDis[pos] = make(map[string][]int)
			VC.ChrDiff[pos] = make(map[string][]int)
			VC.MapProb[pos] = make(map[string][]float64)
			VC.AlnProb[pos] = make(map[string][]float64)
			VC.ChrProb[pos] = make(map[string][]float64)
			VC.StartPos1[pos] = make(map[string][]int)
			VC.StartPos2[pos] = make(map[string][]int)
			VC.Strand1[pos] = make(map[string][]bool)
			VC.Strand2[pos] = make(map[string][]bool)
			VC.VarBQual[pos] = make(map[string][][]byte)
			VC.ReadInfo[pos] = make(map[string][][]byte)
		}
	}
	//if found new variants at existing locations
	var l float64
	var b string
	if _, var_exist := VC.VarProb[pos][a]; !var_exist {
		l = float64(len(VC.VarProb[pos]))
		if t == 0 {
			for b, _ = range VC.VarProb[pos] {
				VC.VarProb[pos][b] = VC.VarProb[pos][b] - (1/l)*NEW_SNP_RATE
			}
			VC.VarProb[pos][a] = NEW_SNP_RATE
		} else {
			for b, _ = range VC.VarProb[pos] {
				VC.VarProb[pos][b] = VC.VarProb[pos][b] - (1/l)*NEW_INDEL_RATE
			}
			VC.VarProb[pos][a] = NEW_INDEL_RATE
		}
	}
	VC.VarType[pos][a] = t
	VC.VarRNum[pos][a] += 1
	if INPUT_INFO.Debug_mode {
		VC.ChrDis[pos][a] = append(VC.ChrDis[pos][a], var_info.CDis)
		VC.ChrDiff[pos][a] = append(VC.ChrDiff[pos][a], var_info.CDiff)
		VC.MapProb[pos][a] = append(VC.MapProb[pos][a], var_info.MProb)
		VC.AlnProb[pos][a] = append(VC.AlnProb[pos][a], var_info.AProb)
		VC.ChrProb[pos][a] = append(VC.ChrProb[pos][a], var_info.IProb)
		VC.StartPos1[pos][a] = append(VC.StartPos1[pos][a], var_info.SPos1)
		VC.StartPos2[pos][a] = append(VC.StartPos2[pos][a], var_info.SPos2)
		VC.Strand1[pos][a] = append(VC.Strand1[pos][a], var_info.Strand1)
		VC.Strand2[pos][a] = append(VC.Strand2[pos][a], var_info.Strand2)
		VC.VarBQual[pos][a] = append(VC.VarBQual[pos][a], var_info.BQual)
		VC.ReadInfo[pos][a] = append(VC.ReadInfo[pos][a], var_info.RInfo)
	}

	var q byte
	p1, p2 := 1.0, 1.0
	for _, q = range var_info.BQual {
		p1 *= Q2P[q]
	}
	for _, q = range var_info.BQual {
		p2 *= Q2E[q]
	}
	var p, p_b float64
	p_a := 0.0
	p_ab := make(map[string]float64)
	_, is_known_var := INDEX.VarProf[int(pos)]
	for b, p_b = range VC.VarProb[pos] {
		if b == a {
			p_ab[b] = p1
		} else {
			if !is_known_var && var_info.Type == 2 {
				p_ab[b] = L2E[len(a)]
			} else {
				p_ab[b] = p2
			}
		}
		p = p_b * p_ab[b]
		p_ab[b] = p
		p_a += p
	}
	i := 0
	for b, p_b = range VC.VarProb[pos] {
		VC.VarProb[pos][b] = p_ab[b] / p_a
		i++
	}
}

//---------------------------------------------------------------------------------------------------
// OutputVarCalls determines variant calls, convert their probabilities to Phred scores, and writes
// all variant calls to file in VCF format.
//---------------------------------------------------------------------------------------------------
func (VC *VarCall) OutputVarCalls() {

	log.Printf("Outputing variant calls...")
	start_time := time.Now()
	f, e := os.Create(INPUT_INFO.Var_call_file)
	if e != nil {
		log.Printf("Error: Create output file, %s, (err: %s)", INPUT_INFO.Var_call_file, e)
		os.Exit(1)
	}
	defer f.Close()

	w := bufio.NewWriter(f)
	if INPUT_INFO.Debug_mode == false {
		w.WriteString("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" +
			"VAR_PROB\tMAP_PROB\tCOM_QUAL\tBASE_NUM\n")
	} else {
		w.WriteString("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" +
			"VAR_PROB\tMAP_PROB\tCOM_QUAL\tBASE_NUM\tBASE_QUAL\tCHR_DIS\tCHR_DIFF\tMAP_PROB\t" +
			"ALN_PROB\tPAIR_PROB\tS_POS1\tBRANCH1\tS_POS2\tBRANCH2\tREAD_HEADER\tALN_BASE\tBASE_NUM\n")
	}
	var var_pos uint32
	Var_Pos := make([]int, 0, len(VC.VarProb))
	for var_pos, _ = range VC.VarProb {
		Var_Pos = append(Var_Pos, int(var_pos))
	}
	sort.Ints(Var_Pos)

	var var_base, var_call, str_qual, str_aln, str_base, str_ivc string
	var line_aln, line_base, line_ivc []string
	var var_prob, var_call_prob, map_prob, p float64
	var is_var bool
	var idx, var_num int
	for _, pos := range Var_Pos {
		var_pos = uint32(pos)
		//Get variant call by considering maximum prob
		var_call_prob = 0
		for var_base, var_prob = range VC.VarProb[var_pos] {
			if var_call_prob < var_prob {
				var_call_prob = var_prob
				var_call = var_base
			}
		}
		//Start getting variant call info
		line_aln = make([]string, 0)
		//#CHROM
		line_aln = append(line_aln, ".")
		//POS
		line_aln = append(line_aln, strconv.Itoa(pos+1))
		//ID
		line_aln = append(line_aln, ".")
		//REF & ALT
		if _, is_var = INDEX.VarProf[pos]; is_var {
			if var_call == string(INDEX.VarProf[pos][0]) { //Do not report known variants which are same with the reference
				continue
			}
			if VC.VarRNum[var_pos][var_call] == 0 { //Do not report known variants at locations without aligned reads
				continue
			}
			line_aln = append(line_aln, string(INDEX.VarProf[pos][0]))
			line_aln = append(line_aln, var_call)
		} else {
			if VC.VarType[var_pos][var_call] >= 0 {
				if VC.VarType[var_pos][var_call] == 2 { //DEL
					line_aln = append(line_aln, var_call)
					line_aln = append(line_aln, string(INDEX.Seq[pos]))
				} else if VC.VarType[var_pos][var_call] == 1 { //INS
					line_aln = append(line_aln, string(INDEX.Seq[pos]))
					line_aln = append(line_aln, var_call)
				} else { //SUB
					//Ignore variants that are identical with ref
					if var_call == string(INDEX.Seq[pos]) {
						continue
					}
					line_aln = append(line_aln, string(INDEX.Seq[pos]))
					line_aln = append(line_aln, var_call)
				}
			} else {
				continue
			}
		}
		//QUAL
		str_qual = strconv.FormatFloat(-10*math.Log10(1-var_call_prob), 'f', 5, 64)
		if str_qual != "+Inf" {
			line_aln = append(line_aln, str_qual)
		} else {
			line_aln = append(line_aln, "1000")
		}
		//FILTER
		line_aln = append(line_aln, ".")
		//INFO
		line_aln = append(line_aln, ".")
		//FORMAT
		line_aln = append(line_aln, ".")
		//IVC-INFO
		line_aln = append(line_aln, strconv.FormatFloat(var_call_prob, 'f', 20, 64))
		map_prob = 1.0
		for _, p = range VC.MapProb[var_pos][var_call] {
			map_prob *= p
		}
		line_aln = append(line_aln, strconv.FormatFloat(map_prob, 'f', 20, 64))
		str_qual = strconv.FormatFloat(-10*math.Log10(1-var_call_prob*map_prob), 'f', 5, 64)
		if str_qual != "+Inf" {
			line_aln = append(line_aln, str_qual)
		} else {
			line_aln = append(line_aln, "1000")
		}
		line_aln = append(line_aln, strconv.Itoa(VC.VarRNum[var_pos][var_call]))
		str_aln = strings.Join(line_aln, "\t")
		if INPUT_INFO.Debug_mode == false {
			w.WriteString(str_aln + "\n")
		} else {
			line_base = make([]string, 0)
			for var_base, var_num = range VC.VarRNum[var_pos] {
				line_base = append(line_base, var_base)
				line_base = append(line_base, strconv.Itoa(var_num))
			}
			str_base = strings.Join(line_base, "\t")
			for idx, _ = range VC.VarBQual[var_pos][var_call] {
				line_ivc = make([]string, 0)
				line_ivc = append(line_ivc, string(VC.VarBQual[var_pos][var_call][idx]))
				line_ivc = append(line_ivc, strconv.Itoa(VC.ChrDis[var_pos][var_call][idx]))
				line_ivc = append(line_ivc, strconv.Itoa(VC.ChrDiff[var_pos][var_call][idx]))
				line_ivc = append(line_ivc, strconv.FormatFloat(VC.MapProb[var_pos][var_call][idx], 'f', 20, 64))
				line_ivc = append(line_ivc, strconv.FormatFloat(VC.AlnProb[var_pos][var_call][idx], 'f', 20, 64))
				line_ivc = append(line_ivc, strconv.FormatFloat(VC.ChrProb[var_pos][var_call][idx], 'f', 20, 64))
				line_ivc = append(line_ivc, strconv.Itoa(VC.StartPos1[var_pos][var_call][idx]))
				line_ivc = append(line_ivc, strconv.FormatBool(VC.Strand1[var_pos][var_call][idx]))
				line_ivc = append(line_ivc, strconv.Itoa(VC.StartPos2[var_pos][var_call][idx]))
				line_ivc = append(line_ivc, strconv.FormatBool(VC.Strand2[var_pos][var_call][idx]))
				line_ivc = append(line_ivc, string(VC.ReadInfo[var_pos][var_call][idx]))
				str_ivc = strings.Join(line_ivc, "\t")
				w.WriteString(str_aln + "\t" + str_ivc + "\t" + str_base + "\n")
			}
		}
	}
	w.Flush()
	PrintProcessMem("Memstats after outputing variant calls")
	output_var_time := time.Since(start_time)
	log.Printf("Time for outputing variant calls:\t%s", output_var_time)
	log.Printf("Finish outputing variant calls.")
	log.Printf("Check results in the file: %s", INPUT_INFO.Var_call_file)
}
