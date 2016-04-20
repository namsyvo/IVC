//---------------------------------------------------------------------------------------------------
// IVC: variantcall.go
// Calling genomic variants based on alignment between reads and multigenomes.
// Known variant information (if available) is taken into account when performing alignment.
// Variants and probability of variant calls to be corect is determined using Bayesian update.
// Copyright 2015 Nam Sy Vo.
//---------------------------------------------------------------------------------------------------

package ivc

import (
	"bufio"
	"bytes"
	"github.com/namsyvo/IVC/fmi"
	"log"
	"math"
	"math/rand"
	"os"
	"runtime/pprof"
	"sort"
	"strconv"
	"strings"
	"sync"
	"time"
)

//---------------------------------------------------------------------------------------------------
// VarCallIndex represents preprocessed information of the reference genome and variant profile,
// includes an FM-index of (reverse of) the multigenome, which is used to speed up variant calling.
// This struct also consists of functions for calling variants.
//---------------------------------------------------------------------------------------------------
type VarCallIndex struct {
	Seq        []byte            // multi-sequence
	SeqLen     int               // length of multi-sequence
	ChrPos     []int             // position (first base) of the chromosome on whole-genome
	ChrName    [][]byte          // chromosome names
	Variants   map[int][][]byte  // variants (position, variants).
	VarAF      map[int][]float32 // allele frequency of variants (position, allele frequency)
	SameLenVar map[int]int       // indicate if variants has same length (SNPs or MNPs)
	DelVar     map[int]int       // length of deletions if variants are deletion
	RevFMI     *fmi.Index        // FM-index of reverse multi-sequence (to do forward search)
}

//--------------------------------------------------------------------------------------------------
// VarProf represents variant profile and related info of the individual genome.
//--------------------------------------------------------------------------------------------------
type VarProf struct {
	// VarProb stores all possible variants at each position and their confident probablilities.
	// Prior probablities will be obtained from reference genomes and variant profiles.
	// Posterior probabilities will be updated during alignment phase based on incomming aligned bases
	VarProb   map[uint32]map[string]float64   // probability of the variant call
	VarType   map[uint32]map[string]int       // pype of variants (0: sub, 1: ins, 2: del; other types will be considered in future)
	VarRNum   map[uint32]map[string]int       // numer of aligned reads corresponding to each variant
	ChrDis    map[uint32]map[string][]int     // chromosomal distance between two aligned read-ends
	ChrDiff   map[uint32]map[string][]int     // chromosomal distance betwwen the aligned postion and true postion (for simulated data)
	MapProb   map[uint32]map[string][]float64 // probability of mapping read to be corect (mapping quality)
	AlnProb   map[uint32]map[string][]float64 // probability of aligning read to be correct (alignment quality)
	ChrProb   map[uint32]map[string][]float64 // probability of insert size to be correct (for pair-end reads)
	StartPos1 map[uint32]map[string][]int     // start position (on read) of alignment of the first end
	StartPos2 map[uint32]map[string][]int     // start position (on read) of alignment of the second end
	Strand1   map[uint32]map[string][]bool    // strand indicator of the first end ("true" if read has same strand with ref, "false" otherwise)
	Strand2   map[uint32]map[string][]bool    // strand indicator of the second end ("true" if read has same strand with ref, "false" otherwise)
	VarBQual  map[uint32]map[string][][]byte  // quality sequences (in FASTQ format) of aligned bases at the variant call position
	ReadInfo  map[uint32]map[string][][]byte  // information sequences (in FASTQ format) of aligned reads (header of reads in FASTQ format)
}

//---------------------------------------------------------------------------------------------------
// VarInfo represents information of detected variants, which serves as temporary variables.
//---------------------------------------------------------------------------------------------------
type VarInfo struct {
	Pos     uint32  // postion of variant (on the reference)
	Bases   []byte  // aligned bases to be the variant
	BQual   []byte  // quality sequences (in FASTQ format) of bases to be the variant
	Type    int     // type of the variant (0: sub, 1: ins, 2: del; other types will be considered in future)
	CDis    int     // chromosomal distance between alignment positions of two read-ends
	CDiff   int     // chromosomal distance between aligned pos and true pos
	MProb   float64 // probability of mapping read corectly (mapping quality)
	AProb   float64 // probability of aligning read correctly (alignment quality)
	IProb   float64 // probability of insert size to be correct (for pair-end reads)
	SPos1   int     // starting position on read1 of exact match (or ending position from backward search with FM-index)
	SPos2   int     // starting position on read2 of exact match (or ending position from backward search with FM-index)
	Strand1 bool    // strand (backward/forward) of read1 of exact match
	Strand2 bool    // strand (backward/forward) of read2 of exact match
	RInfo   []byte  // information sequences (in FASTQ format) of aligned reads (header of reads in FASTQ format)
}

//---------------------------------------------------------------------------------------------------
// UnAlnReadInfo represents information of unaligned-reads, which serves as temporary variables.
//---------------------------------------------------------------------------------------------------
type UnAlnReadInfo struct {
	read_info1 []byte // info of first-end of read
	read_info2 []byte // info of second-end of read
}

//---------------------------------------------------------------------------------------------------
// Set of variant calls, each element cover a certain region on the multigenome.
//---------------------------------------------------------------------------------------------------
var VarCall []*VarProf // number of elements will be set equal to number of cores to run parallel updates

//---------------------------------------------------------------------------------------------------
// NewVariantCaller creates an instance of VarCallIndex and sets up its variables.
// This function will be called from the main program.
//---------------------------------------------------------------------------------------------------
func NewVariantCaller() *VarCallIndex {
	log.Printf("----------------------------------------------------------------------------------------")
	log.Printf("Initializing the variant caller...")
	start_time := time.Now()

	VC := new(VarCallIndex)

	log.Printf("Loading FM-index of the reference...")
	VC.RevFMI = fmi.Load(PARA.Rev_index_file)
	log.Printf("Finish loading FM-index of the reference.")
	if PARA.Debug_mode {
		log.Printf("Memstats (golang name):\tAlloc\tTotalAlloc\tSys\tHeapAlloc\tHeapSys")
		PrintMemStats("Memstats after loading index of multi-sequence")
	}

	log.Printf("Loading the reference...")
	VC.ChrPos, VC.ChrName, VC.Seq = LoadMultiSeq(PARA.Ref_file)
	VC.SeqLen = len(VC.Seq)
	log.Printf("Finish loading the reference.")
	if PARA.Debug_mode {
		PrintMemStats("Memstats after loading multi-sequence")
	}

	log.Printf("Loading the variant profile...")
	VC.Variants, VC.VarAF = LoadVarProf(PARA.Var_prof_file)
	log.Printf("Finish loading the variant profile.")
	if PARA.Debug_mode {
		PrintMemStats("Memstats after loading variant profile")
	}

	log.Printf("Creating auxiliary data structures...")
	VC.SameLenVar = make(map[int]int)
	VC.DelVar = make(map[int]int)
	var same_len_flag, del_flag bool
	var var_len int
	for var_pos, var_bases := range VC.Variants {
		var_len = len(var_bases[0])
		same_len_flag, del_flag = true, true
		for _, val := range var_bases[1:] {
			if var_len != len(val) {
				same_len_flag = false
			}
			if var_len <= len(val) {
				del_flag = false
			}
		}
		if same_len_flag {
			VC.SameLenVar[var_pos] = var_len
		}
		if del_flag {
			VC.DelVar[var_pos] = var_len - 1
		}
	}

	// Set up pre-calculated cost
	// Notice: Phred-encoding factor is set to 33 here. It is better to be determined from input data.
	Q2C = make(map[byte]float64)           // alignment cost based on Phred-scale quality
	Q2E = make(map[byte]float64)           // error probability based on Phred-scale quality
	Q2P = make(map[byte]float64)           // non-error probability based on Phred-scale quality
	L2E = make([]float64, PARA.Read_len+1) // indel-error rate based on indel-length
	var q byte
	for i := 33; i < 105; i++ {
		q = byte(i)
		Q2C[q] = -math.Log10(1.0 - math.Pow(10, -(float64(q)-33)/10.0))
		Q2E[q] = math.Pow(10, -(float64(q)-33)/10.0) / 3.0
		Q2P[q] = 1.0 - math.Pow(10, -(float64(q)-33)/10.0)
	}
	for i := 0; i < PARA.Read_len+1; i++ {
		L2E[i] = math.Pow(INDEL_ERR_RATE, float64(i))
	}

	log.Printf("Finish creating auxiliary data structures.")
	if PARA.Debug_mode {
		PrintMemStats("Memstats after creating auxiliary data structures")
	}

	// Initialize VarCallIndex object for calling variants
	log.Printf("Initializing variant call data structure...")
	VarCall = make([]*VarProf, PARA.Proc_num)
	for rid := 0; rid < PARA.Proc_num; rid++ {
		VarCall[rid] = new(VarProf)
		VarCall[rid].VarProb = make(map[uint32]map[string]float64)
		VarCall[rid].VarType = make(map[uint32]map[string]int)
		VarCall[rid].VarRNum = make(map[uint32]map[string]int)
		if PARA.Debug_mode {
			VarCall[rid].ChrDis = make(map[uint32]map[string][]int)
			VarCall[rid].ChrDiff = make(map[uint32]map[string][]int)
			VarCall[rid].MapProb = make(map[uint32]map[string][]float64)
			VarCall[rid].AlnProb = make(map[uint32]map[string][]float64)
			VarCall[rid].ChrProb = make(map[uint32]map[string][]float64)
			VarCall[rid].StartPos1 = make(map[uint32]map[string][]int)
			VarCall[rid].StartPos2 = make(map[uint32]map[string][]int)
			VarCall[rid].Strand1 = make(map[uint32]map[string][]bool)
			VarCall[rid].Strand2 = make(map[uint32]map[string][]bool)
			VarCall[rid].VarBQual = make(map[uint32]map[string][][]byte)
			VarCall[rid].ReadInfo = make(map[uint32]map[string][][]byte)
		}
	}

	//At this point, assume that all variants are biallelic
	var pos uint32
	var rid int
	c := 0
	for var_pos, var_prof := range VC.Variants {
		pos = uint32(var_pos)
		rid = PARA.Proc_num * var_pos / VC.SeqLen
		VarCall[rid].VarProb[pos] = make(map[string]float64)
		rbase, vbase := string(var_prof[0]), string(var_prof[1])
		VarCall[rid].VarProb[pos][rbase+"|"+rbase] = float64(VC.VarAF[var_pos][0]) * 2.0 / 3.0
		VarCall[rid].VarProb[pos][rbase+"|"+vbase] = float64(VC.VarAF[var_pos][0])/3.0 + float64(VC.VarAF[var_pos][1])/3.0
		VarCall[rid].VarProb[pos][vbase+"|"+vbase] = float64(VC.VarAF[var_pos][1]) * 2.0 / 3.0
		VarCall[rid].VarType[pos] = make(map[string]int)
		VarCall[rid].VarRNum[pos] = make(map[string]int)
		if PARA.Debug_mode {
			VarCall[rid].ChrDis[pos] = make(map[string][]int)
			VarCall[rid].ChrDiff[pos] = make(map[string][]int)
			VarCall[rid].MapProb[pos] = make(map[string][]float64)
			VarCall[rid].AlnProb[pos] = make(map[string][]float64)
			VarCall[rid].ChrProb[pos] = make(map[string][]float64)
			VarCall[rid].StartPos1[pos] = make(map[string][]int)
			VarCall[rid].StartPos2[pos] = make(map[string][]int)
			VarCall[rid].Strand1[pos] = make(map[string][]bool)
			VarCall[rid].Strand2[pos] = make(map[string][]bool)
			VarCall[rid].VarBQual[pos] = make(map[string][][]byte)
			VarCall[rid].ReadInfo[pos] = make(map[string][][]byte)
		}
		if (c+1)%(len(VC.Variants)/10) == 0 {
			log.Println("Finish initializing", (c+1)/(len(VC.Variants)/100), "% of variant call data structure.")
		}
		c++
	}
	log.Printf("Finish initializing variant call data structure.")
	if PARA.Debug_mode {
		PrintMemStats("Memstats after initializing the variant caller")
	}

	index_time := time.Since(start_time)
	log.Printf("Time for initializing the variant caller:\t%s", index_time)
	log.Printf("Finish initializing the variant caller.")
	return VC
}

//---------------------------------------------------------------------------------------------------
// CallVariants searches for variants and updates variant information in VarCallIndex.
// This function will be called from main program.
//---------------------------------------------------------------------------------------------------
func (VC *VarCallIndex) CallVariants() {
	log.Printf("----------------------------------------------------------------------------------------")
	log.Printf("Calling variants...")
	start_time := time.Now()

	read_data := make(chan *ReadInfo, PARA.Proc_num)
	// The channel read_signal is used for signaling between goroutines which run ReadReads and SearchVariants.
	// When a SearchVariants goroutine finish copying a read to its own memory, it signals ReadReads goroutine
	// to scan next reads.
	read_signal := make(chan bool)

	var_info := make([]chan *VarInfo, PARA.Proc_num)
	for i := 0; i < PARA.Proc_num; i++ {
		var_info[i] = make(chan *VarInfo)
	}
	uar_info := make(chan *UnAlnReadInfo)

	// Read input reads
	go VC.ReadReads(read_data, read_signal)

	var wg sync.WaitGroup
	// Search for variants
	for i := 0; i < PARA.Proc_num; i++ {
		wg.Add(1)
		go VC.SearchVariants(read_data, read_signal, var_info, uar_info, &wg)
	}

	//Collect variants from results channel and update variant probabilities
	for i := 0; i < PARA.Proc_num; i++ {
		go func(i int) {
			for vi := range var_info[i] {
				VC.UpdateVariantProb(vi)
			}
		}(i)
	}

	go func() {
		wg.Wait()
		for i := 0; i < PARA.Proc_num; i++ {
			close(var_info[i])
		}
		close(uar_info)
	}()

	// Get unaligned reads and related info
	i := 0
	for uar := range uar_info {
		i++
		if PARA.Debug_mode {
			UNALIGN_READ_INFO = append(UNALIGN_READ_INFO, uar)
		}
	}
	log.Printf("Number of un-aligned reads:\t%d", i)

	if PARA.Debug_mode {
		ProcessNoAlignReadInfo()
		PrintMemStats("Memstats after calling variants")
	}
	call_var_time := time.Since(start_time)
	log.Printf("Time for calling variants:\t%s", call_var_time)
	log.Printf("Finish calling variants.")
}

//---------------------------------------------------------------------------------------------------
// ReadReads reads all reads from input FASTQ files and put them into data channel.
//---------------------------------------------------------------------------------------------------
func (VC *VarCallIndex) ReadReads(read_data chan *ReadInfo, read_signal chan bool) {

	fn1, fn2 := PARA.Read_file_1, PARA.Read_file_2
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
	read_info := InitReadInfo(PARA.Read_len, PARA.Info_len)
	for scanner1.Scan() && scanner2.Scan() {
		read_info.Info1 = read_info.Info1[:len(scanner1.Bytes())]
		read_info.Info2 = read_info.Info2[:len(scanner2.Bytes())]
		copy(read_info.Info1, scanner1.Bytes()) // use 1st line in 1st FASTQ file
		copy(read_info.Info2, scanner2.Bytes()) // use 1st line in 2nd FASTQ file
		scanner1.Scan()
		scanner2.Scan()
		read_info.Read1 = read_info.Read1[:len(scanner1.Bytes())]
		read_info.Read2 = read_info.Read2[:len(scanner2.Bytes())]
		copy(read_info.Read1, scanner1.Bytes()) // use 2nd line in 1st FASTQ file
		copy(read_info.Read2, scanner2.Bytes()) // use 2nd line in 2nd FASTQ file
		scanner1.Scan()                         // ignore 3rd line in 1st FASTQ file
		scanner2.Scan()                         // ignore 3rd line in 2nd FASTQ file
		scanner1.Scan()
		scanner2.Scan()
		read_info.Qual1 = read_info.Qual1[:len(scanner1.Bytes())]
		read_info.Qual2 = read_info.Qual2[:len(scanner2.Bytes())]
		copy(read_info.Qual1, scanner1.Bytes()) // use 4th line in 1st FASTQ file
		copy(read_info.Qual2, scanner2.Bytes()) // use 4th line in 2nd FASTQ file
		if len(read_info.Read1) > 0 && len(read_info.Read2) > 0 {
			read_num++
			read_data <- read_info
			read_signal <- true
		}
		if read_num%100000 == 0 {
			log.Println("Processed " + strconv.Itoa(read_num) + " reads.")
			if PARA.Debug_mode {
				PrintMemStats("Memstats after distributing " + strconv.Itoa(read_num) + " reads")
				pprof.WriteHeapProfile(MEM_FILE)
			}
		}
	}
	log.Printf("Number of reads:\t%d", read_num)
	close(read_data)
}

//---------------------------------------------------------------------------------------------------
// SearchVariants takes data from data channel, searches for variants and put them into results channel.
//---------------------------------------------------------------------------------------------------
func (VC *VarCallIndex) SearchVariants(read_data chan *ReadInfo, read_signal chan bool,
	var_info []chan *VarInfo, uar_info chan *UnAlnReadInfo, wg *sync.WaitGroup) {

	defer wg.Done()

	// Initialize inter-function share variables
	read_info := InitReadInfo(PARA.Read_len, PARA.Info_len)
	edit_aln_info := InitEditAlnInfo(2 * PARA.Read_len)
	seed_pos := make([][]int, 4)
	for i := 0; i < 4; i++ {
		seed_pos[i] = make([]int, PARA.Max_snum)
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

		RevComp(read_info.Read1, read_info.Qual1, read_info.Rev_comp_read1, read_info.Rev_qual1)
		RevComp(read_info.Read2, read_info.Qual2, read_info.Rev_comp_read2, read_info.Rev_qual2)

		VC.SearchVariantsPE(read_info, edit_aln_info, seed_pos, rand_gen, var_info, uar_info)
	}
}

//---------------------------------------------------------------------------------------------------
// SearchVariantsPE searches for variants from alignment between pair-end reads and the multigenome.
// It uses seed-and-extend strategy and looks for the best alignment candidates through several iterations.
//---------------------------------------------------------------------------------------------------
func (VC *VarCallIndex) SearchVariantsPE(read_info *ReadInfo, edit_aln_info *EditAlnInfo, seed_pos [][]int,
	rand_gen *rand.Rand, var_info []chan *VarInfo, uar_info chan *UnAlnReadInfo) {

	//-----------------------------------------------------------------------------------------------
	// in case of simulated reads, get info with specific format of testing dataset
	true_pos1, true_pos2 := 0, 0
	if PARA.Debug_mode {
		read_info1_tokens := bytes.Split(read_info.Info1, []byte{'_'})
		var tmp int64
		var err error
		if read_info1_tokens[0][1] != 'r' && len(read_info1_tokens) >= 4 {
			if tmp, err = strconv.ParseInt(string(read_info1_tokens[2]), 10, 64); err != nil {
				true_pos1 = int(tmp)
			}
			if tmp, err = strconv.ParseInt(string(read_info1_tokens[3]), 10, 64); err != nil {
				true_pos2 = int(tmp)
			}
		} else if len(read_info1_tokens) >= 3 {
			if tmp, err = strconv.ParseInt(string(read_info1_tokens[1]), 10, 64); err != nil {
				true_pos1 = int(tmp)
			}
			if tmp, err = strconv.ParseInt(string(read_info1_tokens[2]), 10, 64); err != nil {
				true_pos2 = int(tmp)
			}
		}
	}
	//-------------------------------------------------------------------------------------------------

	read_info1, read_info2 := make([]byte, len(read_info.Info1)), make([]byte, len(read_info.Info2))
	copy(read_info1, read_info.Info1)
	copy(read_info2, read_info.Info2)
	var vars1, vars2, vars_get1, vars_get2 []*VarInfo
	var l_aln_pos1, l_aln_pos2 int
	var seed_info1, seed_info2 *SeedInfo
	var has_seeds bool
	var aln_dist1, aln_dist2 float64
	var cand_num []int
	var p_idx, s_idx, c_num int

	paired_dist := math.MaxFloat64
	loop_has_cand := 0
	for loop_num := 1; loop_num <= PARA.Iter_num; loop_num++ {
		seed_info1, seed_info2, has_seeds = VC.SearchSeedsPE(read_info, seed_pos, rand_gen)
		if !has_seeds {
			cand_num = append(cand_num, 0)
			continue
		}
		c_num = 0
		for p_idx = 0; p_idx < len(seed_info1.s_pos); p_idx++ {
			// For conventional paired-end sequencing (i.e. Illumina) the directions should be F-R
			// For other kinds of variants (e.g inversions) or other technologies, they can be F-F or R-R
			// For mate-pair, they can be R-F (need to be confirmed)
			if seed_info1.strand[p_idx] == seed_info2.strand[p_idx] {
				continue
			}
			// Search variants for the first end
			if seed_info1.strand[p_idx] == true {
				vars1, _, _, aln_dist1 = VC.ExtendSeeds(seed_info1.s_pos[p_idx], seed_info1.e_pos[p_idx],
					seed_info1.m_pos[p_idx], read_info.Read1, read_info.Qual1, edit_aln_info)
			} else {
				vars1, _, _, aln_dist1 = VC.ExtendSeeds(seed_info1.s_pos[p_idx], seed_info1.e_pos[p_idx],
					seed_info1.m_pos[p_idx], read_info.Rev_comp_read1, read_info.Rev_qual1, edit_aln_info)
			}
			// Search variants for the second end
			if seed_info2.strand[p_idx] == true {
				vars2, _, _, aln_dist2 = VC.ExtendSeeds(seed_info2.s_pos[p_idx], seed_info2.e_pos[p_idx],
					seed_info2.m_pos[p_idx], read_info.Read2, read_info.Qual2, edit_aln_info)
			} else {
				vars2, _, _, aln_dist2 = VC.ExtendSeeds(seed_info2.s_pos[p_idx], seed_info2.e_pos[p_idx],
					seed_info2.m_pos[p_idx], read_info.Rev_comp_read2, read_info.Rev_qual2, edit_aln_info)
			}
			// Currently, variants can be called iff both read-ends can be aligned
			if aln_dist1 != -1 && aln_dist2 != -1 {
				c_num++
				ins_prob := -math.Log10(math.Exp(-math.Pow(math.Abs(float64(l_aln_pos1-l_aln_pos2))-400.0, 2.0) / (2 * 50 * 50)))
				if paired_dist > aln_dist1+aln_dist2 {
					paired_dist = aln_dist1 + aln_dist2
					//PrintGetVariants("Find_min", paired_dist, aln_dist1, aln_dist2, vars1, vars2)
					vars_get1 = make([]*VarInfo, len(vars1)) // need to reset vars_get1 here
					vars_get2 = make([]*VarInfo, len(vars2)) // need to reset vars_get2 here
					loop_has_cand = loop_num
					for s_idx = 0; s_idx < len(vars1); s_idx++ {
						vars_get1[s_idx] = vars1[s_idx]
						if PARA.Debug_mode {
							// Update vars_get1 with other info
							vars_get1[s_idx].CDis = l_aln_pos1 - l_aln_pos2
							vars_get1[s_idx].CDiff = l_aln_pos1 - true_pos1
							vars_get1[s_idx].AProb = aln_dist1
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
						if PARA.Debug_mode {
							// Update vars_get2 with other info
							vars_get2[s_idx].CDis = l_aln_pos1 - l_aln_pos2
							vars_get2[s_idx].CDiff = l_aln_pos2 - true_pos2
							vars_get2[s_idx].AProb = aln_dist2
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
		if paired_dist < PARA.Gap_open { // there are no gaps, in this case, the alignment is likely to be correct
			break
		}
	}
	var rid int
	if loop_has_cand != 0 {
		map_qual := 1.0 / float64(cand_num[loop_has_cand-1]) // a simple mapping quality estimation, might be changed later
		if PARA.Debug_mode {
			PrintGetVariants("Final_var", paired_dist, aln_dist1, aln_dist2, vars_get1, vars_get2)
		}
		for _, var1 := range vars_get1 {
			var1.MProb = map_qual
			rid = PARA.Proc_num * int(var1.Pos) / VC.SeqLen
			var_info[rid] <- var1
		}
		for _, var2 := range vars_get2 {
			var2.MProb = map_qual
			rid = PARA.Proc_num * int(var2.Pos) / VC.SeqLen
			var_info[rid] <- var2
		}
		return
	}
	// Get unaligned paired-end reads
	uar := new(UnAlnReadInfo)
	if PARA.Debug_mode {
		uar.read_info1 = read_info1
		uar.read_info2 = read_info2
	}
	uar_info <- uar
}

//---------------------------------------------------------------------------------------------------
// ExtendSeeds performs alignment between extensions from seeds on reads and multigenomes
// and determines variants from the alignment of both left and right extensions.
//---------------------------------------------------------------------------------------------------
func (VC *VarCallIndex) ExtendSeeds(s_pos, e_pos, m_pos int, read, qual []byte, edit_aln_info *EditAlnInfo) ([]*VarInfo, int, int, float64) {

	var i, j, del_len int
	var is_var, is_del bool

	l_read_flank_len := s_pos + PARA.Seed_backup
	l_read_flank, l_qual_flank := read[:l_read_flank_len], qual[:l_read_flank_len]

	l_ref_flank := make([]byte, 0)
	l_ref_pos_map := make([]int, 0)
	l_aln_e_pos := m_pos - 1 + PARA.Seed_backup
	i = l_aln_e_pos
	j = 0 // to check length of l_ref_flank
	for j < l_read_flank_len+PARA.Indel_backup && i >= 0 {
		if _, is_var = VC.Variants[i]; is_var {
			if del_len, is_del = VC.DelVar[i]; is_del {
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
		l_ref_flank = append(l_ref_flank, VC.Seq[i])
		j++
		i--
	}
	l_aln_s_pos := i + 1

	// Reverse l_ref_pos_map and l_ref_flank to get them in original direction
	for i, j = 0, len(l_ref_pos_map)-1; i < j; i, j = i+1, j-1 {
		l_ref_pos_map[i], l_ref_pos_map[j] = l_ref_pos_map[j], l_ref_pos_map[i]
	}
	for i, j = 0, len(l_ref_flank)-1; i < j; i, j = i+1, j-1 {
		l_ref_flank[i], l_ref_flank[j] = l_ref_flank[j], l_ref_flank[i]
	}

	seed_len := e_pos - s_pos + 1
	r_read_flank_len := len(read) - e_pos - 1 + PARA.Seed_backup
	r_read_flank, r_qual_flank := read[len(read)-r_read_flank_len:], qual[len(read)-r_read_flank_len:]

	r_ref_flank := make([]byte, 0)
	r_ref_pos_map := make([]int, 0)
	r_aln_s_pos := m_pos + seed_len - PARA.Seed_backup
	i = r_aln_s_pos
	j = 0 //to check length of r_ref_flank
	for j < r_read_flank_len+PARA.Indel_backup && i < VC.SeqLen {
		r_ref_pos_map = append(r_ref_pos_map, i)
		r_ref_flank = append(r_ref_flank, VC.Seq[i])
		if _, is_var = VC.Variants[i]; is_var {
			if del_len, is_del = VC.DelVar[i]; is_del {
				if del_len < r_read_flank_len-j && i+del_len < VC.SeqLen {
					i += del_len
				} else {
					//continue to align without remaning part of read and ref
					r_ref_flank = r_ref_flank[:len(r_ref_flank)-1]
					break
				}
			}
		}
		j++
		i++
	}
	if PARA.Debug_mode {
		PrintComparedReadRef(l_read_flank, l_ref_flank, r_read_flank, r_ref_flank)
	}
	l_Ham_dist, l_Edit_dist, l_bt_mat, l_m, l_n, l_var_pos, l_var_base, l_var_qual, l_var_type :=
		VC.LeftAlign(l_read_flank, l_qual_flank, l_ref_flank, l_aln_s_pos, edit_aln_info.l_Dist_D, edit_aln_info.l_Dist_IS,
			edit_aln_info.l_Dist_IT, edit_aln_info.l_Trace_D, edit_aln_info.l_Trace_IS, edit_aln_info.l_Trace_IT, edit_aln_info.l_Trace_K, l_ref_pos_map)
	r_Ham_dist, r_Edit_dist, r_bt_mat, r_m, r_n, r_var_pos, r_var_base, r_var_qual, r_var_type :=
		VC.RightAlign(r_read_flank, r_qual_flank, r_ref_flank, r_aln_s_pos, edit_aln_info.r_Dist_D, edit_aln_info.r_Dist_IS,
			edit_aln_info.r_Dist_IT, edit_aln_info.r_Trace_D, edit_aln_info.r_Trace_IS, edit_aln_info.r_Trace_IT, edit_aln_info.r_Trace_K, r_ref_pos_map)

	aln_dist := l_Ham_dist + l_Edit_dist + r_Ham_dist + r_Edit_dist
	if aln_dist <= PARA.Dist_thres {
		if l_m > 0 && l_n > 0 {
			l_pos, l_base, l_qual, l_type := VC.LeftAlignEditTraceBack(l_read_flank, l_qual_flank, l_ref_flank, l_m, l_n,
				l_aln_s_pos, l_bt_mat, edit_aln_info.l_Trace_D, edit_aln_info.l_Trace_IS, edit_aln_info.l_Trace_IT, edit_aln_info.l_Trace_K, l_ref_pos_map)
			if PARA.Debug_mode {
				PrintVarInfo("LeftAlnitTraceBack, variant info", l_pos, l_base, l_qual)
			}
			l_var_pos = append(l_var_pos, l_pos...)
			l_var_base = append(l_var_base, l_base...)
			l_var_qual = append(l_var_qual, l_qual...)
			l_var_type = append(l_var_type, l_type...)
		}
		if PARA.Debug_mode {
			PrintMatchTraceInfo(m_pos, l_aln_s_pos, aln_dist, l_var_pos, read)
		}
		if r_m > 0 && r_n > 0 {
			r_pos, r_base, r_qual, r_type := VC.RightAlignEditTraceBack(r_read_flank, r_qual_flank, r_ref_flank, r_m, r_n,
				r_aln_s_pos, r_bt_mat, edit_aln_info.r_Trace_D, edit_aln_info.r_Trace_IS, edit_aln_info.r_Trace_IT, edit_aln_info.r_Trace_K, r_ref_pos_map)
			if PARA.Debug_mode {
				PrintVarInfo("RightAlnEditTraceBack, variant info", r_pos, r_base, r_qual)
			}
			r_var_pos = append(r_var_pos, r_pos...)
			r_var_base = append(r_var_base, r_base...)
			r_var_qual = append(r_var_qual, r_qual...)
			r_var_type = append(r_var_type, r_type...)
		}
		if PARA.Debug_mode {
			PrintMatchTraceInfo(m_pos, r_aln_s_pos, aln_dist, r_var_pos, read)
		}
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
		return vars_arr, l_aln_s_pos, r_aln_s_pos, aln_dist
	}
	return nil, -1, -1, -1
}

//---------------------------------------------------------------------------------------------------
// UpdateVariantProb updates probablilities of variants at a variant location using Bayesian update.
//---------------------------------------------------------------------------------------------------
func (VC *VarCallIndex) UpdateVariantProb(var_info *VarInfo) {
	pos := var_info.Pos
	//vtype := var_info.Type
	vbase := strings.Split(string(var_info.Bases), "|")
	rid := PARA.Proc_num * int(pos) / VC.SeqLen
	MUT.Lock()
	// if new variant locations
	if _, var_call_exist := VarCall[rid].VarProb[pos]; !var_call_exist {
		VarCall[rid].VarProb[pos] = make(map[string]float64)
		if len(vbase[0]) == len(vbase[1]) { // SUB
			VarCall[rid].VarProb[pos][vbase[0]+"|"+vbase[0]] = 1 - 1.5*NEW_SNP_RATE
			VarCall[rid].VarProb[pos][vbase[0]+"|"+vbase[1]] = NEW_SNP_RATE
			VarCall[rid].VarProb[pos][vbase[1]+"|"+vbase[1]] = 0.5 * NEW_SNP_RATE
		} else if len(vbase[0]) < len(vbase[1]) { // INS
			VarCall[rid].VarProb[pos][vbase[0]+"|"+vbase[0]] = 1 - 1.5*NEW_INDEL_RATE
			VarCall[rid].VarProb[pos][vbase[0]+"|"+vbase[1]] = NEW_INDEL_RATE
			VarCall[rid].VarProb[pos][vbase[1]+"|"+vbase[1]] = 0.5 * NEW_INDEL_RATE
		} else {
			VarCall[rid].VarProb[pos][vbase[0]+"|"+vbase[0]] = 0.5 * NEW_INDEL_RATE
			VarCall[rid].VarProb[pos][vbase[0]+"|"+vbase[1]] = NEW_INDEL_RATE
			VarCall[rid].VarProb[pos][vbase[1]+"|"+vbase[1]] = 1 - 1.5*NEW_INDEL_RATE
		}
		VarCall[rid].VarType[pos] = make(map[string]int)
		if len(vbase[0]) == len(vbase[1]) { //SUB
			VarCall[rid].VarType[pos][vbase[0]+"|"+vbase[0]] = 0
			VarCall[rid].VarType[pos][vbase[0]+"|"+vbase[1]] = 0
			VarCall[rid].VarType[pos][vbase[1]+"|"+vbase[1]] = 0
		} else if len(vbase[0]) < len(vbase[1]) { //INS
			VarCall[rid].VarType[pos][vbase[0]+"|"+vbase[0]] = 0
			VarCall[rid].VarType[pos][vbase[0]+"|"+vbase[1]] = 1
			VarCall[rid].VarType[pos][vbase[1]+"|"+vbase[1]] = 1
		} else { //DEL
			VarCall[rid].VarType[pos][vbase[0]+"|"+vbase[0]] = 2
			VarCall[rid].VarType[pos][vbase[0]+"|"+vbase[1]] = 2
			VarCall[rid].VarType[pos][vbase[1]+"|"+vbase[1]] = 0
		}
		VarCall[rid].VarRNum[pos] = make(map[string]int)
		VarCall[rid].VarRNum[pos][vbase[0]+"|"+vbase[0]] = 0
		VarCall[rid].VarRNum[pos][vbase[0]+"|"+vbase[1]] = 1
		VarCall[rid].VarRNum[pos][vbase[1]+"|"+vbase[1]] = 1
		if PARA.Debug_mode {
			VarCall[rid].ChrDis[pos] = make(map[string][]int)
			VarCall[rid].ChrDiff[pos] = make(map[string][]int)
			VarCall[rid].MapProb[pos] = make(map[string][]float64)
			VarCall[rid].AlnProb[pos] = make(map[string][]float64)
			VarCall[rid].ChrProb[pos] = make(map[string][]float64)
			VarCall[rid].StartPos1[pos] = make(map[string][]int)
			VarCall[rid].StartPos2[pos] = make(map[string][]int)
			VarCall[rid].Strand1[pos] = make(map[string][]bool)
			VarCall[rid].Strand2[pos] = make(map[string][]bool)
			VarCall[rid].VarBQual[pos] = make(map[string][][]byte)
			VarCall[rid].ReadInfo[pos] = make(map[string][][]byte)
		}
	} else { // if existing variant locations
		var l1, l2 float64
		var b, hap string
		hap_map := make(map[string]bool)
		for b, _ = range VarCall[rid].VarProb[pos] {
			hap_arr := strings.Split(b, "|")
			hap_map[hap_arr[0]], hap_map[hap_arr[1]] = true, true
		}
		// if new variants at existing locations
		if _, var_exist := hap_map[vbase[1]]; !var_exist {
			l1 = float64(len(hap_map) + 1)
			l2 = float64(len(VarCall[rid].VarProb[pos]))
			min_prob := 1.0
			for b, _ = range VarCall[rid].VarProb[pos] {
				if min_prob > VarCall[rid].VarProb[pos][b] {
					min_prob = VarCall[rid].VarProb[pos][b]
				}
			}
			if len(vbase[0]) == len(vbase[1]) {
				for b, _ = range VarCall[rid].VarProb[pos] {
					VarCall[rid].VarProb[pos][b] = VarCall[rid].VarProb[pos][b] - (l1/l2)*min_prob*NEW_SNP_RATE
				}
				for hap, _ = range hap_map {
					VarCall[rid].VarProb[pos][hap+"|"+vbase[1]] = min_prob * NEW_SNP_RATE
					VarCall[rid].VarRNum[pos][hap+"|"+vbase[1]] += 1
				}
				VarCall[rid].VarProb[pos][vbase[1]+"|"+vbase[1]] = min_prob * NEW_SNP_RATE
				VarCall[rid].VarRNum[pos][vbase[1]+"|"+vbase[1]] += 1
			} else if len(vbase[0]) < len(vbase[1]) {
				for b, _ = range VarCall[rid].VarProb[pos] {
					VarCall[rid].VarProb[pos][b] = VarCall[rid].VarProb[pos][b] - (l1/l2)*min_prob*NEW_INDEL_RATE
				}
				for hap, _ = range hap_map {
					VarCall[rid].VarProb[pos][hap+"|"+vbase[1]] = min_prob * NEW_INDEL_RATE
					VarCall[rid].VarRNum[pos][hap+"|"+vbase[1]] += 1
					VarCall[rid].VarType[pos][hap+"|"+vbase[1]] = 1
				}
				VarCall[rid].VarProb[pos][vbase[1]+"|"+vbase[1]] = min_prob * NEW_INDEL_RATE
				VarCall[rid].VarRNum[pos][vbase[1]+"|"+vbase[1]] += 1
				VarCall[rid].VarType[pos][vbase[1]+"|"+vbase[1]] = 1
			} else {
				for b, _ = range VarCall[rid].VarProb[pos] {
					VarCall[rid].VarProb[pos][b] = VarCall[rid].VarProb[pos][b] - (l1/l2)*min_prob*NEW_INDEL_RATE
				}
				for hap, _ = range hap_map {
					VarCall[rid].VarProb[pos][hap+"|"+vbase[1]] = min_prob * NEW_INDEL_RATE
					VarCall[rid].VarRNum[pos][hap+"|"+vbase[1]] += 1
					VarCall[rid].VarType[pos][hap+"|"+vbase[1]] = 2
				}
				VarCall[rid].VarProb[pos][vbase[1]+"|"+vbase[1]] = min_prob * NEW_INDEL_RATE
				VarCall[rid].VarRNum[pos][vbase[1]+"|"+vbase[1]] += 1
				VarCall[rid].VarType[pos][vbase[1]+"|"+vbase[1]] = 2
			}
		} else { //if existing variants
			for b, _ = range VarCall[rid].VarProb[pos] {
				if strings.Contains(b, vbase[1]) {
					VarCall[rid].VarRNum[pos][b] += 1
				}
			}
		}
	}
	if PARA.Debug_mode {
		var_str := string(var_info.Bases)
		VarCall[rid].ChrDis[pos][var_str] = append(VarCall[rid].ChrDis[pos][var_str], var_info.CDis)
		VarCall[rid].ChrDiff[pos][var_str] = append(VarCall[rid].ChrDiff[pos][var_str], var_info.CDiff)
		VarCall[rid].MapProb[pos][var_str] = append(VarCall[rid].MapProb[pos][var_str], var_info.MProb)
		VarCall[rid].AlnProb[pos][var_str] = append(VarCall[rid].AlnProb[pos][var_str], var_info.AProb)
		VarCall[rid].ChrProb[pos][var_str] = append(VarCall[rid].ChrProb[pos][var_str], var_info.IProb)
		VarCall[rid].StartPos1[pos][var_str] = append(VarCall[rid].StartPos1[pos][var_str], var_info.SPos1)
		VarCall[rid].StartPos2[pos][var_str] = append(VarCall[rid].StartPos2[pos][var_str], var_info.SPos2)
		VarCall[rid].Strand1[pos][var_str] = append(VarCall[rid].Strand1[pos][var_str], var_info.Strand1)
		VarCall[rid].Strand2[pos][var_str] = append(VarCall[rid].Strand2[pos][var_str], var_info.Strand2)
		VarCall[rid].VarBQual[pos][var_str] = append(VarCall[rid].VarBQual[pos][var_str], var_info.BQual)
		VarCall[rid].ReadInfo[pos][var_str] = append(VarCall[rid].ReadInfo[pos][var_str], var_info.RInfo)
	}

	pm := 1.0
	for _, q := range var_info.BQual {
		pm *= Q2P[q]
	}
	pi := 1.0
	for _, q := range var_info.BQual {
		pi *= Q2E[q]
	}
	pd := L2E[len(vbase[0])-1]
	p_a := 0.0
	p_ab := make(map[string]float64)
	_, is_known_del := VC.DelVar[int(pos)]
	if PARA.Debug_mode {
		check_pos := uint32(234535446)
		if pos == check_pos {
			log.Println("Before: var_prof, vbase, pm, pi, pd", VarCall[rid].VarProb[pos], vbase, pm, pi, pd, string(var_info.RInfo))
		}
	}
	for b, p_b := range VarCall[rid].VarProb[pos] {
		d := strings.Split(b, "|")
		if len(vbase[0]) > len(vbase[1]) { //DEL
			if vbase[0] == d[0] && vbase[0] == d[1] {
				p_ab[b] = pm
			} else if vbase[0] != d[0] && vbase[0] != d[1] {
				p_ab[b] = pd
			} else {
				p_ab[b] = pm/2.0 + pd/2.0
			}
		} else {
			if is_known_del { //Known DEL
				if len(vbase[0]) == len(vbase[1]) {
					if string(vbase[0][0]) == d[0] && string(vbase[0][0]) == d[1] {
						p_ab[b] = pm
					} else if string(vbase[1][0]) != d[0] && string(vbase[1][0]) != d[1] {
						p_ab[b] = pd
					} else {
						p_ab[b] = pm/2.0 + pd/2.0
					}
				}
			} else {
				if vbase[1] == d[0] && vbase[1] == d[1] {
					p_ab[b] = pm
				} else if vbase[1] != d[0] && vbase[1] != d[1] {
					p_ab[b] = pi
				} else {
					p_ab[b] = pm/2.0 + pi/2.0
				}
			}
		}
		p_a += p_b * p_ab[b]
		if PARA.Debug_mode {
			if pos == check_pos {
				log.Println("Update: b, p_b, p_ab[b], p_a", b, p_b, p_ab[b], p_a)
			}
		}
	}
	for b, p_b := range VarCall[rid].VarProb[pos] {
		VarCall[rid].VarProb[pos][b] = p_b * p_ab[b] / p_a
	}
	if PARA.Debug_mode {
		if pos == check_pos {
			log.Println("After:", VarCall[rid].VarProb[pos])
			log.Println()
		}
	}
	MUT.Unlock()
}

//---------------------------------------------------------------------------------------------------
// OutputVarCalls determines variant calls and writes them to file in VCF format.
//---------------------------------------------------------------------------------------------------
func (VC *VarCallIndex) OutputVarCalls() {
	log.Printf("----------------------------------------------------------------------------------------")
	log.Printf("Outputing variant calls...")
	start_time := time.Now()
	f, e := os.OpenFile(PARA.Var_call_file, os.O_APPEND|os.O_WRONLY, 0666)
	if e != nil {
		log.Panicf("Error: %s", e)
	}
	defer f.Close()

	w := bufio.NewWriter(f)
	var var_pos uint32
	Var_Pos := make([]int, 0)
	for i := 0; i < PARA.Proc_num; i++ {
		for var_pos, _ = range VarCall[i].VarProb {
			Var_Pos = append(Var_Pos, int(var_pos))
		}
	}
	sort.Ints(Var_Pos)
	var var_base, var_call, str_qual, str_aln string
	var line_aln, line_base, line_ivc []string
	var p, var_prob, var_call_prob, map_prob float64
	var i, chr_id, var_num int
	var is_known_var, is_known_del bool
	for _, pos := range Var_Pos {
		var_pos = uint32(pos)
		rid := PARA.Proc_num * pos / VC.SeqLen
		// Get variant call by considering maximum prob
		var_call_prob = 0
		for var_base, var_prob = range VarCall[rid].VarProb[var_pos] {
			if var_call_prob < var_prob {
				var_call_prob = var_prob
				var_call = var_base
			}
		}
		if VarCall[rid].VarRNum[var_pos][var_call] == 0 { // do not report variants without aligned reads (happen at known locations)
			continue
		}
		// Start getting variant call info
		line_aln = make([]string, 0)
		// Get the largest ChrPos that is <= pos
		for chr_id = 0; chr_id < len(VC.ChrPos) && VC.ChrPos[chr_id] <= pos; chr_id++ {
		}
		// #CHROM
		line_aln = append(line_aln, string(VC.ChrName[chr_id-1]))
		// POS
		line_aln = append(line_aln, strconv.Itoa(pos+1-VC.ChrPos[chr_id-1]))
		// ID
		line_aln = append(line_aln, ".")
		// REF & ALT
		hap_arr := strings.Split(var_call, "|")
		if _, is_known_var = VC.Variants[pos]; is_known_var {
			if _, is_known_del = VC.DelVar[pos]; is_known_del {
				//Do not report known variants which are identical with the reference
				if hap_arr[0] == string(VC.Variants[pos][0][0]) && hap_arr[1] == string(VC.Variants[pos][0][0]) {
					continue
				}
				line_aln = append(line_aln, var_call)
				line_aln = append(line_aln, string(VC.Variants[pos][1]))
			} else {
				//Do not report known variants which are identical with the reference
				if hap_arr[0] == string(VC.Variants[pos][0]) && hap_arr[1] == string(VC.Variants[pos][0]) {
					continue
				}
				line_aln = append(line_aln, string(VC.Variants[pos][0]))
				line_aln = append(line_aln, var_call)
			}
		} else {
			//Do not report variants which are identical with the reference
			if hap_arr[0] == string(VC.Seq[pos]) && hap_arr[1] == string(VC.Seq[pos]) {
				continue
			}
			if VarCall[rid].VarType[var_pos][var_call] >= 0 {
				if VarCall[rid].VarType[var_pos][var_call] == 2 { //DEL
					line_aln = append(line_aln, var_call)
					line_aln = append(line_aln, string(VC.Seq[pos]))
				} else { //SUB or INS
					line_aln = append(line_aln, string(VC.Seq[pos]))
					line_aln = append(line_aln, var_call)
				}
			} else {
				continue
			}
		}
		// QUAL
		str_qual = strconv.FormatFloat(-10*math.Log10(1-var_call_prob), 'f', 5, 64)
		if str_qual != "+Inf" {
			line_aln = append(line_aln, str_qual)
		} else {
			line_aln = append(line_aln, "1000")
		}
		// FILTER
		line_aln = append(line_aln, ".")
		// INFO
		if _, is_known_var = VC.Variants[pos]; is_known_var {
			line_aln = append(line_aln, "Known variants")
		} else {
			line_aln = append(line_aln, ".")
		}
		// FORMAT
		line_aln = append(line_aln, ".")
		// IVC-INFO
		line_aln = append(line_aln, strconv.FormatFloat(var_call_prob, 'f', 20, 64))
		map_prob = 1.0
		for _, p = range VarCall[rid].MapProb[var_pos][var_call] {
			map_prob *= p
		}
		line_aln = append(line_aln, strconv.FormatFloat(map_prob, 'f', 20, 64))
		str_qual = strconv.FormatFloat(-10*math.Log10(1-var_call_prob*map_prob), 'f', 5, 64)
		if str_qual != "+Inf" {
			line_aln = append(line_aln, str_qual)
		} else {
			line_aln = append(line_aln, "1000")
		}
		line_aln = append(line_aln, strconv.Itoa(VarCall[rid].VarRNum[var_pos][var_call]))
		read_depth := 0
		for var_base, var_num = range VarCall[rid].VarRNum[var_pos] {
			read_depth += var_num
		}
		line_aln = append(line_aln, strconv.Itoa(read_depth))
		str_aln = strings.Join(line_aln, "\t")
		if PARA.Debug_mode == false {
			w.WriteString(str_aln + "\n")
		} else {
			line_base = make([]string, 0)
			for var_base, var_num = range VarCall[rid].VarRNum[var_pos] {
				line_base = append(line_base, var_base)
				line_base = append(line_base, strconv.Itoa(var_num))
			}
			for i = 0; i < len(VarCall[rid].VarBQual[var_pos][var_call]); i++ {
				line_ivc = make([]string, 0)
				line_ivc = append(line_ivc, string(VarCall[rid].VarBQual[var_pos][var_call][i]))
				line_ivc = append(line_ivc, strconv.Itoa(VarCall[rid].ChrDis[var_pos][var_call][i]))
				line_ivc = append(line_ivc, strconv.Itoa(VarCall[rid].ChrDiff[var_pos][var_call][i]))
				line_ivc = append(line_ivc, strconv.FormatFloat(VarCall[rid].MapProb[var_pos][var_call][i], 'f', 20, 64))
				line_ivc = append(line_ivc, strconv.FormatFloat(VarCall[rid].AlnProb[var_pos][var_call][i], 'f', 20, 64))
				line_ivc = append(line_ivc, strconv.FormatFloat(VarCall[rid].ChrProb[var_pos][var_call][i], 'f', 20, 64))
				line_ivc = append(line_ivc, strconv.Itoa(VarCall[rid].StartPos1[var_pos][var_call][i]))
				line_ivc = append(line_ivc, strconv.FormatBool(VarCall[rid].Strand1[var_pos][var_call][i]))
				line_ivc = append(line_ivc, strconv.Itoa(VarCall[rid].StartPos2[var_pos][var_call][i]))
				line_ivc = append(line_ivc, strconv.FormatBool(VarCall[rid].Strand2[var_pos][var_call][i]))
				line_ivc = append(line_ivc, string(VarCall[rid].ReadInfo[var_pos][var_call][i]))
				w.WriteString(str_aln + "\t" + strings.Join(line_ivc, "\t") + "\t" + strings.Join(line_base, "\t") + "\n")
			}
		}
	}
	w.Flush()
	output_var_time := time.Since(start_time)
	if PARA.Debug_mode {
		PrintMemStats("Memstats after outputing variant calls")
		pprof.StopCPUProfile()
		CPU_FILE.Close()
		MEM_FILE.Close()
	}
	log.Printf("Time for outputing variant calls:\t%s", output_var_time)
	log.Printf("Finish outputing variant calls.")
	log.Printf("------------------------------------------------------")
	log.Printf("Check results in the file: %s", PARA.Var_call_file)
}
