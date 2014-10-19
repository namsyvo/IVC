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
	"bytes"
	"strings"
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
//For debugging
var (
	TRUE_VAR_COMP, TRUE_VAR_PART, TRUE_VAR_NONE map[int][]byte
)
type Debug_info struct {
	read_info1, read_info2 []byte
	align_pos1, align_pos2 int
	snp_num1, snp_num2 int
	snp_pos1, snp_pos2 []uint32
	snp_base1, snp_base2 [][]byte
	snp_baseq1, snp_baseq2 [][]byte
	snp_eval1, snp_eval2 []bool
}
//--------------------------------------------------------------------------------------------------


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
	PARA_INFO = *SetPara(100, 0.001)
	INDEX.Init()

	// Assign all possible SNPs and their prior probabilities from SNP profile.
	S.SNP_Calls = make(map[uint32]map[string]float64)
	
	RAND_GEN = rand.New(rand.NewSource(time.Now().UnixNano()))

	TRUE_VAR_COMP = LoadTrueVar("/data/nsvo/test_data/GRCh37_chr1/refs/mutate-0.3300/variant_comp.txt")
	TRUE_VAR_PART = LoadTrueVar("/data/nsvo/test_data/GRCh37_chr1/refs/mutate-0.3300/variant_part.txt")
	TRUE_VAR_NONE = LoadTrueVar("/data/nsvo/test_data/GRCh37_chr1/refs/mutate-0.3300/variant_none.txt")
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
	debug_info := make(chan Debug_info)
	var wg sync.WaitGroup
	for i := 0; i < INPUT_INFO.Routine_num; i++ {
		go S.FindSNPs(read_data, read_signal, snp_results, &wg, debug_info)
	}
	go func() {
		wg.Wait()
		close(snp_results)
		close(debug_info)
	}()
	//Collect SNPs from results channel and update SNPs and their probabilities
	go func() {
		var snp SNP
		for snp = range snp_results {
			if len(snp.Bases) == 1 {
				S.UpdateSNPProb(snp)
			} else {
				S.UpdateIndelProb(snp)
			}
		}
	}()
	ProcessDebugInfo(debug_info)
}

//Write debug info to files
func ProcessDebugInfo(debug_info chan Debug_info) {
	var i int
	var snp_pos uint32
	for d := range debug_info {
		if d.snp_num1 > 0 {
			for i, snp_pos = range d.snp_pos1 {
				if len(d.snp_base1[i]) == 0 {
					OutputDebugInfo(d, snp_pos, []byte{'.'}, []byte{'I'})
				} else {
					OutputDebugInfo(d, snp_pos, d.snp_base1[i], d.snp_baseq1[i])
				}
			}
		}
		if d.snp_num2 > 0 {
			for i, snp_pos = range d.snp_pos2 {
				if len(d.snp_base2[i]) == 0 {
					OutputDebugInfo(d, snp_pos, []byte{'.'}, []byte{'I'})
				} else {
					OutputDebugInfo(d, snp_pos, d.snp_base2[i], d.snp_baseq2[i])
				}
			}
		}
	}
}

func OutputDebugInfo(d Debug_info, snp_pos uint32, snp_base []byte, snp_baseq []byte) {
	var ok bool
	var val []byte
	if val, ok = TRUE_VAR_COMP[int(snp_pos)]; ok {
		if len(snp_base) == 1 {
			if snp_base[0] == val[0] {
				WriteDebugInfo("tp_snp_comp", d, snp_pos, snp_base, snp_baseq)
			} else {
				WriteDebugInfo("fp_snp_comp", d, snp_pos, snp_base, snp_baseq)
			}
		} else {
			if bytes.Equal(snp_base, val) {
				WriteDebugInfo("tp_indel_comp", d, snp_pos, snp_base, snp_baseq)
			} else {
				WriteDebugInfo("fp_indel_comp", d, snp_pos, snp_base, snp_baseq)
			}
		}
	} else if val, ok = TRUE_VAR_PART[int(snp_pos)]; ok {
		if len(snp_base) == 1 {
			if snp_base[0] == val[0] {
				WriteDebugInfo("tp_snp_part", d, snp_pos, snp_base, snp_baseq)
			} else {
				WriteDebugInfo("fp_snp_part", d, snp_pos, snp_base, snp_baseq)
			}
		} else {
			if bytes.Equal(snp_base, val) {
				WriteDebugInfo("tp_indel_part", d, snp_pos, snp_base, snp_baseq)
			} else {
				WriteDebugInfo("fp_indel_part", d, snp_pos, snp_base, snp_baseq)
			}
		}
	} else if val, ok = TRUE_VAR_NONE[int(snp_pos)]; ok {
		if len(snp_base) == 1 {
			if snp_base[0] == val[0] {
				WriteDebugInfo("tp_snp_none", d, snp_pos, snp_base, snp_baseq)
			} else {
				WriteDebugInfo("fp_snp_none", d, snp_pos, snp_base, snp_baseq)
			}
		} else {
			if bytes.Equal(snp_base, val) {
				WriteDebugInfo("tp_indel_none", d, snp_pos, snp_base, snp_baseq)
			} else {
				WriteDebugInfo("fp_indel_none", d, snp_pos, snp_base, snp_baseq)
			}
		}
	} else {
		if len(snp_base) == 1 {
			WriteDebugInfo("fp_snp_other", d, snp_pos, snp_base, snp_baseq)
		} else {
			WriteDebugInfo("fp_indel_other", d, snp_pos, snp_base, snp_baseq)
		}
	}
}

func WriteDebugInfo(file_name string, d Debug_info, snp_pos uint32, snp_base []byte, snp_baseq []byte) {
	file, _ := os.OpenFile(INPUT_INFO.SNP_call_file + "." + file_name, os.O_APPEND|os.O_WRONLY|os.O_CREATE, 0666)
	defer file.Close()

	if d.align_pos1 != 0 && d.align_pos2 != 0 {
		file.WriteString(strconv.Itoa(d.align_pos1 - d.align_pos2) + "\t")
	} else {
		file.WriteString("None\t")
	}
	file.WriteString(strconv.Itoa(d.align_pos1) + "\t" + strconv.Itoa(d.align_pos2) + "\t")
	file.WriteString(strconv.Itoa(int(snp_pos)) + "\t" + string(snp_base) + "\t")
	file.WriteString(strconv.FormatFloat(math.Pow(10, -(float64(snp_baseq[0]) - 33)/10.0), 'f', 5, 32) + "\t")

	tokens := bytes.Split(d.read_info1, []byte{'_'})
	if len(tokens) >= 11 {
		true_pos1, err1 := strconv.ParseInt(string(tokens[2]), 10, 32)
		true_pos2, err2 := strconv.ParseInt(string(tokens[3]), 10, 32)
		if err1 == nil && err2 == nil {
			file.WriteString(strconv.FormatInt(true_pos1 - true_pos2, 10) + "\t" + strconv.FormatInt(true_pos1, 10) + "\t" + strconv.FormatInt(true_pos2, 10) + "\t")
		} else {
			file.WriteString("None\tNone\tNone\t")
		}
		file.WriteString(string(tokens[10]) + "\n")
	} else {
		file.WriteString("None\tNone\tNone\tNone\n")
	}

	file.Sync()
}

//Load true variants for debugging
func LoadTrueVar(file_name string) map[int][]byte {
	barr := make(map[int][]byte)

	f, err := os.Open(file_name)
	if err != nil {
		fmt.Printf("%v\n", err)
		os.Exit(1)
	}
	defer f.Close()
	br := bufio.NewReader(f)
	for {
		line, err := br.ReadString('\n')
		if err != nil {
			break
		}
		sline := string(line[:len(line)-1])
		split := strings.Split(sline, "\t")
		k, _ := strconv.ParseInt(split[0], 10, 64)
		b := make([]byte, len(split[1]))
		copy(b, split[1])
		barr[int(k)] = b
	}
	return barr
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
	read_info := InitReadInfo(100)
	for scanner1.Scan() && scanner2.Scan() {
		copy(read_info.Info1, scanner1.Bytes()) //use 1st line in input FASTQ file 1
		copy(read_info.Info2, scanner2.Bytes()) //use 1st line in input FASTQ file 2
		read_info.Info1 = read_info.Info1[0:len(scanner1.Bytes())]
		read_info.Info2 = read_info.Info2[0:len(scanner2.Bytes())]
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
		//if read_num%10000 == 0 {
		//	PrintMemStats("Memstats after distributing 10000 reads")
		//}
	}
	close(read_data)
}

//--------------------------------------------------------------------------------------------------
// FindSNPs takes data from data channel, find all possible SNPs and put them into results channel.
//--------------------------------------------------------------------------------------------------
func (S *SNP_Prof) FindSNPs(read_data chan *ReadInfo, read_signal chan bool, snp_results chan SNP, 
	wg *sync.WaitGroup, debug_info chan Debug_info) {

	//Initialize inter-function share variables
	read_info := InitReadInfo(PARA_INFO.Read_len)
	align_info := InitAlignInfo(PARA_INFO.Read_len)
	match_pos := make([]int, PARA_INFO.Max_match)

	wg.Add(1)
	defer wg.Done()
	for read := range read_data {
		//PrintMemStats("Before copying all info from data chan")
		copy(read_info.Info1, read.Info1)
		copy(read_info.Info2, read.Info2)
		read_info.Info1 = read_info.Info1[0:len(read.Info1)]
		read_info.Info2 = read_info.Info2[0:len(read.Info2)]
		copy(read_info.Read1, read.Read1)
		copy(read_info.Read2, read.Read2)
		copy(read_info.Qual1, read.Qual1)
		copy(read_info.Qual2, read.Qual2)
		<- read_signal
		//PrintMemStats("After copying all info from data chan")
		RevComp(read_info.Read1, read_info.Qual1, read_info.Rev_read1, read_info.Rev_comp_read1, read_info.Comp_read1, read_info.Rev_qual1)
		//PrintMemStats("After calculating RevComp for Read1")
		RevComp(read_info.Read2, read_info.Qual2, read_info.Rev_read2, read_info.Rev_comp_read2, read_info.Comp_read2, read_info.Rev_qual2)
		//PrintMemStats("After calculating RevComp for Read2")
		S.FindSNPsFromReads(read_info, snp_results, align_info, match_pos, debug_info)
		//PrintMemStats("After finding all SNPs from reads")
	}
}

//--------------------------------------------------------------------------------------------------
// FindSNPsFromReads returns SNPs found from alignment between pair-end reads and the multigenome.
// This version treats each end of the reads independently.
//--------------------------------------------------------------------------------------------------
func (S *SNP_Prof) FindSNPsFromReads(read_info *ReadInfo, snp_results chan SNP, align_info *AlignInfo, match_pos []int, debug_info chan Debug_info) {

	var snps1, snps2 []SNP
	var left_most_pos1, left_most_pos2 int
	//Find SNPs for the first end
	//PrintMemStats("Before FindSNPsFromEnd1")
	snps1, left_most_pos1 = S.FindSNPsFromEachEnd(read_info.Read1, read_info.Rev_read1, read_info.Rev_comp_read1, 
		read_info.Comp_read1, read_info.Qual1, read_info.Rev_qual1, align_info, match_pos)
	//PrintMemStats("After FindSNPsFromEnd1")

	//Find SNPs for the second end
	//PrintMemStats("Before FindSNPsFromEnd2")
	snps2, left_most_pos2 = S.FindSNPsFromEachEnd(read_info.Read2, read_info.Rev_read2, read_info.Rev_comp_read2, 
		read_info.Comp_read2, read_info.Qual2, read_info.Rev_qual2, align_info, match_pos)
	//PrintMemStats("After FindSNPsFromEnd2")

	//Will process constrants of two ends here
	//...
	var d Debug_info
	d.read_info1 = read_info.Info1
	d.read_info2 = read_info.Info2
	d.align_pos1 = left_most_pos1
	d.align_pos2 = left_most_pos2
	d.snp_num1 = len(snps1)
	d.snp_num2 = len(snps2)

	var snp SNP
	if len(snps1) > 0 {
		for _, snp = range snps1 {
			snp_results <- snp
			d.snp_pos1 = append(d.snp_pos1, snp.Pos)
			d.snp_base1 = append(d.snp_base1, snp.Bases)
			d.snp_baseq1 = append(d.snp_baseq1, snp.BaseQ)
		}
	}
	if len(snps2) > 0 {
		for _, snp = range snps2 {
			snp_results <- snp
			d.snp_pos2 = append(d.snp_pos2, snp.Pos)
			d.snp_base2 = append(d.snp_base2, snp.Bases)
			d.snp_baseq2 = append(d.snp_baseq2, snp.BaseQ)
		}
	}
	debug_info <- d
}

//---------------------------------------------------------------------------------------------------
// FindSNPsFromEachEnd find SNPs from alignment between read (one end) and multigenome.
//---------------------------------------------------------------------------------------------------
func (S *SNP_Prof) FindSNPsFromEachEnd(read, rev_read, rev_comp_read, comp_read, qual, rev_qual []byte, 
	align_info *AlignInfo, match_pos []int) ([]SNP, int) {
	var has_seeds bool
	var p, s_pos, e_pos int
	var loop_num, match_num int
	var snps []SNP

	p = INPUT_INFO.Start_pos
	loop_num = 1
	var left_most_pos int
	for loop_num <= PARA_INFO.Iter_num {
		//fmt.Println(loop_num, "\tread2")
		//PrintMemStats("Before FindSeeds, original_read, loop_num " + strconv.Itoa(loop_num))
		s_pos, e_pos, match_num, has_seeds = INDEX.FindSeeds(read, rev_read, p, match_pos)
		//PrintMemStats("After FindSeeds, original_read, loop_num " + strconv.Itoa(loop_num))
		if has_seeds {
			//fmt.Println("read2, has seed\t", s_pos, "\t", e_pos, "\t", string(read_info.Read2))
			//PrintMemStats("Before FindSNPsFromMatch, original_read, loop_num " + strconv.Itoa(loop_num))
			snps, left_most_pos = S.FindSNPsFromMatch(read, qual, s_pos, e_pos, match_pos, match_num, align_info)
			//PrintMemStats("After FindSeeds, original_read, loop_num " + strconv.Itoa(loop_num))
			if len(snps) > 0 {
				//fmt.Println("read2, has snp\t", s_pos, "\t", e_pos, "\t", string(read_info.Read2))
				return snps, left_most_pos
			}
		}
		//Find SNPs for the reverse complement of the second end
		//PrintMemStats("Before FindSeeds, revcomp_read, loop_num " + strconv.Itoa(loop_num))
		s_pos, e_pos, match_num, has_seeds = INDEX.FindSeeds(rev_comp_read, comp_read, p, match_pos)
		//PrintMemStats("After FindSeeds, revcomp_read, loop_num " + strconv.Itoa(loop_num))
		if has_seeds {
			//fmt.Println("rc_read2, has seed\t", s_pos, "\t", e_pos, "\t", string(read_info.Rev_comp_read2))
			//PrintMemStats("Before FindSNPsFromMatch, revcomp_read, loop_num " + strconv.Itoa(loop_num))
			snps, left_most_pos = S.FindSNPsFromMatch(rev_comp_read, rev_qual, s_pos, e_pos, match_pos, match_num, align_info)
			//PrintMemStats("After FindSNPsFromMatch, revcomp_read, loop_num " + strconv.Itoa(loop_num))
			if len(snps) > 0 {
				//fmt.Println("rc_read2, has snp\t", s_pos, "\t", e_pos, "\t", string(read_info.Rev_comp_read2))
				return snps, left_most_pos
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
	return snps, left_most_pos
}

//---------------------------------------------------------------------------------------------------
// FindSNPsFromMatch finds SNPs from extensions of matches between read (one end) and multigenome.
//---------------------------------------------------------------------------------------------------
func (S *SNP_Prof) FindSNPsFromMatch(read, qual []byte, s_pos, e_pos int, 
	match_pos []int, match_num int, align_info *AlignInfo) ([]SNP, int) {

	var pos, k, dis int
	var left_snp_pos, right_snp_pos, left_snp_idx, right_snp_idx []int
	var left_snp_val, right_snp_val [][]byte
	var snps []SNP
	var snp SNP

	min_dis := INF
	var left_most_pos, min_pos int
	for i := 0; i < match_num; i++ {
		pos = match_pos[i]
		//PrintMemStats("Before FindExtensions, match_num " + strconv.Itoa(i))
		dis, left_snp_pos, left_snp_val, left_snp_idx, right_snp_pos, right_snp_val, right_snp_idx, left_most_pos =
			 INDEX.FindExtensions(read, s_pos, e_pos, pos, align_info)
		//PrintMemStats("After FindExtensions, match_num " + strconv.Itoa(i))
		if dis <= PARA_INFO.Dist_thres {
			if len(left_snp_pos) != 0 || len(right_snp_pos) != 0 {
				if min_dis > dis {
					min_dis = dis
					min_pos = left_most_pos
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
	return snps, min_pos
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