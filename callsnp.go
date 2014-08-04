//---------------------------------------------------------------------------------------------------
// Calling SNPs based on read-multigenome alignment results.
// Copyright 2014 Nam Sy Vo
//---------------------------------------------------------------------------------------------------

package isc

import (
	"fmt"
	"log"
	"math/rand"
	"os"
	"runtime"
	"strconv"
	"time"
	"sort"
	//"bytes"
	"sync"
)

//Global variables
var (
	INDEX       Index      //SNP caller Index
	RAND_GEN    *rand.Rand //pseudo-random number generator
	SEARCH_MODE int        //searching mode for finding seeds
	START_POS   int        //start postion on reads to search
	SEARCH_STEP int        //step for searching in deterministic mode
)

//SNP stores SNP info
type SNP struct {
	SNP_Idx int
	SNP_Val []byte
}

//SNPProf stores SNP Calling information and defines SNP calling functions
type SNPProf struct {
	SNP_Prof map[int][][]byte  // to store all possible SNPs at each position
	SNP_Conf map[int][]float32 // to store quality of all possible SNPS at each position
	SNP_Call map[int][]byte    // to store SNP call at each position
	SNP_Prob map[int][]int     // to store percentage of called SNP among all possilbe SNPs at each position
}

//Initialize parameters
func (S *SNPProf) Init(input_info InputInfo, para_info ParaInfo) {

	memstats := new(runtime.MemStats)

	INDEX.Init(input_info, para_info)
	RAND_GEN = rand.New(rand.NewSource(time.Now().UnixNano()))
	SEARCH_MODE = input_info.Search_mode
	START_POS = input_info.Start_pos
	SEARCH_STEP = input_info.Search_step
	fmt.Println("SEARCH_MODE: ", SEARCH_MODE)
	fmt.Println("START_POS: ", START_POS)
	fmt.Println("SEARCH_STEP: ", SEARCH_STEP)

	runtime.ReadMemStats(memstats)
	log.Printf("callsnp.go: memstats after initializing indexes:\t%d\t%d\t%d\t%d\t%d", memstats.Alloc, memstats.TotalAlloc, memstats.Sys, memstats.HeapAlloc, memstats.HeapSys)

	S.SNP_Prof = make(map[int][][]byte)
	S.SNP_Conf = make(map[int][]float32)
	S.SNP_Call = make(map[int][]byte)
	S.SNP_Prob = make(map[int][]int)
	S.SNP_Conf = make(map[int][]float32)

	for k, v := range INDEX.SNP_AF {
		S.SNP_Conf[k] = v
	}

	runtime.ReadMemStats(memstats)
	log.Printf("callsnp.go: memstats after initializing SNP Prof:\t%d\t%d\t%d\t%d\t%d", memstats.Alloc, memstats.TotalAlloc, memstats.Sys, memstats.HeapAlloc, memstats.HeapSys)

}

//--------------------------------------------------------------------------------------------------
// FindSNPProfile returns SNP profile of new genome based on SNP profile of reference multi-genomes
// and alignment between reads and multi-genomes.
//
//	- Find SNPs for pairend reads, treat each end separately and independently.
//--------------------------------------------------------------------------------------------------
func (S *SNPProf) FindSNP(read_info ReadInfo, align_info AlignInfo, match_pos []int) []SNP {

	var has_seeds bool
	var p, s_pos, e_pos int
	var loop_num, match_num int
	var SNPs, snps []SNP

	//Find SNPs for the first end
	p = START_POS
	loop_num = 1
	for loop_num <= ITER_NUM {
		//fmt.Println(loop_num, "\tread1")
		s_pos, e_pos, match_num, has_seeds = INDEX.FindSeeds(read_info.Read1, read_info.Rev_read1, p, match_pos)
		if has_seeds {
			//fmt.Println("read1, has seed\t", s_pos, "\t", e_pos, "\t", string(read_info.Read1))
			snps = S.FindSNPFromMatch(read_info.Read1, s_pos, e_pos, match_pos, match_num, align_info)
			if len(snps) > 0 {
				//fmt.Println("read1, has snp\t", s_pos, "\t", e_pos, "\t", string(read_info.Read1))
				//fmt.Println(loop_num, "\tori1\t", string(read1))
				SNPs = append(SNPs, snps...)
				break
			}
		}
		//Find SNPs for the reverse complement of the first end
		s_pos, e_pos, match_num, has_seeds = INDEX.FindSeeds(read_info.Rev_comp_read1, read_info.Comp_read1, p, match_pos)
		if has_seeds {
			//fmt.Println("rc_read1, has seed\t", s_pos, "\t", e_pos, "\t", string(read_info.Rev_comp_read1))
			snps = S.FindSNPFromMatch(read_info.Rev_comp_read1, s_pos, e_pos, match_pos, match_num, align_info)
			if len(snps) > 0 {
				//fmt.Println("rc_read1, has snp\t", s_pos, "\t", e_pos, "\t", string(read_info.Rev_comp_read1))
				//fmt.Println(loop_num, "\trev1\t", string(rev_read1))
				SNPs = append(SNPs, snps...)
				break
			}
		}
		//Take a new position to search
		if SEARCH_MODE == 1 {
			p = RAND_GEN.Intn(read_info.Read_len-1) + 1
		} else if SEARCH_MODE == 2 {
			p = p + SEARCH_STEP
		}
		loop_num++
	}

	//Find SNPs for the second end
	p = START_POS
	loop_num = 1
	for loop_num <= ITER_NUM {
		//fmt.Println(loop_num, "\tread2")
		s_pos, e_pos, match_num, has_seeds = INDEX.FindSeeds(read_info.Read2, read_info.Rev_read2, p, match_pos)
		if has_seeds {
			//fmt.Println("read2, has seed\t", s_pos, "\t", e_pos, "\t", string(read_info.Read2))
			snps = S.FindSNPFromMatch(read_info.Read2, s_pos, e_pos, match_pos, match_num, align_info)
			if len(snps) > 0 {
				//fmt.Println("read2, has snp\t", s_pos, "\t", e_pos, "\t", string(read_info.Read2))
				//fmt.Println(loop_num, "\tori2\t", string(read2))
				SNPs = append(SNPs, snps...)
				return SNPs
			}
		}
		//Find SNPs for the reverse complement of the second end
		s_pos, e_pos, match_num, has_seeds = INDEX.FindSeeds(read_info.Rev_comp_read2, read_info.Comp_read2, p, match_pos)
		if has_seeds {
			//fmt.Println("rc_read2, has seed\t", s_pos, "\t", e_pos, "\t", string(read_info.Rev_comp_read2))
			snps = S.FindSNPFromMatch(read_info.Rev_comp_read2, s_pos, e_pos, match_pos, match_num, align_info)
			if len(snps) > 0 {
				//fmt.Println("rc_read2, has snp\t", s_pos, "\t", e_pos, "\t", string(read_info.Rev_comp_read2))
				//fmt.Println(loop_num, "\trev2\t", string(rev_read2))
				SNPs = append(SNPs, snps...)
				return SNPs
			}
		}
		//Take a new position to search
		if SEARCH_MODE == 1 {
			p = RAND_GEN.Intn(read_info.Read_len-1) + 1
		} else if SEARCH_MODE == 2 {
			p = p + SEARCH_STEP
		}
		loop_num++
	}
	return SNPs
}

//---------------------------------------------------------------------------------------------------
// UpdateSNPCall updates SNP Call found from alignment between reads and multi-genomes.
//---------------------------------------------------------------------------------------------------
func (S *SNPProf) FindSNPFromMatch(read []byte, s_pos, e_pos int, match_pos []int, match_num int, align_info AlignInfo) []SNP {

	var pos, dis, k int
	var left_snp_idx, right_snp_idx []int
	var left_snp_val, right_snp_val [][]byte
	var snp SNP
	var snps []SNP

	for i := 0; i < match_num; i++ {
		pos = match_pos[i]
		//Call IntervalHasSNP to determine whether extension is needed
		//if index.IntervalHasSNP(A.SORTED_SNP_POS, pos - e_pos, pos - e_pos + len(read1)) {
		dis, left_snp_idx, left_snp_val, right_snp_idx, right_snp_val = INDEX.FindExtensions(read, s_pos, e_pos, pos, align_info)
		if dis <= DIST_THRES {
			//Determine SNP profile
			if len(left_snp_idx) == 0 && len(right_snp_idx) == 0 {
				continue
			} else {
				for k = 0; k < len(left_snp_idx); k++ {
					snp.SNP_Idx, snp.SNP_Val = left_snp_idx[k], left_snp_val[k]
					snps = append(snps, snp)
					//fmt.Println("left snp\t", idx, "\t", string(val))
				}
				for k = 0; k < len(right_snp_idx); k++ {
					snp.SNP_Idx, snp.SNP_Val = right_snp_idx[k], right_snp_val[k]
					snps = append(snps, snp)
					//fmt.Println("right snp\t", idx, "\t", string(val))
				}
			}
		}
		//}
	}
	return snps
}

/*
//---------------------------------------------------------------------------------------------------
// CalcQual calculates called SNP quality.
//---------------------------------------------------------------------------------------------------
func (S *SNPProf) CalcQual() {

	cond_prob := make([]float32, len(S.SNP_Conf[snp_pos]))
	for idx, base_conf := range S.SNP_Conf[snp_pos] {
		SNP := index.SNP_PROF[snp_pos][idx]
		for _, snp := range snp_arr {
			if bytes.Equal(snp, SNP) {
				//fmt.Println("Diff: ", string(snp), string(SNP), base_conf)
				cond_prob[idx] = base_conf * float32(1 - 0.01)
			} else {
				//fmt.Println("Same: ", string(snp), string(SNP), base_conf)
				cond_prob[idx] = base_conf * float32(0.01)
			}
		}
	}
	base_prob := float32(0)
	for _, value := range(cond_prob) {
		base_prob += value
	}
	for idx, value := range(cond_prob) {
		S.SNP_Conf[snp_pos][idx] = value / base_prob
	}
}
*/

type CallSNP struct {
	SNPPos int
	SNPCall []byte
	SNPProb []int
}

//-----------------------------------------------------------------------------------------------------
// CallSNP returns called SNPs and related information based on SNP profile constructed from
// alignment between reads and multi-genomes.
//-----------------------------------------------------------------------------------------------------
func (S *SNPProf) CallSNP(routine_num int) {

	SNP_Pos := make([]int, 0, len(S.SNP_Prof))
	for snp_pos, _ := range S.SNP_Prof {
		SNP_Pos = append(SNP_Pos, snp_pos)
	}
	sort.Ints(SNP_Pos)

	snp_pos_chan := make(chan int, 32)
	result_chan := make(chan CallSNP)
	go func() {
		for _, snp_pos := range SNP_Pos {
			snp_pos_chan <- snp_pos
		}
	}()
	var wg sync.WaitGroup
	for i := 0; i < routine_num; i++ {
		go func() {
			wg.Add(1)
			defer wg.Done()
			snp_pos := <- snp_pos_chan
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
			var callSNP CallSNP
			callSNP.SNPPos = snp_pos
			callSNP.SNPCall = []byte(major_snp)
			callSNP.SNPProb = []int{major_num, len(S.SNP_Prof[snp_pos])}
			result_chan <- callSNP
		}()
	}
	go func() {
		wg.Wait()
		close(result_chan)
	}()
	for snp_prof := range result_chan {
		S.SNP_Call[snp_prof.SNPPos] = snp_prof.SNPCall
		S.SNP_Prob[snp_prof.SNPPos] = snp_prof.SNPProb
	}		
}

//-------------------------------------------------------------------------------------------------------
// SNPCall_tofile writes called SNPs and related information to given output file in tab-delimited format
//-------------------------------------------------------------------------------------------------------
func (S *SNPProf) SNPCall_tofile(file_name string) {

	file, err := os.Create(file_name)
	if err != nil {
		return
	}
	defer file.Close()
	var nums []int
	var str_snp_pos, str_snp, str_snp_num1, str_snp_num2, str_snp_prob string
/*
	var flag bool
	var idx int
	var value []byte
	var str_snp_conf string
*/
	//log.Printf("New Alleles:\n")
	for snp_pos, snp_val := range S.SNP_Call {
		str_snp_pos = strconv.Itoa(snp_pos)
		str_snp = string(snp_val)
		nums = S.SNP_Prob[snp_pos]
		str_snp_num1, str_snp_num2 = strconv.Itoa(nums[0]), strconv.Itoa(nums[1])
		str_snp_prob = strconv.FormatFloat(float64(nums[0])/float64(nums[1]), 'f', 5, 32)

		if str_snp != "" {
			_, err = file.WriteString(str_snp_pos + "\t" + str_snp + "\t" +
				str_snp_num1 + "\t" + str_snp_num2 + "\t" + str_snp_prob + "\t")
		} else {
			_, err = file.WriteString(str_snp_pos + "\t.\t" +
				str_snp_num1 + "\t" + str_snp_num2 + "\t" + str_snp_prob + "\t")
		}
/*
		//Write SNP Qual - testing////////////////////////////
		flag = false
		for idx, value = range INDEX.SNP_PROF[snp_pos] {
			if bytes.Equal(value, S.SNP_Call[snp_pos]) {
				str_snp_conf = strconv.FormatFloat(float64(S.SNP_Conf[snp_pos][idx]), 'f', 5, 32)
				_, err = file.WriteString(str_snp_conf + "\n")
				flag = true
			}
		}
		if !flag {
			//_, err = file.WriteString(".\n")
			//log.Printf("%s\t%s\n", str_snp_pos, str_snp)
			//for _, val := range INDEX.SNP_PROF[snp_pos] {
			//	log.Printf("%s\t", string(val))
			//}
			//log.Printf("\n")
		}
		//////////////////////////////////////////////////////
*/
		if err != nil {
			fmt.Println(err)
			break
		}
	}
}
