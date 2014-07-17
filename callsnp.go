//---------------------------------------------------------------------------------------------------
// Calling SNPs based on read-multigenome alignment results.
// Copyright 2014 Nam Sy Vo
//---------------------------------------------------------------------------------------------------

package isc

import (
    "fmt"
    "os"
    "math/rand"
    "time"
    "strconv"
	"bytes"
    "sort"
	"runtime"
	"log"
)

//Global variables
var (
	INDEX Index
	RAND_GEN *rand.Rand
	START_POS int //start postion on reads to search
	READ_LEN int //length of reads
)

//SNP caller object with Parameters
type SNPProf struct {
    SNP_Prof map[int][][]byte // to store all possible SNPs at each position
	SNP_Conf map[int][]float32 // to store quality of all possible SNPS at each position
    SNP_Call map[int][]byte // to store SNP call at each position
    SNP_Prob map[int][]int // to store percentage of called SNP among all possilbe SNPs at each position
}

//Initialize parameters
func (S *SNPProf) Init(genome_file, snp_file, index_file, rev_index_file string, read_len int,
	seq_err float32, k, a, n int) {

	memstats := new(runtime.MemStats)

    INDEX.Init(genome_file, snp_file, index_file, rev_index_file, read_len, seq_err, k, a, n)
    RAND_GEN = rand.New(rand.NewSource(time.Now().UnixNano()))
	START_POS = 0
	READ_LEN = read_len
    fmt.Println("START_POS: ", START_POS)
    fmt.Println("READ_LEN: ", READ_LEN)

	runtime.ReadMemStats(memstats)
	log.Printf("callsnp.go: memstats after initializing indexes:\t%d\t%d\t%d\t%d\t%d", memstats.Alloc, memstats.TotalAlloc,
		memstats.Sys, memstats.HeapAlloc, memstats.HeapSys)

    S.SNP_Prof = make(map[int][][]byte)
    S.SNP_Conf = make(map[int][]float32)
    S.SNP_Call = make(map[int][]byte)
    S.SNP_Prob = make(map[int][]int)
	S.SNP_Conf = INDEX.SNP_AF

    runtime.ReadMemStats(memstats)
    log.Printf("callsnp.go: memstats after initializing SNP Prof:\t%d\t%d\t%d\t%d\t%d", memstats.Alloc, memstats.TotalAlloc, memstats.Sys, memstats.HeapAlloc, memstats.HeapSys)

}

//--------------------------------------------------------------------------------------------------
// FindSNPProfile returns SNP profile of new genome based on SNP profile of reference multi-genomes
// and alignment between reads and multi-genomes.
//
//	- Find SNPs for pairend reads, treat each end separately and independently.
//--------------------------------------------------------------------------------------------------
func (S *SNPProf) UpdateSNPCall(read1, read2, rev_read1, rev_read2, rev_comp_read1, rev_comp_read2, comp_read1, comp_read2 []byte, bw_snp_idx, fw_snp_idx []int, bw_snp_val, fw_snp_val [][]byte, bw_D, fw_D [][]int, bw_T, fw_T [][][]byte, match_pos []int) bool {

	var hasSeeds, has_snp_1, has_snp_2 bool
    var p, s_pos, e_pos int
	var loop_num, match_num int

    //Find SNPs for the first end
	p = START_POS
    loop_num = 1
	has_snp_1 = false
    for loop_num <= ITER_NUM {
        s_pos, e_pos, match_num, hasSeeds = INDEX.FindSeeds(read1, rev_read1, p, match_pos)
        if hasSeeds {
			has_snp_1 = S.FindSNPCall(read1, s_pos, e_pos, match_pos, match_num, bw_snp_idx, fw_snp_idx, bw_snp_val, fw_snp_val, bw_D, fw_D, bw_T, fw_T)
            if has_snp_1 {
		        //fmt.Println(loop_num, "\tori1\t", string(read1))
                break
            }
        }
        //Find SNPs for the reverse complement of the first end
        s_pos, e_pos, match_num, hasSeeds = INDEX.FindSeeds(rev_comp_read1, comp_read1, p, match_pos)
        if hasSeeds {
			has_snp_1 = S.FindSNPCall(rev_comp_read1, s_pos, e_pos, match_pos, match_num, bw_snp_idx, fw_snp_idx, bw_snp_val, fw_snp_val, bw_D, fw_D, bw_T, fw_T)
            if has_snp_1 {
		        //fmt.Println(loop_num, "\trev1\t", string(rev_read1))
                break
            }
        }
        //Take a random position to search
        p = RAND_GEN.Intn(len(read1) - 1) + 1
        loop_num++
    }

    //Find SNPs for the second end
    p = START_POS
    loop_num = 1
    has_snp_2 = false
    for loop_num <= ITER_NUM {
        //Call FindSeeds to determine seed
        s_pos, e_pos, match_num, hasSeeds = INDEX.FindSeeds(read2, rev_read2, p, match_pos)
        if hasSeeds {
			has_snp_2 = S.FindSNPCall(read2, s_pos, e_pos, match_pos, match_num, bw_snp_idx, fw_snp_idx, bw_snp_val, fw_snp_val, bw_D, fw_D, bw_T, fw_T)
			if has_snp_2 {
				//fmt.Println(loop_num, "\tori2\t", string(read2))
				return true
			}
		}
		//Find SNPs for the reverse complement of the second end
        s_pos, e_pos, match_num, hasSeeds = INDEX.FindSeeds(rev_comp_read2, comp_read2, p, match_pos)
        if hasSeeds {
			has_snp_2 = S.FindSNPCall(rev_comp_read2, s_pos, e_pos, match_pos, match_num, bw_snp_idx, fw_snp_idx, bw_snp_val, fw_snp_val, bw_D, fw_D, bw_T, fw_T)
			if has_snp_2 {
				//fmt.Println(loop_num, "\trev2\t", string(rev_read2))
				return true
			}
		}
        //Take a random position to search
		p = RAND_GEN.Intn(len(read2) - 1) + 1
        loop_num++
    }

    if has_snp_1 {
        return true
    } else {
        return false
    }
}

//---------------------------------------------------------------------------------------------------
// UpdateSNPCall updates SNP Call found from alignment between reads and multi-genomes.
//---------------------------------------------------------------------------------------------------
func (S *SNPProf) FindSNPCall(read []byte, s_pos, e_pos int, match_pos []int, match_num int, bw_snp_idx, fw_snp_idx []int, bw_snp_val, fw_snp_val [][]byte, bw_D, fw_D [][]int, bw_T, fw_T [][][]byte) bool {

    var k, idx, left_num, right_num, pos int
	var val []byte
	var isExtended, has_snp bool

	has_snp = false
    for i := 0 ; i < match_num ; i++ {
		pos = match_pos[i]
        //Call IntervalHasSNP to determine whether extension is needed
        //if index.IntervalHasSNP(A.SORTED_SNP_POS, pos - e_pos, pos - e_pos + len(read1)) {
        _, left_num, right_num, isExtended = INDEX.FindExtensions(read, s_pos, e_pos, pos, bw_snp_idx, fw_snp_idx, bw_snp_val, fw_snp_val, bw_D, fw_D, bw_T, fw_T)
        if isExtended {
            //Determine SNP profile
			if left_num > 0 {
				has_snp = true
	            for k = 0; k < left_num ; k++ {
    	            idx, val = bw_snp_idx[k], bw_snp_val[k]
					S.SNP_Prof[idx] = append(S.SNP_Prof[idx], val)
    	        }
			}
			if right_num > 0 {
				has_snp = true
				for k = 0; k < right_num ; k++ {
    	            idx, val = fw_snp_idx[k], fw_snp_val[k]
					S.SNP_Prof[idx] = append(S.SNP_Prof[idx], val)
    	        }
			}
		}
		//}
	}
	return has_snp
}

//---------------------------------------------------------------------------------------------------
// UpdateSNPProfile updates SNP profile found from alignment between reads and multi-genomes.
//---------------------------------------------------------------------------------------------------
func (S *SNPProf) UpdateSNPProfile(read1, read2, rev_read1, rev_read2, rev_comp_read1, rev_comp_read2, comp_read1, comp_read2 []byte, bw_snp_idx, fw_snp_idx []int, bw_snp_val, fw_snp_val [][]byte, bw_D, fw_D [][]int, bw_T, fw_T [][][]byte, match_pos []int) bool {

	has_SNP_call := S.UpdateSNPCall(read1, read2, rev_read1, rev_read2, rev_comp_read1, rev_comp_read2, comp_read1, comp_read2, bw_snp_idx, fw_snp_idx, bw_snp_val, fw_snp_val, bw_D, fw_D, bw_T, fw_T, match_pos)

	if has_SNP_call {
		return true
	} else {
		return false
	}
}


/*
func (S *SNPProf) CalculateQual() {

	    for i := 0 ; i < len(bw_snp_idx) := range snp_prof {
	        S.SNP_Prof[snp_pos] = append(S.SNP_Prof[snp_pos], snp_arr...)

			//Update SNP Qual - testing/////////////////////////////////////////////
			cond_prob := make([]float32, len(S.SNP_Conf[snp_pos]))
			for idx, base_conf := range S.SNP_Conf[snp_pos] {
				SNP := index.SNP_PROFILE[snp_pos][idx]
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
			///////////////////////////////////////////////////////////////////////

	    }
}
*/

//-----------------------------------------------------------------------------------------------------
// GenerateSNP returns called SNPs and related information based on SNP profile constructed from
// alignment between reads and multi-genomes.
//-----------------------------------------------------------------------------------------------------
func (S *SNPProf) CallSNP() {
    var snp []byte
    var snp_pos int
    var snp_prof [][]byte
    var major_snp string
    var major_num int
    var SNP_Qlt map[string]int	

    for snp_pos, snp_prof = range S.SNP_Prof {
        SNP_Qlt = make(map[string]int)
        for _, snp = range snp_prof {
            SNP_Qlt[string(snp)] = SNP_Qlt[string(snp)] + 1
        }
        major_num = 0
        for snp_val, snp_num := range SNP_Qlt {
            if snp_num > major_num {
                major_num = snp_num
                major_snp = snp_val
            }
        }
        S.SNP_Call[snp_pos] = []byte(major_snp)
        S.SNP_Prob[snp_pos] = []int{major_num, len(snp_prof)}
    }
}

//-------------------------------------------------------------------------------------------------------
// SNPCall_tofile writes called SNPs and related information to given output file in tab-delimited format
//-------------------------------------------------------------------------------------------------------
func (S *SNPProf) SNPCall_tofile(file_name string) {
    var snp_pos int
    var SNP_Pos = make([]int, 0, len(S.SNP_Call))
    for snp_pos, _ := range S.SNP_Call {
        SNP_Pos = append(SNP_Pos, snp_pos)
    }
    sort.Ints(SNP_Pos)

    file, err := os.Create(file_name)
    if  err != nil {
        return
    }
    defer file.Close()
    var num []int
    var str_snp_pos, str_snp, str_snp_num1, str_snp_num2, str_snp_prob string

	var flag bool
	var idx int
	var value []byte

    for _, snp_pos = range SNP_Pos {
        str_snp_pos = strconv.Itoa(snp_pos)
        str_snp = string(S.SNP_Call[snp_pos])
        num = S.SNP_Prob[snp_pos]
        str_snp_num1, str_snp_num2 = strconv.Itoa(num[0]), strconv.Itoa(num[1])
        str_snp_prob = strconv.FormatFloat(float64(num[0])/float64(num[1]), 'f', 5, 32)
        //fmt.Println(snp_pos, "\t", str_snp)
        _, err := file.WriteString(str_snp_pos + "\t" + str_snp + "\t" + 
					str_snp_num1 + "\t" + str_snp_num2 + "\t" + str_snp_prob + "\t");

		//Write SNP Qual - testing////////////////////////////
		flag = false
		for idx, value = range INDEX.SNP_PROFILE[snp_pos] {
			if bytes.Equal(value, S.SNP_Call[snp_pos]) {
		        str_snp_conf := strconv.FormatFloat(float64(S.SNP_Conf[snp_pos][idx]), 'f', 5, 32)
                _, err = file.WriteString(str_snp_conf + "\n");
				flag = true
			}
		}
		if !flag {
			_, err = file.WriteString("." + "\n");
			fmt.Println("--------")
			fmt.Println(snp_pos, " : ", string(S.SNP_Call[snp_pos]))
			for _, val := range INDEX.SNP_PROFILE[snp_pos] {
				fmt.Print(string(val), " - ")
			}
			fmt.Println()
		}
		//////////////////////////////////////////////////////

        if err != nil {
            fmt.Println(err)
            break
        }
    }   
}
