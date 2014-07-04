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
)

//Global variables
var (
    index Index
    start_pos int = 0
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

    S.SNP_Prof = make(map[int][][]byte)
    S.SNP_Conf = make(map[int][]float32)
    S.SNP_Call = make(map[int][]byte)
    S.SNP_Prob = make(map[int][]int)

    index.Init(genome_file, snp_file, index_file, rev_index_file, read_len, seq_err, k, a, n)
	S.SNP_Conf = index.SNP_AF

    //Print out Allele Freq at each SNP - testing/////////////////
    /*
    for snp_pos, _ := range(index.SNP_PROFILE) {
            for idx, value := range(index.SNP_PROFILE[snp_pos]) {
                fmt.Println(snp_pos, " : ", string(value), " : ", S.SNP_Conf[snp_pos][idx])
            }
            fmt.Println()
    }
    */
    //testing////////////////////////////////////////////////////
}

//--------------------------------------------------------------------------------------------------
// FindSNPProfile returns SNP profile of new genome based on SNP profile of reference multi-genomes
// and alignment between reads and multi-genomes.
//--------------------------------------------------------------------------------------------------
func (S SNPProf) FindSNPProfile(read1, read2 []byte) (map[int][][]byte, bool) {

    snp_profile := make(map[int][][]byte)
    var k int
    var v []byte

    //Consider reverse complement of reads.
    rev_read1 := ReverseComplement(read1)
    rev_read2 := ReverseComplement(read2)
    /*
    fmt.Println("read1", string(read1))
    fmt.Println("rev_read1", string(rev_read1))
    fmt.Println("read2", string(read2))
    fmt.Println("rev_read2", string(rev_read2))
    */

    //Find SNPs for pairend reads, treat each end separately and independently.
    s_pos, e_pos, match_pos, hasExactMatches := -1, -1, []int{}, false
    r := rand.New(rand.NewSource(time.Now().UnixNano()))

    //Find SNPs for the first end
    var p int = start_pos
    loop_num := 1
    for loop_num <= ITER_NUM {
        //fmt.Println(loop_num, "read1", string(read1))
        //Call FindSeeds to determine seed
        s_pos, e_pos, match_pos, hasExactMatches = index.FindSeeds(read1, p)
        if hasExactMatches {
            //fmt.Println(s_pos, e_pos, match_pos)
            for _, pos := range match_pos {
                //Call IntervalHasSNP to determine whether extension is needed
                //if index.IntervalHasSNP(index.SORTED_SNP_POS, pos - e_pos, pos - e_pos + len(read1)) {
                    //Call ApproxSearch to determine extension
                    _, _, _, _, _, left_snp, right_snp, isExtended := index.FindExtensions(read1, s_pos, e_pos, pos)
                    if isExtended {
                        //Determine SNP profile
                        for k, v = range left_snp {
                            snp_profile[k] = append(snp_profile[k], v)
                        }
                        for k, v = range right_snp {
                            snp_profile[k] = append(snp_profile[k], v)
                        }
                    }
                //}
            }
            if len(snp_profile) > 0 {
                break
            }
        }
        //Find SNPs for the reverse complement of second end
        //fmt.Println(loop_num, "rev_read1", string(rev_read1))
        //Call FindSeeds to determine seed
        s_pos, e_pos, match_pos, hasExactMatches = index.FindSeeds(rev_read1, p)
        if hasExactMatches {
            //fmt.Println(s_pos, e_pos, match_pos)
            for _, pos := range match_pos {
                //Call IntervalHasSNP to determine whether extension is needed
                //if index.IntervalHasSNP(A.SORTED_SNP_POS, pos - e_pos, pos - e_pos + len(read1)) {
                    //Call ApproxSearch to determine extension
                    _, _, _, _, _, left_snp, right_snp, isExtended := index.FindExtensions(rev_read1, s_pos, e_pos, pos)
                    if isExtended {
                        //Determine SNP profile
                        for k, v = range left_snp {
                            snp_profile[k] = append(snp_profile[k], v)
                        }
                        for k, v = range right_snp {
                            snp_profile[k] = append(snp_profile[k], v)
                        }
                    }
                //}
            }
            if len(snp_profile) > 0 {
                break
            }
        }
        //Take a random position to search
        p = r.Intn(len(read1) - 1) + 1
        loop_num++
    }

    snp_prof_len_1 := len(snp_profile)

    //Find SNPs for the second end
    loop_num = 1
    snp_found_num := 0
    p = start_pos
    for loop_num <= ITER_NUM {
        //fmt.Println(loop_num, "read2", string(read2))
        //Call FindSeeds to determine seed
        s_pos, e_pos, match_pos, hasExactMatches = index.FindSeeds(read2, p)
        if hasExactMatches {
            //fmt.Println(s_pos, e_pos, match_pos)
            for _, pos := range match_pos {
                //Call IntervalHasSNP to determine whether extension is needed
                //if index.IntervalHasSNP(A.SORTED_SNP_POS, pos - e_pos, pos - e_pos + len(read1)) {
                    //Call ApproxSearch to determine extension
                    _, _, _, _, _, left_snp, right_snp, isExtended := index.FindExtensions(read2, s_pos, e_pos, pos)
                    if isExtended {
                        snp_found_num += 1
                        //Determine SNP profile
                        for k, v = range left_snp {
                            snp_profile[k] = append(snp_profile[k], v)
                        }
                        for k, v = range right_snp {
                            snp_profile[k] = append(snp_profile[k], v)
                        }
                    }
                //}
            }
            if snp_found_num > 0 {
                return snp_profile, true
            }
        }
        //Find SNPs for the reverse complement of second end
        //fmt.Println(loop_num, "rev_read2", string(rev_read2))
        //Call FindSeeds to determine seed
        s_pos, e_pos, match_pos, hasExactMatches = index.FindSeeds(rev_read2, p)
        if hasExactMatches {
            //fmt.Println(s_pos, e_pos, match_pos)
            for _, pos := range match_pos {
                //Call IntervalHasSNP to determine whether extension is needed
                //if index.IntervalHasSNP(A.SORTED_SNP_POS, pos - e_pos, pos - e_pos + len(read1)) {
                    //Call ApproxSearch to determine extension
                    _, _, _, _, _, left_snp, right_snp, isExtended := index.FindExtensions(rev_read2, s_pos, e_pos, pos)
                    if isExtended {
                        snp_found_num += 1
                        //Determine SNP profile
                        for k, v = range left_snp {
                            snp_profile[k] = append(snp_profile[k], v)
                        }
                        for k, v = range right_snp {
                            snp_profile[k] = append(snp_profile[k], v)
                        }
                    }
                //}
            }
            if snp_found_num > 0 {
                return snp_profile, true
            }
        }
        //Take a random position to search
        p = r.Intn(len(read2) - 1) + 1
        loop_num++
    }

    if snp_prof_len_1 > 0 {
        return snp_profile, true
    } else {
        return snp_profile, false
    }
}

//---------------------------------------------------------------------------------------------------
// UpdateSNPProfile updates SNP profile found from alignment between reads and multi-genomes.
//---------------------------------------------------------------------------------------------------
func (S *SNPProf) UpdateSNPProfile(read1, read2 []byte) bool {
	snp_prof, has_SNP_call := S.FindSNPProfile(read1, read2)
	if has_SNP_call {
	    for snp_pos, snp_arr := range snp_prof {
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
            /*
            fmt.Println("Base: ", string(snp_arr[0]))
            for idx, value := range(index.SNP_PROFILE[snp_pos]) {
                fmt.Println(snp_pos, " : ", string(value), " : ", S.SNP_Conf[snp_pos][idx])
            }
            fmt.Println()
            */
			/////////////////////////////////////////////////////////////////////////
	    }
		return true
	} else {
		return false
	}
}

//-----------------------------------------------------------------------------------------------------
// GenerateSNP returns called SNPs and related information based on SNP profile constructed from
// alignment between reads and multi-genomes.
//-----------------------------------------------------------------------------------------------------
func (S SNPProf) CallSNP() {
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
func (S SNPProf) SNPCall_tofile(file_name string) {
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
		for idx, value := range index.SNP_PROFILE[snp_pos] {
			if bytes.Equal(value, S.SNP_Call[snp_pos]) {
		        str_snp_conf := strconv.FormatFloat(float64(S.SNP_Conf[snp_pos][idx]), 'f', 5, 32)
                _, err = file.WriteString(str_snp_conf + "\n");
			}
		}
		//////////////////////////////////////////////////////

        if err != nil {
            fmt.Println(err)
            break
        }
    }   
}
