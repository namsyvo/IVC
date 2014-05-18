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
    "sort"
)

//Global variables
var (
    index Index
)

//SNP caller object with Parameters
type SNPProf struct {
    SNP_Prof map[int][][]byte
    SNP_Call map[int][]byte
    SNP_Prob map[int][]int
}

//Initialize parameters
func (S *SNPProf) Init(genome_file, snp_file, index_file, rev_index_file string, read_len int,
	seq_err float32, k, a, n int) {

    index.Init(genome_file, snp_file, index_file, rev_index_file, read_len, seq_err, k, a, n)
    S.SNP_Prof = make(map[int][][]byte)
    S.SNP_Call = make(map[int][]byte)
    S.SNP_Prob = make(map[int][]int)
}

//--------------------------------------------------------------------------------------------------
// FindSNPProfile returns SNP profile of new genome based on SNP profile of reference multi-genomes
// and alignment between reads and multi-genomes.
//--------------------------------------------------------------------------------------------------
func (I *Index) FindSNPProfile(read1, read2 []byte) (map[int][][]byte, bool) {

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
    var p int
    s_pos, e_pos, match_pos, hasExactMatches := -1, -1, []int{}, false
    loop_num := 1
    r := rand.New(rand.NewSource(time.Now().UnixNano()))

    //Find SNPs for the first end
    for loop_num <= ITER_NUM {
        //fmt.Println(loop_num, "read1", string(read1))
        p = r.Intn(len(read1) - 1) + 1
        //Call FindLCS to determine seed
        s_pos, e_pos, match_pos, hasExactMatches = index.FindLCS(read1, p)
        if hasExactMatches {
            for _, pos := range match_pos {
                //Call IntervalHasSNP to determine whether extension is needed
                //if index.IntervalHasSNP(index.SORTED_SNP_POS, pos - e_pos, pos - e_pos + len(read1)) {
                    //Call ApproxSearch to determine extension
                    _, _, _, _, _, left_snp, right_snp, isExtended := index.FindExtension(read1, s_pos, e_pos, pos)
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
        loop_num++
    }

    if (len(snp_profile) == 0) {
        //Find SNPs for the reverse complement of first end
        for loop_num <= ITER_NUM {
            //fmt.Println(loop_num, "rev_read1", string(rev_read1))
            p = r.Intn(len(rev_read1) - 1) + 1
            //Call FindLCS to determine seed
            s_pos, e_pos, match_pos, hasExactMatches = index.FindLCS(rev_read1, p)
            if hasExactMatches {
                for _, pos := range match_pos {
                    //Call IntervalHasSNP to determine whether extension is needed
                    //if index.IntervalHasSNP(A.SORTED_SNP_POS, pos - e_pos, pos - e_pos + len(read1)) {
                        //Call ApproxSearch to determine extension
                        _, _, _, _, _, left_snp, right_snp, isExtended := index.FindExtension(rev_read1, s_pos, e_pos, pos)
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
            loop_num++
        }
    }

    snp_prof_len_1 := len(snp_profile)

    //Find SNPs for the second end
    loop_num = 1
    snp_found_num := 0
    for loop_num <= ITER_NUM {
        //fmt.Println(loop_num, "read2", string(read2))
        p = r.Intn(len(read2) - 1) + 1
        //Call FindLCS to determine seed
        s_pos, e_pos, match_pos, hasExactMatches = index.FindLCS(read2, p)
        if hasExactMatches {
            for _, pos := range match_pos {
                //Call IntervalHasSNP to determine whether extension is needed
                //if index.IntervalHasSNP(A.SORTED_SNP_POS, pos - e_pos, pos - e_pos + len(read1)) {
                    //Call ApproxSearch to determine extension
                    _, _, _, _, _, left_snp, right_snp, isExtended := index.FindExtension(read2, s_pos, e_pos, pos)
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
        loop_num++
    }

    //Find SNPs for the reverse complement of second end
    loop_num = 1
    snp_found_num = 0
    for loop_num <= ITER_NUM {
        //fmt.Println(loop_num, "rev_read2", string(rev_read2))
        p = r.Intn(len(rev_read2) - 1) + 1
        //Call FindLCS to determine seed
        s_pos, e_pos, match_pos, hasExactMatches = index.FindLCS(rev_read2, p)
        if hasExactMatches {
            for _, pos := range match_pos {
                //Call IntervalHasSNP to determine whether extension is needed
                //if index.IntervalHasSNP(A.SORTED_SNP_POS, pos - e_pos, pos - e_pos + len(read1)) {
                    //Call ApproxSearch to determine extension
                    _, _, _, _, _, left_snp, right_snp, isExtended := index.FindExtension(rev_read2, s_pos, e_pos, pos)
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
	snp_prof, has_SNP_call := index.FindSNPProfile(read1, read2)
	if has_SNP_call {
	    for snp_pos, snp_prof := range snp_prof {
	        S.SNP_Prof[snp_pos] = append(S.SNP_Prof[snp_pos], snp_prof...)
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
func (S SNPProf) CallSNP() (map[int][]byte, map[int][]int) {
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
    return S.SNP_Call, S.SNP_Prob
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
					str_snp_num1 + "\t" + str_snp_num2 + "\t" + str_snp_prob + "\n");
        if err != nil {
            fmt.Println(err)
            break
        }
    }   
}
