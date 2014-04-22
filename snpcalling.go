//---------------------------------------------------------------------------------------------------
// Calling SNPs based on read-multigenome alignment results.
// Copyright 2014 Nam Sy Vo
//---------------------------------------------------------------------------------------------------

package isc

import (
    "fmt"
    "os"
    "strconv"
    "sort"
)

//Global variables
var (
    aligner Aligner
)

//SNP caller object with Parameters
type SNPCaller struct {
    SNP_Prof map[int][][]byte
    SNP_Call map[int][]byte
    SNP_Prob map[int][]int
}

//Initialize parameters
func (S *SNPCaller) Init(genome_file, snp_file, index_file, rev_index_file string, read_len int,
	seq_err float32, k, a, n int) {

    aligner.Init(genome_file, snp_file, index_file, rev_index_file, read_len, seq_err, k, a, n)
    S.SNP_Prof = make(map[int][][]byte)
    S.SNP_Call = make(map[int][]byte)
    S.SNP_Prob = make(map[int][]int)
}

//---------------------------------------------------------------------------------------------------
// UpdateSNPProfile updates SNP profile found from alignment between reads and multi-genomes.
//---------------------------------------------------------------------------------------------------
func (S *SNPCaller) UpdateSNPProfile(read []byte) bool {
	snp_prof, has_SNP_call := aligner.FindSNPProfile_read(read)
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
func (S SNPCaller) CallSNP() (map[int][]byte, map[int][]int) {
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
func (S SNPCaller) SNPCall_tofile(file_name string) {
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
