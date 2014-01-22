//----------------------------------------------------------------------------------------
// Copyright 2013 Nam S. Vo
// SNP calling utility functions for read-multigenome alignment
//----------------------------------------------------------------------------------------

package randalx

import (
    "fmt"
    "os"
    "strconv"
    "sort"
)

// Search object with Parameters
type SNPCaller struct {
    SNP_Profile map[int][][]byte
}


func (S *SNPCaller) Init() {
    S.SNP_Profile = make(map[int][][]byte)
}
//-----------------------------------------------------------------------------------------------------
// UpdateSNPProfile updates SNP profile found from alignment between reads and multi-genomes.
//-----------------------------------------------------------------------------------------------------
func (S *SNPCaller) UpdateSNPProfile(snp_profile map[int][][]byte) {
    for snp_pos, snp_prof := range snp_profile {
        S.SNP_Profile[snp_pos] = append(S.SNP_Profile[snp_pos], snp_prof...)
    }
}

//-----------------------------------------------------------------------------------------------------
// GenerateSNP returns SNP profile of new genome from SNP profile of reference multi-genomes and
// alignment between reads and multi-genomes.
//-----------------------------------------------------------------------------------------------------
func (S SNPCaller) CallSNP() (map[int][]byte, map[int][]int) {

    var snp_pos int
    var snp []byte
    var snp_prof [][]byte
    var major_snp string
    var major_num int
    var SNP_Qlt map[string]int

    var SNP_Calling = make(map[int][]byte)
    var SNP_Prob = make(map[int][]int)

    for snp_pos, snp_prof = range S.SNP_Profile {
        SNP_Qlt = make(map[string]int)
        //fmt.Println("SNP Profiles at ", snp_pos, ": ")
        for _, snp = range snp_prof {
            SNP_Qlt[string(snp)] = SNP_Qlt[string(snp)] + 1
            //fmt.Print(string(snp), "\t")
        }
        //fmt.Println()
        major_num = 0
        for snp_val, snp_num := range SNP_Qlt {
            if snp_num > major_num {
                major_num = snp_num
                major_snp = snp_val
            }
        }
        SNP_Calling[snp_pos] = []byte(major_snp)
        SNP_Prob[snp_pos] = []int{major_num, len(snp_prof)}
        //fmt.Println("Percentage of majority snp: ", float32(major_num)/float32(len(snp_prof)))
    }
    return SNP_Calling, SNP_Prob
}

func (S SNPCaller) SNPCall_tofile(SNP_Calling map[int][]byte, SNP_Prob map[int][]int, file_name string) {
    var snp_pos int
    var SNP_Pos = make([]int, 0, len(SNP_Calling))
    for snp_pos, _ := range SNP_Calling {
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
        str_snp = string(SNP_Calling[snp_pos])
        num = SNP_Prob[snp_pos]
        str_snp_num1, str_snp_num2 = strconv.Itoa(num[0]), strconv.Itoa(num[1])
        str_snp_prob = strconv.FormatFloat(float64(num[0])/float64(num[1]), 'f', 5, 32)
        //fmt.Println(snp_pos, "\t", str_snp)
        _, err := file.WriteString(str_snp_pos + "\t" + str_snp + "\t" + str_snp_num1 + "\t" + str_snp_num2 + "\t" + str_snp_prob + "\n");
        if err != nil {
            fmt.Println(err)
            break
        }
    }   
}