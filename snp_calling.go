//----------------------------------------------------------------------------------------
// Copyright 2013 Nam S. Vo
// Approximate searching of reads on multigenomes based on FM index and Edit distance
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
func (S SNPCaller) CallSNP() (map[int][]byte) {

    var snp_pos int
    var snp []byte
    var snp_prof [][]byte
    var major_snp string
    var major_num int
    var SNP_Qlt map[string]int

    var SNP_Calling = make(map[int][]byte)

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
        //fmt.Println("Percentage of majority snp: ", float32(major_num)/float32(len(snp_prof)))
    }
    return SNP_Calling
}

func (S SNPCaller) SNPCall_tofile(SNP_Calling map[int][]byte, file_name string) {
    var snp_pos int
    var snp []byte
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
    var str_snp_pos, str_snp string
    for _, snp_pos = range SNP_Pos {
        snp = SNP_Calling[snp_pos]
        //fmt.Println(snp_pos, "\t", string(snp))
        str_snp_pos = strconv.Itoa(snp_pos)
        str_snp = string(snp)
        _, err := file.WriteString(str_snp_pos + "\t" + str_snp + "\n"); 
        if err != nil {
            fmt.Println(err)
            break
        }
    }   
}