//----------------------------------------------------------------------------------------
// Test for IVC main program
// Copyright 2015 Nam Sy Vo
//----------------------------------------------------------------------------------------

package ivc_test

import (
	"fmt"
    "runtime"
    "strings"
	//"os"
	//"bufio"
	//"testing"
    //"github.com/namsyvo/IVC"
)

func o_() string {
	pc, _, _, _ := runtime.Caller(1)
	name := runtime.FuncForPC(pc).Name()
	if p := strings.LastIndexAny(name, `./\`); p >= 0 {
		name = name[p+1:]
	} // if
	fmt.Println()
	fmt.Println("== BEGIN", name, "===")
	return name
}

func __(name string) {
	fmt.Println("== END", name, "===")
	fmt.Println()
}
/*
func TestIVC_NC007194(t *testing.T) {
	defer __(o_())

    fmt.Println("@@@-Read-MultiGenome Alignment ----------------------")

    var genome_file = "../test_data/indexes/multigenome.txt"
    var snp_file = "../test_data/indexes/SNPLocation.txt"
    var index_file = "../test_data/indexes/multigenome.txt.index"
    var rev_index_file = "../test_data/indexes/multigenome_rev.txt.index"

    var snpcaller ivc.SNPProf
    snpcaller.Init(genome_file, snp_file, index_file, rev_index_file, 100, 0.02, 4, 1, 32)

    var q_file_1 = "../test_data/reads/test_reads_1.fq"
    var q_file_2 = "../test_data/reads/test_reads_2.fq"
    var read1, read2 []byte
	var read_num int = 0
	var snp_aligned_read_num int = 0
    var read_len int = 100
    var has_SNP_call bool

    f1, err1 := os.Open(q_file_1)
    if err1 != nil {
        panic("Error opening file " + q_file_1)
    }
    f2, err2 := os.Open(q_file_2)
    if err2 != nil {
        panic("Error opening file " + q_file_2)
    }
    data1 := bufio.NewReader(f1)
    data2 := bufio.NewReader(f2)
    var line_f1, line_f2 []byte
    if q_file_1[len(q_file_1)-3:] == ".fq" || q_file_1[len(q_file_1)-6:] == ".fastq"  {
        for {
            data1.ReadBytes('\n') //ignore 1st line in input FASTQ file
            data2.ReadBytes('\n') //ignore 1st line in input FASTQ file
            line_f1, err1 = data1.ReadBytes('\n')
            line_f2, err2 = data2.ReadBytes('\n')
            if err1 != nil {
                break
            }
            if err2 != nil {
                break
            }
            if len(line_f1) >= read_len && len(line_f2) >= read_len{
                read_num++
                read1 = line_f1[0 : read_len]
                read2 = line_f2[0 : read_len]
                has_SNP_call = snpcaller.UpdateSNPProfile(read1, read2)
                if has_SNP_call {
                    snp_aligned_read_num++
                }
            }
            data1.ReadBytes('\n') //ignore 3rd line in input FASTQ file
            data1.ReadBytes('\n') //ignore 4th line in input FASTQ file
            data2.ReadBytes('\n') //ignore 3rd line in input FASTQ file
            data2.ReadBytes('\n') //ignore 4th line in input FASTQ file
        }
    }
    fmt.Println("read_num: ", read_num)
    fmt.Println("snp_aligned_read_num: ", snp_aligned_read_num)

    fmt.Println("@@@-SNP Calling ------------------------------------")
    snpcaller.CallSNP()
    snpcaller.SNPCall_tofile("../test_data/results/test_called_snps.txt")
    fmt.Println("@@@-Finish ----------------------------------------")
}

func TestMem(t *testing.T) {
    var a, b, c, d []byte
    memstats := new(runtime.MemStats)
    a = make([]byte, 100)
    for i:= 0 ; i < 100 ; i++ {
        a[i] = 'A'
    }
    for i:= 0 ; i < 1000000 ; i++ {
        b, c, d = MakeMem(a)
        if len(b) == 0 && len(c) == 0 && len(d) == 0 {
            println("Error!")
        }
        if i%10000 == 0 {
            runtime.ReadMemStats(memstats)
            fmt.Println(memstats.Alloc)
        }
    }
}

func MakeMem(a []byte) ([]byte, []byte, []byte) {
    b, c, d := make([]byte, 100), make([]byte, 100), make([]byte, 100)
    for i, v := range a {
        b[i], c[i], d[i] = v, v, v
    }
    return b, c, d
}
*/