//----------------------------------------------------------------------------------------
// Copyright 2013 Nam S. Vo
// Test for randalx
//----------------------------------------------------------------------------------------

package randalx_test

import (
	"randalx"
	"fmt"
	"os"
    "bufio"
	"runtime"
	"strings"
	"testing"
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

type TestCase struct {
	A []int
	i, j int
	result bool
}

func TestApproxSearch_4M(t *testing.T) {
	defer __(o_())

    fmt.Println("@@@-Read-MultiGenome Alignment ----------------------")

    var genome_file = "../data-test-4M/genome_starred.txt"
    var snp_file = "../data-test-4M/SNPlocation.txt"
    var index_file = "../data-test-4M/genome_starred.txt.index"
    var rev_index_file = "../data-test-4M/genome_starred_rev.txt.index"

    var search randalx.Search
    search.Init_seq(genome_file, snp_file, index_file, rev_index_file)
    search.Init_para(0.02, 4, 1, 100)

    var snpcaller randalx.SNPCaller
    snpcaller.Init()

    var queries_file = "../data-test-4M/reads-test1.txt"
    var read []byte
	var snp_profile map[int][][]byte
	var isSNPCalled bool
	var read_num int = 0
	var snp_aligned_read_num int = 0

    f, err := os.Open(queries_file)
    if err != nil {
        panic("error opening file " + queries_file)
    }
    data := bufio.NewReader(f)
    fmt.Println("@@@-Aligning reads to mutigenome -------------------")
    for {
        line, err := data.ReadBytes('\n')
        if err != nil {
            break
        }
        if len(line) > 1 {
        	read_num++
            read = line[0:100]
            //pos_mate = line[30:len(line) - 1]
            //fmt.Println("\nread to search, pos: ", string(read), string(pos_mate))
            snp_profile, isSNPCalled = search.FindSNPProfile_read(read)
            if isSNPCalled {
            	snpcaller.UpdateSNPProfile(snp_profile)
            	snp_aligned_read_num++
            }
        }
    }
    fmt.Println("read_num: ", read_num)
    fmt.Println("snp_aligned_read_num: ", snp_aligned_read_num)

    fmt.Println("@@@-Construct SNP profiles -------------------------")
    SNP_Calling, SNP_Prob := snpcaller.CallSNP()
    fmt.Println("@@@-SNP Calling ------------------------------------")
    snpcaller.SNPCall_tofile(SNP_Calling, SNP_Prob, "../data-test-4M/snp_calling_test.txt")
    fmt.Println("@@@-Finish ----------------------------------------")
}