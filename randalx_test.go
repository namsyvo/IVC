//----------------------------------------------------------------------------------------
// Copyright 2013 Nam S. Vo
// Test for Approximate searching package
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

/*
func TestIntervalHasSNP(t *testing.T) {
	defer __(o_())

	var test_cases = []TestCase{
		{[]int{1, 4, 7, 15, 20}, 3, 5, true},
		{[]int{1, 4, 7, 15, 20}, 5, 6, false},
		{[]int{1, 4, 7, 15, 20}, 1, 2, true},
		{[]int{1, 4, 7, 15, 20}, 2, 4, true},
		{[]int{1, 4, 7, 15, 20}, 16, 19, false},
		{[]int{1, 4, 7, 15, 20}, 1, 20, true},
		{[]int{1, 4, 7, 15, 20}, 20, 20, true},
		{[]int{1, 4, 7, 15, 20}, 2, 19, true},
		{[]int{1, 4, 7, 15, 20}, 3, 19, true},
		{[]int{1, 4, 7, 15, 20}, 4, 19, true},
		{[]int{1, 4, 7, 15, 20}, 5, 19, true},
		{[]int{1, 4, 7, 15, 20}, 6, 19, true},
		{[]int{1, 4, 7, 15, 20}, 7, 19, true},
		{[]int{1, 4, 7, 15, 20}, 8, 19, true},
		{[]int{1, 4, 7, 15, 20}, 9, 19, true},
		{[]int{1, 4, 7, 15, 20}, 10, 19, true},
		{[]int{1, 4, 7, 15, 20}, 11, 19, true},
		{[]int{1, 4, 7, 15, 20}, 13, 19, true},
		{[]int{1, 4, 7, 15, 20}, 14, 19, true},
		{[]int{1, 4, 7, 15, 20}, 15, 19, true},
		{[]int{1, 4, 7, 15, 20}, 16, 19, false},
		{[]int{1, 4, 7, 15, 20}, 17, 19, false},
		{[]int{1, 4, 7, 15, 20}, 18, 19, false},
		{[]int{1, 4, 7, 15, 20}, 19, 22, true},
		{[]int{1, 4, 7, 15, 20}, 22, 25, false},
		{[]int{1, 4, 7, 15, 20}, 0, 1, true},
		{[]int{3, 4, 7, 15, 20}, 1, 3, true},
		{[]int{3, 4, 7, 15, 20}, 1, 2, false},
		{[]int{3, 4, 7, 15, 20}, 5, 6, false},
		{[]int{3, 4, 7, 15, 20}, 16, 16, false},
		{[]int{3, 3}, 3, 3, true},
		{[]int{3, 4, 7, 15, 20}, 8, 10, false},
		{[]int{3, 4, 7, 15, 20}, 14, 15, true},
		{[]int{3, 4, 7, 15, 20}, 19, 20, true},
		{[]int{3, 4, 7, 15, 20}, 4, 5, true},
		{[]int{3, 4, 7, 15, 20}, 7, 8, true},
	}
	var search randalx.Search
	for k := 0; k < len(test_cases); k++ {
		a := search.IntervalHasSNP(test_cases[k].A, test_cases[k].i, test_cases[k].j)
		if a != test_cases[k].result {
			t.Errorf("Fail searching (A, i, j)", test_cases[k].A, test_cases[k].i, test_cases[k].j)
		} else {
			fmt.Println("Succesful searching (A, i, j)", test_cases[k].A, test_cases[k].i, test_cases[k].j, test_cases[k].result)	
		}
	}
}
*/
/*
func TestApproxSearch_200(t *testing.T) {
	defer __(o_())

    fmt.Println("@@@-Read-MultiGenome Alignment ----------------------")

    var genome_file = "../data-test-200/genome_starred.txt"
    var snp_file = "../data-test-200/SNPlocation.txt"
    var index_file = "../data-test-200/genome_starred.txt.index"
    var rev_index_file = "../data-test-200/genome_starred_rev.txt.index"
    var search randalx.Search
    search.Init_seq(genome_file, snp_file, index_file, rev_index_file)
    search.Init_para(0.02, 4, 1, 30)

    var snpcaller randalx.SNPCaller
    snpcaller.Init()

    var queries_file = "../data-test-200/reads.txt"
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
            read = line[0:30]
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
    SNP_Calling := snpcaller.CallSNP()
    fmt.Println("@@@-SNP Calling ------------------------------------")
    snpcaller.SNPCall_tofile(SNP_Calling, "../data-test-200/snp_calling.txt")
    fmt.Println("@@@-Finish ----------------------------------------")
}
*/

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

    var queries_file = "../data-test-4M/reads.txt"
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
    SNP_Calling := snpcaller.CallSNP()
    fmt.Println("@@@-SNP Calling ------------------------------------")
    snpcaller.SNPCall_tofile(SNP_Calling, "../data-test-4M/snp_calling.txt")
    fmt.Println("@@@-Finish ----------------------------------------")
}