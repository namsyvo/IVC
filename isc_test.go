//----------------------------------------------------------------------------------------
// Copyright 2013 Nam Sy Vo
// Test for ISC
//----------------------------------------------------------------------------------------

package isc

import (
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

func TestISC_NC007194(t *testing.T) {
	defer __(o_())

    fmt.Println("@@@-Read-MultiGenome Alignment ----------------------")

    var genome_file = "test_data/indexes/genome_starred.txt"
    var snp_file = "test_data/indexes/SNPlocation.txt"
    var index_file = "test_data/indexes/genome_starred.txt.index"
    var rev_index_file = "test_data/indexes/genome_starred_rev.txt.index"

    var snpcaller SNPCaller
    snpcaller.Init(genome_file, snp_file, index_file, rev_index_file, 100, 0.02, 4, 1, 32)

    var queries_file = "test_data/reads/test_reads.txt"
    var read []byte
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
			isSNPCalled = snpcaller.UpdateSNPProfile(read)
			if isSNPCalled {
	       	    snp_aligned_read_num++
	        }
        }
    }
    fmt.Println("read_num: ", read_num)
    fmt.Println("snp_aligned_read_num: ", snp_aligned_read_num)

    fmt.Println("@@@-SNP Calling ------------------------------------")
    snpcaller.CallSNP()
    snpcaller.SNPCall_tofile("test_data/results/test_called_snps.txt")
    fmt.Println("@@@-Finish ----------------------------------------")
}
