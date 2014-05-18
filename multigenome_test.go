//----------------------------------------------------------------------------------------
// Test for building multigenome
// Copyright 2014 Nam Sy Vo
//----------------------------------------------------------------------------------------

package multigenome

import (
	"fmt"
	"runtime"
	"strings"
	"testing"
)

//Functions for displaying information
func o_() string {
    pc, _, _, _ := runtime.Caller(1)
    name := runtime.FuncForPC(pc).Name()
    if p := strings.LastIndexAny(name, `./\`); p >= 0 {
        name = name[p+1:]
    } // if
    fmt.Println("== BEGIN", name, "===")
    return name
}

func __(name string) {
    fmt.Println("== END", name, "===")
    fmt.Println()
}

func TestMultiGenomeBuild(t *testing.T) {
    defer __(o_())

	sequence := fastaRead("test_data/chr1.fasta")
	SNP_array := vcfRead("test_data/vcf_chr_1.vcf")
	genome := buildMultigenome2(SNP_array, sequence)

	SaveMulti("test_data/genomestar.txt", genome)
	SaveSNPLocation("test_data/SNPLocation.txt", SNP_array)

	saved_genome := LoadMulti("test_data/genomestar.txt")
	saved_SNP_array, saved_SameLen_SNP := LoadSNPLocation("test_data/SNPLocation.txt")

	fmt.Println(len(saved_genome))
	fmt.Println(len(saved_SNP_array))
	fmt.Println(len(saved_SameLen_SNP))
}
