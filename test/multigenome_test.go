//----------------------------------------------------------------------------------------
// Test for building multigenome
// Copyright 2014 Nam Sy Vo
//----------------------------------------------------------------------------------------

package isc_test

import (
	"fmt"
	"testing"
    "github.com/namsyvo/ISC"
)

func TestMultiGenomeBuild(t *testing.T) {
    defer __(o_())

	var index isc.Index
	//index.Init(genome_file, snp_file, index_file, rev_index_file, read_len, seq_err, k, a, n)
	sequence := isc.ReadFASTA("../test_data/refs/chr1.fasta")
	SNP_array := isc.ReadVCF("../test_data/refs/vcf_chr_1.vcf")
	genome := isc.BuildMultigenome(SNP_array, sequence)

	index.SaveMultigenome("../test_data/refs/genomestar.txt", genome)
	index.SaveSNPLocation("../test_data/refs/SNPLocation.txt", SNP_array)

	saved_genome := index.LoadMultigenome("../test_data/refs/genomestar.txt")
	saved_SNP_array, saved_SameLen_SNP := index.LoadSNPLocation("../test_data/refs/SNPLocation.txt")

	fmt.Println(len(saved_genome))
	fmt.Println(len(saved_SNP_array))
	fmt.Println(len(saved_SameLen_SNP))
}
