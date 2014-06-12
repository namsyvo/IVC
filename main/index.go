//----------------------------------------------------------------------------------------
// ISC - main program
// Copyright 2014 Nam Sy Vo
//----------------------------------------------------------------------------------------

package main

import (
    "fmt"
    "flag"
    "github.com/namsyvo/ISC"
    "github.com/vtphan/fmi"
)

func main() {
    var genome_file = flag.String("g", "", "reference genome file")
    var dbsnp_file = flag.String("p", "", "snp profile file")
    var snp_file = flag.String("s", "", "index file of snp profile")
    var multigenome_file = flag.String("m", "", "multigenome file")
    var index_file = flag.String("i", "", "index file of multigenome")
    var rev_index_file = flag.String("r", "", "index file of reverse of multigenome")
    //var workers = flag.Int("w", 1, "number of workers")
    //flag.BoolVar(&Debug, "debug", false, "Turn on debug mode.")
    flag.Parse()

	fmt.Println("ISC - Integrated SNP Calling based on Read-Multigenome Alignment")

    fmt.Println("Creating multigenome and SNP profile index...")

	sequence := isc.ReadFASTA(*genome_file)
	SNP_array := isc.ReadVCF(*dbsnp_file)
	multigenome := isc.BuildMultigenome(SNP_array, sequence)

	multigenome_len := len(multigenome)
	rev_multigenome := make([]byte, multigenome_len)
	for i := range multigenome {
		rev_multigenome[i] = multigenome[multigenome_len-1-i]
	}
	isc.SaveMultigenome(*multigenome_file, multigenome)
	rev_multigenome_file := *multigenome_file + ".rev"
	isc.SaveMultigenome(rev_multigenome_file, rev_multigenome)

	isc.SaveSNPLocation(*snp_file, SNP_array)

	saved_genome := isc.LoadMultigenome(*multigenome_file)
	saved_SNP_array, saved_SameLen_SNP := isc.LoadSNPLocation(*snp_file)

	fmt.Println("Multigenome length: ", len(saved_genome))
	fmt.Println("SNP profile index size: ", len(saved_SNP_array))
	fmt.Println("Same length SNP profile index size: ", len(saved_SameLen_SNP))

    fmt.Println("Finish creating multigenome and SNP profile index...")

    fmt.Println("Indexing multigenome...")
	
    var idx fmi.Index
	idx = *fmi.New(*multigenome_file)
	idx.Save(*index_file)
	idx = *fmi.New(rev_multigenome_file)
	idx.Save(*rev_index_file)

    fmt.Println("Finish indexing multigenome...")
}
