//----------------------------------------------------------------------------------------
// ISC - main program
// Copyright 2014 Nam Sy Vo
//----------------------------------------------------------------------------------------

package main

import (
	"flag"
	"fmt"
	"time"
	"runtime"
	"log"
	"path"
	"github.com/namsyvo/ISC"
	"github.com/vtphan/fmi"
)

func main() {
	fmt.Println("ISC - Integrated SNP Calling based on Read-Multigenome Alignment")
	fmt.Println("ISC-index: Indexing reference genomes and variant profiles.")

	memstats := new(runtime.MemStats)
	runtime.ReadMemStats(memstats)
	log.Printf("ISC-index: memstats:\tmemstats.Alloc\tmemstats.TotalAlloc\tmemstats.Sys\tmemstats.HeapAlloc\tmemstats.HeapSys")
	isc.PrintProcessMem("ISC-index: memstats at the beginning")

	var genome_file = flag.String("g", "", "reference genome file")
	var var_prof_file = flag.String("v", "", "variant profile file")
	var idx_dir = flag.String("i", "", "index directory")
	//var workers = flag.Int("w", 1, "number of workers")
	//flag.BoolVar(&Debug, "debug", false, "Turn on debug mode.")
	flag.Parse()

	fmt.Println("Creating multigenome and variant profile index...")
	start_time := time.Now()

	sequence := isc.ReadFASTA(*genome_file)
	runtime.ReadMemStats(memstats)
	isc.PrintProcessMem("ISC-index: memstats after reading genome")

	SNP_array := isc.ReadVCF(*var_prof_file)
	runtime.ReadMemStats(memstats)
	isc.PrintProcessMem("ISC-index: memstats after reading variant profile")

	multigenome := isc.BuildMultigenome(SNP_array, sequence)

	runtime.ReadMemStats(memstats)
	isc.PrintProcessMem("ISC-index: memstats after building multigenome")

	multigenome_len := len(multigenome)
	rev_multigenome := make([]byte, multigenome_len)
	for i := range multigenome {
		rev_multigenome[i] = multigenome[multigenome_len - 1 - i]
	}

	_, genome_file_name := path.Split(*genome_file)
	multigenome_file := path.Join(*idx_dir, genome_file_name) + ".mgf"
	rev_multigenome_file := path.Join(*idx_dir, genome_file_name) + "_rev.mgf"
	_, var_prof_file_name := path.Split(*var_prof_file)
	var_prof_file := path.Join(*idx_dir, var_prof_file_name) + ".idx"

	isc.SaveMultigenome(multigenome_file, multigenome)
	//isc.SaveMultigenome(rev_multigenome_file, rev_multigenome)
	isc.SaveVarProf(var_prof_file, SNP_array)

	gen_time := time.Since(start_time)
	log.Printf("ISC-index: time for creating multigenome and variant profile index:\t%s", gen_time)

	runtime.ReadMemStats(memstats)
	isc.PrintProcessMem("ISC-index: memstats after creating multigenome and variant profile index")

	fmt.Println("Multigenome length: ", multigenome_len)
	fmt.Println("Variant profile index size: ", len(SNP_array))
	fmt.Println("Multigenome file: ", multigenome_file)
	//fmt.Println("Reverse multigenome file: ", rev_multigenome_file)
	fmt.Println("Variant profile index: ", var_prof_file)
	fmt.Println("Finish creating multigenome and variant profile index.")

	fmt.Println("Indexing multigenome...")

	var idx fmi.Index

	//start_time = time.Now()
	//idx = *fmi.New(multigenome_file)
	//idx.Save(multigenome_file)

	index_time := time.Since(start_time)
	log.Printf("ISC-index: time for indexing multigenome:\t%s", index_time)

	runtime.ReadMemStats(memstats)
	isc.PrintProcessMem("ISC-index: memstats after indexing multigenome")

	start_time = time.Now()

	idx = *fmi.New(rev_multigenome_file)
	idx.Save(rev_multigenome_file)

	index_time = time.Since(start_time)
	log.Printf("ISC-index: time for indexing reverse multigenome:\t%s", index_time)

	runtime.ReadMemStats(memstats)
	isc.PrintProcessMem("ISC-index: memstats after indexing reverse multigenome")
	
	//fmt.Println("Index directory for multigenome: ", multigenome_file + ".index/")
	fmt.Println("Index directory for reverse multigenome: ", rev_multigenome_file + ".index/")
	fmt.Println("Finish indexing multigenome.")

}
