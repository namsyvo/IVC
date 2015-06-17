//----------------------------------------------------------------------------------------
// IVC: index.go - Indexing program.
// Copyright 2015 Nam Sy Vo.
//----------------------------------------------------------------------------------------

package main

import (
	"flag"
	"fmt"
	"github.com/namsyvo/IVC"
	"github.com/vtphan/fmi"
	"log"
	"path"
	"runtime"
	"time"
)

func main() {
	fmt.Println("IVC - Integrated Variant Caller using Next-generation sequencing data.")
	fmt.Println("IVC-index: Indexing reference genomes and variant profiles.")

	memstats := new(runtime.MemStats)
	runtime.ReadMemStats(memstats)
	log.Printf("IVC-index: memstats:\tmemstats.Alloc\tmemstats.TotalAlloc\tmemstats.Sys\tmemstats.HeapAlloc\tmemstats.HeapSys")
	ivc.PrintProcessMem("IVC-index: memstats at the beginning")

	var genome_file = flag.String("g", "", "reference genome file")
	var var_prof_file = flag.String("v", "", "variant profile file")
	var idx_dir = flag.String("i", "", "index directory")
	//var workers = flag.Int("w", 1, "number of workers")
	//flag.BoolVar(&Debug, "debug", false, "Turn on debug mode.")
	flag.Parse()

	fmt.Println("Creating multigenome and variant profile index...")
	start_time := time.Now()

	genome := ivc.ReadFASTA(*genome_file)
	runtime.ReadMemStats(memstats)
	ivc.PrintProcessMem("IVC-index: memstats after reading genome")

	var_prof := ivc.ReadVCF(*var_prof_file)
	runtime.ReadMemStats(memstats)
	ivc.PrintProcessMem("IVC-index: memstats after reading variant profile")

	multigenome := ivc.BuildMultigenome(var_prof, genome)

	runtime.ReadMemStats(memstats)
	ivc.PrintProcessMem("IVC-index: memstats after building multigenome")

	multigenome_len := len(multigenome)
	rev_multigenome := make([]byte, multigenome_len)
	for i := range multigenome {
		rev_multigenome[i] = multigenome[multigenome_len-1-i]
	}

	_, genome_file_name := path.Split(*genome_file)
	multigenome_file_name := path.Join(*idx_dir, genome_file_name) + ".mgf"
	rev_multigenome_file_name := path.Join(*idx_dir, genome_file_name) + "_rev.mgf"
	_, var_prof_file_name := path.Split(*var_prof_file)
	var_prof_idx_file_name := path.Join(*idx_dir, var_prof_file_name) + ".idx"

	ivc.SaveMultigenome(multigenome_file_name, multigenome)
	ivc.SaveMultigenome(rev_multigenome_file_name, rev_multigenome)
	ivc.SaveVarProf(var_prof_idx_file_name, var_prof)

	gen_time := time.Since(start_time)
	log.Printf("IVC-index: time for creating multigenome and variant profile index:\t%s", gen_time)

	runtime.ReadMemStats(memstats)
	ivc.PrintProcessMem("IVC-index: memstats after creating multigenome and variant profile index")

	fmt.Println("Multigenome length: ", multigenome_len)
	fmt.Println("Variant profile index size: ", len(var_prof))
	fmt.Println("Multigenome file: ", multigenome_file_name)
	//fmt.Println("Reverse multigenome file: ", rev_multigenome_file_name)
	fmt.Println("Variant profile index: ", var_prof_idx_file_name)
	fmt.Println("Finish creating multigenome and variant profile index.")

	fmt.Println("Indexing multigenome...")

	var idx fmi.Index

	//start_time = time.Now()
	//idx = *fmi.New(multigenome_file)
	//idx.Save(multigenome_file)

	index_time := time.Since(start_time)
	log.Printf("IVC-index: time for indexing multigenome:\t%s", index_time)

	runtime.ReadMemStats(memstats)
	ivc.PrintProcessMem("IVC-index: memstats after indexing multigenome")

	start_time = time.Now()

	idx = *fmi.New(rev_multigenome_file_name)
	idx.Save(rev_multigenome_file_name)

	index_time = time.Since(start_time)
	log.Printf("IVC-index: time for indexing reverse multigenome:\t%s", index_time)

	runtime.ReadMemStats(memstats)
	ivc.PrintProcessMem("IVC-index: memstats after indexing reverse multigenome")

	//fmt.Println("Index directory for multigenome: ", multigenome_file + ".index/")
	fmt.Println("Index directory for reverse multigenome: ", rev_multigenome_file_name+".index/")
	fmt.Println("Finish indexing multigenome.")

}
