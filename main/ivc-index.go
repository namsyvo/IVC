//----------------------------------------------------------------------------------------
// IVC: index.go - Indexing program.
// Copyright 2015 Nam Sy Vo.
//----------------------------------------------------------------------------------------

package main

import (
	"flag"
	"github.com/namsyvo/IVC"
	"github.com/vtphan/fmi"
	"log"
	"os"
	"path"
	"runtime"
	"time"
)

func main() {
	//Starting program------------------------------------------------------------------//
	log.Printf("IVC - Integrated Variant Caller using next-generation sequencing data.")
	log.Printf("IVC-index: Indexing reference genomes and variant profiles.")
	//----------------------------------------------------------------------------------//
	var genome_file = flag.String("G", "", "reference genome file")
	var var_prof_file = flag.String("V", "", "variant profile file")
	var idx_dir = flag.String("I", "", "index directory")
	flag.Parse()

	if _, err := os.Stat(*idx_dir); err != nil {
		if os.IsNotExist(err) {
			if err := os.Mkdir(*idx_dir, 0777); err != nil {
				log.Printf("Could not create index directory, possibly due to a path error.")
				os.Exit(1)
			}
		} else {
			os.Exit(1)
		}
	}

	//Creating multigenome and variant profile index---------------------------//
	log.Printf("Creating multigenome and variant profile index...")
	ivc.MEM_STATS = new(runtime.MemStats)
	log.Printf("IVC-index: memstats:\tmemstats.Alloc\tmemstats.TotalAlloc\tmemstats.Sys\tmemstats.HeapAlloc\tmemstats.HeapSys")

	start_time := time.Now()
	genome := ivc.ReadFASTA(*genome_file)
	ivc.PrintProcessMem("Memstats after reading reference genome")
	var_prof := ivc.ReadVCF(*var_prof_file)
	ivc.PrintProcessMem("Memstats after reading variant profile")
	multigenome := ivc.BuildMultigenome(var_prof, genome)
	ivc.PrintProcessMem("Memstats after building multigenome")

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
	log.Printf("Time for creating multigenome and variant profile index:\t%s", gen_time)
	ivc.PrintProcessMem("Memstats after creating multigenome and variant profile index")

	log.Printf("Multigenome length: %d", multigenome_len)
	log.Printf("Variant profile index size: %d", len(var_prof))
	log.Printf("Multigenome file: %s", multigenome_file_name)
	//log.Printf("Reverse multigenome file: %s", rev_multigenome_file_name)
	log.Printf("Variant profile index: %s", var_prof_idx_file_name)
	log.Printf("Finish creating multigenome and variant profile index.")
	//--------------------------------------------------------------------------//

	//Indexing multigenome------------------------------------------------------//
	log.Printf("Indexing multigenome...")
	var idx fmi.Index

	//start_time = time.Now()
	//idx = *fmi.New(multigenome_file)
	//idx.Save(multigenome_file)
	//index_time := time.Since(start_time)
	//log.Printf("Time for indexing multigenome:\t%s", index_time)
	//ivc.PrintProcessMem("Memstats after indexing multigenome")

	start_time = time.Now()
	idx = *fmi.New(rev_multigenome_file_name)
	idx.Save(rev_multigenome_file_name)
	index_time := time.Since(start_time)
	log.Printf("Time for indexing reverse multigenome:\t%s", index_time)
	ivc.PrintProcessMem("Memstats after indexing reverse multigenome")

	//log.Printf("Index directory for multigenome: %s", multigenome_file+".index/")
	log.Printf("Index directory for reverse multigenome: %s", rev_multigenome_file_name+".index/")
	log.Printf("Finish indexing multigenome.")

}
