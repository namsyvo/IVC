//----------------------------------------------------------------------------------------
// IVC: index.go
// Main program for indexing module.
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
	var genome_file = flag.String("R", "", "reference genome file")
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

	//Creating multi_seq and variant profile index---------------------------//
	log.Printf("----------------------------------------------------------------------------------------")
	log.Printf("Creating multi-sequence and variant profile index...")
	log.Printf("Memstats (golang name):\tAlloc\tTotalAlloc\tSys\tHeapAlloc\tHeapSys")
	ivc.MEM_STATS = new(runtime.MemStats)

	start_time := time.Now()
	multi_seq, var_prof := ivc.BuildMultiGenome(*genome_file, *var_prof_file)
	ivc.PrintProcessMem("Memstats after building multi-sequence")

	multi_seq_len := len(multi_seq)
	rev_multi_seq := make([]byte, multi_seq_len)
	for i := range multi_seq {
		rev_multi_seq[i] = multi_seq[multi_seq_len-1-i]
	}
	_, genome_file_name := path.Split(*genome_file)
	multi_seq_file_name := path.Join(*idx_dir, genome_file_name) + ".mgf"
	rev_multi_seq_file_name := path.Join(*idx_dir, genome_file_name) + "_rev.mgf"
	_, var_prof_file_name := path.Split(*var_prof_file)
	var_prof_idx_file_name := path.Join(*idx_dir, var_prof_file_name) + ".idx"

	ivc.SaveMultiSeq(multi_seq_file_name, multi_seq)
	ivc.SaveMultiSeq(rev_multi_seq_file_name, rev_multi_seq)
	ivc.SaveVarProf(var_prof_idx_file_name, var_prof)
	gen_time := time.Since(start_time)

	log.Printf("Multi-sequence length: %d", multi_seq_len)
	log.Printf("Variant profile index size: %d", len(var_prof))
	log.Printf("Multi-sequence file: %s", multi_seq_file_name)
	log.Printf("Variant profile index file: %s", var_prof_idx_file_name)

	log.Printf("Time for creating multi-sequence and variant profile index:\t%s", gen_time)
	ivc.PrintProcessMem("Memstats after creating multi-sequence and variant profile index")
	log.Printf("Finish creating multi-sequence and variant profile index.")
	//--------------------------------------------------------------------------//

	//Indexing multi-sequence------------------------------------------------------//
	log.Printf("----------------------------------------------------------------------------------------")
	log.Printf("Indexing multi-sequence...")
	start_time = time.Now()

	var idx fmi.Index
	idx = *fmi.New(rev_multi_seq_file_name)
	idx.Save(rev_multi_seq_file_name)
	index_time := time.Since(start_time)
	log.Printf("Time for indexing reverse multi-sequence:\t%s", index_time)
	ivc.PrintProcessMem("Memstats after indexing reverse multi-sequence")
	log.Printf("Index directory for reverse multi-sequence: %s", rev_multi_seq_file_name+".index/")
	log.Printf("Finish indexing multi-sequence.")
}
