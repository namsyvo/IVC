//----------------------------------------------------------------------------------------
// ISC - main program
// Copyright 2014 Nam Sy Vo
//----------------------------------------------------------------------------------------

package main

import (
	"flag"
	"fmt"
	"github.com/namsyvo/ISC"
	"log"
	"runtime"
	"time"
	"path"
)

func main() {

	//Starting Program-----------------------------------------------------------//
	fmt.Println("ISC - Integrated SNP Calling based on Read-Multigenome Alignment")

	log.Printf("memstats:\tmemstats.Alloc\tmemstats.TotalAlloc\tmemstats.Sys\tmemstats.HeapAlloc\tmemstats.HeapSys")
	isc.PrintProcessMem("memstats at the beginning")

	input_info := ReadInputInfo()
	runtime.GOMAXPROCS(input_info.Proc_num)
	var snp_prof isc.SNP_Prof
	//--------------------------------------------------------------------------//

	//Initializing Indexes------------------------------------------------------//
	fmt.Println("Initializing indexes and parameters...")
	start_time := time.Now()
	snp_prof.Init(input_info)
	index_time := time.Since(start_time)
	log.Printf("time for initializing SNP caller\t%s", index_time)
	isc.PrintProcessMem("memstats after initializing SNP caller")
	//-------------------------------------------------------------------------//

	//Call SNPs from read-multigenome alignment--------------------------------//
	fmt.Println("Calling SNPs based on aligning reads to the mutigenome...")
	start_time = time.Now()
	snp_prof.CallSNPs()
	align_time := time.Since(start_time)
	log.Printf("time for calling SNPs:\t%s", align_time)
	isc.PrintProcessMem("memstats after calling SNPs")
	//-------------------------------------------------------------------------//

	//Finishing Program--------------------------------------------------------//
	WriteOutputInfo(input_info)
	fmt.Println("Finish SNP calling process.")
	//-------------------------------------------------------------------------//
}

//--------------------------------------------------------------------------------------------------
//Read input information and parameters
//--------------------------------------------------------------------------------------------------
func ReadInputInfo() isc.InputInfo {
	var genome_file = flag.String("g", "", "reference genome file")
	var dbsnp_file = flag.String("s", "", "snp profile file")
	var idx_dir = flag.String("i", "", "index directory")
	var read_file_1 = flag.String("1", "", "pairend read file, first end")
	var read_file_2 = flag.String("2", "", "pairend read file, second end")
	var snp_call_file = flag.String("o", "", "snp calling file")
	var search_mode = flag.Int("m", 1, "searching mode for finding seeds (1: random, 2: deterministic)")
	var start_pos = flag.Int("p", 0, "starting position on reads for finding seeds")
	var search_step = flag.Int("j", 5, "step for searching in deterministic mode")
	var proc_num = flag.Int("w", 0, "maximum number of CPUs using by Go")
	var routine_num = flag.Int("t", 0, "number of goroutines")
	var max_snum = flag.Int("n", 1024, "maximum number of seeds")
	var min_slen = flag.Int("l", 10, "minimum length of seeds")
	//flag.BoolVar(&Debug, "debug", false, "Turn on debug mode.")
	flag.Parse()

	_, genome_file_name := path.Split(*genome_file)
	multigenome_file := path.Join(*idx_dir, genome_file_name) + ".mgf"
	rev_multigenome_file := path.Join(*idx_dir, genome_file_name) + "_rev.mgf"
	_, dbsnp_file_name := path.Split(*dbsnp_file)
	snp_prof_file := path.Join(*idx_dir, dbsnp_file_name) + ".idx"

	input_info := isc.InputInfo{}
	input_info.Genome_file = multigenome_file
	input_info.SNP_file = snp_prof_file
	input_info.Index_file = multigenome_file + ".index/"
	input_info.Rev_index_file = rev_multigenome_file + ".index/"
	input_info.Read_file_1 = *read_file_1
	input_info.Read_file_2 = *read_file_2
	input_info.SNP_call_file = *snp_call_file
	input_info.Search_mode = *search_mode
	input_info.Start_pos = *start_pos
	input_info.Search_step = *search_step
	input_info.Proc_num = *proc_num
	input_info.Routine_num = *routine_num
	if *proc_num <= 0 || *routine_num <= 0 {
		input_info.Proc_num = runtime.NumCPU()
		input_info.Routine_num = runtime.NumCPU()
	}
	input_info.Max_snum = *max_snum
	input_info.Min_slen = *min_slen

	log.Printf("Input files: Genome_file: %s, SNP_file: %s, Index_file: %s, Rev_index_file: %s," + 
		" Read_file_1: %s, Read_file_2: %s, SNP_call_file: %s", 
		input_info.Genome_file, input_info.SNP_file, input_info.Index_file, input_info.Rev_index_file, input_info.Read_file_1, 
		input_info.Read_file_2, input_info.SNP_call_file)

	log.Printf("Input parameters: Search_mode: %d, Start_pos: %d, Search_step: %d, Proc_num: %d," + 
		" Routine_num: %d, Max_snum: %d, Min_slen: %d", 
		input_info.Search_mode, input_info.Start_pos, input_info.Search_step, input_info.Proc_num, input_info.Routine_num, 
		input_info.Max_snum, input_info.Min_slen)

	return input_info
}

//--------------------------------------------------------------------------------------------------
//Write output information and parameters
//--------------------------------------------------------------------------------------------------
func WriteOutputInfo(input_info isc.InputInfo) {
	fmt.Println("Check SNPs in the file", input_info.SNP_call_file)
}
