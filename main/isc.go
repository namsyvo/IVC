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
	"path"
	"runtime"
	"time"
)

func main() {

	//Starting Program----------------------------------------------------------//
	fmt.Println("IVC - Integrated Variant Caller using Next-generation sequencing data.")
	log.Printf("memstats:\tmemstats.Alloc\tmemstats.TotalAlloc\tmemstats.Sys\tmemstats.HeapAlloc\tmemstats.HeapSys")
	//--------------------------------------------------------------------------//

	//Initializing Indexes and parameters---------------------------------------//
	fmt.Println("Initializing indexes and parameters...")
	start_time := time.Now()
	input_info := ReadInputInfo()
	runtime.GOMAXPROCS(input_info.Proc_num)
	variant_caller := isc.NewVariantCaller(input_info)
	index_time := time.Since(start_time)
	log.Printf("Time for initializing the variant caller\t%s", index_time)
	isc.PrintProcessMem("Memstats after initializing the variant caller")
	fmt.Println("Finish initializing indexes and parameters.")
	//-------------------------------------------------------------------------//

	//Call Variants from read-multigenome alignment----------------------------//
	fmt.Println("Calling variants...")
	start_time = time.Now()
	variant_caller.CallVariants()
	call_var_time := time.Since(start_time)
	log.Printf("Time for calling variants:\t%s", call_var_time)
	isc.PrintProcessMem("Memstats after calling variants")
	variant_caller.OutputVarCalls()
	fmt.Println("Check results in the file", input_info.Var_call_file)
	fmt.Println("Finish calling variants.")
	//-------------------------------------------------------------------------//
}

//--------------------------------------------------------------------------------------------------
//Read input information and parameters
//--------------------------------------------------------------------------------------------------
func ReadInputInfo() isc.InputInfo {
	var genome_file = flag.String("g", "", "reference genome file")
	var var_prof_file = flag.String("s", "", "variant profile file")
	var idx_dir = flag.String("i", "", "index directory")
	var read_file_1 = flag.String("1", "", "pairend read file, first end")
	var read_file_2 = flag.String("2", "", "pairend read file, second end")
	var var_call_file = flag.String("o", "", "variant calling file")
	var search_mode = flag.Int("m", 1, "searching mode for finding seeds (1: random, 2: deterministic)")
	var start_pos = flag.Int("p", 0, "starting position on reads for finding seeds")
	var search_step = flag.Int("j", 5, "step for searching in deterministic mode")
	var proc_num = flag.Int("w", 0, "maximum number of CPUs using by Go")
	var routine_num = flag.Int("t", 0, "number of goroutines")
	var max_snum = flag.Int("n", 1024, "maximum number of seeds")
	var min_slen = flag.Int("l", 10, "minimum length of seeds")
	var max_slen = flag.Int("h", 100, "maximum length of seeds")
	var max_psnum = flag.Int("k", 1024, "maximum number of paired-seeds")
	var dist_thres = flag.Int("d", 0, "threshold of alignment distances")
	var iter_num = flag.Int("r", 0, "maximum number of iterations")
	//flag.BoolVar(&Debug, "debug", false, "Turn on debug mode.")
	flag.Parse()

	_, genome_file_name := path.Split(*genome_file)
	multigenome_file_name := path.Join(*idx_dir, genome_file_name) + ".mgf"
	rev_multigenome_file_name := path.Join(*idx_dir, genome_file_name) + "_rev.mgf"
	_, var_prof_file_name := path.Split(*var_prof_file)
	var_prof_index_file_name := path.Join(*idx_dir, var_prof_file_name) + ".idx"

	input_info := isc.InputInfo{}
	input_info.Ref_file = multigenome_file_name
	input_info.Var_prof_file = var_prof_index_file_name
	input_info.Index_file = multigenome_file_name + ".index/"
	input_info.Rev_index_file = rev_multigenome_file_name + ".index/"
	input_info.Read_file_1 = *read_file_1
	input_info.Read_file_2 = *read_file_2
	input_info.Var_call_file = *var_call_file
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
	input_info.Max_slen = *max_slen
	input_info.Max_psnum = *max_psnum
	input_info.Dist_thres = *dist_thres
	input_info.Iter_num = *iter_num

	log.Printf("Input files:\tGenome_file: %s, Var_file: %s, Index_file: %s, Rev_index_file: %s,"+
		" Read_file_1: %s, Read_file_2: %s, Var_call_file: %s",
		input_info.Ref_file, input_info.Var_prof_file, input_info.Index_file, input_info.Rev_index_file,
		input_info.Read_file_1, input_info.Read_file_2, input_info.Var_call_file)

	log.Printf("Input parameters:\tSearch_mode: %d, Start_pos: %d, Search_step: %d, Proc_num: %d,"+
		" Routine_num: %d, Max_snum: %d, Min_slen: %d, Max_slen: %d, Max_psnum: %d, Dist_thres: %d, Iter_num: %d",
		input_info.Search_mode, input_info.Start_pos, input_info.Search_step, input_info.Proc_num, input_info.Routine_num,
		input_info.Max_snum, input_info.Min_slen, input_info.Max_slen, input_info.Max_psnum, input_info.Dist_thres, input_info.Iter_num)

	return input_info
}
