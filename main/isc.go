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

	memstats := new(runtime.MemStats)
	runtime.ReadMemStats(memstats)
	log.Printf("ISC: memstats:\tmemstats.Alloc\tmemstats.TotalAlloc\tmemstats.Sys\tmemstats.HeapAlloc\tmemstats.HeapSys")
	log.Printf("ISC: memstats at the beginning:\t%d\t%d\t%d\t%d\t%d", memstats.Alloc, memstats.TotalAlloc,
		memstats.Sys, memstats.HeapAlloc, memstats.HeapSys)

	para_info := SetPara(100, 0.01)
	input_info := ReadInputInfo()
	runtime.GOMAXPROCS(input_info.Proc_num)
	var snpcaller isc.SNPProf
	//-----------------------------------------------------------------//

	//Initializing Indexes----------------------------------------------//
	fmt.Println("Initializing indexes and parameters...")
	start_time := time.Now()
	snpcaller.InitIndex(input_info, para_info)
	index_time := time.Since(start_time)
	log.Printf("ISC: time for SNP caller init:\t%s", index_time)
	runtime.ReadMemStats(memstats)
	log.Printf("ISC: memstats after SNP caller init:\t%d\t%d\t%d\t%d\t%d", memstats.Alloc, memstats.TotalAlloc,
		memstats.Sys, memstats.HeapAlloc, memstats.HeapSys)
	//-----------------------------------------------------------------//

	//Aligning Reads---------------------------------------------------//
	fmt.Println("Aligning reads to the reference mutigenome...")
	start_time = time.Now()
	snp_aligned_read_num := snpcaller.AlignReads()
	fmt.Println("\tNumber of aligned reads: ", snp_aligned_read_num)
	align_time := time.Since(start_time)
	log.Printf("ISC: time for alignment:\t%s", align_time)
	runtime.ReadMemStats(memstats)
	log.Printf("ISC: memstats after alignment:\t%d\t%d\t%d\t%d\t%d", memstats.Alloc, memstats.TotalAlloc,
		memstats.Sys, memstats.HeapAlloc, memstats.HeapSys)
	//-----------------------------------------------------------------//

	//Calling SNPs-----------------------------------------------------//
	fmt.Println("Calling SNPs from alignment results...")
	start_time = time.Now()
	snpcaller.CallSNP()
	callsnp_time := time.Since(start_time)
	log.Printf("ISC: time for calling SNPs:\t%s", callsnp_time)
	runtime.ReadMemStats(memstats)
	log.Printf("ISC: memstats after calling SNPs:\t%d\t%d\t%d\t%d\t%d", memstats.Alloc, memstats.TotalAlloc,
		memstats.Sys, memstats.HeapAlloc, memstats.HeapSys)
	//----------------------------------------------------------------//

	//Writing SNPs-----------------------------------------------------//
	fmt.Println("Writing SNPs to output file...")
	start_time = time.Now()
	snpcaller.WriteSNPCalls()
	writetofile_time := time.Since(start_time)
	log.Printf("ISC: time for writing SNPs to file:\t%s", writetofile_time)
	runtime.ReadMemStats(memstats)
	log.Printf("ISC: memstats after writing SNPs:\t%d\t%d\t%d\t%d\t%d", memstats.Alloc, memstats.TotalAlloc,
		memstats.Sys, memstats.HeapAlloc, memstats.HeapSys)
	//----------------------------------------------------------------//

	//Finishing Program-----------------------------------------------------------//
	WriteOutputInfo(input_info)
	fmt.Println("ISC - Finish SNP calling process!")

	runtime.ReadMemStats(memstats)
	log.Printf("ISC: memstats:\tmemstats.Alloc\tmemstats.TotalAlloc\tmemstats.Sys\tmemstats.HeapAlloc\tmemstats.HeapSys")
	log.Printf("ISC: memstats at the end:\t%d\t%d\t%d\t%d\t%d", memstats.Alloc, memstats.TotalAlloc,
		memstats.Sys, memstats.HeapAlloc, memstats.HeapSys)
}

//--------------------------------------------------------------------------------------------------
//Set values for parameters
//--------------------------------------------------------------------------------------------------
func SetPara(read_len int, seq_err float32) isc.ParaInfo {
	para_info := isc.ParaInfo{}
	para_info.Max_match = 32
	para_info.Err_var_factor = 4
	para_info.Iter_num_factor = 1
	para_info.Seq_err = seq_err //will be replaced by seq_err estimated from input reads
	para_info.Read_len = read_len //will be replaced by read length taken from input reads
	return para_info
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
	var search_mode = flag.Int("m", 2, "searching mode for finding seeds (1: random, 2: deterministic)")
	var start_pos = flag.Int("p", 0, "starting position on reads for finding seeds")
	var search_step = flag.Int("j", 5, "step for searching in deterministic mode")
	var proc_num = flag.Int("w", 0, "maximum number of CPUs using by Go")
	var routine_num = flag.Int("t", 0, "number of goroutines")
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
	if *proc_num == 0 || *routine_num == 0 {
		input_info.Proc_num = runtime.NumCPU()
		input_info.Routine_num = runtime.NumCPU()
	}

	return input_info
}

//--------------------------------------------------------------------------------------------------
//Write output information and parameters
//--------------------------------------------------------------------------------------------------
func WriteOutputInfo(input_info isc.InputInfo) {
	fmt.Println("Check SNPs in the file", input_info.SNP_call_file)
}
