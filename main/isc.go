//----------------------------------------------------------------------------------------
// ISC - main program
// Copyright 2014 Nam Sy Vo
//----------------------------------------------------------------------------------------

package main

import (
    "fmt"
    "os"
    "bufio"
    "flag"
    "github.com/namsyvo/ISC"
	"runtime"
	"time"
	"log"
	"sync"
	//"runtime/pprof"
)

func main() {

	fmt.Println("ISC - Integrated SNP Calling based on Read-Multigenome Alignment")

		
	//f, err := os.Create(*memprofile)
	//if err != nil {
	//	log.Fatal(err)
	//}
    //defer f.Close()
	//pprof.WriteHeapProfile(f)

	memstats := new(runtime.MemStats)
    log.Printf("ISC: memstats:\tmemstats.Alloc\tmemstats.TotalAlloc\tmemstats.Sys\tmemstats.HeapAlloc\tmemstats.HeapSys")

    runtime.ReadMemStats(memstats)
    log.Printf("isc.go: memstats at the beginning:\t%d\t%d\t%d\t%d\t%d", memstats.Alloc, memstats.TotalAlloc, memstats.Sys, memstats.HeapAlloc, memstats.HeapSys)

    fmt.Println("Initializing indexes and parameters...")
	start_time := time.Now()

    var genome_file = flag.String("g", "", "multi-genome file")
    var snp_file = flag.String("s", "", "snp profile file")
    var index_file = flag.String("i", "", "index file of multigenome")
    var rev_index_file = flag.String("r", "", "index file of reverse of multigenome")
    var read_file_1 = flag.String("1", "", "pairend read file, first end")
    var read_file_2 = flag.String("2", "", "pairend read file, second end")
    var snp_call_file = flag.String("c", "", "snp calling file")
    var read_len = flag.Int("l", 100, "read length")
    var seq_err = flag.Float64("e", 0.01, "sequencing error")
    var search_mode = flag.Int("m", 2, "searching mode for finding seeds(1: random, 2: deterministic)")
    var start_pos = flag.Int("p", 0, "starting position on reads for finding seeds")
    var search_step = flag.Int("j", 5, "step for searching in deterministic mode")
	var proc_num = flag.Int("w", 1, "maximum number of CPUs using by Go")
	var routine_num = flag.Int("t", 1, "number of goroutines")
	//var memprofile = flag.String("memprofile", "", "write memory profile to this file")
    //flag.BoolVar(&Debug, "debug", false, "Turn on debug mode.")
    flag.Parse()

	input_info := isc.InputInfo{}
    input_info.Genome_file = *genome_file
    input_info.SNP_file = *snp_file
	input_info.Index_file = *index_file
    input_info.Rev_index_file = *rev_index_file
    input_info.Read_file_1 = *read_file_1
    input_info.Read_file_2 = *read_file_2
    input_info.SNP_call_file = *snp_call_file
    input_info.Search_mode = *search_mode
    input_info.Start_pos = *start_pos
    input_info.Search_step = *search_step
	input_info.Proc_num = *proc_num
	input_info.Routine_num = *routine_num

	para_info := isc.ParaInfo{}
	para_info.Max_match = 32
	para_info.Err_var_factor = 4
	para_info.Iter_num_factor = 1

	read_info := isc.ReadInfo{}
    read_info.Read_len = *read_len
    read_info.Seq_err = float32(*seq_err)

	runtime.GOMAXPROCS(input_info.Proc_num)

    var snpcaller isc.SNPProf
    snpcaller.Init(input_info, read_info, para_info)

	index_time := time.Since(start_time)
	log.Printf("ISC: time for SNP caller init:\t%s", index_time)

	runtime.ReadMemStats(memstats)
    log.Printf("align.go: memstats after SNP caller init:\t%d\t%d\t%d\t%d\t%d", memstats.Alloc, memstats.TotalAlloc, memstats.Sys, memstats.HeapAlloc, memstats.HeapSys)

	align_mem := make([]isc.AlignMem, input_info.Routine_num)
	for i := 0; i < input_info.Routine_num; i++ {
		align_mem[i].AllocMem(read_info.Read_len)
	}
	match_pos := make([][]int, input_info.Routine_num)
	for i := 0; i < input_info.Routine_num; i++ {
		match_pos[i] = make([]int, isc.MAXIMUM_MATCH)
	}

	runtime.ReadMemStats(memstats)
    log.Printf("isc.go: memstats after loading all global variables and indexes:\t%d\t%d\t%d\t%d\t%d", memstats.Alloc, memstats.TotalAlloc, memstats.Sys, memstats.HeapAlloc, memstats.HeapSys)
		
    fmt.Println("Aligning reads to the mutigenome...")
	start_time = time.Now()

	data := make([]chan isc.ReadInfo, input_info.Routine_num)
	for i := 0; i < input_info.Routine_num; i++ {
		data[i] = make(chan isc.ReadInfo)
	}
	results := make(chan []isc.SNP)
	quit := make(chan int)

	go ReadReads(input_info, data, results, quit)

	var wg sync.WaitGroup
	wg.Add(input_info.Routine_num)
	for i := 0; i < input_info.Routine_num; i++ {
		go ProcessReads(&snpcaller, data[i], results, quit, &wg, align_mem[i], match_pos[i])
	}
	go func() {
		wg.Wait()
		close(results)
	}()

	//Collect SNPS from results channel and update SNPs
    snp_aligned_read_num := 0
	var snp isc.SNP
	for SNPs := range(results) {
		snp_aligned_read_num++
		for _, snp = range SNPs {
			snpcaller.SNP_Prof[snp.SNP_Idx] = append(snpcaller.SNP_Prof[snp.SNP_Idx], snp.SNP_Val)
		}
	}
    fmt.Println("\tNumber of aligned reads: ", snp_aligned_read_num)

	align_time := time.Since(start_time)
	log.Printf("ISC: time for alignment:\t%s", align_time)

	runtime.ReadMemStats(memstats)
    log.Printf("isc.go: memstats after alignment:\t%d\t%d\t%d\t%d\t%d", memstats.Alloc, memstats.TotalAlloc,
        memstats.Sys, memstats.HeapAlloc, memstats.HeapSys)
	
    fmt.Println("Calling SNPs from alignment results...")
	start_time = time.Now()

    snpcaller.CallSNP()
	runtime.ReadMemStats(memstats)
    log.Printf("isc.go: memstats after snp call:\t%d\t%d\t%d\t%d\t%d", memstats.Alloc, memstats.TotalAlloc,
        memstats.Sys, memstats.HeapAlloc, memstats.HeapSys)

    snpcaller.SNPCall_tofile(input_info.SNP_call_file)
	runtime.ReadMemStats(memstats)
    log.Printf("isc.go: memstats at the end:\t%d\t%d\t%d\t%d\t%d", memstats.Alloc, memstats.TotalAlloc,
        memstats.Sys, memstats.HeapAlloc, memstats.HeapSys)

	callsnp_time := time.Since(start_time)
	log.Printf("ISC: time for calling SNPs:\t%s", callsnp_time)

    fmt.Println("Finish, check the file", input_info.SNP_call_file, "for results")
}

func ReadReads(input_info isc.InputInfo, data []chan isc.ReadInfo, results chan []isc.SNP, quit chan int) {
	memstats := new(runtime.MemStats)
	
	var i, j int

	read_info := make([]isc.ReadInfo, input_info.Routine_num)
	for i = 0; i < input_info.Routine_num; i++ {
		read_info[i].Read_len = 100
		read_info[i].Seq_err = float32(0.01)
	}
	
	fn1, fn2 := input_info.Read_file_1, input_info.Read_file_2
	f1, err_f1 := os.Open(fn1)
	if err_f1 != nil {
		panic("Error opening file " + fn1)
	}
	defer f1.Close()
	f2, err_f2 := os.Open(fn2)
	if err_f2 != nil {
    	panic("Error opening file " + fn2)
	}
	defer f2.Close()
	
	var read_num int
	var line_f1, line_f2 []byte
	data1 := bufio.NewReader(f1)
	data2 := bufio.NewReader(f2)
	if fn1[len(fn1)-3: ] == ".fq" || fn1[len(fn1)-6: ] == ".fastq"  {
		for {
			for i = 0; i < input_info.Routine_num; i++ {
				data1.ReadBytes('\n') //ignore 1st line in input FASTQ file
				data2.ReadBytes('\n') //ignore 1st line in input FASTQ file
				line_f1, err_f1 = data1.ReadBytes('\n') //use 2nd line in input FASTQ file
				line_f2, err_f2 = data2.ReadBytes('\n') //use 2nd line in input FASTQ file
				if err_f1 != nil || err_f2 != nil {
					fmt.Println("\tNumber of reads: ", read_num)
					for j = 0; j < input_info.Routine_num; j++ {
						quit <- 0
					}
					return
				}
				if len(line_f1) >= read_info[i].Read_len && len(line_f2) >= read_info[i].Read_len{
	        		read_num++
					read_info[i].Read1 = line_f1[0 : read_info[i].Read_len]
					read_info[i].Read2 = line_f2[0 : read_info[i].Read_len]
					read_info[i].AllocMem()
					isc.RevComp(read_info[i].Read1, read_info[i].Rev_read1, read_info[i].Rev_comp_read1, read_info[i].Comp_read1)
					isc.RevComp(read_info[i].Read2, read_info[i].Rev_read2, read_info[i].Rev_comp_read2, read_info[i].Comp_read2)
					data[i] <- read_info[i]
				}
				//pprof.WriteHeapProfile(f)
				if (read_num % 10000 == 0) {
					runtime.ReadMemStats(memstats)
					log.Printf("isc.go: memstats after aligning each 10,000 reads:\t%d\t%d\t%d\t%d\t%d", memstats.Alloc, memstats.TotalAlloc, memstats.Sys, memstats.HeapAlloc, memstats.HeapSys)
				}
				data1.ReadBytes('\n') //ignore 3rd line in 1st input FASTQ file
				data2.ReadBytes('\n') //ignore 3rd line in 2nd input FASTQ file
				data1.ReadBytes('\n') //ignore 4th line in 1st input FASTQ file
				data2.ReadBytes('\n') //ignore 4th line in 2nd input FASTQ file
			}
		}
	}
}

func ProcessReads(snpcaller *isc.SNPProf, data chan isc.ReadInfo, results chan []isc.SNP, quit chan int, wg *sync.WaitGroup, align_mem isc.AlignMem, match_pos []int) {
	defer wg.Done()
	var read_info isc.ReadInfo
	var SNPs []isc.SNP
	for {
		select {
		case read_info = <- data:
			SNPs = (*snpcaller).FindSNP(read_info, align_mem, match_pos)
			if len(SNPs) > 0 {
				results <- SNPs
			}
		case <- quit:
			return
		}
	}
}