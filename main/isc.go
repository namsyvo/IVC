//----------------------------------------------------------------------------------------
// ISC - main program
// Copyright 2014 Nam Sy Vo
//----------------------------------------------------------------------------------------

package main

import (
	"bufio"
	"flag"
	"fmt"
	"github.com/namsyvo/ISC"
	"log"
	"os"
	"runtime"
	"sync"
	"time"
	"path"
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
	runtime.ReadMemStats(memstats)
	log.Printf("ISC: memstats:\tmemstats.Alloc\tmemstats.TotalAlloc\tmemstats.Sys\tmemstats.HeapAlloc\tmemstats.HeapSys")
	log.Printf("ISC: memstats at the beginning:\t%d\t%d\t%d\t%d\t%d", memstats.Alloc, memstats.TotalAlloc, memstats.Sys, memstats.HeapAlloc, memstats.HeapSys)

	fmt.Println("Initializing indexes and parameters...")
	start_time := time.Now()

	var genome_file = flag.String("g", "", "reference genome file")
	var dbsnp_file = flag.String("s", "", "snp profile file")
	//var index_file = flag.String("i", "", "index file of multigenome")
	//var rev_index_file = flag.String("r", "", "index file of reverse of multigenome")
	var idx_dir = flag.String("i", "", "index directory")
	var read_file_1 = flag.String("1", "", "pairend read file, first end")
	var read_file_2 = flag.String("2", "", "pairend read file, second end")
	var snp_call_file = flag.String("o", "", "snp calling file")
	var read_len = flag.Int("l", 100, "read length")
	var seq_err = flag.Float64("e", 0.01, "sequencing error")
	var search_mode = flag.Int("m", 2, "searching mode for finding seeds (1: random, 2: deterministic)")
	var start_pos = flag.Int("p", 0, "starting position on reads for finding seeds")
	var search_step = flag.Int("j", 5, "step for searching in deterministic mode")
	var proc_num = flag.Int("w", 1, "maximum number of CPUs using by Go")
	var routine_num = flag.Int("t", 1, "number of goroutines")
	//var memprofile = flag.String("memprofile", "", "write memory profile to this file")
	//flag.BoolVar(&Debug, "debug", false, "Turn on debug mode.")
	flag.Parse()

	/*Will validate all input info here ...
	if strings.ToLower(*read_file_1[len(*read_file_1) - 3: ]) != ".fq" && strings.ToLower(*read_file_1[len(*read_file_1) - 6: ]) != ".fastq"  {
		fmt.Println("Improper input FASTA file extension")
		return
	}
	*/

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

	para_info := isc.ParaInfo{}
	para_info.Max_match = 32
	para_info.Err_var_factor = 4
	para_info.Iter_num_factor = 1
	para_info.Seq_err = float32(*seq_err)
	para_info.Read_len = *read_len

	runtime.GOMAXPROCS(input_info.Proc_num)

	var snpcaller isc.SNPProf
	snpcaller.Init(input_info, para_info)

	align_info := make([]isc.AlignInfo, input_info.Routine_num)
	for i := 0; i < input_info.Routine_num; i++ {
		align_info[i].InitAlignInfo(*read_len)
	}
	match_pos := make([][]int, input_info.Routine_num)
	for i := 0; i < input_info.Routine_num; i++ {
		match_pos[i] = make([]int, isc.MAXIMUM_MATCH)
	}

	index_time := time.Since(start_time)
	log.Printf("ISC: time for SNP caller init:\t%s", index_time)

	runtime.ReadMemStats(memstats)
	log.Printf("ISC: memstats after SNP caller init:\t%d\t%d\t%d\t%d\t%d", memstats.Alloc, memstats.TotalAlloc, memstats.Sys, memstats.HeapAlloc, memstats.HeapSys)

	fmt.Println("Aligning reads to the mutigenome...")
	start_time = time.Now()

	data := make(chan isc.ReadInfo, input_info.Routine_num)
	results := make(chan []isc.SNP)

	go ReadReads(input_info, data, results)

	var wg sync.WaitGroup
	for i := 0; i < input_info.Routine_num; i++ {
		go ProcessReads(&snpcaller, data, results, &wg, align_info[i], match_pos[i])
	}
	go func() {
		wg.Wait()
		close(results)
	}()

	//Collect SNPS from results channel and update SNPs
	snp_aligned_read_num := 0
	var snp isc.SNP
	for SNPs := range results {
		snp_aligned_read_num++
		for _, snp = range SNPs {
			snpcaller.SNP_Prof[snp.SNP_Idx] = append(snpcaller.SNP_Prof[snp.SNP_Idx], snp.SNP_Val)
		}
	}
	fmt.Println("\tNumber of aligned reads: ", snp_aligned_read_num)

	align_time := time.Since(start_time)
	log.Printf("ISC: time for alignment:\t%s", align_time)

	runtime.ReadMemStats(memstats)
	log.Printf("ISC: memstats after alignment:\t%d\t%d\t%d\t%d\t%d", memstats.Alloc, memstats.TotalAlloc,
		memstats.Sys, memstats.HeapAlloc, memstats.HeapSys)

	fmt.Println("Calling SNPs from alignment results...")

	start_time = time.Now()
	snpcaller.CallSNP(input_info.Routine_num)
	callsnp_time := time.Since(start_time)
	log.Printf("ISC: time for calling SNPs:\t%s", callsnp_time)

	start_time = time.Now()
	snpcaller.SNPCall_tofile(input_info.SNP_call_file)
	writetofile_time := time.Since(start_time)
	log.Printf("ISC: time for writing SNPs to file:\t%s", writetofile_time)

	runtime.ReadMemStats(memstats)
	log.Printf("ISC: memstats after snp call:\t%d\t%d\t%d\t%d\t%d", memstats.Alloc, memstats.TotalAlloc,
		memstats.Sys, memstats.HeapAlloc, memstats.HeapSys)

	fmt.Println("Finish, check the file", input_info.SNP_call_file, "for results")
}

//--------------------------------------------------------------------------------------------------
//Read input FASTQ files and put data into data channel
//--------------------------------------------------------------------------------------------------
func ReadReads(input_info isc.InputInfo, data chan isc.ReadInfo, results chan []isc.SNP) {
	memstats := new(runtime.MemStats)

	read_info := isc.ReadInfo{}

	fn1, fn2 := input_info.Read_file_1, input_info.Read_file_2
	f1, err_f1 := os.Open(fn1)
	if err_f1 != nil {
		panic("Error opening input read file " + fn1)
	}
	defer f1.Close()
	f2, err_f2 := os.Open(fn2)
	if err_f2 != nil {
		panic("Error opening input read file " + fn2)
	}
	defer f2.Close()

	var read_num int = 0
	var line_f1, line_f2 []byte
	scanner1 := bufio.NewScanner(f1)
	scanner2 := bufio.NewScanner(f2)
	for scanner1.Scan() && scanner2.Scan() {
		//ignore 1st lines in input FASTQ files
		scanner1.Scan() //use 2nd line in input FASTQ file 1
		scanner2.Scan() //use 2nd line in input FASTQ file 2
		line_f1 = scanner1.Bytes()
		line_f2 = scanner2.Bytes()
		if len(line_f1) > 0 && len(line_f2) > 0 {
			read_num++
			read_info.AssignReads(line_f1, line_f2)
			read_info.CalcRevComp()
			data <- read_info
		}
		//pprof.WriteHeapProfile(f)
		if read_num%10000 == 0 {
			runtime.ReadMemStats(memstats)
			log.Printf("isc.go: memstats after aligning each 10,000 reads:\t%d\t%d\t%d\t%d\t%d", memstats.Alloc, memstats.TotalAlloc, memstats.Sys, memstats.HeapAlloc, memstats.HeapSys)
		}
		scanner1.Scan() //ignore 3rd line in 1st input FASTQ file 1
		scanner2.Scan() //ignore 3rd line in 2nd input FASTQ file 2
		scanner1.Scan() //ignore 4th line in 1st input FASTQ file 1
		scanner2.Scan() //ignore 4th line in 2nd input FASTQ file 2
	}
	fmt.Println("\tNumber of inout reads: ", read_num)
	close(data)
}

//--------------------------------------------------------------------------------------------------
//Take data from data channel, process them (find SNPs) and put results (SNPs) into results channel
//--------------------------------------------------------------------------------------------------
func ProcessReads(snpcaller *isc.SNPProf, data chan isc.ReadInfo, results chan []isc.SNP,
	wg *sync.WaitGroup, align_info isc.AlignInfo, match_pos []int) {
	wg.Add(1)
	defer wg.Done()
	var read_info isc.ReadInfo
	var SNPs []isc.SNP
	for read_info = range data {
		SNPs = (*snpcaller).FindSNP(read_info, align_info, match_pos)
		if len(SNPs) > 0 {
			results <- SNPs
		}
	}
}
