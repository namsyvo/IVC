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
	//"runtime/pprof"
	"log"
)

func main() {

	memstats := new(runtime.MemStats)
    log.Printf("ISC: memstats:\tmemstats.Alloc\tmemstats.TotalAlloc\tmemstats.Sys\tmemstats.HeapAlloc\tmemstats.HeapSys")

	fmt.Println("ISC - Integrated SNP Calling based on Read-Multigenome Alignment")

    var genome_file = flag.String("g", "", "multi-genome file")
    var snp_file = flag.String("s", "", "snp profile file")
    var index_file = flag.String("i", "", "index file of multigenome")
    var rev_index_file = flag.String("r", "", "index file of reverse of multigenome")
    var read_file_1 = flag.String("1", "", "pairend read file, first end")
    var read_file_2 = flag.String("2", "", "pairend read file, second end")
    var snp_call_file = flag.String("c", "", "snp calling file")
    var read_len = flag.Int("l", 100, "read length")
    var seq_err = flag.Float64("e", 0.01, "sequencing error")
	//var memprofile = flag.String("memprofile", "", "write memory profile to this file")
    //var workers = flag.Int("w", 1, "number of workers")
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

	read_info := isc.ReadInfo{}
    read_info.Read_len = *read_len
    read_info.Seq_err = float32(*seq_err)
	read_info.AllocMem()

	para_info := isc.ParaInfo{}
	para_info.Max_match = 32
	para_info.Std_dev_factor = 4
	para_info.Iter_num_factor = 1

println("@@@", input_info.SNP_file)

    runtime.ReadMemStats(memstats)
    log.Printf("isc.go: memstats at the beginning:\t%d\t%d\t%d\t%d\t%d", memstats.Alloc, memstats.TotalAlloc, memstats.Sys, memstats.HeapAlloc, memstats.HeapSys)

    fmt.Println("Initializing indexes and parameters...")

    var snpcaller isc.SNPProf
    snpcaller.Init(input_info, read_info, para_info)

	runtime.ReadMemStats(memstats)
    log.Printf("align.go: memstats after SNP caller init:\t%d\t%d\t%d\t%d\t%d", memstats.Alloc, memstats.TotalAlloc, memstats.Sys, memstats.HeapAlloc, memstats.HeapSys)

	align_mem := isc.AlignMem{}
	align_mem.AllocMem(read_info.Read_len)

    var read_num, snp_aligned_read_num int = 0, 0
	var match_pos = make([]int, isc.MAXIMUM_MATCH)

	runtime.ReadMemStats(memstats)
    log.Printf("isc.go: memstats after loading all global variables and indexes:\t%d\t%d\t%d\t%d\t%d", memstats.Alloc, memstats.TotalAlloc, memstats.Sys, memstats.HeapAlloc, memstats.HeapSys)

    fmt.Println("Aligning reads to the mutigenome...")
	
	fn1, fn2 := input_info.Read_file_1, input_info.Read_file_2
    f1, err1 := os.Open(fn1)
    if err1 != nil {
        panic("Error opening file " + fn1)
    }
    f2, err2 := os.Open(fn2)
    if err2 != nil {
        panic("Error opening file " + fn2)
    }
    data1 := bufio.NewReader(f1)
    data2 := bufio.NewReader(f2)
    var line_f1, line_f2 []byte

	//f, err := os.Create(*memprofile)
	//if err != nil {
	//	log.Fatal(err)
	//}
	//pprof.WriteHeapProfile(f)
    
    var has_SNP_call bool
	if fn1[len(fn1) - 3 : ] == ".fq" || fn1[len(fn1) - 6 : ] == ".fastq"  {
		for {
			data1.ReadBytes('\n') //ignore 1st line in input FASTQ file
			data2.ReadBytes('\n') //ignore 1st line in input FASTQ file
			line_f1, err1 = data1.ReadBytes('\n') //use 2nd line in input FASTQ file
			line_f2, err2 = data2.ReadBytes('\n') //use 2nd line in input FASTQ file
			if err1 != nil {
				break
            }
			if err2 != nil {
				break
            }
            if len(line_f1) >= isc.READ_LEN && len(line_f2) >= isc.READ_LEN{
	        	read_num++
                read_info.Read1 = line_f1[0 : isc.READ_LEN]
                read_info.Read2 = line_f2[0 : isc.READ_LEN]
				isc.RevComp(read_info.Read1, read_info.Rev_read1, read_info.Rev_comp_read1, read_info.Comp_read1)
				isc.RevComp(read_info.Read2, read_info.Rev_read2, read_info.Rev_comp_read2, read_info.Comp_read2)

				has_SNP_call = snpcaller.UpdateSNPProfile(read_info, align_mem, match_pos)
				if has_SNP_call {
	        	    snp_aligned_read_num++
		        }
            }

			//pprof.WriteHeapProfile(f)
			if (read_num % 10000 == 0) {
				runtime.ReadMemStats(memstats)
				log.Printf("isc.go: memstats after aligning each 10,000 reads:\t%d\t%d\t%d\t%d\t%d", memstats.Alloc, memstats.TotalAlloc, memstats.Sys, memstats.HeapAlloc, memstats.HeapSys)
			}

			data1.ReadBytes('\n') //ignore 3rd line in input FASTQ file
            data2.ReadBytes('\n') //ignore 3rd line in input FASTQ file

            data1.ReadBytes('\n') //ignore 4th line in input FASTQ file
            data2.ReadBytes('\n') //ignore 4th line in input FASTQ file

        }
    }
/*
 else { // for isc format, each read is in one line
		for {
			line_f1, err1 = data1.ReadBytes('\n')
			if err1 != nil {
				break
	        }
	        if len(line_f1) >= isc.READ_LEN {
				read_num++
				read1 = line_f1[0 : isc.READ_LEN]
				read2 = line_f1[0 : isc.READ_LEN]
				isc.RevComp(read1, rev_read1, rev_comp_read1, comp_read1)
				isc.RevComp(read2, rev_read2, rev_comp_read2, comp_read2)
				has_SNP_call = snpcaller.UpdateSNPProfile(read1, read2, rev_read1, rev_read2, rev_comp_read1, rev_comp_read2, comp_read1, comp_read2, bw_snp_idx, fw_snp_idx, bw_snp_val, fw_snp_val, D, T)
				if has_SNP_call {
    	       	    snp_aligned_read_num++
    	        }
			}
        }
    }
*/
    //f.Close()

	f1.Close()
	f2.Close()

	runtime.ReadMemStats(memstats)
    log.Printf("isc.go: memstats after alignment:\t%d\t%d\t%d\t%d\t%d", memstats.Alloc, memstats.TotalAlloc,
        memstats.Sys, memstats.HeapAlloc, memstats.HeapSys)

    fmt.Println("\tNumber of reads: ", read_num)
    fmt.Println("\tNumber of aligned reads: ", snp_aligned_read_num)
	
    fmt.Println("Calling SNPs from alignment results...")
    snpcaller.CallSNP()
    snpcaller.SNPCall_tofile(input_info.SNP_call_file)

	runtime.ReadMemStats(memstats)
    log.Printf("isc.go: memstats after snp call:\t%d\t%d\t%d\t%d\t%d", memstats.Alloc, memstats.TotalAlloc,
        memstats.Sys, memstats.HeapAlloc, memstats.HeapSys)

    fmt.Println("Finish, check the file", input_info.SNP_call_file, "for results")
}
