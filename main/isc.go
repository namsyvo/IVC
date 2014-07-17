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
    var query_file_1 = flag.String("1", "", "pairend read file, first end")
    var query_file_2 = flag.String("2", "", "pairend read file, second end")
    var snp_call_file = flag.String("c", "", "snp calling file")
    var read_len = flag.Int("l", 100, "read length")
    var seq_err = flag.Float64("e", 0.01, "sequencing error")
	//var memprofile = flag.String("memprofile", "", "write memory profile to this file")
    //var workers = flag.Int("w", 1, "number of workers")
    //flag.BoolVar(&Debug, "debug", false, "Turn on debug mode.")
    flag.Parse()

    fmt.Println("Initializing indexes and parameters...")

    runtime.ReadMemStats(memstats)
    log.Printf("isc.go: memstats before indexing:\t%d\t%d\t%d\t%d\t%d", memstats.Alloc, memstats.TotalAlloc,
        memstats.Sys, memstats.HeapAlloc, memstats.HeapSys)

    var snpcaller isc.SNPProf
    snpcaller.Init(*genome_file, *snp_file, *index_file, *rev_index_file, *read_len, float32(*seq_err), 4, 1, 32)

	var bw_snp_idx, fw_snp_idx []int
	var bw_snp_val, fw_snp_val [][]byte
	var bw_D, fw_D [][]int
	var bw_T, fw_T [][][]byte

	InitMatrix(*read_len, &bw_snp_idx, &bw_snp_val, &bw_D, &bw_T)
	InitMatrix(*read_len, &fw_snp_idx, &fw_snp_val, &fw_D, &fw_T)

    var read1, read2, rev_read1, rev_read2, rev_comp_read1, rev_comp_read2, comp_read1, comp_read2 []byte
	read1 = make([]byte, *read_len)
	read2 = make([]byte, *read_len)
	rev_read1 = make([]byte, *read_len)
	rev_read2 = make([]byte, *read_len)
	rev_comp_read1 = make([]byte, *read_len)
	rev_comp_read2 = make([]byte, *read_len)
	comp_read1 = make([]byte, *read_len)
	comp_read2 = make([]byte, *read_len)

    var read_num, snp_aligned_read_num int = 0, 0
	var match_pos = make([]int, isc.MAXIMUM_MATCH)

	runtime.ReadMemStats(memstats)
    log.Printf("isc.go: memstats after loading all global variables and indexes:\t%d\t%d\t%d\t%d\t%d", memstats.Alloc, memstats.TotalAlloc,
        memstats.Sys, memstats.HeapAlloc, memstats.HeapSys)

    fmt.Println("Aligning reads to the mutigenome...")
	
    var q_file_1 string = *query_file_1
    var q_file_2 string = *query_file_2
    f1, err1 := os.Open(q_file_1)
    if err1 != nil {
        panic("Error opening file " + q_file_1)
    }
    f2, err2 := os.Open(q_file_2)
    if err2 != nil {
        panic("Error opening file " + q_file_2)
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
	if q_file_1[len(q_file_1) - 3 : ] == ".fq" || q_file_1[len(q_file_1) - 6 : ] == ".fastq"  {
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
                read1 = line_f1[0 : isc.READ_LEN]
                read2 = line_f2[0 : isc.READ_LEN]
				isc.RevComp(read1, rev_read1, rev_comp_read1, comp_read1)
				isc.RevComp(read2, rev_read2, rev_comp_read2, comp_read2)

				has_SNP_call = snpcaller.UpdateSNPProfile(read1, read2, rev_read1, rev_read2, rev_comp_read1, rev_comp_read2, comp_read1, comp_read2, bw_snp_idx, fw_snp_idx, bw_snp_val, fw_snp_val, bw_D, fw_D, bw_T, fw_T, match_pos)
				if has_SNP_call {
	        	    snp_aligned_read_num++
		        }
            }

			//pprof.WriteHeapProfile(f)
			if (read_num % 10000 == 0) {
				runtime.ReadMemStats(memstats)
				log.Printf("isc.go: memstats after aligning each 10,000 reads:\t%d\t%d\t%d\t%d\t%d", memstats.Alloc, memstats.TotalAlloc,
					memstats.Sys, memstats.HeapAlloc, memstats.HeapSys)
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
    snpcaller.SNPCall_tofile(*snp_call_file)

	runtime.ReadMemStats(memstats)
    log.Printf("isc.go: memstats after snp call:\t%d\t%d\t%d\t%d\t%d", memstats.Alloc, memstats.TotalAlloc,
        memstats.Sys, memstats.HeapAlloc, memstats.HeapSys)

    fmt.Println("Finish, check the file", *snp_call_file, "for results")
}

//-------------------------------------------------------------------------------------------------
// Initializing variables for computing distance and alignment between reads and multi-genomes.
//-------------------------------------------------------------------------------------------------
func InitMatrix(arr_len int, snp_idx *[]int, snp_val *[][]byte, D *[][]int, T *[][][]byte) {
	*snp_idx = make([]int, arr_len)
	*snp_val = make([][]byte, arr_len)
	*D = make([][]int, arr_len + 1)
	for i:= 0; i <= arr_len; i++ {
		(*D)[i] = make([]int, arr_len + 1)
    }
	*T = make([][][]byte, arr_len)
    for i:= 0; i < arr_len; i++ {
		(*T)[i] = make([][]byte, arr_len)
    }
}
