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
)

func main() {
    var genome_file = flag.String("g", "../test_data/indexes/genome_starred.txt", "multi-genome file")
    var snp_file = flag.String("s", "../test_data/indexes/SNPlocation.txt", "snp profile file")
    var index_file = flag.String("i", "../test_data/indexes/genome_starred.txt.index", "index file of multigenome")
    var rev_index_file = flag.String("r", "../test_data/indexes/genome_starred_rev.txt.index", "index file of reverse of multigenome")
    var query_file_1 = flag.String("1", "../test_data/reads/test_reads_1.fq", "pairend read file, first end")
    var query_file_2 = flag.String("2", "../test_data/reads/test_reads_2.fq", "pairend read file, second end")
    var snp_call_file = flag.String("c", "../test_data/results/test_called_snps.txt", "snp calling file")
    var read_len = flag.Int("l", 100, "read length")
    var seq_err = flag.Float64("e", 0.01, "sequencing error")
    //var workers = flag.Int("w", 1, "number of workers")
    //flag.BoolVar(&Debug, "debug", false, "Turn on debug mode.")
    flag.Parse()

	fmt.Println("ISC - Integrated SNP Calling based on Read-Multigenome Alignment")

    fmt.Println("Initializing indexes and parameters...")
    var snpcaller isc.SNPCaller
    snpcaller.Init(*genome_file, *snp_file, *index_file, *rev_index_file,
					 *read_len, float32(*seq_err), 4, 1, 32)

    fmt.Println("Aligning reads to the mutigenome...")
    var read1, read2 []byte
    var has_SNP_call bool
    var read_num int = 0
    var snp_aligned_read_num int = 0
	
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
    if q_file_1[len(q_file_1)-3:] == ".fq" || q_file_1[len(q_file_1)-6:] == ".fastq"  {
		for {
			data1.ReadBytes('\n') //ignore 1st line in input FASTQ file
			data2.ReadBytes('\n') //ignore 1st line in input FASTQ file
			line_f1, err1 = data1.ReadBytes('\n')
			line_f2, err2 = data2.ReadBytes('\n')
			if err1 != nil {
				break
            }
			if err2 != nil {
				break
            }
            if len(line_f1) >= *read_len && len(line_f2) >= *read_len{
	        	read_num++
                read1 = line_f1[0 : *read_len]
                read2 = line_f2[0 : *read_len]
				has_SNP_call = snpcaller.UpdateSNPProfile(read1, read2)
				if has_SNP_call {
	        	    snp_aligned_read_num++
		        }
            }
            data1.ReadBytes('\n') //ignore 3rd line in input FASTQ file
            data1.ReadBytes('\n') //ignore 4th line in input FASTQ file
            data2.ReadBytes('\n') //ignore 3rd line in input FASTQ file
            data2.ReadBytes('\n') //ignore 4th line in input FASTQ file
        }
    } else { // for isc format, each read is in one line
		for {
			line_f1, err1 = data1.ReadBytes('\n')
			if err1 != nil {
				break
	        }
	        if len(line_f1) >= *read_len {
				read_num++
				read1 = line_f1[0 : *read_len]
				has_SNP_call = snpcaller.UpdateSNPProfile(read1, read2)
				if has_SNP_call {
    	       	    snp_aligned_read_num++
    	        }
			}
        }
    }
    fmt.Println("\tNumber of reads: ", read_num)
    fmt.Println("\tNumber of aligned reads: ", snp_aligned_read_num)
	
    fmt.Println("Calling SNPs from alignment results...")
    snpcaller.CallSNP()
    snpcaller.SNPCall_tofile(*snp_call_file)
    fmt.Println("Finish, check the file", *snp_call_file, "for results")
}
