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
    var genome_file = flag.String("g", "", "multi-genome file")
    var snp_file = flag.String("s", "", "snp profile file")
    var index_file = flag.String("i", "", "index file of multigenome")
    var rev_index_file = flag.String("r", "", "index file of reverse of multigenome")
    var query_file = flag.String("q", "", "query file")
    var snp_call_file = flag.String("c", "", "snp calling file")
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
    var read []byte
    var has_SNP_call bool
    var read_num int = 0
    var snp_aligned_read_num int = 0
	
    var q_file string = *query_file
    f, err := os.Open(q_file)
    if err != nil {
        panic("Error opening file " + q_file)
    }
    data := bufio.NewReader(f)
    var line []byte
    if q_file[len(q_file)-3:] == ".fq" || q_file[len(q_file)-6:] == ".fastq"  {
		for {
			data.ReadBytes('\n') //ignore 1st line in input FASTQ file
			line, err = data.ReadBytes('\n')
			if err != nil {
				break
            }
            if len(line) >= *read_len {
	        	read_num++
                read = line[0 : *read_len]
				has_SNP_call = snpcaller.UpdateSNPProfile(read)
				if has_SNP_call {
	        	    snp_aligned_read_num++
		        }
            }
            data.ReadBytes('\n') //ignore 3rd line in input FASTQ file
            data.ReadBytes('\n') //ignore 4th line in input FASTQ file
        }
    } else { // for simple tab-delimited format
		for {
			line, err = data.ReadBytes('\n')
			if err != nil {
				break
	        }
	        if len(line) >= *read_len {
				read_num++
				read = line[0 : *read_len]
				has_SNP_call = snpcaller.UpdateSNPProfile(read)
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
