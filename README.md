# IVC - An Integrated Variant Caller


## 1. Overview

IVC is a tool for calling genomic variants from next-generation sequencing data. The tool is developed based on a novel approach to variant calling which leverages existing genetic variants to improve the accuracy of called variants, including new variants and hard-to-detect INDELs. By design, IVC integrates read alignment, alignment sorting, and variant calling into a unified process. The simplified workflow eliminates many intermediate steps and consequently reduces human intervention and errors.

IVC is written in Go programming language (see https://golang.org). It currently supports Illumina paired-end reads. Other data formats will be supported soon.


## 2. Installation

### 2.1 Download and run IVC with Go
First you need to install Go and setup its environment as guided in https://golang.org/doc/install   
Then you should create a workspace directory, e.g., $HOME/workspace/gocode, for which you will download IVC source code to, and run the following command to setup $GOPATH:   
```
export GOPATH=$HOME/workspace/gocode
```
Now you can get IVC source code by running the following command:   
```
go get github.com/namsyvo/IVC
```
After these steps, IVC source code should be in the directory $GOPATH/src/github.com/namsyvo/IVC   
Then you can go to the IVC directory, from which IVC can be run as a Go program:   
```
cd $GOPATH/src/github.com/namsyvo/IVC
```
Conventionally, one often creates 3 following sub directories in the Go workspace directory ($HOME/gocode in above example) and then put the source code into to the sub directory src/ (if you use "go get", it might create such directories for you):   
```
$GOPATH/   
    bin/   
    pkg/   
    src/
```

Then you can run IVC directly from source code using commands:   
```
go run main/ivc-index.go ...
go run main/ivc.go ...
```

A detail description of these commands will be described in section 3.1.

You can also make binary executable files for IVC by compiling IVC source code in Go:   
```
go build main/ivc-index.go 
go build main/ivc.go
```
And then run IVC using the following commands:   
```
./ivc-index
./ivc
```

### 2.2 Download and run IVC without Go
You can also get IVC without Go using Git:   
```
git clone https://github.com/namsyvo/IVC.git
cd IVC
```
Or you can download IVC without Git from its releases at https://github.com/namsyvo/IVC/releases   
Then you can run IVC using its pre-compiled binary executable files on several platforms (currently Linux, MacOS and Windows). The binary files can be found in directory binaries/ in IVC root directory.   
```
./binaries/ivc-index_macos-amd64
./binaries/ivc_macos-amd64
```

## 3. Usage

### 3.1 Example commands
IVC comes with a test dataset which includes the following directories:   
./test_data/refs: includes a reference genome and a corresponding variant profile for NC_007194.1 (Aspergillus fumigatus Af293 chromosome 1, whole genome shotgun sequence, see http://www.ncbi.nlm.nih.gov/nuccore/AAHF00000000)   
./test_data/reads: includes a set of 10.000 simulated paired-end reads generated with DWGSIM (see how to install DWGSIM and generate simulated reads at https://github.com/nh13/DWGSIM).

#### 3.1.1. Indexing reference genomes with known variant profiles
Run the following command to index the reference genome with the associated variant profile in our test data.   
```
cd $GOPATH/src/github.com/namsyvo/IVC   
go run main/ivc-index.go -R test_data/refs/chr1_ref.fasta -V test_data/refs/chr1_variant_prof.vcf -I test_data/indexes
```
The command "go run main/ivc-index.go" can be replaced by the command "./ivc-index" if you have the binary file ivc-index (with or without Go). Then you should see the following output:
```
2018/02/18 02:43:16 IVC - Integrated Variant Caller using next-generation sequencing data.   
2018/02/18 02:43:16 IVC-index: Indexing reference genomes and variant profiles.   
2018/02/18 02:43:16 ----------------------------------------------------------------------------------------   
2018/02/18 02:43:16 Creating multi-sequence and variant profile index...   
2018/02/18 02:43:16 Multi-sequence file: test_data/indexes/chr1_ref.fasta.mgf   
2018/02/18 02:43:16 Variant profile index file: test_data/indexes/chr1_variant_prof.vcf.idx   
2018/02/18 02:43:16 Time for creating multi-sequence and variant profile index: 304.72818ms   
2018/02/18 02:43:16 Finish creating multi-sequence and variant profile index.   
2018/02/18 02:43:16 ----------------------------------------------------------------------------------------   
2018/02/18 02:43:16 Indexing multi-sequence...   
2018/02/18 02:43:16 Building suffix array...   
2018/02/18 02:43:20 Finish building suffix array.   
2018/02/18 02:43:20 Building bwt and fm-index...   
2018/02/18 02:43:26 Finish building bwt and fm-index.   
2018/02/18 02:43:26 Time for indexing multi-sequence:   10.602177628s   
2018/02/18 02:43:26 Index directory for multi-sequence: test_data/indexes/chr1_ref.fasta.rev.mgf.index/   
2018/02/18 02:43:26 Finish indexing multi-sequence.   
```
The resulted index will be stored in directory "test_data/indexes".

#### 3.1.2. Calling variants from reads and the reference
Run the following command to call variants from simulated reads in our test data using the index created above.   
```
cd $GOPATH/src/github.com/namsyvo/IVC   
go run main/ivc.go -R test_data/refs/chr1_ref.fasta -V test_data/refs/chr1_variant_prof.vcf -I test_data/indexes -1 test_data/reads/chr1_dwgsim_100_0.001-0.01.bwa.read1.fastq -2 test_data/reads/chr1_dwgsim_100_0.001-0.01.bwa.read2.fastq -O test_data/results/chr1_variant_calls.vcf   
```
The command "go run main/ivc.go" can be replaced by the command "./ivc" of you have the binary file ivc (with or without Go). Then you should see the following output:   
```
2018/02/18 02:46:29 IVC - Integrated Variant Caller using next-generation sequencing data.   
2018/02/18 02:46:29 IVC-main: Calling variants based on alignment between reads and reference multi-genomes.   
2018/02/18 02:46:29 ----------------------------------------------------------------------------------------   
2018/02/18 02:46:29 Checking input information and seting up parameters...   
2018/02/18 02:46:29 No or invalid input for searching mode, use default strategy (randomizaion).   
2018/02/18 02:46:29 No or invalid input for maximum number of seeds, use default value (4096).   
    ...   
2018/02/18 02:46:29 No or invalid input for number of threads, use maximum number of CPUs of the current machine (32).   
2018/02/18 02:46:29 Input files:    Genome_file: test_data/indexes/chr1_ref.fasta.mgf, Var_file: test_data/indexes/chr1_variant_prof.vcf.idx, Index_file=test_data/indexes/chr1_ref.fasta.rev.mgf.index/, Read_file_1=test_data/reads/chr1_dwgsim_100_0.001-0.01.bwa.read1.fastq, Read_file_2=test_data/reads/chr1_dwgsim_100_0.001-0.01.bwa.read2.fastq, Var_call_file=test_data/results/chr1_variant_calls.vcf   
2018/02/18 02:46:29 Input paras:    Search_mode=1, Start_pos=0, Search_step=0, Max_snum=4096, Max_psnum=128, Min_slen=15, Max_slen=25, Dist_thres=36.0, Iter_num=12, Sub_cost=4.0, Gap_open=4.1, Gap_ext=1.0, Proc_num=32, Debug_mode=false   
2018/02/18 02:46:29 Prog paras: Max_ins=1500, Max_err=0.00150, Mut_rate=0.01000, Err_var_factor=4, Mut_var_factor=2, Iter_num_factor=2, Read_len=100, Info_len=62, Seed_backup=10, Ham_backup=15, Indel_backup=30   
2018/02/18 02:46:29 Finish checking input information and seting up parameters.   
2018/02/18 02:46:29 ----------------------------------------------------------------------------------------   
2018/02/18 02:46:29 Initializing the variant caller...   
2018/02/18 02:46:29 Loading FM-index of the reference...   
2018/02/18 02:46:29 Finish loading 10 % of index file occ.A   
2018/02/18 02:46:29 Finish loading 10 % of index file occ.G   
2018/02/18 02:46:29 Finish loading 10 % of index file sa   
    ...   
2018/02/18 02:46:29 Finish loading 100 % of index file occ.G   
2018/02/18 02:46:29 Finish loading 90 % of index file occ.T   
2018/02/18 02:46:29 Finish loading 100 % of index file occ.T   
2018/02/18 02:46:29 Finish loading FM-index of the reference.   
2018/02/18 02:46:29 Loading the reference...   
2018/02/18 02:46:29 Finish loading the reference.   
2018/02/18 02:46:29 Loading the variant profile...   
2018/02/18 02:46:29 Finish loading the variant profile.   
2018/02/18 02:46:29 Creating auxiliary data structures...   
2018/02/18 02:46:29 Finish creating auxiliary data structures.   
2018/02/18 02:46:29 Initializing variant call data structure...   
2018/02/18 02:46:29 Finish initializing 10 % of variant call data structure.   
2018/02/18 02:46:29 Finish initializing 20 % of variant call data structure.   
    ...   
2018/02/18 02:46:29 Finish initializing 100 % of variant call data structure.   
2018/02/18 02:46:29 Finish initializing variant call data structure.   
2018/02/18 02:46:29 Time for initializing the variant caller:   835.573557ms   
2018/02/18 02:46:29 Finish initializing the variant caller.   
2018/02/18 02:46:29 ----------------------------------------------------------------------------------------   
2018/02/18 02:46:29 Calling variants...   
2018/02/18 02:49:26 Processed 100000 reads.   
2018/02/18 02:49:26 Number of reads:    100000   
2018/02/18 02:49:26 Number of un-aligned reads: 4994   
2018/02/18 02:49:26 Time for calling variants:  2m56.917841744s   
2018/02/18 02:49:26 Finish calling variants.   
2018/02/18 02:49:26 ----------------------------------------------------------------------------------------   
2018/02/18 02:49:26 Outputing variant calls...   
2018/02/18 02:49:27 Time for outputing variant calls:   627.28324ms   
2018/02/18 02:49:27 Finish outputing variant calls.   
2018/02/18 02:49:27 ------------------------------------------------------   
2018/02/18 02:49:27 Check results in the file: test_data/results/chr1_variant_calls.vcf   
2018/02/18 02:49:27 Finish whole variant calling process.   
```
The resulted variant calls will be stored in file "test_data/results/chr1_variant_calls.vcf".   

### 3.2 Commands and options

#### 3.2.1. Creating and indexing reference genomes with variant profile:
Required:   
	-R: reference genome (FASTA format).  
	-V: known variant profile (VCF format).  
	-I: directory for storing index.

#### 3.2.2. Calling Variants:
Required:   
	-R: reference genome (FASTA format).  
	-V: known variant profile (VCF format).  
	-I: directory for storing index.  
	-1: the read file (for single-end reads) (FASTQ format).  
	-2: the second end file (for pair-end reads) (FASTQ format).  
	-O: variant call result file (VCF format).  

Options:   
	-d: threshold of alignment distances (float, default: determined by the program).  
	-t: maximum number of CPUs to run (integer, default: number of CPU of running computer).  
	-r: maximum number of iterations for random searching (int, default: determined by the program).  
	-s: substitution cost (float, default: 4). 
	-o: gap open cost (float, default: 4.1).   
	-e: gap extension cost (float, default: 1.0).   
	-mode: searching mode for finding seeds (1: random (default), 2: deterministic).  
	-start: starting position on reads for finding seeds (integer, default: 0).  
	-step: step for searching in deterministic mode (integer, default: 5).  
	-maxs: maximum number of seeds for single-end reads (default: 1024).  
	-maxp: maximum number of paired-seeds for paired-end reads (default: 128).  
	-lmin: minimum length of seeds for each end (default: 15).  
	-lmax: maximum length of seeds for each end (default: 30).  
	-debug: debug mode (boolean, default: false)

## 4. Data preparation

### 4.1 Simulated data
IVC comes with a simulator which simulates mutant genomes based on the reference genome and its associated variant profile. Reads are then can be generated from the mutant genome using other simulators, such as DWGSIM (https://github.com/nh13/DWGSIM).

You can run the following commands to get that simulator:   
```
cd $GOPATH/src/github.com/namsyvo   
git clone https://github.com/namsyvo/varcall-tools.git
cd varcall-tools/ivc-tools/genome-simulator
```
Then you can run the following commands to generate a simulated mutant genome from the reference and its associated variant profile in our test data and you should see the following output:   
```
go run gen_af_sid_mutant.go ../../../IVC/test_data/refs/chr1_ref.fasta ../../../IVC/test_data/refs/chr1_variant_prof.vcf simulated_data
Reading reference...
Reading variant profile...
Generating mutant genome and corresponding variant profile...
Total number of variants: 16272
Total, sub_diff_ref, ins_diff_ref, del_diff_ref, sub_same_ref, ins_same_ref, del_same_ref, SV, OL_VAR:
16272 8117 1 0 8153 0 0 0 0
Saving mutant gennome and its variant profile...
Done!
```
The resulted genome (simulated_genomes/mutant_genome.fasta) can be used by a simulator (such as DWGSIM) to generate simulated reads without mutation introduced. For example you can run the following commands to get DWGSIM and to generate simualted reads and you should see the following output:   
```
cd $GOPATH/src/github.com   
mkdir nh13; cd nh13   
git clone --recursive https://github.com/nh13/DWGSIM.git   
cd DWGSIM   
make
make[1]: Entering directory '/home/nsvo/workspace/goprojects/src/github.com/DWGSIM/samtools'
make[2]: Entering directory '/home/nsvo/workspace/goprojects/src/github.com/DWGSIM/samtools'
gcc -c -g -Wall -O3  -D_FILE_OFFSET_BITS=64 -D_LARGEFILE64_SOURCE -D_USE_KNETFILE -DPACKAGE_VERSION="0.1.11" -DBGZF_CACHE -I. bgzf.c -o bgzf.o
        ...
gcc -g -Wall -O3  -o dwgsim_eval src/dwgsim_eval.o samtools/knetfile.o samtools/bgzf.o samtools/kstring.o samtools/bam_aux.o samtools/bam.o samtools/bam_import.o samtools/sam.o samtools/bam_index.o samtools/bam_pileup.o samtools/bam_lpileup.o samtools/bam_md.o samtools/razf.o samtools/faidx.o samtools/bedidx.o samtools/bam_sort.o samtools/sam_header.o samtools/bam_reheader.o samtools/kprobaln.o samtools/bam_cat.o -Lsamtools -lm -lz -lpthread
make[1]: Leaving directory '/home/nsvo/workspace/goprojects/src/github.com/DWGSIM'

cd ../../namsyvo/varcall-tools/ivc-tools/genome-simulator   
$GOPATH/src/github.com/nh13/DWGSIM/dwgsim -e 0.001 -E 0.01 -N 100000 -1 100 -2 100 -r 0.0 -o 1 simulated_data/mutant_genome.fasta simulated_data/dwgsims   
[dwgsim_core] 1 length: 4918980   
[dwgsim_core] 1 sequences, total length: 4918980   
[dwgsim_core] Currently on:   
[dwgsim_core] 100000   
[dwgsim_core] Complete!   
```

Then you should find two fastq files which can be used as input for IVC (and BWA as well) "dwgsims.bwa.read1.fastq.gz" and "dwgsims.bwa.read2.fastq.gz" (please unzip them before using with IVC). You should also find two mutation files "dwgsims.mutations.txt" and "dwgsims.mutations.vcf" and in our simulation there should be no mutations there (the mutations are already introduced to the mutant genome "simulated_data/mutant_genome.fasta" previously by the IVC genome simulator and we are generating reads from that mutant genome without mutations). More instructions about simulating reads with DWGSIM can be found at https://github.com/nh13/DWGSIM/wiki

### 4.2 Real data
* Human (and other species) reference genomes and variant profiles (1000 Genomes, dbSNP, dbVar...) can be downloaded at http://www.ncbi.nlm.nih.gov   
	* Human reference genome (GRCh37): http://www.ncbi.nlm.nih.gov/assembly/GCF_000001405.25   
	* Human variant profile (1000 Genomes data): http://hgdownload.cse.ucsc.edu/gbdb/hg19/1000Genomes   
* Human (and other species) sequencing reads can be downloaded at http://sra.dnanexus.com


## 5. Contact

Nam Sy Vo   
namsyvo@uchicago.edu   
vosynam@gmail.com