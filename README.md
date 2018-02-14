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

### 3.1 Example command
IVC comes with a test dataset which includes the following directories:   
./test_data/refs: includes a reference genome and a corresponding variant profile for NC_007194.1 (Aspergillus fumigatus Af293 chromosome 1, whole genome shotgun sequence, see http://www.ncbi.nlm.nih.gov/nuccore/AAHF00000000)   
./test_data/reads: includes a set of 10.000 simulated paired-end reads generated with DWGSIM (see how to install DWGSIM and generate simulated reads at https://github.com/nh13/DWGSIM).

#### 3.1.1. Creating and indexing reference genomes with variant profile:
```
go run main/ivc-index.go -R test_data/refs/chr1_ref.fasta -V test_data/refs/chr1_variant_prof.vcf -I test_data/indexes
```
The command "go run main/ivc-index.go" can be replaced by the command "./ivc-index" if you have the binary file ivc-index.

#### 3.1.2. Calling variants from reads and the reference
```
go run main/ivc.go -R test_data/refs/chr1_ref.fasta -V test_data/refs/chr1_variant_prof.vcf -I test_data/indexes -1 test_data/reads/chr1_dwgsim_100_0.001-0.01.bwa.read1.fastq -2 test_data/reads/chr1_dwgsim_100_0.001-0.01.bwa.read2.fastq -O test_data/results/chr1_variant_calls.vcf
```
The command "go run main/ivc.go" can be replaced by the command "./ivc" of you have the binary file ivc.

### 3.2 Commands and options

####3.2.1. Creating and indexing reference genomes with variant profile:
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
IVC comes with a simulator which simulates mutant genomes based on the reference genome and its associated variant profile. Reads are then can be generated from the mutant genome using other simulators, such as DWGSIM.

Get the mutant genome simulator:   
```
git clone https://github.com/namsyvo/ivc-tools.git
cd ivc-tools/genome-simulator
```
Then follow the instructions to generate simulated mutant genomes and evaluate the called variants.


### 4.2 Real data
* Human (and other species) reference genomes and variant profiles (dbSNP, dbVar) can be downloaded at http://www.ncbi.nlm.nih.gov   
	* Human reference genome: http://www.ncbi.nlm.nih.gov/assembly/GCF_000001405.25   
	* Human variant profile: http://hgdownload.cse.ucsc.edu/gbdb/hg19/1000Genomes   

* Human (and other species) sequence reads can be downloaded at http://sra.dnanexus.com


## 5. Contact

Nam Sy Vo   
namsyvo@uchicago.edu   
vosynam@gmail.com