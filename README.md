IVC - An Integrated Variant Caller
==================================


1. Overview
-----------

IVC is a tool for calling genomic variants from next-generation sequencing data. The tool is developed based on a novel approach to variant calling which leverages existing genetic variants to improve the accuracy of called variants, including new variants and hard-to-detect INDELs. By design, IVC integrates read alignment, alignment sorting, and variant calling into a unified process. The simplified workflow eliminates many intermediate steps and consequently reduces human intervention and errors.

IVC is written in Go programming language (see https://golang.org). It currently supports Illumina paired-end reads. Other data formats will be supported soon.


2. Install IVC
--------------

### 2.1 Download IVC source code with Go
Pre-requirement: Go environment is already set up properly.

Get IVC source code:
```
go get github.com/namsyvo/IVC
```
After these steps, IVC source code should be in the directory $GOPATH/github.com/namsyvo/IVC   
Then go to the IVC directory, from which IVC can be run as a Go program:
```
cd $GOPATH/src/github.com/namsyvo/IVC
```

### 2.2 Download IVC source code without Go
Get IVC source code with pre-compiled binary executable files of IVC (compiled on GNU/Linux 3.2.0-4-amd64 #1 SMP Debian 3.2.68-1+deb7u3 x86_64):
```
git clone https://github.com/namsyvo/IVC.git
cd IVC
```
Those binary executable files were obtained by compiling source code with Go:
```
go build main/ivc-index.go 
go build main/ivc.go
```
The source code can be also downloaded from releases of IVC at https://github.com/namsyvo/IVC/releases


3. Usage
--------

### 3.1 Example command
IVC comes with a test dataset which includes the following directories:   
./test_data/refs: includes a reference genome and a corresponding variant profile for NC_007194.1 (Aspergillus fumigatus Af293 chromosome 1, whole genome shotgun sequence, see http://www.ncbi.nlm.nih.gov/nuccore/AAHF00000000)   
./test_data/reads: includes a set of 10.000 simulated paired-end reads generated with DWGSIM (see https://github.com/nh13/DWGSIM).

3.1.1. Creating and indexing reference genomes with variant profile:
```
go run main/ivc-index.go -R test_data/refs/chr1_ref.fasta -V test_data/refs/chr1_variant_prof.vcf -I test_data/indexes
```
The command "go run main/ivc-index.go" can be replaced by the command "./ivc-index.go".

3.1.2. Calling variants from reads and the reference

```
go run main/ivc.go -R test_data/refs/chr1_ref.fasta -V test_data/refs/chr1_variant_prof.vcf -I test_data/indexes -1 test_data/reads/chr1_dwgsim_100_0.001-0.01.bwa.read1.fastq -2 test_data/reads/chr1_dwgsim_100_0.001-0.01.bwa.read2.fastq -O test_data/results/chr1_variant_calls.vcf
```
The command "go run main/ivc.go" can be replaced by the command "./ivc.go".

### 3.2 Commands and options

3.2.1. Creating and indexing reference genomes with variant profile:

Required:

	-R: reference genome (FASTA format).  
	-V: known variant profile (VCF format).  
	-I: directory for storing index.  

Options:


3.2.2. Calling Variants:

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


4. Preparing data and performing experiments
--------------------------------------------

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


5. Contact
----------

Nam Sy Vo  
nsvo1@memphis.edu