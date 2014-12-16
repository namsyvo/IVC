ISC - An Integrated SNP Caller
===


1. Overview
-----------


2. Install ISC
--------------
Pre-requirement: GO environment is already set up properly.  
Check Go path to make sure GOPATH is set up properly. For example:
```
echo $GOPATH
/home/nsvo/workspace/goprojects
```

Get fmi and ISC source code:
```
go get github.com/vtphan/fmi
go get github.com/namsyvo/ISC
```
After these steps, fmi source code and ISC source code should be in the directories $GOPATH/github.com/vtphan/fmi and $GOPATH/github.com/namsyvo/ISC, respectively.  
Then go to the ISC directory:
```
cd $GOPATH/src/github.com/namsyvo/ISC
```

3. Usage
--------

### 3.1 Example command
ISC comes with a test dataset which includes the following directories:  
test_data/refs: includes a reference genome and the corresponding dbsnp (NC_007194.1: Aspergillus fumigatus Af293 chromosome 1, whole genome shotgun sequence).  
test_data/reads: includes two test reads files for above reference.

3.1.1. Creating and indexing reference genomes with SNP profile:
```
mkdir test_data/index
go run main/index.go -g test_data/refs/chr1.fasta -s test_data/refs/vcf_chr_1.vcf -i test_data/index
```

3.1.2. Calling SNPs from reads and the reference

```
mkdir test_data/results
go run main/isc.go -g test_data/refs/chr1.fasta -s test_data/refs/vcf_chr_1.vcf -i test_data/index -1 test_data/reads/test_reads_1.fq -2 test_data/reads/test_reads_2.fq -o test_data/results/test_called_snps.vcf
```

### 3.2 Commands and options

3.2.1. Creating and indexing reference genomes with SNP profile:

Required:

	-g: reference genome (FASTA format).  
	-s: SNP profile (dbSNP with VCF format).  
	-i: directory for storing index.  

Options:


3.2.2. Calling SNPs:

Required:

	-g: reference genome (FASTA format).  
	-s: SNP profile (dbSNP with VCF format).  
	-i: directory for storing index.  
	-1: the read file (for single-end reads) (FASTQ format).  
	-2: the second end file (for pair-end reads) (FASTQ format).  
	-o: called SNP result file (in VCF format).  

Options:  

	-m: searching mode for finding seeds (1: random, 2: deterministic; default: 1).  
	-p: starting position on reads for finding seeds (integer, default: 0).  
	-j: step for searching in deterministic mode (integer, default: 5).  
	-w: maximum number of CPUs using by Go (integer, default: number of CPU of runnign system).  
	-t: number of used goroutines (integer, default: number of CPU of runnign system).  
	-n: maximum number of seeds for each ends.  
	-l: minimum length of seeds for each ends.  
	-k: maximum number of paired-seeds which are considers as proper seeds for paired-end reads.  


### 3.3 Parameters


### 3.4 Data-related problems


3. Notes
--------


4. Contact:
-----------
Nam Sy Vo  
nsvo1@memphis.edu