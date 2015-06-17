IVC - An Integrated Variant Caller.
===================================


1. Overview
-----------


2. Install IVC
--------------
Pre-requirement: GO environment is already set up properly.  
Check Go path to make sure GOPATH is set up properly. For example:
```
echo $GOPATH
/home/nsvo/workspace/goprojects
```

Get fmi and IVC source code:
```
go get github.com/vtphan/fmi
go get github.com/namsyvo/IVC
```
After these steps, fmi source code and IVC source code should be in the directories $GOPATH/github.com/vtphan/fmi and $GOPATH/github.com/namsyvo/IVC, respectively.  
Then go to the IVC directory:
```
cd $GOPATH/src/github.com/namsyvo/IVC
```

3. Usage
--------

### 3.1 Example command
IVC comes with a test dataset which includes the following directories:  
test_data/refs: includes a reference genome and the corresponding dbsnp (NC_007194.1: Aspergillus fumigatus Af293 chromosome 1, whole genome shotgun sequence).  
test_data/reads: includes two test reads files for above reference.

3.1.1. Creating and indexing reference genomes with variant profile:
```
mkdir test_data/index
go run main/index.go -g test_data/refs/chr1.fasta -s test_data/refs/vcf_chr_1.vcf -i test_data/index
```

3.1.2. Calling Variants from reads and the reference

```
mkdir test_data/results
go run main/ivc.go -g test_data/refs/chr1.fasta -v test_data/refs/vcf_chr_1.vcf -i test_data/index -1 test_data/reads/test_reads_1.fq -2 test_data/reads/test_reads_2.fq -o test_data/results/test_called_snps.vcf
```

### 3.2 Commands and options

3.2.1. Creating and indexing reference genomes with Variant profile:

Required:

	-g: reference genome (FASTA format).  
	-v: variant profile, like dbSNP (VCF format).  
	-i: directory for storing index.  

Options:


3.2.2. Calling Variants:

Required:

	-g: reference genome (FASTA format).  
	-v: variant profile, like dbSNP (VCF format).  
	-i: directory for storing index.  
	-1: the read file (for single-end reads) (FASTQ format).  
	-2: the second end file (for pair-end reads) (FASTQ format).  
	-o: variant call result file (in VCF format).  

Options:  

	-m: searching mode for finding seeds (1: random, 2: deterministic; default: 1).  
	-p: starting position on reads for finding seeds (integer, default: 0).  
	-j: step for searching in deterministic mode (integer, default: 5).  
	-w: maximum number of CPUs using by Go (integer, default: number of CPU of running computer).  
	-t: number of used goroutines (integer, default: number of CPU of running computer).  
	-n: maximum number of seeds for each end (default: 1024).  
	-l: minimum length of seeds for each end (default: 10).  
	-h: maximum length of seeds for each end (default: length of input reads).  
	-k: maximum number of paired-seeds which are considers as proper seeds for paired-end reads (default: 1024).  
	-d: threshold of alignment distances (default: determined by the program).  
	-r: maximum number of iterations (default: determined by the program).  



### 3.3 Parameter calculation


### 3.4 Data-related problems


3. Notes
--------


4. Contact:
-----------
Nam Sy Vo  
nsvo1@memphis.edu