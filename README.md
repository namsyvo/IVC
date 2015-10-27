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
Then go to the IVC directory, from which IVC can be run as a Go program:
```
cd $GOPATH/src/github.com/namsyvo/IVC
```

3. Usage
--------

### 3.1 Example command
IVC comes with a test dataset which includes the following directories:  
test_data/refs: includes a reference genome and a corresponding variant profile (NC_007194.1: Aspergillus fumigatus Af293 chromosome 1, whole genome shotgun sequence).  
test_data/reads: includes a set of paired-end reads with 10.000 simulated reads.

3.1.1. Creating and indexing reference genomes with variant profile:
```
mkdir test_data/index
go run main/ivc-index.go -R test_data/refs/chr1_ref.fasta -V test_data/refs/chr1_variant_prof.vcf -I test_data/indexes
```

3.1.2. Calling Variants from reads and the reference

```
mkdir test_data/results
go run main/ivc.go -R test_data/refs/chr1_ref.fasta -V test_data/refs/chr1_variant_prof.vcf -I test_data/indexes/ -1 test_data/reads/chr1_reads_1.fq -2 test_data/reads/chr1_reads_2.fq -O test_data/results/chr1_variant_calls.vcf
```

### 3.2 Commands and options

3.2.1. Creating and indexing reference genomes with Variant profile:

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
	-O: variant call result file (in VCF format).  

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


### 3.3 Parameter calculation


### 3.4 Data-related problems


3. Notes
--------


4. Contact:
-----------
Nam Sy Vo  
nsvo1@memphis.edu