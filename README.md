ISC - An Integrated SNP Caller
===


1. Overview
-----------


2. Install ISC
--------------
Pre-requirement: GO environment is already set up properly.

```
go get github.com/namsyvo/ISC
cd $GOPATH/github.com/namsyvo/ISC
```

3. Usage
--------

### 3.1 Example command

3.1.1. Creating and indexing a reference multi-genome

```
go run main/index.go -g test_data/refs/chr1.fasta -p test_data/refs/vcf_chr_1.vcf -s test_data/indexes/SNPLocation.txt -m test_data/indexes/multigenome.txt -i test_data/indexes/multigenome.txt -r test_data/indexes/multigenome_rev.txt
```

3.1.2. Calling SNPs from reads and the reference

```
go run main/isc.go -g test_data/indexes/multigenome.txt -s test_data/indexes/SNPLocation.txt -i test_data/indexes/multigenome.txt.index -r test_data/indexes/multigenome_rev.txt.index -1 test_data/reads/test_reads_1.fq -2 test_data/reads/test_reads_2.fq -c test_data/results/test_called_snps.vcf
```

### 3.2 Commands and options


### 3.3 Parameters


### 3.4 Data-related problems


3. Notes
--------


4. Contact
----------
+ nsvo1@memphis.edu
+ qmtran@memphis.edu
+ vphan@memphis.edu
