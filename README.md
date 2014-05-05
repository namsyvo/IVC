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
go run main/index.go
 -g test_data/refs/genome.txt -s test_data/refs/ref_dbSNP.vcf
 -i test_data/indexes/genomestar.txt.index -r test_data/indexes/genomestar_rev.txt.index
```

3.1.2. Calling SNPs from reads and the reference

```
go run main/isc.go
 -g test_data/refs/multigenome.txt -s test_data/refs/SNPLocation.txt
 -i test_data/indexes/genomestar.txt.index -r test_data/indexes/genomestar_rev.txt.index
  -q test_data/reads/test_reads.fastq -c test_data/results/test_called_snps.vcf
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
