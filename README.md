ISC - An Integrated SNP Caller
==============================

1. Install ISC
---------------
```
go get github.com/namsyvo/ISC
```
2. Usage
--------

### 2.1 Istall ISC
```
cd $GOPATH/github.com/namsyvo/ISC
```

### 2.2 Commands

2.2.1. Creating and indexing a reference multi-genome:
```
go run main/index.go -g test\_data/refs/genome.txt -s test\_data/refs/ref\_dbSNP.vcf -i test\_data/indexes/genomestar.txt.index -r test\_data/indexes/genomestar_rev.txt.index
```

2.2.3. Calling SNPs from reads and the reference:
```
go run main/isc.go -g test\_test/refs/multigenome.txt -s data-test/refs/SNPLocation.txt -i data-test/indexes/genomestar.txt.index -r data-test/indexes/genomestar_rev.txt.index -q data-test/reads/reads\_test1.fastq -c data-test/results/called\_snp\_test1.vcf
```

### 2.3 Parameters:


3. Data and other problems:


4. Contact
----------
+ nsvo1@memphis.edu
+ qmtran@memphis.edu
+ vphan@memphis.edu
