//----------------------------------------------------------------------------------------
// Test for calculating distance
// Copyright 2014 Nam Sy Vo
//----------------------------------------------------------------------------------------

package isc_test

import (
	"fmt"
	"testing"
    "github.com/namsyvo/ISC"
)

type type_snpprofile map[int][][]byte
type type_samelensnp map[int]int
type TestCase struct {
	Profile type_snpprofile
	SNPlen type_samelensnp
	ref string
	read string
	qual string
	d float64
}
/*
// Test for alignment between reads and "starred" multi-genomes
func TestBackwarddistance2MultiAlignment(t *testing.T) {
	defer __(o_())

	var test_cases = []TestCase{
		{ type_snpprofile{}, type_samelensnp{}, "ACG", "G", 0 },
		{ type_snpprofile{}, type_samelensnp{}, "TTACG", "ACT", 1 },

		{ type_snpprofile{3: {{'A'}, {'C'}} }, type_samelensnp{3: 1}, "ACC*CGT", "ACCACGT", 0 },
		{ type_snpprofile{3: {{'A'}, {'C'}} }, type_samelensnp{3: 1}, "ACC*CGT", "ACCTCGT", isc.INF },
		{ type_snpprofile{3: {{'A'}, {'C'}} }, type_samelensnp{3: 1}, "ACC*CGT", "ACCCCGA", 1 },
		{ type_snpprofile{3: {{'A'}, {'C'}} }, type_samelensnp{3: 1}, "ACC*CGT", "CACGT", 0 },
		{ type_snpprofile{3: {{'A'}, {'C'}} }, type_samelensnp{3: 1}, "ACC*CGT", "CGCGT", isc.INF },
		{ type_snpprofile{3: {{'A'}, {'C'}} }, type_samelensnp{3: 1}, "ACC*CGT", "ACCTCGT", isc.INF },

		{ type_snpprofile{3: {{'A'}, {'C'}, {'.'}} }, 	type_samelensnp{}, "ACC*CGT", "ACCACGT", 0 },
		{ type_snpprofile{3: {{'A'}, {'C'}, {'.'}} }, 	type_samelensnp{}, "ACC*CGT", "ACCCCGT", 0 },
		{ type_snpprofile{3: {{'A'}, {'C'}, {'.'}} }, 	type_samelensnp{}, "ACC*CGT", "ACCCGT", 0 },
		{ type_snpprofile{3: {{'A'}, {'C'}, {'.'}} }, 	type_samelensnp{}, "ACC*CGT", "ACCTCGT", isc.INF },

		{ type_snpprofile{3: {{'A'}, {'A','C'}} }, type_samelensnp{}, "ACC*CGT", "ACCACGT", 0 },
		{ type_snpprofile{3: {{'A'}, {'A','C'}} }, type_samelensnp{}, "ACC*CGT", "ACCACCGT", 0 },
		{ type_snpprofile{3: {{'A'}, {'A','C'}} }, type_samelensnp{}, "ACC*CGT", "ACCATCGT", isc.INF },

		{ type_snpprofile{3: {{'A'}, {'A','C'}} }, type_samelensnp{}, "ACC*CGT", "ACCTCGT", isc.INF },
		{ type_snpprofile{3: {{'A'}, {'A','C'}} }, type_samelensnp{}, "ACC*CGT", "ACCACCGT", 0 },
		{ type_snpprofile{3: {{'A'}, {'A','C'}} }, type_samelensnp{}, "ACC*CGT", "ACCTCCGT", isc.INF },

		{ type_snpprofile{3: {{'T'}, {'T', 'T', 'A'}} }, type_samelensnp{}, "ACC*CGT", "ACCTCGT", 0 },
		{ type_snpprofile{3: {{'T'}, {'T', 'T', 'A'}} }, type_samelensnp{}, "ACC*CGT", "ACCTTACGT", 0 },
		{ type_snpprofile{3: {{'T'}, {'T', 'T', 'A'}} }, type_samelensnp{}, "ACC*CGT", "ACCGTACGT", isc.INF },
		{ type_snpprofile{3: {{'T'}, {'T', 'T', 'A'}} }, type_samelensnp{}, "ACC*CGT", "ACCTACGT", isc.INF },
		{ type_snpprofile{3: {{'T'}, {'T', 'T', 'A'}} }, type_samelensnp{}, "ACC*CGT", "ACCTTCGT", isc.INF },
		{ type_snpprofile{3: {{'T'}, {'T', 'T', 'A'}} }, type_samelensnp{}, "ACC*CGT", "ACCTTCGG", isc.INF },

		{ type_snpprofile{3: {{'A'}, {'C'}} }, type_samelensnp{3: 1}, "ACC*CGT", "TTTACCACGT", isc.INF },
		{ type_snpprofile{3: {{'A'}, {'C'}} }, type_samelensnp{3: 1}, "ACC*CGT", "TTTACCACGT", isc.INF },

		{ type_snpprofile{3: {{'A'}, {'T', 'A'}, {'T', 'T', 'A'}} }, type_samelensnp{}, "ACC*CGT", "ACCGTACGT", isc.INF },
		{ type_snpprofile{3: {{'A'}, {'T', 'A'}, {'T', 'T', 'A'}, {'.'}} }, type_samelensnp{}, "ACC*CGT", "ACCTTACGT", 0 },
		{ type_snpprofile{3: {{'A'}, {'T', 'A'}, {'T', 'T', 'A'}, {'.'}} }, type_samelensnp{}, "ACC*CGT", "ACCGTACGT", isc.INF },

		{ type_snpprofile{0: {{'A'}, {'T', 'A'}, {'T', 'T', 'A'}, {'.'}} }, type_samelensnp{}, "*ACGT", "TTAACGT", 0 },
		{ type_snpprofile{0: {{'A'}, {'T', 'A'}, {'T', 'T', 'A'}, {'.'}} }, type_samelensnp{}, "*ACGT", "GTAACGT", isc.INF },
		{ type_snpprofile{0: {{'A'}, {'T', 'A'}, {'T', 'T', 'A'}, {'.'}} }, type_samelensnp{}, "*ACGT", "ATACGT", isc.INF },

		{ type_snpprofile{5: {{'A'}, {'T', 'A'}, {'T', 'T', 'A'}, {'.'}} }, type_samelensnp{},
		 "TAACC*CGT", "ACCGTACGT", 2},

		{ type_snpprofile{7: {{'A'}, {'T', 'A'}, {'T', 'T', 'A'}, {'.'}} }, type_samelensnp{}, "CCCACGT*", "ACGTA", 0 },

	}
	for i := 0; i < len(test_cases); i++ {
		var index isc.Index
		//index.Init(genome_file, snp_file, index_file, rev_index_file, read_len, seq_err, k, a, n)
		//Init(DIST_THRES, test_cases[i].Profile, test_cases[i].SNPlen, 100)
		read, genome := []byte(test_cases[i].read), []byte(test_cases[i].genome)
		d, D, m, n, S, T, _ := index.BackwardDistance(read, genome, 0)
		if d + D != test_cases[i].d {
			t.Errorf("Fail alignment (case, read, genome, calculated distance2, true distance2, d, m, n):",
			 i, string(read), string(genome), d + D, test_cases[i].d, m, n)
		} else if d + D >= isc.INF {
			fmt.Println("Successful alignment but with infinity distance2 (distance2, read, genome, d, m, n, case):",
			 d + D, string([]byte(test_cases[i].read)), string([]byte(test_cases[i].genome)), m, n, i)
		} else {
			fmt.Println("Successful alignment (distance2, read, genome, d, m, n, case):",
			 d + D, string([]byte(test_cases[i].read)), string([]byte(test_cases[i].genome)), m, n, i)
			fmt.Println(index.BackwardTraceBack(read, genome, m, n, S, T, 0))
		}
	}
}


// Test for alignment between reads and "starred" multi-genomes
func TestForwarddistance2MultiAlignment(t *testing.T) {
	defer __(o_())

	var test_cases = []TestCase{

		{ type_snpprofile{}, type_samelensnp{}, "ACG", "A", 0 },
		{ type_snpprofile{}, type_samelensnp{}, "TTACG", "TTG", 1 },

		{ type_snpprofile{3: {{'A'}, {'C'}} }, type_samelensnp{3: 1}, "ACC*CGT", "ACCACGT", 0 },
		{ type_snpprofile{3: {{'A'}, {'C'}} }, type_samelensnp{3: 1}, "ACC*CGT", "ACCTCGT", isc.INF },
		{ type_snpprofile{3: {{'A'}, {'C'}} }, type_samelensnp{3: 1}, "ACC*CGT", "ACCCCGA", 1 },
		{ type_snpprofile{3: {{'A'}, {'C'}} }, type_samelensnp{3: 1}, "ACC*CGT", "ACCAC", 0 },
		{ type_snpprofile{3: {{'A'}, {'C'}} }, type_samelensnp{3: 1}, "ACC*CGT", "ACCGC", isc.INF },
		{ type_snpprofile{3: {{'A'}, {'C'}} }, type_samelensnp{3: 1}, "ACC*CGT", "ACCTCGT", isc.INF },

		{ type_snpprofile{3: {{'A'}, {'C'}, {'.'}} }, 	type_samelensnp{}, "ACC*CGT", "ACCACGT", 0 },
		{ type_snpprofile{3: {{'A'}, {'C'}, {'.'}} }, 	type_samelensnp{}, "ACC*CGT", "ACCCCGT", 0 },
		{ type_snpprofile{3: {{'A'}, {'C'}, {'.'}} }, 	type_samelensnp{}, "ACC*CGT", "ACCCGT", 0 },
		{ type_snpprofile{3: {{'A'}, {'C'}, {'.'}} }, 	type_samelensnp{}, "ACC*CGT", "ACCTCGT", isc.INF },

		{ type_snpprofile{3: {{'A'}, {'A','C'}} }, type_samelensnp{}, "ACC*CGT", "ACCACGT", 0 },
		{ type_snpprofile{3: {{'A'}, {'A','C'}} }, type_samelensnp{}, "ACC*CGT", "ACCACCGT", 0 },
		{ type_snpprofile{3: {{'A'}, {'A','C'}} }, type_samelensnp{}, "ACC*CGT", "ACCATCGT", isc.INF },

		{ type_snpprofile{3: {{'A'}, {'A','C'}} }, type_samelensnp{}, "ACC*CGT", "ACCTCGT", isc.INF },
		{ type_snpprofile{3: {{'A'}, {'A','C'}} }, type_samelensnp{}, "ACC*CGT", "ACCACCGT", 0 },
		{ type_snpprofile{3: {{'A'}, {'A','C'}} }, type_samelensnp{}, "ACC*CGT", "ACCTCCGT", isc.INF },

		{ type_snpprofile{3: {{'T'}, {'T', 'T', 'A'}} }, type_samelensnp{}, "ACC*CGT", "ACCTCGT", 0 },
		{ type_snpprofile{3: {{'T'}, {'T', 'T', 'A'}} }, type_samelensnp{}, "ACC*CGT", "ACCTTACGT", 0 },
		{ type_snpprofile{3: {{'T'}, {'T', 'T', 'A'}} }, type_samelensnp{}, "ACC*CGT", "ACCGTACGT", isc.INF },
		{ type_snpprofile{3: {{'T'}, {'T', 'T', 'A'}} }, type_samelensnp{}, "ACC*CGT", "ACCTACGT", isc.INF },
		{ type_snpprofile{3: {{'T'}, {'T', 'T', 'A'}} }, type_samelensnp{}, "ACC*CGT", "ACCTTCGT", isc.INF },
		{ type_snpprofile{3: {{'T'}, {'T', 'T', 'A'}} }, type_samelensnp{}, "ACC*CGT", "ACCTTCGG", isc.INF },

		{ type_snpprofile{3: {{'A'}, {'T', 'A'}, {'T', 'T', 'A'}} }, type_samelensnp{}, "ACC*CGT", "ACCGTACGT", isc.INF },
		{ type_snpprofile{3: {{'A'}, {'T', 'A'}, {'T', 'T', 'A'}, {'.'}} }, type_samelensnp{}, "ACC*CGT", "ACCTTACGT", 0 },
		{ type_snpprofile{3: {{'A'}, {'T', 'A'}, {'T', 'T', 'A'}, {'.'}} }, type_samelensnp{}, "ACC*CGT", "ACCGTACGT", isc.INF },

		{ type_snpprofile{0: {{'A'}, {'T', 'A'}, {'T', 'T', 'A'}, {'.'}} }, type_samelensnp{}, "*ACGT", "TTAACGT", 0 },
		{ type_snpprofile{0: {{'A'}, {'T', 'A'}, {'T', 'T', 'A'}, {'.'}} }, type_samelensnp{}, "*ACGT", "GTAACGT", isc.INF },
		{ type_snpprofile{0: {{'A'}, {'T', 'A'}, {'T', 'T', 'A'}, {'.'}} }, type_samelensnp{}, "*ACGT", "ATACGT", isc.INF },

		{ type_snpprofile{3: {{'A'}, {'T', 'A'}, {'T', 'T', 'A'}, {'.'}} }, type_samelensnp{}, "ACC*CGTAC", "ACCGTACGT", isc.INF },
		{ type_snpprofile{4: {{'A'}, {'T', 'A'}, {'T', 'T', 'A'}, {'.'}} }, type_samelensnp{}, "ACGT*GCCC", "ACGTAG", 0 },

	}
	for i := 0; i < len(test_cases); i++ {
		var index isc.Index
		//index.Init(genome_file, snp_file, index_file, rev_index_file, read_len, seq_err, k, a, n)
		//Init(DIST_THRES, test_cases[i].Profile, test_cases[i].SNPlen, 100)
		read, genome := []byte(test_cases[i].read), []byte(test_cases[i].genome)
		d, D, m, n, S, T, _ := index.ForwardDistance(read, genome, 0)
		if d + D != test_cases[i].d {
			t.Errorf("Fail alignment (read, genome, calculated distance2, true distance2, m, n, case):",
			 string(read), string(genome), d + D, test_cases[i].d, m, n, i)
		} else if d + D >= isc.INF {
			fmt.Println("Successful alignment but with infinity distance2 (distance2, read, genome, d, m, n, case):",
			 d + D, string([]byte(test_cases[i].read)), string([]byte(test_cases[i].genome)), m, n, i)
		} else {
			fmt.Println("Successful alignment (distance2, read, genome, d, m, n, case):",
			 d + D, string([]byte(test_cases[i].read)), string([]byte(test_cases[i].genome)), m, n, i)
			fmt.Println(index.ForwardTraceBack(read, genome, m, n, S, T, 0))
		}
	}
}

// Test for alignment between reads and "starred" multi-genomes
// Some more complex cases
func TestBackwarddistance2MultiAlignment2(t *testing.T) {
	defer __(o_())

	var test_cases = []TestCase{
		//test for 2 snp pos
		{ type_snpprofile{26042387: {{'G'}, {'G', 'C'}}, 26042385: {{'A'}, {'A', 'C', 'T'}} },
		 type_samelensnp{}, "CT*A*GGTTAAACAATTT", "AAGGGTTATTCAATTA", 3 },
	}
	for i := 0; i < len(test_cases); i++ {
		var index isc.Index
		//index.Init(genome_file, snp_file, index_file, rev_index_file, read_len, seq_err, k, a, n)
		//Init(DIST_THRES, test_cases[i].Profile, test_cases[i].SNPlen, 100)
		read, genome := []byte(test_cases[i].read), []byte(test_cases[i].genome)
		d, D, m, n, S, T, _ := index.BackwardDistance(read, genome, 26042383)
		if d + D != test_cases[i].d {
			t.Errorf("Fail alignment (read, genome, calculated distance2, true distance2, m, n, case):",
			 string(read), string(genome), d + D, test_cases[i].d, m, n, i)
		} else {
			fmt.Println("Successful alignment (distance2, read, genome, profile, m, n, case):",
				d + D, string([]byte(test_cases[i].read)), string([]byte(test_cases[i].genome)),
				 test_cases[i].Profile, m, n, i)
			snp := index.BackwardTraceBack(read, genome, m, n, S, T, 26042383)
			for k, v := range snp {
				fmt.Println(k, string(v))
			}
		}
	}
}

// Test for alignment between reads and "starred" multi-genomes
// Some more complex cases
func TestForwarddistance2MultiAlignment2(t *testing.T) {
	defer __(o_())

	var test_cases = []TestCase{
		//test for 2 snp pos
		{ type_snpprofile{26042385: {{'G'}, {'G', 'C'}}, 26042387: {{'A'}, {'A', 'C', 'T'}} },
		 type_samelensnp{}, "TTTAACAAATTGG*A*TC", "ATTAACTTATTGGGAA", 3 },
	}
	for i := 0; i < len(test_cases); i++ {
		var index isc.Index
		//index.Init(genome_file, snp_file, index_file, rev_index_file, read_len, seq_err, k, a, n)
		//Init(DIST_THRES, test_cases[i].Profile, test_cases[i].SNPlen, 100)
		read, genome := []byte(test_cases[i].read), []byte(test_cases[i].genome)
		d, D, m, n, S, T, _ := index.ForwardDistance(read, genome, 26042372)
		if d + D != test_cases[i].d {
			t.Errorf("Fail alignment (read, genome, calculated distance2, true distance2, m, n, case):",
			 string(read), string(genome), d + D, test_cases[i].d, m, n, i)
		} else {
			fmt.Println("Successful alignment (distance2, read, genome, profile, m, n, case):",
				d + D, string([]byte(test_cases[i].read)), string([]byte(test_cases[i].genome)),
				 test_cases[i].Profile, m, n, i)
			snp := index.ForwardTraceBack(read, genome, m, n, S, T, 26042372)
			for k, v := range snp {
				fmt.Println(k, string(v))
			}
		}
	}
}
*/
// Test for alignment between reads and "starred" multi-genomes
// Some more complex cases
func TestBackwardAffineGapEditDistance(t *testing.T) {
	defer __(o_())

	var test_cases = []TestCase{
		//test for 2 snp pos
		{ type_snpprofile{3230632: {{'A'}, {'A', 'C'}}, 3230636: {{'G'}, {'G', 'G', 'G', 'G', 'G', 'T'}}, 3230637: {{'G'}, {'C'}} },
		 type_samelensnp{3230637: 1}, "GAATGCCGTCCTTCCCC*CCG**GGGG", "CCTTCCCCACCCGGGGGGTGCGGGG", "@>=?>=>>>???@;@@>?>=>>?=>", 38.0 },
	}
	isc.INDEX.SNP_PROF, _, isc.INDEX.SAME_LEN_SNP = isc.LoadSNPLocation("/data/nsvo/test-data/GRCh37_chr1/indexes/af_mutant_index/index_0.70/isc_snp_prof_0.70.vcf.idx")
	isc.PARA_INFO = *isc.SetPara(100, 500, 700, 0.0015, 0.01)
	align_info := isc.InitAlignInfo(2 * isc.PARA_INFO.Read_len)
	var snp_prof isc.SNP_Prof
	for i := 0; i < len(test_cases); i++ {
		//index.Init(genome_file, snp_file, index_file, rev_index_file, read_len, seq_err, k, a, n)
		//Init(DIST_THRES, test_cases[i].Profile, test_cases[i].SNPlen, 100)
		read, qual, ref := []byte(test_cases[i].read), []byte(test_cases[i].qual), []byte(test_cases[i].ref)
		d, D, m, n, _, _, _ := snp_prof.BackwardDistance(read, qual, ref, 3230615, align_info.Bw_Dis, align_info.Bw_Trace)
		fmt.Println("Successful alignment (distance2, read, ref, profile, m, n, case):", d + D, string([]byte(test_cases[i].read)), string([]byte(test_cases[i].ref)), test_cases[i].Profile, m, n, i)
		l_pos, l_val, l_idx := snp_prof.BackwardTraceBack(read, qual, ref, m, n, 3230615, align_info.Bw_Dis, align_info.Bw_Trace)
		for i := 0; i < len(l_pos); i++ {
				fmt.Println(l_pos[i], string(l_val[i]), l_idx[i])
		}
	}
}
