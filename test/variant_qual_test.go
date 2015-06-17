//----------------------------------------------------------------------------------------
// Test for variant quality calculation
// Copyright 2015 Nam Sy Vo
//----------------------------------------------------------------------------------------

package ivc_test

import (
	//"fmt"
	//"testing"
    //"github.com/namsyvo/IVC"
)
/*
func TestSNPQual(t *testing.T) {
	defer __(o_())

	var S ivc.SNP_Prof
	S.SNP_Calls = make(map[uint32]map[string]float64)
	S.SNP_Calls[100] = make(map[string]float64)
	S.SNP_Calls[100]["A"] = 0.97
	S.SNP_Calls[100]["C"] = 0.01
	S.SNP_Calls[100]["G"] = 0.01
	S.SNP_Calls[100]["T"] = 0.01
	//S.SNP_Calls[100]["ACT"] = 0.1
	
	snps := []ivc.SNP {
		{100, []byte{'C'}, []byte{'I'}},
		//{100, []byte{'C'}, []byte{'1'}},
		//{100, []byte{'A'}, []byte{'5'}},
		//{100, []byte{'A'}, []byte{'I'}},
		//{100, []byte{'C'}, []byte{'5'}},
		//{100, []byte{'T'}, []byte{'I'}},
		//{100, []byte{'T'}, []byte{'I'}},
	}

	snps2 := []ivc.SNP {
		{100, []byte{'A'}, []byte{'I'}},
		{100, []byte{'A'}, []byte{'I'}},
		{100, []byte{'G'}, []byte{'I'}},
		{100, []byte{'G'}, []byte{'I'}},
		{100, []byte{'G'}, []byte{'I'}},
		{100, []byte{'G'}, []byte{'I'}},
		{100, []byte{'G'}, []byte{'I'}},
	}
	snps3 := []ivc.SNP {
		{100, []byte{}, []byte{}},
		{100, []byte{'A', 'G', 'T'}, []byte{'4', '4', '4'}},
		{100, []byte{'A', 'G', 'T'}, []byte{'4', '4', '4'}},
		{100, []byte{'A', 'C', 'T'}, []byte{'4', '4', '4'}},
		{100, []byte{'A', 'C', 'T'}, []byte{'4', '4', '4'}},
		{100, []byte{'A', 'C', 'T'}, []byte{'4', '4', '4'}},
		{100, []byte{'A', 'C', 'T'}, []byte{'4', '4', '4'}},
	}


	for snp, prob := range S.SNP_Calls[100] {
		fmt.Print("SNP: ", snp, "\t")
		fmt.Println("Prob: ", prob, "\tQual: ", ivc.ProbtoQual(prob))
	}
	for _, snp := range snps {
		fmt.Println("a: ", string(snp.Bases))
		fmt.Println("e: ", string(snp.BaseQ))
		S.UpdateSNPProb(snp)
		for snp, prob := range S.SNP_Calls[100] {
			fmt.Print("SNP: ", snp, "\t")
			fmt.Println("Prob: ", prob, "\tQual: ", ivc.ProbtoQual(prob))
		}
	}
}


var bases []byte = []byte{'A', 'C', 'G', 'T'}
var p_b []float64 = []float64{0.7, 0.2, 0.00001, 0.09999}

func TestSNPQual(t *testing.T) {
	defer __(o_())

	fmt.Println("p_b: ")
	for _, p := range p_b {
		fmt.Print(p, " ")
	}
	fmt.Println()

	var a [][]byte
	var e [][]byte
	a = append(a, []byte{'A', 'A', 'A', 'A', 'C', 'T', 'T'})
	e = append(e, []byte{'4', '4', '4', '4', '4', '4', '4'})
	a = append(a, []byte{'A', 'A', 'C', 'C', 'T', 'T'})
	e = append(e, []byte{'4', '4', '4', '4', '4', '4'})
	a = append(a, []byte{'A', 'C', 'C', 'C', 'T', 'T'})
	e = append(e, []byte{'4', '4', '4', '4', '4', '4'})
	a = append(a, []byte{'A', 'C', 'C', 'T', 'T', 'T'})
	e = append(e, []byte{'4', '4', '4', '4', '4', '4'})
	a = append(a, []byte{'A', 'C', 'C', 'T', 'T'})
	e = append(e, []byte{'4', '4', '4', '4', '4'})
	a = append(a, []byte{'A', 'A', 'C', 'T', 'T'})
	e = append(e, []byte{'4', '4', '4', '4', '4'})
	a = append(a, []byte{'A', 'A', 'A', 'A', 'A', 'T'})
	e = append(e, []byte{'4', '4', '4', '4', '4', '4'})
	a = append(a, []byte{'G', 'G', 'G', 'G', 'G', 'T'})
	e = append(e, []byte{'4', '4', '4', '4', '4', '4'})

	fmt.Println("a: ", string(a[0]))
	fmt.Println("e: ", string(e[0]))
	for j, v := range bases {
		fmt.Print("SNP: ", string(v), "\t")
		fmt.Println("Prob: ", p_b[j], "\tQual: ", ivc.ProbtoQual(p_b[j]))
	}
	for i, _ := range(a[0]) {
		fmt.Println("base: ", string(a[0][i]))
		snp_prob := ivc.CalcPosProb(a[0][i], e[0][i], bases, p_b)
		copy(p_b, snp_prob)
		for j, v := range bases {
			fmt.Print("SNP: ", string(v), "\t")
			fmt.Println("Prob: ", p_b[j], "\tQual: ", ivc.ProbtoQual(p_b[j]))
		}
	}
}
*/