//----------------------------------------------------------------------------------------
// Test for variant quality calculation
// Copyright 2015 Nam Sy Vo
//----------------------------------------------------------------------------------------

package ivc_test

import (
	"fmt"
	"github.com/namsyvo/IVC"
	"math"
	"testing"
)

func TestVarQual(t *testing.T) {
	defer __(o_())

	ivc.Q2C = make(map[byte]float64)
	ivc.Q2E = make(map[byte]float64)
	ivc.Q2P = make(map[byte]float64)
	var q byte
	for i := 33; i < 74; i++ {
		q = byte(i)
		//Phred-encoding factor (33) need to be estimated from input data
		ivc.Q2C[q] = -math.Log10(1.0 - math.Pow(10, -(float64(q)-33)/10.0))
		ivc.Q2E[q] = math.Pow(10, -(float64(q)-33)/10.0) / 3.0
		ivc.Q2P[q] = 1.0 - math.Pow(10, -(float64(q)-33)/10.0)
	}

	VC := new(ivc.VarCall)
	VC.VarProb = make(map[uint32]map[string]float64)
	VC.VarBQual = make(map[uint32]map[string][][]byte)
	VC.VarType = make(map[uint32]map[string]int)
	VC.VarRNum = make(map[uint32]map[string]int)
	VC.ChrDis = make(map[uint32]map[string][]int)
	VC.ChrDiff = make(map[uint32]map[string][]int)
	VC.AlnProb = make(map[uint32]map[string][]float64)
	VC.ChrProb = make(map[uint32]map[string][]float64)
	VC.MapProb = make(map[uint32]map[string][]float64)
	VC.StartPos1 = make(map[uint32]map[string][]int)
	VC.StartPos2 = make(map[uint32]map[string][]int)
	VC.Strand1 = make(map[uint32]map[string][]bool)
	VC.Strand2 = make(map[uint32]map[string][]bool)
	VC.ReadInfo = make(map[uint32]map[string][][]byte)

	VC.VarProb[100] = make(map[string]float64)
	VC.VarBQual[100] = make(map[string][][]byte)
	VC.VarType[100] = make(map[string]int)
	VC.VarRNum[100] = make(map[string]int)
	VC.ChrDis[100] = make(map[string][]int)
	VC.ChrDiff[100] = make(map[string][]int)
	VC.AlnProb[100] = make(map[string][]float64)
	VC.ChrProb[100] = make(map[string][]float64)
	VC.MapProb[100] = make(map[string][]float64)
	VC.StartPos1[100] = make(map[string][]int)
	VC.StartPos2[100] = make(map[string][]int)
	VC.Strand1[100] = make(map[string][]bool)
	VC.Strand2[100] = make(map[string][]bool)
	VC.ReadInfo[100] = make(map[string][][]byte)

	//**********************************
	assign_var := func(m, n int) []*ivc.VarInfo {
		var_info := make([]*ivc.VarInfo, m+n)
		for j := 0; j < m; j++ {
			var_info[j] = new(ivc.VarInfo)
			var_info[j].Pos, var_info[j].Bases, var_info[j].BQual, var_info[j].Type = 100, []byte("A"), []byte("I"), 0
		}
		for j := m; j < m+n; j++ {
			var_info[j] = new(ivc.VarInfo)
			var_info[j].Pos, var_info[j].Bases, var_info[j].BQual, var_info[j].Type = 100, []byte("C"), []byte("I"), 0
			//var_info[j].Pos, var_info[j].Bases, var_info[j].BQual, var_info[j].Type = 100, []byte("ACGT"), []byte("IIII"), 1
		}
		return var_info
	}
	//*******************************
	init_prob := func() {
		VC.VarProb[100]["A"] = 0.999969
		VC.VarProb[100]["C"] = 0.00001
		VC.VarProb[100]["G"] = 0.00001
		VC.VarProb[100]["T"] = 0.00001
		VC.VarProb[100]["ACGT"] = 0.000001
	}

	init_prob()
	fmt.Println("Initial")
	for v, p := range VC.VarProb[100] {
		fmt.Println("Var: ", string(v), "\tProb: ", p, "\tQual: ", -10*math.Log10(1-p))
	}
	fmt.Println()
	var_info := assign_var(1, 2)
	for _, var_item := range var_info {
		fmt.Println("Update for var", string(var_item.Bases))
		VC.UpdateVariantProb(var_item)
		for v, p := range VC.VarProb[100] {
			fmt.Println("Post var: ", string(v), "\tProb: ", p, "\tQual: ", -10*math.Log10(1-p))
		}
		fmt.Println()
	}
}

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
