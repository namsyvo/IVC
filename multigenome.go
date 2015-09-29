//-------------------------------------------------------------------------------------------------
// IVC: multigenome.go - Building multigenome by combining SNPs and INDELs from dbSNPs with the reference genome.
// Multigenome includes a multi-sequence (include standard bases and *) and a variant profile for variant locations.
// Copyright 2015 Nam Sy Vo.
//-------------------------------------------------------------------------------------------------

package ivc

import (
	"bufio"
	"bytes"
	"fmt"
	"github.com/vtphan/fmi" //to use FM index
	"io/ioutil"
	"log"
	"os"
	"sort"
	"strconv"
	"strings"
)

//--------------------------------------------------------------------------------------------------
//MultiGenome represents multi-genome.
//--------------------------------------------------------------------------------------------------
type MultiGenome struct {
	Seq        []byte            //store "starred" sequence of the standard reference genome.
	Variants   map[int][][]byte  //store variants (position, variants).
	VarAF      map[int][]float32 //store allele frequency of variants (position, allele frequency).
	SameLenVar map[int]int       //indicate if variants has same length (SNPs or MNPs).
	DelVar     map[int]int       //store length of deletions if variants are deletion.
	RevFMI     fmi.Index         //FM-index of reverse multigenomes (to do forward search).
}

//--------------------------------------------------------------------------------------------------
// MultiGenome creates a multi-genome instance and set up its content.
//--------------------------------------------------------------------------------------------------
func NewMultiGenome() *MultiGenome {
	log.Printf("Memstats (golang name):\tAlloc\tTotalAlloc\tSys\tHeapAlloc\tHeapSys")

	M := new(MultiGenome)
	M.Seq = LoadMultiSeq(PARA_INFO.Ref_file)
	PrintMemStats("Memstats after loading multi-sequence")

	M.Variants, M.VarAF = LoadVarProf(PARA_INFO.Var_prof_file)
	PrintMemStats("Memstats after loading variant profile")

	M.SameLenVar = make(map[int]int)
	M.DelVar = make(map[int]int)
	var same_len_flag, del_flag bool
	var var_len int
	for var_pos, var_bases := range M.Variants {
		var_len = len(var_bases[0])
		same_len_flag, del_flag = true, true
		for _, val := range var_bases[1:] {
			if var_len != len(val) {
				same_len_flag = false
			}
			if var_len <= len(val) {
				del_flag = false
			}
		}
		if same_len_flag {
			M.SameLenVar[var_pos] = var_len
		}
		if del_flag {
			M.DelVar[var_pos] = var_len - 1
		}
	}
	PrintMemStats("Memstats after creating auxiliary data structures for multi-genome")

	M.RevFMI = *fmi.Load(PARA_INFO.Rev_index_file)
	PrintMemStats("Memstats after loading index of reverse multi-sequence")

	return M
}

type VarProf struct {
	Variant [][]byte
	AleFreq []float32
}

//-------------------------------------------------------------------------------------------------
// BuildMultiGenome builds multi-sequence from a standard reference genome and a variant profile.
//-------------------------------------------------------------------------------------------------
func BuildMultiGenome(genome_file, var_prof_file string) ([]byte, map[int]VarProf) {
	genome := ReadFASTA(genome_file)
	PrintProcessMem("Memstats after reading reference genome")
	var_prof := ReadVCF(var_prof_file)
	PrintProcessMem("Memstats after reading variant profile")
	multi_seq := make([]byte, len(genome))
	copy(multi_seq, genome)
	for key, _ := range var_prof {
		multi_seq[key] = '*'
	}
	return multi_seq, var_prof
}

//-------------------------------------------------------------------------------------------------
// LoadMultiSeq loads multi-sequence from file.
//-------------------------------------------------------------------------------------------------
func LoadMultiSeq(file_name string) []byte {
	bs, err := ioutil.ReadFile(file_name)
	if err != nil {
		fmt.Println("Error: Read multigenome file", err)
		os.Exit(1)
	}
	return bs
}

//-------------------------------------------------------------------------------------------------
// SaveMultiSeq saves multi-sequence to file.
//-------------------------------------------------------------------------------------------------
func SaveMultiSeq(file_name string, multi_seq []byte) {
	file, err := os.Create(file_name)
	if err != nil {
		fmt.Println("Error: Create multi-sequence file", err)
		os.Exit(1)
	}
	defer file.Close()
	file.Write(multi_seq)
}

//-------------------------------------------------------------------------------------------------
// LoadVarProf loads variant profile from file and return a map of variants.
//-------------------------------------------------------------------------------------------------
func LoadVarProf(file_name string) (map[int][][]byte, map[int][]float32) {
	Variant := make(map[int][][]byte)
	AleFreq := make(map[int][]float32)

	f, err := os.Open(file_name)
	if err != nil {
		fmt.Println("Error: Open variant profile file", err)
		os.Exit(1)
	}
	defer f.Close()
	br := bufio.NewReader(f)
	for {
		line, err := br.ReadString('\n')
		if err != nil {
			break
		}
		sline := string(line[:len(line)-1])
		split := strings.Split(sline, "\t")
		k, _ := strconv.ParseInt(split[0], 10, 64)
		t := make([]string, (len(split)-1)/2)
		p := make([]float32, (len(split)-1)/2)
		for i := 0; i < len(t); i++ {
			t[i] = split[2*i+1]
			tmp, _ := strconv.ParseFloat(split[2*i+2], 32)
			p[i] = float32(tmp)
		}

		// convert to [][]byte & map[int]int
		b := make([][]byte, len(t))

		for i := range b {
			b[i] = make([]byte, len(t[i]))
			copy(b[i], []byte(t[i]))
		}
		Variant[int(k)] = b
		AleFreq[int(k)] = p
	}
	return Variant, AleFreq
}

//-------------------------------------------------------------------------------------------------
// SaveVarProf saves variant profile to file.
//-------------------------------------------------------------------------------------------------
func SaveVarProf(file_name string, var_prof map[int]VarProf) {
	file, err := os.Create(file_name)
	if err != nil {
		fmt.Println("Error: Create variant profile index file", err)
		os.Exit(1)
	}
	defer file.Close()
	var var_pos []int
	for pos, _ := range var_prof {
		var_pos = append(var_pos, pos)
	}
	sort.Sort(sort.IntSlice(var_pos))
	for _, pos := range var_pos {
		_, err = file.WriteString(strconv.Itoa(pos) + "\t")
		for idx, val := range var_prof[pos].Variant {
			_, err = file.WriteString(string(val) + "\t" +
				strconv.FormatFloat(float64(var_prof[pos].AleFreq[idx]), 'f', 10, 32) + "\t")
		}
		_, err = file.WriteString("\n")
		if err != nil {
			break
		}
	}
}

//--------------------------------------------------------------------------------------------------
// ReadFASTA reads reference genome from FASTA files.
//--------------------------------------------------------------------------------------------------
func ReadFASTA(sequence_file string) []byte {
	f, err := os.Open(sequence_file)
	if err != nil {
		fmt.Println("Error: Open FASTA file", err)
		os.Exit(1)
	}
	defer f.Close()
	br := bufio.NewReader(f)
	byte_array := bytes.Buffer{}

	var isPrefix bool
	_, isPrefix, err = br.ReadLine()
	if err != nil || isPrefix {
		fmt.Printf("%v\n", err)
		os.Exit(1)
	}
	var line []byte
	for {
		line, isPrefix, err = br.ReadLine()
		if err != nil || isPrefix {
			break
		} else {
			byte_array.Write(line)
		}
	}
	return []byte(byte_array.String())
}

//--------------------------------------------------------------------------------------------------
// ReadVCF reads variant profile from VCF files.
//--------------------------------------------------------------------------------------------------
func ReadVCF(file_name string) map[int]VarProf {
	var_prof := make(map[int]VarProf)
	f, err := os.Open(file_name)
	if err != nil {
		fmt.Println("Error: Open VCF file", err)
		os.Exit(1)
	}
	defer f.Close()

	var line, sub_info []byte
	var alt, info, sub_info_part [][]byte
	var pos int
	var p float32

	data := bufio.NewReader(f)
	for {
		line, err = data.ReadBytes('\n')
		if err != nil {
			break
		}
		if bytes.Equal(line[0:1], []byte("#")) {
			continue
		} else {
			var_prof_elem := VarProf{}
			sub_line, _ := SplitN(line, []byte("\t"), 9)
			pos, _ = strconv.Atoi(string(sub_line[1]))
			var_prof_elem.Variant = append(var_prof_elem.Variant, sub_line[3])
			var_prof_elem.AleFreq = append(var_prof_elem.AleFreq, 0)
			alt = bytes.Split(sub_line[4], []byte(","))
			info = bytes.Split(sub_line[7], []byte(";"))
			for _, sub_info = range info {
				sub_info_part = bytes.Split(sub_info, []byte("="))
				if bytes.Equal(sub_info_part[0], []byte("AF")) {
					tmp_p, _ := strconv.ParseFloat(string(sub_info_part[1]), 32)
					p = float32(tmp_p)
					break
				}
			}
			for i := 0; i < len(alt); i++ {
				var_prof_elem.Variant = append(var_prof_elem.Variant, alt[i])
				var_prof_elem.AleFreq = append(var_prof_elem.AleFreq, p)
			}
			var_prof_elem.AleFreq[0] = 1 - p
			var_prof[pos-1] = var_prof_elem
		}
	}
	return var_prof
}
