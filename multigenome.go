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
	"io/ioutil"
	"os"
	"sort"
	"strconv"
	"strings"
)

type VarProf struct {
	Variant [][]byte
	AleFreq []float32
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

//-------------------------------------------------------------------------------------------------
// SaveMultigenome saves multi-sequence to file.
//-------------------------------------------------------------------------------------------------
func SaveMultigenome(file_name string, multigenome []byte) {
	file, err := os.Create(file_name)
	if err != nil {
		fmt.Println("Error: Create multigenome file", err)
		os.Exit(1)
	}
	defer file.Close()
	file.Write(multigenome)
}

//-------------------------------------------------------------------------------------------------
// LoadMultigenome loads multi-sequence from file.
//-------------------------------------------------------------------------------------------------
func LoadMultigenome(file_name string) []byte {
	bs, err := ioutil.ReadFile(file_name)
	if err != nil {
		fmt.Println("Error: Read multigenome file", err)
		os.Exit(1)
	}
	return bs
}

//-------------------------------------------------------------------------------------------------
// BuildMultigenome biulds multi-sequence from variant profile.
//-------------------------------------------------------------------------------------------------
func BuildMultigenome(var_prof map[int]VarProf, seq []byte) []byte {
	multigenome := make([]byte, len(seq))
	copy(multigenome, seq)
	for key, _ := range var_prof {
		multigenome[key] = '*'
	}
	return multigenome
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
			tmp := VarProf{}
			sub_line, _ := SplitN(line, []byte("\t"), 9)
			pos, _ = strconv.Atoi(string(sub_line[1]))
			tmp.Variant = append(tmp.Variant, sub_line[3])
			tmp.AleFreq = append(tmp.AleFreq, 0)
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
				tmp.Variant = append(tmp.Variant, alt[i])
				tmp.AleFreq = append(tmp.AleFreq, p)
			}
			tmp.AleFreq[0] = 1 - p
			var_prof[pos-1] = tmp
		}
	}
	return var_prof
}
