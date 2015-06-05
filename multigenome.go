//-------------------------------------------------------------------------------------------------
// Multigenome package: multigenome module.
// Combine SNPs and INDELs from dbSNPs (vcf file) with a reference genome (fasta file).
// Copyright 2015 by Nam Sy Vo.
//-------------------------------------------------------------------------------------------------

package isc

import (
	"bufio"
	"bytes"
	"fmt"
	"io/ioutil"
	"os"
	"strconv"
	"strings"
	"sort"
	//"log"
)

type VarProf struct {
	Variant    	[][]byte
	AleFreq 	[]float32
}

func LoadVarProf(file_name string) (map[int][][]byte, map[int][]float32) {
	Variant := make(map[int][][]byte)
	AleFreq := make(map[int][]float32)

	f, err := os.Open(file_name)
	if err != nil {
		fmt.Println("Error: Load Variant Profile", err)
		os.Exit(1)
	}
	defer f.Close()
	br := bufio.NewReader(f)
	for {
		line, err := br.ReadString('\n')
		if err != nil {
			fmt.Println("Finish reading variant profile index")
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

func SaveVarProf(file_name string, var_prof map[int]VarProf) {
	file, err := os.Create(file_name)
	if err != nil {
		fmt.Println("Error: Save Variant Profile Index", err)
		return
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

func SaveMultigenome(file_name string, multigenome []byte) {
	file, err := os.Create(file_name)
	if err != nil {
		fmt.Println("Error: Save Multigenome", err)
		return
	}
	defer file.Close()
	file.Write(multigenome)
}

func LoadMultigenome(file_name string) []byte {
	bs, err := ioutil.ReadFile(file_name)
	if err != nil {
		fmt.Println("Error: Load Multigenome", err)
		return nil
	}
	return bs
}

func BuildMultigenome(var_prof map[int]VarProf, seq []byte) []byte {
	multigenome := make([]byte, len(seq))
	copy(multigenome, seq)
	for key, _ := range var_prof {
		multigenome[key] = '*'
	}
	return multigenome
}

//--------------------------------------------------------------------------------------------------
// Read FASTA files
//--------------------------------------------------------------------------------------------------
func ReadFASTA(sequence_file string) []byte {
	f, err := os.Open(sequence_file)
	if err != nil {
		fmt.Println("Error: Read FASTA file", err)
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
	//fmt.Printf("%s",line)
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
// Read VCF files
//--------------------------------------------------------------------------------------------------
func ReadVCF(file_name string) map[int]VarProf {
	var_prof := make(map[int]VarProf)
	f, err := os.Open(file_name)
	if err != nil {
		fmt.Println("Error: Read VCF file", err)
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
			//fmt.Printf("%v\n",err)
			break
		}
		if bytes.Equal(line[0:1], []byte("#")) {
			//fmt.Printf("%s \n",line)
		} else {
			tmp := VarProf{}
			sub_line, _ := SplitN(line, []byte("\t"), 9)
			pos, _ = strconv.Atoi(string(sub_line[1]))
			tmp.Variant = append(tmp.Variant, sub_line[3])
			//Temporal impl for allele freq, need to be changed later
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
			//sort.Strings(tmp.Variant)
			var_prof[pos - 1] = tmp
		}
	}
	return var_prof
}

func SplitN(s, sep []byte, n int) ([][]byte, int) {
	first_idx, sep_idx := 0, 0
	sep_num := 0
	t := make([][]byte, 0)
	for first_idx < len(s) {
		sep_idx = bytes.Index(s[first_idx:], sep)
		if sep_idx != -1 {
			sep_num++
			tmp := make([]byte, first_idx + sep_idx - first_idx)
			copy(tmp, s[first_idx : first_idx + sep_idx])
			t = append(t, tmp)
			if sep_num == n {
				return t, sep_num
			}
			first_idx = first_idx + sep_idx + 1
		} else {
			return t, sep_num
		}
	}
	return t, sep_num
}

func IndexN(s, sep []byte, n int) int {
	first_idx, sep_idx := 0, 0
	sep_num := 0
	for first_idx < len(s) {
		sep_idx = bytes.Index(s[first_idx:], sep)
		if sep_idx != -1 {
			sep_num++
			if sep_num == n {
				return first_idx + sep_idx
			}
			first_idx = first_idx + sep_idx + 1
		} else {
			return -1
		}
	}
	return -1
}
