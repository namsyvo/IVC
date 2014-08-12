//-------------------------------------------------------------------------------------------------
// Multigenome package: multigenome module.
// Combine SNPs and INDELs from dbSNPs (vcf file) with a reference genome (fasta file).
// Copyright 2014 Quang Minh Tran.
// Modified by Nam Sy Vo.
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
	//"sort"
)

type SNPProfile struct {
	Profile    [][]byte
	AlleleFreq []float32
}

func LoadSNPLocation(file_name string) (map[int][][]byte, map[int][]float32, map[int]int) {
	//location := make(map[int]SNPProfile)
	barr := make(map[int][][]byte)
	af := make(map[int][]float32)
	is_equal := make(map[int]int)

	f, err := os.Open(file_name)
	if err != nil {
		fmt.Printf("%v\n", err)
		os.Exit(1)
	}
	defer f.Close()
	br := bufio.NewReader(f)
	for {
		line, err := br.ReadString('\n')
		if err != nil {
			//fmt.Printf("%v\n",err)
			break
		}
		sline := string(line[:len(line)-1])
		split := strings.Split(sline, "\t")
		k, _ := strconv.ParseInt(split[0], 10, 64)
		t := make([]string, len(split)-1)
		for i := 0; i < len(t); i++ {
			t[i] = split[i+1]
		}
		//location[int(k)] = SNPProfile{t}

		// convert to [][]byte & map[int]int
		flag := len(t[0])
		b := make([][]byte, len(t))

		for i := range b {
			b[i] = make([]byte, len(t[i]))
			copy(b[i], []byte(t[i]))
			//b[i] = []byte(t[i])
			if flag != len(b[i]) || t[i] == "." {
				flag = 0
			}
		}
		barr[int(k)] = b

		//Read allele freq - testing///////////////////
		p := make([]float32, len(split)-1)
		for i := range p {
			//tmp, _ := strconv.ParseFloat(split[i + 1 + (len(split)-1)/2], 32)
			//p[i] = float32(tmp)
			p[i] = 1 / float32(len(p))
		}
		af[int(k)] = p
		//////////////////////////////////////////////

		if flag != 0 {
			is_equal[int(k)] = flag
		}
	}
	return barr, af, is_equal
}

func SaveSNPLocation(file_name string, SNP_arr map[int]SNPProfile) {
	file, err := os.Create(file_name)
	if err != nil {
		return
	}
	defer file.Close()
	for i, item := range SNP_arr {
		str := ""
		for _, v := range item.Profile {
			str = "\t" + string(v) + str
		}
		key := strconv.Itoa(i)
		_, err := file.WriteString(key + str + "\n")
		if err != nil {
			fmt.Println(err)
			break
		}
	}
}

func SaveMultigenome(file_name string, multi []byte) {
	file, err := os.Create(file_name)
	if err != nil {
		// handle the error here
		return
	}
	defer file.Close()
	file.Write(multi)
}

func LoadMultigenome(file_name string) []byte {
	bs, err := ioutil.ReadFile(file_name)
	if err != nil {
		return nil
	}
	return bs
}

// string * multi-genome
func BuildMultigenome(SNP_arr map[int]SNPProfile, seq []byte) []byte {
	multi := make([]byte, len(seq))
	copy(multi, seq)
	for key, _ := range SNP_arr {
		multi[key] = '*'
	}
	return multi
}

//--------------------------------------------------------------------------------------------------
// Read FASTA files
//--------------------------------------------------------------------------------------------------
func ReadFASTA(sequence_file string) []byte {
	f, err := os.Open(sequence_file)
	if err != nil {
		fmt.Printf("%v\n", err)
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
func ReadVCF(sequence_file string) map[int]SNPProfile {
	array := make(map[int]SNPProfile)
	f, err := os.Open(sequence_file)
	if err != nil {
		fmt.Printf("%v\n", err)
		os.Exit(1)
	}
	defer f.Close()

	var line []byte
	var alt [][]byte
	var pos int

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
			tmp := SNPProfile{}
			sub_line, _ := SplitN(line, []byte("\t"), 9)
			pos, _ = strconv.Atoi(string(sub_line[1]))
			tmp.Profile = append(tmp.Profile, sub_line[3])
			alt = bytes.Split(sub_line[4], []byte(","))
			for i := 0; i < len(alt); i++ {
				if bytes.Equal(alt[i], []byte("<DEL>")) {
					tmp.Profile = append(tmp.Profile, []byte("."))
				} else {
					tmp.Profile = append(tmp.Profile, alt[i])
				}
			}

			//sort.Strings(tmp.Profile)
			array[pos - 1] = tmp // append SNP at pos
		}
	}
	return array
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
