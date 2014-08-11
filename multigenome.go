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
	"sort"
	"strconv"
	"strings"
)

type SNPProfile struct {
	Profile    []string
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
			str = "\t" + v + str
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
	str := string(bs)
	return ([]byte(str))
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
	br := bufio.NewReader(f)
	var line string
	var split, alt []string
	t := make([]string, 1000)
	var pos int64

	for {
		line, err = br.ReadString('\n')
		if err != nil {
			//fmt.Printf("%v\n",err)
			break
		}
		if line[0] == byte('#') {
			//fmt.Printf("%s \n",line)
		} else {
			split = strings.Split(line, "\t")
			//fmt.Printf("%s %s %s\n", split[1], split[3], split[4])
			pos, _ = strconv.ParseInt(split[1], 10, 64)
			pos = pos - 1
			if len(split[4]) > 1 {
				alt = strings.Split(split[4], ",")
				//t := make([]string, len(alt)+1)
				t[0] = split[3]
				for i := 0; i < len(alt); i++ {
					println("len(alt): ", len(alt))
					if alt[i] == "<DEL>" {
						t[i+1] = "."
					} else {
println("in for: ", i+1)
						t[i+1] = alt[i]
					}
				}
				tmp, ok := array[int(pos)]
				if ok {
					t = append(t[:0], t[1:len(alt)+1]...)
					tmp.Profile = append(tmp.Profile, t...)
				} else {
					tmp.Profile = append(tmp.Profile, t...)
				}
				sort.Strings(tmp.Profile)
				array[int(pos)] = tmp // append SNP at pos
				//fmt.Printf("pos=%d %q \n", pos, alt)
			} else {
				tmp, ok := array[int(pos)]
				if ok {
					if split[4] == "<DEL>" {
						tmp.Profile = append(tmp.Profile, ".")
					} else {
						tmp.Profile = append(tmp.Profile, split[4])
					}
				} else {
					if split[4] == "<DEL>" {
						tmp.Profile = append(tmp.Profile, []string{split[3], "."}...)
					} else {
						tmp.Profile = append(tmp.Profile, []string{split[3], split[4]}...)
					}
				}
				sort.Strings(tmp.Profile)

				array[int(pos)] = tmp // append SNP at pos
				//fmt.Println(pos)
			}
		}
	}
	return array
}
