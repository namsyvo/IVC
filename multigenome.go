//-------------------------------------------------------------------------------------------------
// Multigenome package: multigenome module.
// Combine SNPs and INDELs from dbSNPs (vcf file) with a reference genome (fasta file).
// Copyright 2014 Quang Minh Tran.
// Modified by Nam Sy Vo.
//-------------------------------------------------------------------------------------------------

package isc

import (
	"bufio"
	"io/ioutil"
	"fmt"
	"os"
	"strings"
	"strconv"
)

type SNP struct{
	Profile []string
	AlleleFreq []float32
}

func LoadSNPLocation(file_name string ) (map[int][][]byte, map[int][]float32, map[int]int) {
	//location := make(map[int]SNP)
	barr := make(map[int][][]byte)
	af := make(map[int][]float32)
	is_equal := make(map[int]int)
	
	f,err := os.Open(file_name)
    if err != nil{
        fmt.Printf("%v\n",err)
        os.Exit(1)
    }
	defer f.Close()
    br := bufio.NewReader(f)
	for{
		line , err := br.ReadString('\n')
		if err != nil {
			//fmt.Printf("%v\n",err)
			break
		}
		sline := string(line[:len(line)-1])
		split := strings.Split(sline, "\t");
		k, _ := strconv.ParseInt(split[0], 10, 64)
		t := make([]string, (len(split)-1)/2)
		for i := 0; i < len(t); i++ {
			t[i] = split[i + 1]
		}
		//location[int(k)] = SNP{t} 
		
		// convert to [][]byte & map[int]int
		flag := len(t[0]);
		b := make([][]byte, len(t))

		for i:= range b {
			b[i] = make([]byte, len(t[i]))
			copy(b[i], []byte(t[i]))
			//b[i] = []byte(t[i])
			if (flag != len(b[i]) || t[i] == ".") {
				flag = 0;
			}
		}		
		barr[int(k)] = b

		//Read allele freq - testing///////////////////
		p := make([]float32, (len(split)-1)/2)
		for i:= range p {
			tmp, _ := strconv.ParseFloat(split[i + 1 + (len(split)-1)/2], 32)
			p[i] = float32(tmp)
		}		
		af[int(k)] = p
		//////////////////////////////////////////////

		if flag != 0 {
			is_equal[int(k)] = flag			
		}
	}
	return barr, af, is_equal
}

func SaveSNPLocation(file_name string , SNP_arr map[int]SNP) {
	file, err := os.Create(file_name)
	if  err != nil {
        return
    }
	defer file.Close()
	for i, item := range SNP_arr {
		str := ""
		for _, v := range item.Profile {
			str = "\t" + v + str
		}
		key := strconv.Itoa(i)
		_, err := file.WriteString(key + str + "\n"); 
		if err != nil {
            fmt.Println(err)
            break
        }
	}	
}

func SaveMultigenome(file_name string , multi []byte) {
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
func BuildMultigenome(SNP_arr map[int]SNP, seq []byte) []byte {
	multi := make([]byte, len(seq))
	copy(multi, seq)
	for key, _ := range SNP_arr {
		 multi[key] = '*'
	}
	return multi
}
