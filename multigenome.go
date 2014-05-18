//-------------------------------------------------------------------------------------------------
// Multigenome package: multigenome module.
// Combine SNPs and INDELs from dbSNPs (vcf file) with a reference genome (fasta file).
// Copyright 2014 Quang Minh Tran.
// Modified by Nam Sy Vo.
//-------------------------------------------------------------------------------------------------

package multigenome

import (
	"bufio"
	"io/ioutil"
	"fmt"
	"os"
	"strings"
	"strconv"
	"bytes"
	"sort"
)

type SNP struct{
	profile []string
}

func LoadSNPLocation(file_name string )  (map[int] [][]byte, map[int]int) {
	//location := make(map[int]SNP)
	barr := make(map[int][][]byte)
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
		t := make([]string, len(split)-1)
		for i := 1; i<len(split); i++ {
			t[i-1] = split[i]
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
		if flag != 0 {
			is_equal[int(k)] = flag			
		}
	}
	return barr, is_equal
}

func SaveSNPLocation(file_name string , SNP_arr map[int]SNP) {
	file, err := os.Create(file_name)
	if  err != nil {
        return
    }
	defer file.Close()
	for i, item := range SNP_arr {
		str := ""
		for _, v := range item.profile {
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

func SaveMulti(file_name string , multi []byte) {
	file, err := os.Create(file_name)
    if err != nil {
        // handle the error here
        return
    }
    defer file.Close()
    file.Write(multi)
}

func LoadMulti(file_name string) []byte {
	bs, err := ioutil.ReadFile(file_name)
    if err != nil {
        return nil
    }
    str := string(bs)	
	return ([]byte(str))
}

// string * multi-genome
func buildMultigenome2(SNP_arr map[int]SNP, seq []byte) []byte {
	multi := make([]byte, len(seq))
	copy(multi, seq)
	for key, _ := range SNP_arr {
		 multi[key] = '*'
	}
	return multi
}

func vcfRead(sequence_file string) map[int]SNP {
	array := make(map[int]SNP)
	f,err := os.Open(sequence_file)
    if err != nil{
        fmt.Printf("%v\n",err)
        os.Exit(1)
    }

    defer f.Close()
    br := bufio.NewReader(f)
    //byte_array := bytes.Buffer{}

	for{
		line , err := br.ReadString('\n')
		if err != nil {
			//fmt.Printf("%v\n",err)
			break
		}		
		if line[0]==byte('#') {
			//fmt.Printf("%s \n",line)
		} else {
			sline := string(line)
			split := strings.Split(sline, "\t");
			//fmt.Printf("%s %s %s\n", split[1], split[3], split[4])
			pos, _ := strconv.ParseInt(split[1], 10, 64)
			pos = pos - 1
			if len(split[4])>1 {
				alt := strings.Split(split[4], ",")
				t := make([]string, len(alt)+1)
				t[0] = split[3]				
				for i:=0; i<len(alt); i++ {
					if alt[i] == "<DEL>" {
						t[i+1] = "."
					} else {
						t[i+1] = alt[i]
					}					
				}	
				//sort.Strings(t)
				//array[int(pos)] = SNP{t} // asign SNP at pos
				tmp, ok := array[int(pos)]
				if ok {
					t = append(t[:0], t[1:]...)
					tmp.profile = append(tmp.profile, t...)
				} else {
					tmp.profile = append(tmp.profile, t...)
				}
				sort.Strings(tmp.profile)
				array[int(pos)] = tmp // append SNP at pos
				//fmt.Printf("pos=%d %q \n", pos, alt)
			} else {				
				//array[int(pos)] = SNP{[]string{split[3], split[4]}} // asign SNP at pos
				tmp, ok := array[int(pos)]
				if ok {
					if split[4] == "<DEL>" {
						tmp.profile = append(tmp.profile, ".")
					} else {
						tmp.profile = append(tmp.profile, split[4])
					}					
				} else {
					if split[4] == "<DEL>" {
						tmp.profile = append(tmp.profile, []string{split[3], "."}...)
					} else {
						tmp.profile = append(tmp.profile, []string{split[3], split[4]}...)
					}
				}
				sort.Strings(tmp.profile)
				array[int(pos)]= tmp // append SNP at pos
				//fmt.Println(pos)
			}
		}
	}
    return array
}

func fastaRead(sequence_file string) []byte {
    f,err := os.Open(sequence_file)
    if err != nil{
        fmt.Printf("%v\n",err)
        os.Exit(1)
    }

    defer f.Close()
    br := bufio.NewReader(f)
    byte_array := bytes.Buffer{}

    //line , err := br.ReadString('\n')
	_ , isPrefix, err := br.ReadLine()
	if err != nil || isPrefix{
		fmt.Printf("%v\n",err)
		os.Exit(1)
	}
    //fmt.Printf("%s",line)

    for {
        line , isPrefix, err := br.ReadLine()
        if err != nil || isPrefix{
            break
        } else {
            byte_array.Write([]byte(line))
        }
    }
    //byte_array.Write([]byte("$"))
    input := []byte(byte_array.String())
    return input
}

