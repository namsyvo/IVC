//---------------------------------------------------------------------------------------------------
// Sharing variables and functions for modules in ISC package.
// Copyright 2014 Nam Sy Vo
//---------------------------------------------------------------------------------------------------

package isc

import (
	"bufio"
	"fmt"
	"os"
	"strings"
	"strconv"
	"bytes"
	"sort"
    "github.com/vtphan/fmi"
)

type Index struct {
    SEQ []byte //multigenomes
    SNP_PROFILE map[int][][]byte //hash table of SNP Profile (position, snps)
    SNP_AF map[int][]float32 //allele frequency of SNP Profile (position, af of snps)
    SAME_LEN_SNP map[int]int //hash table to indicate if SNPs has same length
    SORTED_SNP_POS []int //sorted array of SNP positions
    REV_FMI fmi.Index //FM-index of reverse multigenomes
}

//--------------------------------------------------------------------------------------------------
// IntervalHasSNP determines whether [i, j] contains SNP positions which are stores in array A.
// This function impelements interpolation search. The array A must be sorted in increasing order.
//--------------------------------------------------------------------------------------------------
func IntervalHasSNP(A []int, i, j int) bool {
    L := 0
    R := len(A) - 1
    var m int
	for A[L] <= i && i <= A[R] && A[L] != A[R] {
		m = L + (R - L) * ((i - A[L]) / (A[R] - A[L]))  //out of range is possible here		
		if (A[m] < i) {
			L = m + 1;
		} else if A[m] > i {
			R = m - 1
		} else {
			return i <= j
		}
	}
	return i <= j && L < len(A) && i <= A[L] && j >= A[L]
}

//--------------------------------------------------------------------------------------------------
// ReverseComplement returns reverse complement of a read.
//--------------------------------------------------------------------------------------------------

func ReverseComplement(read []byte) []byte {
    l := len(read)
    rev_read := make([]byte, l)
    for idx, elem := range read {
        if elem == 'A' {
            rev_read[l-1-idx] = 'T'
        } else if (elem == 'T') {
            rev_read[l-1-idx] = 'A'
        } else if (elem == 'C') {
            rev_read[l-1-idx] = 'G'
        } else if (elem == 'G') {
            rev_read[l-1-idx] = 'C'
        } else {
            rev_read[l-1-idx] = elem
        }
    }
    return rev_read
}

func ReadVCF(sequence_file string) map[int]SNP {
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
					tmp.Profile = append(tmp.Profile, t...)
				} else {
					tmp.Profile = append(tmp.Profile, t...)
				}
				sort.Strings(tmp.Profile)
				array[int(pos)] = tmp // append SNP at pos
				//fmt.Printf("pos=%d %q \n", pos, alt)
			} else {				
				//array[int(pos)] = SNP{[]string{split[3], split[4]}} // asign SNP at pos
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

				array[int(pos)]= tmp // append SNP at pos
				//fmt.Println(pos)
			}
		}
	}
    return array
}

func ReadFASTA(sequence_file string) []byte {
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

