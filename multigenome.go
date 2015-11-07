//-------------------------------------------------------------------------------------------------
// IVC: multigenome.go
// Building multigenome by combining SNPs and INDELs from dbSNPs with the reference genome.
// Multigenome includes a multi-sequence (include standard bases and *) and a variant profile.
// Copyright 2015 Nam Sy Vo.
//-------------------------------------------------------------------------------------------------

package ivc

import (
	"bufio"
	"bytes"
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
	Seq        []byte            //multi-genome sequence.
	ChrPos     []int             //position (first base) of the chromosome on whole-genome.
	ChrName    [][]byte          //chromosome names
	Variants   map[int][][]byte  //variants (position, variants).
	VarAF      map[int][]float32 //allele frequency of variants (position, allele frequency).
	SameLenVar map[int]int       //indicate if variants has same length (SNPs or MNPs).
	DelVar     map[int]int       //length of deletions if variants are deletion.
	RevFMI     FMIndex           //FM-index of reverse multigenomes (to do forward search).
}

//--------------------------------------------------------------------------------------------------
// MultiGenome creates a multi-genome instance and set up its content.
//--------------------------------------------------------------------------------------------------
func NewMultiGenome(para_info *ParaInfo) *MultiGenome {
	log.Printf("Memstats (golang name):\tAlloc\tTotalAlloc\tSys\tHeapAlloc\tHeapSys")

	M := new(MultiGenome)
	M.ChrPos, M.ChrName, M.Seq = LoadMultiSeq(para_info.Ref_file)
	PrintMemStats("Memstats after loading multi-sequence")
	M.Variants, M.VarAF = LoadVarProf(para_info.Var_prof_file)
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
	PrintMemStats("Memstats after creating auxiliary data structures")

	M.RevFMI = *Load(para_info.Rev_index_file)
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
func BuildMultiGenome(genome_file, var_prof_file string) (chr_pos []int, chr_name [][]byte,
	seq []byte, var_prof map[string]map[int]VarProf) {

	chr_pos, chr_name, seq = GetGenome(genome_file)
	PrintProcessMem("Memstats after reading reference genome")
	var_prof = GetVarProf(var_prof_file)
	PrintProcessMem("Memstats after reading variant profile")
	var contig_name string
	var ok bool
	for i, pos := range chr_pos {
		contig_name = string(chr_name[i])
		if _, ok = var_prof[contig_name]; ok {
			for key, _ := range var_prof[contig_name] {
				seq[pos+key] = '*'
			}
		} else {
			log.Println("Missing chromosome " + contig_name + " in variant profile")
		}
	}
	return chr_pos, chr_name, seq, var_prof
}

//-------------------------------------------------------------------------------------------------
// LoadMultiSeq loads multi-sequence from file.
//-------------------------------------------------------------------------------------------------
func LoadMultiSeq(file_name string) (chr_pos []int, chr_name [][]byte, multi_seq []byte) {
	f, e := os.Open(file_name)
	if e != nil {
		log.Panicf("Error: %s", e)
	}
	defer f.Close()
	chr_pos = make([]int, 0)
	r := bufio.NewReader(f)
	var pos int
	var line, sline []byte
	multi_seq = make([]byte, 0)
	for {
		line, e = r.ReadBytes('\n')
		if e != nil { //reach EOF
			break
		}
		sline = bytes.Trim(line, "\n\r")
		if len(sline) == 0 {
			continue
		}
		if sline[0] == '>' {
			split := bytes.Split(sline, []byte("\t"))
			pos, _ = strconv.Atoi(string(split[1]))
			chr_pos = append(chr_pos, pos)
			chr_name = append(chr_name, split[0][1:])
		} else {
			multi_seq = append(multi_seq, sline...)
		}
	}
	return chr_pos, chr_name, multi_seq
}

//-------------------------------------------------------------------------------------------------
// SaveMultiSeq saves multi-sequence to file.
//-------------------------------------------------------------------------------------------------
func SaveMultiSeq(file_name string, chr_pos []int, chr_name [][]byte, multi_seq []byte) {
	f, e := os.Create(file_name)
	if e != nil {
		log.Panicf("Error: %s", e)
	}
	defer f.Close()
	w := bufio.NewWriter(f)
	for i := 0; i < len(chr_pos); i++ {
		w.WriteString(">" + string(chr_name[i]) + "\t" + strconv.Itoa(chr_pos[i]) + "\n")
	}
	w.Write(multi_seq)
	w.WriteString("\n")
	w.Flush()
}

//-------------------------------------------------------------------------------------------------
// LoadVarProf loads variant profile from file and return a map of variants.
//-------------------------------------------------------------------------------------------------
func LoadVarProf(file_name string) (variant map[int][][]byte, af map[int][]float32) {

	f, e := os.Open(file_name)
	if e != nil {
		log.Panicf("Error: %s", e)
	}
	defer f.Close()

	variant = make(map[int][][]byte)
	af = make(map[int][]float32)
	var line []byte
	var sline string
	var split, t []string
	var i int
	var k int64
	r := bufio.NewReader(f)
	for {
		line, e = r.ReadBytes('\n')
		if e != nil {
			break
		}
		sline = string(bytes.Trim(line, "\n\r"))
		if len(sline) == 0 {
			continue
		}
		split = strings.Split(sline, "\t")
		k, e = strconv.ParseInt(split[0], 10, 64)
		t = make([]string, (len(split)-1)/2)
		p := make([]float32, (len(split)-1)/2)
		for i = 0; i < len(t); i++ {
			t[i] = split[2*i+1]
			tmp, _ := strconv.ParseFloat(split[2*i+2], 32)
			p[i] = float32(tmp)
		}
		b := make([][]byte, len(t))
		for i = 0; i < len(b); i++ {
			b[i] = make([]byte, len(t[i]))
			copy(b[i], []byte(t[i]))
		}
		variant[int(k)] = b
		af[int(k)] = p
	}
	return variant, af
}

//-------------------------------------------------------------------------------------------------
// SaveVarProf saves variant profile to file.
//-------------------------------------------------------------------------------------------------
func SaveVarProf(file_name string, chr_pos []int, chr_name [][]byte, var_prof map[string]map[int]VarProf) {
	f, e := os.Create(file_name)
	if e != nil {
		log.Panicf("Error: %s", e)
	}
	defer f.Close()
	w := bufio.NewWriter(f)
	var var_pos []int
	var var_prof_chr map[int]VarProf
	for i, contig_name := range chr_name {
		var_prof_chr = var_prof[string(contig_name)]
		var_pos = make([]int, 0)
		for pos, _ := range var_prof_chr {
			var_pos = append(var_pos, pos)
		}
		sort.Sort(sort.IntSlice(var_pos))
		for _, pos := range var_pos {
			w.WriteString(strconv.Itoa(chr_pos[i]+pos) + "\t")
			for idx, val := range var_prof_chr[pos].Variant {
				w.WriteString(string(val) + "\t" + strconv.FormatFloat(float64(var_prof_chr[pos].AleFreq[idx]), 'f', 10, 32) + "\t")
			}
			w.WriteString("\n")
		}
	}
	w.Flush()
}

//--------------------------------------------------------------------------------------------------
// GetGenome gets reference genome from FASTA files.
//--------------------------------------------------------------------------------------------------
func GetGenome(file_name string) (chr_pos []int, chr_name [][]byte, seq []byte) {
	f, e := os.Open(file_name)
	if e != nil {
		log.Panicf("Error: %s", e)
	}
	defer f.Close()

	chr_pos = make([]int, 0)
	chr_name = make([][]byte, 0)
	seq = make([]byte, 0)
	var line []byte
	var sub_line [][]byte
	scanner := bufio.NewScanner(f)
	for scanner.Scan() {
		line = scanner.Bytes()
		if len(line) == 0 {
			continue
		}
		if line[0] == '>' {
			sub_line = bytes.Split(line, []byte(" "))
			chr_pos = append(chr_pos, len(seq))
			contig_name := make([]byte, len(sub_line[0][1:]))
			copy(contig_name, sub_line[0][1:])
			chr_name = append(chr_name, contig_name)
		} else {
			seq = append(seq, line...)
		}
	}
	return chr_pos, chr_name, seq
}

//--------------------------------------------------------------------------------------------------
// GetVarProf gets variant profile from VCF files.
//--------------------------------------------------------------------------------------------------
func GetVarProf(file_name string) map[string]map[int]VarProf {

	f, e := os.Open(file_name)
	if e != nil {
		log.Panicf("Error: %s", e)
	}
	defer f.Close()

	var_prof := make(map[string]map[int]VarProf)
	var line, sline, info, sub_info, tmp_af []byte
	var sub_line, sub_info_part, info_arr [][]byte
	var i, var_pos int
	var alt_prob float32
	var tmp_p float64
	var af []float32
	r := bufio.NewReader(f)
	for {
		line, e = r.ReadBytes('\n')
		if e != nil {
			break
		}
		sline = bytes.Trim(line, "\n\r")
		if sline[0] == '#' || len(sline) == 0 {
			continue
		} else {
			sub_line = bytes.SplitN(sline, []byte("\t"), 9)

			var_prof_elem := VarProf{}
			ref := make([]byte, len(sub_line[3]))
			copy(ref, sub_line[3])
			var_prof_elem.Variant = append(var_prof_elem.Variant, ref)

			alt := make([]byte, len(sub_line[4]))
			copy(alt, sub_line[4])
			alt_arr := bytes.Split(alt, []byte(","))
			info = make([]byte, len(sub_line[7]))
			copy(info, sub_line[7])
			info_arr = bytes.Split(sub_line[7], []byte(";"))
			af = make([]float32, 0)
			for _, sub_info = range info_arr {
				sub_info_part = bytes.Split(sub_info, []byte("="))
				if bytes.Equal(sub_info_part[0], []byte("AF")) {
					for _, tmp_af = range bytes.Split(sub_info_part[1], []byte(",")) {
						tmp_p, _ = strconv.ParseFloat(string(tmp_af), 32)
						af = append(af, float32(tmp_p))
					}
					break
				}
			}
			var_prof_elem.AleFreq = append(var_prof_elem.AleFreq, 0)
			if len(af) == len(alt_arr) {
				alt_prob = float32(0)
				for i := 0; i < len(alt_arr); i++ {
					var_prof_elem.Variant = append(var_prof_elem.Variant, alt_arr[i])
					var_prof_elem.AleFreq = append(var_prof_elem.AleFreq, af[i])
					alt_prob += af[i]
				}
				var_prof_elem.AleFreq[0] = 1 - alt_prob
			} else {
				alt_prob = 1 / float32(1+len(alt_arr))
				for i = 0; i < len(alt_arr); i++ {
					var_prof_elem.Variant = append(var_prof_elem.Variant, alt_arr[i])
					var_prof_elem.AleFreq = append(var_prof_elem.AleFreq, alt_prob)
				}
				var_prof_elem.AleFreq[0] = alt_prob
			}
			chr_name := string(sub_line[0])
			if _, ok := var_prof[chr_name]; !ok {
				var_prof[chr_name] = make(map[int]VarProf)
			}
			var_pos, _ = strconv.Atoi(string(sub_line[1]))
			var_prof[chr_name][var_pos-1] = var_prof_elem
		}
	}
	return var_prof
}
