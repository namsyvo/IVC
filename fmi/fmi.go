//----------------------------------------------------------------------------------------
// IVC: fmi.go
// Constructing FM-index.
// Copyright 2013 Vinhthuy Phan.
// Modified 2014 Nam Sy Vo.
//----------------------------------------------------------------------------------------

package fmi

import (
	"bufio"
	"encoding/binary"
	"fmt"
	"log"
	"os"
	"path"
	"sort"
	"sync"
)

//-----------------------------------------------------------------------------
// Global variables: sequence (SEQ), suffix array (SA), BWT, FM index (C, OCC)
//-----------------------------------------------------------------------------

var SEQ []byte

type Index struct {
	SA  []uint32          // suffix array
	OCC map[byte][]uint32 // occurence table
	C   map[byte]uint32   // count table
	EP  map[byte]uint32   // ending row/position of each symbol

	LEN     uint32
	END_POS uint32          // position of "$" in the text
	SYMBOLS []int           // sorted symbols
	Freq    map[byte]uint32 // Frequency of each symbol
}

//-----------------------------------------------------------------------------

func check_for_error(e error) {
	if e != nil {
		panic(e)
	}
}

//-----------------------------------------------------------------------------
// Build FM index given the file storing the text.

func New(seq []byte) *Index {
	I := new(Index)
	GetSeq(seq)
	log.Println("Building suffix array...")
	I.build_suffix_array()
	log.Println("Finish building suffix array.")
	log.Println("Building bwt and fm-index...")
	I.build_bwt_fmindex()
	log.Println("Finish building bwt and fm-index.")
	return I
}

//-----------------------------------------------------------------------------

type Symb_OCC struct {
	Symb int
	OCC  []uint32
}

//-----------------------------------------------------------------------------
// Load FM index. Usage:  idx := Load(index_file)
func Load(dirname string) *Index {

	I := new(Index)

	_load_slice := func(filename string, length uint32) []uint32 {
		f, err := os.Open(filename)
		check_for_error(err)
		defer f.Close()

		seq_len := int(I.LEN / 100)
		_, idx_fn := path.Split(filename)
		v := make([]uint32, length)
		scanner := bufio.NewScanner(f)
		scanner.Split(bufio.ScanBytes)
		for i := 0; scanner.Scan(); i++ {
			// convert 4 consecutive bytes to a uint32 number
			v[i] = uint32(scanner.Bytes()[0])
			scanner.Scan()
			v[i] += uint32(scanner.Bytes()[0]) << 8
			scanner.Scan()
			v[i] += uint32(scanner.Bytes()[0]) << 16
			scanner.Scan()
			v[i] += uint32(scanner.Bytes()[0]) << 24
			if (i+1)%(10*seq_len) == 0 {
				log.Println("Finish loading", (i+1)/seq_len, "% of index file", idx_fn)
			}
		}
		return v
	}

	// First, load "others"
	f, err := os.Open(path.Join(dirname, "others"))
	check_for_error(err)
	defer f.Close()

	var symb byte
	var freq, c, ep uint32
	scanner := bufio.NewScanner(f)
	scanner.Scan()
	fmt.Sscanf(scanner.Text(), "%d%d\n", &I.LEN, &I.END_POS)

	I.Freq = make(map[byte]uint32)
	I.C = make(map[byte]uint32)
	I.EP = make(map[byte]uint32)
	for scanner.Scan() {
		fmt.Sscanf(scanner.Text(), "%c%d%d%d", &symb, &freq, &c, &ep)
		I.SYMBOLS = append(I.SYMBOLS, int(symb))
		I.Freq[symb], I.C[symb], I.EP[symb] = freq, c, ep
	}

	// Second, load Suffix array and OCC
	I.OCC = make(map[byte][]uint32)
	var wg sync.WaitGroup
	wg.Add(5)
	go func() {
		defer wg.Done()
		I.SA = _load_slice(path.Join(dirname, "sa"), I.LEN)
	}()
	Symb_OCC_chan := make(chan Symb_OCC)
	for _, symb := range I.SYMBOLS[0:4] {
		go func(symb int) {
			defer wg.Done()
			Symb_OCC_chan <- Symb_OCC{symb, _load_slice(path.Join(dirname, "occ."+string(symb)), I.LEN)}
		}(symb)
	}
	go func() {
		wg.Wait()
		close(Symb_OCC_chan)
	}()

	for symb_occ := range Symb_OCC_chan {
		I.OCC[byte(symb_occ.Symb)] = symb_occ.OCC
	}
	return I
}

//-----------------------------------------------------------------------------
func (I *Index) Save(dirname string) {

	_save_slice := func(s []uint32, filename string) {
		f, err := os.Create(filename)
		check_for_error(err)
		defer f.Close()
		w := bufio.NewWriter(f)
		binary.Write(w, binary.LittleEndian, s)
		w.Flush()
	}

	dir := dirname + ".index"
	os.Mkdir(dir, 0777)

	var wg sync.WaitGroup
	wg.Add(5)

	go func() {
		defer wg.Done()
		_save_slice(I.SA, path.Join(dir, "sa"))
	}()

	for symb := range I.OCC {
		go func(symb byte) {
			defer wg.Done()
			_save_slice(I.OCC[symb], path.Join(dir, "occ."+string(symb)))
		}(symb)
	}

	f, err := os.Create(path.Join(dir, "others"))
	check_for_error(err)
	defer f.Close()
	w := bufio.NewWriter(f)
	fmt.Fprintf(w, "%d %d\n", I.LEN, I.END_POS)
	for i := 0; i < len(I.SYMBOLS); i++ {
		symb := byte(I.SYMBOLS[i])
		fmt.Fprintf(w, "%s %d %d %d\n", string(symb), I.Freq[symb], I.C[symb], I.EP[symb])
	}
	w.Flush()

	wg.Wait()
}

//-----------------------------------------------------------------------------
// BWT is saved into a separate file
func (I *Index) build_suffix_array() {
	I.LEN = uint32(len(SEQ))
	I.SA = make([]uint32, I.LEN)
	SA := make([]int, I.LEN)
	ws := &WorkSpace{}
	ws.ComputeSuffixArray(SEQ, SA)
	for i := range SA {
		I.SA[i] = uint32(SA[i])
	}
}

//-----------------------------------------------------------------------------
func (I *Index) build_bwt_fmindex() {
	I.Freq = make(map[byte]uint32)
	seq_len := I.LEN
	bwt := make([]byte, seq_len)
	var i uint32
	for i = 0; i < seq_len; i++ {
		I.Freq[SEQ[i]]++
		if I.SA[i] == 0 {
			bwt[i] = SEQ[seq_len-1]
		} else {
			bwt[i] = SEQ[I.SA[i]-1]
		}
		if bwt[i] == '$' {
			I.END_POS = i
		}
	}

	I.C = make(map[byte]uint32)
	I.OCC = make(map[byte][]uint32)
	for c := range I.Freq {
		I.SYMBOLS = append(I.SYMBOLS, int(c))
		I.OCC[c] = make([]uint32, seq_len)
		I.C[c] = 0
	}
	sort.Ints(I.SYMBOLS)
	I.EP = make(map[byte]uint32)
	for j := 1; j < len(I.SYMBOLS); j++ {
		curr_c, prev_c := byte(I.SYMBOLS[j]), byte(I.SYMBOLS[j-1])
		I.C[curr_c] = I.C[prev_c] + I.Freq[prev_c]
		I.EP[curr_c] = I.C[curr_c] + I.Freq[curr_c] - 1
	}

	for j := 0; j < len(bwt); j++ {
		I.OCC[bwt[j]][j] = 1
		if j > 0 {
			for symbol := range I.OCC {
				I.OCC[symbol][j] += I.OCC[symbol][j-1]
			}
		}
	}
	I.SYMBOLS = I.SYMBOLS[1:] // Remove $, which is the first symbol
	delete(I.OCC, '$')
	delete(I.C, '$')
	delete(I.OCC, 'X')
	delete(I.C, 'X')
	delete(I.OCC, 'Y')
	delete(I.C, 'Y')
	delete(I.OCC, 'Z')
	delete(I.C, 'Z')
}

//-----------------------------------------------------------------------------
func GetSeq(seq []byte) {
	SEQ = make([]byte, len(seq))
	copy(SEQ, seq)
	SEQ = append(SEQ, byte('$'))
	// replace N with X, '*' with Y, and other characters with Z (last character is '$')
	for i := 0; i < len(SEQ)-1; i++ {
		if SEQ[i] == 'N' {
			SEQ[i] = 'X'
		} else if SEQ[i] == '*' {
			SEQ[i] = 'Y'
		} else if SEQ[i] != 'A' && SEQ[i] != 'C' && SEQ[i] != 'G' && SEQ[i] != 'T' {
			log.Println("Sequence contains a non-standard base", string(SEQ[i]), "at location", i, "(will be replaced by Z)")
			SEQ[i] = 'Z'
		}
	}
}
