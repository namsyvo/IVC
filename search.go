//--------------------------------------------------------------------------------------------------
// Finding matches betwwen reads and multigenomes using exact search with FM index. Exact search
// is perfomed with regard to a random position on reads.
// Determining whether an interval on multigenomes contains a SNP position using interpolation search.
// Copyright 2014 Nam Sy Vo
//--------------------------------------------------------------------------------------------------

package isc

import (
	"github.com/vtphan/fmi" //to use FM index
)

//Global variables
var (
    MAXIMUM_MATCH int //maximum number of matches
)

// Search object
type Search struct {
    FMI fmi.Index //FM-index of multigenomes
    REV_FMI fmi.Index //FM-index of multigenomes
}

func (S *Search) Init(index_file, rev_index_file string, n int) {
    S.FMI = *fmi.Load(index_file)
    S.REV_FMI = *fmi.Load(rev_index_file)
    MAXIMUM_MATCH = n
}

//--------------------------------------------------------------------------------------------------
// Bachward Search with FM-index, start from any position on the pattern.
//--------------------------------------------------------------------------------------------------

func (S Search) BackwardSearch(index fmi.Index, pattern []byte, start_pos int) []int {
	var sp, ep, offset int
	var ok bool
	
	c := pattern[start_pos]
	sp, ok = index.C[c]
	if ! ok {
		return make([]int, 0)
	}
	ep = index.EP[c]
	var sp0, ep0 int
	// if Debug { fmt.Println("pattern: ", string(pattern), "\n\t", string(c), sp, ep) }
	for i:= start_pos - 1 ; i >= 0 ; i-- {
		//fmt.Println("pos, # candidates: ", i, ep - sp + 1)
  		c = pattern[i]
  		offset, ok = index.C[c]
  		if ok {
			sp0 = offset + index.OCC[c][sp - 1]
			ep0 = offset + index.OCC[c][ep] - 1
			if sp0 <= ep0 {
				sp = sp0
				ep = ep0
			} else {
				return []int{sp, ep, i + 1}
			}
		} else {
			return make([]int, 0)
		}
  		// if Debug { fmt.Println("\t", string(c), sp, ep) }
	}
 	return []int{sp, ep, 0}
}

//--------------------------------------------------------------------------------------------------
// FindLCS function returns positions and distances of LCS between reads and multi-genomes.
// It uses both backward search and forward search (backward search on reverse references).
//--------------------------------------------------------------------------------------------------
func (S Search) FindLCS(read []byte, p int) (int, int, []int, bool) {
    read_len := len(read)
    rev_read := make([]byte, read_len)
    for i := 0; i < read_len; i++ {
        rev_read[i] = read[read_len - i - 1]
    }
    var sp, ep int = 0, MAXIMUM_MATCH
    var rev_s_pos, rev_e_pos, s_pos, e_pos int
    var rev_result, result []int
	
    rev_s_pos = p
    rev_result = S.BackwardSearch(S.REV_FMI, rev_read, rev_s_pos)
    _, _, rev_e_pos = rev_result[0], rev_result[1], rev_result[2]
	
    //convert rev_e_pos in forward search to s_pos in backward search
    s_pos = read_len - rev_e_pos - 1
    result = S.BackwardSearch(S.FMI, read, s_pos)
    sp, ep, e_pos = result[0], result[1], result[2]
	
    if ep - sp + 1 <= MAXIMUM_MATCH {
        match_pos := make([]int, 0, MAXIMUM_MATCH)
        for p := sp; p <= ep; p++ {
            match_pos = append(match_pos, S.FMI.SA[p])
        }
        return s_pos, e_pos, match_pos, true
    }
    return -1, -1, []int{}, false
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
