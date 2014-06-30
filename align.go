//--------------------------------------------------------------------------------------------------
// Aligning reads to multigenomes by extending exact matches based on read-multigenome edit distance.
// Finding exact matches betwwen reads and multigenomes using exact search with FM index.
// Exact search is perfomed with regard to a random position on reads.
// Finding inexact matches betwwen reads and multigenomes by extending FM-index based exact matches
// using edit distance between reads and multigenomes.
// Determining whether an interval on multigenomes contains a SNP position using interpolation search.
// Copyright 2014 Nam Sy Vo
//--------------------------------------------------------------------------------------------------


package isc

import (
    "fmt"
    "math"
    "sort"
    "github.com/vtphan/fmi" //to use FM index
)

//Global variables
var (
    DIST_THRES int = math.MaxInt16 //threshold for distances between reads and multigenomes
    ITER_NUM int //number of random iterations to find proper seeds
    MAXIMUM_MATCH int //maximum number of matches
)

//--------------------------------------------------------------------------------------------------
// Init function sets initial values for global variables and parameters for Index object
//--------------------------------------------------------------------------------------------------
func (I *Index) Init(genome_file, snp_file, index_file, rev_index_file string, read_len int, re float32, k, a, n int) {
    I.SEQ = LoadMultigenome(genome_file)
    I.SNP_PROFILE, I.SNP_AF, I.SAME_LEN_SNP = LoadSNPLocation(snp_file)
    I.SORTED_SNP_POS = make([]int, 0, len(I.SNP_PROFILE))
    for k := range I.SNP_PROFILE {
        I.SORTED_SNP_POS = append(I.SORTED_SNP_POS, k)
    }
    sort.Sort(sort.IntSlice(I.SORTED_SNP_POS))

    //I.FMI = *fmi.Load(index_file)
    I.REV_FMI = *fmi.Load(rev_index_file)

    //Const for computing distance
    DIST_THRES = int(math.Ceil(float64(re) * float64(read_len) +
     float64(k) * math.Sqrt(float64(read_len) * float64(re) * float64((1 - re)))))
    ITER_NUM = a * (DIST_THRES + 1)
    MAXIMUM_MATCH = n
	//testing///
	//ITER_NUM = 1
	////////////
    fmt.Println("DIST_THRES: ", DIST_THRES)
    fmt.Println("ITER_NUM: ", ITER_NUM)
    fmt.Println("MAXIMUM_MATCH: ", MAXIMUM_MATCH)
    fmt.Println("INF: ", INF)
}

//--------------------------------------------------------------------------------------------------
// Bachward Search with FM-index, start from any position on the pattern.
//--------------------------------------------------------------------------------------------------

func (I Index) BackwardSearchFrom(index fmi.Index, pattern []byte, start_pos int) []int {
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
// FindSeeds function returns positions and distances of LCS between reads and multi-genomes.
// It uses both backward search and forward search (backward search on reverse references).
//--------------------------------------------------------------------------------------------------
func (I Index) FindSeeds(read []byte, p int) (int, int, []int, bool) {
    read_len := len(read)
    rev_read := make([]byte, read_len)
    for i := 0; i < read_len; i++ {
        rev_read[i] = read[read_len - i - 1]
    }
    var rev_sp, rev_ep int = 0, MAXIMUM_MATCH
    var rev_s_pos, rev_e_pos, s_pos, e_pos int
    var rev_result []int
	
    rev_s_pos = read_len - 1 - p
    rev_result = I.BackwardSearchFrom(I.REV_FMI, rev_read, rev_s_pos)
    rev_sp, rev_ep, rev_e_pos = rev_result[0], rev_result[1], rev_result[2]

    //convert rev_e_pos in forward search to s_pos in backward search
    s_pos = rev_e_pos
    e_pos = p
    if rev_ep - rev_sp + 1 <= MAXIMUM_MATCH {
        match_pos := make([]int, 0, MAXIMUM_MATCH)
        for p := rev_sp; p <= rev_ep; p++ {
            match_pos = append(match_pos, len(I.SEQ) - 1 - I.REV_FMI.SA[p])
        }
        return s_pos, e_pos, match_pos, true
    }
    /*
    result = I.BackwardSearchFrom(I.FMI, read, s_pos)
    sp, ep, e_pos = result[0], result[1], result[2]

    if ep - sp + 1 <= MAXIMUM_MATCH {
        match_pos := make([]int, 0, MAXIMUM_MATCH)
        for p := sp; p <= ep; p++ {
            match_pos = append(match_pos, I.FMI.SA[p])
        }
        return s_pos, e_pos, match_pos, true
    }
    */
    return -1, -1, []int{}, false
}

//-----------------------------------------------------------------------------------------------------
// FindExtension function returns alignment (snp report) between between reads and multi-genomes.
// The alignment is built within a given threshold of distance.
//-----------------------------------------------------------------------------------------------------
func (I Index) FindExtensions(read []byte, s_pos, e_pos int, match_pos int) (int, []byte, []byte,
	 []byte, []byte, map[int][]byte, map[int][]byte, bool) {

    var ref_left_flank, ref_right_flank []byte
    read_left_flank := read[ : e_pos]
    read_right_flank := read[s_pos + 1 : ]
    var lcs_len int = s_pos - e_pos + 1

    var isSNP, isSameLenSNP bool
    left_ext_add_len, right_ext_add_len := 0, 0
    i := 0
    for i = match_pos - e_pos; i < match_pos; i++ {
        _, isSNP = I.SNP_PROFILE[i]
        _, isSameLenSNP = I.SAME_LEN_SNP[i]
        if isSNP && !isSameLenSNP {
            left_ext_add_len++
        }
    }
    for i = match_pos + lcs_len; i < (match_pos + lcs_len) + (len(read) - s_pos) - 1; i++ {
        _, isSNP = I.SNP_PROFILE[i]
        _, isSameLenSNP = I.SAME_LEN_SNP[i]
        if isSNP && !isSameLenSNP {
            right_ext_add_len++
        }
    }
    left_most_pos := match_pos - e_pos - left_ext_add_len
    if left_most_pos >= 0 {
        ref_left_flank = I.SEQ[left_most_pos : match_pos]
    } else {
        ref_left_flank = I.SEQ[0 : match_pos]
    }
    right_most_pos := (match_pos + lcs_len) + (len(read) - s_pos) - 1 + right_ext_add_len
    if  right_most_pos <= len(I.SEQ) {
        ref_right_flank = I.SEQ[match_pos + lcs_len : right_most_pos]
    } else {
        ref_right_flank = I.SEQ[match_pos + lcs_len : len(I.SEQ)]
    }

    left_d, left_D, left_m, left_n, left_A, left_T, _ :=
     I.BackwardDistanceMulti(read_left_flank, ref_left_flank, left_most_pos)
    right_d, right_D, right_m, right_n, right_A, right_T, _ :=
     I.ForwardDistanceMulti(read_right_flank, ref_right_flank, match_pos + lcs_len)

    dis := left_d + right_d + left_D + right_D
    if dis <= DIST_THRES {
        left_snp := I.BackwardTraceBack(read_left_flank, ref_left_flank,
         left_m, left_n, left_A, left_T, left_most_pos)
        right_snp := I.ForwardTraceBack(read_right_flank, ref_right_flank,
         right_m, right_n, right_A, right_T, match_pos + lcs_len)
        return dis, read_left_flank, read_right_flank, ref_left_flank, ref_right_flank,
		left_snp, right_snp, true
    }
    return dis, read_left_flank, read_right_flank, ref_left_flank, ref_right_flank,
	map[int][]byte{}, map[int][]byte{}, false
}
