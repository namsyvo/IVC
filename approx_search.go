//----------------------------------------------------------------------------------------
// Copyright 2013 Nam S. Vo
// Approximate searching of reads on multigenomes based on FM index and a modified-version
// of edit distance for read-multigenome alignment
//----------------------------------------------------------------------------------------

package randalx

import (
    "fmt"
    "math"
    "math/rand"
    "time"
    "sort"

    //user-defined package
    "distance2"
    "multigenome"
    "fmi"
)

// Search object with Parameters
type Search struct {
    SEQ []byte //multigenomes
    SNP_PROFILE map[int][][]byte //hash table of SNP profile (position, values)
    SAME_LEN_SNP map[int]int //hash table to indicate if SNPs has same length
    SORTED_SNP_POS []int //sorted array of SNP positions

    FMI fmi.Index //FM-index of multigenomes
    REV_FMI fmi.Index //FM-index of multigenomes
}

//Global variables
var (
    DIST_THRES int //threshold for distances between reads and multigenomes
    ITER_NUM int //number of random iterations to find proper seeds
    MIN_LCS int // minimum length of seeds
    MAXIMUM_MATCH int //maximum of number of matches
)

//-----------------------------------------------------------------------------------------------------
// Init function sets initial values for global variables and parameters for Search object
//-----------------------------------------------------------------------------------------------------
func (S *Search) Init_seq(genome_file, snp_file, index_file, rev_index_file string) {
    S.SEQ = multigenome2.LoadMulti(genome_file)
    S.SNP_PROFILE, S.SAME_LEN_SNP = multigenome2.LoadSNPLocation(snp_file)
    S.SORTED_SNP_POS = make([]int, 0, len(S.SNP_PROFILE))
    for k := range S.SNP_PROFILE {
        S.SORTED_SNP_POS = append(S.SORTED_SNP_POS, k)
    }
    sort.Sort(sort.IntSlice(S.SORTED_SNP_POS))
    S.FMI = *fmi.Load(index_file)
    S.REV_FMI = *fmi.Load(rev_index_file)

    /*
    fmt.Println("SEQ: ", string(S.SEQ))
    for _, value := range S.SNP_PROFILE {
        for _, snp := range value {
            fmt.Print(string(snp), "\t")
        }
        fmt.Println()
    }
    fmt.Println("SAME_LEN_SNP: ", S.SAME_LEN_SNP)
    fmt.Println("SORTED_SNP_POS: ", S.SORTED_SNP_POS)
    */
    
}

func (S *Search) Init_para(re float32, k, A int, read_len int) {
    //Const for computing distance
    DIST_THRES = int(math.Ceil(float64(re) * float64(read_len) +
     float64(k) * math.Sqrt(float64(read_len) * float64(re) * float64((1 - re)))))
    ITER_NUM = A * (DIST_THRES + 1)
    MIN_LCS = int(math.Ceil(float64(read_len) / float64(DIST_THRES + 1)))
    MAXIMUM_MATCH = 32
    distance.Init(DIST_THRES, S.SNP_PROFILE, S.SAME_LEN_SNP)
    
    fmt.Println("DIST_THRES: ", DIST_THRES)
    fmt.Println("ITER_NUM: ", ITER_NUM)
    fmt.Println("MIN_LCS: ", MIN_LCS)
    fmt.Println("MAXIMUM_MATCH: ", MAXIMUM_MATCH)
    
}

//-----------------------------------------------------------------------------------------------------
// FindLCS function returns positions and distances of LCS between reads and multi-genomes.
// It will use both BackwardSearch and ForwardSearch on multi-genomes
//-----------------------------------------------------------------------------------------------------
func (S Search) FindLCS(read []byte) (int, int, []int, bool) {
    r := rand.New(rand.NewSource(time.Now().UnixNano()))
    read_len := len(read)
    rev_read := make([]byte, read_len)
    for i := 0; i < read_len; i++ {
        rev_read[i] = read[read_len - i - 1]
    }
    var sp, ep int = 0, MAXIMUM_MATCH
    var rev_s_pos, rev_e_pos, s_pos, e_pos int
    var rev_result, result []int
    var count int = 0

    //re-find exact matches if number of match (= ep - sp + 1) > MAXIMUM_MATCH
    for (ep - sp + 1 > MAXIMUM_MATCH) || (s_pos - e_pos + 1 < MIN_LCS) {
        rev_s_pos = r.Intn(read_len - 1) + 1
        rev_result = S.REV_FMI.BackwardSearch(rev_read, rev_s_pos)
        _, _, rev_e_pos = rev_result[0], rev_result[1], rev_result[2]
        //convert rev_e_pos in forward search to s_pos in backward search
        s_pos = read_len - rev_e_pos - 1
        result = S.FMI.BackwardSearch(read, s_pos)
        sp, ep, e_pos = result[0], result[1], result[2]
        count++
        if count > ITER_NUM {
            return -1, -1, []int{}, false
        }
    }
    match_pos := make([]int, 0, MAXIMUM_MATCH)
    for p := sp; p <= ep; p++ {
        match_pos = append(match_pos, S.FMI.SA[p])
    }
    //fmt.Println(e_pos, s_pos, match_pos, string(read[e_pos:s_pos+1]))
    return s_pos, e_pos, match_pos, true
}

//-----------------------------------------------------------------------------------------------------
// IntervalHasSNP function determins whether [i, j] contains SNP positions which are stores in array A.
// Implement interpolation search. A must be sorted in increasing order.
//-----------------------------------------------------------------------------------------------------
func (S Search) IntervalHasSNP(A []int, i, j int) bool {
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

//-----------------------------------------------------------------------------------------------------
// FindExtension function returns alignment (snp report) between between reads and multi-genomes.
// The alignment is built within a given threshold of distance.
//-----------------------------------------------------------------------------------------------------
func (S Search) FindExtension(read []byte, s_pos, e_pos int, match_pos int) (int, []byte, []byte, []byte, []byte, map[int][]byte, map[int][]byte, bool) {
    var ref_left_flank, ref_right_flank []byte
    read_left_flank := read[ : e_pos]
    read_right_flank := read[s_pos + 1 : ]
    var lcs_len int = s_pos - e_pos + 1

    var isSNP, isSameLenSNP bool
    left_ext_add_len, right_ext_add_len := 0, 0
    i := 0
    for i = match_pos - e_pos; i < match_pos; i++ {
        _, isSNP = S.SNP_PROFILE[i]
        _, isSameLenSNP = S.SAME_LEN_SNP[i]
        if isSNP && !isSameLenSNP {
            left_ext_add_len++
        }
    }
    for i = match_pos + lcs_len; i < (match_pos + lcs_len) + (len(read) - s_pos) - 1; i++ {
        _, isSNP = S.SNP_PROFILE[i]
        _, isSameLenSNP = S.SAME_LEN_SNP[i]
        if isSNP && !isSameLenSNP {
            right_ext_add_len++
        }
    }
    left_most_pos := match_pos - e_pos - left_ext_add_len
    if left_most_pos >= 0 {
        ref_left_flank = S.SEQ[left_most_pos : match_pos]
    } else {
        ref_left_flank = S.SEQ[0 : match_pos]
    }
    right_most_pos := (match_pos + lcs_len) + (len(read) - s_pos) - 1 + right_ext_add_len
    if  right_most_pos <= len(S.SEQ) {
        ref_right_flank = S.SEQ[match_pos + lcs_len : right_most_pos]
    } else {
        ref_right_flank = S.SEQ[match_pos + lcs_len : len(S.SEQ)]
    }

    left_d, left_D, left_m, left_n, left_S, left_T :=
     distance.BackwardDistanceMulti(read_left_flank, ref_left_flank, left_most_pos)

    right_d, right_D, right_m, right_n, right_S, right_T :=
     distance.ForwardDistanceMulti(read_right_flank, ref_right_flank, match_pos + lcs_len)

    if left_d + right_d + left_D + right_D <= DIST_THRES {
        left_snp := distance.BackwardTraceBack(read_left_flank, ref_left_flank,
         left_m, left_n, left_S, left_T, left_most_pos)
        right_snp := distance.ForwardTraceBack(read_right_flank, ref_right_flank,
         right_m, right_n, right_S, right_T, match_pos + lcs_len)
        return left_d + right_d + left_D + right_D, read_left_flank, read_right_flank, ref_left_flank, ref_right_flank, left_snp, right_snp, true
    }
    return left_d + right_d + left_D + right_D, read_left_flank, read_right_flank, ref_left_flank, ref_right_flank,
     map[int][]byte{}, map[int][]byte{}, false
}

//-----------------------------------------------------------------------------------------------------
// FindSNPProfile_read returns SNP profile of new genome based on SNP profile of reference multi-genomes
// and alignment between reads and multi-genomes.
//-----------------------------------------------------------------------------------------------------
func (S *Search) FindSNPProfile_read(read []byte) (map[int][][]byte, bool) {
    snp_profile := make(map[int][][]byte)
    //Call FindLCS to determine seed
    s_pos, e_pos, match_pos, hasExactMatches := S.FindLCS(read)
    if !hasExactMatches {
        return snp_profile, false
    }
    var k int
    var v []byte
    for _, pos := range match_pos {
        //Call IntervalHasSNP to determine whether extension is needed
        if S.IntervalHasSNP(S.SORTED_SNP_POS, pos - e_pos, pos - e_pos + len(read)) {
            //Call ApproxSearch to determine extension
            _, _, _, _, _, left_snp, right_snp, isExtended := S.FindExtension(read, s_pos, e_pos, pos)
            if isExtended {
                //Determine SNP profile
                for k, v = range left_snp {
                    snp_profile[k] = append(snp_profile[k], v)
                }
                for k, v = range right_snp {
                    snp_profile[k] = append(snp_profile[k], v)
                }
            }
        }
    }
    if (len(snp_profile)) > 0 {
        return snp_profile, true
    }
    /*
    fmt.Println("Start------------------------------")
    fmt.Println("read, e_pos, s_pos, match_pos, seed ", string(read), e_pos, s_pos, match_pos, string(read[e_pos: s_pos+1]))
    for _, pos := range match_pos {
        dis, read_left_flank, read_right_flank, ref_left_flank, ref_right_flank, _, _, _ := S.FindExtension(read, s_pos, e_pos, pos)
        fmt.Println("dis ", dis)
        fmt.Println("left read ", string(read_left_flank))
        fmt.Println("left ref  ", string(ref_left_flank))
        fmt.Println("right read ", string(read_right_flank))
        fmt.Println("right ref  ", string(ref_right_flank))
    }
    fmt.Println("End------------------------------")
    */
    return snp_profile, false
}