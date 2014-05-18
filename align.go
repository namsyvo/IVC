//--------------------------------------------------------------------------------------------------
// Aligning reads to multigenomes by extending exact matches based on read-multigenome edit distance.
// Copyright 2014 Nam Sy Vo
//--------------------------------------------------------------------------------------------------

package isc

import (
    "fmt"
    "math"
    "math/rand"
    "time"
    "sort"
    "github.com/namsyvo/multigenome"
)

//Global variables
var (
    search Search //to find exact matches (seeds) between reads and multigenomes
    DIST_THRES int //threshold for distances between reads and multigenomes
    ITER_NUM int //number of random iterations to find proper seeds
)

type Aligner struct {
    SEQ []byte //multigenomes
    SNP_PROFILE map[int][][]byte //hash table of SNP profile (position, values)
    SAME_LEN_SNP map[int]int //hash table to indicate if SNPs has same length
    SORTED_SNP_POS []int //sorted array of SNP positions
}

//--------------------------------------------------------------------------------------------------
// Init function sets initial values for global variables and parameters for Aligner object
//--------------------------------------------------------------------------------------------------
func (A *Aligner) Init(genome_file, snp_file, index_file, rev_index_file string, read_len int, re float32, k, a, n int) {
    A.SEQ = multigenome.LoadMulti(genome_file)
    A.SNP_PROFILE, A.SAME_LEN_SNP = multigenome.LoadSNPLocation(snp_file)
    A.SORTED_SNP_POS = make([]int, 0, len(A.SNP_PROFILE))
    for k := range A.SNP_PROFILE {
        A.SORTED_SNP_POS = append(A.SORTED_SNP_POS, k)
    }
    sort.Sort(sort.IntSlice(A.SORTED_SNP_POS))

    //Const for computing distance
    DIST_THRES = int(math.Ceil(float64(re) * float64(read_len) +
     float64(k) * math.Sqrt(float64(read_len) * float64(re) * float64((1 - re)))))
    ITER_NUM = a * (DIST_THRES + 1)
    multigenome.Init(DIST_THRES, A.SNP_PROFILE, A.SAME_LEN_SNP, read_len)
    fmt.Println("DIST_THRES: ", DIST_THRES)
    fmt.Println("ITER_NUM: ", ITER_NUM)

    search.Init(index_file, rev_index_file, n)
}

//-----------------------------------------------------------------------------------------------------
// FindExtension function returns alignment (snp report) between between reads and multi-genomes.
// The alignment is built within a given threshold of distance.
//-----------------------------------------------------------------------------------------------------
func (A Aligner) FindExtension(read []byte, s_pos, e_pos int, match_pos int) (int, []byte, []byte,
	 []byte, []byte, map[int][]byte, map[int][]byte, bool) {

    var ref_left_flank, ref_right_flank []byte
    read_left_flank := read[ : e_pos]
    read_right_flank := read[s_pos + 1 : ]
    var lcs_len int = s_pos - e_pos + 1

    var isSNP, isSameLenSNP bool
    left_ext_add_len, right_ext_add_len := 0, 0
    i := 0
    for i = match_pos - e_pos; i < match_pos; i++ {
        _, isSNP = A.SNP_PROFILE[i]
        _, isSameLenSNP = A.SAME_LEN_SNP[i]
        if isSNP && !isSameLenSNP {
            left_ext_add_len++
        }
    }
    for i = match_pos + lcs_len; i < (match_pos + lcs_len) + (len(read) - s_pos) - 1; i++ {
        _, isSNP = A.SNP_PROFILE[i]
        _, isSameLenSNP = A.SAME_LEN_SNP[i]
        if isSNP && !isSameLenSNP {
            right_ext_add_len++
        }
    }
    left_most_pos := match_pos - e_pos - left_ext_add_len
    if left_most_pos >= 0 {
        ref_left_flank = A.SEQ[left_most_pos : match_pos]
    } else {
        ref_left_flank = A.SEQ[0 : match_pos]
    }
    right_most_pos := (match_pos + lcs_len) + (len(read) - s_pos) - 1 + right_ext_add_len
    if  right_most_pos <= len(A.SEQ) {
        ref_right_flank = A.SEQ[match_pos + lcs_len : right_most_pos]
    } else {
        ref_right_flank = A.SEQ[match_pos + lcs_len : len(A.SEQ)]
    }

    left_d, left_D, left_m, left_n, left_A, left_T, _ :=
     multigenome.BackwardDistanceMulti(read_left_flank, ref_left_flank, left_most_pos)
    right_d, right_D, right_m, right_n, right_A, right_T, _ :=
     multigenome.ForwardDistanceMulti(read_right_flank, ref_right_flank, match_pos + lcs_len)

    dis := left_d + right_d + left_D + right_D
    if dis <= DIST_THRES {
        left_snp := multigenome.BackwardTraceBack(read_left_flank, ref_left_flank,
         left_m, left_n, left_A, left_T, left_most_pos)
        right_snp := multigenome.ForwardTraceBack(read_right_flank, ref_right_flank,
         right_m, right_n, right_A, right_T, match_pos + lcs_len)
        return dis, read_left_flank, read_right_flank, ref_left_flank, ref_right_flank,
		left_snp, right_snp, true
    }
    return dis, read_left_flank, read_right_flank, ref_left_flank, ref_right_flank,
	map[int][]byte{}, map[int][]byte{}, false
}

//-----------------------------------------------------------------------------------------------------
// FindSNPProfile_read returns SNP profile of new genome based on SNP profile of reference multi-genomes
// and alignment between reads and multi-genomes.
//-----------------------------------------------------------------------------------------------------
func (A *Aligner) FindSNPProfile_read(read1, read2 []byte) (map[int][][]byte, bool) {

    snp_profile := make(map[int][][]byte)
    var k int
    var v []byte

    //Consider reverse complement of reads.
    rev_read1 := reverse_complement(read1)
    rev_read2 := reverse_complement(read2)
    /*
    fmt.Println("read1", string(read1))
    fmt.Println("rev_read1", string(rev_read1))
    fmt.Println("read2", string(read2))
    fmt.Println("rev_read2", string(rev_read2))
    */
    //Find SNPs for pairend reads, treat each end separately and independently.
    var p int
    s_pos, e_pos, match_pos, hasExactMatches := -1, -1, []int{}, false
    loop_num := 1
    r := rand.New(rand.NewSource(time.Now().UnixNano()))

    //Find SNPs for the first end
    for loop_num <= ITER_NUM {
        //fmt.Println(loop_num, "read1", string(read1))
        p = r.Intn(len(read1) - 1) + 1
        //Call FindLCS to determine seed
        s_pos, e_pos, match_pos, hasExactMatches = search.FindLCS(read1, p)
        if hasExactMatches {
            for _, pos := range match_pos {
                //Call IntervalHasSNP to determine whether extension is needed
                //if search.IntervalHasSNP(A.SORTED_SNP_POS, pos - e_pos, pos - e_pos + len(read1)) {
                    //Call ApproxSearch to determine extension
                    _, _, _, _, _, left_snp, right_snp, isExtended := A.FindExtension(read1, s_pos, e_pos, pos)
                    if isExtended {
                        //Determine SNP profile
                        for k, v = range left_snp {
                            snp_profile[k] = append(snp_profile[k], v)
                        }
                        for k, v = range right_snp {
                            snp_profile[k] = append(snp_profile[k], v)
                        }
                    }
                //}
            }
            if len(snp_profile) > 0 {
                break
            }
        }
        loop_num++
    }

    if (len(snp_profile) == 0) {
        //Find SNPs for the reverse complement of first end
        for loop_num <= ITER_NUM {
            //fmt.Println(loop_num, "rev_read1", string(rev_read1))
            p = r.Intn(len(rev_read1) - 1) + 1
            //Call FindLCS to determine seed
            s_pos, e_pos, match_pos, hasExactMatches = search.FindLCS(rev_read1, p)
            if hasExactMatches {
                for _, pos := range match_pos {
                    //Call IntervalHasSNP to determine whether extension is needed
                    //if search.IntervalHasSNP(A.SORTED_SNP_POS, pos - e_pos, pos - e_pos + len(read1)) {
                        //Call ApproxSearch to determine extension
                        _, _, _, _, _, left_snp, right_snp, isExtended := A.FindExtension(rev_read1, s_pos, e_pos, pos)
                        if isExtended {
                            //Determine SNP profile
                            for k, v = range left_snp {
                                snp_profile[k] = append(snp_profile[k], v)
                            }
                            for k, v = range right_snp {
                                snp_profile[k] = append(snp_profile[k], v)
                            }
                        }
                    //}
                }
                if len(snp_profile) > 0 {
                    break
                }
            }
            loop_num++
        }
    }

    snp_prof_len_1 := len(snp_profile)

    //Find SNPs for the second end
    loop_num = 1
    snp_found_num := 0
    for loop_num <= ITER_NUM {
        //fmt.Println(loop_num, "read2", string(read2))
        p = r.Intn(len(read2) - 1) + 1
        //Call FindLCS to determine seed
        s_pos, e_pos, match_pos, hasExactMatches = search.FindLCS(read2, p)
        if hasExactMatches {
            for _, pos := range match_pos {
                //Call IntervalHasSNP to determine whether extension is needed
                //if search.IntervalHasSNP(A.SORTED_SNP_POS, pos - e_pos, pos - e_pos + len(read1)) {
                    //Call ApproxSearch to determine extension
                    _, _, _, _, _, left_snp, right_snp, isExtended := A.FindExtension(read2, s_pos, e_pos, pos)
                    if isExtended {
                        snp_found_num += 1
                        //Determine SNP profile
                        for k, v = range left_snp {
                            snp_profile[k] = append(snp_profile[k], v)
                        }
                        for k, v = range right_snp {
                            snp_profile[k] = append(snp_profile[k], v)
                        }
                    }
                //}
            }
            if snp_found_num > 0 {
                return snp_profile, true
            }
        }
        loop_num++
    }

    //Find SNPs for the reverse complement of second end
    loop_num = 1
    snp_found_num = 0
    for loop_num <= ITER_NUM {
        //fmt.Println(loop_num, "rev_read2", string(rev_read2))
        p = r.Intn(len(rev_read2) - 1) + 1
        //Call FindLCS to determine seed
        s_pos, e_pos, match_pos, hasExactMatches = search.FindLCS(rev_read2, p)
        if hasExactMatches {
            for _, pos := range match_pos {
                //Call IntervalHasSNP to determine whether extension is needed
                //if search.IntervalHasSNP(A.SORTED_SNP_POS, pos - e_pos, pos - e_pos + len(read1)) {
                    //Call ApproxSearch to determine extension
                    _, _, _, _, _, left_snp, right_snp, isExtended := A.FindExtension(rev_read2, s_pos, e_pos, pos)
                    if isExtended {
                        snp_found_num += 1
                        //Determine SNP profile
                        for k, v = range left_snp {
                            snp_profile[k] = append(snp_profile[k], v)
                        }
                        for k, v = range right_snp {
                            snp_profile[k] = append(snp_profile[k], v)
                        }
                    }
                //}
            }
            if snp_found_num > 0 {
                return snp_profile, true
            }
        }
        loop_num++
    }

    if snp_prof_len_1 > 0 {
        return snp_profile, true
    } else {
        return snp_profile, false
    }
}

func reverse_complement(read []byte) []byte {
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