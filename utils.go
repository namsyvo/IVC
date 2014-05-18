//---------------------------------------------------------------------------------------------------
// Utility functions for ISC package.
// Copyright 2014 Nam Sy Vo
//---------------------------------------------------------------------------------------------------

package isc

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