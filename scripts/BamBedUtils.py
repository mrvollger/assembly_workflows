#!/usr/bin/env python
"""
A module for bed and bam operations
"""
from numba import njit
import numpy as np
import sys
M=0 #M  BAM_CMATCH      0
I=1 #I  BAM_CINS        1
D=2 #D  BAM_CDEL        2
N=3 #N  BAM_CREF_SKIP   3
S=4 #S  BAM_CSOFT_CLIP  4
H=5 #H  BAM_CHARD_CLIP  5
P=6 #P  BAM_CPAD        6
E=7 #=  BAM_CEQUAL      7
X=8 #X  BAM_CDIFF       8
B=9 #B  BAM_CBACK       9
NM=10 #NM       NM tag  10
CON_T = [M, D, N, E, X] # these ones "consume" the reference
CON_Q = [M, I, S, E, X] # these ones "consume" the query
CON_A = [M, I, D, N, S, E, X] # these ones "consume" the alignments

INT_TO_CHAR_CIGAR = {0:"M", 1:"I", 2:"D", 3:"N", 4:"S",
                     5:"H", 6:"P", 7:"=", 8:"X", 9:"B"}
CHAR_TO_INT_CIGAR = {"M":0, "I":1, "D":2, "N":3, "S":4,
                     "H":5, "P":6, "=":7, "X":8, "B":9}


@njit
def intersect(a, a0, a1, b, b0, b1, dist = 0, limit=2):
    """
    check if two genomic intervals overlap
    """
    if limit :
        dist = min(dist, limit*(a1-a0 + b1-b0))
    return (a == b) and ((a1 + dist) >= b0) and ((b1 + dist) >= a0)

@njit
def sd_or_intersect(qrow, trow1, trow2):
    """
    does qrow overlap either side of the SD
    """
    return intersect(qrow, trow1) or intersect(qrow, trow2)

def df_or_intersect(qrow, df):
    """
    check if qrow intersect with the first three columns or the next three
    """
    return df.apply(lambda x: sd_or_intersect(qrow, df.iloc[:, 0:3], df.iloc[:, 3:6]), axis=1)

@njit
def intersect_list(a_chrs, a_starts, a_ends, b_chr, b_start, b_end):
    rtn = np.zeros(a_chrs.shape[-1], dtype=bool)
    i = 0
    for ct, st, en in zip(a_chrs, a_starts, a_ends):
        rtn[i] = intersect(ct, st, en, b_chr, b_start, b_end)
        i += 1
    return rtn

def intersect_dfs(df_a, df_b):
    """
    check for elements of df_a that intersect with any of df_w, output like -wa bedtools
    """
    a = np.array(df_a.iloc[:, 0:3])
    b = np.array(df_b.iloc[:, 0:3])




@njit
def _str_to_int(s):
    """
    from https://github.com/numba/numba/issues/5650
    """
    final_index, result = len(s) - 1, 0
    for i,v in enumerate(s):
        result += (ord(v) - 48) * (10 ** (final_index - i))
    return result

@njit
def parse_cigar(cigar_s):
    """
    function to turn a cigar string into a list of tuples
    M=0 #M  BAM_CMATCH      0
    I=1 #I  BAM_CINS        1
    D=2 #D  BAM_CDEL        2
    N=3 #N  BAM_CREF_SKIP   3
    S=4 #S  BAM_CSOFT_CLIP  4
    H=5 #H  BAM_CHARD_CLIP  5
    P=6 #P  BAM_CPAD        6
    E=7 #=  BAM_CEQUAL      7
    X=8 #X  BAM_CDIFF       8
    B=9 #B  BAM_CBACK       9
    NM=10 #NM       NM tag  10
    """
    digits = [str(i) for i in range(10)]
    i=0
    cigar = []
    while i < len(cigar_s):
        digit = ""
        while cigar_s[i] in digits:
            digit += cigar_s[i]
            i+=1
        length = _str_to_int(digit)
        c =  cigar_s[i]
        if c=="M":
            opt = M
        elif c=="I":
            opt = I
        elif c=="D":
            opt = D
        elif c=="N":
            opt = N
        elif c=="S":
            opt = S
        elif c=="H":
            opt = H
        elif c=="P":
            opt = P
        elif c=="E":
            opt = E
        elif c=="X":
            opt = X
        elif c=="B":
            opt = B
        else:
            raise RuntimeError("Cigar opt not known.")
        cigar.append((length,opt))
        #print(length, opt, c)
        i+=1
    return np.array(cigar)

@njit
def infer_query_cigar(cigar, strand):
    """
    function to invert a cigar string
    """
    q_cigar = []
    for length, opt in cigar:
        if opt not in [M, E, X, I, D]:
            raise RuntimeError("Cannot handle opts outside of M,E,X,I,D (TODO).\n")
        if opt == I:
            opt = D
        elif opt == D:
            opt = I
        elif opt in (N, S):
            continue
        q_cigar.append((length, opt))
    if strand == "-":
        q_cigar = q_cigar[::-1]
    return q_cigar

def cigar_tup_to_s(cigar):
    """
    convert cigar tuple to string
    """
    rtn = []
    for length, opt in cigar:
        rtn.append(str(length) + INT_TO_CHAR_CIGAR[opt] )
    return "".join(rtn)

def infer_query_cigar_s(cigar_s, strand):
    """
    convert cigar string to the equivelent query cigar string
    """
    #print(cigar_s)
    q_cigar = infer_query_cigar(parse_cigar(cigar_s), strand)
    return cigar_tup_to_s(q_cigar)

@njit
def _lift(cigar, t_start, t_end):
    """
    helper function for liftover
    returns the qstart and qend (no strand taken into affect)
    """
    INT_TO_CHAR_CIGAR = {0:"M", 1:"I", 2:"D", 3:"N", 4:"S",
                         5:"H", 6:"P", 7:"=", 8:"X", 9:"B"}
    conRef = [M, D, N, E, X] # these ones "consume" the reference
    conQuery = [M, I, S, E, X] # these ones "consume" the query
    conAln = [M, I, D, N, S, E, X] # these ones "consume" the alignments
    first = True
    q_start = -1
    q_end = -1
    t_pos = 0
    q_pos = 0
    for length, opt in cigar:
        i = 0
        #sys.stdout.write(str(length) + ":" + INT_TO_CHAR_CIGAR[opt] + ", ")
        #print(length,INT_TO_CHAR_CIGAR[opt])
        while i < length:
            # increment counters
            if opt in conRef :
                t_pos += 1
            if opt in conQuery :
                q_pos += 1

            # set the query start position
            if q_start == -1 and t_pos >= t_start :
                q_start = max(q_pos - 1, 0)
            # set the end cord if we are there
            if q_end == -1 and t_pos >= t_end :
                q_end = q_pos
                return (q_start, q_end)
            #print(t_pos, t_end, q_pos, q_end, opt)
            i+=1
    return (q_start, q_end)

@njit
def infer_q_length(cigar):
    """
    infer the query length from the cigar
    """
    con_q = [M, I, S, E, X] # these ones "consume" the query
    rtn = 0
    for length, opt in cigar:
        if opt in con_q:
            rtn += length
    return rtn

@njit
def liftover(cigar, t_start, t_end, strand, q_length=None):
    """
    lift reference alignment to query sequence
    """
    q_start, q_end = _lift(cigar, t_start, t_end)
    
    if q_start == -1 or q_end == -1:
        return (q_start, q_end)
    
    if strand == "-":
        if q_length is None:
            q_length = infer_q_length(cigar)
        return (q_length - q_end, q_length - q_start)
    return (q_start, q_end)


