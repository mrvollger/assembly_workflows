#!/usr/bin/env python
import argparse
import os 
import sys
import pandas as pd
import numpy as np
import numba as nb

# global var for inputs
args=None 

@nb.jit(nopython=True)
def dseg(s, S, D, chrm, minlen):
    N = len(s)
    cumul=0; mmax = 0; start = 0;
    for i in range(N):
        cumul = cumul + s[i]
        if (cumul >= mmax):
            mmax = cumul; end = i
        if( cumul <= 0 or cumul <= mmax + D or i == N-1 ):
            if (mmax >= S):
                #print("{}\t{}\t{}\t".format(start, end, mmax) )
                if(end-start > minlen):
                    print(chrm,start, end, mmax)
            mmax = 0 
            cumul = 0 
            start = i 
            end = i 


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("bed", help="positional input")
    parser.add_argument("fai", help="positional input")
    parser.add_argument("-S", "--score", help="numeric option", type=int, default=1000000)
    parser.add_argument("-l", "--minlen", help="numeric option", type=int, default=0)
    parser.add_argument("-D", "--drop", help="numeric option", type=int, default=-500000)
    parser.add_argument("--scale", help="positve regons get + scale per base", type=int, default=3)
    parser.add_argument('-d', help="store args.d as true if -d",  action="store_true", default=False)
    args = parser.parse_args()
    
    assert args.score >= -args.drop and args.drop < 0

    fai=pd.read_csv(args.fai, sep="\s+", names=["chr", "chrlen", "x","y","z"])
    fai.set_index("chr", inplace=True)
    bed=pd.read_csv(args.bed, sep="\s+")
    bed.rename(columns = {bed.columns[0]:"chr", bed.columns[1]:"start",bed.columns[2]:"end"}, inplace=True)
   
    for chrm, region in bed.groupby("chr"):
        s = np.zeros( fai.chrlen[chrm], dtype=np.int8); s[:]=-1
        sys.stderr.write(chrm+"\n")
        for row in region.itertuples():
            s[row.start:row.end]=args.scale
            continue
        #print(chrm, s.sum()/1000000)
        dseg(s, args.score, args.drop, chrm, args.minlen)


