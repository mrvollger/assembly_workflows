#!/usr/bin/env python
import argparse
import os 
import sys
import pandas as pd
from numba import njit
import numpy as np

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    #parser.add_argument("infile", help="positional input")
    #parser.add_argument("output", help="positional input")
    parser.add_argument("-s", "--symetric", help="make sedef output symetric, set to 0 to disable.", default = 1)
    parser.add_argument("-l", "--minlength", help=" ", type=int, default=1000)
    parser.add_argument("-i", "--minidentity", help=" ", type=float, default=0.9)
    parser.add_argument("--minindelidentity", help=" ", type=float, default=0.5)
    parser.add_argument("--sat", help="Remove dups that are this fraction of sat or more", type=float, default=0.70)
    parser.add_argument('-d', help="store args.d as true if -d",  action="store_true", default=False)
    args = parser.parse_args()

    df = pd.read_csv(sys.stdin, sep="\t", header=None, comment="#")

    flip = (df[0] > df[9]) | ( (df[0] == df[9]) & (df[1] >= df[10]) ) | ( (df[0] == df[9]) & (df[1] == df[10]) & (df[2] >= df[11]) )
    df["c1"] = df[0]
    df["s1"] = df[1]
    df["e1"] = df[2]
    df["c2"] = df[9]
    df["s2"] = df[10]
    df["e2"] = df[11]
    df.c1.loc[flip] = df[9].loc[flip]
    df.s1.loc[flip] = df[10].loc[flip]
    df.e1.loc[flip] = df[11].loc[flip]
    df.c2.loc[flip] = df[0].loc[flip]
    df.s2.loc[flip] = df[1].loc[flip]
    df.e2.loc[flip] = df[2].loc[flip]
   
    newcols = ["c1", "s1", "e1", "c2", "s2", "e2"]
    # sort by frac_match
    df.sort_values(by = newcols + [23], inplace=True)
    sys.stderr.write(f"{df.shape}\n") 
    df.drop_duplicates(subset=newcols, inplace=True, keep="last")
    df.drop(newcols, inplace=True, axis=1)
    sys.stderr.write(f"{df.shape}\n")

    df.sort_values(by = [0,1,2], inplace=True)
    df.to_csv(sys.stdout, index=False, sep="\t", header=False)

 
