#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Author: Mitchell R. Vollger
import os
import sys
import argparse
import pandas as pd
from numba import njit
import numpy as np
import BamBedUtils as bb

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("infile", help="positional input")
    parser.add_argument("-d", help="store args.d as true if -d",  action="store_true", default=False)
    args = parser.parse_args()

    f="/net/eichler/vol26/projects/chm13_t2t/nobackups/FES_stuff/acro_cn/probe.locations.bed
    pd.read_csv("")



"""
    df = pd.read_csv(args.infile, sep="\t")
    df.dropna(inplace=True)
    df["passed"]=True
    for idx, row in df.iterrows():
        cigar = bb.parse_cigar(row.cigar)
        qstart, qend = bb.liftover(cigar, 0, row.end1-row.start1, row.strand2)
        #print(row, qstart, qend)
        if qstart + row.start2 == row.start2 and qend + row.start2 == row.end2:
            continue
        else:
            print(f"{row}\n{qstart},{0},{qend},{row.end2-row.start2}")
            df["passed"][idx]=False
    
    df[~df.passed].to_csv("tmp.failed.bed", sep="\t", index=False)
"""
