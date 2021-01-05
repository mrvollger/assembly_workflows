#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Author: Mitchell R. Vollger
"""
Module for adding region annotations to the sedef results
"""
import argparse
import BamBedUtils as bb
import pandas as pd
import sys


def is_peri(sedef, cens, slop):
    peri = cens.copy()
    peri.start = peri.start - slop
    peri.end = peri.end + slop
    return bb.df_df_intersect(peri, sedef)

def is_telo(sedef, fai, slop):
    telo1 = pd.DataFrame({"chr":fai.chr, "start":0, "end":slop})
    telo2 = pd.DataFrame({"chr":fai.chr, "start":fai.length-slop, "end":fai.length})
    telo = telo1.append(telo2)
    return bb.df_df_intersect(telo, sedef)


def is_acro(sedef, cens):
    acros = ["chr13", "chr14", "chr15", "chr21", "chr22"]
    acro = pd.DataFrame({"chr":cens.chr, "start":0, "end":cens.end})
    acro = acro[acro.chr.isin(acros)]  
    return bb.df_df_intersect(acro, sedef)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("bed", help="sedef bed")
    parser.add_argument("--cens", help="centromeric locations")
    parser.add_argument("--fai", help="lengths of the genome")
    parser.add_argument("--peri", 
                        help="size to slop to edges of cen for peri",
                        type=int, default=5000000)
    parser.add_argument("--telo", 
                        help="size to slop from ends for telomeric",
                        type=int, default=500000)
    parser.add_argument("-d",
                        help="store args.d as true if -d",
                        action="store_true", default=False)
    ARGS = parser.parse_args()

    sedef = pd.read_csv(ARGS.bed, sep="\t")
    cens = pd.read_csv(ARGS.cens,
                       sep="\t",
                       names=["chr", "start", "end"],
                       header=None)
    fai = pd.read_csv(ARGS.fai,
                      sep="\t",
                      names=["chr", "length", "x", "y", "z"],
                      header=None)
    
    sedef["telo"] = is_telo(sedef, fai, ARGS.telo)
    sedef["peri"] = is_peri(sedef, cens, ARGS.peri)
    sedef["acro"] = is_acro(sedef, cens)
    sedef["telo2"] = is_telo(sedef[["chr2","start2","end2"]], fai, ARGS.telo)
    sedef["peri2"] = is_peri(sedef[["chr2","start2","end2"]], cens, ARGS.peri)
    sedef["acro2"] = is_acro(sedef[["chr2","start2","end2"]], cens)
    
    sedef.to_csv(sys.stdout, sep="\t", index=False)





