#!/usr/bin/env python
import argparse
import os 
import sys
from pybedtools import BedTool
import pandas as pd

# global var for inputs
args=None 
ACRO = ["chr"+str(x) for x in [13,14,15,21,22] ]


def parm(rec):
    return(rec.start < CEN[rec.chrom][0])

def qarm(rec):
    return(rec.end > CEN[rec.chrom][1])

def acrofilt(rec):
    cond = (rec[chr1] in ACRO) and (rec[chr2] in ACRO) and parm(rec)
    return(cond)

def pairfilt(rec, chrm, other):
    cond = (rec[chr1] in chrm) and (rec[chr2] in other)
    return(cond)

def intra(rec):
    return(rec[chr1] == rec[chr2])

def inter(rec):
    return(rec[chr1] != rec[chr2])

def chrfilt(rec, CHR):
    return(rec[0]==CHR)

def bedcount(bed):
    bp = 0
    count = 0
    for rec in bed:
        bp += int(rec[2]) - int(rec[1])
        count += 1
    return(bp, count)

def count_up(bed, prefix="", header=False):
    if(header):
        print("#label\tnon-redundant\t\t\t\t\t\tredundant")
        print("#label\ttotal\t\tintra\t\tinter\t\ttotal\t\tintra\t\tinter")
        print("#label"+"\tbp\tcount"*6)

    bed=bed.saveas()
    out = [] 
    for bp, count in [
            bedcount(bed.merge()),
            bedcount(bed.filter(intra).merge()),
            bedcount(bed.filter(inter).merge()),
            bedcount(bed),
            bedcount(bed.filter(intra)),
            bedcount(bed.filter(inter))
            ]:
        #out += "\t{:,}\t{:,}".format(bp, count)
        out.append(bp)
        out.append(count)
    return(out)

def make_cen(rmf):
    rm = BedTool(rmf)
    cen = rm.filter(
            lambda x: x[3]=="ALR/Alpha" and x.end - x.start > 50 ).merge(d=100).filter(
                    lambda x: x.end - x.start > 10000
                    ).merge(
                            d=200000
                    )
    global CEN
    CEN={}
    for rec in cen:
        if(rec.chrom not in CEN):
            CEN[rec.chrom]=(rec.start, rec.end)
        elif(CEN[rec.chrom][1] - CEN[rec.chrom][0] < rec.end - rec.start):
            CEN[rec.chrom]=(rec.start, rec.end)



def unique(sequence):
    seen = set()
    return [x for x in sequence if not (x in seen or seen.add(x))]

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("sedef", help="positional input")
    parser.add_argument("-x", "--excel", help="place to write excel", default=None)
    parser.add_argument("-r", "--rm", help="RepeatMasker Bed", default=None)
    parser.add_argument("-n", "--number", help="numeric option", type=int, default=5)
    parser.add_argument("-l", "--list", nargs="*", help="list with zero or more entries")
    parser.add_argument("-l2", "--list2", nargs="+", help="list one or more entries")
    parser.add_argument('-d', help="store args.d as true if -d",  action="store_true", default=False)
    args = parser.parse_args()
    
    # make a cen track so I can get p and q arms
    if(args.rm):
       make_cen(args.rm)

    # readin file and make headers global varables so we can acess them
    header = open(args.sedef).readline().strip("#").split()
    for idx, col in enumerate(header):
        exec("global " + col)
        exec(col + " = " + str(idx))
    sedef = BedTool(args.sedef)
    CHRS=unique([rec[0] for rec in sedef])
    CHRS= ["chr"+str(x) for x in range(1,23) ] + ["chrX", "chrY"]

    #
    # column names for the table
    #
    columns = pd.MultiIndex.from_product([['non-redundant', 'redundant'], ['total', 'inter', 'intra'], ['bp','total']])
    
    #
    # 
    #
    index = pd.MultiIndex.from_product( [[""], ["all"]+CHRS])
    out = [count_up(sedef, prefix="all")]
    for chrm in CHRS:
        out.append(count_up(sedef.filter(chrfilt, chrm), prefix=chrm))
    chrtbl=pd.DataFrame(out, columns=columns, index=index)
    chrtbl.label="hello"
    print(chrtbl)

    #
    # Acrocentric table 
    #
    acro = sedef.filter(acrofilt)
    acro = acro.saveas()
    index2 = pd.MultiIndex.from_product( [ACRO, ACRO])
    a = []#[count_up(acro, prefix="Acrocentric", header=True)]
    for idx, chrm in enumerate(ACRO):
        for other in ACRO:
            a.append(count_up(acro.filter(pairfilt, chrm, other)) )
    acrotbl=pd.DataFrame(a, columns=columns, index=index2)
    print(acrotbl)

    #
    # paired table 
    #
    index3 = pd.MultiIndex.from_product( [CHRS, CHRS])
    a = []
    for idx, chrm in enumerate(CHRS):
        for other in CHRS:
            a.append(count_up(sedef.filter(pairfilt, chrm, other)) )
    pairtbl=pd.DataFrame(a, columns=columns, index=index3)
    print(pairtbl)



    #
    # write output
    #
    if(args.excel):
        with pd.ExcelWriter(args.excel) as writer:  
            chrtbl.to_excel(writer, sheet_name='Chromosomes')
            acrotbl.to_excel(writer, sheet_name='Acrocentric')
            pairtbl.to_excel(writer, sheet_name='PairedChromosomes')



