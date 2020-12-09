#!/usr/bin/env python
import argparse
import os 
import sys
import pandas as pd 
import re
global args
first = True

def describe(row):
    size = "{:.0f}kbp".format( (row.qe - row.qs)/1000.0 )
    sm = row.r
    if (row.strand == "+" and row.direction == ">") or (row.strand == "-" and row.direction == "<"):
        t = "INS"
    else:
        t = "INV"
    if(sm == args.ref and t == "INS"):
        return(f"Syntenic") 
    else:
        return("{}_{}".format(t, size, sm))

    #f1 = "fc"
    #f2 = "fc"
    #if row.strand == "-" : f1 = "rc"
    #if row.direction == "<": f2 = "rc"
    #return f"the {f1} of {row.q}:{row.qs}-{row.qe} aligns to the {f2} of {row.r}:{row.rs}-{row.re}"
    




def parse_path(row):
    global first
    row = row.copy()
    paths = re.findall(r"(>|<)([^<>:]+):(\d+)-(\d+)", row.path) 
    # 100% syntenicA
    #out = {"direction":[], "q":[], "qs":[], "qe":[], "path":[],"paths":[], "pathe":[]}
    out = []
    if(len(paths) == 0 ):
        #print(f"{row.q}\t{row.qs}\t{row.qe}\t{row.path}\t{row.paths}\t{row.pathe}\tsyntenic")
        out.append((row.q, row.qs, row.qe, row.ql, row.path, row.paths, row.pathe, row.strand, ">"))
    
    #print()
    #print(paths) 
    # path through the graph
    first = True
    for idx, path in enumerate(paths):
        direction, ref, start, end = path[0], path[1], int(path[2]),int(path[3])
        if(idx==0):
            start = start + row.paths
        if(idx == len(paths)-1 ):
            end = end - (row.pathl - row.pathe)

        qs = row.qs
        qe = qs + end-start
        row.qs += end - start
        #if(row.reference == ref ):
        #    status = "syntenic"
        #else:
        #    status = "insertion"
           
        out.append((row.q, qs, qe, row.ql, ref, start, end, row.strand, direction))

    out = pd.DataFrame(out, columns=["q", "qs", "qe", "ql", "r", "rs", "re", "strand", "direction"])
    out["description"] = out.apply(describe, axis=1)
    #out.sort_values(by=["q","qs"]).to_csv(sys.stdout, sep="\t", header=first, index=False)
    out.sort_values(by=["q","qs"], inplace=True)
    first=False
    return(out)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("infile", help="positional input")
    parser.add_argument("-r", "--ref", help="Main reference sequence to compare to", default = "CHM13.pri__1")
    parser.add_argument("-n", "--number", help="numeric option", type=int, default=5)
    parser.add_argument("-l", "--list", nargs="*", help="list with zero or more entries")
    parser.add_argument("-l2", "--list2", nargs="+", help="list one or more entries")
    parser.add_argument('-d', help="store args.d as true if -d",  action="store_true", default=False)
    args = parser.parse_args()

    colnames = ["q", "ql", "qs", "qe","strand", "path", "pathl", "paths", "pathe", "matches", "block", "qual"]
    gaf = pd.read_csv(args.infile, header=None, sep="\t").loc[: ,0:(len(colnames)-1)]; gaf.columns = colnames
    gaf["reference"] = args.ref
    #print(gaf[colnames[0:7]])

    out = pd.concat( list(gaf.apply(parse_path, axis = 1)) )
    out.sort_values(by=["q","qs"]).to_csv(sys.stdout, sep="\t", header=True, index=False)
    for (q, r), group in out.groupby(by = ["q", "r"]):
        syn = group[group.description == "Syntenic"]
        if(r != args.ref or syn.shape[0] <= 1):
            continue
        shift = syn.shift(periods=1)
        deletions = ((syn.qs - shift.qe).abs() < 5000)
        for idx, d in enumerate(deletions):
            if(not d):
                continue
            cur = syn.iloc[idx -1]
            nxt = syn.iloc[idx]
            size = "{:.0f}kbp".format( (nxt.rs - cur.re)/1000 )
            #print(cur, nxt, size)
            print(f"{cur.q}\t{cur.qe}\t{cur.qe}\t{cur.ql}\t{cur.r}\t{cur.re}\t{nxt.rs}\t+\t>\tDEL_{size}")
    






