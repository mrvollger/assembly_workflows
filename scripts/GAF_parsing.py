#!/usr/bin/env python
import argparse
import os 
import sys
import pandas as pd 
import re

first = True

def describe(row):
    f1 = "fc"
    f2 = "fc"
    if row.strand == "-" : f1 = "rc"
    if row.direction == "<": f2 = "rc"
    return f"the {f1} of {row.q}:{row.qs}-{row.qe} aligns to the {f2} of {row.r}:{row.rs}-{row.re}"




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
    for path in paths:
        direction, ref, start, end = path[0], path[1], int(path[2]),int(path[3])

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
    out.sort_values(by=["q","qs"]).to_csv(sys.stdout, sep="\t", header=first, index=False)
    first=False

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

    gaf.apply(parse_path, axis = 1)