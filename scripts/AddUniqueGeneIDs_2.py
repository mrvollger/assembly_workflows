#!/usr/bin/env python
import argparse
import os 
import sys
import re
import pandas 

# global var for inputs
args=None 

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Adds unique IDs for duplicated entries from liftoff so the output will work with UCSC tools", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("gff", help="gff3 from liftoff with --extra-copies ")
    parser.add_argument('-d', help="store args.d as true if -d",  action="store_true", default=False)
    args = parser.parse_args()
    

    dupatts = "(gene_name|ID|Parent|transcript_id|gene_id|protein_id|transcript_name)"
    attpat = re.compile(f"({dupatts}=[^;]+);")
    #sys.stderr.write(attpat)
    copy_num_ID = ""
    for line in open(args.gff):
        t = line.strip().split()
        if(line[0]=="#"): 
            sys.stdout.write(line)
            continue
        
        typ = t[2]
        if(typ == "gene"):
            match = re.findall("copy_num_ID=([^;]+)(_\d+)", t[8] )
            #sys.stderr.write(str(match)+"\n")
            if(len(match)==0):
                continue
            assert len(match) == 1
            geneid, copy_num_ID = match[0]
            if(copy_num_ID == "_0"):
                copy_num_ID = ""
            else:
                copy_num_ID = "_"+str(int(copy_num_ID.strip("_")) + 1)

        # replace all instances of attributes with attributes that have the dup number
        line, count = re.subn( attpat, r"\1"+copy_num_ID+";", line)
        sys.stdout.write(line)
