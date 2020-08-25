#!/usr/bin/env python
import argparse
import os 
import sys
import pandas as pd

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("liftoff", help="positional input, sedef outputs")
    parser.add_argument("-x", "--excel", help="place to write excel", default=None)
    parser.add_argument("-n", "--number", help="numeric option", type=int, default=5)
    parser.add_argument('-d', help="store args.d as true if -d",  action="store_true", default=False)
    args = parser.parse_args()

    df = pd.read_csv(args.liftoff, sep="\t", low_memory=False)
    df["extra_copy_num"] = df.copy_num_ID.str.extract(".*(_\d+)")[0].str.lstrip("_").astype(int)
    df["sequence_ID"] = df.sequence_ID.str.replace("_\d+","").astype(float)
    df.insert(6,"origin gene",df.gene_name.str.replace(r"_\d+",""))
    df.sort_values(["origin gene","geneid","extra_copy_num","cdslen"], ascending=[True,True,True,False],  inplace=True)

    df["ORF"] = False
    for geneid, group in df.groupby("geneid"):
        if(group.cdslen.sum()>0):
            df.loc[group.index, "ORF"]=True


    df.drop(["copy_num_ID"], inplace=True, axis=1)
    new = df.extra_copy_num > 0
    longest = ~ df.duplicated(subset=["geneid"], keep="first")

    #
    # write output
    #
    if(args.excel):
        with pd.ExcelWriter(args.excel) as writer:  
            df[new & longest & df.ORF].to_excel(writer, sheet_name='NewCopiesLongestCDS', index=False)
            df[new & longest].to_excel(writer, sheet_name='NewCopiesLongestCDSAll', index=False)
            df[df.ORF].to_excel(writer, sheet_name='TranscriptsORF', index=False)
            df.to_excel(writer, sheet_name='TranscriptsAll', index=False)



