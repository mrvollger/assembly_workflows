#!/usr/bin/env python
"""
Convert sedef output to a bedfile and add the symetric duplications
"""
import argparse
import pandas as pd
import BamBedUtils as bb


def get_color(frac_id):
    """
    convert segdup identity to color
    """
    if  frac_id >= .99 :
        return "255,103,0"
    if frac_id >= .98 :
        return "204,204,0"
    if frac_id >= 0.9 :
        mingray = 105
        maxgray = 220
        gray = int( mingray + (maxgray-mingray) * (  1  - (frac_id-.9)*10  ) )
        return f"{gray},{gray},{gray}"
    # return purple
    return "147,112,219"


SEDEF_HEADER = """chr1   start1  end1    chr2    start2  end2
name    score   strand1 strand2 max_len aln_len
comment indel_a indel_b alnB    matchB  mismatchB
transitionsB     transversions   fracMatch       fracMatchIndel  jck     k2K
aln_gaps        uppercaseA      uppercaseB      uppercaseMatches        aln_matches  aln_mismatches
aln_gaps.1        aln_gap_bases   cigar   filter_score    count_ovls      sat_bases
total_bases     sat_coverage""".strip().split()
DROP = ["aln_gaps.1", "comment",     "count_ovls", "total_bases", "sat_coverage"]


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("infile", help="positional input")
    parser.add_argument("output", help="positional input")
    parser.add_argument("filt", help="positional input")
    parser.add_argument("-s", "--symetric",
                        help="make sedef output symetric, set to 0 to disable.", default = 1)
    parser.add_argument("-l", "--minlength", help=" ", type=int, default=1000)
    parser.add_argument("-i", "--minidentity", help=" ", type=float, default=0.9)
    parser.add_argument("--minindelidentity", help=" ", type=float, default=0.5)
    parser.add_argument("--sat",
                        help="Remove dups that are this fraction of sat or more",
                        type=float, default=0.70)
    parser.add_argument('-d',
                        help="store args.d as true if -d",
                        action="store_true", default=False)
    args = parser.parse_args()

    # read in the segdups
    df =pd.read_csv(args.infile, sep="\t", names=SEDEF_HEADER, header=None, comment="#")

    # filter out nan. Sometimes sedef with report an alignment without a cigar
    # usually these are just alignments between the same coordiantes
    df.dropna(inplace=True)

    # filter out high sat regions
    df = df[df.sat_coverage <= args.sat]

    # remove extra columns
    df.drop(DROP, axis=1, inplace=True)

    # remove duplicate SDs and keep only the best alignment
    dup_cols = ["chr1", "start1", "end1", "chr2", "start2", "end2"]
    df.sort_values(by=dup_cols+["matchB"], inplace=True)
    df.drop_duplicates(subset=dup_cols, inplace=True, keep="last")

    # add a color for the browser
    df["color"] = df.fracMatch.map(get_color)

    # add a unique duplication identifier
    df["unique_id"] = list(range(df.shape[0]))
    df["original"] = True


    # make the symetric SDs so that there is evidence at both locations
    if args.symetric :
        df2 = df.copy()
        df2[["chr1", "start1", "end1"]] = df[["chr2", "start2", "end2"]]
        df2[["chr2", "start2", "end2"]] = df[["chr1", "start1", "end1"]]
        df2["cigar"] = df.apply(lambda x: bb.infer_query_cigar_s(x['cigar'], x['strand2']), axis=1)
        # fliping the strand is actually the wrong idea, since I flip the cigar to be query based
        #df2["strand2"] = df["strand1"]
        #df2["strand1"] = df["strand2"]
        df2["original"] = False
        df = pd.concat([df,df2], ignore_index=True)

    # make bed 9 format for the browser
    df.sort_values(by=["chr1", "start1"], inplace=True)
    bed9 = ["chr1", "start1", "end1", "name", "fakeScore", "strand1", "start1", "end1", "color"]
    df["name"] = df.chr2 + ":" + df.start2.astype(str) + "-" + df.end2.astype(str)
    df["fakeScore"] = 0
    extra = [ col for col in SEDEF_HEADER if col not in bed9 and col not in DROP]
    extra += ["unique_id", "original"]
    df = df[bed9 + extra]
    df.rename(columns={"chr1":"#chr1"}, inplace=True)

    # write output
    cond = ((df.aln_len >= args.minlength) &
            (df.fracMatch >= args.minidentity) &
            (df.fracMatchIndel >= args.minindelidentity))
    sd = df.loc[cond]
    filt = df.loc[~cond]
    sd.to_csv(args.output, index=False, sep="\t")
    filt.to_csv(args.filt, index=False, sep="\t")
