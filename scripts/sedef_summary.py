#!/usr/bin/env python
import argparse
import os 
import sys
from pybedtools import BedTool
import pandas as pd
import warnings
import numba as nb
warnings.simplefilter(action='ignore', category=pd.errors.PerformanceWarning)

# global var for inputs
args=None 
ACRO = ["chr"+str(x) for x in [13,14,15,21,22] ]
subsets = {}


@nb.jit
def sd_info(rec):
    st, en = rec['start1'], rec['end1']
    st_o, en_o = rec['start2'], rec['end2']
    

    cst, cen = CEN[rec["#chr1"]]
    mn = cst - 5e6; mx = cen + 5e6
    cst_o, cen_o = CEN[rec['chr2']]
    mn_o = cst_o - 5e6; mx_o = cen_o + 5e6
    
    p = st < cst
    q = en > cen
    pi = en >= mn and st <= mx 

    p_o = st_o < cst_o
    q_o = en_o > cen_o
    pi_o = en_o >= mn_o and st_o <= mx_o



    telo = (st < 0 + 500000) or (en > FAI.loc[rec["#chr1"]]["length"] - 500000)
    telo_o = (st_o < 0 + 500000) or (en_o > FAI.loc[rec["chr2"]]["length"] - 500000)
    
    acro = (rec["#chr1"] in ACRO)
    acro_o =  (rec["chr2"] in ACRO)
    if(acro or acro_o):
        acro = True


    inter = rec["#chr1"] != rec["chr2"]
    intra = rec["#chr1"] == rec["chr2"]
    return(p, pi, q, p_o, pi_o, q_o, telo, telo_o, acro, acro, inter, inter, intra, intra)


@nb.jit
def parm(rec):
    return(rec.start < CEN[rec.chrom][0])

@nb.jit
def qarm(rec):
    return(rec.end > CEN[rec.chrom][1])

@nb.jit
def peri(rec):
    mn = CEN[rec.chrom][0] - 5e6
    mx = CEN[rec.chrom][1] + 5e6
    return( rec.end >= mn and rec.start <= mx )

@nb.jit
def acrofilt(rec):
    cond = (rec[chr1] in ACRO) and (rec[chr2] in ACRO) and parm(rec)
    return(cond)

@nb.jit
def pairfilt(rec, chrm, other):
    cond = (rec[chr1] == chrm) and (rec[chr2] == other)
    return(cond)

@nb.jit
def intra(rec):
    return(rec[chr1] == rec[chr2])

@nb.jit
def inter(rec):
    return(rec[chr1] != rec[chr2])

@nb.jit
def chrfilt(rec, CHR):
    return(rec[0]==CHR)

@nb.jit
def chrfilt2(rec, CHR):
    return(rec[chr2]==CHR)

@nb.jit
def bedcount(bed):
    bp = 0
    count = 0
    for rec in bed:
        bp += int(rec[2]) - int(rec[1])
        count += 1
    return(bp, count)

@nb.jit
def bases(bed):
    return(bedcount(bed)[0])

@nb.jit
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

def read_cen(cenf):
    global CEN
    cen = BedTool(cenf)
    CEN={}
    for rec in cen:
        if(rec.chrom not in CEN):
            CEN[rec.chrom]=(rec.start, rec.end)
        elif(CEN[rec.chrom][1] - CEN[rec.chrom][0] < rec.end - rec.start):
            CEN[rec.chrom]=(rec.start, rec.end)

def sd_stats(bed):
    nalns = bed.count()  
    bp = bases(bed.merge())
    rbp = bases(bed)
    return(nalns, bp, rbp)

def unique(sequence):
    seen = set()
    return [x for x in sequence if not (x in seen or seen.add(x))]

def add_summ_stats(summ, chr1, chr2, df):
    pair=df
    for rgn1 in RGNS:
        for rgn2 in RGNS:
            if(rgn1 == "All" and rgn2 == "All"):
                rgn = df
            elif(rgn1 == "All"):
                rgn = pair[pair[rgn2+"2"]]
            elif(rgn2 == "All"):
                rgn = pair[pair[rgn1]]
            else:
                rgn = pair[pair[rgn1] & pair[rgn2+"2"]]
            rgnbed = BedTool.from_dataframe(rgn)
            summ.loc[(chr1, rgn1), (chr2, rgn2)] = sd_stats(rgnbed)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("sedef", help="positional input, sedef outputs")
    parser.add_argument("-x", "--excel", help="place to write excel")
    parser.add_argument("-r", "--rm", help="RepeatMasker Bed", default=None)
    parser.add_argument("--cen", help="File with cen locations, Bed")
    parser.add_argument("--fai", help="Fai of the assembly", default=None)
    parser.add_argument("-n", "--number", help="numeric option", type=int, default=5)
    parser.add_argument('-d', help="store args.d as true if -d",  action="store_true", default=False)
    args = parser.parse_args()

    global FAI
    FAI = pd.read_csv(args.fai, sep="\s+", index_col=0, names=["chr","length","x","y","z"])
    #.loc["chr1"]["length"]

    # make a cen track so I can get p and q arms
    if(args.cen):
       read_cen(args.cen)

    # readin file and make headers global varables so we can acess them
    header = open(args.sedef).readline().strip("#").split()
    for idx, col in enumerate(header):
        exec("global " + col)
        exec(col + " = " + str(idx))
    CHRS= ["All"] + ["chr"+str(x) for x in range(1,23) ] + ["chrX", "chrY"]
    RGNS = ["All","Interchromosomal", "Intrachromosomal", "pericentromeric", "telomeric", "acrocentric"]
    STATS = ['#SDs',"bp", "AlignedBases"]

    df = pd.read_csv(args.sedef, sep="\t")
    sys.stderr.write(f"\r{df.shape}")
    df = df[ df["#chr1"].isin(CHRS) & df["chr2"].isin(CHRS) ]
    sys.stderr.write(f"\r{df.shape}")
    # add info to the df
    df["p-arm"],df["pericentromeric"],df["q-arm"], df["p-arm2"],df["pericentromeric2"],df["q-arm2"],df["telomeric"],df["telomeric2"], df["acrocentric"], df["acrocentric2"], df["Interchromosomal"], df["Interchromosomal2"], df["Intrachromosomal"], df["Intrachromosomal2"] = zip(*df.apply(sd_info, axis=1))
    
    #df["info"] = df.apply(sd_info, axis=1)
    
    index = pd.MultiIndex.from_product([CHRS, RGNS, STATS])
    columns = pd.MultiIndex.from_product([CHRS, RGNS])
    summ = pd.DataFrame(index=index, columns=columns)
   
    # add stats for all chromosomes vs a chromosome
    add_summ_stats(summ, "All", "All", df)
    for chrm in CHRS[1:]:
        sys.stderr.write(f"\r{chrm}")
        add_summ_stats(summ, "All", chrm, df[df["chr2"] == chrm])
        add_summ_stats(summ, chrm, "All", df[df["#chr1"]==chrm])

    # add stats for pairs of chromosomes
    for (chr1, chr2), pair in df.groupby(by=["#chr1", "chr2"]):
        sys.stderr.write(f"\r{chr1}\t{chr2}")
        add_summ_stats(summ, chr1, chr2, pair)

    with pd.ExcelWriter(args.excel) as writer:  
        #first = pd.DataFrame(summ.loc[(slice(None), slice(None)), ("All","All")]).reset_index(level=[0,1,2])
        #first = first.pivot_table(index="level_0", columns=["level_1","level_2"])
        first = pd.DataFrame(summ.loc[(slice(None), slice(None)), ("All","All")]).reset_index(level=[0,1,2]).pivot_table(index=["level_0", "level_2"], columns=["level_1"])
        first.columns = first.columns.droplevel(level=[0,1])
        first.to_excel(writer, sheet_name="SD summary",  float_format="{:,}")

        for rgn in RGNS:
            for stat in STATS:
                name = (rgn + " " + stat).lstrip("All ")
                tmp = summ.loc[(slice(None),rgn,stat), (slice(None),"All")]
                tmp.columns = tmp.columns.droplevel(1)
                tmp.index = tmp.index.droplevel([1,2])
                tmp.to_excel(writer, sheet_name=name)
        summ.to_excel(writer, sheet_name='Matrix')
    exit(0)






    sedef = BedTool(args.sedef)
    info = []
    for rec in sedef:
        info.append(sd_info(rec))

    pp   
    #CHRS=unique([rec[0] for rec in sedef])

    # calculate all subsets 
    for cr in CHRS:
        sys.stderr.write(f"\r{cr}")
        crbed = sedef.filter(chrfilt, cr).saveas()
        p = crbed.filter(parm).saveas()
        q = crbed.filter(qarm).saveas()
        pericen = crbed.filter(peri).saveas()

        subsets[cr]={"all":crbed,
                "p-arm":p,
                "pericen":pericen,
                "q-arm":q}
    sys.stderr.write("\n")

    index = pd.MultiIndex.from_product([CHRS, RGNS, ['#SDs',"bp", "AlignedBases"]])
    #columns = pd.MultiIndex.from_product([CHRS, RGNS])
    df = pd.DataFrame(index=index, columns=CHRS)

    for cr1 in subsets:
        for cr2 in CHRS:
            sys.stderr.write(f"\r{cr1}\t{cr2}")
            for rgn in subsets[cr1]:
                data = subsets[cr1][rgn]
                nalns = data.filter(chrfilt2, cr2).count()  
                bp = bases(data.filter(chrfilt2, cr2).merge())
                rbp = bases(data.filter(chrfilt2, cr2))

                df.loc[(cr1, rgn, "bp"), cr2] = bp
                df.loc[(cr1, rgn, "AlignedBases"), cr2] = rbp
                df.loc[(cr1, rgn, "#SDs"), cr2] = nalns

    with pd.ExcelWriter(args.excel) as writer:  
        df.to_excel(writer, sheet_name='Matrix')
    exit()
    


    #
    # column names for the table
    #
    columns = pd.MultiIndex.from_product([['non-redundant', 'redundant'], ['total', 'intra', 'inter'], ['bp','count']])
    
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
    index2 = pd.MultiIndex.froiim_product( [ACRO, ACRO])
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
            chrtbl.to_excel(writer, sheet_name='Chromosomes', float_format="{:,}")
            acrotbl.to_excel(writer, sheet_name='Acrocentric', float_format="{:,}")
            pairtbl.to_excel(writer, sheet_name='PairedChromosomes', float_format="{:,}")



