import os 
import sys
import re
import re
import pysam 
import pandas as pd
from datetime import date
from snakemake.remote.HTTP import RemoteProvider as HTTPRemoteProvider
from snakemake.remote.FTP import RemoteProvider as FTPRemoteProvider
FTP = FTPRemoteProvider()
HTTP = HTTPRemoteProvider()

today = date.today()
DATE =  today.strftime("%Y/%m/%d")

SDIR=os.path.realpath(os.path.dirname(srcdir("env.cfg"))+"/..")
shell.prefix(f"source {SDIR}/env.cfg ; set -eo pipefail; ")

# delete if not debug
DEBUG=True
def tempd(fname):
	if(DEBUG):
		return(fname)
	return(temp(fname))


FASTA = os.path.abspath( config["fasta"] )
FAI = FASTA + ".fai"
assert os.path.exists(FAI), f"Index must exist. Try: samtools faidx {FASTA}"

# WILDCARDS
NIDS = min(200, len(open(FAI).readlines()) )
IDS = [ "{:08}".format(ID+1) for ID in range(NIDS) ]

SM = "asm"
if("sample" in config): SM = config["sample"]
SPECIES = "human"
if("species" in config): SPECIES = config["species"]
THREADS = 16
if("threads" in config): THREADS = config["threads"]


SMS = [SM]

wildcard_constraints:
	SM="|".join(SMS),
	ID="\d+",

# temp output
FASTA_FMT = f"Masked/temp/{SM}_{{ID}}.fasta"

#
# OUTPUTS
#
DUP = f"Masked/{SM}_dupmasker.tbl"
COLOR = f"Masked/{SM}_dupmasker_colors.tbl"
BED = f"Masked/{SM}_dupmasker_colors.bed"
DUPLICONS = f"Masked/{SM}.duplicons"
DUPBED = f"Masked/{SM}.duplicons.bed"


DMHTML = f"Masked/{SM}_dupmasker_colors.html"
RM = os.path.abspath(f"Masked/{SM}_repeatmasker.out")
RMBED = os.path.abspath(f"Masked/{SM}_repeatmasker.out.bed")
TRFBED = os.path.abspath(f"Masked/{SM}_trf.bed")


rule split_fasta:
	input:
		fasta = FASTA,
	output:
		fastas = tempd(expand(FASTA_FMT, ID=IDS)),
	threads: 1
	resources:
		mem=8
	run:
		fasta = pysam.FastaFile(input["fasta"])
		outs = [open(f,"w+") for f in output.fastas]
		outidx = 0
		for name in fasta.references:
			seq = fasta.fetch(name)
			outs[outidx].write( ">{}\n{}\n".format(name, seq) )
			outidx += 1
			if(outidx == NIDS): outidx = 0 

		for out in outs: 
			out.close()


####################################################################
#################### REPEAT MASKER #################################
####################################################################

rule RunRepeatMasker:
	input:
		fasta = FASTA_FMT,
	output:
		out = tempd(FASTA_FMT + ".out"),
		msk = tempd(FASTA_FMT + ".masked"),
	resources:
		mem=8,
	threads: THREADS 
	shell:"""
echo "RM on {input.fasta}"
RepeatMasker \
	-s \
	-xsmall \
	-e ncbi \
	-species {SPECIES} \
	-dir $(dirname {input.fasta}) \
	-pa {threads} \
	{input.fasta} 
	
if [ -f "{output.msk}" ]; then
    echo "masked fasta exists"
else 
    echo "No repeats found, copying unmasked fasta to masked fasta"
	cp {input.fasta} {output.msk}
fi
"""

rule splitRepeatMaskerBed:
	input:
		out = rules.RunRepeatMasker.output.out,
	output:
		bed = tempd(FASTA_FMT + "_rm.bed"),
	resources:
		mem=8,
	threads: 1
	shell:"""
RM2Bed.py -d $(dirname {output.bed}) {input.out}
"""
	
#
# RepeatMasker merge output
#
rule mergeRM:
	input:
		outs = expand( rules.RunRepeatMasker.output.out, ID=IDS, SM=SM),
	output:
		out=RM,
	resources:
		mem=4,
	shell:"""
head -n 3 {input.outs[0]} > {output.out}  && \
	tail -q -n +4 {input.outs} >> {output.out}
"""

rule mergeRMbed:
	input:
		beds = expand( rules.splitRepeatMaskerBed.output.bed, ID=IDS, SM=SM),
	output:
		bed=RMBED,
		fofn=temp(f"Masked/{SM}.rm.fofn"),
	resources:
		mem=4,
	threads: 1
	run:
		open(output.fofn, "w+").write("\n".join(input.beds) + "\n")
		shell("while IFS= read -r line; do cat $line; done < {output.fofn} | bedtools sort -i - > {output.bed}")


rule RepeatMasker:
    input:
        out = rules.mergeRM.output.out,
        bed = rules.mergeRMbed.output.bed,

####################################################################
####################### DUP MASKER #################################
####################################################################

rule DupMaskerRM:
	input:
		fasta = FASTA_FMT,
		out = rules.RunRepeatMasker.output.out,
	output:
		dupout = tempd(FASTA_FMT + ".dupout"),
	resources:
		mem=8,
	threads: THREADS
	shell:"""
DupMaskerParallel \
	-pa {threads} -dupout \
	-engine ncbi \
	{input.fasta}
"""
"""
if [ -f "{output.dupout}" ]; then
    echo "exists"
else 
    echo "No repeats found, touching output"
	cp {input.fasta} {output.msk}
fi
"
"""

rule RunDupMasker:
	input:
		fasta = FASTA_FMT,
		out = rules.RunRepeatMasker.output.out,
		dupout = rules.DupMaskerRM.output.dupout,
	output:
		dups = tempd(FASTA_FMT + ".duplicons"),
	resources:
		mem=8,
	threads:1
	shell:"""
DupMaskerParallel \
	-engine ncbi \
	{input.fasta}
"""


rule DupMaskerColor:
	input:
		dups = rules.RunDupMasker.output.dups,
	output:
		dupcolor = tempd(FASTA_FMT + ".duplicons.extra"),
	resources:
		mem=4,
	shell:"""
{SDIR}/scripts/DupMask_parserV6.pl -i {input.dups} -E -o {output.dupcolor}
"""




"""
SW score = smith-waterman score of the match (complexity-adjusted )
perc div. = %substitutions in matching region.
perc del. = %deletions (in query seq rel to subject) in matching region.
perc ins. = %insertions (in query seq rel to subject) in matching region.
qry seq = id of query sequence.
qry begin = starting position of match in query sequence.
qry end = ending position of match in query sequence.
qry (left) = no. of bases in query sequence past the ending position of match (so 0 means that the match extended all the way to the end of the query sequence).
C = "C" match is found on the reverse strand
subj seq = id of the duplicon.
subj (left) = The remaining bases in (complement of) subject sequence prior to beginning of the match.
subj end = starting position of match in subject sequence (using top-strand numbering).
subj begin = ending position of match in subject sequence.
"""
rule duplicons:
    input:
        dups = expand( rules.RunDupMasker.output.dups, ID=IDS, SM=SM),
    output:
        dups=DUPLICONS,
        bed=DUPBED,
    resources:
        mem=4,
    threads: 1 
    run:
        shell(" cat {input.dups} | sort -k 5,5 -k6,6n > {output.dups}") 
        bed = open(output.bed,"w+")
        bed.write("#contig\tstart\tend\tduplicon\tscore\tstrand\tancestral_position\tchr_band\tsub_rate\tdel_rate\tins_rate\n")
        for line in open(output.dups):
            line=line.strip().split()
            if(line[8]=="C"):
                line[8]="-"
            else:
                line.insert(8, "+")
           
            score, sub, D, I, q_nm, q_st, q_en, q_left, strand, r_nm, r_left, r_st, r_en = line[0:13]
            duplicon, anc, band = r_nm.split("|")
            band = anc.split(":")[0].strip("chr") + band
            if(band=="NANA"): band="NA"

            bed.write(
                    ("{}"+"\t{}"*10 + "\n").format(
                        q_nm, q_st, q_en, duplicon,
                        score, strand, anc, band,
                        sub, D, I)
                    )
        bed.close()






#
# DupMakser merge output
#
rule mergeDM:
	input:
		dups = expand(rules.RunDupMasker.output.dups,ID=IDS, SM=SM),
		colors = expand(rules.DupMaskerColor.output.dupcolor,ID=IDS, SM=SM),
	output:
		dup = DUP,
		color = COLOR,
	resources:
		mem=4,
	threads: 1 
	run:
		f=open(output.dup, "w+")
		for dup in input.dups:
			f.write(open(dup).read())
		f.close()
	
		f=open(output.color, "w+")
		for idx, color in enumerate(input.colors):
			if(idx == 0):
				f.write(open(color).read())
			else:
				# all but the first line
				f.write( "".join( open(color).readlines()[1:] ) )
		f.close()
	


def hex_to_rgb(h):
	h = h.lstrip('#')
	return( ",".join( tuple( str(int(h[i:i+2], 16)) for i in (0, 2, 4))  )  )

rule DupMaskerBed:
  input:
    #color = rules.mergeDM.output.color,
    colors = expand(rules.DupMaskerColor.output.dupcolor,ID=IDS, SM=SM),
  output:
    bed =  BED
  resources:
    mem=4,
  threads: 1 
  run:
    colors = []
    for f in input.colors:
      #chr    chrStart	chrEnd  orient  Repeat  color   width   offset
      color = pd.read_csv(f, sep="\s+")
      colors.append(color)

    color = pd.concat(colors, ignore_index=True)
    print(color)
    
    # these rows are no good, they come from contigs that have messed up results. fix TODO
    bad= (color["chr"] == 0 ) & ( color["chrEnd"] == "#BEBEBE")
    color.drop(color[bad].index, inplace = True)

    color.sort_values(by=["chr", "chrStart"], inplace=True)

    color["strand"] = "+"
    color.loc[color["orient"] == "R", "strand"] = "-"

    color["rgb"] = color["color"].map(hex_to_rgb)
    color["score"] = 0

    color.sort_values(by=["chr", "chrStart"], inplace=True)
    out = color[ ["chr", "chrStart", "chrEnd", "Repeat", "score", "strand", "chrStart", "chrEnd", "rgb"] ]
    out.to_csv(output["bed"], sep="\t", header=False, index=False)			


rule DupMaskerHTML:
	input:
		bed=BED,
	output:
		html = DMHTML,
	resources:
		mem=4,
	threads: 1 
	run:
		html = open(f"{SDIR}/templates/dupmasker.html").read()
		open(output["html"], "w+").write(html.format(DATE=DATE, SM=SM))

rule DupMaskerSummary:
    input:
        bed=rules.duplicons.output.bed,
    output:
        tbl = f"Masked/{SM}.duplicons.tbl",
        excel = f"Masked/{SM}_duplicons.xlsx",
    resources:
        mem=4,
    threads: 1 
    run:
        shell("""sort -k4,4 {input.bed} | \
                awk '{{print $0"\t"$3-$2}}' | \
                bedtools groupby -g 4,8,7 -c 4,12 -o count,sum -i - > {output.tbl}""")
        df = pd.read_csv(output.tbl, sep="\t", names=["Duplicon", "Chr band", "Ancestral position", "Count", "bp"])
        with pd.ExcelWriter(output.excel) as writer:
            df.to_excel(writer, sheet_name='Duplicons', index=False)


rule DupMasker:
    input:
        bed=BED,
        color=COLOR,
        dup=DUP,

 
####################################################################
####################### TRF MASKER #################################
####################################################################

rule run_trf:
    input:
        fasta = FASTA_FMT,
    output:
        dat = tempd(FASTA_FMT + ".dat")
    benchmark:
        FASTA_FMT + ".bench"
    resources:
        mem=24,
    threads: 1
    shell:"""
trf {input.fasta} 2 7 7 80 10 50 15 -l 25 -h -ngs > {output.dat}
"""

rule trf_bed:
	input:
		dats = expand(rules.run_trf.output.dat, ID=IDS, SM=SM),
	output:
		bed = TRFBED,
	resources:
		mem=8,
	threads: 1
	run:
		trf = []
		header = '#chr start end PeriodSize CopyNumber ConsensusSize PercentMatches PercentIndels Score A C G T Entropy Motif Sequence'.split()
		for datf in input.dats:
			chrom = None
			sys.stderr.write( "\r" + datf )			
			with open(datf, 'r') as dat:
				for line in dat:
					splitline = line.split()
					if( line.startswith("Sequence:") ):
						chrom = int(line.split()[1].strip())
						#sys.stderr.write(chrom + "\n")
					elif( line.startswith("@") ):
						chrom = splitline[0][1:].strip() # grab everything after the @ in the first word
					else:
						# Catch index errors when line is blank
						try:
							# Check if in header sequence (all non-header lines start with an int: start pos)
							try:
								int(splitline[0])
							except ValueError:
								continue
							trf.append([chrom] + splitline[ 0: (len(header)-1) ] )
						except IndexError:
							pass
		trf = pd.DataFrame(trf, columns=header)
		print(trf.shape)
		
		trf["start"] = trf["start"].astype(int)
		trf.sort_values(by=["#chr", "start"], inplace=True)
		print("done sorting trf")

		trf.to_csv(output.bed, sep="\t", index=False)

rule trf:
    input:
        bed = rules.trf_bed.output.bed


####################################################################
####################### MASKED FASTA ###############################
####################################################################
# TODO


