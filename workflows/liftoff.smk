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



#
# hg38 information
#
GRCH38=HTTP.remote("https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.39_GRCh38.p13/GRCh38_major_release_seqs_for_alignment_pipelines/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz")
RGFF=FTP.remote("ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_34/gencode.v34.annotation.gff3.gz")

#
# Sample config
#
FASTA = config["fasta"] 
assert os.path.isabs(FASTA), f"Must specify absolute path for {FASTA}"
FAI = FASTA + ".fai"
assert os.path.exists(FAI), f"Index must exist. Try: samtools faidx {FASTA}"

SM = "asm"
if("sample" in config): SM = config["sample"]
THREADS=16
if("threads" in config): THREADS = config["threads"]
SMS = [SM]

if("ref" in config):
    REF = config["ref"]
    assert os.path.isabs(REF), f"Must specify absolute path for {REF}"
    assert os.path.exists(REF), f"Index must exist. Try: samtools faidx {REF}"
else:
    REF = 'Liftoff/ref/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna'

if("gff" in config):
    GFF = config["gff"]
    assert os.path.isabs(GFF), f"Must specify absolute path for {GFF}"
else:
	GFF = 'Liftoff/ref/gencode.v34.primary_assembly.annotation.gff3',

if("regions" in config):
    RGN = config["regions"]
    assert os.path.isabs(RGN), f"Must specify absolute path for {RGN}"
else:
    RGN = ""

wildcard_constraints:
	SM="|".join(SMS),
	ID="\d+",


localrules: run_liftoff


rule get_hg38:
	input:
		fasta = GRCH38,
	output:
		fasta = REF,
		fai = REF+'.fai',
	threads: 1
	resources:
		mem=8
	shell:'''
seqtk seq -l 50 {input.fasta} > {output.fasta}
samtools faidx {output.fasta}
'''

rule get_gff:
	input:
		gff = RGFF,
	output:
		gff = GFF,
	threads: 1
	resources:
		mem=8
	shell:'''
gunzip -c {input.gff} > {output.gff}
'''


rule subset_gff:
    input:
        gff = rules.get_gff.output.gff,
    output:
        gff = temp("Liftoff/ref/{SM}.subset.gff"),
    threads: 1
    resources:
        mem=8
    run:
        if("regions" in config):
            shell("bedtools intersect -wa -a {input.gff} -b {RGN} > {output.gff} ")
        else:
            shell("ln -s {input.gff} {output.gff}")


#################################################################3
# LIFTOFF
#################################################################3

rule run_liftoff:
	input:
		r = REF,
		gff = rules.subset_gff.output.gff,
		t = FASTA,
	output:
		gff = "Liftoff/{SM}.gff3",
		unmapped = "Liftoff/{SM}.unmapped.gff3",
		temp = directory('Liftoff/temp.{SM}')
	threads: 32
	resources:
		mem=8,
	shell:"""
~mvollger/.local/bin/liftoff -dir {output.temp} -sc 0.95 -copies -p {threads} -r {input.r} -t {input.t} -g {input.gff} -o {output.gff} -u {output.unmapped}
"""


rule orf_gff:
    input:
        gff = rules.run_liftoff.output.gff,
        fasta =  FASTA,
    output:
        gff="Liftoff/{SM}.orf_only.gff3",
        all="Liftoff/{SM}.all.gff3",
    threads: 1
    resources:
        mem=8,
    shell:"""
{SDIR}/scripts/AddUniqueGeneIDs_2.py {input.gff} \
    | gffread --adj-stop -C -F -g {input.fasta} /dev/stdin \
    | gffread --keep-comments -F -J -g {input.fasta} /dev/stdin \
    | gffread --keep-comments -F -M -K /dev/stdin \
    | sed 's/CDStopAdjusted/cDStopAdjusted/g' > {output.gff}


{SDIR}/scripts/AddUniqueGeneIDs_2.py {input.gff} \
    | gffread --adj-stop -F -g {input.fasta} /dev/stdin \
    | gffread --keep-comments -F -M -K /dev/stdin \
    | sed 's/CDStopAdjusted/cDStopAdjusted/g' > {output.all}
"""

rule orf_bb:
    input:
        gff = rules.orf_gff.output.gff,
        fasta =  FASTA,
        fai = FASTA + ".fai",
        all = rules.orf_gff.output.all,
    output:
        bb="Liftoff/{SM}.orf_only.bb",
        bed="Liftoff/{SM}.orf_only.bed",
        all="Liftoff/{SM}.all.bed",
        allbb="Liftoff/{SM}.all.bb",
    threads: 1
    resources:
        mem=8,
    shell:"""
gff3ToGenePred -geneNameAttr=gene_name -useName {input.gff} /dev/stdout | \
    genePredToBigGenePred /dev/stdin /dev/stdout | \
    awk -F $'\t' '{{ t = $4; $4 = $13; $13 = t; print; }}' OFS=$'\t' | \
    bedtools sort -i - > {output.bed} 

bedToBigBed -extraIndex=name,name2 -type=bed12+7 -tab -as={SDIR}/templates/bigGenePred.as {output.bed} {input.fai} {output.bb}

gff3ToGenePred -geneNameAttr=gene_name -warnAndContinue -useName {input.all}  /dev/stdout | \
    genePredToBigGenePred /dev/stdin /dev/stdout | \
    awk -F $'\t' '{{ t = $4; $4 = $13; $13 = t; print; }}' OFS=$'\t' | \
    bedtools sort -i - > {output.all}

bedToBigBed -extraIndex=name,name2 -type=bed12+7 -tab -as={SDIR}/templates/bigGenePred.as {output.all} {input.fai} {output.allbb}

"""

rule liftoff:
	input:
		orf = expand(rules.orf_gff.output, SM=[SM]),
		bbs = expand(rules.orf_bb.output, SM=[SM]),


