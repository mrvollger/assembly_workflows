import os 
import sys
from Bio import SeqIO
import itertools
import re
import tempfile
import glob
import pandas as pd
import numpy as np 

ENV = os.path.dirname(workflow.snakefile) + "/../env.cfg"
shell.prefix("source {ENV}; set -eo pipefail; ")


SMS=list(config.keys())


df = pd.read_csv("samples_trio.tbl", sep="\t")
df.replace(np.nan, '', regex=True, inplace=True)
PARENTS = ["pat", "mat"]
config={}

NOT_SMS=[]
TRIO_SMS=[]
for idx, row in df.iterrows():
    if(row["pat"] == "" or row["mat"] == ""):
        config[row["sample"]] = {"fofn":row["fofn"]}
        NOT_SMS.append(row["sample"])
    else:
        config[row["sample"]] = {"pat":row["pat"],"mat":row["mat"], "fofn":row["fofn"]}
        TRIO_SMS.append(row["sample"])
SMS = list(config.keys())

print(SMS)

wildcard_constraints:
	SM= "|".join(SMS),
	PAR = "|".join(PARENTS),

workdir: "hifiasm_out"

MAXT=32
TMP = "/tmp/mvollger/hifi_asm"
shell("mkdir -p {TMP} logs")

rule all:
	input:
		mats = expand("{SM}.pat.fa.gz", SM=TRIO_SMS),
		pats = expand("{SM}.mat.fa.gz", SM=TRIO_SMS),
		pris = expand("{SM}.pri.fa.gz", SM=NOT_SMS),
		alts = expand("{SM}.alt.fa.gz", SM=NOT_SMS),


def get_illumina(wc):
    #print(config[wc.SM][wc.PAR])
    return(config[wc.SM][wc.PAR] )

rule yac:
    input:
        illumina = get_illumina, 
    output:
        yak = temp(TMP+"/{SM}.{PAR}.yak"),
        fasta = temp(TMP+"/{SM}.{PAR}.ill.fasta"),
    threads: MAXT
    log: "logs/yak_{SM}_{PAR}.log"
    benchmark: "logs/yak_{SM}_{PAR}.b"
    run:
        if(input.illumina[-5:]==".cram"):
            shell("samtools fasta -@ {threads} {input.illumina} > {output.fasta}")
            shell("yak count -k31 -b37 -t {threads} -o {output.yak} {output.fasta} &> {log}")
        else:
            shell("touch {output.fasta}")
            shell("yak count -k31 -b37 -t {threads} -o {output.yak} {input.illumina} &> {log}")


def get_fofn(wc):
    return(config[wc.SM]["fofn"])

rule reads:
    input:
        fofn=get_fofn,
    output:
        reads = temp(TMP + "/{SM}.reads"),
    threads: MAXT
    priority: 50
    run:
        shell("> {output.reads}")
        for read in open(input.fofn):
            read=read.strip()
            if(read[-4:]==".bam"):
                shell("samtools fasta -@ {threads} {read} >> {output.reads}")
            else:
                shell("cat {read} >> {output.reads}")


rule trio_hifiasm:
    input:
        reads=rules.reads.output.reads,
        pat_yak = TMP+"/{SM}.pat.yak",
        mat_yak = TMP+"/{SM}.mat.yak",
        mat_ill = TMP+"/{SM}.mat.ill.fasta",
        pat_ill = TMP+"/{SM}.pat.ill.fasta",
    output:
        #dipy = "all_out/{SM}.asm.dip.r_utg.gfa",
        #dipn = "all_out/{SM}.asm.dip.r_utg.noseq.gfa",
        hap1 = "all_out/{SM}.asm.hap1.p_ctg.gfa",
        #hap1n = "all_out/{SM}.asm.hap1.p_ctg.noseq.gfa",
        hap2 = "all_out/{SM}.asm.hap2.p_ctg.gfa",
        #hap2n = "all_out/{SM}.asm.hap2.p_ctg.noseq.gfa",
        #ec = "all_out/{SM}.asm.ec.bin",
        #reverse = "all_out/{SM}.asm.ovlp.reverse.bin",
        #source = "all_out{SM}.asm.ovlp.source.bin",
    threads:MAXT
    priority: 100
    log: "logs/hifiasm_{SM}.log"
    benchmark: "logs/hifiasm_{SM}.b"
    shell:"""
LOG=$(readlink -f {log})
pushd all_out
hifiasm -o {wildcards.SM}.asm -t {threads} \
	-1 {input.pat_yak} -2 {input.mat_yak} {input.reads} &> $LOG
popd
"""

rule trio_make_fasta:	
	input:
		hap1 = rules.trio_hifiasm.output.hap1,
		hap2 = rules.trio_hifiasm.output.hap2,
	output:
	    pat = "{SM}.pat.fa.gz",
		mat = "{SM}.mat.fa.gz",
	    fpat = "{SM}.pat.fa.gz.fai",
		fmat = "{SM}.mat.fa.gz.fai",
	threads:4 
	shell:"""
awk '/^S/{{print ">"$2"_pat\\n"$3}}' {input.hap1} | seqtk seq -l 80 | bgzip -@ {threads} > {output.pat}
awk '/^S/{{print ">"$2"_mat\\n"$3}}' {input.hap2} | seqtk seq -l 80 | bgzip -@ {threads} > {output.mat}
samtools faidx {output.pat}
samtools faidx {output.mat}
"""


rule hifiasm:
    input:
        reads=rules.reads.output.reads,
    output:
        pri = "all_out/{SM}.asm.p_ctg.gfa",
        alt = "all_out/{SM}.asm.a_ctg.gfa",
        #hap1n = "all_out/{SM}.asm.hap1.p_ctg.noseq.gfa",
        #dipn = "all_out/{SM}.asm.dip.r_utg.noseq.gfa",
        #dipy = "all_out/{SM}.asm.dip.r_utg.gfa",
        #hap2n = "all_out/{SM}.asm.hap2.p_ctg.noseq.gfa",
        #ec = "all_out/{SM}.asm.ec.bin",
        #reverse = "all_out/{SM}.asm.ovlp.reverse.bin",
        #source = "all_out{SM}.asm.ovlp.source.bin",
    threads:MAXT
    priority: 100
    log: "logs/hifiasm_{SM}.log"
    benchmark: "logs/hifiasm_{SM}.b"
    shell:"""
LOG=$(readlink -f {log})
pushd all_out
hifiasm -o {wildcards.SM}.asm -t {threads} {input.reads} &> $LOG
popd
"""

rule make_fasta:	
	input:
		pri = rules.hifiasm.output.pri,
		alt = rules.hifiasm.output.alt,
	output:
		pri = "{SM}.pri.fa.gz",
		alt = "{SM}.alt.fa.gz",
		fpri = "{SM}.pri.fa.gz.fai",
		falt = "{SM}.alt.fa.gz.fai",
	threads:4
	shell:"""
awk '/^S/{{print ">"$2"_pri\\n"$3}}' {input.pri} | seqtk seq -l 80 | bgzip -@ {threads} > {output.pri} 
awk '/^S/{{print ">"$2"_alt\\n"$3}}' {input.alt} | seqtk seq -l 80 | bgzip -@ {threads} > {output.alt} 
samtools faidx {output.pri}
samtools faidx {output.alt}
"""




