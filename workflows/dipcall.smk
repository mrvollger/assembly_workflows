import os 
import sys
import re
import re
import pysam 
import pandas as pd
from datetime import date

today = date.today()
DATE =  today.strftime("%Y/%m/%d")

SDIR=os.path.realpath(os.path.dirname(srcdir("env.cfg"))+"/..")
shell.prefix(f"source {SDIR}/env.cfg ; set -eo pipefail; ")

#
#
#
configfile: "dipcall.yaml"
REF=os.path.abspath(config["ref"])
GFF=os.path.abspath(config["gff"])
TBL=os.path.abspath(config["tbl"])
df = pd.read_csv(TBL, sep="\s+")

SMS=df["sample"]
SM_ASM = {}
for idx, row in df.iterrows():
    SM_ASM[row["sample"]] = {"pat": os.path.abspath(row.pat), "mat":os.path.abspath(row.mat)}

HAPS=["hap1","hap2"]

#
# RULES
#
workdir: "dipcall"

wildcard_constraints:
    SM="|".join(SMS),
    HAP="|".join(HAPS),

rule all:
    input:
        vcf = expand("{SM}.dip.vcf.gz", SM=SMS),
        vcf_a = expand("{SM}.dip.anno.vcf.gz", SM=SMS),

rule get_dipcall:
    input:
    output:
        run="dipcall.kit/run-dipcall",
    threads: 1
    shell:"""
wget https://github.com/lh3/dipcall/releases/download/v0.1/dipcall-0.1_x64-linux.tar.bz2
tar -jxf dipcall-0.1_x64-linux.tar.bz2
rm dipcall-0.1_x64-linux.tar.bz2
"""

THREADS=32
MM_OPT=f"-x asm20 -r50k --cs -t {THREADS} -I 8G -2K 1500m"
def get_pat(wc):
    return(SM_ASM[wc.SM]["pat"])
def get_mat(wc):
    return(SM_ASM[wc.SM]["mat"])
def get_hap(wc):
  if(wc.HAP == "hap1"):
      return(SM_ASM[wc.SM]["pat"])
  if(wc.HAP == "hap2"):
      return(SM_ASM[wc.SM]["mat"])

rule dipcall_paf:
    input:
        ref=REF,
        hap=get_hap,
        dipcall=rules.get_dipcall.output,
    output:
        paf = temp("{SM}.{HAP}.paf.gz"),
    log:
        log=temp("{SM}.{HAP}.paf.gz.log"),
    threads:THREADS
    shell:"""
dipcall.kit/minimap2 -c --paf-no-hit {MM_OPT} {input.ref} {input.hap} \
    2> {log.log} | pigz -p {threads} > {output.paf}
"""
 
rule dipcall_sam:
    input:
        ref=REF,
        hap=get_hap,
        dipcall=rules.get_dipcall.output,
    output:
        sam = temp("{SM}.{HAP}.sam.gz"),
    log:
        log=temp("{SM}.{HAP}.sam.gz.log"),
    threads:THREADS
    shell:"""
dipcall.kit/minimap2 -R @RG\\tID:{wildcards.SM}\\tSM:{wildcards.SM} -a {MM_OPT} {input.ref} {input.hap} \
    2> {log.log} | pigz -p {threads} > {output.sam}
"""

rule dipcall_bam:
    input:
        sam = rules.dipcall_sam.output.sam,
    output:
        bam = temp("{SM}.{HAP}.bam"),
    threads:8
    shell:"""
 dipcall.kit/k8 dipcall.kit/dipcall-aux.js samflt {input.sam} | \
   dipcall.kit/samtools sort -m4G --threads {threads} -o {output.bam} -
"""
               
rule dipcall_pair:
    input:
        ref = REF,
        bam1 = "{SM}.hap1.bam",
        bam2 = "{SM}.hap2.bam",
    output:
        vcf = "{SM}.pair.vcf.gz"
    threads: 4
    shell:"""
dipcall.kit/htsbox pileup -q5 -evcf {input.ref} {input.bam1} {input.bam2} | dipcall.kit/htsbox bgzip > {output.vcf}
"""



rule dipcall_vcf:
    input:
        vcf = rules.dipcall_pair.output.vcf,
    output:
        vcf = "{SM}.dip.vcf.gz"
    threads: 4
    shell:"""
dipcall.kit/k8 dipcall.kit/dipcall-aux.js vcfpair {input.vcf} | \
    bcftools reheader -s <(printf "syndip {wildcards.SM}\n") | \
    dipcall.kit/htsbox bgzip > {output.vcf}
"""

rule dipcall_hap_bed:
    input:
        paf=rules.dipcall_paf.output.paf,
    output:
        var="{SM}.{HAP}.var.gz",
        vst="{SM}.{HAP}.vst",
        bed="{SM}.{HAP}.bed",
    threads: 1
    shell: """
gzip -dc {input.paf} | sort -k6,6 -k8,8n | dipkall.kit/k8 dipcall.kit/paftools.js call - \
       2> {output.vst} | gzip > {output.var}
gzip -dc {output.var} | grep ^R | cut -f2- > {output.bed}
"""

rule dipcall_bed:
    input:
        bed1="{SM}.hap1.bed",
        bed2="{SM}.hap2.bed",
    output:
        bed="{SM}.dip.bed"
    threads:1
    shell:"""
dipcall.kit/bedtk isec -m {input.bed1} {input.bed2} > {output.bed}
"""

'''
rule run_dipcall:
    input:
        ref=REF,
        pat=get_pat,
        mat=get_mat,
        dipcall=rules.get_dipcall.output,
    output:
        log1=temp("{SM}.hap1.paf.gz.log"),
        log2=temp("{SM}.hap2.paf.gz.log"),
        log3=temp("{SM}.hap1.sam.gz.log"),
        log4=temp("{SM}.hap2.sam.gz.log"),
        sam1=temp("{SM}.hap1.sam.gz"),
        sam2=temp("{SM}.hap2.sam.gz"),
        paf1=temp("{SM}.hap1.paf.gz"),
        paf2=temp("{SM}.hap2.paf.gz"),
        var1=temp("{SM}.hap1.var.gz"),
        var2=temp("{SM}.hap2.var.gz"),
        bed1=temp("{SM}.hap1.bed"),
        bed2=temp("{SM}.hap2.bed"),
        mak=temp("{SM}.mak"),
        pair=temp("{SM}.pair.vcf.gz"),
        vcf="{SM}.dip.vcf.gz",
        bed="{SM}.dip.bed",
        bam1="{SM}.hap1.bam",
        bam2="{SM}.hap2.bam",
    threads: 32
    shell:"""
dipcall.kit/run-dipcall -t {threads} {wildcards.SM} {input.ref} {input.pat} {input.mat} > {output.mak} 

make -j2 -f {output.mak}

# or for male, requiring PAR regions in BED
# dipcall.kit/run-dipcall -x dipcall.kit/hs38.PAR.bed prefix hs38.fa pat.fa.gz mat.fa.gz > prefix.mak
"""
'''

rule snpeff_setup:
    input:
        ref=REF,
        gff=GFF,
    output:
        fa="data/REF/sequences.fa",
        gff="data/REF/genes.gff",
        jar="snpEff/snpEff.jar",
        config="snpEff.config",
    threads:1
    shell:"""
wget https://snpeff.blob.core.windows.net/versions/snpEff_latest_core.zip
# Unzip file
unzip snpEff_latest_core.zip
rm snpEff_latest_core.zip

mkdir -p data/REF
ln -s {input.ref} {output.fa}
ln -s {input.gff} {output.gff}

cp {SDIR}/templates/snpEff.config.template {output.config} 
"""

rule snpeff_build:
    input:
        config=rules.snpeff_setup.output.config,
    output:
        db="data/REF/snpEffectPredictor.bin"
    threads: 1 
    shell:"""
java -jar snpEff/snpEff.jar build -c {input.config} -gff3 -v REF
"""

rule snpeff_annotate:
    input:
        db=rules.snpeff_build.output.db,
        vcf=rules.dipcall_vcf.output.vcf,
        bed=rules.dipcall_bed.output.bed,
    output:
        vcf="{SM}.dip.anno.vcf.gz",
        html="{SM}.anno.html",
        csv="{SM}.anno.csv",
        genes="{SM}.anno.genes.txt",
    threads: 4
    shell:"""
java -jar snpEff/snpEff.jar \
       -c snpEff.config -s {output.html} -csvStats {output.csv} -fi {input.bed} \
       REF {input.vcf} | bgzip -@ 4 > {output.vcf}
"""




