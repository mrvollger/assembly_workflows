import os 
import sys
import re
import re
import pysam
import pandas as pd

SDIR=os.path.realpath(os.path.dirname(srcdir("env.cfg"))+"/..")
shell.prefix(f"source {SDIR}/env.cfg ; set -eo pipefail; ")

#include: "liftoff.smk"
#include: "sedef.smk"




#
# TODO old remove this 
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
    TO_SPLIT_REF = config["ref"]
    assert os.path.isabs(TO_SPLIT_REF), f"Must specify absolute path for {REF}"
    assert os.path.exists(TO_SPLIT_REF), f"Index must exist. Try: samtools faidx {REF}"
else:
    #TO_SPLIT_REF = 'Liftoff/tmp/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna'
    TO_SPLIT_REF = '/net/eichler/vol26/projects/chm13_t2t/nobackups/assemblies/hg38.chr_only.fa'

#
# 
#
T2T="/net/eichler/vol26/projects/chm13_t2t/nobackups/assemblies/chm13.draft_v1.0_plus38Y.fasta"
T2TSDs="/net/eichler/vol26/projects/chm13_t2t/nobackups/Assembly_analysis/SEDEF/chm13.draft_v1.0_plus38Y.SDs.bed"
HG38='/net/eichler/vol26/projects/chm13_t2t/nobackups/assemblies/hg38.chr_only.fa'

wildcard_constraints:
	SM="|".join(SMS),

localrules: run_aln

#
# minimap2 rules
#
rule align_to_ref:
	input:
		q = FASTA,
		r = TO_SPLIT_REF,
	output:
		paf = "Align/{SM}.paf",
	threads: 128
	resources:
		mem=4,
	shell:"""
minimap2 -I 8G -2K 2000m -t {threads} \
	--secondary=no -c --eqx \
	-x asm20 -s 200000 \
	 -z 10000 -r 50000 \
	--paf-no-hit \
	{input.r}  {input.q} > {output.paf}
"""

rule run_aln:
	input:
		q = TO_SPLIT_REF,
		r = FASTA,
	output:
    #bed = temp("Align/{SM}.bed"),
		split = temp("Align/{SM}.split"),
		bam = "Align/{SM}.bam",
		bai = "Align/{SM}.bam.bai",
	threads: 160
	resources:
		mem=4,
	shell:"""
module load unimap/0.1

bedtools getfasta -fi {input.q} \
    -bed /net/eichler/vol26/projects/chm13_t2t/nobackups/Assembly_analysis/snyteny/GRCh38.snyteny_1Mbp.bed \
    > {output.split}

unimap \
  -R @RG'\\tID':{wildcards.SM}'\\t'SM:{wildcards.SM} \
  -2K 2000m -t {threads} \
	--secondary=no -c --eqx -Y \
	-ax asm20  \
	-r 500000 \
	{input.r}  {output.split} \
  | samtools view -@ 8 -u - \
  | samtools sort -@ 8 -m 16G - > {output.bam}

samtools index {output.bam}
"""
#minimap2 -I 8G 
#bedtools makewindows -g {input.q}.fai -w 5000 > {output.bed}
#bedtools getfasta -fi {input.q} -bed {output.bed} > {output.split}

rule mpileup:
  input:
    bam = rules.run_aln.output.bam,
    r = FASTA,
  output:
    bcf = "Align/{SM}.bcf",
	threads: 24
	resources:
		mem=4,
	shell:"""
parallel -j {threads} --will-cite --colsep '\\t' \
  "bcftools mpileup -r {{1}} -f {input.r} {input.bam} > Align/temp.{wildcards.SM}.{{1}}.bcf" \
  :::: {input.r}.fai 
bcftools concat --threads {threads} Align/temp.{wildcards.SM}.*.bcf > {output.bcf}
rm Align/temp.{wildcards.SM}.*.bcf
"""

rule vcf:
  input:
    bcf = "Align/{SM}.bcf",
    r = FASTA,
  output:
    vcf = "Align/{SM}.vcf",
	threads: 24
	resources:
		mem=4,
	shell:"""
cat {input.bcf} \
  | bcftools call -mv > {output.vcf}
"""

rule vcf_bed:
  input:
    bam = rules.run_aln.output.bam,
    r = FASTA,
    vcf = rules.vcf.output.vcf,
  output:
    snvbed = "Align/{SM}.snv.bed",
    insbed = "Align/{SM}.ins.bed",
    delbed = "Align/{SM}.del.bed",
	threads: 4
	resources:
		mem=4,
	shell:"""
module load bedops/2.4.35
vcf2bed --snvs < {input.vcf} > {output.snvbed}
vcf2bed --insertions < {input.vcf} > {output.insbed}
vcf2bed --deletions < {input.vcf} > {output.delbed}
"""


#
#
# filp reference and query to match previous anlysis 
# TODO, remove and replace with above analysis 
#
#
rule run_aln_flip:
	input:
		r = TO_SPLIT_REF,
		q = FASTA,
	output:
		bed = temp("Align/{SM}.flip.bed"),
		split = temp("Align/{SM}.flip.split"),
		bam = "Align/{SM}.flip.bam",
		bai = "Align/{SM}.flip.bam.bai",
	threads: 128
	resources:
		mem=4,
	shell:"""
bedtools makewindows -g {input.q}.fai -w 5000 > {output.bed}
bedtools getfasta -fi {input.q} -bed {output.bed} > {output.split}
module load unimap/0.1
#minimap2 -I 8G 
unimap \
  -R @RG'\\tID':{wildcards.SM}'\\t'SM:{wildcards.SM} \
  -2K 2000m -t {threads} \
	--secondary=no -c --eqx -Y \
	-ax asm20  \
	-r 50000 \
	{input.r}  {output.split} \
  | samtools view -b - \
  | samtools sort -m 16G - > {output.bam}

samtools index {output.bam}
"""


rule split_bed:
	input:
		bam = rules.run_aln_flip.output.bam,
	output:
		bed = "Align/{SM}.split.bed",
	threads: 1
	resources:
		mem=8,
	shell:"""
~mvollger/projects/utility/samIdentity.py --qbed --split {input.bam} > {output.bed}
"""

FLANK=50000
MERGE_DIST = max(50000, FLANK+1)
MINSIZE=MERGE_DIST
rule inter_sd:
  input:
    #sd = rules.sedef_browser.output.bed,
    sd = T2TSDs,
    bed = rules.split_bed.output.bed,
    fai = TO_SPLIT_REF + ".fai",
  output:
    bed = "Align/{SM}.split.sd.bed",
    nosd = "Align/{SM}.split.nosd.bed",
    flank = "Align/{SM}.flank.bed",
    sdflank = "Align/{SM}.split.sdflank.bed",
  threads: 1
  resources:
    mem=8,
  shell:"""
bedtools intersect -header -f 0.5 -u -a {input.bed} -b {input.sd} > {output.bed}
bedtools intersect -header -v -a {input.bed} -b {input.sd} > {output.nosd}

bedtools sort -i {input.sd} | bedtools merge -d {MERGE_DIST} -i - | \
	awk '$3-$2 > {MINSIZE}{{print $0}}' | \
	bedtools flank -g {input.fai} -b {FLANK} -i - > {output.flank}

bedtools intersect -header -f 0.5 -u -a {input.bed} -b {output.flank} > {output.sdflank}

"""

rule sd_genes:
	input:	
		bed = rules.inter_sd.output.bed,
		genes = T2TSDs, # TODO fix    rules.orf_bb.output.bed,
	output:
		bed = "Align/{SM}.split.sd.genes.bed",
	resources:
		mem=8,
	threads: 1
	shell:"""
bedtools intersect -wa -wb -a {input.genes} -b {input.bed} | bedtools groupby -g 1,2,3,4 -c 12 -o mean > {output.bed}
"""

rule align:
	input:
#		genes = expand(rules.sd_genes.output, SM = [SM]),
		vcfbed = expand(rules.vcf_bed.output.snvbed, SM = [SM]),


