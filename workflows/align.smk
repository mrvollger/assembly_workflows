import os 
import sys
import re
import re
import pysam 
import pandas as pd

SDIR=os.path.realpath(os.path.dirname(srcdir("env.cfg"))+"/..")
shell.prefix(f"source {SDIR}/env.cfg ; set -eo pipefail; ")

include: "liftoff.smk"
include: "sedef.smk"

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
    REF = 'Liftoff/tmp/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna'

wildcard_constraints:
	SM="|".join(SMS),

localrules: run_aln

#
# minimap2 rules
#
rule align_to_ref:
	input:
		q = FASTA,
		r = REF,
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
		q = FASTA,
		r = REF,
	output:
		bed = temp("Align/{SM}.bed"),
		split = temp("Align/{SM}.split"),
		bam = "Align/{SM}.bam",
		bai = "Align/{SM}.bam.bai",
	threads: 128
	resources:
		mem=4,
	shell:"""
bedtools makewindows -g {input.q}.fai -w 5000 > {output.bed}
bedtools getfasta -fi {input.q} -bed {output.bed} > {output.split}

minimap2 -I 8G -2K 2000m -t {threads} \
	--secondary=no -c --eqx -Y \
	-ax asm20  \
	-r 50000 \
	{input.r}  {output.split} | samtools view -b - | samtools sort - > {output.bam}

samtools index {output.bam}
"""

rule split_bed:
	input:
		bam = rules.run_aln.output.bam,
		bai = rules.run_aln.output.bai,
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
		sd=rules.sedef_browser.output.bed,
		bed = rules.split_bed.output.bed,
		fai = FASTA + ".fai",
	output:
		bed = "Align/{SM}.split.sd.bed",
		nosd = "Align/{SM}.split.nosd.bed",
		flank = "Align/{SM}.flank.bed",
		sdflank = "Align/{SM}.split.sdflank.bed",
	threads: 1
	resources:
		mem=8,
	shell:"""
bedtools intersect -header -f 0.5 -wa -a {input.bed} -b {input.sd} > {output.bed}
bedtools intersect -header -v  -a {input.bed} -b {input.sd} > {output.nosd}

bedtools sort -i {input.sd} | bedtools merge -d {MERGE_DIST} -i - | \
	awk '$3-$2 > {MINSIZE}{{print $0}}' | \
	bedtools flank -g {input.fai} -b {FLANK} -i - > {output.flank}

bedtools intersect -header -f 0.5 -wa -a {input.bed} -b {output.flank} > {output.sdflank}

"""

rule sd_genes:
	input:	
		bed = rules.inter_sd.output.bed,
		genes = rules.orf_bb.output.bed,
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
		genes = expand(rules.sd_genes.output, SM = [SM]),


