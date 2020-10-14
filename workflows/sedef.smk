import os 
import sys
import re
import re
import pysam 
import pandas as pd
from datetime import date
from pybedtools import BedTool

today = date.today()
DATE =  today.strftime("%Y/%m/%d")

SDIR=os.path.realpath(os.path.dirname(srcdir("env.cfg"))+"/..")
shell.prefix(f"source {SDIR}/env.cfg ; set -eo pipefail; ")

FASTA = config["fasta"] 
assert os.path.isabs(FASTA), f"Must specify absolute path for {FASTA}"
FAI = FASTA + ".fai"
assert os.path.exists(FAI), f"Index must exist. Try: samtools faidx {FASTA}"

SM = "asm"
if("sample" in config): SM = config["sample"]
THREADS=16
if("threads" in config): THREADS = config["threads"]
SMS = [SM]

include: "mask.smk"
RMBED = rules.RepeatMasker.input.bed
TRFBED = rules.trf.input.bed

wildcard_constraints:
	SM="|".join(SMS),

localrules: sedef, run_sedef


#
# Make masked version of the fasta for sedef
#
MASKED = os.path.abspath(f"SEDEF/{SM}_masked.fasta")
rule sedef_masked_fasta:
	input:
		fasta = FASTA, 
		trf = TRFBED,
		rm = RMBED,
	output:
		bed = temp(f"SEDEF/{SM}.tmp.msk.bed"),
		bed2 = temp(f"SEDEF/{SM}.tmp2.msk.bed"),
		fasta = MASKED,
		fai = MASKED+".fai",
	resources:
		mem=8,
	threads:1
	shell:"""
# very it is very common to find a 27 base pair gap between alpha sat annotations
# similarly it is commone to find a 49 bp gap in anotations in HSAT arrays
# therefor I merge some of the sat features from rm
# additionally there are some huge HSAT arrays that are not annotated as  HSAT by rm and instead are
# Simple_repeat


grep Alpha {input.rm} | bedtools merge -d 35 -i - > {output.bed}  
grep HSAT {input.rm} | bedtools merge -d 75 -i - >> {output.bed}  
# merge large simple_repeats
grep Simple_repeat {input.rm} | \
        awk '$3-$2 > 1000 {{print $0}}' | \
        bedtools merge -d 100 -i - | \
        bedtools slop -b 100 -g {input.fasta}.fai -i - >> {output.bed}  

# combine custome merges with trf and rm
cat {input.trf} {input.rm} | cut -f 1-3 >> {output.bed}  

# make large merge where large entries are allowed to merge 
# together further
cut -f 1-3 {output.bed} | bedtools sort -i - | bedtools merge -i - | \
        awk '$3-$2 > 2000 {{print $0}}' | \
        bedtools merge -d 100 -i - > {output.bed2}  


cut -f 1-3 {output.bed} {output.bed2} | bedtools sort -i - | bedtools merge -i - | \
	seqtk seq -l 50 -M /dev/stdin {input.fasta} > {output.fasta}

samtools faidx {output.fasta}
"""
		

#
# sedef rules
#
rule run_sedef:
	input:
		fasta = rules.sedef_masked_fasta.output.fasta,
	output:
		tmpf = temp(f"/tmp/sedef/sedef_{SM}.fasta"),
		fai = temp(f"/tmp/sedef/sedef_{SM}.fasta.fai"),
		bed = f"SEDEF/{SM}/final.bed"
	resources:
		mem=8
	threads: 120
	shell:"""
cp {input.fasta} {output.tmpf} 
samtools faidx {output.tmpf}
sedef.sh -f -o $(dirname {output.bed}) -j {threads} {output.tmpf}
"""

# unnessisary can calculate from uppercase matches.
rule count_sat_sedef:
    input:
        rm = RMBED,
        trf = TRFBED,
        sedef = rules.run_sedef.output.bed,
    output:
        bed = f"SEDEF/{SM}/final.sat.count.bed",
        tmp = f"SEDEF/{SM}/tmp_final.sat.count.bed",
    resources:
        mem=8,
    threads: 1
	shell:"""
cat {input.rm} | grep Satellite | cut -f 1-3 | bedtools sort -i - | bedtools merge -i - | \
		bedtools coverage -header -a {input.sedef} -b - > {output.tmp}
grep "^#" {input.sedef} > {output.bed}
cat {output.tmp} >> {output.bed}
sed -i '1{{s/$/\tcount_ovls\tsat_bases\ttotal_bases\tsat_coverage/}}' {output.bed}
"""

rule sedef_browser:
	input:
		bed = rules.count_sat_sedef.output.bed,
	output:
		bed = f"SEDEF/{SM}.SDs.bed",
		lowid = f"SEDEF/{SM}.SDs.lowid.bed",
		html = f"SEDEF/{SM}.SDs.html",
	resources:
		mem=8,
	threads: 1
	run:
		html = open(f"{SDIR}/templates/sedef.html").read()
		open(output["html"], "w+").write(html.format(DATE=DATE, SM=SM))
		shell("{SDIR}/scripts/sedef_to_bed.py --sat 0.70 {input.bed} {output.bed} {output.lowid}")

rule sedef_bb:
    input:
        bed = rules.sedef_browser.output.bed,
        lowid = rules.sedef_browser.output.lowid,
        fai = FASTA + ".fai",
    output:
        bb = rules.sedef_browser.output.bed + ".bb",
        lowid = rules.sedef_browser.output.lowid+".bb",
    resources:
        mem=8,
    threads: 1
    shell:"""
bedToBigBed -type=bed9+27 -tab -as={SDIR}/templates/sedef.as {input.bed} {input.fai} {output.bb}
bedToBigBed -type=bed9+27 -tab -as={SDIR}/templates/sedef.as {input.lowid} {input.fai} {output.lowid}
"""


rule find_cen:
    input:
        rm = RMBED,
    output:
        tmp = temp(f"SEDEF/tmp.{SM}.cen.bed"),
        bed = f"SEDEF/{SM}.cen.bed",
    resources:
        mem=8,
    threads:1
    shell:"""
grep 'ALR/Alpha' {input.rm} | \
        bedtools merge -d 50 -i - | \
        awk '($3-$2)>50' | \
        bedtools merge -d 100 -i - | \
        awk '($3-$2)>10000' | \
        bedtools merge -d 200000 -i - | \
        awk '{{print $0"\t"$3-$2}}' | \
        sort -k 1,1 -k4,4n > {output.tmp} 

bedtools groupby -g 1 -c 2,3 -o last,last -i {output.tmp} > {output.bed} 
"""

rule sum_sedef:
    input:
        bed = rules.sedef_browser.output.bed,
        lowid = rules.sedef_browser.output.lowid,
        cen = rules.find_cen.output.bed,
        fai = MASKED+".fai",
    output:
        xlsx = f"SEDEF/{SM}.sedef.summary.xlsx",
        lowid = f"SEDEF/{SM}.sedef.lowid.summary.xlsx",
    resources:
        mem=8,
    threads:1
    shell:"""
{SDIR}/scripts/sedef_summary.py --fai {input.fai} --cen {input.cen} --excel {output.xlsx} {input.bed}
{SDIR}/scripts/sedef_summary.py --fai {input.fai} --cen {input.cen} --excel {output.lowid} {input.lowid}
"""

rule sedef:
	input:
		rules.sedef_browser.output,
		rules.sedef_bb.output,
        rules.find_cen.output.bed,
        rules.sum_sedef.output,







