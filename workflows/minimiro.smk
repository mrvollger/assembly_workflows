import os 
import sys
import re
import re
import pandas as pd 



configfile: "minimiro.yaml"
#SDIR=os.path.dirname(workflow.snakefile)
SDIR=os.path.realpath(os.path.dirname(srcdir("env.cfg"))+"/..")

DEBUG=True

shell.prefix(f"source {SDIR}/env.cfg ; set -eo pipefail; ")

SMS = list(config.keys())
SEQS=["ref", "query"]
SCORES=config["scores"];  SMS.remove("scores") # minimum peak dp socre 

RS = {}
RGNS = {}
QS = {}
QRGNS = {}
RCS = {}
GENES = {}
GFF = {}
BAMS = {}
for SM in SMS:
	RS[SM] = os.path.abspath( config[SM]["ref"] ) # reference seqeunce 
	assert os.path.exists(RS[SM]+".fai")
	RGNS[SM]	=	config[SM]["regions"] # region(s) to pull from reference 
	QS[SM]		=	os.path.abspath( config[SM]["query"] )# query sequence
	assert os.path.exists(QS[SM]+".fai")
	QRGNS[SM]	=	config[SM]["queryregions"] # region(s) to pull from query
	if("rc" in config[SM]):
		RCS[SM] = config[SM]["rc"] 
	else: 
		RCS[SM] = False
	GFF[SM] = os.path.abspath( config[SM]["gff"] )

	if("bam" in config[SM]):
		BAMS[SM] = config[SM]["bam"]

wildcard_constraints:
	SEQ="|".join(SEQS),
	SM="|".join(SMS),
	SCORE="\d+",

rule all:
	input:
		pdf	= expand("minimiro_smk_out/{SM}_{SCORE}_aln.pdf", SM=SMS, SCORE=SCORES),
		png	= expand("minimiro_smk_out/{SM}_coverage.png", SM=list(BAMS.keys())),

# delete if not debug
def tempd(fname):
	if(DEBUG):
		return(fname)
	return(temp(fname))

def get_ref(wildcards):
	SM = str(wildcards.SM)
	return(RS[SM])

def get_gff(wildcards):
	SM = str(wildcards.SM)
	return(GFF[SM])

def get_query(wildcards):
	SM = str(wildcards.SM)
	return(QS[SM])

def get_ref_rgns(wildcards):
	SM = str(wildcards.SM)
	return( " ".join( RGNS[SM] ) )
def get_query_rgns(wildcards):
	SM = str(wildcards.SM)
	return( " ".join( QRGNS[SM] ) )

def get_rc(wildcards):
	SM = str(wildcards.SM)
	return(RCS[SM])

rule get_rgns:
	input:
		ref=get_ref, 
		query=get_query,
	output:
		ref = tempd("temp/{SM}_ref.fasta"),
		query = tempd("temp/{SM}_query.fasta"),
	params:
		rgns = get_ref_rgns,
		qrgns = get_query_rgns,
		rc = get_rc,
	threads:1
	run:
		shell("samtools faidx {input.ref} {params.rgns} > {output.ref}")	
		if(params["rc"]):
			shell("samtools faidx {input.query} {params.qrgns} | seqtk seq -r - > {output.query}")	
		else:	
			shell("samtools faidx {input.query} {params.qrgns} > {output.query}")	

rule RepeatMasker:
	input:
		fasta = "temp/{SM}_{SEQ}.fasta",
	output:
		out = tempd("temp/{SM}_{SEQ}.fasta.out"),
		cat = tempd("temp/{SM}_{SEQ}.fasta.cat"),
		tbl = tempd("temp/{SM}_{SEQ}.fasta.tbl"),
		msk = tempd("temp/{SM}_{SEQ}.fasta.masked"),
	resources:
		mem=8,
	threads:8
	shell:"""
RepeatMasker \
	-e ncbi \
	-species human \
	-dir $(dirname {input.fasta}) \
	-pa {threads} \
	{input.fasta}
"""


rule DupMasker:
	input:
		fasta = "temp/{SM}_{SEQ}.fasta",
		out = rules.RepeatMasker.output.out,
	output:
		dups = "temp/{SM}_{SEQ}.fasta.duplicons",
	threads:8
	shell:"""
DupMaskerParallel -pa {threads} -engine ncbi \
	{input.fasta}
"""
#-pa {threads} \

rule DupMaskerColor:
	input:
		dups = rules.DupMasker.output.dups,
	output:
		dupcolor = "temp/{SM}_{SEQ}.fasta.duplicons.extra",
	shell:"""
{SDIR}/scripts/DupMask_parserV6.pl -i {input.dups} -E -o {output.dupcolor}
"""

#
# rules to get genes
#
def get_genes(wildcards):
    SM = str(wildcards.SM)
    return(GENES[SM])

def get_ref_bed(wildcards):
    SM = str(wildcards.SM)
    rtn = []
    for rgn in RGNS[SM]:
        match = re.match("(.+):(\d+)-(\d+)", rgn.strip().replace(",",""))
        if(match):
            rtn.append('{}\t{}\t{}\n'.format(*match.groups()))
        else:
            rtn.append('{}\t{}\t{}\n'.format(rgn, 0, 1000000000))
    return( rtn )


rule clean_gff:
    input:
        ref = get_ref,
        gff = get_gff,
    output:
        gff = temp("temp/{SM}.small.gff"),
        tmpgff = temp("temp/{SM}.tmpsmall.gff"),
        rgns = temp("temp/{SM}.rgns.bed"),
        bed = temp("temp/{SM}.small.bed"),
    params:
        bed = get_ref_bed,
    threads: 8
    run:
        open(output.rgns, "w+").write("\n".join(params.bed) )
        shell("""
                #echo '##gff-version 3' > {output.gff} 
                bedtools intersect -header -f 1.0 -a {input.gff} -b {output.rgns} | \
                    {SDIR}/bin/gffread-0.12.3.Linux_x86_64/gffread --keep-genes -F > {output.gff} """)
        shell("""
            {SDIR}/scripts/AddUniqueGeneIDs_2.py {output.gff} \
                    | {SDIR}/bin/gffread-0.12.3.Linux_x86_64/gffread --keep-comments --adj-stop -C -F -g {input.ref} /dev/stdin \
                    | {SDIR}/bin/gffread-0.12.3.Linux_x86_64/gffread --keep-comments -F -J -g {input.ref} /dev/stdin \
                    | {SDIR}/bin/gffread-0.12.3.Linux_x86_64/gffread --keep-comments -F -M -K /dev/stdin \
                    | sed 's/CDStopAdjusted/cDStopAdjusted/g' > {output.tmpgff}""")
        
        gff=open(output.tmpgff).readlines()
        count=0
        for line in gff:
            if(line[0]=="#"): continue
            count+=1
        if(count==0):
            shell("touch {output.bed}")
        else:
            shell("""
            gff3ToGenePred -geneNameAttr=gene_name -useName {output.tmpgff} /dev/stdout | \
                    genePredToBigGenePred /dev/stdin /dev/stdout | \
                    awk -F $'\t' '{{ t = $4; $4 = $13; $13 = t; print; }}' OFS=$'\t' | \
                    bedtools sort -i - > {output.bed} 
                """)

rule get_ref_genes:
    input:
        bed = rules.clean_gff.output.bed,
    output:
        bed12 = "temp/{SM}.ref.genes.12.bed",
        tmp = temp("temp/tmp.{SM}.ref.genes.bed"),
    params:
        bed = get_ref_bed,
    run:
        # make a bed file with all the corrdiantes in the correct space
        rtn = ""
        for bed in params["bed"]:
            # subset gene bed, and then fix coords
            shell("""bedtools intersect -header -f 1.0 -a {input.bed} -b <(printf "{bed}") | bedtools sort -i - > {output.tmp} """)

            chrm, start, end = bed.split(); start = int(start);
            name = f"{chrm}:{start}-{end}"
            for line in open(output["tmp"]).readlines():
                t = line.strip().split()
                t[0] = name
                t[1] = int(t[1]) - start	
                t[2] = int(t[2]) - start	
                rtn += (11*"{}\t" + "{}\n").format(*t)

        open(output["bed12"], "w+").write(rtn)
            

rule query_genes:
    input:
        fasta = rules.get_rgns.output.query,
        ref = get_ref,
        #gff = get_gff,
        rgns = rules.clean_gff.output.rgns,
        gff = rules.clean_gff.output.gff,
    output:
        bed12 = "Liftoff/{SM}.orf_only.bed",
    params:
        bed = get_ref_bed,
    threads: 8
    run:
        shell("samtools faidx {input.fasta}")
        shell("""snakemake -s {SDIR}/workflows/liftoff.smk \
                -j {threads} -p \
                {output.bed12} \
                --config \
                    fasta=$(readlink -f {input.fasta}) \
                    ref={input.ref} \
                    gff=$(readlink -f {input.gff}) \
                    sample={wildcards.SM} \
                --nolock""")


#
# make the alignments and the miropeats 	
#
def get_score(wildcards):
	return( int(str(wildcards.SCORE)))

rule minimap2:
	input:
		ref = rules.get_rgns.output.ref,
		query = rules.get_rgns.output.query,
	output:
		paf = tempd("temp/{SM}_{SCORE}_aln.paf"),
	params:
		score = get_score,
	threads: 16
	shell:"""
# YOU HAVE TO INCLUDE --cs FOR MINIMIRO TO WORK
minimap2 -x asm20 -r 200000 -s {params.score} -p 0.01 -N 1000 --cs {input.ref} {input.query} > {output.paf}
"""


rule minimiro:
    input:
        paf = rules.minimap2.output.paf,
        rmout = expand("temp/{{SM}}_{SEQ}.fasta.out", SEQ=SEQS),
        dmout = expand("temp/{{SM}}_{SEQ}.fasta.duplicons.extra", SEQ=SEQS),
        genes = rules.get_ref_genes.output.bed12,
        query_genes = rules.query_genes.output.bed12,
    output:
        ps	= "temp/{SM}_{SCORE}_aln.ps",
        pdf	= "minimiro_smk_out/{SM}_{SCORE}_aln.pdf",
    threads: 1
	shell:"""
{SDIR}/scripts/minimiro.py --paf {input.paf} \
	--rm {input.rmout} \
	--dm {input.dmout} \
    --bed <(cut -f 1-12 {input.genes}) <(cut -f 1-12 {input.query_genes}) \
	--bestn 1000 \
	-o {output.ps} && \
	ps2pdf {output.ps} {output.pdf}
"""

def get_bam(wc):
	SM = str(wc.SM)
	return(BAMS[SM])

rule coverage:
	input:
		bam = get_bam,
	output:
		png	= "minimiro_smk_out/{SM}_coverage.png",
	params:
		rgn = get_query_rgns, 
	threads: 1
	shell:"""
NucPlot.py {input.bam} {output.png} --regions {params.rgn} --height 4 --width 16
"""






