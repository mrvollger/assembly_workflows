import os
import sys
import re
import re
import pandas as pd


configfile: "minimiro.yaml"


# SDIR=os.path.dirname(workflow.snakefile)
SDIR = os.path.realpath(os.path.dirname(srcdir("env.cfg")) + "/..")
MAXT = 16
DEBUG = True

shell.prefix(f"source {SDIR}/env.cfg ; set -eo pipefail; ")

SEQS = ["ref", "query"]
SCORES = config.pop("scores", [1000])
SIM = config.pop("simple", False)
GENES = config.pop("genes", True)
GREP = config.pop("grep", ".")
SECONDARY = config.pop("secondary", "yes")
INFER = config.pop("infer", True)
MAX_INDEL = config.pop("max_indel", 1e6)
SMS = list(config.keys())
print(SMS)

simple_param = ""
if SIM:
    simple_param = " --simple "

RS = {}
RGNS = {}
QS = {}
QRGNS = {}
RCS = {}
GENEBED = {}
GFF = {}
BAMS = {}
for SM in SMS:
    RS[SM] = os.path.abspath(config[SM]["ref"])  # reference seqeunce
    assert os.path.exists(RS[SM] + ".fai")
    RGNS[SM] = config[SM]["regions"]  # region(s) to pull from reference
    QS[SM] = os.path.abspath(config[SM]["query"])  # query sequence
    print(QS[SM])
    assert os.path.exists(QS[SM] + ".fai")
    QRGNS[SM] = config[SM]["queryregions"]  # region(s) to pull from query
    if "rc" in config[SM]:
        RCS[SM] = config[SM]["rc"]
    else:
        RCS[SM] = False

    if "gff" in config[SM]:
        GFF[SM] = os.path.abspath(config[SM]["gff"])

    if "bam" in config[SM]:
        BAMS[SM] = config[SM]["bam"]

    if "genebed" in config[SM]:
        GENEBED[SM] = config[SM]["genebed"]


workdir: "minimiro"


wildcard_constraints:
    SEQ="|".join(SEQS),
    SM="|".join(SMS),
    SCORE="\d+",


localrules:
    DupMaskerColor,
    get_ref_genes,
    all,


rule all:
    input:
        pdf=expand("miro_figures/{SM}_{SCORE}_aln.pdf", SM=SMS, SCORE=SCORES),
        png=expand("miro_figures/{SM}_coverage.png", SM=list(BAMS.keys())),


# delete if not debug
def tempd(fname):
    if DEBUG:
        return fname
    return temp(fname)


def get_ref(wildcards):
    SM = str(wildcards.SM)
    return RS[SM]


def get_gff(wildcards):
    SM = str(wildcards.SM)
    if SM in GFF:
        return GFF[SM]
    else:
        return get_ref(wildcards)  # this is a fake git so it can be skipped


def get_query(wildcards):
    SM = str(wildcards.SM)
    return QS[SM]


def get_ref_rgns(wildcards):
    SM = str(wildcards.SM)
    return " ".join(RGNS[SM])


def get_query_rgns(wildcards):
    SM = str(wildcards.SM)
    return " ".join(QRGNS[SM])


def get_rc(wildcards):
    SM = str(wildcards.SM)
    return RCS[SM]


rule get_rgns:
    input:
        ref=get_ref,
        query=get_query,
    output:
        ref=tempd("temp_minimiro/{SM}_ref.fasta"),
        query=tempd("temp_minimiro/{SM}_query.fasta"),
    params:
        rgns=get_ref_rgns,
        qrgns=get_query_rgns,
        rc=get_rc,
    threads: 8
    resources:
        mem=4,
        hrs=24,
    run:
        shell("samtools faidx {input.ref} {params.rgns} > {output.ref}")
        if params["rc"]:
            shell(
                """samtools faidx {input.query} {params.qrgns} | seqtk seq -r - | awk '{{if (NR == 1) print ">{wildcards.SM}"; else print}}' > {output.query}"""
            )
        elif INFER:
            shell(
                """
        minimap2 -t {threads} -ax asm20 -r 200000 --eqx -Y \
            {output.ref} <(samtools faidx {input.query} {params.qrgns}) | \
            samtools view -F 2304 | \
            awk '{{OFS="\\t"; print ">{SM}_{SEQ}""\\n"$10}}' - | \
        seqtk seq -l 60 > {output.query} """
            )
        else:
            shell(
                """samtools faidx {input.query} {params.qrgns} | awk '{{if (NR == 1) print ">{wildcards.SM}"; else print}}' > {output.query}"""
            )

        shell("samtools faidx {output.query}")
        shell("samtools faidx {output.ref}")


rule RepeatMasker:
    input:
        fasta="temp_minimiro/{SM}_{SEQ}.fasta",
    output:
        #cat = tempd("temp_minimiro/{SM}_{SEQ}.fasta.cat"),
        #tbl = tempd("temp_minimiro/{SM}_{SEQ}.fasta.tbl"),
        #msk = tempd("temp_minimiro/{SM}_{SEQ}.fasta.masked"),
        out=tempd("temp_minimiro/{SM}_{SEQ}.fasta.out"),
    resources:
        mem=8,
        hrs=24,
    threads: MAXT
    shell:
        """
        RepeatMasker \
            -e ncbi \
            -species human \
            -dir $(dirname {input.fasta}) \
            -pa {threads} \
            {input.fasta}
        """


rule DupMasker:
    input:
        fasta="temp_minimiro/{SM}_{SEQ}.fasta",
        out=rules.RepeatMasker.output.out,
    output:
        dups="temp_minimiro/{SM}_{SEQ}.fasta.duplicons",
    threads: MAXT
    resources:
        mem=8,
        hrs=24,
    shell:
        """
        samtools faidx {input.fasta}
        DupMaskerParallel -pa {threads} -engine ncbi \
            {input.fasta}
        """


# -pa {threads} \


rule DupMaskerColor:
    input:
        dups=rules.DupMasker.output.dups,
    output:
        dupcolor="temp_minimiro/{SM}_{SEQ}.fasta.duplicons.extra",
    shell:
        """
        {SDIR}/scripts/DupMask_parserV6.pl -i {input.dups} -E -o {output.dupcolor}
        """


#
# rules to get genes
#
# def get_genes(wildcards):
#    SM = str(wildcards.SM)
#    return(GENES[SM])


def get_ref_bed(wildcards):
    SM = str(wildcards.SM)
    rtn = []
    for rgn in RGNS[SM]:
        match = re.match("(.+):(\d+)-(\d+)", rgn.strip().replace(",", ""))
        if match:
            rtn.append("{}\t{}\t{}\n".format(*match.groups()))
        else:
            rtn.append("{}\t{}\t{}\n".format(rgn, 0, 1000000000))
    return rtn


rule clean_gff:
    input:
        ref=get_ref,
        gff=get_gff,
    output:
        gff=temp("temp_minimiro/{SM}.small.gff"),
        tmpgff=temp("temp_minimiro/{SM}.tmpsmall.gff"),
        rgns=temp("temp_minimiro/{SM}.rgns.bed"),
        bed=temp("temp_minimiro/{SM}.small.bed"),
    params:
        bed=get_ref_bed,
    threads: 8
    resources:
        mem=4,
        hrs=24,
    run:
        if not GENES:
            shell("touch {output}")
        else:
            open(output.rgns, "w+").write("\n".join(params.bed))
            shell(
                """
        #echo '##gff-version 3' > {output.gff} 
        bedtools intersect -header -f 1.0 -a {input.gff} -b {output.rgns} | \
        {SDIR}/bin/gffread-0.12.3.Linux_x86_64/gffread --keep-genes -F > {output.gff} """
            )
            shell(
                """
        {SDIR}/scripts/AddUniqueGeneIDs_2.py {output.gff} \
                | {SDIR}/bin/gffread-0.12.3.Linux_x86_64/gffread --keep-comments --adj-stop -C -F -g {input.ref} /dev/stdin \
                | {SDIR}/bin/gffread-0.12.3.Linux_x86_64/gffread --keep-comments -F -J -g {input.ref} /dev/stdin \
                | {SDIR}/bin/gffread-0.12.3.Linux_x86_64/gffread --keep-comments -F -M -K /dev/stdin \
        | sed 's/CDStopAdjusted/cDStopAdjusted/g' > {output.tmpgff}"""
            )

            gff = open(output.tmpgff).readlines()
            count = 0
            for line in gff:
                if line[0] == "#":
                    continue
                count += 1
            if count == 0:
                shell("touch {output.bed}")
            else:
                shell(
                    """
        gff3ToGenePred -warnAndContinue -geneNameAttr=gene_name -useName {output.tmpgff} /dev/stdout | \
                genePredToBigGenePred /dev/stdin /dev/stdout | \
                awk -F $'\t' '{{ t = $4; $4 = $13; $13 = t; print; }}' OFS=$'\t' | \
                bedtools sort -i - > {output.bed} 
        """
                )


rule get_ref_genes:
    input:
        bed=rules.clean_gff.output.bed,
    output:
        bed12="temp_minimiro/{SM}.ref.genes.12.bed",
        tmp=temp("temp_minimiro/tmp.{SM}.ref.genes.bed"),
    params:
        bed=get_ref_bed,
    run:
        if not GENES:
            shell("touch {output}")
        else:
            # make a bed file with all the corrdiantes in the correct space
            rtn = ""
            for bed in params["bed"]:
                # subset gene bed, and then fix coords
                shell(
                    """bedtools intersect -header -f 1.0 -a {input.bed} -b <(printf "{bed}") \
        | bedtools sort -i - > {output.tmp} """
                )

                chrm, start, end = bed.split()
                start = int(start)
                name = f"{chrm}:{start}-{end}"
                for line in open(output["tmp"]).readlines():
                    t = line.strip().split()
                    t[0] = name
                    t[1] = int(t[1]) - start
                    t[2] = int(t[2]) - start
                    rtn += (11 * "{}\t" + "{}\n").format(*t)

            open(output["bed12"], "w+").write(rtn)


rule liftoff_genes:
    input:
        fa1=rules.get_rgns.output.query,
        fa2=rules.get_rgns.output.ref,
        rgns=rules.clean_gff.output.rgns,
        ref=get_ref,
        gff=rules.clean_gff.output.gff,
    output:
        bed12="Liftoff/{SM}.all.bed",
        fasta="Liftoff/tmp.combined.{SM}.fasta",
    params:
        bed=get_ref_bed,
    threads: 8
    resources:
        mem=4,
        hrs=24,
    run:
        if not GENES:
            shell("touch {output}")
        else:
            shell(
                "cat {input.fa1} {input.fa2} | seqtk seq -l 60 > {output.fasta}; sleep 5"
            )
            shell("samtools faidx {output.fasta}")

            shell(
                """snakemake -s {SDIR}/workflows/liftoff.smk \
        -j {threads} -p \
        {output.bed12} \
        --config \
            fasta=$(readlink -f {output.fasta}) \
            ref={input.ref} \
            gff=$(readlink -f {input.gff}) \
            sample={wildcards.SM} \
        --nolock || touch {output.bed12}"""
            )


#
# make the alignments and the miropeats
#
def get_score(wildcards):
    return int(str(wildcards.SCORE))


rule minimap2:
    input:
        ref=rules.get_rgns.output.ref,
        query=rules.get_rgns.output.query,
    output:
        paf=tempd("temp_minimiro/{SM}_{SCORE}_aln.paf"),
    params:
        score=get_score,
    threads: 16
    resources:
        mem=4,
        hrs=24,
    shell:
        """
        # YOU HAVE TO INCLUDE --cs FOR MINIMIRO TO WORK
        minimap2 -x asm20 --eqx -s {params.score} \
                   -p 0.01 -N 1000 --secondary={SECONDARY} \
                   --cs \
                   {input.ref} {input.query} 
            | rb break-paf --max-size {MAX_INDEL} \
            > {output.paf}
        """


def get_in_gene(wildcards):
    if GENES:
        out = [(rules.liftoff_genes.output.bed12).format(SM=wildcards.SM)]
    elif wildcards.SM in GENEBED:
        out = [GENEBED[wildcards.SM]]
    else:
        shell("touch tmp.fake.bed")
        out = ["tmp.fake.bed"]
    return out


rule minimiro:
    input:
        paf=rules.minimap2.output.paf,
        rmout=expand("temp_minimiro/{{SM}}_{SEQ}.fasta.out", SEQ=SEQS),
        dmout=expand("temp_minimiro/{{SM}}_{SEQ}.fasta.duplicons.extra", SEQ=SEQS),
        genes=get_in_gene,
    output:
        ps="temp_minimiro/{SM}_{SCORE}_aln.ps",
        pdf="miro_figures/{SM}_{SCORE}_aln.pdf",
    threads: 1
    resources:
        mem=16,
        hrs=24,
    shell:
        """
        {SDIR}/scripts/minimiro.py --paf {input.paf} \
          --rm {input.rmout} \
          --dm {input.dmout} \
          {simple_param} \
          --bed <(cat {input.genes} | cut -f 1-12 | grep -i '{GREP}' ) \
          --bestn 1000 \
          -o {output.ps} && \
        ps2pdf {output.ps} {output.pdf}
        """


def get_bam(wc):
    SM = str(wc.SM)
    return BAMS[SM]


rule coverage:
    input:
        bam=get_bam,
    output:
        png="miro_figures/{SM}_coverage.png",
    params:
        rgn=get_query_rgns,
    threads: 1
    resources:
        mem=8,
        hrs=24,
    shell:
        """
        NucPlot.py {input.bam} {output.png} --regions {params.rgn} --height 4 --width 16
        """
