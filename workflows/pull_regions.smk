import os
import pandas as pd
import networkx as nx
import pysam

SDIR=os.path.realpath(os.path.dirname(srcdir("env.cfg"))+"/..")
shell.prefix(f"source {SDIR}/env.cfg ; set -eo pipefail; ")

if("table" in config.keys()):
    table = config["table"]
else:
    table = "Master.tbl"

if os.path.exists("Master_SD_freeze.tbl"):
  table = "Master_SD_freeze.tbl"
df = pd.read_csv(table, sep="\t", comment="#")
CHM13D="/net/eichler/vol26/projects/chm13_t2t/nobackups"
REF=os.path.abspath(f"{CHM13D}/assemblies/chm13_v1.1_plus38Y.fasta")
FAI=os.path.abspath(f"{REF}.fai")
GFF=os.path.abspath(f"{CHM13D}/Assembly_analysis/Liftoff/chm13_v1.1_plus38Y.all.gff3")
FLANK=50000
RMFLANK=FLANK-2000
RMFLANK=0
#df = df[(df.hap != "pri") & (df.hap != "alt")]
#df = df[(df.hap != "alt")]
df.set_index(["sample","hap"], inplace=True)
regions = pd.read_csv("regions.bed", header=None,sep="\t", names=["chr", "start","end","name"], comment="#")
BED = os.path.abspath("regions.bed")
regions.set_index("name", inplace=True)
sms, haps = zip(*df.index.values)
rgns = list(regions.index.values)
haps=set(haps)

workdir: "pull_sd_regions"

def get_bed(wc):
  x=df.loc[(wc.sm, wc.h)].all_bed
  return(x)

def get_fasta(wc):
  return(df.loc[(wc.sm,wc.h)].fasta)


def get_all(wc):
  f="pulled_fastas/{r}.{sm}.{h}.fasta"
  o = []
  for sm, hap in df.index.values:
    for rgn in rgns:
      o.append(f.format(sm=sm,h=hap,r=rgn))
  return(o)

def get_rgn_fasta(wc):
  f="pulled_fastas/{r}.{sm}.{h}.fasta"
  o = []
  for sm, hap in df.index.values:
      rgn=wc.r
      o.append(f.format(sm=sm,h=hap,r=rgn))
  return(o)

def get_liftoff_choice(wc):
    ''' give user option of to generate liftoff results for all haplotypes vs. only simple.fasta result.'''
    print("config_keys" + str( config.keys()) )
    if("liftoff_samples" in config.keys()):
        if(config['liftoff_samples'] == 'all'):
            return( rules.simple_fasta.output.allfasta )
    return(rules.simple_fasta.output.fasta)

wildcard_constraints:
  sm="|".join(sms),
  h="|".join(haps),
  r="|".join(rgns)

rule all:
  input:
    f=get_all,
    simple = expand("combined/{r}.simple.fasta", r=rgns),
    yaml = "minimiro.yaml"
#run:
#     shell("mkdir -p combined")
#     for rgn in rgns:
#       shell("cat pulled_fastas/{rgn}*fasta > combined/{rgn}.fasta")
#       shell("samtools faidx combined/{rgn}.fasta")

rule get_gene:
  input:
    bed=get_bed,
    fasta=get_fasta,
  output:
    bed="temp/{sm}.{h}.bed",
    exon="temp/{sm}.{h}.exons",
    fasta="temp/{sm}.{h}.fasta",
    tbl="temp/{sm}.{h}.tbl",
  shell:"""
module load bedtools/2.29.2
grep -P "LPA\s+" {input.bed} | head -n 1 > {output.bed}
bedtools getfasta -s -split -fi {input.fasta} -bed <(cut -f 1-12 {output.bed}) > {output.exon}

bedtools getfasta -s -name -fi {input.fasta} -bed <(cut -f 1-12 {output.bed} |  bedtools slop -b 400000 -i - -g {input.fasta}.fai) > {output.fasta}

cds=/net/eichler/vol26/home/mvollger/targeted_hifi_asm/nobackups/all_ccs_data_aln/Targeted_trio/old_scripts/gene_map/CDS

minimap2 -ax splice --eqx -s 500 -p .3 -f 0 -N 1000 -r 150000 {output.fasta} {output.exon} | \
           samtools view -b - | \
           bedtools bamtobed -split -i - | \
           sed "s/$/\t{wildcards.sm}\t{wildcards.h}/" \
           > {output.tbl}
"""

#print(regions)
def get_rgn(wc):
  x = regions.loc[wc.r]
  f1 = "{}:{}-{}".format(x.chr, x.start - FLANK, x.start)
  f2 = "{}:{}-{}".format(x.chr, x.end, x.end + FLANK)
  return(f1 + " " + f2)

def get_whole_rgn(wc):
  x = regions.loc[wc.r]
  rgn = "{}:{}-{}".format(x.chr, x.start - FLANK, x.end + FLANK)
  return(rgn)

MM_OPTS=f" -r 50000 -x asm20 -z 10000 -s {int(FLANK/2)} "
rule index:
  input:
    fasta=get_fasta,
  output:
    mmi="temp/index/{sm}.{h}.mmi",
  threads: 8
  shell:"minimap2 {MM_OPTS} -t {threads} -d {output.mmi} {input.fasta}"

rule get_rgn_query:
  input:
    ref = REF,
    genome = FAI,
  output:
    query = temp("temp/{sm}.{h}.query.flank.fasta"),
  threads: 4
  shell:"""
bedtools getfasta -name+ \
  -bed <(bedtools flank -b {FLANK} -g {input.genome} -i {BED}) \
  -fi {input.ref} \
  > {output.query}
"""

rule get_rgn_paf:
  input:
    mmi=rules.index.output.mmi,
    query="temp/{sm}.{h}.query.flank.fasta",
  output:
    all_paf=temp("temp/{sm}.{h}.paf"),
  threads: 4
  shell:"""
minimap2 {MM_OPTS} -2 --secondary=no --eqx -Y -t {threads} \
  {input.mmi} {input.query} > {output.all_paf}
"""

rule get_rgn:
  input:
    all_paf=ancient("temp/{sm}.{h}.paf"),
  output:
    paf="temp/{sm}.{h}.{r}.paf",
  threads: 1
  run:
    prefix = wildcards.r + "::"
    out = open(output.paf, "w+")
    for line in open(input.all_paf):    
      t = line.strip().split() 
      # check for gene of interest
      if t[0].startswith(prefix):
        # remove prefix
        t[0] = t[0].replace(prefix, "")
        line = "\t".join(t) + "\n"
        out.write(line)
    out.close()

def get_dist(wc):
  x = regions.loc[wc.r]
  dist = 5*(x.end - x.start)
  return(dist)

rule pull_fasta:
  input:
    # adding ancient because we don't want to rerun the pipeline when one new region is added
    paf = ancient(rules.get_rgn.output.paf), 
    fasta=get_fasta,
    ref=REF,
  output:
    bed="temp/{sm}.{h}.{r}.bed",
    fasta="temp/{sm}.{h}.{r}.fasta",
  params:
    dist=get_dist,
    rgn=get_whole_rgn,
  shell:"""
awk  '{{print $6"\t"$8"\t"$9"\t"$1"_"$6"\t"$3"\t"$4}}' {input.paf} | \
      bedtools sort -i - | bedtools merge -d {params.dist} -i - | \
      bedtools slop -i - -g {input.fasta}.fai -b -{RMFLANK} \
      > {output.bed}

# make fasta be in the same orientation as the reference 
minimap2 -t {threads} -m 10000 -ax asm20 -r 200000 --eqx -Y \
  <(samtools faidx {input.ref} {params.rgn})\
  <(bedtools getfasta -fi {input.fasta} -bed {output.bed}) | \
  samtools view -F 2308 | \
  awk '{{OFS="\\t"; print ">"$1"\\n"$10}}' - | \
  seqtk seq -l 60 > {output.fasta} 

"""


rule good_fasta:
  input:
      fasta = rules.pull_fasta.output.fasta
  output:
      good = "pulled_fastas/{r}.{sm}.{h}.fasta"
  run:
    import pysam
    o = open(output.good, "w+")
    if os.stat(input.fasta).st_size > 0:
      shell("samtools faidx {input.fasta}")
      f = pysam.FastaFile(input.fasta)
      if(len(f.references) == 1):
        for idx, name in enumerate(f.references):
          seq = f.fetch(name)
          o.write(f">{wildcards.r}__{wildcards.sm}__{wildcards.h}__ID.{idx+1}\n")
          o.write(seq+"\n")
      else:
        o.write("")
    o.close()

def get_minscore(wc):
  x = regions.loc[wc.r]
  total = x.end - x.start 
  return(max(FLANK-RMFLANK, int(total/3)))

rule region_fasta:
  input:
      fasta = get_rgn_fasta,
  output:
      tbl = "combined/{r}.tbl",
      tmp = "combined/{r}.tmp",
      bam = "combined/{r}.unimap.bam",
      utbl = "combined/{r}.unimap.tbl",
  params:
      minscore = get_minscore,
  threads: 80
  shell:"""
cat {input.fasta} > {output.tmp}
samtools faidx {output.tmp}
minimap2 -r 50000 -ax asm20 -s {params.minscore} --eqx -Y -t {threads} \
          {output.tmp} {output.tmp} \
         | samtools view -F 4 -b - | samtools sort -@ {threads} - \
         | ~mvollger/projects/utility/samIdentity.py --header /dev/stdin > {output.tbl}

module load unimap/0.1
unimap -ax asm20 -s {params.minscore} --eqx -Y -t {threads} \
          {output.tmp} {output.tmp} \
          | samtools view -F 4 -b - \
          | samtools sort \
          > {output.bam}

~mvollger/projects/utility/samIdentity.py \
         --header {output.bam} \
         > {output.utbl}
"""

def sort_ctgs(ctgs, extra=True):
  first = []
  human = []
  other = []
  for c in ctgs:
    if("CHM13" in c.upper()):
      first.insert(0, c)
    elif("GRCH38" in c.upper() or "hg38" in c.lower()):
      first.insert(1, c)
    elif(c[:2].upper() in ["NA", "HG"] or "CHM1" in c.upper()):
      human.append(c)
    else:
      other.append(c)
  if(extra):
    human = sorted(human, key=lambda x: ( -int(x.split("__")[-1]), x.split("__")[0] ))
  else:
    human.sort()
  other.sort()
  return(first+human+other)


rule simple_fasta:
  input:
      tbl = "combined/{r}.tbl",
      tmp = "combined/{r}.tmp",
  output:      
      allfasta = "combined/{r}.all.fasta",
      allfastafai = "combined/{r}.all.fasta.fai",
      fasta = "combined/{r}.simple.fasta",
      fastafai = "combined/{r}.simple.fasta.fai",
      yaml = "combined/{r}.simple.yaml",
  run:
    # read in the fasta data
    fasta=pysam.FastaFile(input.tmp)
    print(fasta.references)
    
    # make a table of shared connections 
    df = pd.read_csv(input.tbl, sep="\t")
    e = df[(df.perID_by_all > 99) 
            & ((df.query_end-df.query_start)/df.query_length > 0.9) ]
#& ~df.query_name.str.contains("CHM13|GRCh38")
#          & ~df.reference_name.str.contains("CHM13|GRCh38")]
   
    # make a graph of the connections
    g = nx.Graph()
    g.add_nodes_from(fasta.references)
    #g.add_nodes_from(set(list(df.query_name) + list(df.reference_name)))
    for i, x, y in e[["query_name", "reference_name"]].itertuples():
      g.add_edge(x,y)
      g.add_edge(y,x)
   
    # write the connection groups to file
    out = open(output.fasta, "w+")
    ctgs=[]
    seen = set()
    for x in nx.connected_components(g):
      sorted_names = sort_ctgs(x, extra=False)
      clean_names = []
      # get all the names going into this set  
      for n in sorted_names:
        m = re.match(".+__(.+)__(.+)__ID.\d+", n)
        clean_names.append( ".".join(m.groups()) )
      cmt = ":".join(clean_names)
      
      # make sure we process each contig only once
      s = False 
      for n in sorted_names:
        if n in seen:
          s = True
        seen.add(n)
      if s: break 

      # write the sequence to file 
      name = clean_names[0] + "__" + str(len(clean_names)) + ""
      sys.stderr.write(name + "\t" + cmt + "\n")
      ctg=sorted_names[0]
      seq = fasta.fetch(ctg)
      ctgs.append(name)
      out.write(f">{name}\t{cmt}\n{seq}\n")
    
    out.close() 
    shell("samtools faidx {output.fasta}")
    
    row = regions.loc[wildcards.r]
    r = "{}:{}-{}".format(row.chr, row.start-FLANK, row.end+FLANK)
    ctgs = sort_ctgs(ctgs)
    q="\n    - ".join(ctgs)
    z=f"""
{wildcards.r}:
  ref: {REF}
  gff: {GFF}
  regions:
    - {r}
  query: {os.path.abspath(output.fasta)}
  queryregions:
    - {q}

""" 
    open(output.yaml,"w+").write(z)
    #
    # write all contigs to output
    #
    allo = open(output.allfasta, "w+")
    for name in sort_ctgs(set(fasta.references), extra=False):
      seq = fasta.fetch(name)
      allo.write(f">{name}\n{seq}\n")
    allo.close() 
    shell("samtools faidx {output.allfasta}")
       
rule yaml:
    input:
      yaml = expand("combined/{r}.simple.yaml", r=rgns),
    output:
      yaml = "minimiro.yaml"
    shell: """
if [ -f {output.yaml} ]; then
  mv {output.yaml} {output.yaml}.bck
fi
printf 'scores:\\n  - 5000\\nsimple: True\\n\\n' > {output.yaml}
cat {input.yaml} >> {output.yaml}
"""

#
#
# STARTING THE MINIGRAPH SECTION
#
#
rule get_duplicons:
  input:
    fasta = rules.simple_fasta.output.fasta,
  output:
    bed = "Masked/{r}_dupmasker_colors.bed"
  threads: 12
  shell:"""
snakemake -s {SDIR}/workflows/mask.smk \
  -j {threads} -p -k \
  {output.bed} \
  --config \
      threads=4 \
      fasta=$(readlink -f {input.fasta}) \
      sample={wildcards.r} \
  --nolock 
"""

rule get_genes:
    input:
      fasta = get_liftoff_choice, #rules.simple_fasta.output.fasta
      bed="temp/GRCh38chrOnly.pri.{r}.bed",
      #ref = REF,
      #gff = GFF,
    output:
      #subset = "Liftoff/{r}.subset.bed",
      bed12 = "Liftoff/{r}.orf_only.bed",
      bedall = "Liftoff/{r}.all.bed",
    threads: 12
    run:

      # liftoff does weird things with reruns so I am going to clean the dir first
      shell("rm -rf Liftoff/temp.{wildcards.r}/")
      shell("rm -f {input.fasta}.mmi")
      shell("rm -f Liftoff/{wildcards.r}.all.*")
      shell("rm -f Liftoff/{wildcards.r}.orf_only.*")
      shell("rm -f Liftoff/{wildcards.r}.gff3")

      shell("""snakemake -s {SDIR}/workflows/liftoff.smk \
                -j {threads} -p \
                {output.bed12} {output.bedall} \
                --config \
                    fasta=$(readlink -f {input.fasta}) \
                    regions=$(readlink -f {input.bed}) \
                    sample={wildcards.r} \
                --nolock  """)#|| touch {output.bed12} {output.bedall}""")



      if False: # old pipeline that uses CHM13 as the reference for genes
        x = regions.loc[wildcards.r]
        open(output.subset, "w+").write(f"{x.chr}\t{x.start}\t{x.end}\n")
        shell("""snakemake -s {SDIR}/workflows/liftoff.smk \
                -j {threads} -p \
                {output.bed12} {output.bedall} \
                --config \
                    fasta=$(readlink -f {input.fasta}) \
                    ref={input.ref} \
                    gff=$(readlink -f {input.gff}) \
                    regions=$(readlink -f {output.subset}) \
                    sample={wildcards.r} \
                --nolock  || touch {output.bed12} {output.bedall}""")



rule mg_make_gfa:
  input:
    fasta = rules.simple_fasta.output.fasta,
  output:
    gfa = "Minigraph/{r}.gfa",
    fastas = directory("Minigraph/temp.{r}/"),
    fasta = "Minigraph/{r}.all.fasta",
  threads: 16
  run:
      shell("samtools faidx {input.fasta}")
      pairs = { line.split()[0]:int(line.split()[1]) for line in open(input.fasta + ".fai") }
      names = sorted(list(pairs), key = lambda x: pairs[x]) 
      ordered = [None, None]
      for name in names:
        print(name)
        path=f"Minigraph/temp.{wildcards.r}/{name}.fasta"
        if name.startswith("GRCh38chrOnly"):
          ordered.insert(1, path)
        elif name.startswith("CHM13"):
          ordered.insert(0, path)
        elif name.startswith("CHM1"):
          ordered.append(path)
        elif ".pri" in name or ".alt" in name:
          continue
        else:
          ordered.append(path)
        shell(f"samtools faidx {input.fasta} {name} > {path}") 
        ordered = [i for i in ordered if i != None]
        shell("""
cat {ordered} | seqtk seq -l 80 > {output.fasta}
minigraph -xggs -L 5000 -r 100000 -t {threads} {ordered} > {output.gfa} """)
    
rule mg_map:
  input:
    gfa = "Minigraph/{r}.gfa",
    fasta = "Minigraph/{r}.all.fasta",
  output:
    gaf = "Minigraph/{r}.gaf",
    bed = "Minigraph/{r}.bed",
    tmp = temp("Minigraph/{r}.bed.tmp.fa"),
  threads: 16
  run:
    shell("minigraph -x asm -N 0 -t {threads} {input.gfa} {input.fasta} > {output.gaf}")
    shell("samtools faidx {input.fasta} && > {output.bed}")
    for line in open(input.fasta + ".fai"):
      name = line.split()[0]
      shell("""
          samtools faidx {input.fasta} {name} > {output.bed}.tmp.fa 
          minigraph -x asm -N 0 -t {threads} \
              --call {input.gfa} \
              {output.bed}.tmp.fa \
              >> {output.bed} || echo failed minigraph call on {name} """)


rule mg_parse:
  input:
    gaf = "Minigraph/{r}.gaf",
    gfa = "Minigraph/{r}.gfa",
  output:
    tbl = "Minigraph/{r}.tbl",
    csv = "Minigraph/{r}.csv",
  threads: 1
  shell:"""
{SDIR}/scripts/GAF_parsing.py \
    --ref $(cut -f 1 {input.gaf} | grep "CHM13.pri" | head -n 1 ) \
    {input.gaf} > {output.tbl}
{SDIR}/scripts/make_csv_from_gfa.sh {input.gfa} > {output.csv}
"""


rule sub_table:
  input:
    fasta = rules.simple_fasta.output.fasta,
    tbl = rules.mg_parse.output.tbl, 
  output:
    tbl = "Tables/{r}.tbl",
  threads:1
  run:
    fai = pd.read_csv(input.fasta+".fai",
        sep="\t", 
        names=["contig", "length", "x","y","z"])[["contig", "length"]]
    fai = fai.set_index('contig')
    fai["region"] = wildcards.r 
    fai["haplotypes"] = ""
    for line in open(input.fasta).readlines():
      if line[0] is not ">":
        continue 
      ts = line[1:].strip().split()
      sup_hap, all_haps = ts[0], ts[1]
      fai.at[sup_hap, "haplotypes"] = all_haps
   

    mg = pd.read_csv(input.tbl, sep="\t")
    df = mg.merge(fai, left_on="q", right_index=True)
    print(df)
    df.to_csv(output.tbl, index=False, sep = "\t")


rule table:
  input:
    fastas = expand(rules.simple_fasta.output.allfasta, r=rgns),
    tbls = expand(rules.sub_table.output.tbl, r=rgns),
  output:
    tbl="results.tbl",
    svtbl="sv.results.tbl",
  run:
    dfs=[]
    for fasta in input.fastas:
      df = pd.read_csv(fasta+".fai", sep="\t", names=["contig", "length", "x","y","z"])
      df[['Region','sm','hap','drop']] = df.contig.str.split("__", expand=True)
      df[['sm','Species']] = df.sm.str.split("_", expand=True)
      df.Species.replace({None:"Human"}, inplace=True)
      df["fasta"] = os.path.abspath(fasta)
      dfs.append(df)
    df = pd.concat(dfs, ignore_index=True)   
    df.drop(['x','y','z','drop'], axis=1, inplace=True)
    print(df)
    df.to_csv(output.tbl, sep="\t", index=False)

    # the other table
    shell("head -n 1 {input.tbls[0]} > {output.svtbl}" )
    shell("tail -n +2 -q {input.tbls} >> {output.svtbl}" )
    
rule minigraph: 
  input:
    tbl = expand(rules.mg_parse.output.tbl, r=rgns),
    bed = expand(rules.mg_map.output.bed, r=rgns),
    genes = expand(rules.get_genes.output.bed12, r=rgns),
    duplicons = expand(rules.get_duplicons.output.bed, r=rgns),
    fulltbl=rules.table.output.tbl,
    svtbl=rules.table.output.svtbl,


