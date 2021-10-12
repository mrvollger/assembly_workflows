[![DOI](https://zenodo.org/badge/288580731.svg)](https://zenodo.org/badge/latestdoi/288580731)
# assembly_workflows
Interconnected snakemake workflows for annotation and analysis of assemblies.

# Running minimiro
To run `minimiro` please see the example and directions in `examples/minimiro/`

# including a workflow
To include a workflow (e.g. mask.smk) just add the following to you snakemake near the top:
```
include: "/path/to/workflows/mask.smk"
```
All the workflows here expect two things to exist in the config before they are included. 
```
config["fasta"]="/path/to/fasta/or/assembly/something.fasta"
config["sample"]="SomeSampleName"
```
Note the path to the fasta must be a full path, becuase the relative path may change within the included workflow.

These options can be within the snakemake directly, or you can specify them via command line:
```
snakemake --config fasta=/path/to/fasta/or/assembly/something.fasta sample=SomeSampleName
```
if you want to do a cluster submission try something like:
```
snakemake --config fasta=/path/to/fasta/or/assembly/something.fasta sample=SomeSampleName \
    --drmaa " -l centos=7 -l h_rt=48:00:00 -l mfree={resources.mem}G -pe serial {threads} -V -cwd -S /bin/bash -w n" \
    --drmaa-log-dir Assembly_analysis/logs

```

A minimal example of how to include can be found in `workflows/run.smk`.
