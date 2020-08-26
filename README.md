# assembly_workflows
Interconnected snakemake workflows for annotation and analysis of assemblies.

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
A minimal example of how to include can be found in `workflows/run.smk`.
