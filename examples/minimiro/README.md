To run minimiro you first need to setup your config file called `minimiro.yaml`. 
You can find a simple example with some instructions about inputs in the `minimiro.yaml` 
included in this directory. Remember to index both your reference and query fasta files with `samtools faidx` before running.

Once you have that set up you can run minimiro using `snakemake`. 

To run the test case in this dir you would type:
```bash
snakemake -s ../../workflows/minimiro.smk --configfile minimiro.yaml --cores 60 
```
