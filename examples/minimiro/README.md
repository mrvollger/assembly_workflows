# Running minimiro
To run minimiro you first need to setup your config file called `minimiro.yaml`. 
You can find a simple example with some instructions about inputs in the `minimiro.yaml` 
included in this directory. Remember to index both your reference and query fasta files with `samtools faidx` before running.

Once you have that set up you can run minimiro using `snakemake`. 

To run the test case in this directory you would type:
```bash
snakemake -s ../../workflows/minimiro.smk --configfile minimiro.yaml --cores 60 
```
This `snakemake` is not set up for cluster distribution, please use on a single node like `ocelot`, `lynx`, or `liger`.

## Running into errors

Try these things to before making an issue: 
  - Read the `minimiro.yaml` example configuration carefully.
  - Make sure this test case runs.
  - Make a new and unique prefix name and try rerunning.
  - Check to make sure your ref and query are indexed.
  - Remove all minimiro output and try again.

If you do need to contact me with an issue please be ready to share a minimal
test case as well as a log file with the output of your run. 
You can capture the output of your run like this for example:
```
snakemake ... &> minimiro.log
```

