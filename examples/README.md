## SEDEF
```
snakemake --config fasta={target.fasta} sample={sample.name} -s ../workflows/liftoff.smk -j 120 sedef
```

## Liftoff
### single sample
```
snakemake --config fasta={target.fasta} sample={sample.name} -s ../workflows/liftoff.smk -j 120 liftoff
```
### multi sample
```
snakemake --configfile lift.yaml -s ../workflows/liftoff.smk -j 130 -p liftoff
```
