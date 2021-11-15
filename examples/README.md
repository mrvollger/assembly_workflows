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

### Pull Regions
## 6 samples
#  using liftoff command to apply liftOff to ALL samples
```
snakemake -s ../workflows/pull_regions.smk --config liftoff_samples=all table=../test_data/test_Master.tbl -j 32

```

