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
### 6 samples
### using liftoff and minigraph applied to ALL samples
```

snakemake -s {path_to_assembly_workflows}/assembly_workflows/workflows/pull_regions.smk minigraph --config minigraph_samples=all liftoff_samples=all table=test_Master.tbl -p -j 32


```

