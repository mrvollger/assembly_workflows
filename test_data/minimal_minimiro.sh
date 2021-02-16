#!/usr/bin/env bash
set -euo pipefail

mkdir -p temp

# this will show primary alignments,
# but if you want to see all homology like calissic miropeats
# change your -N and -p as desiered
### --cs MUST BE INCLUDED IN MM2 CMD
if [ ! -f temp/test.paf ]; then
  minimap2 -x asm20 -t 2 \
    --cs \
    -s 1000 --secondary=no \
    ./CHM13.pri.NOTCH2NL.fasta.gz \
    ./Clint_PTR.pri.NOTCH2NL.fasta.gz \
    > temp/test.paf 
fi

../scripts/minimiro.py --paf temp/test.paf -o temp/test.ps 


