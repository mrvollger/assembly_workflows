import os 
import sys
import re
import re
import pysam 
import pandas as pd
from datetime import date

today = date.today()
DATE =  today.strftime("%Y/%m/%d")

SDIR=os.path.dirname(workflow.snakefile)
CWD=os.getcwd()
shell.prefix(f"source {SDIR}/env.cfg ; set -eo pipefail; ")

workdir: "Assembly_analysis"

include: "liftoff.smk"
include: "mask.smk"
# sedef also includes mask.smk so it can run standalone need be
include: "sedef.smk"
include: "align.smk"

rule all:
    input:
        rules.liftoff.input,
        rules.sedef.input,
        rules.DupMaskerBed.output.bed,
        rules.DupMaskerSummary.output,
#rules.align.input,


