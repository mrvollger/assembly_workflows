#!/usr/bin/env bash
set -euo pipefail

printf "Name,Label,Offset\n" 
grep "^S" $1 \
  | sed -r 's/^S\t(\S+).*(SN:Z:)(\S+).*(SO:i:)(\S+).*$/\1,\3,\5/' 

