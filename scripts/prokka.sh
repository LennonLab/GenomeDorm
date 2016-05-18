#!/bin/bash
####Note#####
# This script runs on my personal mac. Prokka is not installed on dc2


mkdir -p "/Users/WRShoemaker/github/GenomeDorm/data/prokka/Bacillus_subtilis_subsp_subtilis_str_168/"

cd "/Users/WRShoemaker/github/GenomeDorm/data/fa/Bacillus_subtilis_subsp_subtilis_str_168/"

for i in *.fa
do
  iType="$(echo "$i" | cut -d "." -f1-1)"
  echo $iType
  OUT="/Users/WRShoemaker/github/GenomeDorm/data/prokka/Bacillus_subtilis_subsp_subtilis_str_168/${iType}"
  prokka --compliant --centre I --outdir $OUT --locustag G --prefix G-Chr1 $i --force
done
