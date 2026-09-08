#!/bin/bash
# Usage: bash 2_run_multiinter.sh data/filtered data/multiinter

in_dir=$1
out_dir=$2
mkdir -p "$out_dir"

declare -A groups
groups[AFR]="ESN GWD LWK MSL YRI ACB ASW"
groups[EUR]="CEU FIN GBR IBS TSI"
groups[AMR]="CLM PEL MXL PUR"
groups[EAS]="CDX CHB CHS JPT KHV"
groups[SAS]="BEB GIH ITU PJL STU"

for continent in "${!groups[@]}"; do
  samples=""
  for pop in ${groups[$continent]}; do
    samples+=" $in_dir/${pop}.all_0.05.bed"
  done

  bedtools multiinter -i $samples \
    | awk -v OFS='\t' '$4 >= 1 {print $1, $2, $3, $4}' \
    > "$out_dir/${continent}.multi.bed"
