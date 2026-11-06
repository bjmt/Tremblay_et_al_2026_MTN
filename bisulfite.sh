#!/usr/bin/env bash

# Usage: bash bisulfite.sh align|extract|profiles
# Run dmr_calling.R after extract to call differentially methylated regions.

set -euo pipefail

root=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
stage=${1:-}

case "$stage" in
  align|extract|profiles)
    ;;
  *)
    echo 'Stage must be align, extract, or profiles' >&2
    exit 2
    ;;
esac

threads=${THREADS:-2}
data=${MTN_DATA:-$root/../mtn-data}
output=${MTN_OUTPUT:-$root/../mtn-output}

data=$(cd "$data" && pwd)
mkdir -p "$output"
output=$(cd "$output" && pwd)

reference=${REFERENCE:-$data/reference/TAIR10.fa}
annotation=${ANNOTATION:-$data/reference/Araport11_GFF3_genes_transposons.gff3}
work=$output/preprocessing/bs_bud

mkdir -p "$work"
cd "$work"

need() {
  if ! command -v "$1" >/dev/null; then
    echo "Missing executable: $1" >&2
    exit 1
  fi
}

if [[ "$stage" == align ]]; then

  for tool in trim_galore bismark bismark_genome_preparation deduplicate_bismark; do
    need "$tool"
  done

  if [[ ! -s "$reference" ]]; then
    echo "Missing reference: $reference" >&2
    exit 1
  fi

  mkdir -p genome trimmed bams
  ln -sf "$reference" genome/TAIR10.fa

  if [[ ! -d genome/Bisulfite_Genome ]]; then
    bismark_genome_preparation genome/
  fi

  for prefix in gm gs mm ms wm ws; do
    for replicate in 1 2 3; do

      sample=${prefix}${replicate}
      r1=$data/raw/bs_bud/${sample}_R1.fastq.gz
      r2=$data/raw/bs_bud/${sample}_R2.fastq.gz

      if [[ ! -s "$r1" || ! -s "$r2" ]]; then
        echo "Missing FASTQs for $sample" >&2
        exit 1
      fi

      trim_galore \
        --illumina --gzip --paired \
        --cores "$threads" \
        -o trimmed \
        "$r1" "$r2"

      bismark \
        --parallel "$threads" \
        -o bams \
        genome/ \
        -1 "trimmed/${sample}_R1_val_1.fq.gz" \
        -2 "trimmed/${sample}_R2_val_2.fq.gz"

      deduplicate_bismark \
        --paired --bam \
        --output_dir bams \
        "bams/${sample}_R1_val_1_bismark_bt2_pe.bam"

    done
  done

elif [[ "$stage" == extract ]]; then

  need bismark_methylation_extractor
  need makeTagDirectory

  mkdir -p tagDirs methcalls

  for prefix in gm gs mm ms wm ws; do
    for replicate in 1 2 3; do

      sample=${prefix}${replicate}
      bam=bams/${sample}_R1_val_1_bismark_bt2_pe.deduplicated.bam

      if [[ ! -s "$bam" ]]; then
        echo "Missing $bam; run align first" >&2
        exit 1
      fi

      mkdir -p "methcalls/$sample"

      bismark_methylation_extractor \
        --gzip --comprehensive --bedGraph \
        --ignore 2 --ignore_r2 2 \
        --cytosine_report --CX_context --ample_memory \
        --genome_folder genome \
        --output "methcalls/$sample" \
        "$bam"

      reports=("methcalls/$sample/"*CX_report.txt.gz)

      if [[ ${#reports[@]} != 1 || ! -s "${reports[0]}" ]]; then
        echo "Expected one CX report for $sample" >&2
        exit 1
      fi

      for context in CG CHG CHH; do

        makeTagDirectory \
          "tagDirs/${sample}_${context}-tagDir" \
          "${reports[0]}" \
          -format bismark \
          -minCounts 4 \
          -mCcontext "$context"

      done
    done
  done

elif [[ "$stage" == profiles ]]; then

  need gff2pos.pl
  need makeMetaGeneProfile.pl

  if [[ ! -s "$annotation" ]]; then
    echo "Missing annotation: $annotation" >&2
    exit 1
  fi

  mkdir -p profiles
  tagdirs=()

  # Use individual replicates, excluding merged tag directories.
  for group in gm gs mm ms wm ws; do
    for replicate in 1 2 3; do
      for context in CG CHG CHH; do

        directory=tagDirs/${group}${replicate}_${context}-tagDir

        if [[ ! -d "$directory" ]]; then
          echo "Missing $directory; run extract first" >&2
          exit 1
        fi

        tagdirs+=("$directory")

      done
    done
  done

  for feature in gene transposable_element; do

    awk -F '\t' -v feature="$feature" '
      $0 !~ /^#/ && $3 == feature
    ' "$annotation" > "profiles/$feature.gff3"

    gff2pos.pl "profiles/$feature.gff3" > "profiles/$feature.pos"

    makeMetaGeneProfile.pl "profiles/$feature.pos" tair10 \
      -d "${tagdirs[@]}" \
      -len 1 -strand both -ratio \
      -min 0 -gbin 200 -size 2000 -bin 100 \
      > "profiles/$feature.txt"

  done

fi
