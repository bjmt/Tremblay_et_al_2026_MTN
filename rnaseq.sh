#!/usr/bin/env bash

# Usage: bash rnaseq.sh rna_bud|rna_seedling align|counts|de

set -euo pipefail

root=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
assay=${1:-}
stage=${2:-}

case "$assay" in
  rna_bud)
    read_length=150
    threshold=16
    ;;
  rna_seedling)
    read_length=126
    threshold=14
    ;;
  *)
    echo 'Assay must be rna_bud or rna_seedling' >&2
    exit 2
    ;;
esac

case "$stage" in
  align|counts|de)
    ;;
  *)
    echo 'Stage must be align, counts, or de' >&2
    exit 2
    ;;
esac

threads=${THREADS:-4}
data=${MTN_DATA:-$root/../mtn-data}
output=${MTN_OUTPUT:-$root/../mtn-output}

data=$(cd "$data" && pwd)
mkdir -p "$output"
output=$(cd "$output" && pwd)

reference=${REFERENCE:-$data/reference/TAIR10.fa}
annotation=${ANNOTATION:-$data/reference/Araport11_GFF3_genes_transposons.gff3}
work=$output/preprocessing/$assay

mkdir -p "$work"
cd "$work"

need() {
  if ! command -v "$1" >/dev/null; then
    echo "Missing executable: $1" >&2
    exit 1
  fi
}

if [[ "$stage" == align ]]; then

  for tool in cutadapt STAR stringtie fastqc; do
    need "$tool"
  done

  if [[ ! -s "$reference" || ! -s "$annotation" ]]; then
    echo 'Supply TAIR10 FASTA and Araport11 GFF3; see README.md' >&2
    exit 1
  fi

  mkdir -p trimmed aligned stringtie qc star_index logs

  if [[ ! -s star_index/Genome ]]; then

    STAR \
      --runMode genomeGenerate \
      --runThreadN "$threads" \
      --genomeDir star_index \
      --genomeFastaFiles "$reference" \
      --sjdbGTFfile "$annotation" \
      --sjdbGTFtagExonParentTranscript Parent \
      --sjdbGTFtagExonParentGene ID \
      --sjdbOverhang "$((read_length - 1))" \
      --genomeSAindexNbases 12

  fi

  if [[ "$assay" == rna_bud ]]; then
    prefixes=(gm gs mm ms wm ws)
    replicates=4
  else
    prefixes=(5m 5s mm ms wm ws)
    replicates=3
  fi

  for prefix in "${prefixes[@]}"; do
    for ((replicate=1; replicate<=replicates; replicate++)); do

      sample=${prefix}${replicate}
      r1=$data/raw/$assay/${sample}_R1.fastq.gz
      r2=$data/raw/$assay/${sample}_R2.fastq.gz

      if [[ ! -s "$r1" || ! -s "$r2" ]]; then
        echo "Missing FASTQs for $sample" >&2
        exit 1
      fi

      fastqc \
        -t "$threads" \
        -o qc \
        "$r1" "$r2" \
        > "logs/$sample.fastqc.log" 2>&1

      cutadapt \
        -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA \
        -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
        -j "$threads" -m 15 \
        -o "trimmed/${sample}_1.fastq.gz" \
        -p "trimmed/${sample}_2.fastq.gz" \
        "$r1" "$r2" \
        > "logs/$sample.cutadapt.log" 2>&1

      STAR \
        --genomeDir star_index \
        --runThreadN "$threads" \
        --readFilesCommand gunzip -c \
        --readFilesIn "trimmed/${sample}_1.fastq.gz" "trimmed/${sample}_2.fastq.gz" \
        --outFileNamePrefix "aligned/$sample." \
        --outSAMstrandField intronMotif \
        --outMultimapperOrder Random \
        --outSAMmultNmax 1 \
        --outSAMtype BAM SortedByCoordinate \
        --outFilterMultimapNmax 10000 \
        --limitOutSAMoneReadBytes 10000000 \
        > "logs/$sample.STAR.log" 2>&1

      stringtie "aligned/$sample.Aligned.sortedByCoord.out.bam" \
        -G "$annotation" \
        -A "stringtie/$sample.gene_abund.txt" \
        -e -p "$threads" \
        -o "stringtie/$sample.gtf"

    done
  done

elif [[ "$stage" == counts ]]; then

  # prepDE.py3 is distributed with StringTie.
  prepde=${PREPDE:-prepDE.py3}
  need "$prepde"

  find stringtie -name '*.gtf' -type f | LC_ALL=C sort | \
    awk '{
      sample=$0
      sub(/^.*\//, "", sample)
      sub(/[.]gtf$/, "", sample)
      print sample, $0
    }' > gtfs.txt

  if [[ ! -s gtfs.txt ]]; then
    echo 'No StringTie GTF files found; run align first' >&2
    exit 1
  fi

  "$prepde" -l "$read_length" -i gtfs.txt

elif [[ "$stage" == de ]]; then

  need getDiffExpression.pl

  if [[ ! -s gene_count_matrix.csv ]]; then
    echo 'Run counts first' >&2
    exit 1
  fi

  # Exclude rRNA genes and apply the total-count threshold.
  # The seedling threshold uses the first 16 samples.
  awk -F ',' -v limit="$threshold" -v assay="$assay" '
    BEGIN {
      OFS="\t"
      split("AT2G01010 AT2G01020 AT3G41768 AT3G41979 ATCG00920 ATCG00950 ATCG00960 ATCG00970 ATCG01160 ATCG01170 ATCG01180 ATCG01210", ids, " ")
      for (i in ids) rrna[ids[i]]=1
    }
    FNR==1 {
      gsub(/\r/, "")
      gsub(/,/, "\t")
      print
      next
    }
    {
      gsub(/\r/, "")
      total=0
      last=(assay=="rna_seedling" ? 17 : NF)
      for (i=2; i<=last; i++) total+=$i
      if (total>limit && !($1 in rrna)) {
        gsub(/,/, "\t")
        print
      }
    }
  ' gene_count_matrix.csv > gene_count_matrix.noRRNA.tsv

  # Assign groups in count-matrix column order.
  groups=()
  IFS=$'\t' read -r -a header < gene_count_matrix.noRRNA.tsv

  for sample in "${header[@]:1}"; do

    prefix=${sample%[0-9]}

    if [[ "$assay" == rna_seedling ]]; then

      case "$prefix" in
        5m|5s|mm|ms|wm|ws)
          group=$prefix
          ;;
        *)
          echo "Unknown sample: $sample" >&2
          exit 1
          ;;
      esac

    else

      case "$prefix" in
        gm)
          group=g3_ms
          ;;
        gs)
          group=g3_spd
          ;;
        mm)
          group=mtn_ms
          ;;
        ms)
          group=mtn_spd
          ;;
        wm)
          group=wt_ms
          ;;
        ws)
          group=wt_spd
          ;;
        *)
          echo "Unknown sample: $sample" >&2
          exit 1
          ;;
      esac

    fi

    groups+=("$group")

  done

  getDiffExpression.pl gene_count_matrix.noRRNA.tsv "${groups[@]}" \
    -basic -AvsA -fdr 0.1 -DESeq2 \
    > gene_count_diff.tsv \
    2> gene_count_diff.log.txt

  Rscript --vanilla "$root/annotate_rnaseq.R" "$assay" "$work" "$annotation"

fi
