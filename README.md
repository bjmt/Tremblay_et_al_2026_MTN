# MTN paper analysis

Shell and R scripts for [Methylthioadenosine-driven rewiring of sulfur assimilation lowers the methyl index and fertility in Arabidopsis](https://doi.org/10.1093/jxb/erag135).

## Data and dependencies

Keep this folder as code only. By default the scripts read `../mtn-data/` and write `../mtn-output/`. Set `MTN_DATA` and `MTN_OUTPUT` to use other locations.

Download the supplementary workbook from the paper as `../mtn-data/supplementary_datasets.xlsx`. Download the files from [Zenodo](https://zenodo.org/records/18282656) into `../mtn-data/zenodo/`, keeping their published filenames. Optional GEO count matrices (`flowerbuds_counts.csv`, `seedlings_counts.csv`) go directly in `../mtn-data/`.

In R, install the packages needed for the analyses you want to run:

```r
install.packages(c("ggplot2", "dplyr", "readr", "readxl", "tidyr", "rstatix",
  "multcompView", "cowplot", "eulerr", "matrixStats", "UpSetR", "circlize", "viridis"))
install.packages("BiocManager")
BiocManager::install(c("ComplexHeatmap", "rtracklayer", "DESeq2"))
```

## Analyse the downloaded tables

Run these commands from this folder; choose whichever analyses you need:

```sh
Rscript --vanilla prepare_data.R
Rscript --vanilla phenotypes.R
Rscript --vanilla phenotypes_supplement.R
Rscript --vanilla metabolites.R
Rscript --vanilla tor.R
Rscript --vanilla methyl_index.R
Rscript --vanilla rnaseq.R
Rscript --vanilla rnaseq_go.R
```

`prepare_data.R` reads supplementary Datasets S1–S4 and the Zenodo phenotype workbook. RNA-seq plotting uses the supplied contrasts and TPMs; count correlations additionally need the GEO count matrices. TE volcano plots need Araport11 GFF3 at `../mtn-data/reference/Araport11_GFF3_genes_transposons.gff3`. The GO script plots the published term P-values. Calculations use the values in the public downloads. These scripts generate analysis components, not final figure layouts or PCA plots.

## Process sequencing reads

Get paired reads through the SRA links in [GEO GSE288809](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE288809). Match samples by genotype, treatment and replicate. Use these filename prefixes:

| Genotype | MS | Spd |
| --- | --- | --- |
| WT | `wm` | `ws` |
| mtn1-1mtn2-1 | `mm` | `ms` |
| G3 (buds) | `gm` | `gs` |
| mtn1-1mtn2-5 (seedlings) | `5m` | `5s` |

Append the replicate number and mate, e.g. `wm1_R1.fastq.gz` and `wm1_R2.fastq.gz`. Put bud RNA-seq reads in `../mtn-data/raw/rna_bud/` (replicates 1–4), seedling RNA-seq in `raw/rna_seedling/` (1–3), and bud BS-seq in `raw/bs_bud/` (1–3). Combine runs belonging to the same biological sample before assigning these names.

Put TAIR10 FASTA and Araport11 GFF3 in `../mtn-data/reference/`, named `TAIR10.fa` and `Araport11_GFF3_genes_transposons.gff3`. Their chromosome names must agree. Install FastQC, Cutadapt, STAR, StringTie with `prepDE.py3`, Trim Galore, Bismark/Bowtie2 and HOMER as needed, with executables on `PATH`. HOMER also needs its `tair10` resources.

```sh
bash rnaseq.sh rna_bud align
bash rnaseq.sh rna_bud counts
bash rnaseq.sh rna_bud de
# Repeat with rna_seedling instead of rna_bud.

bash bisulfite.sh align
bash bisulfite.sh extract
bash bisulfite.sh profiles
Rscript --vanilla methylation_profiles.R
```

The RNA processing script explicitly uses HOMER/DESeq2. `THREADS`, `REFERENCE`, `ANNOTATION` and `PREPDE` can override defaults. The raw-read commands have not been tested by rerunning the full sequencing experiment. DMR calling is not included.
