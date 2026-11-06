# MTN paper analysis

Shell and R scripts for [Methylthioadenosine-driven rewiring of sulfur assimilation lowers the methyl index and fertility in Arabidopsis](https://doi.org/10.1093/jxb/erag135).

## Data and dependencies

Keep this folder as code only. By default the scripts read `../mtn-data/` and write `../mtn-output/`. Set `MTN_DATA` and `MTN_OUTPUT` to use other locations.

Download the supplementary workbook from the paper as `../mtn-data/supplementary_datasets.xlsx`. Download the files from [Zenodo](https://zenodo.org/records/18282656) into `../mtn-data/zenodo/`, keeping their published filenames.

Download the [reference, PCA and GO inputs](https://github.com/bjmt/Tremblay_et_al_2026_MTN/releases/download/rnaseq-counts/mtn-analysis-inputs.tar.gz) into the same data folder:

```sh
mkdir -p ../mtn-data
curl -fL -o ../mtn-data/mtn-analysis-inputs.tar.gz \
  https://github.com/bjmt/Tremblay_et_al_2026_MTN/releases/download/rnaseq-counts/mtn-analysis-inputs.tar.gz
tar -xzf ../mtn-data/mtn-analysis-inputs.tar.gz -C ../mtn-data
```

This supplies the exact TAIR10 FASTA and Araport11 GFF3, full internode measurements, metabolite PCA coordinates and GO enrichment exports used for plotting.

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
Rscript --vanilla development.R
Rscript --vanilla pca.R
Rscript --vanilla phenotypes_supplement.R
Rscript --vanilla metabolites.R
Rscript --vanilla tor.R
Rscript --vanilla methyl_index.R
Rscript --vanilla rnaseq.R
Rscript --vanilla rnaseq_go.R
Rscript --vanilla dmrs.R
```

`prepare_data.R` reads supplementary Datasets S1–S5 and the Zenodo phenotype workbook. RNA-seq plotting uses the supplied contrasts and TPMs; count correlations additionally use the original count matrices described below. TE volcano plots and DMR/gene overlaps need Araport11 GFF3 at `../mtn-data/reference/Araport11_GFF3_genes_transposons.gff3`.

## Scripts and figures

Figure numbers refer to the published paper.

| Figure | Script(s) |
| --- | --- |
| 1B–E | `phenotypes.R` |
| 2B, D, G, H | `development.R` |
| 2C, 3B–C | `pca.R` |
| 3D–H | `metabolites.R` |
| 4A–E | `metabolites.R` |
| 5A | `metabolites.R` |
| 5C–D | `tor.R` |
| 6A–B, D–I | `rnaseq.R` |
| 6C | `rnaseq_go.R` |
| 7A–B | `methyl_index.R` |
| 7C–F | `bisulfite.sh`, `methylation_profiles.R` |
| S1A–B | `phenotypes.R` |
| S2A–B | `phenotypes_supplement.R` |
| S3A–C | `metabolites.R` |
| S5C–D | `tor.R` |
| S6–S7 | `rnaseq.R` |
| S8A–C, S10A–B | `dmrs.R` |
| S9 | `dmrs.R --go` |

`prepare_data.R` supplies inputs for the table-based analyses. `rnaseq.sh` and `annotate_rnaseq.R` provide the raw-read RNA-seq processing route; `common.R` supplies shared paths and package checks.

To recalculate the metabolite PCA from measurements, install `pcaMethods` with `BiocManager::install("pcaMethods")` and run `Rscript --vanilla pca.R --recalculate`. The internode PCA is calculated from measurements in both modes.

To run fresh GO enrichment through g:Profiler, install `gprofiler2` with `install.packages("gprofiler2")`, then use `Rscript --vanilla rnaseq_go.R --recalculate` or `Rscript --vanilla dmrs.R --go`. The RNA-seq GO script also accepts a contrast folder after `--recalculate`, such as `../mtn-output/rnaseq_de`.

## Recalculate RNA-seq differential expression

Download the [original StringTie/prepDE count matrices](https://github.com/bjmt/Tremblay_et_al_2026_MTN/releases/download/rnaseq-counts/mtn-rnaseq-counts.tar.gz) and extract them into `../mtn-data/rnaseq/`:

```sh
mkdir -p ../mtn-data/rnaseq
curl -fL -o ../mtn-data/rnaseq/mtn-rnaseq-counts.tar.gz \
  https://github.com/bjmt/Tremblay_et_al_2026_MTN/releases/download/rnaseq-counts/mtn-rnaseq-counts.tar.gz
tar -xzf ../mtn-data/rnaseq/mtn-rnaseq-counts.tar.gz -C ../mtn-data/rnaseq
```

Alternatively, generate the matrices from reads with the `align` and `counts` stages below; the R scripts also read those outputs directly.

The following applies the original count/rRNA filters and fits each comparison with DESeq2. Run `prepare_data.R` first to supply supplementary-table TPMs for plotting.

```sh
Rscript --vanilla rnaseq_de.R
Rscript --vanilla rnaseq.R ../mtn-output/rnaseq_de
```

Use the actual output path in the second command if `MTN_OUTPUT` is set. `rnaseq_de.R bud` or `rnaseq_de.R seedling` runs one tissue only; `rnaseq.R` expects both tissues.

## Process sequencing reads

Get paired reads through the SRA links in [GEO GSE288809](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE288809). Match samples by genotype, treatment and replicate. Use these filename prefixes:

| Genotype | MS | Spd |
| --- | --- | --- |
| WT | `wm` | `ws` |
| mtn1-1mtn2-1 | `mm` | `ms` |
| G3 (buds) | `gm` | `gs` |
| mtn1-1mtn2-5 (seedlings) | `5m` | `5s` |

Append the replicate number and mate, e.g. `wm1_R1.fastq.gz` and `wm1_R2.fastq.gz`. Put bud RNA-seq reads in `../mtn-data/raw/rna_bud/` (replicates 1–4), seedling RNA-seq in `raw/rna_seedling/` (1–3), and bud BS-seq in `raw/bs_bud/` (1–3). Combine runs belonging to the same biological sample before assigning these names.

The reference download above places `TAIR10.fa` and `Araport11_GFF3_genes_transposons.gff3` in `../mtn-data/reference/`. Install FastQC, Cutadapt, STAR, StringTie with `prepDE.py3`, Trim Galore, Bismark/Bowtie2 and HOMER as needed, with executables on `PATH`. HOMER also needs its `tair10` resources.

```sh
bash rnaseq.sh rna_bud align
bash rnaseq.sh rna_bud counts
bash rnaseq.sh rna_bud de
# Repeat with rna_seedling instead of rna_bud.

# After processing both RNA-seq tissues:
Rscript --vanilla rnaseq.R ../mtn-output/preprocessing/rnaseq

bash bisulfite.sh align
bash bisulfite.sh extract
bash bisulfite.sh profiles
Rscript --vanilla methylation_profiles.R
```

The RNA processing script uses HOMER/DESeq2. `THREADS`, `REFERENCE`, `ANNOTATION` and `PREPDE` can override defaults.
