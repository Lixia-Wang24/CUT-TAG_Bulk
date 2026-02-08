# CUT&Tag Analysis Pipeline

This repository contains a comprehensive analysis pipeline for **CUT&Tag (Cleavage Under Targets and Tagmentation)** sequencing data, focusing on histone modification profiling (**H3K27me3** and **H3K4me3**). The pipeline includes data processing, quality control, peak calling, visualization, and differential analysis with *E. coli* spike-in normalization.

---

## Dependencies

### Command-line Tools
* **Data download and preprocessing**: 
    * SRA Toolkit (`prefetch`, `fasterq-dump`)
    * pigz (parallel gzip)
* **Quality control**: 
    * FastQC (v0.11.9 or higher)
    * MultiQC (v1.9 or higher)
    * fastp (optional, for adapter trimming)
* **Alignment**: 
    * Bowtie2 (v2.4.0 or higher)
    * SAMtools (v1.10 or higher)
    * BEDTools (v2.29 or higher)
* **Duplicate marking**: 
    * Picard Tools (v2.23 or higher)
* **Peak calling**: 
    * SEACR (v1.3)
* **Visualization**: 
    * deepTools (`bamCoverage`, `computeMatrix`, `plotHeatmap`)

### R Packages
* **Core packages**: `library(dplyr)`, `library(stringr)`, `library(ggplot2)`, `library(viridis)`
* **Genomic analysis**: `library(GenomicRanges)`, `library(chromVAR)`, `library(DESeq2)`
* **Visualization**: `library(ggpubr)`, `library(corrplot)`

### Reference Genomes
* **Human genome**: GRCh38/hg38
* **E. coli genome**: K12 MG1655 (for spike-in normalization)

---

## Run Analysis Pipeline

### Step 1: Download data from GEO
`bash 1_Download.sh`
* Downloads SRA files from GEO database.
* Converts SRA to FASTQ format.
* Merges technical replicates.
* Renames files according to sample names.

### Step 2: Quality control
`bash 2_Quality_cutadaptor.sh`
* Runs FastQC on all FASTQ files.
* Aggregates QC reports with MultiQC.
* Optional adapter trimming with `fastp`.

### Step 3: Alignment to hg38 and E. coli
`bash 3_Alignment.sh`

#### 3.1 Alignment to Human Genome (hg38)
**Parameters:**
* `--end-to-end`: End-to-end alignment.
* `--very-sensitive`: High sensitivity mode.
* `--no-mixed`: Suppress unpaired alignments.
* `--no-discordant`: Suppress discordant alignments.
* `-I 10`: Minimum fragment length.
* `-X 700`: Maximum fragment length.

#### 3.2 Alignment to E. coli (Spike-in)
* Calculates spike-in fragment counts for normalization.
* Generates scaling factors.

#### 3.3 Quality Filtering
* Removes multimapped reads (MAPQ < 2).
* Filters unmapped reads.
* Keeps properly paired reads.
* Removes fragments > 1000 bp.

#### 3.4 Duplication Analysis
* Uses **Picard MarkDuplicates**.
* Calculates duplication rates and estimates library complexity.

**Outputs:** SAM/BAM files, fragment length distributions, alignment statistics, and duplication metrics.

### Step 4: Spike-in calibration
`bash 4_Spikein_calibration.sh`
* Calculates scaling factors based on *E. coli* spike-in reads.
* **Formula**: `scaling_factor = 10000 / spike_in_depth`
* Generates normalized **bedGraph** files.
* *Purpose*: Normalizes for differences in antibody efficiency and cell number across samples.

### Step 5: Peak calling with SEACR
`bash 5_Peak_calling.sh`

**Two modes supported:**
1.  **Control mode** (Uses IgG control):
    ```bash
    bash SEACR_1.3.sh sample.bedgraph IgG.bedgraph non stringent output
    ```
2.  **Top percentile mode** (Calls top 1% peaks):
    ```bash
    bash SEACR_1.3.sh sample.bedgraph 0.01 non stringent output
    ```
* **non**: Non-normalized mode (bedGraph already normalized).
* **stringent**: Stringent threshold.

### Step 6: Generate heatmaps
`bash 6_Heatmap_regions.sh`

* **6.1 BigWig Generation**: Converts BAM to bigWig for genome browser visualization (IGV, UCSC).
* **6.2 Heatmaps**:
    * **Over Gene Bodies**: 3 kb upstream of TSS, 5 kb gene body, 3 kb downstream of TES.
    * **Over Peak Summits**: ±3 kb around peak center using K-means clustering (k=2).

### Step 7: R analysis and visualization

#### 1. Alignment Summary
`Rscript Alignment_summary.R`
* **Metrics**: Sequencing depth, alignment rates, mapped fragment counts, duplication rates, library complexity, fragment size distributions, and scaling factors.
* **Outputs**: `Alignment_summary_figure.pdf`, `Fragment_size_figure.pdf`, `Alignment_summary_table.csv`.

#### 2. Replicate Reproducibility
`Rscript Replicate_reproducibility.R`
* Assesses consistency between replicates.
* **Outputs**: `SpikeIn_Normalized_count.pdf`, `Correlation_plots.pdf`.

#### 3. Peaks Summary
`Rscript Peaks_summary.R`
* **Metrics**: Peak counts, width distributions, and **FRiP (Fraction of Reads in Peaks)** scores.
* **Quality Standard**: FRiP Score should be > 1% for high quality; Peak reproducibility aim is > 70%.
* **Outputs**: `Peaks_summary.pdf`, `Peaks_summary_table.csv`.

#### 4. Differential Analysis
`Rscript Differential_analysis.R`
* Uses **DESeq2** for differential binding analysis.
* Creates a master peak list (union) and counts fragments.
* Performs statistical testing (e.g., H3K27me3 vs H3K4me3).
* **Outputs**: `DESeq2_results_all_counts.csv`.

---

## Expected Results

### Quality Metrics Reference Table
| Metric | Expected Value |
| :--- | :--- |
| **Alignment Rate (hg38)** | > 80% |
| **Duplication Rate** | < 20% |
| **FRiP Score** | > 1-5% |
| **Replicate Correlation** | > 0.9 |
| **Peak Reproducibility** | > 70% |

### Biological Expectations
* **H3K4me3**: Narrow peaks (~1-2 kb) concentrated at promoters/TSSs.
* **H3K27me3**: Broad peaks (~10-50 kb) covering repressed genomic domains.
* **IgG**: Few to no peaks (background noise).
