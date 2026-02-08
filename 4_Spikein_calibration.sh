#!/bin/bash
:<<!
### =================== chrom.sizes =========================================== ###
cd /Users/lixia/Data/database/ref_genome/GRCh38.UCSC/bowtie2.index/GRCh38.index

bowtie2-inspect -s GRCh38_noalt_as | grep "^Sequence" | awk '{print $2"\t"$NF}' > hg38.chrom.sizes

head hg38.chrom.sizes

### =================== Spike-in calibration ================================== ###
cd /Users/lixia/Data/data/Cut_tag
chromSize="/Users/lixia/Data/database/ref_genome/GRCh38.UCSC/bowtie2.index/GRCh38.index/hg38.chrom.sizes"

sample="K27me3_rep1 K27me3_rep2 K4me3_rep1 K4me3_rep2 IgG_rep1 IgG_rep2"
for s in $sample
do

# Read spike in sequence depth
seqDepthFile="./alignment/sam/bowtie2_summary/${s}_bowtie2_spikeIn.seqDepth"

if [ -f "$seqDepthFile" ]; then
    seqDepth=$(cat $seqDepthFile)
else
    echo "Error: Spike-in depth file not found for ${s}"
    exit 1
fi

# Scaling factor
if [[ "$seqDepth" -gt "1" ]]; then
    mkdir -p ./alignment/bedgraph

    scale_factor=`echo "10000 / $seqDepth" | bc -l`

    echo "Scaling factor for ${s} is: ${scale_factor} (E. coli reads: ${seqDepth})"

# BedGraph
bedtools genomecov -bg -scale ${scale_factor} \
        -i ./alignment/bed/${s}_bowtie2.fragments.bed \
        -g ${chromSize} \
        > ./alignment/bedgraph/${s}_bowtie2.fragments.normalized.bedgraph
fi

done
!
### =================== Add scaling factor to summary.table ===================== ###
source /Users/lixia/Data/Ubin/miniforge3/etc/profile.d/conda.sh
conda activate base
Rscript Alignment_summary.R
conda activate RNAexpdiff
