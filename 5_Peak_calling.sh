#!/bin/bash
seacr="/Users/lixia/Data/Ubin/miniforge3/envs/RNAexpdiff/bin/SEACR_1.3.sh"
mkdir -p ./peakCalling/SEACR

sample="K27me3_rep1 K27me3_rep2 K4me3_rep1 K4me3_rep2"
for s in $sample
do
# Set control ssample
    if [[ "$s" == *"rep1"* ]]; then
        histControl="IgG_rep1"
    elif [[ "$s" == *"rep2"* ]]; then
        histControl="IgG_rep2"
    else
        histControl="IgG_rep1"
    fi

    echo "Calling peaks for $s vs $histControl ..."

# Run SEACR
    bash $seacr \
       ./alignment/bedgraph/${s}_bowtie2.fragments.normalized.bedgraph \
       ./alignment/bedgraph/${histControl}_bowtie2.fragments.normalized.bedgraph \
       non stringent \
       ./peakCalling/SEACR/${s}_seacr_control.peaks

    bash $seacr \
        ./alignment/bedgraph/${s}_bowtie2.fragments.normalized.bedgraph \
        0.01 non stringent \
        ./peakCalling/SEACR/${s}_seacr_top0.01.peaks

done
