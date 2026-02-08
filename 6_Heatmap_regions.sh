#!/bin/bash
:<<!
### ================ convert bam into Bigwig .bw file format ================== ###
mkdir -p ./alignment/bigwig

sample="K27me3_rep1 K27me3_rep2 K4me3_rep1 K4me3_rep2 IgG_rep1 IgG_rep2"
for s in $sample
do

samtools sort \
    ./alignment/bam/${s}_bowtie2.mapped.bam \
    -o ./alignment/bam/${s}.sorted.bam

samtools index \
    ./alignment/bam/${s}.sorted.bam

bamCoverage \
    -b ./alignment/bam/${s}.sorted.bam \
    -o ./alignment/bigwig/${s}_raw.bw

done
!
:<<!
### ================ Heatmap over transcription units ========================== ###
Threads=8
hg38_gene="/Users/lixia/Data/database/ref_genome/GRCh38.UCSC/bowtie2.index/hg38_gene.tsv"
mkdir -p ./heatmap/trans_units

computeMatrix scale-regions \
    -S ./alignment/bigwig/K27me3_rep1_raw.bw \
       ./alignment/bigwig/K27me3_rep2_raw.bw \
       ./alignment/bigwig/K4me3_rep1_raw.bw \
       ./alignment/bigwig/K4me3_rep2_raw.bw \
       ./alignment/bigwig/IgG_rep1_raw.bw \
       ./alignment/bigwig/IgG_rep2_raw.bw \
    -R ${hg38_gene} \
    --beforeRegionStartLength 3000 \
    --regionBodyLength 5000 \
    --afterRegionStartLength 3000 \
    --skipZeros \
    -o ./heatmap/trans_units/matrix_gene.mat.gz \
    -p ${Threads}

plotHeatmap \
    -m ./heatmap/trans_units/matrix_gene.mat.gz \
    -out ./heatmap/trans_units/Histone_gene.png \
    --sortUsing sum
!
:<<!
### ================ Heatmap on CUT&Tag peaks ============= ###
Threads=8
mkdir -p ./heatmap/peaks

sample="K27me3_rep1 K4me3_rep1"
for s in $sample
do

# Convert "chr:start-end" into "chr start end" format
awk '{split($6, summit, ":"); split(summit[2], region, "-"); print summit[1]"\t"region[1]"\t"region[2]}' \
    ./peakCalling/SEACR/${s}_seacr_control.peaks.stringent.bed > \
    ./peakCalling/SEACR/${s}_seacr_control.peaks.summitRegion.bed

# Peak matrix
computeMatrix reference-point \
    -S ./alignment/bigwig/${s}_raw.bw \
    -R ./peakCalling/SEACR/${s}_seacr_control.peaks.summitRegion.bed \
    --skipZeros \
    -o ./peakCalling/SEACR/${s}_SEACR.mat.gz \
    -p ${Threads} \
    -a 3000 \
    -b 3000 \
    --referencePoint center

# Figure
plotHeatmap \
    -m ./peakCalling/SEACR/${s}_SEACR.mat.gz \
    -out ./heatmap/peaks/${s}_SEACR_heatmap.png \
    --sortUsing sum \
    --startLabel "Peak Start" \
    --endLabel "Peak End" \
    --xAxisLabel "" \
    --regionsLabel "Peaks" \
    --samplesLabel "${s}"
done
!

### ================ Heatmap on CUT&Tag peaks (combined)============= ###
Threads=8
mkdir -p ./heatmap/peaks

sample="K27me3_rep1 K27me3_rep2 K4me3_rep1 K4me3_rep2"
for s in $sample
do
# Convert "chr:start-end" into "chr start end" format
awk '{split($6, summit, ":"); split(summit[2], region, "-"); print summit[1]"\t"region[1]"\t"region[2]}' \
    ./peakCalling/SEACR/${s}_seacr_control.peaks.stringent.bed > \
    ./peakCalling/SEACR/${s}_seacr_control.peaks.summitRegion.bed
done

# Combine summitRegion.bed
cat ./peakCalling/SEACR/*_seacr_control.peaks.summitRegion.bed | sort -k1,1 -k2,2n | bedtools merge > ./peakCalling/SEACR/All_Merged_Summits.bed

# Combine peak matrix
computeMatrix reference-point \
    -S ./alignment/bigwig/K27me3_rep1_raw.bw \
       ./alignment/bigwig/K27me3_rep2_raw.bw \
       ./alignment/bigwig/K4me3_rep1_raw.bw \
       ./alignment/bigwig/K4me3_rep2_raw.bw \
    -R ./peakCalling/SEACR/All_Merged_Summits.bed \
    --skipZeros \
    -o ./peakCalling/SEACR/metrix_allpeaks_SEACR.mat.gz \
    -p ${Threads} \
    -a 3000 \
    -b 3000 \
    --referencePoint center

# Combined figure
plotHeatmap \
    -m ./peakCalling/SEACR/metrix_allpeaks_SEACR.mat.gz \
    -out ./heatmap/peaks/Heatmap_allpeaks_SEACR.png \
    --sortUsing sum \
    --colorMap Reds \
    --kmeans 2 \
    --samplesLabel "K27me3_rep1" "K27me3_rep2" "K4me3_rep1" "K4me3_rep2" \
    --zMax 15 15 50 50
