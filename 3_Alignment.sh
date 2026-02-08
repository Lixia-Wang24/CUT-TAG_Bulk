#!/bin/bash
:<<!
### ===================== build bowtie2 index ==================================== ###
##### download pre-builded human index ###
hg38_in_dir="/Users/lixia/Data/database/ref_genome/GRCh38.UCSC/bowtie2.index"
cd "$hg38_in_dir"
wget https://genome-idx.s3.amazonaws.com/bt/GRCh38_noalt_as.zip
unzip GRCh38_noalt_as.zip
##### bulid E.coli index ###
cd "/Users/lixia/Data/database/ref_genome/E.coli.K12.MG1655"

wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/GCF_000005845.2_ASM584v2/GCF_000005845.2_ASM584v2_genomic.fna.gz

gunzip GCF_000005845.2_ASM584v2_genomic.fna.gz

mkdir -p bowtie2.index
bowtie2-build GCF_000005845.2_ASM584v2_genomic.fna bowtie2.index/Ecoli
!
:<<!
### ===================== Allign to human hg38 genome =========================== ###
cd "/Users/lixia/Data/data/Cut_tag"

hg38_in_dir="/Users/lixia/Data/database/ref_genome/GRCh38.UCSC/bowtie2.index/GRCh38.index/GRCh38_noalt_as"
Threads=8
mkdir -p ./alignment/sam
mkdir -p ./alignment/sam/bowtie2_summary

sample="K27me3_rep1 K27me3_rep2 K4me3_rep1 K4me3_rep2 IgG_rep1 IgG_rep2"
for s in $sample
do
bowtie2 \
    --end-to-end \
    --very-sensitive \
    --no-mixed \
    --no-discordant \
    --phred33 \
    -I 10 \
    -X 700 \
    -p ${Threads} \
    -x ${hg38_in_dir} \
    -1 ${s}_1.fastq.gz \
    -2 ${s}_2.fastq.gz \
    -S ./alignment/sam/${s}_bowtie2.sam \
    &> ./alignment/sam/bowtie2_summary/${s}_bowtie2.txt

done
!
:<<!
### ===================== Allign to E.coli genome ================================ ###
Ecoli_in_dir="/Users/lixia/Data/database/ref_genome/E.coli.K12.MG1655/bowtie2.index/Ecoli"
Threads=8

sample="K27me3_rep1 K27me3_rep2 K4me3_rep1 K4me3_rep2 IgG_rep1 IgG_rep2"
for s in $sample
do
bowtie2 \
    --end-to-end \
    --very-sensitive \
    --no-overlap \
    --no-dovetail \
    --no-mixed \
    --no-discordant \
    --phred33 \
    -I 10 \
    -X 700 \
    -p ${Threads} \
    -x ${Ecoli_in_dir} \
    -1 ${s}_1.fastq.gz \
    -2 ${s}_2.fastq.gz \
    -S ./alignment/sam/${s}_bowtie2_spikeIn.sam \
    &> ./alignment/sam/bowtie2_summary/${s}_bowtie2_spikeIn.txt

#### Calculate the E. coli fragemnts, which can be used to normalize epitope abundance
seqDepthDouble=$(samtools view -F 0x04 ./alignment/sam/${s}_bowtie2_spikeIn.sam | wc -l)
seqDepth=$((seqDepthDouble/2))
echo ${seqDepth} > ./alignment/sam/bowtie2_summary/${s}_bowtie2_spikeIn.seqDepth

done
!
:<<!
### ===================== Check the duplication rate =========================== ###
Threads=8
picardCMD="java -jar /Users/lixia/Data/Ubin/picard.jar"
mkdir -p ./alignment/removeDuplicate/picard_summary

sample="K27me3_rep1 K27me3_rep2 K4me3_rep1 K4me3_rep2 IgG_rep1 IgG_rep2"
for s in $sample
do

### Add Read Group
$picardCMD AddOrReplaceReadGroups \
       I=./alignment/sam/${s}_bowtie2.sam \
       O=./alignment/sam/${s}_bowtie2.sorted.sam \
       SORT_ORDER=coordinate \
       RGID=${s} \
       RGLB=lib1 \
       RGPL=illumina \
       RGPU=unit1 \
       RGSM=${s}

### Mark duplicates
$picardCMD MarkDuplicates \
    I=./alignment/sam/${s}_bowtie2.sorted.sam \
    O=./alignment/removeDuplicate/${s}_bowtie2.sorted.dupMarked.sam \
    METRICS_FILE=./alignment/removeDuplicate/picard_summary/${s}_picard.dupMark.txt \
    REMOVE_DUPLICATES=false

### Remove duplicates
$picardCMD MarkDuplicates \
    I=./alignment/sam/${s}_bowtie2.sorted.sam \
    O=./alignment/removeDuplicate/${s}_bowtie2.sorted.rmDup.sam \
    METRICS_FILE=./alignment/removeDuplicate/picard_summary/${s}_picard.rmDup.txt \
    REMOVE_DUPLICATES=true
done
!
:<<!
### ===================== Fragment size distribution =========================== ###
mkdir -p ./alignment/sam/fragmentLen

sample="K27me3_rep1 K27me3_rep2 K4me3_rep1 K4me3_rep2 IgG_rep1 IgG_rep2"
for s in $sample
do
## Extract the 9th column from the alignment sam file which is the fragment length
samtools view -F 0x04 ./alignment/sam/${s}_bowtie2.sam | awk -F'\t' 'function abs(x){return ((x < 0.0) ? -x : x)} {print abs($9)}' | sort | uniq -c | awk -v OFS="\t" '{print $2, $1/2}' > ./alignment/sam/fragmentLen/${s}_fragmentLen.txt
done
!
:<<!
### ===================== Report sequencing mapping summary ==================== ###
source /Users/lixia/Data/Ubin/miniforge3/etc/profile.d/conda.sh
conda activate base
Rscript Alignment_summary.R
conda activate RNAexpdiff
!
:<<!
### ===================== Filter low quality reads =========================== ###
mkdir -p ./alignment/bam
mkdir -p ./alignment/bed

sample="K27me3_rep1 K27me3_rep2 K4me3_rep1 K4me3_rep2 IgG_rep1 IgG_rep2"
for s in $sample
do

### Filter multimapped reads and keep unique mapped reads
minQualityScore=2
samtools view -h -q ${minQualityScore} ./alignment/sam/${s}_bowtie2.sam > ./alignment/sam/${s}_bowtie2.qualityScore${minQualityScore}.sam

## Filter and keep the mapped read pairs (convert sam into bam)
samtools view -bS -F 0x04 ./alignment/sam/${s}_bowtie2.qualityScore${minQualityScore}.sam  > ./alignment/bam/${s}_bowtie2.mapped.bam

## Convert bam into bed file format
bedtools bamtobed -i ./alignment/bam/${s}_bowtie2.mapped.bam -bedpe > ./alignment/bed/${s}_bowtie2.bed

## Keep the read pairs that are on the same chromosome and fragment length less than 1000bp.
awk '$1==$4 && $6-$2 < 1000 {print $0}' ./alignment/bed/${s}_bowtie2.bed > ./alignment/bed/${s}_bowtie2.clean.bed

## Only extract the fragment related columns
cut -f 1,2,6 ./alignment/bed/${s}_bowtie2.clean.bed | sort -k1,1 -k2,2n -k3,3n  > ./alignment/bed/${s}_bowtie2.fragments.bed

done
!
:<<!
### =====================  Assess replicate reproducibility ====================== ###
binLen=500
sample="K27me3_rep1 K27me3_rep2 K4me3_rep1 K4me3_rep2 IgG_rep1 IgG_rep2"
for s in $sample
do
awk -v w=${binLen} '{print $1, int(($2 + $3)/(2*w))*w + w/2}' ./alignment/bed/${s}_bowtie2.fragments.bed | sort -k1,1V -k2,2n | uniq -c | awk -v OFS="\t" '{print $2, $3, $1}' |  sort -k1,1V -k2,2n  > ./alignment/bed/${s}_bowtie2.fragmentsCount.bin${binLen}.bed

done
!
### ===================== Visualize replicate reproducibility==================== ###
source /Users/lixia/Data/Ubin/miniforge3/etc/profile.d/conda.sh
conda activate base
Rscript Replicate_reproducibility.R
conda activate RNAexpdiff
