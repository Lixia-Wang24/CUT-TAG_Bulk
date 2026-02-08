#!/bin/bash
###======================  download data from GEO ====================###
 prefetch --option-file GEO.download.txt -v -p

###====================== convert sra to fastq =======================###
sample="SRR12246717 SRR11074240 SRR11074254 SRR11074258 SRR11923224 SRR8754611 SRR8754612"
for s in $sample
do
fasterq-dump -e 8 --split-files -O . "$s"
pigz -p 8 "$s"*.fastq
done

###====================== Merge lanes ================================###
cat SRR8754611_1.fastq.gz SRR8754612_1.fastq.gz > IgG_rep2_1.fastq.gz
cat SRR8754611_2.fastq.gz SRR8754612_2.fastq.gz > IgG_rep2_2.fastq.gz

rm SRR8754611_1.fastq.gz SRR8754612_1.fastq.gz SRR8754611_2.fastq.gz SRR8754612_2.fastq.gz

###====================== Rename samples =============================###
mv "SRR12246717_1.fastq.gz" "K27me3_rep1_1.fastq.gz"
mv "SRR12246717_2.fastq.gz" "K27me3_rep1_2.fastq.gz"

mv "SRR11074240_1.fastq.gz" "K27me3_rep2_1.fastq.gz"
mv "SRR11074240_2.fastq.gz" "K27me3_rep2_2.fastq.gz"

mv "SRR11074254_1.fastq.gz" "K4me3_rep1_1.fastq.gz"
mv "SRR11074254_2.fastq.gz" "K4me3_rep1_2.fastq.gz"

mv "SRR11074258_1.fastq.gz" "K4me3_rep2_1.fastq.gz"
mv "SRR11074258_2.fastq.gz" "K4me3_rep2_2.fastq.gz"

mv "SRR11923224_1.fastq.gz" "IgG_rep1_1.fastq.gz"
mv "SRR11923224_2.fastq.gz" "IgG_rep1_2.fastq.gz"

