#!/bin/bash
### ==================== Quality Control =======================###
qc_outdir="./qc_reports"
mkdir -p ${qc_outdir}

fastqc -o ${qc_outdir} -t 8 *.fastq.gz

multiqc ${qc_outdir} -o ${qc_outdir}/mutiqc

### ==================== Cut adaptor =========================###
# There is no need to trim reads from out standard 25x25 PE sequencing
# as adapter sequences will not be included in reads of inserts >25 bp

:<<!
### ==================== Cut adaptor =========================###
trim_outdir="./trimmed.reads"
mkdir -p ${trim_outdir}

sample="K27me3_rep1 K27me3_rep2 K4me3_rep1 K4me3_rep2 IgG_rep1 IgG_rep2"
for s in $sample
do
fastp \
        -i "${s}_1.fastq.gz"
        -I "${s}_2.fastq.gz" \
        -o "${s}.1.trimmed.fastq.gz" \
        -O "${s}.2.trimmed.fastq.gz" \
        --detect_adapter_for_pe \
        --trim_poly_g \
        --correction \
        --length_required 30 \
        --qualified_quality_phred 20 \
        --unqualified_percent_limit 30 \
        --thread 8 \
        --html "${s}.fastp.html" \
        --json "${s}.fastp.json"
done
!
