### Trim, align, clean, mark PCR duplicates. Generate QC, alignment and coverage reports 
### GNU parallel by sample

module load fastp
module load fastqc
module load bwa
module load samtools

# create directory variables 
DATDIR="/g/data/ht96/McLay_UQ/avneet_paper/D01H01"
QCDIR="/g/data/ht96/McLay_UQ/avneet_paper/qc"
ALDIR="/g/data/ht96/McLay_UQ/avneet_paper/align"
CLDIR="/g/data/ht96/McLay_UQ/avneet_paper/clean"

# fastp - trim etc. 
fastp --in1 ${DATDIR}/${1}_1.fq.gz \
    --in2 ${DATDIR}/${1}_2.fq.gz \
    --out1 ${QCDIR}/${1}_R1_trimmed.fastq.gz --out2 ${QCDIR}/${1}_R2_trimmed.fastq.gz \
    --unpaired1 ${QCDIR}/${1}_R1_unpaired.fastq.gz --unpaired2 ${QCDIR}/${1}_R2_unpaired.fastq.gz \
    -q 10 -u 50 -l 50 -h ${QCDIR}/${1}.html &> ${QCDIR}/${1}.log 

# fastQC - quality report. 
fastqc ${QCDIR}/${1}_R1_trimmed.fastq.gz ${QCDIR}/${1}_R2_trimmed.fastq.gz -o ${QCDIR}

# align and add read groups with BWA, sort the bam file with Samtools  
bwa mem -t 12 -M -R "@RG\tSM:${1}\tID:${1}\tLB:${1}\tPL:ILLUMINA\tPU:${1}" /scratch/ht96/km6006/SLv141Asm_Ch20RN.fasta \
    ${QCDIR}/${1}_R1_trimmed.fastq.gz ${QCDIR}/${1}_R2_trimmed.fastq.gz | samtools sort -@ 12 -T ${ALDIR}/${1} -o ${ALDIR}/${1}_trm_srt.bam 

# alignment stats report.
samtools flagstat ${ALDIR}/${1}_trm_srt.bam &> ${ALDIR}/${1}_trm_srt_stats.txt

# clean the bam file with Picard 
java -Xmx8g -jar /home/564/km6006/picard.jar CleanSam INPUT=${ALDIR}/${1}_trm_srt.bam OUTPUT=${CLDIR}/${1}_cln_srt.bam

# marking PCR duplicates with Picard MarkDuplicates
java -XX:ConcGCThreads=1 -XX:ParallelGCThreads=1 -Xmx4g -jar /home/564/km6006/picard.jar MarkDuplicates \
        INPUT=${CLDIR}/${1}_cln_srt.bam \
        OUTPUT=${CLDIR}/${1}_PCRm_cln_srt.bam \
        ASSUME_SORTED=true \
        VALIDATION_STRINGENCY=SILENT \
        READ_NAME_REGEX=null \
        MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=900 \
        METRICS_FILE=${CLDIR}/${1}_PCRm_cln_srt_metrics

# index the .bam file
samtools index ${CLDIR}/${1}_PCRm_cln_srt.bam

# check coverage
samtools coverage ${CLDIR}/${1}_PCRm_cln_srt.bam --output ${CLDIR}/${1}_coverage.txt


