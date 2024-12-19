#!/bin/bash
#author:g. ozan bozdag
#Description:
#This pipeline processes paired-end sequencing data to identify variants.
#After obtaining variant calls (VCF), we first run a quality-based filtering. Next, we use 'bcftools isec' to compare variants
#between the evolved and ancestral VCF files, retaining
#mutations present in the evolved strain but not in the ancestor. Additional filtering
#steps, including cross-comparison filtering of indels that are located in hard-to-map regions
#and appear as false positives across strains from independent lineages, are performed to ensure high-quality variant calls.

# Define the reference genome paths for alignment and variant calling
REFERENCE_GENOME="./S288C_reference_sequence_R64-3-1_20210421.fsa"
REFERENCE_GENOME_WITH_DICTIONARY="./S288C_ref.fa"
GATK_JAR="/usr/local/bin/gatk-4.2.4.1/gatk-package-4.2.4.1-local.jar"
PICARD_JAR="$WORK/usr/local/bin/picard/build/libs/picard.jar"
TRIMMOMATIC_JAR="$WORK/usr/local/bin/trimmomatic-0.39-hdfd78af_2/share/trimmomatic-0.39-2/trimmomatic.jar"

# Check if the reference genome files are present
if [[ ! -f "$REFERENCE_GENOME" ]] || [[ ! -f "$REFERENCE_GENOME_WITH_DICTIONARY" ]]; then
    echo "ERROR: One or both reference genome files not found."
    echo "Please ensure the following files are in the current folder:"
    echo "1) $REFERENCE_GENOME"
    echo "2) $REFERENCE_GENOME_WITH_DICTIONARY"
    exit 1
fi

# Ask user to confirm that the reference genomes are correct
echo "Please confirm that the reference genome files are correct and in the current directory."
echo "Press 'Y' to continue or any other key to exit."
read -n 1 -r user_input
echo # move to a new line
if [[ ! $user_input =~ ^[Yy]$ ]]; then
    echo "User did not confirm. Exiting script."
    exit 1
fi

echo "Proceeding with pipeline..."

# Step 0: Index the main reference genome
bwa index -a bwtsw $REFERENCE_GENOME

# Step 1: Trim and filter out low quality reads using Trimmomatic
for R1 in *R1*.fastq.gz; do
    R2="${R1//R1/R2}"
    R1paired="${R1//.fastq.gz/_paired.fastq.gz}"
    R1unpaired="${R1//.fastq.gz/_unpaired.fastq.gz}"
    R2paired="${R2//.fastq.gz/_paired.fastq.gz}"
    R2unpaired="${R2//.fastq.gz/_unpaired.fastq.gz}"
    java -jar $TRIMMOMATIC_JAR PE -phred33 $R1 $R2 $R1paired $R1unpaired $R2paired $R2unpaired LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:36
done

# Step 2: Align trimmed reads to the reference genome using BWA MEM
for R1 in *_paired.fastq.gz; do
    R2="${R1//R1/R2}"
    output="${R1//_R1_paired.fastq.gz/.sam}"
    bwa mem $REFERENCE_GENOME $R1 $R2 > $output
done

# Step 3: Sort and convert SAM to BAM files using Picard's SortSam
for samfile in *.sam; do
    bamfile="${samfile//.sam/.bam}"
    java -jar $PICARD_JAR SortSam INPUT=$samfile OUTPUT=$bamfile SORT_ORDER=coordinate
done

# Step 4: Index BAM files using Picard's BuildBamIndex
for bamfile in *.bam; do
    java -jar $PICARD_JAR BuildBamIndex INPUT=$bamfile
done

# Step 5: Mark duplicates in BAM files using Picard's MarkDuplicates
for bamfile in *.bam; do
    dedup_bamfile="${bamfile%.*}_dedup.bam"
    metrics_file="${bamfile%.*}_metrics.txt"
    java -jar $PICARD_JAR MarkDuplicates VALIDATION_STRINGENCY=LENIENT INPUT=$bamfile OUTPUT=$dedup_bamfile METRICS_FILE=$metrics_file
done

# Step 6: Fix read group errors in BAM files using Picard's AddOrReplaceReadGroups
for bamfile in *_dedup.bam; do
    fixed_bamfile="${bamfile%.*}_fixed.bam"
    java -jar $PICARD_JAR AddOrReplaceReadGroups VALIDATION_STRINGENCY=LENIENT I=$bamfile O=$fixed_bamfile RGID=4 RGLB=lib1 RGPL=illumina RGPU=unit1 RGSM=20
done

# Step 7: Index final BAM files with fixed read groups using Samtools
for bamfile in *_fixed.bam; do
    samtools index $bamfile
done

# Step 8: Calculate coverage for indexed BAM files
echo "Calculating coverage for BAM files..."
for file in *fixed.bam; do
    samtools depth -aa $file > "${file%.*}.coverage"
done

# Step 9: Index the alternative reference genome
samtools faidx $REFERENCE_GENOME_WITH_DICTIONARY

# Step 10: Generate a sequence dictionary for the alternative reference genome
java -jar $PICARD_JAR CreateSequenceDictionary R=$REFERENCE_GENOME_WITH_DICTIONARY O=${REFERENCE_GENOME_WITH_DICTIONARY%.*}.dict

# Step 11: Call variants using GATK's HaplotypeCaller
for bamfile in *_fixed.bam; do
    vcf_output="${bamfile%.*}.vcf"
    java -jar $GATK_JAR HaplotypeCaller --native-pair-hmm-threads 20 -I $bamfile -O $vcf_output -R $REFERENCE_GENOME_WITH_DICTIONARY
done

# Echo completion message
echo "Pipeline completed. VCF files generated."
# Note: After this pipeline, use bcftools isec to compare evolved and ancestral VCFs
# Then perform additional filtering steps including manual curation