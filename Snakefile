SRA = "SRR1972739"
REF_ID = "AF086833.2"
RESULTS_FOLDER = "results"
RAW_DIR=f"{RESULTS_FOLDER}/raw"
ALIGNED_DIR=f"{RESULTS_FOLDER}/aligned"
VARIANT_DIR=f"{RESULTS_FOLDER}/variants"
ANNOTATED_DIR=f"{RESULTS_FOLDER}/annotated"
QC_DIR=f"{RESULTS_FOLDER}/qc"
SNPEFF_DIR=f"{RESULTS_FOLDER}/snpEff"
SNPEFF_DATA_DIR=f"{SNPEFF_DIR}/data/reference_db"
SNAKEMAKE_DIR=f"{RESULTS_FOLDER}/snakemake"

        
rule create_dirs:
    output:
        marker = f"{SNAKEMAKE_DIR}/.dirs_created"
    shell:
        """
        mkdir -p {RESULTS_FOLDER} {RAW_DIR} {ALIGNED_DIR} {VARIANT_DIR} {ANNOTATED_DIR} {QC_DIR} {SNPEFF_DATA_DIR} {SNAKEMAKE_DIR}
        touch {output.marker}
        """

rule download_reference:
    input:
        marker = rules.create_dirs.output.marker
    output:
        reference_fasta = f"{RAW_DIR}/reference.fasta"
    shell:
        """
        echo Downloading reference genome...
        efetch -db nucleotide -id {REF_ID} -format fasta > {RAW_DIR}/reference.fasta
        echo Downloaded reference genome!
        """

rule download_sra:
    input:
        marker = rules.create_dirs.output.marker
    output:
        sequence_sra = f"{RAW_DIR}/{SRA}/{SRA}.sra"
    shell:
        """
        echo Downloading sequencing data...
        prefetch {SRA} -O {RAW_DIR}
        echo Downloaded sequencing data!
        """

rule extract_sequence:
    input:
        marker = rules.create_dirs.output.marker,
        sequence_sra = rules.download_sra.output.sequence_sra
    output:
        sequence_fastq = f"{RAW_DIR}/{SRA}.fastq"
    shell:
        """
        echo Extracting sequencing data...
        fastq-dump -X 10000 {RAW_DIR}/{SRA}/{SRA}.sra -O {RAW_DIR}
        echo Extracted sequencing data!
        """


rule check_files:
	input:
	    fasta = f"{RAW_DIR}/reference.fasta",
	    fastq = f"{RAW_DIR}/{SRA}.fastq"
	output:
	    touch(f"{RESULTS_FOLDER}/.checked_files")
	shell:
	    """
	    if [ ! -s {input.fasta} ]; then
	        echo "Error: Reference genome file is missing or empty." >&2
	        exit 1
	    fi
	
	    if [ ! -s {input.fastq} ]; then
            echo "Error: FASTQ file is missing or empty." >&2
	        exit 1
	    fi
	
	    touch {output}
	    """

rule fastqc_raw:
    input:
        marker = rules.create_dirs.output.marker,
        sequence_fastq = rules.extract_sequence.output.sequence_fastq
    output:
        html = f"{QC_DIR}/{SRA}_fastqc.html",
        zip = f"{QC_DIR}/{SRA}_fastqc.zip"
    shell:
        """
        echo Running FastQC on raw reads...
        fastqc -o {QC_DIR} {input.sequence_fastq}
        echo FastQC complete!
        """

rule samtools_faidx:
    input:
        marker = rules.create_dirs.output.marker,
        reference_fasta = rules.download_reference.output.reference_fasta
    output:
        fai = f"{RAW_DIR}/reference.fasta.fai"
    shell:
        """
        echo Indexing reference genome with samtools...
        samtools faidx {input.reference_fasta}
        echo Reference genome indexed!
        """

rule bwa_index:
    input:
        marker = rules.create_dirs.output.marker,
        reference_fasta = rules.download_reference.output.reference_fasta
    output:
        [
            f"{RAW_DIR}/reference.fasta.amb",
            f"{RAW_DIR}/reference.fasta.ann",
            f"{RAW_DIR}/reference.fasta.bwt",
            f"{RAW_DIR}/reference.fasta.pac",
            f"{RAW_DIR}/reference.fasta.sa"
        ]
    shell:
        """
        echo Building BWA index...
        bwa index {input.reference_fasta}
        echo BWA index built!
        """

rule gatk_dict:
    input:
        marker = rules.create_dirs.output.marker,
        reference_fasta = rules.download_reference.output.reference_fasta
    output:
        dict = f"{RAW_DIR}/reference.dict"
    shell:
        """
        echo Creating FASTA dictionary using GATK...
        gatk CreateSequenceDictionary -R {input.reference_fasta} -O {output.dict}
        echo FASTA dictionary created!
        """

rule bwa_mem:
    input:
        marker = rules.create_dirs.output.marker,
        reference_fasta = rules.download_reference.output.reference_fasta,
        fastq = rules.extract_sequence.output.sequence_fastq
    output:
        sam = f"{ALIGNED_DIR}/aligned.sam"
    shell:
        """
        echo Aligning reads with read groups...
        bwa mem -R '@RG\\tID:1\\tLB:lib1\\tPL:illumina\\tPU:unit1\\tSM:sample1' {input.reference_fasta} {input.fastq} > {output.sam}
        echo Aligned reads!
        """

rule sam_to_sorted_bam:
    input:
        marker = rules.create_dirs.output.marker,
        sam = rules.bwa_mem.output.sam
    output:
        bam = f"{ALIGNED_DIR}/aligned.sorted.bam"
    shell:
        """
        echo Converting SAM to sorted BAM...
        samtools view -b {input.sam} | samtools sort -o {output.bam}
        echo Conversion complete!
        """

rule validate_bam:
    input:
        marker = rules.create_dirs.output.marker,
        bam = rules.sam_to_sorted_bam.output.bam
    output:
        touch(f"{ALIGNED_DIR}/validated_bam.txt")
    shell:
        """
        echo Validating BAM file...
        gatk ValidateSamFile -I {input.bam} -MODE SUMMARY
        echo Validation complete!
        touch {output}
        """

rule mark_duplicates:
    input:
        marker = rules.create_dirs.output.marker,
        bam = rules.sam_to_sorted_bam.output.bam
    output:
        dedup_bam = f"{ALIGNED_DIR}/dedup.bam",
        metrics = f"{ALIGNED_DIR}/dup_metrics.txt"
    shell:
        """
        echo Marking duplicates...
        gatk MarkDuplicates -I {input.bam} -O {output.dedup_bam} -M {output.metrics}
        echo Duplicates marked!
        """

rule index_dedup_bam:
    input:
        marker = rules.create_dirs.output.marker,
        dedup_bam = rules.mark_duplicates.output.dedup_bam
    output:
        bai = f"{ALIGNED_DIR}/dedup.bam.bai"
    shell:
        """
        echo Indexing deduplicated BAM file...
        samtools index {input.dedup_bam}
        echo BAM indexing complete!
        """

rule haplotype_caller:
    input:
        marker = rules.create_dirs.output.marker,
        reference_fasta = rules.download_reference.output.reference_fasta,
        dedup_bam = rules.mark_duplicates.output.dedup_bam
    output:
        vcf = f"{VARIANT_DIR}/raw_variants.vcf"
    shell:
        """
        echo Calling variants...
        gatk HaplotypeCaller -R {input.reference_fasta} -I {input.dedup_bam} -O {output.vcf}
        echo Called variants!
        """

rule variant_filtration:
    input:
        marker = rules.create_dirs.output.marker,
        reference_fasta = rules.download_reference.output.reference_fasta,
        raw_vcf = rules.haplotype_caller.output.vcf
    output:
        filtered_vcf = f"{VARIANT_DIR}/filtered_variants.vcf"
    shell:
        """
        echo Filtering variants...
        gatk VariantFiltration -R {input.reference_fasta} -V {input.raw_vcf} -O {output.filtered_vcf} --filter-expression "QD < 2.0 || FS > 60.0" --filter-name FILTER
        echo Variant filtration complete!
        """
