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

rule all:
    input: 
        f"{SNAKEMAKE_DIR}/.dirs_created",
        f"{RAW_DIR}/reference.fasta",
        f"{RAW_DIR}/{SRA}/{SRA}.sra",
        f"{RAW_DIR}/{SRA}.fastq",
        f"{RESULTS_FOLDER}/.checked_files",
        f"{QC_DIR}/{SRA}_fastqc.html",
        f"{QC_DIR}/{SRA}_fastqc.zip",
        f"{RAW_DIR}/reference.fasta.fai"
        
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