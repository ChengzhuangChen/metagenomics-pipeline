#!/bin/bash
# ============================================================
# Metagenomic QC and Co-assembly Pipeline
# For Nature article supplementary materials
#
# Version: 1.0
# Dependencies: Trimmomatic, FastQC, MultiQC, metaSPAdes,
#               MEGAHIT, Bowtie2, samtools, QUAST, seqkit
#
# Usage:
#   Single sample:  bash 01_QC_Co_assembly.sh sample_name R1.fastq.gz R2.fastq.gz
#   Multiple samples: bash 01_QC_Co_assembly.sh samples.txt
#
# Input samples.txt format (one line per sample):
#   sample_name,reads_R1.fastq.gz,reads_R2.fastq.gz
# ============================================================

set -e  # Exit on error

# ================== CONFIGURATION ==================
# -------------------- EDIT THESE VALUES --------------------
THREADS=32                      # Number of threads
MIN_CONTIG_LENGTH=1000          # Minimum contig length (bp)
ADAPTERS="TruSeq3-PE.fa"        # Trimmomatic adapter file path
WORK_DIR="results"              # Output directory
# ---------------------------------------------------

# ================== HELPER FUNCTIONS ==================

log() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $1"
}

error_exit() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] ERROR: $1" >&2
    exit 1
}

check_software() {
    local software=$1
    if ! command -v ${software} &> /dev/null; then
        error_exit "${software} not found. Please install it first."
    fi
}

check_dependencies() {
    log "Checking dependencies..."
    local deps=("trimmomatic" "fastqc" "multiqc" "spades.py" "megahit" "bowtie2" "samtools" "quast.py" "seqkit")
    for dep in "${deps[@]}"; do
        check_software ${dep}
    done
    log "All dependencies found."
}

# ================== PIPELINE STEPS ==================

# Step 1: Quality Control with Trimmomatic
# -----------------------------------------
run_trimmomatic() {
    local r1=$1
    local r2=$2
    local output_prefix=$3
    local sample_dir=$4

    log "[Step 1/5] Quality filtering with Trimmomatic v0.39..."

    mkdir -p ${sample_dir}/clean

    trimmomatic PE \
        -threads ${THREADS} \
        -phred33 \
        ${r1} ${r2} \
        ${output_prefix}_R1_paired.fastq.gz \
        ${output_prefix}_R1_unpaired.fastq.gz \
        ${output_prefix}_R2_paired.fastq.gz \
        ${output_prefix}_R2_unpaired.fastq.gz \
        ILLUMINACLIP:${ADAPTERS}:2:30:10:2:keepBothReads \
        LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:65

    log "  -> Clean reads: ${output_prefix}_R*_paired.fastq.gz"
}

# Step 2: Quality Assessment with FastQC & MultiQC
# -------------------------------------------------
run_qc_reports() {
    local reads_dir=$1
    local sample_dir=$2
    local sample_name=$3

    log "[Step 2/5] Quality assessment with FastQC v0.12.1 and MultiQC v1.30..."

    mkdir -p ${sample_dir}/qc/fastqc
    fastqc ${reads_dir}/*_paired.fastq.gz -o ${sample_dir}/qc/fastqc -t ${THREADS}

    mkdir -p ${sample_dir}/qc/multiqc
    multiqc ${sample_dir}/qc/fastqc -o ${sample_dir}/qc/multiqc -n ${sample_name}_qc_report

    log "  -> QC report: ${sample_dir}/qc/multiqc/${sample_name}_qc_report.html"
}

# Step 3: Two-step Hybrid Assembly
# ---------------------------------
run_hybrid_assembly() {
    local r1=$1
    local r2=$2
    local sample_dir=$3

    log "[Step 3/5] Two-step hybrid assembly (metaSPAdes v4.2.0 -> MEGAHIT v1.2.9)..."

    # 3a: Primary assembly with metaSPAdes
    # ------------------------------------
    log "  [3a] Primary assembly with metaSPAdes..."
    mkdir -p ${sample_dir}/metaspades

    spades.py \
        --meta \
        --careful \
        -1 ${r1} \
        -2 ${r2} \
        -o ${sample_dir}/metaspades \
        -t ${THREADS} \
        -m 500

    local primary_contigs="${sample_dir}/metaspades/scaffolds.fasta"
    log "  -> Primary contigs: ${primary_contigs}"

    # 3b: Extract unmapped reads
    # --------------------------
    log "  [3b] Extracting unmapped reads..."
    mkdir -p ${sample_dir}/unmapped
    local index_dir="${sample_dir}/index"
    mkdir -p ${index_dir}

    bowtie2-build ${primary_contigs} ${index_dir}/contigs

    bowtie2 -p ${THREADS} \
        -x ${index_dir}/contigs \
        -1 ${r1} \
        -2 ${r2} \
        -S ${sample_dir}/unmapped/aligned.sam

    samtools view -bS -@ ${THREADS} ${sample_dir}/unmapped/aligned.sam | \
        samtools sort -@ ${THREADS} -o ${sample_dir}/unmapped/aligned_sorted.bam

    samtools view -b -f 12 -F 256 ${sample_dir}/unmapped/aligned_sorted.bam > \
        ${sample_dir}/unmapped/both_unmapped.bam

    samtools fastq -1 ${sample_dir}/unmapped/unmapped_R1.fastq.gz \
                   -2 ${sample_dir}/unmapped/unmapped_R2.fastq.gz \
                   ${sample_dir}/unmapped/both_unmapped.bam

    log "  -> Unmapped reads: ${sample_dir}/unmapped/unmapped_R*.fastq.gz"

    # 3c: Secondary assembly with MEGAHIT
    # ------------------------------------
    log "  [3c] Secondary assembly with MEGAHIT..."
    mkdir -p ${sample_dir}/megahit

    megahit \
        -1 ${sample_dir}/unmapped/unmapped_R1.fastq.gz \
        -2 ${sample_dir}/unmapped/unmapped_R2.fastq.gz \
        -o ${sample_dir}/megahit \
        -t ${THREADS} \
        --min-count 2 \
        --k-min 27 \
        --k-max 127 \
        --k-step 10

    local secondary_contigs="${sample_dir}/megahit/final.contigs.fa"
    log "  -> Secondary contigs: ${secondary_contigs}"
}

# Step 4: Merge and Filter Contigs
# ---------------------------------
merge_and_filter() {
    local metaspades_contigs=$1
    local megahit_contigs=$2
    local sample_dir=$3
    local min_len=$4

    log "[Step 4/5] Merging assemblies and filtering (≥${min_len} bp)..."

    mkdir -p ${sample_dir}/final

    # Merge contigs from both assemblies
    cat ${metaspades_contigs} ${megahit_contigs} > ${sample_dir}/final/merged_contigs.fa

    # Filter by length using seqkit
    seqkit seq -m ${min_len} \
        ${sample_dir}/final/merged_contigs.fa \
        -o ${sample_dir}/final/contigs_${min_len}bp.fa

    # Statistics
    local total_contigs=$(grep -c "^>" ${sample_dir}/final/merged_contigs.fa)
    local filtered_contigs=$(grep -c "^>" ${sample_dir}/final/contigs_${min_len}bp.fa)
    local total_bp=$(seqkit stats -b ${sample_dir}/final/merged_contigs.fa | tail -1 | awk '{print $5}')
    local filtered_bp=$(seqkit stats -b ${sample_dir}/final/contigs_${min_len}bp.fa | tail -1 | awk '{print $5}')

    log "  -> Merged contigs: ${total_contigs} (${total_bp} bp)"
    log "  -> Filtered contigs (≥${min_len} bp): ${filtered_contigs} (${filtered_bp} bp)"
    log "  -> Final contigs: ${sample_dir}/final/contigs_${min_len}bp.fa"
}

# Step 5: Assembly Quality Assessment with QUAST
# ------------------------------------------------
run_quast() {
    local contigs=$1
    local sample_dir=$2

    log "[Step 5/5] Quality assessment with QUAST v5.3.0..."

    mkdir -p ${sample_dir}/quast

    quast.py ${contigs} \
        -o ${sample_dir}/quast \
        -t ${THREADS} \
        -l "hybrid_assembly"

    log "  -> QUAST report: ${sample_dir}/quast/report.html"
}

# ================== SINGLE SAMPLE PROCESSING ==================

process_single_sample() {
    local sample_name=$1
    local reads_r1=$2
    local reads_r2=$3

    log "=========================================="
    log "Processing sample: ${sample_name}"
    log "=========================================="
    log "R1: ${reads_r1}"
    log "R2: ${reads_r2}"

    local sample_dir="${WORK_DIR}/${sample_name}"

    # Run pipeline steps
    run_trimmomatic ${reads_r1} ${reads_r2} ${sample_dir}/clean/${sample_name} ${sample_dir}
    run_qc_reports ${sample_dir}/clean ${sample_dir} ${sample_name}
    run_hybrid_assembly \
        ${sample_dir}/clean/${sample_name}_R1_paired.fastq.gz \
        ${sample_dir}/clean/${sample_name}_R2_paired.fastq.gz \
        ${sample_dir}
    merge_and_filter \
        ${sample_dir}/metaspades/scaffolds.fasta \
        ${sample_dir}/megahit/final.contigs.fa \
        ${sample_dir} ${MIN_CONTIG_LENGTH}
    run_quast ${sample_dir}/final/contigs_${MIN_CONTIG_LENGTH}bp.fa ${sample_dir}

    log "=========================================="
    log "Sample ${sample_name} complete!"
    log "=========================================="
}

# ================== CO-ASSEMBLY PROCESSING ==================

process_coassembly() {
    local sample_list=$1

    log "=========================================="
    log "Processing co-assembly from: ${sample_list}"
    log "=========================================="

    local coassembly_dir="${WORK_DIR}/coassembly"
    mkdir -p ${coassembly_dir}/clean

    # Combine all cleaned reads for co-assembly
    log "[Co-assembly Step 1/4] Combining reads from all samples..."

    local first_sample=true
    local r1_files=""
    local r2_files=""

    while IFS=, read -r sample r1 r2; do
        # Skip header if present
        [[ "$sample" == "sample_name" ]] && continue
        [[ -z "$sample" ]] && continue

        if [ "$first_sample" = true ]; then
            r1_files="${r1}"
            r2_files="${r2}"
            first_sample=false
        else
            r1_files="${r1_files} ${r1}"
            r2_files="${r2_files} ${r2}"
        fi

        log "  + ${sample}: ${r1}, ${r2}"
    done < "${sample_list}"

    # Step 1: Quality control for all samples
    log "[Co-assembly Step 2/4] Quality control for all samples..."

    local r1_paired=""
    local r2_paired=""

    while IFS=, read -r sample r1 r2; do
        [[ "$sample" == "sample_name" ]] && continue
        [[ -z "$sample" ]] && continue

        local clean_prefix="${coassembly_dir}/clean/${sample}"

        trimmomatic PE \
            -threads ${THREADS} \
            -phred33 \
            ${r1} ${r2} \
            ${clean_prefix}_R1_paired.fastq.gz \
            ${clean_prefix}_R1_unpaired.fastq.gz \
            ${clean_prefix}_R2_paired.fastq.gz \
            ${clean_prefix}_R2_unpaired.fastq.gz \
            ILLUMINACLIP:${ADAPTERS}:2:30:10:2:keepBothReads \
            LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:65

        r1_paired="${r1_paired} ${clean_prefix}_R1_paired.fastq.gz"
        r2_paired="${r2_paired} ${clean_prefix}_R2_paired.fastq.gz"

    done < "${sample_list}"

    # Step 2: Run FastQC/MultiQC
    log "[Co-assembly Step 3/4] Quality assessment..."
    mkdir -p ${coassembly_dir}/qc/fastqc
    fastqc ${coassembly_dir}/clean/*_paired.fastq.gz -o ${coassembly_dir}/qc/fastqc -t ${THREADS}

    mkdir -p ${coassembly_dir}/qc/multiqc
    multiqc ${coassembly_dir}/qc/fastqc -o ${coassembly_dir}/qc/multiqc -n coassembly_qc_report

    # Step 3: Co-assembly with MEGAHIT
    log "[Co-assembly Step 4/4] Co-assembly with MEGAHIT..."

    megahit \
        -1 $(echo $r1_paired | tr ' ' ',') \
        -2 $(echo $r2_paired | tr ' ' ',') \
        -o ${coassembly_dir}/megahit \
        -t ${THREADS} \
        --min-count 2 \
        --k-min 27 \
        --k-max 127 \
        --k-step 10

    # Step 4: Filter and assess
    log "[Post-assembly] Filtering and quality assessment..."

    mkdir -p ${coassembly_dir}/final

    seqkit seq -m ${MIN_CONTIG_LENGTH} \
        ${coassembly_dir}/megahit/final.contigs.fa \
        -o ${coassembly_dir}/final/contigs_${MIN_CONTIG_LENGTH}bp.fa

    mkdir -p ${coassembly_dir}/quast
    quast.py ${coassembly_dir}/final/contigs_${MIN_CONTIG_LENGTH}bp.fa \
        -o ${coassembly_dir}/quast -t ${THREADS}

    log "=========================================="
    log "Co-assembly complete!"
    log "=========================================="
    log "  -> Final contigs: ${coassembly_dir}/final/contigs_${MIN_CONTIG_LENGTH}bp.fa"
    log "  -> QUAST report: ${coassembly_dir}/quast/report.html"
}

# ================== MAIN ENTRY POINT ==================

print_usage() {
    cat << EOF
Usage: bash 01_QC_Co_assembly.sh <mode> [options]

Modes:
  1. Single sample assembly:
     bash 01_QC_Co_assembly.sh -s <sample_name> <R1.fastq.gz> <R2.fastq.gz>

  2. Multiple samples (batch):
     bash 01_QC_Co_assembly.sh -m <samples.txt>

  3. Co-assembly (combine all samples):
     bash 01_QC_Co_assembly.sh -c <samples.txt>

Examples:
  # Single sample
  bash 01_QC_Co_assembly.sh -s sample1 raw/sample1_R1.fastq.gz raw/sample1_R2.fastq.gz

  # Multiple samples
  bash 01_QC_Co_assembly.sh -m samples.txt

  # Co-assembly (combine all samples)
  bash 01_QC_Co_assembly.sh -c samples.txt

samples.txt format (comma-separated, one sample per line):
  sample_name,reads_R1.fastq.gz,reads_R2.fastq.gz

Configuration (edit in script header):
  THREADS           - Number of threads (default: 32)
  MIN_CONTIG_LENGTH - Minimum contig length in bp (default: 1000)
  ADAPTERS          - Path to Trimmomatic adapter file
  WORK_DIR          - Output directory (default: results)
EOF
}

main() {
    # Check dependencies
    check_dependencies

    if [ $# -eq 0 ]; then
        print_usage
        exit 1
    fi

    local mode=$1

    case ${mode} in
        -s|--single)
            # Single sample mode
            if [ $# -ne 4 ]; then
                echo "Error: Single sample mode requires <sample_name> <R1> <R2>"
                print_usage
                exit 1
            fi
            process_single_sample $2 $3 $4
            ;;
        -m|--multiple)
            # Multiple samples mode (individual assembly)
            if [ $# -ne 2 ]; then
                echo "Error: Multiple samples mode requires <samples.txt>"
                print_usage
                exit 1
            fi
            while IFS=, read -r sample r1 r2; do
                [[ "$sample" == "sample_name" ]] && continue
                [[ -z "$sample" ]] && continue
                process_single_sample ${sample} ${r1} ${r2}
            done < $2
            ;;
        -c|--coassembly)
            # Co-assembly mode
            if [ $# -ne 2 ]; then
                echo "Error: Co-assembly mode requires <samples.txt>"
                print_usage
                exit 1
            fi
            process_coassembly $2
            ;;
        -h|--help)
            print_usage
            ;;
        *)
            echo "Error: Unknown mode '${mode}'"
            print_usage
            exit 1
            ;;
    esac

    log "=========================================="
    log "Pipeline completed successfully!"
    log "=========================================="
}

main "$@"