#!/bin/bash
# ============================================================
# ARG Annotation Pipeline using OAP (SARG)
# For Nature article supplementary materials
#
# Method: Antibiotic resistance genes (ARGs) were identified using
# OAP (Online Antibiotic Resistance platform) against the SARG
# (Structured Antibiotic Resistance Gene) database. The pipeline
# consists of two stages: Stage 1 identifies ARG-like reads using
# DIAMOND alignment, and Stage 2 assembles ARG fragments and
# maps reads back for accurate quantification.
#
# Version: 1.0
# Dependencies: args_oap (OAP pipeline)
#
# Usage:
#   bash 02_ARGs_OAP.sh -i <fastq_dir> -o <output_dir> -d <sarg_db>
# ============================================================

set -e

# ================== CONFIGURATION ==================
THREADS=32
FASTQ_EXT="fq.gz"                # FASTQ file extension
SARG_DB=""                       # SARG database path (required)
# ====================================================

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
    local deps=("args_oap")
    for dep in "${deps[@]}"; do
        check_software ${dep}
    done
    log "All dependencies found."
}

# ================== STEP 1: OAP Stage One ==================
run_oap_stage1() {
    local fastq_dir=$1
    local output_dir=$2

    log "[Step 1/2] OAP Stage 1: Identifying ARG-like reads..."

    local stage1_dir="${output_dir}/stage_one"
    mkdir -p ${stage1_dir}

    # Run OAP stage 1
    # -i: Input FASTQ directory
    # -o: Output directory
    # -f: FASTQ file extension
    # -t: Number of threads
    # --database: SARG database path
    args_oap stage_one \
        -i ${fastq_dir} \
        -o ${stage1_dir} \
        -f ${FASTQ_EXT} \
        -t ${THREADS} \
        --database ${SARG_DB}

    log "  -> Stage 1 complete!"
    log "  -> Results: ${stage1_dir}/"

    echo "${stage1_dir}"
}

# ================== STEP 2: OAP Stage Two ==================
run_oap_stage2() {
    local stage1_output=$1
    local output_dir=$2

    log "[Step 2/2] OAP Stage 2: Assembling and quantifying ARGs..."

    local stage2_dir="${output_dir}/stage_two"
    mkdir -p ${stage2_dir}

    # Run OAP stage 2
    # -i: Input from stage 1
    # -t: Number of threads
    # -o: Output directory
    args_oap stage_two \
        -i ${stage1_output} \
        -t ${THREADS} \
        -o ${stage2_dir}

    log "  -> Stage 2 complete!"
    log "  -> Results: ${stage2_dir}/"

    echo "${stage2_dir}"
}

# ================== STEP 3: Generate Summary ==================
generate_summary() {
    local output_dir=$1

    log "[Step 3/3] Generating summary report..."

    local summary_file="${output_dir}/summary_report.txt"

    {
        echo "=========================================="
        echo "ARG Annotation Pipeline (OAP/SARG)"
        echo "=========================================="
        echo ""
        echo "Parameters:"
        echo "  Threads: ${THREADS}"
        echo "  FASTQ extension: ${FASTQ_EXT}"
        echo "  SARG database: ${SARG_DB}"
        echo ""
        echo "Output Directories:"
        echo "  - ${output_dir}/stage_one/ (ARG-like reads)"
        echo "  - ${output_dir}/stage_two/ (Assembled ARGs)"
        echo ""
        echo "Stage 1 Output Files:"
        echo "  - *_ORF (nucleotide sequences)"
        echo "  - *_diamond (DIAMOND alignment results)"
        echo "  - *_extended_diamond (Extended alignment)"
        echo "  - *normalized (Normalized counts)"
        echo ""
        echo "Stage 2 Output Files:"
        echo "  - * ARG abundance by subtype.xls"
        echo "  - * ARG abundance by type.xls"
        echo "  - * ARG abundance by module.xls"
        echo "  - * ARG abundance by class.xls"
        echo "=========================================="
    } > ${summary_file}

    log "  -> Summary report: ${summary_file}"
}

# ================== MAIN ENTRY POINT ==================

print_usage() {
    cat << EOF
Usage: bash 02_ARGs_OAP.sh -i <fastq_dir> -o <output_dir> -d <sarg_db>

Description:
    Identify and quantify antibiotic resistance genes (ARGs) using
    OAP (Online Antibiotic Resistance Platform) against the SARG
    (Structured Antibiotic Resistance Gene) database.

    Stage 1: Identifies ARG-like reads via DIAMOND alignment
    Stage 2: Assembles ARG fragments and quantifies abundance

Required Arguments:
  -i, --input     Directory containing FASTQ files
  -o, --output    Output directory
  -d, --database  Path to SARG database (FASTA format)

Optional Arguments:
  -t, --threads   Number of threads (default: 32)
  -f, --ext       FASTQ file extension (default: fq.gz)

Example:
  bash 02_ARGs_OAP.sh \\
      -i results/clean_data \\
      -o results/ARGs \\
      -d /path/to/SARG database.fasta

  # With custom parameters
  bash 02_ARGs_OAP.sh -i fastq/ -o args_output -d sarg_db.fasta -t 64 -f fastq.gz

Input:
  - fastq_dir: Directory with paired-end FASTQ files
  - SARG database: Structured ARG Reference database (FASTA)

Output:
  Stage 1:
    - ORF files: Nucleotide sequences of ARG-like reads
    - Diamond results: DIAMOND alignment files
    - Normalized counts: Per-sample ARG abundances

  Stage 2:
    - ARG abundance tables at different levels:
      - Subtype level (*ARG_abundance_by_subtype.xls)
      - Type level (*ARG_abundance_by_type.xls)
      - Class level (*ARG_abundance_by_class.xls)

Installation:
  # Create conda environment
  conda create -n args_oap -c bioconda -c conda-forge args_oap

  conda activate args_oap

  # Download SARG database
  # The SARG database should be downloaded from:
  # https://smileassistant1.github.io/OAP/#!/download

Database:
  - SARG v3.2 (Short subdatabase): ~50 MB
  - Full SARG database: ~500 MB

  Database files should include:
  - *database.fasta: Reference ARG sequences
  - *database.map: Sequence to ARG mapping
EOF
}

main() {
    local fastq_dir=""
    local output_dir=""
    local sarg_db=""

    # Parse arguments
    while [[ $# -gt 0 ]]; do
        case $1 in
            -i|--input) fastq_dir="$2"; shift 2 ;;
            -o|--output) output_dir="$2"; shift 2 ;;
            -d|--database) sarg_db="$2"; shift 2 ;;
            -t|--threads) THREADS="$2"; shift 2 ;;
            -f|--ext) FASTQ_EXT="$2"; shift 2 ;;
            -h|--help) print_usage; exit 0 ;;
            *) echo "Unknown option: $1"; print_usage; exit 1 ;;
        esac
    done

    # Validate required arguments
    if [ -z "$fastq_dir" ] || [ -z "$output_dir" ] || [ -z "$sarg_db" ]; then
        echo "Error: Missing required arguments!"
        print_usage
        exit 1
    fi

    # Check dependencies
    check_dependencies

    log "=========================================="
    log "ARG Annotation Pipeline (OAP/SARG)"
    log "=========================================="
    log "FASTQ directory: ${fastq_dir}"
    log "Output directory: ${output_dir}"
    log "SARG database: ${sarg_db}"
    log "Threads: ${THREADS}"
    log "FASTQ extension: ${FASTQ_EXT}"
    log "=========================================="

    # Step 1: OAP Stage 1
    local stage1_output=$(run_oap_stage1 ${fastq_dir} ${output_dir})

    # Step 2: OAP Stage 2
    local stage2_output=$(run_oap_stage2 ${stage1_output} ${output_dir})

    # Step 3: Generate summary
    generate_summary ${output_dir}

    log "=========================================="
    log "Pipeline completed successfully!"
    log "=========================================="
    log "Results: ${output_dir}/"
    log "Abundance tables: ${stage2_output}/"
}

main "$@"