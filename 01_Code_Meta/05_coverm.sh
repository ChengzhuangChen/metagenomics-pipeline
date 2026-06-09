#!/bin/bash
# ============================================================
# MAG Abundance Profiling with CoverM
# For Nature article supplementary materials
#
# Method: MAG abundance was calculated using CoverM with RPKM
# and TPM normalization methods. Mapping was performed with
# Bowtie2 and reads with <95% identity or <10% coverage were
# excluded to ensure accurate quantification.
#
# Version: 1.0
# Dependencies: CoverM, Bowtie2
#
# Usage:
#   bash 05_coverm.sh -m <mag_dir> -f <fastq_dir> -o <output_dir>
# ============================================================

set -e

# ================== CONFIGURATION ==================
THREADS=32
MIN_IDENTITY=95                  # Minimum read percent identity (%)
MIN_COVERAGE=10                  # Minimum covered fraction (%)
EXTENSION="fa"                   # MAG file extension
# Normalization methods: rpkm, tpm, counts, coverage, CPM
NORMALIZATION="rpkm"
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
    local deps=("coverm")
    for dep in "${deps[@]}"; do
        check_software ${dep}
    done
    log "All dependencies found."
}

# ================== STEP 1: Prepare Sample List ==================
prepare_samples() {
    local fastq_dir=$1
    local output_dir=$2

    log "[Step 1/3] Preparing sample list..."

    mkdir -p ${output_dir}

    local coupled_reads=""
    local sample_count=0

    for r1 in ${fastq_dir}/*_R1.fq.gz; do
        local sample_name=$(basename ${r1} _R1.fq.gz)
        local r2="${fastq_dir}/${sample_name}_R2.fq.gz"

        if [[ -f "${r2}" ]]; then
            coupled_reads="${coupled_reads} ${r1} ${r2}"
            sample_count=$((sample_count + 1))
            log "  -> Found sample: ${sample_name}"
        fi
    done

    if [ ${sample_count} -eq 0 ]; then
        error_exit "No valid paired-end read files found in ${fastq_dir}"
    fi

    log "  -> Total samples: ${sample_count} pairs"
    echo "${coupled_reads}"
    echo "${sample_count}"
}

# ================== STEP 2: Run CoverM ==================
run_coverm() {
    local mag_dir=$1
    local coupled_reads=$2
    local output_dir=$3

    log "[Step 2/3] Running CoverM genome abundance..."

    local output_file="${output_dir}/abundance_matrix.tsv"

    coverm genome \
        --genome-fasta-directory ${mag_dir} \
        --genome-fasta-extension ${EXTENSION} \
        -c ${coupled_reads} \
        -m ${NORMALIZATION} \
        -t ${THREADS} \
        --min-read-percent-identity ${MIN_IDENTITY} \
        --min-covered-fraction ${MIN_COVERAGE} \
        -o ${output_file}

    log "  -> Abundance matrix: ${output_file}"
    echo "${output_file}"
}

# ================== STEP 3: Generate Additional Matrices ==================
generate_matrices() {
    local mag_dir=$1
    local coupled_reads=$2
    local output_dir=$3

    log "[Step 3/3] Generating additional normalization matrices..."

    local methods=("tpm" "counts" "coverage" "cpm")

    for method in "${methods[@]}"; do
        local output_file="${output_dir}/abundance_matrix_${method}.tsv"
        log "  -> Generating ${method} matrix..."

        coverm genome \
            --genome-fasta-directory ${mag_dir} \
            --genome-fasta-extension ${EXTENSION} \
            -c ${coupled_reads} \
            -m ${method} \
            -t ${THREADS} \
            --min-read-percent-identity ${MIN_IDENTITY} \
            --min-covered-fraction ${MIN_COVERAGE} \
            -o ${output_file}
    done

    log "  -> Additional matrices generated"
}

# ================== STEP 4: Summary Report ==================
generate_summary() {
    local output_dir=$1
    local sample_count=$2

    log "[Step 4/4] Generating summary report..."

    local summary_file="${output_dir}/summary_report.txt"

    {
        echo "=========================================="
        echo "CoverM MAG Abundance Profiling Summary"
        echo "=========================================="
        echo ""
        echo "Parameters:"
        echo "  Threads: ${THREADS}"
        echo "  Min identity: ${MIN_IDENTITY}%"
        echo "  Min coverage: ${MIN_COVERAGE}%"
        echo "  Normalization: ${NORMALIZATION}"
        echo "  File extension: ${EXTENSION}"
        echo "  Total samples: ${sample_count}"
        echo ""
        echo "Output files:"
        echo "  - ${output_dir}/abundance_matrix.tsv (${NORMALIZATION})"
        echo "  - ${output_dir}/abundance_matrix_tpm.tsv"
        echo "  - ${output_dir}/abundance_matrix_counts.tsv"
        echo "  - ${output_dir}/abundance_matrix_coverage.tsv"
        echo "  - ${output_dir}/abundance_matrix_cpm.tsv"
        echo ""
        echo "Matrix format:"
        echo "  Column 1: Genome (MAG identifier)"
        echo "  Column 2+: Sample abundances"
        echo "=========================================="
    } > ${summary_file}

    log "  -> Summary report: ${summary_file}"

    # Count MAGs
    if [ -f "${output_dir}/abundance_matrix.tsv" ]; then
        local mag_count=$(tail -n +1 ${output_dir}/abundance_matrix.tsv | wc -l)
        mag_count=$((mag_count - 1))  # Subtract header
        log "  === Abundance Summary ==="
        log "  -> MAGs in matrix: ${mag_count}"
        log "  -> Samples: ${sample_count}"
    fi
}

# ================== MAIN ENTRY POINT ==================

print_usage() {
    cat << EOF
Usage: bash 05_coverm.sh -m <mag_dir> -f <fastq_dir> -o <output_dir>

Description:
    Calculate MAG abundance across samples using CoverM.
    Generates RPKM, TPM, counts, coverage, and CPM matrices.

Required Arguments:
  -m, --mag-dir     Directory containing MAG FASTA files
  -f, --fastq-dir   Directory containing FASTQ files
  -o, --output      Output directory

Optional Arguments:
  -t, --threads     Number of threads (default: 32)
  -i, --min-id      Minimum read identity (default: 95)
  -c, --min-cov     Minimum covered fraction (default: 10)
  -e, --ext         MAG file extension (default: fa)
  -n, --norm        Normalization method: rpkm, tpm, counts, coverage, cpm (default: rpkm)

Example:
  bash 05_coverm.sh \\
      -m results/mag_quality/drep/dereplicated_genomes \\
      -f results/clean_data \\
      -o results/abundance

  # With custom parameters
  bash 05_coverm.sh -m bins/ -f fastq/ -o output -t 64 -i 97 -n tpm

Input:
  - mag_dir: Directory with MAGs as FASTA files (.fa or .fasta)
  - fastq_dir: Directory with paired-end FASTQ files (*_R1.fq.gz, *_R2.fq.gz)

Output:
  - abundance_matrix.tsv: Primary abundance matrix (default: RPKM)
  - abundance_matrix_tpm.tsv: TPM normalized matrix
  - abundance_matrix_counts.tsv: Raw read counts
  - abundance_matrix_coverage.tsv: Coverage values
  - abundance_matrix_cpm.tsv: Counts per million

Normalization Methods:
  - rpkm: Reads Per Kilobase per Million mapped reads
  - tpm: Transcripts Per Million
  - counts: Raw mapped read counts
  - coverage: Mean covered bases per position
  - cpm: Counts per million reads
EOF
}

main() {
    local mag_dir=""
    local fastq_dir=""
    local output_dir=""

    # Parse arguments
    while [[ $# -gt 0 ]]; do
        case $1 in
            -m|--mag-dir) mag_dir="$2"; shift 2 ;;
            -f|--fastq-dir) fastq_dir="$2"; shift 2 ;;
            -o|--output) output_dir="$2"; shift 2 ;;
            -t|--threads) THREADS="$2"; shift 2 ;;
            -i|--min-id) MIN_IDENTITY="$2"; shift 2 ;;
            -c|--min-cov) MIN_COVERAGE="$2"; shift 2 ;;
            -e|--ext) EXTENSION="$2"; shift 2 ;;
            -n|--norm) NORMALIZATION="$2"; shift 2 ;;
            -h|--help) print_usage; exit 0 ;;
            *) echo "Unknown option: $1"; print_usage; exit 1 ;;
        esac
    done

    # Validate required arguments
    if [ -z "$mag_dir" ] || [ -z "$fastq_dir" ] || [ -z "$output_dir" ]; then
        echo "Error: Missing required arguments!"
        print_usage
        exit 1
    fi

    # Check dependencies
    check_dependencies

    log "=========================================="
    log "CoverM MAG Abundance Profiling"
    log "=========================================="
    log "MAG directory: ${mag_dir}"
    log "FASTQ directory: ${fastq_dir}"
    log "Output directory: ${output_dir}"
    log "Threads: ${THREADS}"
    log "Min identity: ${MIN_IDENTITY}%"
    log "Min coverage: ${MIN_COVERAGE}%"
    log "Normalization: ${NORMALIZATION}"
    log "=========================================="

    # Step 1: Prepare sample list
    local result=$(prepare_samples ${fastq_dir} ${output_dir})
    local coupled_reads=$(echo "$result" | head -n1)
    local sample_count=$(echo "$result" | tail -n1)

    # Step 2: Run CoverM with primary normalization
    local output_file=$(run_coverm ${mag_dir} "${coupled_reads}" ${output_dir})

    # Step 3: Generate additional matrices
    generate_matrices ${mag_dir} "${coupled_reads}" ${output_dir}

    # Step 4: Generate summary
    generate_summary ${output_dir} ${sample_count}

    log "=========================================="
    log "Pipeline completed successfully!"
    log "=========================================="
    log "Output: ${output_dir}/"
}

main "$@"