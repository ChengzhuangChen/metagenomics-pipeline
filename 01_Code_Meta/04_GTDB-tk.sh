#!/bin/bash
# ============================================================
# GTDB-Tk Taxonomy Classification Pipeline
# For Nature article supplementary materials
#
# Method: MAG taxonomy was determined using GTDB-Tk (v2.4.0)
# classify_wf workflow with the GTDB r223 database. Taxonomic
# assignment includes domain, phylum, class, order, family, genus,
# and species levels. Phylogenetic trees were inferred using
# FastTree for both bacterial and archaeal MAGs.
#
# Version: 1.0
# Dependencies: GTDB-Tk, FastTree, pplacer
#
# Usage:
#   bash 04_GTDB-tk.sh -i <bins_dir> -o <output_dir>
# ============================================================

set -e

# ================== CONFIGURATION ==================
THREADS=32
EXTENSION="fa"                    # Bin file extension
PREFIX="bin"                      # Output prefix
SKIP_ANI=false                    # Skip ANI screening (for large datasets)
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
    local deps=("gtdbtk")
    for dep in "${deps[@]}"; do
        check_software ${dep}
    done
    log "All dependencies found."
}

# ================== STEP 1: GTDB-Tk Classify ==================
run_gtdbtk_classify() {
    local genome_dir=$1
    local output_dir=$2

    log "[Step 1/2] Running GTDB-Tk classify_wf..."

    mkdir -p ${output_dir}

    # Build command
    local cmd="gtdbtk classify_wf \
        --genome_dir ${genome_dir} \
        --out_dir ${output_dir} \
        --extension ${EXTENSION} \
        --prefix ${PREFIX} \
        --cpus ${THREADS}"

    # Add optional flags
    if [ "${SKIP_ANI}" = true ]; then
        cmd="${cmd} --skip_ani_screen"
        log "  -> ANI screening skipped"
    fi

    # Run GTDB-Tk classify
    eval ${cmd}

    log "  -> GTDB-Tk classification complete"
}

# ================== STEP 2: Phylogenetic Tree Inference ==================
run_gtdbtk_infer() {
    local output_dir=$1

    log "[Step 2/2] Inferring phylogenetic trees..."

    # Decompress alignment files
    if [ -f "${output_dir}/align/bin.ar53.user_msa.fasta.gz" ]; then
        gunzip -k "${output_dir}/align/bin.ar53.user_msa.fasta.gz"
        local ar53_msa="${output_dir}/align/bin.ar53.user_msa.fasta"
    else
        local ar53_msa=""
    fi

    if [ -f "${output_dir}/align/bin.bac120.user_msa.fasta.gz" ]; then
        gunzip -k "${output_dir}/align/bin.bac120.user_msa.fasta.gz"
        local bac120_msa="${output_dir}/align/bin.bac120.user_msa.fasta"
    else
        local bac120_msa=""
    fi

    # Infer archaeal tree (if there are archaeal MAGs)
    if [ -n "${ar53_msa}" ] && [ -f "${ar53_msa}" ]; then
        log "  [2a] Inferring archaeal tree..."
        mkdir -p ${output_dir}/tree_ar53
        gtdbtk infer \
            --msa_file ${ar53_msa} \
            --out_dir ${output_dir}/tree_ar53 \
            --cpus ${THREADS} \
            --prefix ${PREFIX}
        log "  -> Archaeal tree: ${output_dir}/tree_ar53/${PREFIX}.ar53.user_msa.tree"
    fi

    # Infer bacterial tree
    if [ -n "${bac120_msa}" ] && [ -f "${bac120_msa}" ]; then
        log "  [2b] Inferring bacterial tree..."
        mkdir -p ${output_dir}/tree_bac120
        gtdbtk infer \
            --msa_file ${bac120_msa} \
            --out_dir ${output_dir}/tree_bac120 \
            --cpus ${THREADS} \
            --prefix ${PREFIX}
        log "  -> Bacterial tree: ${output_dir}/tree_bac120/${PREFIX}.bac120.user_msa.tree"
    fi

    log "  -> Phylogenetic tree inference complete"
}

# ================== STEP 3: Generate Summary ==================
generate_summary() {
    local output_dir=$1

    log "[Step 3/3] Generating summary..."

    local summary_file="${output_dir}/summary_report.txt"

    {
        echo "=========================================="
        echo "GTDB-Tk Taxonomy Classification Summary"
        echo "=========================================="
        echo ""
        echo "GTDB-Tk version: $(gtdbtk --version 2>/dev/null || echo 'unknown')"
        echo "Database: GTDB r223 (or later)"
        echo ""
        echo "Results:"
        echo "  - Classification: ${output_dir}/classify/${PREFIX}.bac120.summary.tsv"
        echo "  - Archaea: ${output_dir}/classify/${PREFIX}.ar53.summary.tsv"
        echo "  - Bacterial tree: ${output_dir}/tree_bac120/${PREFIX}.bac120.user_msa.tree"
        echo "  - Archaeal tree: ${output_dir}/tree_ar53/${PREFIX}.ar53.user_msa.tree"
        echo ""
        echo "Taxonomic levels:"
        echo "  - d__: Domain (Bacteria/Archaea)"
        echo "  - p__: Phylum"
        echo "  - c__: Class"
        echo "  - o__: Order"
        echo "  - f__: Family"
        echo "  - g__: Genus"
        echo "  - s__: Species"
        echo "=========================================="
    } > ${summary_file}

    log "  -> Summary report: ${summary_file}"

    # Count taxonomy distribution if summary file exists
    if [ -f "${output_dir}/classify/${PREFIX}.bac120.summary.tsv" ]; then
        local total=$(tail -n +2 ${output_dir}/classify/${PREFIX}.bac120.summary.tsv | wc -l)
        log "  === Taxonomy Summary ==="
        log "  -> Total bacterial MAGs classified: ${total}"
        log "  -> Top phyla:"
        tail -n +2 ${output_dir}/classify/${PREFIX}.bac120.summary.tsv | \
            cut -f2 | cut -d';' -f1 | sed 's/p__//' | sort | uniq -c | sort -rn | head -5
    fi
}

# ================== MAIN ENTRY POINT ==================

print_usage() {
    cat << EOF
Usage: bash 04_GTDB-tk.sh -i <bins_dir> -o <output_dir>

Description:
    Taxonomic classification of MAGs using GTDB-Tk with the classify_wf
    workflow. Generates phylogenetic trees for both bacterial and archaeal MAGs.

Required Arguments:
  -i, --input     Directory containing MAG bins (FASTA format)
  -o, --output    Output directory

Optional Arguments:
  -t, --threads   Number of threads (default: 32)
  -e, --ext       File extension for bins (default: fa)
  -p, --prefix    Output prefix (default: bin)
  -s, --skip-ani  Skip ANI screening (faster, less accurate)

Example:
  bash 04_GTDB-tk.sh \\
      -i results/mag_quality/drep/dereplicated_genomes \\
      -o results/taxonomy

  # With custom parameters
  bash 04_GTDB-tk.sh -i bins/ -o taxonomy -t 64 -p mygenomes -s

Input:
  - bins_dir: Directory containing MAG bins as FASTA files
  - Each bin should be a separate file with contigs belonging to one MAG
  - File extension should match: .fa or .fasta

Output:
  - classify/: GTDB-Tk classification results
    - *.bac120.summary.tsv: Bacterial classification
    - *.ar53.summary.tsv: Archaeal classification
    - *.bac120.markers_summary.tsv: Bacterial marker gene summary
    - *.ar53.markers_summary.tsv: Archaeal marker gene summary
  - align/: Multiple sequence alignments
  - tree_bac120/: Bacterial phylogenetic tree
  - tree_ar53/: Archaeal phylogenetic tree

Database Setup:
  # Download GTDB-Tk reference data
  download-gtdbtk-data.sh --full  # ~50 GB

  # Or use conda
  conda install -c bioconda gtdbtk
  gtdbtk download --db-version r223
EOF
}

main() {
    local genome_dir=""
    local output_dir=""

    # Parse arguments
    while [[ $# -gt 0 ]]; do
        case $1 in
            -i|--input) genome_dir="$2"; shift 2 ;;
            -o|--output) output_dir="$2"; shift 2 ;;
            -t|--threads) THREADS="$2"; shift 2 ;;
            -e|--ext) EXTENSION="$2"; shift 2 ;;
            -p|--prefix) PREFIX="$2"; shift 2 ;;
            -s|--skip-ani) SKIP_ANI=true; shift ;;
            -h|--help) print_usage; exit 0 ;;
            *) echo "Unknown option: $1"; print_usage; exit 1 ;;
        esac
    done

    # Validate required arguments
    if [ -z "$genome_dir" ] || [ -z "$output_dir" ]; then
        echo "Error: Missing required arguments!"
        print_usage
        exit 1
    fi

    # Check dependencies
    check_dependencies

    log "=========================================="
    log "GTDB-Tk Taxonomy Classification Pipeline"
    log "=========================================="
    log "Genome directory: ${genome_dir}"
    log "Output directory: ${output_dir}"
    log "Threads: ${THREADS}"
    log "Extension: ${EXTENSION}"
    log "Prefix: ${PREFIX}"
    log "Skip ANI screening: ${SKIP_ANI}"
    log "=========================================="

    # Run pipeline
    run_gtdbtk_classify ${genome_dir} ${output_dir}/classify
    run_gtdbtk_infer ${output_dir}/classify
    generate_summary ${output_dir}

    log "=========================================="
    log "Pipeline completed successfully!"
    log "=========================================="
    log "Results: ${output_dir}/classify/"
}

main "$@"