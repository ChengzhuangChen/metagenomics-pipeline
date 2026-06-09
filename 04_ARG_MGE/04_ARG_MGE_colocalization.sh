#!/bin/bash
# ============================================================
# ARG-MGE Co-localization Analysis Pipeline (10 kb window)
# For Nature article supplementary materials
#
# Method: To investigate potential horizontal gene transfer (HGT)
# of antibiotic resistance genes, ARG-MGE co-localization analysis
# was performed. ARGs and MGEs within 10 kb on the same contig
# were regarded as co-occurring. Contigs meeting this criterion
# were identified, and their functional annotations were extracted.
#
# Version: 1.0
# Dependencies: Python 3.8+, Biopython
#
# Usage:
#   bash 04_ARG_MGE_colocalization.sh -r <arg_csv> -m <mge_csv> -g <gff_dir> -f <faa_dir> -o <output>
# ============================================================

set -e

# ================== CONFIGURATION ==================
THREADS=32
DISTANCE_THRESHOLD=10000              # 10 kb = 10000 bp
# ====================================================

log() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $1"
}

error_exit() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] ERROR: $1" >&2
    exit 1
}

check_dependencies() {
    log "Checking dependencies..."
    if ! command -v python3 &> /dev/null; then
        error_exit "Python3 not found"
    fi
    log "All dependencies found."
}

# ================== STEP 1: Identify ARG-MGE Co-localized ORFs ==================
identify_cocol_orfs() {
    local arg_csv=$1
    local mge_csv=$2
    local gff_dir=$3
    local output_dir=$4

    log "[Step 1/4] Identifying ARGs and MGEs within ${DISTANCE_THRESHOLD} bp (10 kb)..."

    local cocol_dir="${output_dir}/01_cocol_pairs"
    mkdir -p ${cocol_dir}

    python3 identify_cocol_contigs.py \
        --arg-csv ${arg_csv} \
        --mge-csv ${mge_csv} \
        --gff-dir ${gff_dir} \
        --output ${cocol_dir}/cocol_pairs_within_10kb.tsv \
        --stats ${cocol_dir}/cocol_stats.txt \
        --distance ${DISTANCE_THRESHOLD}

    local pair_count=$(tail -n +2 ${cocol_dir}/cocol_pairs_within_10kb.tsv 2>/dev/null | wc -l || echo 0)
    local cocol_contigs=$(cat ${cocol_dir}/cocol_pairs_within_10kb_contigs.txt 2>/dev/null | wc -l || echo 0)

    log "  -> Found ${pair_count} ARG-MGE pairs within ${DISTANCE_THRESHOLD} bp"
    log "  -> Identified ${cocol_contigs} co-localized contigs"
    log "  -> Results: ${cocol_dir}/"

    echo "${cocol_dir}/cocol_pairs_within_10kb_contigs.txt"
}

# ================== STEP 2: Extract Contig Sequences ==================
extract_sequences() {
    local cocol_contigs=$1
    local all_contigs=$2
    local output_dir=$3

    log "[Step 2/4] Extracting co-localized contig sequences..."

    local seq_dir="${output_dir}/02_cocol_sequences"
    mkdir -p ${seq_dir}

    seqtk subseq ${all_contigs} ${cocol_contigs} > ${seq_dir}/cocol_contigs.fa

    local seq_count=$(grep -c "^>" ${seq_dir}/cocol_contigs.fa)
    log "  -> Extracted ${seq_count} contig sequences"
    log "  -> Sequences: ${seq_dir}/cocol_contigs.fa"

    echo "${seq_dir}/cocol_contigs.fa"
}

# ================== STEP 3: Extract Functional Annotations ==================
extract_annotations() {
    local cocol_contigs=$1
    local faa_dir=$2
    local gff_dir=$3
    local output_dir=$4

    log "[Step 3/4] Extracting functional annotations for co-localized contigs..."

    local annot_dir="${output_dir}/03_annotations"
    mkdir -p ${annot_dir}

    python3 extract_cocol_annotations.py \
        --contigs ${cocol_contigs} \
        --faa-dir ${faa_dir} \
        --gff-dir ${gff_dir} \
        --output-dir ${annot_dir}

    log "  -> Annotations: ${annot_dir}/"
}

# ================== STEP 4: Generate Summary Report ==================
generate_summary() {
    local output_dir=$1
    local pair_count=$2
    local cocol_count=$3
    local arg_contig_count=$4
    local mge_contig_count=$5

    log "[Step 4/4] Generating summary report..."

    local summary_file="${output_dir}/summary_report.txt"

    {
        echo "=========================================="
        echo "ARG-MGE Co-localization Analysis Summary"
        echo "=========================================="
        echo ""
        echo "Criterion: ARGs and MGEs within ${DISTANCE_THRESHOLD} bp (10 kb)"
        echo "           on the same contig were regarded as co-occurring."
        echo ""
        echo "Analysis Parameters:"
        echo "  Distance threshold: ${DISTANCE_THRESHOLD} bp (10 kb)"
        echo ""
        echo "Basic Statistics:"
        echo "  Contigs with ARGs: ${arg_contig_count}"
        echo "  Contigs with MGEs: ${mge_contig_count}"
        echo ""
        echo "Co-localization Results:"
        echo "  ARG-MGE pairs within 10 kb: ${pair_count}"
        echo "  Co-localized contigs: ${cocol_count}"
        local rate=$(echo "scale=2; ${cocol_count}/${arg_contig_count}*100" | bc 2>/dev/null || echo "N/A")
        echo "  Co-localization rate: ${rate}%"
        echo ""
        echo "Output Directories:"
        echo "  - ${output_dir}/01_cocol_pairs/"
        echo "  - ${output_dir}/02_cocol_sequences/"
        echo "  - ${output_dir}/03_annotations/"
        echo ""
        echo "Key Output Files:"
        echo "  - 01_cocol_pairs/cocol_pairs_within_10kb.tsv"
        echo "  - 01_cocol_pairs/cocol_stats.txt"
        echo "  - 02_cocol_sequences/cocol_contigs.fa"
        echo "  - 03_annotations/cocol_proteins.faa"
        echo "  - 03_annotations/cocol_annotation_summary.tsv"
        echo "=========================================="
    } > ${summary_file}

    log "  -> Summary: ${summary_file}"
}

# ================== MAIN ENTRY POINT ==================

print_usage() {
    cat << EOF
Usage: bash 04_ARG_MGE_colocalization.sh -r <arg_csv> -m <mge_csv> -g <gff_dir> -f <faa_dir> -o <output>

Description:
    Analyze co-localization of antibiotic resistance genes (ARGs)
    and mobile genetic elements (MGEs) on the same contigs.

    Criterion: ARGs and MGEs within 10 kb on the same contig
    were regarded as co-occurring.

Required Arguments:
  -r, --arg-csv       Merged ARG DIAMOND results (all_samples_merged_arg.csv)
  -m, --mge-csv       Merged MGE DIAMOND results (all_samples_merged_mge.csv)
  -g, --gff-dir       Directory with Prodigal GFF files
  -f, --faa-dir       Directory with Prodigal FAA files
  -o, --output        Output directory

Optional Arguments:
  -t, --threads       Number of threads (default: 32)
  -d, --distance      Distance threshold in bp (default: 10000 = 10 kb)

Example:
  bash 04_ARG_MGE_colocalization.sh \\
      -r results/annotation/00_ARG_ALL/all_samples_merged_arg.csv \\
      -m results/annotation/00_MGE_ALL/all_samples_merged_mge.csv \\
      -g results/prodigal \\
      -f results/prodigal \\
      -o results/arg_mge_analysis

Input:
  - all_samples_merged_arg.csv: ARG DIAMOND results (from 03_Contigs)
  - all_samples_merged_mge.csv: MGE DIAMOND results (from 03_Contigs)
  - gff/*.gff: Prodigal gene annotations (for position info)
  - faa/*.faa: Prodigal-predicted protein sequences

Output:
  - 01_cocol_pairs/cocol_pairs_within_10kb.tsv: ARG-MGE pairs within 10 kb
  - 01_cocol_pairs/cocol_stats.txt: Statistics
  - 02_cocol_sequences/cocol_contigs.fa: Extracted sequences
  - 03_annotations/cocol_proteins.faa: All proteins
  - 03_annotations/cocol_annotation_summary.tsv: Annotation summary
EOF
}

main() {
    local arg_csv=""
    local mge_csv=""
    local faa_dir=""
    local gff_dir=""
    local output_dir=""

    while [[ $# -gt 0 ]]; do
        case $1 in
            -r|--arg-csv) arg_csv="$2"; shift 2 ;;
            -m|--mge-csv) mge_csv="$2"; shift 2 ;;
            -g|--gff-dir) gff_dir="$2"; shift 2 ;;
            -f|--faa-dir) faa_dir="$2"; shift 2 ;;
            -o|--output) output_dir="$2"; shift 2 ;;
            -t|--threads) THREADS="$2"; shift 2 ;;
            -d|--distance) DISTANCE_THRESHOLD="$2"; shift 2 ;;
            -h|--help) print_usage; exit 0 ;;
            *) echo "Unknown option: $1"; print_usage; exit 1 ;;
        esac
    done

    # Validate arguments
    if [ -z "$arg_csv" ] || [ -z "$mge_csv" ] || [ -z "$gff_dir" ] || \
       [ -z "$faa_dir" ] || [ -z "$output_dir" ]; then
        echo "Error: Missing required arguments!"
        print_usage
        exit 1
    fi

    check_dependencies

    log "=========================================="
    log "ARG-MGE Co-localization Analysis"
    log "Distance threshold: ${DISTANCE_THRESHOLD} bp (10 kb)"
    log "=========================================="
    log "ARG CSV: ${arg_csv}"
    log "MGE CSV: ${mge_csv}"
    log "GFF directory: ${gff_dir}"
    log "FAA directory: ${faa_dir}"
    log "Output directory: ${output_dir}"
    log "=========================================="

    # Step 1: Identify co-localized ORFs (using GFF positions)
    local cocol_contigs=$(identify_cocol_orfs ${arg_csv} ${mge_csv} ${gff_dir} ${output_dir})

    # Get statistics for summary
    local pair_count=$(tail -n +2 ${output_dir}/01_cocol_pairs/cocol_pairs_within_10kb.tsv 2>/dev/null | wc -l || echo 0)
    local cocol_count=$(wc -l < ${cocol_contigs})

    # Extract all contigs with ARGs and MGEs for rate calculation
    local arg_contig_count=$(tail -n +2 ${output_dir}/01_cocol_pairs/cocol_stats.txt 2>/dev/null | grep "Contigs with ARGs:" | awk '{print $4}')
    local mge_contig_count=$(tail -n +2 ${output_dir}/01_cocol_pairs/cocol_stats.txt 2>/dev/null | grep "Contigs with MGEs:" | awk '{print $4}')

    # Step 2: Extract sequences (if we have the combined contigs file)
    # Note: This requires all_arg_contigs.fa from 03_Contigs step
    if [ -f "results/annotation/03_ARG_total/all_arg_contigs.fa" ]; then
        extract_sequences ${cocol_contigs} "results/annotation/03_ARG_total/all_arg_contigs.fa" ${output_dir}
    else
        log "Note: all_arg_contigs.fa not found, skipping sequence extraction"
        mkdir -p ${output_dir}/02_cocol_sequences
    fi

    # Step 3: Extract annotations
    extract_annotations ${cocol_contigs} ${faa_dir} ${gff_dir} ${output_dir}

    # Step 4: Generate summary
    generate_summary ${output_dir} ${pair_count} ${cocol_count} ${arg_contig_count:-0} ${mge_contig_count:-0}

    log "=========================================="
    log "Pipeline completed successfully!"
    log "=========================================="
}

main "$@"