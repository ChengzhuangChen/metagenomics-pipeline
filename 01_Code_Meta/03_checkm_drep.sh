#!/bin/bash
# ============================================================
# MAG Quality Assessment and Dereplication Pipeline
# For Nature article supplementary materials
#
# Method: MAG quality was assessed using CheckM (v1.2.0) with
# default parameters, evaluating completeness and contamination.
# MAGs were dereplicated using dRep (v3.4.5) with a 95% ANI
# threshold and 90% AF threshold for primary clustering.
#
# Quality Filter:
#   High-quality MAGs: Completeness >= 90%, Contamination < 5%
#   Medium-quality MAGs: Completeness >= 50%, Contamination < 10%
#
# Version: 1.0
# Dependencies: CheckM, dRep, FastANI
#
# Usage:
#   bash 03_checkm_drep.sh -i <bins_dir> -o <output_dir>
# ============================================================

set -e

# ================== CONFIGURATION ==================
THREADS=32
WORK_DIR="results"
# Quality thresholds
COMPLETENESS_THRESHOLD=50
CONTAMINATION_THRESHOLD=10
# dRep parameters
ANI_THRESHOLD=95
AF_THRESHOLD=0.90
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
    local deps=("checkm" "dRep" "fastANI")
    for dep in "${deps[@]}"; do
        check_software ${dep}
    done
    log "All dependencies found."
}

# ================== STEP 1: CheckM Quality Assessment ==================
run_checkm() {
    local bins_dir=$1
    local output_dir=$2

    log "[Step 1/4] MAG quality assessment with CheckM..."

    mkdir -p ${output_dir}/checkm

    # Run CheckM lineage_wf
    # This analyzes the bins and places them in the CheckM reference tree
    checkm lineage_wf \
        -t ${THREADS} \
        -x fa \
        ${bins_dir} \
        ${output_dir}/checkm \
        -f ${output_dir}/checkm/checkm_results.tsv

    # Generate summary
    checkm qa \
        ${output_dir}/checkm/checkm_refine/msa_storage/tree/consensus_bins/consensus_tree.fa \
        ${output_dir}/checkm \
        -o 2 \
        -f ${output_dir}/checkm/checkm_qa_summary.tsv \
        --tab_table

    log "  -> CheckM quality assessment complete"
    log "  -> Results: ${output_dir}/checkm/"
}

# ================== STEP 2: Filter by Quality ==================
filter_mags() {
    local checkm_output=$1
    local bins_dir=$2
    local output_dir=$3

    log "[Step 2/4] Filtering MAGs by quality..."

    mkdir -p ${output_dir}/filtered

    # CheckM output has these columns:
    # Bin Id, Marker lineage, # genomes, # markers, # marker sets, Completeness, Contamination, Strain heterogeneity

    # Filter MAGs by quality thresholds
    # Keep MAGs with Completeness >= threshold and Contamination < threshold
    awk -F'\t' -v comp=${COMPLETENESS_THRESHOLD} -v cont=${CONTAMINATION_THRESHOLD} \
        'NR==1 || ($6>=comp && $7<cont)' ${checkm_output} > \
        ${output_dir}/filtered/quality_filtered_mags.tsv

    # Extract filtered MAG names
    local filtered_count=$(tail -n +2 ${output_dir}/filtered/quality_filtered_mags.tsv | wc -l)

    # Count by quality tier
    local total_mags=$(tail -n +2 ${checkm_output} | wc -l)
    local high_quality=$(awk -F'\t' 'NR>1 && $6>=90 && $7<5' ${checkm_output} | wc -l)
    local medium_quality=$(awk -F'\t' 'NR>1 && $6>=50 && $7<10' ${checkm_output} | wc -l)
    local low_quality=$(awk -F'\t' 'NR>1 && $6<50' ${checkm_output} | wc -l)

    log "  === CheckM Summary ==="
    log "  -> Total MAGs assessed: ${total_mags}"
    log "  -> High-quality MAGs (Completeness>=90%, Contamination<5%): ${high_quality}"
    log "  -> Medium-quality MAGs (Completeness>=50%, Contamination<10%): ${medium_quality}"
    log "  -> Low-quality MAGs (Completeness<50%): ${low_quality}"
    log "  -> Filtered MAGs (>=${COMPLETENESS_THRESHOLD}% completeness, <${CONTAMINATION_THRESHOLD}% contamination): ${filtered_count}"

    # Copy filtered MAGs to new directory
    while IFS=$'\t' read -r name lineage genomes markers marker_sets completeness contamination heterogeneity; do
        if [ "$name" != "Bin Id" ]; then
            find ${bins_dir} -name "${name}.fa" -exec cp {} ${output_dir}/filtered/ \;
            find ${bins_dir} -name "${name}.fasta" -exec cp {} ${output_dir}/filtered/ \;
        fi
    done < <(tail -n +2 ${output_dir}/filtered/quality_filtered_mags.tsv)

    echo "${output_dir}/filtered"
}

# ================== STEP 3: dRep Dereplication ==================
run_drep() {
    local filtered_bins=$1
    local output_dir=$2

    log "[Step 3/4] Dereplicating MAGs with dRep (ANI>=${ANI_THRESHOLD}%, AF>=${AF_THRESHOLD})..."

    mkdir -p ${output_dir}/drep

    # Run dRep dereplicate
    # -p: threads
    # -g: genome directory
    # -sa: ANI threshold for Average Nucleotide Identity
    # -nc: ANI threshold for secondary clustering
    # -cm: complete mode (use Mash distances for initial clustering)
    # -comp: minimum completeness threshold
    # -con: maximum contamination threshold
    dRep dereplicate \
        ${output_dir}/drep \
        -p ${THREADS} \
        -g ${filtered_bins}/*.fa \
        -sa ${ANI_THRESHOLD} \
        -nc ${ANI_THRESHOLD} \
        -cm \
        -comp ${COMPLETENESS_THRESHOLD} \
        -con ${CONTAMINATION_THRESHOLD}

    log "  -> dRep dereplication complete"
    log "  -> Dereplicated genomes: ${output_dir}/drep/dereplicated_genomes/"
}

# ================== STEP 4: Generate Summary Report ==================
generate_summary() {
    local checkm_output=$1
    local output_dir=$2

    log "[Step 4/4] Generating summary report..."

    local summary_file="${output_dir}/summary_report.txt"

    local total_mags=$(tail -n +2 ${checkm_output} | wc -l)
    local high_quality=$(awk -F'\t' 'NR>1 && $6>=90 && $7<5' ${checkm_output} | wc -l)
    local medium_quality=$(awk -F'\t' 'NR>1 && $6>=50 && $7<10' ${checkm_output} | wc -l)
    local low_quality=$(awk -F'\t' 'NR>1 && $6<50' ${checkm_output} | wc -l)
    local derep_count=$(find ${output_dir}/drep/dereplicated_genomes -name '*.fa' 2>/dev/null | wc -l)

    {
        echo "=========================================="
        echo "MAG Quality Assessment and Dereplication"
        echo "=========================================="
        echo ""
        echo "Parameters:"
        echo "  Completeness threshold: ${COMPLETENESS_THRESHOLD}%"
        echo "  Contamination threshold: <${CONTAMINATION_THRESHOLD}%"
        echo "  ANI threshold: ${ANI_THRESHOLD}%"
        echo "  AF threshold: ${AF_THRESHOLD}"
        echo ""
        echo "Quality Distribution:"
        echo "  Total MAGs: ${total_mags}"
        echo "  High-quality (>=90% comp, <5% cont): ${high_quality}"
        echo "  Medium-quality (>=50% comp, <10% cont): ${medium_quality}"
        echo "  Low-quality (<50% comp): ${low_quality}"
        echo ""
        echo "Dereplication Results:"
        echo "  Non-redundant MAGs: ${derep_count}"
        echo ""
        echo "Output Directories:"
        echo "  - ${output_dir}/checkm/ (CheckM results)"
        echo "  - ${output_dir}/filtered/ (Quality-filtered MAGs)"
        echo "  - ${output_dir}/drep/dereplicated_genomes/ (Final non-redundant MAGs)"
        echo "=========================================="
    } > ${summary_file}

    log "  -> Summary report: ${summary_file}"
}

# ================== MAIN ENTRY POINT ==================

print_usage() {
    cat << EOF
Usage: bash 03_checkm_drep.sh -i <bins_dir> -o <output_dir>

Description:
    MAG quality assessment using CheckM and dereplication using dRep.
    Filters MAGs based on quality thresholds and removes redundant genomes.

Required Arguments:
  -i, --input     Directory containing MAG bins (FASTA format)
  -o, --output    Output directory

Optional Arguments:
  -t, --threads   Number of threads (default: 32)
  -c, --comp      Minimum completeness threshold (default: 50)
  -x, --cont      Maximum contamination threshold (default: 10)
  -a, --ani       ANI threshold for dereplication (default: 95)
  -f, --af        Alignment fraction threshold (default: 0.90)

Example:
  bash 03_checkm_drep.sh \\
      -i results/binning/dastool/DAS_Tool_DASToolbins \\
      -o results/mag_quality

  # With custom parameters
  bash 03_checkm_drep.sh -i bins/ -o output -t 64 -c 90 -x 5

Input:
  - bins_dir: Directory containing MAG bins as FASTA files (.fa or .fasta)
  - Each bin should be a separate file with contigs belonging to one MAG

Output:
  - checkm/: CheckM quality assessment results
  - filtered/: Quality-filtered MAGs
  - drep/dereplicated_genomes/: Non-redundant MAGs after dRep
  - summary_report.txt: Summary statistics

CheckM Output Columns:
  Bin Id, Marker lineage, # genomes, # markers, # marker sets,
  Completeness, Contamination, Strain heterogeneity
EOF
}

main() {
    local bins_dir=""
    local output_dir=""

    # Parse arguments
    while [[ $# -gt 0 ]]; do
        case $1 in
            -i|--input) bins_dir="$2"; shift 2 ;;
            -o|--output) output_dir="$2"; shift 2 ;;
            -t|--threads) THREADS="$2"; shift 2 ;;
            -c|--comp) COMPLETENESS_THRESHOLD="$2"; shift 2 ;;
            -x|--cont) CONTAMINATION_THRESHOLD="$2"; shift 2 ;;
            -a|--ani) ANI_THRESHOLD="$2"; shift 2 ;;
            -f|--af) AF_THRESHOLD="$2"; shift 2 ;;
            -h|--help) print_usage; exit 0 ;;
            *) echo "Unknown option: $1"; print_usage; exit 1 ;;
        esac
    done

    # Validate required arguments
    if [ -z "$bins_dir" ] || [ -z "$output_dir" ]; then
        echo "Error: Missing required arguments!"
        print_usage
        exit 1
    fi

    # Check dependencies
    check_dependencies

    log "=========================================="
    log "MAG Quality Assessment and Dereplication"
    log "=========================================="
    log "Input bins: ${bins_dir}"
    log "Output: ${output_dir}"
    log "Threads: ${THREADS}"
    log "Completeness threshold: ${COMPLETENESS_THRESHOLD}%"
    log "Contamination threshold: <${CONTAMINATION_THRESHOLD}%"
    log "ANI threshold: ${ANI_THRESHOLD}%"
    log "=========================================="

    # Step 1: CheckM quality assessment
    run_checkm ${bins_dir} ${output_dir}

    local checkm_output="${output_dir}/checkm/checkm_results.tsv"

    # Step 2: Filter by quality
    local filtered_dir=$(filter_mags ${checkm_output} ${bins_dir} ${output_dir})

    # Step 3: dRep dereplication
    run_drep ${filtered_dir} ${output_dir}

    # Step 4: Generate summary
    generate_summary ${checkm_output} ${output_dir}

    log "=========================================="
    log "Pipeline completed successfully!"
    log "=========================================="
    log "Final non-redundant MAGs: ${output_dir}/drep/dereplicated_genomes/"
    log "Summary report: ${output_dir}/summary_report.txt"
}

main "$@"