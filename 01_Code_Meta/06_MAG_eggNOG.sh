#!/bin/bash
# ============================================================
# MAG Functional Annotation Pipeline
# For Nature article supplementary materials
#
# Method: Protein-coding genes were predicted using Prodigal (v2.6.3)
# for each MAG. Predicted proteins were clustered using MMseqs2
# (v13.10711) with 95% identity and 80% coverage thresholds to
# remove redundancy. Functional annotation was performed using
# eggNOG-mapper (v2.1.12) with DIAMOND against the eggNOG
# database v2.1.12.
#
# Version: 1.0
# Dependencies: Prodigal, MMseqs2, eggNOG-mapper
#
# Usage:
#   bash 06_MAG_eggNOG.sh -i <mags_dir> -o <output_dir>
# ============================================================

set -e

# ================== CONFIGURATION ==================
THREADS=32
# Prodigal parameters
TRANSLATION_TABLE=11              # Bacterial genetic code
# MMseqs2 clustering parameters
MMSEQS_MIN_ID=0.95               # 95% sequence identity
MMSEQS_COV=0.80                  # 80% coverage
# eggNOG parameters
EGGNOG_EVALUE=1e-5
EGGNOG_SENSMODE="sensitive"
EGGNOG_DB_DIR=""                 # Set if custom database path needed
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
    local deps=("prodigal" "mmseqs" "emapper.py")
    for dep in "${deps[@]}"; do
        check_software ${dep}
    done
    log "All dependencies found."
}

# ================== STEP 1: Gene Prediction with Prodigal ==================
run_prodigal() {
    local mags_dir=$1
    local output_dir=$2

    log "[Step 1/4] Gene prediction with Prodigal v2.6.3..."

    local faa_dir="${output_dir}/01_prodigal_faa"
    mkdir -p ${faa_dir}

    local total_mags=0
    local total_genes=0

    for mag in ${mags_dir}/*.fa ${mags_dir}/*.fasta; do
        if [ -f "$mag" ]; then
            local mag_name=$(basename ${mag} .fa)
            mag_name=$(basename ${mag_name} .fasta)

            log "  -> Processing: ${mag_name}"

            # Predict genes (translated proteins)
            prodigal \
                -i ${mag} \
                -o ${faa_dir}/${mag_name}.gff \
                -a ${faa_dir}/${mag_name}.faa \
                -p meta \
                -g ${TRANSLATION_TABLE} \
                -q

            local gene_count=$(grep -c "^>" ${faa_dir}/${mag_name}.faa)
            total_genes=$((total_genes + gene_count))
            total_mags=$((total_mags + 1))
        fi
    done

    log "  === Prodigal Summary ==="
    log "  -> MAGs processed: ${total_mags}"
    log "  -> Total genes predicted: ${total_genes}"
    log "  -> FAA files: ${faa_dir}/"

    echo "${faa_dir}"
    echo "${total_genes}"
}

# ================== STEP 2: Merge FAA Files ==================
merge_faa() {
    local faa_dir=$1
    local output_dir=$2

    log "[Step 2/4] Merging FAA files..."

    local merged_faa="${output_dir}/02_merge_faa/all_genes.faa"
    mkdir -p ${output_dir}/02_merge_faa

    cat ${faa_dir}/*.faa > ${merged_faa}

    local total_sequences=$(grep -c "^>" ${merged_faa})
    log "  -> Merged FAA: ${merged_faa}"
    log "  -> Total sequences: ${total_sequences}"

    echo "${merged_faa}"
}

# ================== STEP 3: MMseqs2 Clustering ==================
run_mmseqs2() {
    local merged_faa=$1
    local output_dir=$2

    log "[Step 3/4] Clustering with MMseqs2..."

    local cluster_dir="${output_dir}/03_mmseqs2_result"
    mkdir -p ${cluster_dir}

    local tmp_dir="${cluster_dir}/mmseqs_tmp"
    mkdir -p ${tmp_dir}

    local out_prefix="${cluster_dir}/Lakesed_clustered"

    # Run MMseqs2 easy-cluster
    # --min-seq-id: Minimum sequence identity (0.95 = 95%)
    # -c: Coverage threshold (0.80 = 80%)
    # --cov-mode 1: Coverage of query and target
    # --cluster-mode 2: Greedy set cover algorithm
    mmseqs easy-cluster \
        ${merged_faa} \
        ${out_prefix} \
        ${tmp_dir} \
        --min-seq-id ${MMSEQS_MIN_ID} \
        -c ${MMSEQS_COV} \
        --cov-mode 1 \
        --threads ${THREADS} \
        --cluster-mode 2

    # Count clusters (representative sequences)
    local cluster_count=$(grep -c "^>" ${out_prefix}_rep_seq.faa)

    log "  === MMseqs2 Clustering Summary ==="
    log "  -> Clustered genes: $(grep -c "^>" ${merged_faa})"
    log "  -> Non-redundant clusters (unigenes): ${cluster_count}"
    log "  -> Representative sequences: ${out_prefix}_rep_seq.faa"
    log "  -> All cluster members: ${out_prefix}_all_seqs.fasta"

    echo "${out_prefix}_rep_seq.faa"
    echo "${cluster_count}"
}

# ================== STEP 4: eggNOG Annotation ==================
run_eggnog() {
    local rep_faa=$1
    local output_dir=$2

    log "[Step 4/4] Functional annotation with eggNOG-mapper..."

    local eggnog_dir="${output_dir}/04_eggnog"
    mkdir -p ${eggnog_dir}

    # Build command
    local cmd="emapper.py \
        -i ${rep_faa} \
        -o ${eggnog_dir}/annotation \
        -m diamond \
        --cpu ${THREADS} \
        --evalue ${EGGNOG_EVALUE} \
        --sensmode ${EGGNOG_SENSMODE} \
        --report_orthologs"

    # Add database path if specified
    if [ -n "${EGGNOG_DB_DIR}" ]; then
        cmd="${cmd} --data_dir ${EGGNOG_DB_DIR}"
    fi

    # Run eggNOG-mapper
    eval ${cmd}

    log "  === eggNOG Annotation Summary ==="

    # Count annotated genes
    if [ -f "${eggnog_dir}/annotation.emapper.annotations" ]; then
        local annotated=$(tail -n +2 ${eggnog_dir}/annotation.emapper.annotations | wc -l)
        log "  -> Annotated unigenes: ${annotated}"
        log "  -> Full annotations: ${eggnog_dir}/annotation.emapper.annotations"
    fi

    if [ -f "${eggnog_dir}/annotation.emapper.seed_orthologs" ]; then
        log "  -> Seed orthologs: ${eggnog_dir}/annotation.emapper.seed_orthologs"
    fi

    echo "${eggnog_dir}"
}

# ================== STEP 5: Generate Summary ==================
generate_summary() {
    local output_dir=$1
    local total_genes=$2
    local cluster_count=$3

    log "[Step 5/5] Generating summary report..."

    local summary_file="${output_dir}/summary_report.txt"

    {
        echo "=========================================="
        echo "MAG Functional Annotation Pipeline"
        echo "=========================================="
        echo ""
        echo "Parameters:"
        echo "  Threads: ${THREADS}"
        echo "  Translation table: ${TRANSLATION_TABLE} (Bacterial)"
        echo "  MMseqs2 min identity: ${MMSEQS_MIN_ID} (${MMSEQS_MIN_ID}%)"
        echo "  MMseqs2 coverage: ${MMSEQS_COV} (${MMSEQS_COV}%)"
        echo "  eggNOG evalue: ${EGGNOG_EVALUE}"
        echo "  eggNOG sensmode: ${EGGNOG_SENSMODE}"
        echo ""
        echo "Pipeline Results:"
        echo "  Total genes predicted: ${total_genes}"
        echo "  Non-redundant unigenes: ${cluster_count}"
        echo "  Clustering reduction: $(echo "scale=1; (1-${cluster_count}/${total_genes})*100" | bc)%"
        echo ""
        echo "Output Directories:"
        echo "  - ${output_dir}/01_prodigal_faa/"
        echo "  - ${output_dir}/02_merge_faa/"
        echo "  - ${output_dir}/03_mmseqs2_result/"
        echo "  - ${output_dir}/04_eggnog/"
        echo ""
        echo "Key Output Files:"
        echo "  - ${output_dir}/03_mmseqs2_result/Lakesed_clustered_rep_seq.faa"
        echo "  - ${output_dir}/04_eggnog/annotation.emapper.annotations"
        echo "=========================================="
    } > ${summary_file}

    log "  -> Summary report: ${summary_file}"
}

# ================== MAIN ENTRY POINT ==================

print_usage() {
    cat << EOF
Usage: bash 06_MAG_eggNOG.sh -i <mags_dir> -o <output_dir>

Description:
    MAG functional annotation pipeline including:
    1. Gene prediction with Prodigal
    2. Merging protein sequences
    3. Clustering with MMseqs2 (95% identity, 80% coverage)
    4. Functional annotation with eggNOG-mapper

Required Arguments:
  -i, --input     Directory containing MAG FASTA files
  -o, --output    Output directory

Optional Arguments:
  -t, --threads   Number of threads (default: 32)
  --min-id        MMseqs2 minimum identity (default: 0.95)
  --cov           MMseqs2 coverage threshold (default: 0.80)
  --evalue        eggNOG evalue threshold (default: 1e-5)
  --db-dir        eggNOG database directory (if custom location)

Example:
  bash 06_MAG_eggNOG.sh \\
      -i results/mag_quality/drep/dereplicated_genomes \\
      -o results/functional_annotation

  # With custom parameters
  bash 06_MAG_eggNOG.sh -i bins/ -o output -t 64 --min-id 0.98

Input:
  - mags_dir: Directory containing MAG FASTA files (.fa or .fasta)

Output:
  - 01_prodigal_faa/: Individual MAG protein predictions
  - 02_merge_faa/all_genes.faa: Combined protein sequences
  - 03_mmseqs2_result/
      - Lakesed_clustered_rep_seq.faa: Non-redundant proteins
      - Lakesed_clustered_all_seqs.fasta: All cluster members
  - 04_eggnog/
      - annotation.emapper.annotations: Full annotation results
      - annotation.emapper.seed_orthologs: Best hits

Installation:
  # Create conda environment
  conda create -n annotation -c conda-forge -c bioconda \
      prodigal=2.6.3 \
      mmseqs2=13.10711 \
      eggnog-mapper=2.1.12

  conda activate annotation

  # Download eggNOG database (first time only)
  download_eggnog_data.py
EOF
}

main() {
    local mags_dir=""
    local output_dir=""

    # Parse arguments
    while [[ $# -gt 0 ]]; do
        case $1 in
            -i|--input) mags_dir="$2"; shift 2 ;;
            -o|--output) output_dir="$2"; shift 2 ;;
            -t|--threads) THREADS="$2"; shift 2 ;;
            --min-id) MMSEQS_MIN_ID="$2"; shift 2 ;;
            --cov) MMSEQS_COV="$2"; shift 2 ;;
            --evalue) EGGNOG_EVALUE="$2"; shift 2 ;;
            --db-dir) EGGNOG_DB_DIR="$2"; shift 2 ;;
            -h|--help) print_usage; exit 0 ;;
            *) echo "Unknown option: $1"; print_usage; exit 1 ;;
        esac
    done

    # Validate required arguments
    if [ -z "$mags_dir" ] || [ -z "$output_dir" ]; then
        echo "Error: Missing required arguments!"
        print_usage
        exit 1
    fi

    # Check dependencies
    check_dependencies

    log "=========================================="
    log "MAG Functional Annotation Pipeline"
    log "Prodigal + MMseqs2 + eggNOG-mapper"
    log "=========================================="
    log "MAG directory: ${mags_dir}"
    log "Output directory: ${output_dir}"
    log "Threads: ${THREADS}"
    log "MMseqs2 min identity: ${MMSEQS_MIN_ID}"
    log "MMseqs2 coverage: ${MMSEQS_COV}"
    log "=========================================="

    # Step 1: Prodigal gene prediction
    local prodigal_result=$(run_prodigal ${mags_dir} ${output_dir})
    local faa_dir=$(echo "$prodigal_result" | head -n1)
    local total_genes=$(echo "$prodigal_result" | tail -n1)

    # Step 2: Merge FAA files
    local merged_faa=$(merge_faa ${faa_dir} ${output_dir})

    # Step 3: MMseqs2 clustering
    local cluster_result=$(run_mmseqs2 ${merged_faa} ${output_dir})
    local rep_faa=$(echo "$cluster_result" | head -n1)
    local cluster_count=$(echo "$cluster_result" | tail -n1)

    # Step 4: eggNOG annotation
    run_eggnog ${rep_faa} ${output_dir}

    # Step 5: Generate summary
    generate_summary ${output_dir} ${total_genes} ${cluster_count}

    log "=========================================="
    log "Pipeline completed successfully!"
    log "=========================================="
    log "Final annotations: ${output_dir}/04_eggnog/"
}

main "$@"