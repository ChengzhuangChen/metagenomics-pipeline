#!/bin/bash
# ============================================================
# Contig-level Annotation Pipeline (ARG, MGE, VF)
# For Nature article supplementary materials
#
# Method: Protein sequences predicted by Prodigal were aligned
# against three databases using DIAMOND: SARG (antibiotic resistance
# genes), MobileOG (mobile genetic elements), and VFDB (virulence
# factors). Contigs harboring annotated ORFs were extracted for
# downstream analysis.
#
# Parameters:
#   Identity: >= 80%
#   Query coverage: >= 70%
#   E-value: 1e-7
#
# Version: 1.0
# Dependencies: DIAMOND, Prodigal, seqtk
#
# Usage:
#   bash 03_Contigs_annotation.sh -i <protein_dir> -c <contig_dir> -o <output_dir> -a <sarg_db> -m <mge_db> -v <vf_db>
# ============================================================

set -e

# ================== CONFIGURATION ==================
THREADS=32
MIN_IDENTITY=80                  # Minimum identity (%)
MIN_QUERY_COVER=70               # Minimum query coverage (%)
EVALUE=1e-7                      # E-value threshold
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
    local deps=("diamond" "prodigal" "seqtk")
    for dep in "${deps[@]}"; do
        if ! command -v ${dep} &> /dev/null; then
            error_exit "${dep} not found. Please install it first."
        fi
    done
    log "All dependencies found."
}

# ================== STEP 1: ARG Annotation ==================
run_arg_annotation() {
    local protein_dir=$1
    local contig_dir=$2
    local output_dir=$3
    local sarg_db=$4

    log "[Step 1/3] ARG annotation against SARG database..."

    local arg_out="${output_dir}/01_ARG_out"
    mkdir -p ${arg_out} ${output_dir}/02_ARG_contigs ${output_dir}/03_ARG_total

    # Define output files
    local arg_result="${arg_out}/${SAMPLE_ID}_diamond_blastp_arg.txt"

    # Run DIAMOND blastp
    log "  -> Running DIAMOND against SARG..."
    diamond blastp \
        --query ${protein_dir}/${SAMPLE_ID}.faa \
        --db ${sarg_db} \
        --out ${arg_result} \
        --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore \
        --evalue ${EVALUE} \
        --id ${MIN_IDENTITY} \
        --query-cover ${MIN_QUERY_COVER} \
        --max-target-seqs 1 \
        --sensitive \
        --max-hsps 1 \
        --threads ${THREADS}

    # Add header to result file
    add_header "${arg_result}"

    # Extract ARG-containing contigs
    extract_contigs "${arg_result}" "${contig_dir}/${SAMPLE_ID}_contig.fa" \
                    "${output_dir}/02_ARG_contigs/${SAMPLE_ID}_arg_contigs.fa" "arg"

    log "  -> ARG annotation complete"
}

# ================== STEP 2: MGE Annotation ==================
run_mge_annotation() {
    local protein_dir=$1
    local contig_dir=$2
    local output_dir=$3
    local mge_db=$4

    log "[Step 2/3] MGE annotation against MobileOG database..."

    local mge_out="${output_dir}/04_MGE_out"
    mkdir -p ${mge_out} ${output_dir}/05_MGE_contigs ${output_dir}/06_MGE_total

    local mge_result="${mge_out}/${SAMPLE_ID}_diamond_blastp_mge.txt"

    # Run DIAMOND blastp
    log "  -> Running DIAMOND against MobileOG..."
    diamond blastp \
        --query ${protein_dir}/${SAMPLE_ID}.faa \
        --db ${mge_db} \
        --out ${mge_result} \
        --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore \
        --evalue ${EVALUE} \
        --id ${MIN_IDENTITY} \
        --query-cover ${MIN_QUERY_COVER} \
        --max-target-seqs 1 \
        --sensitive \
        --max-hsps 1 \
        --threads ${THREADS}

    # Add header
    add_header "${mge_result}"

    # Extract MGE-containing contigs
    extract_contigs "${mge_result}" "${contig_dir}/${SAMPLE_ID}_contig.fa" \
                    "${output_dir}/05_MGE_contigs/${SAMPLE_ID}_mge_contigs.fa" "mge"

    log "  -> MGE annotation complete"
}

# ================== STEP 3: VF Annotation ==================
run_vf_annotation() {
    local protein_dir=$1
    local contig_dir=$2
    local output_dir=$3
    local vf_db=$4

    log "[Step 3/3] VF annotation against VFDB..."

    local vf_out="${output_dir}/07_VF_out"
    mkdir -p ${vf_out} ${output_dir}/08_VF_contigs ${output_dir}/09_VF_total

    local vf_result="${vf_out}/${SAMPLE_ID}_diamond_blastp_vf.txt"

    # Run DIAMOND blastp
    log "  -> Running DIAMOND against VFDB..."
    diamond blastp \
        --query ${protein_dir}/${SAMPLE_ID}.faa \
        --db ${vf_db} \
        --out ${vf_result} \
        --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore \
        --evalue ${EVALUE} \
        --id ${MIN_IDENTITY} \
        --query-cover ${MIN_QUERY_COVER} \
        --max-target-seqs 1 \
        --sensitive \
        --max-hsps 1 \
        --threads ${THREADS}

    # Add header
    add_header "${vf_result}"

    # Extract VF-containing contigs
    extract_contigs "${vf_result}" "${contig_dir}/${SAMPLE_ID}_contig.fa" \
                    "${output_dir}/08_VF_contigs/${SAMPLE_ID}_vf_contigs.fa" "vf"

    log "  -> VF annotation complete"
}

# ================== HELPER FUNCTIONS ==================

add_header() {
    local result_file=$1
    local header_file="${result_file}.tmp"
    local header="qseqid\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore"

    if [ -s "${result_file}" ]; then
        echo -e "${header}" > "${header_file}"
        cat "${result_file}" >> "${header_file}"
        mv "${header_file}" "${result_file}"
    else
        echo -e "${header}" > "${result_file}"
    fi
}

extract_contigs() {
    local result_file=$1
    local contig_file=$2
    local output_file=$3
    local type=$4

    # Check if there are valid results
    if [ $(tail -n +2 "${result_file}" | wc -l) -eq 0 ]; then
        log "  -> No ${type} contigs to extract for ${SAMPLE_ID}"
        return
    fi

    # Extract ORF IDs and convert to contig IDs
    tail -n +2 "${result_file}" | \
        cut -f1 | \
        sed -E 's/_[^_]+$//' | \
        sort -u > "${SAMPLE_ID}_${type}_ctg.tmp"

    if [ -s "${SAMPLE_ID}_${type}_ctg.tmp" ]; then
        seqtk subseq "${contig_file}" "${SAMPLE_ID}_${type}_ctg.tmp" > "${output_file}"
        log "  -> Extracted $(wc -l < ${SAMPLE_ID}_${type}_ctg.tmp) ${type} contigs"
    fi

    rm -f "${SAMPLE_ID}_${type}_ctg.tmp"
}

merge_contigs() {
    local type=$1
    local output_dir=$2

    if [ "$type" = "ARG" ]; then
        local in_dir="${output_dir}/02_ARG_contigs"
        local out_dir="${output_dir}/03_ARG_total"
    elif [ "$type" = "MGE" ]; then
        local in_dir="${output_dir}/05_MGE_contigs"
        local out_dir="${output_dir}/06_MGE_total"
    elif [ "$type" = "VF" ]; then
        local in_dir="${output_dir}/08_VF_contigs"
        local out_dir="${output_dir}/09_VF_total"
    fi

    mkdir -p "${out_dir}"
    local out_file="${out_dir}/all_${type,,}_contigs.fa"
    > "${out_file}"

    for fa in "${in_dir}"/*.fa; do
        [ -f "${fa}" ] && [ -s "${fa}" ] || continue
        local fa_basename=$(basename "${fa}" .fa)
        local sid=$(echo "${fa_basename}" | sed -E "s/_${type,,}_contigs$//")
        sed "s/^>/>${sid}_/" "${fa}" >> "${out_file}"
    done

    log "  -> Merged ${type} contigs: ${out_file}"
}

merge_all_results() {
    local output_dir=$1

    log "[Merging] Combining all sample results..."

    local types=("ARG" "MGE" "VF")
    local pattern_suffixes=("diamond_blastp_arg.txt" "diamond_blastp_mge.txt" "diamond_blastp_vf.txt")

    for i in "${!types[@]}"; do
        local type="${types[$i]}"
        local suffix="${pattern_suffixes[$i]}"
        local out_dir="${output_dir}/00_${type}_ALL"

        mkdir -p "${out_dir}"

        case "$type" in
            "ARG") local src_dir="${output_dir}/01_ARG_out" ;;
            "MGE") local src_dir="${output_dir}/04_MGE_out" ;;
            "VF") local src_dir="${output_dir}/07_VF_out" ;;
        esac

        # Merge blast results
        python3 merge_blast_results.py \
            --source-dir "${src_dir}" \
            --output-dir "${out_dir}" \
            --suffix "${suffix}" \
            --sample-prefix "${SAMPLE_ID}_"

        merge_contigs "${type}" "${output_dir}"
    done

    log "  -> All results merged"
}

# ================== MAIN ENTRY POINT ==================

print_usage() {
    cat << EOF
Usage: bash 03_Contigs_annotation.sh -i <protein_dir> -c <contig_dir> -o <output_dir> -a <sarg_db> -m <mge_db> -v <vf_db> -s <sample_id>

Description:
    Annotate contigs for ARGs, MGEs, and VFs using DIAMOND against
    reference databases. Extract annotated contigs for downstream analysis.

Required Arguments:
  -i, --protein-dir   Directory with Prodigal protein files (*.faa)
  -c, --contig-dir    Directory with contig FASTA files (*_contig.fa)
  -o, --output        Output directory
  -a, --sarg-db       Path to SARG database (DIAMOND format)
  -m, --mge-db        Path to MobileOG database (DIAMOND format)
  -v, --vf-db         Path to VFDB database (DIAMOND format)
  -s, --sample        Sample ID

Optional Arguments:
  -t, --threads       Number of threads (default: 32)
  --min-id            Minimum identity (default: 80)
  --min-cover         Minimum query coverage (default: 70)
  --evalue            E-value threshold (default: 1e-7)

Example:
  bash 03_Contigs_annotation.sh \\
      -i results/prodigal \\
      -c results/contigs \\
      -o results/annotation \\
      -a /path/to/SARG.dmnd \\
      -m /path/to/MobileOG.dmnd \\
      -v /path/to/VFDB.dmnd \\
      -s sample1

Input:
  - protein_dir/*.faa: Prodigal-predicted protein sequences
  - contig_dir/*_contig.fa: Assembled contigs

Output:
  - 01_ARG_out/: ARG DIAMOND results
  - 02_ARG_contigs/: Extracted ARG-containing contigs
  - 04_MGE_out/: MGE DIAMOND results
  - 05_MGE_contigs/: Extracted MGE-containing contigs
  - 07_VF_out/: VF DIAMOND results
  - 08_VF_contigs/: Extracted VF-containing contigs
  - 00_*_ALL/: Merged results across samples

Installation:
  conda create -n annotation -c bioconda -c conda-forge \\
      diamond=2.0.15 prodigal=2.6.3 seqtk=1.3 python=3.8

Databases:
  - SARG: https://smile.hku.hk/ARGs/Indexing/download
  - MobileOG: https://mobileogdb.flsi.cloud.vt.edu/entries/database_download
  - VFDB: http://www.mgc.ac.cn/VFs/
EOF
}

main() {
    local protein_dir=""
    local contig_dir=""
    local output_dir=""
    local sarg_db=""
    local mge_db=""
    local vf_db=""
    local sample_id=""

    # Parse arguments
    while [[ $# -gt 0 ]]; do
        case $1 in
            -i|--protein-dir) protein_dir="$2"; shift 2 ;;
            -c|--contig-dir) contig_dir="$2"; shift 2 ;;
            -o|--output) output_dir="$2"; shift 2 ;;
            -a|--sarg-db) sarg_db="$2"; shift 2 ;;
            -m|--mge-db) mge_db="$2"; shift 2 ;;
            -v|--vf-db) vf_db="$2"; shift 2 ;;
            -s|--sample) sample_id="$2"; shift 2 ;;
            -t|--threads) THREADS="$2"; shift 2 ;;
            --min-id) MIN_IDENTITY="$2"; shift 2 ;;
            --min-cover) MIN_QUERY_COVER="$2"; shift 2 ;;
            --evalue) EVALUE="$2"; shift 2 ;;
            -h|--help) print_usage; exit 0 ;;
            *) echo "Unknown option: $1"; print_usage; exit 1 ;;
        esac
    done

    # Validate required arguments
    if [ -z "$protein_dir" ] || [ -z "$contig_dir" ] || [ -z "$output_dir" ] || \
       [ -z "$sarg_db" ] || [ -z "$mge_db" ] || [ -z "$vf_db" ] || [ -z "$sample_id" ]; then
        echo "Error: Missing required arguments!"
        print_usage
        exit 1
    fi

    export SAMPLE_ID="${sample_id}"

    log "=========================================="
    log "Contig Annotation Pipeline"
    log "ARG, MGE, VF Annotation"
    log "=========================================="
    log "Sample ID: ${SAMPLE_ID}"
    log "Protein directory: ${protein_dir}"
    log "Contig directory: ${contig_dir}"
    log "Output directory: ${output_dir}"
    log "Threads: ${THREADS}"
    log "Min identity: ${MIN_IDENTITY}%"
    log "Min query coverage: ${MIN_QUERY_COVER}%"
    log "E-value: ${EVALUE}"
    log "=========================================="

    # Run annotations
    run_arg_annotation ${protein_dir} ${contig_dir} ${output_dir} ${sarg_db}
    run_mge_annotation ${protein_dir} ${contig_dir} ${output_dir} ${mge_db}
    run_vf_annotation ${protein_dir} ${contig_dir} ${output_dir} ${vf_db}

    log "=========================================="
    log "Pipeline completed successfully!"
    log "=========================================="
}

main "$@"