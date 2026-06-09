#!/bin/bash
# ============================================================
# Batch Contig Annotation Pipeline
# For processing multiple samples in parallel
# ============================================================

set -e

# ================== CONFIGURATION ==================
# Edit these paths
PROTEIN_DIR="results/prodigal"          # Prodigal *.faa files
CONTIG_DIR="results/contigs"            # Contig *_contig.fa files
OUTPUT_DIR="results/annotation"         # Output directory
SARG_DB="/path/to/SARG.dmnd"            # SARG database
MGE_DB="/path/to/MobileOG.dmnd"         # MobileOG database
VF_DB="/path/to/VFDB.dmnd"              # VFDB database
SAMPLE_LIST="samples.txt"               # List of sample IDs
THREADS_PER_SAMPLE=16                   # Threads per sample
# ====================================================

log() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $1"
}

# Check if sample list exists, otherwise use all files in protein_dir
if [ ! -f "${SAMPLE_LIST}" ]; then
    log "Sample list not found, using all protein files..."
    SAMPLES=$(ls ${PROTEIN_DIR}/*.faa 2>/dev/null | xargs -n1 basename | sed 's/.faa//' | sort)
else
    SAMPLES=$(cat ${SAMPLE_LIST} | grep -v "^#" | grep -v "^$" | tr -d '\r')
fi

TOTAL=$(echo "$SAMPLES" | wc -l)
COUNT=0

log "=========================================="
log "Batch Contig Annotation Pipeline"
log "=========================================="
log "Total samples: ${TOTAL}"
log "Output directory: ${OUTPUT_DIR}"
log "=========================================="

# Create output directories
mkdir -p ${OUTPUT_DIR}

# Process each sample
for SAMPLE_ID in $SAMPLES; do
    COUNT=$((COUNT + 1))
    log "[${COUNT}/${TOTAL}] Processing sample: ${SAMPLE_ID}"

    bash 03_Contigs_annotation.sh \
        -i ${PROTEIN_DIR} \
        -c ${CONTIG_DIR} \
        -o ${OUTPUT_DIR} \
        -a ${SARG_DB} \
        -m ${MGE_DB} \
        -v ${VF_DB} \
        -s ${SAMPLE_ID} \
        -t ${THREADS_PER_SAMPLE}

    if [ $? -eq 0 ]; then
        log "[${COUNT}/${TOTAL}] Sample ${SAMPLE_ID} completed"
    else
        log "ERROR: Sample ${SAMPLE_ID} failed"
    fi
done

# Merge all results after all samples complete
log "=========================================="
log "Merging all sample results..."
log "=========================================="

for type in ARG MGE VF; do
    case "$type" in
        "ARG") src_dir="${OUTPUT_DIR}/01_ARG_out"; out_dir="${OUTPUT_DIR}/00_ARG_ALL" ;;
        "MGE") src_dir="${OUTPUT_DIR}/04_MGE_out"; out_dir="${OUTPUT_DIR}/00_MGE_ALL" ;;
        "VF")  src_dir="${OUTPUT_DIR}/07_VF_out";  out_dir="${OUTPUT_DIR}/00_VF_ALL" ;;
    esac

    mkdir -p "${out_dir}"

    python3 merge_blast_results.py \
        --source-dir "${src_dir}" \
        --output-dir "${out_dir}" \
        --suffix "_diamond_blastp_${type,,}.txt" \
        --sample-prefix ""

    # Merge contigs
    for type2 in ARG MGE VF; do
        case "$type2" in
            "ARG") in_dir="${OUTPUT_DIR}/02_ARG_contigs"; out_d="${OUTPUT_DIR}/03_ARG_total" ;;
            "MGE") in_dir="${OUTPUT_DIR}/05_MGE_contigs"; out_d="${OUTPUT_DIR}/06_MGE_total" ;;
            "VF")  in_dir="${OUTPUT_DIR}/08_VF_contigs";  out_d="${OUTPUT_DIR}/09_VF_total" ;;
        esac
    done
done

log "=========================================="
log "Batch processing complete!"
log "=========================================="
log "Results:"
log "  ARG: ${OUTPUT_DIR}/00_ARG_ALL/all_samples_merged_arg.csv"
log "  MGE: ${OUTPUT_DIR}/00_MGE_ALL/all_samples_merged_mge.csv"
log "  VF:  ${OUTPUT_DIR}/00_VF_ALL/all_samples_merged_vf.csv"
log "  ARG contigs: ${OUTPUT_DIR}/03_ARG_total/all_arg_contigs.fa"
log "  MGE contigs: ${OUTPUT_DIR}/06_MGE_total/all_mge_contigs.fa"
log "  VF contigs:  ${OUTPUT_DIR}/09_VF_total/all_vf_contigs.fa"