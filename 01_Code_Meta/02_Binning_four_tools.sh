#!/bin/bash
# ============================================================
# Metagenome Binning Pipeline - Four Tools + DAS_Tool
# For Nature article supplementary materials
#
# Method: Metagenome binning was conducted using four independent
# tools: MetaBAT2 (v2.12.1), MaxBin2 (v2.2.6), CONCOCT (v1.1.0),
# and SemiBin2 (v2.2). Only contigs with lengths >= 2000 bp were
# subjected to binning, generating preliminary MAGs. These MAGs
# were integrated and dereplicated via DAS_Tool (v1.1.6) to obtain
# non-redundant MAGs.
#
# Note: MetaBAT2, MaxBin2, and CONCOCT are executed via metaWRAP.
#
# Version: 1.0
# Dependencies: metaWRAP, SemiBin2, DAS_Tool, Bowtie2, SAMtools, seqkit
#
# Usage:
#   bash 02_Binning_four_tools.sh -a <assembly.fa> -1 <R1.fq.gz> -2 <R2.fq.gz> -o <output_dir>
# ============================================================

set -e

# ================== CONFIGURATION ==================
THREADS=32
MIN_CONTIG_LENGTH=2000
WORK_DIR="results"
# ============================================================

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
    local deps=("metawrap" "SemiBin2" "DAS_Tool" "bowtie2" "samtools" "seqkit")
    for dep in "${deps[@]}"; do
        check_software ${dep}
    done
    log "All dependencies found."
}

# ================== PREPROCESSING ==================
prepare_input() {
    local assembly=$1
    local r1=$2
    local r2=$3
    local work_dir=$4

    log "[Preprocessing] Preparing input files..."

    local preprocess_dir="${work_dir}/preprocess"
    mkdir -p ${preprocess_dir}

    # Filter contigs by length (>= MIN_CONTIG_LENGTH bp)
    log "  [1] Filtering contigs by length (>=${MIN_CONTIG_LENGTH} bp)..."
    seqkit seq -m ${MIN_CONTIG_LENGTH} ${assembly} -o ${preprocess_dir}/contigs_${MIN_CONTIG_LENGTH}bp.fa

    local contig_count=$(grep -c "^>" ${preprocess_dir}/contigs_${MIN_CONTIG_LENGTH}bp.fa)
    log "  -> Retained ${contig_count} contigs (>=${MIN_CONTIG_LENGTH} bp)"

    # Build Bowtie2 index
    log "  [2] Building Bowtie2 index..."
    mkdir -p ${preprocess_dir}/bt2_index
    bowtie2-build ${preprocess_dir}/contigs_${MIN_CONTIG_LENGTH}bp.fa \
        ${preprocess_dir}/bt2_index/contigs

    # Map reads to contigs
    log "  [3] Mapping reads to contigs..."
    bowtie2 -p ${THREADS} \
        -x ${preprocess_dir}/bt2_index/contigs \
        -1 ${r1} -2 ${r2} \
        --very-sensitive \
        -S ${preprocess_dir}/aligned.sam

    samtools view -bS -@ ${THREADS} ${preprocess_dir}/aligned.sam | \
        samtools sort -@ ${THREADS} -o ${preprocess_dir}/aligned_sorted.bam

    samtools index ${preprocess_dir}/aligned_sorted.bam

    log "  -> Preprocessing complete!"
    log "  -> Contigs: ${preprocess_dir}/contigs_${MIN_CONTIG_LENGTH}bp.fa"
    log "  -> BAM: ${preprocess_dir}/aligned_sorted.bam"

    echo "${preprocess_dir}/contigs_${MIN_CONTIG_LENGTH}bp.fa"
    echo "${preprocess_dir}/aligned_sorted.bam"
}

# ================== STEP 1: MetaWRAP Binning (MetaBAT2, MaxBin2, CONCOCT) ==================
run_metawrap_binning() {
    local contigs=$1
    local r1=$2
    local r2=$3
    local work_dir=$4

    log "[Step 1/3] MetaWRAP binning (MetaBAT2, MaxBin2, CONCOCT)..."

    local metawrap_dir="${work_dir}/metawrap_binning"
    mkdir -p ${metawrap_dir}

    # Run metaWRAP binning with all three tools
    # metaWRAP handles read mapping, coverage calculation, and binning internally
    metawrap binning \
        -t ${THREADS} \
        -a ${contigs} \
        -o ${metawrap_dir}/raw_bins \
        --metabat2 \
        --maxbin2 \
        --concoct \
        -1 ${r1} \
        -2 ${r2}

    local metabat_count=$(find ${metawrap_dir}/raw_bins/metabat2 -name '*.fa' 2>/dev/null | wc -l)
    local maxbin_count=$(find ${metawrap_dir}/raw_bins/maxbin2 -name '*.fa' 2>/dev/null | wc -l)
    local concoct_count=$(find ${metawrap_dir}/raw_bins/concoct -name '*.fa' 2>/dev/null | wc -l)

    log "  -> MetaBAT2 bins: ${metabat_count}"
    log "  -> MaxBin2 bins: ${maxbin_count}"
    log "  -> CONCOCT bins: ${concoct_count}"
}

# ================== STEP 2: SemiBin2 (Independent) ==================
run_semibin2() {
    local contigs=$1
    local r1=$2
    local r2=$3
    local bam=$4
    local work_dir=$5

    log "[Step 2/3] Binning with SemiBin2 v2.2..."

    local semibin_dir="${work_dir}/semibin2"
    mkdir -p ${semibin_dir}

    # Run SemiBin2 in single-sample mode
    SemiBin2 single \
        -i ${contigs} \
        -b ${bam} \
        -o ${semibin_dir} \
        -p ${THREADS} \
        --min-contig-length ${MIN_CONTIG_LENGTH}

    local bin_count=$(find ${semibin_dir} -name "bin_*.fa" 2>/dev/null | wc -l)
    log "  -> SemiBin2 bins: ${bin_count}"
}

# ================== STEP 3: DAS_Tool Integration ==================
run_dastool() {
    local work_dir=$1
    local contigs=$2

    log "[Step 3/3] Integrating and dereplicating with DAS_Tool v1.1.6..."

    local dastool_dir="${work_dir}/dastool"
    mkdir -p ${dastool_dir}

    # Rename SemiBin2 bins to match DAS_Tool expected format
    local semibin_renamed="${work_dir}/semibin2_renamed"
    mkdir -p ${semibin_renamed}

    local i=1
    for bin in ${work_dir}/semibin2/bin_*.fa; do
        if [ -f "$bin" ]; then
            cp "$bin" "${semibin_renamed}/bin.${i}.fa"
            i=$((i+1))
        fi
    done

    # Create bin collection file (format: tool<tab>path_to_bins)
    cat > ${dastool_dir}/bin_collection.tsv << EOF
MetaBAT2	${work_dir}/metawrap_binning/raw_bins/metabat2
MaxBin2	${work_dir}/metawrap_binning/raw_bins/maxbin2
CONCOCT	${work_dir}/metawrap_binning/raw_bins/concoct
SemiBin2	${semibin_renamed}
EOF

    # Run DAS_Tool
    DAS_Tool \
        -i ${dastool_dir}/bin_collection.tsv \
        -l MetaBAT2,MaxBin2,CONCOCT,SemiBin2 \
        -c ${contigs} \
        -o ${dastool_dir}/DAS_Tool \
        --write_bins \
        -p ${THREADS}

    # Summary statistics
    local total_preliminary=0

    local metabat_count=$(find ${work_dir}/metawrap_binning/raw_bins/metabat2 -name '*.fa' 2>/dev/null | wc -l)
    local maxbin_count=$(find ${work_dir}/metawrap_binning/raw_bins/maxbin2 -name '*.fa' 2>/dev/null | wc -l)
    local concoct_count=$(find ${work_dir}/metawrap_binning/raw_bins/concoct -name '*.fa' 2>/dev/null | wc -l)
    local semibin_count=$(find ${work_dir}/semibin2 -name 'bin_*.fa' 2>/dev/null | wc -l)

    total_preliminary=$((metabat_count + maxbin_count + concoct_count + semibin_count))

    local derep_bins=$(find ${dastool_dir}/DAS_Tool_DASToolbins -name '*.fa' 2>/dev/null | wc -l)

    log "  === Binning Summary ==="
    log "  -> MetaBAT2: ${metabat_count} bins"
    log "  -> MaxBin2: ${maxbin_count} bins"
    log "  -> CONCOCT: ${concoct_count} bins"
    log "  -> SemiBin2: ${semibin_count} bins"
    log "  -> Total preliminary MAGs: ${total_preliminary}"
    log "  -> Non-redundant MAGs (after DAS_Tool): ${derep_bins}"
    log "  -> Final bins: ${dastool_dir}/DAS_Tool_DASToolbins/"
}

# ================== GENERATE SUMMARY ==================
generate_summary() {
    local work_dir=$1

    log "=========================================="
    log "Binning Pipeline Complete!"
    log "=========================================="
    log ""
    log "Configuration:"
    log "  Min contig length: >=${MIN_CONTIG_LENGTH} bp"
    log "  Threads: ${THREADS}"
    log ""
    log "Individual tool results:"
    log "  - MetaBAT2:  $(find ${work_dir}/metawrap_binning/raw_bins/metabat2 -name '*.fa' 2>/dev/null | wc -l) bins"
    log "  - MaxBin2:   $(find ${work_dir}/metawrap_binning/raw_bins/maxbin2 -name '*.fa' 2>/dev/null | wc -l) bins"
    log "  - CONCOCT:   $(find ${work_dir}/metawrap_binning/raw_bins/concoct -name '*.fa' 2>/dev/null | wc -l) bins"
    log "  - SemiBin2:  $(find ${work_dir}/semibin2 -name 'bin_*.fa' 2>/dev/null | wc -l) bins"
    log ""
    log "Dereplicated MAGs:"
    log "  - DAS_Tool:  $(find ${work_dir}/dastool/DAS_Tool_DASToolbins -name '*.fa' 2>/dev/null | wc -l) non-redundant MAGs"
    log ""
    log "Output directories:"
    log "  - ${work_dir}/metawrap_binning/raw_bins/metabat2/"
    log "  - ${work_dir}/metawrap_binning/raw_bins/maxbin2/"
    log "  - ${work_dir}/metawrap_binning/raw_bins/concoct/"
    log "  - ${work_dir}/semibin2/"
    log "  - ${work_dir}/dastool/DAS_Tool_DASToolbins/"
    log "=========================================="
}

# ================== MAIN ENTRY POINT ==================

print_usage() {
    cat << EOF
Usage: bash 02_Binning_four_tools.sh -a <assembly.fa> -1 <R1.fq.gz> -2 <R2.fq.gz> -o <output_dir>

Description:
    Metagenome binning using four independent tools (MetaBAT2, MaxBin2,
    CONCOCT via metaWRAP, and SemiBin2), integrated and dereplicated
    via DAS_Tool.

Required Arguments:
  -a, --assembly    Assembly contigs file (FASTA)
  -1, --read1       Forward reads (FASTQ)
  -2, --read2       Reverse reads (FASTQ)
  -o, --output      Output directory

Optional Arguments:
  -t, --threads     Number of threads (default: 32)
  -m, --min-length  Minimum contig length for binning (default: 2000)

Example:
  bash 02_Binning_four_tools.sh \\
      -a results/assembly/contigs_1000bp.fa \\
      -1 results/clean/sample_R1.fastq.gz \\
      -2 results/clean/sample_R2.fastq.gz \\
      -o results/binning

  # With custom parameters
  bash 02_Binning_four_tools.sh -a assembly/contigs.fasta \\
                                -1 clean/R1.fq.gz \\
                                -2 clean/R2.fq.gz \\
                                -o binning \\
                                -t 64 \\
                                -m 2500

Input:
  - assembly.fa: Contigs from metagenomic assembly (>=1000 bp recommended)
  - R1/R2.fastq.gz: Quality-filtered, cleaned reads

Output:
  - Four binning methods (MetaBAT2, MaxBin2, CONCOCT, SemiBin2)
  - DAS_Tool dereplicated non-redundant MAGs
EOF
}

main() {
    local assembly=""
    local read1=""
    local read2=""
    local output_dir=""

    # Parse arguments
    while [[ $# -gt 0 ]]; do
        case $1 in
            -a|--assembly) assembly="$2"; shift 2 ;;
            -1|--read1) read1="$2"; shift 2 ;;
            -2|--read2) read2="$2"; shift 2 ;;
            -o|--output) output_dir="$2"; shift 2 ;;
            -t|--threads) THREADS="$2"; shift 2 ;;
            -m|--min-length) MIN_CONTIG_LENGTH="$2"; shift 2 ;;
            -h|--help) print_usage; exit 0 ;;
            *) echo "Unknown option: $1"; print_usage; exit 1 ;;
        esac
    done

    # Validate required arguments
    if [ -z "$assembly" ] || [ -z "$read1" ] || [ -z "$read2" ] || [ -z "$output_dir" ]; then
        echo "Error: Missing required arguments!"
        print_usage
        exit 1
    fi

    log "=========================================="
    log "Metagenome Binning Pipeline"
    log "MetaBAT2, MaxBin2, CONCOCT, SemiBin2 + DAS_Tool"
    log "=========================================="
    log "Assembly: ${assembly}"
    log "Read 1: ${read1}"
    log "Read 2: ${read2}"
    log "Output: ${output_dir}"
    log "Threads: ${THREADS}"
    log "Min contig length: >=${MIN_CONTIG_LENGTH} bp"
    log "=========================================="

    # Preprocessing
    local result=$(prepare_input ${assembly} ${read1} ${read2} ${output_dir})
    local contigs=$(echo "$result" | head -n1)
    local bam=$(echo "$result" | tail -n1)

    # Step 1: MetaWRAP binning
    run_metawrap_binning ${contigs} ${read1} ${read2} ${output_dir}

    # Step 2: SemiBin2
    run_semibin2 ${contigs} ${read1} ${read2} ${bam} ${output_dir}

    # Step 3: DAS_Tool integration
    run_dastool ${output_dir} ${contigs}

    generate_summary ${output_dir}

    log "Pipeline completed successfully!"
}

main "$@"