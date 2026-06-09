# MAG Abundance Profiling with CoverM

## Method Summary

MAG abundance was calculated using CoverM with multiple normalization methods (RPKM, TPM, counts, coverage, CPM) to quantify relative abundances across samples.

> MAG abundance was calculated using CoverM with RPKM normalization. Mapping was performed with Bowtie2, and reads with less than 95% identity or less than 10% covered fraction were excluded to ensure accurate quantification.

---

## Dependencies

### Software Requirements

| Tool | Version | Purpose |
|------|---------|---------|
| CoverM | 0.6.1+ | MAG abundance calculation |
| Bowtie2 | 2.4+ | Read mapping |
| samtools | 1.14+ | BAM processing |

### Installation via Conda

```bash
conda create -n coverm -c bioconda -c conda-forge \
    coverm=0.6.1 \
    bowtie2=2.4.5 \
    samtools=1.14

conda activate coverm
```

---

## Pipeline Steps

### Step 1: Prepare Sample List

Identify paired-end FASTQ files in the input directory.

```bash
#!/bin/bash
# prepare_samples.sh

FASTQ_DIR="results/clean_data"
OUTPUT_DIR="results/abundance"

mkdir -p ${OUTPUT_DIR}

# Find all R1 files and pair with R2
COUPLED_READS=()
for r1 in ${FASTQ_DIR}/*_R1.fq.gz; do
    sample_name=$(basename ${r1} _R1.fq.gz)
    r2="${FASTQ_DIR}/${sample_name}_R2.fq.gz"

    if [[ -f "${r2}" ]]; then
        COUPLED_READS+=("${r1}")
        COUPLED_READS+=("${r2}")
        echo "Found: ${sample_name}"
    fi
done

echo "Total sample pairs: ${#COUPLED_READS[@]}"
```

---

### Step 2: Run CoverM

Calculate MAG abundance with normalization.

```bash
#!/bin/bash
# coverm_run.sh

MAG_DIR="results/mag_quality/drep/dereplicated_genomes"
OUTPUT_DIR="results/abundance"
THREADS=32

mkdir -p ${OUTPUT_DIR}

# Build coupled reads list
COUPLED_READS=()
for r1 in results/clean_data/*_R1.fq.gz; do
    sample=$(basename ${r1} _R1.fq.gz)
    COUPLED_READS+=("${r1}")
    COUPLED_READS+=("results/clean_data/${sample}_R2.fq.gz")
done

# Run CoverM with RPKM normalization
coverm genome \
    --genome-fasta-directory ${MAG_DIR} \
    --genome-fasta-extension fa \
    -c "${COUPLED_READS[@]}" \
    -m rpkm \
    -t ${THREADS} \
    --min-read-percent-identity 95 \
    --min-covered-fraction 10 \
    -o ${OUTPUT_DIR}/abundance_matrix.tsv

echo "Done!"
```

---

## Complete Workflow Script

```bash
#!/bin/bash
# 05_coverm.sh

set -e

MAG_DIR="results/mag_quality/drep/dereplicated_genomes"
FASTQ_DIR="results/clean_data"
OUTPUT_DIR="results/abundance"
THREADS=32

echo "=========================================="
echo "CoverM MAG Abundance Profiling"
echo "=========================================="

# Prepare sample list
echo "[1/3] Preparing sample list..."
COUPLED_READS=()
for r1 in ${FASTQ_DIR}/*_R1.fq.gz; do
    sample=$(basename ${r1} _R1.fq.gz)
    r2="${FASTQ_DIR}/${sample}_R2.fq.gz"
    if [[ -f "${r2}" ]]; then
        COUPLED_READS+=("${r1}" "${r2}")
        echo "  + ${sample}"
    fi
done

# Run CoverM
echo "[2/3] Running CoverM..."
coverm genome \
    --genome-fasta-directory ${MAG_DIR} \
    --genome-fasta-extension fa \
    -c "${COUPLED_READS[@]}" \
    -m rpkm \
    -t ${THREADS} \
    --min-read-percent-identity 95 \
    --min-covered-fraction 10 \
    -o ${OUTPUT_DIR}/abundance_matrix.tsv

# Generate additional matrices
echo "[3/3] Generating additional matrices..."
for method in tpm counts coverage cpm; do
    coverm genome \
        --genome-fasta-directory ${MAG_DIR} \
        --genome-fasta-extension fa \
        -c "${COUPLED_READS[@]}" \
        -m ${method} \
        -t ${THREADS} \
        --min-read-percent-identity 95 \
        --min-covered-fraction 10 \
        -o ${OUTPUT_DIR}/abundance_matrix_${method}.tsv
done

echo "=========================================="
echo "Pipeline complete!"
echo "=========================================="
```

---

## Usage

```bash
# Standard usage
bash 05_coverm.sh \
    -m results/mag_quality/drep/dereplicated_genomes \
    -f results/clean_data \
    -o results/abundance

# With custom parameters
bash 05_coverm.sh \
    -m bins/ \
    -f fastq/ \
    -o output \
    -t 64 \
    -i 97 \
    -c 20 \
    -n tpm
```

---

## Output Structure

```
results/abundance/
├── abundance_matrix.tsv          # RPKM normalized
├── abundance_matrix_tpm.tsv      # TPM normalized
├── abundance_matrix_counts.tsv   # Raw counts
├── abundance_matrix_coverage.tsv # Coverage values
├── abundance_matrix_cpm.tsv      # Counts per million
└── summary_report.txt
```

---

## Output Format

### Abundance Matrix (TSV)

```tsv
Genome	DB1	DB2	DB3	LAKE1	LAKE2
MAG.001	125.45	98.32	156.78	45.23	67.89
MAG.002	78.34	112.56	89.12	234.56	198.34
MAG.003	456.23	389.45	423.67	123.45	156.78
```

**Columns:**
| Column | Description |
|--------|-------------|
| Genome | MAG identifier |
| Sample1-N | Abundance values for each sample |

---

## Normalization Methods

| Method | Full Name | Description |
|--------|-----------|-------------|
| `rpkm` | Reads Per Kilobase per Million | Normalized by gene length and total reads |
| `tpm` | Transcripts Per Million | Similar to RPKM, different normalization order |
| `counts` | Raw Counts | Total mapped reads per MAG |
| `coverage` | Mean Coverage | Average covered bases per position |
| `cpm` | Counts Per Million | Counts normalized to per million total reads |

**When to use which:**
- **RPKM/TPM**: Cross-study comparison, gene expression-like analysis
- **Counts**: Differential abundance analysis (DESeq2, edgeR)
- **Coverage**: Assembly quality consideration
- **CPM**: Quick relative comparison

---

## Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--min-read-percent-identity` | 95 | Minimum Bowtie2 alignment identity (%) |
| `--min-covered-fraction` | 10 | Minimum genome coverage (%) |
| `--genome-fasta-extension` | fa | File extension for MAGs |

---

## Expected Runtime

| Samples x MAGs | Estimated Time |
|----------------|----------------|
| 10 x 500 | ~30 minutes |
| 50 x 1000 | ~2-3 hours |
| 100 x 5000 | ~4-6 hours |

---

## Downstream Analysis

### Import into R for analysis

```R
library(readr)
library(vegan)
library(ggplot2)

# Load abundance matrix
abund <- read_tsv("abundance_matrix.tsv")
rownames(abund) <- abund$Genome
abund$Genome <- NULL

# Transpose if samples are in columns
abund_t <- t(abund)

# Calculate beta diversity
bray <- vegdist(abund_t, method = "bray")

# PCoA visualization
pcoa <- cmdscale(bray, k = 2)
pcoa_df <- data.frame(PC1 = pcoa[,1], PC2 = pcoa[,2], Sample = rownames(pcoa))
ggplot(pcoa_df, aes(x = PC1, y = PC2, label = Sample)) + geom_point() + geom_text()
```

---

## Troubleshooting

### No reads mapped
```bash
# Check genome format
head -1 results/mag_quality/drep/dereplicated_genomes/*.fa

# Verify FASTQ pairing
ls results/clean_data/*_R1.fq.gz | wc -l
ls results/clean_data/*_R2.fq.gz | wc -l
```

### Slow processing
```bash
# Reduce mapping stringency
# (maps more reads, less accurate)
--min-read-percent-identity 90

# Skip secondary alignments
coverm genome ... --discard-unmapped
```

---

## Citation

If you use this pipeline, please cite:

1. Woodcroft BJ et al. CoverM: fast and accurate calculation of genome coverages. https://github.com/wwood/CoverM