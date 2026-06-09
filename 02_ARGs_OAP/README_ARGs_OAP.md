# ARG Annotation Pipeline using OAP/SARG

## Method Summary

Antibiotic resistance genes (ARGs) were identified using OAP (Online Antibiotic Resistance Platform) against the SARG (Structured Antibiotic Resistance Gene) database.

> Antibiotic resistance genes (ARGs) were identified using OAP (v1.0) against the SARG database (v3.2). ARG-like reads were first identified through DIAMOND alignment, followed by assembly and read mapping for accurate quantification at subtype, type, and class levels.

```
FASTQ Files -> OAP Stage 1 -> OAP Stage 2 -> ARG Abundance Tables
                          |                  |
                   ARG-like reads      Assembled ARGs + Quantification
```

---

## Dependencies

### Software Requirements

| Tool | Version | Purpose |
|------|---------|---------|
| OAP (args_oap) | 1.0+ | ARG annotation pipeline |
| DIAMOND | 2.0+ | Sequence alignment (included with OAP) |

### Installation via Conda

```bash
conda create -n args_oap -c bioconda -c conda-forge \
    args_oap=1.0 \
    diamond=2.0.15

conda activate args_oap
```

### Database Download

```bash
# Download SARG database
# Visit: https://smileassistant1.github.io/OAP/#!/download

# Or use wget (example for v3.2)
wget https://sourceforge.net/projects/oap/files/SARG-3.2.zip
unzip SARG-3.2.zip

# Required files:
# - SARG_v3.2_Short_subdatabase.fasta (for short reads)
# - SARG_v3.2_Short_subdatabase.map
```

**Database versions:**
| Version | Size | Description |
|---------|------|-------------|
| SARG v3.2 (Full) | ~500 MB | Complete ARG reference database |
| SARG v3.2 (Short) | ~50 MB | Optimized for short reads |

---

## Pipeline Steps

### Stage 1: Identify ARG-like Reads

Screen all reads against SARG database using DIAMOND.

```bash
#!/bin/bash
# oap_stage1.sh

FASTQ_DIR="results/clean_data"
OUTPUT_DIR="results/ARGs/stage_one"
SARG_DB="/path/to/SARG_v3.2_Short_subdatabase.fasta"
THREADS=32

mkdir -p ${OUTPUT_DIR}

args_oap stage_one \
    -i ${FASTQ_DIR} \
    -o ${OUTPUT_DIR} \
    -f fq.gz \
    -t ${THREADS} \
    --database ${SARG_DB}
```

**Parameters:**
| Parameter | Description |
|-----------|-------------|
| `-i` | Input directory with FASTQ files |
| `-o` | Output directory |
| `-f` | FASTQ file extension (fq.gz, fastq.gz) |
| `-t` | Number of threads |
| `--database` | Path to SARG database FASTA |

**Stage 1 Output:**
- `*_ORF`: Nucleotide sequences of ARG-like reads
- `*_diamond`: Raw DIAMOND alignment results
- `*_extended_diamond`: Extended alignment with homologs
- `*_normalized`: Normalized ARG counts per sample

---

### Stage 2: Assemble and Quantify

Assemble ARG fragments and calculate abundance.

```bash
#!/bin/bash
# oap_stage2.sh

STAGE1_OUTPUT="results/ARGs/stage_one"
OUTPUT_DIR="results/ARGs/stage_two"
THREADS=32

mkdir -p ${OUTPUT_DIR}

args_oap stage_two \
    -i ${STAGE1_OUTPUT} \
    -t ${THREADS} \
    -o ${OUTPUT_DIR}
```

**Stage 2 Output:**
- `*ARG_abundance_by_subtype.xls`: ARG subtypes (e.g., tetQ, strB)
- `*ARG_abundance_by_type.xls`: ARG types (e.g., tetracycline, aminoglycoside)
- `*ARG_abundance_by_module.xls`: ARG modules
- `*ARG_abundance_by_class.xls`: ARG classes (e.g., antibiotic resistance)

---

## Complete Workflow Script

```bash
#!/bin/bash
# 02_ARGs_OAP.sh

set -e

FASTQ_DIR="results/clean_data"
OUTPUT_DIR="results/ARGs"
SARG_DB="/path/to/SARG_v3.2_Short_subdatabase.fasta"
THREADS=32
FASTQ_EXT="fq.gz"

echo "=========================================="
echo "ARG Annotation Pipeline (OAP/SARG)"
echo "=========================================="

# Stage 1: Identify ARG-like reads
echo "[1/2] OAP Stage 1: Identifying ARG-like reads..."
mkdir -p ${OUTPUT_DIR}/stage_one

args_oap stage_one \
    -i ${FASTQ_DIR} \
    -o ${OUTPUT_DIR}/stage_one \
    -f ${FASTQ_EXT} \
    -t ${THREADS} \
    --database ${SARG_DB}

# Stage 2: Assemble and quantify
echo "[2/2] OAP Stage 2: Assembling and quantifying ARGs..."
mkdir -p ${OUTPUT_DIR}/stage_two

args_oap stage_two \
    -i ${OUTPUT_DIR}/stage_one \
    -t ${THREADS} \
    -o ${OUTPUT_DIR}/stage_two

echo "=========================================="
echo "Pipeline complete!"
echo "=========================================="
echo "Results: ${OUTPUT_DIR}/stage_two/"
```

---

## Usage

```bash
# Standard usage
bash 02_ARGs_OAP.sh \
    -i results/clean_data \
    -o results/ARGs \
    -d /path/to/SARG_v3.2_Short_subdatabase.fasta

# With custom parameters
bash 02_ARGs_OAP.sh -i fastq/ -o args_output -d sarg_db.fasta -t 64 -f fastq.gz
```

---

## Output Structure

```
results/ARGs/
├── stage_one/
│   ├── sample1_ORF.fasta           # ARG-like nucleotide sequences
│   ├── sample1_diamond.m8          # DIAMOND alignments
│   ├── sample1_extended_diamond.m8 # Extended alignments
│   └── sample1_normalized.txt      # Normalized counts
├── stage_two/
│   ├── ARG_abundance_by_subtype.xls  # Subtype level
│   ├── ARG_abundance_by_type.xls     # Type level
│   ├── ARG_abundance_by_module.xls   # Module level
│   └── ARG_abundance_by_class.xls    # Class level
└── summary_report.txt
```

---

## Output Format

### Subtype Abundance Table

```tsv
ARG_subtype	Sample1	Sample2	Sample3
tetQ	125.45	98.32	156.78
strB	78.34	112.56	89.12
sul1	456.23	389.45	423.67
```

### ARG Classification Hierarchy

```
Class (e.g., "Antibiotic Resistance")
    │
    ├── Module (e.g., "Tetracycline")
    │       │
    │       └── Type (e.g., "Tetracycline resistance")
    │               │
    │               └── Subtype (e.g., "tetQ", "tetM")
```

---

## ARG Categories in SARG

| Class | Examples |
|-------|----------|
| Aminoglycoside | strA, strB, aadA |
| Beta-lactam | blaTEM, blaCTX-M, blaOXA |
| Tetracycline | tetQ, tetM, tetW |
| Sulfonamide | sul1, sul2, sul3 |
| Macrolide | ermB, ermF, mefA |
| Glycopeptide | vanA, vanB |
| Quinolone | qnrA, qnrB, qepA |
| Colistin | mcr-1, mcr-2 |

---

## Expected Runtime

| Samples | Estimated Time |
|---------|----------------|
| 10 samples | ~2-4 hours |
| 50 samples | ~8-12 hours |
| 100 samples | ~24-48 hours |

---

## Downstream Analysis

### Visualization in R

```R
library(ggplot2)
library(pheatmap)

# Load abundance data
abund <- read.delim("ARG_abundance_by_type.xls", row.names=1)

# Barplot of top ARGs
top_args <- sort(rowSums(abund), decreasing=TRUE)[1:20]
barplot(top_args, las=2, cex.names=0.5)

# Heatmap of samples vs ARG types
pheatmap(log10(abund + 1),
         clustering_distance_rows="bray",
         clustering_distance_cols="bray",
         cutree_rows=4,
         cutree_cols=4)
```

### Statistical Analysis

```R
library(vegan)

# Beta diversity
bray <- vegdist(t(abund), method="bray")
pcoa <- cmdscale(bray, k=2)

# Ordination with environmental variables
plot(pcoa, pch=21, bg=group, cex=2)
```

---

## Troubleshooting

### OAP not found
```bash
# Check installation
which args_oap

# Reinstall if needed
conda install -c bioconda args_oap
```

### Database errors
```bash
# Verify database files exist
ls -la /path/to/SARG*.fasta
ls -la /path/to/SARG*.map

# Check file permissions
chmod 644 /path/to/SARG*.fasta
```

### Memory issues
```bash
# Reduce thread count
# DIAMOND uses ~8GB per thread

# Or use smaller SARG database
# Short subdatabase is optimized for memory
```

---

## Citation

If you use this pipeline, please cite:

1. Yin X et al. OAP: A pipeline for detecting and quantifying antibiotic resistance genes in metagenomes. Preprint.
2. Arango-Argoty G et al. The structured antibiotic resistance genes (SARG) database. *Nucleic Acids Res*. 2018.

---

## Data Availability

```
SARG database: https://smileassistant1.github.io/OAP/#!/download
```