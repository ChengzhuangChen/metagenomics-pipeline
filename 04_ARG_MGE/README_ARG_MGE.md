# ARG-MGE Co-localization Analysis Pipeline

## Method Summary

To investigate potential horizontal gene transfer (HGT) of antibiotic resistance genes, ARG-MGE co-localization analysis was performed. **ARGs and MGEs within 10 kb on the same contig were regarded as co-occurring.**

> ARG-MGE co-localization analysis was performed to assess potential horizontal gene transfer. Only ARGs and MGEs located within 10 kb of each other on the same contig were considered as co-occurring.

```
all_samples_merged_arg.csv + all_samples_merged_mge.csv
              |
              + GFF files (for position information)
                        |
                        v
              identify_cocol_contigs.py
              (calculates distance between ORFs)
                        |
              ARGs and MGEs within 10 kb = co-occurring
                        |
                        v
              Extract sequences + annotations
```

---

## Key Criterion

**"ARGs and MGEs within 10 kb on the same contig were regarded as co-occurring"**

This means:
1. Both an ARG and an MGE must be on the **same contig**
2. The **distance between them** (measured from ORF midpoints) must be **≤ 10,000 bp (10 kb)**

---

## Dependencies

### Software Requirements

| Tool | Version | Purpose |
|------|---------|---------|
| Python | 3.8+ | Script execution |
| Biopython | 1.79+ | Sequence handling |
| seqtk | 1.3+ | Sequence extraction |

### Installation via Conda

```bash
conda create -n arg_mge -c conda-forge -c bioconda \
    python=3.8 \
    biopython=1.79 \
    seqtk=1.3

conda activate arg_mge
```

---

## Pipeline Steps

### Step 1: Identify Co-localized ORFs (within 10 kb)

```bash
#!/bin/bash
# identify_cocol.sh

python3 identify_cocol_contigs.py \
    --arg-csv results/annotation/00_ARG_ALL/all_samples_merged_arg.csv \
    --mge-csv results/annotation/00_MGE_ALL/all_samples_merged_mge.csv \
    --gff-dir results/prodigal \
    --output cocol_pairs_within_10kb.tsv \
    --stats cocol_stats.txt \
    --distance 10000
```

**How it works:**
1. Load ARG and MGE ORF IDs from DIAMOND CSV results
2. Load ORF positions from GFF files (start, end coordinates)
3. For each contig with both ARGs and MGEs:
   - Calculate midpoint of each ORF
   - Check if any ARG-MGE pair has distance ≤ 10 kb
4. Output all ARG-MGE pairs within the threshold

---

### Step 2: Extract Sequences

```bash
#!/bin/bash
# extract_seqs.sh

seqtk subseq all_arg_contigs.fa cocol_contigs.txt > cocol_contigs.fa
```

---

### Step 3: Extract Annotations

```bash
#!/bin/bash
# extract_annot.sh

python3 extract_cocol_annotations.py \
    --contigs cocol_contigs.txt \
    --faa-dir results/prodigal \
    --gff-dir results/prodigal \
    --output-dir annotations
```

---

## Complete Workflow Script

```bash
#!/bin/bash
# 04_ARG_MGE_colocalization.sh

set -e

ARG_CSV="results/annotation/00_ARG_ALL/all_samples_merged_arg.csv"
MGE_CSV="results/annotation/00_MGE_ALL/all_samples_merged_mge.csv"
GFF_DIR="results/prodigal"
FAA_DIR="results/prodigal"
OUTPUT_DIR="results/arg_mge_analysis"
DISTANCE=10000  # 10 kb

echo "=========================================="
echo "ARG-MGE Co-localization Analysis"
echo "Criterion: ARGs and MGEs within 10 kb"
echo "=========================================="

# Step 1: Identify co-localized ORFs
echo "[1/4] Identifying ARGs and MGEs within ${DISTANCE} bp (10 kb)..."
mkdir -p ${OUTPUT_DIR}/01_cocol_pairs

python3 identify_cocol_contigs.py \
    --arg-csv ${ARG_CSV} \
    --mge-csv ${MGE_CSV} \
    --gff-dir ${GFF_DIR} \
    --output ${OUTPUT_DIR}/01_cocol_pairs/cocol_pairs.tsv \
    --stats ${OUTPUT_DIR}/01_cocol_pairs/cocol_stats.txt \
    --distance ${DISTANCE}

# Step 2: Extract sequences
echo "[2/4] Extracting co-localized contig sequences..."
mkdir -p ${OUTPUT_DIR}/02_cocol_sequences

seqtk subseq all_arg_contigs.fa \
    ${OUTPUT_DIR}/01_cocol_pairs/cocol_pairs_contigs.txt \
    > ${OUTPUT_DIR}/02_cocol_sequences/cocol_contigs.fa

# Step 3: Extract annotations
echo "[3/4] Extracting annotations..."
mkdir -p ${OUTPUT_DIR}/03_annotations

python3 extract_cocol_annotations.py \
    --contigs ${OUTPUT_DIR}/01_cocol_pairs/cocol_pairs_contigs.txt \
    --faa-dir ${FAA_DIR} \
    --gff-dir ${GFF_DIR} \
    --output-dir ${OUTPUT_DIR}/03_annotations

# Step 4: Summary
echo "[4/4] Summary..."
echo "Pipeline complete!"

echo "=========================================="
echo "Key Results:"
echo "  Criterion: ARGs and MGEs within 10 kb"
echo "  Co-localized contigs: $(wc -l < ${OUTPUT_DIR}/01_cocol_pairs/cocol_pairs_contigs.txt)"
echo "=========================================="
```

---

## Usage

```bash
bash 04_ARG_MGE_colocalization.sh \
    -r results/annotation/00_ARG_ALL/all_samples_merged_arg.csv \
    -m results/annotation/00_MGE_ALL/all_samples_merged_mge.csv \
    -g results/prodigal \
    -f results/prodigal \
    -o results/arg_mge_analysis
```

---

## Input Files

| File | Source | Description |
|------|--------|-------------|
| `all_samples_merged_arg.csv` | 03_Contigs | Merged ARG DIAMOND results |
| `all_samples_merged_mge.csv` | 03_Contigs | Merged MGE DIAMOND results |
| `*.gff` | Prodigal | Gene coordinates (start, end) |
| `*.faa` | Prodigal | Protein sequences |

---

## Output Structure

```
results/arg_mge_analysis/
├── 01_cocol_pairs/
│   ├── cocol_pairs_within_10kb.tsv    # All ARG-MGE pairs within 10 kb
│   ├── cocol_pairs_contigs.txt        # List of co-localized contigs
│   └── cocol_stats.txt                # Statistics
├── 02_cocol_sequences/
│   └── cocol_contigs.fa               # Extracted sequences
├── 03_annotations/
│   ├── cocol_proteins.faa
│   ├── cocol_genes.gff
│   └── cocol_annotation_summary.tsv
└── summary_report.txt
```

---

## Output Format

### Co-localization Pairs (TSV)

```tsv
Contig	ARG_ORF	ARG_Start	ARG_End	MGE_ORF	MGE_Start	MGE_End	Distance_bp	Distance_kb
sample1_contig_001	sample1_contig_001_orf3	1000	1800	sample1_contig_001_orf7	6500	7200	4000	4.00
sample1_contig_001	sample1_contig_001_orf3	1000	1800	sample1_contig_001_orf9	8500	9100	7000	7.00
sample1_contig_002	sample1_contig_002_orf5	500	1200	sample1_contig_002_orf8	9000	9800	8000	8.00
```

### Statistics Report

```
ARG-MGE Co-localization Analysis (10 kb window)
============================================================

Distance threshold: 10000 bp (10 kb)

Basic Statistics:
----------------------------------------
Total ARG ORFs detected: 450
Total MGE ORFs detected: 1200
Contigs with ARGs: 150
Contigs with MGEs: 280
Contigs with both ARGs and MGEs: 85

Co-localization Results (within 10 kb):
----------------------------------------
ARG-MGE pairs within 10 kb: 120
Co-localized contigs: 45
Co-localization rate (of ARG contigs): 30.00%
Co-localization rate (of MGE contigs): 16.07%

Distance Statistics:
----------------------------------------
Average distance: 4500 bp (4.50 kb)
Median distance: 3800 bp (3.80 kb)
Min distance: 200 bp
Max distance: 9800 bp
```

---

## Distance Calculation Method

```
Distance = |ARG_midpoint - MGE_midpoint|

Where:
  ARG_midpoint = (ARG_start + ARG_end) / 2
  MGE_midpoint = (MGE_start + MGE_end) / 2
```

**Why midpoint?**
- More robust than edge-to-edge distance
- Better represents the "center" of each ORF
- Standard practice in genomic co-localization analysis

---

## Expected Results

| Metric | Typical Range |
|--------|---------------|
| Distance threshold | 10,000 bp (10 kb) |
| Contigs with both ARGs and MGEs | 50-100 |
| ARG-MGE pairs within 10 kb | 30-150 |
| Co-localization rate (of ARG contigs) | 20-50% |

---

## Downstream Analysis

### ARG Class Distribution

```R
library(ggplot2)

# Load pairs
pairs <- read.delim("cocol_pairs_within_10kb.tsv")

# Parse ARG subtype from ORF ID
pairs$ARG_subtype <- gsub(".*_", "", pairs$ARG_ORF)

# Barplot
ggplot(pairs, aes(x = ARG_subtype)) +
    geom_bar() +
    coord_flip() +
    labs(x = "ARG Subtype", y = "Count",
         title = "ARG Subtypes in Co-localized ARG-MGE Pairs")
```

### Distance Distribution

```R
# Histogram of distances
ggplot(pairs, aes(x = Distance_kb)) +
    geom_histogram(binwidth = 1) +
    labs(x = "Distance (kb)", y = "Count",
         title = "Distance Distribution of ARG-MGE Pairs")
```

---

## Citation

If you use this pipeline, please cite:

1. Partridge SR et al. Mobile genetic elements associated with antibiotic resistance. *Int J Med Microbiol*. 2018.
2. Gillings MR. Integrons: past, present, and future. *Microbiol Mol Biol Rev*. 2014.
3. Zhang AN et al. Host quotas and cooperative hitches shift the evolutionary dynamics of antibiotic resistance on conjugative plasmids. *Mol Biol Evol*. 2022.