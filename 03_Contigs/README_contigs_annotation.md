# Contig-level Annotation Pipeline (ARG, MGE, VF)

## Method Summary

Protein sequences predicted by Prodigal were aligned against three reference databases using DIAMOND: SARG (antibiotic resistance genes), MobileOG (mobile genetic elements), and VFDB (virulence factors). Contigs harboring annotated ORFs were extracted for downstream analysis.

> Contig-level annotation was performed by aligning Prodigal-predicted proteins against the SARG database (antibiotic resistance genes), MobileOG database (mobile genetic elements), and VFDB (virulence factors) using DIAMOND with 80% identity, 70% query coverage, and 1e-7 E-value thresholds.

```
Contigs -> Prodigal -> Protein sequences -> DIAMOND alignment -> Annotated contigs
                                                          |
                    SARG database -> ARG annotation -------
                    MobileOG -> MGE annotation ------------
                    VFDB -> VF annotation -----------------
```

---

## Dependencies

### Software Requirements

| Tool | Version | Purpose |
|------|---------|---------|
| DIAMOND | 2.0+ | Protein sequence alignment |
| Prodigal | 2.6.3 | Gene prediction |
| seqtk | 1.3+ | Sequence manipulation |
| Python | 3.8+ | Result processing |

### Installation via Conda

```bash
conda create -n annotation -c bioconda -c conda-forge \
    diamond=2.0.15 \
    prodigal=2.6.3 \
    seqtk=1.3 \
    python=3.8

conda activate annotation
```

### Database Download

```bash
# SARG database
# Download from: https://smileassistant1.github.io/OAP/#!/download
wget https://sourceforge.net/projects/oap/files/SARG-3.2.zip
unzip SARG-3.2.zip
diamond makedb --in SARG_v3.2_Short_subdatabase.fasta -d SARG_v3.2_Short_subdatabase

# MobileOG database
wget https://mobileogdb.fresh桔life.org/download
diamond makedb --in mobileOG-db.fasta -d mobileOG-db

# VFDB
wget http://www.mgc.ac.cn/VFs/Down/VFDB_setA_pro.fas
diamond makedb --in VFDB_setA_pro.fas -d VFDB
```

---

## Pipeline Steps

### Step 1: ARG Annotation

```bash
#!/bin/bash
# ARG annotation

PROTEIN_FILE="results/prodigal/sample1.faa"
CONTIG_FILE="results/contigs/sample1_contig.fa"
OUTPUT="results/annotation"
SARG_DB="/path/to/SARG.dmnd"

mkdir -p ${OUTPUT}/01_ARG_out ${OUTPUT}/02_ARG_contigs

# DIAMOND blastp against SARG
diamond blastp \
    --query ${PROTEIN_FILE} \
    --db ${SARG_DB} \
    --out ${OUTPUT}/01_ARG_out/sample1_diamond_blastp_arg.txt \
    --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore \
    --evalue 1e-7 \
    --id 80 \
    --query-cover 70 \
    --max-target-seqs 1 \
    --sensitive \
    --max-hsps 1 \
    --threads 32
```

### Step 2: Extract Annotated Contigs

```bash
#!/bin/bash
# Extract contigs with ARG annotations

RESULT_FILE="${OUTPUT}/01_ARG_out/sample1_diamond_blastp_arg.txt"
CONTIG_FILE="results/contigs/sample1_contig.fa"

# Extract contig IDs from ORF IDs
tail -n +2 ${RESULT_FILE} | \
    cut -f1 | \
    sed -E 's/_[^_]+$//' | \
    sort -u > arg_contig_ids.txt

# Extract sequences
seqtk subseq ${CONTIG_FILE} arg_contig_ids.txt > ${OUTPUT}/02_ARG_contigs/sample1_arg_contigs.fa
```

---

## Complete Workflow Script

```bash
#!/bin/bash
# 03_Contigs_batch.sh - Batch processing for multiple samples

PROTEIN_DIR="results/prodigal"
CONTIG_DIR="results/contigs"
OUTPUT_DIR="results/annotation"
SARG_DB="/path/to/SARG.dmnd"
MGE_DB="/path/to/MobileOG.dmnd"
VF_DB="/path/to/VFDB.dmnd"

# Process samples
for sample in sample1 sample2 sample3; do
    echo "Processing: ${sample}"

    bash 03_Contigs_annotation.sh \
        -i ${PROTEIN_DIR} \
        -c ${CONTIG_DIR} \
        -o ${OUTPUT_DIR} \
        -a ${SARG_DB} \
        -m ${MGE_DB} \
        -v ${VF_DB} \
        -s ${sample}
done

# Merge results
python3 merge_blast_results.py \
    --source-dir ${OUTPUT_DIR}/01_ARG_out \
    --output-dir ${OUTPUT_DIR}/00_ARG_ALL \
    --suffix "_diamond_blastp_arg.txt"

python3 merge_blast_results.py \
    --source-dir ${OUTPUT_DIR}/04_MGE_out \
    --output-dir ${OUTPUT_DIR}/00_MGE_ALL \
    --suffix "_diamond_blastp_mge.txt"

python3 merge_blast_results.py \
    --source-dir ${OUTPUT_DIR}/07_VF_out \
    --output-dir ${OUTPUT_DIR}/00_VF_ALL \
    --suffix "_diamond_blastp_vf.txt"
```

---

## Usage

```bash
# Single sample
bash 03_Contigs_annotation.sh \
    -i results/prodigal \
    -c results/contigs \
    -o results/annotation \
    -a /path/to/SARG.dmnd \
    -m /path/to/MobileOG.dmnd \
    -v /path/to/VFDB.dmnd \
    -s sample1

# Batch processing
bash 03_Contigs_batch.sh
```

---

## Output Structure

```
results/annotation/
├── 01_ARG_out/
│   └── sample1_diamond_blastp_arg.txt    # ARG DIAMOND results
├── 02_ARG_contigs/
│   └── sample1_arg_contigs.fa            # Extracted ARG contigs
├── 03_ARG_total/
│   └── all_arg_contigs.fa                # Merged ARG contigs
├── 04_MGE_out/
│   └── sample1_diamond_blastp_mge.txt    # MGE DIAMOND results
├── 05_MGE_contigs/
│   └── sample1_mge_contigs.fa            # Extracted MGE contigs
├── 06_MGE_total/
│   └── all_mge_contigs.fa                # Merged MGE contigs
├── 07_VF_out/
│   └── sample1_diamond_blastp_vf.txt     # VF DIAMOND results
├── 08_VF_contigs/
│   └── sample1_vf_contigs.fa             # Extracted VF contigs
├── 09_VF_total/
│   └── all_vf_contigs.fa                 # Merged VF contigs
├── 00_ARG_ALL/
│   └── all_samples_merged_arg.csv        # Merged ARG results
├── 00_MGE_ALL/
│   └── all_samples_merged_mge.csv        # Merged MGE results
└── 00_VF_ALL/
    └── all_samples_merged_vf.csv         # Merged VF results
```

---

## Output Format (DIAMOND results)

| Column | Description |
|--------|-------------|
| qseqid | Query sequence ID (ORF ID) |
| sseqid | Subject sequence ID (database entry) |
| pident | Percentage identity |
| length | Alignment length |
| mismatch | Number of mismatches |
| gapopen | Number of gap openings |
| qstart | Query start position |
| qend | Query end position |
| sstart | Subject start position |
| send | Subject end position |
| evalue | E-value |
| bitscore | Bit score |

---

## Annotation Parameters

| Parameter | Value | Description |
|-----------|-------|-------------|
| Identity | >= 80% | Minimum sequence identity |
| Query coverage | >= 70% | Minimum query coverage |
| E-value | 1e-7 | Statistical significance threshold |
| Max targets | 1 | Best hit only |
| Sensitivity | sensitive | Sensitive alignment mode |

---

## Databases

### SARG (Structured ARG Database)
- Purpose: Antibiotic resistance gene annotation
- Version: v3.2 (Short subdatabase recommended for short reads)
- Categories: Tetracycline, Beta-lactam, Aminoglycoside, Sulfonamide, Macrolide, etc.

### MobileOG Database
- Purpose: Mobile genetic element annotation
- Version: beatrix-1-6
- Categories: Integrases, Transposases, Plasmids, Phages, etc.

### VFDB (Virulence Factor Database)
- Purpose: Virulence factor annotation
- Version: SetA (curated virulence factors)
- Categories: Adhesins, Toxins, Secretion systems, etc.

---

## Expected Results

| Database | Typical Hits | Description |
|----------|-------------|-------------|
| SARG | 50-500 per sample | ARG-carrying contigs |
| MobileOG | 100-1000 per sample | MGE-carrying contigs |
| VFDB | 20-200 per sample | VF-carrying contigs |

---

## Downstream Analysis

### ARG-MGE co-occurrence

```R
library(ggplot2)
library(ggvenn)

# Load merged results
arg_df <- read.csv("00_ARG_ALL/all_samples_merged_arg.csv")
mge_df <- read.csv("00_MGE_ALL/all_samples_merged_mge.csv")

# Extract contig IDs
arg_contigs <- unique(gsub("_[^_]+$", "", arg_df$qseqid))
mge_contigs <- unique(gsub("_[^_]+$", "", mge_df$qseqid))

# Find overlap (potential co-selection)
shared <- intersect(arg_contigs, mge_contigs)
cat("ARG-MGE co-occurring contigs:", length(shared))
```

---

## Citation

If you use this pipeline, please cite:

1. Buchfink B et al. Fast and sensitive protein alignment using DIAMOND. *Nat Methods*. 2015.
2. Arango-Argoty G et al. The structured antibiotic resistance genes (SARG) database. *Nucleic Acids Res*. 2018.
3. Jia B et al. CARD 2020: antibiotic resistome surveillance with the comprehensive antibiotic resistance database. *Nucleic Acids Res*. 2020.