# GTDB-Tk Taxonomy Classification Pipeline

## Method Summary

MAG taxonomy was determined using GTDB-Tk with the classify_wf workflow, assigning MAGs to taxonomic levels from domain to species based on the Genome Taxonomy Database.

> MAG taxonomy was determined using GTDB-Tk (v2.4.0) with the GTDB r223 database. Taxonomic assignment included domain, phylum, class, order, family, genus, and species levels based on 120 bacterial and 53 archaeal marker genes.

---

## Dependencies

### Software Requirements

| Tool | Version | Purpose |
|------|---------|---------|
| GTDB-Tk | 2.4.0+ | Taxonomic classification |
| FastTree | 2.1+ | Phylogenetic tree inference |
| pplacer | 1.1+ | Phylogenetic placement |

### Installation via Conda

```bash
conda create -n gtdbtk -c bioconda -c conda-forge \
    gtdbtk=2.4.0 \
    fasttree=2.1.10 \
    pplacer=1.1

conda activate gtdbtk

# Download GTDB reference data (~50 GB)
gtdbtk download --db-version r223
```

**Database sizes:**
| Version | Approximate Size |
|---------|-----------------|
| r202 | ~40 GB |
| r214 | ~50 GB |
| r223 | ~60 GB |

---

## Pipeline Steps

### Step 1: GTDB-Tk Classify

Taxonomic classification using the classify_wf workflow.

```bash
#!/bin/bash
# gtdbtk_classify.sh

BINS_DIR="results/mag_quality/drep/dereplicated_genomes"
OUTPUT="results/taxonomy"
THREADS=32

mkdir -p ${OUTPUT}

gtdbtk classify_wf \
    --genome_dir ${BINS_DIR} \
    --out_dir ${OUTPUT} \
    --extension fa \
    --prefix bin \
    --cpus ${THREADS}
```

**Parameters:**
| Parameter | Description |
|-----------|-------------|
| `--genome_dir` | Directory containing MAG FASTA files |
| `--out_dir` | Output directory |
| `--extension` | File extension (fa, fasta) |
| `--prefix` | Output file prefix |
| `--cpus` | Number of threads |
| `--skip_ani_screen` | Skip ANI screening (faster) |
| `--full_tree` | Use full reference tree |

**Workflow steps:**
1. Identify marker genes (120 bacterial, 53 archaeal)
2. Align marker genes with reference proteins
3. Identify and filter chimeric alignments
4. Place MAGs in reference tree (pplacer)
5. Derive taxonomy from tree placement

---

### Step 2: Decompress Alignments

Decompress the MSA files for tree inference.

```bash
# Decompress alignment files
gunzip -k ${OUTPUT}/align/bin.ar53.user_msa.fasta.gz
gunzip -k ${OUTPUT}/align/bin.bac120.user_msa.fasta.gz
```

---

### Step 3: Phylogenetic Tree Inference

Infer phylogenetic trees for bacterial and archaeal MAGs.

```bash
#!/bin/bash
# gtdbtk_infer.sh

OUTPUT="results/taxonomy"
THREADS=32

# Infer archaeal tree
if [ -f "${OUTPUT}/align/bin.ar53.user_msa.fasta" ]; then
    mkdir -p ${OUTPUT}/tree_ar53
    gtdbtk infer \
        --msa_file ${OUTPUT}/align/bin.ar53.user_msa.fasta \
        --out_dir ${OUTPUT}/tree_ar53 \
        --cpus ${THREADS} \
        --prefix bin
fi

# Infer bacterial tree
if [ -f "${OUTPUT}/align/bin.bac120.user_msa.fasta" ]; then
    mkdir -p ${OUTPUT}/tree_bac120
    gtdbtk infer \
        --msa_file ${OUTPUT}/align/bin.bac120.user_msa.fasta \
        --out_dir ${OUTPUT}/tree_bac120 \
        --cpus ${THREADS} \
        --prefix bin
fi
```

---

## Complete Workflow Script

```bash
#!/bin/bash
# 04_GTDB-tk.sh

set -e

BINS_DIR="results/mag_quality/drep/dereplicated_genomes"
OUTPUT="results/taxonomy"
THREADS=32
EXTENSION="fa"
PREFIX="bin"

echo "=========================================="
echo "GTDB-Tk Taxonomy Classification"
echo "=========================================="

# Step 1: Classify
echo "[1/2] Running GTDB-Tk classify_wf..."
mkdir -p ${OUTPUT}

gtdbtk classify_wf \
    --genome_dir ${BINS_DIR} \
    --out_dir ${OUTPUT} \
    --extension ${EXTENSION} \
    --prefix ${PREFIX} \
    --cpus ${THREADS}

# Step 2: Decompress and infer trees
echo "[2/2] Inferring phylogenetic trees..."

gunzip -k ${OUTPUT}/align/bin.ar53.user_msa.fasta.gz 2>/dev/null || true
gunzip -k ${OUTPUT}/align/bin.bac120.user_msa.fasta.gz 2>/dev/null || true

if [ -f "${OUTPUT}/align/bin.ar53.user_msa.fasta" ]; then
    mkdir -p ${OUTPUT}/tree_ar53
    gtdbtk infer --msa_file ${OUTPUT}/align/bin.ar53.user_msa.fasta \
        --out_dir ${OUTPUT}/tree_ar53 --cpus ${THREADS} --prefix ${PREFIX}
fi

if [ -f "${OUTPUT}/align/bin.bac120.user_msa.fasta" ]; then
    mkdir -p ${OUTPUT}/tree_bac120
    gtdbtk infer --msa_file ${OUTPUT}/align/bin.bac120.user_msa.fasta \
        --out_dir ${OUTPUT}/tree_bac120 --cpus ${THREADS} --prefix ${PREFIX}
fi

echo "=========================================="
echo "Pipeline complete!"
echo "=========================================="
```

---

## Usage

```bash
# Standard usage
bash 04_GTDB-tk.sh \
    -i results/mag_quality/drep/dereplicated_genomes \
    -o results/taxonomy

# With custom parameters
bash 04_GTDB-tk.sh \
    -i bins/ \
    -o taxonomy \
    -t 64 \
    -e fa \
    -p mygenomes

# Skip ANI screening (faster, for large datasets)
bash 04_GTDB-tk.sh -i bins/ -o taxonomy -s
```

---

## Output Structure

```
results/taxonomy/
├── classify/
│   ├── bin.bac120.summary.tsv         # Bacterial classification
│   ├── bin.ar53.summary.tsv           # Archaeal classification
│   ├── bin.bac120.markers_summary.tsv # Bacterial markers
│   ├── bin.ar53.markers_summary.tsv   # Archaeal markers
│   ├── bin.bac120.metadata.tsv        # Bacterial metadata
│   └── bin.ar53.metadata.tsv          # Archaeal metadata
├── align/
│   ├── bin.bac120.user_msa.fasta      # Bacterial alignment
│   └── bin.ar53.user_msa.fasta        # Archaeal alignment
├── tree_bac120/
│   └── bin.bac120.user_msa.tree       # Bacterial tree (Newick)
└── tree_ar53/
    └── bin.ar53.user_msa.tree         # Archaeal tree (Newick)
```

---

## Output Format

### Classification Summary (TSV)

```tsv
user_genome	classification	fastani_reference	similarity	ANI	AF
bin.001	d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Enterobacterales;f__Enterobacteriaceae;g__Escherichia;s__Escherichia_coli	RS_GCF_000005845.2	-	-	-
bin.002	d__Bacteria;p__Bacteroidota;c__Bacteroidia;o__Bacteroidales;f__Bacteroidaceae;g__Bacteroides;s__Bacteroides_thetaiotaomicron	RS_GCF_000009125.1	-	-	-
```

**Columns:**
| Column | Description |
|--------|-------------|
| user_genome | MAG identifier |
| classification | GTDB taxonomy string |
| fastani_reference | Closest reference genome |
| similarity | Similarity to reference |
| ANI | Average Nucleotide Identity |
| AF | Alignment Fraction |

### Taxonomy String Format

```
d__Domain;p__Phylum;c__Class;o__Order;f__Family;g__Genus;s__Species
```

Example: `d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Enterobacterales;f__Enterobacteriaceae;g__Escherichia;s__Escherichia_coli`

---

## Expected Runtime

| Number of MAGs | Estimated Time |
|----------------|----------------|
| 100 MAGs | ~2-4 hours |
| 500 MAGs | ~8-12 hours |
| 1000 MAGs | ~24-48 hours |

---

## Troubleshooting

### Out of memory
```bash
# Increase memory allocation
# GTDB-Tk requires ~8 GB per 100 MAGs

# Use skip_ani_screen for faster processing
gtdbtk classify_wf ... --skip_ani_screen
```

### Missing GTDB database
```bash
# Download database
gtdbtk download --db-version r223

# Check database location
echo $GTDBTK_DATA_PATH
```

### Low classification rate
- Check MAG quality (CheckM completeness)
- Ensure proper file naming and format
- Verify genome files are complete FASTA

---

## Downstream Analysis

### Extract taxonomy summary
```bash
# Get taxonomy at phylum level
cut -f2 ${OUTPUT}/classify/bin.bac120.summary.tsv | \
    tail -n +2 | cut -d';' -f2 | sed 's/p__//' | \
    sort | uniq -c | sort -rn

# Get species-level classification
cut -f2 ${OUTPUT}/classify/bin.bac120.summary.tsv | \
    tail -n +2 | cut -d';' -f7 | sed 's/s__//' | \
    sort | uniq -c | sort -rn | head -20
```

### Visualize phylogenetic tree
```R
# Using ggtree in R
library(ggtree)
tree <- read.tree("tree_bac120/bin.bac120.user_msa.tree")
ggtree(tree) + geom_tiplab()
```

---

## Citation

If you use this pipeline, please cite:

1. Chaumeil PA et al. GTDB-Tk: a toolkit to classify genomes with the Genome Taxonomy Database. *Bioinformatics*. 2020.
2. Parks DH et al. A standardized bacterial taxonomy based on genome phylogeny substantially revises the tree of life. *Nat Biotechnol*. 2018.

---

## Data Availability

```
GTDB database versions available at: https://gtdb.ecogenomic.org/
```