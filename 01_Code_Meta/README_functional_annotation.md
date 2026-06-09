# MAG Functional Annotation Pipeline

## Method Summary

Protein-coding genes were predicted using Prodigal, clustered with MMseqs2 to remove redundancy, and annotated using eggNOG-mapper.

> Protein-coding genes were predicted using Prodigal (v2.6.3), and predicted proteins were clustered using MMseqs2 (v13.10711) with 95% identity and 80% coverage thresholds to remove redundancy. Functional annotation was performed using eggNOG-mapper (v2.1.12) with DIAMOND against the eggNOG database.

```
MAGs -> Prodigal -> Merged FAA -> MMseqs2 Clustering -> eggNOG Annotation
                                                              |
                                          COG, KEGG, GO, Pfam, CAZy annotations
```

---

## Dependencies

### Software Requirements

| Tool | Version | Purpose |
|------|---------|---------|
| Prodigal | 2.6.3 | Gene prediction |
| MMseqs2 | 13.10711 | Sequence clustering |
| eggNOG-mapper | 2.1.12 | Functional annotation |
| DIAMOND | 2.0.15+ | Sequence alignment |

### Installation via Conda

```bash
conda create -n annotation -c conda-forge -c bioconda \
    prodigal=2.6.3 \
    mmseqs2=13.10711 \
    eggnog-mapper=2.1.12 \
    diamond=2.0.15

conda activate annotation

# Download eggNOG database (first time, ~40 GB)
download_eggnog_data.py
```

---

## Pipeline Steps

### Step 1: Gene Prediction with Prodigal

Predict protein-coding genes for each MAG.

```bash
#!/bin/bash
# prodigal_prediction.sh

MAG_DIR="results/mag_quality/drep/dereplicated_genomes"
OUTPUT="results/functional_annotation/01_prodigal_faa"

mkdir -p ${OUTPUT}

for mag in ${MAG_DIR}/*.fa; do
    mag_name=$(basename ${mag} .fa)

    prodigal \
        -i ${mag} \
        -o ${OUTPUT}/${mag_name}.gff \
        -a ${OUTPUT}/${mag_name}.faa \
        -p meta \
        -g 11 \
        -q
done

echo "Total MAGs: $(ls ${OUTPUT}/*.faa | wc -l)"
echo "Total genes: $(cat ${OUTPUT}/*.faa | grep '^>' | wc -l)"
```

**Parameters:**
| Parameter | Description |
|-----------|-------------|
| `-p meta` | Metagenomic mode (shorter genes acceptable) |
| `-g 11` | Translation table 11 (Bacterial) |
| `-a` | Output protein translations |
| `-o` | Output GFF3 format |

---

### Step 2: Merge FAA Files

Combine all protein sequences.

```bash
#!/bin/bash
# merge_faa.sh

FAA_DIR="results/functional_annotation/01_prodigal_faa"
OUTPUT_DIR="results/functional_annotation/02_merge_faa"

mkdir -p ${OUTPUT_DIR}

cat ${FAA_DIR}/*.faa > ${OUTPUT_DIR}/all_genes.faa

echo "Total sequences: $(grep -c '^>' ${OUTPUT_DIR}/all_genes.faa)"
```

---

### Step 3: MMseqs2 Clustering

Remove redundant proteins.

```bash
#!/bin/bash
# mmseqs2_clustering.sh

INPUT_FAA="results/functional_annotation/02_merge_faa/all_genes.faa"
OUTPUT_DIR="results/functional_annotation/03_mmseqs2_result"
THREADS=32

mkdir -p ${OUTPUT_DIR}/mmseqs_tmp

# Run MMseqs2 easy-cluster
# --min-seq-id: 95% identity
# -c: 80% coverage
# --cov-mode 1: Query and target coverage
# --cluster-mode 2: Greedy set cover
mmseqs easy-cluster \
    ${INPUT_FAA} \
    ${OUTPUT_DIR}/clustered \
    ${OUTPUT_DIR}/mmseqs_tmp \
    --min-seq-id 0.95 \
    -c 0.80 \
    --cov-mode 1 \
    --cluster-mode 2 \
    --threads ${THREADS}

echo "Input sequences: $(grep -c '^>' ${INPUT_FAA})"
echo "Non-redundant clusters: $(grep -c '^>' ${OUTPUT_DIR}/clustered_rep_seq.faa)"
```

**Clustering Parameters:**
| Parameter | Value | Description |
|-----------|-------|-------------|
| `--min-seq-id` | 0.95 | 95% sequence identity |
| `-c` | 0.80 | 80% coverage |
| `--cov-mode` | 1 | Qcov + Tcov (both query and target) |
| `--cluster-mode` | 2 | Greedy set cover (fastest) |

---

### Step 4: eggNOG Annotation

Functional annotation of non-redundant proteins.

```bash
#!/bin/bash
# eggnog_annotation.sh

INPUT_FAA="results/functional_annotation/03_mmseqs2_result/clustered_rep_seq.faa"
OUTPUT_DIR="results/functional_annotation/04_eggnog"
THREADS=32

mkdir -p ${OUTPUT_DIR}

emapper.py \
    -i ${INPUT_FAA} \
    -o ${OUTPUT_DIR}/annotation \
    -m diamond \
    --cpu ${THREADS} \
    --evalue 1e-5 \
    --sensmode sensitive \
    --report_orthologs
```

**Parameters:**
| Parameter | Description |
|-----------|-------------|
| `-m diamond` | Use DIAMOND for alignment |
| `--evalue` | E-value threshold |
| `--sensmode` | Sensitivity mode (sensitive, ultra-sensitive) |
| `--report_orthologs` | Report ortholog information |

---

## Complete Workflow Script

```bash
#!/bin/bash
# 06_MAG_eggNOG.sh

set -e

MAGS_DIR="results/mag_quality/drep/dereplicated_genomes"
OUTPUT_DIR="results/functional_annotation"
THREADS=32

echo "=========================================="
echo "MAG Functional Annotation Pipeline"
echo "=========================================="

# Step 1: Prodigal
echo "[1/4] Gene prediction with Prodigal..."
mkdir -p ${OUTPUT_DIR}/01_prodigal_faa

TOTAL_GENES=0
for mag in ${MAGS_DIR}/*.fa; do
    mag_name=$(basename ${mag} .fa)
    prodigal -i ${mag} -a ${OUTPUT_DIR}/01_prodigal_faa/${mag_name}.faa -p meta -q
    TOTAL_GENES=$((TOTAL_GENES + $(grep -c '^>' ${OUTPUT_DIR}/01_prodigal_faa/${mag_name}.faa)))
done
echo "  -> Total genes: ${TOTAL_GENES}"

# Step 2: Merge
echo "[2/4] Merging FAA files..."
mkdir -p ${OUTPUT_DIR}/02_merge_faa
cat ${OUTPUT_DIR}/01_prodigal_faa/*.faa > ${OUTPUT_DIR}/02_merge_faa/all_genes.faa

# Step 3: MMseqs2
echo "[3/4] Clustering with MMseqs2..."
mkdir -p ${OUTPUT_DIR}/03_mmseqs2_result/mmseqs_tmp
mmseqs easy-cluster \
    ${OUTPUT_DIR}/02_merge_faa/all_genes.faa \
    ${OUTPUT_DIR}/03_mmseqs2_result/clustered \
    ${OUTPUT_DIR}/03_mmseqs2_result/mmseqs_tmp \
    --min-seq-id 0.95 -c 0.80 --cov-mode 1 --threads ${THREADS}

CLUSTER_COUNT=$(grep -c '^>' ${OUTPUT_DIR}/03_mmseqs2_result/clustered_rep_seq.faa)
echo "  -> Non-redundant clusters: ${CLUSTER_COUNT}"

# Step 4: eggNOG
echo "[4/4] Annotating with eggNOG..."
mkdir -p ${OUTPUT_DIR}/04_eggnog
emapper.py -i ${OUTPUT_DIR}/03_mmseqs2_result/clustered_rep_seq.faa \
    -o ${OUTPUT_DIR}/04_eggnog/annotation \
    -m diamond --cpu ${THREADS} --evalue 1e-5 --sensmode sensitive

echo "=========================================="
echo "Pipeline complete!"
echo "=========================================="
```

---

## Usage

```bash
# Standard usage
bash 06_MAG_eggNOG.sh \
    -i results/mag_quality/drep/dereplicated_genomes \
    -o results/functional_annotation

# With custom parameters
bash 06_MAG_eggNOG.sh -i bins/ -o output -t 64 --min-id 0.98 --evalue 1e-10
```

---

## Output Structure

```
results/functional_annotation/
├── 01_prodigal_faa/
│   ├── MAG.001.faa
│   ├── MAG.002.faa
│   └── ...
├── 02_merge_faa/
│   └── all_genes.faa              # All merged proteins
├── 03_mmseqs2_result/
│   ├── clustered_rep_seq.faa      # Non-redundant proteins
│   ├── clustered_all_seqs.fasta   # All cluster members
│   └── clustered_cluster.tsv      # Cluster assignments
├── 04_eggnog/
│   ├── annotation.emapper.annotations   # Full annotations
│   └── annotation.emapper.seed_orthologs
└── summary_report.txt
```

---

## eggNOG Output Format

### Annotations File (TSV)

```tsv
query	seed_ortholog	evalue	score	eggNOG OGs	COG functional categories	KEGG KOs	PFAMs
gene_001	COG0001@1@root	1.2e-45	245.5	COG0001@2@root	[H]	ko:K00001	PF00001
gene_002	COG0002@1@root	5.3e-32	198.3	COG0002@2@root	[E]	ko:K00002	PF00002
```

**Columns:**
| Column | Description |
|--------|-------------|
| query | Protein query ID |
| seed_ortholog | Best matching ortholog |
| evalue | E-value of alignment |
| score | Alignment score |
| eggNOG OGs | Orthologous group IDs |
| COG functional categories | COG category codes |
| KEGG KOs | KEGG Orthology IDs |
| PFAMs | Pfam domain IDs |

---

## Expected Results

| Step | Expected Value |
|------|----------------|
| Input MAGs | 500-1000 |
| Predicted genes | 50,000-500,000 |
| Non-redundant clusters (95% ID) | 20,000-200,000 |
| Annotated unigenes | 15,000-150,000 (75-80%) |

---

## Downstream Analysis

### Extract specific annotations

```bash
# Get KO annotations
cut -f8 results/functional_annotation/04_eggnog/annotation.emapper.annotations | \
    tail -n +2 | grep -v '-' | sort | uniq -c | sort -rn > ko_counts.txt

# Get COG category distribution
cut -f7 results/functional_annotation/04_eggnog/annotation.emapper.annotations | \
    tail -n +2 | grep -v '-' | grep -o '.' | sort | uniq -c | sort -rn > cog_categories.txt
```

### Pathway enrichment analysis

```R
library(clusterProfiler)

# Load KO annotations
kodata <- read.delim("annotation.emapper.annotations", header=TRUE)
kos <- gsub("ko:", "", strsplit(kodata$KEGG.KOs, ","))

# KEGG pathway enrichment
kk <- enrichKEGG(gene = kos, organism = 'ko', pvalueCutoff = 0.05)
barplot(kk, showCategory = 20)
```

---

## Troubleshooting

### Prodigal fails on some MAGs
```bash
# Try meta mode with single genome mode
prodigal -i mag.fa -a mag.faa -p single -q
```

### MMseqs2 out of memory
```bash
# Use more threads to reduce memory per thread
# Or use greedy set cover mode
mmseqs easy-cluster ... --cluster-mode 2
```

### eggNOG database missing
```bash
# Download database
download_eggnog_data.py

# Check database location
echo $EGGNOG_DATA_DIR
```

---

## Citation

If you use this pipeline, please cite:

1. Hyatt D et al. Prodigal: prokaryotic gene recognition and translation initiation site identification. *BMC Bioinformatics*. 2010.
2. Steinegger M, Soding J. MMseqs2 enables sensitive protein sequence searching. *Nat Methods*. 2017.
3. Cantalapiedra CP et al. eggNOG-mapper v2: functional annotation, orthology assignments, and domain prediction at the metagenomic scale. *Mol Biol Evol*. 2021.