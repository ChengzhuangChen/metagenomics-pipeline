# MAG Quality Assessment and Dereplication Pipeline

## Method Summary

MAG quality was assessed using CheckM and dereplicated using dRep to obtain non-redundant, high-quality genomes.

```
MAGs (from DAS_Tool) -> CheckM Quality Check -> Quality Filter -> dRep Dereplication -> Non-redundant MAGs
```

> MAG quality was assessed using CheckM (v1.2.0), and MAGs were dereplicated using dRep (v3.4.5) with a 95% ANI threshold.

---

## Dependencies

### Software Requirements

| Tool | Version | Purpose |
|------|---------|---------|
| CheckM | 1.2.0 | MAG quality assessment (completeness, contamination) |
| dRep | 3.4.5 | Genome dereplication |
| FastANI | 1.1+ | Average Nucleotide Identity calculation |

### Installation via Conda

```bash
conda create -n mag_quality -c bioconda -c conda-forge \
    checkm-genome=1.2.0 \
    drep=3.4.5 \
    fastani=1.33 \
    pandas=1.3.5

conda activate mag_quality
```

---

## Pipeline Steps

### Step 1: CheckM Quality Assessment

Assess completeness and contamination for all MAGs using marker gene analysis.

```bash
#!/bin/bash
# quality_assessment.sh

BINS_DIR="results/binning/dastool/DAS_Tool_DASToolbins"
OUTPUT="results/mag_quality/checkm"

mkdir -p ${OUTPUT}

# Run CheckM lineage_wf
# This places bins in the CheckM reference tree and calculates quality metrics
checkm lineage_wf \
    -t 32 \
    -x fa \
    ${BINS_DIR} \
    ${OUTPUT} \
    -f ${OUTPUT}/checkm_results.tsv

# Generate QA summary table
checkm qa \
    ${OUTPUT}/checkm_refine/msa_storage/tree/consensus_bins/consensus_tree.fa \
    ${OUTPUT} \
    -o 2 \
    -f ${OUTPUT}/checkm_qa_summary.tsv \
    --tab_table
```

**CheckM Output Columns:**
| Column | Description |
|--------|-------------|
| Bin Id | MAG identifier |
| Marker lineage | Taxonomy lineage |
| # genomes | Number of reference genomes |
| # markers | Number of marker genes |
| # marker sets | Number of marker sets |
| Completeness | Estimated genome completeness (%) |
| Contamination | Estimated contamination (%) |
| Strain heterogeneity | Strain heterogeneity estimate |

**Quality tiers (MIMAG standard):**
| Quality | Completeness | Contamination |
|---------|-------------|---------------|
| High-quality | >= 90% | < 5% |
| Medium-quality | >= 50% | < 10% |
| Low-quality | < 50% | - |

---

### Step 2: Quality Filtering

Filter MAGs based on quality thresholds.

```bash
#!/bin/bash
# filter_quality.sh

CHECKM_OUTPUT="results/mag_quality/checkm/checkm_results.tsv"
BINS_DIR="results/binning/dastool/DAS_Tool_DASToolbins"
OUTPUT="results/mag_quality/filtered"

mkdir -p ${OUTPUT}

# Filter: Completeness >= 50% AND Contamination < 10%
# CheckM columns: Bin Id (1), Marker lineage (2), # genomes (3), # markers (4),
#                 # marker sets (5), Completeness (6), Contamination (7)
awk -F'\t' -v comp=50 -v cont=10 \
    'NR==1 || ($6>=comp && $7<cont)' ${CHECKM_OUTPUT} > \
    ${OUTPUT}/quality_filtered_mags.tsv

# Copy filtered MAGs
while IFS=$'\t' read -r bin_id lineage genomes markers marker_sets completeness contamination heterogeneity; do
    [ "$bin_id" = "Bin Id" ] && continue
    find ${BINS_DIR} -name "${bin_id}.fa" -exec cp {} ${OUTPUT}/ \;
done < <(tail -n +2 ${OUTPUT}/quality_filtered_mags.tsv)

echo "Filtered MAGs: $(ls ${OUTPUT}/*.fa | wc -l)"
```

---

### Step 3: dRep Dereplication

Remove redundant MAGs based on ANI and alignment fraction.

```bash
#!/bin/bash
# dereplicate.sh

FILTERED_BINS="results/mag_quality/filtered/*.fa"
OUTPUT="results/mag_quality/drep"

mkdir -p ${OUTPUT}

# dRep dereplicate
# -sa: ANI threshold (default 95%)
# -nc: ANI for secondary clustering
# -cm: use Mash for coarse clustering
# -comp: minimum completeness
# -con: maximum contamination
dRep dereplicate \
    ${OUTPUT} \
    -p 32 \
    -g ${FILTERED_BINS} \
    -sa 95 \
    -nc 95 \
    -comp 50 \
    -con 10
```

**Parameters:**
| Parameter | Default | Description |
|-----------|---------|-------------|
| `-sa` | 95 | ANI threshold for primary clustering (%) |
| `-nc` | 95 | ANI threshold for secondary clustering (%) |
| `-cm` | - | Use Mash distances for coarse clustering |
| `-comp` | 50 | Minimum completeness threshold (%) |
| `-con` | 10 | Maximum contamination threshold (%) |

---

### Complete Workflow Script

```bash
#!/bin/bash
# 03_checkm_drep.sh

set -e

THREADS=32
BINS_DIR="results/binning/dastool/DAS_Tool_DASToolbins"
OUTPUT="results/mag_quality"

# Quality thresholds
COMPLETENESS=50
CONTAMINATION=10
ANI_THRESHOLD=95

echo "=========================================="
echo "MAG Quality and Dereplication Pipeline"
echo "=========================================="

# Step 1: CheckM
echo "[1/4] CheckM quality assessment..."
mkdir -p ${OUTPUT}/checkm
checkm lineage_wf -t ${THREADS} -x fa ${BINS_DIR} ${OUTPUT}/checkm \
    -f ${OUTPUT}/checkm/checkm_results.tsv

# Step 2: Filter
echo "[2/4] Filtering by quality..."
mkdir -p ${OUTPUT}/filtered
awk -F'\t' -v comp=${COMPLETENESS} -v cont=${CONTAMINATION} \
    'NR==1 || ($6>=comp && $7<cont)' \
    ${OUTPUT}/checkm/checkm_results.tsv > \
    ${OUTPUT}/filtered/quality_filtered_mags.tsv

while IFS=$'\t' read -r bin_id lineage genomes markers marker_sets completeness contamination heterogeneity; do
    [ "$bin_id" = "Bin Id" ] && continue
    find ${BINS_DIR} -name "${bin_id}.fa" -exec cp {} ${OUTPUT}/filtered/ \;
done < <(tail -n +2 ${OUTPUT}/filtered/quality_filtered_mags.tsv)

# Step 3: dRep
echo "[3/4] Dereplication with dRep..."
mkdir -p ${OUTPUT}/drep
dRep dereplicate ${OUTPUT}/drep \
    -p ${THREADS} \
    -g ${OUTPUT}/filtered/*.fa \
    -sa ${ANI_THRESHOLD} \
    -nc ${ANI_THRESHOLD} \
    -comp ${COMPLETENESS} \
    -con ${CONTAMINATION}

echo "=========================================="
echo "Pipeline complete!"
echo "=========================================="
echo "CheckM results: ${OUTPUT}/checkm/"
echo "Filtered MAGs: ${OUTPUT}/filtered/"
echo "Non-redundant MAGs: ${OUTPUT}/drep/dereplicated_genomes/"
```

---

## Usage

```bash
# Standard usage
bash 03_checkm_drep.sh \
    -i results/binning/dastool/DAS_Tool_DASToolbins \
    -o results/mag_quality

# With custom parameters
bash 03_checkm_drep.sh \
    -i bins/ \
    -o output \
    -t 64 \
    -c 90 \
    -x 5 \
    -a 98
```

---

## Output Structure

```
results/mag_quality/
├── checkm/
│   ├── checkm_results.tsv         # CheckM quality results
│   ├── checkm_refine/             # CheckM refine output
│   └── checkm_qa_summary.tsv      # QA summary table
├── filtered/
│   ├── quality_filtered_mags.tsv
│   └── MAG.*.fa                   # Quality-filtered MAGs
├── drep/
│   └── dereplicated_genomes/      # Final non-redundant MAGs
│       └── *.fa
└── summary_report.txt             # Summary statistics
```

---

## Expected Results

| Metric | Expected Range |
|--------|----------------|
| Input MAGs | 500-5000 |
| After quality filter | 200-2000 |
| After dRep dereplication | 100-1000 |

---

## Citation

If you use this pipeline, please cite:

1. Parks DH et al. CheckM: assessing the quality of microbial genomes recovered from isolates, single cells, and metagenomes. *Genome Res*. 2015.
2. Olm MR et al. dRep: a tool for fast and accurate genomic comparisons that enables improved genome recovery from metagenomes. *ISME J*. 2017.