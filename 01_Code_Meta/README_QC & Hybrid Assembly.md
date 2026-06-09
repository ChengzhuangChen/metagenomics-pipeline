# Metagenomic Data Processing Pipeline

## Overview

Complete pipeline for metagenomic data analysis including quality control, hybrid assembly, and genome binning.

**Pipeline workflow:**
```
Raw FASTQ → QC & Assembly → Binning → Quality-filtered MAGs
```

---

## Part I: QC & Hybrid Assembly (`01_QC_Co_assembly.sh`)

### Method Summary

Two-step hybrid assembly combining metaSPAdes (for primary assembly) and MEGAHIT (for unmapped reads) to maximize genome recovery.

```
Raw FASTQ → Trimmomatic → metaSPAdes → Bowtie2 mapping → MEGAHIT → Merge & Filter → QUAST
```

### Quick Start

```bash
# Create environment
conda env create -f environment.yml
conda activate metagenomics

# Single sample assembly
bash 01_QC_Co_assembly.sh -s sample1 raw/sample1_R1.fastq.gz raw/sample1_R2.fastq.gz

# Multiple samples
bash 01_QC_Co_assembly.sh -m samples.txt

# Co-assembly (combine all samples)
bash 01_QC_Co_assembly.sh -c samples.txt
```

### Output Structure

```
results/
├── sample1/
│   ├── clean/              # Trimmed reads
│   ├── qc/                 # FastQC/MultiQC reports
│   ├── metaspades/         # Primary assembly
│   ├── unmapped/           # Unmapped reads
│   ├── megahit/            # Secondary assembly
│   ├── final/
│   │   └── contigs_1000bp.fa   # Final contigs (≥1kb)
│   └── quast/              # Assembly quality
├── sample2/...
└── coassembly/
    └── final/contigs_1000bp.fa
```

---

## Part II: Metagenome Binning (`02_Binning_four_tools.sh`)

### Method Summary

Metagenome-assembled genomes (MAGs) recovered using four independent binning tools and integrated via DAS_Tool.

> Metagenome binning was conducted using four independent tools: MetaBAT2 (v2.12.1), MaxBin2 (v2.2.6), CONCOCT (v1.1.0), and SemiBin2 (v2.2). Only contigs with lengths ≥ 2000 bp were subjected to binning, generating preliminary MAGs. These MAGs were integrated and dereplicated via DAS_Tool to obtain non-redundant MAGs.

### Pipeline Workflow

```
Contigs (≥2000bp)
    │
    ├─→ MetaBAT2  ─┐
    ├─→ MaxBin2  ──┼──→ DAS_Tool ─→ Non-redundant MAGs
    ├─→ CONCOCT  ──┤
    └─→ SemiBin2 ──┘
```

### Quick Start

```bash
# Run complete binning pipeline
bash 02_Binning_four_tools.sh \
    -a results/assembly/contigs_1000bp.fa \
    -1 results/clean/sample_R1_paired.fastq.gz \
    -2 results/clean/sample_R2_paired.fastq.gz \
    -o results/binning

# With custom parameters
bash 02_Binning_four_tools.sh \
    -a assembly/contigs.fasta \
    -1 clean/R1.fq.gz \
    -2 clean/R2.fq.gz \
    -o binning \
    -t 64 \
    -m 2500
```

### Output Structure

```
results/binning/
├── contigs_2000bp.fa           # Filtered contigs
├── aligned_sorted.bam          # Read alignments
├── metabat2/                   # MetaBAT2 bins
│   ├── depth.txt
│   └── bin.1.fa, bin.2.fa...
├── maxbin2/                    # MaxBin2 bins
│   └── bin.001.fasta...
├── concoct/                    # CONCOCT bins
│   └── bins/bin_1.fa...
├── semibin2/                   # SemiBin2 bins
│   └── bin_0.fa...
└── dastool/                    # Dereplicated MAGs
    ├── DAS_Tool_DASToolbins/   # <-- Final MAGs
    │   └── MAG.1.fa, MAG.2.fa...
    └── DAS_Tool_summary.txt
```

### Expected Results

| Metric | Expected Range |
|--------|----------------|
| MetaBAT2 bins | 200-2000 |
| MaxBin2 bins | 100-800 |
| CONCOCT bins | 150-1500 |
| SemiBin2 bins | 200-2500 |
| **Total preliminary MAGs** | **500-5000** |
| **DAS_Tool non-redundant MAGs** | **100-1000** |

---

## Installation

### Prerequisites

- Linux/macOS with bash
- Conda or Mamba
- 32+ CPU cores recommended
- 256+ GB RAM recommended for large datasets

### Setup

```bash
# Clone repository
git clone https://github.com/ChengzhuangChen/metagenomics-pipeline.git
cd metagenomics-pipeline

# Create conda environment
conda env create -f environment.yml
conda activate metagenomics

# Download Trimmomatic adapter files (if needed)
# Adapter files are included in Trimmomatic distribution
```

---

## File Descriptions

| File | Description |
|------|-------------|
| `01_QC_Co_assembly.sh` | QC, trimming, and hybrid assembly pipeline |
| `02_Binning_four_tools.sh` | Four binning tools + DAS_Tool integration |
| `environment.yml` | Conda environment specification |
| `samples.txt` | Sample list template for batch processing |
| `QUICKSTART.md` | Quick start guide for assembly |
| `README_binning.md` | Detailed binning documentation |

---

## Citation

If you use this pipeline in your research, please cite:

### Pipeline
```
Metagenomic Data Processing Pipeline.
[GitHub URL]
```

### Underlying Tools

1. Bolger AM et al. Trimmomatic. *Bioinformatics*. 2014.
2. Nurk S et al. metaSPAdes. *Genome Res*. 2017.
3. Li D et al. MEGAHIT. *BMC Bioinformatics*. 2015.
4. Kang DD et al. MetaBAT2. *PeerJ*. 2019.
5. Wu YW et al. MaxBin2. *Nucleic Acids Res*. 2016.
6. Alneberg J et al. CONCOCT. *Nat Methods*. 2014.
7. Pan S et al. SemiBin2. *Bioinformatics*. 2022.
8. Sieber CMK et al. DAS_Tool. *Nat Microbiol*. 2018.

---

## Data Availability Statement (Template)

```
Raw sequencing data are available at NCBI SRA (PRJNAxxxxxx).
Assembled contigs are available at Zenodo (doi:10.xxxx/zenodo.xxxxxx).
Non-redundant MAGs are available at Zenodo (doi:10.xxxx/zenodo.xxxxxx).
Code is available at: https://github.com/ChengzhuangChen/metagenomics-pipeline
```

---

## License

MIT License