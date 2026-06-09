# Metagenomic Data Processing Pipeline

## Overview

Complete pipeline for metagenomic data analysis including quality control, hybrid assembly, genome binning, MAG quality assessment, taxonomy classification, abundance profiling, and functional annotation.

**Pipeline workflow:**
```
Raw FASTQ → QC & Assembly → Binning → Quality-filtered MAGs → Taxonomy → Abundance → Functional Annotation
```

---

## Quick Start

### 1. Installation

```bash
# Clone repository
git clone https://github.com/ChengzhuangChen/metagenomics-pipeline.git
cd metagenomics-pipeline

# Create conda environment
conda env create -f environment.yml
conda activate metagenomics
```

### 2. Run Pipeline

```bash
# Step 1: QC and Assembly
bash 01_QC_Co_assembly.sh -s sample1 raw/sample1_R1.fastq.gz raw/sample1_R2.fastq.gz

# Step 2: Binning
bash 02_Binning_four_tools.sh -a results/assembly/contigs_1000bp.fa \
    -1 results/clean/sample_R1.fastq.gz -2 results/clean/sample_R2.fastq.gz \
    -o results/binning

# Step 3: Quality Assessment and Dereplication
bash 03_checkm_drep.sh -i results/binning/dastool/DAS_Tool_DASToolbins \
    -o results/mag_quality

# Step 4: Taxonomy Classification
bash 04_GTDB-tk.sh -i results/mag_quality/drep/dereplicated_genomes \
    -o results/taxonomy

# Step 5: Abundance Profiling
bash 05_coverm.sh -m results/mag_quality/drep/dereplicated_genomes \
    -f results/clean_data -o results/abundance

# Step 6: Functional Annotation
bash 06_MAG_eggNOG.sh -i results/mag_quality/drep/dereplicated_genomes \
    -o results/functional_annotation
```

---

## Pipeline Scripts

| Script | Description |
|--------|-------------|
| `01_QC_Co_assembly.sh` | Quality control, trimming, and hybrid assembly |
| `02_Binning_four_tools.sh` | MAG binning (MetaBAT2, MaxBin2, CONCOCT, SemiBin2) + DAS_Tool |
| `03_checkm_drep.sh` | MAG quality assessment (CheckM) and dereplication (dRep) |
| `04_GTDB-tk.sh` | Taxonomic classification using GTDB-Tk |
| `05_coverm.sh` | MAG abundance profiling using CoverM |
| `06_MAG_eggNOG.sh` | Functional annotation (Prodigal + MMseqs2 + eggNOG) |

---

## Documentation

Detailed documentation for each step is available in the README files:

- `README_QC & Hybrid Assembly.md` - Assembly pipeline details
- `README_binning.md` - Binning pipeline details
- `README_MAG_quality.md` - Quality assessment and dereplication
- `README_GTDB-tk.md` - Taxonomy classification
- `README_coverm.md` - Abundance profiling
- `README_functional_annotation.md` - Functional annotation

---

## File Descriptions

| File | Description |
|------|-------------|
| `environment.yml` | Conda environment specification |
| `samples.txt` | Sample list template for batch processing |

---

## Directory Structure

```
results/
├── sample1/
│   ├── clean/              # Trimmed reads
│   ├── metaspades/         # Primary assembly
│   ├── unmapped/           # Unmapped reads
│   ├── megahit/            # Secondary assembly
│   └── final/              # Final contigs
├── binning/
│   ├── metabat2/           # MetaBAT2 bins
│   ├── maxbin2/            # MaxBin2 bins
│   ├── concoct/            # CONCOCT bins
│   ├── semibin2/           # SemiBin2 bins
│   └── dastool/            # Dereplicated MAGs
├── mag_quality/
│   ├── checkm/              # CheckM results
│   └── drep/                # Dereplicated MAGs
├── taxonomy/
│   └── classify/           # GTDB-Tk results
├── abundance/
│   └── abundance_matrix.tsv
└── functional_annotation/
    ├── 01_prodigal_faa/
    ├── 03_mmseqs2_result/
    └── 04_eggnog/
```

---

## Dependencies

All dependencies are specified in `environment.yml`:

- Trimmomatic 0.39
- FastQC 0.12.1
- MultiQC 1.30
- SPAdes 4.2.0
- MEGAHIT 1.2.9
- MetaBAT2 2.12.1
- MaxBin2 2.2.6
- CONCOCT 1.1.0
- SemiBin2 2.2.0
- DAS_Tool 1.1.6
- CheckM 1.2.0
- dRep 3.4.5
- GTDB-Tk 2.4.0
- CoverM 0.6.1
- Prodigal 2.6.3
- MMseqs2 13.10711
- eggNOG-mapper 2.1.12

---

## Citation

If you use this pipeline in your research, please cite:

### Pipeline
```
Metagenomic Data Processing Pipeline.
https://github.com/ChengzhuangChen/metagenomics-pipeline
```

### Key References

1. Bolger AM et al. Trimmomatic: a flexible trimmer for Illumina sequence data. *Bioinformatics*. 2014;30(15):2114-2120.
2. Nurk S et al. metaSPAdes: a new versatile metagenomic assembler. *Genome Res*. 2017;27(5):824-834.
3. Li D et al. MEGAHIT: an ultra-fast single-node solution for large and complex metagenomics assembly. *BMC Bioinformatics*. 2015;16:1-12.
4. Kang DD et al. MetaBAT2: an adaptive binning algorithm for robust and efficient genome reconstruction. *PeerJ*. 2019.
5. Wu YW et al. MaxBin2: automated binning of metagenomic contigs. *Nucleic Acids Res*. 2016.
6. Alneberg J et al. Binning metagenomic contigs by coverage and composition. *Nat Methods*. 2014.
7. Chaumeil PA et al. GTDB-Tk: a toolkit to classify genomes with GTDB. *Bioinformatics*. 2020.
8. Olm MR et al. dRep: a tool for fast and accurate genomic comparisons. *ISME J*. 2017.
9. Parks DH et al. CheckM: assessing the quality of microbial genomes recovered from isolates, single cells, and metagenomes. *Genome Res*. 2015.

---

## License

MIT License. See [LICENSE](LICENSE) file for details.

---

## Data Availability Statement (Template)

```
Raw sequencing data are available at NCBI SRA (PRJNAxxxxxx).
Assembled contigs are available at Zenodo (doi:10.xxxx/zenodo.xxxxxx).
Non-redundant MAGs are available at Zenodo (doi:10.xxxx/zenodo.xxxxxx).
Code is available at: https://github.com/ChengzhuangChen/metagenomics-pipeline
```