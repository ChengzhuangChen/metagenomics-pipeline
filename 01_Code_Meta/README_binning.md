# Metagenome Binning Pipeline - Tutorial

## Method Summary

Metagenome-assembled genomes (MAGs) were recovered using four independent binning tools and integrated via DAS_Tool.

```
Contigs (≥2000bp) → [MetaBAT2, MaxBin2, CONCOCT, SemiBin2] → DAS_Tool → Non-redundant MAGs
```

---

## Dependencies

### Software Requirements

| Tool | Version | Purpose |
|------|---------|---------|
| MetaBAT2 | 2.12.1 | Coverage-based binning |
| MaxBin2 | 2.2.6 | Probabilistic binning |
| CONCOCT | 1.1.0 | Coverage + composition-based |
| SemiBin2 | 2.2.0 | Semi-supervised binning |
| DAS_Tool | 1.1.6 | Binning integration & dereplication |
| Bowtie2 | 2.4+ | Read mapping |
| SAMtools | 1.14+ | BAM processing |
| seqkit | 2.0+ | Sequence manipulation |

### Installation via Conda

```bash
conda create -n binning -c bioconda -c conda-forge \
    metabat2=2.12.1 \
    maxbin2=2.2.6 \
    concoct=1.1.0 \
    semibin2=2.2.0 \
    das_tool=1.1.6 \
    bowtie2=2.4.5 \
    samtools=1.14 \
    seqkit=2.0.0 \
    hmmer=3.3.2 \
    prodigal=2.6.3 \
    python=3.8

conda activate binning
```

---

## Pipeline Steps

### Step 0: Prepare Input Files

Filter contigs by length and map reads to contigs.

```bash
# Filter contigs (≥2000 bp)
seqkit seq -m 2000 assembly/contigs_1000bp.fa -o contigs_2000bp.fa

# Build Bowtie2 index
bowtie2-build contigs_2000bp.fa bt2_index/contigs

# Map reads
bowtie2 -p 32 -x bt2_index/contigs \
    -1 clean/sample_R1.fastq.gz \
    -2 clean/sample_R2.fastq.gz \
    --very-sensitive \
    -S aligned.sam

# Convert to sorted BAM
samtools view -bS aligned.sam | samtools sort -o aligned_sorted.bam
samtools index aligned_sorted.bam
```

---

### Step 1: MetaBAT2

Coverage-based binning using depth information.

```bash
#!/bin/bash
# binning_metabat2.sh

CONTIGS="contigs_2000bp.fa"
BAM="aligned_sorted.bam"
OUTPUT="bins_metabat2"

mkdir -p ${OUTPUT}

# Calculate coverage depth
jgi_summarize_bam_contig_depths --outputDepth ${OUTPUT}/depth.txt ${BAM}

# Run MetaBAT2
metabat1 -t 32 -m 1500 -o ${OUTPUT}/bin ${CONTIGS} ${OUTPUT}/depth.txt

echo "MetaBAT2 complete: $(find ${OUTPUT} -name 'bin.*.fa' | wc -l) bins"
```

**Parameters:**
- `-t`: Number of threads
- `-m 1500`: Minimum contig length (bp)
- `-o`: Output prefix

---

### Step 2: MaxBin2

Probabilistic binning using tetranucleotide frequencies and read coverage.

```bash
#!/bin/bash
# binning_maxbin2.sh

CONTIGS="contigs_2000bp.fa"
R1="clean/sample_R1.fastq.gz"
R2="clean/sample_R2.fastq.gz"
OUTPUT="bins_maxbin2"

mkdir -p ${OUTPUT}

# Run MaxBin2
run_MaxBin.pl \
    -thread 32 \
    -contig ${CONTIGS} \
    -reads1 ${R1} \
    -reads2 ${R2} \
    -out ${OUTPUT}/bin \
    -min_contig_length 2000

echo "MaxBin2 complete: $(find ${OUTPUT} -name 'bin.*.fasta' | wc -l) bins"
```

**Parameters:**
- `-thread`: Number of threads
- `-contig`: Input contigs
- `-reads1/2`: Paired-end reads
- `-min_contig_length`: Minimum contig length

---

### Step 3: CONCOCT

Bin clustering using both coverage and composition data.

```bash
#!/bin/bash
# binning_concoct.sh

CONTIGS="contigs_2000bp.fa"
BAM="aligned_sorted.bam"
OUTPUT="bins_concoct"

mkdir -p ${OUTPUT}
mkdir -p ${OUTPUT}/bins

# Step 1: Generate coverage table
concoct_coverage_table.py ${CONTIGS} ${BAM} > ${OUTPUT}/coverage.tsv

# Step 2: Cut contigs into pieces (CONCOCT requirement: ≥1000 bp)
cut_up_fasta.py ${CONTIGS} -c 10000 -o 0 > ${OUTPUT}/contigs_cut.fa

# Step 3: Run CONCOCT clustering
concut --composition_file ${OUTPUT}/contigs_cut.fa \
    --coverage_file ${OUTPUT}/coverage.tsv \
    -p 32 \
    -o ${OUTPUT}/

# Step 4: Merge cut contigs back to original lengths
merge_cutup_clustering.py ${OUTPUT}/clustering_gt1000.csv > \
    ${OUTPUT}/clustering_merged.csv

# Step 5: Extract bins
python3 /path/to/extract_fasta_bins.py \
    ${CONTIGS} \
    ${OUTPUT}/clustering_merged.csv \
    --output_path ${OUTPUT}/bins/

echo "CONCOCT complete: $(find ${OUTPUT}/bins -name '*.fa' | wc -l) bins"
```

**Note:** `extract_fasta_bins.py` is included with CONCOCT installation.

---

### Step 4: SemiBin2

Semi-supervised binning using self-training with reference genomes.

```bash
#!/bin/bash
# binning_semibin2.sh

CONTIGS="contigs_2000bp.fa"
BAM="aligned_sorted.bam"
OUTPUT="bins_semibin2"

mkdir -p ${OUTPUT}

# Run SemiBin2 (single-sample mode)
SemiBin2 single \
    -i ${CONTIGS} \
    -b ${BAM} \
    -o ${OUTPUT} \
    -p 32 \
    --min-contig-length 2000

echo "SemiBin2 complete: $(find ${OUTPUT} -name 'bin_*.fa' | wc -l) bins"
```

**Parameters:**
- `single`: Single-sample binning mode
- `-i`: Input contigs
- `-b`: BAM file
- `-o`: Output directory
- `--min-contig-length`: Minimum contig length (bp)

---

### Step 5: DAS_Tool Integration & Dereplication

Integrate bins from all four tools and dereplicate to non-redundant MAGs.

```bash
#!/bin/bash
# dastool_integration.sh

CONTIGS="contigs_2000bp.fa"
OUTPUT="bins_dastool"

mkdir -p ${OUTPUT}

# Step 1: Create bin collection file
# Format: tool_name<TAB>path_to_bins
cat > ${OUTPUT}/bin_collection.tsv << EOF
MetaBAT2	bins_metabat2
MaxBin2	bins_maxbin2
CONCOCT	bins_concoct/bins
SemiBin2	bins_semibin2
EOF

# Step 2: Run DAS_Tool
DAS_Tool \
    -i ${OUTPUT}/bin_collection.tsv \
    -l MetaBAT2,MaxBin2,CONCOCT,SemiBin2 \
    -c ${CONTIGS} \
    -o ${OUTPUT}/DAS_Tool \
    --write_bins \
    -p 32

# Count final MAGs
echo "Non-redundant MAGs: $(find ${OUTPUT}/DAS_Tool_DASToolbins -name '*.fa' | wc -l)"
```

**Parameters:**
- `-i`: Bin collection file
- `-l`: Tool names (comma-separated)
- `-c`: Contigs file
- `-o`: Output prefix
- `--write_bins`: Write final bins to file

---

## Complete Workflow Script

```bash
#!/bin/bash
# 02_Binning_four_tools.sh

set -e

THREADS=32
MIN_CONTIG_LENGTH=2000
WORK_DIR="results/binning"
ASSEMBLY="results/assembly/contigs_1000bp.fa"
R1="results/clean/sample_R1_paired.fastq.gz"
R2="results/clean/sample_R2_paired.fastq.gz"

mkdir -p ${WORK_DIR}

echo "=========================================="
echo "Metagenome Binning Pipeline"
echo "MetaBAT2, MaxBin2, CONCOCT, SemiBin2"
echo "=========================================="

# Step 0: Prepare input
echo "[0/5] Preparing input files..."
seqkit seq -m ${MIN_CONTIG_LENGTH} ${ASSEMBLY} -o ${WORK_DIR}/contigs_2000bp.fa
bowtie2-build ${WORK_DIR}/contigs_2000bp.fa ${WORK_DIR}/bt2_index/contigs
bowtie2 -p ${THREADS} -x ${WORK_DIR}/bt2_index/contigs -1 ${R1} -2 ${R2} --very-sensitive -S ${WORK_DIR}/aligned.sam
samtools view -bS ${WORK_DIR}/aligned.sam | samtools sort -o ${WORK_DIR}/aligned_sorted.bam
samtools index ${WORK_DIR}/aligned_sorted.bam

CONTIGS=${WORK_DIR}/contigs_2000bp.fa
BAM=${WORK_DIR}/aligned_sorted.bam

# Step 1: MetaBAT2
echo "[1/5] MetaBAT2..."
mkdir -p ${WORK_DIR}/metabat2
jgi_summarize_bam_contig_depths --outputDepth ${WORK_DIR}/metabat2/depth.txt ${BAM}
metabat1 -t ${THREADS} -m 1500 -o ${WORK_DIR}/metabat2/bin ${CONTIGS} ${WORK_DIR}/metabat2/depth.txt

# Step 2: MaxBin2
echo "[2/5] MaxBin2..."
mkdir -p ${WORK_DIR}/maxbin2
run_MaxBin.pl -thread ${THREADS} -contig ${CONTIGS} -reads1 ${R1} -reads2 ${R2} -out ${WORK_DIR}/maxbin2/bin -min_contig_length ${MIN_CONTIG_LENGTH}

# Step 3: CONCOCT
echo "[3/5] CONCOCT..."
mkdir -p ${WORK_DIR}/concoct/bins
concoct_coverage_table.py ${CONTIGS} ${BAM} > ${WORK_DIR}/concoct/coverage.tsv
cut_up_fasta.py ${CONTIGS} -c 10000 -o 0 > ${WORK_DIR}/concoct/contigs_cut.fa
concoct --composition_file ${WORK_DIR}/concoct/contigs_cut.fa --coverage_file ${WORK_DIR}/concoct/coverage.tsv -p ${THREADS} -o ${WORK_DIR}/concoct/
merge_cutup_clustering.py ${WORK_DIR}/concoct/clustering_gt1000.csv > ${WORK_DIR}/concoct/clustering_merged.csv
python3 $(which extract_fasta_bins.py) ${CONTIGS} ${WORK_DIR}/concoct/clustering_merged.csv --output_path ${WORK_DIR}/concoct/bins/

# Step 4: SemiBin2
echo "[4/5] SemiBin2..."
mkdir -p ${WORK_DIR}/semibin2
SemiBin2 single -i ${CONTIGS} -b ${BAM} -o ${WORK_DIR}/semibin2 -p ${THREADS} --min-contig-length ${MIN_CONTIG_LENGTH}

# Step 5: DAS_Tool
echo "[5/5] DAS_Tool integration..."
mkdir -p ${WORK_DIR}/dastool
cat > ${WORK_DIR}/dastool/bin_collection.tsv << EOF
MetaBAT2	${WORK_DIR}/metabat2
MaxBin2	${WORK_DIR}/maxbin2
CONCOCT	${WORK_DIR}/concoct/bins
SemiBin2	${WORK_DIR}/semibin2
EOF
DAS_Tool -i ${WORK_DIR}/dastool/bin_collection.tsv -l MetaBAT2,MaxBin2,CONCOCT,SemiBin2 -c ${CONTIGS} -o ${WORK_DIR}/dastool/DAS_Tool --write_bins -p ${THREADS}

# Summary
echo ""
echo "=========================================="
echo "Pipeline Complete!"
echo "=========================================="
echo "MetaBAT2 bins:    $(find ${WORK_DIR}/metabat2 -name 'bin.*.fa' | wc -l)"
echo "MaxBin2 bins:     $(find ${WORK_DIR}/maxbin2 -name 'bin.*.fasta' | wc -l)"
echo "CONCOCT bins:     $(find ${WORK_DIR}/concoct/bins -name '*.fa' | wc -l)"
echo "SemiBin2 bins:    $(find ${WORK_DIR}/semibin2 -name 'bin_*.fa' | wc -l)"
echo "DAS_Tool MAGs:    $(find ${WORK_DIR}/dastool/DAS_Tool_DASToolbins -name '*.fa' | wc -l)"
echo "=========================================="
```

---

## Usage

```bash
# Run complete pipeline
bash 02_Binning_four_tools.sh -a assembly/contigs_1000bp.fa \
                              -1 clean/sample_R1.fastq.gz \
                              -2 clean/sample_R2.fastq.gz \
                              -o results/binning

# With custom parameters
bash 02_Binning_four_tools.sh -a assembly/contigs.fasta \
                              -1 clean/R1.fq.gz \
                              -2 clean/R2.fq.gz \
                              -o binning \
                              -t 64 \
                              -m 2500
```

---

## Output Structure

```
results/binning/
├── contigs_2000bp.fa           # Filtered contigs
├── aligned_sorted.bam          # Read alignments
├── metabat2/                   # MetaBAT2 bins
│   ├── depth.txt
│   └── bin.1.fa, bin.2.fa...
├── maxbin2/                    # MaxBin2 bins
│   ├── bin.001.fasta...
├── concoct/                    # CONCOCT bins
│   └── bins/bin_1.fa, bin_2.fa...
├── semibin2/                   # SemiBin2 bins
│   └── bin_0.fa, bin_1.fa...
└── dastool/                    # Dereplicated MAGs
    ├── DAS_Tool_DASToolbins/   # <-- Final MAGs
    │   ├── MAG.1.fa, MAG.2.fa...
    └── DAS_Tool_summary.txt
```

---

## Expected Results

Based on typical metagenomic datasets:

| Metric | Expected Range |
|--------|----------------|
| MetaBAT2 bins | 200-2000 |
| MaxBin2 bins | 100-800 |
| CONCOCT bins | 150-1500 |
| SemiBin2 bins | 200-2500 |
| **Total preliminary MAGs** | **500-5000** |
| **DAS_Tool non-redundant MAGs** | **100-1000** |

---

## Citation

If you use this pipeline, please cite:

1. Kang DD et al. MetaBAT2: an adaptive binning algorithm for robust and efficient genome reconstruction from metagenome assemblies. *PeerJ*. 2019.
2. Wu YW et al. MaxBin2: automated binning of metagenomic contigs with evidence-based recruitment and cross-validation. *Nucleic Acids Res*. 2016.
3. Alneberg J et al. Binning metagenomic contigs by coverage and composition. *Nat Methods*. 2014.
4. Pan S et al. SemiBin2: improving metagenome binning with semi-supervised classification. *Bioinformatics*. 2022.
5. Sieber CMK et al. Recovery of genomes from metagenomes via a dereplication, annotation and standardization workflow. *Nat Microbiol*. 2018.