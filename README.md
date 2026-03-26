``````md
# Microbiome Bioinformatics Learning Series
## Project 2: IBD Metagenomics Analysis

Whole-genome shotgun (WGS) metagenomics analysis of Inflammatory Bowel Disease (IBD) patients using real NIH Human Microbiome Project 2 (HMP2) data.

## Overview

This project applies a full metagenomics analysis pipeline to 5 HMP2 samples spanning three clinical groups: healthy controls (nonIBD), Ulcerative Colitis (UC), and Crohn's Disease (CD). The analysis covers taxonomic profiling, functional pathway abundance, strain-level tracking, gene family annotation, differential analysis, species-function correlations, and machine learning classification.

> **Data source:** [NIH Human Microbiome Project 2 (iHMP)](https://ibdmdb.org/) — real patient WGS metagenomics data

---

## Samples

| Sample ID | Diagnosis |
|-----------|-----------|
| HSM7CZ2A | nonIBD (healthy control) |
| MSM5LLHV | UC (Ulcerative Colitis) |
| HSM6XRQE | UC (Ulcerative Colitis) |
| CSM9X1ZO | UC (Ulcerative Colitis) |
| CSM5FZ4C | CD (Crohn's Disease) |

---

## Key Findings

- **109 species** detected across 5 samples via taxonomic profiling
- **Bacteroides vulgatus** — strongest positive correlation with IBD-associated metabolic pathways
- **Faecalibacterium prausnitzii** — negative correlations across IBD pathways (consistent with its known anti-inflammatory role and depletion in IBD)
- **Top ML-predictive features:** tRNA charging pathway + peptidoglycan biosynthesis pathway
- **ML classifier (Random Forest, LOOCV) AUC = 0.000** — expected given n=5 with 4:1 class imbalance; documented as a methodological limitation
- **16S rRNA baseline AUC = 0.894** (from companion Project 1, 1,233-sample dataset) — contrast highlights the sample size requirement for ML in microbiome studies

---

## Repository Structure

```
IBD-metagenomics-analysis/
├── scripts/
│   ├── 01_explore_taxonomic_profiles.py      # Alpha diversity, top species per sample
│   ├── 02_functional_pathways.py             # Metabolic pathway abundance analysis
│   ├── 03_compare_16S_vs_metagenomics.py     # 16S vs WGS capability comparison
│   ├── 04_strain_level_analysis.py           # Strain-level tracking
│   ├── 05_gene_family_analysis.py            # Gene family abundance profiles
│   ├── 06_decode_gene_functions.py           # EC enzyme annotation
│   ├── 07_differential_analysis.py           # Statistical differential analysis
│   ├── 08_species_function_correlation.py    # Spearman correlation analysis
│   ├── 09_machine_learning.py                # Random Forest LOOCV classifier
│   └── 10_visualization_dashboard.py         # 6-panel comprehensive dashboard
├── data/                                     # HMP2 data (not tracked — see DATA_DOWNLOAD.md)
├── figures/                                  # Output figures
├── results/                                  # Output tables and CSVs
├── environment.yml
├── DATA_DOWNLOAD.md
├── PROJECT_SUMMARY.md
└── README.md
```

---

## Setup

### 1. Clone the repository

```bash
git clone https://github.com/YOUR_USERNAME/IBD-metagenomics-analysis.git
cd IBD-metagenomics-analysis
```

### 2. Create the conda environment

```bash
conda env create -f environment.yml
conda activate metagenomics
```

### 3. Download the data

See DATA_DOWNLOAD.md for instructions on obtaining the HMP2 WGS data files.

### 4. Run the analysis pipeline

Scripts are numbered and designed to run in order:

```bash
python scripts/01_explore_taxonomic_profiles.py
python scripts/02_functional_pathways.py
python scripts/03_compare_16S_vs_metagenomics.py
python scripts/04_strain_level_analysis.py
python scripts/05_gene_family_analysis.py
python scripts/06_decode_gene_functions.py
python scripts/07_differential_analysis.py
python scripts/08_species_function_correlation.py
python scripts/09_machine_learning.py
python scripts/10_visualization_dashboard.py
```

---

## Analysis Pipeline

| Script | Analysis | Output |
|--------|----------|--------|
| `01` | Taxonomic profiling, alpha diversity | Species tables, diversity plots |
| `02` | Metabolic pathway abundance | Pathway heatmaps |
| `03` | 16S vs WGS resolution comparison | Comparison figures |
| `04` | Strain-level tracking | Strain profiles |
| `05` | Gene family abundance | Gene family tables |
| `06` | EC enzyme annotation | Decoded function tables |
| `07` | Differential analysis (IBD vs nonIBD) | Statistical results, volcano plots |
| `08` | Species–function correlation | Spearman correlation heatmaps |
| `09` | Random Forest LOOCV classification | Feature importances, ROC curve |
| `10` | Visualization dashboard | 6-panel summary figure |

---

## Environment

- **OS:** Ubuntu 22.04 (WSL2)
- **Conda environment:** `metagenomics`
- **Python:** 3.10

---

## Data

Raw data is from the NIH HMP2 (iHMP) study. Due to file size and data use policies, raw data files are not included in this repository. See DATA_DOWNLOAD.md.

---

## Limitations

- **n=5 samples** — insufficient for robust statistical inference; all ML results reflect methodology demonstration only
- **4:1 class imbalance** (IBD:nonIBD) — affects classifier training; AUC = 0.000 is an expected artifact, not a model failure
- Pre-computed HUMAnN/MetaPhlAn profile files were used (not raw FASTQ processing) due to computational constraints

---

## License

MIT License. Data derived from the NIH HMP2 study — see [ibdmdb.org](https://ibdmdb.org/) for data use terms.
```

 ```
