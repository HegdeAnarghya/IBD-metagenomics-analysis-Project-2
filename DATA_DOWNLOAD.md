```
# Data Download Instructions

The raw data for this project comes from the NIH Human Microbiome Project 2 (HMP2), also known as the Integrative Human Microbiome Project (iHMP). Due to file size and data use policies, raw data files are not included in this repository.

---

## Data Source

**Portal:** NIH Human Microbiome Project 2 — Inflammatory Bowel Disease Multi-omics Database (ibdmdb.org)  
**Study:** HMP2 / iHMP  
**Data type:** Whole-genome shotgun (WGS) metagenomics, pre-computed profiles

---

## Samples Used

| Sample ID | Diagnosis |
|-----------|-----------|
| HSM7CZ2A  | nonIBD    |
| MSM5LLHV  | UC        |
| HSM6XRQE  | UC        |
| CSM9X1ZO  | UC        |
| CSM5FZ4C  | CD        |

---

## Option 1: Download Pre-Computed Profiles (Recommended)

The analysis scripts use pre-computed MetaPhlAn taxonomic profiles and HUMAnN functional profiles. These are much smaller than raw FASTQ files and sufficient to reproduce all results.

### Step 1: Accept the Data Use Agreement

Go to the HMP2 / iHMP data portal (search "ibdmdb" or "iHMP portal") and accept the data use terms before downloading any files.

### Step 2: Navigate to the Files

On the portal, navigate to:

```
Products → WGS → Taxonomic Profiles (MetaPhlAn)
Products → WGS → Functional Profiles (HUMAnN)
Products → Metadata
```

### Step 3: Download the Following Files

| File | Description |
|------|-------------|
| `hmp2_metaphlan_species_abundances.tsv` | Merged MetaPhlAn species abundance table |
| `hmp2_pathabundances_relab.tsv.gz` | HUMAnN pathway abundances (relative) |
| `hmp2_genefamilies.tsv.gz` | HUMAnN gene family abundances (UniRef90) |
| `hmp2_metadata.csv` | Sample metadata with diagnosis labels |

Decompress the `.gz` files after downloading:

```bash
gunzip hmp2_pathabundances_relab.tsv.gz
gunzip hmp2_genefamilies.tsv.gz
```

### Step 4: Place Files in the Data Directory

```
data/
├── hmp2_metaphlan_species_abundances.tsv
├── hmp2_pathabundances_relab.tsv
├── hmp2_genefamilies.tsv
└── hmp2_metadata.csv
```

---

## Option 2: Download Raw FASTQ Files

If you want to rerun the full pipeline from raw reads (requires significant compute and storage).

### Requirements

- 50–200 GB storage per sample
- 16+ GB RAM for HUMAnN
- Several hours of compute per sample

### Download via SRA Toolkit

Sample-to-SRA accession mappings are available in `hmp2_metadata.csv`. Use the SRA toolkit to download by accession:

```bash
# Install SRA toolkit
conda install -c bioconda sra-tools

# Download by SRA accession
prefetch SRR_ACCESSION
fastq-dump --split-files SRR_ACCESSION
```

---

## Option 3: Bulk Download via the Portal File Browser

The HMP2 portal provides a file browser and bulk download tool. Filter by:
- **Data type:** WGS
- **File type:** `taxonomic_profile` or `pathabundance`
- **Sample IDs:** HSM7CZ2A, MSM5LLHV, HSM6XRQE, CSM9X1ZO, CSM5FZ4C

---

## Running the Full Pipeline from FASTQ

If starting from raw FASTQ files, run the following tools before the analysis scripts:

### 1. Quality Control

```bash
trimmomatic PE -threads 8 \
    sample_R1.fastq.gz sample_R2.fastq.gz \
    sample_R1_trimmed.fastq.gz sample_R1_unpaired.fastq.gz \
    sample_R2_trimmed.fastq.gz sample_R2_unpaired.fastq.gz \
    ILLUMINACLIP:adapters.fa:2:30:10 LEADING:3 TRAILING:3 MINLEN:36
```

### 2. Host Read Removal

```bash
bowtie2 -x human_genome_index \
    -1 sample_R1_trimmed.fastq.gz \
    -2 sample_R2_trimmed.fastq.gz \
    --un-conc-gz sample_host_removed.fastq.gz \
    -S /dev/null
```

### 3. Taxonomic Profiling (MetaPhlAn)

```bash
metaphlan sample_host_removed_1.fastq.gz,sample_host_removed_2.fastq.gz \
    --input_type fastq \
    --nproc 8 \
    -o sample_metaphlan_profile.txt
```

### 4. Functional Profiling (HUMAnN)

```bash
humann --input sample_host_removed_1.fastq.gz \
    --output humann_output/ \
    --threads 8 \
    --taxonomic-profile sample_metaphlan_profile.txt
```

### 5. Merge Profiles Across Samples

```bash
# Merge HUMAnN pathway tables
humann_join_tables \
    --input humann_output/ \
    --output merged_pathabundances.tsv \
    --file_name pathabundance

# Merge MetaPhlAn tables
merge_metaphlan_tables.py *_metaphlan_profile.txt > merged_metaphlan.tsv
```

---

## File Size Reference

| File | Approximate Size |
|------|-----------------|
| Raw FASTQ per sample (paired) | 2–10 GB |
| MetaPhlAn merged profile | < 5 MB |
| HUMAnN pathway abundances (merged) | 10–100 MB |
| HUMAnN gene families (merged) | 500 MB – 2 GB |
| Sample metadata | < 1 MB |

---

## Citation

Lloyd-Price, J., Arze, C., Ananthakrishnan, A.N. et al. **Multi-omics of the gut microbial ecosystem in inflammatory bowel diseases.** *Nature* 569, 655–662 (2019). DOI: 10.1038/s41586-019-1237-9
```

