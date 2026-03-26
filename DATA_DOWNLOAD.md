# Data Download Instructions

The raw data for this project comes from the NIH Human Microbiome Project 2 (HMP2), also known as the Integrative Human Microbiome Project (iHMP). Due to file size and data use policies, raw data files are not included in this repository.

---

## Data Source

**Portal:** [https://ibdmdb.org](https://ibdmdb.org)
**Study:** HMP2 — Inflammatory Bowel Disease Multi-omics Database
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

The analysis scripts use pre-computed MetaPhlAn taxonomic profiles and HUMAnN functional profiles provided by the HMP2 consortium. These are much smaller than raw FASTQ files and sufficient to reproduce all results in this project.

### Step 1: Register and Accept Data Use Agreement

Visit [https://ibdmdb.org](https://ibdmdb.org) and accept the data use terms.

### Step 2: Download Files

Navigate to the HMP2 data portal and download the following for each sample:

- `*_genefamilies_cpm.tsv` — Gene family abundances (HUMAnN output)
- `*_level4ec.tsv` — EC enzyme classifications (HUMAnN output)
- `*_pathabundance_cpm.tsv` — Pathway abundances (HUMAnN output)
- `hmp2_metadata.csv` — Sample metadata

### Step 3: Place Files

Place all downloaded files in the `data/` directory:

```
data/
├── MSM5LLHV_genefamilies_cpm.tsv
├── MSM5LLHV_level4ec.tsv
├── MSM5LLHV_pathabundance_cpm.tsv
├── HSM7CZ2A_genefamilies_cpm.tsv
├── HSM7CZ2A_level4ec.tsv
├── HSM7CZ2A_pathabundance_cpm.tsv
├── ... (repeat for all 5 samples)
└── hmp2_metadata.csv
```

---

## Option 2: Download Raw FASTQ Files

If you want to rerun the full pipeline from raw reads (requires significant compute and storage):

**Requirements:**
- 50–200 GB storage per sample
- 16+ GB RAM for HUMAnN
- Hours of compute per sample

**Download via SRA Toolkit:**

```bash
conda install -c bioconda sra-tools
prefetch SRR_ACCESSION
fastq-dump --split-files SRR_ACCESSION
```

Sample-to-SRA accession mappings are available in the HMP2 metadata file at ibdmdb.org.

---

## Running the Full Pipeline from FASTQ

If starting from raw FASTQ files, run these tools before the analysis scripts:

**Quality Control (Trimmomatic):**

```bash
trimmomatic PE -threads 8 \
    sample_R1.fastq.gz sample_R2.fastq.gz \
    sample_R1_trimmed.fastq.gz sample_R1_unpaired.fastq.gz \
    sample_R2_trimmed.fastq.gz sample_R2_unpaired.fastq.gz \
    ILLUMINACLIP:adapters.fa:2:30:10 LEADING:3 TRAILING:3 MINLEN:36
```

**Host Read Removal (Bowtie2):**

```bash
bowtie2 -x human_genome_index \
    -1 sample_R1_trimmed.fastq.gz \
    -2 sample_R2_trimmed.fastq.gz \
    --un-conc-gz sample_host_removed.fastq.gz \
    -S /dev/null
```

**Taxonomic Profiling (MetaPhlAn):**

```bash
metaphlan sample_host_removed_1.fastq.gz,sample_host_removed_2.fastq.gz \
    --input_type fastq \
    --nproc 8 \
    -o sample_metaphlan_profile.txt
```

**Functional Profiling (HUMAnN):**

```bash
humann --input sample_host_removed_1.fastq.gz \
    --output humann_output/ \
    --threads 8 \
    --taxonomic-profile sample_metaphlan_profile.txt
```

**Merge Profiles Across Samples:**

```bash
humann_join_tables --input humann_output/ \
    --output merged_pathabundances.tsv \
    --file_name pathabundance

merge_metaphlan_tables.py *_metaphlan_profile.txt > merged_metaphlan.tsv
```

---

## File Size Reference

| File | Approximate Size |
|------|-----------------|
| Raw FASTQ (per sample, paired) | 2–10 GB |
| MetaPhlAn merged profile | < 5 MB |
| HUMAnN pathway abundances (merged) | 10–100 MB |
| HUMAnN gene families (merged) | 500 MB – 2 GB |
| Sample metadata | < 1 MB |

---

## Citation

Lloyd-Price, J., Arze, C., Ananthakrishnan, A.N. et al. **Multi-omics of the gut microbial ecosystem in inflammatory bowel diseases.** *Nature* 569, 655–662 (2019). [https://doi.org/10.1038/s41586-019-1237-9](https://doi.org/10.1038/s41586-019-1237-9)
