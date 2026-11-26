# Metagenomic Taxonomic–Functional Redundancy and Machine Learning Reveal Temporal Succession and Regional Signatures of the Koji Microbiome

[![Journal](https://img.shields.io/badge/Journal-Food%20Chemistry-blue.svg)](https://www.sciencedirect.com/journal/food-chemistry)
[![Repository](https://img.shields.io/badge/GitHub-Tang2026FC-181717.svg?logo=github)](https://github.com/hbtang25/Tang2026FC)
[![Data Accession](https://img.shields.io/badge/GSA-PRJCA035900-green.svg)](https://ngdc.cncb.ac.cn/gsa/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

## Overview

This repository contains the bioinformatic pipelines, statistical analysis scripts, and visualization workflows supporting the manuscript submitted to ***Food Chemistry***:

> **Title:** Metagenomic taxonomic–functional redundancy and machine learning reveal temporal succession and regional signatures of the koji microbiome  
> **Authors:** Hongbiao Tang, Mouyong Zou, GuoBang Zheng, Xingui Zhu, Liqiong Guo\*, Junfang Lin\*  
> **Repository:** `https://github.com/hbtang25/Tang2026FC`

The codebase integrates **16S rRNA gene and ITS1 amplicon profiling**, **conditional-independence ecological network inference**, **Random Forest geographical discrimination**, and **whole-genome shotgun (WGS) metagenomic functional stratification** to characterize microbial succession (0–40 h) and regional microbial signatures across solid-state koji fermentations.

---

## Data Availability

All raw sequencing data generated in this study have been deposited in the **Genome Sequence Archive (GSA)** at the National Genomics Data Center (NGDC) / Beijing Institute of Genomics (BIG) Data Center, Chinese Academy of Sciences, under BioProject accession **`PRJCA035900`**:

* **Amplicon Sequencing (16S rRNA V3–V4 & Fungal ITS1):** Accession **`CRA022829`**
* **Whole-Genome Shotgun (WGS) Metagenomics:** Accession **`CRA022828`**
* **Public Access Portal:** [https://ngdc.cncb.ac.cn/gsa/](https://ngdc.cncb.ac.cn/gsa/) (or [http://bigd.big.ac.cn/gsa](http://bigd.big.ac.cn/gsa))

> **Note on Reproducibility:** Large raw FASTQ files, intermediate alignment files, and reference databases (`RDP`, `UNITE`, `Homo_sapiens`, `Glycine_max`, `UniRef90`, `MetaPhlAn2`) are excluded from this repository due to file size constraints and should be retrieved from GSA (`PRJCA035900`) or the respective database repositories prior to execution.

---

## Repository Structure

```text
Tang2026FC/
├── README.md                     # Project overview, pipeline instructions, and reproducibility guide
├── amplicon/
│   ├── 16S_moromi.sh             # Bacterial 16S rRNA (V3–V4) QC, VSEARCH/USEARCH denoising & rarefaction
│   └── ITS_moromi.sh             # Fungal ITS1 QC, VSEARCH/USEARCH denoising & rarefaction
├── diversity_phylogeny/
│   ├── Tree.qmd                  # Phylogenetic cladogram & abundance ring preparation (GraPhlAn)
│   ├── tree.sh                   # GraPhlAn cladogram rendering pipeline
│   ├── Heatmap.qmd               # Differential ASV hierarchical clustering heatmaps (Wilcoxon + FDR)
│   └── Volcano.qmd               # Pairwise differential species volcano plots
├── network/
│   └── Network.qmd               # SPIEC-EASI (MB + StARS) network inference & igraph topology
├── machine_learning/
│   └── RandomForest.qmd          # Regional Random Forest classification, 5-fold rfcv & ROC evaluation
└── metagenomics/
│    ├── WGS.sh                   # KneadData host depletion, MetaPhlAn2 & HUMAnN2 EC stratification
│    └── WGS.qmd                  # KEGG Level 1–3 Welch's t-tests & species-enzyme contribution pies
├── data/
    ├── otutab_bac.txt            # ASV abundance matrix for bacteria
    ├── otutab_fun.txt            # ASV abundance matrix for fungi
    ├── sum_g_bac.txt             # genus abundance matrix for bacteria
    ├── sum_g_fun.txt             # genus abundance matrix for fungi
    ├── taxonomy_bac.txt          # taxonomy for bacterial ASV
    ├── taxonomy_fun.txt          # taxonomy for fungal ASV
    ├── metadata_bac.tsv          # Sample metadata for 16S
    ├── metadata_fun.tsv          # Sample metadata for ITS
    ├── metadata_WGS.tsv          # Sample metadata for metagenomics
```

---

## Software & Dependencies

### 1. Command-Line Bioinformatic Tools
* **Amplicon Processing:** `fastp` (v0.23.2), `FastQC` (v0.12.1), `VSEARCH` (v2.15.1), `USEARCH` (v10.0.240)
* **Phylogenetic Visualization:** `GraPhlAn` (v1.1.3)
* **Shotgun Metagenomics:** `FastQC` (v0.11.9), `MultiQC` (v1.11), `KneadData` (v0.6.1), `Trimmomatic` (v0.39), `Bowtie2` (v2.3.5), `HUMAnN2` (v2.8.1), `MetaPhlAn2` (v2.7.5), `DIAMOND` (v0.8.2.84)
* **Differential Biomarker Discovery:** `LEfSe` (v1.1.2)

### 2. R Environment (`R >= 4.0.3`)
Key R packages required across the Quarto (`.qmd`) analytical notebooks:
* **Community Ecology & Networks:** `SpiecEasi` (v1.99.0), `igraph` (v1.6.0), `vegan` (v2.6-10), `phyloseq` (v1.52.0)
* **Machine Learning & Evaluation:** `randomForest` (v4.6-14), `caret` (v7.0-1), `pROC` (v1.18.5)
* **Data Manipulation & Visualization:** `tidyverse` (v2.0.0), `dplyr` (v1.1.4), `ggplot2` (v3.5.4), `pheatmap` (v1.0.12), `ggrepel` (v0.9.6), `patchwork` (v1.3.1)

---

## Analytical Workflow

### Step 1: Amplicon Sequence Processing (`16S_moromi.sh`, `ITS_moromi.sh`)
1. Adapter/quality trimming via `fastp` and paired-end read merging (`vsearch --fastq_mergepairs`).
2. Primer stripping, quality filtering (`--fastq_filter`), and sequence dereplication (`--derep_fulllength`) with minimum abundance thresholds (`minuniquesize` > 25 for 16S; > 10 for ITS).
3. Denoising and chimera removal via `usearch -unoise3` to generate high-resolution Amplicon Sequence Variants (ASVs).
4. Taxonomic classification (`sintax` against RDP for 16S and UNITE for ITS), removal of host chloroplast/mitochondrial sequences, and rarefaction to uniform sequencing depth (`14,446` reads for 16S; `71,766` reads for ITS; random seed = `325`).

### Step 2: Co-Occurrence Network Analysis (`Network.qmd`)
1. Compositionally robust ecological association inference using **SPIEC-EASI** (`method = "mb"`, `sel.criterion = "stars"`, `pulsar.params = list(rep.num = 100)`, stability threshold = `0.05`).
2. Calculation of global and node-level topological metrics (node/edge counts, connectance, average degree, average path length, modularity, and betweenness centrality) via `igraph`.

### Step 3: Phylogenetic & Differential Abundance Profiling (`Tree.qmd`, `tree.sh`, `Heatmap.qmd`, `Volcano.qmd`)
1. Identification of stage- and region-specific biomarkers via **LEfSe** (`TSS` normalization to $1 \times 10^6$, Kruskal–Wallis $\alpha < 0.05$) and **Wilcoxon rank-sum tests** with Benjamini–Hochberg FDR adjustment ($p_{\text{adj}} < 0.05$).
2. Generation of taxonomic cladograms (`GraPhlAn`), clustered ASV heatmaps, and pairwise differential volcano plots.

### Step 4: Regional Random Forest Modeling (`RandomForest.qmd`)
1. Genus-level relative abundance matrices partitioned into training (`Train`: 106 16S samples; 98 ITS samples) and held-out testing (`Test`: 56 16S samples; 61 ITS samples) cohorts across four regions (`Jiangmen`, `Sichuan`, `Shandong`, `Korea`).
2. Random Forest training (`ntree = 1000`, with `classwt` applied in the fungal model to account for regional sample imbalance) and 5-fold cross-validation (`rfcv` across five seeds) executed exclusively on the training partition.
3. Evaluation on the held-out testing set using overall accuracy, 95% confidence intervals, Cohen's Kappa, per-class confusion matrix metrics, and multi-class ROC/AUC curves.

### Step 5: Shotgun Metagenomics & Functional Stratification (`WGS.sh`, `WGS.qmd`)
1. Host decontamination using `KneadData` (`Trimmomatic` + `Bowtie2` in `--very-sensitive` mode) against both human (`Homo sapiens`) and soybean (`Glycine max`) reference genomes.
2. Species-level taxonomic profiling via `MetaPhlAn2` and functional gene family quantification against `UniRef90` via `HUMAnN2`.
3. Regrouping of `UniRef90` families to Enzyme Commission categories (`uniref90_level4ec`, CPM-normalized) and partitioning into community-wide (`unstratified`) and species-resolved (`stratified`) enzymatic profiles (`humann2_split_stratified_table`) to quantify potential species-level contributions to predicted metabolic pathways.
4. Pairwise regional comparisons of KEGG Level 3 pathways using two-sided Welch's *t*-tests.

---

## Methodological & Interpretive Scope

To ensure transparent interpretation in accordance with peer-review standards:
* **Taxonomic Resolution:** V3–V4 16S rRNA, ITS1, and standard `MetaPhlAn2` profiling resolve taxa at the genus or species level and do not differentiate the inoculated *Aspergillus oryzae* 3.042 starter strain from indigenous *Aspergillus* strains.
* **Cohort Scope:** Random Forest classification metrics reflect discriminative performance within the sampled cohorts and require external validation across larger, batch- and stage-balanced multi-facility datasets prior to broad authentication deployment.
* **Functional Potential:** Metagenomic KEGG and `HUMAnN2` EC-stratified outputs represent **predicted genomic/enzymatic potential** inferred from DNA sequences rather than direct measurements of *in situ* enzyme activity or volatile flavor compound formation.

---

## Citation

If you use the scripts or analytical workflows in this repository, please cite:

> Tang, H., Zou, M., Zheng, G., Zhu, X., Guo, L., & Lin, J. (2026). Metagenomic taxonomic–functional redundancy and machine learning reveal the temporal succession and regional signatures of the koji microbiome. *Food Chemistry* (Under Revision).

---

## Contact

For questions regarding the datasets or computational scripts, please contact:
* **Hongbiao Tang** (First Author) — College of Food Science, South China Agricultural University
* **Liqiong Guo** (`guolq@scau.edu.cn`) & **Junfang Lin** (`linjf@scau.edu.cn`) — Corresponding Authors, South China Agricultural University
