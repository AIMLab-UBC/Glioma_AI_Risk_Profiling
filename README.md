Integrative AI profiling redefines risk in infiltrating gliomas
===========
<details>
<summary>
  <b>Integrative AI profiling redefines risk in infiltrating gliomas</b>, JOURNAL.
  <a href="link" target="blank">[HTML]</a>
  <br><em>Katherine Rich, Karina C. Martin, Thomas Roetzer-Pejrimovsky, Kira Tosefsky, Verena Goebeler, Agustina Conrerro, Lesley Hill, Pouya Ahmadvand, Crystal Ma, Hossein Farahani, Steven JM Jones, Michael Underhill, Adelheid Woehrer, Philipp F. Lange, Stephen Yip, Ali Bashashati</em></br>
</summary>
</details>
---

## 🔍 Overview

**Integrative AI profiling redefines risk in infiltrating gliomas** is a novel study for evaluating whether histological features from Whole Slide Images (WSIs) can predict risk in IDH-mutant gliomas. Evaluated on two external datasets, we identified a high-risk subset of IDH-mutant gliomas characterized by distinct molecular patterns in bulk transcriptomic data, and further corroborated by spatial transcriptomics and proteomic analyses. Moreover, this study establishes a framework for leveraging deep learning in clinical discovery, with broad potential applications across other cancers and complex diseases.

---

## 💾 Dataset & Features

We use WSI images and clinical data from the TCGA-LGG and TCGA-GBM studies, covering more than 800 patients. We use a number of publically available pathology foundation models (including CTransPath, SimCLR, UNI, and Prov-Gigapath) to develop a simple MIL risk prediction framework.

Features should be structured in an h5 file. Sample survival data is included, and survival data from other publically available studies can be easily accessed from cBioPortal (https://www.cbioportal.org)

---

## 🗂 Repository Structure

The directory is structured as follows:
<pre lang="markdown">
    ├── r_scripts/   # R code for downstream genomic and survival analysis
    ├── ├── evaluate_mil_model/
    │       ├── calculate_hr.R # R code for computing hazard ratio for risk group, grade, and subtype
    │       ├── c_index.R # R code for computing 5 fold c-index w/ STD and singular c-index based on avg hazard score
    │       ├── findMaxStat.R # R code for computing maxStat threshold
    │       ├── forest_plot.R # R code for plotting multivariate forest plot with risk scores + clinical variables
    │       ├── km_plot.R # R code for plotting survival curves (KM-curves)  
    ├── ├── bulk_analysis/
    │       ├── TCGA_DE.R              # R code for differential expression analysis with TCGA data
    │       ├── enrichment_analysis.R  # R code for enrichment/pathway analysis with TCGA data
    │       ├── volcano_plot.R         # R code for plotting volcano plot with TCGA data    
    ├── ├── spatial_transcriptomics/
    │       ├── rctd.R              # R code for running rctd on xenium data
    │       ├── xenium_clustering.R   # R code for clustering xenium data
    │       ├── xenium_de.R   # R code for running differential expression analysis on xenium data
    │       ├── xenium_haz_roi.R   # R code for adding cell-level hazard scores to xenium data and doing basic plots + analysis
    ├── survival_mil/   # python code for training MIL risk model
    └── README.md                  # Project overview and documentation
</pre>

---

## ⚙️ Getting Started
### 0. Create the Conda environment
To run the code with the correct dependencies, use the provided requirements.txt file to create a conda environment:

```bash
conda create --name <env_name> --file requirements.txt
```
### 1. Clone the Repository

Clone the repository:

```bash
git clone https://github.com/AIMLab-UBC/Glioma_AI_Risk_Profiling
```
### 2. Training an MIL Risk Prediction Model

To train an MIL risk prediction model simply extract features using any available foundation model into h5 files. Then generate splits (see survival_mil/make_splits.py) and ensure survival data is formated like the example in survival_mil/csv_files/tcga_lgggbm_surv.csv. Following this, simply update the run.sh script and:

```bash
cd survival_mil
./run.sh
```

### 3. Evaluating MIL Risk Prediction Model

R-scripts used to evaluate the performance of MIL model using several metrics. Assumes risk scores and survival information are aggregated in a patient-level csv file. R-scripts are provided to calculate c-index and hazard ratio, as well as generate several useful plots such as Kaplan-Meier curves or forest plots.

### 4. Bulk RNAseq Analysis

Code to run bulk analysis on TCGA data (combined LGG + GBM). For analysis, need a csv file with the case-ids of the VALIDATION cases along with their risk groups (High Risk vs Low Risk). Should first run TCGA_DE.R to download the bulk data and run differential expression analysis, data is downloaded from web so need an internet connect. Then run enrichment_analysis.R to get enriched pathways in either risk group and volcano_plot.R  to plot differentially expressed genes.

For further information on using clusterProfiler for enrichment analysis can refer to https://bioinformatics.ccr.cancer.gov/docs/btep-coding-club/CC2023/FunctionalEnrich_clusterProfiler/

### 5. Spatial Transcriptomic Analysis

First need to run Baysor on Xenium data to improve segmentation results. If using Xenium Prime, the default segmentation should be sufficient. Refer to guide: https://www.10xgenomics.com/analysis-guides/using-baysor-to-perform-xenium-cell-segmentation

R scripts to run analysis of spatial transcriptomic analysis can be found in 'r_scripts/spatial_transcriptomics/'. First need to run rctd.R to load in Xenium data, threshold based on transcript count and get cell-type predictions using RCTD. Single-cell reference data can be created using method here: https://raw.githack.com/dmcable/RCTD/master/vignettes/spatial-transcriptomics.html with an applicable single-cell dataset.

xenium_clustering.R runs the clustering pipeline, then projects to full dataset. Cluster markers are computed from full dataset. Cluster cell-types were identified using provided annotations from 10x website. Other tools such as enrichr (https://maayanlab.cloud/Enrichr/) or cellxgene (https://cellxgene.cziscience.com) can also be used.


## Acknowledgements

The survival MIL code is adapted from the [MCAT](https://github.com/mahmoodlab/MCAT) model.

## 📜 Citation

If you use this work, please cite:

```bibtex
@inproceedings{}
```
