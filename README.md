Official repository for:

**Integrative AI profiling redefines risk in infiltrating gliomas**  
_AUTHORS_  
_JOURNAL_

[📄 Paper ]()  


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
    │   ├── coxph/   
    │       ├── multicox.R # R code for plotting the multivariate coxph forest plot
    │       ├── tcga_val_multicoxph_cleaned.csv # example csv data
    ├── ├── EnrichmentAnalysis.R           # Code for pathway enrichment analysis
    ├── ├── km_curves/
    │       ├── c_index_5fold.R # R code for computing 5 fold c-index w/ STD
    │       ├── c_index.R # R code for computing singular c-index
    │       ├── findMaxStat.R # R code for computing maxStat threshold
    │       ├── surv_plot.R # R code for plotting survival curve
    ├── ├── proteomics.R           # R code for proteomic analysis
    ├── ├── TCGA_DE.R              # R code for differential expression analysis with TCGA data
    │   └── XeniumSeurat.R         # R code for analysis of Xenium data
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
### 2. Training an MIL risk prediction model

To train an MIL risk prediction model simply extract features using any available foundation model into h5 files. Then generate splits (see survival_mil/make_splits.py) and ensure survival data is formated like the example in survival_mil/csv_files/tcga_lgggbm_surv.csv. Following this, simply update the run.sh script and:

```bash
cd survival_mil
./run.sh
```


## 📜 Citation

If you use this work, please cite:

```bibtex
@inproceedings{}
```
