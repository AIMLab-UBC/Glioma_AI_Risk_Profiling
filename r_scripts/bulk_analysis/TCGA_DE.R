library("TCGAbiolinks")
library("limma")
library("edgeR")
library("glmnet")
library("factoextra")
library("FactoMineR")
library("caret")
library("SummarizedExperiment")
library("gplots")
library("survival")
library("survminer")
library("RColorBrewer")
library("gProfileR")
library("genefilter")
library(clusterProfiler)
library(org.Hs.eg.db)
library(dplyr)

GDCprojects = getGDCprojects()

#get TCGA project summaries
head(GDCprojects[c("project_id", "name")])
TCGAbiolinks:::getProjectSummary("TCGA-LGG")
TCGAbiolinks:::getProjectSummary("TCGA-GBM")

query_TCGA = GDCquery(
  project = "TCGA-LGG",
  data.category = "Transcriptome Profiling", # parameter enforced by GDCquery
  experimental.strategy = "RNA-Seq",
  workflow.type = "STAR - Counts")

lgg_res = getResults(query_TCGA) # make results as table
# head(lihc_res) # data of the first 6 patients.
colnames(lgg_res) # columns present in the table

query_TCGA = GDCquery(
  project = "TCGA-GBM",
  data.category = "Transcriptome Profiling", # parameter enforced by GDCquery
  experimental.strategy = "RNA-Seq",
  workflow.type = "STAR - Counts")

gbm_res = getResults(query_TCGA) # make results as table
colnames(gbm_res) # columns present in the table

#should have case_ids in one column and the relevant groups for comparison in another
val_cases <- read.csv("val_groups.csv",stringsAsFactors = FALSE)


lgg_filtered <- lgg_res %>% filter(substr(lgg_res$cases, 1, 12) %in% val_cases$CASEID)
#exclude recurrent tumour samples
lgg_filtered <-lgg_filtered %>% filter(lgg_filtered$sample_type!="Recurrent Tumor")

query_TCGA = GDCquery(
  project = "TCGA-LGG",
  data.category = "Transcriptome Profiling", # parameter enforced by GDCquery
  experimental.strategy = "RNA-Seq",
  data.type="Gene Expression Quantification",
  workflow.type = "STAR - Counts",
  barcode = c(lgg_filtered$cases))
GDCdownload(query = query_TCGA)


gbm_filtered <- gbm_res %>% filter(substr(gbm_res$cases, 1, 12) %in% val_cases$CASEID)

query_TCGA = GDCquery(
  project = "TCGA-GBM",
  data.category = "Transcriptome Profiling", # parameter enforced by GDCquery
  experimental.strategy = "RNA-Seq",
  data.type="Gene Expression Quantification",
  workflow.type = "STAR - Counts",
  barcode = c(gbm_filtered$cases))
GDCdownload(query = query_TCGA)

gbm_data = GDCprepare(query_TCGA)
lgg_data = GDCprepare(query_TCGA)

saveRDS(object =gbm_data,
        file = "gbm_data.RDS",
        compress = FALSE)
saveRDS(object =lgg_data,
        file = "lgg_data.RDS",
        compress = FALSE)

gbm_data=readRDS("gbm_data.RDS")
lgg_data=readRDS("lgg_data.RDS")
exp_lgg_preprocessed <- TCGAanalyze_Preprocessing(
  object = lgg_data,
  filename = "LGG_IlluminaHiSeq_RNASeqV2.png"
)

exp_gbm_preprocessed <- TCGAanalyze_Preprocessing(
  object = gbm_data,
  filename = "GBM_IlluminaHiSeq_RNASeqV2.png"
)
exp_preprocessed <- cbind(
  exp_lgg_preprocessed,
  exp_gbm_preprocessed
)


exp_normalized <- TCGAanalyze_Normalization(
  tabDF = exp_preprocessed,
  geneInfo = TCGAbiolinks::geneInfoHT,
  method = "geneLength"
) # 60513   40

exp_filtered <- TCGAanalyze_Filtering(
  tabDF = exp_normalized,
  method = "quantile",
  qnt.cut =  0.25
)  # 44630   40

exp_filtered_hr <- exp_filtered[
  ,substr(colnames(exp_filtered),1,12) %in% val_cases$CASEID[val_cases$GROUP == "HighRisk"]
]

exp_filtered_lr <- exp_filtered[
  ,substr(colnames(exp_filtered),1,12) %in% val_cases$CASEID[val_cases$GROUP == "LowRisk"]
]

diff_expressed_genes <- TCGAanalyze_DEA(
  mat1 = exp_filtered_lr,
  mat2 = exp_filtered_hr,
  pipeline="limma",
  Cond1type = "LR",
  Cond2type = "HR",
  fdr.cut = 0.05 ,
  #logFC.cut = 0.5,
  method = "glmLRT",
  log.trans = TRUE
)

write.csv(diff_expressed_genes,"tcga_val_mrna_ensembl.csv")
genes <- clusterProfiler::bitr(rownames(diff_expressed_genes),fromType = "ENSEMBL",toType = "SYMBOL",OrgDb = org.Hs.eg.db)
write.csv(genes,"gene_names.csv")

diff_expressed_genes$genenames <- genes$SYMBOL[match(rownames(diff_expressed_genes), genes$ENSEMBL)]

write.csv(diff_expressed_genes,"tcga_val_mrna.csv")
