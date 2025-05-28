#https://www.bioconductor.org/packages/devel/bioc/vignettes/rWikiPathways/inst/doc/Pathway-Analysis.html
if(!"rWikiPathways" %in% installed.packages()){
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  BiocManager::install("rWikiPathways", update = FALSE)
}
library(rWikiPathways)


load.libs <- c(
  "DOSE",
  "GO.db",
  "GSEABase",
  "org.Hs.eg.db",
  "clusterProfiler",
  "dplyr",
  "tidyr",
  "ggplot2",
  "stringr",
  "RColorBrewer",
  "rWikiPathways",
  "RCy3")
options(install.packages.check.source = "no")
options(install.packages.compile.from.source = "never")
if (!require("pacman")) install.packages("pacman"); library(pacman)
p_load(load.libs, update = TRUE, character.only = TRUE)
status <- sapply(load.libs,require,character.only = TRUE)
if(all(status)){
  print("SUCCESS: You have successfully installed and loaded all required libraries.")
} else{
  cat("ERROR: One or more libraries failed to install correctly. Check the following list for FALSE cases and try again...\n\n")
  status
}

library(org.Hs.eg.db)
library(clusterProfiler)

glioma.mrna <- read.csv("tcga_val_mrna.csv",stringsAsFactors = FALSE)
nrow(glioma.mrna)
head(glioma.mrna)

up.genes <- glioma.mrna[glioma.mrna$Log2.Ratio > 1.0 & glioma.mrna$q.Value < 0.05, 8] 
dn.genes <- glioma.mrna[glioma.mrna$Log2.Ratio < -1.0 & glioma.mrna$q.Value < 0.05, 8]
bkgd.genes <- glioma.mrna[,8]

up.genes.entrez <- clusterProfiler::bitr(up.genes,fromType = "SYMBOL",toType = "ENTREZID",OrgDb = org.Hs.eg.db)
cat("\n\nWhich column contains my new Entrez IDs?\n")
head(up.genes.entrez)

keytypes(org.Hs.eg.db)

dn.genes.entrez <- bitr(dn.genes,fromType = "SYMBOL",toType = "ENTREZID",OrgDb = org.Hs.eg.db)
bkgd.genes.entrez <- bitr(bkgd.genes,fromType = "SYMBOL",toType = "ENTREZID",OrgDb = org.Hs.eg.db)

egobp <- clusterProfiler::enrichGO(
  gene     = dn.genes.entrez[[2]],
  universe = bkgd.genes.entrez[[2]],
  OrgDb    = org.Hs.eg.db,
  ont      = "BP",
  pAdjustMethod = "fdr",
  pvalueCutoff = 0.05, #p.adjust cutoff (https://github.com/GuangchuangYu/clusterProfiler/issues/104)
  readable = TRUE)

head(egobp,15)

barplot(egobp, showCategory = 15,order=TRUE)
#dotplot(egobp, showCategory = 15)
#goplot(egobp)


egokegg <- clusterProfiler::enrichKEGG(
  gene     = up.genes.entrez[[2]],
  universe = bkgd.genes.entrez[[2]],
  organism = "hsa",
  pAdjustMethod = "fdr",
  pvalueCutoff = 0.05
  )

barplot(egokegg, showCategory = 15,order=TRUE)