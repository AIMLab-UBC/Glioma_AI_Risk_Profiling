library(DOSE)
library(GO.db)
library(GSEABase)
library(org.Hs.eg.db)
library(clusterProfiler)
library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)
library(RColorBrewer)
library(rWikiPathways)
library(RCy3)

glioma.mrna <- read.csv("mrna_val.tsv",sep="\t",stringsAsFactors = FALSE)
nrow(glioma.mrna)
head(glioma.mrna)

up.genes <- glioma.mrna[glioma.mrna$Log2.Ratio > 0.5 & glioma.mrna$p.Value < 0.05, 1]
dn.genes <- glioma.mrna[glioma.mrna$Log2.Ratio < -0.5 & glioma.mrna$p.Value < 0.05, 1]
bkgd.genes <- glioma.mrna[,1]

up.genes.entrez <- clusterProfiler::bitr(up.genes,fromType = "SYMBOL",toType = "ENTREZID",OrgDb = org.Hs.eg.db)
cat("\n\nWhich column contains my new Entrez IDs?\n")
head(up.genes.entrez)

keytypes(org.Hs.eg.db)

dn.genes.entrez <- bitr(dn.genes,fromType = "SYMBOL",toType = "ENTREZID",OrgDb = org.Hs.eg.db)
bkgd.genes.entrez <- bitr(bkgd.genes,fromType = "SYMBOL",toType = "ENTREZID",OrgDb = org.Hs.eg.db)
#  ont      = "BP",
egobp <- clusterProfiler::enrichGO(
  gene     = dn.genes.entrez[[2]],
  universe = bkgd.genes.entrez[[2]],
  OrgDb    = org.Hs.eg.db,
  ont      = "BP",
  pAdjustMethod = "fdr",
  pvalueCutoff = 0.05, #p.adjust cutoff (https://github.com/GuangchuangYu/clusterProfiler/issues/104)
  readable = TRUE)

head(egobp,15)

barplot(egobp, showCategory = 15)
dotplot(egobp, showCategory = 15)
goplot(egobp)
