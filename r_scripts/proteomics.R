####################
########SET UP######
####################

###
#Set working directory and load libraries
###
setwd("~/")

library(data.table)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(ggrepel)
library(pheatmap)
library(reshape2)
library(Rtsne)
library(tibble)
library(tidyr)
library(umap)
library(VennDiagram)
set.seed(1)
####
#Read Data either imputed or not imputed + metadata of all samples
####

data <- read.csv(".\\LGGStudy.renamedColumns-all-good-samples.csv") #Not imputed

# data <- read.csv(".\\imputedData__kNNk5TOGETHER-only-good-samples.csv") #Imputed
# data <- data %>% select(-matches("_imp"))

metadata <-  read.csv(".\\metadataV2WithSampleColumns.csv")

imputed=F #Change this according to the data you are using. It is used for saving the plots. 
####
#Define Colors, regex and titles 
####

colorCondition1 <- "#2c7da0"
colorCondition2 <- "#023047"

#Side Note: The '^' symbol indicates that the string starts with the following characters. For example, '^Fibro' matches any string that begins with 'Fibro'.
regexCondition1 <- "^LR" 
regexCondition2 <- "^HR"
regexIDCols <- "^PG"
regexOR <- paste(regexCondition1, regexCondition2, sep = "|")

titleCondition1 <- "LowRisk" 
titleCondition2 <- "HighRisk"

print("IMPORTANT: Please ensure that Condition 1 and Condition 2 are consistently defined to maintain reproducibility and accuracy in the analysis.")

#Export results
export=TRUE #Select TRUE or FALSE. If FALSE, plots will not be saved in folders.

#If export == TRUE, Create folder to export results
if(export == TRUE){
  dir.create(".//results//")
  dir.create(".//results//pdf//")
  dir.create(".//results//png//")
} #If folders already exist, a warning will appear. Please ignore it.

data <- as.data.frame(data)
data <- data %>%
  mutate(across(matches(regexOR), as.numeric)) 

####################
########RESULTS#####
####################

samplesToDelete <- c(
  "HR_Slide.1.1_33", "HR_Slide.1.1_34", "HR_Slide.1.1_35", "HR_Slide.1.1_36",
  "HR_Slide.1.1_37", "HR_Slide.1.1_38", "HR_Slide.1.1_39", "HR_Slide.1.1_41",
  "HR_Slide.1.1_42", "HR_Slide.1.1_43", "HR_Slide.1.1_44", "HR_Slide.1.1_45",
  "LR_Slide.2.1_49", "LR_Slide.2.1_37", "LR_Slide.2.1_39", "LR_Slide.2.1_41",
  "LR_Slide.2.1_42", "LR_Slide.2.1_43", "LR_Slide.2.1_44", "LR_Slide.2.1_46",
  "LR_Slide.2.1_47", "LR_Slide.2.1_61", "HR_Slide.2.1_74", "HR_Slide.2.1_77",
  "HR_Slide.2.1_81", "HR_Slide.2.1_84", "HR_Slide.2.1_86", "HR_Slide.2.1_88",
  "HR_Slide.2.1_89", "HR_Slide.2.1_90", "HR_Slide.2.1_91", "HR_Slide.2.1_93",
  "HR_Slide.1.2_44", "HR_Slide.1.1_46", "HR_Slide.1.1_52", "HR_Slide.1.1_54",
  "HR_Slide.1.1_56", "HR_Slide.1.2_53", "LR_Slide.1.2_58", "LR_Slide.1.2_66",
  "LR_Slide.1.2_70", "LR_Slide.1.2_73", "LR_Slide.1.2_77", "LR_Slide.1.2_79",
  "LR_Slide.1.2_85", "LR_Slide.1.2_87", "LR_Slide.2.2_11", "LR_Slide.2.2_16",
  "LR_Slide.2.2_01", "LR_Slide.2.2_21", "LR_Slide.2.2_22", "LR_Slide.2.2_23",
  "LR_Slide.2.2_04", "LR_Slide.2.2_05", "LR_Slide.2.2_06", "LR_Slide.2.2_07",
  "LR_Slide.2.2_25", "LR_Slide.2.2_27", "LR_Slide.2.2_37", "LR_Slide.2.2_38",
  "LR_Slide.2.2_40", "LR_Slide.2.2_41", "LR_Slide.2.2_42", "LR_Slide.2.2_46",
  "LR_Slide.2.2_49", "LR_Slide.2.2_64", "HR_Slide.2.2_65", "HR_Slide.2.2_69",
  "HR_Slide.2.2_73", "HR_Slide.2.1_100", "HR_Slide.2.1_102", "HR_Slide.2.1_97",
  "HR_Slide.2.1_99", "HR_Slide.2.2_100", "HR_Slide.2.2_101", "HR_Slide.2.2_85",
  "HR_Slide.2.2_86", "HR_Slide.2.2_87", "HR_Slide.2.2_90", "HR_Slide.2.2_91",
  "HR_Slide.2.2_92", "HR_Slide.2.2_95", "HR_Slide.2.2_97", "HR_Slide.2.2_99",
  "LR_Slide.2.1_10", "LR_Slide.2.1_11", "LR_Slide.2.1_12", "LR_Slide.2.1_13",
  "LR_Slide.2.1_14", "LR_Slide.2.1_15", "LR_Slide.2.1_17", "LR_Slide.2.1_18",
  "LR_Slide.2.1_19", "LR_Slide.2.1_20", "LR_Slide.2.1_21", "LR_Slide.2.1_22",
  "LR_Slide.2.1_23", "LR_Slide.2.1_24", "LR_Slide.2.1_25", "LR_Slide.2.1_26",
  "LR_Slide.2.1_27", "LR_Slide.2.1_02", "LR_Slide.2.1_05", "LR_Slide.2.1_06",
  "LR_Slide.2.1_07", "LR_Slide.2.1_08", "LR_Slide.2.1_09", "LR_Slide.1.1_61",
  "LR_Slide.2.1_28", "LR_Slide.1.1_78", "LR_Slide.1.1_81", "LR_Slide.1.1_88",
  "LR_Slide.2.2_40.1"
)
metadata <- metadata[!metadata$Sample %in% samplesToDelete, ]
# data <- data[, !colnames(data) %in% samplesToDelete]

table(metadata$Condition)
table(metadata$ROI)












####
#Barplot Protein Counts ###It only makes sense for not imputed
####

#Convert the data to a data frame and Set every numeric column as numeric
data <- as.data.frame(data)
data <- data %>%
  mutate(across(matches(regexOR), as.numeric)) 
#Convert into a binary data (0 for NA values and 1 for actual values) for numeric columns only.
dataBinary <- data
numericCols <- grep(regexOR, colnames(dataBinary)) 
dataBinary[, numericCols] <- apply(dataBinary[, numericCols], 2, function(x) ifelse(x >= 1, 1, 0))
dataBinary <- dataBinary %>%
  mutate(across(matches(regexOR), ~ replace_na(., 0)))

#Count how many proteins there are per sample
proteinCounts <- colSums(dataBinary[,numericCols] == 1)
dataCount <- data.frame(
  Sample = names(proteinCounts),  
  Count = proteinCounts             
)

#Add colors and conditions according to regular expressions
dataCount <- dataCount %>%
  mutate(Color = case_when(
    grepl(regexCondition1, Sample) ~ colorCondition1,
    grepl(regexCondition2, Sample) ~ colorCondition2,
    TRUE ~ "#000000"  #If not regex is found
  ))
dataCount <- dataCount %>%
  mutate(Condition = case_when(
    grepl(regexCondition1, Sample) ~ titleCondition1,
    grepl(regexCondition2, Sample) ~ titleCondition2,
    TRUE ~ "No Match"
  ))

#Calculate median and plot
medianLine <- median(dataCount$Count, na.rm = TRUE)

# Plot
barplot <- ggplot(dataCount, aes(x = Sample, y = Count, fill = Color)) +
  geom_bar(stat = "identity") +
  scale_fill_identity() + 
  labs(title = "Identified Proteins by Sample",
       x = "Sample",
       y = "Number of Identified Proteins") +
  theme_minimal() +
  geom_hline(yintercept = medianLine, linetype = "dashed", color = "black") +  # Comment this line to delete median line
  theme(axis.text.x = element_text(size = 6, angle = 90))

barplot


#Saving
if(export){
  ggsave(paste0(".//results//png//barplot-ProteinsBySampleFiltered_imputed", imputed, ".png"), barplot,
         width = 12, height = 4)
  dev.off()
  pdf(paste0(".//results//pdf//barplot-ProteinsBySampleFiltered_imputed", imputed, ".pdf"), #If pdf file has trouble, delete the pdf file from folder and run these 3 lines again (from pdf to devoff)
      width = 12, height = 4)
  barplot
  dev.off()
}

####
#Violin Plot for Intensities
####

#Create a long format data for ggplotting this
idCols <- grep("^PG", colnames(data), value = TRUE)
longData <- melt(data, id.vars = idCols, 
                 variable.name = "Sample", value.name = "Intensity")

longData$Intensity <- as.numeric(longData$Intensity)
longData$Intensity[longData$Intensity < 1] <- NA

# Assign condition and color
longData <- longData %>%
  mutate(Condition = case_when(
    grepl(regexCondition1, Sample) ~ titleCondition1,
    grepl(regexCondition2, Sample) ~ titleCondition2,
    TRUE ~ "No Match"#If not regex is found
  ))
longData <- longData %>%
  mutate(Color = case_when(
    grepl(regexCondition1, Sample) ~ colorCondition1,
    grepl(regexCondition2, Sample) ~ colorCondition2,
    TRUE ~ "#000000"  #If not regex is found
  ))

dataCondition1 <- subset(longData, Condition == titleCondition1)
dataCondition2 <- subset(longData, Condition == titleCondition2)

#Plot
violinPlot1 <- ggplot(dataCondition1, aes(x = Sample, y = log10(Intensity), fill = Color)) +
  geom_violin(alpha = 0.5, color = "black") +
  scale_fill_identity() +
  ylim(0,6)+
  labs(title = paste0("Log10(Intensity) by Sample - ", titleCondition1),
       x = "Sample",
       y = "Log10(Intensity)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 8))

violinPlot1


#Saving
if(export){
  ggsave(paste0(".//results//png//violinPlot-IntensityDistributionBySample1_imputed", imputed,".png"), violinPlot1,
         width = 10, height = 4)
  dev.off()
  pdf(paste0(".//results//pdf//violinPlot-IntensityDistributionBySample1_imputed", imputed,".pdf"),
      width = 10, height = 4)
  violinPlot1
  dev.off()
}


violinPlot2 <- ggplot(dataCondition2, aes(x = Sample, y = log10(Intensity), fill = Color)) +
  geom_violin(alpha = 0.5, color = "black") +
  scale_fill_identity() +
  ylim(0,6)+
  labs(title = paste0("Log10(Intensity) by Sample - ", titleCondition2),
       x = "Sample",
       y = "Log10(Intensity)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 8))

violinPlot2


if(export){
  ggsave(paste0(".//results//png//violinPlot-IntensityDistributionBySample2_imputed", imputed,".png"), violinPlot2,
         width = 10, height = 4)
  dev.off()
  pdf(paste0(".//results//pdf//violinPlot-IntensityDistributionBySample2_imputed", imputed,".pdf"),
      width = 10, height = 4)
  violinPlot2
  dev.off()
}


####
#CV #This only makes sense for not imputed 
#### 

#Calculate CVs per condition
cvResults <- data.frame(
  PG.ProteinGroups = data$PG.ProteinGroups
)
cvResults[[titleCondition1]] <- apply(
  data[grep(regexCondition1, colnames(data))], 
  1, 
  function(x) sd(x, na.rm = TRUE) / mean(x, na.rm = TRUE) * 100
)
cvResults[[titleCondition2]] <- apply(
  data[grep(regexCondition2, colnames(data))], 
  1, 
  function(x) sd(x, na.rm = TRUE) / mean(x, na.rm = TRUE) * 100
)
cvResults <- na.omit(cvResults)

#Create a long format data
cvLong <- melt(cvResults, id.vars = "PG.ProteinGroups", 
               variable.name = "Condition", value.name = "CV")


#Plot
median_cv_values <- aggregate(CV ~ Condition, data = na.omit(cvLong), FUN = median)
mean_cv_values <- aggregate(CV ~ Condition, data = na.omit(cvLong), FUN = mean)

g <- ggplot(data = na.omit(cvLong), aes(x = Condition, y = CV)) +
  geom_violin(draw_quantiles = c(0.25, 0.75), linetype = "dashed", adjust = 1.5) +
  geom_violin(fill = "transparent", draw_quantiles = 0.5, adjust = 1.5) +
  labs(x = "", y = "%CV") +
  coord_flip() +
  theme_pubr() +
  theme(
    axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
    axis.text.y = element_text(size = 10)
  )

g
# add medians
for (i in 1:nrow(median_cv_values)) {
  median_cv <- median_cv_values[i, "CV"]   
  group <- median_cv_values[i, "Condition"]  
  
  g <- g +
    geom_text(data = data.frame(x = group, y = median_cv, label = round(median_cv, digits = 1)), 
              aes(x = x, y = y, label = label),
              vjust = -0.5, hjust =0, color = "black", size = 3.5)
}

g
#add mean
for (i in 1:nrow(mean_cv_values)) {
  mean_cv <- mean_cv_values[i, "CV"]     
  group <- mean_cv_values[i, "Condition"]  
  g <- g +
    geom_point(data = data.frame(x = group, y = mean_cv), aes(x = x, y = y), 
               color = "red", size = 3, shape = 16) }
g


#Saving
if(export){
  ggsave(paste0(".//results//png//CVplot-allSamples_imputed",imputed,".png"), g,
         width = 8, height = 4)
  dev.off()
  pdf(paste0(".//results//pdf//CVplot-allSamples_imputed",imputed,".pdf"),
      width = 8, height = 4)
  g
  dev.off()
}




###
#PCA
###
data <- data[, !grepl("_imp", colnames(data))]
#Select numeric data and log10 it
data <- as.data.frame(data)
rownames(data) <- data$PG.ProteinGroups
numericData <- data[, grepl(regexOR, colnames(data))]
numericData <- na.omit(numericData)
numericData[numericData <= 0] <- 1e-10
numericDataLog10 <- log10(numericData)

#Run pca
pcaResult <- prcomp(t(numericDataLog10), scale. = TRUE)  
pcaData <- as.data.frame(pcaResult$x)

#Add conditions and sample names
pcaData$Condition <- ""
pcaData$Condition <- ifelse(grepl(regexCondition1, rownames(pcaData)), titleCondition1, pcaData$Condition)
pcaData$Condition <- ifelse(grepl(regexCondition2,  rownames(pcaData)), titleCondition2, pcaData$Condition)
pcaData$Sample <- rownames(pcaData)

#calculate variances in order to add them to the plot
variances <- pcaResult$sdev^2
explainedVariance <- round(100 * (pcaResult$sdev^2 / sum(pcaResult$sdev^2)), 2)

pca1 <- ggplot(pcaData, aes(x = PC1, y = PC2, color = Condition)) +
  geom_point(size = 3) +
  labs(
    title = paste0("PCA log10 Intensities (n=", nrow(numericData), ")"),
    x = paste0("PC1 (", explainedVariance[1], "%)"),
    y = paste0("PC2 (", explainedVariance[2], "%)")
  ) +
  scale_color_manual(values = setNames(c(colorCondition1, colorCondition2), 
                                       c(titleCondition1, titleCondition2))) +
  geom_text_repel(aes(label = Sample), size = 3)+  #Comment this line if you want to delete labeling
  theme_minimal()

pca1


#Saving
if(export){
  ggsave(paste0(".//results//png//PCAlog10_imputed", imputed, ".png"), pca1,
         width = 6, height = 4)
  dev.off()
  pdf(paste0(".//results//pdf//PCAlog10_imputed", imputed, ".pdf"),
      width = 6, height = 4)
  pca1
  dev.off()
}

setdiff(metadata$Sample,colnames(data))
pcaDataAll <- merge(pcaData, metadata, by="Sample", all.x = T)
pcaDataAll <- merge(pcaDataAll, dataCount, by="Sample", all.x = T)

darkColors <- c(
  "#2E86C1", 
  "#34499E", 
  "red4", 
  "orange", 
  "green4", 
  "lightgreen", 
  "purple4", 
  "#8E99AD", 
  "pink", 
  "yellow2"  
)

pca2 <- ggplot(pcaDataAll, aes(x = PC1, y = PC2, color = Patient)) +
  geom_point(size = 3) +
  labs(
    title = paste0("PCA log10 Intensities (n=", nrow(numericData), ")"),
    x = paste0("PC1 (", explainedVariance[1], "%)"),
    y = paste0("PC2 (", explainedVariance[2], "%)")
  ) +
  scale_color_manual(values = darkColors) +
  # geom_text_repel(aes(label = Sample), size = 3)+  #Comment this line if you want to delete labeling
  theme_minimal()

pca2

if(export){
  ggsave(paste0(".//results//png//PCAlog10_byPatients_imputed", imputed, ".png"), pca2,
         width = 6, height = 4)
  dev.off()
  pdf(paste0(".//results//pdf//PCAlog10_byPatients_imputed", imputed, ".pdf"),
      width = 6, height = 4)
  pca2
  dev.off()
}

pca3 <- ggplot(pcaDataAll, aes(x = PC1, y = PC2, color = as.factor(Slide))) +
  geom_point(size = 3) +
  labs(
    title = paste0("PCA log10 Intensities (n=", nrow(numericData), ")"),
    x = paste0("PC1 (", explainedVariance[1], "%)"),
    y = paste0("PC2 (", explainedVariance[2], "%)")
  ) +
  # scale_color_manual(values = darkColors) +
  # geom_text_repel(aes(label = Sample), size = 3)+  #Comment this line if you want to delete labeling
  theme_minimal()

pca3

if(export){
  ggsave(paste0(".//results//png//PCAlog10_bySlide_imputed", imputed, ".png"), pca3,
         width = 6, height = 4)
  dev.off()
  pdf(paste0(".//results//pdf//PCAlog10_bySlide_imputed", imputed, ".pdf"),
      width = 6, height = 4)
  pca3
  dev.off()
}

###
##Heatmap
###

rownames(data) <- data$PG.ProteinGroups
numericData <- data[, grepl(regexOR, colnames(data))]
numericData <- na.omit(numericData)
numericData
numericDataScaled <- t(scale(t(numericData)))

#create clusters for annotation 
classLabels <- ifelse(grepl(regexCondition1, colnames(numericData)), titleCondition1,
                      ifelse(grepl(regexCondition2, colnames(numericData)), titleCondition2, ""))


metadata <- merge(metadata,dataCount, by="Sample")
metadata$Slide <- as.character(metadata$Slide)


annotationData <- data.frame(
  Condition = classLabels,
  Slide = metadata$Slide[match(colnames(numericData), metadata$Sample)],
  Patient = metadata$Patient[match(colnames(numericData), metadata$Sample)]
)

colnames(annotationData) <- c("Condition", "Slide", "Patient")
rownames(annotationData) <- colnames(numericDataScaled)


unique(metadata$CountCategory)
annotationColors <- list(
  Condition = setNames(c(colorCondition1, colorCondition2), c(titleCondition1, titleCondition2)),
  Slide = setNames(c("#F8766D", "#7CAE00", "#00BFC4" ,"#C77CFF"), sort(unique(metadata$Slide))),
  Patient = setNames(c(darkColors[1:5]), sort(unique(annotationData$Patient)))
  # Customize as needed
)

heatmap <- pheatmap(
  as.matrix(numericDataScaled),
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  main = paste0("Z-Scaled Intensity (n=", nrow(numericData), ")"),
  show_rownames = FALSE,
  show_colnames = TRUE,
  color = colorRampPalette(c("blue4", "white", "red4"))(100),
  fontsize_col = 7,
  annotation_col = annotationData, #Comment this and the line below to delete annotation
  annotation_colors = annotationColors,#Comment this and the line above to delete annotation
  breaks = c(seq(min(numericDataScaled), 0, length.out = 50), seq(0.1, max(numericDataScaled), length.out = 50))
)

heatmap


#Saving
if(export){
  ggsave(paste0(".//results//png//heatmap-ZScaled_imputed", imputed, ".png"), heatmap,
         width = 15, height =8)
  dev.off()
  pdf(paste0(".//results//pdf//heatmap-ZScaled_imputed", imputed, ".pdf"),
      width = 15, height = 8)
  heatmap
  dev.off()
}

#the end