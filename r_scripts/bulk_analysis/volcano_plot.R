library(ggplot2)
library(ggrepel)


data <- read.csv("tcga_val_mrna.csv")
data$diffexpressed <- "Insignificant"
data$diffexpressed[data$logFC > 0.5 & data$P.Value < 0.05] <- "Upregulated"
data$diffexpressed[data$logFC < -0.5 & data$P.Value < 0.05] <- "Downregulated"

higher_expression_in <- "(A) TCGA_IDHMUT_HR"
lower_expression_in <- "(B) TCGA_IDHMUT_LR"

sizes <- c("Upregulated" = 1.5, "Downregulated" = 1.5, "Insignificant" = 1)
alphas <- c("Upregulated" = 1, "Downregulated" = 1, "Insignificant" = 0.5)

up_count <- sum(data$diffexpressed == "Upregulated")
down_count <- sum(data$diffexpressed == "Downregulated")
# Volcano plot with customized features
volc_plot <- ggplot(data, aes(x = logFC, y = -log10(P.Value), color = diffexpressed,size =diffexpressed,
                 alpha = diffexpressed)) +
  geom_point() +
  geom_text_repel(data = subset(data, diffexpressed != "Insignificant"),
                  aes(label = genenames), size = 2, max.overlaps = 9) +
  scale_color_manual(
    values = c("Upregulated" = "#FF8C00", "Downregulated" = "#003566", "Insignificant" = "gray30"),
    labels = c("Downregulated", "Insignificant", "Upregulated")
  ) +
  labs(
    x = "Log2 Fold Change",
    y = "-log10(p-Value)"
  ) +
  geom_vline(xintercept = c(-0.5, 0.5), col = "red", linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), col = "red", linetype = "dashed") +
  theme_minimal() +
  scale_size_manual(values = sizes) + # Modify point size
  scale_alpha_manual(values = alphas) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),  # Center title and style it
    legend.title = element_blank()  # Remove legend title
  ) +
  annotate("text", x = -max(data$logFC) * 2, y = -1,
           label = paste("←", lower_expression_in), hjust = 0, size = 2) +
  annotate("text", x = max(data$logFC) * 2, y = -1,
           label = paste(higher_expression_in, "→"), hjust = 1, size = 2) +
  xlim(-2, 2)  +
  # Add annotations for up and downregulated gene counts
  annotate("text", x = -2, y = 1,
           label = paste0("Downregulated: ", down_count), size = 3, hjust = 0) +
  annotate("text", x = 2, y = 1,
           label = paste0("Upregulated: ", up_count), size = 3, hjust = 1)
print(volc_plot)
ggsave(filename = "volc_plot.pdf", plot = volc_plot,
       width = 16, height = 12, dpi = 1000, units = "cm")
