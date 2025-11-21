# Cran packages
cran.pkgs <- c("tidyverse", "ashr", "ggplot2", "gplots", "RColorBrewer", "readxl", "pheatmap", "gridExtra", "grid", "extrafont", "ggalt")
options(repos = c(CRAN = "https://cloud.r-project.org"))
for (pkg in cran.pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
}

# Bioconductor packages
bioc.pkgs <- c("DESeq2", "EnhancedVolcano")
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

for (pkg in bioc.pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    BiocManager::install(pkg, ask = FALSE)
  }
}

# Load required libraries
library(tidyverse)
library(DESeq2)
library(EnhancedVolcano)
library(ggplot2)
library(gplots)
library(RColorBrewer)
library(readxl)
library(pheatmap)
library(gridExtra)
library(grid)
library(extrafont)
library(ggalt)

# font_import()
loadfonts()
loadfonts(device = "pdf")
police <- "Arial"
fontsize <- 20


# Load data----

metadata <- read.csv("HRG1_K562_sampleInfo.csv", row.names = 1)
counts <- read.csv("HRG1_K562_counts.csv", row.names = 1)
MitoCartaHuman <- read_xls("Human.MitoCarta3.0.xls", sheet = "A Human MitoCarta3.0")
MitoCartaHumanGenes <- MitoCartaHuman$Symbol

# KO_U vs WT_U----

# Construct a DESeqDataSet object
metadata$Genotype_treatment <- relevel(factor(metadata$Genotype_treatment), ref = "WT_Ctl") # set reference 
dds <- DESeqDataSetFromMatrix(countData = counts,  
                              colData = metadata,
                              design = ~ Genotype_treatment)
#  Filter Out Low-Count Genes
dds<-dds[rowSums(counts(dds))>1,]

# Run DESeq2
dds <- DESeq(dds)

# Extract Differential Expression Results
res_U <- results(dds, contrast = c("Genotype_treatment", "KO_Ctl", "WT_Ctl"), pAdjustMethod="BH")

# lfcShrink ashr helps stabilize the fold changes for genes with low counts
res_U <- lfcShrink(dds, coef = "Genotype_treatment_KO_Ctl_vs_WT_Ctl", type = "ashr")
res_U.mito <- as.data.frame(res_U)[rownames(res_U) %in% MitoCartaHumanGenes, ]

# MA plots display log fold changes (LFC) vs. mean expression levels
pdf("Plots/MAplot_KO_Ctl_vs_WT_Ctl.pdf", width = 8, height = 6)
# par(family = police, font=2)
plotMA(res_U, main = "WT vs KO in Control (Ctl)", ylim = c(-5, 5))
dev.off()


# Volcano plot
pdf("Plots/volcano_plot_KO_Ctl_vs_WT_Ctl.pdf", width = 8, height = 6)#, family = police)
print(
EnhancedVolcano(res_U,
                lab = rownames(res_U),
                x = 'log2FoldChange',
                y = 'padj',
                pCutoff = 0.05,
                FCcutoff = 1, 
                title = "KO_Ctl vs WT_Ctl",
                subtitle = "Ref. WT_Ctl",
                caption = "Differential Expression",
                colAlpha = 0.7,
                legendPosition = 'top',
                legendLabSize = 10)
)
dev.off()


# KO_SB vs WT_SB----

# Construct a DESeqDataSet object
metadata$Genotype_treatment <- relevel(factor(metadata$Genotype_treatment), ref = "WT_SB")
dds <- DESeqDataSetFromMatrix(countData = counts,  
                              colData = metadata,
                              design = ~ Genotype_treatment)
#  Filter Out Low-Count Genes
dds<-dds[rowSums(counts(dds))>1,]

# Run DESeq2
dds <- DESeq(dds)
res_SB <- results(dds, contrast = c("Genotype_treatment", "KO_SB", "WT_SB"), pAdjustMethod="BH")
res_SB <- lfcShrink(dds, coef = "Genotype_treatment_KO_SB_vs_WT_SB", type = "ashr")
res_SB.mito <- as.data.frame(res_SB)[rownames(res_SB) %in% MitoCartaHumanGenes, ]


pdf("Plots/MAplot_KO_SB_vs_WT_SB.pdf", width = 8, height = 6)
# par(family = police, font=2)
plotMA(res_SB, main = "WT vs KO in Treated (SB) Condition", ylim = c(-5, 5))
dev.off()

# Volcano plot
pdf("Plots/volcano_plot_KO_SB_vs_WT_SB.pdf", width = 8, height = 6)#, family = police)
print(
  EnhancedVolcano(res_SB,
                lab = rownames(res_SB),
                x = 'log2FoldChange',
                y = 'padj',
                pCutoff = 0.05,
                FCcutoff = 1,
                title = "KO_SB vs WT_SB",
                subtitle = "Ref. WT_SB",
                caption = "Differential Expression",
                colAlpha = 0.7,
                legendPosition = 'top',
                legendLabSize = 10)
)
dev.off()


# KO_SB vs KO_U----

# Construct a DESeqDataSet object
metadata$Genotype_treatment <- relevel(factor(metadata$Genotype_treatment), ref = "KO_Ctl") # set reference 

dds <- DESeqDataSetFromMatrix(countData = counts, #change count matrix if needed 
                              colData = metadata,
                              design = ~ Genotype_treatment)
#  Filter Out Low-Count Genes
dds<-dds[rowSums(counts(dds))>1,]

# Run DESeq2
dds <- DESeq(dds)

# Extract Differential Expression Results
res_koSB_koU <- results(dds, contrast = c("Genotype_treatment", "KO_SB", "KO_Ctl"), pAdjustMethod="BH")


# lfcShrink ashr helps stabilize the fold changes for genes with low counts

res_koSB_koU <- lfcShrink(dds, coef = "Genotype_treatment_KO_SB_vs_KO_Ctl", type = "ashr")
res_koSB_koU.mito <- as.data.frame(res_koSB_koU)[rownames(res_koSB_koU) %in% MitoCartaHumanGenes, ]

# MA plots display log fold changes (LFC) vs. mean expression levels
pdf("Plots/MAplot_KO_SB_vs_KO_Ctl.pdf", width = 8, height = 6)
# par(family = police, font=2)
plotMA(res_koSB_koU, main = "Treated vs Untreated in KO Condition", ylim = c(-5, 5))
dev.off()


# Volcano plot
pdf("Plots/volcano_plot_KO_SB_vs_KO_Ctl.pdf", width = 8, height = 6)#, family = police)
print(
EnhancedVolcano(res_koSB_koU,
                lab = rownames(res_koSB_koU),
                x = 'log2FoldChange',
                y = 'padj',
                pCutoff = 0.05,
                FCcutoff = 1, 
                title = "KO_SB vs KO_Ctl",
                subtitle = "Ref. KO_Ctl",
                caption = "Differential Expression",
                colAlpha = 0.7,
                legendPosition = 'top',
                legendLabSize = 10)
)
dev.off()


# Save FC tables
KO_SB_vs_KO_U_all <- as.data.frame(res_koSB_koU)

KO_SB_vs_KO_U_filter <- KO_SB_vs_KO_U_all[KO_SB_vs_KO_U_all$log2FoldChange > 3 | KO_SB_vs_KO_U_all$log2FoldChange < -3 | (-log10(KO_SB_vs_KO_U_all$padj) > 10), ]

KO_SB_vs_KO_U_filter = KO_SB_vs_KO_U_filter[!is.na(KO_SB_vs_KO_U_filter),]
# check some genes to see if the filtering is ok
# KO_U_vs_WT_U_filter[rownames(KO_U_vs_WT_U_filter ) == 'FZD3',]

write.csv(KO_SB_vs_KO_U_all, "FC_table/KO_SB_vs_KO_U_all.csv", row.names = TRUE)
write.csv(KO_SB_vs_KO_U_filter, "FC_table/KO_SB_vs_KO_U_filter.csv", row.names = TRUE)

# WT_SB vs WT_U----

# Construct a DESeqDataSet object
metadata$Genotype_treatment <- relevel(factor(metadata$Genotype_treatment), ref = "WT_Ctl") # set reference 

dds <- DESeqDataSetFromMatrix(countData = counts, #change count matrix if needed 
                              colData = metadata,
                              design = ~ Genotype_treatment)
#  Filter Out Low-Count Genes
dds<-dds[rowSums(counts(dds))>1,]

# Run DESeq2
dds <- DESeq(dds)

# Extract Differential Expression Results
res_wtSB_wtU <- results(dds, contrast = c("Genotype_treatment", "WT_SB", "WT_Ctl"), pAdjustMethod="BH")


# lfcShrink ashr helps stabilize the fold changes for genes with low counts

res_wtSB_wtU <- lfcShrink(dds, coef = "Genotype_treatment_WT_SB_vs_WT_Ctl", type = "ashr")
res_wtSB_wtU.mito <- as.data.frame(res_wtSB_wtU)[rownames(res_wtSB_wtU) %in% MitoCartaHumanGenes, ]


# MA plots display log fold changes (LFC) vs. mean expression levels
pdf("Plots/MAplot_WT_SB_vs_WT_Ctl.pdf", width = 8, height = 6)
# par(family = police, font=2)
plotMA(res_wtSB_wtU, main = "Treated vs Untreated in WT Condition", ylim = c(-5, 5))
dev.off()


# Volcano plot
pdf("Plots/volcano_plot_WT_SB_vs_WT_Ctl.pdf", width = 8, height = 6)#, family = police)
print(
  EnhancedVolcano(res_wtSB_wtU,
                lab = rownames(res_wtSB_wtU),
                x = 'log2FoldChange',
                y = 'padj',
                pCutoff = 0.05,
                FCcutoff = 1, 
                title = "WT_SB vs WT_Ctl",
                subtitle = "Ref. WT_Ctl",
                caption = "Differential Expression",
                colAlpha = 0.7,
                legendPosition = 'top',
                legendLabSize = 10))
dev.off()


# Save FC tables
WT_SB_vs_WT_U_all <- as.data.frame(res_wtSB_wtU)
WT_SB_vs_WT_U_filter <- WT_SB_vs_WT_U_all[WT_SB_vs_WT_U_all$log2FoldChange > 3 | WT_SB_vs_WT_U_all$log2FoldChange < -3 | (-log10(WT_SB_vs_WT_U_all$padj) > 10), ]
WT_SB_vs_WT_U_filter = WT_SB_vs_WT_U_filter[!is.na(WT_SB_vs_WT_U_filter),]
# check some genes to see if the filtering is ok
# KO_U_vs_WT_U_filter[rownames(KO_U_vs_WT_U_filter ) == 'FZD3',]
write.csv(WT_SB_vs_WT_U_all, "FC_table/WT_SB_vs_WT_U_all.csv", row.names = TRUE)
write.csv(WT_SB_vs_WT_U_filter, "FC_table/WT_SB_vs_WT_U_filter.csv", row.names = TRUE)



# Plot KO_SB/KO_U vs WT_SB/ST_U ----
data.to.plot <- dplyr::full_join(KO_SB_vs_KO_U_all %>% rownames_to_column("Gene"), 
                                WT_SB_vs_WT_U_all %>% rownames_to_column("Gene"), 
                                by = "Gene", 
                                suffix = c("_KO_SB_vs_KO_Ctl", "_WT_SB_vs_WT_Ctl"))

up_down_treshlod <- 1

# Create subdataframe for each corner of the future plots
data_up_up <- subset(data.to.plot,log2FoldChange_KO_SB_vs_KO_Ctl > up_down_treshlod & log2FoldChange_WT_SB_vs_WT_Ctl > up_down_treshlod)
data_up_down <- subset(data.to.plot,log2FoldChange_KO_SB_vs_KO_Ctl > up_down_treshlod & log2FoldChange_WT_SB_vs_WT_Ctl < -up_down_treshlod)
data_down_up <- subset(data.to.plot,log2FoldChange_KO_SB_vs_KO_Ctl < -up_down_treshlod & log2FoldChange_WT_SB_vs_WT_Ctl > up_down_treshlod)
data_down_down <- subset(data.to.plot,log2FoldChange_KO_SB_vs_KO_Ctl < -up_down_treshlod & log2FoldChange_WT_SB_vs_WT_Ctl < -up_down_treshlod)
data_middle <- data.to.plot[!(data.to.plot$Gene %in% c(data_up_up$Gene, data_up_down$Gene, data_down_up$Gene, data_down_down$Gene)), ]


# Create a matrix with the 4 values (counts) in a 2x2 table
table_data <- matrix(
  c(nrow(data_down_up), nrow(data_down_down), nrow(data_up_up), nrow(data_up_down)),
  nrow = 2, ncol = 2,
  dimnames = list(NULL, NULL)  # No row and column names
)
# Color for text in table (based on groups)
text_colors <- c(
  'orangered', # for data_up_up
  'turquoise', # for data_down_down
  'turquoise', # for data_down_up
  'orangered'  # for data_up_down
)

# Create a table theme for styling (no rownames or colnames)
table_theme <- ttheme_minimal(
  core = list(
    fg_params = list(fontsize = fontsize, fontface = "bold", col = text_colors),  # Set colors for text
    bg_params = list(fill = "white"),
    lwd = 2,
    lty = 2,
    col = "black"  # Border color
  ),
  colhead = list(
    fg_params = list(fontsize = fontsize, fontface = "bold", col = "black"),
    bg_params = list(fill = "white")
  ),
  rowhead = list(
    fg_params = list(fontsize = fontsize, fontface = "bold", col = "black"),
    bg_params = list(fill = "white")
  )
)


pdf("Plots/WT_SB_WT_Ctl_vs_KO_SB_KO_Ctl.pdf", width = 8, height = 6)#, family = police)
# Create the plot
print(
  ggplot(data_middle, aes(x = log2FoldChange_WT_SB_vs_WT_Ctl, y = log2FoldChange_KO_SB_vs_KO_Ctl)) +
  geom_point(color="black", alpha=0.5, size = 5) +
  geom_point(data = data_up_up,
             aes(x = log2FoldChange_WT_SB_vs_WT_Ctl, y = log2FoldChange_KO_SB_vs_KO_Ctl),
             color = 'turquoise', alpha=0.5, size=5) +
  geom_point(data = data_down_down,
             aes(x = log2FoldChange_WT_SB_vs_WT_Ctl, y = log2FoldChange_KO_SB_vs_KO_Ctl), 
             color = 'turquoise', alpha=0.5, size=5) +
  geom_point(data = data_down_up,
             aes(x = log2FoldChange_WT_SB_vs_WT_Ctl, y = log2FoldChange_KO_SB_vs_KO_Ctl), 
             color = 'orangered', alpha=0.5, size=5) +
  geom_point(data = data_up_down, 
             aes(x = log2FoldChange_WT_SB_vs_WT_Ctl, y = log2FoldChange_KO_SB_vs_KO_Ctl), 
             color = 'orangered', alpha=0.5, size=5) +
  
  geom_hline(yintercept = up_down_treshlod, linetype = "dashed", color = "black") +
  geom_hline(yintercept = -up_down_treshlod, linetype = "dashed", color = "black") +
  geom_vline(xintercept = -up_down_treshlod, linetype = "dashed", color = "black") +
  geom_vline(xintercept = up_down_treshlod, linetype = "dashed", color = "black") +
  
  labs(x = bquote(Log[2] * "FC (KO SB vs KO Ctl)"), 
       y = bquote(Log[2] * "WT SB vs WT Ctl")) +
  
  scale_x_continuous() +
  scale_y_continuous() +
  
  theme_minimal() +
  theme(
    legend.position = "none",  
    panel.grid.major = element_blank(),  
    panel.grid.minor = element_blank(), 
    panel.border = element_rect(color = "black", fill = NA, linewidth = 2),
    axis.line = element_blank(),
    axis.text = element_text(size = fontsize, color = "black", face = "bold"), #family = police,
    axis.title = element_text(size = fontsize, color = "black", face = "bold"), #family = police,
    axis.ticks = element_line(color = "black", linewidth = 1)  # Adds tick marks
  ) +
  
  # Add the 2x2 table to the plot, placed in the top-left corner
  annotation_custom(
    grob = tableGrob(
      table_data, 
      theme = table_theme, 
      rows = NULL, # No rownames
      cols = NULL  # No colnames
    ), 
    xmin = -Inf, xmax = -1, ymin = Inf, ymax = 4
  )
)
dev.off()



# PCA----

# Variance Stabilizing Transformation
vsd <- vst(dds,blind=FALSE)
vsd$SampleID <- factor(as.character(vsd$Genotype_treatment), levels=c("WT_Ctl", "KO_Ctl", "WT_SB", "KO_SB"))

# Regularized Log Transformation
rld<-rlog(dds,blind=FALSE)
rld$SampleID <- factor(as.character(rld$Genotype_treatment), levels=c("WT_Ctl", "KO_Ctl", "WT_SB", "KO_SB"))

# Use rld instead of vsd because relatively small dataset = less than 10–15 samples per condition

# Personalization of PCA plot
cbPalette <- c( "darkgrey", "black","pink", "red")

# Extract PCA coordinates
pca_data <- as.data.frame(plotPCA(rld, intgroup = c("Genotype_treatment"), returnData = TRUE)) %>%
  arrange(Genotype_treatment) %>%
  mutate(Genotype_treatment = str_replace_all(Genotype_treatment, "_", " "),  # Replace underscores with spaces
         Genotype_treatment = str_replace_all(Genotype_treatment, "U", "Ctl"))

percentVar <- round(100 * attr(pca_data, "percentVar"))  # Round to whole numbers

# Define custom offsets for each point
x_nudge <- 15
pca_data$nudge_x <- c(x_nudge, x_nudge, x_nudge,
                      -x_nudge,-x_nudge,-x_nudge,
                      x_nudge,x_nudge,x_nudge,
                      -x_nudge,-x_nudge,-x_nudge)
y_nudge <- 0
pca_data$nudge_y <- c(y_nudge+1, y_nudge , y_nudge,# WT Ctl
                      y_nudge,y_nudge-1,y_nudge, # KO Ctl
                      y_nudge-1,y_nudge-1,y_nudge+2, # WT SB
                      y_nudge,y_nudge+1,y_nudge+2) # KO SB

# Create PCA plot with variance in axis labels
pdf("Plots/PCA_plot_review.pdf", width = 6, height = 4)#, family = police)
print(
  ggplot(pca_data, aes(x = PC1, y = PC2)) +
  
  # Plot PCA points with specified aesthetics
  geom_point(aes(fill = Genotype_treatment), 
             size = 8, shape = 21, color = "black", stroke = 1) +
  # Set the colors for Genotype_treatment
  scale_fill_manual(values = cbPalette) +
  # Add text labels for Genotype_treatment with adjusted position
  # geom_text(aes(label = Genotype_treatment, color = Genotype_treatment,
  #               x = PC1 + nudge_x, y = PC2 + nudge_y), 
  #           size = 8, fontface = "bold") +
  # Manually set the text color for labels to match the fill color
  scale_color_manual(values = cbPalette) +
  # Add dashed lines at x = 0 and y = 0 for visual reference
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  # Add % variance explained to x and y axis labels
  labs(x = paste0("PC1 (", percentVar[1], "%)"), 
       y = paste0("PC2 (", percentVar[2], "%)")) +
  # Apply a clean theme for the plot
  theme_bw() +
  theme(
    # Remove the legend
    legend.position = "none",  
    # Remove grid lines (both major and minor)
    panel.grid.major = element_blank(),  
    panel.grid.minor = element_blank(), 
    panel.border = element_rect(color = "black", fill = NA, linewidth = 2),#                           Add a thick border around the plot
    axis.line = element_blank(), #                                                                Remove the default x and y axis lines
    axis.text = element_text(size = fontsize, color = "black",  face = "bold"),#  Format axis text, #family = police,
    axis.title = element_text(size = fontsize, color = "black", face = "bold")#  Format axis titles, family = police,
  )
)
dev.off()

# Heatmap Top 35 ----

# Use rld data to plot heatmap of most variable genes
topVarGenes <- head(order(rowVars(assay(rld)), decreasing = TRUE), 35)
topVarData <- as.data.frame(assay(rld)[topVarGenes, ])

annotation_col <- metadata[, c('Genotype_treatment'), drop = FALSE]
rownames(annotation_col) <- rownames(metadata)

annotation_colors <- list(
  Genotype_treatment = c(
    "WT_Ctl" = "pink",
    "WT_SB" = "red",
    "KO_Ctl" = "grey",
    "KO_SB" = "black"
  )
)


pdf("Plots/heatmap_top35_genes.pdf", width = 8, height = 6)#, family = police)
print(
  pheatmap(
  topVarData,
  scale = "row",             # Scale data by rows (z-score normalization)
  cluster_rows = FALSE,       # Cluster genes by similarity
  cluster_cols = TRUE,       # Cluster samples by similarity
  show_rownames = TRUE,      # Show gene names
  show_colnames = TRUE,      # Show sample names
  color = colorRampPalette(rev(brewer.pal(9, "RdBu")))(100),  # Color palette
  main = "Top 35 Variable Genes",
  annotation_col = annotation_col,  # Add column annotations
  annotation_colors = annotation_colors  # Apply custom colors to annotations
)
)
dev.off()

# Heatmap of AB genes list----

# Use list of genes from AB and plot heatmap
AB_genes_list <- read_excel("list_of_genes.xlsx")
AB_genes_list <- AB_genes_list$both
AB_genes_list <- unique(AB_genes_list)
AB_genes_list <- unlist(strsplit(AB_genes_list, "/"))

all_genes_capitalized <- toupper(rownames(assay(rld)))
AB_genes_capitalized <- toupper(AB_genes_list)

my_genes <- c()
for(i in AB_genes_capitalized) {
  if (i %in% all_genes_capitalized) {
    my_genes <- c(my_genes, i)
  }
}

matching_indices <- which(rownames(assay(rld)) %in% my_genes)
data_for_heatmap <- as.data.frame(assay(rld)[matching_indices, ])
data_for_heatmap <- data_for_heatmap[rowSums(data_for_heatmap != 0) > 0, ]

columns.order <- c("A01", "A05", "A09", "A02", "A06", "A10", "A03", "A07", "A11", "A04", "A08", "A12")
data_for_heatmap <- data_for_heatmap[, columns.order]

# Plot heatmap
pdf("Plots/heatmap_AB_genes.pdf", width = 8, height = 6)
print(
  pheatmap(
    data_for_heatmap,
    scale = "row",          # Scale data by rows (z-score normalization)
    cluster_rows = F,       # Cluster genes by similarity
    cluster_cols = F,       # Cluster samples by similarity
    show_rownames = T,      # Show gene names
    show_colnames = T,      # Show sample names
    color = colorRampPalette(rev(brewer.pal(9, "RdBu")))(100),  # Color palette
    main = "Genes from AB list",
    annotation_col = annotation_col,  # Add column annotations
    annotation_colors = annotation_colors,  # Apply custom colors to annotations
    fontsize_row = 5
  )
)
dev.off()

# Heatmap Heme Biosynthesis ----

# Plot heatmap for Heme biosynthesis + hbb genes
heme.biosynthesis.genes <- read_excel("list_of_genes.xlsx", sheet = "Heme biosynthesis + hbb", col_names = 'gene')
heme.biosynthesis.genes$gene <- toupper(heme.biosynthesis.genes$gene) # capitalized
data.heme.biosynthesis <- as.data.frame(assay(rld)[ which(rownames(assay(rld)) %in% heme.biosynthesis.genes$gene), ])

columns.order <- c("A01", "A05", "A09", "A02", "A06", "A10", "A03", "A07", "A11", "A04", "A08", "A12")
data.heme.biosynthesis <- data.heme.biosynthesis[, columns.order]

rows.order <- c("SLC48A1","HBA1","HBA2","HBB","HBD","HBE1","HBG1","HBG2","ALAS1","ALAS2","ALAD","HMBS","UROD","UROS","CPOX","PPOX","FECH")
data.heme.biosynthesis <- data.heme.biosynthesis[rows.order,]
rownames(data.heme.biosynthesis)[rownames(data.heme.biosynthesis) == "SLC48A1"] <- "HRG1"


pdf("Plots/heatmap_heme_biosynthesis_hbb.pdf", width = 8, height = 6)
print(
  pheatmap(
    data.heme.biosynthesis,
    scale = "row",          # Scale data by rows (z-score normalization)
    cluster_rows = F,       # Cluster genes by similarity
    cluster_cols = F,       # Cluster samples by similarity
    show_rownames = T,      # Show gene names
    show_colnames = T,      # Show sample names
    color = colorRampPalette(rev(brewer.pal(9, "RdBu")))(100),  # Color palette
    main = "Genes from AB list",
    annotation_col = annotation_col,  # Add column annotations
    annotation_names_col = F,
    annotation_colors = annotation_colors,  # Apply custom colors to annotations
    fontsize_row = 10)
)
dev.off()


# Check if all the genes from rownames(assay(rld))[matching_indices] exist in res_SB
matching_genes <- rownames(assay(rld))[matching_indices]
valid_genes <- matching_genes[matching_genes %in% rownames(res_SB)]
invalid_genes <- matching_genes[!matching_genes %in% rownames(res_SB)]
print(invalid_genes)

res_SB_filtered <- res_SB[matching_genes, ]

pCutoff <- 0.05

# Apply the conditions for coloring
keyvals.colour <- ifelse(
  res_SB_filtered$log2FoldChange < -1 & res_SB_filtered$padj < pCutoff, 'red',
  ifelse(res_SB_filtered$log2FoldChange > 1 & res_SB_filtered$padj < pCutoff, 'red',
         ifelse(abs(res_SB_filtered$log2FoldChange) <= 1 & res_SB_filtered$padj < pCutoff, 'blue', 'grey'))
)
keyvals.colour[is.na(keyvals.colour)] <- 'grey'
names(keyvals.colour)[keyvals.colour == 'red'] <- expression("Pval & LogFC")
names(keyvals.colour)[keyvals.colour == 'blue'] <- 'Pval'
names(keyvals.colour)[keyvals.colour == 'grey'] <- 'NS'


# # Filter genes for text labels based on p-value and log2FoldChange conditions
# label_data <- res_SB_filtered[which((res_SB_filtered$padj < pCutoff) & 
#                                       (abs(res_SB_filtered$log2FoldChange) > 1)), ]


# Volcano plot
pdf("Plots/volcano_plot_AB_genes.pdf", width = 8, height = 6)#, family = police)
print(
EnhancedVolcano(res_SB_filtered,
                lab = rownames(res_SB_filtered),
                x = 'log2FoldChange',
                y = 'padj',
                pCutoff = pCutoff,
                FCcutoff = 1,
                title = "",      # "KO_SB vs WT_SB",
                subtitle = "",   #,Ref. WT_SB",
                caption ="",     # "Differential Expression",
                colAlpha = 0.7,
                legendPosition = 'top',
                labSize = 10,
                labFace = "italic",
                labCol = 'black',
                drawConnectors = T,
                widthConnectors = 0.1,
                colConnectors = 'red',
                boxedLabels = TRUE,
                encircle = FALSE,
                ylim = c(0, 18),  # Adjust y-axis max (change 10 to desired value)
                xlim = c(-3.5, 4.5),
                gridlines.major = TRUE,
                gridlines.minor = FALSE,
                border = 'full',
                borderWidth = 1.0,
                borderColour = 'black',
                xlab = "",
                ylab = "",
                pointSize = 5,
                colCustom = keyvals.colour
) +
  theme_classic() +
  theme(
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1.5),  # Add border
    axis.title.x = element_text(size = fontsize, face = "bold", colour = 'black'),  # Bold x-axis title
    axis.title.y = element_text(size = fontsize, face = "bold", colour = 'black'),  # Bold y-axis title
    axis.text.x = element_text(size = fontsize, face = "bold", colour = 'black'),   # Bold x-axis tick labels
    axis.text.y = element_text(size = fontsize, face = "bold", colour = 'black'),   # Bold y-axis tick labels
    legend.text = element_text(size = fontsize, colour = 'black')
  ) +
  guides(color = guide_legend(title = NULL))  # Remove legend title
)
dev.off()



# HBE1 plot ----

data_for_plot_HBE1 <- as.data.frame(assay(rld))
data_for_plot_HBE1 <- data_for_plot_HBE1["HBE1", , drop = FALSE]
data_for_plot_HBE1 <- as.data.frame(t(data_for_plot_HBE1))
SampleName <- rownames(data_for_plot_HBE1)
data_for_plot_HBE1$SampleName <- SampleName
data_for_plot_HBE1 <- merge(data_for_plot_HBE1, metadata)

# Define custom colors
genotype_colors <- c(
  "KO_Ctl" = "pink",
  "KO_SB" = "red",
  "WT_Ctl" = "grey",
  "WT_SB" = "black"
)


pdf("Plots/HBE1.pdf", width = 8, height = 6)
print(
  ggplot(data_for_plot_HBE1, aes(x = Genotype_treatment, y = HBE1, color = Genotype_treatment)) +
    geom_point(size = 4)+
    scale_color_manual(values = genotype_colors) +
    theme_minimal() +
    labs(x = "Group", y = "HBE1 (rld)", color = "Group") +  # Rename legend
    theme(axis.text = element_text(size = 14, face = "bold"),
          axis.title = element_text(size = 16, face = "bold"),
          legend.title = element_text(size = 14, face = "bold"),
          legend.text = element_text(size = 12))
)
dev.off()


# MitoCarta ----

MitoCartaHuman <- read_xls("Human.MitoCarta3.0.xls", sheet = "A Human MitoCarta3.0")
MitoCartaHumanGenes <- MitoCartaHuman$Symbol

mito.data <- list()
mito.data[["KO_Ctl vs WT_Ctl"]] <- res_U.mito
mito.data[["KO_SB vs WT_SB"]] <- res_SB.mito
mito.data[["KO_SB vs KO_Ctl"]] <- res_koSB_koU.mito
mito.data[["WT_SB vs WT_Ctl"]] <- res_wtSB_wtU.mito

heme.biosynthesis.genes.list <- heme.biosynthesis.genes$gene[1:9]

for (i in names(mito.data)) {
  parts <- unlist(strsplit(i, " vs "))
  comparison <- parts[1]
  ref <- parts[2]
  
  df <- mito.data[[i]] %>%
    mutate(
      Significant = case_when(
        padj < pCutoff & abs(log2FoldChange) > 1 ~ "Significant",
        TRUE ~ "Not Significant"
      ),
      HemeBioLabel = ifelse(
        rownames(.) %in% heme.biosynthesis.genes.list & Significant == "Significant",
        rownames(.),
        NA_character_
      ),
      Label = ifelse(
        rownames(.) %in% MitoCartaHumanGenes & Significant == "Significant",
        rownames(.),
        NA_character_
    )
    )
  
  pdf(file.path("Plots", "MitoCarta", paste0("MitoCarta_volcano_plot_", comparison, "_vs_", ref, "_AllLabels.pdf")), width = 8, height = 6)
  print(
    ggplot(df, aes(x = log2FoldChange, y = -log10(padj), color = Significant)) +
      geom_point(alpha = 0.7) +
      geom_text_repel(
        aes(label = Label),
        color = "black",       # labels in black
        size = 3,
        box.padding = 0.5,
        max.overlaps = 20
      )+
      scale_color_manual(values = c("Significant" = "red", "Not Significant" = "grey")) +
      geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black") +
      geom_hline(yintercept = -log10(pCutoff), linetype = "dashed", color = "black") +
      labs(
        title = str_replace_all(i, "_", " "),
        subtitle = "",
        caption = "",
        x = expression(Log[2]*"FC"),
        y = "-log10(p.adj)"
      ) +
      theme_minimal() +
      theme(
        legend.position = "top",
        legend.title = element_blank(),
        legend.text = element_text(size = 12),
        plot.title = element_text(size = 12),
        plot.subtitle = element_text(size = 12),
        plot.caption = element_text(size = 12),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 12),
        panel.border = element_rect(color = "black", fill = NA, size = 1)
      )
  )
  dev.off()
  
  pdf(file.path("Plots", "MitoCarta", paste0("MitoCarta_volcano_plot_", comparison, "_vs_", ref, ".pdf")), width = 8, height = 6)
  print(
    ggplot(df, aes(x = log2FoldChange, y = -log10(padj), color = Significant)) +
      geom_point(alpha = 0.7) +
      geom_text_repel(
        aes(label = HemeBioLabel),
        color = "black",       # labels in black
        size = 3,
        box.padding = 0.5,
        max.overlaps = 20
      )+
    scale_color_manual(values = c("Significant" = "red", "Not Significant" = "grey")) +
      geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black") +
      geom_hline(yintercept = -log10(pCutoff), linetype = "dashed", color = "black") +
      labs(
        title = str_replace_all(i, "_", " "),
        subtitle = "",
        caption = "",
        x = expression(Log[2]*"FC"),
        y = "-log10(p.adj)"
      ) +
      theme_minimal() +
      theme(
        legend.position = "top",
        legend.title = element_blank(),
        legend.text = element_text(size = 12),
        plot.title = element_text(size = 12),
        plot.subtitle = element_text(size = 12),
        plot.caption = element_text(size = 12),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 12),
        panel.border = element_rect(color = "black", fill = NA, size = 1)
      )
  )
  dev.off()

}









# Create subdataframe for each corner of the future plots
data.to.plot.mito <- data.to.plot[data.to.plot$Gene %in% MitoCartaHumanGenes, ]

mito_data_up_up <- subset(data.to.plot.mito,log2FoldChange_KO_SB_vs_KO_Ctl > up_down_treshlod & log2FoldChange_WT_SB_vs_WT_Ctl > up_down_treshlod)
mito_data_up_down <- subset(data.to.plot.mito,log2FoldChange_KO_SB_vs_KO_Ctl > up_down_treshlod & log2FoldChange_WT_SB_vs_WT_Ctl < -up_down_treshlod)
mito_data_down_up <- subset(data.to.plot.mito,log2FoldChange_KO_SB_vs_KO_Ctl < -up_down_treshlod & log2FoldChange_WT_SB_vs_WT_Ctl > up_down_treshlod)
mito_data_down_down <- subset(data.to.plot.mito,log2FoldChange_KO_SB_vs_KO_Ctl < -up_down_treshlod & log2FoldChange_WT_SB_vs_WT_Ctl < -up_down_treshlod)
mito_data_middle <- data.to.plot.mito[!(data.to.plot.mito$Gene %in% c(data_up_up$Gene, data_up_down$Gene, data_down_up$Gene, data_down_down$Gene)), ]

# Create a matrix with the 4 values (counts) in a 2x2 table
mito_table_data <- matrix(
  c(nrow(mito_data_down_up), nrow(mito_data_down_down), nrow(mito_data_up_up), nrow(mito_data_up_down)),
  nrow = 2, ncol = 2,
  dimnames = list(NULL, NULL)  # No row and column names
)


pdf("Plots/MitoCarta/MitoCarta_WT_SB_WT_Ctl_vs_KO_SB_KO_Ctl.pdf", width = 8, height = 6)#, family = police)
# Create the plot
ggplot(mito_data_middle, aes(x = log2FoldChange_WT_SB_vs_WT_Ctl, y = log2FoldChange_KO_SB_vs_KO_Ctl)) +
  geom_point(color="black", alpha=0.5, size = 5) +
  geom_point(data = mito_data_up_up,aes(x = log2FoldChange_WT_SB_vs_WT_Ctl, y = log2FoldChange_KO_SB_vs_KO_Ctl),color = 'turquoise', alpha=0.5, size=5) +
  geom_text_repel(
    data = mito_data_up_up,
    aes(label = Gene),      # adjust 'gene' to your actual column name
    size = 3.5,
    color = "black",
    box.padding = 0.4,
    point.padding = 0.2,
    max.overlaps = 20)+
    
  geom_point(data = mito_data_down_down,
             aes(x = log2FoldChange_WT_SB_vs_WT_Ctl, y = log2FoldChange_KO_SB_vs_KO_Ctl), 
             color = 'turquoise', alpha=0.5, size=5) +
  geom_text_repel(
    data = mito_data_down_down,
    aes(label = Gene),      # adjust 'gene' to your actual column name
    size = 3.5,
    color = "black",
    box.padding = 0.4,
    point.padding = 0.2,
    max.overlaps = 20)+
  geom_point(data = mito_data_down_up,
             aes(x = log2FoldChange_WT_SB_vs_WT_Ctl, y = log2FoldChange_KO_SB_vs_KO_Ctl), 
             color = 'orangered', alpha=0.5, size=5) +
  geom_point(data = mito_data_up_down, 
             aes(x = log2FoldChange_WT_SB_vs_WT_Ctl, y = log2FoldChange_KO_SB_vs_KO_Ctl), 
             color = 'orangered', alpha=0.5, size=5) +
  
  geom_hline(yintercept = up_down_treshlod, linetype = "dashed", color = "black") +
  geom_hline(yintercept = -up_down_treshlod, linetype = "dashed", color = "black") +
  geom_vline(xintercept = -up_down_treshlod, linetype = "dashed", color = "black") +
  geom_vline(xintercept = up_down_treshlod, linetype = "dashed", color = "black") +
  
  
  
  labs(x = bquote(Log[2] * "FC (KO SB vs KO Ctl)"), 
       y = bquote(Log[2] * "WT SB vs WT Ctl")) +
  
  scale_x_continuous() +
  scale_y_continuous() +
  
  theme_minimal() +
  theme(
    legend.position = "none",  
    panel.grid.major = element_blank(),  
    panel.grid.minor = element_blank(), 
    panel.border = element_rect(color = "black", fill = NA, linewidth = 2),
    axis.line = element_blank(),
    axis.text = element_text(size = fontsize, color = "black", face = "bold"), #family = police,
    axis.title = element_text(size = fontsize, color = "black", face = "bold"), #family = police,
    axis.ticks = element_line(color = "black", linewidth = 1)  # Adds tick marks
  ) +
  
  # Add the 2x2 table to the plot, placed in the top-left corner
  annotation_custom(
    grob = tableGrob(
      mito_table_data, 
      theme = table_theme, 
      rows = NULL, # No rownames
      cols = NULL  # No colnames
    ), 
    xmin = -Inf, xmax = -1, ymin = Inf, ymax = 4
  )

dev.off()



# Individual Volcano Plot - create function

selected.ref <- "WT_Ctl"
to.compare <- "KO_Ctl"

run_DESeq2_comparison <- function(count_matrix, metadata, selected.ref, to.compare, output_dir = "Plots",
                                  font_family = "Helvetica", pCutoff = 0.05, FCcutoff = 1) {
  # Load required libraries
  require(DESeq2)
  require(ashr)
  require(EnhancedVolcano)
  
  # Ensure output directory exists
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  
  # Set reference level
  metadata$Genotype_treatment <- relevel(factor(metadata$Genotype_treatment), ref = selected.ref)
  
  # Construct DESeqDataSet
  dds <- DESeqDataSetFromMatrix(countData = count_matrix,
                                colData = metadata,
                                design = ~ Genotype_treatment)
  
  # Filter low counts
  dds <- dds[rowSums(counts(dds)) > 1, ]
  
  # Run DESeq
  dds <- DESeq(dds)
  
  # Get results
  res <- results(dds,
                 contrast = c("Genotype_treatment", to.compare, selected.ref),
                 pAdjustMethod = "BH")
  
  # Perform shrinkage
  # Identify coefficient name dynamically
  coef_name <- grep(paste0("Genotype_treatment_", to.compare, "_vs_", selected.ref),
                    resultsNames(dds), value = TRUE)
  if (length(coef_name) == 0) {
    warning("No coefficient found for this contrast; skipping shrinkage.")
  } else {
    res <- lfcShrink(dds, coef = coef_name, type = "ashr")
  }
  
  # # Save MA plot
  # pdf(file.path(output_dir, paste0("Mito_MAplot_", to.compare, "_vs_", selected.ref, ".pdf")),
  #     width = 8, height = 6)
  # par(family = font_family, font = 2)
  # plotMA(res, main = paste(to.compare, "vs", selected.ref), ylim = c(-5, 5))
  # dev.off()
  
  # # Save Volcano plot
  # pdf(file.path(output_dir, paste0("Mito_VolcanoPlot_", to.compare, "_vs_", selected.ref, ".pdf")),width = 8, height = 6)
  # EnhancedVolcano(res,
  #                 lab = rownames(res),
  #                 x = 'log2FoldChange',
  #                 y = 'padj',
  #                 pCutoff = pCutoff,
  #                 FCcutoff = FCcutoff,
  #                 title = paste(to.compare, "vs", selected.ref),
  #                 subtitle = paste("Ref.", selected.ref),
  #                 caption = "Differential Expression",
  #                 colAlpha = 0.7,
  #                 legendPosition = 'top',
  #                 legendLabSize = 10)
  # dev.off()
  
  # Return DESeq2 results
  return(res)
}


res_KO_Ctl_vs_WT_Ctl <- run_DESeq2_comparison(
  count_matrix = counts,
  metadata = metadata,
  selected.ref = "WT_Ctl",
  to.compare = "KO_Ctl"
)

df <- as.data.frame(res_KO_Ctl_vs_WT_Ctl)

# Keep only genes present in your list
df_subset <- df[rownames(df) %in% MitoCartaHumanGenes, ]

EnhancedVolcano(df_subset,
                lab = rownames(df_subset),
                x = 'log2FoldChange',
                y = 'padj',
                pCutoff = 0.05,
                FCcutoff = 1,
                title = paste(to.compare, "vs", selected.ref),subtitle = paste("Ref.", selected.ref),
                caption = "Differential Expression",
                colAlpha = 0.7,
                legendPosition = 'top',
                legendLabSize = 10)


upregulted.genes <- df_subset[(df_subset$log2FoldChange > 1) & (df_subset$padj < 0.05), ] %>% rownames_to_column("Gene")
downregulted.genes <- df_subset[(df_subset$log2FoldChange < -1) & (df_subset$padj < 0.05), ] %>% rownames_to_column("Gene")

BiocManager::install(c("clusterProfiler", "org.Hs.eg.db", "org.Mm.eg.db", "enrichplot"))


library(clusterProfiler)
library(enrichplot)
library(org.Hs.eg.db)  # Use org.Mm.eg.db for mouse
library(forcats)

up_genes <- bitr(upregulted.genes$Gene, 
                 fromType = "SYMBOL",
                 toType = "ENTREZID",
                 OrgDb = org.Hs.eg.db)

down_genes <- bitr(downregulted.genes$Gene, 
                   fromType = "SYMBOL",
                   toType = "ENTREZID",
                   OrgDb = org.Hs.eg.db)

ego_up <- enrichGO(gene          = up_genes$ENTREZID,
                   OrgDb         = org.Hs.eg.db,
                   keyType       = "ENTREZID",
                   ont           = "BP",     # Biological Process
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.05,
                   qvalueCutoff  = 0.05)

### ADD NODE_SIZE = 10 ?

ego_down <- enrichGO(gene          = down_genes$ENTREZID,
                     OrgDb         = org.Hs.eg.db,
                     keyType       = "ENTREZID",
                     ont           = "BP",
                     pAdjustMethod = "BH",
                     pvalueCutoff  = 0.05,
                     qvalueCutoff  = 0.05)


# Dotplot for enriched GO terms
dotplot(ego_up, showCategory = 20, title = "Upregulated genes - GO BP")
dotplot(ego_down, showCategory = 20, title = "Downregulated genes - GO BP")




ego_down_df <- ego_down@result

# Extract numerator and denominator from GeneRatio and BgRatio
ego_down_df <- ego_down_df %>%
  mutate(
    GeneRatio_num = as.numeric(sub("/.*", "", GeneRatio)),
    GeneRatio_den = as.numeric(sub(".*/", "", GeneRatio)),
    BgRatio_num   = as.numeric(sub("/.*", "", BgRatio)),
    BgRatio_den   = as.numeric(sub(".*/", "", BgRatio))
  )

# Calculate contingency table components
ego_down_df <- ego_down_df %>%
  mutate(
    a = GeneRatio_num,
    b = GeneRatio_den - GeneRatio_num,
    c = BgRatio_num - a,
    d = BgRatio_den - BgRatio_num - b,
    OddsRatio = (a * d) / (b * c)
  )


ego_down_df <- ego_down_df %>%
  arrange(pvalue) %>%
  head(20) %>%
  mutate(Description = fct_rev(fct_inorder(Description)))

library(ggtext)

ggplot(ego_down_df, 
       aes(x = -log10(p.adjust), 
           y = Description,
           size = OddsRatio, 
           color = zScore)) +
  geom_point(alpha = 0.9) +
  scale_size_continuous(range = c(3, 13)) +
  scale_color_gradient(low = "darkgreen", high = "limegreen") +
  labs(
    title = "<span style='color:darkgreen;'>Down</span>regulated genes<br>in KO SB vs WT SB",
    x = expression(-log[10](p~value)),
    y = NULL,
    size = "Odds Ratio",
    color = "z-score"
  ) +
  theme_bw(base_size = 14) +
  theme(
    plot.title = element_markdown(size = 18, face = "bold", hjust = 0.5),
    axis.text.y = element_text(size = 10, color = "black"),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank()
  )

