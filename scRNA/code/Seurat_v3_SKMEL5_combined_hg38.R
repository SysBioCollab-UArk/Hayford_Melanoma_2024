library(Seurat)
library(svglite)
library(ggplot2)
library(dplyr)
library(purrr)
library(ggplot2)
library(cowplot)
library(clusterProfiler)
library(org.Hs.eg.db)
library(reshape2)
library(Hmisc)
library(ggpubr)
library(transport)

# Pull this from the data I shared on Box (data_large)
untreated <- Read10X(data.dir = "../data/Untreated/filtered_feature_bc_matrix/")
idling <- Read10X(data.dir = "../data/Idling/filtered_feature_bc_matrix/")

# Initialize the Seurat object with the raw (non-normalized data).
untreated <- CreateSeuratObject(counts = untreated, project = "Untreated", min.cells = 3, min.features = 200)
idling <- CreateSeuratObject(counts = idling, project = "Idling", min.cells = 3, min.features = 200)

# Add lineage information
## Untreated
cellMeta_untreated <- untreated@meta.data
barcode_cell_df_UT <- read.csv("../data/Untreated_LineageBC_cellBC.csv")
barcode_cell_df_UT <- as.data.frame(t(barcode_cell_df_UT))
colnames(barcode_cell_df_UT) <- c("Barcode", "Cell Barcode")
barcode_cell_UT <- barcode_cell_df_UT[-1,]
barcode_cell_UT$`Cell Barcode` <- paste(barcode_cell_UT$`Cell Barcode`, "-1", sep="")
cellMeta_untreated$lineage <- as.character(barcode_cell_UT$Barcode[match(rownames(cellMeta_untreated),
                                                                         barcode_cell_UT$`Cell Barcode`)])
cellMeta_untreated$lineage[nchar(as.character(cellMeta_untreated$lineage)) != 20] <- "XXX"
lineageMeta_UT <- subset(cellMeta_untreated, select = "lineage")
untreated <- AddMetaData(untreated, lineageMeta_UT, col.name = "lineage")

## Idling
cellMeta_idling <- idling@meta.data
barcode_cell_df_I <- as.data.frame(t(read.csv("~/Documents/QuarantaLab/Treated_LineageBC_cellBC.csv")))
colnames(barcode_cell_df_I) <- c("Barcode", "Cell Barcode")
barcode_cell_I <- barcode_cell_df_I[-1,]
barcode_cell_I$`Cell Barcode` <- paste(barcode_cell_I$`Cell Barcode`, "-1", sep="")
cellMeta_idling$lineage <- as.character(barcode_cell_I$Barcode[match(rownames(cellMeta_idling),
                                                                     barcode_cell_I$`Cell Barcode`)])
cellMeta_idling$lineage[nchar(as.character(cellMeta_idling$lineage)) != 20] <- "XXX"
lineageMeta_I <- subset(cellMeta_idling, select = "lineage")
idling <- AddMetaData(idling, lineageMeta_I, col.name = "lineage")

combined <- merge(untreated, y = idling, add.cell.ids = c("Untreated", "Idling"), project = "Combined")

combined[["percent.mt"]] <- PercentageFeatureSet(combined, pattern = "^MT-")

combined <- SCTransform(combined, assay = 'RNA',
                        new.assay.name = 'SCT',
                        vars.to.regress = c('percent.mt', 'nFeature_RNA', 'nCount_RNA'))

VlnPlot(combined, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# combined <- subset(combined, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
combined <- subset(combined, subset = lineage != "NA")

combined <- NormalizeData(combined)

combined <- CellCycleScoring(object = combined, s.features = cc.genes.updated.2019$s.genes,
                             g2m.features = cc.genes.updated.2019$g2m.genes, assay = 'SCT',
                             set.ident = TRUE)

combined <- SCTransform(combined, assay = 'RNA', new.assay.name = 'SCT',
  vars.to.regress = c('percent.mt', 'nFeature_RNA', 'nCount_RNA', 'S.Score', 'G2M.Score'))

combined <- FindVariableFeatures(combined, assay = 'SCT', selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
# top10 <- head(VariableFeatures(combined), 10)
# 
# # plot variable features with and without labels
# plot1 <- VariableFeaturePlot(combined)
# plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
# plot1 + plot2

all.genes <- rownames(combined)
combined <- ScaleData(combined, assay = 'SCT', features = all.genes)

combined <- RunPCA(combined, assay = 'SCT', features = VariableFeatures(object = combined))
VizDimLoadings(combined, dims = 1:2, reduction = "pca")
DimPlot(combined, reduction = "pca")

combined <- FindNeighbors(combined, assay = 'SCT', dims = 1:10)
combined <- FindClusters(combined, assay = 'SCT', resolution = 0.5)

combined <- RunUMAP(combined, assay = 'SCT', dims = 1:10)

DimPlot(combined, reduction = "umap", group.by = "orig.ident") +
  theme_bw() + scale_color_manual(values = c("blue", "red")) +
  theme(axis.text=element_text(size=14), legend.text = element_text(size=14),
        axis.title=element_text(size=14), legend.position = "none",
        panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  xlim(-11.5, 11.5) + ylim(-11.5, 11.5) +
  ggsave("UMAP_combined_SKMEL5_hg38_qcCCReg.svg", width = 4, height = 3)
  
# DimPlot(combined, reduction = "umap", group.by = "Phase") +
#   theme_bw() + #scale_color_manual(values = c("blue", "red")) +
#   theme(axis.text = element_text(size = 14), legend.position = "right",
#         panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
#   ggsave("UMAP_combined_SKMEL5_hg38_qcCCReg_CCPhase_leg.pdf", width = 5, height = 3)
# 
# DimPlot(combined, reduction = "umap", group.by = "State") +
#   theme_bw() + scale_color_manual(values = c("green3", "gold")) +
#   theme(axis.text = element_text(size = 14), legend.position = "none",
#         panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
#   xlim(-11.5, 11.5) + ylim(-11.5, 11.5) +
#   ggsave("UMAP_combined_SKMEL5_hg38_qcCCReg_CCState_leg.svg", width = 4, height = 3)


# DimPlot(combined, reduction = "umap", group.by = "seurat_clusters") +
#   theme_bw() + #scale_color_manual(values = c("blue", "red")) +
#   theme(axis.text = element_text(size = 14), legend.position = "right",
#         panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
#   ggsave("UMAP_combined_SKMEL5_hg38_qcCCReg_clusters_leg.pdf", width = 5, height = 3)


######

# Determine metrics to plot present in seurat_control@meta.data
metrics <-  c("nCount_RNA", "nFeature_RNA", "S.Score", "G2M.Score", "percent.mt")

# Extract the UMAP coordinates for each cell and include information about the metrics to plot
qc_data <- FetchData(combined, 
                     vars = c(metrics, "ident", "UMAP_1", "UMAP_2"))

# Adding cluster label to center of cluster on UMAP
umap_label <- FetchData(combined, 
                        vars = c("ident", "UMAP_1", "UMAP_2"))  %>%
  as.data.frame() %>% 
  group_by(ident) %>%
  summarise(x=mean(UMAP_1), y=mean(UMAP_2))


# Plot a UMAP plot for each metric
map(metrics, function(qc){
  ggplot(qc_data,
         aes(UMAP_1, UMAP_2)) +
    theme_bw() +
    geom_point(aes_string(color=qc), 
               alpha = 0.7) +
    scale_color_gradient(guide = FALSE, 
                         low = "grey90", 
                         high = "blue")  +
    geom_text(data=umap_label, 
              aes(label=ident, x, y)) +
    ggtitle(qc)

}) %>%
  plot_grid(plotlist = .) +
  ggsave("UMAP_combined_SKMEL5_hg38_QCmetrics_qcCCReg.svg",
         width = 9, height = 6)

########

umap_vals <- as.data.frame(Embeddings(combined, reduction = "umap"))
umap_vals_df <- dplyr::bind_rows(umap_vals)
colnames(umap_vals_df) <- c("UMAP_1", "UMAP_2")
pops_all <- subset(combined@meta.data, select = "orig.ident")
umap_vals_df$Population <- unlist(pops_all)
all_Kmeans <- kmeans(umap_vals_df[,1:2,], centers = 3)
umap_vals_df$Kmeans <- all_Kmeans$cluster

ggplot(umap_vals_df, aes(x=UMAP_1, y=UMAP_2, color = factor(Kmeans))) +
  geom_point() + theme_bw() + 
  scale_color_manual(values = c("magenta", "gold", "skyblue"), name = "Clusters") +
  theme(axis.text = element_text(size = 14), legend.position = "right",
        panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  ggsave("combined_clustersSKMEL5.svg", width = 5, height = 3)

# save(combined, file = "SKMEL5_scRNAseq_hg38_combined.RData")
# load("SKMEL5_scRNAseq_hg38_combined.RData")

### Get proportion of CC cells in each condition
phase_UT <- FetchData(comb_UT, vars = "Phase")
phase_UT$Condition = "Untreated"
phase_I <- FetchData(comb_I, vars = "Phase")
phase_I$Condition = "Idling"
phase_all <- rbind(phase_UT, phase_I)
phase_all$Condition <- factor(phase_all$Condition, levels = c("Untreated", "Idling"))


# ggplot(phase_all, aes(x = Phase, group = Condition, color = Condition)) +
#   theme_bw() + geom_bar(aes(fill=as.factor(Condition), y = (..count..)), 
#                         position = "dodge", color = "black") +
#   scale_fill_manual(values = c("red", "blue")) +
#   labs(x = "Cell Cycle Phase", y = "Number of Cells") +
#   theme(axis.text.y = element_text(size = 14),
#         legend.position = "right", legend.text = element_text(size = 14),
#         plot.title = element_text(size = 16, hjust = 0.5, face = "bold"), 
#         axis.text=element_text(size=14),
#         legend.title = element_blank(), axis.title=element_text(size=14),
#         panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
#   ggsave("SKMEL5_cellCyclePhaseBreakdown_byCondition.pdf", width = 6, height = 4)

# Percentages
phase_UT_table <- as.data.frame.matrix(table(phase_UT))
phase_UT_table <- transform(phase_UT_table, percent = Untreated / sum(phase_UT_table$Untreated))
phase_UT_table$Phase <- rownames(phase_UT_table)
phase_UT_table$Condition <- "Untreated"
names(phase_UT_table) <- c("Count", "Percent", "Phase", "Condition")

phase_I_table <- as.data.frame.matrix(table(phase_I))
phase_I_table <- transform(phase_I_table, percent = Idling / sum(phase_I_table$Idling))
phase_I_table$Phase <- rownames(phase_I_table)
phase_I_table$Condition <- "Idling"
names(phase_I_table) <- c("Count", "Percent", "Phase", "Condition")

phase_all_table <- rbind(phase_UT_table, phase_I_table)
phase_all_table_melt <- melt(phase_all_table, id.vars = c("Phase", "Condition"),
                             measure.vars = "Percent")
phase_all_table_melt$Condition <- factor(phase_all_table_melt$Condition, 
                                         levels = c("Untreated", "Idling"))

# ggplot(phase_all_table_melt, aes(x = Phase, y = value,
#                                  group = Condition, fill = Condition)) +
#   theme_bw() + geom_bar(stat = "identity", position = "dodge", color = "black") +
#   scale_fill_manual(values = c("red", "blue")) +
#   scale_y_continuous(labels = scales::percent) +
#   labs(x = "Cell Cycle Phase", y = "Percentage of Cells") +
#   theme(axis.text.y = element_text(size = 14),
#         legend.position = "right", legend.text = element_text(size = 14),
#         plot.title = element_text(size = 16, hjust = 0.5, face = "bold"), 
#         axis.text=element_text(size=14),
#         legend.title = element_blank(), axis.title=element_text(size=14),
#         panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
#   ggsave("SKMEL5_cellCyclePhaseBreakdown_byCondition_percent.pdf", width = 6, height = 4)

# ggplot(phase_all_table_melt, aes(x = Condition, y = value,
#                                  group = Phase, fill = Phase)) +
#   theme_bw() + geom_bar(stat = "identity", position = "fill", color = "black") +
#   # scale_fill_manual(values = c("red", "blue")) +
#   scale_y_continuous(labels = scales::percent) +
#   labs(x = "Condition", y = "Percentage of Cells") +
#   theme(axis.text.y = element_text(size = 14),
#         legend.position = "right", legend.text = element_text(size = 14),
#         plot.title = element_text(size = 16, hjust = 0.5, face = "bold"), 
#         axis.text=element_text(size=14),
#         legend.title = element_blank(), axis.title=element_text(size=14),
#         panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
#   ggsave("SKMEL5_cellCyclePhaseBreakdown_byPhase_percent.pdf", width = 6, height = 4)

### Get proportion of CC cells in each state
state_UT <- FetchData(comb_UT, vars = "State")
state_UT$Condition = "Untreated"
state_I <- FetchData(comb_I, vars = "State")
state_I$Condition = "Idling"
state_all <- rbind(state_UT, state_I)
state_all$Condition <- factor(state_all$Condition, levels = c("Untreated", "Idling"))


# ggplot(state_all, aes(x = State, group = Condition, color = Condition)) +
#   theme_bw() + geom_bar(aes(fill=as.factor(Condition), y = (..count..)), 
#                         position = "dodge", color = "black") +
#   scale_fill_manual(values = c("red", "blue")) +
#   labs(x = "Cell Cycle State", y = "Number of Cells") +
#   theme(axis.text.y = element_text(size = 14),
#         legend.position = "right", legend.text = element_text(size = 14),
#         plot.title = element_text(size = 16, hjust = 0.5, face = "bold"), 
#         axis.text=element_text(size=14),
#         legend.title = element_blank(), axis.title=element_text(size=14),
#         panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
#   ggsave("SKMEL5_cellCycleStateBreakdown_byCondition.pdf", width = 6, height = 4)

# Percentages
state_UT_table <- as.data.frame.matrix(table(state_UT))
state_UT_table <- transform(state_UT_table, percent = Untreated / sum(state_UT_table$Untreated))
state_UT_table$State <- rownames(state_UT_table)
state_UT_table$Condition <- "Untreated"
names(state_UT_table) <- c("Count", "Percent", "State", "Condition")

state_I_table <- as.data.frame.matrix(table(state_I))
state_I_table <- transform(state_I_table, percent = Idling / sum(state_I_table$Idling))
state_I_table$State <- rownames(state_I_table)
state_I_table$Condition <- "Idling"
names(state_I_table) <- c("Count", "Percent", "State", "Condition")

state_all_table <- rbind(state_UT_table, state_I_table)
state_all_table_melt <- melt(state_all_table, id.vars = c("State", "Condition"),
                             measure.vars = "Percent")
state_all_table_melt$Condition <- factor(state_all_table_melt$Condition, 
                                         levels = c("Untreated", "Idling"))

# ggplot(state_all_table_melt, aes(x = State, y = value,
#                                  group = Condition, fill = Condition)) +
#   theme_bw() + geom_bar(stat = "identity", position = "dodge", color = "black") +
#   scale_fill_manual(values = c("red", "blue")) +
#   scale_y_continuous(labels = scales::percent) +
#   labs(x = "Cell Cycle State", y = "Percentage of Cells") +
#   theme(axis.text.y = element_text(size = 14),
#         legend.position = "right", legend.text = element_text(size = 14),
#         plot.title = element_text(size = 16, hjust = 0.5, face = "bold"), 
#         axis.text=element_text(size=14),
#         legend.title = element_blank(), axis.title=element_text(size=14),
#         panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
#   ggsave("SKMEL5_cellCycleStateBreakdown_byCondition_percent.pdf", width = 6, height = 4)

# ggplot(state_all_table_melt, aes(x = Condition, y = value,
#                                  group = State, fill = State)) +
#   theme_bw() + geom_bar(stat = "identity", position = "fill", color = "black") +
#   scale_fill_manual(values = c("green3", "gold")) +
#   scale_y_continuous(labels = scales::percent) +
#   labs(x = "Condition", y = "Percentage of Cells") +
#   theme(axis.text.y = element_text(size = 14),
#         legend.position = "none", legend.text = element_text(size = 14),
#         plot.title = element_text(size = 16, hjust = 0.5, face = "bold"), 
#         axis.text=element_text(size=14),
#         legend.title = element_blank(), axis.title=element_text(size=14),
#         panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
#   ggsave("SKMEL5_cellCycleStateBreakdown_byState_percent.pdf", width = 2.5, height = 4)


### do findmarkers (3 comparisons)
Idents(combined) <- "seurat_clusters"

# Idling - cycling vs non-cycling
I_cycling.markers <- FindMarkers(object = combined, ident.1 = 6,
                                  ident.2 = c(0,2,4), min.pct = 0.25)
I_noncycling.markers <- FindMarkers(object = combined, ident.1 = c(0,2,4),
                                 ident.2 = 6, min.pct = 0.25)

I_cycling_GO <- enrichGO(gene = rownames(subset(I_cycling.markers, avg_logFC > 0.5)),
                         universe = row.names(combined@assays$RNA@data),
                         OrgDb = org.Hs.eg.db,
                         keyType = 'SYMBOL',
                         ont = "BP",
                         pAdjustMethod = "BH",
                         pvalueCutoff  = 0.01,
                         qvalueCutoff  = 0.05)

dotplot(I_cycling_GO, showCategory=10) + ggsave("Idling_smallClusterGO_BP.pdf", width = 6, height = 4)

I_noncycling_GO <- enrichGO(gene = rownames(subset(I_noncycling.markers, avg_logFC > 0.5)),
                         universe = row.names(combined@assays$RNA@data),
                         OrgDb = org.Hs.eg.db,
                         keyType = 'SYMBOL',
                         ont = "CC",
                         pAdjustMethod = "BH",
                         pvalueCutoff  = 0.01,
                         qvalueCutoff  = 0.05)

dotplot(I_noncycling_GO, showCategory=20)

# Idling vs Untreated
I_markers <- FindMarkers(object = combined, ident.1 = c(0,2,4,6),
                                 ident.2 = c(1,3,5,7,8), min.pct = 0.25)
UT_markers <- FindMarkers(object = combined, ident.1 = c(1,3,5,7,8),
                                  ident.2 = c(0,2,4,6), min.pct = 0.25)

I_GO <- enrichGO(gene = rownames(subset(I_markers, avg_logFC > 0.5)),
                            universe = row.names(combined@assays$RNA@data),
                            OrgDb = org.Hs.eg.db,
                            keyType = 'SYMBOL',
                            ont = "BP",
                            pAdjustMethod = "BH",
                            pvalueCutoff  = 0.01,
                            qvalueCutoff  = 0.05)

dotplot(I_GO, showCategory=10) + ggsave("IdlingGO_BP.pdf", width = 8, height = 4)

UT_GO <- enrichGO(gene = rownames(subset(UT_markers, avg_logFC > 0.5)),
                 universe = row.names(combined@assays$RNA@data),
                 OrgDb = org.Hs.eg.db,
                 keyType = 'SYMBOL',
                 ont = "CC",
                 pAdjustMethod = "BH",
                 pvalueCutoff  = 0.01,
                 qvalueCutoff  = 0.05)

dotplot(UT_GO, showCategory=20)

# Untreated - outcast vs rest
UT_outcast.markers <- FindMarkers(object = combined, ident.1 = 7,
                         ident.2 = c(1,3,5,8), min.pct = 0.25)
UT_rest.markers <- FindMarkers(object = combined, ident.1 = c(1,3,5,8),
                          ident.2 = 7, min.pct = 0.25)

UT_outcast_GO <- enrichGO(gene = rownames(subset(UT_outcast.markers, avg_logFC > 0.5)),
                  universe = row.names(combined@assays$RNA@data),
                  OrgDb = org.Hs.eg.db,
                  keyType = 'SYMBOL',
                  ont = "BP",
                  pAdjustMethod = "BH",
                  pvalueCutoff  = 0.01,
                  qvalueCutoff  = 0.05)

dotplot(UT_outcast_GO, showCategory=10) + ggsave("Untreated_smallClusterGO_BP.pdf", width = 9, height = 4)

UT_rest_GO <- enrichGO(gene = rownames(subset(UT_rest.markers, avg_logFC > 0.5)),
                          universe = row.names(combined@assays$RNA@data),
                          OrgDb = org.Hs.eg.db,
                          keyType = 'SYMBOL',
                          ont = "CC",
                          pAdjustMethod = "BH",
                          pvalueCutoff  = 0.01,
                          qvalueCutoff  = 0.05)

dotplot(UT_rest_GO, showCategory=20)

# Untreated vs Idling - Small
I_small.markers <- FindMarkers(object = combined, ident.1 = 6,
                                  ident.2 = 7, min.pct = 0.25)
UT_small.markers <- FindMarkers(object = combined, ident.1 = 7,
                               ident.2 = 6, min.pct = 0.25)

I_small_GO <- enrichGO(gene = rownames(subset(I_small.markers, avg_logFC > 0.5)),
                          universe = row.names(combined@assays$RNA@data),
                          OrgDb = org.Hs.eg.db,
                          keyType = 'SYMBOL',
                          ont = "CC",
                          pAdjustMethod = "BH",
                          pvalueCutoff  = 0.01,
                          qvalueCutoff  = 0.05)

dotplot(I_small_GO, showCategory=20)

UT_small_GO <- enrichGO(gene = rownames(subset(UT_small.markers, avg_logFC > 0.5)),
                       universe = row.names(combined@assays$RNA@data),
                       OrgDb = org.Hs.eg.db,
                       keyType = 'SYMBOL',
                       ont = "CC",
                       pAdjustMethod = "BH",
                       pvalueCutoff  = 0.01,
                       qvalueCutoff  = 0.05)

dotplot(UT_small_GO, showCategory=20)

# Untreated vs Idling - Large
I_large.markers <- FindMarkers(object = combined, ident.1 = c(0,2,4),
                               ident.2 = c(1,3,5,8), min.pct = 0.25)
UT_large.markers <- FindMarkers(object = combined, ident.1 = c(1,3,5,8),
                                ident.2 = c(0,2,4), min.pct = 0.25)

I_large_GO <- enrichGO(gene = rownames(subset(I_large.markers, avg_logFC > 0.5)),
                       universe = row.names(combined@assays$RNA@data),
                       OrgDb = org.Hs.eg.db,
                       keyType = 'SYMBOL',
                       ont = "CC",
                       pAdjustMethod = "BH",
                       pvalueCutoff  = 0.01,
                       qvalueCutoff  = 0.05)

dotplot(I_large_GO, showCategory=20)

UT_large_GO <- enrichGO(gene = rownames(subset(UT_large.markers, avg_logFC > 0.5)),
                        universe = row.names(combined@assays$RNA@data),
                        OrgDb = org.Hs.eg.db,
                        keyType = 'SYMBOL',
                        ont = "CC",
                        pAdjustMethod = "BH",
                        pvalueCutoff  = 0.01,
                        qvalueCutoff  = 0.05)

dotplot(UT_large_GO, showCategory=20)


### use lineages
# keep_UT <- as.data.frame(sort(table(untreated@meta.data$lineage), decreasing = TRUE))
# keep_UT <- subset(keep_UT, Var1 != "XXX")
# keep_UT_30 <- subset(keep_UT, Freq >= 30)
# keep_I <- as.data.frame(sort(table(idling@meta.data$lineage), decreasing = TRUE))
# keep_I <- subset(keep_I, Var1 != "XXX")
# keep_I_30 <- subset(keep_I, Freq >= 30)
# bc_30 <- Reduce(intersect, list(keep_UT_30$Var1, keep_I_30$Var1))

# Use top25 untreated lineages
load("top25BCs.RData")
bcs_top25 <- as.character(unique(bcNum_prop_compare_sub$Barcode))
bc_not <- setdiff(combined@meta.data$lineage, bcs_top25)
test <- combined@meta.data
test$lineageColored[test$lineage %in% bcs_top25] <- test$lineage[test$lineage %in% bcs_top25]
test$lineageColored[test$lineage %in% bc_not] <- "XXX"
lineageColored <- subset(test, select = c("lineageColored"))
combined <- AddMetaData(combined, lineageColored, col.name = "lineageColored")
# combined <- SetAllIdent(object = combined, id = "lineageColored")

cols <- c("CTGAGTCAGAGTGACACACT" = rainbow(25)[1],
          "CTGAGAGTGAGTCTGTCAGT" = rainbow(25)[2],
          "CTGAGTGTGAGAGTGTGTGT" = rainbow(25)[3],
          "CTGACAGTGTCACACAGTGA" = rainbow(25)[4],
          "CTGAGTCACTCACTGAGTGT" = rainbow(25)[5],
          "CTGAGACTCAGACAGACACT" = rainbow(25)[6],
          "CTGAGAGACTCTGTGACTGA" = rainbow(25)[7],
          "CTGACAGACAGTGACTGTCT" = rainbow(25)[8],
          "CTGACTGTCAGACAGAGTGA" = rainbow(25)[9],
          "CTGAGTCAGTCACACTCTGT" = rainbow(25)[10],
          "CTGACAGTGTGTCAGTCTCT" = rainbow(25)[11],
          "CTGAGTGTGACTGTGTGTGA" = rainbow(25)[12],
          "CTGACAGACACTCTCAGTCT" = rainbow(25)[13],
          "CTGACTGTCTGTCAGTGTGT" = rainbow(25)[14],
          "CTGACTGTGTGTCAGTGTGA" = rainbow(25)[15],
          "CTGAGAGTCACTGAGTGTGT" = rainbow(25)[16],
          "CTGAGTGTCACTCTCTCAGA" = rainbow(25)[17],
          "CTGACACACTGTGACTGTGT" = rainbow(25)[18],
          "CTGAGAGTGACAGACTCAGT" = rainbow(25)[19],
          "CTGAGTGTGAGTCTGTGACA" = rainbow(25)[20],
          "CTGACAGTCACACTGACTCA" = rainbow(25)[21],
          "CTGACACACTCTCACTGACA" = rainbow(25)[22],
          "CTGACTCAGTCTGTCTGTCA" = rainbow(25)[23],
          "CTGACTGAGAGTGAGTCACA" = rainbow(25)[24],
          "CTGAGTGACTGTGAGACTGA" = rainbow(25)[25],
          "XXX" = "grey80")

cols1 <- c("CTGAGTCAGAGTGACACACT" = rainbow(25)[1],
          "CTGAGAGTGAGTCTGTCAGT" = "grey80",
          "CTGAGTGTGAGAGTGTGTGT" = "grey80",
          "CTGACAGTGTCACACAGTGA" = "grey80",
          "CTGAGTCACTCACTGAGTGT" = "grey80",
          "CTGAGACTCAGACAGACACT" = "grey80",
          "CTGAGAGACTCTGTGACTGA" = "grey80",
          "CTGACAGACAGTGACTGTCT" = "grey80",
          "CTGACTGTCAGACAGAGTGA" = "grey80",
          "CTGAGTCAGTCACACTCTGT" = "grey80",
          "CTGACAGTGTGTCAGTCTCT" = "grey80",
          "CTGAGTGTGACTGTGTGTGA" = "grey80",
          "CTGACAGACACTCTCAGTCT" = "grey80",
          "CTGACTGTCTGTCAGTGTGT" = "grey80",
          "CTGACTGTGTGTCAGTGTGA" = "grey80",
          "CTGAGAGTCACTGAGTGTGT" = "grey80",
          "CTGAGTGTCACTCTCTCAGA" = "grey80",
          "CTGACACACTGTGACTGTGT" = "grey80",
          "CTGAGAGTGACAGACTCAGT" = "grey80",
          "CTGAGTGTGAGTCTGTGACA" = "grey80",
          "CTGACAGTCACACTGACTCA" = "grey80",
          "CTGACACACTCTCACTGACA" = "grey80",
          "CTGACTCAGTCTGTCTGTCA" = "grey80",
          "CTGACTGAGAGTGAGTCACA" = "grey80",
          "CTGAGTGACTGTGAGACTGA" = "grey80",
          "XXX" = "grey80")

labels = c("1","2","3","4","5","6","7","8","9",
           "10","11","12","13","14","15","16",
           "17","18","19","20","21","22",
           "23","24","25","others")
# DimPlot(combined, reduction = "umap", group.by = "lineageColored") +
#   # geom_text(x=7, y=6, label=sprintf("%d cells", nrow(combined@meta.data))) +
#   scale_color_manual(
#     values = cols,
#     labels = labels
#   ) + theme_bw() +
#   theme(legend.position = "none", legend.title = element_blank(),
#         plot.title = element_text(size = 14, hjust = 0.5, face = "bold"),
#         axis.text=element_text(size=14), legend.text = element_text(size=10),
#         axis.title=element_text(size=14), panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank()) +
#   ggsave("umap_combined_lineageID_greater30_legendSpecial.pdf", width = 8, height = 5)

tint <- data.frame(ifelse(combined@meta.data$lineageColored == "XXX", 0.1, 1),
                   row.names = rownames(combined@meta.data))
names(tint) <- "Tint"
combined <- AddMetaData(combined, tint, 'Tint')
# tint1 <- data.frame(ifelse(combined@meta.data$lineageColored == "CTGACAGACACTCTCAGTCT", 1, 0.1),
#                     row.names = rownames(combined@meta.data))
# names(tint1) <- "Tint1"
# combined <- AddMetaData(combined, tint1, 'Tint1')
test_subset <- FetchData(combined, vars = c("UMAP_1", "UMAP_2", "lineage", "S.Score", "G2M.Score",
                                                 "Phase", "old.ident", "lineageColored", "Tint", "Tint1"))

test_subset$lineageColored <- factor(test_subset$lineageColored, levels = unique(bcNum_prop_compare_sub$Barcode))
plt <- ggplot(data = test_subset, aes(x = UMAP_1, y = UMAP_2, color = lineageColored, alpha = Tint)) +
  geom_point() + scale_alpha_continuous(range = c(0.1,1)) + 
  theme_bw() +
  scale_color_manual(values = cols, labels = labels, name = "Barcode") +
  guides(alpha = F) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 14, colour = "black"),
        axis.text.y = element_text(size = 14, colour = "black"),
        legend.position = "right", legend.text = element_text(size = 14),
        plot.title = element_text(size = 16, hjust = 0.5, face = "bold"), axis.text=element_text(size=14),
        legend.title = element_text(size=14), axis.title=element_text(size=14)) +
  guides(col = guide_legend(nrow = 13)) #+ ggtitle("CTGACAGACACTCTCAGTCT") +

plt_leg <- ggpubr::get_legend(plt)
as_ggplot(plt_leg) + ggsave("umap_combined_lineageID_tinted_legendOnly.svg", width = 2.5, height = 4)

# ggplot(data = test_subset, aes(x = UMAP_1, y = UMAP_2, color = lineageColored, alpha = Tint)) +
#   geom_point(size = 0.8) + scale_alpha_continuous(range = c(0.1,1)) + 
#   theme_bw() +
#   scale_color_manual(values = cols, labels = labels, name = "Barcode") +
#   guides(alpha = F) +
#   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#         axis.text.x = element_text(size = 14, colour = "black"),
#         axis.text.y = element_text(size = 14, colour = "black"),
#         legend.position = "none", legend.text = element_text(size = 14),
#         plot.title = element_text(size = 16, hjust = 0.5, face = "bold"), axis.text=element_text(size=14),
#         legend.title = element_text(size=14), axis.title=element_text(size=14)) 
#   ggsave("../figures/umap_combined_lineageID_tinted_legendSpecial.svg", width = 3.5, height = 3)


# ggplot(data = test_subset, aes(x = UMAP_1, y = UMAP_2, color = lineageColored, alpha = Tint1)) +
#   geom_point(size = 0.8) + scale_alpha_continuous(range = c(0.1,1)) + 
#   theme_bw() +
#   scale_color_manual(values = cols1, labels = labels, name = "Barcode") +
#   guides(alpha = F) +
#   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#         axis.text.x = element_text(size = 14, colour = "black"),
#         axis.text.y = element_text(size = 14, colour = "black"),
#         legend.position = "none", legend.text = element_text(size = 14),
#         plot.title = element_text(size = 16, hjust = 0.5, face = "bold"), axis.text=element_text(size=14),
#         legend.title = element_text(size=14), axis.title=element_text(size=14)) 
#   ggsave("../figures/umap_combined_lineageID_tinted_legendSpecial.svg", width = 4, height = 3)

# bcs <- c("CTGACAGACACTCTCAGTCT",
#          "CTGACAGTCACACTGACTCA",
#          "CTGACAGTGTCACACAGTGA",
#          "CTGACAGTGTGTCAGTCTCT",
#          "CTGACTGTCAGACAGAGTGA",
#          "CTGACTGTCTGTCAGTGTGT",
#          "CTGACTGTGTGTCAGTGTGA",
#          "CTGAGACTCAGACAGACACT",
#          "CTGAGAGACTCTGTGACTGA",
#          "CTGAGAGTCACTGAGTGTGT",
#          "CTGAGAGTGACAGACTCAGT",
#          "CTGAGAGTGAGTCTGTCAGT",
#          "CTGAGTCACTCACTGAGTGT",
#          "CTGAGTCAGAGTGACACACT",
#          "CTGAGTCAGTCACACTCTGT",
#          "CTGAGTGTCACTCTCTCAGA",
#          "CTGAGTGTGACTGTGTGTGA",
#          "CTGAGTGTGAGAGTGTGTGT",
#          "CTGAGTGTGACTGTGTGTGA",
#          "CTGAGTGTGAGAGTGTGTGT")

t2 <- as.data.frame(combined@meta.data)
t2$BCnum <- as.integer(factor(t2$lineageColored, levels = unique(bcNum_prop_compare_sub$Barcode)))
BCnum <- subset(t2, select = c("BCnum"))
combined <- AddMetaData(combined, BCnum, col.name = "BCnum")

for (i in seq(25)){
  # print(i)
  tint1 <- data.frame(ifelse(combined@meta.data$BCnum == i, 1, 0.1),
                      row.names = rownames(combined@meta.data))
  names(tint1) <- "Tint1"
  combined <- AddMetaData(combined, tint1, 'Tint1')
  test_subset1 <- FetchData(combined, vars = c("UMAP_1", "UMAP_2", "lineage", "S.Score", "G2M.Score",
                                              "Phase", "old.ident", "lineageColored", "Tint", "Tint1"))
  ggplot() +
    theme_bw() +
    geom_density_2d(data = test_subset1, aes(x = UMAP_1, y = UMAP_2, 
                                             color = old.ident), alpha = 0.4) +
    scale_color_manual(values = c("blue", "red")) +
    geom_point(data = subset(test_subset1, BCnum == i), shape = 21,
               aes(x = UMAP_1, y = UMAP_2, fill = lineageColored),
               size = 0.8, stroke = 0.1) + 
    scale_fill_manual(values = cols, labels = labels) +
    # scale_alpha_continuous(range = c(0.1,1)) +
    # scale_color_manual(values = cols, labels = labels, name = "Barcode") +
    # guides(alpha = "none") +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          axis.text.x = element_text(size = 12, colour = "black"),
          axis.text.y = element_text(size = 12, colour = "black"),
          legend.position = "none", legend.text = element_text(size = 12),
          plot.title = element_text(size = 12, hjust = 0.5, face = "bold"), axis.text=element_text(size=14),
          legend.title = element_text(size=12), axis.title=element_text(size=12)) +
    ggtitle(paste("Barcode",i)) +
    ggsave(paste0("../figures/UMAP_top25BCs_byNum/umap_combined_lineageID_tinted_legendSpecial_top25_Barcode", i, ".svg"), width = 4, height = 3)
}

####

combined_untreated <- subset(combined, subset = orig.ident == "Untreated")
combined_idling <- subset(combined, subset = orig.ident == "Idling")

UT_UMAP <- data.frame(combined_untreated@reductions$umap@cell.embeddings)
UT_UMAP$UMAP1_dist <- (UT_UMAP$UMAP_1 - colMeans(UT_UMAP)[1])**2
UT_UMAP$UMAP2_dist <- (UT_UMAP$UMAP_2 - colMeans(UT_UMAP)[2])**2
UT_UMAP$dist <- sqrt(UT_UMAP$UMAP1_dist + UT_UMAP$UMAP2_dist)

I_UMAP <- data.frame(combined_idling@reductions$umap@cell.embeddings)
I_UMAP$UMAP1_dist <- (I_UMAP$UMAP_1 - colMeans(I_UMAP)[1])**2
I_UMAP$UMAP2_dist <- (I_UMAP$UMAP_2 - colMeans(I_UMAP)[2])**2
I_UMAP$dist <- sqrt(I_UMAP$UMAP1_dist + I_UMAP$UMAP2_dist)

dist_mat <- data.frame(Condition = rep(c("Untreated", "Idling"), 
                                       times = c(nrow(UT_UMAP), 
                                                 nrow(I_UMAP))),
                       Distance = c(UT_UMAP$dist, I_UMAP$dist))


dist_mat_melt <- melt(dist_mat, id.vars = "Condition", measure.vars = "Distance")
dist_mat_melt.summary <- aggregate(dist_mat_melt, by = list(dist_mat_melt$Condition),
                                   FUN = mean)[c(1,4)]
names(dist_mat_melt.summary) <- c("Condition", "Distance")

dist_mat$Condition <- factor(dist_mat$Condition, levels = c("Untreated", "Idling"))
dist_mat_melt.summary$Condition <- factor(dist_mat_melt.summary$Condition, 
                                          levels = c("Untreated", "Idling"))
# ggplot(dist_mat, aes(x=Condition, y=Distance, group = Condition,
#                      color = Condition)) + theme_bw() +
#   geom_jitter(width = 0.2, alpha = 0.1) +
#   geom_crossbar(data = dist_mat_melt.summary, aes(ymin=Distance, ymax=Distance, color = Condition),
#                 size=0.5, width = 0.5, color = "black") +
#   # stat_summary(fun.data=mean_sdl, fun.args = list(mult=1), 
#   #              geom="errorbar", width=0.2) +
#   # stat_summary(fun=mean, geom="point", shape=20, size=5) +
#   # ylim(0,5) +
#   scale_color_manual(values = c("red", "blue")) +
#   labs(x = "Condition", y = "Distance from Centroid") +
#   theme(axis.text.y = element_text(size = 14),
#         legend.position = "none", legend.text = element_text(size = 14),
#         plot.title = element_text(size = 16, hjust = 0.5, face = "bold"), 
#         axis.text=element_text(size=14),
#         legend.title = element_text(size=14), axis.title=element_text(size=14),
#         panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
#   ggsave("SKMEL5_conditions_distFromCentroid_byCondition.svg", width = 2.5, height = 5)


###

df1 <- as.data.frame(combined@meta.data)
df1$State <-  with(df1, ifelse(Phase %in% c("G2M", "S"), yes="Dividing", no="Nondividing"))
State <- subset(df1, select = c("State"))
combined <- AddMetaData(combined, State, col.name = "State")
save(combined, file = "combined_includingState.RData")

# BC-by-BC Number/Proportion of cells in each CC state
stateBC_df <- FetchData(combined, vars = c("orig.ident", "lineageColored", "State"))
stateBC_df <- subset(stateBC_df, lineageColored != "XXX")
t1 <- stateBC_df
t1$BCnum <- as.integer(factor(t1$lineageColored, levels = unique(bcNum_prop_compare_sub$Barcode)))
t1$lineageColored <- factor(t1$lineageColored, levels = unique(bcNum_prop_compare_sub$Barcode))

ggplot(subset(t1, orig.ident == "Untreated"), 
       aes(x=BCnum, group = State, fill = State)) +
  theme_bw() + geom_bar(stat = "count", position = "dodge", color = "black") +
  scale_fill_manual(values = c("green3", "gold")) +
  labs(x = "Barcode", y = "Number of Cells") +
  theme(axis.text.y = element_text(size = 14),
        legend.position = "right", legend.text = element_text(size = 14),
        plot.title = element_text(size = 16, hjust = 0.5, face = "bold"), 
        axis.text=element_text(size=14),
        legend.title = element_blank(), axis.title=element_text(size=14),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 

ggplot(subset(t1, orig.ident == "Idling"), 
       aes(x=BCnum, group = State, fill = State)) +
  theme_bw() + geom_bar(stat = "count", position = "dodge", color = "black") +
  scale_fill_manual(values = c("green3", "gold")) +
  labs(x = "Barcode", y = "Number of Cells") +
  theme(axis.text.y = element_text(size = 14),
        legend.position = "right", legend.text = element_text(size = 14),
        plot.title = element_text(size = 16, hjust = 0.5, face = "bold"), 
        axis.text=element_text(size=14),
        legend.title = element_blank(), axis.title=element_text(size=14),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 

# Percentages
state_UT_table_byBC <- subset(t1, orig.ident == "Untreated")[,c("State", "BCnum")]
state_UT_byBC <- state_UT_table_byBC %>% 
  dplyr::group_by(State, BCnum) %>% 
  dplyr::summarise (n = n()) %>%
  dplyr::group_by(BCnum) %>%
  dplyr::mutate(freq = n / sum(n))
state_UT_byBC$Condition <- "Untreated"

state_I_table_byBC <- subset(t1, orig.ident == "Idling")[,c("State", "BCnum")]
state_I_byBC <- state_I_table_byBC %>% 
  dplyr::group_by(State, BCnum) %>% 
  dplyr::summarise (n = n()) %>%
  dplyr::group_by(BCnum) %>%
  dplyr::mutate(freq = n / sum(n))
state_I_byBC$Condition <- "Idling"

state_byBC <- as.data.frame(rbind(state_UT_byBC, state_I_byBC))

# ggplot(subset(state_byBC, Condition == "Untreated"), aes(x = as.factor(BCnum), y = freq,
#                        group = State, fill = State)) +
#   theme_bw() + geom_bar(position = position_fill(reverse = TRUE),
#                         stat = "identity", color = "black") +
#   scale_fill_manual(values = c("green3", "gold")) +
#   scale_y_continuous(labels = scales::percent) +
#   labs(x = "Barcode", y = "Percentage of Cells") +
#   ggtitle("Untreated Cell Cycle State") +
#   theme(axis.text.y = element_text(size = 14),
#         legend.position = "bottom", legend.text = element_text(size = 14),
#         plot.title = element_text(size = 16, hjust = 0.5), 
#         axis.text=element_text(size=14),
#         legend.title = element_blank(), axis.title=element_text(size=14),
#         panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#         axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
#   ggsave("SKMEL5_UT_CellCycleState_proportion.pdf", width = 8, height = 5)

a <- subset(combined@meta.data, orig.ident == "Idling")
a %>%
  dplyr::group_by(State) %>% 
  dplyr::summarise (n = n()) %>%
  dplyr::mutate(prop = n / sum(n))
  
b <- subset(combined@meta.data, orig.ident == "Untreated")
b %>%
  dplyr::group_by(State) %>% 
  dplyr::summarise (n = n()) %>%
  dplyr::mutate(prop = n / sum(n))

# ggplot(subset(state_byBC, Condition == "Idling"), aes(x = as.factor(BCnum), y = freq,
#                                                          group = State, fill = State)) +
#   theme_bw() + geom_bar(position = position_fill(reverse = TRUE),
#                         stat = "identity", color = "black") +
#   scale_fill_manual(values = c("green3", "gold")) +
#   scale_y_continuous(labels = scales::percent) +
#   labs(x = "Barcode", y = "Percentage of Cells") +
#   ggtitle("Idling Cell Cycle State") +
#   theme(axis.text.y = element_text(size = 14),
#         legend.position = "right", legend.text = element_text(size = 14),
#         plot.title = element_text(size = 16, hjust = 0.5), 
#         axis.text=element_text(size=14),
#         legend.title = element_blank(), axis.title=element_text(size=14),
#         panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#         axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
#   ggsave("SKMEL5_I_CellCycleState_proportion.pdf", width = 8, height = 5)


ggplot(subset(state_byBC, State == "Dividing" & Condition == "Untreated"), aes(x = as.factor(BCnum), y = freq)) +
  theme_bw() + geom_bar(stat = "identity", position = "dodge", color = "black", fill = "red") +
  geom_hline(yintercept = 0.670, linetype = 2) +
  scale_y_continuous(labels = scales::percent) +
  labs(x = "Barcode", y = "Percentage of Dividing Cells") +
  theme(axis.text.y = element_text(size = 14),
        legend.position = "right", legend.text = element_text(size = 14),
        plot.title = element_text(size = 16, hjust = 0.5, face = "bold"), 
        axis.text=element_text(size=14),
        legend.title = element_blank(), axis.title=element_text(size=14),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  ggsave("SKMEL5_UT_CellCycleState_onlyDividing_proportionWAverage.pdf", width = 6, height = 4)

ggplot(subset(state_byBC, State == "Dividing" & Condition == "Idling"), aes(x = as.factor(BCnum), y = freq)) +
  theme_bw() + geom_bar(stat = "identity", position = "dodge", color = "black", fill = "blue") +
  geom_hline(yintercept = 0.166, linetype = 2) +
  scale_y_continuous(labels = scales::percent) +
  labs(x = "Barcode", y = "Percentage of Dividing Cells") +
  theme(axis.text.y = element_text(size = 14),
        legend.position = "right", legend.text = element_text(size = 14),
        plot.title = element_text(size = 16, hjust = 0.5, face = "bold"), 
        axis.text=element_text(size=14),
        legend.title = element_blank(), axis.title=element_text(size=14),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  ggsave("SKMEL5_I_CellCycleState_onlyDividing_proportionWAverage.pdf", width = 6, height = 4)


# ggplot(state_all_table_byBC_melt, aes(x = Condition, y = value,
#                                  group = State, fill = State)) +
#   theme_bw() + geom_bar(stat = "identity", position = "fill", color = "black") +
#   scale_fill_manual(values = c("green3", "orange1")) +
#   scale_y_continuous(labels = scales::percent) +
#   labs(x = "Condition", y = "Percentage of Cells") +
#   theme(axis.text.y = element_text(size = 14),
#         legend.position = "right", legend.text = element_text(size = 14),
#         plot.title = element_text(size = 16, hjust = 0.5, face = "bold"), 
#         axis.text=element_text(size=14),
#         legend.title = element_blank(), axis.title=element_text(size=14),
#         panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 

####
# K-means / Clusters

test_KM <- test_twoState
set.seed(10)
kc <- kmeans(test_KM[,c("UMAP_1", "UMAP_2")], centers = 6)
test_KM$Kclusters <- kc$cluster
ggplot(test_KM, aes(x=UMAP_1, y=UMAP_2, color = Kclusters)) +
  geom_point()
# Kclusters <- subset(test_KM, select = c("Kclusters"))
# combined <- AddMetaData(combined, Kclusters, col.name = "Kclusters")
# DimPlot(combined, reduction = "pca", group.by = "Kclusters")




DimPlot(combined, reduction = "umap", group.by = "seurat_clusters")
test_twoState <- FetchData(combined, vars = c("UMAP_1", "UMAP_2", "Phase", "State", "old.ident", "seurat_clusters"))

plt_treat <- ggplot(test_twoState, aes(x = UMAP_1, y = UMAP_2)) +
  geom_point(aes(color = factor(old.ident, levels = c("Untreated", "Idling"))), size = 1) + 
  scale_color_manual(name = "", values = c("red", "blue")) +
  theme_bw() + xlim(-11.5, 11.5) + ylim(-11.5, 11.5) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 12, colour = "black"),
        axis.text.y = element_text(size = 12, colour = "black"),
        legend.position = "right", legend.text = element_text(size = 12),
        plot.title = element_text(size = 12, hjust = 0.5, face = "bold"), axis.text=element_text(size=14),
        legend.title = element_text(size=12), axis.title=element_text(size=12)) 
plt_treat + ggsave("UMAP_combined_SKMEL5_hg38_qcCCReg_treatmentPoint.svg", width = 4, height = 3)
plt_treat_leg <- ggpubr::get_legend(plt_treat)
as_ggplot(plt_treat_leg) + ggsave("UMAP_treatment_legend.svg", width = 2.5, height = 4)

plt_state <- ggplot(test_twoState, aes(x = UMAP_1, y = UMAP_2)) +
  geom_density_2d(aes(color = factor(old.ident, levels = c("Untreated", "Idling"))), alpha = 0.4) +
  scale_color_manual(name = "", values = c("red", "blue")) +
  geom_point(shape = 21,
             aes(fill = State),
             size = 1, stroke = 0.1) + 
  scale_fill_manual(name = "", values = c("green3", "gold")) +
  theme_bw() + xlim(-11.5, 11.5) + ylim(-11.5, 11.5) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 12, colour = "black"),
        axis.text.y = element_text(size = 12, colour = "black"),
        legend.position = "right", legend.text = element_text(size = 12),
        plot.title = element_text(size = 12, hjust = 0.5, face = "bold"), axis.text=element_text(size=14),
        legend.title = element_text(size=12), axis.title=element_text(size=12)) 
plt_state + ggsave("UMAP_combined_SKMEL5_hg38_qcCCReg_treatmentDensity_CCStatePoint.svg", width = 4, height = 3)
plt_state_leg <- ggpubr::get_legend(plt_state)
as_ggplot(plt_state_leg) + ggsave("UMAP_state_legend.svg", width = 2.5, height = 4)



test_twoState_UT <- subset(test_twoState, old.ident == "Untreated")
test_twoState_UT$Region <- with(test_twoState_UT, 
                                ifelse(seurat_clusters %in% as.factor(c("1","3","5","8")), 
                                       yes="large", no="small"))

test_twoState_I <- subset(test_twoState, old.ident == "Idling")
test_twoState_I$Region <- with(test_twoState_I, 
                             ifelse(seurat_clusters %in% as.factor(c("0","2","4")), 
                                    yes="large", no="small"))


region_UT_table <- test_twoState_UT[,c("State", "Region")]
region_UT_table <- region_UT_table %>% 
  dplyr::group_by(State, Region) %>% 
  dplyr::summarise (n = n()) %>%
  dplyr::group_by(Region) %>%
  dplyr::mutate(freq = n / sum(n))
region_UT_table$Condition <- "Untreated"

region_I_table <- test_twoState_I[,c("State", "Region")]
region_I_table <- region_I_table %>%
  dplyr::group_by(State, Region) %>%
  dplyr::summarise (n = n()) %>%
  dplyr::group_by(Region) %>%
  dplyr::mutate(freq = n / sum(n))
region_I_table$Condition <- "Idling"

# region_all_table_n <- region_all_table %>%
#   dplyr::group_by(name) %>%
#   dplyr::summarise (num = sum(n))

region_all_table <- rbind(region_I_table, region_UT_table)
region_all_table$name <- c("I_L", "I_S", "I_L", "I_S",
                           "UT_L", "UT_S", "UT_L", "UT_S")
region_all_table$name <- factor(region_all_table$name,
                                levels = c("UT_L", "UT_S", "I_L", "I_S"))
# region_all_table$tot <- c()


region_all_table_n <- region_all_table %>%
  dplyr::group_by(name) %>%
  dplyr::summarise (num = sum(n))

region_all_n <- data.frame(name = unique(region_all_table_n$name),
                           label = c(paste("n =", as.character(region_all_table_n$num[1])),
                                     paste("n =", as.character(region_all_table_n$num[2])),
                                     paste("n =", as.character(region_all_table_n$num[3])),
                                     paste("n =", as.character(region_all_table_n$num[4]))))

# ggplot(region_all_table, aes(x = name, y = freq,
#                              group = State, fill = State)) +
#   theme_classic() + geom_bar(stat = "identity", color = "black") +
#   scale_fill_manual(values = c("green3", "gold")) +
#   scale_y_continuous(name = "Percentage of Cells", labels = scales::percent) +
#   scale_x_discrete(name = "Condition",
#                    labels = c(expression(Untreated[large]), 
#                               expression(Untreated[small]), 
#                               expression(Idling[large]),
#                               expression(Idling[small]))) +
#   # labs(y = "Percentage of Cells") +
#   # ggtitle("Cluster Cell Cycle State") +
#   # geom_text(data = region_all_n,
#   #           aes(name, label, label = label)) +
#   theme(axis.text.y = element_text(size = 14),
#         legend.position = "bottom", legend.text = element_text(size = 14),
#         plot.title = element_text(size = 16, hjust = 0.5), 
#         axis.text=element_text(size=14),
#         legend.title = element_blank(), axis.title=element_text(size=14),
#         panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
#   ggsave("SKMEL5_allClusters_CellCycleState_proportion.pdf", width = 6, height = 5)


####

EMD_data_UT <- subset(FetchData(combined, vars = c("UMAP_1", "UMAP_2", "old.ident", "BCnum", "lineageColored")), 
                      old.ident == "Untreated")
D_UT <- dist(cbind(EMD_data_UT$UMAP_1, EMD_data_UT$UMAP_2), 
          diag=TRUE, upper=TRUE) 
m_UT <- as.matrix(D_UT) # coerce dist object to a matrix
dimnames(m_UT) <- list(1:length(EMD_data_UT$UMAP_1),
                    1:length(EMD_data_UT$UMAP_2))
xy_UT <- t(combn(colnames(m_UT), 2))
df_UT_EMD <- data.frame(xy_UT, dist=m_UT[xy_UT])
df_UT_EMD$Condition <- "Untreated"

EMD_data_I <- subset(FetchData(combined, vars = c("UMAP_1", "UMAP_2", "old.ident", "BCnum", "lineageColored")), 
                      old.ident == "Idling")
D_I <- dist(cbind(EMD_data_I$UMAP_1, EMD_data_I$UMAP_2), 
             diag=TRUE, upper=TRUE) 
m_I <- as.matrix(D_I) # coerce dist object to a matrix
dimnames(m_I) <- list(1:length(EMD_data_I$UMAP_1),
                       1:length(EMD_data_I$UMAP_2))
xy_I <- t(combn(colnames(m_I), 2))
df_I_EMD <- data.frame(xy_I, dist=m_I[xy_I])
df_I_EMD$Condition <- "Idling"

df_all_EMD <- rbind(df_UT_EMD, df_I_EMD)

df_all_EMD_UT <- df_UT_EMD[sample(nrow(df_UT_EMD), 15000), ]
df_all_EMD_I <- df_I_EMD[sample(nrow(df_I_EMD), 15000), ]
df_all_EMD_sample <- rbind(df_all_EMD_UT, df_all_EMD_I)

EMD_compare <- wasserstein1d(df_all_EMD_UT$dist, df_all_EMD_I$dist)

df_all_EMD_sample$Condition <- factor(df_all_EMD_sample$Condition,
                                      levels = c("Untreated", "Idling"))
ggplot(df_all_EMD_sample, aes(dist, color = Condition)) + 
  theme_bw() + stat_ecdf(size = 1.25) +
  scale_color_manual(values = c("red", "blue")) +
  labs(x = "Distance", y = "Fraction") +
  annotate("text", x = 17.5, y = 0.2, label = paste("EMD =", round(EMD_compare,2)), size = 5) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 12, colour = "black"),
        axis.text.y = element_text(size = 12, colour = "black"),
        legend.position = "none", legend.text = element_text(size = 12),
        plot.title = element_text(size = 12, hjust = 0.5, face = "bold"), axis.text=element_text(size=14),
        legend.title = element_text(size=12), axis.title=element_text(size=12)) +
  ggsave("ECDF_UT-I_withEMD.pdf", width = 4, height = 3)

# ## BC-by-BC ECDF compared to Treatment condition
# # Untreated
# for (i in seq(25)[-8]){
#   # print(i)
#   bcSub <- subset(EMD_data_UT, BCnum == i)
#   D_bcSub <- dist(cbind(bcSub$UMAP_1, bcSub$UMAP_2), 
#                diag=TRUE, upper=TRUE) 
#   m_bcSub <- as.matrix(D_bcSub) # coerce dist object to a matrix
#   dimnames(m_bcSub) <- list(1:length(bcSub$UMAP_1),
#                          1:length(bcSub$UMAP_2))
#   xy_bcSub <- t(combn(colnames(m_bcSub), 2))
#   df_bcSub_EMD <- data.frame(xy_bcSub, dist=m_bcSub[xy_bcSub])
#   df_bcSub_EMD$BCnum <- i
#   df_bcSub_EMD$lineageColored <- unique(bcSub$lineageColored)
#   ggplot() +
#     theme_bw() +
#     stat_ecdf(data = df_all_EMD_UT, aes(dist), color = "black", size = 1.25) +
#     stat_ecdf(data = df_bcSub_EMD, aes(dist, color = lineageColored)) +
#     scale_color_manual(values = cols, labels = labels) +
#     labs(x = "Distance", y = "Fraction") +
#     theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#           axis.text.x = element_text(size = 12, colour = "black"),
#           axis.text.y = element_text(size = 12, colour = "black"),
#           legend.position = "none", legend.text = element_text(size = 12),
#           plot.title = element_text(size = 12, hjust = 0.5, face = "bold"), axis.text=element_text(size=14),
#           legend.title = element_text(size=12), axis.title=element_text(size=12)) +
#     ggtitle(paste("Barcode",i)) +
#     ggsave(paste0("ECDF_top25BCs_byNum/ecdf_combined_lineageID_top25_Barcode", i, ".svg"), width = 4, height = 3)
# }
# 
# 
# bc_list_UT = list()
# 
# for (i in seq(25)[-8]){
#   # print(i)
#   bcSub <- subset(EMD_data_UT, BCnum == i)
#   D_bcSub <- dist(cbind(bcSub$UMAP_1, bcSub$UMAP_2), 
#                   diag=TRUE, upper=TRUE) 
#   m_bcSub <- as.matrix(D_bcSub) # coerce dist object to a matrix
#   dimnames(m_bcSub) <- list(1:length(bcSub$UMAP_1),
#                             1:length(bcSub$UMAP_2))
#   xy_bcSub <- t(combn(colnames(m_bcSub), 2))
#   df_bcSub_EMD <- data.frame(xy_bcSub, dist=m_bcSub[xy_bcSub])
#   df_bcSub_EMD$BCnum <- i
#   df_bcSub_EMD$lineageColored <- unique(bcSub$lineageColored)
#   bc_list_UT[[i]] <- df_bcSub_EMD
# }
# 
# all_bc_list_UT = do.call(rbind, bc_list_UT)
# all_bc_list_UT$lineageColored <- factor(all_bc_list_UT$lineageColored, 
#                                         levels = unique(bcNum_prop_compare_sub$Barcode))
# ggplot() +
#   theme_bw() + stat_ecdf(data = all_bc_list_UT, aes(dist, color = lineageColored)) +
#   stat_ecdf(data = df_all_EMD_UT, aes(dist), color = "black", size = 1.25) +
#   scale_color_manual(values = cols, labels = as.character(seq(25)[-8])) +
#   labs(x = "Distance", y = "Fraction") +
#   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#         axis.text.x = element_text(size = 12, colour = "black"),
#         axis.text.y = element_text(size = 12, colour = "black"),
#         legend.position = "none", legend.text = element_text(size = 12),
#         plot.title = element_text(size = 12, hjust = 0.5, face = "bold"), axis.text=element_text(size=14),
#         legend.title = element_text(size=12), axis.title=element_text(size=12)) +
#   ggtitle("Untreated Barcode Distribution") +
#   ggsave(paste0("ECDF_top25BCs_plottedTogether_UT.svg"), width = 4, height = 3)
# 
# 
# # Idling
# for (i in seq(25)[-8]){
#   # print(i)
#   bcSub <- subset(EMD_data_I, BCnum == i)
#   D_bcSub <- dist(cbind(bcSub$UMAP_1, bcSub$UMAP_2), 
#                   diag=TRUE, upper=TRUE) 
#   m_bcSub <- as.matrix(D_bcSub) # coerce dist object to a matrix
#   dimnames(m_bcSub) <- list(1:length(bcSub$UMAP_1),
#                             1:length(bcSub$UMAP_2))
#   xy_bcSub <- t(combn(colnames(m_bcSub), 2))
#   df_bcSub_EMD <- data.frame(xy_bcSub, dist=m_bcSub[xy_bcSub])
#   df_bcSub_EMD$BCnum <- i
#   df_bcSub_EMD$lineageColored <- unique(bcSub$lineageColored)
#   ggplot() +
#     theme_bw() +
#     stat_ecdf(data = df_all_EMD_I, aes(dist), color = "black", size = 1.25) +
#     stat_ecdf(data = df_bcSub_EMD, aes(dist, color = lineageColored)) +
#     scale_color_manual(values = cols, labels = labels) +
#     labs(x = "Distance", y = "Fraction") +
#     theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#           axis.text.x = element_text(size = 12, colour = "black"),
#           axis.text.y = element_text(size = 12, colour = "black"),
#           legend.position = "none", legend.text = element_text(size = 12),
#           plot.title = element_text(size = 12, hjust = 0.5, face = "bold"), axis.text=element_text(size=14),
#           legend.title = element_text(size=12), axis.title=element_text(size=12)) +
#     ggtitle(paste("Barcode",i)) +
#     ggsave(paste0("ECDF_top25BCs_byNum/ecdf_combined_lineageID_top25_Barcode", i, ".svg"), width = 4, height = 3)
# }
# 
# 
# bc_list_I = list()
# 
# for (i in seq(25)[-8]){
#   # print(i)
#   bcSub <- subset(EMD_data_I, BCnum == i)
#   D_bcSub <- dist(cbind(bcSub$UMAP_1, bcSub$UMAP_2), 
#                   diag=TRUE, upper=TRUE) 
#   m_bcSub <- as.matrix(D_bcSub) # coerce dist object to a matrix
#   dimnames(m_bcSub) <- list(1:length(bcSub$UMAP_1),
#                             1:length(bcSub$UMAP_2))
#   xy_bcSub <- t(combn(colnames(m_bcSub), 2))
#   df_bcSub_EMD <- data.frame(xy_bcSub, dist=m_bcSub[xy_bcSub])
#   df_bcSub_EMD$BCnum <- i
#   df_bcSub_EMD$lineageColored <- unique(bcSub$lineageColored)
#   bc_list_I[[i]] <- df_bcSub_EMD
# }
# 
# all_bc_list_I = do.call(rbind, bc_list_I)
# all_bc_list_I$lineageColored <- factor(all_bc_list_I$lineageColored, 
#                                        levels = unique(bcNum_prop_compare_sub$Barcode))
# 
# ggplot() +
#   theme_bw() + stat_ecdf(data = all_bc_list_I, aes(dist, color = lineageColored)) +
#   stat_ecdf(data = df_all_EMD_I, aes(dist), color = "black", size = 1.25) +
#   scale_color_manual(values = cols, labels = as.character(seq(25)[-8])) +
#   labs(x = "Distance", y = "Fraction") +
#   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#         axis.text.x = element_text(size = 12, colour = "black"),
#         axis.text.y = element_text(size = 12, colour = "black"),
#         legend.position = "none", legend.text = element_text(size = 12),
#         plot.title = element_text(size = 12, hjust = 0.5, face = "bold"), axis.text=element_text(size=14),
#         legend.title = element_text(size=12), axis.title=element_text(size=12)) +
#   ggtitle("Idling Barcode Distribution") +
#   ggsave(paste0("ECDF_top25BCs_plottedTogether_I.svg"), width = 4, height = 3)
# 
# 
# ### EMD
# 
# was_list_UT = list()
# was_list_I = list()
# for (i in seq(25)[-8]){
#   was_BCsub_UT <- subset(all_bc_list_UT, BCnum == i)
#   was_BCsub_I <- subset(all_bc_list_I, BCnum == i)
#   was_list_UT[[i]] <- wasserstein1d(df_all_EMD_UT$dist, was_BCsub_UT$dist)
#   was_list_I[[i]] <- wasserstein1d(df_all_EMD_I$dist, was_BCsub_I$dist)
# }
# 
# was_df <- data.frame(BCnum = seq(25)[-8],
#                      Untreated = unlist(was_list_UT),
#                      Idling = unlist(was_list_I))
# 
# was_df_melt <- melt(was_df, id.vars = "BCnum", 
#                     measure.vars = c("Untreated", "Idling"))
# names(was_df_melt) <- c("BCnum", "Condition", "EMD")
# 
# ggplot(data = was_df_melt, aes(x=factor(BCnum), y=EMD, group = Condition, fill = Condition)) +
#   geom_bar(stat = "identity", position = "dodge", color = "black") +
#   scale_fill_manual(name = "", values = c("red", "blue")) +
#   theme_bw() +
#   labs(x = "Barcode", y = "EMD (from Condition Null)") +
#   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#         axis.text.x = element_text(size = 12, colour = "black"),
#         axis.text.y = element_text(size = 12, colour = "black"),
#         legend.position = "none", legend.text = element_text(size = 12),
#         plot.title = element_text(size = 12, hjust = 0.5, face = "bold"), axis.text=element_text(size=14),
#         legend.title = element_text(size=12), axis.title=element_text(size=12)) +
#   # ggtitle("Idling Barcode Distribution") +
#   ggsave(paste0("EMD_top25BCs_UT-I_plottedTogether.pdf"), width = 8, height = 4)
# 
# ## Do EMD fold change plot
# was_df$FC <- log2(was_df$Idling / was_df$Untreated)
# ggplot(data = was_df, aes(x=factor(BCnum), y=FC)) +
#   geom_bar(stat = "identity", position = "dodge", color = "black", fill = "grey") +
#   theme_bw() +
#   labs(x = "Barcode", y = "Log2 EMD Difference") +
#   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#         axis.text.x = element_text(size = 12, colour = "black"),
#         axis.text.y = element_text(size = 12, colour = "black"),
#         legend.position = "none", legend.text = element_text(size = 12),
#         plot.title = element_text(size = 12, hjust = 0.5, face = "bold"), axis.text=element_text(size=14),
#         legend.title = element_text(size=12), axis.title=element_text(size=12)) +
#   # ggtitle("Idling Barcode Distribution") +
#   ggsave(paste0("EMD_top25BCs_UT-I_FC_plottedTogether.pdf"), width = 8, height = 4)
# 
# ##
# baseline <- wasserstein1d(df_all_EMD_UT$dist, df_all_EMD_I$dist)
# was_list_t = list()
# for (i in seq(25)[-8]){
#   was_BCsub_UT <- subset(all_bc_list_UT, BCnum == i)
#   was_BCsub_I <- subset(all_bc_list_I, BCnum == i)
#   was_list_t[[i]] <- wasserstein1d(was_BCsub_UT$dist, was_BCsub_I$dist)
# }
# 
# was_df_comp <- data.frame(BCnum = seq(25)[-8],
#                           EMD = unlist(was_list_t),
#                           baseline = baseline)
# was_df_comp$EMD_c <- was_df_comp$EMD - was_df_comp$baseline
# 
# ggplot(data = was_df_comp, aes(x=factor(BCnum), y=EMD_c)) +
#   geom_bar(stat = "identity", color = "black")
# 
# DimPlot(combined, reduction = "umap", group.by = "lineage") +
#   # geom_text(x=7, y=6, label=sprintf("%d cells", nrow(combined@meta.data))) +
#   theme_bw() +
#   theme(legend.position = "none", legend.title = element_blank(),
#         plot.title = element_text(size = 14, hjust = 0.5, face = "bold"),
#         axis.text=element_text(size=14), legend.text = element_text(size=10),
#         axis.title=element_text(size=14), panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank()) +
#   ggsave("umap_combined_lineageID.svg", width = 8, height = 6)
