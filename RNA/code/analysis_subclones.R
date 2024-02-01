library(biomaRt)
library(DESeq2)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(ggpubr)  
library(org.Hs.eg.db)
library(clusterProfiler)

# set working directory to RNA folder
if (Sys.getenv("RSTUDIO") == "1") {
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
}

## Get the data in the right format
d <- read.csv("../data/featureCounts_matrix_all.csv", header=T, sep=",")

#Rename columns
cols <- c("ensembl_gene_id", "SC01_day0_rep1", "SC01_day0_rep2", "SC01_day0_rep3",
          "SC01_day3_rep1", "SC01_day3_rep2", "SC01_day3_rep3",
          "SC01_day8_rep1", "SC01_day8_rep2", "SC01_day8_rep3",
          "SC07_day0_rep1", "SC07_day0_rep2", "SC07_day0_rep3",
          "SC07_day3_rep1", "SC07_day3_rep2", "SC07_day3_rep3",
          "SC07_day8_rep1", "SC07_day8_rep2", "SC07_day8_rep3",
          "SC10_day0_rep1", "SC10_day0_rep2", "SC10_day0_rep3",
          "SC10_day3_rep1", "SC10_day3_rep2", "SC10_day3_rep3",
          "SC10_day8_rep1", "SC10_day8_rep2", "SC10_day8_rep3")
names(d) <- cols
ensembl <- useMart("ensembl")
mart <- useDataset("hsapiens_gene_ensembl", mart = ensembl)
genes <- d$ensembl_gene_id
G_list <- getBM(attributes= c("ensembl_gene_id","hgnc_symbol"),
                filters= "ensembl_gene_id",
                values=genes,
                mart=mart)

GE_data <- merge(d, G_list, by = "ensembl_gene_id")
d <- GE_data[, -1]
d <- d[c(28, seq(1:27))]
rownames(d) <- make.names(d$hgnc_symbol, unique = T)
d <- d[, 2:28]

countdata <- d
# baseline <- c(1,2,3,10,11,12,19,20,21)
# treat3d  <- c(4,5,6,13,14,15,22,23,24)
# treat8d  <- c(7,8,9,16,17,18,25,26,27)
# # define the groups by subclones
# sc01 <- c(baseline[1:3], treat3d[1:3], treat8d[1:3])
# sc07 <- c(baseline[4:6], treat3d[4:6], treat8d[4:6])
# sc10 <- c(baseline[7:9], treat3d[7:9], treat8d[7:9])
# # Get the countdata specific to conditions: 
# # countdata <- countdata[,c(baseline)] 
# rownames(countdata) <- d[,"ensembl_gene_id"]

# save(countdata, file = "merged_countdata.RData")

## Creating a DEseq2 object
condition <- c("0", "3", "8")
treatment <- rep(condition, each=3) # Three biological replicates
unique(treatment)
cell <- c("SC01", "SC07","SC10") #sublines used for the analysis
cellName <- rep(cell, each=3)

coldata <- data.frame(cell=rep(cellName), treatment=rep(treatment, each=3))
group = factor(paste(coldata$cell, coldata$treatment, sep="."))
coldata$group = group

dds <- DESeqDataSetFromMatrix(countData = countdata,
                              colData = coldata,
                              design = ~ cell + treatment + cell:treatment)

# dds2 <- dds[rowSums(counts(dds)) > 18, ] # remove rows with minimum of 2 read per condition
# save(dds2, file = "DDS_SC-1,7,10_cell-treat-int.RData")

dds2 <- dds[rowSums(counts(dds)) > 18, ] # remove rows with minimum of 2 read per condition
nrow(dds2)
# save(dds2, file = "DDS_SC-1,7,10_cell-treat-int.RData")

## Plot normalized data
rld <- rlog(dds2, blind = FALSE)
# save(rld, file = "RLD_SC-1,7,10_0,3,8d_20180701.RData")
plotPCA(rld, intgroup = c("cell", "treatment"), ntop=5000)

## Use prcomp function
# Colored by cell line, shape by time point, lines connecting time
pca_DEseq <- prcomp(t(assay(rld)))
pca_DEseq_perc <- round(100*pca_DEseq$sdev^2/sum(pca_DEseq$sdev^2),1)
pca_DEseq_df <- data.frame(PC1 = pca_DEseq$x[,1], 
                           PC2 = pca_DEseq$x[,2], 
                           sample = colnames(assay(rld)),
                           cell.line = rep(c("SC01", "SC07", "SC10"), each = 9),
                           day = rep(c("Day0", "Day3", "Day8"), each = 3),
                           replicate = rep(c("Rep1", "Rep2", "Rep3"), times=9))

library(plyr)
pca_DEseq_means <- ddply(pca_DEseq_df, .(cell.line, day), summarise, meanPC1 = mean(PC1), meanPC2 = mean(PC2))

plt_pca <- ggplot(pca_DEseq_df, aes(PC1,PC2, color = cell.line))+
  geom_point(aes(shape = day), size=2) +
  geom_path(data = pca_DEseq_means, 
            aes(x=meanPC1, y=meanPC2,
                color=cell.line), arrow = arrow(),
            size = 0.8) +
  labs(x=paste0("PC1 (",pca_DEseq_perc[1],"% variance)"), y=paste0("PC2 (",pca_DEseq_perc[2],"% variance)")) +
  theme_bw() + scale_color_manual(values = c("orange", "cyan", "purple")) +
  # ggtitle("PCA - Subclones in Time") +
  theme(legend.text = element_text(size = 12), 
        plot.title = element_text(size = 12, 
                                  hjust = 0.5), 
        axis.text=element_text(size=12),
        legend.title = element_text(size=12),
        legend.position = "right",
        axis.title=element_text(size=12),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())

library(svglite)
plt_pca 
# ggsave("SKMEL5_sublines_timeSeriesRNA_rld_Fig.svg") #, width = 4, height = 3)
ggsave("SKMEL5_sublines_timeSeriesRNA_rld_Fig.pdf") #, width = 4, height = 3)

plt_pca_leg <- get_legend(plt_pca)
as_ggplot(plt_pca_leg) 
# ggsave("SKMEL5_sublines_timeSeriesRNA_rld_Figleg.svg")
ggsave("SKMEL5_sublines_timeSeriesRNA_rld_Figleg.pdf")


## Differential Expression Analysis
dds <- DESeq(dds)
res_0to8d <- results(dds, name="treatment_8_vs_0", cooksCutoff = 0.99, 
                     independentFiltering = TRUE, alpha = 0.05, pAdjustMethod = "BH")
summary(res_0to8d)
# order results table by the smallest adjusted p value:
res_0to8d <- res_0to8d[order(res_0to8d$padj),]
results_0to8d <- as.data.frame(res_0to8d)

results_0to8d <- mutate(results_0to8d, sig=ifelse(results_0to8d$padj<0.05 & results_0to8d$log2FoldChange > 2, "Upregulated", ifelse(results_0to8d$padj<0.05 & results_0to8d$log2FoldChange < -2, "Downregulated", "Not Significant")))

row.names(results_0to8d) <- row.names(res_0to8d)

DEgenes_0to8d <- results_0to8d[which(abs(results_0to8d$log2FoldChange) > log2(1.5) & results_0to8d$padj < 0.05),]

volcano <- ggplot(results_0to8d, aes(log2FoldChange, -log10(pvalue))) +
  geom_point(aes(col = sig)) + theme_bw() +
  scale_color_manual(values = c("red", "grey", "green3")) +
  # ggtitle("Volcano Plot of Untreated vs Idling") +
  labs(x="log2(Fold Change)", y="Log(Odds Ratio)") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.text = element_text(size = 12),
        plot.title = element_text(size = 14, hjust = 0.5, face = "bold"), 
        axis.text=element_text(size=12),
        legend.title = element_text(size=12), 
        axis.title=element_text(size=12),
        legend.position = "none")

# save(results_0to8d, file="untreatedIdling_DEA.RData")

## GO Overrepresentation Analysis
OrgDB <- org.Hs.eg.db
upreg_genes <- subset(results_0to8d, padj<0.05 & log2FoldChange>2)
downreg_genes <-subset(results_0to8d, padj<0.05 & log2FoldChange<(-2))

geneList_up <- as.vector(upreg_genes$log2FoldChange)
names(geneList_up) <- rownames(upreg_genes)
geneList_down <- as.vector(downreg_genes$log2FoldChange)
names(geneList_down) <- rownames(downreg_genes)

genes_up <- as.vector(rownames(upreg_genes))
genes_down <- as.vector(rownames(downreg_genes))
# names(geneList) <- rownames(results_0to8d)
genes_up_ENTREZID <- bitr(genes_up, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = OrgDB)$ENTREZID
genes_down_ENTREZID <- bitr(genes_down, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = OrgDB)$ENTREZID

# Group GO
ggo_up <- clusterProfiler::groupGO(gene     = genes_up_ENTREZID,
                                   OrgDb    = OrgDB,
                                   ont      = "BP",
                                   level    = 3,
                                   readable = TRUE)
ggo_up_df <- as.data.frame(ggo_up)
ggo_up_df <- ggo_up_df[order(-ggo_up_df$Count),] 

ggo_down <- clusterProfiler::groupGO(gene = genes_down_ENTREZID,
                                     OrgDb    = OrgDB,
                                     ont      = "BP",
                                     level    = 3,
                                     readable = TRUE)

# GO over-representation test
ego_genesUp_BP <- clusterProfiler::enrichGO(gene  = genes_up_ENTREZID,
                                            OrgDb         = OrgDB,
                                            ont           = "BP",
                                            pAdjustMethod = "BH",
                                            pvalueCutoff  = 0.05,
                                            qvalueCutoff  = 0.05, 
                                            readable      = TRUE)
ego_genesUp_MF <- clusterProfiler::enrichGO(gene  = genes_up_ENTREZID,
                                            OrgDb         = OrgDB,
                                            ont           = "MF",
                                            pAdjustMethod = "BH",
                                            pvalueCutoff  = 0.05,
                                            qvalueCutoff  = 0.05, 
                                            readable      = TRUE)
ego_genesUp_CC <- clusterProfiler::enrichGO(gene  = genes_up_ENTREZID,
                                            OrgDb         = OrgDB,
                                            ont           = "CC",
                                            pAdjustMethod = "BH",
                                            pvalueCutoff  = 0.05,
                                            qvalueCutoff  = 0.05, 
                                            readable      = TRUE)


ego_genesDown <- clusterProfiler::enrichGO(gene  = genes_down_ENTREZID,
                                           OrgDb         = OrgDB,
                                           ont           = "BP",
                                           pAdjustMethod = "BH",
                                           pvalueCutoff  = 0.05,
                                           qvalueCutoff  = 0.05, 
                                           readable      = TRUE)

## This creates FIG 3A
dotplot(ego_genesUp_MF, font.size = 14, label_format = 40) + 
  # ggtitle("GO Over-representation Upregulated Genes for RNA-seq") +
  labs(x="Gene Ratio") +
  theme(legend.text = element_text(size = 12),
        plot.title = element_text(size = 14, hjust = 0.5, face = "bold"), 
        axis.text=element_text(size=12),
        legend.title = element_text(size=12,face="bold"), 
        axis.title=element_text(size=12, face="bold"))
ggsave("GOenrichment_genesUp_MF.pdf")

# this is not FIG 3A because it calculated BP-type GO, not MF to compare w/ ATAC
dotplot(ego_genesDown, font.size = 14, label_format = 40) + 
  # ggtitle("GO Over-representation Downregulated Genes") +
  labs(x="Gene Ratio") +
  theme(legend.text = element_text(size = 12),
        plot.title = element_text(size = 14, hjust = 0.5, face = "bold"), 
        axis.text=element_text(size=12),
        legend.title = element_text(size=12,face="bold"), 
        axis.title=element_text(size=12, face="bold"))

