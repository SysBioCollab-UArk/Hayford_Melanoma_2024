library(ggplot2)
library(ggpubr)
library(ggrepel)
library(org.Hs.eg.db)
library(clusterProfiler)
library(Rsamtools)
library(BSgenome.Hsapiens.UCSC.hg38)
library(ChIPseeker)
library(soGGi)
library(MotifDb)
library(Biostrings)
library(BSgenome.Hsapiens.UCSC.hg38)
library(GenomicAlignments)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(svglite)
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

# if (Sys.getenv("RSTUDIO") == "1") {
#   setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
# }

## RNAseq ##
load(file=file.path("..", "data", "untreatedIdling_DEA.RData"))

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
## ##

## ATAC-seq ##
## Non-downsampled idling
# load("/Volumes/Transcend/ATACseq/SKMEL5_ATAC_annotatedPeaks_uniqueShared.RData")
## Downsampled idling
load(file.path("..", "..", "ATAC", "data", 
               "SKMEL5_ATACsub25_annotatedPeaks_uniqueShared.RData"))

go_UT_BP <- enrichGO(as.data.frame(as.GRanges(anno_UT))$geneId, 
                     OrgDb = "org.Hs.eg.db", ont = "BP", maxGSSize = 5000)
go_I_BP <- enrichGO(as.data.frame(as.GRanges(anno_I))$geneId, 
                    OrgDb = "org.Hs.eg.db", ont = "BP", maxGSSize = 5000)
go_con_BP <- enrichGO(as.data.frame(as.GRanges(anno_con))$geneId, 
                      OrgDb = "org.Hs.eg.db", ont = "BP", maxGSSize = 5000)

go_UT_MF <- enrichGO(as.data.frame(as.GRanges(anno_UT))$geneId, 
                     OrgDb = "org.Hs.eg.db", ont = "MF", maxGSSize = 5000)
go_I_MF <- enrichGO(as.data.frame(as.GRanges(anno_I))$geneId, 
                    OrgDb = "org.Hs.eg.db", ont = "MF", maxGSSize = 5000)
go_con_MF <- enrichGO(as.data.frame(as.GRanges(anno_con))$geneId, 
                      OrgDb = "org.Hs.eg.db", ont = "MF", maxGSSize = 5000)

go_UT_CC <- enrichGO(as.data.frame(as.GRanges(anno_UT))$geneId, 
                     OrgDb = "org.Hs.eg.db", ont = "CC", maxGSSize = 5000)
go_I_CC <- enrichGO(as.data.frame(as.GRanges(anno_I))$geneId, 
                    OrgDb = "org.Hs.eg.db", ont = "CC", maxGSSize = 5000)
go_con_CC <- enrichGO(as.data.frame(as.GRanges(anno_con))$geneId, 
                      OrgDb = "org.Hs.eg.db", ont = "CC", maxGSSize = 5000)
## ##

## Pull top GO terms from each data modality
# BP
RNA_Up_BP <- 
  subset(ego_genesUp_BP@result, qvalue < 0.05)[,c("ID", "Description", "GeneRatio", "qvalue")]
RNA_Up_BP$logp <- -log(RNA_Up_BP$qvalue)
ATAC_Up_BP <- 
  subset(go_I_BP@result, qvalue < 0.05)[,c("ID", "Description", "GeneRatio", "qvalue")]
ATAC_Up_BP$logp <- -log(ATAC_Up_BP$qvalue)

Up_all_BP <- merge(x=RNA_Up_BP, y=ATAC_Up_BP, by="Description", all.x = T, all.y = T)
Up_all_BP_complete <- Up_all_BP[complete.cases(Up_all_BP),]
Up_all_BP_complete$sort <- Up_all_BP_complete$logp.x + Up_all_BP_complete$logp.y
Up_all_BP_complete_sorted <- Up_all_BP_complete[order(-Up_all_BP_complete$sort),]

ggplot(Up_all_BP_complete_sorted, aes(x=logp.x, y=logp.y, label=Description)) +
  geom_point() + geom_smooth(method = "lm", color = "black") +
  geom_label_repel(data = Up_all_BP_complete_sorted[1:10,], size = 4) +
  theme_bw() +
  labs(x="-Log(q-value) RNA GO Terms", y="-Log(q-value) ATAC GO Terms") +
  theme(axis.text=element_text(size=14), axis.title=element_text(size=14),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 
ggsave("GO-RNA-ATACsub25_pvalComparison_BP.pdf", width = 8, height = 6)

# MF
RNA_Up_MF <- 
  subset(ego_genesUp_MF@result, qvalue < 0.05)[,c("ID", "Description", "GeneRatio", "qvalue")]
RNA_Up_MF$logp <- -log(RNA_Up_MF$qvalue)
ATAC_Up_MF <- 
  subset(go_I_MF@result, qvalue < 0.05)[,c("ID", "Description", "GeneRatio", "qvalue")]
ATAC_Up_MF$logp <- -log(ATAC_Up_MF$qvalue)

Up_all_MF <- merge(x=RNA_Up_MF, y=ATAC_Up_MF, by="Description", all.x = T, all.y = T)
Up_all_MF_complete <- Up_all_MF[complete.cases(Up_all_MF),]
Up_all_MF_complete$sort <- Up_all_MF_complete$logp.x + Up_all_MF_complete$logp.y
Up_all_MF_complete_sorted <- Up_all_MF_complete[order(-Up_all_MF_complete$sort),]

ggplot(Up_all_MF_complete_sorted, aes(x=logp.x, y=logp.y, label=Description)) +
  geom_point() + geom_smooth(method = "lm", color = "black") +
  geom_label_repel(data = Up_all_MF_complete_sorted[1:10,], size = 4) +
  theme_bw() +
  labs(x="-Log(q-value) RNA GO Terms", y="-Log(q-value) ATAC GO Terms") +
  theme(axis.text=element_text(size=14), axis.title=element_text(size=14),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 
ggsave("GO-RNA-ATACsub25_pvalComparison_MF.pdf", width = 8, height = 6)

# CC
RNA_Up_CC <- 
  subset(ego_genesUp_CC@result, qvalue < 0.05)[,c("ID", "Description", "GeneRatio", "qvalue")]
RNA_Up_CC$logp <- -log(RNA_Up_CC$qvalue)
ATAC_Up_CC <- 
  subset(go_I_CC@result, qvalue < 0.05)[,c("ID", "Description", "GeneRatio", "qvalue")]
ATAC_Up_CC$logp <- -log(ATAC_Up_CC$qvalue)

Up_all_CC <- merge(x=RNA_Up_CC, y=ATAC_Up_CC, by="Description", all.x = T, all.y = T)
Up_all_CC_complete <- Up_all_CC[complete.cases(Up_all_CC),]
Up_all_CC_complete$sort <- Up_all_CC_complete$logp.x + Up_all_CC_complete$logp.y
Up_all_CC_complete_sorted <- Up_all_CC_complete[order(-Up_all_CC_complete$sort),]

ggplot(Up_all_CC_complete_sorted, aes(x=logp.x, y=logp.y, label=Description)) +
  geom_point() + geom_smooth(method = "lm", color = "black") +
  geom_label_repel(data = Up_all_CC_complete_sorted[1:10,], size = 4) +
  theme_bw() +
  labs(x="-Log(q-value) RNA GO Terms", y="-Log(q-value) ATAC GO Terms") +
  theme(axis.text=element_text(size=14), axis.title=element_text(size=14),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 
ggsave("GO-RNA-ATACsub25_pvalComparison_CC.pdf", width = 8, height = 6)

#### PUT ALL TOGETHER ####
Up_all_BP_complete_sorted$GOtype <- "BP"
Up_all_MF_complete_sorted$GOtype <- "MF"
Up_all_CC_complete_sorted$GOtype <- "CC"

Up_all_GO <- rbind(Up_all_BP_complete_sorted, Up_all_MF_complete_sorted, 
                   Up_all_CC_complete_sorted)
Up_all_GO$GOtype <- factor(Up_all_GO$GOtype, levels = c("BP", "MF", "CC"))

# ggscatter(Up_all_GO, 
#           x = "logp.x", y = "logp.y", color = "GOtype",
#           add = "reg.line", size = 1) + #, label = "Term", repel = "True") +
#   stat_cor(method = "spearman", aes(label = ..r.label..), size = 7) +
#   ggrepel::geom_label_repel(data = subset(Up_all_GO, logp.x > 15 | logp.y > 30),
#                             aes(label = Description), size = 5) +
#   facet_wrap(~GOtype, ncol = 1, scales = "free_x") +
#   theme_bw() +
#   scale_color_manual(values = c("grey70", "grey40", "grey10")) +
#   labs(x=expression("-"~Log[10]~"(p)"~GO~Terms[Transcriptomics]),
#        y=expression("-"~Log[10]~"(p)"~GO~Terms[Epigenomics])) +
#   theme(axis.text=element_text(size=25), axis.title=element_text(size=25),
#         strip.text.x = element_text(size=25), strip.text.y = element_text(size=25),
#         panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#         legend.position = "none") + 
#   ggsave("RNA-ATACsub25_allGOtypes.svg", width = 12, height = 14)

########### FIGURE S3F ###########
ggscatter(Up_all_BP_complete_sorted, x = "logp.x", y = "logp.y", 
          add = "reg.line", size = 1) +
  geom_smooth(method = "lm", color = "black") +
  stat_cor(method = "spearman", aes(label = after_stat(r.label)), size = 6,
           label.x.npc=0.87, label.y.npc=0.95) +
  ggrepel::geom_label_repel(data = Up_all_BP_complete_sorted[1:10,],
                            aes(label = Description), size = 4) +
  theme_bw() +
  labs(x=expression("-"~Log[10]~"(p)"~GO~Terms[Transcriptomics]),
       y=expression("-"~Log[10]~"(p)"~GO~Terms[Epigenomics])) +
  theme(axis.text=element_text(size=20), axis.title=element_text(size=20),
        strip.text.x = element_text(size=20), strip.text.y = element_text(size=20),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        legend.position = "none") 
ggsave("RNA-ATACsub25_BP.svg", width = 7, height = 5.5)
#################################

########### FIGURE 3D ###########
ggscatter(Up_all_MF_complete_sorted, x = "logp.x", y = "logp.y",
          add = "reg.line", size = 1) +
  geom_smooth(method = "lm", color = "black") +
  stat_cor(method = "spearman", aes(label = after_stat(r.label)), size = 6,
           label.x.npc=0, label.y.npc=1) +
  ggrepel::geom_label_repel(data = Up_all_MF_complete_sorted[1:10,],
                            aes(label = Description), size = 4) +
  theme_bw() +
  labs(x=expression("-"~Log[10]~"(p)"~GO~Terms[Transcriptomics]),
       y=expression("-"~Log[10]~"(p)"~GO~Terms[Epigenomics])) +
  theme(axis.text=element_text(size=20), axis.title=element_text(size=20),
        strip.text.x = element_text(size=20), strip.text.y = element_text(size=20),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        legend.position = "none")
ggsave("RNA-ATACsub25_MF.svg", width = 7, height = 5.5)
#################################

########### FIGURE S3G ###########
ggscatter(Up_all_CC_complete_sorted, x = "logp.x", y = "logp.y",
          add = "reg.line", size = 1) +
  geom_smooth(method = "lm", color = "black") +
  stat_cor(method = "spearman", aes(label = after_stat(r.label)), size = 6,
           label.x.npc=0.1, label.y.npc=0.95) +
  ggrepel::geom_label_repel(data = Up_all_CC_complete_sorted[1:10,],
                            aes(label = Description), size = 4) +
  theme_bw() +
  labs(x=expression("-"~Log[10]~"(p)"~GO~Terms[Transcriptomics]),
       y=expression("-"~Log[10]~"(p)"~GO~Terms[Epigenomics])) +
  theme(axis.text=element_text(size=20), axis.title=element_text(size=20),
        strip.text.x = element_text(size=20), strip.text.y = element_text(size=20),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        legend.position = "none") 
ggsave("RNA-ATACsub25_CC.svg", width = 7, height = 5.5)
#################################
