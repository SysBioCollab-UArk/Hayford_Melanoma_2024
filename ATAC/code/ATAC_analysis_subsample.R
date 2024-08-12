library(Rsamtools)
library(BSgenome.Hsapiens.UCSC.hg38)
library(ggrepel)
library(ChIPseeker)
library(soGGi)
library(MotifDb)
library(Biostrings)
library(BSgenome.Hsapiens.UCSC.hg38)
library(GenomicAlignments)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

library(ChIPQC)
library(rtracklayer)
library(DT)
library(dplyr)
library(tidyr)

# if (Sys.getenv("RSTUDIO") == "1") {
#   setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
# }

blkList <- import.bed(file.path("..", "data", "ENCFF356LFX.bed.gz"))
chrs <- c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8",
          "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16",
          "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX", "chrY")

openRegionPeaks_UT <- file.path("..", "data", "untreated_all_peaks.narrowPeak_peaks.narrowPeak")
qcRes_UT <- ChIPQCsample(
  file.path("..", "data", "bam", "trimmed_3334-CH-1-TAGGCATG-ACTGCATA_S121_aligned_sorted_dedup_unique_fullClean.bam"),
  peaks = openRegionPeaks_UT, annotation ="hg38", chromosomes = chrs, 
  blacklist = blkList, verboseT = FALSE)

openRegionPeaks_I <- file.path("..", "data", "idling_peaks.narrowPeak")
qcRes_I <- ChIPQCsample(
  file.path("..", "data", "bam", "idling_subsample_25p_dedup_unique_fullClean.bam"), 
  peaks = openRegionPeaks_I, annotation ="hg38", chromosomes = chrs, 
  blacklist = blkList, verboseT = FALSE)

# save(qcRes_UT, qcRes_I, file = "allPeaks_narrow_idlingSub25.RData")

### Peak Annotation
MacsCalls_UT <- granges(qcRes_UT)
`%over%` <- IRanges::`%over%` # Added this to avoid a conflict with %over% defined in another package --LAH
data.frame(Blacklisted=sum(MacsCalls_UT %over% blkList),
           Not_Blacklisted=sum(!MacsCalls_UT %over% blkList))
MacsCalls_UT <- MacsCalls_UT[!MacsCalls_UT %over% blkList]

MacsCalls_I <- granges(qcRes_I)
data.frame(Blacklisted=sum(MacsCalls_I %over% blkList),
           Not_Blacklisted=sum(!MacsCalls_I %over% blkList))
MacsCalls_I <- MacsCalls_I[!MacsCalls_I %over% blkList]

promoter <- getPromoters(TxDb=txdb, upstream=3000, downstream=3000)
tagMatrix_UT <- getTagMatrix(MacsCalls_UT, windows=promoter)
tagMatrix_I <- getTagMatrix(MacsCalls_I, windows=promoter)

plotAvgProf(tagMatrix_UT, xlim=c(-3000, 3000),
            xlab="Genomic Region (5'->3')", ylab = "Read Count Frequency")
plotAvgProf(tagMatrix_I, xlim=c(-3000, 3000),
            xlab="Genomic Region (5'->3')", ylab = "Read Count Frequency")

MacsCalls_UT_Anno <-  annotatePeak(MacsCalls_UT, 
                                   TxDb = TxDb.Hsapiens.UCSC.hg38.knownGene)
MacsCalls_I_Anno <-  annotatePeak(MacsCalls_I, 
                                  TxDb = TxDb.Hsapiens.UCSC.hg38.knownGene)

plotAnnoPie(MacsCalls_UT_Anno)
plotAnnoPie(MacsCalls_I_Anno)

vennpie(MacsCalls_UT_Anno)
vennpie(MacsCalls_I_Anno)

library(ggupset)

upsetplot(MacsCalls_UT_Anno)
upsetplot(MacsCalls_I_Anno)

# Limma venn diagram
peaks <- dir("../data/venn/", pattern = "*.narrowPeak", 
             full.names = TRUE)
myPeaks <- lapply(peaks, ChIPQC:::GetGRanges, simple = TRUE)

names(myPeaks) <- c("Idling", "Untreated")
Group <- factor(c("Idling", "Untreated"))

# library(devtools)
# install_github("https://github.com/ColeWunderlich/soGGi.git")

source("soGGi_runConsensusRegions_fixed.R")
consensusToCount <- runConsensusRegions(GRangesList(myPeaks), "none")
# consensusToCount

library(limma)

as.data.frame(elementMetadata(consensusToCount)) %>% dplyr::select(Untreated, Idling) %>% 
  vennDiagram(main = "Overlap for ATAC Open regions")

vd <- as.data.frame(elementMetadata(consensusToCount))[,c("Untreated", "Idling")]
vd_counts <- plyr::count(vd)

library(eulerr)
library(UpSetR)

shared_peaks <- c("Untreated" = vd_counts$freq[2],
                  "Idling" = vd_counts$freq[1],
                  "Untreated&Idling" = vd_counts$freq[3])

########### FIGURE S3C ###########
pdf('venn_sharedPeaks.pdf', width=4, height=4)
venn_sharedPeaks <- euler(shared_peaks)
p <- plot(venn_sharedPeaks, fills = c("red", "lightblue"), shape = "ellipse", 
          quantities = TRUE)
plot(p)
dev.off()
##################################

upset(fromExpression(shared_peaks), order.by = "freq",
      sets.bar.color = c("blue", "red"), 
      mb.ratio = c(0.55, 0.45), nintersects = NA,
      show.numbers = "no", point.size = 5, line.size = 0.5, mainbar.y.max = 120000,
      mainbar.y.label = "Peak Intersections", sets.x.label = "Set Size", 
      text.scale = c(2, 1.5, 2, 1.5, 2, 1))

# Annotate unique regions
library(clusterProfiler)
library(ChIPseeker)

CTC_UT <- consensusToCount[(elementMetadata(consensusToCount)[,"Untreated"] == 1) &
                             (elementMetadata(consensusToCount)[,"Idling"] == 0)]
CTC_I <- consensusToCount[(elementMetadata(consensusToCount)[,"Untreated"] == 0) &
                            (elementMetadata(consensusToCount)[,"Idling"] == 1)]
CTC_con <- consensusToCount[(elementMetadata(consensusToCount)[,"Untreated"] == 1) &
                              (elementMetadata(consensusToCount)[,"Idling"] == 1)]


anno_UT <- annotatePeak(CTC_UT, TxDb = TxDb.Hsapiens.UCSC.hg38.knownGene)
anno_I <- annotatePeak(CTC_I, TxDb = TxDb.Hsapiens.UCSC.hg38.knownGene)
anno_con <- annotatePeak(CTC_con, TxDb = TxDb.Hsapiens.UCSC.hg38.knownGene)

# save(anno_UT, anno_I, anno_con, 
#      file = "SKMEL5_ATACsub25_annotatedPeaks_uniqueShared.RData")

##

# Unique - submit to great, GO, cluster profiler
library(org.Hs.eg.db)
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

# dotplot(go_UT_BP) + ggsave("GOenrichment_sub25_UT_BP.pdf", width = 8, height = 5)
# dotplot(go_I_BP)
# ggsave("GOenrichment_sub25_I_BP.pdf") #, width = 8, height = 5)
# dotplot(go_con_BP) + ggsave("GOenrichment_sub25_con_BP.pdf", width = 8, height = 5)
# 
# dotplot(go_UT_MF) + ggsave("GOenrichment_sub25_UT_MF.pdf", width = 8, height = 5)

########### FIGURE 3C ###########
library(scales)
dotplot(go_I_MF, font.size = 14, label_format = 40) + 
  labs(x="Gene Ratio") +
  theme(legend.text = element_text(size = 12),
        plot.title = element_text(size = 14, hjust = 0.5, face = "bold"), 
        axis.text=element_text(size=6),
        legend.title = element_text(size=12, face="bold"), 
        legend.position = c(0.8, 0.32), legend.box = 'vertical',
        axis.title=element_text(size=14, face="bold")) +
  scale_y_discrete(labels = label_wrap(35))
ggsave("GOenrichment_sub25_I_MF.pdf", width = 6, height = 6)
#################################

# ggsave("GOenrichment_sub25_I_MF.svg") #, width = 8, height = 5)
# dotplot(go_con_MF) + ggsave("GOenrichment_sub25_con_MF.pdf", width = 8, height = 5)
# 
# dotplot(go_UT_CC) + ggsave("GOenrichment_sub25_UT_CC.pdf", width = 8, height = 5)
# dotplot(go_I_CC)
# ggsave("GOenrichment_sub25_I_CC.pdf") #, width = 8, height = 5)
# dotplot(go_con_CC) + ggsave("GOenrichment_sub25_con_CC.pdf", width = 8, height = 5)

df <- data.frame(UT_MF = go_UT_MF@result$Description[1:25],
                 I_MF = go_I_MF@result$Description[1:25],
                 con_MF = go_con_MF@result$Description[1:25])

files <- list(CTC_UT, CTC_I, CTC_con)
names(files) <- c("Untreated", "Idling", "Shared")
peakAnnoList <- lapply(files, annotatePeak, TxDb=txdb,
                       tssRegion=c(-3000, 3000), verbose=FALSE)

########### FIGURE S3E ###########
plotAnnoBar(peakAnnoList) +
  guides(fill=guide_legend(ncol=2)) +
  theme(plot.title = element_blank(), #element_text(size = 12, hjust = 0.5), 
        plot.margin = unit(c(0.1,1,0.1,-0.5), 'cm'), # (t,r,b,l)
        axis.title=element_text(size = 20),
        axis.text=element_text(size = 20),
        legend.text = element_text(size = 16), legend.position = "bottom", 
        legend.title = element_blank(), #element_text(size= 16), 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())
ggsave("ATAC_sub25_annotationDistribution_UniqueShared.pdf", width = 6, height = 8)
##################################

########### FIGURE S3D ###########
plotDistToTSS(peakAnnoList) +
  theme(plot.title = element_blank(), #element_text(size = 12, hjust = 0.5), 
        plot.margin = unit(c(0.1,1,0.1,-0.5), 'cm'), # (t,r,b,l)
        axis.title=element_text(size=24),
        axis.text.x=element_text(size=20),
        axis.text.y=element_text(size=24),
        legend.text = element_text(size = 22), legend.position = "bottom", 
        legend.title = element_blank(), #element_text(size=12), 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 
ggsave("ATAC_sub25_distanceToTSS_UniqueShared.pdf", width = 7.5, height = 8)
##################################

tagMatrixList <- lapply(files, getTagMatrix, windows=promoter)
plotAvgProf(tagMatrixList, xlim=c(-3000, 3000)) +
  theme(legend.text = element_text(size = 12), legend.position = "right",
        plot.title = element_text(size = 14, hjust = 0.5), axis.text=element_text(size=12),
        legend.title = element_text(size=12), axis.title=element_text(size=12),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())
ggsave("ATAC_sub25_averageBindingProfile_UniqueShared.pdf") #, width = 6, height = 4)
