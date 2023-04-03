library(reshape2)
library(ggplot2)
library(scales)

load("../data/merged_countdata.RData")
Ferr = read.delim("../data/KEGGFerroptosis_hsa04216_06-25-18.txt", header=T, stringsAsFactors = F)
Ferr = Ferr[rowSums(is.na(Ferr)) == 0, ]

countdata_Fer <- countdata[Ferr$GeneName,]

Fer_match = colsplit(colnames(countdata_Fer), pattern = "_", names = c("Population", "Time", "Replicate"))
Fer_plot = cbind(Fer_match, t(countdata_Fer))
Fer_melt = melt(data = Fer_plot, id.vars = c("Population", "Time", "Replicate"), measure.vars = Ferr$GeneName)
Fer_dat = summarySE(Fer_melt, measurevar = "value", groupvars = c("Population", "Time", "variable"))
Fer_dat$Time = as.numeric(gsub("[^[:digit:]]","",Fer_dat$Time))
# Fer_dat_sub = subset(Fer_dat, variable %in% c("ITPR2","ITPR3"))
Fer_ggploted <- ggplot(Fer_dat, aes(x=Time, y=value, group = interaction(variable, Population))) + 
  geom_line(size=1.5, aes(color = Population)) + 
  geom_point(size = 1.5, aes(color = Population)) + facet_wrap(~variable, ncol = 5, scales = "free") +
  geom_errorbar(aes(ymin=value-sd, ymax=value+sd, color = Population), width=.2, size=1.5) +
  theme_bw() + xlab("Time (days)") + ylab("Gene Counts") +
  ggtitle("Ferroptosis gene signature") +
  theme(legend.text = element_text(size = 10), legend.position = "right", 
        plot.title = element_text(size = 14, hjust = 0.5, face = "bold"), axis.text=element_text(size=12),
        legend.title = element_text(size=12,face="bold"),
        axis.title=element_text(size=12,face="bold"))

Fer_ggploted + ggsave("FerroptosisGeneSignature_rawCounts_SKMEL5sublines+treatment.pdf", width = 20, height = 25)

load("../data/RLD_SC-1,7,10_0,3,8d_20180701.RData")
normdata_Fer <- as.data.frame(assay(rld))[Ferr$GeneName,]

Fer_match = colsplit(colnames(normdata_Fer), pattern = "_", names = c("Population", "Time", "Replicate"))
Fer_plot = cbind(Fer_match, t(normdata_Fer))
Fer_melt = melt(data = Fer_plot, id.vars = c("Population", "Time", "Replicate"), measure.vars = unique(colnames(Fer_plot))[4:42])
Fer_dat = summarySE(Fer_melt, measurevar = "value", groupvars = c("Population", "Time", "variable"))
Fer_dat$Time = as.numeric(gsub("[^[:digit:]]","",Fer_dat$Time))
# Fer_dat_sub = subset(Fer_dat, variable %in% c("ITPR2","ITPR3"))
# Fer_ggploted <- ggplot(Fer_dat, aes(x=Time, y=value, group = interaction(variable, Population))) + 
#   geom_line(size=1.5, aes(color = Population)) + 
#   geom_point(size = 1.5, aes(color = Population)) + facet_wrap(~variable, ncol = 5, scales = "free") +
#   geom_errorbar(aes(ymin=value-sd, ymax=value+sd, color = Population), width=.2, size=1.5) +
#   theme_bw() + xlab("Time (days)") + ylab("Gene Counts") +
#   ggtitle("Ferroptosis gene signature") +
#   theme(legend.text = element_text(size = 10), legend.position = "right", 
#         plot.title = element_text(size = 14, hjust = 0.5, face = "bold"), axis.text=element_text(size=12),
#         legend.title = element_text(size=12,face="bold"),
#         axis.title=element_text(size=12,face="bold"))
# 
# Fer_ggploted + ggsave("FerroptosisGeneSignature_normCounts_SKMEL5sublines+treatment.pdf", width = 20, height = 25)


####

library(pheatmap)
# normdata_Fer
samples <- c("SC01_day0", "SC01_day3", "SC01_day8", "SC07_day0", "SC07_day3", "SC07_day8", "SC10_day0", "SC10_day3", "SC10_day8")
test <- normdata_Fer
test1 <- sapply(samples, function(x) rowMeans(test[, grep(x, colnames(test))]))
test2 <- as.data.frame(test1[complete.cases(test1),])

####
# Plot pheatmap (z-score) of gene list 
test3 <- test2
test3$Gene <- rownames(test2)
test3_sub <- subset(test3, !Gene %in% c("ALOX15", "ACSL1", "ACSL3",
                                        "ACSL5", "ACSL6", "TP53", "TF",
                                         "CP", "MAP1LC3A", "MAP1LC3C",
                                         "CYBB"))
test3_sub <- test3_sub[,1:9]
pheatmap(test3_sub,cluster_cols=FALSE, cluster_rows = F, scale = "row")

test4_special <- subset(test3, Gene %in% c("SLC3A2", "SLC7A11", "GCLC", "GSS", "GPX4",
                                           "ACSL1", "ACSL3", "ACSL4", "LPCAT3",
                                           "TFRC", "STEAP3", "SLC11A2", "SLC39A8", "SLC39A14", 
                                           "PCBP2", "SLC40A1", 
                                           "HMOX1",
                                           "PRNP",
                                           "PCBP1", "FTH1", "MAP1LC3B", "ATG5", "ATG7", "NCOA4",
                                           "VDAC2", "VDAC3"))
test4_special <- test4_special[,1:9]
pheatmap(test4_special,cluster_cols=FALSE, cluster_rows = F, scale = "row")




allFC <- function(DEProc,startcol,endcol){ 
  GE_fold = DEProc[,-c(startcol:endcol)]
  colvec = colnames(DEProc)[startcol:endcol]
  
  #Last index is a self comparison and is removed
  for(k in 1:(length(colvec)-1)){
    #Start with column that is 1 away from index 
    for(j in (k+1):length(colvec)){
      compnam = paste0(colvec[j],"/",colvec[k])
      #Loop through each gene/row  
      for(i in 1:nrow(DEProc)){
        f = DEProc[i,colvec[j]]
        h = DEProc[i,colvec[k]]
        GE_fold[i, compnam] = log2(f/h)

        #Capture upregulation and down regulation
        # if(f>h){
        #   # GE_fold[i,compnam] = 2^(f-h)
        #   GE_fold[i,compnam] = log2(f/h)
        # }else{
        #   # GE_fold[i,compnam] = -2^(h-f)
        #   GE_fold[i,compnam] = log2(f/h)
      #   }
      #   
     }
    }
  }
  
  return(GE_fold)
  
}

GE_fold <- allFC(test2, 1,9)
ImpRat = c("SC01_day3/SC01_day0", "SC01_day8/SC01_day0", 
           "SC07_day3/SC07_day0", "SC07_day8/SC07_day0", 
           "SC10_day3/SC10_day0", "SC10_day8/SC10_day0")
Imp_fold = GE_fold[,ImpRat]

Ferr = read.delim("~/Documents/QuarantaLab/Ferroptosis/RNAseq/KEGGFerroptosis_hsa04216_06-25-18.txt", header=T, stringsAsFactors = F)
Ferr_fold <- Imp_fold[rownames(Imp_fold) %in% Ferr$GeneName,]
colnames(Ferr_fold) <- c("SC01_short", "SC01_long", "SC07_short", "SC07_long", "SC10_short", "SC10_long")
Ferr_fold$SC01_base <- 0
Ferr_fold$SC07_base <- 0
Ferr_fold$SC10_base <- 0
Ferr_fold <- Ferr_fold[,c(7,1,2,8,3,4,9,5,6)]
Ferr_fold$Gene <- rownames(Ferr_fold)
Ferr_fold_melt <- melt(Ferr_fold, id.vars = "Gene")
Ferr_fold_melt$Gene <- factor(Ferr_fold_melt$Gene, levels = rev(Ferr$GeneName))

Ferr_fold_melt <- subset(Ferr_fold_melt, !Gene %in% c("ALOX15", "ACSL1", "ACSL3",
                                                     "ACSL5", "ACSL6", "TP53", "TF",
                                                     "CP", "MAP1LC3A", "MAP1LC3C",
                                                     "CYBB"))
ggplot(Ferr_fold_melt, aes(variable, Gene, fill = value)) + 
  geom_tile(color = "black") + theme_bw() +
  scale_fill_gradientn(
    colors=c("blue","white","red","red4"),
    values=rescale(c(min(Ferr_fold_melt$value), 0, max(Ferr_fold_melt$value)/2, max(Ferr_fold_melt$value))),
    limits=c(min(Ferr_fold_melt$value),max(Ferr_fold_melt$value)),
    name = "Log2 Fold Change"
  ) + xlab("Population") +
  theme(axis.text=element_text(size=14),
        axis.text.x=element_text(angle = 90, hjust = 0),
        axis.title=element_text(size=14),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  ggsave("Ferroptosis_FCHM_selected_wide.pdf", width = 10, height = 6)
