library(reshape2)
library(ggplot2)

if (Sys.getenv("RSTUDIO") == "1") {
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
}

bcDat <- read.csv('../data/bcPlot.csv')
bcNum <- read.csv('../data/bcCount_allTechReps.csv')

numReads <- data.frame(Untreated = c(5340832,5750058,6238661),
                       Idling = c(6153559,5634086,6160163))
numReads_CPM_SF <- data.frame(apply(numReads, 2, function(x){x/1000000}))

numReads_all <- unlist(numReads_CPM_SF)

names(bcNum) <- c("Barcode", "U1R1", "U2R1", "U2R2", "T1R1", "T2R1", "T2R2")
bcNames <- unlist(bcNum[c(1)])
bcNum1 <- bcNum[c(2:7)]
bcNum1 <- sweep(bcNum1, 2, numReads_all, FUN = '/')
bcNum1$Barcode <- bcNames
bcNum1 <- bcNum1[,c(7,1,2,3,4,5,6)]
bcNum <- bcNum1

bcNum_UT <- bcNum[,c(1:4)]
names(bcNum_UT) <- c("Barcode", "Rep1", "Rep2", "Rep3")
bcNum_I <- bcNum[,c(1, 5:7)]
names(bcNum_I) <- c("Barcode", "Rep1", "Rep2", "Rep3")

# Rank order line plot
bcNum_UT_order <- bcNum_UT[order(bcNum_UT$Rep1, decreasing=TRUE),]
bcNum_UT_melt <- melt(bcNum_UT_order, id = "Barcode")
bcNum_UT_melt$Barcode <- factor(bcNum_UT_melt$Barcode,
                                levels=as.character(bcNum_UT_order$Barcode))
bcNum_UT_melt$BarcodeNumber <- rep(seq(nrow(subset(bcNum_UT_melt, variable == "Rep1"))),
                                   times = 3)
bcNum_UT_melt_sub <- subset(bcNum_UT_melt, BarcodeNumber %in% c(1:50))
# bcNum_UT_melt_sub <- subset(bcNum_UT_melt, value > 1000)
# ggplot(bcNum_UT_melt_sub, aes(x=Barcode, y=value, group = variable,
#                           color = variable)) + geom_path() + theme_bw() +
#   labs(x = "Barcode (rank order)", y = "Reads Per Million (RPM)") +
#   theme(axis.text.x=element_blank(),
#         axis.ticks.x=element_blank(), axis.text.y = element_text(size = 14),
#         legend.position = "right", legend.text = element_text(size = 14),
#         plot.title = element_text(size = 16, hjust = 0.5, face = "bold"),
#         axis.text=element_text(size=14),
#         legend.title = element_text(size=14), axis.title=element_text(size=14),
#         panel.grid.major = element_blank(), panel.grid.minor = element_blank())
#   ggsave("SKMEL5_barcodeUT_rankAbundance.pdf", width = 8, height = 5)


bcNum_I_order <- bcNum_I[order(bcNum_I$Rep1, decreasing=TRUE),]
bcNum_I_melt <- melt(bcNum_I_order, id = "Barcode")
bcNum_I_melt$Barcode <- factor(bcNum_I_melt$Barcode,
                                levels=as.character(bcNum_I_order$Barcode))
bcNum_I_melt$BarcodeNumber <- rep(seq(nrow(subset(bcNum_I_melt, variable == "Rep1"))),
                                   times = 3)
bcNum_I_melt_sub <- subset(bcNum_I_melt, BarcodeNumber %in% c(1:50))
# bcNum_I_melt_sub <- subset(bcNum_I_melt, value > 1000)
# ggplot(bcNum_I_melt_sub, aes(x=Barcode, y=value, group = variable,
#                               color = variable)) + geom_path() + theme_bw() +
#   labs(x = "Barcode (rank order)", y = "Reads Per Million (RPM)") +
#   theme(axis.text.x=element_blank(),
#         axis.ticks.x=element_blank(), axis.text.y = element_text(size = 14),
#         legend.position = "right", legend.text = element_text(size = 14),
#         plot.title = element_text(size = 16, hjust = 0.5, face = "bold"),
#         axis.text=element_text(size=14),
#         legend.title = element_text(size=14), axis.title=element_text(size=14),
#         panel.grid.major = element_blank(), panel.grid.minor = element_blank())
#   ggsave("SKMEL5_barcodeI_rankAbundance.pdf", width = 8, height = 5)



##
UT_R1_totalBC <- subset(bcNum_UT[c(1,2)], Rep1 > 100)
UT_R2_totalBC <- subset(bcNum_UT[c(1,3)], Rep2 > 100)
UT_R3_totalBC <- subset(bcNum_UT[c(1,4)], Rep3 > 100)

I_R1_totalBC <- subset(bcNum_I[c(1,2)], Rep1 > 100)
I_R2_totalBC <- subset(bcNum_I[c(1,3)], Rep2 > 100)
I_R3_totalBC <- subset(bcNum_I[c(1,4)], Rep3 > 100)

shared_R1_totalBC <- sum(I_R1_totalBC$Barcode %in% UT_R1_totalBC$Barcode)
shared_R2_totalBC <- sum(I_R2_totalBC$Barcode %in% UT_R2_totalBC$Barcode)
shared_R3_totalBC <- sum(I_R3_totalBC$Barcode %in% UT_R3_totalBC$Barcode)

bcCount <- data.frame(Untreated = c(nrow(UT_R1_totalBC),
                                    nrow(UT_R2_totalBC),
                                    nrow(UT_R3_totalBC)),
                      Idling = c(nrow(I_R1_totalBC),
                                 nrow(I_R2_totalBC),
                                 nrow(I_R3_totalBC)),
                      Shared = c(shared_R1_totalBC,
                                 shared_R2_totalBC,
                                 shared_R3_totalBC))
bcCount_melt <- melt(bcCount)
bcCount_melt.summary <- aggregate(bcCount_melt, by=list(bcCount_melt$variable), FUN=mean)[c(1,3)]
names(bcCount_melt.summary) <- c("variable", "value")

g1 <- ggplot(bcCount_melt, aes(x=variable,y=value,color=variable)) +
  geom_jitter(width = 0.2)+ theme_bw() +
  geom_crossbar(data = bcCount_melt.summary, aes(ymin=value, ymax=value, color = variable),
                linewidth=0.5, width = 0.5) +
  ylim(0,500) + scale_color_manual(values = c("red", "blue", "green")) +
  labs(x = "Condition", y = "Number of Unique Barcodes") +
  theme(axis.text.y = element_text(size = 14),
        legend.position = "none", legend.text = element_text(size = 14),
        plot.title = element_text(size = 16, hjust = 0.5, face = "bold"),
        axis.text=element_text(size=14),
        legend.title = element_text(size=14), axis.title=element_text(size=14),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  ylim(300, 500)

ggsave("SKMEL5_numUniqueBC_byCondition.pdf") #, width = 4, height = 5)

# ggplot(subset(bcCount_melt, variable %in% c("Untreated", "Idling")),
#               aes(x=variable,y=value,color=variable)) +
#   geom_jitter(width = 0.2)+ theme_bw() +
#   geom_crossbar(data = subset(bcCount_melt.summary, variable %in% c("Untreated", "Idling")),
#                 aes(ymin=value, ymax=value, color = variable),
#                 size=0.5, width = 0.5) +
#   ylim(0,500) + scale_color_manual(values = c("red", "blue")) +
#   labs(x = "Condition", y = "Number of Unique Barcodes") +
#   theme(axis.text.y = element_text(size = 12),
#         legend.position = "none", legend.text = element_text(size = 12),
#         plot.title = element_text(size = 16, hjust = 0.5, face = "bold"),
#         axis.text=element_text(size=12),
#         legend.title = element_text(size=12), axis.title=element_text(size=12),
#         panel.grid.major = element_blank(), panel.grid.minor = element_blank())
#   ggsave("SKMEL5_numUniqueBC_byCondition_noShared.pdf", width = 2.5, height = 5)

###
# Num barcodes shared between replicates
## Untreated
UT_bcShareList <- list(UT_R1_totalBC$Barcode,
                       UT_R2_totalBC$Barcode,
                       UT_R3_totalBC$Barcode)

UT_all3 <- Reduce(intersect, UT_bcShareList)
UT_onlyR1 <- setdiff(UT_bcShareList[[1]],
                     unique(unlist(list(UT_bcShareList[[2]],
                                        UT_bcShareList[[3]]))))
UT_onlyR2 <- setdiff(UT_bcShareList[[2]],
                     unique(unlist(list(UT_bcShareList[[1]],
                                        UT_bcShareList[[3]]))))
UT_onlyR3 <- setdiff(UT_bcShareList[[3]],
                     unique(unlist(list(UT_bcShareList[[1]],
                                        UT_bcShareList[[2]]))))

UT_R1R2 <- intersect(UT_bcShareList[[1]],UT_bcShareList[[2]])
UT_R1R2_t <- UT_R1R2[which(!(UT_R1R2 %in% unique(unlist(list(UT_all3, UT_onlyR1, UT_onlyR2)))))]

UT_R1R3 <- intersect(UT_bcShareList[[1]],UT_bcShareList[[3]])
UT_R1R3_t <- UT_R1R3[which(!(UT_R1R3 %in% unique(unlist(list(UT_all3, UT_onlyR1, UT_onlyR3)))))]

UT_R2R3 <- intersect(UT_bcShareList[[2]],UT_bcShareList[[3]])
UT_R2R3_t <- UT_R2R3[which(!(UT_R2R3 %in% unique(unlist(list(UT_all3, UT_onlyR2, UT_onlyR3)))))]

## Idling
I_bcShareList <- list(I_R1_totalBC$Barcode,
                       I_R2_totalBC$Barcode,
                       I_R3_totalBC$Barcode)

I_all3 <- Reduce(intersect, I_bcShareList)
I_onlyR1 <- setdiff(I_bcShareList[[1]],
                     unique(unlist(list(I_bcShareList[[2]],
                                        I_bcShareList[[3]]))))
I_onlyR2 <- setdiff(I_bcShareList[[2]],
                     unique(unlist(list(I_bcShareList[[1]],
                                        I_bcShareList[[3]]))))
I_onlyR3 <- setdiff(I_bcShareList[[3]],
                     unique(unlist(list(I_bcShareList[[1]],
                                        I_bcShareList[[2]]))))

I_R1R2 <- intersect(I_bcShareList[[1]],I_bcShareList[[2]])
I_R1R2_t <- I_R1R2[which(!(I_R1R2 %in% unique(unlist(list(I_all3, I_onlyR1, I_onlyR2)))))]

I_R1R3 <- intersect(I_bcShareList[[1]],I_bcShareList[[3]])
I_R1R3_t <- I_R1R3[which(!(I_R1R3 %in% unique(unlist(list(I_all3, I_onlyR1, I_onlyR3)))))]

I_R2R3 <- intersect(I_bcShareList[[2]],I_bcShareList[[3]])
I_R2R3_t <- I_R2R3[which(!(I_R2R3 %in% unique(unlist(list(I_all3, I_onlyR2, I_onlyR3)))))]

###
sharingBCs <- data.frame(Shared = c("3", "2", "1"),
                            Untreated = c(length(UT_all3),
                                          length(unique(unlist(list(UT_R1R2_t, UT_R1R3_t, UT_R2R3_t)))),
                                          length(unique(unlist(list(UT_onlyR1, UT_onlyR2, UT_onlyR3))))),
                            Idling = c(length(I_all3),
                                       length(unique(unlist(list(I_R1R2_t, I_R1R3_t, I_R2R3_t)))),
                                       length(unique(unlist(list(I_onlyR1, I_onlyR2, I_onlyR3))))))

library(reshape2)
library(dplyr)
sharingBCs_melt <- melt(sharingBCs, id.vars = "Shared",
                        measure.vars = c("Untreated", "Idling"))
sharingBCs_melt_prop <-sharingBCs_melt %>% group_by(variable) %>%
  mutate(prop = value/sum(value))

# ggplot(sharingBCs_melt_prop, aes(x=variable, y=prop, group = Shared, fill = Shared)) +
#   geom_bar(position = "fill", stat = "identity", color = "black") +
#   xlab("Condition") + theme_classic() +
#   scale_y_continuous(name = "Percent of Unique Barcodes", labels = scales::percent) +
#   theme(axis.text = element_text(size = 12), axis.title=element_text(size=12),
#         legend.position = "none", legend.text = element_text(size = 12),
#         legend.title = element_text(size=14),
#         panel.grid.major = element_blank(), panel.grid.minor = element_blank())
#   ggsave("SKMEL5_propShared_byCondition.pdf", width = 3, height = 5)

sharingBCs_withinRep <- data.frame(Shared = c("3", "2", "1"),
                                   UTR1 = c(length(UT_all3),
                                            length(unique(unlist(list(UT_R1R2_t, UT_R1R3_t)))),
                                            length(UT_onlyR1)),
                                   UTR2 = c(length(UT_all3),
                                            length(unique(unlist(list(UT_R1R2_t, UT_R2R3_t)))),
                                            length(UT_onlyR2)),
                                   UTR3 = c(length(UT_all3),
                                            length(unique(unlist(list(UT_R1R3_t, UT_R2R3_t)))),
                                            length(UT_onlyR3)),
                                   IR1 = c(length(I_all3),
                                            length(unique(unlist(list(I_R1R2_t, I_R1R3_t)))),
                                            length(I_onlyR1)),
                                   IR2 = c(length(I_all3),
                                            length(unique(unlist(list(I_R1R2_t, I_R2R3_t)))),
                                            length(I_onlyR2)),
                                   IR3 = c(length(I_all3),
                                            length(unique(unlist(list(I_R1R3_t, I_R2R3_t)))),
                                            length(I_onlyR3)))

sharingBCs_withinRep_melt <- melt(sharingBCs_withinRep, id.vars = "Shared",
                        measure.vars = c("UTR1","UTR2","UTR3",
                                         "IR1", "IR2", "IR3"))
sharingBCs_withinRep_melt_prop <-sharingBCs_withinRep_melt %>% group_by(variable) %>%
  mutate(prop = value/sum(value))

plot_propShared_byRep <- ggplot(sharingBCs_withinRep_melt_prop, aes(x=variable, y=prop, group = Shared, fill = Shared)) +
  geom_bar(position = "fill", stat = "identity", color = "black") +
  xlab("Condition") + theme_classic() +
  scale_y_continuous(name = "Percent of Unique Barcodes", labels = scales::percent) +
  theme(axis.text = element_text(size = 12), axis.title=element_text(size=12),
        legend.position = "right", legend.text = element_text(size = 12),
        legend.title = element_text(size=14),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())

plot_propShared_byRep 
ggsave("SKMEL5_propShared_byReplicate.pdf") #, width = 4, height = 5)
plot_propShared_byRep_leg <- ggpubr::get_legend(plot_propShared_byRep)

ggpubr::as_ggplot(plot_propShared_byRep_leg)
ggsave("SKMEL5_propShared_byReplicate_legend.pdf")



###
bcNum_order <- bcNum[order(bcNum$U1R1, decreasing=FALSE),]
bcNum_melt <- melt(bcNum_order, id = "Barcode")
bcNum_melt$Barcode <- factor(bcNum_melt$Barcode,
                             levels=as.character(bcNum_order$Barcode))

ggplot(bcNum_melt, aes(x=variable, y=Barcode, fill = value)) +
  geom_tile() + theme_bw() +
  scale_fill_gradient(name = expression(Log[10]~RPM), trans = "log10",
                      low = "grey90", high = "black") +
  labs(x = "Condition", y = "Barcode") +
  # scale_x_discrete(name = "Untreated                      Treated",
  #                  limits = c("Rep1", "Rep2", "Rep3",
  #                             "Rep1", "Rep2", "Rep3"))+
  theme(axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        legend.position = "right", legend.text = element_text(size = 14),
        plot.title = element_text(size = 16, hjust = 0.5, face = "bold"),
        axis.text=element_text(size=14),
        legend.title = element_text(size=14), axis.title=element_text(size=14),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())
ggsave("SKMEL5_barcode_RPM_rank.pdf") #, width = 4, height = 5)

#####

bcNum_UT_test <- bcNum_order[,c(1:4)]
names(bcNum_UT_test) <- c("Barcode", "Rep1", "Rep2", "Rep3")
bcNum_UT_test$Sample <- "Untreated"

bcNum_UT_testMeans <- bcNum_UT_test[,c(2:4)]
bcNum_UT_testMeans$mean <- rowMeans(bcNum_UT_testMeans)
bcNum_UT_testMeans$Barcode <- unlist(bcNum_UT_test[,c(1)])
bcNum_UT_testMeans_order <- bcNum_UT_testMeans[order(bcNum_UT_testMeans$mean, decreasing=TRUE),]

bcNum_I_test <- bcNum_order[,c(1,5:7)]
names(bcNum_I_test) <- c("Barcode", "Rep1", "Rep2", "Rep3")
bcNum_I_test$Sample <- "Idling"
bcNum_test <- rbind(bcNum_UT_test, bcNum_I_test)
bcNum_test_melt <- melt(bcNum_test, id = c("Barcode", "Sample"))
bcNum_test_melt$Barcode <- factor(bcNum_test_melt$Barcode,
                           levels=as.character(bcNum_UT_testMeans_order$Barcode))

# This is the same function as in the Rmisc package
library(Rmisc)
# source('~/Documents/QuarantaLab/SummarySE.R')

bcNum_compare <- summarySE(bcNum_test_melt, measurevar = "value",
                           groupvars = c("Barcode", "Sample"))

bcNum_compare_sub <- bcNum_compare[c(1:100),]
bcNum_compare_sub <- subset(bcNum_compare, value > 1000)
# ggplot(bcNum_compare_sub, aes(x=Barcode, y=value, group = Sample,
#                              color = Sample)) +
#   theme_bw() +
#   geom_line() + geom_ribbon(aes(ymin=value-sd, ymax=value+sd, fill = Sample),
#                             alpha=0.1) +
#   labs(x = "Barcode (rank order)", y = "Abundance") +
#   theme(axis.text.x=element_blank(),
#         axis.ticks.x=element_blank(), axis.text.y = element_text(size = 14),
#         legend.position = "right", legend.text = element_text(size = 14),
#         plot.title = element_text(size = 16, hjust = 0.5, face = "bold"),
#         axis.text=element_text(size=14),
#         legend.title = element_text(size=14), axis.title=element_text(size=14),
#         panel.grid.major = element_blank(), panel.grid.minor = element_blank())
#   ggsave("SKMEL5_barcode_rankAbundance_comparison.pdf", width = 8, height = 5)


#######

bcNum_prop <- data.frame(lapply(bcNum[,-1], function(x){x/sum(x, na.rm=TRUE)}))
bcNum_prop$Barcode <- bcNum$Barcode
bcNum_prop <- bcNum_prop[,c(7,1,2,3,4,5,6)]

bcNum_prop_UT_test <- bcNum_prop[,c(1:4)]
names(bcNum_prop_UT_test) <- c("Barcode", "Rep1", "Rep2", "Rep3")
bcNum_prop_UT_test$Sample <- "Untreated"

bcNum_prop_UT_testMeans <- bcNum_prop_UT_test[,c(2:4)]
bcNum_prop_UT_testMeans$mean <- rowMeans(bcNum_prop_UT_testMeans)
bcNum_prop_UT_testMeans$Barcode <- unlist(bcNum_prop_UT_test[,c(1)])
bcNum_prop_UT_testMeans_order <- bcNum_prop_UT_testMeans[order(bcNum_prop_UT_testMeans$mean, decreasing=TRUE),]

bcNum_prop_I_test <- bcNum_prop[,c(1,5:7)]
names(bcNum_prop_I_test) <- c("Barcode", "Rep1", "Rep2", "Rep3")
bcNum_prop_I_test$Sample <- "Idling"

bcNum_prop_I_testMeans <- bcNum_prop_I_test[,c(2:4)]
bcNum_prop_I_testMeans$mean <- rowMeans(bcNum_prop_I_testMeans)
bcNum_prop_I_testMeans$Barcode <- unlist(bcNum_prop_I_test[,c(1)])
bcNum_prop_I_testMeans_order <- bcNum_prop_I_testMeans[order(bcNum_prop_I_testMeans$mean, decreasing=TRUE),]

bcNum_prop_test <- rbind(bcNum_prop_UT_test, bcNum_prop_I_test)
bcNum_prop_test_melt <- melt(bcNum_prop_test, id = c("Barcode", "Sample"))
bcNum_prop_test_melt$Barcode <- factor(bcNum_prop_test_melt$Barcode,
                                  levels=as.character(bcNum_prop_UT_testMeans_order$Barcode))

bcNum_prop_compare <- Rmisc::summarySE(bcNum_prop_test_melt, measurevar = "value",
                           groupvars = c("Barcode", "Sample"))

bcNum_prop_compare$Sample <- factor(bcNum_prop_compare$Sample,
                                    levels = c('Untreated', 'Idling'))

bcNum_prop_compare_sub <- bcNum_prop_compare[c(1:48),]
ggplot(bcNum_prop_compare_sub, aes(x=Barcode, y=value, group = Sample,
                              fill = Sample)) +
  geom_bar(stat = "identity", position = "dodge", color = "black", linewidth = 0.2) +
  geom_errorbar(aes(ymin=value-sd, ymax=value+sd), width=.2,
                position=position_dodge(.9)) +
  # geom_ribbon(data = bcNum_model_prop_compare_sub,
  #             aes(x = Barcode, y = value, ymin=value-sd, ymax=value+sd),
  #             alpha=0.1, color = "blue", fill = "blue") +
  # coord_flip() + ylab("Barcode Fraction") + xlab("Population") +
  xlab("Barcode") + ylab("Barcode Fraction") +
  theme_classic() + scale_fill_manual(values = c('red', 'blue')) +
  scale_x_discrete(labels= factor(seq(nrow(bcNum_prop_compare_sub)/2))) +
  theme(
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    axis.text.x = element_text(size = 14, colour = "black", angle = 90, vjust = 0.5, hjust = 1),
    axis.text.y = element_text(size = 14, colour = "black"),
    legend.position = "none", legend.text = element_text(size = 14),
    plot.title = element_text(size = 16, hjust = 0.5, face = "bold"), axis.text=element_text(size=14),
    legend.title = element_text(size=14), axis.title=element_text(size=14))

ggsave("SKMEL5_barcode_propRankAbundance_comparison.pdf") #, width = 6, height = 4)

### Fold Change Plot
bcFC <- merge(x = bcNum_prop_UT_testMeans, y = bcNum_prop_I_testMeans,
              by = "Barcode", all.x = T)
bcFC_order <- bcFC[order(bcFC$mean.x, decreasing=TRUE),][,c("Barcode", "mean.x", "mean.y")]
bcFC_order$l2FC <- log2(bcFC_order$mean.y/bcFC_order$mean.x)

# ggplot(bcFC_order, aes(l2FC)) + theme_bw() +
#   geom_density(color = "black", fill = "grey") +
#   # geom_density(data = bcFC_order[c(1:25),],
#   #              aes(l2FC), color = "black", fill = "blue")
#   xlab("Log2 Fold Change") + ylab("Density") +
#   xlim(-3.5,3.5) +
#   theme(
#     panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#     axis.text = element_text(size = 12),
#     legend.position = "none", legend.text = element_text(size = 12),
#     plot.title = element_text(size = 14, hjust = 0.5, face = "bold"),
#     legend.title = element_text(size=12), axis.title=element_text(size=12)) 
#   ggsave("SKMEL5_barcode_FCdensity.pdf", width = 4, height = 3)

# test_hist <- data.frame(x = c(bcFC_order$l2FC, bcFC_order$l2FC[1:25]),
#                         Class = as.factor(c(rep("Whole", length(bcFC_order$l2FC)),
#                                             rep("Sub", length(bcFC_order$l2FC[1:25])))),
#                         BCnum = as.factor(c(seq(length(bcFC_order$l2FC)),
#                                             seq(length(bcFC_order$l2FC[1:25])))))

test_hist1 <- data.frame(x = bcFC_order$l2FC[1:24],
                        Class = as.factor(rep("Sub", length(bcFC_order$l2FC[1:24]))),
                        BCnum = as.factor(seq(length(bcFC_order$l2FC[1:24]))),
                        Barcode = bcFC_order$Barcode[1:24])

n <- dim(test_hist1)[1]
# bw <- 3.49 * sd(test_hist1[, "x"]) * dim(test_hist1)[1]^(-1/3)
bw = 1

plt_hist_BCoverlay <- ggplot(test_hist1, aes(x=x)) + theme_bw() +
  geom_histogram(aes(y= (after_stat(count))/(n*bw), fill = BCnum),
                 binwidth = 1, color = "black", linewidth = 0.1) +
  scale_fill_manual(values = rainbow(24)) +
  labs(fill = "Barcode") +
  geom_density(data = bcFC_order, aes(l2FC), fill = "grey",
               color = "black", alpha = 0.3) +
  xlab("Log2 Fold Change") + ylab("Density") +
  xlim(-3.5,3.5) +
  theme(
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    axis.text = element_text(size = 12),
    legend.position = "right", legend.text = element_text(size = 12),
    plot.title = element_text(size = 14, hjust = 0.5, face = "bold"),
    legend.title = element_text(size=12), axis.title=element_text(size=12))

plt_hist_BCoverlay 
ggsave("SKMEL5_barcode_FCdensity_bcOverlay.pdf") #, width = 3.5, height = 3)
plt_hist_BCoverlay_leg <- ggpubr::get_legend(plt_hist_BCoverlay)
ggpubr::as_ggplot(plt_hist_BCoverlay_leg)
ggsave("SKMEL5_barcode_FCdensity_bcOverlay_leg.pdf")
