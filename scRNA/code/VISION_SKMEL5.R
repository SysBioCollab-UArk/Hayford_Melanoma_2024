library(VISION)
library(Seurat)
library(limma)
library(stringr)
library(ggplot2)
library(ggpubr)

# if (Sys.getenv("RSTUDIO") == "1") {
#   setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
# }

if (!file.exists("combined_includingState.RData")){
  source("Seurat_v5_SKMEL5_combined_hg38.R")
}
load("combined_includingState.RData")

## Run VISION Analysis
counts <- as.data.frame(combined@assays$SCT@counts)
meta <- as.data.frame(combined@meta.data)
meta$lineage <- factor(meta$lineage)
meta$lineageColored <- factor(meta$lineageColored)

# Scale counts within a sample
n.umi <- colSums(counts)
scaled_counts <- t(t(counts) / n.umi) * median(n.umi)

# HALLMARK GENES
gmt_files <- list.files(path = file.path('..', 'data', 'VISION_hallmark'), 
                        pattern = "\\.gmt$")
# gmt_files1 <- paste0("../data/VISION_hallmark/", gmt_files, sep = "")
gmt_files1 <- file.path("..", "data", "VISION_hallmark", 
                        paste0(gmt_files, sep=""))
vis <- Vision(data = counts, signatures = gmt_files1, meta = meta)
vis <- analyze(vis)
visScores <- as.data.frame(getSignatureScores(vis))

## Plot overlay UMAP
combined_VISION <- combined
combined_VISION <- AddMetaData(object = combined, 
                                          metadata = visScores)

hallmark_metrics <- paste0("HALLMARK_", removeExt(gmt_files))

combined_plotDat <- 
  FetchData(
    combined_VISION, 
    vars = c("umap_1", "umap_2", "lineage", "Phase", "State", "old.ident", 
             hallmark_metrics, "HALLMARK_KRAS_SIGNALING", "HALLMARK_UV_RESPONSE"))

umap_vals <- as.data.frame(Embeddings(combined, reduction = "umap"))
umap_vals_df <- dplyr::bind_rows(umap_vals)
colnames(umap_vals_df) <- c("UMAP_1", "UMAP_2")
pops_all <- subset(combined@meta.data, select = "orig.ident")
umap_vals_df$Population <- unlist(pops_all)
all_Kmeans <- kmeans(umap_vals_df[,1:2,], centers = 3)
combined_plotDat$Kmeans <- all_Kmeans$cluster

plt_all_overlay <- ggplot() + theme_bw() +
  geom_density_2d(data = combined_plotDat, 
                  aes(x = UMAP_1, y = UMAP_2,
                      color = old.ident), alpha = 0.8) +
  scale_color_manual(values = c("blue", "red")) +
  geom_point(data = combined_plotDat, shape = 21,
             aes(x = UMAP_1, y = UMAP_2, fill = HALLMARK_GLYCOLYSIS),
             size = 0.8, stroke = 0.01) + 
  scale_fill_gradient(guide = FALSE, 
                      low = "white", 
                      high = "black") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text = element_text(size = 12), legend.position = "none", 
        legend.text = element_text(size = 12), legend.title = element_text(size=12), 
        axis.title=element_text(size=12)) 

## Plot distribution of all 50 hallmark sets for each population
combined_dat_hallmarks <- combined_plotDat[7:ncol(combined_plotDat)]
combined_dat_hallmarks_melt <- reshape2::melt(combined_dat_hallmarks, id.vars = "Kmeans")

combined_dat_hallmarks_melt$variable <- 
  str_remove(as.character(combined_dat_hallmarks_melt$variable), "HALLMARK_")
combined_dat_hallmarks_melt$Kmeans <- factor(combined_dat_hallmarks_melt$Kmeans,
                                            levels = c(1,2,3))

plt_hallmarks <- ggplot(
  combined_dat_hallmarks_melt, aes(x=value, color = Kmeans, fill = Kmeans)) +
  geom_density(alpha = 0.5) + facet_wrap(~variable, ncol = 6, scales = "free") + 
  theme_bw() +
  scale_color_manual(values = c("magenta", "gold", "skyblue"), 
                     labels = c(1,2,3), name = "Clusters") +
  scale_fill_manual(values = c("magenta", "gold", "skyblue"), 
                     labels = c(1,2,3), name = "Clusters") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        legend.position = "right", axis.text = element_text(size = 18), 
        axis.title = element_text(size = 18)) +
  xlab("Signature Score") + ylab("Density") 

plt_hallmarks 
ggsave("allHallmarks_acrossSKMEL5clusters.svg", width = 18, height = 24)
plt_hallmarks_leg <- ggpubr::get_legend(plt_hallmarks)
as_ggplot(plt_hallmarks_leg) 
ggsave("allHallmarks_acrossSKMEL5clusters_legend.svg", width = 6, height = 4)

# combined_dat_hallmarks_melt_subset <- subset(combined_dat_hallmarks_melt, variable == c("GLYCOLYSIS"))
# for (var in c("OXIDATIVE_PHOSPHORYLATION", "DNA_REPAIR", "FATTY_ACID_METABOLISM", 
#               "P53_PATHWAY", "REACTIVE_OXYGEN_SPECIES_PATHWAY")){
#   combined_dat_hallmarks_melt_subset <- rbind(combined_dat_hallmarks_melt_subset, 
#                                               subset(combined_dat_hallmarks_melt, 
#                                                      variable == var))
# }

combined_dat_hallmarks_melt_subset <- 
  subset(combined_dat_hallmarks_melt, 
         variable == c("GLYCOLYSIS", "OXIDATIVE_PHOSPHORYLATION", "DNA_REPAIR", 
                       "FATTY_ACID_METABOLISM", "P53_PATHWAY", 
                       "REACTIVE_OXYGEN_SPECIES_PATHWAY"))

# rename variables
combined_dat_hallmarks_melt_subset$variable <- 
  gsub("OXIDATIVE_PHOSPHORYLATION", "OXIDATIVE PHOSPHORYLATION", 
       combined_dat_hallmarks_melt_subset$variable)
combined_dat_hallmarks_melt_subset$variable <- 
  gsub("DNA_REPAIR", "DNA REPAIR", 
       combined_dat_hallmarks_melt_subset$variable)
combined_dat_hallmarks_melt_subset$variable <- 
  gsub("FATTY_ACID_METABOLISM", "FATTY ACID METABOLISM", 
       combined_dat_hallmarks_melt_subset$variable)
combined_dat_hallmarks_melt_subset$variable <- 
  gsub("P53_PATHWAY", "P53 PATHWAY", 
       combined_dat_hallmarks_melt_subset$variable)
combined_dat_hallmarks_melt_subset$variable <- 
  gsub("REACTIVE_OXYGEN_SPECIES_PATHWAY", "ROS PATHWAY", 
       combined_dat_hallmarks_melt_subset$variable)

# define order of variables for plotting
combined_dat_hallmarks_melt_subset$variable <- 
  factor(combined_dat_hallmarks_melt_subset$variable, 
         levels = c("GLYCOLYSIS", "OXIDATIVE PHOSPHORYLATION", "DNA REPAIR", 
                    "FATTY ACID METABOLISM", "P53 PATHWAY", "ROS PATHWAY"))

########### FIGURE 1C ###########
plt_hallmarks_sub <- ggplot(
  combined_dat_hallmarks_melt_subset, aes(x=value, color = Kmeans, fill = Kmeans)) +
  geom_density(alpha = 0.5) + 
  facet_wrap(~variable, ncol = 6, scales = "free", labeller = label_wrap_gen(width=20)) + 
  theme_bw() +
  scale_color_manual(values = c("magenta", "gold", "skyblue"), 
                     labels = c("UTS","UTL","Idling"), name = NULL) +
  scale_fill_manual(values = c("magenta", "gold", "skyblue"), 
                    labels = c("UTS","UTL","Idling"), name = NULL) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        legend.position = "top", legend.text = element_text(size = 20),
        axis.text = element_text(size = 20), axis.title = element_text(size = 24),
        axis.text.y = element_blank(), panel.spacing = unit(2, "lines"),
        strip.text = element_text(size = 16)) +
  xlab("Signature Score") + ylab("Density") 

# plt_hallmarks_sub
ggsave("subsetHallmarks_acrossSKMEL5clusters.svg", width = 18, height = 4.5)
#################################
# plt_hallmarks_leg_sub <- ggpubr::get_legend(plt_hallmarks_sub)
# as_ggplot(plt_hallmarks_leg_sub) 
# ggsave("subsetHallmarks_acrossSKMEL5clusters_legend.svg",
#        width = 6, height = 4)
