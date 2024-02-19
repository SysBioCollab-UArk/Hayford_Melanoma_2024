library(dplyr)
library(ggplot2)
library(Seurat)

if (Sys.getenv("RSTUDIO") == "1") {
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
}

load("combined_includingState.RData")
DimPlot(combined, reduction = "umap", group.by = "seurat_clusters")

load("../data/top25BCs.RData")
bcs_top25 <- as.character(unique(bcNum_prop_compare_sub$Barcode))
bc_not <- setdiff(combined@meta.data$lineage, bcs_top25)
test <- combined@meta.data
test$lineageColored[test$lineage %in% bcs_top25] <- test$lineage[test$lineage %in% bcs_top25]
test$lineageColored[test$lineage %in% bc_not] <- "XXX"
lineageColored <- subset(test, select = c("lineageColored"))
combined <- AddMetaData(combined, lineageColored, col.name = "lineageColored")

t2 <- as.data.frame(combined@meta.data)
t2$BCnum <- as.integer(factor(t2$lineageColored, levels = unique(bcNum_prop_compare_sub$Barcode)))
BCnum <- subset(t2, select = c("BCnum"))
combined <- AddMetaData(combined, BCnum, col.name = "BCnum")

df <- FetchData(combined, vars = c("umap_1", "umap_2", "orig.ident", "lineage", "lineageColored", "seurat_clusters", "BCnum"))
# Convert Seurat clusters to meta-clusters
df1 <- df %>% 
  mutate(meta_cluster = case_when(
    seurat_clusters %in% c("8") ~ 'UT_S',
    seurat_clusters %in% c("6") ~ 'I_S',
    seurat_clusters %in% c("0","2","5","7") ~ 'I_L',
    TRUE ~ 'UT_L' ) )

# Add large vs small information
df2 <- df1 %>% 
  mutate(cluster_size = case_when(
    seurat_clusters %in% c("8","6") ~ 'Small',
    TRUE ~ 'Large' ) )

# # Percentages
meta_cluster_UT_table_byBC <- subset(df1, orig.ident == "Untreated")[,c("meta_cluster", "BCnum")]
meta_cluster_UT_byBC <- meta_cluster_UT_table_byBC %>%
  dplyr::group_by(meta_cluster, BCnum) %>%
  dplyr::summarise (n = n()) %>%
  dplyr::group_by(BCnum) %>%
  dplyr::mutate(freq = n / sum(n))
meta_cluster_UT_byBC$Condition <- "Untreated"

meta_cluster_I_table_byBC <- subset(df1, orig.ident == "Idling")[,c("meta_cluster", "BCnum")]
meta_cluster_I_byBC <- meta_cluster_I_table_byBC %>%
  dplyr::group_by(meta_cluster, BCnum) %>%
  dplyr::summarise (n = n()) %>%
  dplyr::group_by(BCnum) %>%
  dplyr::mutate(freq = n / sum(n))
meta_cluster_I_byBC$Condition <- "Idling"

meta_cluster_byBC <- as.data.frame(rbind(meta_cluster_UT_byBC, meta_cluster_I_byBC))

ggplot(subset(meta_cluster_byBC, Condition == "Untreated"), aes(x = as.factor(BCnum), y = freq,
                                                         group = meta_cluster, fill = meta_cluster)) +
  theme_bw() + geom_bar(position = position_fill(reverse = TRUE),
                        stat = "identity", color = "black") +
  scale_fill_manual(values = c("pink", "grey", "purple", "brown")) +
  scale_y_continuous(labels = scales::percent) +
  labs(x = "Barcode", y = "Percentage of Cells") +
  ggtitle("Untreated Meta-cluster") +
  theme(axis.text.y = element_text(size = 14),
        legend.position = "bottom", legend.text = element_text(size = 14),
        plot.title = element_text(size = 16, hjust = 0.5),
        axis.text=element_text(size=14),
        legend.title = element_blank(), axis.title=element_text(size=14),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

ggplot(subset(meta_cluster_byBC, Condition == "Idling"), aes(x = as.factor(BCnum), y = freq,
                                                                group = meta_cluster, fill = meta_cluster)) +
  theme_bw() + geom_bar(position = position_fill(reverse = TRUE),
                        stat = "identity", color = "black") +
  scale_fill_manual(values = c("pink", "grey", "purple", "brown")) +
  scale_y_continuous(labels = scales::percent) +
  labs(x = "Barcode", y = "Percentage of Cells") +
  ggtitle("Idling Meta-cluster") +
  theme(axis.text.y = element_text(size = 14),
        legend.position = "bottom", legend.text = element_text(size = 14),
        plot.title = element_text(size = 16, hjust = 0.5),
        axis.text=element_text(size=14),
        legend.title = element_blank(), axis.title=element_text(size=14),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


######

# # Percentages
meta_cluster_UT_table_byBC_x <- subset(df2, orig.ident == "Untreated")[,c("orig.ident", "BCnum", "cluster_size")]
meta_cluster_UT_byBC_x <- meta_cluster_UT_table_byBC_x %>%
  dplyr::group_by(orig.ident, cluster_size, BCnum) %>%
  dplyr::summarise (n = n()) %>%
  dplyr::group_by(BCnum) %>%
  dplyr::mutate(freq = n / sum(n))

meta_cluster_I_table_byBC_x <- subset(df2, orig.ident == "Idling")[,c("orig.ident", "BCnum", "cluster_size")]
meta_cluster_I_byBC_x <- meta_cluster_I_table_byBC_x %>%
  dplyr::group_by(orig.ident, cluster_size, BCnum) %>%
  dplyr::summarise (n = n()) %>%
  dplyr::group_by(BCnum) %>%
  dplyr::mutate(freq = n / sum(n))

meta_cluster_byBC_x <- as.data.frame(rbind(meta_cluster_UT_byBC_x, meta_cluster_I_byBC_x))
meta_cluster_byBC_x <- subset(meta_cluster_byBC_x, BCnum != "8")

meta_cluster_byBC_x$orig.ident <- factor(meta_cluster_byBC_x$orig.ident, levels = c("Untreated", "Idling"))

ggplot(subset(meta_cluster_byBC_x, cluster_size == "Large"), aes(x = as.factor(BCnum), y = freq,
                                                                group = orig.ident, fill = orig.ident)) +
  theme_bw() + geom_bar(width=.5, position = "dodge",
                        stat = "identity", color = "black") +
  scale_fill_manual(values = c("purple", "pink")) +
  scale_y_continuous(labels = scales::percent) +
  coord_cartesian(ylim = c(0.65,0.95)) +
  labs(x = "Barcode", y = "Percentage of Cells") +
  ggtitle("Large Cluster in BRAFi") +
  theme(axis.text.y = element_text(size = 14),
        legend.position = "bottom", legend.text = element_text(size = 14),
        plot.title = element_text(size = 16, hjust = 0.5),
        axis.text=element_text(size=14),
        legend.title = element_blank(), axis.title=element_text(size=14),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

ggplot(subset(meta_cluster_byBC_x, cluster_size == "Small"), aes(x = as.factor(BCnum), y = freq,
                                                                 group = orig.ident, fill = orig.ident)) +
  theme_bw() + geom_bar(width=.5, position = "dodge",
                        stat = "identity", color = "black") +
  scale_fill_manual(values = c("brown", "grey")) +
  scale_y_continuous(labels = scales::percent) +
  # coord_cartesian(ylim = c(0.65,0.95)) +
  labs(x = "Barcode", y = "Percentage of Cells") +
  ggtitle("Small Cluster in BRAFi") +
  theme(axis.text.y = element_text(size = 14),
        legend.position = "bottom", legend.text = element_text(size = 14),
        plot.title = element_text(size = 16, hjust = 0.5),
        axis.text=element_text(size=14),
        legend.title = element_blank(), axis.title=element_text(size=14),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

##
# ggplot(meta_cluster_byBC_x, aes(x = orig.ident, y = freq, fill = cluster_size)) +
#   theme_bw() + geom_bar(width=.5,
#                         stat = "identity", color = "black") +
#   facet_wrap(vars(as.factor(BCnum)), ncol = 5) +
#   scale_fill_manual(values = c("brown", "grey")) +
#   scale_y_continuous(labels = scales::percent) +
#   # coord_cartesian(ylim = c(0.65,0.95)) +
#   labs(x = "", y = "Percentage of Cells") +
#   # ggtitle("Cluster in BRAFi") +
#   theme(axis.text.y = element_text(size = 12),
#         legend.position = "none", legend.text = element_text(size = 12),
#         plot.title = element_text(size = 16, hjust = 0.5),
#         axis.text=element_text(size=12),
#         legend.title = element_blank(), axis.title=element_text(size=12),
#         panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#         axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
#   ggsave(".", width = 10, height =10)

####### Correlation

load("../data/barcodeFC_forCorPlot.RData")
names(test_hist1) <- c("FC", "Class", "BCnum", "lineage")
bcFC <- test_hist1[,c("lineage", "BCnum", "FC")]

##
library(ggpubr)
library(ggrepel)
# BC-by-BC Number/Proportion of cells in each CC state
stateBC_df <- FetchData(combined, vars = c("orig.ident", "lineage", "lineageColored", "State"))
# stateBC_df <- subset(stateBC_df, lineageColored != "XXX")
t1 <- stateBC_df
t1$BCnum <- as.integer(factor(t1$lineageColored, levels = unique(bcNum_prop_compare_sub$Barcode)))
t1$lineageColored <- factor(t1$lineageColored, levels = unique(bcNum_prop_compare_sub$Barcode))
state_I_table_byBC <- subset(t1, orig.ident == "Idling")[,c("State", "BCnum", "lineage", "lineageColored")]
state_I_byBC <- state_I_table_byBC %>% 
  dplyr::group_by(State, BCnum) %>% 
  dplyr::summarise (n = n()) %>%
  dplyr::group_by(BCnum) %>%
  dplyr::mutate(freq = n / sum(n))
state_I_byBC$Condition <- "Idling"
state_bc_I <- as.data.frame(subset(state_I_byBC, State == "Fast_Dividing"))
state_bc_I <- subset(state_bc_I, BCnum != "NA")
state_bc_I_sub <- state_bc_I[,c("BCnum", "freq")]

compareBCs_top25 <- merge(bcFC, state_bc_I_sub, by = "BCnum")
names(compareBCs_top25) <- c("BCnum", "lineage", "FC", "prop_I")

ggscatter(compareBCs_top25, x="prop_I",y="FC", add="reg.line", size=1) +
  stat_cor(method="pearson", aes(label = ..r.label..), size = 5) +
  ggrepel::geom_label_repel(aes(label=BCnum)) +
  scale_x_continuous(labels = scales::percent) +
  labs(x="Percentage of Fast-Diving Idling DTPs", 
       y=expression(Log[2]~"Barcode Fold Change")) +
  theme(axis.text = element_text(size = 14)) 
ggsave("Idling_bcFC_correlation.pdf", width = 6, height = 4)
