library(tidyverse)
library(ggpubr)

# UT_lc_res <- read_tsv('/Volumes/Transcend/ATACseq/preseq_results/trimmed_3334-CH-1-TAGGCATG-ACTGCATA_S121_aligned_sorted_lcExtrapResults.txt') %>%
#   mutate(Library="Untreated")
# 
# I_lc_res_25 <- read_tsv('idling_subsample_25p_lcExtrapResults.txt') %>%
#   mutate(Library="Idling_25")
# I_lc_res_33 <- read_tsv('idling_subsample_33p_lcExtrapResults.txt') %>%
#   mutate(Library="Idling_33")
# I_lc_res_50 <- read_tsv('idling_subsample_33p_lcExtrapResults.txt') %>%
#   mutate(Library="Idling_50")
# 
# I_lc_res_25_dedup <- read_tsv('idling_subsample_25p_dedup_lcExtrapResults.txt') %>%
#   mutate(Library="Idling_25_dedup")
# I_lc_res_33_dedup <- read_tsv('idling_subsample_33p_dedup_lcExtrapResults.txt') %>%
#   mutate(Library="Idling_33_dedup")
# I_lc_res_50_dedup <- read_tsv('idling_subsample_50p_dedup_lcExtrapResults.txt') %>%
#   mutate(Library="Idling_50_dedup")
# 
# lc_res <- bind_rows(UT_lc_res, I_lc_res_25, I_lc_res_33, I_lc_res_50,
#                     I_lc_res_25_dedup, I_lc_res_33_dedup, I_lc_res_50_dedup)
# 
# lc_res <- as.data.frame(lc_res)
# 
# ggplot(lc_res) + theme_bw() +
#   geom_ribbon(aes(x=TOTAL_READS/10^6, 
#                   ymin=LOWER_0.95CI/10^6, 
#                   ymax=UPPER_0.95CI/10^6, 
#                   fill=Library), alpha = 0.4) +
#   geom_line(aes(x=TOTAL_READS/10^6, 
#                 y=EXPECTED_DISTINCT/10^6, 
#                 color=Library), 
#             size=1, linetype="dashed") 

#####


c_res_UT <- read_tsv('/Volumes/Transcend/ATACseq/preseq_results/trimmed_3334-CH-1-TAGGCATG-ACTGCATA_S121_aligned_sorted_cCurveResults.txt') %>%
  mutate(Library="Untreated")

c_res_I <- read_tsv('/Volumes/Transcend/ATACseq/preseq_results/trimmed_3334-CH-2-CGTACTAG-CTCTCTAT_S122_aligned_sorted_cCurveResults.txt') %>%
  mutate(Library="Idling")

c_res_I_25 <- read_tsv('idling_subsample_25p_dedup_cCurveResults.txt') %>%
  mutate(Library="Idling_25")
c_res_I_33 <- read_tsv('idling_subsample_33p_dedup_cCurveResults.txt') %>%
  mutate(Library="Idling_33")
c_res_I_50 <- read_tsv('idling_subsample_50p_dedup_cCurveResults.txt') %>%
  mutate(Library="Idling_50")

c_res_I_25_dedup <- read_tsv('idling_subsample_25p_dedup_cCurveResults.txt') %>%
  mutate(Library="Idling_25_dedup")
c_res_I_33_dedup <- read_tsv('idling_subsample_33p_dedup_cCurveResults.txt') %>%
  mutate(Library="Idling_33_dedup")
c_res_I_50_dedup <- read_tsv('idling_subsample_50p_dedup_cCurveResults.txt') %>%
  mutate(Library="Idling_50_dedup")

c_res <- bind_rows(c_res_UT, c_res_I_25, c_res_I_25_dedup)

c_res <- as.data.frame(c_res)

ggplot(c_res) + theme_bw() +
  geom_line(aes(x=total_reads/10^6, 
                y=distinct_reads/10^6, 
                color=Library, linetype = Library), 
            size=1) + 
  # coord_cartesian(x=c(0,75), y=c(0,37.5)) +
  xlab("Sequenced reads (M)") + 
  ylab("Distinct reads (M)") +
  ggtitle("Observed Complexity of ATAC-seq Libraries") +
  # scale_x_continuous(breaks = seq(0, 75, by = 5)) +
  # scale_y_continuous(breaks = seq(0, 37.5, by = 2.5)) +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        legend.position = "right",
        axis.text=element_text(size=14),
        axis.title=element_text(size=14),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())


# Final plots comparing libraries
## Not modified
c_res_full <- bind_rows(c_res_UT, c_res_I)
c_res_full <- as.data.frame(c_res_full)
c_res_full$Library <- factor(c_res_full$Library, levels = c("Untreated", "Idling"))

## Idling subsampled (25%)
c_res_sub <- bind_rows(c_res_UT, c_res_I_25)
c_res_sub <- as.data.frame(c_res_sub)
c_res_sub$Library <- factor(c_res_sub$Library, levels = c("Untreated", "Idling_25"))

# ## Not modified plot
# ggplot(c_res_full) + theme_bw() +
#   geom_line(aes(x=total_reads/10^6, 
#                 y=distinct_reads/10^6, 
#                 color=Library), 
#             size=1) + 
#   xlab("Sequenced reads (M)") + 
#   ylab("Distinct reads (M)") +
#   # ggtitle("Observed Complexity of ATAC-seq Libraries") +
#   scale_color_manual(values = c("red", "blue")) +
#   theme(plot.title = element_text(hjust = 0.5),
#         plot.subtitle = element_text(hjust = 0.5),
#         legend.position = "none",
#         axis.text=element_text(size=12),
#         axis.title=element_text(size=12),
#         panel.grid.major = element_blank(), 
#         panel.grid.minor = element_blank()) +
#   ggsave("cCurve_Observed_SKMEL5.pdf", width = 4, height = 3)

# Idling subsampled plot
ggplot(c_res_sub) + theme_bw() +
  geom_line(aes(x=total_reads/10^6, 
                y=distinct_reads/10^6, 
                color=Library), 
            size=1) + 
  xlab("Sequenced reads (M)") + 
  ylab("Distinct reads (M)") +
  # ggtitle("Observed Complexity of ATAC-seq Libraries") +
  scale_color_manual(values = c("red", "cyan")) +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        legend.position = "none",
        axis.text=element_text(size=12),
        axis.title=element_text(size=12),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  ggsave("cCurve_IdlingSubsampled_SKMEL5.pdf", width = 4, height = 3)


#### Insert Size Metrics (ISM)
# # Not deduplicate
# ## Untreated
# insert_UT <- read_tsv('/Volumes/Transcend/ATACseq/preseq_results/trimmed_3334-CH-1-TAGGCATG-ACTGCATA_S121_aligned_sorted_ISM.txt', 
#                       skip = 10) %>%
#   mutate(Library="Untreated")
# insert_UT <- mutate(insert_UT, cpm=(All_Reads.fr_count/sum(insert_UT$All_Reads.fr_count))*1e6)
# 
# ## Idling
# insert_I <- read_tsv('/Volumes/Transcend/ATACseq/preseq_results/trimmed_3334-CH-2-CGTACTAG-CTCTCTAT_S122_aligned_sorted_ISM.txt', 
#                      skip = 10) %>%
#   mutate(Library="Idling")
# insert_I <- mutate(insert_I, cpm=(All_Reads.fr_count/sum(insert_I$All_Reads.fr_count))*1e6)
# 
# insert_all <- bind_rows(insert_UT, insert_I)
# insert_all$Library <- factor(insert_all$Library, 
#                              levels = c("Untreated", "Idling"))
# 
# ggplot(insert_all, aes(x=insert_size, y=cpm, color=Library)) +
#   geom_line(size = 0.7) +
#   facet_grid(.~Library) +
#   scale_color_manual(values = c("red", "blue", "purple", "green", "orange")) +
#   coord_cartesian(x=c(0,800)) +
#   xlab("Insert Size (bp)") + 
#   ylab("Counts per Million") +
#   ggtitle("Insert Size Distribution of ATAC-seq Libraries", subtitle = "0-800bp") +
#   scale_x_continuous(breaks = seq(0, 900, by = 150)) +
#   theme_bw() +
#   theme(plot.title = element_text(hjust = 0.5),
#         plot.subtitle = element_text(hjust = 0.5),
#         legend.position = "none",
#         axis.text=element_text(size=14),
#         axis.title=element_text(size=14),
#         panel.grid.major = element_blank(), 
#         panel.grid.minor = element_blank()) 

# Deduplicated
## Normal
## Untreated
insert_UT_dedup <- read_tsv('/Volumes/Transcend/ATACseq/subsample/untreated_dedup_ISM.txt', 
                      skip = 10) %>%
  mutate(Library="Untreated")
insert_UT_dedup <- mutate(insert_UT_dedup, cpm=(All_Reads.fr_count/sum(insert_UT_dedup$All_Reads.fr_count))*1e6)

## Idling
insert_I_dedup <- read_tsv('/Volumes/Transcend/ATACseq/subsample/idling_dedup_ISM.txt', 
                     skip = 10) %>%
  mutate(Library="Idling")
insert_I_dedup <- mutate(insert_I_dedup, cpm=(All_Reads.fr_count/sum(insert_I_dedup$All_Reads.fr_count))*1e6)

## Subsampled Idling
## 25 percent
insert_I_25_dedup <- read_tsv('idling_subsample_25p_dedup_ISM.txt', skip = 10) %>%
  mutate(Library="Idling_25")
insert_I_25_dedup <- mutate(insert_I_25_dedup, cpm=(All_Reads.fr_count/sum(insert_I_25_dedup$All_Reads.fr_count))*1e6)

## 33 percent
insert_I_33_dedup <- read_tsv('idling_subsample_33p_dedup_ISM.txt', skip = 10) %>%
  mutate(Library="Idling_33")
insert_I_33_dedup <- mutate(insert_I_33_dedup, cpm=(All_Reads.fr_count/sum(insert_I_33_dedup$All_Reads.fr_count))*1e6)

## 50 percent
insert_I_50_dedup <- read_tsv('idling_subsample_50p_dedup_ISM.txt', skip = 10) %>%
  mutate(Library="Idling_50")
insert_I_50_dedup <- mutate(insert_I_50_dedup, cpm=(All_Reads.fr_count/sum(insert_I_50_dedup$All_Reads.fr_count))*1e6)


insert_all_dedup <- bind_rows(insert_UT_dedup, insert_I_dedup, insert_I_25_dedup, 
                              insert_I_33_dedup, insert_I_50_dedup)
insert_all_dedup$Library <- factor(insert_all_dedup$Library, 
                             levels = c("Untreated", "Idling", "Idling_25", "Idling_33", "Idling_50"))

# ggplot(insert_all_dedup, aes(x=insert_size, y=All_Reads.fr_count, color=Library)) +
#   geom_line(size = 0.7) +
#   facet_grid(.~Library) +
#   scale_color_manual(values = c("red", "blue", "purple", "green", "orange")) +
#   coord_cartesian(x=c(0,800)) +
#   xlab("Insert Size (bp)") + 
#   ylab("Read Count") +
#   ggtitle("Insert Size Distribution of ATAC-seq Libraries", subtitle = "0-800bp") +
#   scale_x_continuous(breaks = seq(0, 900, by = 150)) +
#   theme_bw() +
#   theme(plot.title = element_text(hjust = 0.5),
#         plot.subtitle = element_text(hjust = 0.5),
#         legend.position = "none",
#         axis.text=element_text(size=14),
#         axis.title=element_text(size=14),
#         panel.grid.major = element_blank(), 
#         panel.grid.minor = element_blank()) 


### Final Plots from subsampled libraries
insert_subsampled_dedup <- bind_rows(insert_UT_dedup, insert_I_25_dedup)
insert_subsampled_dedup$Library <- factor(insert_subsampled_dedup$Library, 
                                        levels = c("Untreated", "Idling_25"))

# ggplot(insert_observed_dedup, aes(x=insert_size, y=All_Reads.fr_count, color=Library)) +
#   geom_line(size = 0.7) +
#   scale_color_manual(values = c("red", "blue")) +
#   coord_cartesian(x=c(0,800)) +
#   xlab("Insert Size (bp)") +
#   ylab("Read Count") +
#   # ggtitle("Insert Size Distribution of ATAC-seq Libraries", subtitle = "0-800bp") +
#   scale_x_continuous(breaks = seq(0, 900, by = 150)) +
#   theme_bw() +
#   theme(plot.title = element_text(hjust = 0.5),
#         plot.subtitle = element_text(hjust = 0.5),
#         legend.position = "none",
#         axis.text=element_text(size=14),
#         axis.title=element_text(size=14),
#         panel.grid.major = element_blank(), 
#         panel.grid.minor = element_blank()) +
#   ggsave("ISM_Observed_SKMEL5.pdf", width = 4, height = 3)

# ggplot(insert_subsampled_dedup, aes(x=insert_size, y=All_Reads.fr_count, color=Library)) +
#   geom_line(size = 0.7) +
#   scale_color_manual(values = c("red", "cyan")) +
#   coord_cartesian(x=c(0,800)) +
#   xlab("Insert Size (bp)") +
#   ylab("Read Count") +
#   # ggtitle("Insert Size Distribution of ATAC-seq Libraries", subtitle = "0-800bp") +
#   scale_x_continuous(breaks = seq(0, 900, by = 150)) +
#   theme_bw() +
#   theme(plot.title = element_text(hjust = 0.5),
#         plot.subtitle = element_text(hjust = 0.5),
#         legend.position = "none",
#         axis.text=element_text(size=14),
#         axis.title=element_text(size=14),
#         panel.grid.major = element_blank(), 
#         panel.grid.minor = element_blank()) +
#   ggsave("ISM_IdlingSubsampled_SKMEL5.pdf", width = 4, height = 3)

# Create legend for colors
c_res_colors <- bind_rows(c_res_UT, c_res_I, c_res_I_25)
c_res_colors <- as.data.frame(c_res_colors)
c_res_colors$Library <- factor(c_res_colors$Library, levels = c("Untreated", "Idling", "Idling_25"))

ggplot(c_res_colors) + theme_bw() +
  geom_line(aes(x=total_reads/10^6, 
                y=distinct_reads/10^6, 
                color=Library), 
            size=1) + 
  xlab("Sequenced reads (M)") + 
  ylab("Distinct reads (M)") +
  # ggtitle("Observed Complexity of ATAC-seq Libraries") +
  scale_color_manual(values = c("red", "blue", "cyan")) +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        legend.position = "none",
        axis.text=element_text(size=12),
        axis.title=element_text(size=12),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  ggsave("cCurve_colors_SKMEL5.pdf", width = 4, height = 3)

insert_colors_dedup <- bind_rows(insert_UT_dedup, insert_I_dedup, insert_I_25_dedup)
insert_colors_dedup$Library <- factor(insert_colors_dedup$Library, 
                                        levels = c("Untreated", "Idling", "Idling_25"))

plt_color <- ggplot(insert_colors_dedup, aes(x=insert_size, y=All_Reads.fr_count, color=Library)) +
  geom_line(size = 0.7) +
  scale_color_manual(values = c("red", "blue", "cyan")) +
  coord_cartesian(x=c(0,800)) +
  xlab("Insert Size (bp)") +
  ylab("Read Count") +
  scale_x_continuous(breaks = seq(0, 900, by = 150)) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        legend.position = "right",
        axis.text=element_text(size=14),
        axis.title=element_text(size=14),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) 

plt_color +  ggsave("ISM_colors_SKMEL5.pdf", width = 4, height = 3)
leg <- get_legend(plt_color)
as_ggplot(leg) + ggsave("ATAClib_colors_legend.pdf")
