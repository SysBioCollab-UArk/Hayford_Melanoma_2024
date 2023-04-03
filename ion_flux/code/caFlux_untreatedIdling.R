library(ggplot2)
library(reshape2)

caFlux_try2_dat <- read.csv('./ion_flux/data/SingleAdd_IP3_Fluo-8_May_01_19.txt', 
                       sep = "\t")
names(caFlux_try2_dat) <- names(caFlux_try2_dat[,-1])
caFlux_try2_dat <- caFlux_try2_dat[,1:384]
caFlux_try2_dat$Time <- rownames(caFlux_try2_dat)
caFlux_try2_melt <- melt(caFlux_try2_dat, 
                    value.name = "Intensity", 
                    variable.name = "Well")
caFlux_try2_melt$time <- as.numeric(caFlux_try2_melt$Time)/1000
caFlux_try2_melt$row <- substr(caFlux_try2_melt$Well, 1, 1)
caFlux_try2_melt$column <- substr(caFlux_try2_melt$Well, 2, 3)
caFlux_try2_melt<-caFlux_try2_melt[!(caFlux_try2_melt$column %in% c("1","2","23","24")),]
caFlux_try2_melt$Agonist <- rep(c("Ionomycin", "Thapsigargin"), each = 322*5)
caFlux_try2_melt$concentration <- rep(c("High", "Medium-High", "Medium", "Medium-Low", "Low"), each = 322)
caFlux_try2_melt$pop <- rep(c("Untreated", "Idling"), each = 322*10)
caFlux_try2_test <- caFlux_try2_melt
dn <- data.frame()
for(i in unique(caFlux_try2_test$Well)){
  # print(i)
  temp <- subset(caFlux_try2_test, caFlux_try2_test$Well==i)
  temp$normIntensity <- temp$Intensity - temp$Intensity[temp$Time == temp$Time[which.min(temp$Intensity)]]
  temp1 <- temp[which.min(temp$Intensity):length(temp$Intensity),]
  dn <- rbind(dn, temp1)
}

caFlux_try2_iono <- subset(dn, Agonist == "Ionomycin")
caFlux_try2_iono_high <- subset(dn, Agonist == "Ionomycin" & concentration == "High")
ggplot_iono <- ggplot(caFlux_try2_iono, aes(x=time,y=normIntensity,color=concentration,fill=concentration))
ggplot_iono_high <- ggplot(caFlux_try2_iono_high, aes(x=time,y=normIntensity,color=pop,fill=pop))

# color=interaction(concentration, column), 
# fill=interaction(concentration, column))) # ,
# group = interaction(concentration, Well)))
# ggplot_iono + stat_summary(fun.data = "mean_cl_normal", geom = "smooth", se=T, aes(linetype=pop)) +
#   theme_bw() + theme(legend.text = element_text(size = 8), 
#                      plot.title = element_text(size = 16, hjust = 0.5, face = "bold"), axis.text=element_text(size=8),
#                      legend.title = element_text(size=14,face="bold"), axis.title=element_text(size=14,face="bold"),
#                      legend.position = "bottom") + 
#   labs(x = "Time (seconds)", y = "Normalized Intensity") + 
#   ggtitle("Calcium flux with Ionomycin agonist") + 
#   ggsave("CalciumFlux/CaFlux_try2_ionomycin.pdf", width = 6, height = 4)

# ggplot_iono_high + stat_summary(fun.data = "mean_cl_normal", geom = "smooth", se=T) +
#   theme_bw() + 
#   theme(legend.text = element_text(size = 12), legend.position = "bottom", 
#         plot.title = element_text(size = 12, hjust = 0.5), axis.text=element_text(size=12),
#         legend.title = element_blank(), axis.title=element_text(size=12),
#         panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
#   labs(x = "Time (seconds)", y = "Normalized Intensity") + scale_color_manual(values = c("blue", "red")) +
#   scale_fill_manual(values = c("blue", "red")) +
#   ggtitle("Calcium flux - Ionomycin Agonist (8uM)") + 
#   ggsave("CalciumFlux/CaFlux_try2_ionomycin_high.pdf", width = 6, height = 4)

caFlux_try2_thap <- subset(dn, Agonist == "Thapsigargin")
caFlux_try2_thap_high <- subset(dn, Agonist == "Thapsigargin" & concentration == "High")
ggplot_thap <- ggplot(caFlux_try2_thap, aes(x=time,y=normIntensity,color=concentration, fill=concentration)) # ,
ggplot_thap_high <- ggplot(caFlux_try2_thap_high, aes(x=time,y=normIntensity,color=pop,fill=pop))
# group = interaction(concentration, Well)))
ggplot_thap + stat_summary(fun.data = "mean_cl_normal", geom = "smooth", se=T, aes(linetype=pop)) 
# ggplot_thap_high + stat_summary(fun.data = "mean_cl_normal", geom = "smooth", se=T) +
#   theme_bw() + 
#   theme(legend.text = element_text(size = 12), legend.position = "bottom", 
#         plot.title = element_text(size = 12, hjust = 0.5), axis.text=element_text(size=12),
#         legend.title = element_blank(), axis.title=element_text(size=12),
#         panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
#   labs(x = "Time (seconds)", y = "Normalized Intensity") + scale_color_manual(values = c("blue", "red")) +
#   scale_fill_manual(values = c("blue", "red")) +
#   ggtitle("Calcium flux - Thapsigargin Agonist (20uM)") + 
#   ggsave("CalciumFlux/CaFlux_try2_thapsigargin_high.pdf", width = 6, height = 4)

### 
caFlux_data_high <- rbind(caFlux_try2_iono_high, caFlux_try2_thap_high)
ggplot(caFlux_data_high, aes(x=time,y=normIntensity,color=pop,fill=pop)) +
  facet_wrap(~Agonist) +
  stat_summary(fun.data = "mean_cl_normal", geom = "smooth", se=T) +
  theme_bw() + 
  theme(legend.text = element_text(size = 12), legend.position = "none", 
        plot.title = element_text(size = 12, hjust = 0.5), axis.text=element_text(size=12),
        legend.title = element_blank(), axis.title=element_text(size=12),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  labs(x = "Time (seconds)", y = "Normalized Intensity") + scale_color_manual(values = c("blue", "red")) +
  scale_fill_manual(values = c("blue", "red")) 
  ggsave("caFlux_high_Iono-Thap_UTvI_2021vis.pdf", width = 8, height = 3)
