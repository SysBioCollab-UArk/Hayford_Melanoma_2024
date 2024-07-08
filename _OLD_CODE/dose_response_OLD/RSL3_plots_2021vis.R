library(lubridate)
library(gplots)
library(diprate)
library(ggplot2)
require(ggplot2)
require(Hmisc)
# source("~/Documents/QuarantaLab/SummarySE.R") # Summarize function to get the summary statistics;
# ============================================================================================================
# ============================================================================================================

# set working directory
if (Sys.getenv("RSTUDIO") == "1") {
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
}

listDir = list.files("../data", pattern = "*.txt")

listDir = paste("../data/", listDir, sep="")

countIJ <- function(platefile){ 
  
  plate1 = read.delim(platefile, header = T, sep="\t")
  plate1$Datetime = substr(plate1$Slice, 1, 14)
  plate1$Datetime = ymd_hms(plate1$Datetime)
  times = difftime(plate1$Datetime, plate1$Datetime[1], units = "hours")
  plate1$Time = as.numeric(times)
  plate1$Well = substr(plate1$Slice, (nchar(as.character(plate1$Slice)) - 6), nchar(as.character(plate1$Slice)))
  
  cellcount = aggregate(Count~ Time + Well, data = plate1, FUN = 'sum')
  cellcount$Row    = substr(cellcount$Well, 1, 3)
  cellcount$Column = substr(cellcount$Well, 5, 7)
  
  return(cellcount)
  
}


data = list()
test = list()
for(i in 1:length(listDir)){
  data[i] = list(countIJ(listDir[i]))
}

#Load platemaps

Eras = read.delim("../Platemaps/ErastinDR plate map.txt", header=TRUE, sep="\t")
RSL3 = read.delim("../Platemaps/RSL3DR plate map.txt", header=TRUE, sep="\t")

platDrug <- function(platemap, countdata){
  
  rows = as.character(platemap$X)
  cols = colnames(platemap)[2:12]
  rownames(platemap) = rows
  platemap$X = NULL
  
  for(i in 2:length(rows)){
    for(j in 2:length(cols)){
      Well = paste0(rows[i],"-",cols[j])
      countdata[countdata$Well == Well,"Drug"] = platemap[i,"Drug"]
      countdata[countdata$Well == Well,"Conc"] = platemap[i,j]
      countdata[countdata$Well == Well,"Units"] = platemap[i,"Conc"]
    }
  }
  
  return(countdata)
}

# ============================================================================================================
# ============================================================================================================

#Apply platemap and find the ln2 growth curve

iEras = platDrug(Eras, data.frame(data[1]))
iRSL3 = platDrug(RSL3, as.data.frame(data[2]))
pEras = platDrug(Eras, as.data.frame(data[3]))
pRSL3 = platDrug(RSL3, as.data.frame(data[4]))
comb = list(iEras,iRSL3,pEras,pRSL3)

getWellRates <- function(raw, time.range=c(70,120))
{
  #' Determine rates of proliferation of cells in each well over specified time range
  #' @param raw Raw \emph{data.frame}; expecting colnames of \code{time, well, date, cell.line, nl2}
  #'  Requires normalized log2(cell.count) \code{nl2} values as a column in \code{raw}
  #'  Uses linear model fit of  nl2 ~ time within time.range
  timeName <- colnames(raw)[grep('[Tt]ime', colnames(raw))]
  wellName <- colnames(raw)[grep('[Ww]ell',colnames(raw))]
  dateName <- colnames(raw)[grep('[Dd]ate',colnames(raw))]
  cellLineName <- colnames(raw)[grepl('[cC]ell',colnames(raw)) & grepl('[lL]ine',colnames(raw))]
  if(length(wellName)>1)    wellName <- wellName[nchar(wellName)==4]
  f <- formula(paste('nl2 ~ ',timeName,' * ',wellName))
  m <- lm(f, data=raw[raw[,timeName] > time.range[1] & raw[,timeName] < time.range[2],])
  wells <- unique(raw[,wellName])
  rates <- coef(m)[grep(timeName,names(coef(m)))]
  rates <- c(rates[1],rates[-1]+rates[1])
  cl     <- unique(raw[,cellLineName])
  expt <- ifelse(is.null(unique(raw[,dateName])), 'unknown date',unique(raw[,dateName]))
  out     <- data.frame(well=wells, DIP.rate=rates, cell.line=cl, expt.date=expt)
  rownames(out) <- NULL
  out
}

source('log2norm.R')
convertWell <- function(dat) {
  well.temp <- unlist(strsplit(dat$Well,'-'))
  row.temp <- LETTERS[as.numeric(substr(well.temp[grep('R',well.temp)],2,3))]
  col.temp <- substr(well.temp[grep('C',well.temp)],2,3)
  dat$Well     <- paste0(row.temp,col.temp)
  dat
}

getDIPs <- function(dat, maxconc) {
  # dat_c <- convertWell(dat)
  dat_c <- dat
  dat_c$nl2 <- log2norm(dat_c$Count, dat_c$Well)
  dat_c$cell.line <- "SKMEL5"
  dat_c$expt.date <- "07-30-2018"
  dat_c$Row <- NULL
  dat_c$Column <- NULL
  names(dat_c) <- c("time", "well", "cell.count", 
                    "drug1", "drug1.conc", "drug1.units",
                    "nl2", "cell.line", "expt.date")
  dat_c_dips <- getWellRates(dat_c, time.range = c(50, 140))
  dat_c_dips$drug.conc <- rep(c(maxconc,maxconc/2,maxconc/4,maxconc/8,maxconc/16,
                                maxconc/32,maxconc/64,maxconc/128,maxconc/256,0),
                              times = 6)
  dat_c_dips$rep <- rep(c(1,2,3), each = 10, times = 2)
  dat_c_dips
}

getDRdata <- function(dat, maxconc, drug2, State, dilution) {
  data <- getDIPs(dat, maxconc)
  if (drug2 == "Y") {
    data <- data[31:60,]
  } else {
    data <- data[1:30,]
  }
  data$State <- State
  dataQC <- subset(data, data[,"DIP.rate"] < 0.05)
  dataQC$drug.conc0 <- dataQC$drug.conc
  n <- dataQC$drug.conc0
  dataQC$drug.conc0[dataQC$drug.conc0 == 0] <- min(n[n!=min(n)])/dilution
  dataQC
}

getDRC <- function(dat, State, drug2) {
  require(drc)
  model <- drm(DIP.rate ~ drug.conc0, data = dat, fct = LL.4())
  newdata <- expand.grid(conc=exp(seq(log(min(dat$drug.conc0)), 
                                      log(max(dat$drug.conc0)),
                                      length = 100)))
  pm <- predict(model, newdata=newdata, interval="confidence")
  newdata$p <- pm[,1]
  newdata$pmin <- pm[,2]
  newdata$pmax <- pm[,3]
  newdata$State <- State
  newdata$drug2 <- drug2
  newdata
}

getIC50 <- function(dat, State, drug2) {
  require(drc)
  model <- drm(DIP.rate ~ drug.conc0, 
               data = dat, 
               fct = LL.4(names = c("hill", "min_value", "max_value", "ec_50")))
  
  coefs <- setNames(
    model$coefficients,
    c("hill", "min_value", "max_value", "ec_50")
  )
  
  ic_50 <- with(
    as.list(coefs),
    exp(
      log(ec_50) + (1 / hill) * log(max_value / (max_value - 2 * min_value))
    )
  )
  ic_50
}

### RSL3
RSL3_DR_Untreated_noFer <- getDRdata(dat = pRSL3, maxconc=10, drug2 = "N", 
                                     State = "Untreated", dilution = 2)
RSL3_DR_Untreated_Fer <- getDRdata(dat = pRSL3, maxconc=10, drug2 = "Y", 
                                   State = "Untreated", dilution = 2)
RSL3_DR_Idling_noFer <- getDRdata(dat = iRSL3, maxconc=10, drug2 = "N", 
                                  State = "Idling", dilution = 2)
RSL3_DR_Idling_Fer <- getDRdata(dat = iRSL3, maxconc=10, drug2 = "Y", 
                                State = "Idling", dilution = 2)
####
RSL3_fit_Untreated_noFer <- getDRC(dat = RSL3_DR_Untreated_noFer, 
                                   State = "Untreated", drug2 = "No")
RSL3_fit_Untreated_Fer <- getDRC(dat = RSL3_DR_Untreated_Fer, 
                                 State = "Untreated", drug2 = "Yes")
RSL3_fit_Idling_noFer <- getDRC(dat = RSL3_DR_Idling_noFer, 
                                State = "Idling", drug2 = "No")
RSL3_fit_Idling_Fer <- getDRC(dat = RSL3_DR_Idling_Fer, 
                              State = "Idling", drug2 = "Yes")
###
IC50_Untreated_noFer <- getIC50(dat = RSL3_DR_Untreated_noFer, 
                                   State = "Untreated", drug2 = "No")
IC50_Untreated_Fer <- getIC50(dat = RSL3_DR_Untreated_Fer, 
                                 State = "Untreated", drug2 = "Yes")
IC50_Idling_noFer <- getIC50(dat = RSL3_DR_Idling_noFer, 
                                State = "Idling", drug2 = "No")
IC50_Idling_Fer <- getIC50(dat = RSL3_DR_Idling_Fer, 
                              State = "Idling", drug2 = "Yes")

IC50_df <- data.frame(IC50 = c(IC50_Untreated_noFer, IC50_Untreated_Fer,
                               IC50_Idling_noFer, IC50_Idling_Fer),
                      State = c("Untreated", "Untreated", "Idling", "Idling"),
                      drug2 = c("No", "Yes", "No", "Yes"))



####
RSL3_UTvI <- rbind(RSL3_fit_Untreated_noFer, RSL3_fit_Idling_noFer)
RSL3_UTvI_Fer <- rbind(RSL3_fit_Untreated_Fer, RSL3_fit_Idling_Fer)
RSL3_UTvI_all <- rbind(RSL3_fit_Untreated_noFer, RSL3_fit_Idling_noFer,
                   RSL3_fit_Untreated_Fer, RSL3_fit_Idling_Fer)


RSL3_UTvI$State <- factor(RSL3_UTvI$State, levels = c("Untreated", "Idling"))
# ggplot(RSL3_UTvI, aes(x=conc, y=p, color=State, fill=State)) +
#   geom_ribbon(aes(ymin=pmin, ymax=pmax), alpha=0.2, color = NA) +
#   geom_line(size = 1) +
#   geom_vline(aes(xintercept=IC50_df$IC50[1]), color = "red", linetype = "dashed", size = 1) +
#   geom_vline(aes(xintercept=IC50_df$IC50[3]), color = "blue", linetype = "dashed", size = 1) +
#   # coord_trans(x="log10") +
#   scale_x_log10() +
#   xlab("[RSL3] (uM)") + ylab("DIP Rate") +
#   theme_bw() + ggtitle("RSL3") +
#   scale_color_manual(values = c("red", "blue")) +
#   scale_fill_manual(values = c("red", "blue")) +
#   theme(legend.text = element_text(size = 12), 
#         plot.title = element_text(size = 14, hjust = 0.5), 
#         axis.text=element_text(size=12),
#         legend.title = element_text(size=12),
#         legend.position = "none",
#         axis.title=element_text(size=12),
#         panel.grid = element_blank()) +
#   ggsave("RSL3_DRC_UTvI_2021vis.pdf", width = 4, height = 3)

RSL3_UTvI_Fer$State <- factor(RSL3_UTvI$State, levels = c("Untreated", "Idling"))
# ggplot(RSL3_UTvI_Fer, aes(x=conc, y=p, color=State, fill=State)) +
#   geom_ribbon(aes(ymin=pmin, ymax=pmax), alpha=0.2, color = NA) +
#   geom_line(size = 1) +
#   geom_vline(aes(xintercept=IC50_df$IC50[2]), color = "red", linetype = "dashed", size = 1) +
#   geom_vline(aes(xintercept=IC50_df$IC50[4]), color = "blue", linetype = "dashed", size = 1) +
#   scale_x_log10() +
#   xlab("[RSL3] (uM)") + ylab("DIP Rate") +
#   theme_bw() + ggtitle("RSL3 + Fer-1") +
#   scale_color_manual(values = c("red", "blue")) +
#   scale_fill_manual(values = c("red", "blue")) +
#   theme(legend.text = element_text(size = 12), 
#         plot.title = element_text(size = 14, hjust = 0.5), 
#         axis.text=element_text(size=12),
#         legend.title = element_text(size=12),
#         legend.position = "none",
#         axis.title=element_text(size=12),
#         panel.grid = element_blank()) +
#   ggsave("RSL3-Fer1_DRC_UTvI_2021vis.pdf", width = 4, height = 3)

##
RSL3_UTvI_all$State <- factor(RSL3_UTvI_all$State, levels = c("Untreated", "Idling"))
ggplot(RSL3_UTvI_all, aes(x=conc, y=p, color=State, fill=State)) +
  facet_wrap(~drug2) +
  geom_ribbon(aes(ymin=pmin, ymax=pmax), alpha=0.2, color = NA) +
  geom_line(size = 1) +
  geom_vline(data = IC50_df, aes(xintercept = IC50, color = State), linetype = "dashed", size = 0.75) +
  scale_x_log10() +
  xlab("[RSL3] (uM)") + ylab("DIP Rate") +
  theme_bw() + 
  scale_color_manual(values = c("red", "blue")) +
  scale_fill_manual(values = c("red", "blue")) +
  theme(legend.text = element_text(size = 12), 
        plot.title = element_text(size = 14, hjust = 0.5), 
        axis.text=element_text(size=12),
        legend.title = element_text(size=12),
        legend.position = "none",
        axis.title=element_text(size=12),
        panel.grid = element_blank()) 
  ggsave("RSL3-all_DRC_UTvI_2021vis.pdf", width = 8, height = 3)

