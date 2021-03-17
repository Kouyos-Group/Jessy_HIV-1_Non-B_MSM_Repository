
### Clear R's memory
rm(list=ls())


### get the packages
if (!requireNamespace("pacman", quietly = TRUE))
  install.packages("pacman")
pacman::p_load("dplyr","tidyr", "data.table", "gmodels", "tidyverse", "boot", "table1","tableone", "devtools", "survival", "vcd")




#FastTree
setwd("~")
setwd("Downloads/MSMHET_MTCS_above50") #NEWNEW: with MSM>50% and GD: 4.5%, BS: NO,>80%
MSMBoottable <- read.csv("MSMcluster_allsubtypes_above50_withboot.csv")
HETBoottable <- read.csv("HETcluster_allsubtypes_above50_withboot.csv")
setwd("~/Desktop/SHCS")
data<- read.csv("Output/fulltable_study_population_update.csv")
data$X <-  NULL #does the same: data <- data[,-1]
setwd("~")
setwd("Downloads/MSMHET_MTCS_above50")
for (i in c("A","C","G","F","CRF01_AE","CRF02_AG")){
  assign(paste0("dataframe_",i),read.csv(file=paste0(paste0("MSMHETcluster_GD45B80_",i),".csv"))[,-1])
  dd <-eval(parse(text=paste0("dataframe_",i)))
  dd$Subtype <- i
  assign(paste0("MSM_",i),dd[dd$MSMQ == T,])
  MSM <-eval(parse(text=paste0("MSM_",i)))
  assign(paste0("largeCL_",i),MSM[ MSM$length >10,])
  LargeCluster <-  eval(parse(text=paste0("largeCL_",i)))
  #assign(paste0("largeCL_",i),LargeCluster["Subtype"==i,])
  assign(paste0("largeCL_",i),LargeCluster[with(LargeCluster, order(length, decreasing = F)),])
  # dataframe <- eval(parse(text=paste0("largeCL_",i)))
 
}

ALL <- rbind(largeCL_A,largeCL_C,largeCL_F,largeCL_G,largeCL_CRF01_AE,largeCL_CRF02_AG)
ALL_N1 <- ALL[ALL$Subtype =="F",]; ALL_N1$clusterNR <- 3
ALL_N3 <- ALL[ALL$Subtype =="CRF02_AG",]; ALL_N3$clusterNR <- 4
ALL_N2_N4 <- ALL[ALL$Subtype =="CRF01_AE",]
ALL_new <- rbind(ALL_N1,ALL_N3, ALL_N2_N4)

data$var_eth
ALL_new <- merge(ALL_new, data[,c('id', 'diag_year', 'var_eth','regroupedregion')], all.x = T)

colnames(ALL_new)

rndr <- function(x, name, ...) {
  if (!is.numeric(x)) return(render.categorical.default(x))
  what <- switch(name,
                 gender ="",
                 riskNEW = "",
                 ethNEW = "",
                 var_eth = "",
                 regionNEW = "",
                 diag_year = "Median [Q1,Q3]",
                 edu = "")
  #age_sampledate = "Median [Min,Max]")
  parse.abbrev.render.code(c("", what))(x)
}


## better labels and units:

label(ALL_new$riskNEW) <- "Risk Group"
label(ALL_new$clusterType) <- "Cluster type"
label(ALL_new$regionNEW) <-"Geographic region of origin"
label(ALL_new$ethNEW) <- "Ethnicity"
label(ALL_new$var_eth) <- "Ethnicity"
ALL_new$diag_year <- as.integer(ALL_new$diag_year)
label(ALL_new$diag_year) <- "Year of diagnosis"



table1(~ 
         Subtype+
         clusterType+
         riskNEW+
         #ethNEW+
         var_eth+
         regionNEW +
         regroupedregion+
         diag_year
       | clusterNR , 
       data = ALL_new,
       render=rndr, 
       rowlabelhead = "Cluster",
       overall = T,
       footnote = "<b>NOTE</b> All numbers are no. (%) unless otherwise stated<br/>",
       caption = "<h3><b> Table 1. D45 BS80 </b></h3>")




