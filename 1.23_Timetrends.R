### MSM Time series analysis

## Clear R's memory
rm(list=ls())

# Libraries
library(haven)
library(dygraphs) # allows to represent time series: a chart where the X axis represent time, and the Y axis the evolution of one or several variables.
library(xts) # To make the convertion data-frame / xts format
library(lubridate)
library(tidyverse)

## load data
setwd("~/Desktop/SHCS")
data<- read.csv("Output/fulltable_study_population.csv")
data$X <-  NULL #does the same: data <- data[,-1]


## Dates as dates
data$sampledate<- as.Date(as.character(data$sampledate), tryFormats = c("%Y-%m-%d"))
data$seq_dt<- as.Date(as.character(data$seq_dt), tryFormats = c("%Y-%m-%d"))
data$exitdate<- as.Date(as.character(data$exitdate), tryFormats = c("%Y-%m-%d"))
data$regdate<- as.Date(as.character(data$regdate), tryFormats = c("%Y-%m-%d"))
data$hiv_posdate<- as.Date(as.character(data$hiv_posdate), tryFormats = c("%Y-%m-%d"))
data$hiv_posdocdate<- as.Date(as.character(data$hiv_posdocdate), tryFormats = c("%Y-%m-%d"))
data$hiv_negdate<- as.Date(as.character(data$hiv_negdate), tryFormats = c("%Y-%m-%d"))
data$born<- as.Date(as.character(data$born), tryFormats = c("%Y"))
data$diag_date<- as.Date(as.character(data$diag_date), tryFormats = c("%Y-%m-%d"))
data$diag_year  <- substr(data$diag_date,1,4)
data$diag_year<- as.Date(as.character(data$diag_year), tryFormats = c("%Y"))



### check if time is indeed recognized as a date:
str(data$diag_year)
str(data$diag_date)


data<- data %>%
  mutate(copy_subtype=case_when(
    data$subtype=="B"  ~ "B",
    #~ "non-B"
  )) #20,802
select <- which(is.na(data$copy_subtype))
data$copy_subtype[select] <- "Non-B"
FUN_nonB <- function(x){
  if (is.na(x)){return (0)}
  else if (x == "B" ){return(0)}
  else if (x == "Non-B"){return(1)}
  else{return(0)}
}
data$nonB <- sapply(data$copy_subtype, FUN_nonB)



### need to use another variable, which is the same for all the entries! 
data$all <- "1"
data$all <- as.numeric(data$all)


# MSM!
MSM <- data[data$risk == "1",] #4751



summary(MSM$diag_date)
# Min.      1st Qu.       Median         Mean      3rd Qu.         Max.         NA's 
# "1982-07-01" "1993-07-15" "2001-07-01" "2001-02-23" "2008-06-21" "2020-03-18"        "400" 



MSM_All <- setNames(aggregate(MSM$all~MSM$diag_year, FUN = sum),c("diag_year","incidence_all"))
MSM_All$prev_all <- MSM_All$incidence_all
## to get the prevelance: 
for (i in 2:nrow(MSM_All)) {
  MSM_All[i,3] <- as.numeric(MSM_All[i-1,3]) + as.numeric(MSM_All[i,2])
}
# Compute % Returns
MSM_All$returns_perc_p <- (MSM_All$prev_all)/max(MSM_All$prev_all)
MSM_All$returns_perc_i <- (MSM_All$incidence_all)/max(MSM_All$incidence_all)


# Format 3: Several variables for each date
dd3 <- data.frame(
  time= MSM_All$diag_year, 
  value1=MSM_All$incidence_all
)

# Then you can create the xts format:
don3=xts( x=dd3[,-1], order.by=dd3$time)
p3 <- dygraph(don3);p3# Chart


##################

MSM_B_NB <- setNames(aggregate(MSM$all~MSM$diag_year + MSM$copy_subtype, FUN = sum),c("diag_year","Subtype","incidence_all"))
MSM_B_NB$prev_all <- MSM_B_NB$incidence_all
MSM_B_NB <- MSM_B_NB[MSM_B_NB$diag_year != "1981-05-13" &MSM_B_NB$diag_year != "1982-05-13",] 

## to get the prevelance: 
for (i in 2:nrow(MSM_B_NB)) {
  if (B_NB[i,2] == MSM_B_NB[i-1,2]) {
    B_NB[i,4] <- as.numeric(MSM_B_NB[i-1,4]) + as.numeric(MSM_B_NB[i,3])
  }}


MSM_B_NB$returns_perc_p <- (MSM_B_NB$prev_all)/max(MSM_B_NB$prev_all)
MSM_B_NB$returns_perc_i <- (MSM_B_NB$incidence_all)/max(MSM_B_NB$incidence_all)



MSM_B <- MSM_B_NB[MSM_B_NB$Subtype =="B",]
MSM_NB <- MSM_B_NB[MSM_B_NB$Subtype =="Non-B",]
MSM_DD <- merge(MSM_B,MSM_NB, by="diag_year")

MSM_DD$returns_perc_p_B <- (MSM_DD$prev_all.x)/(MSM_DD$prev_all.x+MSM_DD$prev_all.y)
MSM_DD$returns_perc_p_NB <- (MSM_DD$prev_all.y)/(MSM_DD$prev_all.x+MSM_DD$prev_all.y)

dm <- data.frame(
  time= MSM_DD$diag_year, 
  "B"=MSM_DD$returns_perc_p_B,
  "Non-B"=MSM_DD$returns_perc_p_NB
)


dm_nB <- data.frame(
  time= MSM_DD$diag_year, 
  value1=MSM_DD$returns_perc_p_NB
)

MSM_don2<-xts(x=dm[,-1], order.by=dm$time)
MSM_don3<-xts(x=dm_nB[,-1], order.by=dm_nB$time)
dygraph(MSM_don2) 
dygraph(MSM_don3) 

dyBarChart <- function(dygraph) {
  dyPlotter(dygraph = dygraph,
            name = "BarChart",
            path = system.file("plotters/barchart.js",
                               package = "dygraphs"))
}


dyMultiColumn <- function(dygraph) {
  dyPlotter(dygraph = dygraph,
            name = "MultiColumn",
            path = system.file("plotters/multicolumn.js",
                               package = "dygraphs"))
}

MSM_BETTER<-dygraph(MSM_don2, main= "Time trends in distribution of subtypes in the SHCS among MSM",
                    ylab = "Percent",
                    xlab= "Year of HIV Diagnosis")  %>%
  dySeries("B", label = "B",  strokeWidth = 2) %>%
  dySeries("Non.B", label = "Non-B",  strokeWidth = 2) %>%
  dyOptions(stackedGraph = TRUE, colors= c( "blue","magenta"),
            fillGraph = TRUE, fillAlpha = 0.4) %>%
  dyRangeSelector(height = 40)  %>%
  dyHighlight(highlightSeriesOpts = list(strokeWidth = 5))%>%
  dyLegend(width = 450);MSM_BETTER




MSM_BETTER2<-dygraph(MSM_don3, main= "Ratio of Non-B Subtype in the SHCS among MSM",
                     ylab = "Percent",
                     xlab= "Year of HIV Diagnosis")  %>%
  dySeries("V1", label = "Non-B",  strokeWidth = 2) %>%
  dyOptions( colors=  "magenta",
             fillGraph = TRUE, fillAlpha = 0.4) %>%
  dyLegend(width = 450);MSM_BETTER2


MSM_BETTER_bar<-dygraph(MSM_don2, main= "Time trends in distribution of subtypes in the SHCS among MSM",
                     ylab = "Percent",
                     xlab= "Year of HIV Diagnosis")  %>%
  dyStackedBarChart() %>%
  dyOptions(stackedGraph = TRUE, colors= c( "blue","magenta"),
            fillGraph = TRUE, fillAlpha = 0.4) %>%
  dyRangeSelector(height = 40, dateWindow = c("1988-06-30", "2020-06-30"))%>%
  dyLegend(width = 450);MSM_BETTER_bar



library(ggplot2)
library(ggthemes)
library(extrafont)
library(plyr)
library(scales)
dd <- ddply(MSM_B_NB, .(diag_year),
                     transform, pos = cumsum(perc) - (0.5 *perc))

dd$perc <- round(dd$returns_perc_p*100,1)
p4 <- ggplot() + 
  geom_bar(aes(y = perc, x = diag_year, fill = Subtype), data = dd,
           stat="identity")  
p4 <- p4 + geom_text(data=dd, aes(x = diag_year, y = pos, label = paste0(perc,"%")),
                     size=4)
p4




########
### same for non-Bs
NB_MSM<- filter(MSM, subtype=="A" | subtype=="02_AG" 
                   | subtype=="C" | subtype=="01_AE"
                | subtype=="F"| subtype=="G") ;unique(NB_MSM$subtype) #310

NB <- setNames(aggregate(NB_MSM$all~NB_MSM$diag_year + NB_MSM$subtype, FUN = sum),c("diag_year","Subtype","incidence_all"))
NB$prev_all <- NB$incidence_all
#NB <- NB[B_NB$diag_year != "1981-05-13" &NB$diag_year != "1982-05-13",] 

## to get the prevelance: 
for (i in 2:nrow(NB)) {
  if (NB[i,2] == NB[i-1,2]) {
    NB[i,4] <- as.numeric(NB[i-1,4]) + as.numeric(NB[i,3])
  }}

# Compute % Returns
NB$returns_perc_p <- (NB$prev_all)/max(NB$prev_all)
NB$returns_perc_i <- (NB$incidence_all)/max(NB$incidence_all)


NB2 <- NB






### percentage 
A <- NB[NB$Subtype =="A",]
C <-  NB[NB$Subtype =="C",]
AE <- NB[NB$Subtype =="01_AE",]
AG <-  NB[NB$Subtype =="02_AG",]
F1 <-  NB[NB$Subtype =="F",]
G <-  NB[NB$Subtype =="G",]
all_NB <- Reduce(function(x, y) merge(x, y, all=TRUE, by="diag_year"), list(A, C, AE,AG, F1,G))
colnames(all_NB) <- c("diag_year", "Subtype_A","inc_A","prev_A",
                      "perc_p_A", "perc_i_A",
                      "Subtype_C","inc_C","prev_C",
                      "perc_p_C", "perc_i_C",
                      "Subtype_AE","inc_AE","prev_AE",
                      "perc_p_AE", "perc_i_AE",
                      "Subtype_AG","inc_AG","prev_AG",
                      "perc_p_AG", "perc_i_AG",
                      "Subtype_F","inc_F","prev_F",
                      "perc_p_F", "perc_i_F",
                      "Subtype_G","inc_G","prev_G",
                      "perc_p_G", "perc_i_G")
all_NB$Subtype_A <- "A"; all_NB$Subtype_C <- "C"; all_NB$Subtype_AG <- "AG"; all_NB$Subtype_AE <- "AE"; all_NB$Subtype_AG <- "F"; all_NB$Subtype_AE <- "G"
all_NB[is.na(all_NB)] <- 0; head(all_NB)

sum_all_NB <-all_NB$perc_p_A+all_NB$perc_p_AE+all_NB$perc_p_C+all_NB$perc_p_AG+all_NB$perc_p_F+all_NB$perc_p_G
all_NB$return_perc_p_A <- all_NB$perc_p_A/sum_all_NB;
all_NB$return_perc_p_C <- all_NB$perc_p_C/sum_all_NB;
all_NB$return_perc_p_AG <- all_NB$perc_p_AG/sum_all_NB;
all_NB$return_perc_p_AE <- all_NB$perc_p_AE/sum_all_NB;
all_NB$return_perc_p_F <- all_NB$perc_p_F/sum_all_NB;
all_NB$return_perc_p_G <- all_NB$perc_p_G/sum_all_NB;

df_NB <- data.frame(
  time= all_NB$diag_year, 
  "A"=all_NB$return_perc_p_A ,
  "AG"=all_NB$return_perc_p_AG ,
  "C"=all_NB$return_perc_p_C ,
  "F"=all_NB$return_perc_p_F,
  "G"=all_NB$return_perc_p_G,
  "AE"=all_NB$return_perc_p_AE,
)
# Then you can create the xts format:
don_NB<-xts(x=df_NB[,-1], order.by=df_NB$time)
dygraph(don_NB) # Chart


gyg_NB<-dygraph(don_NB, main= "Time trends in distribution of non-B subtypes among MSM",
                ylab = "Percent",
                xlab= "Year of HIV Diagnosis")  %>%
  dySeries("AG", label = "AG",  strokeWidth = 2) %>%
  dySeries("C", label = "C",  strokeWidth = 2) %>%
  dySeries("AE", label = "AE",  strokeWidth = 2) %>%
  dySeries("A", label = "A",  strokeWidth = 2) %>%
  dySeries("G", label = "G",  strokeWidth = 2) %>%
  dyOptions(stackedGraph = TRUE,
            fillGraph = TRUE, fillAlpha = 0.4) %>%
  dyRangeSelector(height = 40)  %>%
  dyHighlight(highlightSeriesOpts = list(strokeWidth = 5))%>%
  dyLegend(width = 450); gyg_NB

MSM_NONB_bar<-dygraph(don_NB, main= "Time trends in distribution of non-B subtypes in the SHCS among MSM",
                        ylab = "Percent",
                        xlab= "Year of HIV Diagnosis")  %>%
  dySeries("A", label = "A",  strokeWidth = 2) %>%
  dySeries("AE", label = "AE",  strokeWidth = 2) %>%
  dyStackedBarChart() %>%
  dyOptions(stackedGraph = T,colors= c( "orange","blue","brown","darkblue", "gray", "green")) %>%
  dyRangeSelector(height = 40, dateWindow = c("1992-06-30", "2020-06-30"))%>%
  dyLegend(width = 450);MSM_NONB_bar


