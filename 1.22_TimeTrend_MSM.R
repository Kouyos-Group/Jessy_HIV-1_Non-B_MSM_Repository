
### MSM Time series analysis

## Clear R's memory
rm(list=ls())

# Libraries
library(haven)
library(dygraphs) # allows to represent time series: a chart where the X axis represent time, and the Y axis the evolution of one or several variables.
library(xts) # To make the convertion data-frame / xts format
library(lubridate)
library(tidyverse)
library(ggplot2)

library(grid)
library(ggpubr)

## load data
setwd("~/Desktop/SHCS")
# data<- read.csv("Output/fulltable_study_population.csv")
data<- read.csv("Output/fulltable_study_population_update.csv")
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
MSM <- data[which(data$risk == "1"),] #4351
  # HET <- data[which(data$risk == "2"),] #HET
# IDU <- data[which(data$risk == "3"),] #IDU
MSM <- MSM[!(is.na(MSM$id)),] #4351

summary(MSM$diag_date)
# Min.      1st Qu.       Median         Mean      3rd Qu.         Max.         NA's 
# "1982-07-01" "1993-07-15" "2001-07-01" "2001-02-23" "2008-06-21" "2020-03-18"        "400" 



##################

MSM_B_NB <- setNames(aggregate(MSM$all~MSM$diag_year + MSM$copy_subtype, FUN = sum),c("diag_year","Subtype","incidence_all"))
MSM_B_NB$prev_all <- MSM_B_NB$incidence_all
MSM_B_NB <- MSM_B_NB[MSM_B_NB$diag_year != "1981-05-13" &MSM_B_NB$diag_year != "1982-05-13",] 

## to get the prevelance: 
for (i in 2:nrow(MSM_B_NB)) {
  if (MSM_B_NB[i,2] == MSM_B_NB[i-1,2]) {
    MSM_B_NB[i,4] <- as.numeric(MSM_B_NB[i-1,4]) + as.numeric(MSM_B_NB[i,3])
  }}


# MSM_B_NB$returns_perc_p <- (MSM_B_NB$prev_all)/max(MSM_B_NB$prev_all)
# MSM_B_NB$returns_perc_i <- (MSM_B_NB$incidence_all)/max(MSM_B_NB$incidence_all)



MSM_B <- MSM_B_NB[MSM_B_NB$Subtype =="B",]
MSM_NB <- MSM_B_NB[MSM_B_NB$Subtype =="Non-B",]
MSM_DD <- merge(MSM_B,MSM_NB, by="diag_year")

MSM_DD$returns_perc_p_B <- (MSM_DD$prev_all.x)/(MSM_DD$prev_all.x+MSM_DD$prev_all.y)
MSM_DD$returns_perc_p_NB <- (MSM_DD$prev_all.y)/(MSM_DD$prev_all.x+MSM_DD$prev_all.y)

MSM_DD$returns_perc_inc_B <- (MSM_DD$incidence_all.x)/(MSM_DD$incidence_all.x+MSM_DD$incidence_all.y)
MSM_DD$returns_perc_inc_NB <- (MSM_DD$incidence_all.y)/(MSM_DD$incidence_all.x+MSM_DD$incidence_all.y)

MSM_6NB<- filter(MSM, subtype=="A" | subtype=="02_AG" 
                 | subtype=="C" | subtype=="01_AE"
                 | subtype=="F"| subtype=="G"| subtype=="B") ;unique(MSM_6NB$subtype) #5114




MSM_B_6NB <- setNames(aggregate(MSM_6NB$all~MSM_6NB$diag_year + MSM_6NB$copy_subtype, FUN = sum),c("diag_year","Subtype","incidence_all"))
MSM_B_6NB$prev_all <- MSM_B_6NB$incidence_all
MSM_B_6NB <- MSM_B_6NB[MSM_B_6NB$diag_year != "1981-05-13" &MSM_B_6NB$diag_year != "1982-05-13",] 

## to get the prevelance: 
for (i in 2:nrow(MSM_B_6NB)) {
  if (MSM_B_6NB[i,2] == MSM_B_6NB[i-1,2]) {
    MSM_B_6NB[i,4] <- as.numeric(MSM_B_6NB[i-1,4]) + as.numeric(MSM_B_6NB[i,3])
  }}


# MSM_B_6NB$returns_perc_p <- (MSM_B_6NB$prev_all)/max(MSM_B_6NB$prev_all)
# MSM_B_6NB$returns_perc_i <- (MSM_B_6NB$incidence_all)/max(MSM_B_6NB$incidence_all)
# 


MSM1_B <- MSM_B_6NB[MSM_B_6NB$Subtype =="B",]
MSM1_NB <- MSM_B_6NB[MSM_B_6NB$Subtype =="Non-B",]
MSM1_DD <- merge(MSM1_B,MSM1_NB, by="diag_year", all.x = T)
library("imputeTS")
MSM1_DD<- na.replace(MSM1_DD, 0); MSM1_DD$Subtype.y <- "Non-B"

MSM1_DD$returns_perc_p_B <- (MSM1_DD$prev_all.x)/(MSM1_DD$prev_all.x+MSM1_DD$prev_all.y)
MSM1_DD$returns_perc_p_NB <- (MSM1_DD$prev_all.y)/(MSM1_DD$prev_all.x+MSM1_DD$prev_all.y)

MSM1_DD$returns_perc_inc_B <- (MSM1_DD$incidence_all.x)/(MSM1_DD$incidence_all.x+MSM1_DD$incidence_all.y)
MSM1_DD$returns_perc_inc_NB <- (MSM1_DD$incidence_all.y)/(MSM1_DD$incidence_all.x+MSM1_DD$incidence_all.y)



### For the 6 Non-B subtypes: 

################################
#### Prevalence
################################

#Get the dataframes:
MSM_6PREV_N <- rbind(
  data.frame("diag_year"= MSM1_DD$diag_year, "prevalence" =MSM1_DD$prev_all.x, "type"="B"),
  data.frame("diag_year"= MSM1_DD$diag_year, "prevalence" =MSM1_DD$prev_all.y, "type"="Non-B")
)


MSM_6PREV_FRACTION <- rbind(
  data.frame("diag_year"= MSM1_DD$diag_year, "perc_prev" =MSM1_DD$returns_perc_p_B, "type"="B"),
  data.frame("diag_year"= MSM1_DD$diag_year, "perc_prev" =MSM1_DD$returns_perc_p_NB, "type"="Non-B")
)




################################
## Stacked Barplots
## Plot the prevalence in absolute numbers:
P_6N <- ggplot(MSM_6PREV_N, aes(x=diag_year,y=prevalence, fill= type)) +
  geom_bar(position="stack", stat="identity") +
  scale_x_date(limits = as.Date(c("1990-01-01","2020-01-01"))) +
  scale_fill_manual(values=c("lightblue","darkred")) +
  labs(x="Year of HIV Diagnosis",  y="Frequency")+
  guides(fill=guide_legend(title = "Subtype", reverse = F, size= 14))+
  theme_bw() + 
  theme(
    plot.title = element_text(color="black", size=14, face="bold.italic"),
    axis.title.x = element_text(color="black", size=14, face="bold"),
    axis.title.y = element_text(color="black", size=14, face="bold"),
    #legend.position = c(0.9, 0.8),
    legend.text = element_text(colour="black", size=14, face="plain"),
    legend.title = element_text(colour="black", size=14, face="bold")) +
  geom_hline(yintercept=0, linetype="dashed", color = "black")
P_6N

## Plot the prevalence as Fraction:
P_6F <-ggplot(MSM_6PREV_FRACTION, aes(x=diag_year,y=perc_prev, fill= type)) +
  geom_bar(position="stack", stat="identity") +
  scale_x_date(limits = as.Date(c("1990-01-01","2020-01-01"))) +
  scale_y_continuous(breaks = c(seq(0,1,0.1)))+
  scale_fill_manual(values=c("lightblue","darkred")) +
  labs(x="Year of HIV Diagnosis",  y="Fraction")+
  guides(fill=guide_legend(title = "Subtype", reverse = F, size= 14))+
  theme_bw() + 
  theme(
    plot.title = element_text(color="black", size=14, face="bold.italic"),
    axis.title.x = element_text(color="black", size=14, face="bold"),
    axis.title.y = element_text(color="black", size=14, face="bold"),
    #legend.position = c(0.9, 0.8),
    legend.text = element_text(colour="black", size=14, face="plain"),
    legend.title = element_text(colour="black", size=14, face="bold")) +
  geom_hline(yintercept=0, linetype="dashed", color = "black")
P_6F




################################
#### INCIDENCE
################################

#Get the dataframes:
MSM_6INC_N <- rbind(
  data.frame("diag_year"= MSM1_DD$diag_year, "incidence" =MSM1_DD$incidence_all.x, "type"="B"),
  data.frame("diag_year"= MSM1_DD$diag_year, "incidence" =MSM1_DD$incidence_all.y, "type"="Non-B")
)


MSM_6INC_FRACTION <- rbind(
  data.frame("diag_year"= MSM1_DD$diag_year, "perc_inc" =MSM1_DD$returns_perc_inc_B, "type"="B"),
  data.frame("diag_year"= MSM1_DD$diag_year, "perc_inc" =MSM1_DD$returns_perc_inc_NB, "type"="Non-B")
)


################################
## Stacked Barplots

## Plot the Incidence in absolute numbers:
I_6N <- ggplot(MSM_6INC_N, aes(x=diag_year,y=incidence, fill= type)) +
  geom_bar(position="stack", stat="identity") +
  scale_x_date(limits = as.Date(c("1990-01-26","2020-01-26"))) +
  scale_fill_manual(values=c("lightblue","darkred")) +
  labs(x="Year of HIV Diagnosis",  y="Frequency")+
  guides(fill=guide_legend(title = "Subtype", reverse = F, size= 14))+
  theme_bw() + 
  theme(
    plot.title = element_text(color="black", size=14, face="bold.italic"),
    axis.title.x = element_text(color="black", size=14, face="bold"),
    axis.title.y = element_text(color="black", size=14, face="bold"),
    #legend.position = c(0.9, 0.8),
    legend.text = element_text(colour="black", size=14, face="plain"),
    legend.title = element_text(colour="black", size=14, face="bold")) +
  geom_hline(yintercept=0, linetype="dashed", color = "black")
I_6N

## Plot the Incidence as Fraction:
I_6F <-ggplot(MSM_6INC_FRACTION, aes(x=diag_year,y=perc_inc, fill= type)) +
  geom_bar(position="stack", stat="identity") +
  scale_x_date(limits = as.Date(c("1990-01-01","2020-01-01"))) +
  scale_y_continuous(breaks = c(seq(0,1,0.1)))+
  scale_fill_manual(values=c("lightblue","darkred")) +
  labs(x="Year of HIV Diagnosis",  y="Fraction")+
  guides(fill=guide_legend(title = "Subtype", reverse = F, size= 14))+
  theme_bw() + 
  theme(
    plot.title = element_text(color="black", size=14, face="bold.italic"),
    axis.title.x = element_text(color="black", size=14, face="bold"),
    axis.title.y = element_text(color="black", size=14, face="bold"),
    #legend.position = c(0.9, 0.8),
    legend.text = element_text(colour="black", size=14, face="plain"),
    legend.title = element_text(colour="black", size=14, face="bold")) +
  geom_hline(yintercept=0, linetype="dashed", color = "black")
I_6F





ggarrange(I_6N, P_6N, I_6F, P_6F, ncol = 2, nrow = 2,
          labels = c("A", "B","C","D"), 
          common.legend = TRUE, legend = "bottom",
          font.label = list(size = 18, color = "black"))


setEPS()
postscript("Trend_MSM_6NonB.eps")
ggarrange(I_6N, I_6F, 
          ncol = 1, nrow = 2,
          labels = c("A", "B"), 
          common.legend = TRUE, legend = "bottom",
          font.label = list(size = 18, color = "black"))
dev.off()



################################
################################
################################




### For the ALL  Non-B subtypes: 


#### Prevalence


#Get the dataframes:
MSM_PREV_N <- rbind(
  data.frame("diag_year"= MSM_DD$diag_year, "prevalence" =MSM_DD$prev_all.x, "type"="B"),
  data.frame("diag_year"= MSM_DD$diag_year, "prevalence" =MSM_DD$prev_all.y, "type"="Non-B")
)


MSM_PREV_FRACTION <- rbind(
  data.frame("diag_year"= MSM_DD$diag_year, "perc_prev" =MSM_DD$returns_perc_p_B, "type"="B"),
  data.frame("diag_year"= MSM_DD$diag_year, "perc_prev" =MSM_DD$returns_perc_p_NB, "type"="Non-B")
)




################################
## Stacked Barplots
## Plot the prevalence in absolute numbers:
P_N <- ggplot(MSM_PREV_N, aes(x=diag_year,y=prevalence, fill= type)) +
  geom_bar(position="stack", stat="identity") +
  scale_x_date(limits = as.Date(c("1990-01-01","2020-01-01"))) +
  scale_fill_manual(values=c("lightblue","red")) +
  labs(x="Year of HIV Diagnosis",  y="Frequency")+
  guides(fill=guide_legend(title = "Subtype", reverse = F, size= 14))+
  theme_bw() + 
  theme(
    plot.title = element_text(color="black", size=14, face="bold.italic"),
    axis.title.x = element_text(color="black", size=14, face="bold"),
    axis.title.y = element_text(color="black", size=14, face="bold"),
    #legend.position = c(0.9, 0.8),
    legend.text = element_text(colour="black", size=14, face="plain"),
    legend.title = element_text(colour="black", size=14, face="bold")) +
  geom_hline(yintercept=0, linetype="dashed", color = "black")
P_N

## Plot the prevalence as Fraction:
P_F <-ggplot(MSM_PREV_FRACTION, aes(x=diag_year,y=perc_prev, fill= type)) +
  geom_bar(position="stack", stat="identity") +
  scale_x_date(limits = as.Date(c("1990-01-01","2020-01-01"))) +
  scale_y_continuous(breaks = c(seq(0,1,0.1)))+
  scale_fill_manual(values=c("lightblue","red")) +
  labs(x="Year of HIV Diagnosis",  y="Fraction")+
  guides(fill=guide_legend(title = "Subtype", reverse = F, size= 14))+
  theme_bw() + 
  theme(
    plot.title = element_text(color="black", size=14, face="bold.italic"),
    axis.title.x = element_text(color="black", size=14, face="bold"),
    axis.title.y = element_text(color="black", size=14, face="bold"),
    #legend.position = c(0.9, 0.8),
    legend.text = element_text(colour="black", size=14, face="plain"),
    legend.title = element_text(colour="black", size=14, face="bold")) +
  geom_hline(yintercept=0, linetype="dashed", color = "black")
P_F




################################
#### INCIDENCE
################################

#Get the dataframes:
MSM_INC_N <- rbind(
  data.frame("diag_year"= MSM_DD$diag_year, "incidence" =MSM_DD$incidence_all.x, "type"="B"),
  data.frame("diag_year"= MSM_DD$diag_year, "incidence" =MSM_DD$incidence_all.y, "type"="Non-B")
)


MSM_INC_FRACTION <- rbind(
  data.frame("diag_year"= MSM_DD$diag_year, "perc_inc" =MSM_DD$returns_perc_inc_B, "type"="B"),
  data.frame("diag_year"= MSM_DD$diag_year, "perc_inc" =MSM_DD$returns_perc_inc_NB, "type"="Non-B")
)


################################
## Stacked Barplots

## Plot the Incidence in absolute numbers:
I_N <- ggplot(MSM_INC_N, aes(x=diag_year,y=incidence, fill= type)) +
  geom_bar(position="stack", stat="identity") +
  scale_x_date(limits = as.Date(c("1990-01-01","2020-01-01"))) +
  scale_fill_manual(values=c("lightblue","red")) +
  labs(x="Year of HIV Diagnosis",  y="Frequency")+
  guides(fill=guide_legend(title = "Subtype", reverse = F, size= 14))+
  theme_bw() + 
  theme(
    plot.title = element_text(color="black", size=14, face="bold.italic"),
    axis.title.x = element_text(color="black", size=14, face="bold"),
    axis.title.y = element_text(color="black", size=14, face="bold"),
    #legend.position = c(0.9, 0.8),
    legend.text = element_text(colour="black", size=14, face="plain"),
    legend.title = element_text(colour="black", size=14, face="bold")) +
  geom_hline(yintercept=0, linetype="dashed", color = "black")
I_N

## Plot the Incidence as Fraction:
I_F <-ggplot(MSM_INC_FRACTION, aes(x=diag_year,y=perc_inc, fill= type)) +
  geom_bar(position="stack", stat="identity") +
  scale_x_date(limits = as.Date(c("1990-01-01","2020-01-01"))) +
  scale_y_continuous(breaks = c(seq(0,1,0.1)))+
  scale_fill_manual(values=c("lightblue","red")) +
  labs(x="Year of HIV Diagnosis",  y="Fraction")+
  guides(fill=guide_legend(title = "Subtype", reverse = F, size= 14))+
  theme_bw() + 
  theme(
    plot.title = element_text(color="black", size=14, face="bold.italic"),
    axis.title.x = element_text(color="black", size=14, face="bold"),
    axis.title.y = element_text(color="black", size=14, face="bold"),
    #legend.position = c(0.9, 0.8),
    legend.text = element_text(colour="black", size=14, face="plain"),
    legend.title = element_text(colour="black", size=14, face="bold")) +
  geom_hline(yintercept=0, linetype="dashed", color = "black")
I_F



# 
# # Move to a new page
# grid.newpage()
# # Create layout : nrow = 2, ncol = 2
# pushViewport(viewport(layout = grid.layout(nrow = 2, ncol = 2)))
# 
# # A helper function to define a region on the layout
# define_region <- function(row, col){
#   viewport(layout.pos.row = row, layout.pos.col = col)
# } 
# 
# # Arrange the plots
# print(I_N, vp = define_region(row = 1, col = 1))   # Span over two columns
# print(I_F, vp = define_region(row = 2, col = 1))
# print(P_N, vp = define_region(row = 1, col = 2))
# print(P_F, vp = define_region(row = 2, col = 2))
# 

ggarrange(I_N, P_N, I_F, P_F, ncol = 2, nrow = 2,
          labels = c("A", "B","C","D"), 
          common.legend = TRUE, legend = "bottom",
          font.label = list(size = 18, color = "black"))



ggarrange(I_F,I_N, ncol = 2, nrow = 1,
          labels = c("C", "D"), 
          common.legend = TRUE, legend = "bottom",
          font.label = list(size = 18, color = "black"))







############
#Fraction of MSM and HETs combined:
ggarrange(I_F,I_N, 
          HET_I_F,HET_I_N,
          ncol = 2, nrow = 2,
          labels = c("A", "B","C","D"), 
          common.legend = TRUE, legend = "bottom",
          font.label = list(size = 18, color = "black"))

