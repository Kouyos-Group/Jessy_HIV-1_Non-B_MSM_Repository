
###Jan2021
### Non-B subtype Distrubtion of HETs and MSM


###############


### Clear R's memory
rm(list=ls())


### get the packages
if (!requireNamespace("pacman", quietly = TRUE))
  install.packages("pacman")
pacman::p_load("dplyr","tidyr", "data.table", "gmodels", "tidyverse", "boot", "table1", "devtools", "survival")

### Load the data
setwd("~/Desktop/SHCS")
data <- read.csv("Output/fulltable_study_population_update.csv")
data <- data[,-1] #12 852 


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
# data$diag_year<- as.Date(as.character(data$diag_year), tryFormats = c("%Y"))


#######################################
##change the levels:
NB <- data
NB$sex <- factor(NB$sex, levels=c(1,2),
                 labels=c("Male", 
                          "Female"))

NB$risk <- factor(NB$risk, levels=c(1,2,3,4,5,6), #7,9,0 
                  labels=c("MSM", 
                           "HET",
                           "IDU",
                           "IDU",
                           "Blood Transfusion/ Products",
                           "Blood Transfusion/ Products"))
# "Other/ Unknown", #"MTCT",
# "Other/ Unknown",
# "Other/ Unknown"))


NB$ethnicity <- factor(NB$ethnicity, levels=c(1,2,3,4),
                       labels=c("White", 
                                "Black", #rest non-whites!
                                "Hispanic",
                                "Asian"))


# NB$copyregion <- factor(NB$copyregion, levels=c(39,
#                                                 11,
#                                                 14,
#                                                 5,
#                                                 30,
#                                                 15,
#                                                 151,
#                                                 9),
#                         #1),
#                         labels=c("Western, Northern and Southern Europe",
#                                  "Western and Central Africa",
#                                  "Eastern and Southern Africa",
#                                  "Latin America",
#                                  "Eastern and Southern Asia",
#                                  "Western Asia and North Africa",
#                                  "Eastern Europe",
#                                  "North America and Oceania"))
#"Unknown"))



# NB$groupedregion <- factor(NB$groupedregion, levels=c(2,
#                                                       3,
#                                                       4,
#                                                       5,
#                                                       9),
#                            #1),
#                            labels=c("Europe",
#                                     "America",
#                                     "Asia",
#                                     "Africa",
#                                     "Oceania"))
#"Unknown"))


FUN_Switzerland <- function(x){
  if (is.na(x)){return (NA)}
  else if (x == 9){return(NA)}
  else if (x == 1){return("Swiss")}
  else if (x == 2){return("Not-Swiss")}
  else if (x == 3){return("Not-Swiss")}
}
NB$known_swiss_infect <- sapply(NB$infect_place, FUN_Switzerland)
#MSM <- within(MSM, known_swiss_infect <- relevel(known_swiss_infect , ref = "Swiss"))



FUN_edu <- function (x){
  if (is.na(x)){return (NA)}
  else if (x == 1) {return("no completed school")}
  else if (x == 2) {return("mandatory school")}
  else if (x == 3) {return("apprenticeship")}
  else if (x == 4 | x==5| x==6| x==7) {return("higher education")}
  else {return(NA)}
} ## values: not completed  school, mandatory school, apprentiship, higher education,  unkown

NB$edu <- sapply(NB$education, FUN_edu) 


# 	Cohort-center: document PAT, 
#   	values: recruitment center (10 Zurich,20 Basel, 30 Bern, 40 Geneva, 50 Lausanne, 60 Lugano, 70 St.Gallen) 
FUN_center_cat <- function(x){
  if (is.na(x)){return (NA)}
  else if (x == 10){return("Zurich")}
  else if (x == 20){return("Basel")}
  else if (x == 30){return("Bern")}
  else if (x == 40){return("Geneva")}
  else if (x == 50){return("Lausanne")}
  else if (x == 60){return("Lugano")}
  else if (x == 70){return("St.Gallen")}
}

NB$center_cat <- sapply(NB$center1, FUN_center_cat)


FUN_center_cat1 <- function(x){
  if (is.na(x)){return (NA)}
  else if (x == 10){return("Zurich")}
  else if (x != 10){return("not-Zurich")}
}

NB$center_cat1 <- sapply(NB$center1, FUN_center_cat1)



FUN_CD4_diag <- function(x){
  if (is.na(x)){return(NA)}
  else if (x <200){return("<200")} 
  else if (x >= 200 & x < 500){return("200-499")} 
  else if (x >= 500){return("≥500")}
  else return(NA)
}
NB$cd4_diag_cat <- sapply(NB$cd4_diag, FUN_CD4_diag)


FUN_age_diag_cat <- function(x){
  if (is.na(x)){return(NA)}
  else if (x <25){return("<25")} 
  else if (x >= 25 & x < 35 ){return("25-34")} 
  else if (x >= 35){return("≥35")}
  else return(NA)
}
NB$age_diag_cat <- sapply(NB$age_at_diag, FUN_CD4_diag)


## better labels and units:
label(NB$age)       <- "Age"; units(NB$age)       <- "years"
label(NB$age_at_diag) <- "Age at diagnosis"; units(NB$age_at_diag)     <- "years"
label(NB$copyregion) <- "Geographic region of origin" 
label(NB$groupedregion) <-"Geographic region of origin"
label(NB$regroupedregion) <-"Geographic region of origin"
label(NB$ethnicity) <- "Ethnicity"
label(NB$center_cat) <- "Cohort Center "
NB$diag_year <- as.integer(NB$diag_year)
label(NB$diag_year) <- "Year of diagnosis"
label(NB$cd4_diag) <- "CD4 cell counts at diagnosis"
label(NB$known_swiss_infect) <- "Infection place";
label(NB$edu) <- "Education"


NonB <- filter(NB, subtype!="B")
NonB_MSMHET <- filter(NonB, risk_all=="1"| risk_all=="2")
table1(~subtype , data = NonB_MSMHET)


#by Risk group
HET <- data[data$risk == "2",]  
MSM <- data[data$risk == "1",]  
library(tidyr)
HET <- HET %>% drop_na(id)
MSM <- MSM %>% drop_na(id)

HET_NonB <- filter(HET, subtype!="B")
table1(~subtype , data = HET_NonB)

MSM_NonB <- filter(MSM, subtype!="B")
table1(~subtype , data = MSM_NonB)




