# June 2020, JDR

### Calculate and Transform Variables for Characteristics
# Check the subtype, population and demographics
# Datasets: PAT, ADMIN and table_subtypes

############################

### Clear R's memory
rm(list=ls())


### get the packages
if (!requireNamespace("pacman", quietly = TRUE))
  install.packages("pacman")
pacman::p_load("dplyr","tidyr", "data.table", "gmodels", "tidyverse", "boot", "table1", "devtools", "survival")

### Load the data
setwd("~/Desktop/SHCS")
data <- read.csv("Output/fulltable_study_population.csv")
data <- data[,-1] #11 168 


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




# 
# #### regrouping region
# data$copyregion <- data$region
# 
# data$copyregion[data$copyregion ==419] <- "Latin America and the Carribean"
# data$copyregion[data$copyregion ==5] <- "South America"
# data$copyregion[data$copyregion ==13] <- "Central America"
# data$copyregion[data$copyregion ==11] <- "Western Africa"
# data$copyregion[data$copyregion ==17] <- "Middle Africa"
# data$copyregion[data$copyregion ==14] <- "Eastern Africa"
# data$copyregion[data$copyregion ==18] <- "Southern Africa"
# data$copyregion[data$copyregion ==145] <- "Western Asia"
# data$copyregion[data$copyregion ==15] <- "Northern Africa"
# data$copyregion[data$copyregion ==30] <- "Eastern Asia"
# data$copyregion[data$copyregion ==34] <- "Southern Asia"
# data$copyregion[data$copyregion ==35] <- "South-Eastern Asia"
# data$copyregion[data$copyregion ==21] <- "Northern America"
# data$copyregion[data$copyregion ==9] <- "Oceania"
# data$copyregion[data$copyregion ==151] <- "Eastern Europe"
# data$copyregion[data$copyregion ==154] <- "Northern Europe"
# data$copyregion[data$copyregion ==155] <- "Western Europe"
# data$copyregion[data$copyregion ==39] <- "Southern Europe"
# data$copyregion[data$copyregion ==1] <- "Unknown"



# 
# ## make supper grouped regions: 
# data$groupedregion <- data$region
# data$groupedregion[data$groupedregion ==419] <- 5
# data$groupedregion[data$groupedregion ==5] <- 3
# data$groupedregion[data$groupedregion ==13] <- 3
# data$groupedregion[data$groupedregion ==11] <- 5
# data$groupedregion[data$groupedregion ==17] <- 5
# data$groupedregion[data$groupedregion ==14] <- 5
# data$groupedregion[data$groupedregion ==18] <- 5
# data$groupedregion[data$groupedregion==145] <- 4
# data$groupedregion[data$groupedregion==15] <- 5
# data$groupedregion[data$groupedregion ==30] <- 4
# data$groupedregion[data$groupedregion==34] <- 4
# data$groupedregion[data$groupedregion ==35] <- 4
# data$groupedregion[data$groupedregion ==21] <- 3
# data$groupedregion[data$groupedregion ==9] <- 9
# data$groupedregion[data$groupedregion ==151] <- 2
# data$groupedregion[data$groupedregion ==154] <- 2
# data$groupedregion[data$groupedregion ==155] <- 2
# data$groupedregion[data$groupedregion ==39] <- 2
# data$groupedregion[data$groupedregion==1] <- 1




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


FUN_center_cat <- function(x){
  if (is.na(x)){return (NA)}
  else if (x == 10){return("Zurich")}
  else if (x != 10){return("not-Zurich")}
}

NB$center_cat <- sapply(NB$center1, FUN_center_cat)



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
label(NB$cd4_diag) <- "CD4 cell count at diagnosis"
label(NB$known_swiss_infect) <- "Infection place";
label(NB$edu) <- "Education"


# ## better labels and units:
# label(NB$sex)       <- "Sex"; units(NB$sex)       <- "%"
# label(NB$age)       <- "Age"; units(NB$age)       <- "years"
# label(NB$age_sampledate) <- "Age at sample date"; units(NB$age_sampledate)     <- "years"
# label(NB$age_at_diag) <- "Age at diagnosis date"; units(NB$age_at_diag)     <- "years"
# label(NB$risk)     <- "Tranmission group"; units(NB$risk)       <- "%"
# label(NB$copyregion) <- "Geographic region of origin" ;units(NB$copyregion)       <- "%"
# label(NB$groupedregion) <- "Geographic region of origin";units(NB$groupedregion)       <- "%"
# label(NB$ethnicity) <- "Ethnicity";units(NB$ethnicity)       <- "%"





rndr <- function(x, name, ...) {
  if (!is.numeric(x)) return(render.categorical.default(x))
  what <- switch(name,
                 risk = "",
                 ethnicity = "",
                 copyregion = "",
                 regroupedregion = "",
                 age = "Median [Q1,Q3]",
                 age_sampledate = "Mean ± SD",
                 age_at_diag = "Median [Q1,Q3]",
                 center_cat = "",
                 diag_year = "Median [Q1,Q3]",
                 cd4_diag = "Median [Q1,Q3]",
                 known_swiss_infect = "",
                 edu = "")
  #age_sampledate = "Median [Min,Max]")
  parse.abbrev.render.code(c("", what))(x)
}

tab <- table1(~ age + age_at_diag + 
                center_cat + ethnicity + regroupedregion +
                known_swiss_infect +
                edu+
                cd4_diag 
              #diag_year
              | copy_subtype , 
              data = NB, overall="Total",
              render=rndr, 
              rowlabelhead = "Variables", 
              footnote = "<b>NOTE</b> All numbers are no. (%) unless otherwise stated, All modes of HIV transmission (n = 11168)<br/>
                Non-B: A, CRF01_AE, CRF02_AG, C, D, F, G, CRFs and unrecognized recombinants",
              caption = "<h3><b> Table 1. Patient characteristics and associations with subtype </b></h3>")
print(tab)



#######################################
#MSM Table
MSM <- data[data$risk == "1",]  #4751
library(tidyr)
MSM <- MSM %>% drop_na(id)


#######################################
##change the levels:
MSM$ethnicity <- factor(MSM$ethnicity, levels=c(1,2,3,4),
                       labels=c("White", 
                                "Black", #rest non-whites!
                                "Hispanic",
                                "Asian"))

MSM$copyregion <- factor(MSM$region, levels=c(39,
                                                11,
                                                14,
                                                5,
                                                30,
                                                15,
                                                151,
                                                9),
                        #1),
                        labels=c("Western, Northern and Southern Europe",
                                 "Western and Central Africa",
                                 "Eastern and Southern Africa",
                                 "Latin America",
                                 "Eastern and Southern Asia",
                                 "Western Asia and North Africa",
                                 "Eastern Europe",
                                 "North America and Oceania"))


# MSM$groupedregion <- factor(MSM$groupedregion, levels=c(2,
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


FUN_Switzerland <- function(x){
  if (is.na(x)){return (NA)}
  else if (x == 9){return(NA)}
  else if (x == 1){return("Swiss")}
  else if (x == 2){return("Not-Swiss")}
  else if (x == 3){return("Not-Swiss")}
}
MSM$known_swiss_infect <- sapply(MSM$infect_place, FUN_Switzerland)
#MSM <- within(MSM, known_swiss_infect <- relevel(known_swiss_infect , ref = "Swiss"))



FUN_edu <- function (x){
  if (is.na(x)){return (NA)}
  else if (x == 1) {return("no completed school")}
  else if (x == 2) {return("mandatory school")}
  else if (x == 3) {return("apprenticeship")}
  else if (x == 4 | x==5| x==6| x==7) {return("higher education")}
  else {return(NA)}
} ## values: not completed  school, mandatory school, apprentiship, higher education,  unkown

MSM$edu <- sapply(MSM$education, FUN_edu) 


FUN_center_cat <- function(x){
  if (is.na(x)){return (NA)}
  else if (x == 10){return("Zurich")}
  else if (x != 10){return("not-Zurich")}
}

MSM$center_cat <- sapply(MSM$center1, FUN_center_cat)



FUN_CD4_diag <- function(x){
  if (is.na(x)){return(NA)}
  else if (x <200){return("<200")} 
  else if (x >= 200 & x < 500){return("200-499")} 
  else if (x >= 500){return("≥500")}
  else return(NA)
}
MSM$cd4_diag_cat <- sapply(MSM$cd4_diag, FUN_CD4_diag)


FUN_age_diag_cat <- function(x){
  if (is.na(x)){return(NA)}
  else if (x <25){return("<25")} 
  else if (x >= 25 & x < 35 ){return("25-34")} 
  else if (x >= 35){return("≥35")}
  else return(NA)
}
MSM$age_diag_cat <- sapply(MSM$age_at_diag, FUN_CD4_diag)




## better labels and units:
label(MSM$age)       <- "Age"; units(MSM$age)       <- "years"
label(MSM$age_at_diag) <- "Age at diagnosis"; units(MSM$age_at_diag)     <- "years"
label(MSM$copyregion) <- "Geographic region of origin" 
label(MSM$groupedregion) <-"Geographic region of origin"
label(MSM$regroupedregion) <-"Geographic region of origin"
label(MSM$ethnicity) <- "Ethnicity"
label(MSM$center_cat) <- "Cohort Center "
MSM$diag_year <- as.integer(MSM$diag_year)
label(MSM$diag_year) <- "Year of diagnosis"
label(MSM$cd4_diag) <- "CD4 cell count at diagnosis"
label(MSM$known_swiss_infect) <- "Infection place";
label(MSM$edu) <- "Education"



# label(MSM$age)       <- "Age"; units(MSM$age)       <- "years"
# label(MSM$age_at_diag) <- "Age at diagnosis"; units(MSM$age_at_diag)     <- "years"
# label(MSM$copyregion) <- "Geographic region of origin" ;units(MSM$copyregion)       <- "%"
# label(MSM$groupedregion) <-"Geographic region of origin"; units(MSM$groupedregion)  <- "%"
# label(MSM$regroupedregion) <-"Geographic region of origin"; units(MSM$regroupedregion)  <- "%"
# label(MSM$ethnicity) <- "Ethnicity";units(MSM$ethnicity)       <- "%"
# label(MSM$center_cat) <- "Cohort Center ";units(MSM$center_cat)       <- "%"
# MSM$diag_year <- as.integer(MSM$diag_year)
# label(MSM$diag_year) <- "Year of diagnosis"
# label(MSM$cd4_diag) <- "CD4 cell count at diagnosis"
# label(MSM$known_swiss_infect) <- "Infection place";units(MSM$known_swiss_infect)       <- "%"
# label(MSM$edu) <- "Education";units(MSM$edu)       <- "%"


rndr <- function(x, name, ...) {
  if (!is.numeric(x)) return(render.categorical.default(x))
  what <- switch(name,
                 risk = "",
                 ethnicity = "",
                 copyregion = "",
                 regroupedregion = "",
                 age = "Median [Q1,Q3]",
                 age_sampledate = "Mean ± SD",
                 age_at_diag = "Median [Q1,Q3]",
                 center_cat = "",
                 diag_year = "Median [Q1,Q3]",
                 cd4_diag = "Median [Q1,Q3]",
                 known_swiss_infect = "",
                 edu = "")
  #age_sampledate = "Median [Min,Max]")
  parse.abbrev.render.code(c("", what))(x)
}

tab2 <- table1(~ age + age_at_diag + 
                 center_cat + ethnicity + regroupedregion +
                 known_swiss_infect +
                 edu+
                 cd4_diag 
                 # diag_year
               | copy_subtype , 
               data = MSM, overall="Total",
               render=rndr, 
               rowlabelhead = "Variables", 
               footnote = "<b>NOTE</b> All numbers are no. (%) unless otherwise stated, Homosexual contacts only (n = 4351)<br/>
               Non-B: A, CRF01_AE, CRF02_AG, C, D, F, G, CRFs and unrecognized recombinants",
               caption = "<h3><b> Table 1. Patient characteristics and associations with subtype</b></h3>")
print(tab2)




