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
data <- read.csv("Output/fulltable_study_population_update.csv")
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

NB$cd4_diag_per200

#######################################

table(data$risk)
# 1    2    3    5 
# 4351 3923 2392  102 

#HETTable
HET <- data[data$risk == "2",]  #3923
library(tidyr)
HET <- HET %>% drop_na(id)


#######################################
##change the levels:
HET$ethnicity <- factor(HET$ethnicity, levels=c(1,2,3,4),
                       labels=c("White", 
                                "Black", #rest non-whites!
                                "Hispanic",
                                "Asian"))

HET$copyregion <- factor(HET$region, levels=c(39,
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
HET$known_swiss_infect <- sapply(HET$infect_place, FUN_Switzerland)
#MSM <- within(MSM, known_swiss_infect <- relevel(known_swiss_infect , ref = "Swiss"))



FUN_edu <- function (x){
  if (is.na(x)){return (NA)}
  else if (x == 1) {return("no completed school")}
  else if (x == 2) {return("mandatory school")}
  else if (x == 3) {return("apprenticeship")}
  else if (x == 4 | x==5| x==6| x==7) {return("higher education")}
  else {return(NA)}
} ## values: not completed  school, mandatory school, apprentiship, higher education,  unkown

HET$edu <- sapply(HET$education, FUN_edu) 


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

HET$center_cat <- sapply(HET$center1, FUN_center_cat)

FUN_center_cat1 <- function(x){
  if (is.na(x)){return (NA)}
  else if (x == 10){return("Zurich")}
  else if (x != 10){return("not-Zurich")}
}
HET$center_cat1 <- sapply(HET$center1, FUN_center_cat1)



FUN_CD4_diag <- function(x){
  if (is.na(x)){return(NA)}
  else if (x <200){return("<200")} 
  else if (x >= 200 & x < 500){return("200-499")} 
  else if (x >= 500){return("≥500")}
  else return(NA)
}
HET$cd4_diag_cat <- sapply(HET$cd4_diag, FUN_CD4_diag)


FUN_age_diag_cat <- function(x){
  if (is.na(x)){return(NA)}
  else if (x <25){return("<25")} 
  else if (x >= 25 & x < 35 ){return("25-34")} 
  else if (x >= 35){return("≥35")}
  else return(NA)
}
HET$age_diag_cat <- sapply(HET$age_at_diag, FUN_age_diag_cat)




## better labels and units:
label(HET$age)       <- "Age"; units(HET$age)       <- "years"
label(HET$age_at_diag) <- "Age at diagnosis"; units(HET$age_at_diag)     <- "years"
label(HET$copyregion) <- "Geographic region of origin" 
label(HET$groupedregion) <-"Geographic region of origin"
label(HET$regroupedregion) <-"Geographic region of origin"
label(HET$ethnicity) <- "Ethnicity"
label(HET$center_cat) <- "Cohort Center "
HET$diag_year <- as.integer(HET$diag_year)
label(HET$diag_year) <- "Year of diagnosis"
label(HET$cd4_diag) <- "CD4 count at diagnosis";  units(HET$cd4_diag)       <- "cells per µl blood"
label(HET$cd4_diag_per200) <- "CD4 count at diagnosis";  units(HET$cd4_diag_per200)       <- "cells per 200 µl blood"
label(HET$log_mean_rna) <- "Viral Load";  units(HET$log_mean_rna)       <- "Log10 copies/ml"
label(HET$log_mean_rna_corrected) <- "Viral Load";  units(HET$log_mean_rna_corrected)       <- "Log10 copies/ml"

label(HET$known_swiss_infect) <- "Infection place";
label(HET$edu) <- "Education"



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
                 gender ="",
                 risk = "",
                 ethnicity = "",
                 copyregion = "",
                 regroupedregion = "",
                 age = "Median [Q1,Q3]",
                 age_sampledate = "Mean ± SD",
                 age_at_diag = "Median [Q1,Q3]",
                 center_cat = "",
                 diag_year = "Median [Q1,Q3]",
                 age_diag_cat ="",
                 cd4_diag = "Median [Q1,Q3]",
                 cd4_diag_per200= "Median [Q1,Q3]",
                 log_mean_rna =  "Median [Q1,Q3]",
                 log_mean_rna_corrected = "Median [Q1,Q3]",
                 known_swiss_infect = "",
                 edu = "")
  #age_sampledate = "Median [Min,Max]")
  parse.abbrev.render.code(c("", what))(x)
}

tab2 <- table1(~ 
                 # age + 
                 age_at_diag + 
                 center_cat + ethnicity + regroupedregion +
                 # known_swiss_infect +
                 # edu+
                 cd4_diag +
                 log_mean_rna
                 # diag_year
               | copy_subtype , 
               data = HET, overall="Total",
               render=rndr, 
               rowlabelhead = "Variables", 
               footnote = "<b>NOTE</b> All numbers are no. (%) unless otherwise stated, HETs only (n = 3923)<br/>
               Non-B: A, C, D, F, G, CRF01_AE, CRF02_AG, other CRFs and unrecognized recombinants",
               caption = "<h3><b> Table 1. Patient characteristics and associations with subtype</b></h3>")
print(tab2)


HET2 <- HET %>%
  mutate(subtypesNEW=case_when(
    HET$subtype=="B"  ~ "B",
    HET$subtype %in% c("A","C","F","G","01_AE","02_AG") ~ "Non-B")) #11 168
select <- which(is.na(HET2$subtypesNEW))
HET2$subtypesNEW[select] <- "Others"


tab2NEW <- table1(~ 
                    # age + 
                    age_at_diag + 
                    center_cat + ethnicity + regroupedregion +
                    # known_swiss_infect +
                    # edu+
                    cd4_diag_per200+
                    log_mean_rna
                  # diag_year
                  | subtypesNEW , 
                  data = HET2, overall="Total",
                  render=rndr, 
                  rowlabelhead = "Variables", 
                  footnote = "<b>NOTE</b> All numbers are no. (%) unless otherwise stated, HETs only (n = 3923)<br/>
               Non-B: A, C, D, F, G, CRF01_AE, CRF02_AG<br/>
               Others: Other subtypes, CRFs and unrecognized recombinants",
                  caption = "<h3><b> Table 1. HETs characteristics and associations with subtype</b></h3>")
print(tab2NEW)



HET_BnonB <- filter(HET2, subtype=="A" | subtype=="02_AG" 
                    | subtype=="C" | subtype=="01_AE"
                    | subtype=="F" | subtype=="G"| subtype =="B") 

tab3NEW <- table1(~ 
                    # age + 
                    age_at_diag + 
                    center_cat1 + ethnicity + regroupedregion +
                    # known_swiss_infect +
                    # edu+
                    cd4_diag+
                    #cd4_diag_per200 +
                    log_mean_rna
                  # diag_year
                  | subtypesNEW , 
                  data = HET_BnonB, overall="Total",
                  render=rndr, 
                  rowlabelhead = "Variables", 
                  footnote = "<b>NOTE</b> All numbers are no. (%) unless otherwise stated, HETs only (n = 4076)<br/>
               Non-B: A, C, D, F, G, CRF01_AE, CRF02_AG<br/>
               Others: Other subtypes, CRFs and unrecognized recombinants",
                  caption = "<h3><b> Table 1. HETs characteristics and associations with subtype</b></h3>")
print(tab3NEW)


HET_BnonB$center_cat1 <- sapply(HET_BnonB$center1, FUN_center_cat)
tab3NEW_cond <- table1(~ 
                         # age + 
                         age_at_diag + 
                         center_cat1 + ethnicity + regroupedregion +
                         # known_swiss_infect +
                         # edu+
                         cd4_diag_per200 +
                         log_mean_rna
                       # diag_year
                       | subtypesNEW , 
                       data = HET_BnonB, overall="Total",
                       render=rndr, 
                       rowlabelhead = "Variables", 
                       footnote = "<b>NOTE</b> All numbers are no. (%) unless otherwise stated, HETs only (n = 4076)<br/>
               Non-B: A, C, D, F, G, CRF01_AE, CRF02_AG<br/>
               Others: Other subtypes, CRFs and unrecognized recombinants",
                       caption = "<h3><b> Table 1. HETs characteristics and associations with subtype</b></h3>")
print(tab3NEW_cond)


table1(~subtype , data = HET_BnonB)

################################

#HET Table 2
HET_NonB <- filter(HET, subtype=="A" | subtype=="02_AG" 
                   | subtype=="C" | subtype=="01_AE"
                   | subtype=="F" | subtype=="G" )  ;unique(HET_NonB$subtype) #310
table1(~subtype , data = HET_NonB)



tab <- table1(~   # age + 
                age_at_diag + 
                center_cat1 + ethnicity + regroupedregion +
                # known_swiss_infect +
                # edu+
                cd4_diag +
                #cd4_diag_per200+
                log_mean_rna
              # log_mean_rna_corrected
              # diag_year
              | subtype , 
              data = HET_NonB, overall="Total",
              render=rndr, 
              rowlabelhead = "Variables", 
              footnote = "<b>NOTE</b> All numbers are no. (%) unless otherwise stated,  HETs only (n = 4076)",
              caption = "<h3><b> Table 2. HETs Characteristics Infected With the 6 Most Frequently Observed Non-B Subtypes </b></h3>")
print(tab)






