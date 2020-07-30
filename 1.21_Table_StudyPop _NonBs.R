# June 2020, JDR

# Table1
# Of Non-Bs subtypes: A, AG, C,AE, F and G

############################

### Clear R's memory
rm(list=ls())


### get the packages
if (!requireNamespace("pacman", quietly = TRUE))
  install.packages("pacman")
pacman::p_load("dplyr","tidyr", "data.table", "gmodels", 
               "tidyverse", "boot", "table1", "devtools", "survival",
               "plyr", "lubricate")


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


### need to use another variable, which is the same for all the entries! 
data$all <- "1"
data$all <- as.numeric(data$all)


table1(~subtype , data = data )
### A, AG, C,AE, G and F most prominant after B!



#### regrouping region


data$copyregion <- data$region
data$copyregion[data$copyregion ==419] <- 5
data$copyregion[data$copyregion ==5] <- 5
data$copyregion[data$copyregion ==13] <- 5
data$copyregion[data$copyregion ==11] <- 11
data$copyregion[data$copyregion ==17] <- 11
data$copyregion[data$copyregion ==14] <- 14
data$copyregion[data$copyregion ==18] <- 14
data$copyregion[data$copyregion ==145] <- 15
data$copyregion[data$copyregion ==15] <- 15
data$copyregion[data$copyregion ==30] <- 30
data$copyregion[data$copyregion ==34] <- 30
data$copyregion[data$copyregion ==35] <- 30
data$copyregion[data$copyregion ==21] <- 9
data$copyregion[data$copyregion ==9] <- 9
data$copyregion[data$copyregion ==151] <- 151
data$copyregion[data$copyregion ==154] <- 39
data$copyregion[data$copyregion ==155] <- 39
data$copyregion[data$copyregion ==39] <- 39
data$copyregion[data$copyregion ==1] <- 1


## make supper grouped regions: 
data$groupedregion <- data$region
data$groupedregion[data$groupedregion ==419] <- 5
data$groupedregion[data$groupedregion ==5] <- 3
data$groupedregion[data$groupedregion ==13] <- 3
data$groupedregion[data$groupedregion ==11] <- 5
data$groupedregion[data$groupedregion ==17] <- 5
data$groupedregion[data$groupedregion ==14] <- 5
data$groupedregion[data$groupedregion ==18] <- 5
data$groupedregion[data$groupedregion==145] <- 4
data$groupedregion[data$groupedregion==15] <- 5
data$groupedregion[data$groupedregion ==30] <- 4
data$groupedregion[data$groupedregion==34] <- 4
data$groupedregion[data$groupedregion ==35] <- 4
data$groupedregion[data$groupedregion ==21] <- 3
data$groupedregion[data$groupedregion ==9] <- 9
data$groupedregion[data$groupedregion ==151] <- 2
data$groupedregion[data$groupedregion ==154] <- 2
data$groupedregion[data$groupedregion ==155] <- 2
data$groupedregion[data$groupedregion ==39] <- 2
data$groupedregion[data$groupedregion==1] <- 1



######################################################
##get dataset only with subtypes A, AE, C, AG, F and G!
NB <- filter(data, subtype=="A" | subtype=="02_AG" 
             | subtype=="C" | subtype=="01_AE"
             | subtype=="F" | subtype=="G" )  ;unique(NB$subtype) #310

table1(~subtype , data = NB)





#######################################
##change the levels:
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
               | subtype , 
               data = NB, overall="Total",
               render=rndr, 
               rowlabelhead = "Variables", 
               footnote = "<b>NOTE</b> All numbers are no. (%) unless otherwise stated, All modes of HIV transmission (n = 11168)",
               caption = "<h3><b> Table 2. Characteristics of SHCS Study Participants Infected With the 6 Most Frequently Observed Non-B Subtypes </b></h3>")
print(tab)







################################################
#MSM Table
MSM <- data[data$risk == "1",]  #4751

table1(~subtype , data = MSM)

library(tidyr)
MSM <- MSM %>% drop_na(id)


######################################################
##get dataset only with subtypes A, AE, C, AG, F and G!
NB <- filter(MSM, subtype=="A" | subtype=="02_AG" 
             | subtype=="C" | subtype=="01_AE"
             | subtype=="F" | subtype=="G" )  ;unique(NB$subtype) #310

table1(~subtype , data = NB)





#######################################
##change the levels:
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



tab2 <- table1(~ age + age_at_diag + 
                center_cat + ethnicity + regroupedregion +
                known_swiss_infect +
                edu+
                cd4_diag 
                # diag_year
              | subtype , 
              data = NB, overall="Total",
              render=rndr, 
              rowlabelhead = "Variables", 
              footnote = "<b>NOTE</b> All numbers are no. (%) unless otherwise stated, Homosexual contacts only (n = 4351)",
              caption = "<h3><b> Table 2. Characteristics of SHCS Study Participants Infected With the 6 Most Frequently Observed Non-B Subtypes </b></h3>")
print(tab2)

# 
# table1(~ sex +age+ age_sampledate + age_at_diag +  risk +ethnicity| subtype , data = NB, overall="Total",
#        render.continuous="Mean ± SD")
# 
# table1(~ sex +age+ age_sampledate + age_at_diag +  risk +ethnicity + groupedregion| subtype , data = NB, overall="Total",
#        render.continuous="Mean ± SD")
# 
# table1(~ sex +age+ age_sampledate + age_at_diag +  risk +ethnicity + copyregion| subtype , data = NB, overall="Total",
#        render.continuous="Mean ± SD")
# 
# 
# rndr <- function(x, name, ...) {
#   if (!is.numeric(x)) return(render.categorical.default(x))
#   what <- switch(name,
#                  risk = "",
#                  ethnicity = "",
#                  copyregion = "",
#                  groupedregion = "",
#                  age = "Mean ± SD",
#                  age_sampledate = "Mean ± SD",
#                  age_at_diag = "Mean ± SD")
#   #age_sampledate = "Median [Min,Max]")
#   parse.abbrev.render.code(c("", what))(x)
# }
# 
# 
# tab1 <- table1(~ sex +age + age_sampledate + age_at_diag + risk +ethnicity + groupedregion| subtype , 
#                data = NB, overall="Total",
#                render=rndr, 
#                rowlabelhead = "Variables", 
#                footnote = "Characteristics of the Non-B participants") ; print(tab1)
# 
# tab2 <- table1(~ sex +age + age_sampledate + age_at_diag + risk +ethnicity + regroupedregion| subtype , 
#                data = NB, overall="Total",
#                render=rndr, 
#                rowlabelhead = "Variables", 
#                footnote = "Characteristics of the Non-B participants"); print(tab2)
# 
# 


table.with.row.percentages <- function(formula, data) {
  formula[[2]][[1]] <- quote(`+`)
  m1 <- model.frame(formula, data=data, na.action=na.pass)
  lab <- table1::label(m1[[1]]) #"label of the variable"
  if (is.null(lab)) lab <- names(m1)[1]
  
  t1 <- table(m1[[1]], m1[[2]])
  t2 <- margin.table(t1, 2) #contingency table in array form, compute the sum of table entries for a given index. 2= col, 1= row
  t3 <- prop.table(t1, 1) #Express Table Entries as Fraction of Marginal Table, row percentage
  
  
  tab <- paste0(t1, " (", formatC(100*t3, format="f", digits=1), "%)")
  attributes(tab) <- attributes(t1)
  dimnames(tab)[[2]] <- paste0(dimnames(tab)[[2]], "<br/>(n=", t2, ")")
  
  tab <- rbind(rep("", ncol(tab)), tab)
  dimnames(tab)[[1]][1] <- lab
  
  x <- paste0(
    '<table>\n<thead>\n',
    table.rows(colnames(tab), th=T),
    '</thead>\n<tbody>\n',
    table.rows(tab),
    '</tbody>\n</table>\n')
  
  structure(x, class=c("table1", "html", "character"), html=TRUE)
}


# Example:
table.with.row.percentages(~ copyregion | subtype, data=NB)



