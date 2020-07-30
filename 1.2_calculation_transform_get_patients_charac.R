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
pacman::p_load("dplyr","haven","tidyr", "data.table", "gmodels", "tidyverse")

### Load the data
setwd("~/Desktop/SHCS")
DATA <- read.csv("Output/table_studypop.csv")
DATA <- DATA[,-1]

VAR_REGION <- read_dta("Input/var_region0520.dta")  #24
ADMI <- read_dta("Input/admi0520.dta") #20,904
TEST_INFO <- read_dta("Input/test_info0620.dta") #24 564
PAT <- read_dta("Input/pat0520.dta") #20 904
LAB <- read_dta("Input/lab0520.dta") #667 245

###################################################
#Getting the dataset to work with:
DD <- merge(DATA,ADMI[c("id","born", "exitdate", "region")], all.y=F, by = "id")
DD <- merge(DD,VAR_REGION , all.y=F, by = "region")
DD$sex <- NULL
DD <- merge(DD, PAT[c("id","sex","ethnicity","risk","regdate", 
                      "hiv_posdate", "hiv_posdocdate", "hiv_negdate",
                      "infect_place", "education", "center1")],all.y=F, by = "id")
dim(DD) #11 248 x 52 



######################

## NA's?
sum(is.na(DD$id))#0
sum(is.na(DD$subtype)) #0
sum(is.na(DD$seq_dt)) #0
sum(is.na(DD$sampledate)) #0
sum(is.na(DD$pol_len)) #0
sum(is.na(DD$rega_subtype)) #0
sum(is.na(DD$comet_subtype)) #0
sum(is.na(DD$lab)) #0
sum(is.na(DD$age)) #0

sum(is.na(DD$exitdate)) #9064
sum(is.na(DD$born)) #0


sum(is.na(DD$sex))#0
sum(is.na(DD$ethnicity))#4
sum(is.na(DD$risk)) #0
sum(is.na(DD$region)) #0
sum(is.na(DD$regdate)) #0
sum(is.na(DD$hiv_posdate)) #1351
sum(is.na(DD$hiv_posdocdate)) #3148
sum(is.na(DD$hiv_negdate)) #6414
sum(is.na(DD$cat)) #11,142
sum(is.na(DD$rule)) #0


sum(is.na(DD$education)) #430
sum(is.na(DD$infect_place)) #8369
sum(is.na(DD$center1)) #0


################################

### Output:
# ged rid of Subtype "unknown"!!
DD <- DD[DD$subtype != "UNKNOWN" & DD$subtype != "UNDETERMINED (SHORT SEQUENCE)", ] #11,171
DD <- DD[DD$subtype != "HIV 1 O GROUP" & DD$subtype != "HIV 1 N GROUP", ] #11, 168

### Output:
DD<- DD %>%
  mutate(copy_subtype=case_when(
    DD$subtype=="B"  ~ "B",
    #~ "non-B"
  )) #11 168
select <- which(is.na(DD$copy_subtype))
DD$copy_subtype[select] <- "Non-B"




#   Gender: document PAT, values: gender (sex)
FUN_gender <- function (x){
  if (is.na(x)){return (NA)}
  else if (x == 1) {return("male")}
  else if (x == 2) {return("female")}
  else {return(NA)}
}
DD$gender <- sapply(DD$sex, FUN_gender) 


#  Geographical origin: 
# 	Ethnicity: document PAT, values: (0 other), 1 white,2 black, 3 hispano-american, 4 Asian, 9 unknown
DD$var_eth <- factor(DD$ethnicity, levels=c(1,2,3,4,0,9),
                             labels=c("White", 
                                      "Black",
                                      "Hispano-American",
                                      "Asian",
                                      "Others",
                                      "Unknown"))
FUN_eth <- function (x){
  if (is.na(x)){return (NA)}
  else if (x == 0){return (NA)}
  else if (x == 9){return (NA)}
  else {return (as.character(x))}
}
DD$eth <- sapply(DD$ethnicity, FUN_eth)

FUN_eth_cat <- function (x){
  if (is.na(x)){return (NA)}
  else if (x == 9){return (NA)}
  else if (x == 1){return ("White")}
  else if (x != 1){return ("Not-white")}
  else {return (as.character(x))}
}
DD$ethnicity_cat <- sapply(DD$ethnicity, FUN_eth_cat)


# region 
DD$regroupedregion  <- factor(DD$region, levels=c(155,154,39,
                                              11,17,
                                              14,18,
                                              5,13,419,
                                              30,34,35,
                                              15,145,
                                              151,
                                              21,9,1),
                          labels=c("Western, Northern and Southern Europe",
                                   "Western, Northern and Southern Europe",
                                   "Western, Northern and Southern Europe",
                                   "Western and Central Africa",
                                   "Western and Central Africa",
                                   "Eastern and Southern Africa",
                                   "Eastern and Southern Africa",
                                   "Latin America",
                                   "Latin America",
                                   "Latin America",
                                   "Eastern and Southern Asia",
                                   "Eastern and Southern Asia",
                                   "Eastern and Southern Asia",
                                   "Western Asia and North Africa",
                                   "Western Asia and North Africa",
                                   "Eastern Europe",
                                   "North America and Oceania",
                                   "North America and Oceania",
                                   "Unknown"))







############
# 	source of infection (risk): document PAT, 
#values: homosexual contact, heterosexual contact, IDU, clotting factors against haemophilia, blood transfusion, other blood products, perinatal transmission (MTCT), other sources or unknown
DD$risk_all <- DD$risk
FUN_risk <- function (x){
  if (is.na(x)){return (NA)}
  else if (x == 0){return (NA)}
  else if (x == 4){return (3)}
  else if (x == 6){return (5)}
  else if (x == 7){return (NA)} #MTCT
  else if (x == 9){return (NA)}
  else {return (as.character(x))}
}
DD$risk<- sapply(DD$risk, FUN_risk)


FUN_risk_cat <- function (x){
  if (is.na(x)){return (NA)}
  else if (x == 1){return ("MSM")}
  else if (x == 2){return ("HET")}
  else if (x == 3){return ("IDU")}
  else if (x == 4){return ("IDU")}
  else if (x == 5){return ("Blood Transfusion/ Products")}
  else if (x == 6){return ("Blood Transfusion/ Products")}
  else if (x == 7){return ("Other/ Unknown")}
  else if (x == 9){return ("Other/ Unknown")}
  else if (x == 0){return ("Other/ Unknown")}
  else {return (as.character(x))}
}
DD$risk_cat <- sapply(DD$risk_all, FUN_risk_cat)


#### ADMI Table: Exitdate, born
##Calculate the age
DD$exitdate <- as.Date(as.character(DD$exitdate), tryFormats = c("%Y-%m-%d"))
DD <- DD %>%
  mutate(age=case_when(
    is.na(DD$exitdate)==TRUE  ~ year(Sys.Date())-DD$born,
    is.na(DD$exitdate)==FALSE ~ year(DD$exitdate)-DD$born))


## Calculate the age at sample date 
DD$sampledate<- as.Date(as.character(DD$sampledate), tryFormats = c("%Y-%m-%d"))
DD <- DD %>%
  mutate(age_sampledate=case_when(
    is.na(DD$sampledate)==FALSE  ~ year(DD$sampledate)-DD$born)) 



## get the diagnosis date:
FUN_diag_date <- function(ID){
  id_test_set <- DD[DD$id == ID,]
  dates <- data.frame(matrix(0, ncol = 1))
  dates[1,1] <- as.numeric(as.Date(id_test_set$hiv_posdate))
  dates[2,1] <- as.numeric(as.Date(id_test_set$hiv_posdocdate))
  dates[3,1] <- as.numeric(as.Date(id_test_set$regdate))
  diag_date <- min(dates,na.rm = TRUE)
  diag_date<- format(as.Date(diag_date,origin="1970-01-01"))
  return(paste(diag_date))}

DD$diag_date <- sapply(DD$id, FUN_diag_date)
DD$diag_year  <- substr(DD$diag_date,1,4)


#Age at Diagnosis
DD <- DD %>%
  mutate(age_at_diag=case_when(
    is.na(DD$diag_date)==FALSE  ~ year(DD$diag_date) - DD$born)) 




# 	Infected in Switzerland: document PAT, value: if infect place is Switzerland (others unknown or no)
## 1=in Switzerland, 2=while temporarily abroad, 3=as a resident abroad, 9=unknown
FUN_Switzerland <- function(x){
  if (is.na(x)){return (NA)}
  else if (x == 9){return(NA)}
  else if (x == 1){return("Swiss")}
  else if (x == 2){return("Not-Swiss")}
  else if (x == 3){return("Not-Swiss")}
}
DD$known_swiss_infect <- sapply(DD$infect_place, FUN_Switzerland)




#   education

# 1=no completed school N(1) yes
# (check only one). or professional education
# 2=mandatory school(9 years
#                    in Switzerland)
# 3=finished apprenticeship
# 4=bachelor
# 5=higher professional education
# 6=higher technical or commercial school
# 7=university
# 0=other
# 9=no information
FUN_edu <- function (x){
  if (is.na(x)){return (NA)}
  else if (x == 1) {return(1)}
  else if (x == 2) {return(2)}
  else if (x == 3) {return(3)}
  else if (x == 4 | x==5| x==6| x==7) {return(4)}
  else {return(NA)}
} ## values: not complied school, mandatory school, apprenticeship, higher education,  unknown

DD$edu <- sapply(DD$education, FUN_edu) 




# 	Cohort-center: document PAT, 
#   	values: recruitment center (10 Zurich,20 Basel, 30 Bern, 40 Geneva, 50 Lausanne, 60 Lugano, 70 St.Gallen) 
FUN_center_cat <- function(x){
  if (is.na(x)){return (NA)}
  else if (x == 10){return("Zurich")}
  else if (x != 10){return("not-Zurich")}
}

DD$center_cat <- sapply(DD$center1, FUN_center_cat)



# 	CD4: Document LAB, value: CD4 cells per ul blood
LAB$cd4date <- as.Date(LAB$cd4date, tryFormats  = c("%Y-%m-%d"))
lab <- LAB[!is.na(LAB$cd4date), ] #629 589
lab <- lab[order(lab$cd4date), ]
lab <- lab[names(lab) %in% c('id', 'labdate', 'cd4date', 'cd4', 'rna')]

#   	Mean: of all CD4 values of a patient
agg_ratioCD4<- aggregate(cd4 ~ id, data = lab, FUN = mean) #mean of all tests
names(agg_ratioCD4) <- c("id", "mean_cd4")
DD <- merge(DD, agg_ratioCD4, by="id", all.x = T) #mean in cd4 of all tests, 11 163

#   	Number below 500: of all CD4 values how many times value below 500 
CD4below500<- aggregate(cd4 ~ id, data = lab, FUN = function(x){x<500}) #mean of all tests
names(CD4below500) <- c("id", "num_cd4_below500")

#   	Nadir: lowest CD4 value ever measured
agg_minCD4<- aggregate(cd4 ~ id, data = lab, FUN = min)
names(agg_minCD4) <- c("id", "CD4_nadir")
DD<- merge(DD, agg_minCD4, by="id", all.x = T) 
# need to divide by 200 for each point:
agg_minCD4per200<-aggregate(CD4_nadir ~ id, data= DD, FUN = function(x){x/200})
names(agg_minCD4per200) <- c("id", "CD4_nadir_per200")
DD <- merge(DD, agg_minCD4per200, by="id", all.x = T) 



# cd4 diagnosis / earliest cd4
lab1 <- lab[order(lab$cd4date),]
lab2 <- lab1[!is.na(lab1$cd4),] #need to have a cd4 value!
lab_cd4_diag <- lab2[!duplicated(lab2$id, formatLast=TRUE),]
lab_cd4_diag <- lab_cd4_diag[names(lab_cd4_diag) %in% c('id', "cd4", 'labdate')]; names(lab_cd4_diag) <- c("id","cd4_diag", "labdate")
DD <- merge(DD, lab_cd4_diag, by="id", all.x = T)
# need to divide by 200 for each point:
agg_cd4_diagper200<-aggregate(cd4_diag ~ id, data= DD, FUN = function(x){x/200})
names(agg_cd4_diagper200) <- c("id", "cd4_diag_per200")
DD <- merge(DD, agg_cd4_diagper200, by="id", all.x = T)


FUN_CD4_diag <- function(x){
  if (is.na(x)){return(NA)}
  else if (x <200){return("<200")} 
  else if (x >= 200 & x <= 500){return("200-500")} 
  else if (x >= 500){return(">500")}
  else return(NA)
}
DD$cd4_diag_cat <- sapply(DD$cd4_diag, FUN_CD4_diag)




## 	viral load: Document LAB, value: number of copies of HIV-1 RNA in a mL of blood.
#need to transform variable: log10 of viral RNA
#Logarithms allow exponential behavior to be visualized in a linear manner!
lab$labdate <-  as.Date(as.character(lab$labdate), tryFormats  = c("%Y-%m-%d"))
lab_noNA <- lab[!is.na(lab$rna),] #need to have a cd4 value!
lab$logRNA <- log10(lab$rna)

#   	Mean: of all viral loads of a patient (includes 0)
agg_ratioRNA_log <- aggregate(logRNA~ id, data = lab, FUN = mean) #mean of all tests
names(agg_ratioRNA_log) <- c("id", "log_mean_rna")
agg_ratioRNA <- aggregate(rna~ id, data = lab, FUN = mean) #mean of all tests
names(agg_ratioRNA) <- c("id", "mean_rna")
DD <- merge(DD, agg_ratioRNA, by="id", all.x = T)
DD <- merge(DD, agg_ratioRNA_log, by="id", all.x = T)


dim(DD) # 11 168 76
############################
### Save the dataset:
#write.csv(DD,"Output/fulltable_study_population.csv")


