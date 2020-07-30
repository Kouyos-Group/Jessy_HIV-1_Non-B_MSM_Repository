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




#### regrouping region
data$copyregion <- data$region

data$copyregion[data$copyregion ==419] <- "Latin America and the Carribean"
data$copyregion[data$copyregion ==5] <- "South America"
data$copyregion[data$copyregion ==13] <- "Central America"
data$copyregion[data$copyregion ==11] <- "Western Africa"
data$copyregion[data$copyregion ==17] <- "Middle Africa"
data$copyregion[data$copyregion ==14] <- "Eastern Africa"
data$copyregion[data$copyregion ==18] <- "Southern Africa"
data$copyregion[data$copyregion ==145] <- "Western Asia"
data$copyregion[data$copyregion ==15] <- "Northern Africa"
data$copyregion[data$copyregion ==30] <- "Eastern Asia"
data$copyregion[data$copyregion ==34] <- "Southern Asia"
data$copyregion[data$copyregion ==35] <- "South-Eastern Asia"
data$copyregion[data$copyregion ==21] <- "Northern America"
data$copyregion[data$copyregion ==9] <- "Oceania"
data$copyregion[data$copyregion ==151] <- "Eastern Europe"
data$copyregion[data$copyregion ==154] <- "Northern Europe"
data$copyregion[data$copyregion ==155] <- "Western Europe"
data$copyregion[data$copyregion ==39] <- "Southern Europe"
data$copyregion[data$copyregion ==1] <- "Unknown"




##### Table1
dd <- data


table1(~factor(sex)+ age + factor(risk)+ factor(region)+ factor(eth)| copy_subtype, data = dd )


##change the levels:
dd$sex <- factor(dd$sex, levels=c(1,2),
                 labels=c("Male", 
                          "Female"))

dd$risk <- factor(dd$risk, levels=c(1,2,3,4,5,6),
                  labels=c("MSM", 
                           "HET",
                           "IDU",
                           "IDU",
                           "Blood Transfusion/ Products",
                           "Blood Transfusion/ Products"))
                           # "Other/ Unknown", #MTCT
                           # "Other/ Unknown",
                           # "Other/ Unknown"))


dd$ethnicity <- factor(dd$ethnicity, levels=c(1,2,3,4),
                       labels=c("White", 
                                "Black",
                                "Hispano-American",
                                "Asian"))

dd$region <- factor(dd$region, levels=c(155,154,39,
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


## better labels and units:
label(dd$sex)       <- "Sex"
label(dd$age)       <- "Age"
label(dd$risk)     <- "Tranmission group"
label(dd$region) <- "Geographic region of origin"
label(dd$ethnicity) <- "Ethnicity"
label(dd$age_sampledate) <- "Age at sample date"
label(dd$age_at_diag) <- "Age at diagnosis date"

units(dd$sex)       <- "%"
units(dd$age)       <- "years"
units(dd$age_sampledate)     <- "years"
units(dd$age_at_diag)     <- "years"
units(dd$risk)       <- "%"
units(dd$ethnicity)       <- "%"
units(dd$region)       <- "%"

table1(~ sex +age + risk +ethnicity| copy_subtype , data = dd, overall="Total" )



table1(~ sex +age  + risk +ethnicity| copy_subtype , data = dd, overall="Total",
       render.continuous="Mean (SD)")

table1(~copy_subtype| risk+sex, data=dd)

tabled <-table1(~ sex +age+ age_sampledate + age_at_diag +  risk +ethnicity| copy_subtype , data = dd, overall="Total",
                render.continuous="Mean (SD)")


dd[dd== -Inf] <- 1
table1(~ sex +age + risk +ethnicity + log_mean_rna| copy_subtype , data = dd, overall="Total" )


rndr <- function(x, name, ...) {
  if (!is.numeric(x)) return(render.categorical.default(x))
  what <- switch(name,
                 ethnicity = "",
                 risk = "",
                 age = "Mean (Q1,Q3)",
                 age_sampledate = "Mean (SD)",
                 age_at_diag = "Mean ± SD",
                 log_mean_rna = "Mean ± S ")
  #age_sampledate = "Median [Min,Max]")
  parse.abbrev.render.code(c("", what))(x)
}



tab1 <- table1(~ sex +age + age_sampledate + age_at_diag + risk +ethnicity + region| copy_subtype , 
               data = dd, overall="Total",
               render=rndr, rowlabelhead = "Variables")

print(tab1)





#################################
##### Table one
## Vector of variables to summarize
library("tableone")
col_order <- c("id",  "sampledate","seq_dt","rega_subtype", "comet_subtype","subtype",
               "lab","pol_len","pr_len", "rt_len", "cat","rule","nseq","total_nseq",
               "sex","born", "age", 
               "exitdate", "regdate", "hiv_posdate",  "hiv_posdocdate" ,"hiv_negdate",
               "region", "var_desc", "ethnicity", "risk")

myVars <-  c("sex", "age","risk","region", "ethnicity")

## Vector of categorical variables that need transformation
catVars <- c("sex","risk","region", "ethnicity")

tab2 <- CreateTableOne(vars = myVars, data = data, strata = c("copy_subtype"), factorVars = catVars)
tab2
#Binary categorical variables are summarized as counts and percentages of the second level. 
#For example, if it is codsed as 0 and 1, the “1” level is summarized. For 3+ category variable all levels are summarized. 
#Please bear in mind, the percentages are calculated after excluding missing values.


##Showing all levels for categorical variables
print(tab2, showAllLevels = TRUE)


##Detailed information including missingness
summary(tab2)

####################################


#MSM Table
dd3 <- dd[dd$risk == "MSM",]  #4751

# better labels and units:
label(dd3$sex)       <- "Sex"
label(dd3$age)       <- "Age"
label(dd3$region) <- "Geographic region of origin"
label(dd3$ethnicity) <- "Ethnicity"
label(dd3$age_sampledate) <- "Age at sample date"
label(dd3$age_at_diag) <- "Age at diagnosis"
label(dd3$cd4_diag_cat) <- "Cd4 at diagnosis"
units(dd3$sex)       <- "%"
units(dd3$age)       <- "years"
units(dd3$age_sampledate)     <- "years"
units(dd3$age_at_diag)     <- "years"
units(dd3$ethnicity)       <- "%"
units(dd3$region)       <- "%"

table1(~ sex +age + age_sampledate + age_at_diag +ethnicity + region| copy_subtype , data = dd3, overall="Total",
       render=rndr, rowlabelhead = "Variables")


#better to remove the gender! 
table1(~age + age_sampledate + age_at_diag +ethnicity + region| copy_subtype , data = dd3, overall="Total",
       render=rndr, rowlabelhead = "Variables")


table1(~  age_at_diag +ethnicity + region| copy_subtype , data = dd3, overall="Total",
       render=rndr, rowlabelhead = "Variables")



tab2 <- table1(~ age_at_diag +ethnicity + region + known_swiss_infect + edu
               + cd4_diag_cat + log_mean_rna | copy_subtype , 
               data = dd3, overall="Total",
               render=rndr, 
               rowlabelhead = "Variables", 
               footnote = "Non-B: A, AE, AG, C, D, F, G, CRFs and unrecognized recombinants",
               caption = "<h3><b> Table 1. Patient characteristics and associations with subtype</b></h3>"
               ); print(tab2)





#### make the same table for White MSM:
w_MSM <- dd3[dd3$ethnicity == "White",] #3996 obs


# better labels and units:
label(w_MSM$age)       <- "Age"
label(w_MSM$region) <- "Geographic region of origin"
label(w_MSM$ethnicity) <- "Ethnicity"
label(w_MSM$age_sampledate) <- "Age at sample date"
label(w_MSM$age_at_diag) <- "Age at diagnosis date"
units(w_MSM$sex)       <- "%"
units(w_MSM$age)       <- "years"
units(w_MSM$age_sampledate)     <- "years"
units(w_MSM$age_at_diag)     <- "years"
units(w_MSM$ethnicity)       <- "%"
units(w_MSM$region)       <- "%"

table1(~ age + age_sampledate + age_at_diag + region| copy_subtype , data = w_MSM, overall="Total",
       render=rndr, rowlabelhead = "Variables", caption = 'MSM Whites')






#### ROW Percentage Table:


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
table.with.row.percentages(~ risk | copy_subtype, data=dd)
table.with.row.percentages(~ ethnicity | copy_subtype, data=dd)
table.with.row.percentages(~ region | copy_subtype, data=dd)



###### Combination of the tables 
m1 <- model.frame(risk ~ copy_subtype, data=dd, na.action=na.pass)
m2 <- model.frame(ethnicity ~ copy_subtype, data=dd, na.action=na.pass)
m3 <- model.frame(region ~ copy_subtype, data=dd, na.action=na.pass)
lab1 <- table1::label(m1[[1]]) #"label of the variable"
lab2 <- table1::label(m2[[1]])
lab3 <- table1::label(m3[[1]])
t1 <- table(m1[[1]], m1[[2]])
t2 <- table(m2[[1]], m2[[2]])
t3 <- table(m3[[1]], m3[[2]])

print(t3)

t1_1<- margin.table(t1, 2) #contingency table in array form, compute the sum of table entries for a given index. 2= col, 1= row
t2_1<- margin.table(t2, 2)
t3_1<- margin.table(t3, 2)

t1_2 <- prop.table(t1, 1) #Express Table Entries as Fraction of Marginal Table, row percentage
t2_2 <- prop.table(t2, 1) 
t3_2 <- prop.table(t3, 1) 

tab1 <- paste0(t1, " (", formatC(100*t1_2, format="f", digits=1), "%)")
attributes(tab1) <- attributes(t1)
tab2 <- paste0(t2, " (", formatC(100*t2_2, format="f", digits=1), "%)")
attributes(tab2) <- attributes(t2)
tab3 <- paste0(t3, " (", formatC(100*t3_2, format="f", digits=1), "%)")
attributes(tab3) <- attributes(t3)


dimnames(tab1)[[2]] <- paste0(dimnames(tab1)[[2]], "<br/>(n=", t1_1, ")")
tab1 <- rbind(rep("", ncol(tab1)), tab1)
dimnames(tab1)[[1]][1] <- lab1
tab2 <- rbind(rep("", ncol(tab2)), tab2)
dimnames(tab2)[[1]][1] <- lab2
tab3 <- rbind(rep("", ncol(tab3)), tab3)
dimnames(tab3)[[1]][1] <- lab3
#print(tab3) 


l <- list(tab1,tab2, tab3)
tab <- do.call('rbind', l)
#htmlTable::htmlTable(tab)


x <- paste0(
  '<table>\n<thead>\n',
  table.rows(colnames(tab), th=T),
  '</thead>\n<tbody>\n',
  table.rows(tab),
  '</tbody>\n</table>\n')
print(x)
structure(x, class=c("table1", "html", "character"), html=TRUE)




