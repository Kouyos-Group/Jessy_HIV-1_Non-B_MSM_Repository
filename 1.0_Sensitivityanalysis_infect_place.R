
#### JDR
#### FEB 21
### Sensitivity analysis: Information on the likely location of transmission 

### Clear R's memory
rm(list=ls())


### get the packages
if (!requireNamespace("pacman", quietly = TRUE))
  install.packages("pacman")
pacman::p_load("dplyr","tidyr","plyr", "data.table", "gmodels", "tidyverse", "boot", "table1","tableone", "devtools", "survival", "vcd")

### Load the data
setwd("~/Desktop/SHCS")
data <- read.csv("Output/fulltable_study_population_update.csv")
data <- data[,-1] #12852

### infect_place:
# 1: Switzerland 
# 2: while temporarily abroad
# 3: as a resident abroad 
# 9: unknown


FUN_Switzerland <- function(x){
  if (is.na(x)){return (NA)}
  else if (x == 9){return(NA)}
  else if (x == 1){return("Switzerland")}
  else if (x == 2){return("Temporarily abroad")}
  else if (x == 3){return("Resident abroad")}
}
data$known_swiss_infect <- sapply(data$infect_place, FUN_Switzerland)
data$known_swiss_infect <- factor(data$known_swiss_infect, levels = c("Switzerland","Temporarily abroad","Resident abroad"))
#data <- within(data, known_swiss_infect <- relevel(known_swiss_infect , ref = "Switzerland"))
data$abroad_infect <- ifelse(data$known_swiss_infect == "Temporarily abroad", "Abroad",
                             ifelse(data$known_swiss_infect == "Resident abroad", "Abroad",0))

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


## subset the data with information on the likely location of transmission
subdata <- subset(data, infect_place >0) #4320
table(subdata$infect_place)
# 1    2    3    9 
# 1966  636  895  823 
CrossTable(subdata$known_swiss_infect)
# Not-Swiss     Swiss 
# 1531(43.8%)      1966  (56.2%)


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


label(data$copyregion) <- "Geographic region of origin" 
label(data$regroupedregion) <-"Geographic region of origin"

tab <- table1(~ regroupedregion 
              | known_swiss_infect , 
              data = data, 
              overall=F,
              render=rndr, 
              rowlabelhead = "Variables", 
              footnote = "<b>NOTE</b> All numbers are no. (%) unless otherwise stated",
              caption = "<h3><b> Sensitivity analysis: Information on the likely location of transmission</b></h3>")
print(tab)


tab2 <- table1(~ copyregion 
               | known_swiss_infect , 
               data = data, 
               overall=F,
               render=rndr, 
               rowlabelhead = "Variables", 
               footnote = "<b>NOTE</b> All numbers are no. (%) unless otherwise stated",
               caption = "<h3><b> Sensitivity analysis: Information on the likely location of transmission</b></h3>")
print(tab2)


labels1 <- list(
  variables=list(copyregion="Geographic region of origin"),
  groups=list("","Abroad")
)

labels2 <- list(
  variables=list(regroupedregion="Geographic region of origin"),
  groups=list("","Abroad")
)
levels(data$known_swiss_infect) <- c("Switzerland","Temporarily","Resident") #remove the word "abroad"
strata <- c(split(data,data$known_swiss_infect), list("All abroad"=subset(data, abroad_infect=="Abroad"))) #list(Overall=data))



table1(strata, 
       labels1,
       groupspan=c(1,3),
       data = data, 
       overall=F,
       render=rndr, 
       rowlabelhead = "Variables", 
       footnote = "<b>NOTE</b> All numbers are no. (%) unless otherwise stated",
       caption = "<h3><b> Sensitivity analysis: Geographic region of origin stratified by the information on the likely location of transmission</b></h3>")


table1(strata, 
       labels2,
       groupspan=c(1,3),
       data = data, 
       overall=F,
       render=rndr, 
       rowlabelhead = "Variables", 
       footnote = "<b>NOTE</b> All numbers are no. (%) unless otherwise stated",
       caption = "<h3><b> Sensitivity analysis: Geographic region of origin stratified by the information on the likely location of transmission</b></h3>")

#### MSM
MSM <- data[data$risk == "1",]  
library(tidyr)
MSM <- MSM %>% drop_na(id) #5215

MSM_6NB<- filter(MSM, subtype=="A" | subtype=="02_AG" 
                 | subtype=="C" | subtype=="01_AE"
                 | subtype=="F"| subtype=="G"| subtype=="B") ;unique(MSM_6NB$subtype) #5114


strata_MSM <- c(split(MSM_6NB,MSM_6NB$known_swiss_infect), list("All abroad"=subset(MSM_6NB, abroad_infect=="Abroad"))) #list(Overall=data))

table1(strata_MSM, 
       labels1,
       groupspan=c(1,3),
       data = MSM, 
       overall=T,
       render=rndr, 
       rowlabelhead = "Variables", 
       footnote = "<b>NOTE</b> All numbers are no. (%) unless otherwise stated",
       caption = "<h3><b> Sensitivity analysis: Geographic region of origin stratified by the information on the likely location of transmission</b></h3>")



#### HET
HET <- data[data$risk == "2",]  
HET <- HET %>% drop_na(id) #4563

HET_6NB<- filter(HET, subtype=="A" | subtype=="02_AG" 
                 | subtype=="C" | subtype=="01_AE"
                 | subtype=="F"| subtype=="G"| subtype=="B") ;unique(HET_6NB$subtype) #5114



strata_HET <- c(split(HET_6NB,HET_6NB$known_swiss_infect), list("All abroad"=subset(HET_6NB, abroad_infect=="Abroad"))) #list(Overall=data))

table1(strata_HET, 
       labels1,
       groupspan=c(1,3),
       data = HET, 
       overall=T,
       render=rndr, 
       rowlabelhead = "Variables", 
       footnote = "<b>NOTE</b> All numbers are no. (%) unless otherwise stated",
       caption = "<h3><b> Sensitivity analysis: Geographic region of origin stratified by the information on the likely location of transmission</b></h3>")



############################################################################################
#### correlation between infection place and region of origin by subtype$?
ALL_NEW <- data 
ALL_NEW<-filter(ALL_NEW, subtype=="A" | subtype=="02_AG" 
       | subtype=="C" | subtype=="01_AE"
       | subtype=="F"| subtype=="G")

ALL_NEW<-filter(ALL_NEW, risk=="1" | risk=="2")
levels(ALL_NEW$known_swiss_infect)
ALL_NEW$Infection_place <- mapvalues(ALL_NEW$known_swiss_infect, from = c("Switzerland" ,"Temporarily" ,"Resident" ), to = c("Switzerland",rep("Abroad",2)))

#Region Europe or Not-Europe:
levels(ALL_NEW$copyregion)
ALL_NEW <- subset(ALL_NEW, copyregion !=  "Unknown")
# ALL_NEW$regionclass <- mapvalues(ALL_NEW$copyregion, 
#                                  from = c("Western Europe", 
#                                           "Southern Europe","Northern Europe","Eastern Europe",
#                                           "Northern America","South America" ,"Latin America and the Carribean", "Central America"  ,"Oceania" ,
#                                           "South-Eastern Asia" , "Western Asia"  ,"Southern Asia"  ,"Eastern Asia"  ,
#                                           "Western Africa" ,"Northern Africa" ,"Eastern Africa","Southern Africa" ,"Middle Africa"),
#                                  to = c("W_EU",rep("R_EU",3),rep("AM",5),rep("AS",4),rep("AF",5)))
# ALL_NEW$regionclass <- factor(ALL_NEW$regionclass, levels=c("AF","AM","AS","R_EU","W_EU"))
ALL_NEW$Region <- mapvalues(ALL_NEW$copyregion, 
                                 from = c("Western Europe", 
                                          "Southern Europe","Northern Europe","Eastern Europe",
                                          "Northern America","South America" ,"Latin America and the Carribean", "Central America"  ,"Oceania" ,
                                          "South-Eastern Asia" , "Western Asia"  ,"Southern Asia"  ,"Eastern Asia"  ,
                                          "Western Africa" ,"Northern Africa" ,"Eastern Africa","Southern Africa" ,"Middle Africa"),
                                 to = c("W_EU",rep("Others",17)))
ALL_NEW$Region <- factor(ALL_NEW$Region, levels=c("W_EU","Others"))
levels(ALL_NEW$Region)


ALL_NEW$Subtype <- factor(ALL_NEW$subtype, levels=c("A","C","G","F","CRF01_AE","CRF02_AG"))
levels(ALL_NEW$Subtype)

labels_ALL <- list(
  variables=list(Region="Geographic region of origin"),
  groups=list("","Abroad")
)

strata_ALL <- c(split(ALL_NEW,ALL_NEW$Infection_place))#, list("All abroad"=subset(ALL_NEW, abroad_infect=="Abroad"))) #list(Overall=data))

table1(strata_ALL, 
       labels_ALL,
       #groupspan=c(1,1),
       data = ALL_NEW, 
       overall=T,
       render=rndr, 
       rowlabelhead = "Variables", 
       footnote = "<b>NOTE</b> All numbers are no. (%) unless otherwise stated",
       caption = "<h3><b> Sensitivity analysis: Geographic region of origin stratified by the information on the likely location of transmission</b></h3>")




##changing to factors
ALL <- ALL_NEW  %>% select(Region,Infection_place, Subtype)
must_convert <- sapply(ALL, is.factor)
DATAFRAME <- sapply(ALL[must_convert], unclass)
str(DATAFRAME)

library(ggcorrplot)
#install.packages("ggcorrplot")
corr <- round(cor(DATAFRAME),1)
head(corr)
p.mat <- cor_pmat(DATAFRAME)
head(p.mat)

ggcorrplot(p.mat)


#Infection<-ftable(ALL)

table_reduced <- xtabs(~ Region+ Infection_place + Subtype , data = ALL) #to get the value for infection place by  region 

mantelhaen.test(table_reduced)
oddsratio(table_reduced, log = F) ##violation of no 3 way association
lor <- oddsratio(table_reduced) #capture log odds ratios
summary(lor)

woolf_test(table_reduced) #examine each 2x2 table 

fourfoldplot(table_reduced,mfrow = c(2,3))


vars<- c("Region","Infection_place", "Subtype")
assocstats(table(ALL))
catcorrm <- sapply(vars, function(y) sapply(vars, function(x) assocstats(table(ALL[,x], ALL[,y]))$cramer))
ggcorrplot(catcorrm)


#table_reduced <- xtabs(~ Region+ Infection_place + Subtype , data = ALL) #to get the value for infection place by  region 
prop.table(table_reduced,1) #express Table Entries as Fraction of Marginal Table
cotabplot(~ Region+ Infection_place |Subtype, data = table_reduced, split_vertical = TRUE)
mycol_reduced <- c("darkred","darkred", "gold","gold")
long.labels <- list(set_varnames = c(Region="Region", Infection_place="Infection Place"))
cotabplot(~   Infection_place + Region |Subtype,  data = table_reduced, 
          split_vertical = TRUE,  #to get the infection place above
          labeling_args = long.labels, offset_varnames = c(left= 0.2, top= 0.2),
          margins = c(left = 0.1),
          gp = gpar(fill =mycol_reduced), col = mycol_reduced,
          spacing = spacing_equal(unit(0.6, "lines")),
          labeling = labeling_values, value_type = "observed", gp_text = gpar(fontface = 2),
          newpage = TRUE) #to get the absolute value printed in the cell



table_all <- xtabs(~  Infection_place + Region, data = ALL) #to get the value for infection place by  region 
prop.table(table_all,1) #express Table Entries as Fraction of Marginal Table
cotabplot(~ Region+ Infection_place , data = table_all, split_vertical = TRUE)
mycol_reduced <- c("darkred","darkred", "gold","gold")
cotabplot(~  Infection_place+ Region,  data = table_all, 
          split_vertical = TRUE,  #to get the infection place above
          labeling_args = long.labels, offset_varnames = c(left= 0.2, top= 0.2),
          offset_labels = c(left = 0.0),
          margins = c(left = 0.2, top= 0.2),
          gp = gpar(fill =mycol_reduced), col = mycol_reduced,
          spacing = spacing_equal(unit(0.4, "lines")),
          labeling = labeling_values, value_type = "observed", gp_text = gpar(fontface = 4),
          newpage = TRUE) #to get the absolute value printed in the cell



#### two sided Fisher extact test
mosaicplot(table_all, color = c("darkred", "gold"), xlab ="Region", ylab = "Infection place")
Infection_Region_fisher =  fisher.test(table_all)
Infection_Region_fisher$estimate #odd ratio
Infection_Region_fisher$conf.int #OR CI
Infection_Region_fisher$p.value#assocation between Infection and Region!



for (i in c("A","C","G","F","CRF01_AE","CRF02_AG")){
  assign(paste0("ALL_",i), ALL[ALL$Subtype==i,]) 
  assign(paste0("ALL_",i),  eval(parse(text=paste0("ALL_",i)))[!(is.na(eval(parse(text=paste0("ALL_",i)))$Subtype)),]) 
  dataframe <- eval(parse(text=paste0("ALL_",i)))
  assign(paste0("Infection_Region_",i), xtabs(~  Infection_place + Region, data = dataframe)) 
  datafisher <- eval(parse(text=paste0("Infection_Region_",i)))
  assign(paste0("Infection_Region_fisher_",i), fisher.test(datafisher))
  #print(eval(parse(text=paste0("Infection_Region_fisher_",i))))
}
round(Infection_Region_fisher_A$p.value, 3)
round(Infection_Region_fisher_C$p.value, 3); Infection_Region_fisher_C$p.value
round(Infection_Region_fisher_G$p.value, 3); Infection_Region_fisher_G$p.value
round(Infection_Region_fisher_F$p.value, 3)
round(Infection_Region_fisher_CRF01_AE$p.value, 3)
round(Infection_Region_fisher_CRF02_AG$p.value, 3)

# 
# # example dataframe
# df <- data.frame(
#   group = c('A', 'A', 'A', 'A', 'A', 'B', 'C'),
#   student = c('01', '01', '01', '02', '02', '01', '02'),
#   exam_pass = c('Y', 'N', 'Y', 'N', 'Y', 'Y', 'N'),
#   subject = c('Math', 'Science', 'Japanese', 'Math', 'Science', 'Japanese', 'Math')
# )
# 
# library(tidyverse)
# library(lsr)
# 
# # function to get chi square p value and Cramers V
# f = function(x,y) {
#   tbl = ALL %>% select(x,y) %>% table()
#   chisq_pval = round(chisq.test(tbl)$p.value, 4)
#   cramV = round(cramersV(tbl), 4)
#   data.frame(x, y, chisq_pval, cramV) }
# 
# # create unique combinations of column names
# # sorting will help getting a better plot (upper triangular)
# df_comb = data.frame(t(combn(sort(names(ALL)), 2)), stringsAsFactors = F)
# 
# # apply function to each variable combination
# df_res = map2_df(df_comb$X1, df_comb$X2, f)
# 
# 
# # plot results
# df_res %>%
#   ggplot(aes(x,y,fill=chisq_pval))+
#   geom_tile()+
#   geom_text(aes(x,y,label=cramV))+
#   scale_fill_gradient(low="red", high="yellow")+
#   theme_classic()
# 
