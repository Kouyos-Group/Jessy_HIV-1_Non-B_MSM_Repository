# October 2020, JDR

### Compare cluster vs. non cluster members

############################

### Clear R's memory
rm(list=ls())


### get the packages
if (!requireNamespace("pacman", quietly = TRUE))
  install.packages("pacman")
pacman::p_load("dplyr","tidyr", "data.table", "gmodels", "tidyverse", "boot", "table1","tableone", "devtools", "survival", "vcd")


### Load in all the datas
 # setwd("~/Downloads/compareclusternoncluster")
# setwd("~/Downloads/MTC_nonMTC") #NEW: with MSM>50% and GD: <1.5%, BS: >80%
# setwd("~/Downloads/MTC_NonMTC_D45BS80") #WRONGRESULTS: with MSM>50% and GD: <4.5%, BS: >80%
# 
# for (i in c("A","C","G","F","CRF01_AE","CRF02_AG")){
#   assign(paste0("compareMSM_",i),read.csv(file=paste0(paste0("MSMclusterNONcluster_D45B80_subtype",i),".csv"))[,-1])
#   assign(paste0("compareHET_",i),read.csv(file=paste0(paste0("HETclusterNONcluster_D45B80_subtype",i),".csv"))[,-1])
#   assign(paste0(i,"s"),rbind(eval(parse(text=paste0("compareHET_",i))),eval(parse(text=paste0("compareMSM_",i)))))
#   dataframe <- eval(parse(text=paste0(i,"s")))
#   dataframe$Subtype <- i; assign(paste0(i,"s"),dataframe)
# }


setwd("~")
setwd("Downloads/NEW_MTC_non_MTC") #NEWNEW: with MSM>50% and GD: NO,1.5%,4.5%, BS: NO,>80%


for (i in c("A","C","G","F","CRF01_AE","CRF02_AG")){
  assign(paste0("compareMSM_",i),read.csv(file=paste0(paste0("NEW_MSMclusterNONcluster_",i),".csv"))[,-1])
  assign(paste0("compareHET_",i),read.csv(file=paste0(paste0("NEW_HETclusterNONcluster_",i),".csv"))[,-1])
  MSM <-eval(parse(text=paste0("compareMSM_",i)))
  assign(paste0("compareMSM_",i),MSM %>%
           dplyr::rename( GD45B80 = MSMCluster_GD45B80,
                          GD15B80 = MSMCluster_GD15B80,
                          NoDBCutoff = MSMCluster_NoDBCutoff))
 
  HET <-eval(parse(text=paste0("compareHET_",i)))
  assign(paste0("compareHET_",i),HET %>%
           dplyr::rename( GD45B80 = HETCluster_GD45B80,
                          GD15B80 = HETCluster_GD15B80,
                          NoDBCutoff = HETCluster_NoDBCutoff))
  assign(paste0(i,"s"),rbind(eval(parse(text=paste0("compareHET_",i))),eval(parse(text=paste0("compareMSM_",i)))))
  dataframe <- eval(parse(text=paste0(i,"s")))
  dataframe$Subtype <- i; assign(paste0(i,"s"),dataframe)
}




#### Functions
rndr <- function(x, name, ...) {
  if (!is.numeric(x)) return(render.categorical.default(x))
  what <- switch(name,
                 risk = "",
                 riskNEW = "",
                 ethnicity = "",
                 ethNEW = "",
                 copyregion = "",
                 regroupedregion = "",
                 regionNEW = "",
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



#### get the Tables
for (i in c("A","C","G","F","CRF01_AE","CRF02_AG")){
  assign(paste0("tab_",i)
         ,table1(~  ethNEW + regionNEW +
                   diag_year
                 | riskNEW * Cluster , 
                 data = eval(parse(text=paste0(i,"s"))), overall="Total",
                 render=rndr, 
                 rowlabelhead = "Variables", 
                 footnote = "<b>NOTE</b> All numbers are no. (%) unless otherwise stated <br/>
               ",
                 caption = paste0("<h3><b> Table 1.",i," </b></h3>")))
  print(eval(parse(text=paste0("tab_",i))))
  
  ## Make categorical variables factors
  varsToFactor <- c("ethNEW","regionNEW") 
  # assign(eval(parse(text=paste0(i,"s")))[varsToFactor], 
  #        lapply(eval(parse(text=paste0(i,"s")))[varsToFactor], factor))
  # 
  
  ## Create a variable list
  vars <- c("ethNEW", "regionNEW", 
            "diag_year")
  
  
  ## Create Table 1 stratified by risk and cluster/phylogeny
  dataframe <- eval(parse(text=paste0(i,"s")))
  assign(paste0("tableOne_",i),
         CreateTableOne(vars = vars, strata = c("riskNEW","Cluster"), data = dataframe))
  # assign(paste0("tableOne_",i),
  #        CreateTableOne(vars = vars, strata = c("riskNEW","GD45B80"), data = dataframe))
  
  
  ## Specifying nonnormal variables will show the variables appropriately,
  ## and show nonparametric test p-values. Specify variables in the exact
  ## argument to obtain the exact test p-values. cramVars can be used to
  ## show both levels for a 2-level categorical variables.
  ## If your work flow includes copying to Excel and Word when writing manuscripts,
  ## you may benefit from the quote argument. This will quote everything so that
  ## Excel does not mess up the cells.
  ## If you want to center-align values in Word, use noSpaces option.
  
  print(eval(parse(text=paste0("tableOne_",i))), 
        nonnormal = c("diag_year"),
        exact = c("regionNEW"), cramVars = "ethNEW", 
        smd = FALSE, quote = FALSE, noSpaces = TRUE)
  
}
  
# ## Use the summary.TableOne method for detailed summary
# summary(tableOne_A)
# 
# ## See the categorical part only using $ operator
# tableOne_A$CatTable
# summary(tableOne_A$CatTable)
# 
# ## See the continuous part only using $ operator
# tableOne_A$ContTable
# summary(tableOne_A$ContTable)
# 
# ## If SMDs are needed as numericals, use ExtractSmd()
# ExtractSmd(tableOne_A)


###### Show the tables as mosaicplots: 
ALL <- rbind(As,Cs, Fs,Gs,CRF01_AEs,CRF02_AGs)
#Creating a frequency table:
#table <- xtabs(~ regionNEW+ Cluster + riskNEW + Subtype , data = ALL)
table <- xtabs(~ regionNEW+  GD45B80 + riskNEW  + Subtype , data = ALL) #to get the value for noncluster's vs. cluster's by risk group and continent
prop.table(table,1) #express Table Entries as Fraction of Marginal Table
# cotabplot(~ riskNEW+ regionNEW + Cluster |Subtype, data = table, split_vertical = TRUE)
cotabplot(~ riskNEW + regionNEW + GD45B80 |Subtype, data = table, split_vertical = TRUE)

mycol <- c("#0073C2BF","#0073C2BF","#EFC000CC","#EFC000CC", "#868686FF","#868686FF","#CD534C99","#CD534C99","white","white")
cotabplot(~  riskNEW + regionNEW + Cluster | Subtype, data = table, 
          split_vertical = TRUE,  #to get the risk above
          labeling_args = list(rot_labels = c( left = 0),
                               offset_varnames = c(left= 0.8),
                               offset_labels = c(left = 0.3)),
          margins = c(left = 0.1),
          gp = gpar(fill =mycol), col = mycol,
          spacing = spacing_equal(unit(0.6, "lines")),
          labeling = labeling_values, value_type = "observed", gp_text = gpar(fontface = 2),
          newpage = TRUE) #to get the absolute value printed in the cell

cotabplot(~  riskNEW+ regionNEW + GD45B80 | Subtype, data = table, 
          split_vertical = TRUE,  #to get the risk above
          labeling_args = list(rot_labels = c( left = 0),
                               offset_varnames = c(left= 0.8),
                               offset_labels = c(left = 0.3)),
          margins = c(left = 0.1),
          gp = gpar(fill =mycol), col = mycol,
          spacing = spacing_equal(unit(0.6, "lines")),
          labeling = labeling_values, value_type = "observed", gp_text = gpar(fontface = 2),
          newpage = TRUE) #to get the absolute value printed in the cell


###make a mosIC PLOT   RENAMED
ALL2 <- ALL
ALL2 <- ALL2 %>% #Rename by taking out NEW
  dplyr::rename_all(
    funs(stringr::str_replace_all(., "NEW", ""))) #want to “delete” the variable names’ suffix.


#Region:
levels(ALL2$region)
ALL2$region <- plyr::mapvalues	(ALL2$region, from = c("Africa", "America","Asia","Europe","Unknown"), to =c("AF","AM","AS","EU","NA"))
levels(ALL2$region)


ALL2 <- ALL2 %>% #Rename 
  dplyr::rename(., Cluster = GD45B80)  ##change to D15B80
ALL2$Cluster <- plyr::mapvalues	(ALL2$Cluster, from = c(0,1), to = c("N","Y"))

ALL2 <- ALL2 %>% #Rename by taking out NEW
  dplyr::rename(., Risk = risk)
ALL2 <- ALL2 %>% #Rename by taking out NEW
  dplyr::rename(., Region = region)
test <- subset(ALL2, Region != "NA")
test$Region <- plyr::mapvalues	(test$Region, from = c("AF", "AM","AS","EU","NA"), to =c("AF","AM","AS","EU"," "))
test$Subtype <- factor(test$Subtype, levels=c("A","C","G","F","CRF01_AE","CRF02_AG"))


table2 <- xtabs(~ Region+ Cluster + Risk + Subtype , data = test) #to get the value for noncluster's vs. cluster's by risk group and continent
prop.table(table2,1) #express Table Entries as Fraction of Marginal Table
cotabplot(~ Risk+ Region + Cluster |Subtype, data = table2, split_vertical = TRUE)
mycol2 <- c("#0073C299","#0073C299","#EFC000CC","#EFC000CC", "#868686FF","#868686FF","#CD534CBF","#CD534CBF","white","white") #with All2
# mycol2 <- c("#0073C299","#0073C299","#EFC000CC","#EFC000CC", "#868686FF","#868686FF","#CD534CBF","#CD534CBF")

display.brewer.pal(n = 8, name = 'Dark2')
brewer.pal(n = 8, name = "Dark2")

display.brewer.pal(n = 8, name = 'Set2')
brewer.pal(n = 8, name = "Set2")

mycol2 <- c("#0073C299","#0073C299",
            "#EFC000CC","#EFC000CC",
            "#868686FF","#868686FF",
            "#CD534CBF","#CD534CBF",
            "white","white",
            "#0073C2BF","#0073C2BF",
            "#E6AB02","#E6AB02", 
            "#666666","#666666",
            "darkred","darkred",
            "white","white")


mycol2 <- c("#92c5de","#92c5de", #AF,N
            "#dfc27d","#dfc27d",
            "#80cdc1","#80cdc1",
            "#f4a582","#f4a582", #EU
            "white","white",
            "#0571b0","#0571b0", #AF,Y
            "#a6611a","#a6611a", 
            "#018571","#018571",
            "#ca0020","#ca0020", #EU
            "white","white")



setwd("~/OneDrive/Master Thesis/Figures/ClustervsRisk")
# setEPS()
# postscript("NEWRiskvsCluster_panelB.eps")
cotabplot(~ Risk+ Region + Cluster | Subtype, data = table2, 
          split_vertical = TRUE,  #to get the risk above
          labeling_args = list( rot_labels = c( left = 0),
                                offset_varnames = c(left= 0.2),
                                offset_labels = c(left = 0.0)),
          margins = c(left = 0.1),
          gp = gpar(fill =mycol2), col = mycol2,
          spacing = spacing_equal(unit(0.6, "lines")),
          labeling = labeling_values, value_type = "observed", gp_text = gpar(fontface = 2),
          newpage = TRUE) #to get the absolute value printed in the cell
# dev.off()


##############################################################
#Cluster vs. Risk:

# mycolRiskCluster <- c("#af8dc3","#7fbf7b",
#                       "#762a83","#1b7837")
mycolRiskCluster <- c("#a6dba0", "#a6dba0",
                      "#af8dc3","#af8dc3")

table3 <- xtabs(~  Cluster + Risk + Subtype , data = test) #to get the value for noncluster's vs. cluster's by risk group 
prop.table(table3,1) #express Table Entries as Fraction of Marginal Table
cotabplot(~ Risk + Cluster |Subtype, data = table3, split_vertical = TRUE)
# setEPS()
# postscript("NEWRiskvsCluster_panelA2.eps")
cotabplot(~ Risk + Cluster | Subtype, data = table3, 
          split_vertical = T,  #to get the risk above
          labeling_args = list( rot_labels = c( left = 0),
                                offset_varnames = c(left= 0.2),
                                offset_labels = c(left = 0.0)),
          margins = c(left = 0.1),
          gp = gpar(fill =mycolRiskCluster), col = mycolRiskCluster,
          spacing = spacing_equal(unit(0.15, "lines")),
          labeling = labeling_values, value_type = "observed", gp_text = gpar(fontface = 2),
          newpage = T) #to get the absolute value printed in the cell
# dev.off()
###save plot 500x700



#### two sided Fisher extact test
Risk_Cluster_table <- table(test$Risk, test$Cluster)
prop.table(Risk_Cluster_table,1) #express Table Entries as Fraction of Marginal Table
mosaicplot(Risk_Cluster_table, color = c("darkred", "gold"), xlab ="Risk", ylab = "Cluster")
Risk_Cluster_fisher =  fisher.test(Risk_Cluster_table)
Risk_Cluster_fisher$estimate #odd ratio
Risk_Cluster_fisher$conf.int #OR CI
round(Risk_Cluster_fisher$p.value, 3)#assocation between risk and cluster!



Risk_Cluster_Subtype_table <- table(test$Risk, test$Cluster, test$Subtype)
mosaicplot(Risk_Cluster_Subtype_table, color = c("darkred", "gold"), xlab ="Risk", ylab = "Cluster")

for (i in c("A","C","G","F","CRF01_AE","CRF02_AG")){
  dataframe <- eval(parse(text=paste0(i,"s")))
  assign(paste0("Risk_Cluster_",i), table(dataframe$riskNEW, dataframe$GD45B80)) ##change to D15B80
  datafisher <- eval(parse(text=paste0("Risk_Cluster_",i)))
  assign(paste0("Risk_Cluster_fisher_",i), fisher.test(datafisher))
  print(eval(parse(text=paste0("Risk_Cluster_fisher_",i))))
}
round(Risk_Cluster_fisher_A$p.value, 3)
round(Risk_Cluster_fisher_C$p.value, 3)
round(Risk_Cluster_fisher_G$p.value, 3)
round(Risk_Cluster_fisher_F$p.value, 3)
round(Risk_Cluster_fisher_CRF01_AE$p.value, 3)
round(Risk_Cluster_fisher_CRF02_AG$p.value, 3)

############################################################################################
#### make the mosaik plot with renamed variables and only stratified by EU vs. not-EU
ALL_NEW <- ALL
ALL_NEW <- ALL_NEW %>% #Rename by taking out NEW
  dplyr::rename_all(
    funs(stringr::str_replace_all(., "NEW", ""))) #want to “delete” the variable names’ suffix.

#Region Europe or Not-Europe:
levels(ALL_NEW$region)
ALL_NEW$region <- mapvalues(ALL_NEW$region, from = c("Africa", "America","Asia","Unknown"), to = rep("Not-Europe",4))
levels(ALL_NEW$region)



ALL_NEW$Cluster <- mapvalues(ALL_NEW$Cluster, from = c(0,1), to = c("N","Y"))


table_reduced <- xtabs(~ region+ Cluster + risk + Subtype , data = ALL_NEW) #to get the value for noncluster's vs. cluster's by risk group and continent
prop.table(table_reduced,1) #express Table Entries as Fraction of Marginal Table
cotabplot(~ risk+ region + Cluster |Subtype, data = table_reduced, split_vertical = TRUE)
mycol_reduced <- c("#0073C299","#0073C299","#CD534CBF","#CD534CBF")
cotabplot(~  risk+ region + Cluster | Subtype, data = table_reduced, 
          split_vertical = TRUE,  #to get the risk above
          labeling_args = list( offset_varnames = c(left= 0.2),
                               offset_labels = c(left = 0.0)),
          margins = c(left = 0.1),
          gp = gpar(fill =mycol_reduced), col = mycol_reduced,
          spacing = spacing_equal(unit(0.6, "lines")),
          labeling = labeling_values, value_type = "observed", gp_text = gpar(fontface = 2),
          newpage = TRUE) #to get the absolute value printed in the cell


lor_reduced <- oddsratio(table_reduced)  #log odds ratios for regionNEW and Cluster by riskNEW, Subtype 
summary(lor_reduced)
plot(lor_reduced, xlab="Continent", 
     ylab="Log Odds Ratio (Continent | Cluster)",
     legend("bottomleft")) #Log odds ratio plot- by default, it shows the 95% confidence interval for the log odds ratio.

