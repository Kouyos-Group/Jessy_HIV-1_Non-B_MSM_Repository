# June 2020, JDR
# Part B - Phylogeny 

### Make Tree
##########################################
# clear R's brain
rm(list=ls())
## packages
library(dplyr)
library(readstata13)
library(ape) #Analyses of Phylogenetics and Evolution
library(phytools)

# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("ggtree")
library(ape)
library(readstata13)
library(ggplot2)
library(ggtree)
library(tidyverse)
library("RColorBrewer")
library("cowplot")


### Load the SHCS sequences, the test_info file and the pat.dta:
setwd("~/Desktop/SHCS")
TEST_INFO <- read.dta13("Input/test_info0620.dta") #24 564
PAT <- read.dta13("Input/pat0620.dta") #20 925
ADMI <- read.dta13("Input/admi0620.dta") #20 925
#data <- read.csv("Output/fulltable_study_population.csv")[,-c(1)]#11'168
data <- read.csv("Output/fulltable_study_population_update.csv")[,-c(1)]#12'852
#setwd("~/Desktop")
#Clusters <-  read.csv("MSMHETcluster_distance15boot80_subtypeA.csv")[,-c(1)]
setwd("~")
Clusters <- read.csv("Downloads/MSMHET_MTCS_above50/MSMHETcluster_GD45B80_A.csv")[,-c(1)]
setwd("~/Desktop/SHCS")




#### new risk group
PAT$riskNEW <- ifelse(is.na(PAT$risk), 'Other',
                      ifelse(PAT$risk == 1, 'MSM', 
                             ifelse(PAT$risk == 2,'HET','Other')))
#ifelse(PAT$risk %in% c(3,4), 'IDU','Other'))))


#### new ethnicity group
sum(is.na(PAT$ethnicity)) #9
PAT$ethNEW <- ifelse(PAT$ethnicity == 1, 'white', 'non-white')


### new region group: 
ADMI$regionNEW <- ifelse(is.na(ADMI$region), 'Unknown',
                         ifelse(ADMI$region %in% c(151,154,155,39), 'Europe', 
                                ifelse(ADMI$region %in% c(21,5,13),'America',
                                       ifelse(ADMI$region %in% c(35,34,30,145),'Asia',
                                              ifelse(ADMI$region %in% c(419,15,18,14,17,11),'Africa',
                                                     ifelse(ADMI$region %in% c(9), 'Oceania','Unknown'))))))



#### load the tree (tree generated with FastTree)
MyTree <- read.tree("Input/A_tree.tre")


## Table with all tree tips and node number
dataTipsandLabel <- data.frame(MyTree$tip.label) 
dataTipsandLabel$tip <- row.names(dataTipsandLabel)
names(dataTipsandLabel) <- c('header','tip')


# Merge SHCS id with header ids
TestInfo_need <- TEST_INFO[,c('header_id', 'id')]
names(TestInfo_need) <- c('header', 'id')
dataTipsandLabel <- merge(dataTipsandLabel, TestInfo_need, all.x = T)
dataTipsandLabel$header <- as.character(dataTipsandLabel$header)
dataTipsandLabel$id <- ifelse(is.na(dataTipsandLabel$id), "LA", as.character(dataTipsandLabel$id))
dataTipsandLabel <- merge(dataTipsandLabel, PAT[,c('id', 'riskNEW', 'ethNEW','sex')], all.x = T)
dataTipsandLabel <- merge(dataTipsandLabel, ADMI[,c('id', 'regionNEW')], all.x = T)
dataTipsandLabel <- merge(dataTipsandLabel, data[,c('id', 'diag_year')], all.x = T)
dataTipsandLabel$sex <- ifelse(dataTipsandLabel$sex == 1, 'male', 'female')
dataTipsandLabel$riskNEW <- ifelse(is.na(dataTipsandLabel$riskNEW), "Other", as.character(dataTipsandLabel$riskNEW))
dataTipsandLabel$ethNEW <- ifelse(is.na(dataTipsandLabel$ethNEW), "Unknown", as.character(dataTipsandLabel$ethNEW))

## We make a second tree (subtree), were we only consider the SHCS sequences.
SHCStips <- which(grepl("[[:alpha:]]",MyTree$tip.label)==FALSE)
AlamosTips <- which(grepl("[[:alpha:]]",MyTree$tip.label))
#LAtips <-  as.numeric(dataTipsandLabel[grepl("[[:alpha:]]", dataTipsandLabel$header),'tip'])
SHCStree <- drop.tip(MyTree,c(AlamosTips,272))



## Table with all tips and node number for SHCS tree
dataTipsandLabelSHCS <- data.frame(SHCStree$tip.label)
dataTipsandLabelSHCS$tip <- row.names(dataTipsandLabelSHCS)
names(dataTipsandLabelSHCS) <- c('header','tipSHCS')
dataTipsandLabelSHCS<- merge(dataTipsandLabelSHCS, TestInfo_need, all.x = T)
dataTipsandLabelSHCS <- merge(dataTipsandLabelSHCS, PAT[,c('id', 'riskNEW', 'ethNEW', 'sex')], all.x = T)
dataTipsandLabelSHCS <- merge(dataTipsandLabelSHCS, ADMI[,c('id', 'regionNEW')], all.x = T)
dataTipsandLabelSHCS<- merge(dataTipsandLabelSHCS, data[,c('id', 'diag_year')], all.x = T)
dataTipsandLabelSHCS$id <- ifelse(is.na(dataTipsandLabelSHCS$id), "LA", as.character(dataTipsandLabelSHCS$id))
dataTipsandLabelSHCS$sex <- ifelse(dataTipsandLabelSHCS$sex == 1, 'male', 'female')
dataTipsandLabelSHCS$riskNEW <- ifelse(is.na(dataTipsandLabelSHCS$riskNEW), "Other", as.character(dataTipsandLabelSHCS$riskNEW))
dataTipsandLabelSHCS$ethNEW <- ifelse(is.na(dataTipsandLabelSHCS$ethNEW), "Unknown", as.character(dataTipsandLabelSHCS$ethNEW))



#change the clusters dataframe:
MSMClusters <- Clusters[Clusters$MSMQ == TRUE,] #19
unique(MSMClusters$clusterNR) 
length(unique(MSMClusters$clusterNR)) #8
MSMClusters$Cl_NR <- ifelse(MSMClusters$clusterNR == 8, 7,
                            ifelse(MSMClusters$clusterNR == 9, 4,
                                   ifelse(MSMClusters$clusterNR == 16, 5,
                                          ifelse(MSMClusters$clusterNR == 44, 1,
                                                 ifelse(MSMClusters$clusterNR == 46, 2,
                                                        ifelse(MSMClusters$clusterNR == 48, 3,
                                                               ifelse(MSMClusters$clusterNR == 54, 6,
                                                                      ifelse(MSMClusters$clusterNR == 66, 8,0))))))))
MSMClusters$Cl_NR <- paste0("Cl_",MSMClusters$Cl_NR)

## preparing dataframe:
dd <- dataTipsandLabelSHCS
dd <- dd[,c(2,1,3:8)]
dd <- na.omit(dd)
dd <- merge(dd, MSMClusters[,c('id', 'Cl_NR', 'clusterType', 'length')], all.x = T)
dd <- dd[,c(2,1,3:11)]
row.names(dd) <- NULL
rownames(dd) <- dd[,1]
dd$riskNEW<-as.factor(dd$riskNEW)
dd$ethNEW<-as.factor(dd$ethNEW)
dd$regionNEW<-as.factor(dd$regionNEW)
dd$Cl_NR <-as.factor(dd$Cl_NR)



#################Preparing color vectors
#For tips
USZcol<-c("black","blue") 
SHCScol<-c("#af8dc3","#7fbf7b", "darkgrey")
#For the linked dataframe(I show three variables, but I have all colors in one vector)
ncol<-c("darkblue","red", "darkgrey",
        "#0073C2BF","#EFC000CC", "#868686FF","#CD534C99","white",
        brewer.pal(8, "Set2"),"white")
# ncol<-c("lightblue","darkred", "darkgrey",
#         "#92c5de","#a6611a", "#80cdc1","#ca0020","white",
#         brewer.pal(8, "Set2"),"white")

names(ncol) = c("HET","MSM", "Other",
                "Africa", "America", "Asia", "Europe","Unknown",
                paste0("Cl_",seq(1,8,1)), "NA")




# brewer.pal(12, "Set3"),brewer.pal(7, "Set2"),"white"
# paste0("Cl",seq(1,19,1))


#ggtree:
##plotting/adjusting the phylogenetic tree
no_legend <- function() theme(legend.position="none")

# SHCStree$tip.label <- paste(dd$diag_year,dd$riskNEW,sep="_")
# gzoom(SHCStree, grep("MSM", SHCStree$tip.label))
plotTree <-ggtree(SHCStree)
plotTree
plotTree +geom_text2(aes(subset=!isTip, label=node), hjust=-.3) + geom_tiplab(size=1)
#plotTree %<+% dd + geom_tiplab()


##Adding tips and cluster labels
p <- plotTree %<+% dd +
  geom_point2(aes(subset=isTip,color=Cl_NR))+ #geom_tiplab(aes(color=Cl_NR))+
  # geom_tippoint(aes(color=Clade),size=2, alpha=.75) +
  #geom_tiplab(aes(label = paste0(diag_year,"_P_",sample(1:100, nrow(dd), replace=TRUE))),size=2.0,align=TRUE,linetype='dashed', linesize=.1, offset = 0.01) +
  geom_treescale(fontsize=5.5,x=-0.08, y=80) +
  scale_color_manual(breaks =  c(paste0("Cl_",seq(1,8,1))) , values = c( brewer.pal(8, "Set2")))+
  no_legend()


p 


##Tree/heatmap plots for each variable separately
pfirst<-gheatmap(p,dd[, "riskNEW", drop=F], offset=0.015, width=.15,colnames=FALSE)+ no_legend()
psecond<-gheatmap(p,dd[, "regionNEW", drop=F], offset=0.005, width=.15,colnames=FALSE)+ no_legend()
pthird<-gheatmap(p,dd[, "Cl_NR", drop=F], offset=0.015, width=.15,colnames=FALSE)+ no_legend()



##Tree/heatmap plots for each variable separately with appropriate colors
p1x<-pfirst + 
  scale_fill_manual(breaks = c("HET","MSM", "Other"), values = c("darkblue","red", "darkgrey"))+ 
  guides(color = FALSE)+
  theme(legend.title = element_text(size=20, face = "bold"),legend.title.align=0,legend.text = element_text(size=20))+ 
  labs(fill = "Risk Group")
p2x<-psecond+
  scale_fill_manual(breaks = c("Africa", "America", "Asia", "Europe", "NA"), values = c("#0073C2BF","#EFC000CC", "#868686FF","#CD534C99", "white"))+ 
  guides(color = FALSE)+
  theme(legend.title = element_text(size=20, face = "bold"),legend.title.align=0,legend.text = element_text(size=20))+ 
  labs(fill = "Continent of origin")
p3x<-pthird+
  scale_fill_manual(breaks =  c(paste0("Cl_",seq(1,8,1))) , values = c( brewer.pal(8, "Set2")))+ 
  guides(color = FALSE)+
  theme(legend.title = element_text(size=20, face = "bold"),legend.title.align=0,legend.text = element_text(size=20))+ 
  labs(fill = "MSM Cluster")




##Getting legends from each of these plots (with cowplot package)
leg1 <- get_legend(p1x)
leg2 <- get_legend(p2x)
leg3 <- get_legend(p3x)


##Combining tree with three variables heatmaps
# pfirstsecond <- gheatmap(pfirst,dd[, "regionNEW", drop=F], offset=0.040, width=.15,colnames=FALSE)
# #scale_fill_manual(values=ncol)+ theme(legend.position="none")
pthirdfirst <- gheatmap(pthird,dd[, "riskNEW", drop=F], offset=0.085, width=.15,colnames=FALSE)+ no_legend()
#scale_fill_manual(values=ncol)+ theme(legend.position="none")
pfirstthirdsecond<- gheatmap(pthirdfirst, dd[, "regionNEW", drop=F], offset=0.155, width=.15,colnames=FALSE)+ no_legend()


##Adding the appropriate colors for a final plot
# all<-pfirstsecond+scale_fill_manual(values=ncol)+ theme(legend.position="none", plot.title = element_text(size = 40, face = "bold"))+ ggtitle("Subtype F")

allA<-pfirstthirdsecond+
  scale_fill_manual(values=ncol)+
  theme(legend.position="none", plot.title = element_text(size = 20, face = "bold"))+
  ggtitle("Subtype A")


##Adding the legends to this final plot
# plot_grid(all, leg1, leg2, ncol=4, rel_widths=c(1.5, .15,.55,.2)) 
Aplot<-plot_grid(allA, leg3, leg1,leg2, ncol=4, rel_widths=c(1.5, .2,.3,.5)) 
ATREE<-plot_grid(allA,  ncol=4, rel_widths=c(4.5, .2,.3,.5)) 



## Save the plot:
# setwd("~/OneDrive/Master Thesis/Figures/Trees")
# ggsave("A_Tree.tiff", width = 50, height = 50, units = "cm", limitsize = FALSE)

