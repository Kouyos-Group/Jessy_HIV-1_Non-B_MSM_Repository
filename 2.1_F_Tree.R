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
#Clusters <-  read.csv("MSMHETcluster_distance15boot80_subtypeF.csv")[,-c(1)]
setwd("~")
Clusters <- read.csv("Downloads/MSMHET_MTCS_above50/MSMHETcluster_GD45B80_F.csv")[,-c(1)]
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
MyTree <- read.tree("Input/F_tree.tre")


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
MSMClusters <- Clusters[Clusters$MSMQ == TRUE,] #20
unique(MSMClusters$clusterNR) #5
length(unique(MSMClusters$clusterNR))
MSMClusters$Cl_NR <- ifelse(MSMClusters$clusterNR == 1, 1,
                            ifelse(MSMClusters$clusterNR == 4, 5, 
                                   ifelse(MSMClusters$clusterNR == 11, 2, 
                                          ifelse(MSMClusters$clusterNR == 13, 4, 
                                                 ifelse(MSMClusters$clusterNR == 14, 3, 0)))))
# MSMClusters$Cl_NR <- ifelse(MSMClusters$clusterNR == 1, 1,
#                             ifelse(MSMClusters$clusterNR == 7, 5, 
#                                    ifelse(MSMClusters$clusterNR == 8, 3, 
#                                           ifelse(MSMClusters$clusterNR == 9, 2, 
#                                                  ifelse(MSMClusters$clusterNR == 10, 4, 
#                                                         ifelse(MSMClusters$clusterNR == 12, 7,
#                                                                ifelse(MSMClusters$clusterNR == 13, 6, 
#                                                                       ifelse(MSMClusters$clusterNR == 14, 8,0))))))))
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
  geom_treescale(fontsize=5.5,x=-0.05, y=20) +
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
pthirdfirst <- gheatmap(pthird,dd[, "riskNEW", drop=F], offset=0.038, width=.15,colnames=FALSE)+ no_legend()
#scale_fill_manual(values=ncol)+ theme(legend.position="none")
pfirstthirdsecond<- gheatmap(pthirdfirst, dd[, "regionNEW", drop=F], offset=0.062, width=.15,colnames=FALSE)+ no_legend()


##Adding the appropriate colors for a final plot
# all<-pfirstsecond+scale_fill_manual(values=ncol)+ theme(legend.position="none", plot.title = element_text(size = 40, face = "bold"))+ ggtitle("Subtype F")

allF<-pfirstthirdsecond+
  scale_fill_manual(values=ncol)+
  theme(legend.position="none", plot.title = element_text(size = 20, face = "bold"))+
  ggtitle("Subtype F")


##Adding the legends to this final plot
# plot_grid(all, leg1, leg2, ncol=4, rel_widths=c(1.5, .15,.55,.2)) 
Fplot <- plot_grid(allF, leg3, leg1,leg2, ncol=4, rel_widths=c(1.5, .2,.3,.5)) 
FTREE <- plot_grid(allF, ncol=4, rel_widths = c(4.5,.2,.3,.5))

## Save the plot:
# setwd("~/OneDrive/Master Thesis/Figures/Trees")
# ggsave("F_Tree.tiff", width = 50, height = 50, units = "cm", limitsize = FALSE)
# 
# 

##############################
##############################
##############################

gzoom(pfirstsecond, grep("MSM", SHCStree$tip.label), xmax_adjust=2)

#8) Get inner-node names (branch or clade names):
clades <- MyTree$node.label
clades <- clades[clades != ""]  # remove empty node-names 




###########################
## clear R's brain
rm(list=ls())






### Load the SHCS sequences, the test_info file and the pat.dta:
setwd("~/Desktop/SHCS")
# SHCS_SEQ <- read.fasta(file="final0620.fas") #29 837 
TEST_INFO <- read.dta13("Input/test_info0620.dta") #24 564
PAT <- read.dta13("Input/pat0520.dta") #20 904
F_info <- read.csv("Input/OLD/F_info.csv")
F_info$X <- NULL
data <- read.csv("Output/fulltable_study_population.csv")



# plot tree for visualization
MyTree <- read.tree("Input/OLD/FirstTree.tre")
plot(MyTree)
int_nodes <- unique(MyTree$edge[,1])
nr_ext_nodes <- length(MyTree$tip.label) #127

## Drop all Los Alamos tips
SHCStips <- which(grepl("[[:alpha:]]",MyTree$tip.label)==FALSE)
AlamosTips <- which(grepl("[[:alpha:]]",MyTree$tip.label))
MyTreeSHCS <- drop.tip(MyTree,AlamosTips)
# write.tree(MyTreeSHCS,"MyTreeSHCS_F.tre")
##Phylogenetic tree with 127 tips and 125 internal nodes:
F_Tree <- read.tree("Input/OLD/MyTreeSHCS_F.tre")


par(mfrow=c(1,1))
plot(F_Tree, type= "cladogram")
plot(F_Tree, type= "fan")
plot(F_Tree, type= "unrooted")
summary(F_Tree)


x <- as_tibble(F_Tree)




information <- merge(PAT[c("id","sex","ethnicity","risk","hiv_posdate", "hiv_posdocdate", "hiv_negdate")],F_info, all.x=F, by = "id")
FUN_risk_cat <- function (x){
  if (is.na(x)){return (NA)}
  else if (x == 1){return ("MSM")}
  else if (x == 2){return ("HET")}
  else if (x == 3){return ("IDU")}
  else if (x == 4){return ("IDU")}
  else if (x == 5){return ("Blood Transfusion/ Products")}
  else if (x == 6){return ("Blood Transfusion/ Products")}
  else if (x == 7){return ("Others")}
  else if (x == 9){return ("Others")}
  else if (x == 0){return ("Others")}
  else {return (as.character(x))}
}
information$risk_cat <- sapply(information$risk, FUN_risk_cat)

information <- merge(information,data[c("id","diag_year","age_at_diag","region")], all.y=F,all.x = T, by.x = "id")

information$eth_cat <- ifelse(information$ethnicity == 1, 'white', 'non-white')


x_with_info <- merge(x,information[c("header_id","risk_cat","diag_year","eth_cat")], by.x = "label", by.y = "header_id", all.y = F, all.x = F)
x_with_info <- as_tibble(x_with_info)
x_with_info$diag_year[46] <- 1986
x_with_info$diag_year[47] <- 2014

# as.treedata(x_with_info)
# plot(x_with_info)
# 
# setwd("~/Desktop/SHCS/R_Sequences/Sequences")
# source("my.plot.phylo.r")
# my.plot.phylo(F_Tree, align.tip.label= T, label.offset = 1, 
#               no.margin = T, cex= .8, font=2, type = "fan")
# 
# my.plot.phylo(F_Tree)

#Set up margins
par(mar = c(2,2,2,2))
plot(F_Tree, align.tip.label= F, edge.width = 2, label.offset = 0, cex = 0.4)

F_Tree1 <- F_Tree
F_Tree2 <- F_Tree
F_Tree3 <- F_Tree
F_Tree2$tip.label <- paste(x_with_info$risk_cat,x_with_info$diag_year,sep=" | ")
F_Tree1$tip.label <- x_with_info$risk_cat
F_Tree3$tip.label <- x_with_info$eth_cat
plot(F_Tree1, align.tip.label= F, edge.width = 2, label.offset = 0, cex = 0.4)


# Change edge color
# generate edge colors:
# mycols <- c("orange", "purple", "blue")
# risk_col <- edge.color(F_Tree1, label, col = mycols)

# color=c('red','blue','green')
# nrisk= length(unique(F_Tree1$tip.label))
# color_list=rep(color,nrisk/length(color))
# # plot tree:
# plot(F_Tree1, type = "fan",
#      no.margin = TRUE, edge.width = 2,label.offset = 0,cex = 0.6,
#      tip.color = color_list[nrisk], front=2)

## make unique rownames (equal rownames are not allowed)
rownames(F_Tree1) <- unique(F_Tree1$tip.label)

colorCodes <- c(M="red", H="violet", O="darkgrey")
lineType <- c(w=1,n=2)

plot(F_Tree1, tip.color=colorCodes, type="fan")
plot(F_Tree2, tip.color=colorCodes[substr(F_Tree1$tip.label, 1, 1)],cex = 0.6, edge.width = 2,type = "fan")
plot(F_Tree2, tip.color=colorCodes[substr(F_Tree1$tip.label, 1, 1)],cex = 0.6, edge.width = 2)

plot(F_Tree2, tip.color=colorCodes[substr(F_Tree1$tip.label, 1, 1)],
     edge.lty=lineType[substr(F_Tree3$tip.label, 1, 1)],
     cex = 0.6, edge.width = 2, type="fan")






######################

## Descripte the Tree
DataTipsandLabel <- data.frame(MyTree$tip.label) # Table with all tips and node number
DataTipsandLabel$tip <- row.names(DataTipsandLabel)
names(DataTipsandLabel) <- c('header_id','tip')
pat_info <- PAT[,c('id', 'born', 'risk', 'sex')]
pat_info$risk <- ifelse(is.na(pat_info$risk), 'other',
                       ifelse(pat_info$risk == 1, 'MSM', 
                              ifelse(pat_info$risk == 2,'HET',
                                     ifelse(pat_info$risk %in% c(3,4), 'IDU','other'))))
names(pat_info) <- c('id', 'born', 'risk', 'sex')
DataTipsandLabel <- merge(DataTipsandLabel, TEST_INFO[,c('id', 'header_id')])
DataTipsandLabel <- merge(DataTipsandLabel, pat_info)
DataTipsandLabel$sex <- ifelse(DataTipsandLabel$sex == 1, 'male', 'female')


print(paste0("In the final tree, there are ", length(MyTree$tip.label), " tips with ",
             length(DataTipsandLabel$id), " SHCS sequences."))

plot(table(DataTipsandLabel$sex, DataTipsandLabel$risk), 
     col = c('violet', 'red', 'darkgrey', 'white'),las = 2,cex = 1.2,
     main = paste0("Sex and Risk Groups, Total tips: ",length(DataTipsandLabel$id)))

#############################################################################
## Study Population
###
dataTipsandLabel <- DataTipsandLabel
studyPop <- data.frame()
studyPop[1,1] <- "Total"
studyPop[1,2] <- length(dataTipsandLabel$id)
studyPop[1,3] <- " "
studyPop[2,1] <- "Male"
studyPop[2,2] <- length(dataTipsandLabel[dataTipsandLabel$sex == 'male','id'])
studyPop[2,3] <- round(100*length(dataTipsandLabel[dataTipsandLabel$sex == 'male','id'])/length(dataTipsandLabel$id),2)
studyPop[3,1] <- "Female"
studyPop[3,2] <- length(dataTipsandLabel[dataTipsandLabel$sex == 'female','id'])
studyPop[3,3] <- round(100*length(dataTipsandLabel[dataTipsandLabel$sex == 'female','id'])/length(dataTipsandLabel$id),2)
studyPop[4,1] <- "MSM"
studyPop[4,2] <- length(dataTipsandLabel[dataTipsandLabel$risk == 'MSM','id'])
studyPop[4,3] <- round(100*length(dataTipsandLabel[dataTipsandLabel$risk == 'MSM','id'])/length(dataTipsandLabel$id),2)
studyPop[5,1] <- "HET"
studyPop[5,2] <- length(dataTipsandLabel[dataTipsandLabel$risk == 'HET','id'])
studyPop[5,3] <- round(100*length(dataTipsandLabel[dataTipsandLabel$risk == 'HET','id'])/length(dataTipsandLabel$id),2)
studyPop[6,1] <- "HET male"
studyPop[6,2] <- length(dataTipsandLabel[dataTipsandLabel$risk == 'HET' & dataTipsandLabel$sex == 'male','id'])
studyPop[6,3] <- round(100*length(dataTipsandLabel[dataTipsandLabel$risk == 'HET'  & dataTipsandLabel$sex == 'male','id'])/length(dataTipsandLabel$id),2)
studyPop[7,1] <- "HET female"
studyPop[7,2] <- length(dataTipsandLabel[dataTipsandLabel$risk == 'HET'  & dataTipsandLabel$sex == 'female','id'])
studyPop[7,3] <- round(100*length(dataTipsandLabel[dataTipsandLabel$risk == 'HET'  & dataTipsandLabel$sex == 'female','id'])/length(dataTipsandLabel$id),2)
studyPop[8,1] <- "Other"
studyPop[8,2] <- length(dataTipsandLabel[dataTipsandLabel$risk == 'other','id'])
studyPop[8,3] <- round(100*length(dataTipsandLabel[dataTipsandLabel$risk == 'other','id'])/length(dataTipsandLabel$id),2)
write.csv(studyPop,"F_StudyPop.csv")



######
x_with_info$all <- 1

incidence_F <- setNames(aggregate(x_with_info$all~x_with_info$diag_year, FUN = sum),c("diag_year","Incidence_F_all"))
incidence_F$prev_F_all <- incidence_F$Incidence_F_all
## to get the prevelance: 
for (i in 2:nrow(incidence_F)) {
  incidence_F[i,3] <- as.numeric(incidence_F[i-1,3]) + as.numeric(incidence_F[i,2])
}


ggplot(data = incidence_F, aes(x = diag_year, y = Incidence_F_all))+
  #geom_line(color = "#00AFBB", size = 2)+
  geom_bar(stat="identity", fill = "#00AFBB", colour="black") +
  scale_x_date(date_labels = "%y", date_breaks = "1 years",
               limits = c(min(incidence_F$diag_year),NA))+
  
  theme_bw()+
  scale_y_continuous( limits=c(0,600)) +
  xlab("\nYear of diagnosis [1985-2018]") +
  ylab("Incidence")+
  theme(axis.text.x = element_text(size=10, angle = 90, hjust=1),
        axis.text.y=element_text(size=14),
        axis.title=element_text(size=15,face="bold"))+
  coord_cartesian(
    xlim = (as.Date(c("1986-01-01","2019-01-01"))))
