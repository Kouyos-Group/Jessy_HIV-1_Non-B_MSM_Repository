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
setwd("~/Desktop")
data <- read.csv("one_seq_B_recombis_test_info.csv")[,-c(1)]#11'168


CrossTable(data$subtype)
CrossTable(data$rega_subtype)
CrossTable(data$comet_subtype)



#### load the tree (tree generated with FastTree)
setwd("~/Downloads")
MyTree <- read.tree("BRecombinants_withLA_Tree.tre")
 Tree <- read.tree("BRecombinants_Tree.tre")


## Table with all tree tips and node number
dataTipsandLabel <- data.frame(MyTree$tip.label) 
dataTipsandLabel$tip <- row.names(dataTipsandLabel)
names(dataTipsandLabel) <- c('header','tip')



# Merge SHCS id with header ids
TestInfo_need <- TEST_INFO[,c('header_id', 'id','subtype','rega_subtype','comet_subtype')]
names(TestInfo_need) <- c('header', 'id','subtype','rega_subtype','comet_subtype')
dataTipsandLabel <- merge(dataTipsandLabel, TestInfo_need, all.x = T)
dataTipsandLabel$header <- as.character(dataTipsandLabel$header)
dataTipsandLabel$id <- ifelse(is.na(dataTipsandLabel$id), "LA", as.character(dataTipsandLabel$id))


## We make a second tree (subtree), were we only consider the SHCS sequences.
SHCStips <- which(grepl("[[:alpha:]]",MyTree$tip.label)==FALSE)
AlamosTips <- which(grepl("[[:alpha:]]",MyTree$tip.label))
#LAtips <-  as.numeric(dataTipsandLabel[grepl("[[:alpha:]]", dataTipsandLabel$header),'tip'])
SHCStree <- drop.tip(MyTree,c(AlamosTips))



## Table with all tips and node number for SHCS tree
dataTipsandLabelSHCS <- data.frame(SHCStree$tip.label)
dataTipsandLabelSHCS$tip <- row.names(dataTipsandLabelSHCS)
names(dataTipsandLabelSHCS) <- c('header','tipSHCS')
dataTipsandLabelSHCS<- merge(dataTipsandLabelSHCS, TestInfo_need, all.x = T)
dataTipsandLabelSHCS$id <- ifelse(is.na(dataTipsandLabelSHCS$id), "LA", as.character(dataTipsandLabelSHCS$id))

## preparing dataframe: header id must be in first row.names 
dd <- dataTipsandLabelSHCS
dd <- na.omit(dd)
row.names(dd) <- NULL
rownames(dd) <- dd[,1]


#################Preparing color vectors
#For tips
display.brewer.pal(3, "Set2")
ncol<-c(brewer.pal(3, "Set2"),"white")
names(ncol) = c("B", "RECOMBINANT OF B, D","RECOMBINANT OF B, F1","Unknown")
# levels(unique(data$rega_subtype))




#ggtree:
##plotting/adjusting the phylogenetic tree
plotTree <-ggtree(SHCStree) + coord_flip()
plotTree
plotTree +geom_text2(aes(subset=!isTip, label=node), hjust=-.3) #+ geom_tiplab(size=1)
plotTree %<+% dd + geom_tiplab()




##Adding tips and cluster labels

p <- plotTree %<+% dd +
  geom_tree(aes(color=subtype), size=1.5) +
  # geom_tiplab(aes(label = subtype)) +
  scale_color_manual(breaks =  c("B","Recombinant", "Unknown") , values = c("grey","#66C2A5","red"),
                     guide=guide_legend(title="Subtype",keywidth =0.5, keyheight = 0.5, order =1, 
                                        override.aes = list(size=5)))+
  theme(legend.position = c(0.83,0.9),
        legend.background = element_rect(fill=NA),
        legend.text = element_text(size=18.5),
        legend.title = element_text(size = 30.5),
        legend.spacing.y = unit(0.02,"cm"));p


prega <- plotTree %<+% dd +
  geom_tree(aes(color=rega_subtype), size=1.5) +
  # geom_tiplab(aes(label = subtype)) +
  scale_color_manual(breaks =  c("B","Recombinant of B, D","Recombinant of B, F1") , values = c("grey","#66C2A5","#8DA0CB"));prega

pcomet <- plotTree %<+% dd +
  geom_tree(aes(color=comet_subtype), size=1.5) +
  # geom_tiplab(aes(label = subtype)) +
  scale_color_manual(breaks =  c(unique(dd$comet_subtype)) , values = c("grey",brewer.pal(6, "Set2")));pcomet


headers <- dd[which(dd$subtype=="Recombinant"),"header"] #3019
index <- which(MyTree$tip.label %in% headers)


subtree<- treeio::tree_subset(MyTree, node=12079, levels_back=3)
sub_tree <-ggtree(subtree) + coord_flip() #+ geom_tiplab() #+ geom_nodelab() 
sub_tree+geom_text2(aes(subset=!isTip, label=node), hjust=-.3) #+ geom_tiplab(size=1)

psub <- sub_tree %<+% dd +
  geom_tree(aes(color=subtype), size=1) +
  # geom_tiplab(aes(label = subtype)) +
  scale_color_manual(breaks =  c("B","Recombinant", "Unknown") , values = c("grey","#66C2A5","red"),
                     guide=guide_legend(title="Subtype",keywidth =0.5, keyheight = 0.5, order =1, 
                                        override.aes = list(size=5)))+
  theme(legend.position = c(0.83,0.9),
        legend.background = element_rect(fill=NA),
        legend.text = element_text(size=18.5),
        legend.title = element_text(size = 30.5),
        legend.spacing.y = unit(0.02,"cm"));psub

subsubtree<- treeio::tree_subset(MyTree, node=12076, levels_back=7)
subsub_tree <-ggtree(subsubtree) + coord_flip() #+ geom_tiplab() #+ geom_nodelab() 
psubsub <- subsub_tree %<+% dd +
  geom_tree(aes(color=subtype), size=1.0) +
  # geom_tiplab(aes(label = subtype)) +
  scale_color_manual(breaks =  c("B","Recombinant", "Unknown") , values = c("grey","#66C2A5","red"),
                     guide=guide_legend(title="Subtype",keywidth =0.5, keyheight = 0.5, order =1, 
                                        override.aes = list(size=5)))+
  theme(legend.position = c(0.93,0.9),
        legend.background = element_rect(fill=NA),
        legend.text = element_text(size=18.5),
        legend.title = element_text(size = 30.5),
        legend.spacing.y = unit(0.02,"cm"));psubsub


pregarega <- subsub_tree %<+% dd +
  geom_tree(aes(color=rega_subtype), size=1.0) +
  # geom_tiplab(aes(label = subtype)) +
  scale_color_manual(breaks =  c("B","Recombinant of B, D","Recombinant of B, F1") , values = c("grey","purple","orange"),
                     guide=guide_legend(title="REGA Subtype",keywidth =0.5, keyheight = 0.5, order =1, 
                                        override.aes = list(size=5)))+
  theme(legend.position = c(0.87,0.9),
        legend.background = element_rect(fill=NA),
        legend.text = element_text(size=18.5),
        legend.title = element_text(size = 30.5),
        legend.spacing.y = unit(0.02,"cm"));pregarega


pcometcomet <- subsub_tree %<+% dd +
  geom_tree(aes(color=comet_subtype), size=1.0) +
  # geom_tiplab(aes(label = subtype)) +
  scale_color_manual(breaks =  c(unique(dd$comet_subtype)) , values = c("grey",brewer.pal(6, "Dark2")),
                     guide=guide_legend(title="COMET Subtype",keywidth =0.5, keyheight = 0.5, order =1, 
                                        override.aes = list(size=5)))+
  theme(legend.position = c(0.87,0.9),
        legend.background = element_rect(fill=NA),
        legend.text = element_text(size=18.5),
        legend.title = element_text(size = 30.5),
        legend.spacing.y = unit(0.02,"cm"));pcometcomet

