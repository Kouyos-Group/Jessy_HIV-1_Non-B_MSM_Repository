##### Figure MSM and HET Cluster Size Analysis
### WITHOUT Bootstrapping! 
# clear R's brain
rm(list=ls())


### load the libraries:
library("scales"); library("ggsci"); library(ggplot2);library(RColorBrewer)

#FastTree
setwd("~")
MSM <- read.csv("Downloads/MSMcluster_allsubtypes_cutoff50.csv")
HET<- read.csv("Downloads/HETcluster_allsubtypes.csv")


### get the colors: 
display.brewer.pal(n = 12, name = 'Paired')
brewer.pal(n = 12, name = "Paired")



#######
MSM$sizeMSM <- as.numeric(MSM$sizeMSM)
HET$sizeHET <- as.numeric(HET$sizeHET)
MSMBoottable$sizeMSM <- as.numeric(MSMBoottable$sizeMSM)
HETBoottable$sizeMSM <- as.numeric(HETBoottable$sizeMSM)

M_d15 <-MSM[MSM$distance=="d15",]; H_d15 <-HET[HET$distance=="d15",]
M_d45 <-MSM[MSM$distance=="d45",]; H_d45 <-HET[HET$distance=="d45",]


# 
# par(mfrow=c(2,2))
# barplot(table(M_d15$sizeMSM), ylim = range(pretty(c(0, (table(M_d15$sizeMSM))))),col ="red",cex.names=1.)
# title(xlab="Cluster size",ylab="N of MSM clusters",
#       main = "Distance = 1.5 % \n Non-B Subtypes",  cex.lab = 1.4, cex.main = 1.8)
# 
# barplot(table(M_d45$sizeMSM), ylim = range(pretty(c(0, (table(M_d45$sizeMSM))))),col ="red",cex.names=1.)
# title(xlab="Cluster size",ylab="N of MSM clusters",
#       main = "Distance = 4.5 % \n Non-B Subtypes",  cex.lab = 1.4, cex.main = 1.8)
# 
# barplot(table(H_d15$sizeHET), ylim = range(pretty(c(0, (table(H_d15$sizeHET))))),col ="darkblue",cex.names=1.)
# title(xlab="Cluster size",ylab="N of HET clusters",
#       main = "Distance = 1.5 % \n Non-B Subtypes",  cex.lab = 1.4, cex.main = 1.8)
# 
# barplot(table(H_d45$sizeHET), ylim = range(pretty(c(0, (table(H_d45$sizeHET))))),col ="darkblue",cex.names=1.)
# title(xlab="Cluster size",ylab="N of HET clusters",
#       main = "Distance = 4.5 % \n Non-B Subtypes",  cex.lab = 1.4, cex.main = 1.8)
# 
# 
# 


### changing the visualization to a log log plot 
A1<-table(M_d15$sizeMSM) ; A1 <- as.data.frame(A1); A1$type <- "MSM_d15"
B1<-table(H_d15$sizeHET) ; B1 <- as.data.frame(B1); B1$type <- "HET_d15"
C1<-table(M_d45$sizeMSM) ; C1 <- as.data.frame(C1); C1$type <- "MSM_d45" #top
D1<-table(H_d45$sizeHET) ; D1 <- as.data.frame(D1); D1$type <- "HET_d45" #top
scalelogFASTTREE<- rbind(A1,B1,C1,D1)
scalelogFASTTREE$Var1 <- as.numeric(paste(scalelogFASTTREE$Var1))
scalelogFASTTREE <-scalelogFASTTREE[order(scalelogFASTTREE$Var1),]




#### 1. Graph FastTree, 1.5% and 4.% GD, no bootstrapping
ggplot(scalelogFASTTREE,  aes(x = Var1, y = Freq, color=type) ) +
  geom_point(size=5) + 
  geom_path(aes(group = type)) + 
  scale_y_log10(breaks = c(1,5,25,50,100,150,200,250))+
  scale_x_log10(breaks = c(2,4,6,8,10,15,20))+
  theme_bw() + 
  scale_color_manual(values=c("darkviolet","blue","orange","red") , name = "Cluster Type", labels = c("HET Cluster 1.5 Distance",
                                                                                          "HET Cluster 4.5 Distance",
                                                                                          "MSM Cluster 1.5 Distance ",
                                                                                          "MSM Cluster 4.5 Distance"))+
  # # scale_color_brewer(palette="Paired",name = "Cluster Type", labels = c("HET Cluster 1.5 Distance", "HET Cluster 4.5 Distance", 
  #                                                                       "MSM Cluster 1.5 Distance ", "MSM CLuster 4.5 Distance")) +
  labs(x="Cluster size",  y="Frequency", title ="FastTree")+
  guides(fill=guide_legend(title = "Cluster", reverse = TRUE, size= 14))+
  theme(
    plot.title = element_text(color="black", size=14, face="bold.italic"),
    axis.title.x = element_text(color="black", size=14, face="bold"),
    axis.title.y = element_text(color="black", size=14, face="bold"),
    legend.position = c(0.8, 0.8),
    legend.text = element_text(colour="black", size=10, face="bold")) +
  geom_hline(yintercept=1, linetype="dashed", color = "black")




#seperate Graphs for GD 1.5% and GD 4.5%
scalelogFASTTREE_d15<- rbind(A1,B1)
scalelogFASTTREE_d15$Var1 <- as.numeric(paste(scalelogFASTTREE_d15$Var1))
scalelogFASTTREE_d15 <-scalelogFASTTREE_d15[order(scalelogFASTTREE_d15$Var1),]


scalelogFASTTREE_d45<- rbind(C1,D1)
scalelogFASTTREE_d45$Var1 <- as.numeric(paste(scalelogFASTTREE_d45$Var1))
scalelogFASTTREE_d45 <-scalelogFASTTREE_d45[order(scalelogFASTTREE_d45$Var1),]




# Graph GD 1.5% 
ggplot(scalelogFASTTREE_d15,  aes(x = Var1, y = Freq, color=type) ) +
  geom_point(size=5) + 
  geom_path(aes(group = type)) + 
  scale_y_log10(breaks = c(1,5,25,50,100,150,200,250))+
  scale_x_log10(breaks = c(2,4,6,8,10,15,20))+
  theme_bw() + 
  scale_color_manual(values=c("#A6CEE3","#B2DF8A") , name = "Cluster Type", labels = c("HET Cluster 1.5 Distance", 
                                                                           "MSM Cluster 1.5 Distance "))+
  # scale_color_brewer(palette="Paired",name = "Cluster Type", labels = c("HET Cluster 1.5 Distance", "HET Cluster 4.5 Distance", 
  #                                                                       "MSM Cluster 1.5 Distance ", "MSM CLuster 4.5 Distance")) +
  labs(x="Cluster size",  y="Frequency", title ="")+
  guides(fill=guide_legend(title = "Cluster", reverse = TRUE, size= 14))+
  theme(
    plot.title = element_text(color="black", size=14, face="bold.italic"),
    axis.title.x = element_text(color="black", size=14, face="bold"),
    axis.title.y = element_text(color="black", size=14, face="bold"),
    legend.position = c(0.8, 0.8),
    legend.text = element_text(colour="black", size=10, face="bold")) +
  geom_hline(yintercept=1, linetype="dashed", color = "black")


# Graph GD 4.5% 
ggplot(scalelogFASTTREE_d45,  aes(x = Var1, y = Freq, color=type) ) +
  geom_point(size=5) + 
  geom_path(aes(group = type)) + 
  scale_y_log10(breaks = c(1,5,25,50,100,150,200,250))+
  scale_x_log10(breaks = c(2,4,6,8,10,15,20))+
  theme_bw() + 
  scale_color_manual(values=c("#1F78B4","#33A02C") , name = "Cluster Type", labels = c("HET Cluster 4.5 Distance", 
                                                                           "MSM Cluster 4.5 Distance"))+
  # scale_color_brewer(palette="Paired",name = "Cluster Type", labels = c("HET Cluster 1.5 Distance", "HET Cluster 4.5 Distance", 
  #                                                                       "MSM Cluster 1.5 Distance ", "MSM CLuster 4.5 Distance")) +
  labs(x="Cluster size",  y="Frequency", title ="")+
  guides(fill=guide_legend(title = "Cluster", reverse = TRUE, size= 14))+
  theme(
    plot.title = element_text(color="black", size=14, face="bold.italic"),
    axis.title.x = element_text(color="black", size=14, face="bold"),
    axis.title.y = element_text(color="black", size=14, face="bold"),
    legend.position = c(0.8, 0.8),
    legend.text = element_text(colour="black", size=10, face="bold")) +
  geom_hline(yintercept=1, linetype="dashed", color = "black")






################################################################################
################################################################################

#RaxML
tableMSM <- read.csv("Downloads/raxml_MSMcluster_allsubtypes_cutoff50.csv")
tableHET <- read.csv("Downloads/raxml_HETcluster_allsubtypes_cutoff50.csv")

tableMSM$sizeMSM <- as.numeric(tableMSM$sizeMSM)
tableHET$sizeHET <- as.numeric(tableHET$sizeHET)


MSM_d15 <-tableMSM[tableMSM$distance=="d15",]; HET_d15 <-tableHET[tableHET$distance=="d15",]
MSM_d45 <-tableMSM[tableMSM$distance=="d45",]; HET_d45 <-tableHET[tableHET$distance=="d45",]



# par(mfrow=c(2,2))
# barplot(table(MSM_d15$sizeMSM), ylim = range(pretty(c(0, (table(MSM_d15$sizeMSM))))),col ="red",cex.names=1.)
# title(xlab="Cluster size",ylab="N of MSM clusters",
#       main = "Distance = 1.5 % \n Non-B Subtypes",  cex.lab = 1.4, cex.main = 1.8)
# 
# barplot(table(MSM_d45$sizeMSM), ylim = range(pretty(c(0, (table(MSM_d45$sizeMSM))))),col ="red",cex.names=1.)
# title(xlab="Cluster size",ylab="N of MSM clusters",
#       main = "Distance = 4.5 % \n Non-B Subtypes",  cex.lab = 1.4, cex.main = 1.8)
# 
# barplot(table(HET_d15$sizeHET), ylim = range(pretty(c(0, (table(HET_d15$sizeHET))))),col ="darkblue",cex.names=1.)
# title(xlab="Cluster size",ylab="N of HET clusters",
#       main = "Distance = 1.5 % \n Non-B Subtypes",  cex.lab = 1.4, cex.main = 1.8)
# 
# barplot(table(HET_d45$sizeHET), ylim = range(pretty(c(0, (table(HET_d45$sizeHET))))),col ="darkblue",cex.names=1.)
# title(xlab="Cluster size",ylab="N of HET clusters",
#       main = "Distance = 4.5 % \n Non-B Subtypes",  cex.lab = 1.4, cex.main = 1.8)
# 



### changing the visualization to a log log plot 

A<-table(MSM_d15$sizeMSM) ; names(A)<- as.numeric(names(A)) ; A <- as.data.frame(A); A$type <- "R_MSM_d15"
B<-table(HET_d15$sizeHET) ; B <- as.data.frame(B); B$type <- "R_HET_d15"
C<-table(MSM_d45$sizeMSM) ; C <- as.data.frame(C); C$type <- "R_MSM_d45"
D<-table(HET_d45$sizeHET) ; D <- as.data.frame(D); D$type <- "R_HET_d45"
scalelog<- rbind(A,B,C,D)
scalelog$Var1 <- as.numeric(paste(scalelog$Var1))
scalelog <-scalelog[order(scalelog$Var1),]



#### 2. Graph, RAxML, 1.5% and 4.% GD, no bootstrapping
ggplot(scalelog,  aes(x = Var1, y = Freq, color=type) ) +
  geom_point(size=5) + 
  geom_path(aes(group = type)) + 
  scale_y_log10(breaks = c(1,5,10,25,50,100,150))+
  scale_x_log10(breaks = scales::pretty_breaks())+
  theme_bw() + 
  scale_color_manual(values=c("darkviolet","blue","orange","red") , name = "Cluster Type", labels = c("HET Cluster 1.5 Distance", "HET Cluster 4.5 Distance", 
                                                                                          "MSM Cluster 1.5 Distance ", "MSM Cluster 4.5 Distance"))+
  # scale_color_brewer(palette="Paired",name = "Cluster:", labels = c("HET Cluster 1.5 Distance", "HET Cluster 4.5 Distance", 
  # "MSM Cluster 1.5 Distance ", "MSM CLuster 4.5 Distance"))  +
  labs(x="Cluster size",  y="Frequency", title ="RAxML")+
  guides(fill=guide_legend(title = "Cluster", reverse = TRUE, size= 14))+
  theme(
    plot.title = element_text(color="black", size=14, face="bold.italic"),
    axis.title.x = element_text(color="black", size=14, face="bold"),
    axis.title.y = element_text(color="black", size=14, face="bold"),
    legend.position = c(0.8, 0.8),
    legend.text = element_text(colour="black", size=10, face="bold"))+
  geom_hline(yintercept=1, linetype="dashed", color = "black")  





################################################################################

#both Approaches in one Graph:
scaleboth <- rbind(scalelog,scalelogFASTTREE)
scaleboth$Var1 <- as.numeric(paste(scaleboth$Var1))
scaleboth <-scaleboth[order(scaleboth$Var1),]



ggplot(scaleboth,  aes(x = Var1, y = Freq, color=type) ) +
  geom_point(size=5) + 
  geom_path(aes(group = type)) + 
  scale_y_log10(breaks = c(1,5,25,50,100,150,200,250))+
  scale_x_log10(breaks = c(2,4,6,8,10,15,20))+
  # scale_y_log10(breaks = scales::pretty_breaks())+
  # scale_x_log10(breaks = scales::pretty_breaks())+
  theme_bw() + 
  scale_color_manual(values=c("#A6CEE3", "#1F78B4", "#CAB2D6" ,"#6A3D9A", "#FB9A99", "#E31A1C", "#FDBF6F" ,"#FF7F00" )
                     ,name = "Cluster Type",labels = c("HET 1.5% D", "HET 4.5% D", 
                                                                     "MSM 1.5% D", "MSM 4.5% D",
                                                                     "RAxML - HET 1.5% D", "RAxML - HET 4.5% D", 
                                                                     "RAxML - MSM 1.5% D", "RAxML - MSM 4.5% D")) +
  # scale_color_brewer(palette="Paired",name = "Cluster Type",labels = c("HET 1.5% D", "HET 4.5% D", 
  #                                                                      "MSM 1.5% D", "MSM 4.5% D",
  #                                                                      "RAxML - HET 1.5% D", "RAxML - HET 4.5% D", 
  #                                                                      "RAxML - MSM 1.5% D", "RAxML - MSM 4.5% D")) +
  labs(x="Cluster size",  y="Frequency")+
  guides(fill=guide_legend(title = "Cluster", reverse = TRUE, size= 14))+
  theme(
    plot.title = element_text(color="black", size=14, face="bold.italic"),
    axis.title.x = element_text(color="black", size=14, face="bold"),
    axis.title.y = element_text(color="black", size=14, face="bold"),
    # legend.position = c(0.8, 0.6),
    legend.text = element_text(colour="black", size=10, face="bold")) +
  geom_hline(yintercept=1, linetype="dashed", color = "black")




################################################################################
################################################################################
