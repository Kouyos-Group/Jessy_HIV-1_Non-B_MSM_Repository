##### Figure MSM and HET Cluster Size Analysis
### WITH Bootstraping! 
# clear R's brain
rm(list=ls())



library("scales"); library("ggsci"); library(ggplot2);library(RColorBrewer)

#FastTree
setwd("~")
# MSMBoottable <- read.csv("Downloads/sizeDistribution/MSMcluster_allsubtypes_above50_withboot.csv")
# HETBoottable <- read.csv("Downloads/sizeDistribution/HETcluster_allsubtypes_above50_withboot.csv")
MSMBoottable <- read.csv("Downloads/MSMHET_MTCS_above50/MSMcluster_allsubtypes_above50_withboot.csv")
HETBoottable <- read.csv("Downloads/MSMHET_MTCS_above50/HETcluster_allsubtypes_above50_withboot.csv")


### get the colors: 
display.brewer.pal(n = 12, name = 'Paired')
brewer.pal(n = 12, name = "Paired")


#####
MSMBoottable$sizeMSM <- as.numeric(MSMBoottable$sizeMSM)
HETBoottable$sizeMSM <- as.numeric(HETBoottable$sizeMSM)

MBOOT_dNO <-MSMBoottable[MSMBoottable$distance=="NO",]; HBOOT_dNO <-HETBoottable[HETBoottable$distance=="NO",]
MBOOT_d15 <-MSMBoottable[MSMBoottable$distance=="d15b80",]; HBOOT_d15 <-HETBoottable[HETBoottable$distance=="d15b80",]
MBOOT_d45 <-MSMBoottable[MSMBoottable$distance=="d45b80",]; HBOOT_d45 <-HETBoottable[HETBoottable$distance=="d45b80",]


#MSM in red, HETs in blue
par(mfrow=c(2,3))
barplot(table(MBOOT_dNO$sizeMSM), ylim = range(pretty(c(0, (table(MBOOT_dNO$sizeMSM))))),col ="red",cex.names=1.)
title(xlab="Cluster size",ylab="N of MSM clusters",
      main = "Distance = No Cut off \n Non-B Subtypes",  cex.lab = 1.4, cex.main = 1.8)
barplot(table(MBOOT_d15$sizeMSM), ylim = range(pretty(c(0, (table(MBOOT_d15$sizeMSM))))),col ="red",cex.names=1.)
title(xlab="Cluster size",ylab="N of MSM clusters",
      main = "Distance = 1.5% Bootstrap = 80% \n Non-B Subtypes",  cex.lab = 1.4, cex.main = 1.8)
barplot(table(MBOOT_d45$sizeMSM), ylim = range(pretty(c(0, (table(MBOOT_d45$sizeMSM))))),col ="red",cex.names=1.)
title(xlab="Cluster size",ylab="N of MSM clusters",
      main = "Distance = 4.5% Bootstrap = 80%\n Non-B Subtypes",  cex.lab = 1.4, cex.main = 1.8)
barplot(table(HBOOT_dNO$sizeHET), ylim = range(pretty(c(0, (table(HBOOT_dNO$sizeHET))))),col ="darkblue",cex.names=1.)
title(xlab="Cluster size",ylab="N of MSM clusters",
      main = "Distance = No Cut off \n Non-B Subtypes",  cex.lab = 1.4, cex.main = 1.8)

barplot(table(HBOOT_d15$sizeHET), ylim = range(pretty(c(0, (table(HBOOT_d15$sizeHET))))),col ="darkblue",cex.names=1.)
title(xlab="Cluster size",ylab="N of HET clusters",
      main = "Distance = 1.5% Bootstrap = 80%\n Non-B Subtypes",  cex.lab = 1.4, cex.main = 1.8)

barplot(table(HBOOT_d45$sizeHET), ylim = range(pretty(c(0, (table(HBOOT_d45$sizeHET))))),col ="darkblue",cex.names=1.)
title(xlab="Cluster size",ylab="N of HET clusters",
      main = "Distance = 4.5% Bootstrap = 80% \n Non-B Subtypes",  cex.lab = 1.4, cex.main = 1.8)
par(mfrow=c(1,1))




### changing the visualization to a log log plot 
XBOOT<-table(MBOOT_dNO$sizeMSM) ; XBOOT <- as.data.frame(XBOOT); XBOOT$type <- "MSM_NO" #
YBOOT<-table(HBOOT_dNO$sizeHET) ; YBOOT <- as.data.frame(YBOOT); YBOOT$type <- "HET_NO" #
ABOOT<-table(MBOOT_d15$sizeMSM) ; ABOOT <- as.data.frame(ABOOT); ABOOT$type <- "MSM_d15b80" #
BBOOT<-table(HBOOT_d15$sizeHET) ; BBOOT <- as.data.frame(BBOOT); BBOOT$type <- "HET_d15b80" #
CBOOT<-table(MBOOT_d45$sizeMSM) ; CBOOT <- as.data.frame(CBOOT); CBOOT$type <- "MSM_d45b80" #top
DBOOT<-table(HBOOT_d45$sizeHET) ; DBOOT <- as.data.frame(DBOOT); DBOOT$type <- "HET_d45b80" #top
BOOTscalelogFASTTREE<- rbind(ABOOT,BBOOT,CBOOT,DBOOT,XBOOT,YBOOT)
BOOTscalelogFASTTREE$Var1 <- as.numeric(paste(BOOTscalelogFASTTREE$Var1))
BOOTscalelogFASTTREE <-BOOTscalelogFASTTREE[order(BOOTscalelogFASTTREE$Var1),]




#### 1. Graph FastTree, no cut off, 1.5/ 4.5% GD & bootstrapping 80%
ggplot(BOOTscalelogFASTTREE,  aes(x = Var1, y = Freq, color=type) ) +
  geom_point(size=5) + 
  geom_path(aes(group = type)) + 
  scale_y_log10(breaks = c(1,5,25,50,100,150,200,250))+
  scale_x_log10(breaks = c(2,4,6,8,10,15,20,25,35))+
  theme_bw() + 
  scale_color_manual(values=c("#deebf7" ,"#9ecae1" ,"#3182bd", 
                              "#e5f5e0", "#a1d99b" ,"#31a354") , name = "Cluster Type", labels = c(
                                                                                          "HET 1.5% D 80% BS",
                                                                                          "HET 4.5% D 80% BS",
                                                                                          "HET no threshold",
                                                                                          "MSM 1.5% D 80% BS ",
                                                                                          "MSM 4.5% D 80% BS",
                                                                                          "MSM no threshold"))+
  # scale_color_brewer(palette="Paired",name = "Cluster Type", labels = c("HET Cluster 1.5 Distance", "HET Cluster 4.5 Distance",
  #                                                                       "MSM Cluster 1.5 Distance ", "MSM CLuster 4.5 Distance")) +
  labs(x="Cluster size",  y="Frequency")+
       #title ="FastTree")+
  guides(fill=guide_legend(title = "Cluster", reverse = TRUE, size= 14))+
  theme(
    plot.title = element_text(color="black", size=14, face="bold.italic"),
    axis.title.x = element_text(color="black", size=14, face="bold"),
    axis.title.y = element_text(color="black", size=14, face="bold"),
    #legend.position = c(0.8, 0.6),
    legend.text = element_text(colour="black", size=10, face="bold")) +
  geom_hline(yintercept=1, linetype="dashed", color = "black")



#seperate Graphs for no cut off, GD 1.5%, and GD 4.5%
BOOTscalelogFASTTREE_dNO<- rbind(XBOOT,YBOOT)
BOOTscalelogFASTTREE_dNO$Var1 <- as.numeric(paste(BOOTscalelogFASTTREE_dNO$Var1))
BOOTscalelogFASTTREE_dNO <- BOOTscalelogFASTTREE_dNO[order(BOOTscalelogFASTTREE_dNO$Var1),]


BOOTscalelogFASTTREE_d15<- rbind(ABOOT,BBOOT)
BOOTscalelogFASTTREE_d15$Var1 <- as.numeric(paste(BOOTscalelogFASTTREE_d15$Var1))
BOOTscalelogFASTTREE_d15 <- BOOTscalelogFASTTREE_d15[order(BOOTscalelogFASTTREE_d15$Var1),]


BOOTscalelogFASTTREE_d45<- rbind(CBOOT,DBOOT)
BOOTscalelogFASTTREE_d45$Var1 <- as.numeric(paste(BOOTscalelogFASTTREE_d45$Var1))
BOOTscalelogFASTTREE_d45 <-BOOTscalelogFASTTREE_d45[order(BOOTscalelogFASTTREE_d45$Var1),]


# Graph no cut off
ggplot(BOOTscalelogFASTTREE_dNO,  aes(x = Var1, y = Freq, color=type) ) +
  geom_point(size=5) + 
  geom_path(aes(group = type)) + 
  scale_y_log10(breaks = c(1,5,25,50,100,150,200,250))+
  scale_x_log10(breaks = c(2,4,6,8,10,15,20))+
  theme_bw() + 
  scale_color_manual(values=c("#3182bd" ,"#31a354") , name = "Cluster Type", labels = c("HET no threshold", 
                                                                                       "MSM no threshold "))+
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




# Graph GD 1.5% BS 80& 
ggplot(BOOTscalelogFASTTREE_d15,  aes(x = Var1, y = Freq, color=type) ) +
  geom_point(size=5) + 
  geom_path(aes(group = type)) + 
  scale_y_log10(breaks = c(1,5,25,50,100,150,200,250))+
  scale_x_log10(breaks = c(2,4,6,8,10,15,20))+
  theme_bw() + 
  scale_color_manual(values=c("#deebf7","#e5f5e0") , name = "Cluster Type", labels = c("HET 1.5% D 80% BS", 
                                                                           "MSM 1.5% D 80% BS "))+
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


# Graph GD 4.5% BS 80& 
ggplot(BOOTscalelogFASTTREE_d45,  aes(x = Var1, y = Freq, color=type) ) +
  geom_point(size=5) + 
  geom_path(aes(group = type)) + 
  scale_y_log10(breaks = c(1,5,25,50,100,150,200,250))+
  scale_x_log10(breaks = c(2,4,6,8,10,15,20))+
  theme_bw() + 
  scale_color_manual(values=c("#9ecae1","#a1d99b") , name = "Cluster Type", labels = c("HET 4.5% D 80% BS", 
                                                                           "MSM 4.5% D 80% BS"))+
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



##### Get the Cluster size Table:
MSMBoottable$sizeMSM_factor <- factor(MSMBoottable$sizeMSM) 
D45_boot <- MSMBoottable[MSMBoottable$distance == "d45b80",]

table1(~ 
         sizeMSM_factor
       | subtype , 
       data = D45_boot,
       rowlabelhead = "Sequences", 
       footnote = "<b>NOTE</b> All numbers are no. (%) unless otherwise stated<br/>",
       caption = "<h3><b> Table 1. D45 BS80 </b></h3>")




################################################################################
################################################################################

#### Here place to include HIV Trace!
BOOTscaleFandH <- rbind(BOOTscalelogFASTTREE,BOOTscalelogHIVTRACE)

BOOTscaleFandH$Var1 <- as.numeric(paste(BOOTscaleFandH$Var1))
BOOTscaleFandH <-BOOTscaleFandH[order(BOOTscaleFandH$Var1),]


ggplot(BOOTscaleFandH,  aes(x = Var1, y = Freq, color=type) ) +
  geom_point(size=5) + 
  geom_path(aes(group = type)) + 
  scale_y_log10(breaks = c(1,5,25,50,100,150,200,250))+
  scale_x_log10(breaks = c(2,4,6,8,10,15,20,25,35))+
  # scale_y_log10(breaks = scales::pretty_breaks())+
  # scale_x_log10(breaks = scales::pretty_breaks())+
  theme_bw() +
  scale_color_manual(values=c("#deebf7" ,"#9ecae1" ,"#3182bd", 
                              "#e5f5e0", "#a1d99b" ,"#31a354",
                              "#bcbddc","#fdae6b")
                     ,name = "Cluster Type",labels = c("HET 1.5% D 80% BS",
                                                       "HET 4.5% D 80% BS",
                                                       "HET no threshold",
                                                       "MSM 1.5% D 80% BS ",
                                                       "MSM 4.5% D 80% BS",
                                                       "MSM no threshold",
                                                       "HIV-TRACE - HET","HIV-TRACE - MSM")) +
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



BOOTscaleFd45andH <- rbind(BOOTscalelogFASTTREE_d45,BOOTscalelogHIVTRACE)

ggplot(BOOTscaleFd45andH,  aes(x = Var1, y = Freq, color=type) ) +
  geom_point(size=5) + 
  geom_path(aes(group = type)) + 
  scale_y_log10(breaks = c(1,5,25,50,100,150,200,250))+
  scale_x_log10(breaks = c(2,4,6,8,10,15,20,25,35))+
  # scale_y_log10(breaks = scales::pretty_breaks())+
  # scale_x_log10(breaks = scales::pretty_breaks())+
  theme_bw() +
  scale_color_manual(values=c("#9ecae1" ,
                               "#a1d99b" ,
                              "#bcbddc","#fdae6b")
                     ,name = "Cluster Type",labels = c(
                                                       "HET 4.5% D 80% BS",
                                                       "MSM 4.5% D 80% BS",
                                                       "HIV-TRACE - HET","HIV-TRACE - MSM")) +
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








