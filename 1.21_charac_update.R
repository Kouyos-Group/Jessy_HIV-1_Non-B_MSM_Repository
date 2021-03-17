# DEZ 2020, JDR

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
data <- read.csv("Output/fulltable_study_population_update.csv")
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
data$diag_year <- as.Date(as.character(data$diag_year), tryFormats = c("%Y"))
data$diag_year  <- substr(data$diag_date,1,4)
data$diag_year<- as.Date(as.character(data$diag_year), tryFormats = c("%Y"))


###SUBTYPE
ta<-table(data$subtype)
cross_table <-CrossTable(data$subtype)


mosaicplot(t(ta),col=TRUE, main="HIV Subtypes distribution within the SHCS")

## Study Population
###
Others <- data.frame()
Others[1,1] <- "Sequence"
Others[1,2] <- "03_AB"
Others[1,3] <- "05_DF"
Others[1,4] <- "06_CPX"
Others[1,5] <- "07_BC"
Others[1,6] <- "09_CPX"
Others[1,7] <- "10_CD"
Others[1,8] <- "11_CPX"
Others[1,9] <- "12_BF"
Others[1,10] <- "13_CPX"
Others[1,11] <- "14_BG"
Others[1,12] <- "18_CPX"
Others[1,13] <- "19_CPX"
Others[1,14] <- "20_BG"
Others[1,15] <- "24_BG"
Others[1,16] <- "25_CPX"
Others[1,17] <- "27_CPX"
Others[1,18] <- "29_BF"
Others[1,19] <- "31_BC"
Others[1,20] <- "35_AD"
Others[1,21] <- "37_CPX"
Others[1,22] <- "40_BF"
Others[1,23] <- "43_02G"
Others[1,24] <- "D"
Others[1,25] <- "H"
Others[1,26] <- "J"
Others[1,27] <- "K"
Others[1,28] <- "Recombinant"



Others[2,1] <- "N"
Others[2,2]<- cross_table$t[3]
Others[2,3]<- cross_table$t[4]
Others[2,4]<- cross_table$t[5]
Others[2,5]<- cross_table$t[6]
Others[2,6]<- cross_table$t[7]
Others[2,7]<- cross_table$t[8]
Others[2,8]<- cross_table$t[9]
Others[2,9]<- cross_table$t[10]
Others[2,10] <- cross_table$t[11]
Others[2,11]<- cross_table$t[12]
Others[2,12]<- cross_table$t[13]
Others[2,13]<- cross_table$t[14]
Others[2,14]<- cross_table$t[15]
Others[2,15]<- cross_table$t[16]
Others[2,16]<- cross_table$t[17]
Others[2,17]<- cross_table$t[18]
Others[2,18]<- cross_table$t[19]
Others[2,19]<- cross_table$t[20]
Others[2,20] <- cross_table$t[21]
Others[2,21]<- cross_table$t[22]
Others[2,22]<- cross_table$t[23]
Others[2,23]<- cross_table$t[24]
Others[2,24]<- cross_table$t[28]
Others[2,25]<- cross_table$t[31]
Others[2,26]<- cross_table$t[32]
Others[2,27]<- cross_table$t[33]
Others[2,28]<- cross_table$t[34]



Others[3,1] <- "Proportion"
Others[3,2]<- round(cross_table$prop.row[3],3)*100
Others[3,3]<- round(cross_table$prop.row[4],3)*100
Others[3,4]<- round(cross_table$prop.row[5],3)*100
Others[3,5]<- round(cross_table$prop.row[6],3)*100
Others[3,6]<- round(cross_table$prop.row[7],3)*100
Others[3,7]<- round(cross_table$prop.row[8],3)*100
Others[3,8]<- round(cross_table$prop.row[9],3)*100
Others[3,9]<- round(cross_table$prop.row[10],3)*100
Others[3,10] <- round(cross_table$prop.row[11],3)*100
Others[3,11]<- round(cross_table$prop.row[12],3)*100
Others[3,12]<- round(cross_table$prop.row[13],3)*100
Others[3,13]<- round(cross_table$prop.row[14],3)*100
Others[3,14]<- round(cross_table$prop.row[15],3)*100
Others[3,15]<- round(cross_table$prop.row[16],3)*100
Others[3,16]<- round(cross_table$prop.row[17],3)*100
Others[3,17]<- round(cross_table$prop.row[18],3)*100
Others[3,18]<- round(cross_table$prop.row[19],3)*100
Others[3,19]<- round(cross_table$prop.row[20],3)*100
Others[3,20] <- round(cross_table$prop.row[21],3)*100
Others[3,21]<- round(cross_table$prop.row[22],3)*100
Others[3,22]<- round(cross_table$prop.row[23],3)*100
Others[3,23]<- round(cross_table$prop.row[24],3)*100
Others[3,24]<- round(cross_table$prop.row[28],3)*100
Others[3,25]<- round(cross_table$prop.row[32],3)*100
Others[3,26]<- round(cross_table$prop.row[32],3)*100
Others[3,27]<- round(cross_table$prop.row[33],3)*100
Others[3,28]<- round(cross_table$prop.row[34],3)*100

Others <- Others[-c(1),]
colnames(Others) <- c("Sequences","03_AB","05_DF","06_CPX","07_BC","09_CPX",
                      "10_CD","11_CPX","12_BF","13_CPX",
                      "14_BG","18_CPX","19_CPX","20_BG",
                      "24_BG","25_CPX","27_CPX","29_BF",
                      "31_BC","35_AD","37_CPX","40_BF",
                      "43_02G","D","H","J","K","Recombinant")

library(grid)
library(gridExtra)
library(gtable)
grid.newpage()
grid.draw(tableGrob(format(Others, big.mark=","), rows = NULL))
tt3 <- ttheme_minimal(
  core=list(bg_params = list(fill = "white", col="lightgray"),
            fg_params=list(fontface=1)))#,
#colhead=list(fg_params=list(col="black", fontface=1)))#,
#rowhead=list(fg_params=list(col="black", fontface=3L)))
g <- tableGrob(Others, theme = tt3, rows = NULL)
grid.draw(g)
g <- gtable_add_grob(g,
                     grobs = rectGrob(gp = gpar(fill = NA, lwd = 2)),
                     t = 2, b = nrow(g), l = 1, r = ncol(g))
g <- gtable_add_grob(g,
                     grobs = rectGrob(gp = gpar(fill = NA, lwd = 2)),
                     t = 1, l = 1, r = ncol(g))
grid.draw(g)
h = grid::convertHeight(sum(g$heights), "in", TRUE)
w = grid::convertWidth(sum(g$widths), "in", TRUE)
# ggplot2::ggsave("raxml_table_cluster_distribution.eps", g, width=w, height=h)

recombinants <- data[data$subtype == "RECOMBINANT",]
a<-table(droplevels(recombinants$rega_subtype))
a <- data.frame(a)
b <-table(droplevels(recombinants$comet_subtype))
b <- data.frame(b)

grid.newpage()
grid.draw(tableGrob(format(a, big.mark=","), rows = NULL)) #rega
grid.newpage()
grid.draw(tableGrob(format(b, big.mark=","), rows = NULL)) #comet
###################




ta2<-table(data$copy_subtype)
CrossTable(data$copy_subtype)
#B: 9512 (74.0%), non-B: 3340 (26.0%)
mosaicplot(t(ta2),col=TRUE, main="HIV Subtypes distribution within the SHCS")



n_Bs <-data.frame(matrix(NA, ncol =3, nrow = 2))
colnames(n_Bs ) <- c("Subtype","N","Percentage")
n_Bs$Subtype <- c("B", "non-B")
n_Bs$N <- c( t(ta2)[1], t(ta2)[2])
n_Bs$Percentage <- c( t(ta2)[1]/(sum(t(ta2))),  t(ta2)[2]/(sum(t(ta2))))
n_Bs$Percentage <- round(n_Bs$Percentage,3)

ggplot(n_Bs, aes(x = Subtype, y = Percentage, fill=Subtype)) +
  geom_bar( stat = "identity", position=position_dodge(), colour="black" )+
  scale_fill_manual(values=c("#999999", "lightblue"))+
  guides(fill=FALSE)+
  geom_text(aes(label = paste0(Percentage*100, "%")), vjust = c(+1.5,+1.5), size=10)+
  theme_bw()+
  scale_y_continuous(labels = scales::percent, limits=c(0,1)) +
  #ggtitle("Subtype distribution in the SHCS DB") +
  xlab("Subtype") +
  ylab("Frequency [%]")+
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=15,face="bold"))


ggplot(n_Bs, aes(x = Subtype, y = N, fill=Subtype)) +
  geom_bar( stat = "identity", position=position_dodge(), colour="black" )+
  scale_fill_manual(values=c("#999999", "lightblue"))+
  guides(fill=FALSE)+
  geom_text(aes(label = paste0(N)), vjust = c(+1.5,+1.5), size=10)+
  theme_bw()+
  scale_y_continuous(limits=c(0,10000)) +
  #ggtitle("Subtype distribution in the SHCS DB") +
  xlab("Subtype") +
  ylab("Number of Sequences")+
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=15,face="bold"))


ggplot(n_Bs, aes(x = Subtype, y = N, fill=Subtype)) +
  geom_bar( stat = "identity", position=position_dodge(), colour="black" )+
  scale_fill_manual(values=c("#999999", "lightblue"))+
  guides(fill=FALSE)+
  geom_text(aes(label = paste0(N)), vjust = c(+1.2,+1.2), size=10)+
  geom_text(aes(label = paste0("(",Percentage*100, "%)")), vjust = c(+3.0,+3.0), size=7)+
  theme_bw()+
  scale_y_continuous(limits=c(0,10000)) +
  #ggtitle("Subtype distribution in the SHCS DB") +
  xlab("Subtype") +
  ylab("Number of Sequences")+
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=15,face="bold"))


### LAB
ta4<-table(data$lab)
CrossTable(data$lab)



n_labs <-data.frame(matrix(NA, ncol = 3, nrow = 4))
colnames(n_labs ) <- c("Lab","N","Percentage")
n_labs$Lab <- c("Basel",  "Geneva","Lausanne","Zurich")
n_labs$N <- c(t(ta4)[1],t(ta4)[2], t(ta4)[3], t(ta4)[4])
total <- sum(t(ta4))
n_labs$Percentage <- c(t(ta4)[1]/total, t(ta4)[2]/total, 
                       t(ta4)[3]/total, t(ta4)[4]/total)
n_labs$Percentage <- round(n_labs$Percentage,3)



ggplot(n_labs, aes(x = reorder(Lab,-Percentage), y = Percentage, fill=Lab)) +
  geom_bar( stat = "identity", position=position_dodge(), colour="black" )+
  scale_fill_manual(values=c(  "#999999", "#999999", "#999999","#999999" ))+
  guides(fill=FALSE)+
  geom_text(aes(label = paste0(Percentage*100, "%")), vjust = +1.2, size=7)+
  theme_bw()+
  scale_y_continuous(labels = scales::percent, limits=c(0,1)) +
  #ggtitle("Lab distribution in the SHCS DB") +
  xlab("Lab") +
  ylab("Frequency [%]")+
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=15,face="bold"))


ggplot(n_labs, aes(x = reorder(Lab,-N), y = N, fill=Lab)) +
  geom_bar( stat = "identity", position=position_dodge(), colour="black" )+
  scale_fill_manual(values=c(  "#999999", "#999999", "#999999","#999999" ))+
  guides(fill=FALSE)+
  geom_text(aes(label = paste0(N)), vjust = +1.2, size=10)+
  geom_text(aes(label = paste0("(",Percentage*100, "%)")), vjust = +3.2, size=6, colour ="blue")+
  theme_bw()+
  scale_y_continuous( limits=c(0,6000),
                      sec.axis = sec_axis(~./6000 , name = "Frequency [%]",
                                          labels = scales::percent)) +
  #ggtitle("Lab distribution in the SHCS DB") +
  xlab("Lab") +
  ylab("Number of Sequences")+
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=15,face="bold"),
        axis.title.y.right = element_text(colour = "blue"))





##Pol_len
summary(data$pol_len)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 858     887    1302    1180    1302    1965 

boxplot(data$pol_len,
        main = "Mean partial pol sequence length",
        xlab = "Length [nt]",
        ylab = "Pol",
        col = "orange",
        border = "brown",
        horizontal = TRUE,
        notch = TRUE)

boxplot(pol_len~lab,
        data= data,
        main = "Mean partial pol sequence length by lab",
        xlab = "Length [nt]",
        ylab = "Lab",
        cex.axis=1.2, cex.lab=1.5, cex.main=2,
        face="bold",
        col = "orange",
        border = "brown",
        horizontal = TRUE,
        notch = F)

ggplot(data, aes(x=lab, y=pol_len)) + 
  geom_boxplot(outlier.colour="black", outlier.shape=16,
               outlier.size=2, notch=FALSE, col="orange")+
  labs(title="Mean partial pol sequence length by lab",y="Length [nt]", x = "Lab")+
  coord_flip()+
  theme_bw() + 
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=15,face="bold"),
        plot.title = element_text(size=20,face="bold.italic"))


ggplot(data, aes(x=lab, y=pol_len)) + 
  geom_boxplot(outlier.colour="black", outlier.shape=16,
               outlier.size=2, notch=FALSE, col="orange")+
  #labs(title="Mean partial pol sequence length by lab",y="Length [nt]", x = "Lab")+
  labs(y="Length [nt]", x = "Lab")+
  theme_bw() + 
  theme(
    plot.title = element_text(color="black", size=14, face="bold.italic"),
    axis.title.x = element_text(color="black", size=14, face="bold"),
    axis.title.y = element_text(color="black", size=14, face="bold"),
    #legend.position = c(0.9, 0.8),
    legend.text = element_text(colour="black", size=14, face="plain"),
    legend.title = element_text(colour="black", size=14, face="bold")) +
  coord_flip()+
  geom_hline(yintercept=887, linetype="dashed", color = "blue")


aggregate(data$pol_len, list(data$lab), mean)
# Group.1        x
# 1    Basel 1286.6730
# 2   Geneva  889.8317
# 3 Lausanne 1432.9896
# 4   Zurich  1302.2978




### GENDER:
ta3<-table(data$sex) 
CrossTable(data$sex)
#1 = male: 9,276 (72,2%), 2= female: 3,576 (27.8%)

n_sex <-data.frame(matrix(NA, ncol = 3, nrow = 2))
colnames(n_sex ) <- c("Sex","N","Percentage")
n_sex$Sex <- c("Male",  "Female")
n_sex$N <- c(t(ta3)[1],t(ta3)[2])
n_sex$Percentage <- c(t(ta3)[1]/total, t(ta3)[2]/total)
n_sex$Percentage <- round(n_sex$Percentage,3)


mosaicplot(t(ta3),col=TRUE,set_labels=list(sex = c("M", "F")),
           main="Gender distribution within the SHCS")

data$copysex <- factor(data$sex, levels=c(1,2),
                       labels=c("Male", 
                                "Female"))


barplot(table(data$copysex), ylim=c(0,10000), 
        col    = c("gray90", "gray60"),
        xlab   = "Gender",
        ylab   = "Count")

pie(table(data$copysex), pch = 19, cex = 1.5,
    init.angle = 90,
    col = c("lightblue", "gray60"))


ggplot(n_sex, aes(x = reorder(Sex,-Percentage), y = Percentage, fill=Sex)) +
  geom_bar( stat = "identity", position=position_dodge(), colour="black" )+
  scale_fill_manual(values=c(  "#999999", "lightblue" ))+
  guides(fill=FALSE)+
  geom_text(aes(label = paste0(Percentage*100, "%")), vjust = +1.2, size=10)+
  theme_bw()+
  scale_y_continuous(labels = scales::percent, limits=c(0,1)) +
  #ggtitle("Sex distribution in the SHCS DB") +
  xlab("Sex") +
  ylab("Frequency [%]")+
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=15,face="bold"))


ggplot(n_sex, aes(x = reorder(Sex,-N), y = N, fill=Sex)) +
  geom_bar( stat = "identity", position=position_dodge(), colour="black" )+
  scale_fill_manual(values=c(  "#999999", "lightblue" ))+
  guides(fill=FALSE)+
  geom_text(aes(label = paste0(N)), vjust = +1.2, size=10)+
  theme_bw()+
  scale_y_continuous( limits=c(0,10000)) +
  #ggtitle("Sex distribution in the SHCS DB") +
  xlab("Sex") +
  ylab("Count")+
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=15,face="bold"))


ggplot(n_sex, aes(x = reorder(Sex,-N), y = N, fill=Sex)) +
  geom_bar( stat = "identity", position=position_dodge(), colour="black" )+
  scale_fill_manual(values=c(  "#999999", "lightblue" ))+
  guides(fill=FALSE)+
  geom_text(aes(label = paste0(N)), vjust = +1.2, size=10)+
  geom_text(aes(label = paste0("(",Percentage*100, "%)")), vjust = +3.5, size=6, colour = "blue")+
  theme_bw()+
  scale_y_continuous( limits=c(0,10000),
                      sec.axis = sec_axis(~./10000 , name = "Frequency [%]",
                                          labels = scales::percent)) +
  #ggtitle("Sex distribution in the SHCS DB") +
  xlab("Sex") +
  ylab("Count")+
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=15,face="bold"),
        axis.title.y.right = element_text(colour="blue"))



### Age
summary(data$age)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 18.00   44.00   52.00   52.08   59.00  101.00 


#check for normality:
hist(data$age, xlab="Age (years)", ylab="Freuency", main= "Histogram of Age")

ggplot(data, aes(age)) +
  geom_histogram(aes(x=age, y=..density..), bins=50, fill="#d3d3d3", color="black") +
  stat_function(aes(colour = "Normal"),fun=dnorm, args = list(mean=mean(data$age), sd=sd(data$age)), size=1.5) +
  geom_density(aes(color = "Density")) +
  scale_colour_manual("", values = c("red", "blue")) +
  labs(title="Histogram of Age", x="Age (years)",y="Density") +
  theme_bw() +
  coord_cartesian(xlim=c(15,100), ylim=c(0,0.06))

ggpubr::ggqqplot(data$age)
shapiro.test(data$age[1:5000]) #W = 0.99197, p-value = 3.243e-16
#normal distributed so able to take the mean




### Age at Sample date 
summary(data$age_sampledate)

# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#    4.00   31.00   37.00   38.67   45.00   82.00
# 


#check for normality:
hist(data$age_sampledate, xlab="Age (years)", ylab="Freuency", main= "Histogram of Age at Sample Date")


ggplot(data, aes(age_sampledate)) +
  geom_histogram(aes(x=age_sampledate, y=..density..), bins=50, fill="#d3d3d3", color="black") +
  stat_function(aes(colour = "Normal"),fun=dnorm, args = list(mean=mean(data$age_sampledate), sd=sd(data$age_sampledate)), size=1.5) +
  geom_density(aes(color = "Density")) +
  scale_colour_manual("", values = c("red", "blue")) +
  labs(title="Histogram of Age at Sample Date", x="Age (years)",y="Density") +
  theme_bw() +
  coord_cartesian(xlim=c(15,100), ylim=c(0,0.06))


ggpubr::ggqqplot(data$age_sampledate)
shapiro.test(data$age_sampledate[1:5000]) #W = 0.9617, p-value < 2.2e-16
#also (more or less) normal distributed, can take the mean as well



### Age at Diagnosis date 
summary(data$age_at_diag)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.00   26.00   32.00   34.54   41.00   81.00 



#check for normality:
hist(data$age_at_diag, xlab="Age (years)", ylab="Freuency", main= "Histogram of Age at Diagnosis Date")


ggplot(data, aes(age_at_diag)) +
  geom_histogram(aes(x=age_at_diag, y=..density..), bins=50, fill="#d3d3d3", color="black") +
  stat_function(aes(colour = "Normal"),fun=dnorm, args = list(mean=mean(data$age_at_diag), sd=sd(data$age_at_diag)), size=1.5) +
  geom_density(aes(colour = "Density")) +
  scale_colour_manual("", values = c("red", "blue")) +
  ggtitle("Histogram of Age at Diagnosis Date") +
  xlab("Age at diagnosis date") +
  ylab("Density") +
  theme_bw() +
  coord_cartesian(xlim=c(15,100), ylim=c(0,0.06))


ggpubr::ggqqplot(data$age_at_diag)
shapiro.test(data$age_at_diag[1:5000]) #W = 0.93837, p-value < 2.2e-16
#also (more of less) normal distributed, can take the mean as well



####### Pearson correlation test:
ggpubr::ggscatter(data, x = "age_at_diag", y = "age_sampledate", 
                  add = "reg.line", conf.int = TRUE, 
                  cor.coef = TRUE, cor.method = "pearson",
                  xlab = "Age at Diagnosis date", ylab = "Age at Sample date")

res <- cor.test(data$age_at_diag, data$age_sampledate, 
                method = "pearson")
res

res$p.value # Extract the p.value
res$estimate # Extract the correlation coefficient

#We can conclude that age and age_sampledate are significantly correlated with a 
#correlation coefficient of 0.86 and p-value of  < 2.2e-16.




## ethnicity
data$copyethnicity <- factor(data$ethnicity, levels=c(1,2,3,4,0,9),
                             labels=c("White", 
                                      "Black",
                                      "Hispano-American",
                                      "Asian",
                                      "Others",
                                      "Unknown"))
ta5 <- table(data$copyethnicity)
CrossTable(data$copyethnicity)
# White            Black Hispano-American            Asian 
# 9874            1603             400            429
# Others          Unknown 
# 27              515



n_ethnicity <-data.frame(matrix(NA, ncol = 3, nrow = 6))
colnames(n_ethnicity ) <- c("Ethnicity","N","Percentage")
n_ethnicity$Ethnicity <- c("White", 
                           "Black",
                           "Hispano-American",
                           "Asian",
                           "Others",
                           "Unknown")

n_ethnicity$N <- c(t(ta5)[1],t(ta5)[2],t(ta5)[3],t(ta5)[4],t(ta5)[5], t(ta5)[6])
n_ethnicity$Percentage <- c(t(ta5)[1]/total, t(ta5)[2]/total,
                            t(ta5)[3]/total, t(ta5)[4]/total,
                            t(ta5)[5]/total, t(ta5)[6]/total)
n_ethnicity$Percentage <- round(n_ethnicity$Percentage,3)


ggplot(n_ethnicity, aes(x = reorder(Ethnicity,-N), y = N, fill=Ethnicity)) +
  geom_bar( stat = "identity", position=position_dodge(), colour="black" )+
  scale_fill_manual(values=c(  rep("#999999",6 )))+
  guides(fill=FALSE)+
  geom_text(aes(label = paste0(N)), vjust = c(+1.2,+1.2,-1.5,-1.5,-1.5,-1.5), size=7)+
  geom_text(aes(label = paste0("(",Percentage*100, "%)")), vjust = c(+3.0,+2.8,-0.5,-0.5,-0.5,-0.5), size=5, colour= "blue")+
  theme_bw()+
  scale_y_continuous( limits=c(0,10000),
                      sec.axis = sec_axis(~./10000 , name = "Frequency [%]",
                                          labels = scales::percent)) +
  #ggtitle("Ethnicity distribution in the SHCS DB") +
  xlab("Ethnicity") +
  ylab("Count")+
  theme(axis.text=element_text(size=14),
        axis.text.x=element_text(angle=45, hjust=1),
        axis.title=element_text(size=15,face="bold"),
        axis.title.y.right =  element_text(colour="blue"))


ggplot(n_ethnicity, aes(x = reorder(Ethnicity,-Percentage), y = Percentage, fill=Ethnicity)) +
  geom_bar( stat = "identity", position=position_dodge(), colour="black" )+
  scale_fill_manual(values=c(  rep("#999999",6 )))+
  guides(fill=FALSE)+
  geom_text(aes(label = paste0(Percentage*100, "%")), vjust = c(+1.2,+1.2,-0.5,-1.0,-0.5,-0.5), size=7)+
  theme_bw()+
  scale_y_continuous(labels = scales::percent, limits=c(0,1)) +
  #ggtitle("Ethnicity distribution in the SHCS DB") +
  xlab("Ethnicity") +
  ylab("Frequency [%]")+
  theme(axis.text=element_text(size=14),
        axis.text.x=element_text(angle=45, hjust=1),
        axis.title=element_text(size=15,face="bold"))



## region 
data$copyregion <- factor(data$region, levels=c(155,154,39,
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


ta6 <- table(data$copyregion)
# Western, Northern and Southern Europe 
# 10157
# Western and Central Africa 
# 856
# Eastern and Southern Africa 
# 540
# Latin America 
# 531
# Eastern and Southern Asia 
# 361 
# Western Asia and North Africa 
# 168 
# Eastern Europe 
# 151 
# North America and Oceania 
# 79 
# Unknown 
# 9 


n_region <-data.frame(matrix(NA, ncol = 3, nrow = 9))
colnames(n_region ) <- c("Region","N","Percentage")
n_region$Region <- c("Western, Northern and Southern Europe",
                     "Western and Central Africa",
                     "Eastern and Southern Africa",
                     "Latin America",
                     "Eastern and Southern Asia",
                     "Western Asia and North Africa",
                     "Eastern Europe",
                     "North America and Oceania",
                     "Unknown")

n_region$N <- c(t(ta6)[1],t(ta6)[2],t(ta6)[3],t(ta6)[4],t(ta6)[5], t(ta6)[6],
                t(ta6)[7],t(ta6)[8],t(ta6)[9])
n_region$Percentage <- c(t(ta6)[1]/total, t(ta6)[2]/total,
                         t(ta6)[3]/total, t(ta6)[4]/total,
                         t(ta6)[5]/total, t(ta6)[6]/total,
                         t(ta6)[7]/total,
                         t(ta6)[8]/total, t(ta6)[9]/total)
n_region$Percentage <- round(n_region$Percentage,3)


ggplot(n_region, aes(x = reorder(Region,-N), y = N, fill=Region)) +
  geom_bar( stat = "identity", position=position_dodge(), colour="black" )+
  scale_fill_manual(values=c(  rep("#999999",9 )))+
  guides(fill=FALSE)+
  geom_text(aes(label = paste0(N)), vjust = c(+1.2, rep(-1.2,8)), size=9)+
  geom_text(aes(label = paste0("(",Percentage*100, "%)")), vjust = c(+3.5, rep(-0.5,8)), size=5, colour="blue")+
  theme_bw()+
  scale_y_continuous( limits=c(0,10000),
                      sec.axis = sec_axis(~./10000 , name = "Frequency [%]",
                                          labels = scales::percent)) +
  #ggtitle("Region distribution in the SHCS DB") +
  xlab("Region") +
  ylab("Count")+
  theme(axis.text.x = element_text(size=10, angle = 45, hjust = 1),
        axis.text.y=element_text(size=14),
        axis.title=element_text(size=15,face="bold"),
        axis.title.y.right =element_text(colour = "blue"))


ggplot(n_region, aes(x =reorder(Region,-Percentage), y = Percentage, fill=Region)) +
  geom_bar( stat = "identity", position=position_dodge(), colour="black" )+
  scale_fill_manual(values=c(  rep("#999999",9 )))+
  guides(fill=FALSE)+
  geom_text(aes(label = paste0(Percentage*100, "%")), vjust = c(+1.2, rep(-0.5,8)), size=6)+
  theme_bw()+
  scale_y_continuous(labels = scales::percent, limits=c(0,1)) +
  #ggtitle("Region distribution in the SHCS DB") +
  xlab("Region") +
  ylab("Frequency [%]")+
  theme(axis.text.x = element_text(size=9.5, angle = 45, hjust = 1),
        axis.text.y=element_text(size=14),
        axis.title=element_text(size=15,face="bold"))


## risk
data$copyrisk <- factor(data$risk_all, levels=c(1,2,3,4,5,6,7,9,0),
                        labels=c("MSM", 
                                 "HET",
                                 "IDU",
                                 "IDU",
                                 "Blood Transfusion/ Products",
                                 "Blood Transfusion/ Products",
                                 "Other/ Unknown", #MTCT
                                 "Other/ Unknown",
                                 "Other/ Unknown"))
ta7 <- table(data$copyrisk)
# MSM                         HET 
# 5215                        4563 
# IDU Blood Transfusion/ Products 
# 2487                         108
# Other/ Unknown 
# 479



n_risk <-data.frame(matrix(NA, ncol = 3, nrow = 5))
colnames(n_risk) <- c("Risk","N","Percentage")
n_risk$Risk <- c("MSM", "HET","IDU",
                 "Blood Transfusion/ Products",
                 "Other/ Unknown") #MTCT

n_risk$N <- c(t(ta7)[1],t(ta7)[2],t(ta7)[3],t(ta7)[4], t(ta7)[5])
n_risk$Percentage <- c(t(ta7)[1]/total, t(ta7)[2]/total,
                       t(ta7)[3]/total, t(ta7)[4]/total, t(ta7)[5]/total)
n_risk$Percentage <- round(n_risk$Percentage,3)


ggplot(n_risk, aes(x = reorder(Risk,-N), y = N, fill=Risk)) +
  geom_bar( stat = "identity", position=position_dodge(), colour="black" )+
  scale_fill_manual(values=c( rep("#999999",3), "lightblue", "#999999"))+
  guides(fill=FALSE)+
  geom_text(aes(label = paste0(N)), vjust = c(rep(+1.2,3),rep(-1.5,2)), size=7)+
  geom_text(aes(label = paste0("(",Percentage*100, "%)")), vjust =c(rep(+3.0,3),rep(-0.5,2)), size=5, colour="blue")+
  theme_bw()+
  scale_y_continuous( limits=c(0,10000),
                      sec.axis = sec_axis(~./10000 , name = "Frequency [%]",
                                          labels = scales::percent)) +
  #ggtitle("Mode of HIV Tranmission in the SHCS DB") +
  xlab("Risk") +
  ylab("Count")+
  theme(axis.text=element_text(size=14),
        axis.text.x=element_text(angle=45, hjust=1),
        axis.title=element_text(size=15,face="bold"),
        axis.title.y.right =element_text(colour ="blue"))


ggplot(n_risk, aes(x = reorder(Risk,-Percentage), y = Percentage, fill=Risk)) +
  geom_bar( stat = "identity", position=position_dodge(), colour="black" )+
  scale_fill_manual(values=c( rep("#999999",3), "lightblue", "#999999"))+
  guides(fill=FALSE)+
  geom_text(aes(label = paste0(Percentage*100, "%")), vjust =c(rep(+1.2,3),rep(-0.5,2)), size=7)+
  theme_bw()+
  scale_y_continuous(labels = scales::percent, limits=c(0,1)) +
  #ggtitle("Mode of HIV Tranmission in the SHCS DB") +
  xlab("Risk") +
  ylab("Frequency [%]")+
  theme(axis.text=element_text(size=14),
        axis.text.x=element_text(angle=45, hjust=1),
        axis.title=element_text(size=15,face="bold"))


## diagnosis date
summary(data$diag_date)
# Min.      1st Qu.       Median         Mean      3rd Qu.         Max. 
# "1981-07-01" "1991-07-14" "1998-07-08" "1999-01-07" "2006-01-04" "2020-03-18" 

boxplot(data$diag_date)

boxplot(data$diag_date,
        xlab = "Diagnosis year",
        cex.axis=1.2, cex.lab=1.5, cex.main=2,
        face="bold",
        col = "orange",
        border = "brown",
        horizontal = TRUE,
        notch = F)





### SANPLEDATE:
summary(data$sampledate)
# Min.      1st Qu.       Median         Mean      3rd Qu.         Max. 
# "1989-08-31" "1997-10-05" "2002-07-08" "2003-07-16" "2008-02-04" "2020-05-26" 

hist(data$sampledate, breaks = 15)

ggplot(data, aes(sampledate)) +
  geom_histogram(aes(x=sampledate, y=..density..), bins=50, fill="#d3d3d3", color="black") +
  #stat_function(fun=dnorm, args = list(mean=mean(data$sampledate), sd=sd(data$sampledate)), color="red", size=1.5) +
  geom_density(color="blue") +
  ggtitle("Histogram of Sample date") +
  xlab("Sample date") +
  ylab("Density") +
  theme_bw() #coord_cartesian(xlim=c(1989,2020), ylim=c(0,0.1))



###SEQUENCEDATE:
summary(data$seq_dt)
# Min.      1st Qu.       Median         Mean      3rd Qu.         Max. 
# "1996-05-02" "2007-06-25" "2008-09-15" "2009-04-02" "2012-07-02" "2020-06-04"





#######################
## non-B's 
#######################


nonB <- data[data$copy_subtype == "Non-B",] #2'859


summary(nonB$age_at_diag)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.0    28.0    34.0    36.19    42.0    81.0 
hist(nonB$age_at_diag)
ggplot(nonB, aes(age_at_diag)) +
  geom_histogram(aes(x=age_at_diag, y=..density..), bins=50, fill="#d3d3d3", color="black") +
  stat_function(fun=dnorm, args = list(mean=mean(nonB$age_at_diag), sd=sd(nonB$age_at_diag)), color="red", size=1.5) +
  geom_density(color="blue") +
  ggtitle("Histogram of Age at Diagnosis Date for Non-B Subtype") +
  xlab("Age at diagnosis date") +
  ylab("Density") +
  theme_bw() +
  coord_cartesian(xlim=c(15,95), ylim=c(0,0.06))
#not normally distributed 

ggpubr::ggqqplot(nonB$age_at_diag)
shapiro.test(nonB$age_at_diag[1:5000]) #W = 0.94876, p-value < 2.2e-16
#also (more of less) normal distributed, can take the mean as well






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




