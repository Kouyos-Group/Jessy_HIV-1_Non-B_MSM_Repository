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
data$diag_year <- as.Date(as.character(data$diag_year), tryFormats = c("%Y"))
data$diag_year  <- substr(data$diag_date,1,4)
data$diag_year<- as.Date(as.character(data$diag_year), tryFormats = c("%Y"))


###SUBTYPE
ta<-table(data$subtype)
CrossTable(data$subtype)


mosaicplot(t(ta),col=TRUE, main="HIV Subtypes distribution within the SHCS")



ta2<-table(data$copy_subtype)
CrossTable(data$copy_subtype)
#B: 8309 (73.9%), non-B: 2859 (25.7%)
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

aggregate(data$pol_len, list(data$lab), mean)
# Group.1        x
# 1    Basel 1286.6730
# 2   Geneva  889.8317
# 3 Lausanne 1432.9896
# 4   Zurich  1302.2978




### GENDER:
ta3<-table(data$sex) 
CrossTable(data$sex)
#1 = male: 7998 (71,6%), 2= female: 3,170 (28.5%)

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
# 18.00   45.00   53.00   52.66   60.00  101.00 


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
shapiro.test(data$age[1:5000]) #W = 0.99255, p-value = 1.589e-15
#normal distributed so able to take the mean




### Age at Sample date 
summary(data$age_sampledate)

# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 5.00   31.00   37.00   38.54   44.00   82.00 
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
shapiro.test(data$age_sampledate[1:5000]) #W = 0.95919, p-value < 2.2e-16
#also (more or less) normal distributed, can take the mean as well



### Age at Diagnosis date 
summary(data$age_at_diag)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.00   26.00   32.00   34.02   40.00   81.00 



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
shapiro.test(data$age_at_diag[1:5000]) #W = 0.93753, p-value < 2.2e-16
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
# 8662            1327             3097             330
# Others          Unknown 
# 25              513



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
# 8972
# Western and Central Africa 
# 710
# Eastern and Southern Africa 
# 445
# Latin America 
# 432
# Eastern and Southern Asia 
# 278 
# Western Asia and North Africa 
# 142 
# Eastern Europe 
# 113 
# North America and Oceania 
# 70 
# Unknown 
# 6 


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
# 4351                        3923 
# IDU Blood Transfusion/ Products 
# 2392                         102
# Other/ Unknown 
# 400



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




