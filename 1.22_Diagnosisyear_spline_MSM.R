#### MSM - Spline Regression

## Clear R's memory
rm(list=ls())


# packages
library(tidyverse)
library(ggplot2)
library(cowplot)
library('forestplot')
library(CIplot)
library(haven)

library(grid)
library(vcd)
library(rms)
library(epiDisplay)

## load data
setwd("~/Desktop/SHCS")
# data<- read.csv("Output/fulltable_study_population.csv")
data <- read.csv("Output/fulltable_study_population_update.csv")
data$X <-  NULL #does the same: data <- data[,-1]


################################
table(data$subtype)
# 01_AE       02_AG       03_AB       05_DF      06_CPX       07_BC      09_CPX      11_CPX       12_BF 
# 137          44           0           0           1           1           0           1           3 
# 13_CPX       14_BG      18_CPX      19_CPX       20_BG       24_BG      25_CPX      27_CPX       29_BF 
# 0           0           2           4           5           3           0           0           0 
# 31_BC       35_AD      37_CPX       40_BF      43_02G           A           B           C           D 
# 1           0           1           0           1          45        3894          28           5 
# F           G           H           J           K RECOMBINANT 
# 45          11           0           0           0         119 


## Same for the 6 non-B subtypes with diagnosis year:
#get non-B subpopulation: "A", "C", "F", "G", "CRF01_AE","CRFO2_AG"
# nonB_data <- data[data$subtype != "B",] # 457 sequences
nonB_data <-data[data$subtype %in% c("A", "C","F","G","01_AE","02_AG"),] # 310 in total

################################
FUN_nonB <- function(x){
  if (is.na(x)){return (0)}
  else if (x == "B" ){return(0)}
  else if (x == "Non-B"){return(1)}
  else{return(0)}
}
data$nonB <- sapply(data$copy_subtype, FUN_nonB)
#data$nonB <- ifelse(test=data$nonB == 0, yes="B", no="Non-B")
#data$nonB <- as.factor(data$nonB-1 )
#data$nonB <- factor(data$nonB, levels = c("B","Non-B"))


# MSM!
data <- data[data$risk == "1" ,] 
data <- data[!(is.na(data$id)),] #4351




##	year of diagnosis 
### 2. Visualize the data
# data$diag_year <-as.numeric(data$diag_year)
# #diag_year as numerics!
# data$diag_year_num <- factor(data$diag_year); data$diag_year_num <- as.numeric(data$diag_year_num); data$diag_year_num
# 
# 
# plot(data$diag_year, data$nonB)
# mod <- glm(data$nonB~ as.numeric(data$diag_year), family = "binomial"(link="logit")); summary(model)
# summary(data$diag_year)
# # Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# # 1982    1993    2001    2001    2008    2020
# 
# 
# #AIC remains equal!
# plot(data$diag_year_num, data$nonB)
# model_num <- glm(data$nonB~ data$diag_year_num, family = "binomial"(link="logit")); summary(model_num)
# summary(data$diag_year_num)
# # Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# # 1.0    12.0    20.0    19.6    27.0    39.0



require(splines)
s3 <- glm(nonB ~ bs(diag_year, degree = 1, df=4), family = "binomial",data = data); summary(s3)
attr(terms(s3), "predvars")
# list(nonB, bs(diag_year, degree = 1L, 
#knots = c(`25%` = 1995,  `50%` = 20010 `75%` = 2005), Boundary.knots = c(1982, 2020), 
#               intercept = FALSE))
termplot(s3, se=T)


s4 <- glm(nonB ~ bs(diag_year, degree = 2, df=4), family = "binomial",data = data); summary(s4)
attr(terms(s4), "predvars")
# list(nonB, bs(diag_year, degree = 2L, knots = c(`33.33333%` = 1996, 
#                                                 `66.66667%` = 2003), Boundary.knots = c(1982, 2020), intercept = FALSE))
termplot(s4, se=T)


s5 <- glm(nonB ~ bs(diag_year, degree = 3, df=4), family = "binomial",data = data); summary(s5)
attr(terms(s5), "predvars")
# list(nonB, bs(diag_year, degree = 3L, knots = c(`50%` = 2001), 
#               Boundary.knots = c(1981, 2019), intercept = FALSE))
termplot(s5, se=T)


mod0 <- glm(nonB ~  1  , family = "binomial",data = data)
mod1 <- glm(nonB ~  ethnicity + bs(diag_year, degree = 2, df=4) , family = "binomial",data = data)
mod2 <- glm(nonB ~  ethnicity , family = "binomial",data = data)

mod3 <- glm(nonB ~  regroupedregion , family = "binomial",data = data)

logistic<- glm(nonB ~  diag_year  , family = "binomial",data = data)
# round(exp(coef(test)[2]),2) # 1.134967 OR 
logistic.display(logistic) #1.13 (1.12,1.15) , < 0.001 Wald's or LR-test

# get the likelihood ratio
lmtest::lrtest (mod0,logistic)
A <- logLik(mod0)
B <- logLik(logistic)
teststat1 <- 2 * (as.numeric(B)-as.numeric(A))
p.val1 <- pchisq(teststat1, df = 1, lower.tail = FALSE)

with(logistic, null.deviance - deviance)
with(logistic, df.null - df.residual)
with(logistic, pchisq(null.deviance - deviance, df.null - df.residual, lower.tail = F))


lmtest::lrtest (mod0,s4)
A <- logLik(mod0)
C<- logLik(s4)
teststat2 <- -2 * (as.numeric(A)-as.numeric(C))
p.val2 <- pchisq(teststat2, df = 4, lower.tail = FALSE)


lmtest::lrtest (logistic,s4)
teststat3 <- 2 * (as.numeric(C)-as.numeric(B))
df3 <- abs(print(anova(logistic)[,"Df"][2]) - print(anova(s4)[,"Df"][2]))
p.val3 <- pchisq(teststat3, df = df3, lower.tail = FALSE)



#lrtest (mod01,mod02)
lmtest::lrtest (mod2,mod1)
D <- logLik(mod2)
E<- logLik(mod1)
teststat3 <- -2 * (as.numeric(D)-as.numeric(E))
p.val3 <- pchisq(teststat3, df = 4, lower.tail = FALSE)


lmtest::lrtest (mod0,mod1)
A <- logLik(mod0) #null
E<- logLik(mod1)
teststat4 <- -2 * (as.numeric(A)-as.numeric(E))
p.val4 <- pchisq(teststat4, df = 1, lower.tail = FALSE)


#lrtest (mod0,mod2)
#lrtest (mod0,mod3)

################################
###### For the publication
pal <- function(SplineModel, LogisticModel)
{
  par(mar=c(5.1, 5.1, 4.1, 16.1), xpd=TRUE)
  termplot(SplineModel,se = F, lwd.term=2, las=1, ylim= c(-4,4),xlabs = " ", ylabs = " ", col.term ="firebrick4",col.se = "navy", lwd.se = 0.4)
  par(new=TRUE)
  termplot(LogisticModel, lty.se =0.8, las=1, ylim= c(-4,4),xlabs = " ", ylabs = " ",col.term ="indianred",lwd.term=2)
  segments(x0=1981.6,y0=0,x1=2020.4,y1=0,col="black", lwd=2,lty=2 )
  mtext(side = 1, text = c("Year of HIV Diagnosis"), line = 3, cex=1.5)
  mtext(side = 2,  text = c("Log-Odds"),
        #, "Having HIV-1 Non-B Subtype "), 
        line = c(3), cex= 1.5)
  # legend(par("usr")[2]+0.2 ,par("usr")[4],
  #        c("",
  #          "univariable (Year of Diagnosis)", "logistic regression, df = 1",
  #          "univariable (Year of Diagnosis)","B-spline regression, df = 4",
  #          " ",
  #          "compare to null model","LR > 10, p-value < 2.2e-16 ***",
  #          "compare to null model","LR > 10 , p-value < 2.2e-16 ***",
  #          "compare models","LR > 10, p-value = 0.002 **"),
  #        fill= c("white", "indianred",
  #                "white", "firebrick4",
  #                "white", "white","indianred",
  #                "white","firebrick4",
  #                "white","black","white"),
  #        border = c("white", "indianred",
  #                   "white", "firebrick4",
  #                   "white", "white","indianred",
  #                   "white","firebrick4",
  #                   "white","black","white"),
  #        cex = 1.1, bty="n")
}
overall <-pal(s4,logistic)






###########

# par(mar=c(5.1, 5.1, 4.1, 16.1), xpd=TRUE)
# termplot(test, lty.se =0.8, las=1, ylim= c(-4,4),xlabs = " ", ylabs = " ",col.term ="indianred",lwd.term=2)
# mtext(side = 1, text = c("Year of HIV Diagnosis"), line = 3, cex=1.5)
# mtext(side = 2,  text = c("Log-Odds of", "Having HIV-1 Non-B Subtype "), line = c(3.5,2.5), cex= 1.5)
# segments(x0=1981.6,y0=0,x1=2020.4,y1=0,col="black", lwd=2,lty=2 )
# legend(par("usr")[2],par("usr")[4], 
#        c("",
#          "univariable (Year of Diagnosis)", "logistic regression, df = 1", "",
#          "", "",
#          "","",
#          " "," "," ",
#          "",
#          "", 
#          ""),
#        fill= c("white", "indianred",
#                "white", "white",
#                "white", "white",
#                "white", "white","white","white","white",
#                "white","white",
#                "white","white","white"),
#        border = "white",
#        cex = 1, bty="n")
# par(new=TRUE)
# termplot(s4,se = F, lwd.term=2, las=1, ylim= c(-4,4),xlabs = " ", ylabs = " ", col.term ="firebrick4",col.se = "navy", lwd.se = 0.4)
# legend(par("usr")[2],par("usr")[4], 
#        c("",
#          " ", "", "",
#          "univariable (Year of Diagnosis)","B-spline regression, df = 4", 
#          " "," "," ",
#          "",
#          "", 
#          ""),
#        fill= c("white", "indianred", "white",
#                "white", "firebrick4",
#                "white", "white",
#                "white", "white","white","white","white",
#                "white","white",
#                "white","white","white"),
#        border = "white",
#        cex = 1, bty="n")
# legend(par("usr")[2],par("usr")[4], 
#        c("",
#          "", "", 
#          "",
#          "", "",
#          " "," "," ",
#          "compare to null model","LR > 10, p-value < 2.2e-16 ***", "",
#          "compare to null model","LR > 10 , p-value < 2.2e-16 ***", "",
#          "compare models","LR = 9.90, p-value = 0.02 *"),
#        fill= c("white", "indianred", "white",
#                "white", "firebrick4",
#                "white", "white","white","white","indianred","white",
#                "white","firebrick4",
#                "white","white","white"),
#        border = c("white", "indianred", "white",
#                   "white", "firebrick4",
#                   "white", "white","white","white","indianred","white",
#                   "white","firebrick4",
#                   "white","white","black"),
#        cex = 1, bty="n")
# 
# 


###############################################################
library(RColorBrewer)
display.brewer.all(colorblindFriendly = TRUE)
display.brewer.pal(n=6, name = "Dark2"); brewer.pal(n=6, name = "Dark2")
library(plotrix)
sapply(brewer.pal(n=6, name = "Dark2"), color.id)


data$As <- ifelse(data$subtype == "A",1,0)
data$AEs <- ifelse(data$subtype == "01_AE",1,0)
data$Cs <- ifelse(data$subtype == "C",1,0)
data$Fs <- ifelse(data$subtype == "F",1,0)
data$Gs <- ifelse(data$subtype == "G",1,0)
data$AGs <- ifelse(data$subtype == "02_AG",1,0)
regAE <- glm(AEs ~ diag_year , family = "binomial",data = data); summary(regAE)
# exp(regAE$coefficients[2]) #1.114426 
# exp(confint(regAE)) # 1.089497e+00 1.141273e+00
# regAE$coefficients
# (Intercept)    diag_year 
# -304.0886388    0.1493071 
logistic.display(regAE) #1.11 (1.08,1.13)
splineAE <- glm(AEs ~ bs(diag_year, degree = 2, df=2), family = "binomial",data = data); summary(splineAE)
attr(terms(splineAE), "predvars"); termplot(splineAE, se=T)
# list(subtype, bs(diag_year, degree = 2L, knots = c(`33.33333%` = 2006, 
#                                                    `66.66667%` = 2013), Boundary.knots = c(1984, 2020), intercept = FALSE))

regAG <- glm(AGs ~ diag_year , family = "binomial",data = data); summary(regAG)
# exp(regAG$coefficients[2]) # 1.161029 
# exp(regAE$coefficients)
logistic.display(regAG) #1.15 (1.10,1.20) 
# exp(confint(regAG))
# regAG$coefficients
# (Intercept)    diag_year 
# -304.0886388    0.1493071 THE HIGHEST INCREASE!
splineAG <- glm(AGs ~ bs(diag_year, degree = 2, df=2), family = "binomial",data = data); summary(splineAG)
attr(terms(splineAG), "predvars"); termplot(splineAG, se=T) ##no increase

regG <- glm(Gs ~ diag_year , family = "binomial",data = data); summary(regG)
# exp(coef(regG)) # 1.08499 
# regG$coefficients 
# (Intercept)     diag_year 
# -169.43452268    0.08157072  
logistic.display(regG) #1.08 (1.02,1.15)
splineG <- glm(Gs ~ bs(diag_year, degree = 2, df=2), family = "binomial",data = data); summary(splineG)
attr(terms(splineG), "predvars"); termplot(splineG, se=T)

regF <- glm(Fs ~ diag_year , family = "binomial",data = data); summary(regF)
# exp(regF$coefficients[2]) #1.109073 
# regF$coefficients
# (Intercept)    diag_year 
# -212.0825502    0.1035243 
logistic.display(regF) #1.10 (1.06,1.14)
splineF <- glm(Fs ~ bs(diag_year, degree = 2, df=2), family = "binomial",data = data); summary(splineF)
attr(terms(splineF), "predvars"); termplot(splineF, se=T) #best increase

regC <- glm(Cs ~ diag_year , family = "binomial",data = data); summary(regC)
# exp(regC$coefficients[2]) #1.143739 
# regC$coefficients
# (Intercept)    diag_year 
# -274.3927423    0.1343025 
logistic.display(regC) #1.14 (1.09,1.20)
splineC <- glm(Cs ~ bs(diag_year, degree = 2, df=2), family = "binomial",data = data); summary(splineC)
attr(terms(splineC), "predvars"); termplot(splineC, se=T) #quiet staty 

regA <- glm(As ~ diag_year , family = "binomial",data = data); summary(regA)
# exp(regA$coefficients[2]) #1.160308 
# regA$coefficients
# (Intercept)    diag_year 
# -302.8159460    0.1486854 
logistic.display(regA) #1.16 (1.11,1.21)
splineA <- glm(As ~ bs(diag_year, degree = 2, df=2), family = "binomial",data = data); summary(splineA)
attr(terms(splineA), "predvars"); termplot(splineA, se=T) #increase



#### spline regression:
fig_6_splines <- function(SplineModelAE, SplineModelAG,SplineModelA,SplineModelC,SplineModelF,SplineModelG)
{
  par(mar=c(5.1, 5.1, 4.1, 16.1), xpd=TRUE)
  termplot(SplineModelAE,lwd.term=3, las=1, ylim= c(-10,4), xlabs = " ", ylabs = " ", col.term ="darkorange3", se = F, col.se= "darkorange3", lwd.se =  0.4)
  par(new=TRUE)
  termplot(SplineModelAG,lwd.term=3, las=1,ylim= c(-10,4), xlabs = " ", ylabs = " ", col.term ="darkgoldenrod2",col.se = "darkgoldenrod2", se = F, lwd.se =  0.4)
  par(new=TRUE)
  termplot(SplineModelA,lwd.term=3, las=1, ylim= c(-10,4), xlabs = " ", ylabs = " ", col.term ="violetred2",col.se = "violetred2", se = F, lwd.se =  0.4)
  par(new=TRUE)
  termplot(SplineModelC,lwd.term=3, las=1, ylim= c(-10,4), xlabs = " ", ylabs = " ", col.term ="cyan4",col.se = "cyan4", se = F, lwd.se =  0.4)
  par(new=TRUE)
  termplot(SplineModelF,lwd.term=3, las=1, ylim= c(-10,4), xlabs = " ", ylabs = " ", col.term ="olivedrab",col.se = "olivedrab", se = F, lwd.se =  0.4)
  par(new=TRUE)
  termplot(SplineModelG,lwd.term=3, las=1, ylim= c(-10,4), xlabs = " ", ylabs = " ", col.term ="mediumpurple3",col.se = "mediumpurple3", se = F, lwd.se =  0.4)
  mtext(side = 1, text = c("Year of HIV Diagnosis"), line = 3, cex=1.5)
  mtext(side = 2,  text = c("Log-Odds"),
        #, "Having HIV-1 Non-B Subtype "), 
        line = c(3), cex= 1.5)
  segments(x0=1981.6,y0=0,x1=2020.4,y1=0,col="black", lwd=2,lty=2 )
  legend(par("usr")[2]+2.2 ,par("usr")[4]+1.2, 
         c("",
           "CRF01_AE", "",
           "CRF02_AG", "",
           "A","",
           "C","",
           "F","",
           "G"),
         fill= c("white", "darkorange3",
                 "white", "darkgoldenrod2",
                 "white", "violetred2",
                 "white", "cyan4",
                 "white","olivedrab",
                 "white","mediumpurple3"),
         border = c("white", "darkorange3",
                    "white", "darkgoldenrod2",
                    "white", "violetred2",
                    "white", "cyan4",
                    "white","olivedrab",
                    "white","mediumpurple3"),
         cex = 1.3, bty="n")

}

splines6 <-fig_6_splines(splineAE, splineAG,splineA,splineC,splineF,splineG)


#### logistic regression:
fig_6_logistic <- function(LogisticModelAE, LogisticModelAG,LogisticModelA,LogisticModelC,LogisticModelF,LogisticModelG)
{
  par(mar=c(5.1, 5.1, 4.1, 16.1), xpd=TRUE)
  termplot(LogisticModelAE,lwd.term=2, las=1, ylim= c(-4,4), xlabs = " ", ylabs = " ", col.term ="darkorange3", se = F, col.se= "darkorange3", lwd.se =  0.4)
  par(new=TRUE)
  termplot(LogisticModelAG,lwd.term=3, las=1, ylim= c(-4,4), xlabs = " ", ylabs = " ", col.term ="darkgoldenrod2",col.se = "darkgoldenrod2", se = F, lwd.se =  0.4)
  par(new=TRUE)
  termplot(LogisticModelA,lwd.term=2, pch=10,las=1, ylim= c(-4,4), xlabs = " ", ylabs = " ", col.term ="violetred2",col.se = "violetred2", se = F, lwd.se =  0.4)
  par(new=TRUE)
  termplot(LogisticModelC,lwd.term=2, las=1, ylim= c(-4,4), xlabs = " ", ylabs = " ", col.term ="cyan4",col.se = "cyan4", se = F, lwd.se =  0.4)
  par(new=TRUE)
  termplot(LogisticModelF,lwd.term=2, las=1, ylim= c(-4,4), xlabs = " ", ylabs = " ", col.term ="olivedrab",col.se = "olivedrab", se = F, lwd.se =  0.4)
  par(new=TRUE)
  termplot(LogisticModelG,lwd.term=2, las=1, ylim= c(-4,4), xlabs = " ", ylabs = " ", col.term ="mediumpurple3",col.se = "mediumpurple3", se = F, lwd.se =  0.4)
  mtext(side = 1, text = c("Year of HIV Diagnosis"), line = 3, cex=1.5)
  mtext(side = 2,  text = c("Log-Odds"),
                            #, "Having HIV-1 Non-B Subtype "), 
        line = c(3), cex= 1.5)
  segments(x0=1981.6, y0=0, x1=2020.4, y1=0, col="black", lwd=2, lty=2)
  legend(par("usr")[2]+2.2 ,par("usr")[4]+0.5, 
         c("",
           "CRF01_AE", "",
           "CRF02_AG", "",
           "A","",
           "C","",
           "F","",
           "G"),
         fill= c("white", "darkorange3",
                 "white", "darkgoldenrod2",
                 "white", "violetred2",
                 "white", "cyan4",
                 "white","olivedrab",
                 "white","mediumpurple3"),
         border = c("white", "darkorange3",
                    "white", "darkgoldenrod2",
                    "white", "violetred2",
                    "white", "cyan4",
                    "white","olivedrab",
                    "white","mediumpurple3"),
         cex = 1.3, bty="n")
}
logistic6 <-fig_6_logistic(regAE, regAG,regA,regC,regF,regG)






# library(grid)
# library(ggpubr)
# # Move to a new page
# grid.newpage()
# # Create layout : nrow = 2, ncol = 2
# pushViewport(viewport(layout = grid.layout(nrow = 2, ncol = 2)))
# 
# # A helper function to define a region on the layout
# define_region <- function(row, col){
#   viewport(layout.pos.row = row, layout.pos.col = col)
# } 
# 
# # Arrange the plots
# print(p1, vp = define_region(row = 1, col = 1:2))   # Span over two columns
# print(p2, vp = define_region(row = 2, col = 1))
# print(p3, vp = define_region(row = 2, col = 2))


# ggarrange(pal(s4,logistic), 
#           logistic6, splines6, ncol = 2, nrow = 1,
#           labels = c("A", "B","C"), 
#           common.legend = TRUE, legend = "bottom",
#           font.label = list(size = 18, color = "black"))
library(broom)
summary(logistic)
summary(s4)

#augmented data

logistic.df <- augment(logistic)
glimpse(logistic.df)

ModelPredictions <- predict(logistic , type="response")

mydata <- data %>%
  dplyr::select_if(is.numeric)
predictors <- colnames(mydata)

mydata <- mydata %>%
  mutate(logit = log(ModelProb/ (1-ModelProb))) %>%
  gather( key= "predictors", value = "predictors.value", -logit)


ggplot(mydata, aes(x=logit , y=predictors.value )) +
  geom_point()  +  
  stat_smooth(size=1.5, linetype = "11",
              method="glm", formula = y ~ x,
              se=FALSE, method.args = list(family=binomial)) +
  theme_bw()+
  facet_wrap(~ predictors, scales ="free_y")


ggplot(data, aes(x=data$diag_year , y=ModelPredictions)) +
  geom_line( color = "red")  +  
    # stat_smooth(size=1.5, linetype = "11",
    #         se=FALSE) +
    theme_bw()


ggplot(logistic.df, aes(diag_year, nonB)) + 
  #geom_point() + 
  stat_smooth(size=1.5, linetype = "11",
              method="glm", formula = y ~ x,
              se=FALSE, method.args = list(family=binomial)) +
  



logistic.df$nonB

ggplot(data, aes(x= diag_year, y= log()))


library(cowplot)
library(gridGraphics)

## A
pal(s4,logistic)
####record previous plot
p1 <- recordPlot()

fig_6_splines(splineAE, splineAG,splineA,splineC,splineF,splineG)
p2 <- recordPlot()



fig_6_logistic(regAE, regAG,regA,regC,regF,regG)
p3 <- recordPlot()
ggarrange(p1,p3,
          ncol=1,
          nrow=2,
          labels = "AUTO", 
          font.label = list(size = 18, color = "black"))

ggarrange(p1,p3,
          nrow=2,
          labels = "AUTO",
          font.label = list(size = 18, color = "black"))


ggarrange(p1,
          ggarrange(p2,p3, ncol = 2,labels = c("B", "C"), 
                    font.label = list(size = 18, color = "black")),
          nrow=2,
          labels = "A",
          font.label = list(size = 18, color = "black"))

# 
# ggarrange(p2,p3, ncol = 2, nrow = 1,
#           labels = c("A", "B"), 
#           common.legend = F, legend = "right",
#           font.label = list(size = 15, color = "black"))


# 
# bottom_row<-plot_grid(p2,p3,
#                       labels = c("B","C"),
#                       hjust = 0, vjust = 1)
# 
# plot_grid(p1,bottom_row,
#           labels = c("A",""),
#           label_size = 12, ncol=1,
#           hjust = 0, vjust = 1)

