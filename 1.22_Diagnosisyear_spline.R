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

## load data
setwd("~/Desktop/SHCS")
data<- read.csv("Output/fulltable_study_population.csv")
data$X <-  NULL #does the same: data <- data[,-1]




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
data <- data[data$risk == "1",] #4751


##	year of diagnosis 
### 2. Visualize the data
data$diag_year <-as.numeric(data$diag_year)
#diag_year as numerics!
data$diag_year_num <- factor(data$diag_year); data$diag_year_num <- as.numeric(data$diag_year_num); data$diag_year_num 


plot(data$diag_year, data$nonB)
mod <- glm(data$nonB~ as.numeric(data$diag_year), family = "binomial"(link="logit")); summary(model)
assummary(data$diag_year)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 1982    1993    2001    2001    2008    2020


#AIC remains equal!
plot(data$diag_year_num, data$nonB)
model_num <- glm(data$nonB~ data$diag_year_num, family = "binomial"(link="logit")); summary(model_num)
summary(data$diag_year_num)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 1.0    12.0    20.0    19.6    27.0    39.0 



require(splines)
s3 <- glm(nonB ~ bs(diag_year, degree = 1, df=4), family = "binomial",data = data); summary(s3)
attr(terms(s3), "predvars")
# list(nonB, bs(diag_year, degree = 1L, 
#knots = c(`25%` = 1993,  `50%` = 2001, `75%` = 2008), Boundary.knots = c(1982, 2020), 
#               intercept = FALSE))
termplot(s3, se=T)


s4 <- glm(nonB ~ bs(diag_year, degree = 2, df=4), family = "binomial",data = data); summary(s3)
attr(terms(s4), "predvars")
# list(nonB, bs(diag_year, degree = 2L, knots = c(`33.33333%` = 1995, 
#                                                 `66.66667%` = 2006), Boundary.knots = c(1982, 2020), intercept = FALSE))
termplot(s4, se=T)


s5 <- glm(nonB ~ bs(diag_year, degree = 3, df=4), family = "binomial",data = data); summary(s3)
attr(terms(s5), "predvars")
# list(nonB, bs(diag_year, degree = 3L, knots = c(`50%` = 2001), 
#               Boundary.knots = c(1982, 2020), intercept = FALSE))
termplot(s5, se=T)


mod0 <- glm(nonB ~  1  , family = "binomial",data = data)
mod1 <- glm(nonB ~  ethnicity + bs(diag_year, degree = 2, df=4) , family = "binomial",data = data)
mod2 <- glm(nonB ~  ethnicity , family = "binomial",data = data)

mod3 <- glm(nonB ~  regroupedregion , family = "binomial",data = data)

test<- glm(nonB ~  diag_year  , family = "binomial",data = data)


# get the likelihood ratio
lmtest::lrtest (mod0,test)
A <- logLik(mod0)
B <- logLik(test)
teststat1 <- -2 * (as.numeric(A)-as.numeric(B))
p.val1 <- pchisq(teststat1, df = 1, lower.tail = FALSE)

lmtest::lrtest (mod0,s4)
A <- logLik(mod0)
C<- logLik(s4)
teststat2 <- -2 * (as.numeric(A)-as.numeric(C))
p.val2 <- pchisq(teststat2, df = 4, lower.tail = FALSE)

lmtest::lrtest (test,s4)
teststat3 <- -2 * (as.numeric(B)-as.numeric(C))
p.val3 <- pchisq(teststat3, df = 4, lower.tail = FALSE)



#lrtest (mod01,mod02)
lmtest::lrtest (mod2,mod1)
D <- logLik(mod2)
E<- logLik(mod1)
teststat3 <- -2 * (as.numeric(D)-as.numeric(E))
p.val3 <- pchisq(teststat3, df = 4, lower.tail = FALSE)

#lrtest (mod0,mod2)
#lrtest (mod0,mod3)

################################
# logistic and spline regression:
## use regression output: Non-B
# visualize
par(mar=c(5.1, 5.1, 4.1, 16.1), xpd=TRUE)
termplot(mod1,lwd.term=2, las=1, ylim= c(-4,4),xlabs = " ", ylabs = " ", col.term ="green", se = F, col.se= "lightgreen", lwd.se =  0.4)
par(new=TRUE)
termplot(s4,se = F, lwd.term=2, las=1, ylim= c(-4,4),xlabs = " ", ylabs = " ", col.term ="navyblue",col.se = "navy", lwd.se = 0.4)
par(new=TRUE)
termplot(test, lty.se =0.8, las=1, ylim= c(-4,4),xlabs = " ", ylabs = " ",col.term ="steelblue",lwd.term=2)
# abline(model, col = "orange", lwd = 1, ylim= c(-4,4), xlim = c(1,15))
# segments(x0=1987,y0= -4.3,x1=2020,y1=-0.15,col="red", lwd = 2)

mtext(side = 1, text = c("Year of HIV Diagnosis"), line = 3, cex=1.5)
mtext(side = 2,  text = c("Odds of", "Having HIV-1 Non-B Subtype "), line = c(3.5,2.5), cex= 1.5)
legend(par("usr")[2],par("usr")[4], 
       c("",
         "univariable (Year of Diagnosis)", "logistic regression, df = 1", 
         "univariable (Year of Diagnosis)","B-spline regression, df = 4", 
         "multivariable (Ethnicity and Year of Diagnosis)", "B-spline regression, df = 4",
         " "," "," ",
         "compare to null model","LR > 10, p-value < 2.2e-16 ***",
         "compare to null model","LR > 10 , p-value < 2.2e-16 ***", 
         "compare to ethnicity model","LR > 10, p-value < 2.2e-16 ***"),
       fill= c("white", "steelblue",
               "white", "navyblue",
               "white", "green",
               "white", "white","white","white","steelblue",
               "white","navyblue",
               "white","green","white"),
       border = c("white", "steelblue",
                  "white", "navyblue",
                  "white", "green",
                  "white", "white","white","white","steelblue",
                  "white","navyblue",
                  "white","green","white"),
       cex = 0.8, bty="n")










par(mar=c(5.1, 5.1, 4.1, 16.1), xpd=TRUE)
termplot(test, lty.se =0.8, las=1, ylim= c(-4,4),xlabs = " ", ylabs = " ",col.term ="steelblue",lwd.term=2)
mtext(side = 1, text = c("Year of HIV Diagnosis"), line = 3, cex=1.5)
mtext(side = 2,  text = c("Odds of", "Having HIV-1 Non-B Subtype "), line = c(3.5,2.5), cex= 1.5)
legend(par("usr")[2],par("usr")[4], 
       c("",
         "univariable (Year of Diagnosis)", "logistic regression, df = 1", "",
         "", "",
         "","",
         " "," "," ",
         "",
         "", 
         ""),
       fill= c("white", "steelblue",
               "white", "white",
               "white", "white",
               "white", "white","white","white","white",
               "white","white",
               "white","white","white"),
       border = "white",
       cex = 1, bty="n")
par(new=TRUE)
termplot(s4,se = F, lwd.term=2, las=1, ylim= c(-4,4),xlabs = " ", ylabs = " ", col.term ="navyblue",col.se = "navy", lwd.se = 0.4)
legend(par("usr")[2],par("usr")[4], 
       c("",
         " ", "", "",
         "univariable (Year of Diagnosis)","B-spline regression, df = 4", 
         " ", "",
         " "," "," ",
         "",
         "", 
         ""),
       fill= c("white", "steelblue", "white",
               "white", "navyblue",
               "white", "white",
               "white", "white","white","white","white",
               "white","white",
               "white","white","white"),
       border = "white",
       cex = 1, bty="n")
legend(par("usr")[2],par("usr")[4], 
       c("",
         "", "", 
         "",
         "","", 
         "", "",
         " "," "," ",
         "compare to null model","LR > 10, p-value < 2.2e-16 ***", "",
         "compare to null model","LR  10 , p-value < 2.2e-16 ***", "",
         "compare models","LR = 9.90, p-value = 0.02 *"),
       fill= c("white", "steelblue", "white",
               "white", "navyblue",
               "white", "white",
               "white", "white","white","white","steelblue","white",
               "white","navyblue",
               "white","white","white"),
       border = c("white", "steelblue", "white",
                  "white", "navyblue",
                  "white", "white",
                  "white", "white","white","white","steelblue","white",
                  "white","navyblue",
                  "white","white","black"),
       cex = 1, bty="n")



################################

data$diag_year_cat <- NA
data$diag_year_cat[data$diag_year <1993] <- 1
data$diag_year_cat[data$diag_year >= 1993 & data$diag_year<2008] <- 2
data$diag_year_cat[data$diag_year >= 2008] <- 3
data$diag_year_cat <- as.factor(data$diag_year_cat)
#data$diag_year_cat <- factor(data$diag_year_cat, levels = c("<1993","1993-2008",">2008"))
log_reg_diag_year_cat <- glm(data$nonB ~ data$diag_year_cat, family = "binomial"); summary(log_reg_diag_year_cat)
lines(data$nonB, fitted(log_reg_diag_year_cat), col = "red", lwd=3)


diag_year <-data$diag_year
nonB <- data$nonB
#Distribution Summaries for Predictor Variables:
dd <- datadist(nonB,diag_year)
options(datadist="dd"); dd
## used by summary, plot, survplot, sometimes predict:
# nonB diag_year
# Low:effect         0      1993
# Adjust to          0      2001
# High:effect        1      2008
# Low:prediction     0      1983
# High:prediction    1      2019
# Low                0      1982
# High               1      2019
# 
# Values:
#   
#   nonB : 0 1 

splines3<- lrm(nonB ~  rcs(diag_year,3), data = data, x=TRUE,y=TRUE); splines3
#lrm: Fit binary and proportional odds ordinal logistic regression models using maximum likelihood estimation or penalized maximum likelihood estimation
#rcs:set up special attributes (such as knots and nonlinear term indicators) that are carried through to fits 
# attr(,"parms")
# [1] 1988 2001 2013

summary(splines3, diag_year = c( 2001,1988))
# Factor      Low  High Diff. Effect  S.E.    Lower 0.95 Upper 0.95
# diag_year   2001 1988 -13   -1.5451 0.24467 -2.02460   -1.06550  
# Odds Ratio 2001 1988 -13    0.2133      NA  0.13205    0.34455  
summary(splines3, diag_year = c(2001, 2013))
# Factor      Low  High Diff. Effect S.E.    Lower 0.95 Upper 0.95
# diag_year   2001 2013 12    1.5672 0.11102 1.3496     1.7848    
# Odds Ratio 2001 2013 12    4.7932      NA 3.8559     5.9584 




################################
### OR of having non-B with splines
par(mfrow=c(1,1))
par(new=FALSE)
FUN_diag_splines <- function(i){
  sum_year<- summary(splines3, diag_year = c(2001, i))
  OR_year <- sum_year[2,4] #get the effect: OR
  plot(i, OR_year, ylim= c(0.1,8),log = "y",pch = 16, cex = .9,
       las=1, xlim = c(1982,2019), ylab = "", xlab = "")
  par(new=TRUE)
}
sapply(data$diag_year[!is.na(data$diag_year)], FUN_diag_splines)
title( ylab = c("Having Non-B Subtype [OR]"), xlab = "Year of HIV Diagnosis", line = 3)
segments(x0=1981.6,y0=1,x1=2020.4,y1=1,col="lightsteelblue3", lwd=2)

################################





#replot
#37/ 4334 = 0.008537148
xv <- seq(min(data$diag_year), max(data$diag_year),0.008539)
yv <- predict(model,list(diag_year=xv), type="response")
lines(xv,yv, col="red")
detach(data)

plot(data$diag_year, data$nonB,
     ylim=c(-0.25,1.25))
noise = 0.03*rnorm(data$diag_year)
yobs = data$nonB + noise
points(data$diag_year,yobs, col="red")


### 3. Fit a linear regression model and plot the data with the fitted result
fit <- lm(nonB ~ diag_year, data = data)
plot(data$diag_year, data$nonB,
     xlab="Diagnosis Year",
     ylab= "Subtype", 
     main= "HIV-1 B or Non-B Subtype")
abline(fit, col = "red", lwd = 3)
plot(fit)
## assumtions not met! 


### 4. Fit a polynomial regression models with degrees 2, 3, and 10 and visualize as before
m1 <- lm(nonB ~ poly(diag_year, degree = 2, row =T), data = data); summary(m1)
plot(data$diag_year, data$nonB,
     xlab="Diagnosis Year",
     ylab= "Subtype", 
     main= "HIV-1 B or Non-B Subtype")
lines(data$diag_year, fitted(m1), col = "red", lwd=3)


m2 <- lm(nonB ~ poly(diag_year, degree = 3, row = T), data = data); summary(m2)
plot(data$diag_year, data$nonB,
     xlab="Diagnosis Year",
     ylab= "Subtype", 
     main= "HIV-1 B or Non-B Subtype")
lines(data$diag_year,fitted(m2), col = "red", lwd=3)
#slightly better

log_reg_diagnosed_year <- glm(data$nonB~ data$diag_year, family = "binomial");  summary(log_reg_diagnosed_year)



