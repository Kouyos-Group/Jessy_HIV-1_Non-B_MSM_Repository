#MSM:

par(mfrow=c(2,2))
par(cex = 1)
par(mar=c(4.1,4.1,4.1,2.1))
#par(mar=c(5.1, 5.1, 4.1, 16.1), xpd=TRUE)
termplot(s4,se = F, lwd.term=2, las=1, ylim= c(-4,4),xlabs = " ", ylabs = " ", col.term ="firebrick4",col.se = "navy", lwd.se = 0.4)
par(new=TRUE)
termplot(logistic, lty.se =0.8, las=1, ylim= c(-4,4),xlabs = " ", ylabs = " ",col.term ="indianred",lwd.term=2)
segments(x0=1981.6,y0=0,x1=2020.4,y1=0,col="black", lwd=2,lty=2 )
mtext(side = 1, text = c("Year of HIV Diagnosis"), line = 2.5, cex=1.5)
mtext(side = 2,  text = c("Log-Odds"),
      #, "Having HIV-1 Non-B Subtype "), 
      line = c(2.5), cex= 1.5)
par(xpd=NA)
text(par("usr")[1]-10.2, par("usr")[4]+6.2, "A", cex=1.5, font = 2)
par(cex = 1.5)
legend(par("usr")[2]+10 ,par("usr")[4]+8.2,
       c("",
         "univariable (Year of Diagnosis)", "logistic regression, df = 1",
         "univariable (Year of Diagnosis)","B-spline regression, df = 4",
         " ",
         "compare to null model","LR > 10, p-value < 2.2e-16 ***",
         "compare to null model","LR > 10 , p-value < 2.2e-16 ***",
         "compare models","LR > 10, p-value = 0.002 **"),
       fill= c("white", "indianred",
               "white", "firebrick4",
               "white", "white","indianred",
               "white","firebrick4",
               "white","black","white"),
       border = c("white", "indianred",
                  "white", "firebrick4",
                  "white", "white","indianred",
                  "white","firebrick4",
                  "white","black","white"),
       cex = 0.6, bty="n")

par(cex = 1)
plot.new()
par(cex = 1)
par(mar=c(4.1,4.1,4.1,2.1))
termplot(regAE,lwd.term=1, las=1, ylim= c(-4,4), xlabs = " ", ylabs = " ", col.term ="darkorange3", se = F, col.se= "darkorange3", lwd.se =  0.4)
par(new=TRUE)
termplot(regAG,lwd.term=1, las=1,ylim= c(-4,4), xlabs = " ", ylabs = " ", col.term ="darkgoldenrod1",col.se = "darkgoldenrod1", se = F, lwd.se =  0.4)
par(new=TRUE)
termplot(regA,lwd.term=1, las=1, ylim= c(-4,4), xlabs = " ", ylabs = " ", col.term ="violetred2",col.se = "violetred2", se = F, lwd.se =  0.4)
par(new=TRUE)
termplot(regC,lwd.term=1, las=1, ylim= c(-4,4), xlabs = " ", ylabs = " ", col.term ="cyan4",col.se = "cyan4", se = F, lwd.se =  0.4)
par(new=TRUE)
termplot(regF,lwd.term=1, las=1, ylim= c(-4,4), xlabs = " ", ylabs = " ", col.term ="olivedrab",col.se = "olivedrab", se = F, lwd.se =  0.4)
par(new=TRUE)
termplot(regF,lwd.term=1, las=1, ylim= c(-4,4), xlabs = " ", ylabs = " ", col.term ="mediumpurple3",col.se = "mediumpurple3", se = F, lwd.se =  0.4)
mtext(side = 1, text = c("Year of HIV Diagnosis"), line = 2.5, cex=1.5)
mtext(side = 2,  text = c("Log-Odds"),
      #, "Having HIV-1 Non-B Subtype "), 
      line = c(2.5), cex= 1.5)
segments(x0=1981.6,y0=0,x1=2020.4,y1=0,col="black", lwd=2,lty=2 )
par(xpd=NA)
text(par("usr")[1]-10.2, par("usr")[4]+6.2, "B", cex=1.5, font = 2)
par(cex = 1.5)
legend(par("usr")[2]+10 ,par("usr")[4]+8.2,
       c("",
         "CRF01_AE", 
         "CRF02_AG", 
         "A",
         "C",
         "F",
         "G"),
       fill= c("white", "darkorange3",
                "darkgoldenrod2",
                "violetred2",
                "cyan4",
               "olivedrab",
               "mediumpurple3"),
       border = c("white", "darkorange3",
                   "darkgoldenrod2",
                   "violetred2",
                   "cyan4",
                  "olivedrab",
                  "mediumpurple3"),
       cex = 0.8, bty="n")






#HET:
par(mfrow=c(2,2))
par(cex = 1)
par(mar=c(4.1,4.1,4.1,2.1))
#par(mar=c(5.1, 5.1, 4.1, 16.1), xpd=TRUE)
termplot(s4,se = F, lwd.term=2, las=1, ylim= c(-4,4),xlabs = " ", ylabs = " ", col.term ="navyblue",col.se = "navy", lwd.se = 0.4)
par(new=TRUE)
termplot(logistic, lty.se =0.8, las=1, ylim= c(-4,4),xlabs = " ", ylabs = " ",col.term ="steelblue",lwd.term=2)
segments(x0=1981.6,y0=0,x1=2020.4,y1=0,col="black", lwd=2,lty=2 )
mtext(side = 1, text = c("Year of HIV Diagnosis"), line = 2.5, cex=1.5)
mtext(side = 2,  text = c("Log-Odds"),
      #, "Having HIV-1 Non-B Subtype "), 
      line = c(2.5), cex= 1.5)
par(xpd=NA)
text(par("usr")[1]-10.2, par("usr")[4]+6.2, "A", cex=1.5, font = 2)
par(cex = 1.5)
legend(par("usr")[2]+10 ,par("usr")[4]+8.2,
       c("",
         "univariable (Year of Diagnosis)", "logistic regression, df = 1", 
         "univariable (Year of Diagnosis)","B-spline regression, df = 4", 
         " ",
         "compare to null model","LR > 10, p-value < 2.2e-16 ***",
         "compare to null model","LR > 10, p-value < 2.2e-16 ***", 
         "compare models","LR > 10, p-value < 2.2e-16 ***"),
       fill= c("white", "steelblue",
               "white", "navyblue",
               "white", "white","steelblue",
               "white","navyblue",
               "white","black","white"),
       border = c("white", "steelblue",
                  "white", "navyblue",
                  "white", "white","steelblue",
                  "white","navyblue",
                  "white","black","white"),
       cex = 0.6, bty="n")

par(cex = 1)
plot.new()
par(cex = 1)
par(mar=c(4.1,4.1,4.1,2.1))
termplot(regAE,lwd.term=1, las=1, ylim= c(-4,4), xlabs = " ", ylabs = " ", col.term ="darkorange3", se = F, col.se= "darkorange3", lwd.se =  0.4)
par(new=TRUE)
termplot(regAG,lwd.term=1, las=1,ylim= c(-4,4), xlabs = " ", ylabs = " ", col.term ="darkgoldenrod1",col.se = "darkgoldenrod1", se = F, lwd.se =  0.4)
par(new=TRUE)
termplot(regA,lwd.term=1, las=1, ylim= c(-4,4), xlabs = " ", ylabs = " ", col.term ="violetred2",col.se = "violetred2", se = F, lwd.se =  0.4)
par(new=TRUE)
termplot(regC,lwd.term=1, las=1, ylim= c(-4,4), xlabs = " ", ylabs = " ", col.term ="cyan4",col.se = "cyan4", se = F, lwd.se =  0.4)
par(new=TRUE)
termplot(regF,lwd.term=1, las=1, ylim= c(-4,4), xlabs = " ", ylabs = " ", col.term ="olivedrab",col.se = "olivedrab", se = F, lwd.se =  0.4)
par(new=TRUE)
termplot(regF,lwd.term=1, las=1, ylim= c(-4,4), xlabs = " ", ylabs = " ", col.term ="mediumpurple3",col.se = "mediumpurple3", se = F, lwd.se =  0.4)
mtext(side = 1, text = c("Year of HIV Diagnosis"), line = 2.5, cex=1.5)
mtext(side = 2,  text = c("Log-Odds"),
      #, "Having HIV-1 Non-B Subtype "), 
      line = c(2.5), cex= 1.5)
segments(x0=1981.6,y0=0,x1=2020.4,y1=0,col="black", lwd=2,lty=2 )
par(xpd=NA)
text(par("usr")[1]-10.2, par("usr")[4]+6.2, "B", cex=1.5, font = 2)
par(cex = 1.5)
legend(par("usr")[2]+10 ,par("usr")[4]+8.2,
       c("",
         "CRF01_AE", 
         "CRF02_AG", 
         "A",
         "C",
         "F",
         "G"),
       fill= c("white", "darkorange3",
               "darkgoldenrod2",
               "violetred2",
               "cyan4",
               "olivedrab",
               "mediumpurple3"),
       border = c("white", "darkorange3",
                  "darkgoldenrod2",
                  "violetred2",
                  "cyan4",
                  "olivedrab",
                  "mediumpurple3"),
       cex = 0.8, bty="n")







# 
# par(mfrow=c(2,3))
# par(cex = 1)
# par(mar=c(4.1,4.1,4.1,2.1))
# #par(mar=c(5.1, 5.1, 4.1, 16.1), xpd=TRUE)
# termplot(s4,se = F, lwd.term=2, las=1, ylim= c(-4,4),xlabs = " ", ylabs = " ", col.term ="firebrick4",col.se = "navy", lwd.se = 0.4)
# par(new=TRUE)
# termplot(logistic, lty.se =0.8, las=1, ylim= c(-4,4),xlabs = " ", ylabs = " ",col.term ="indianred",lwd.term=2)
# segments(x0=1981.6,y0=0,x1=2020.4,y1=0,col="black", lwd=2,lty=2 )
# mtext(side = 1, text = c("Year of HIV Diagnosis"), line = 2.5, cex=1.5)
# mtext(side = 2,  text = c("Log-Odds"),
#       #, "Having HIV-1 Non-B Subtype "), 
#       line = c(2.5), cex= 1.5)
# par(xpd=NA)
# text(-0, 10.2, "A", cex=1.5, font = 2)
# plot.new()
# par(cex = 1.5)
# legend(-0.5 ,2,
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
#        cex = 0.6, bty="n")
# par(cex = 1)
# par(mar=c(4.1,4.1,4.1,2.1))
# plot.new()
# text(-0, 1.75, "B", cex=1.5, font = 2)
# par(cex = 1)
# par(mar=c(4.1,4.1,4.1,2.1))
# termplot(splineAE,lwd.term=2, las=1, ylim= c(-8,4), xlabs = " ", ylabs = " ", col.term ="darkorange3", se = F, col.se= "darkorange3", lwd.se =  0.4)
# par(new=TRUE)
# termplot(splineAG,lwd.term=2, las=1,ylim= c(-8,4), xlabs = " ", ylabs = " ", col.term ="darkgoldenrod2",col.se = "darkgoldenrod2", se = F, lwd.se =  0.4)
# par(new=TRUE)
# termplot(splineA,lwd.term=2, las=1, ylim= c(-8,4), xlabs = " ", ylabs = " ", col.term ="violetred2",col.se = "violetred2", se = F, lwd.se =  0.4)
# par(new=TRUE)
# termplot(splineC,lwd.term=2, las=1, ylim= c(-8,4), xlabs = " ", ylabs = " ", col.term ="cyan4",col.se = "cyan4", se = F, lwd.se =  0.4)
# par(new=TRUE)
# termplot(splineF,lwd.term=2, las=1, ylim= c(-8,4), xlabs = " ", ylabs = " ", col.term ="olivedrab",col.se = "olivedrab", se = F, lwd.se =  0.4)
# par(new=TRUE)
# termplot(splineF,lwd.term=2, las=1, ylim= c(-8,4), xlabs = " ", ylabs = " ", col.term ="mediumpurple3",col.se = "mediumpurple3", se = F, lwd.se =  0.4)
# mtext(side = 1, text = c("Year of HIV Diagnosis"), line = 2.5, cex=1.5)
# mtext(side = 2,  text = c("Log-Odds"),
#       #, "Having HIV-1 Non-B Subtype "), 
#       line = c(2.5), cex= 1.5)
# segments(x0=1981.6,y0=0,x1=2020.4,y1=0,col="black", lwd=2,lty=2 )
# par(xpd=NA)
# text(-0, 10.2, "C", cex=1.5, font = 2)
# plot.new()
# par(cex = 1.5)
# legend(-0.5 ,2,
#        c("",
#          "CRF01_AE", 
#          "CRF02_AG", 
#          "A",
#          "C",
#          "F",
#          "G"),
#        fill= c("white", "darkorange3",
#                "darkgoldenrod2",
#                "violetred2",
#                "cyan4",
#                "olivedrab",
#                "mediumpurple3"),
#        border = c("white", "darkorange3",
#                   "darkgoldenrod2",
#                   "violetred2",
#                   "cyan4",
#                   "olivedrab",
#                   "mediumpurple3"),
#        cex = 0.8, bty="n")
# par(cex = 1)
# par(mar=c(4.1,4.1,4.1,2.1))
# plot.new()
# text(-0, 1.75, "D", cex=1.5, font = 2)
# 
