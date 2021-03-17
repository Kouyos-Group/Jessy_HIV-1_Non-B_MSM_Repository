
### Clear R's memory
rm(list=ls())


### GET THE PACKAGES AND FUNCTIONS
if (!requireNamespace("pacman", quietly = TRUE))
  install.packages("pacman")
pacman::p_load("dplyr","haven","tidyr", "data.table","GGally", "gmodels","readr","ape","seqinr")
`%notin%` <- Negate(`%in%`)

library(readstata13)

### IMPORT DATASET
### Load the SHCS sequences, the test_info file and the pat.dta:
setwd("~/Desktop/SHCS")  
SHCS_SEQ <- read.fasta(file="Input/final0620.fas") #29 837 
TEST_INFO <- read.dta13("Input/test_info0620.dta") #24 564
PAT <- read.dta13("Input/pat0520.dta") #20 904
ids_for_rules<- read.csv("Output/ids_for_rules.csv")
admin <- as.data.table(read_dta("Input/admi0520.dta")) #20,904
#hiv_subtype <- as.data.table(read_dta("Input/hiv_subtype.dta")) #10773 obs

### Sandra's Analysis: 
# cat = 1 2 sequences, cat = 2&3 >2 sequences
#from cat 1 to 3 ==> less to most likely HIV-1 Superinfection
ids_inSI_325potential<- read.csv("OLD Output/ids_inSI_325potential.csv") #108 potential HIV-1 superinfection patients, 48 belonging to cat=1
ids_inSI_325potential <- ids_inSI_325potential[,-1]
#{ids_inSI_71Sequenced <- read.csv("Output/ids_inSI_71Sequenced.csv") #16 out of 108 NGS re-sequenced 
# ids_inSI_44Confirmed <- read.csv("Output/ids_inSI_44Confirmed.csv") #9 out of 16 re-sequenced are part of the 44 confirmed HIV-1 superinfections
# ids_inSI_44Confirmed <- ids_inSI_44Confirmed[,-1]
# ids_inSI_71Sequenced <- ids_inSI_71Sequenced[,-1]


length(unique(TEST_INFO$id))#13278

### A) QUALITY CHECK:
# partial pol ( #PR > 250nt and RT > 500 nt):
# sum(is.na(TEST_INFO$pr_len)) #553
# sum(is.na(TEST_INFO$rt_len)) #546
# sum(is.na(TEST_INFO$int_len)) #21,139
# length(which(TEST_INFO$pr_len < 250)) #6
# length(which(TEST_INFO$rt_len < 500)) #1
# sum(TEST_INFO$subtype == "HIV2 subtype A") #3
# sum(TEST_INFO$subtype  == "") #0
# sum(TEST_INFO$subtype == "HIV 1 N group") #1
# sum(TEST_INFO$subtype == "HIV 1 O group") #5
# boxplot(TEST_INFO$rt_len)
# boxplot(TEST_INFO$pr_len)

## get rid of NA's in  RTs & PRs:
TEST_INFO <- TEST_INFO[!is.na(TEST_INFO$pr_len),] #24,011
TEST_INFO <- TEST_INFO[!is.na(TEST_INFO$rt_len),] #23,706
length(unique(TEST_INFO$id))#12990


##patients, which are in pat.dta & test_info.dta
# sum(duplicated(TEST_INFO$header_id)) #0
# sum(!TEST_INFO$id %in% PAT$id) #13
TEST_INFO <- TEST_INFO[TEST_INFO$id %in% PAT$id,] #20,693
length(unique(TEST_INFO$id))#12978



###### ###########
#Add some variables
###########

TEST_INFO<- TEST_INFO %>%
  mutate(pol_len = pr_len+rt_len) ## add the total_nseq per unique id to the table data


## order by id, sampledate, seq_dt
#But first make sure that dates are dates & change all the subtypes toupper()
TEST_INFO$sampledate <- as.Date(as.character(TEST_INFO$sampledate), tryFormats = c("%Y-%m-%d"))
TEST_INFO$seq_dt <- as.Date(as.character(TEST_INFO$seq_dt), tryFormats = c("%Y-%m-%d"))
TEST_INFO$subtype <- toupper(TEST_INFO$subtype)
TEST_INFO$rega_subtype <- toupper(TEST_INFO$rega_subtype)
TEST_INFO$comet_subtype <- toupper(TEST_INFO$comet_subtype)
TEST_INFO <- TEST_INFO[order(TEST_INFO$id, TEST_INFO$sampledate,TEST_INFO$seq_dt),]


## add Sandra's Analysis cat & patients that requieres a rule
TEST_INFO <- merge(TEST_INFO,ids_inSI_325potential,all.x=T,by="id")
TEST_INFO <- TEST_INFO %>%
  mutate(rule=case_when(
    id %in% ids_for_rules$x  ~ "1",
    id %notin% ids_for_rules$x ~ "0",
  )) #20,603


names(TEST_INFO)

## continue to work with the subset:
data <- TEST_INFO %>%
  select(id, header_id, sampledate, seq_dt, rega_subtype, comet_subtype, subtype, lab, pol_len, cat, rule, linked) %>% #include the lab!
  group_by(id) #12603 obs. ordered by id  #data <- data[order(data$id),] 


# Numbering the sequencies per unique id:
data$nseq<- unlist(tapply(data$id, data$id,
                          function(x) seq(1,length(x),1)))


# Total Number of sequencies per unique id:
total_nseq<-cbind(unique(data$id),total_nseq=sapply(X=unique(data$id), function(x) sum(data$id==x)))
colnames(total_nseq)<-c("id","total_nseq")
data <- merge(data, total_nseq, bi="id")  ## add the total_nseq per unique id to the table data


### From here part B) SUBTYPE CLASSIFICATION ALGORITHM:



linked <- TEST_INFO[which(TEST_INFO$linked=="1"),] #3019
length(unique(linked$id)) #2671
length(unique(data$id)) #12978


select <- unique(na.omit(linked$id));length(unique(na.omit(linked$id))) #ids: 2671: 
newdata <- data[data$id %in% select,] #4951 obs.
length(unique(newdata$id))#2671



#create a vector to store id with changing subtype
rules <- rep(NA, nrow(newdata))
#selecting same ID and different subtypes, if fulfilled the statement, it is true and we output the ID. Otherwise it saves NA
#save it in select
for (i in 2:nrow(newdata)){
  rules[i] <- ifelse(newdata$id[i]==newdata$id[i-1] & newdata$subtype[i]!=newdata$subtype[i-1],newdata$id[i],NA)
}
rules <- unique(na.omit(rules)) #76
rules_dd <- newdata[newdata$id %in% rules,] #255 obs.
rules_dd <- rules_dd[order(rules_dd$id, rules_dd$linked, rules_dd$pol_len, rules_dd$sampledate,rules_dd$seq_dt),] #255




linkeddata <- data[which(data$linked=="1"),] #3023
length(unique(linkeddata$id))#2671
### 2)  ABOVE THRESHOLD
threshold<-  summary(data$pol_len)[[2]] #887bp , threshold = 1st Qu.
linkedabove <- linkeddata[sapply(linkeddata$pol_len, 
                  function(x) any(x>=threshold) )==T, ] #3017 above the threshold obs.
length(unique(linkedabove$id))#2665





select2 <- unique(na.omit(linkedabove$id));length(unique(na.omit(linkedabove$id))) #ids: 2665: 
newdata2 <- data[data$id %in% select2,] #4937 obs.
length(unique(newdata2$id))#2665

#create a vector to store id with changing subtype
rulesabove <- rep(NA, nrow(newdata2))
#selecting same ID and different subtypes, if fulfilled the statement, it is true and we output the ID. Otherwise it saves NA
#save it in select
for (i in 2:nrow(newdata2)){
  rulesabove[i] <- ifelse(newdata2$id[i]==newdata2$id[i-1] & newdata2$subtype[i]!=newdata2$subtype[i-1],newdata2$id[i],NA)
}
rulesabove <- unique(na.omit(rulesabove)) #76
rulesabove_dd <- linkedabove[linkedabove$id %in% rulesabove,] #101 obs.
rulesabove_dd <- rulesabove_dd[order(rulesabove_dd$id, rulesabove_dd$linked, rulesabove_dd$pol_len, rulesabove_dd$sampledate,rulesabove_dd$seq_dt),] 




### 3)  EARLIEST SAMPLE DATE
linkedsample <- merge(aggregate(sampledate ~ id, min, data = linkedabove), linkedabove) #2675
linkedsample1 <- linkedsample
ids<- rep(NA, nrow(linkedsample)); newsubtype <-  rep(NA, nrow(linkedsample))
for (i in 1:nrow(linkedsample)){
    if (linkedsample$subtype[i] == "RECOMBINANT" & linkedsample$total_nseq[i] > 1){
      ids[i] <- linkedsample$id[i] 
      # ids <- unique(na.omit(ids))
      newsubtype[i] <- ifelse(linkedsample$comet_subtype[i] == data[data$id== ids[i] & data$nseq != linkedsample$nseq[i], "subtype"],linkedsample$comet_subtype[i],"RECOMBINANT")
      linkedsample$subtype[i] <- ifelse(linkedsample$subtype[i]== "RECOMBINANT" ,newsubtype[i],linkedsample$subtype[i])
    }
}
unique(na.omit(newsubtype))
index <- which(linkedsample$id %in% ids)
# i =2083 #STAYS RECOMBI
# i =1970 #change
# i =570 #change
# i =1828 #change

length(unique(linkedsample$id))#2665

### rules for unknown subtype impossible!
# linkedsample2 <- linkedsample1
# idsU<- rep(NA, nrow(linkedsample1)); newsubtypeU <-  rep(NA, nrow(linkedsample1))
# for (i in 1:nrow(linkedsample1)){
#   if (linkedsample1$subtype[i] == "UNKNOWN" & linkedsample1$total_nseq[i] > 1){
#     idsU[i] <- linkedsample1$id[i] 
#     # ids <- unique(na.omit(ids))
#     newsubtypeU[i] <- ifelse(linkedsample1$comet_subtype[i] == data[data$id== idsU[i] & data$nseq != linkedsample1$nseq[i], "subtype"],linkedsample1$comet_subtype[i],"UNKNOWN")
#     linkedsample2$subtype[i] <- ifelse(linkedsample1$subtype[i]== "UNKNOWN" ,newsubtypeU[i],linkedsample1$subtype[i])
#   }
# }
# unique(na.omit(newsubtypeU));iiindex <-which(linkedsample$id %in% idsU)
# idsU[iiindex]
# i =1246


# duplicates<- linkedsample1[duplicated(linkedsample1$id),]#10



### 4)   EARLIEST SEQUENCE DATE
# have only 2 sequences, taking the most frequent one difficult, 
linkedsequence <- merge(aggregate(seq_dt ~ id, min, data = linkedsample), linkedsample) #2666
sum(duplicated(linkedsequence$id)) # have  1 case where have duplicates ##ids: 52130
linkedduplicates<- linkedsequence[!duplicated(linkedsequence$id),]#2665


##########################



NOTlinked <- TEST_INFO[which(is.na(TEST_INFO$linked)),] #20670
NOT_select <- data[which(data$id %notin% select),] #18742
length(unique(NOT_select$id))#10307

#create a vector to store id with changing subtype
rulesWOlink <- rep(NA, nrow(NOT_select))
#selecting same ID and different subtypes, if fulfilled the statement, it is true and we output the ID. Otherwise it saves NA
#save it in select
for (i in 2:nrow(NOT_select)){
  rulesWOlink[i] <- ifelse(NOT_select$id[i]==NOT_select$id[i-1] & NOT_select$subtype[i]!=NOT_select$subtype[i-1],NOT_select$id[i],NA)
}
rulesWOlink <- unique(na.omit(rulesWOlink)) #230
rulesWOdd <- NOT_select[NOT_select$id %in% rulesWOlink,] #853 obs.
rulesWOdd <- rulesWOdd[order(rulesWOdd$id, rulesWOdd$linked, rulesWOdd$pol_len, rulesWOdd$sampledate,rulesWOdd$seq_dt),]

# normal <- data[data$id %notin% c(rulesWOlink,select),] #17886 obs.
# length(unique(normal$id)) #10,077

### 2)  ABOVE THRESHOLD
threshold<-  summary(data$pol_len)[[2]] #887bp , threshold = 1st Qu.
notlinkedabove <- NOT_select[sapply(NOT_select$pol_len, 
                                 function(x) any(x>=threshold) )==T, ] #18668 above the threshold obs.


length(unique(notlinkedabove$id))#10297


#create a vector to store id with changing subtype
NOTrulesabove <- rep(NA, nrow(notlinkedabove))
#selecting same ID and different subtypes, if fulfilled the statement, it is true and we output the ID. Otherwise it saves NA
#save it in select
for (i in 2:nrow(notlinkedabove)){
  NOTrulesabove[i] <- ifelse(notlinkedabove$id[i]==notlinkedabove$id[i-1] & notlinkedabove$subtype[i]!=notlinkedabove$subtype[i-1],notlinkedabove$id[i],NA)
}
NOTrulesabove <- unique(na.omit(NOTrulesabove)) #224
NOTrulesabove_dd <- notlinkedabove[notlinkedabove$id %in% NOTrulesabove,] #833 obs.
NOTrulesabove_dd <- NOTrulesabove_dd[order(NOTrulesabove_dd$id, NOTrulesabove_dd$linked, NOTrulesabove_dd$pol_len, NOTrulesabove_dd$sampledate,NOTrulesabove_dd$seq_dt),] 




### 3)  EARLIEST SAMPLE DATE
NOTlinkedsample <- merge(aggregate(sampledate ~ id, min, data = notlinkedabove), notlinkedabove) #10414
NOTlinkedsample1 <- NOTlinkedsample
NOTids<- rep(NA, nrow(NOTlinkedsample)); NOTnewsubtype <-  rep(NA, nrow(NOTlinkedsample))
for (i in 1:nrow(NOTlinkedsample)){
  if (NOTlinkedsample$subtype[i] == "RECOMBINANT" & NOTlinkedsample$total_nseq[i] > 1){
    NOTids[i] <- NOTlinkedsample$id[i] 
    # ids <- unique(na.omit(ids))
    NOTnewsubtype[i] <- ifelse(NOTlinkedsample$comet_subtype[i] == data[data$id== NOTids[i] & data$nseq != NOTlinkedsample$nseq[i], "subtype"],NOTlinkedsample$comet_subtype[i],"RECOMBINANT")
    NOTlinkedsample$subtype[i] <- ifelse(NOTlinkedsample$subtype[i]== "RECOMBINANT" ,newsubtype[i],NOTlinkedsample$subtype[i])
  }
}
unique(na.omit(NOTnewsubtype))
length(unique(NOTlinkedsample$id))#10297


### 4)   EARLIEST SEQUENCE DATE
# have only 2 sequences, taking the most frequent one difficult, 
NOTlinkedsequence <- merge(aggregate(seq_dt ~ id, min, data = NOTlinkedsample), NOTlinkedsample) #10299
sum(duplicated(NOTlinkedsequence$id)) # have 2 cases where have duplicates ##ids: 90289 and 57333
NOTlinkedduplicates<- NOTlinkedsequence[!duplicated(NOTlinkedsequence$id),]#10297



#### merge them:
final <- rbind(linkedduplicates,NOTlinkedduplicates) #12962
finaldataset <- final[final$subtype %notin% c("HIV2 SUBTYPE A", "HIV 1 N GROUP","HIV 1 O GROUP"),] #12954
finaldatasetTEST <- finaldataset
finaldatasetTEST$subtype<- ifelse(finaldatasetTEST$subtype== "RECOMBINANT" & finaldatasetTEST$comet_subtype=="B" & finaldataset$rega_subtype=="RECOMBINANT OF B, D"  ,"B",finaldatasetTEST$subtype)
finaldatasetTEST$subtype<- ifelse(finaldatasetTEST$subtype== "RECOMBINANT" & finaldatasetTEST$comet_subtype=="B" & finaldataset$rega_subtype=="RECOMBINANT OF B, F1"  ,"B",finaldatasetTEST$subtype)



table(finaldatasetTEST$subtype)
table(finaldataset$subtype)
###
#name the sequences
names(SHCS_SEQ)<-sapply(names(SHCS_SEQ),
                          FUN=function(a){unlist(strsplit(a,split="\\|"))[1]}) #29837
#remove duplicates of the sequences
SHCS_SEQ <-SHCS_SEQ[!duplicated(names(SHCS_SEQ))] #29837
# only take sequences who are in SHCS (header_ids from test_info)
SHCS_SEQ_FINAL <-SHCS_SEQ[names(SHCS_SEQ) %in% finaldataset$header_id] #12954


############################
### Save the dataset:
#write.csv(finaldatasetTEST,"Output/study_population_update.csv")
# write.csv(finaldatasetTEST,"Output/real_studypopulation_FEB21.csv")

############################





###Make a Table with Recombi B, Comet =B
recombi <- finaldataset[finaldataset$rega_subtype == "RECOMBINANT OF B, D" & finaldataset$comet_subtype == "B", c("id","comet_subtype","rega_subtype","subtype","linked")] #155
r  <- finaldataset[finaldataset$rega_subtype == "RECOMBINANT OF B, F1", c("id","comet_subtype","rega_subtype","subtype","linked")] #155
recombisub  <- finaldataset[finaldataset$subtype == "RECOMBINANT", c("id","comet_subtype","rega_subtype","subtype","linked")] #493
#### new risk group
PAT$riskNEW <- ifelse(is.na(PAT$risk), 'Other',
                      ifelse(PAT$risk == 1, 'MSM', 
                             ifelse(PAT$risk == 2,'HET','Other')))

#ifelse(PAT$risk %in% c(3,4), 'IDU','Other'))))
recombi <- merge(recombi, PAT[,c('id','risk', 'riskNEW', 'ethnicity', 'sex')], all.x = T, by="id")
r <- merge(r, PAT[,c('id','risk', 'riskNEW', 'ethnicity', 'sex')], all.x = T, by="id")

recombi$copyethnicity <- factor(recombi$ethnicity, levels=c(1,2,3,4,0,9),
                             labels=c("White", 
                                      "Black",
                                      "Hispanic",
                                      "Asian",
                                      "Others",
                                      "Others"))

recombi$copyrisk <- factor(recombi$risk, levels=c(1,2,3,4,5,6,7,9,0),
                        labels=c("MSM", 
                                 "HET",
                                 "IDU",
                                 "IDU",
                                 "Blood Transfusion/ Products",
                                 "Blood Transfusion/ Products",
                                 "Other/ Unknown", #MTCT
                                 "Other/ Unknown",
                                 "Other/ Unknown"))
mosaic(recombi$copyethnicity~recombi$riskNEW, direction = "v", rot_labels=c(0,90,0,0))
spineplot(recombi$riskNEW~recombi$copyethnicity)
require(vcd)

a<-table(recombi$copyrisk,recombi$copyethnicity)
a1 <- prop.table(a, margin = 1) %>% round(3)
a1[is.na(a1)] <- 0
library(gridExtra)
library(grid)
grid.newpage()
grid.table(a)
grid.newpage()
grid.table(a1)

##### Taking earlist seq per patient: 
#remove duplicates of the sequences
test_info <- read.dta13("Input/test_info0620.dta") #24 564
test_info <- test_info[!is.na(test_info$pr_len),] 
test_info <- test_info[!is.na(test_info$rt_len),] 
test_info <- test_info[!duplicated(test_info$header_id),]
test_info <- test_info[test_info$id %in% PAT$id,]
test_info <- test_info[order(test_info$sampledate),]
test_info <- test_info[!duplicated(test_info$id),] #12978
test_info <- test_info[test_info$id %in% finaldataset$id,] #12954




###################################################################################################################
###################################################################################################################

#save old settings
op <- par(no.readonly = TRUE)
#change settings
par(mar=c(8, 4, 2, 2) + 0.1)
##linked:
barplot(table(linkeddata[linkeddata$id %in% rules,"subtype"]),ylim = c(0,35),main="Without applied Sutype Classification Algorithm",ylab="Frequency",xlab="Subtype",las=2)
barplot(table(finaldataset[finaldataset$id %in% rules,"subtype"]),ylim = c(0,35),main="With applied Sutype Classification Algorithm",ylab="Frequency",xlab="Subtype",las=2)

## without link:
barplot(table(NOTlinked[NOTlinked$id %in% rulesWOlink,"subtype"]),ylim = c(0,350),main="Without applied Sutype Classification Algorithm",ylab="Frequency",xlab="Subtype",las=2)
barplot(table(finaldataset[finaldataset$id %in% rulesWOlink,"subtype"]),ylim = c(0,350),main="With applied Sutype Classification Algorithm",ylab="Frequency",xlab="Subtype",las=2)
#reset settings
par(op)





#save old settings
op <- par(no.readonly = TRUE)
#change settings
par(mar=c(8, 4, 2, 2) + 0.1)
##linked:
withalgo <-as.data.frame(table(finaldataset[finaldataset$id %in% rules,"subtype"]))
withoutalgo <-as.data.frame(table(test_info[test_info$id %in% rules,"subtype"]))
withoutalgo$Var1 <- toupper(withoutalgo$Var1)
withoutalgo <- rbind(withoutalgo,c("09_CPX",0)); withoutalgo <- rbind(withoutalgo,c("12_BF",0)); withoutalgo <- rbind(withoutalgo,c("13_CPX",0)); withoutalgo <- rbind(withoutalgo,c("C",0))
withoutalgo <-withoutalgo[order(withoutalgo$Var1),]
withoutalgo$Freq <- as.integer(withoutalgo$Freq)
### diff:
withalgo<- withalgo %>%
  mutate(diff=case_when(
    withalgo$Freq>withoutalgo$Freq  ~ 1,
    withalgo$Freq<withoutalgo$Freq  ~ -1,
    withalgo$Freq==withoutalgo$Freq  ~ 0,
  )) #11 168


ggplot(withoutalgo, aes(x=Var1, y=Freq))+
  geom_bar(stat='identity', position='identity') +
  theme_bw() +
  labs( x='Subtype', y='Count', title ="Without applied Subtype Classification Algorithm" ) +
  theme(axis.text.x=element_text(angle=-90, vjust=0.5))
  
ggplot(withalgo, aes(x=Var1, y=Freq, fill=factor(diff)))+
  geom_bar(stat='identity', position='identity') +
  scale_fill_manual(values=c("darkviolet","grey70","darkgreen"), 
                    labels=c("Decrease","No Change","Increase"),
                    name="") +
  theme_bw() +
  labs( x='Subtype', y='Count',title = "With applied Subtype Classification Algorithm") +
  theme(axis.text.x=element_text(angle=-90, vjust=0.5))
                          

##not linked:
withalgo1 <-as.data.frame(table(finaldataset[finaldataset$id %in% rulesWOlink,"subtype"]))
withoutalgo1 <-as.data.frame(table(test_info[test_info$id %in% rulesWOlink,"subtype"]))
withoutalgo1$Var1 <- toupper(withoutalgo1$Var1)
### diff:
withalgo1<- withalgo1 %>%
  mutate(diff=case_when(
    withalgo1$Freq>withoutalgo1$Freq  ~ 1,
    withalgo1$Freq<withoutalgo1$Freq  ~ -1,
    withalgo1$Freq==withoutalgo1$Freq  ~ 0,
  )) #11 168


ggplot(withoutalgo1, aes(x=Var1, y=Freq))+
  geom_bar(stat='identity', position='identity') +
  theme_bw() +
  labs( x='Subtype', y='Count', title ="Without applied Subtype Classification Algorithm" ) +
  theme(axis.text.x=element_text(angle=-90, vjust=0.5))

ggplot(withalgo1, aes(x=Var1, y=Freq, fill=factor(diff)))+
  geom_bar(stat='identity', position='identity') +
  scale_fill_manual(values=c("darkviolet","grey70","darkgreen"), 
                    labels=c("Decrease","No Change","Increase"),
                    name="") +
  theme_bw() +
  labs( x='Subtype', y='Count',title = "With applied Subtype Classification Algorithm") +
  theme(axis.text.x=element_text(angle=-90, vjust=0.5))



withalgoALL <-as.data.frame(table(finaldatasetTEST$subtype))
withoutalgoALL <-as.data.frame(table(test_info$subtype))
withoutalgoALL$Var1 <- toupper(withoutalgoALL$Var1)
### diff:
withalgoALL<- withalgoALL %>%
  mutate(diff=case_when(
    withalgoALL$Freq>withoutalgoALL$Freq  ~ 1,
    withalgoALL$Freq<withoutalgoALL$Freq  ~ -1,
    withalgoALL$Freq==withoutalgoALL$Freq  ~ 0,
  )) #11 168

# withalgoALL[withalgoALL$Var1 == "UNDETERMINED (SHORT SEQUENCE)","Var1"] <- "SHORT"
# withalgoALL$Var1[c("UNDETERMINED (SHORT SEQUENCE)")] <- "SHORT"


ggplot(withoutalgoALL, aes(x=Var1, y=Freq))+
  geom_bar(stat='identity', position='identity') +
  theme_bw() +
  labs( x='Subtype', y='Count', title ="Without applied Subtype Classification Algorithm" ) +
  # theme(axis.text.x=element_text(angle=-90, vjust=0.5))+
  coord_cartesian(ylim= c(0, 800))+
  theme(axis.text=element_text(size=14),
        axis.text.x=element_text(angle=-90, vjust=0.5),
        axis.title=element_text(size=15,face="bold"),
        axis.title.y.right =element_text(colour ="blue"))

ggplot(withoutalgoALL, aes(x=Var1, y=Freq))+
  geom_bar(stat='identity', position='identity') +
  theme_bw() +
  labs( x='Subtype', y='Count', title ="Without applied Subtype Classification Algorithm" ) +
  theme(axis.text.x=element_text(angle=-90, vjust=0.5))+
  coord_cartesian(ylim= c(1000, 9500))+
  theme(axis.text=element_text(size=14),
        axis.text.x=element_text(angle=-90, vjust=0.5),
        axis.title=element_text(size=15,face="bold"),
        axis.title.y.right =element_text(colour ="blue"))

ggplot(withalgoALL, aes(x=Var1, y=Freq, fill=factor(diff)))+
  geom_bar(stat='identity', position='identity') +
  scale_fill_manual(values=c("darkviolet","grey70","darkgreen"), 
                    labels=c("Decrease","No Change","Increase"),
                    name="") +
  theme_bw() +
  labs( x='Subtype', y='Count',title = "With applied Subtype Classification Algorithm") +
  theme(axis.text.x=element_text(angle=-90, vjust=0.5))+
  coord_cartesian(ylim= c(0, 800))+
  theme(axis.text=element_text(size=14),
        axis.text.x=element_text(angle=-90, vjust=0.5),
        axis.title=element_text(size=15,face="bold"),
        axis.title.y.right =element_text(colour ="blue"))
  # ylim(0, 600)


#reset settings
par(op)







###SUBTYPE
ta<-table(finaldatasetTEST$subtype)
cross_table <-CrossTable(finaldatasetTEST$subtype)


mosaicplot(t(ta),col=TRUE, main="HIV Subtypes distribution within the SHCS")

### Output:
# ged rid of Subtype "unknown"!!
DD <- finaldataset[finaldataset$subtype != "UNKNOWN" & finaldataset$subtype != "UNDETERMINED (SHORT SEQUENCE)", ] #11,171


### Output:
DD<- DD %>%
  mutate(copy_subtype=case_when(
    DD$subtype=="B"  ~ "B",
    #~ "non-B"
  )) #11 168
select <- which(is.na(DD$copy_subtype))
DD$copy_subtype[select] <- "Non-B"

CrossTable(DD$subtype)



### LAB
ta4<-table(finaldataset$lab)
CrossTable(finaldataset$lab)
CrossTable(data$lab) #no big changes observable



n_labs <-data.frame(matrix(NA, ncol = 3, nrow = 4))
colnames(n_labs ) <- c("Lab","N","Percentage")
n_labs$Lab <- c("Basel",  "Geneva","Lausanne","Zurich")
n_labs$N <- c(t(ta4)[1],t(ta4)[2], t(ta4)[3], t(ta4)[4])
total <- sum(t(ta4))
n_labs$Percentage <- c(t(ta4)[1]/total, t(ta4)[2]/total, 
                       t(ta4)[3]/total, t(ta4)[4]/total)
n_labs$Percentage <- round(n_labs$Percentage,3)



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
summary(finaldataset$pol_len)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 858     887    1302    1180    1302    1965 

boxplot(finaldataset$pol_len,
        main = "Mean partial pol sequence length",
        xlab = "Length [nt]",
        ylab = "Pol",
        col = "orange",
        border = "brown",
        horizontal = TRUE,
        notch = TRUE)


ggplot(finaldataset, aes(x=lab, y=pol_len)) + 
  geom_boxplot(outlier.colour="black", outlier.shape=16,
               outlier.size=2, notch=FALSE, col="orange")+
  labs(title="Mean partial pol sequence length by lab",y="Length [nt]", x = "Lab")+
  coord_flip()+
  theme_bw() + 
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=15,face="bold"),
        plot.title = element_text(size=20,face="bold.italic"))

aggregate(finaldataset$pol_len, list(finaldataset$lab), mean)
# Group.1        x
# 1    Basel 1279.808
# 2   Geneva  889.205
# 3 Lausanne 1455.252
# 4   Zurich 1302.744

