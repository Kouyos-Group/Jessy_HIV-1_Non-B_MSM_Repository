# non-B subtypes 
# June 2020, Jessy Duran Ramirez


### Clear R's memory
rm(list=ls())

### get the packages
if (!requireNamespace("pacman", quietly = TRUE))
  install.packages("pacman")
pacman::p_load("dplyr","haven","tidyr", "data.table")


## IMPORT DATASET
setwd("~/Desktop/SHCS")
test_info <- as.data.table(read_dta("Input/test_info0620.dta")) #24,564
pat <- as.data.table(read_dta("Input/pat0520.dta")) #20,904
#admin <- read_dta("admi.dta") #20,802




### Quality check:
# partial pol ( #PR > 250nt and RT > 500 nt):
sum(is.na(test_info$pr_len)) #553
sum(is.na(test_info$rt_len)) #546
sum(is.na(test_info$int_len)) #21,139
sum(test_info$pr_len < 250) #3
sum(test_info$rt_len < 500) #1
sum(test_info$subtype == "HIV2 subtype A") #3
sum(test_info$subtype  == "") #0


# get rid of NA's and too short RTs & PRs:
test_info <- test_info[!is.na(test_info$pr_len),] #24,259
test_info <- test_info[!is.na(test_info$rt_len),] #23,706
test_info <- test_info[!test_info$pr_len < 250,] # 23,703
test_info <- test_info[!test_info$rt_len < 500,] #23,702


## get rid of HIV-2 and empty entries in the Subtype
test_info <- test_info[!test_info$subtype  == "",] #23,702
test_info <- test_info[!test_info$subtype == "HIV2 subtype A",] #23,699


#Include patients, which are in pat.dta & test_info.dta
sum(duplicated(test_info$header_id)) #0
sum(!test_info$id %in% pat$id) #13
test_info <- test_info[test_info$id %in% pat$id,] #23,686



#order by id, sampledate, seq_dt
#But first make sure that dates are dates & change all the subtypes toupper()
test_info$sampledate <- as.Date(as.character(test_info$sampledate), tryFormats = c("%Y-%m-%d"))
test_info$seq_dt <- as.Date(as.character(test_info$seq_dt), tryFormats = c("%Y-%m-%d"))
test_info$subtype <- toupper(test_info$subtype)
test_info$rega_subtype <- toupper(test_info$rega_subtype)
test_info$comet_subtype <- toupper(test_info$comet_subtype)
test_info <- test_info[order(test_info$id, test_info$sampledate,test_info$seq_dt),]


# # get variables list:
# # names(test_info)
#  varlist <- c("header_id", "lab", "sex", "id", "sampledate", "seq_dt", "rega_subtype", "comet_subtype", "subtype")
# # varlist2 <- c( "id", "sampledate", "seq_dt", "rega_subtype", "comet_subtype", "subtype")
#  sub_info <- test_info[varlist]

# 
# data <- test_info %>%
#   select(id, sampledate, seq_dt, rega_subtype, comet_subtype, subtype,  lab, linked) %>% #include the lab!
#   group_by(id) #12945 ordered by id  #data <- data[order(data$id),] 
# 

data <- test_info

### Numbering the sequencies per unique id:
data$nseq<- unlist(tapply(data$id, data$id,
                  function(x) seq(1,length(x),1)))

## Total Number of sequencies per unique id:
total_nseq<-cbind(unique(data$id),total_nseq=sapply(X=unique(data$id), function(x) sum(data$id==x)))
colnames(total_nseq)<-c("id","total_nseq")

data <- merge(data, total_nseq, by="id") ## add the total_nseq per unique id to the table data



### Number of Subtypes in  subtypes, Rega and Comet:
unique(test_info$subtype)
length(unique(test_info$subtype)) #39
n_subtype<-cbind(unique(test_info$subtype),n_subtype=sapply(X=unique(test_info$subtype), function(x) sum(test_info$subtype==x)))
colnames(n_subtype)<-c("subtype","n_subtype"); n_subtype <- as.data.frame(n_subtype)


unique(test_info$rega_subtype)
length(unique(test_info$rega_subtype)) #142
n_rega_subtype<-cbind(unique(test_info$rega_subtype),n_rega_subtype=sapply(X=unique(test_info$rega_subtype), function(x) sum(test_info$rega_subtype==x)))
colnames(n_rega_subtype)<-c("rega_subtype","n_rega_subtype"); n_rega_subtype <- as.data.frame(n_rega_subtype)


unique(test_info$comet_subtype)
length(unique(test_info$comet_subtype)) #76
n_comet_subtype<-cbind(unique(test_info$comet_subtype),n_comet_subtype=sapply(X=unique(test_info$comet_subtype), function(x) sum(test_info$comet_subtype==x)))
colnames(n_comet_subtype)<-c("comet_subtype","n_comet_subtype"); n_comet_subtype <- as.data.frame(n_comet_subtype)



### Number of UNDETERMINED (SHORT SEQUENCE): 
undef <-data[data$subtype=="UNDETERMINED (SHORT SEQUENCE)",] #toal 48 undertermined sequences, from 32 patients
#length(unique(undef$id)) #32
undef_data <- data[data$id %in% undef$id,] #72



### Search and Delete total_nseq==1 &&  having undertermined sequences:
# =====>  "HOW TO DO THIS MORE ELEGANT?? "
x<-0
for (i in 1:nrow(undef_data)){
  #new_undef_data[-i,] <- if(new_undef_data$total_nseq[i]==1 & new_undef_data$subtype == "UNDETERMINED (SHORT SEQUENCE)")}
  if (undef_data$total_nseq[i]==1 && undef_data$subtype[i] == "UNDETERMINED (SHORT SEQUENCE)") {
  x <- c(x,i) 
  }
}
new_undef_data<- undef_data[-x,] #56
length(unique(new_undef_data$id)) #16, so 18 patients were deleted







## rega and subtype different: 1049
subset_diff<- data[data$subtype!=data$rega_subtype,]#1,111 cases changed to 1045, after setting all toupper()
subset_diff_without_recombi <- data[data$subtype!=data$rega_subtype & data$subtype!="RECOMBINANT",] #171


### keeping the earlist per patient:
#first <- mydata[!duplicated(mydata$id),] #gives the same! 
dat_1 <- data %>%
  filter(nseq == 1)  #12,974

length(unique(dat_1$id)) #12,974

## >1 sampledate per patient
#dat_2 <- mydata[duplicated(mydata$id),] #gives the same
dat_2 <- data %>%
  filter(nseq >= 2)  ##10 712

length(unique(dat_2$id)) #4785

dat_3 <- data %>%
  filter(nseq <=2)  ##17 759


#which(mydata$subtype!= dplyr::lag(mydata$subtype))
changes <- data[which(data$subtype!= dplyr::lag(data$subtype)),] #6 008


same <- data[which(data$subtype!= dplyr::lag(data$subtype) & data$id== dplyr::lag(data$id)),] #427


#create a vector to store id with changing subtype
select <- rep(NA, nrow(data))
#selecting same ID and different subtypes, if fulfilled the statement, it is true and we output the ID. Otherwise it saves NA
#save it in select
for (i in 2:nrow(data)){
  select[i] <- ifelse(data$id[i]==data$id[i-1] & data$subtype[i]!=data$subtype[i-1],data$id[i],NA)
}
select <- unique(na.omit(select)) #306
newdata <- data[data$id %in% select,] #1 108 obs.



##### SaME if n == 2
#create a vector to store id with changing subtype
select2 <- rep(NA, nrow(dat_3))
#selecting same ID and different subtypes, if fulfilled the statement, it is true and we output the ID. Otherwise it saves NA
#save it in select
for (i in 2:nrow(dat_3)){
  select2[i] <- ifelse(dat_3$id[i]==dat_3$id[i-1] & dat_3$subtype[i]!=dat_3$subtype[i-1],dat_3$id[i],NA)
}
select2 <- unique(na.omit(select2)) #226
newdata_2cases <- dat_3[dat_3$id %in% select2,] #452


# ### Save the dataset:
# write.csv(select,"ids_for_rules.csv")
# write.csv(newdata,"subset_for_rules.csv")





#### WHAT RULE??
newdata[newdata$id==13898,]
newdata[newdata$id==	15248,]
newdata[newdata$id==	15716,]


newdata[newdata$id==	16099,]
newdata[newdata$id==16727	,]


newdata[newdata$id==	16742,]


newdata[newdata$id==16843	,]
newdata[newdata$id==17123	,]
newdata[newdata$id==	18278,]
newdata[newdata$id==18625	,]
newdata[newdata$id==18889	,]

# earliest?
newdata[newdata$id==18995	,]
newdata[newdata$id==19216	,]
newdata[newdata$id==19461	,]
### non-B and B: 2 cases

newdata[newdata$id==19284	,]
newdata[newdata$id==25259	,]
# Delete:

newdata[newdata$id==10555	,]
newdata[newdata$id==	18803,]

#### looked untill row 240

newdata[newdata$id==	16843,]
