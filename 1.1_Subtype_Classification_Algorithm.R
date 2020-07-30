# Algorithm for Non-B and B Subtypes Classification
# Jun 2020, Jessy Duran Ramirez


### Clear R's memory
rm(list=ls())


### GET THE PACKAGES AND FUNCTIONS
if (!requireNamespace("pacman", quietly = TRUE))
  install.packages("pacman")
pacman::p_load("dplyr","haven","tidyr", "data.table","GGally", "gmodels")
`%notin%` <- Negate(`%in%`)


### IMPORT DATASET
### Load the SHCS sequences, the test_info file and the pat.dta:
setwd("~/Desktop/SHCS")
SHCS_SEQ <- read.fasta(file="Input/final0620.fas") #29 837 
TEST_INFO <- read.dta13("Input/test_info0620.dta") #24 564
PAT <- read.dta13("Input/pat0520.dta") #20 904
ids_for_rules<- read.csv("Output/ids_for_rules.csv")
admin <- as.data.table(read_dta("Input/admi0520.dta")) #20,802
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


### A) QUALITY CHECK:
# partial pol ( #PR > 250nt and RT > 500 nt):
sum(is.na(test_info$pr_len)) #552
sum(is.na(test_info$rt_len)) #545
sum(is.na(test_info$int_len)) #21,140
length(which(test_info$pr_len < 250)) #6
length(which(test_info$rt_len < 500)) #1
sum(test_info$subtype == "HIV2 subtype A") #3
sum(test_info$subtype  == "") #1
# sum(test_info$subtype == "HIV 1 N group") #1
# sum(test_info$subtype == "HIV 1 O group") #5


## get rid of NA's and too short RTs & PRs:
test_info <- test_info[!is.na(test_info$pr_len),] #23,948
test_info <- test_info[!is.na(test_info$rt_len),] #23,643
test_info <- test_info[!test_info$pr_len < 250,] # 23,640
test_info <- test_info[!test_info$rt_len < 500,] #23,639


## get rid of HIV-2 and empty entries in the Subtype
test_info <- test_info[!test_info$subtype  == "",] #23,638
test_info <- test_info[!test_info$subtype == "HIV2 subtype A",] #23,635
# test_info <- test_info[!test_info$subtype == "HIV 1 N group",] #23,634
# test_info <- test_info[!test_info$subtype == "HIV 1 O group",] #23,629


## Exclude patients, which are  NOT in pat.dta & test_info.dta
sum(duplicated(test_info$header_id)) #0
sum(!test_info$id %in% pat$id) #13
test_info <- test_info[test_info$id %in% pat$id,] #20,622

test_info<- test_info %>%
  mutate(pol_len = pr_len+rt_len) ## add the total_nseq per unique id to the table data


### From here part B) SUBTYPE CLASSIFICATION ALGORITHM:

## 1) If retrospectively linked get rid of them: 
#linked = 1 
length(which(test_info$linked=="1")) #3019
test_info <- test_info[is.na(test_info$linked),] #20,603


## order by id, sampledate, seq_dt
#But first make sure that dates are dates & change all the subtypes toupper()
test_info$sampledate <- as.Date(as.character(test_info$sampledate), tryFormats = c("%Y-%m-%d"))
test_info$seq_dt <- as.Date(as.character(test_info$seq_dt), tryFormats = c("%Y-%m-%d"))
test_info$subtype <- toupper(test_info$subtype)
test_info$rega_subtype <- toupper(test_info$rega_subtype)
test_info$comet_subtype <- toupper(test_info$comet_subtype)
test_info <- test_info[order(test_info$id, test_info$sampledate,test_info$seq_dt),]


## add Sandra's Analysis cat & patients that requieres a rule
test_info <- merge(test_info,ids_inSI_325potential,all.x=T,by="id")
test_info <- test_info %>%
  mutate(rule=case_when(
    id %in% ids_for_rules$x  ~ "1",
    id %notin% ids_for_rules$x ~ "0",
  )) #20,603


## continue to work with the subset:
data <- test_info %>%
  select(id, sampledate, seq_dt, rega_subtype, comet_subtype, subtype,  lab, pr_len, rt_len,cat, rule) %>% #include the lab!
  group_by(id) #12603 obs. ordered by id  #data <- data[order(data$id),] 


# Getting the total sequence length:
data <- data %>%
  mutate(pol_len = pr_len+rt_len) ## add the total_nseq per unique id to the table data

# Numbering the sequencies per unique id:
data$nseq<- unlist(tapply(data$id, data$id,
                          function(x) seq(1,length(x),1)))


# Total Number of sequencies per unique id:
total_nseq<-cbind(unique(data$id),total_nseq=sapply(X=unique(data$id), function(x) sum(data$id==x)))
colnames(total_nseq)<-c("id","total_nseq")

data <- data %>%
  mutate(total_nseq = max(nseq)) ## add the total_nseq per unique id to the table data

###
# below_threshold <- data[!data$pol_len > 886,] #83ob.
# length(unique(below_threshold$id)) #74
# check <- below_threshold[below_threshold$id %in% ids_for_rules$x,] #13ob.


### 2) find the longest sequence per patient and those > threshold => LONGEST OR ABOVE THRESHOLD
d1 <- merge(aggregate(pol_len ~ id, max, data = data), data) #15988 longest obs.
# below_thresholdd <- d1[d1$pol_len < 887,] #11ob.
# length(unique(below_thresholdd$id)) #11
summary(data$pol_len)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 858     887    1302    1167    1302    1977 

threshold<- 887 #threshold = 1st Qu.
d2 <- data[sapply(data$pol_len, 
                     function(x) any(x>threshold) )==T, ] #13539 above the threshold obs.
#below_thresholddd <- d2[d2$pol_len < 887,] #0ob.
data <-full_join(d2,d1) #18769 obs. 
length(unique(data$id)) #11216
# dd2<- data[data$id %in% ids_for_rules$x,] #812 obs.
# length(unique(dd2$id)) #296
below_threshold2 <- data[data$pol_len < 887,] #11ob.
length(unique(below_threshold2$id)) #11

### 3) find the earlieast among the longest & above the treshhold => EARLIEST SAMPLE DATE
# take only the first sequence of the patients:
data <- merge(aggregate(sampledate ~ id, min, data = data), data) #11291
length(unique(data$id)) #11216, meaning that 75 patients have multiple tests (n=2)


dd2<- data[data$id %in% ids_for_rules$x,] #301 obs.
length(unique(dd2$id)) #296

 #==> ONLY 2 tests per id
# Numbering the sequencies per unique id:
 data$nseq2<- unlist(tapply(data$id, data$id,
                            function(x) seq(1,length(x),1)))

dat_2 <- data%>%
    filter(nseq2 >= 2)  ##75

info<- data[data$id %in% dat_2$id,] #150 id's with multiple tests
 
 
 #==>2 CHANGING SUBTYPES!
 select <- rep(NA, nrow(data))
#selecting same ID and different subtypes, if fulfilled the statement, it is true and we output the ID.
#Otherwise it saves NA save it in select
for (i in 2:nrow(data)){
  select[i] <- ifelse(data$id[i]==data$id[i-1] & data$subtype[i]!=data$subtype[i-1],data$id[i],NA)
}
select <- unique(na.omit(select)) #2 ids: 31102 and 52425
newdata <- data[data$id %in% select,] #4 obs. 



### 3)  if have multiple tests, take most frequent subtype => EARLIEST SEQUENCE DATE
# have only 2 sequences, taking the most frequent one difficult
data <- merge(aggregate(seq_dt ~ id, min, data = data), data) #11218!! 
sum(duplicated(data$id)) # have 2 cases where have duplicates ##ids: 90289 and 57333
data<- data[!duplicated(data$id),]#11216 ## ==> FINAL, SAVE THIS
data <- data[order(data$id, data$sampledate),]

####--------------------

### Save the dataset:
#write.csv(data,"Output/table_subtypes.csv") 

####--------------------

newdata2 <- data[data$id %in% ids_for_rules$x,]
newdata3 <- pat[pat$id %in% data$id ,]
sel <-c(which(data$subtype=="UNKNOWN"), which(data$subtype=="UNDETERMINED (SHORT SEQUENCE)"))
newdata4 <- data[data$id %notin% sel,] #11214

