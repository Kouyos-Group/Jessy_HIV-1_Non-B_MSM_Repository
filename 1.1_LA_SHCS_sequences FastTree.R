# June 2020, JDR
# Part B - Phylogeny 

### Make Tree
##########################################

## clear R's brain
rm(list=ls())


## packages
library(dplyr)
library(readstata13)
library(ape) #Analyses of Phylogenetics and Evolution
library(seqinr) #Biological Sequences Retrieval and Analysis
library(ggtree)



### Load the SHCS sequences, the test_info file and the pat.dta:
setwd("~/Desktop/SHCS/Input")
SHCS_SEQ <- read.fasta(file="final0620.fas") #29 837 
TEST_INFO <- read.dta13("test_info0620.dta") #24 564
PAT <- read.dta13("pat0520.dta") #20 904


#################################### 
length(unique(TEST_INFO$id))# 13278


### Subtype Classification Algorithm,  quality check (cleaning/filtering SHCS sequences) and one sequence per patient:

Sub_Classification_SHCS <- function(SHCS_SEQ,min_pr,min_rt,TEST_INFO,PAT, subtype){
  
  #rename the sequences to only header id
  names(SHCS_SEQ)<-sapply(names(SHCS_SEQ),
                          FUN=function(a){unlist(strsplit(a,split="\\|"))[1]})
  
  ### Quality Check:
  #get rid of NA's and too short RTs & PRs:
  TEST_INFO <- TEST_INFO[!is.na(TEST_INFO$pr_len),] #24,259
  TEST_INFO <- TEST_INFO[!is.na(TEST_INFO$rt_len),] #23,706
  TEST_INFO <- TEST_INFO[!TEST_INFO$pr_len < min_pr,] # 23,703
  TEST_INFO <- TEST_INFO[!TEST_INFO$rt_len < min_rt,] #23,702
  
  #check if they are in TEST_INFO and PAT tables
  TEST_INFO <- TEST_INFO[TEST_INFO$id %in% PAT$id,] #23 689
  #length(unique(TEST_INFO$id))# 12 977
  

  ### SUBTYPE CLASSIFICATION ALGORITHM:
  ##1) Exclude If retrospectively linked: 
  TEST_INFO  <- TEST_INFO [is.na(TEST_INFO$linked),] #20,667
  
  #order by id, sampledate, seq_dt
  #But first make sure that dates are dates & change all the subtypes toupper()
  TEST_INFO$sampledate <- as.Date(as.character(TEST_INFO$sampledate), tryFormats = c("%Y-%m-%d"))
  TEST_INFO$seq_dt <- as.Date(as.character(TEST_INFO$seq_dt), tryFormats = c("%Y-%m-%d"))
  TEST_INFO$subtype <- toupper(TEST_INFO$subtype)
  TEST_INFO$rega_subtype <- toupper(TEST_INFO$rega_subtype)
  TEST_INFO$comet_subtype <- toupper(TEST_INFO$comet_subtype)
  TEST_INFO <- TEST_INFO[order(TEST_INFO$id, TEST_INFO$sampledate,TEST_INFO$seq_dt),]
  
  
  ### Adding some variables:
  # Total Partial Pol Length:
  TEST_INFO<- TEST_INFO %>%
    mutate(pol_len = pr_len+rt_len) ## add the total_nseq per unique id to TEST_INFO
  # Numbering the sequencies per unique id:
  TEST_INFO$nseq<- unlist(tapply(TEST_INFO$id, TEST_INFO$id,
                              function(x) seq(1,length(x),1)))
  # Total Number of sequencies per unique id:
  total_nseq<-cbind(unique(TEST_INFO$id),total_nseq=sapply(X=unique(TEST_INFO$id), function(x) sum(TEST_INFO$id==x)));colnames(total_nseq)<-c("id","total_nseq")
  TEST_INFO <- merge(TEST_INFO, total_nseq, by="id") ## add the total_nseq per unique id to the table data
  
  
  
  ##2) find the longest sequence per patient and those > threshold => LONGEST OR ABOVE THRESHOLD
  LONGEST  <- merge(aggregate(pol_len ~ id, max, data = TEST_INFO), TEST_INFO) #15 995
  threshold <- summary(TEST_INFO$pol_len)[[2]] #threshold = 1st Qu., 887
  ABOVE <- TEST_INFO[sapply(TEST_INFO$pol_len, 
                    function(x) any(x>threshold) )==T, ] #13 604 above the threshold obs
  TEST_INFO <-full_join(ABOVE,LONGEST) #18 821
  
  ### 3) find the earlieast among the longest & above the treshhold => EARLIEST SAMPLE DATE
  TEST_INFO <- merge(aggregate(sampledate ~ id, min, data = TEST_INFO), TEST_INFO) #11 324
  #length(unique(TEST_INFO$id)) #11248
  
  
  ### 4)  if have multiple tests (n=2), take most frequent subtype => EARLIEST SEQUENCE DATE
  TEST_INFO <- merge(aggregate(seq_dt ~ id, min, data = TEST_INFO), TEST_INFO) #11 251
  #sum(duplicated(data$id)) # have 3 cases where have duplicates ##ids: 90289, 57333 91477
  
  ### 5)  if have duplicates, remove them
  TEST_INFO <- TEST_INFO[!duplicated(TEST_INFO$id),]#11248  
  #length(unique(TEST_INFO$id)) #11248  
  
  
  ###
  #NONB_INFO <- TEST_INFO[TEST_INFO$subtype != "B",] #2940
  #NONB_INFO <- TEST_INFO[TEST_INFO$subtype %in% c("A","C","F","01_AE", "02_AG"),] #1986
  SUBTYPE_INFO <- TEST_INFO[TEST_INFO$subtype %in% subtype,] #A: 492
  
  
  # SHCS headers from NONB_INFO
  keep.header <- SUBTYPE_INFO$header_id
  SHCS_SEQ<-SHCS_SEQ[names(SHCS_SEQ) %in% keep.header]

  
  return(SHCS_SEQ)
}


# clean_seq_A <- Sub_Classification_SHCS(SHCS_SEQ = SHCS_SEQ ,min_pr = 250, min_rt = 500,TEST_INFO = TEST_INFO,PAT = PAT, subtype =  "A")
# write.dna(clean_seq_A,"clean_seq.fas",format = "fasta")
## -> 492 sequences

# clean_seq_C <- Sub_Classification_SHCS(SHCS_SEQ = SHCS_SEQ ,min_pr = 250, min_rt = 500,TEST_INFO = TEST_INFO,PAT = PAT, subtype =  "C")
# write.dna(clean_seq_C,"clean_seq_C.fas",format = "fasta")
# ## ->  sequences
# 
# clean_seq_AE <- Sub_Classification_SHCS(SHCS_SEQ = SHCS_SEQ ,min_pr = 250, min_rt = 500,TEST_INFO = TEST_INFO,PAT = PAT, subtype =  "01_AE")
# write.dna(clean_seq_AE,"clean_seq_AE.fas",format = "fasta")
# ## ->  sequences
# 
# clean_seq_AG <- Sub_Classification_SHCS(SHCS_SEQ = SHCS_SEQ ,min_pr = 250, min_rt = 500,TEST_INFO = TEST_INFO,PAT = PAT, subtype =  "02_AG")
# write.dna(clean_seq_AG,"clean_seq_AG.fas",format = "fasta")
# ## ->  sequences
# 
c(clean_seq_F, F_info) <- Sub_Classification_SHCS(SHCS_SEQ = SHCS_SEQ ,min_pr = 250, min_rt = 500,TEST_INFO = TEST_INFO,PAT = PAT, subtype =  "F")
### Save the dataset:
write.csv(F_info,"F_info.csv")
write.dna(clean_seq_F,"clean_seq.fas",format = "fasta")

## ->  127 sequences

#howMany <- merge(TEST_INFO[,c('header_id','id','sampledate','subtype')],PAT) #11 248
#howMany.MSM <- howMany[howMany$risk == 1,] #4 368
# howMany.A <- howMany[howMany$subtype == "A",] #A: 492
# howMany.C <- howMany[howMany$subtype == "C",] #C: 448
# howMany.AE <- howMany[howMany$subtype == "01_AE",] #AE: 431
# howMany.AG <- howMany[howMany$subtype == "02_AG",] #AG: 488
# howMany.F <- howMany[howMany$subtype == "F",] #A: 127
# table(TEST_INFO$subtype)
#table(howMany.MSM$subtype) #AE: 137, AG: 44, C:28, A:45, F:45


#################################### 
# use LA sequences
#myLosAlamosSequences <- read.fasta(file="LA_seq0620.fasta") #28796

## 89841 sequences in LAallSTuntil 2007.fasta
myLosAlamosSequences1 <- read.fasta(file="LAallSTuntil 2007.fasta")
##  148774 sequences in LAallST2008to2020.fasta
myLosAlamosSequences2 <- read.fasta(file="LAallST2008to2017.fasta")
## -> 239,535 background sequences


# quality check for the Los Alamos background sequences
qual_check_LA <- function(LA_file, removeSwiss){
  ## Input: 
  #LA_file: The file containing the background sequences
  #removeSwiss: boolean variable to indicate whether Swiss sequences should be removed: 0 = no, 1 = yes.
  #selectSubtyp: boolean variable to indicate whether only want specifc subtpyes: 0 = no, 1 = yes.
  # subtype: select the subtype wanted, e.g. "B"
  #----------------------------------------------------------------
  ## Output:
  # Cleaned file containing the background sequences
  
  if (removeSwiss == 1){ # removing the Swiss sequences if 1
    LA_file <- LA_file[!sapply(names(LA_file),FUN = function(a){unlist(strsplit(a,split = "[.]"))[2]}) == "CH"]
  }
  # 
  # ##select for subtype
  # if (selectSubtype == 1){ # if 1, specific subtypes will be selected
  #   LA_file <- LA_file[sapply(names(LA_file),FUN = function(a){unlist(strsplit(a,split = "[.]"))[1]}) == subtype]
  # }
  
  #rename to only accession number/ [5] = last number ("A07867")
  names(LA_file)<-sapply(names(LA_file),FUN=function(a){unlist(strsplit(a,split="[.]"))[5]})
  LA_file <- LA_file[!duplicated(names(LA_file))] #remove duplicates

  return(LA_file)
}




## quality check LA:
# clean_LA1 <- qual_check_LA(myLosAlamosSequences1,1)
# clean_LA2 <- qual_check_LA(myLosAlamosSequences2,1)
# clean_LA <- append(clean_LA1, clean_LA2)
# write.fasta(clean_LA,names(clean_LA),"clean_LA.fasta")



# clean_LA1_A <- qual_check_LA(myLosAlamosSequences1,1,1,c("A","A1","A2","A3","A4","A5","A6")) #494
# clean_LA2_A <- qual_check_LA(myLosAlamosSequences2,1,1,c("A","A1","A2","A3","A4","A5","A6")) #1567
# clean_LA_A <- append(clean_LA1_A, clean_LA2_A) #2061
# write.fasta(clean_LA_A,names(clean_LA_A),"clean_LA_A.fasta")

# clean_LA1_C <- qual_check_LA(myLosAlamosSequences1,1,1,"C") #4130
# clean_LA2_C <- qual_check_LA(myLosAlamosSequences2,1,1,"C") #4130
# clean_LA_C <- append(clean_LA1_C, clean_LA2_C) #
# write.fasta(clean_LA_C,names(clean_LA_C),"clean_LA_C.fasta")
# 
# 
# clean_LA1_AE <- qual_check_LA(myLosAlamosSequences1,1,1,"01_AE") #3547
# clean_LA2_AE <- qual_check_LA(myLosAlamosSequences2,1,1,"01_AE") #3547
# clean_LA_AE <- append(clean_LA1_AE, clean_LA2_AE) #
# write.fasta(clean_LA_AE,names(clean_LA_AE),"clean_LA_AE.fasta")
# 
# 
# clean_LA1_AG <- qual_check_LA(myLosAlamosSequences1,1,1,"02_AG") #270
# clean_LA2_AG <- qual_check_LA(myLosAlamosSequences2,1,1,"02_AG") #270
# clean_LA_AG <- append(clean_LA1_AG, clean_LA2_AG) #
# write.fasta(clean_LA_AG,names(clean_LA_AG),"clean_LA_AG.fasta")
# 
# 


clean_LA1_F <- qual_check_LA(myLosAlamosSequences1,0,1,c("F","F1","F2")) #1
clean_LA2_F <- qual_check_LA(myLosAlamosSequences2,0,1,c("F","F1","F2")) #1
clean_LA_F <- append(clean_LA1_F, clean_LA2_F) #
write.fasta(clean_LA_F,names(clean_LA_F),"clean_LA.fasta")


clean_LA1_F <- qual_check_LA(myLosAlamosSequences1,1) #1
clean_LA2_F <- qual_check_LA(myLosAlamosSequences2,1) #1
clean_LA_F <- append(clean_LA1_F, clean_LA2_F) #
write.fasta(clean_LA_F,names(clean_LA_F),"clean_LA.fasta")


#################################### 
## BLAST
# Basic Local Alignment Search Tool: finds regions of similarity between biological sequences. 
#The program compares nucleotide or protein	sequences to sequence databases and calculates the statistical significance.

####################################
# makeblastdb
makedatabase <- function(fastadirectory, fastafilename, blastdirectory){
  ## Input:
  #fastadirectory The directory of the fasta file containing the Cleaned file with the background sequences
  #fastafilename: The name of the fasta file containing the Cleaned file with the background sequences
  #blastdirectory: The directory of the makeblastdb executable file
  #----------------------------------------------------------------
  ## Output:
  # produces BLAST databases from FASTA files named 'clean_LA' 
  stringcommand<-paste(blastdirectory ,"/makeblastdb -in ", 
                       fastadirectory, "/", fastafilename, " -input_type fasta -dbtype nucl -title clean_LA -out ", 
                       fastadirectory, "/","clean_LA", sep = "")
  system(stringcommand)
  #return(syscommand)
}

#system("/Users/jessyduran/Desktop/Blast/bin/makeblastdb -help") #working
directory <- getwd()
makedatabase(directory, "clean_LA.fasta", directory)
LA <- read.fasta("clean_LA.fasta")


# blasting
blasting <- function(seqdirectory, fastaname, blastdirectory, maxtar, percid, dbname){
  ## Input:
  # seqdirectory: The path to the place where the SHCS sequences are stored
  # fastaname: The name with the file containing the (cleaned) SHCS sequences
  # blastdirectory: The directory of the blastn executable file
  # maxtar: An integer indicating the maximum number of hits per sequence
  # percid: An integer between 0 and 100 indicating the percentage of identity necessary to include a sequence
  # dbname: The name of the database which will be blasted
  #----------------------------------------------------------------
  ## Output:
  # The SHCS sequences are BLASTed against the background sequences
  # A new fasta file 'sequencesplushitsbutunaligned.fas' is created
  fastafulldir<-paste(seqdirectory,"/", fastaname, sep = "")
  jkl<-read.FASTA(fastafulldir)
  mno<-del.gaps(jkl)
  write.dna(mno, "sequencesplushitsbutunaligned.fas", "fasta")
  syscommand <- paste(blastdirectory, " -query ", seqdirectory, 
                      "/sequencesplushitsbutunaligned.fas -out ", seqdirectory, 
                      "/seq.txt -outfmt 6 -db ", seqdirectory, 
                      "/clean_LA -evalue 0.05 -max_target_seqs ", 
                      maxtar," -max_hsps 1 -perc_identity ", percid, sep = "")
  system(syscommand)
  #return(syscommand)
  dbfulldir<- paste(seqdirectory,"/", dbname, sep = "")
  #return(dbfulldir)
  database<-read.dna(dbfulldir, format = "fasta", as.matrix = FALSE)
  blasthitsdir<-paste(seqdirectory,"/", "seq.txt", sep = "") 
  blasthits<-read.table(blasthitsdir)
  blasthits<-as.matrix(blasthits)
  hitnames<-blasthits[,2]
  for (q in 1:length(hitnames)) {
    outputdir<- paste(seqdirectory,"/","sequencesplushitsbutunaligned.fas", sep = "") 
    write.dna(database[toString(hitnames[q])], outputdir, format = "fasta", append = TRUE)
  }
  duptemp<-read.FASTA(outputdir)
  duptemp<-duptemp[!duplicated(names(duptemp))]
  write.dna(duptemp, outputdir, format = "fasta", append = FALSE)
}

maxtar <- 10
percid <- 90
seqdirectory <- directory
blastdirectory <- paste(directory,"/blastn",sep = "")
fastaname <- "clean_seq.fas"
dbname <- "clean_LA.fasta"
blasting(seqdirectory, fastaname, blastdirectory, maxtar, percid, dbname)



#################################### 
## MUSCLE


#VA containing the fasta Reference sequence = HXB2 Pol Sequence from position 2253-3870
# HXB2_full <- read.GenBank("K03455.1")
# write.dna(HXB2_full, file ="HXB2_full.fas", format = "fasta", colsep = "")
# HXB2_full_trans <- read.dna(file = "HXB2_full.fas",as.matrix = TRUE,format = "fasta")
# write.dna(HXB2_full_trans[1,2253:3870], format = "fasta", "HXB2.fas") 
Reference <- read.fasta(file = "HXB2.fas")

# VA containing the fasta Sequences (SHCS + Background) output of the BLAST function
Sequences <- read.fasta(file = "sequencesplushitsbutunaligned.fas")

## Only first time:
# download.file("https://www.drive5.com/muscle/downloads3.8.31/muscle3.8.31_i86darwin32.tar.gz", "muscle3.8.31_i86darwin32.tar.gz")
# untar("muscle3.8.31_i86darwin32.tar.gz")
# Rename it to muscle! 




Pairwise_Alignment <- function(blastdirectory,removeFile,Ref, Seqs){
  ## Input:
  # blastdirectory:
  # removeFile
  # Ref: The fasta file conaining the reference sequence
  # Seqs: The fasta file containing the sequences which should be aligned,  i.e. the file 'sequencesplushitsbutunaligned.fas' created with 'blasting'
  #----------------------------------------------------------------
  ## Output:
  # A new fasta file 'Seqs_Aligned_gaps_removed.fas' with the aligned sequences is created.
  # The sequences which were blasted before are aligned against a reference sequence using the MUSCLE executable file
  # if (removeFile == 1){ # removing the existing files
  #   file.remove("Seqs_Aligned.fas")
  #   file.remove("Seqs_Aligned_gaps_removed.fas")
  # }
 
  for(i in 1:length(Seqs)){
    if (i %% 100 ==0){
      Sys.sleep(0.01)
      print(i)
    }
    write.fasta(sequences = c(Ref,Seqs[i]), names = names(c(Ref,Seqs[i])) , file.out = "To_be_Aligned.fasta")
    syscommand <- paste(blastdirectory, "/muscle -in To_be_Aligned.fasta -out Seqs_Aligned.afa", sep = "")
    system(syscommand,  intern=TRUE, wait=TRUE)
    #system("/Users/jessyduran/Desktop/SHCS/Input/muscle -in To_be_Aligned.fasta -out Seqs_Aligned.afa", intern=TRUE, wait=TRUE)
    Seqs_Aligned <- read.dna(file = "Seqs_Aligned.afa",as.character = TRUE,as.matrix = TRUE,format = "fasta")
    keepCol <- which(Seqs_Aligned[1,]!="-") ## check where the gaps in HXB2 are
    Seqs_Aligned_gaps_removed <- t(Seqs_Aligned[2,keepCol])# w
    Seqs_Aligned_gaps_removed <- as.DNAbin(Seqs_Aligned_gaps_removed)
    write.dna(Seqs_Aligned_gaps_removed,append = TRUE,"Seqs_Aligned_gaps_removed.fasta",format = "fasta")
  }
  Seqs_Aligned_gaps_removed <- as.list(read.dna("Seqs_Aligned_gaps_removed.fasta",format = "fasta"))
  names(Seqs_Aligned_gaps_removed) <- names(Seqs)
  write.dna(Seqs_Aligned_gaps_removed,"Seqs_Aligned_gaps_removed.fasta",format = "fasta")
}

blastdirectory_align <- directory
Pairwise_Alignment(blastdirectory_align,0,Reference, Sequences)


deleteMutationPositions <- function(fas_name){
  ## Input:
  # fas_name: The fasta file 'Seqs_Aligned_gaps_removed.fas'  with the aligned sequences
  #----------------------------------------------------------------
  ## Output:
  # A new fasta file 'shcs_seq,Seqs_Aligned_gaps_removed.fas' with the deleted mutation poisitions is created.
  
  tempseq <-read.dna(fas_name, format = "fasta", as.character = TRUE, as.matrix = TRUE)
  #deletions are done on columns of matrix of alignment  
  tempseq<-tempseq[, -c(89,90,91,95,96,97,98,99,100,137,138,139,140,141,142,143,144,145,149,
                        150,151,161,162,163,173,174,175,221,222,223,227,228,229,245,246,247,
                        248,249,250,251,252,253,263,264,265,269,270,271,419,420,421,482,483,
                        484,491,492,493,497,498,499,506,507,508,518,519,520,521,522,523,527,
                        528,529,596,597,598,599,600,601,605,606,607,614,615,616,620,621,622,
                        641,642,643,644,645,646,710,711,712,749,750,751,833,834,835,839,840,
                        841,848,849,850,860,861,862,866,867,868,926,927,928,941,942,943,953,
                        954,955,959,960,961,971,972,973,977,978,979,986,987,988)]   
  shcs_seq<-as.DNAbin(tempseq)
  write.dna(shcs_seq, fas_name,format="fasta",nbcol=ncol(shcs_seq))
}


deleteMutationPositions("Seqs_Aligned_gaps_removed.fasta")




# Trimming
## Only first time:
# download.file("http://trimal.cgenomics.org/_media/trimal.v1.2rev59.tar.gz", "trimal.v1.2rev59.tar.gz")
# untar("trimal.v1.2rev59.tar.gz")
# system("O:/Groups/Epi/duraj/Phylo/trimAl/source make sudo") #make trimal in trimAl/source
# in command line: make the package --> sudo make

## Only first time Windows:
# download.file("http://trimal.cgenomics.org/_media/trimal.v1.2rev59.zip", "trimal.v1.2rev59.zip")
# unzip("trimal.v1.2rev59.zip")
# trimal in bin folder!


trimming <- function(seqdirectory,fastaname,blastdirectory, GTvalue, CONSvalue){
  ## Input:
  # seqdirectory:
  # blastdirectory:
  # GTvalue: 0.7
  # CONSvalue: 0.5
  #----------------------------------------------------------------
  ## Output:
  # The sequences end are trimmed
  # A new fasta file 'trimmed.fas' is created
  
  commandstring<- paste(blastdirectory, "/trimal -in ", seqdirectory,"/",  fastaname, " -gt " , GTvalue,
                        " -cons ", CONSvalue,
                        " -out ", seqdirectory, "/trimmed.fas", sep = "")  
  system(commandstring)
}

#system("/Users/jessyduran/Desktop/SHCS/Input/trimAl/source/trimal -in Seqs_Aligned_gaps_removed.fasta -gt 0.7 -cons 0.5 -out trimmed.fas")

GTvalue <- 0.7
CONSvalue <- 0.5
seqdirectory <- directory
blastdirectory_trim <- paste(directory,"/trimAl/source", sep = "")
fastaname_trim <- "Seqs_Aligned_gaps_removed.fasta"
trimming(seqdirectory, fastaname_trim, blastdirectory_trim, GTvalue, CONSvalue)




#################################### 
## FastTree

## FastTree uses the Jukes-Cantor or generalized time-reversible (GTR) models of nucleotide evolution
# FastTree uses a single rate for each site (the "CAT" approximation). To quickly estimate the reliability 
# of each split in the tree, FastTree computes local support values with the Shimodaira-Hasegawa test (these are the same as PhyML 3's "SH-like local supports").



# Only first time: 
# go to http://www.microbesonline.org/fasttree/
# and install: http://www.microbesonline.org/fasttree/FastTree.c
# system("gcc -O3 -finline-functions -funroll-loops -Wall -o FastTree FastTree.c -lm")

buildTree <- function(blastdirectory,Seq_Name, Tree_Name){
  ## Input:
  # Seq_Name: The fasta file containing the sequences which should be trimmed i.e. the file 'trimmed.fas' created with 'trimming'
  #----------------------------------------------------------------
  ## Output:
  # Maximum likelihood tree, with the name: Tree_Name
  
  h <- paste(blastdirectory,"/FastTree -nt -gtr < ",Seq_Name," > ",Tree_Name, sep="") #infer a tree for a nucleotide alignment with the GTR+CAT model
  system(h, intern=TRUE, wait=TRUE)
}

buildTree(blastdirectory=directory,"trimmed_A.fas","A_Tree.tre")


# plot tree for visualization
MyTree <- read.tree("FirstTree.tre")
plot(MyTree)
int_nodes <- unique(MyTree$edge[,1])
nr_ext_nodes <- length(MyTree$tip.label) #127

## Drop all Los Alamos tips
SHCStips <- which(grepl("[[:alpha:]]",MyTree$tip.label)==FALSE)
AlamosTips <- which(grepl("[[:alpha:]]",MyTree$tip.label))
MyTreeSHCS <- drop.tip(MyTree,AlamosTips)
write.tree(MyTreeSHCS,"MyTreeSHCS_F.tre")


par(mfrow=c(1,1))
plot(MyTreeSHCS, type= "cladogram")
plot(MyTreeSHCS, type= "fan")
plot(MyTreeSHCS, type= "unrooted")
summary(MyTreeSHCS)

#####
ATree<- read.tree("A_Tree.tre")
plot(ATree)
int_nodes_A <- unique(ATree$edge[,1]) #125
nr_ext_nodes_A <- length(ATree$tip.label) #960

## Drop all Los Alamos tips
SHCStips_A <- which(grepl("[[:alpha:]]",ATree$tip.label)==FALSE) #497
AlamosTips_A <- which(grepl("[[:alpha:]]",ATree$tip.label)) #463
A_SHCS_Tree <- drop.tip(ATree,AlamosTips_A)
write.tree(A_SHCS_Tree,"A_SHCS_Tree.tre")
plot(A_SHCS_Tree, type= "fan")
plot(A_SHCS_Tree)


# Build bootstrap alignments
bootstrapTree <- function(Tree_Name, BTree_Name){
  ## Input:
  #  Tree_Name: The tre file containing the phylogenetic tree  i.e. the file 'FirstTree.tre' created with 'buildTree'
  #----------------------------------------------------------------
  ## Output:
  # 
  
  h <- paste("fseqboot ",Seq_Path," -gb",BTree_Name," -gb")
  system(h)
}

bootstrapTree("FirstTree.tre", "Bootstrap.fseqboot")


