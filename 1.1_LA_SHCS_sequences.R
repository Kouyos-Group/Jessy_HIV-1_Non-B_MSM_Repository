# June 2020, JDR
# Part B - Phylogeny 

## First Time to load the sequences
### Make Tree
##########################################

## clear R's brain
rm(list=ls())


### Load the SHCS sequences, the test_info file and the pat.dta:
setwd("~/Desktop/SHCS/Input/db_0620")
SHCS_sequences <- read.FASTA(file="final0620.fas")
TEST_INFO <- read.dta13("test_info0620.dta")
setwd("~/Desktop/SHCS/Input/NEW")
PAT <- read.dta13("pat0320.dta")


#################################### 
### Subtype Classification Algorithm,  quality check and one sequence per patient:










#howMany <- merge(TEST_INFO[,c('header_id','id','sampledate')],PAT)
#howMany.MSM <- howMany[howMany$risk == 1,]
#length(unique(howMany.MSM$id))


#################################### 
# use LA sequences
myLosAlamosSequences1 <- read.fasta(file="LAallST.fasta") #


# quality check for the Los Alamos background sequences
qual_check_LA <- function(LA_file, removeSwiss){
  ## Input: 
  #LA_file: The file containing the background sequences
  #removeSwiss: boolean variable to indicate whether Swiss sequences should be removed: 0 = no, 1 = yes.
  #----------------------------------------------------------------
  ## Output:
  # Cleaned file containing the background sequences
  
  if (removeSwiss == 1){ # removing the Swiss sequences if 1
    LA_file <- LA_file[!sapply(names(LA_file),FUN = function(a){unlist(strsplit(a,split = "[.]"))[2]}) == "CH"]
  }
  #rename to only accession number/ [5] = last number ("A07867")
  names(LA_file)<-sapply(names(LA_file),FUN=function(a){unlist(strsplit(a,split="[.]"))[5]})
  LA_file <- LA_file[!duplicated(names(LA_file))] #remove duplicates
  return(LA_file)
}


# quality checked LA:
clean_LA1 <- qual_check_LA(LosAlamosSequences1,1) # removing the Swiss sequences if 1
#clean_LA <- append(clean_LA1, clean_LA2)

write.fasta(clean_LA,names(clean_LA),"clean_LA.fasta")





#################################### 
## BLAST
# Basic Local Alignment Search Tool: finds regions of similarity between biological sequences. 
#The program compares nucleotide or protein	sequences to sequence databases and calculates the statistical significance.

####################################
# makeblastdb
makedatabase <- function(fastadirectory, fastafilename, blastdirectory){
  ## Output:
  #  produces BLAST databases from FASTA files
  stringcommand<-paste(blastdirectory ,"makeblastdb -in ", 
                       fastadirectory, fastafilename, " -input_type fasta -dbtype nucl -title testdatabase -out ", 
                       fastadirectory, "testdatabase", sep = "")
  system(stringcommand)
}


system("makeblastdb -in clean_LA.fasta -input_type fasta -dbtype nucl -title clean_LA -out clean_LA")


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
  fastafulldir<-paste(seqdirectory, fastaname, sep = "")
  jkl<-read.FASTA(fastafulldir)
  mno<-del.gaps(jkl)
  write.dna(mno, "sequencesplushitsbutunaligned.fas", "fasta")
  syscommand <- paste(blastdirectory, " -query ", seqdirectory, 
                      "sequencesplushitsbutunaligned.fas -out ", seqdirectory, 
                      "seq.txt -outfmt 6 -db ", seqdirectory, 
                      "clean_LA -evalue 0.05 -max_target_seqs ", 
                      maxtar," -max_hsps 1 -perc_identity ", percid, sep = "")
  system(syscommand)
  dbfulldir<- paste(seqdirectory, dbname, sep = "")
  database<-read.dna(dbfulldir, format = "fasta", as.matrix = FALSE)
  blasthitsdir<-paste(seqdirectory, "seq.txt", sep = "")
  blasthits<-read.table(blasthitsdir)
  blasthits<-as.matrix(blasthits)
  hitnames<-blasthits[,2]
  for (q in 1:length(hitnames)) {
    outputdir<- paste(seqdirectory, "sequencesplushitsbutunaligned.fas", sep = "") 
    write.dna(database[toString(hitnames[q])], outputdir, format = "fasta", append = TRUE)
  }
  duptemp<-read.FASTA(outputdir)
  duptemp<-duptemp[!duplicated(names(duptemp))]
  write.dna(duptemp, outputdir, format = "fasta", append = FALSE)
}

numtargets <- 10 
percentid <- 90
blasting("", "clean_seq.fas", "blastn", numtargets, percentid, "clean_LA.fasta")



#################################### 
## MUSCLE


#VA containing the fasta Reference sequence = HXB2 Pol Sequence from position 2253-3870
HXB2_full <- read.GenBank("K03455.1")
write.dna(HXB2_full, file ="HXB2_full.fas", format = "fasta", colsep = "")
HXB2_full_trans <- read.dna(file = "HXB2_full.fas",as.matrix = TRUE,format = "fasta")
write.dna(HXB2_full_trans[1,2253:3870], format = "fasta", "HXB2.fas") 
Reference <- read.fasta(file = "HXB2.fas")

# VA containing the fasta Sequences (SHCS + Background) output of the BLAST function
Sequences <- read.fasta(file = "sequencesplushitsbutunaligned.fas")


Pairwise_Alignment <- function(Ref, Seqs){
  ## Input:
  # Ref: The fasta file conaining the reference sequence
  # Seqs: The fasta file containing the sequences which should be aligned,  i.e. the file 'sequencesplushitsbutunaligned.fas' created with 'blasting'
  #----------------------------------------------------------------
  ## Output:
  # A new fasta file 'Seqs_Aligned_gaps_removed.fas' with the aligned sequences is created.
  # The sequences which were blasted before are aligned against a reference sequence using the MUSCLE executable file
  
  for(i in 1:length(Seqs)){
    write.fasta(sequences = c(Ref,Seqs[i]), names = names(c(Ref,Seqs[i])) , file.out = "To_be_Aligned.fasta")
    system("muscle -in To_be_Aligned.fasta -out Seqs_Aligned.fasta", intern=TRUE, wait=TRUE)
    Seqs_Aligned <- read.dna(file = "Seqs_Aligned.fasta",as.character = TRUE,as.matrix = TRUE,format = "fasta")
    keepCol <- which(Seqs_Aligned[1,]!="-") ## check where the gaps in HXB2 are
    Seqs_Aligned_gaps_removed <- t(Seqs_Aligned[2,keepCol])# w
    Seqs_Aligned_gaps_removed <- as.DNAbin(Seqs_Aligned_gaps_removed)
    write.dna(Seqs_Aligned_gaps_removed,append = TRUE,"Seqs_Aligned_gaps_removed.fasta",format = "fasta")
  }
  Seqs_Aligned_gaps_removed <- as.list(read.dna("Seqs_Aligned_gaps_removed.fasta",format = "fasta"))
  names(Seqs_Aligned_gaps_removed) <- names(Seqs)
  write.dna(Seqs_Aligned_gaps_removed,"Seqs_Aligned_gaps_removed.fasta",format = "fasta")
}


Pairwise_Alignment(Reference, Sequences)

deleteMutationPositions <- function(fas_name){
  ## Input:
  # fas_name: The fasta file 'Seqs_Aligned_gaps_removed.fas'  with the aligned sequences
  #----------------------------------------------------------------
  ## Output:
  # A new fasta file 'shcs_seq,Seqs_Aligned_gaps_removed.fas' with the deleted mutation poisitions is created.
  
  tempseq <-read.dna(fas_name, format = "fasta", as.character = TRUE, as.matrix = TRUE)
  #deletions are done on columns of matrix of alignment  
  tempseq<-tempseq[, -c(89,90,91,95,96,97,98,99,100,137,138,139,140,141,142,143,144,145,149,150,151,161,162,163,173,174,175,221,222,223,227,228,229,245,246,247,248,249,250,251,252,253,263,264,265,269,270,271,419,420,421,482,483,484,491,492,493,497,498,499,506,507,508,518,519,520,521,522,523,527,528,529,596,597,598,599,600,601,605,606,607,614,615,616,620,621,622,641,642,643,644,645,646,710,711,712,749,750,751,833,834,835,839,840,841,848,849,850,860,861,862,866,867,868,926,927,928,941,942,943,953,954,955,959,960,961,971,972,973,977,978,979,986,987,988)]   
  shcs_seq<-as.DNAbin(tempseq)
  write.dna(shcs_seq, fas_name,format="fasta",nbcol=ncol(shcs_seq))
}


deleteMutationPositions("Seqs_Aligned_gaps_removed.fasta")




# Trimming
## Only first time:
# download.file("http://trimal.cgenomics.org/_media/trimal.v1.2rev59.tar.gz", "trimal.v1.2rev59.tar.gz")
# untar("trimal.v1.2rev59.tar.gz")
# in command line: make the package

trimming <- function(wd, GTvalue, CONSvalue){
  ## Input:
  # wd: working directory
  # GTvalue: 0.7
  # CONSvalue: 0.5
  #----------------------------------------------------------------
  ## Output:
  # The sequences end are trimmed
  # A new fasta file 'trimmed.fas' is created
  
  commandstring<- paste(wd, "\\trimal.exe -in ", wd, "\\pretrimmed.fas -gt ", GTvalue, " -cons ", CONSvalue, " -out ", wd, "\\trimmed.fas", sep = "")  
  system(commandstring)
}

system("trimAl/source/trimal -in Seqs_Aligned_gaps_removed.fas -gt 0.7 -cons 0.5 -out trimmed.fas")




#################################### 
## FastTree
# Only first time: 
# system("gcc -O3 -finline-functions -funroll-loops -Wall -o FastTree FastTree.c -lm")

buildTree <- function(Seq_Name, Tree_Name){
  ## Input:
  # Seq_Name: The fasta file containing the sequences which should be trimmed i.e. the file 'trimmed.fas' created with 'trimming'
  #----------------------------------------------------------------
  ## Output:
  # Maximum likelihood tree, with the name: Tree_Name
  
  h <- paste("./FastTree -nt -gtr < ",Seq_Name," > ",Tree_Name)
  system(h, intern=TRUE, wait=TRUE)
}

buildTree("trimmed.fas","FirstTree.tre")


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


