---
# title: "Geoduck_TagSeq_DESeq2"
# author: "Samuel Gurr"
# date: "December 29, 2020"
---
  
# LOAD PACKAGES
library(DESeq2) # note: this was previously installed with the command `BiocManager::install("DESeq2")`
library(pasilla) # note: this was previously installed with the command `BiocManager::install("pasilla")`
library(dplyr)
library(GenomicFeatures)
library(pheatmap)
library(data.table)
library(calibrate)
library(gplots)
library(RColorBrewer)
library(affycoretools) # note: this was previously installed with the BiocManager::install("affycoretools")
library(edgeR)
library(data.table)
library(EnhancedVolcano)  # note: this was previously installed with the command `BiocManager::install("EnhancedVolcano")`
library(pcaExplorer) 
library(vsn)
# devtools::install_github('kevinblighe/EnhancedVolcano') # requires ggplot2 verssion 3.3.3 # install and load package "remotes " to run 'install_version("ggplot2", version = "3.3.3", repos = "http://cran.us.r-project.org")'


# SET WORKING DIRECTORY AND LOAD DATA
setwd("C:/Users/samjg/Documents/My_Projects/Pgenerosa_TagSeq_Metabolomics/TagSeq/")

# LOAD DATA
cts <- read.csv(file="HPC_Bioinf/outputs/transcript_count_matrix.csv", sep=',', header=TRUE) # read the output count matrix from prepDE.py
UT_seq_map <- read.csv(file="20201020_Gurr_TagSeq_UTAustin.csv", sep=',', header=TRUE)
smpl_ref <- read.csv(file="Sample_reference.csv", sep=',', header=TRUE)
treatment_ref <- read.csv(file="Extraction_checklist.csv", sep=',', header=TRUE)

# =====================================================================================================================
#                                                                                                                       
#                             FORMAT EXPERIMENT DESIGN DATAFRAME
# ======================================================================================================================

# MASTER REFERENCE DATA.FRAME
# format and merge to buld master reference dataframe
smpl_ref$Seq_Pos <- paste(smpl_ref$ï..TagSeq_Plate, smpl_ref$TagSeq_Well, sep="_")
smpl_ref <- smpl_ref[,-c(1:2)]
UT_seq_map$Seq_Pos  <- paste(UT_seq_map$ï..Plate, UT_seq_map$Well, sep="_")
UT_seq_map <- UT_seq_map[-c(1:2)]
Seq.Ref <- merge(smpl_ref, UT_seq_map, by = "Seq_Pos")
Mstr.Ref <- merge(Seq.Ref, treatment_ref, by = "Geoduck_ID")

#====================================================================================================================== #
# ALL TIMEPOINTS:
# call all experiment design treatments as 'exp.data'
exp.data <- Mstr.Ref[,c("Sample.Name","All_Treatment", "Primary_Treatment", "Second_Treament", "Third_Treatment", "Time")]


# NOTE: the following datasets expore the pairwise effects through time to explore the persistant DEGs 
# we found that primary treatment (first 110 days) affected diff expression regardless of subseq encounters
# thus, I want to investigate these persistant DEGs in Ambient vs. Moderate through time
# Ambient control (ambient throughout)
exp.data.AmbControl <- exp.data %>%  dplyr::filter(All_Treatment %in% c('A', 'AA', 'AAA')) # call all controls under ambient
exp.data.AmbControl$Treat.day <- paste(exp.data.AmbControl$All_Treatment, (substr(exp.data.AmbControl$Time, 4,5)), sep='_')

# Moderate throughout 
exp.data.Moderate <- exp.data %>%  dplyr::filter(All_Treatment %in% c('M', 'MM', 'MMM')) # call all controls under ambient
exp.data.Moderate$Treat.day <- paste(exp.data.Moderate$All_Treatment, (substr(exp.data.Moderate$Time, 4,5)), sep='_')


#====================================================================================================================== #
# DAY 0 
exp.data.d0 <- exp.data %>% dplyr::filter(Time %in% 'Day0')
nrow(exp.data.d0) # 8 total samples on Day 0
#====================================================================================================================== #
# DAY 7 
exp.data.d7 <- exp.data %>% dplyr::filter(Time %in% 'Day7')
nrow(exp.data.d7) # 36 total samples on Day 0
#====================================================================================================================== #
# DAY 14 
exp.data.d14 <- exp.data %>% dplyr::filter(Time %in% 'DAY14')
nrow(exp.data.d14) # 36 total samples on Day 0
#====================================================================================================================== #
# DAY 21 
exp.data.d21 <- exp.data %>% dplyr::filter(Time %in% 'DAY21')
nrow(exp.data.d21) # 62 total sampels on day 21


#======================================================================================================================
#                                                                                                                       
#                             Write treatment data for WGCNA analysis
#  
#====================================================================================================================== 

# write csv
path.wgcna = 'C:/Users/samjg/Documents/My_Projects/Pgenerosa_TagSeq_Metabolomics/TagSeq/RAnalysis/WGCNA/'
write.csv(exp.data, paste(path.wgcna,"all.treatment.data.csv"))
write.csv(exp.data.d0, paste(path.wgcna,"d0.treatment.data.csv"))
write.csv(exp.data.d7, paste(path.wgcna,"d7.treatment.data.csv"))
write.csv(exp.data.d14, paste(path.wgcna,"d14.treatment.data.csv"))
write.csv(exp.data.d21, paste(path.wgcna,"d21.treatment.data.csv"))

#======================================================================================================================
#                                                                                                                       
#                             FORMAT COUNT MATRIX
#                                                                                                                       
#====================================================================================================================== 

# About: need to call the cts.matrix and call only samples that match the IDs at day 0 and 21 
# Why? DESeq2 can only apply to factors with >=2 levels, so we cannot address 'Day' at just "Day21" for example
# alternatively, we can call subset matrices for the timepoints of interest

# NOTE: this cts matrix was NOT merged by lanes, run'rowsum' on unique delimiters to merge 'cts.merged'
#  Format count matrix (cts has 2x columns per sample on different lanes - sum counts from unique columns + build matrix)
ncol(cts) # 282 samples (not counting gene ID column) - should be 141 samples, need to sum columns by unique ID
cts.merged <- data.frame(cts[,-1], row.names=cts[,1]) # call new dataframe with first column now as row names, now all row values are numeric
names(cts.merged) <- sapply(strsplit(names(cts.merged), "_"), '[', 1) # split the column names by "_" delimiter  and call the first field SG##
cts.merged <- t(rowsum(t(cts.merged), group = colnames(cts.merged), na.rm = TRUE)) # merge all unique columns and sum counts 
ncol(cts.merged) # now 141 samples
# path for outputting all .csv filtered count files
path = 'C:/Users/samjg/Documents/My_Projects/Pgenerosa_TagSeq_Metabolomics/TagSeq/RAnalysis/DESeq2/' # run this for all count matrix outputs!!!
# ========================================================== 
#
# ALL TIMEPOINTS (3 CPM in 50% samples using edgeR)
# ========================================================== 
cts.matrix <-as.matrix(cts.merged, row.names="transcript_id") # call dataframe as matrix
ncol(cts.matrix) # 141 samples
nrow(cts.matrix) # 34947 total genes

cts.merged.as.table <- data.frame(transcript_id = row.names(cts.merged), cts.merged) # add back the rownames 'transcript_ID'
rownames(cts.merged.as.table) <- NULL # ommit the rownames
ncol(cts.merged.as.table) # 142 counting the transcript.ID that we want to keep! 
# pre-filtering; genes ommitted if < 3 counts per million reads in 50% of samples
cts.matrix.all <- cts.matrix
colSums(cts.matrix.all) # view the colSums of our all samples  - notice the read sums are around 1 million
CPM.all <- cpm(cts.matrix.all) # Obtain CPMs (counts oer million) using egdeR
head(CPM.all) # Have a look at the output
thresh.all <- CPM.all > 3 # Which values in myCPM are greater than 3?
head(thresh.all) # This produces a logical matrix with TRUEs and FALSES
rowSums(head(thresh.all)) # Summary of how many TRUEs there are in each row
table(rowSums(thresh.all)) # 2618 genes with TRUE in all 141 samples 
keep.all <- rowSums(thresh.all) >= (ncol(thresh.all)/2) # we would like to keep genes that have at least 50% TRUES in each row of thresh
summary(keep.all) # FALSE 26598 & TRUE  8349 -- more than 2/3 of the genes did not pass
cts.matrix.all.filtered <- cts.matrix.all[keep.all,] # Subset the rows of countdata to keep the more highly expressed genes
dim(cts.matrix.all.filtered) # 13874 genes & 141 samples
head(cts.matrix.all.filtered) # FINAL ALL DATASET
# We will look at a cople samples to check our defined threshold.. 
plot(CPM.all[,1], cts.matrix.all[,1], xlab="CPM", ylab="Raw Count", ylim=c(0,50), xlim=c(0,50), main=colnames(CPM.all)[1]) # check if our threshold of 3 CPM in 50% samples corresponds to a count of about 10-15
abline(v=3) 
abline(h=3)
plot(CPM.all[,70], cts.matrix.all[,70], xlab="CPM", ylab="Raw Count", ylim=c(0,50), xlim=c(0,50), main=colnames(CPM.all)[1]) # check if our threshold of 3 CPM in 50% samples corresponds to a count of about 10-15
abline(v=2.37) 
abline(h=3)

# write csv
all.counts.filtered <- cbind(rownames(cts.matrix.all.filtered), data.frame(cts.matrix.all.filtered, row.names=NULL))
colnames(all.counts.filtered)[1] <- "Gene.ID"
write.csv(all.counts.filtered, paste(path,"all.counts.filtered.csv"))


# ========================================================== 
#
# DAY 0  (3 CPM in 50% samples using edgeR)
# ========================================================== 
# About: run dyplr 'antijoin' to call cts columns that match 'Sample.Name' in the data frame 'exp.data.d0'
cts.merged.d0 <- cts.merged.as.table[,c(1,na.omit(match(exp.data.d0$Sample.Name, colnames(cts.merged.as.table))))]
cts.merged.d0 <- data.frame(cts.merged.d0[,-1], row.names=cts.merged.d0[,1])
cts.matrix.d0  <-as.matrix(cts.merged.d0, row.names="transcript_id")
ncol(cts.matrix.d0) # 8  samples from just Day 0
nrow(cts.matrix.d0) # 34947 total genes
colnames(cts.matrix.d0) == exp.data.d0$Sample.Name # check if TRUE, means the same as the exp/design dataframe exp.data.d0
# pre-filtering; genes ommitted if < 3 counts per million reads in 50% of samples
colSums(cts.matrix.d0) # view the colSums of our Day0 samples  - notice the read sums are around 1 million
CPM.d0 <- cpm(cts.matrix.d0) # Obtain CPMs (counts oer million) using egdeR
head(CPM.d0) # Have a look at the output
thresh.d0 <- CPM.d0 > 3 # Which values in myCPM are greater than 3?
head(thresh.d0) # This produces a logical matrix with TRUEs and FALSES
rowSums(head(thresh.d0)) # Summary of how many TRUEs there are in each row
table(rowSums(thresh.d0)) # 9631 genes with TRUE in all 8 samples 
keep.d0 <- rowSums(thresh.d0) >= (ncol(thresh.d0)/2) # we would like to keep genes that have at least 50% TRUES in each row of thresh
summary(keep.d0) # FALSE 19911 & TRUE 15036 -- more than half of the genes did not pass
cts.matrix.d0.filtered <- cts.matrix.d0[keep.d0,] # Subset the rows of countdata to keep the more highly expressed genes
dim(cts.matrix.d0.filtered) # 15036 genes & 8 samples
head(cts.matrix.d0.filtered) # FINAL DAY 0 DATASET
# We will look at a cople samples to check our defined threshold.. 
plot(CPM.d0[,1], cts.matrix.d0[,1], xlab="CPM", ylab="Raw Count", ylim=c(0,50), xlim=c(0,50), main=colnames(CPM.d0)[1]) # check if our threshold of 3 CPM in 50% samples corresponds to a count of about 10-15
abline(v=10) 
abline(h=10)
plot(CPM.d0[,5], cts.matrix.d0[,5], xlab="CPM", ylab="Raw Count", ylim=c(0,50), xlim=c(0,50), main=colnames(CPM.d0)[1]) # check if our threshold of 3 CPM in 50% samples corresponds to a count of about 10-15
abline(v=10) 
abline(h=10)

# write csv
day0.counts.filtered <- cbind(rownames(cts.matrix.d0.filtered), data.frame(cts.matrix.d0.filtered, row.names=NULL))
colnames(day0.counts.filtered)[1] <- "Gene.ID"
write.csv(day0.counts.filtered, paste(path,"day0.counts.filtered.csv")) # 'path' called in previous # write .csv section

# ========================================================== 
#
# DAY 7 (3 CPM in 50% samples using edgeR)
# ========================================================== 
cts.merged.d7 <- cts.merged.as.table[,c(1,na.omit(match(exp.data.d7$Sample.Name, colnames(cts.merged.as.table))))]
cts.merged.d7 <- data.frame(cts.merged.d7[,-1], row.names=cts.merged.d7[,1])
cts.matrix.d7  <-as.matrix(cts.merged.d7, row.names="transcript_id")
ncol(cts.matrix.d7) # 36 samples from just Day 7
colnames(cts.matrix.d7) == exp.data.d7$Sample.Name # check if TRUE, means the same as the exp/design dataframe exp.data.d7
# pre-filtering; genes ommitted if < 3 counts per million reads in 50% of samples
colSums(cts.matrix.d7) # view the colSums of our Day7 samples 
CPM.d7 <- cpm(cts.matrix.d7) # Obtain CPMs (counts oer million) using egdeR
head(CPM.d7) # Have a look at the output
thresh.d7 <- CPM.d7 > 3 # Which values in myCPM are greater than 3?
head(thresh.d7) # This produces a logical matrix with TRUEs and FALSES
rowSums(head(thresh.d7)) # Summary of how many TRUEs there are in each row
table(rowSums(thresh.d7)) # 6718 genes with TRUE in all 36 samples 
keep.d7 <- rowSums(thresh.d7) >= (ncol(thresh.d7)/2) # we would like to keep genes that have at least 50% TRUES in each row of thresh
summary(keep.d7) # FALSE 21024 & TRUE 13923 -- more than half of the genes did not pass
cts.matrix.d7.filtered <- cts.matrix.d7[keep.d7,] # Subset the rows of countdata to keep the more highly expressed genes
dim(cts.matrix.d7.filtered) # 13923 genes & 36 samples
head(cts.matrix.d7.filtered) # FINAL DAY 7 DATASET
# We will look at a cople samples to check our defined threshold.. 
plot(CPM.d7[,1], cts.matrix.d7[,1], xlab="CPM", ylab="Raw Count", ylim=c(0,50), xlim=c(0,50), main=colnames(CPM.d7)[1]) # check if our threshold of 3 CPM in 50% samples corresponds to a count of about 10-15
abline(v=10) 
abline(h=10)
plot(CPM.d7[,5], cts.matrix.d7[,5], xlab="CPM", ylab="Raw Count", ylim=c(0,50), xlim=c(0,50), main=colnames(CPM.d7)[1]) # check if our threshold of 3 CPM in 50% samples corresponds to a count of about 10-15
abline(v=10) 
abline(h=10)

# write csv
day7.counts.filtered <- cbind(rownames(cts.matrix.d7.filtered), data.frame(cts.matrix.d7.filtered, row.names=NULL))
colnames(day7.counts.filtered)[1] <- "Gene.ID"
write.csv(day7.counts.filtered, paste(path,"day7.counts.filtered.csv"))

# ========================================================== 
#
# DAY 14 (3 CPM in 50% samples using edgeR)
# ========================================================== 
cts.merged.d14 <- cts.merged.as.table[,c(1,na.omit(match(exp.data.d14$Sample.Name, colnames(cts.merged.as.table))))]
cts.merged.d14 <- data.frame(cts.merged.d14[,-1], row.names=cts.merged.d14[,1])
cts.matrix.d14  <-as.matrix(cts.merged.d14, row.names="transcript_id")
ncol(cts.matrix.d14) # 35 samples from just Day 14
colnames(cts.matrix.d14) == exp.data.d14$Sample.Name # SG92 in exp.data.d14$Sample.Name but not colnames(cts.matrix.d14)
UT_seq_map %>% dplyr::filter(Sample.Name == "SG92") # there was no sample in SG92 for TagSeq; 35 total is correct!
# pre-filtering; genes ommitted if < 3 counts per million reads in 50% of samples
colSums(cts.matrix.d14) # view the colSums of our Day14 samples - notice counts are near 1 million
CPM.d14 <- cpm(cts.matrix.d14) # Obtain CPMs (counts oer million) using egdeR
head(CPM.d14) # Have a look at the output
thresh.d14 <- CPM.d14 > 3 # Which values in myCPM are greater than 3?
head(thresh.d14) # This produces a logical matrix with TRUEs and FALSES
rowSums(head(thresh.d14)) # Summary of how many TRUEs there are in each row
table(rowSums(thresh.d14)) # 3451 genes with TRUE in all 35 samples 
keep.d14 <- rowSums(thresh.d14) >= (ncol(thresh.d14)/2) # we would like to keep genes that have at least 50% TRUES in each row of thresh
summary(keep.d14) # FALSE 26406 & TRUE 8541 -- more than half of the genes did not pass
cts.matrix.d14.filtered <- cts.matrix.d14[keep.d14,] # Subset the rows of countdata to keep the more highly expressed genes
dim(cts.matrix.d14.filtered) # 8541 genes & 35 samples
head(cts.matrix.d14.filtered) # FINAL DAY 14 DATASET
# We will look at a cople samples to check our defined threshold.. 
plot(CPM.d14[,1], cts.matrix.d14[,1], xlab="CPM", ylab="Raw Count", ylim=c(0,50), xlim=c(0,50), main=colnames(CPM.d14)[1]) # check if our threshold of 3 CPM in 50% samples corresponds to a count of about 10-15
abline(v=10) 
abline(h=10)
plot(CPM.d14[,5], cts.matrix.d14[,5], xlab="CPM", ylab="Raw Count", ylim=c(0,50), xlim=c(0,50), main=colnames(CPM.d14)[1]) # check if our threshold of 3 CPM in 50% samples corresponds to a count of about 10-15
abline(v=10) 
abline(h=10)

# write csv
day14.counts.filtered <- cbind(rownames(cts.matrix.d14.filtered), data.frame(cts.matrix.d14.filtered, row.names=NULL))
colnames(day14.counts.filtered)[1] <- "Gene.ID"
write.csv(day14.counts.filtered, paste(path,"day14.counts.filtered.csv"))

# ========================================================== 
#
# DAY 21  (3 CPM in 50% samples using edgeR)
# ========================================================== 
cts.merged.d21 <- cts.merged.as.table[,c(1,na.omit(match(exp.data.d21$Sample.Name, colnames(cts.merged.as.table))))]
cts.merged.d21 <- data.frame(cts.merged.d21[,-1], row.names=cts.merged.d21[,1])
cts.matrix.d21  <-as.matrix(cts.merged.d21, row.names="transcript_id")
ncol(cts.matrix.d21) # # 62 total sampels on day 21
colnames(cts.matrix.d21) == exp.data.d21$Sample.Name # check if TRUE, means the same as the exp/design dataframe exp.data.d21
# pre-filtering; genes ommitted if < 3 counts per million reads in 50% of samples
colSums(cts.matrix.d21) # view the colSums of our Day21 samples - notice SOME counts are near 1 million
CPM.d21 <- cpm(cts.matrix.d21) # Obtain CPMs (counts oer million) using egdeR
head(CPM.d21) # Have a look at the output
thresh.d21 <- CPM.d21 > 3 # filter CPM by threshold
head(thresh.d21) # This produces a logical matrix with TRUEs and FALSES
rowSums(head(thresh.d21)) # Summary of how many TRUEs there are in each row
table(rowSums(thresh.d21)) # 3451 genes with TRUE in all 35 samples 
keep.d21 <- rowSums(thresh.d21) >= (ncol(thresh.d21)/2) # we would like to keep genes that have at least 50% TRUES in each row of thresh
summary(keep.d21) # FALSE 26622 & TRUE 8325 -- more than three quarters of the genes did not pass
cts.matrix.d21.filtered <- cts.matrix.d21[keep.d21,] # Subset the rows of countdata to keep the more highly expressed genes
dim(cts.matrix.d21.filtered) # 13747 genes &  62 samples
head(cts.matrix.d21.filtered) # FINAL DAY 21 DATASET
# We will look at a cople samples to check our defined threshold.. 
plot(CPM.d21[,60], cts.matrix.d21[,60], xlab="CPM", ylab="Raw Count", ylim=c(0,50), xlim=c(0,50), main=colnames(CPM.d21)[1]) # check if our threshold of 3 CPM in 50% samples corresponds to a count of about 10-15
abline(v=10) 
abline(h=10)
plot(CPM.d21[,5], cts.matrix.d21[,5], xlab="CPM", ylab="Raw Count", ylim=c(0,50), xlim=c(0,50), main=colnames(CPM.d21)[1]) # check if our threshold of 3 CPM in 50% samples corresponds to a count of about 10-15
abline(v=10) 
abline(h=10)

# write csv
day21.counts.filtered <- cbind(rownames(cts.matrix.d21.filtered), data.frame(cts.matrix.d21.filtered, row.names=NULL))
colnames(day21.counts.filtered)[1] <- "Gene.ID"
write.csv(day21.counts.filtered, paste(path, "day21.counts.filtered.csv"))


#====================================================================================================================== 
#
#                             PREP dds files for DESeq2    
#
#====================================================================================================================== 


# ========================================================== 
# All TIMEPOINTS (test effect of time & Ambient_ALL vs. Moderate_ALL)
#---- INPUT: Experiment/design matrix == 'exp.data.AmbControl' and 'exp.data.Moderate'
#---- INPUT: TagSeq count matrix == 'cts.matrix.all.filtered' (3 CPM in 50% samples)
# ========================================================== 

# format Experiment/design dataframe into matrix
exp.data.Amb_all <- exp.data.AmbControl %>% dplyr::select(c('Sample.Name', 'Treat.day', 'All_Treatment')) # coondense dataset to build target matrix
exp.data.Amb_all$Treat.day <- factor(exp.data.Amb_all$Treat.day) # change Primary_Treatment to factor
# exp.data.d0.PRIMARY$Time <- factor(exp.data.d0.PRIMARY$Time) # change Time to factor
exp.data.Amb_all <- data.frame(exp.data.Amb_all[,-1], row.names=exp.data.Amb_all[,1]) # move Sample.Name column as row names  
exp.data.Amb_all.mtx <- as.matrix(exp.data.Amb_all, row.names="Geoduck.ID") # create matrix 

cts.merged.Ambient.all <- cts.matrix.all.filtered[,c(1,na.omit(match(exp.data.AmbControl$Sample.Name, colnames(cts.matrix.all.filtered))))]
dim(cts.merged.Ambient.all) # genes 13874   number samples 23

# CKECK THE cts.matrix AND THE exp.data.mtx
exp.data.Amb_all.mtx <- exp.data.Amb_all.mtx[match(colnames(cts.merged.Ambient.all),rownames(exp.data.Amb_all.mtx)), ]
all(rownames(exp.data.Amb_all.mtx) %in% colnames(cts.merged.Ambient.all)) # should be TRUE
all(rownames(exp.data.Amb_all.mtx) == colnames(cts.merged.Ambient.all)) # should be TRUE
all(rownames(exp.data.Amb_all.mtx) == colnames(cts.merged.Ambient.all))  # should be TRUE

# build dds
dds.Amb_all <- DESeqDataSetFromMatrix(countData = cts.merged.Ambient.all,
                                 colData = exp.data.Amb_all.mtx,
                                 design = ~ Treat.day) # DESeq Data Set (dds) - design as ~Primary_Treatment
dds.Amb_all # view dds

# prep for DESeq2
dds.Amb_all$Treat.day <- relevel(dds.Amb_all$Treat.day, ref = "A_0") # specify the reference level for count analysis - A = the control treatment
dds.Amb_all$All_Treatment <- droplevels(dds.Amb_all$All_Treatment)
levels(dds.Amb_all$All_Treatment) # NULL
levels(dds.Amb_all$Treat.day) #  levels for condition 'Treatment': "A_0"    "AA_14"  "AA_7"   "AAA_21"
design(dds.Amb_all) # view the design we have specified 
ncol(dds.Amb_all) # 23 cols - Good
nrow(dds.Amb_all) # 13874 (cut-off of 3 CPM)

# DEG ANALYSIS READY!!

# ========================================================== 
# DAY 0 
#---- INPUT: Experiment/design matrix == 'exp.data.d0' 
#---- INPUT: TagSeq count matrix == 'cts.matrix.d0.filtered' (3 CPM in 50% samples)
# ========================================================== 

# format Experiment/design dataframe into matrix
exp.data.d0.PRIMARY <- exp.data.d0 %>% dplyr::select(c('Sample.Name', 'Primary_Treatment', 'Time')) # coondense dataset to build target matrix
exp.data.d0.PRIMARY$Primary_Treatment <- factor(exp.data.d0.PRIMARY$Primary_Treatment) # change Primary_Treatment to factor
exp.data.d0.PRIMARY$Time <- factor(exp.data.d0.PRIMARY$Time) # change Time to factor
exp.data.d0.PRIMARY <- data.frame(exp.data.d0.PRIMARY[,-1], row.names=exp.data.d0.PRIMARY[,1]) # move Sample.Name column as row names  
exp.data.d0.PRIMARY.mtx <- as.matrix(exp.data.d0.PRIMARY, row.names="Geoduck.ID") # create matrix 
# fix(exp.data.mtx) # view data - remove # to open
# CKECK THE cts.matrix AND THE exp.data.mtx
exp.data.d0.PRIMARY.mtx <- exp.data.d0.PRIMARY.mtx[match(colnames(cts.merged.d0),rownames(exp.data.d0.PRIMARY.mtx)), ]
all(rownames(exp.data.d0.PRIMARY.mtx) %in% colnames(cts.matrix.d0.filtered)) # should be TRUE
all(rownames(exp.data.d0.PRIMARY.mtx) == colnames(cts.matrix.d0.filtered)) # should be TRUE
all(rownames(exp.data.d0.PRIMARY.mtx) == colnames(cts.matrix.d0.filtered))  # should be TRUE

# build dds
dds.d0 <- DESeqDataSetFromMatrix(countData = cts.matrix.d0.filtered,
                                      colData = exp.data.d0.PRIMARY.mtx,
                                      design = ~ Primary_Treatment) # DESeq Data Set (dds) - design as ~Primary_Treatment
dds.d0 # view dds

# prep for DESeq2
dds.d0$Primary_Treatment <- relevel(dds.d0$Primary_Treatment, ref = "A") # specify the reference level for count analysis - A = the control treatment
dds.d0$Time <- droplevels(dds.d0$Time)
levels(dds.d0$Time) # NULL
levels(dds.d0$Primary_Treatment) #  levels for condition 'Treatment': "A" "M"

design(dds.d0) # view the design we have specified 
ncol(dds.d0) # 8 cols - Good
nrow(dds.d0) # 9091 cols (cutoff at 10 CPM) 15036 (cut-off of 3 CPM)

# DEG ANALYSIS READY!!

# ========================================================== 
# DAY 7
#---- INPUT: Experiment/design matrix == 'exp.data.d7' 
#---- INPUT: TagSeq count matrix == 'cts.matrix.d7.filtered' (3 CPM in 50% samples)
# ==========================================================
exp.data.d7$Primary_Treatment <- factor(exp.data.d7$Primary_Treatment) # change Primary_Treatment to factor
exp.data.d7$Second_Treament <- factor(exp.data.d7$Second_Treament) # change Second_Treament to factor
exp.data.d7$Time <- factor(exp.data.d7$Time) # change Time to factor

#==================== #
# Primary*Second (ALL) Treatment
# format Experiment/design dataframe into matrix
exp.data.d7.PRIM_SEC <- exp.data.d7 %>% dplyr::select(c('Sample.Name', 'Primary_Treatment', 'Second_Treament', 'Time')) # coondense dataset to build target matrix
exp.data.d7.PRIM_SEC <- data.frame(exp.data.d7.PRIM_SEC[,-1], row.names=exp.data.d7.PRIM_SEC[,1]) # move Sample.Name column as row names  
exp.data.d7.PRIM_SEC.mtx <- as.matrix(exp.data.d7.PRIM_SEC, row.names="Geoduck.ID") # create matrix 
# CKECK THE cts.matrix AND THE exp.data.mtx
exp.data.d7.PRIM_SEC.mtx <- exp.data.d7.PRIM_SEC.mtx[match(colnames(cts.merged.d7),rownames(exp.data.d7.PRIM_SEC.mtx)), ]
all(rownames(exp.data.d7.PRIM_SEC.mtx) %in% colnames(cts.matrix.d7.filtered)) # should be TRUE
all(rownames(exp.data.d7.PRIM_SEC.mtx) == colnames(cts.matrix.d7.filtered)) # should be TRUE
all(rownames(exp.data.d7.PRIM_SEC.mtx) == colnames(cts.matrix.d7.filtered))  # should be TRUE
# build dds
dds.d7 <- DESeqDataSetFromMatrix(countData = cts.matrix.d7.filtered,
                                         colData = exp.data.d7.PRIM_SEC.mtx,
                                         design = ~ Primary_Treatment+Second_Treament) # DESeq Data Set (dds) - design as ~Primary_Treatment
# dds.primary.d0 # view dds
# prep for DESeq2
dds.d7$Primary_Treatment <- relevel(dds.d7$Primary_Treatment, ref = "A") # specify the reference level for count analysis - A = the control treatment
dds.d7$Second_Treament <- relevel(dds.d7$Second_Treament, ref = "A") # specify the reference level for count analysis - A = the control treatment
dds.d7$Time <- droplevels(dds.d7$Time) # drop levels of time
levels(dds.d7$Time) # NULL
levels(dds.d7$Primary_Treatment) #  levels for condition 'Treatment':  "A" "M"
levels(dds.d7$Second_Treament) #  levels for condition 'Treatment': "A" "M" "S"
design(dds.d7) # view the design we have specified 
ncol(dds.d7) # 36 cols - Good
nrow(dds.d7) # 8456 cols (10 CPM) 13923 (3 CPM)

# DEG ANALYSIS READY!!

# ------------------------------------------------------------------------------------------------------------------------------------------------- #

# Interaction term '...INT': All Treatment ( AA, AM, AS, MM, MA, MS)
# format Experiment/design dataframe into matrix
exp.data.d7.INT <- exp.data.d7 %>% dplyr::select(c('Sample.Name', 'All_Treatment', 'Time')) # coondense dataset to build target matrix
exp.data.d7.INT <- data.frame(exp.data.d7.INT[,-1], row.names=exp.data.d7.INT[,1]) # move Sample.Name column as row names  
exp.data.d7.INT.mtx <- as.matrix(exp.data.d7.INT, row.names="Geoduck.ID") # create matrix 
# CKECK THE cts.matrix AND THE exp.data.mtx
exp.data.d7.INT.mtx <- exp.data.d7.INT.mtx[match(colnames(cts.merged.d7),rownames(exp.data.d7.INT.mtx)), ]
all(rownames(exp.data.d7.INT.mtx) %in% colnames(cts.matrix.d7.filtered)) # should be TRUE
all(rownames(exp.data.d7.INT.mtx) == colnames(cts.matrix.d7.filtered)) # should be TRUE
all(rownames(exp.data.d7.INT.mtx) == colnames(cts.matrix.d7.filtered))  # should be TRUE
# build dds
dds.d7.INT <- DESeqDataSetFromMatrix(countData = cts.matrix.d7.filtered,
                                      colData = exp.data.d7.INT.mtx,
                                      design = ~ All_Treatment) # DESeq Data Set (dds) - design as ~All_Treatment
# dds.d7 # view dds
# prep for DESeq2
dds.d7.INT$All_Treatment <- relevel(dds.d7.INT$All_Treatment, ref = "AA") # specify the reference level for count analysis - A = the control treatment
dds.d7.INT$Time <- droplevels(dds.d7.INT$Time) # drop levels of time
levels(dds.d7.INT$Time) # NULL
levels(dds.d7.INT$All_Treatment) #  levels for condition 'Treatment': "AA" "AM" "AS" "MA" "MM" "MS"
design(dds.d7.INT) # view the design we have specified 
ncol(dds.d7.INT) # 36 cols - Good
nrow(dds.d7.INT)  # 8456 (10 CPM) 13923 (3 CPM)

# DEG ANALYSIS READY!!

# ========================================================== 
# DAY 14
#---- INPUT: Experiment/design matrix == 'exp.data.d14' 
#---- INPUT: TagSeq count matrix == 'cts.matrix.d14.filtered' (3 CPM in 50% samples)
# ========================================================== 
exp.data.d14$Primary_Treatment <- factor(exp.data.d14$Primary_Treatment) # change Primary_Treatment to factor
exp.data.d14$Second_Treament <- factor(exp.data.d14$Second_Treament) # change Second_Treament to factor
exp.data.d14$Time <- factor(exp.data.d14$Time) # change Time to factor

#==================== #
# Primary + Second; additive effect
# format Experiment/design dataframe into matrix
exp.data.d14.PRIM_SEC <- exp.data.d14 %>% dplyr::select(c('Sample.Name', 'Primary_Treatment', 'Second_Treament', 'Time')) # coondense dataset to build target matrix
exp.data.d14.PRIM_SEC <- data.frame(exp.data.d14.PRIM_SEC[,-1], row.names=exp.data.d14.PRIM_SEC[,1]) # move Sample.Name column as row names  
exp.data.d14.PRIM_SEC.mtx <- as.matrix(exp.data.d14.PRIM_SEC, row.names="Geoduck.ID") # create matrix 
# CKECK THE cts.matrix AND THE exp.data.mtx
exp.data.d14.PRIM_SEC.mtx <- exp.data.d14.PRIM_SEC.mtx[match(colnames(cts.merged.d14),rownames(exp.data.d14.PRIM_SEC.mtx)), ]
all(rownames(exp.data.d14.PRIM_SEC.mtx) %in% colnames(cts.matrix.d14.filtered)) # should be TRUE
all(rownames(exp.data.d14.PRIM_SEC.mtx) == colnames(cts.matrix.d14.filtered)) # should be TRUE
all(rownames(exp.data.d14.PRIM_SEC.mtx) == colnames(cts.matrix.d14.filtered))  # should be TRUE
# build dds
dds.d14 <- DESeqDataSetFromMatrix(countData = cts.matrix.d14.filtered,
                                     colData = exp.data.d14.PRIM_SEC.mtx,
                                     design = ~ Primary_Treatment+Second_Treament) # DESeq Data Set (dds) - design as ~Primary_Treatment
# dds.d14 # view dds
# prep for DESeq2
dds.d14$Primary_Treatment <- relevel(dds.d14$Primary_Treatment, ref = "A") # specify the reference level for count analysis - A = the control treatment
dds.d14$Second_Treament <- relevel(dds.d14$Second_Treament, ref = "A") # specify the reference level for count analysis - A = the control treatment
dds.d14$Time <- droplevels(dds.d14$Time) # drop levels of time
levels(dds.d14$Time) # NULL
levels(dds.d14$Primary_Treatment) #  levels for condition 'Treatment': "A" "M"
levels(dds.d14$Second_Treament) #  levels for condition 'Treatment': "A" "M" "S"
design(dds.d14) # view the design we have specified 
ncol(dds.d14) # 35 cols - Good
nrow(dds.d14)  # 8541 (10 CPM) 14241 (3 CPM)

# DEG ANALYSIS READY!!

# ------------------------------------------------------------------------------------------------------------------------------------------------- #

# Interaction term '...INT': All Treatment ( AA, AM, AS, MM, MA, MS)
# format Experiment/design dataframe into matrix
exp.data.d14.INT <- exp.data.d14 %>% dplyr::select(c('Sample.Name', 'All_Treatment', 'Time')) # coondense dataset to build target matrix
exp.data.d14.INT <- data.frame(exp.data.d14.INT[,-1], row.names=exp.data.d14.INT[,1]) # move Sample.Name column as row names  
exp.data.d14.INT.mtx <- as.matrix(exp.data.d14.INT, row.names="Geoduck.ID") # create matrix 
# CKECK THE cts.matrix AND THE exp.data.mtx
exp.data.d14.INT.mtx <- exp.data.d14.INT.mtx[match(colnames(cts.merged.d14),rownames(exp.data.d14.INT.mtx)), ]
all(rownames(exp.data.d14.INT.mtx) %in% colnames(cts.matrix.d14.filtered)) # should be TRUE
all(rownames(exp.data.d14.INT.mtx) == colnames(cts.matrix.d14.filtered)) # should be TRUE
all(rownames(exp.data.d14.INT.mtx) == colnames(cts.matrix.d14.filtered))  # should be TRUE
# build dds
dds.d14.INT <- DESeqDataSetFromMatrix(countData = cts.matrix.d14.filtered,
                                  colData = exp.data.d14.INT.mtx,
                                  design = ~ All_Treatment) # DESeq Data Set (dds) - design as ~All_Treatment
# dds.d14 # view dds
# prep for DESeq2
dds.d14.INT$All_Treatment <- relevel(dds.d14.INT$All_Treatment, ref = "AA") # specify the reference level for count analysis - A = the control treatment
dds.d14.INT$Time <- droplevels(dds.d14.INT$Time) # drop levels of time
levels(dds.d14$Time) # NULL
levels(dds.d14.INT$All_Treatment) #  levels for condition 'Treatment': "AA" "AM" "AS" "MA" "MM" "MS"
design(dds.d14.INT) # view the design we have specified 
ncol(dds.d14.INT) # 35 cols - Good
nrow(dds.d14.INT)  # 8541 (10 CPM) 14241 (3 CPM)

# DEG ANALYSIS READY!!

# ========================================================== 
# DAY 21 
#---- INPUT: Experiment/design matrix == 'exp.data.d21'
#---- INPUT: TagSeq count matrix == 'cts.matrix.d21.filtered' (3 CPM in 50% samples)
# ========================================================== 
exp.data.d21$Primary_Treatment <- factor(exp.data.d21$Primary_Treatment) # change Primary_Treatment to factor
exp.data.d21$Second_Treament <- factor(exp.data.d21$Second_Treament) # change Second_Treament to factor
exp.data.d21$Third_Treatment <- factor(exp.data.d21$Third_Treatment) # change Third_Treatment to factor
exp.data.d21$All_Treatment <- factor(exp.data.d21$All_Treatment) # change All_Treatment to factor
exp.data.d21$Time <- factor(exp.data.d21$Time) # change Time to factor

#==================== #
# Primary + Second +Third (ALL) Treatment
# format Experiment/design dataframe into matrix
exp.data.d21.PRIM_SEC_THIRD <- exp.data.d21 %>% dplyr::select(c('Sample.Name', 'Primary_Treatment','Second_Treament', 'Third_Treatment', 'Time')) # coondense dataset to build target matrix
exp.data.d21.PRIM_SEC_THIRD <- data.frame(exp.data.d21.PRIM_SEC_THIRD[,-1], row.names=exp.data.d21.PRIM_SEC_THIRD[,1]) # move Sample.Name column as row names  
exp.data.d21.PRIM_SEC_THIRD.mtx <- as.matrix(exp.data.d21.PRIM_SEC_THIRD, row.names="Geoduck.ID") # create matrix 
# CKECK THE cts.matrix AND THE exp.data.mtx
exp.data.d21.PRIM_SEC_THIRD.mtx <- exp.data.d21.PRIM_SEC.mtx[match(colnames(cts.merged.d21),rownames(exp.data.d21.PRIM_SEC_THIRD.mtx)), ]
all(rownames(exp.data.d21.PRIM_SEC_THIRD.mtx) %in% colnames(cts.matrix.d21.filtered)) # should be TRUE
all(rownames(exp.data.d21.PRIM_SEC_THIRD.mtx) == colnames(cts.matrix.d21.filtered)) # should be TRUE
all(rownames(exp.data.d21.PRIM_SEC_THIRD.mtx) == colnames(cts.matrix.d21.filtered))  # should be TRUE
# build dds
dds.d21 <- DESeqDataSetFromMatrix(countData = cts.matrix.d21.filtered,
                                      colData = exp.data.d21.PRIM_SEC_THIRD.mtx,
                                      design = ~ Primary_Treatment+Second_Treament+Third_Treatment) # DESeq Data Set (dds) - design as ~Primary_Treatment
# dds.d21 # view dds
# prep for DESeq2
dds.d21$Primary_Treatment <- relevel(dds.d21$Primary_Treatment, ref = "A") # specify the reference level for count analysis - A = the control treatment
dds.d21$Second_Treament <- relevel(dds.d21$Second_Treament, ref = "A") 
dds.d21$Third_Treatment <- relevel(dds.d21$Third_Treatment, ref = "A") 

dds.d21$Time <- droplevels(dds.d21$Time) # drop levels of time
levels(dds.d21$Time) # NULL
levels(dds.d21$Primary_Treatment) #  levels for condition 'Treatment': "A" "M"
design(dds.d21) # view the design we have specified 
ncol(dds.d21) #  62 cols - Good
nrow(dds.d21) # 8325 (10 CPM) 13506 (3 CPM)

# DEG ANALYSIS READY!!

# ------------------------------------------------------------------------------------------------------------------------------------------------- #

# Interaction term '...INT': All Treatment ( AA, AM, AS, MM, MA, MS)
# format Experiment/design dataframe into matrix
exp.data.d21.INT <- exp.data.d21 %>% dplyr::select(c('Sample.Name', 'All_Treatment', 'Time')) # coondense dataset to build target matrix
exp.data.d21.INT <- data.frame(exp.data.d21.INT[,-1], row.names=exp.data.d21.INT[,1]) # move Sample.Name column as row names  
exp.data.d21.INT.mtx <- as.matrix(exp.data.d21.INT, row.names="Geoduck.ID") # create matrix 
# CKECK THE cts.matrix AND THE exp.data.mtx
exp.data.d21.INT.mtx <- exp.data.d21.INT.mtx[match(colnames(cts.merged.d21),rownames(exp.data.d21.INT.mtx)), ]
all(rownames(exp.data.d21.INT.mtx) %in% colnames(cts.matrix.d21.filtered)) # should be TRUE
all(rownames(exp.data.d21.INT.mtx) == colnames(cts.matrix.d21.filtered)) # should be TRUE
all(rownames(exp.data.d21.INT.mtx) == colnames(cts.matrix.d21.filtered))  # should be TRUE
# build dds
dds.d21.INT <- DESeqDataSetFromMatrix(countData = cts.matrix.d21.filtered,
                                      colData = exp.data.d21.INT.mtx,
                                      design = ~ All_Treatment) # DESeq Data Set (dds) - design as ~All_Treatment
# dds.d21 # view dds
# prep for DESeq2
dds.d21.INT$All_Treatment <- relevel(dds.d21.INT$All_Treatment, ref = "AAA") # specify the reference level for count analysis - A = the control treatment
dds.d21.INT$Time <- droplevels(dds.d21.INT$Time) # drop levels of time
levels(dds.d21.INT$Time) # NULL
levels(dds.d21.INT$All_Treatment) #  levels for condition 'Treatment': "AAA" "AAM" "AMA" "AMM" "ASA" "ASM" "MAA" "MAM" "MMA" "MMM" "MSA" "MSM"
design(dds.d21.INT) # view the design we have specified 
ncol(dds.d21.INT) #  62 cols - Good
nrow(dds.d21.INT)  # 13747 (3 CPM)

# DEG ANALYSIS READY!!

#====================================================================================================================== 
#
#                          Differential expression analysis (DESeq2)
#
#====================================================================================================================== 

# ==================
# About and notes...
# ================== 
# dds.d0 (not filtered by edgeR = dds.primary.d0.NoFilt)
#
# dds.d7
# 
# dds.d14
# 
# dds.d21

# About: Log fold change shrinkage for visualization and ranking
# Run DESeq : Modeling counts with the specified design(dds) treatment effects 

# Note about NA values! =============================================== #
# p-values are set to NA by DESeq2 for the following reasions (review URL: https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#pvaluesNA)
# 1.) all samples in a row have zero counts 
# 2.) Cooks distance computed in a row finds an extreme count outlier - both the pvalue and the adjusted p value will be set to NA
# 3.) Row is filtered by automatic independent filtering (for having low mean normalized count) 
# the adjusted pvalue will be set to NA

# Note about filtering! =============================================== #
# DESeq2 runs independent filtering autopmatically unless specified otherwise
# filter is passed on a FDR cutoff  - default if alpha = 0.1
# can set the FDR cut-off manually as res <- results(dds, alpha=0.05)
# use 'independentFiltering=FALSE' to turn OFF this automatic filter
# example: res <- results(dds, independentFiltering=FALSE)
# in this R script (as of 20200108) we used edgeR to filter reads and design dds prior to DESeq2
# FDR == False Discovery Rates: think about the raw read counts, what is the liklihood that 
# a treatment effect (DEG) can arise due to group diffs in read count distribution? (false positive DEG)
# false positives are rare - 5% of the time is a common rule of thumb that sample distribtions do not overlap
# note we have ~35,000 genes so 5% is 1,750 genes! 
# VIDEO EXPLAINING FDR https://www.youtube.com/watch?v=K8LQSvtjcEo&feature=emb_logo



# ========================================================== 
#
# ALL TIMEPOINTS - Ambient controls    (3 CPM in 50% samples)
# ========================================================== 

nrow(dds.Amb_all) # 13874 total genes pre filtered; 

#  pre-filtered in edgeR - use independentFiltering=FALSE
dds.Amb_all <- DESeq(dds.Amb_all) # wait for this to complete....
resultsNames(dds.Amb_all) # view the names of your results model "Treat.day_AA_14_vs_A_0"  "Treat.day_AA_7_vs_A_0"   "Treat.day_AAA_21_vs_A_0"

resAmb.d0vd7 <- results(dds.Amb_all, name="Treat.day_AA_7_vs_A_0", alpha = 0.05)
hist(resAmb.d0vd7$pvalue, breaks=20, col="grey") # view histogram,  about 425 genes are likely false positives, alpha = 0.0
abline(h=(nrow(dds.Amb_all)*0.05),col="red") # add line at expected 5% false positive
table(resAmb.d0vd7$padj<0.05) # 0 DEGs


resAmb.d0vd14 <- results(dds.Amb_all, name="Treat.day_AA_14_vs_A_0", alpha = 0.05)
hist(resAmb.d0vd14$pvalue, breaks=20, col="grey") # view histogram,  about 425 genes are likely false positives, alpha = 0.0
abline(h=(nrow(dds.Amb_all)*0.05),col="red") # add line at expected 5% false positive
table(resAmb.d0vd14$padj<0.05) # 11 DEGs


resAmb.d0vd21 <- results(dds.Amb_all, name="Treat.day_AAA_21_vs_A_0", alpha = 0.05)
hist(resAmb.d0vd21$pvalue, breaks=20, col="grey") # view histogram,  about 425 genes are likely false positives, alpha = 0.0
abline(h=(nrow(dds.Amb_all)*0.05),col="red") # add line at expected 5% false positive
table(resAmb.d0vd21$padj<0.05) # 94 DEGs


# Write results
# path_out.d0 = 'C:/Users/samjg/Documents/My_Projects/Pgenerosa_TagSeq_Metabolomics/TagSeq/RAnalysis/DESeq2/output/Day0/'
# write.csv(resdata.d0.primary, paste(path_out.d0, "Day0.PrimaryTreatment_diffexpr-results.csv"))

# volcano plot ------------------------------------------------------------------------------------------------------ #

# png("RAnalysis/DESeq2/output/Day0/Day0.PrimaryTreatment-VolcanoPlot.png", 1000, 1000, pointsize=20)
EnhancedVolcano(resAmb.d0vd21,
                lab = rownames(resAmb.d0vd21),
                x = 'log2FoldChange',
                y = 'pvalue',
                title = 'Day 21 v. Day 0 (Ambient Control)',
                subtitle = "DESeq2 - Differential expression",
                FCcutoff = 2.0,
                pointSize = 4.0,
                labSize = 6.0,
                colAlpha = 1,
                legendPosition = 'right',
                legendLabSize = 12,
                legendIconSize = 4.0,
                widthConnectors = 0.75)
# dev.off()


# ========================================================== 
#
# DAY 0    (3 CPM in 50% samples)
# ========================================================== 
nrow(dds.d0) # 15036 total genes pre filtered; FDR cutoff  0.05 = 751.8 false positives

#  pre-filtered in edgeR - use independentFiltering=FALSE
dds.d0 <- DESeq(dds.d0) # wait for this to complete....
resultsNames(dds.d0) # view the names of your results model 'Primary_Treatment_M_vs_A'
resd0.primary <- results(dds.d0, alpha = 0.5) # already filtered in edgeR; FDR of 0.05 (~ 450 in 9091 genes)
hist(resd0.primary$pvalue, breaks=20, col="grey") # view histogram,  about 425 genes are likely false positives, alpha = 0.05 FDR is good!
abline(h=(nrow(dds.d0)*0.05),col="red") # add line at expected 5% false positive
table(resd0.primary$padj<0.05) #  13  DEGs; 11200 FALSE
resd0.primary <- resd0.primary[order(resd0.primary$padj), ] ## Order by adjusted p-value
resdata.d0.primary  <- merge(as.data.frame(resd0.primary), as.data.frame(counts(dds.d0, normalized=TRUE)), by="row.names", sort=FALSE) ## Merge with normalized count data
names(resdata.d0.primary)[1] <- "Gene"
head(resdata.d0.primary)

# Write results
path_out.d0 = 'C:/Users/samjg/Documents/My_Projects/Pgenerosa_TagSeq_Metabolomics/TagSeq/RAnalysis/DESeq2/output/Day0/'
write.csv(resdata.d0.primary, paste(path_out.d0, "Day0.PrimaryTreatment_diffexpr-results.csv"))

# volcano plot ------------------------------------------------------------------------------------------------------ #

png("RAnalysis/DESeq2/output/Day0/Day0.PrimaryTreatment-VolcanoPlot.png", 1000, 1000, pointsize=20)
EnhancedVolcano(resd0.primary,
                lab = rownames(resd0.primary),
                x = 'log2FoldChange',
                y = 'pvalue',
                title = 'Moderate versus Ambient Treatment (Primary)',
                subtitle = "DESeq2 - Differential expression",
                FCcutoff = 2.0,
                pointSize = 4.0,
                labSize = 6.0,
                colAlpha = 1,
                legendPosition = 'right',
                legendLabSize = 12,
                legendIconSize = 4.0,
                widthConnectors = 0.75)
dev.off()

# Plot dispersions ------------------------------------------------------------------------------------------------------ #

png("RAnalysis/DESeq2/output/Day0/Day0.dispersions.png", 1000, 1000, pointsize=20)
plotDispEsts(dds.d0, main="Day 0 dispersions")
dev.off()

# Data transformations  ------------------------------------------------------------------------------------------------------ #

# rlog - regularized log transformation of origin count data to log2 scale - fit for each sample and dist. of coefficients in the data
rld.d0<- rlogTransformation(dds.d0) # rlog transform (regularized log)
head(assays(rld.d0)) # view first few rows
hist(assay(rld.d0)) # view histogram 
meanSdPlot(assay(rld.d0)) # plot using "vsn" bioconductor package - shows the sd y axis (sq root of varaince in all samples) - flat curve may seem like a goals, BUT may be unreasonable in cases with MANY true DEGs from experimental conditions
# this gives log2(n + 1)
ntd.d0 <- normTransform(dds.d0)
head(assay(ntd.d0)) # view first few rows
hist(assay(ntd.d0)) # view histogram 
meanSdPlot(assay(ntd.d0)) # plot using "vsn" bioconductor package - shows the sd y axis (sq root of varaince in all samples) - flat curve may seem like a goals, BUT may be unreasonable in cases with MANY true DEGs from experimental conditions

# Plot heat map ------------------------------------------------------------------------------------------------------ #

##################  heat map of count matrix using "pheatmap" 

save_pheatmap <- function(x, filename, width=1000, height=960) { # template for saving pheatmap outputs
  stopifnot(!missing(x))
  stopifnot(!missing(filename))
  png(filename,width = width, height=height)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}

select <- order(rowMeans(counts(dds.d0,normalized=TRUE)), decreasing=TRUE)[1:20] 
df <- as.data.frame(colData(dds.d0)[c("Primary_Treatment")])
annotation_colors = list(Treatment = c(A="Blue", M="Orange"))

# rlog heatmap
d0.rlog.heatmap<- pheatmap(assay(rld.d0)[select,], cluster_rows=T, show_rownames=F, # rlog heatmap
         cluster_cols=T, annotation_col=df, annotation_colors = annotation_colors,
         main = "rlog_Day0_MvsA")
save_pheatmap(d0.rlog.heatmap, filename = "RAnalysis/DESeq2/output/Day0/Day0.PrimaryTreatment-rlog_heatmap.png")

# ntd heatmap
d0.ntd.heatmap<- pheatmap(assay(ntd.d0)[select,], cluster_rows=T, show_rownames=F, # rlog heatmap
         cluster_cols=T, annotation_col=df, annotation_colors = annotation_colors,
         main = "ntd_Day0_MvsA")
save_pheatmap(d0.ntd.heatmap, filename = "RAnalysis/DESeq2/output/Day0/Day0.PrimaryTreatment-ntd_heatmap.png")

##################  heatmap of the sample-to-sample distances

# Colors for plots below; Use RColorBrewer, better and assign mycols variable
mycols.d0 <- brewer.pal(8, "Dark2")[1:length(unique(rld.d0$Primary_Treatment))]
# Sample distance heatmap
sampleDists.d0 <- as.matrix(dist(t(assay(rld.d0))))
png("RAnalysis/DESeq2/output/Day0/Day0.PrimaryTreatment-heatmap.png", 1000, 1000, pointsize=20)
heatmap.2(as.matrix(sampleDists.d0), key=F, trace="none",
          col=colorpanel(100, "black", "white"),
          ColSideColors=mycols.d0[rld.d0$Primary_Treatment], RowSideColors=mycols.d0[rld.d0$Primary_Treatment],
          margin=c(10, 10), main="Moderate versus Ambient Treatment (Primary)")
dev.off()

# Principal components analysis ------------------------------------------------------------------------------------------ #
pcaExplorer(dds = dds.d0)
# DESeq2::plotPCA(rld.d0, intgroup="Primary_Treatment") ## Could do with built-in DESeq2 function:
# rld_pca.do.primary <- function (rld.d0.primary, intgroup = "Primary_Treatment", ntop = 500, colors=NULL, legendpos="bottomleft", main="PCA Biplot", textcx=1, ...) {
#   require(genefilter)
#   require(calibrate)
#   require(RColorBrewer)
#   rv = rowVars(assay(rld.d0.primary))
#   select = order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
#   pca = prcomp(t(assay(rld.d0.primary)[select, ]))
#   fac = factor(apply(as.data.frame(colData(rld.d0.primary)[, intgroup, drop = FALSE]), 1, paste, collapse = " : "))
#   if (is.null(colors)) {
#     if (nlevels(fac) >= 3) {
#       colors = brewer.pal(nlevels(fac), "Paired")
#     }   else {
#       colors = c("black", "red")
#     }
#   }
#   pc1var <- round(summary(pca)$importance[2,1]*100, digits=1)
#   pc2var <- round(summary(pca)$importance[2,2]*100, digits=1)
#   pc1lab <- paste0("PC1 (",as.character(pc1var),"%)")
#   pc2lab <- paste0("PC1 (",as.character(pc2var),"%)")
#   plot(PC2~PC1, data=as.data.frame(pca$x), bg=colors[fac], pch=21, xlab=pc1lab, ylab=pc2lab, main=main, ...)
#   with(as.data.frame(pca$x), textxy(PC1, PC2, labs=rownames(as.data.frame(pca$x)), cex=textcx))
#   legend(legendpos, legend=levels(fac), col=colors, pch=20)
# }
# png("Day0.Primary.qc-pca.png", 1000, 1000, pointsize=20)
# rld_pca.do.primary(rld.d0.primary, colors=mycols, intgroup="Primary_Treatment", xlim=c(-75, 35))
# dev.off()

# ========================================================== 
#
# DAY 7    (3 CPM in 50% samples)
# ========================================================== 
# Primary + Secondary Treatment  =========================================================================================== #
nrow(dds.d7) # 8456 total genes pre filtered; FDR cutoff of 0.1 = 845.6 false positives; 0.05 = 422.8 flase positives

dds.d7 <- DESeq(dds.d7) # wait for this to complete....
resultsNames(dds.d7) # view the names of your results model 

resd7.all.primary_M_vs_A<- results(dds.d7, name="Primary_Treatment_M_vs_A", alpha = 0.05)
hist(resd7.all.primary_M_vs_A$pvalue, breaks=20, col="grey") # view histogram,  about 425 genes are likely false positives, alpha = 0.0
abline(h=(nrow(dds.d7)*0.05),col="red") # add line at expected 5% false positive
table(resd7.all.primary_M_vs_A$padj<0.05) # 94 DEGs
resd7.all.primary_M_vs_A <- resd7.all.primary_M_vs_A[order(resd7.all.primary_M_vs_A$padj), ] ## Order by adjusted p-value
resdata.d7.all.prim.M_vs_A  <- merge(as.data.frame(resd7.all.primary_M_vs_A), as.data.frame(counts(dds.d7, normalized=TRUE)), by="row.names", sort=FALSE) ## Merge with normalized count data
names(resdata.d7.all.prim.M_vs_A)[1] <- "Gene" # assign col 1

resd7.all.second_M_vs_A<- results(dds.d7, name="Second_Treament_M_vs_A", alpha = 0.05) # FDR 5% 
hist(resd7.all.second_M_vs_A$pvalue, breaks=20, col="grey") # view histogram,  about 425 genes are likely false positives, alpha = 0.05 FDR is good!
abline(h=(nrow(dds.d7)*0.05),col="red") # add line at expected 5% false positive
table(resd7.all.second_M_vs_A$padj<0.05) # 105 DEGs
resd7.all.second_M_vs_A <- resd7.all.second_M_vs_A[order(resd7.all.second_M_vs_A$padj), ] ## Order by adjusted p-value
resdata.d7.second_M_vs_A  <- merge(as.data.frame(resd7.all.second_M_vs_A), as.data.frame(counts(dds.d7, normalized=TRUE)), by="row.names", sort=FALSE) ## Merge with normalized count data
names(resdata.d7.second_M_vs_A)[1] <- "Gene" # assign col 1

resd7.all.second_S_vs_A<- results(dds.d7, name="Second_Treament_S_vs_A", alpha = 0.05) # FDR 5% 
hist(resd7.all.second_S_vs_A$pvalue, breaks=20, col="grey") # view histogram,  about 425 genes are likely false positives, alpha = 0.05 FDR is good!
abline(h=(nrow(dds.d7)*0.05),col="red") # add line at expected 5% false positive
table(resd7.all.second_S_vs_A$padj<0.05) # 1 DEGs
resd7.all.second_S_vs_A <- resd7.all.second_S_vs_A[order(resd7.all.second_S_vs_A$padj), ] ## Order by adjusted p-value
resdata.d7.all.second_S_vs_A    <- merge(as.data.frame(resd7.all.second_S_vs_A), as.data.frame(counts(dds.d7, normalized=TRUE)), by="row.names", sort=FALSE) ## Merge with normalized count data
names(resdata.d7.all.second_S_vs_A)[1] <- "Gene" # assign col 1


## Write results
path_out.d7 = 'C:/Users/samjg/Documents/My_Projects/Pgenerosa_TagSeq_Metabolomics/TagSeq/RAnalysis/DESeq2/output/Day7/'
write.csv(resdata.d7.all.prim.M_vs_A, paste(path_out.d7, "Day7.PrimaryTreament_diffexpr-results.csv"))
write.csv(resdata.d7.second_M_vs_A, paste(path_out.d7, "Day7.SecondTreament_MvsA_diffexpr-results.csv"))
write.csv(resdata.d7.all.second_S_vs_A, paste(path_out.d7, "Day7.SecondTreament_SvsA_diffexpr-results.csv"))


# volcano plot ------------------------------------------------------------------------------------------------------ #

# Effect of primary treatment 'resd7.all.primary_M_vs_A' - 94 DEGs
png("RAnalysis/DESeq2/output/Day7/Day7.PrimaryTreatment-VolcanoPlot.png", 1000, 1000, pointsize=20)
EnhancedVolcano(resd7.all.primary_M_vs_A,
                lab = rownames(resd7.all.primary_M_vs_A),
                x = 'log2FoldChange',
                y = 'pvalue',
                title = 'Moderate versus Ambient Treatment (Primary)',
                subtitle = "DESeq2 - Differential expression",
                FCcutoff = 2.0,
                pointSize = 4.0,
                labSize = 6.0,
                colAlpha = 1,
                legendPosition = 'right',
                legendLabSize = 12,
                legendIconSize = 4.0,
                widthConnectors = 0.75)
dev.off()


# Effect of second treatment M v. A 'resdata.d7.second_M_vs_A' - 105 DEGs
png("RAnalysis/DESeq2/output/Day7/Day7.SecondTreatment_MvA-VolcanoPlot.png", 1000, 1000, pointsize=20)
EnhancedVolcano(resdata.d7.second_M_vs_A,
                lab = rownames(resdata.d7.second_M_vs_A),
                x = 'log2FoldChange',
                y = 'pvalue',
                title = 'Moderate versus Ambient Treatment (Second)',
                subtitle = "DESeq2 - Differential expression",
                FCcutoff = 2.0,
                pointSize = 4.0,
                labSize = 6.0,
                colAlpha = 1,
                legendPosition = 'right',
                legendLabSize = 12,
                legendIconSize = 4.0,
                widthConnectors = 0.75)
dev.off()

# Effect of second treatment S v. A 'resd7.all.second_S_vs_A' - 1 DEG
png("RAnalysis/DESeq2/output/Day7/Day7.SecondTreatment_SvA-VolcanoPlot.png", 1000, 1000, pointsize=20)
EnhancedVolcano(resd7.all.second_S_vs_A,
                lab = rownames(resd7.all.second_S_vs_A),
                x = 'log2FoldChange',
                y = 'pvalue',
                title = 'Severe versus Ambient Treatment (Second)',
                subtitle = "DESeq2 - Differential expression",
                FCcutoff = 2.0,
                pointSize = 4.0,
                labSize = 6.0,
                colAlpha = 1,
                legendPosition = 'right',
                legendLabSize = 12,
                legendIconSize = 4.0,
                widthConnectors = 0.75)
dev.off()

# Plot dispersions ------------------------------------------------------------------------------------------------------ #

png("RAnalysis/DESeq2/output/Day7/Day7-dispersions.png", 1000, 1000, pointsize=20)
plotDispEsts(dds.d7, main="Day 7 data dispersions")
dev.off()

# Data transformations  ------------------------------------------------------------------------------------------------------ #

# rlog - regularized log transformation of origin count data to log2 scale - fit for each sample and dist. of coefficients in the data
rld.d7<- rlogTransformation(dds.d7) # rlog transform (regularized log)
head(assays(rld.d7)) # view first few rows
hist(assay(rld.d7)) # view histogram 
meanSdPlot(assay(rld.d7)) # plot using "vsn" bioconductor package - shows the sd y axis (sq root of varaince in all samples) - flat curve may seem like a goals, BUT may be unreasonable in cases with MANY true DEGs from experimental conditions
# this gives log2(n + 1)
ntd.d7 <- normTransform(dds.d7)
head(assay(ntd.d7)) # view first few rows
hist(assay(ntd.d7)) # view histogram 
meanSdPlot(assay(ntd.d7)) # plot using "vsn" bioconductor package - shows the sd y axis (sq root of varaince in all samples) - flat curve may seem like a goals, BUT may be unreasonable in cases with MANY true DEGs from experimental conditions

# Plot heat map ------------------------------------------------------------------------------------------------------ #

##################  heat map of count matrix using "pheatmap" 

save_pheatmap <- function(x, filename, width=1000, height=960) { # template for saving pheatmap outputs
  stopifnot(!missing(x))
  stopifnot(!missing(filename))
  png(filename,width = width, height=height)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}

# settings for the heatmaps
select <- order(rowMeans(counts(dds.d7,normalized=TRUE)), decreasing=TRUE)[1:20] 
df <- as.data.frame(colData(dds.d7)[,c("Primary_Treatment","Second_Treament")])
annotation_colors = list(Treatment = c(A="Blue", M="Orange"))

# rlog heatmap
d7.rlog.heatmap<- pheatmap(assay(rld.d7)[select,], cluster_rows=T, show_rownames=F, # rlog heatmap
                           cluster_cols=T, annotation_col=df, annotation_colors = annotation_colors,
                           main = "rlog_Day7")
save_pheatmap(d7.rlog.heatmap, filename = "RAnalysis/DESeq2/output/Day7/Day7-rlog_heatmap.png")

# ntd heatmap
d7.ntd.heatmap<- pheatmap(assay(ntd.d7)[select,], cluster_rows=T, show_rownames=F, # rlog heatmap
                          cluster_cols=T, annotation_col=df, annotation_colors = annotation_colors,
                          main = "ntd_Day7")
save_pheatmap(d7.ntd.heatmap, filename = "RAnalysis/DESeq2/output/Day7/Day7-ntd_heatmap.png")

##################  heatmap of the sample-to-sample distances

# Colors for plots below; Use RColorBrewer, better and assign mycols variable
mycols.d7 <- brewer.pal(8, "Dark2")[1:length(unique(rld.d7$Primary_Treatment))]
mycols.d7_2 <- brewer.pal(8, "Dark2")[1:length(unique(rld.d7$Second_Treament))]
# Sample distance heatmap
sampleDists.d7 <- as.matrix(dist(t(assay(rld.d7))))
png("RAnalysis/DESeq2/output/Day7/Day7-heatmap.png", 1000, 1000, pointsize=20)
heatmap.2(as.matrix(sampleDists.d7), key=F, trace="none",
          col=colorpanel(100, "black", "white"),
          ColSideColors=mycols.d7_2[rld.d7$Second_Treament], RowSideColors=mycols.d7[rld.d7$Primary_Treatment],
          margin=c(10, 10), main="Day 7 Sample to Sample heatmap")
dev.off()

# Principal components analysis ------------------------------------------------------------------------------------------ #
pcaExplorer(dds = dds.d7)

# ========================================================== 
#
# DAY 14   (3 CPM in 50% samples)
# ========================================================== 

path_out.d14 = 'C:/Users/samjg/Documents/My_Projects/Pgenerosa_TagSeq_Metabolomics/TagSeq/RAnalysis/DESeq2/output/Day14/'

# Primary + Secondary Treatment  =========================================================================================== #
nrow(dds.d14) # 8456 total genes pre filtered; FDR cutoff of 0.1 = 845.6 false positives; 0.05 = 422.8 flase positives

dds.d14 <- DESeq(dds.d14) # wait for this to complete....
resultsNames(dds.d14) # view the names of your results model 

resd14.all.primary_M_vs_A<- results(dds.d14, name="Primary_Treatment_M_vs_A", alpha = 0.05)
hist(resd14.all.primary_M_vs_A$pvalue, breaks=20, col="grey") # view histogram,  about 425 genes are likely false positives, alpha = 0.0
abline(h=(nrow(dds.d14)*0.05),col="red") # add line at expected 5% false positive
table(resd14.all.primary_M_vs_A$padj<0.05) #  502 DEGs
resd14.all.primary_M_vs_A <- resd14.all.primary_M_vs_A[order(resd14.all.primary_M_vs_A$padj), ] ## Order by adjusted p-value
resdata.d14.all.prim.M_vs_A  <- merge(as.data.frame(resd14.all.primary_M_vs_A), as.data.frame(counts(dds.d14, normalized=TRUE)), by="row.names", sort=FALSE) ## Merge with normalized count data
names(resdata.d14.all.prim.M_vs_A)[1] <- "Gene" # assign col 1

resd14.all.second_M_vs_A<- results(dds.d14, name="Second_Treament_M_vs_A", alpha = 0.05) # FDR 5% 
hist(resd14.all.second_M_vs_A$pvalue, breaks=20, col="grey") # view histogram,  about 425 genes are likely false positives, alpha = 0.05 FDR is good!
abline(h=(nrow(dds.d14)*0.05),col="red") # add line at expected 5% false positive
table(resd14.all.second_M_vs_A$padj<0.05) # 0 DEGs
resd14.all.second_M_vs_A <- resd14.all.second_M_vs_A[order(resd14.all.second_M_vs_A$padj), ] ## Order by adjusted p-value
resdata.d14.second_M_vs_A  <- merge(as.data.frame(resd14.all.second_M_vs_A), as.data.frame(counts(dds.d14, normalized=TRUE)), by="row.names", sort=FALSE) ## Merge with normalized count data
names(resdata.d14.second_M_vs_A)[1] <- "Gene" # assign col 1

resd14.all.second_S_vs_A<- results(dds.d14, name="Second_Treament_S_vs_A", alpha = 0.05) # FDR 5% 
hist(resd14.all.second_S_vs_A$pvalue, breaks=20, col="grey") # view histogram,  about 425 genes are likely false positives, alpha = 0.05 FDR is good!
abline(h=(nrow(dds.d14)*0.05),col="red") # add line at expected 5% false positive
table(resd14.all.second_S_vs_A$padj<0.05) # 1 DEGs
resd14.all.second_S_vs_A <- resd14.all.second_S_vs_A[order(resd14.all.second_S_vs_A$padj), ] ## Order by adjusted p-value
resdata.d14.all.second_S_vs_A    <- merge(as.data.frame(resd14.all.second_S_vs_A), as.data.frame(counts(dds.d14, normalized=TRUE)), by="row.names", sort=FALSE) ## Merge with normalized count data
names(resdata.d14.all.second_S_vs_A)[1] <- "Gene" # assign col 1

## Write results
path_out.d14 = 'C:/Users/samjg/Documents/My_Projects/Pgenerosa_TagSeq_Metabolomics/TagSeq/RAnalysis/DESeq2/output/Day14/'
write.csv(resdata.d14.all.prim.M_vs_A, paste(path_out.d14, "Day14.PrimaryTreament_diffexpr-results.csv"))
write.csv(resdata.d14.second_M_vs_A, paste(path_out.d14, "Day14.SecondTreament_MvsA_diffexpr-results.csv"))
write.csv(resdata.d14.all.second_S_vs_A, paste(path_out.d14, "Day14.SecondTreament_SvsA_diffexpr-results.csv"))

# volcano plot ------------------------------------------------------------------------------------------------------ #

# Effect of primary treatment 'resd14.all.primary_M_vs_A' - 502 DEGs
png("RAnalysis/DESeq2/output/Day14/Day14.PrimaryTreatment-VolcanoPlot.png", 1000, 1000, pointsize=20)
EnhancedVolcano(resd14.all.primary_M_vs_A,
                lab = rownames(resd14.all.primary_M_vs_A),
                x = 'log2FoldChange',
                y = 'pvalue',
                title = 'Moderate versus Ambient Treatment (Primary)',
                subtitle = "DESeq2 - Differential expression",
                FCcutoff = 2.0,
                pointSize = 4.0,
                labSize = 6.0,
                colAlpha = 1,
                legendPosition = 'right',
                legendLabSize = 12,
                legendIconSize = 4.0,
                widthConnectors = 0.75)
dev.off()


# Effect of second treatment M v. A 'resdata.d14.second_M_vs_A' - 0 DEGs
png("RAnalysis/DESeq2/output/Day14/Day14.SecondTreatment_MvA-VolcanoPlot.png", 1000, 1000, pointsize=20)
EnhancedVolcano(resdata.d14.second_M_vs_A,
                lab = rownames(resdata.d14.second_M_vs_A),
                x = 'log2FoldChange',
                y = 'pvalue',
                title = 'Moderate versus Ambient Treatment (Second)',
                subtitle = "DESeq2 - Differential expression",
                FCcutoff = 2.0,
                pointSize = 4.0,
                labSize = 6.0,
                colAlpha = 1,
                legendPosition = 'right',
                legendLabSize = 12,
                legendIconSize = 4.0,
                widthConnectors = 0.75)
dev.off()

# Effect of second treatment S v. A 'resd14.all.second_S_vs_A' - 1 DEG
png("RAnalysis/DESeq2/output/Day14/Day14.SecondTreatment_SvA-VolcanoPlot.png", 1000, 1000, pointsize=20)
EnhancedVolcano(resd14.all.second_S_vs_A,
                lab = rownames(resd14.all.second_S_vs_A),
                x = 'log2FoldChange',
                y = 'pvalue',
                title = 'Severe versus Ambient Treatment (Second)',
                subtitle = "DESeq2 - Differential expression",
                FCcutoff = 2.0,
                pointSize = 4.0,
                labSize = 6.0,
                colAlpha = 1,
                legendPosition = 'right',
                legendLabSize = 12,
                legendIconSize = 4.0,
                widthConnectors = 0.75)
dev.off()

# Plot dispersions ------------------------------------------------------------------------------------------------------ #

png("RAnalysis/DESeq2/output/Day14/Day14-dispersions.png", 1000, 1000, pointsize=20)
plotDispEsts(dds.d14, main="Day 14 data dispersions")
dev.off()

# Data transformations  ------------------------------------------------------------------------------------------------------ #

# rlog - regularized log transformation of origin count data to log2 scale - fit for each sample and dist. of coefficients in the data
rld.d14<- rlogTransformation(dds.d14) # rlog transform (regularized log)
head(assays(rld.d14)) # view first few rows
hist(assay(rld.d14)) # view histogram 
meanSdPlot(assay(rld.d14)) # plot using "vsn" bioconductor package - shows the sd y axis (sq root of varaince in all samples) - flat curve may seem like a goals, BUT may be unreasonable in cases with MANY true DEGs from experimental conditions
# this gives log2(n + 1)
ntd.d14 <- normTransform(dds.d14)
head(assay(ntd.d14)) # view first few rows
hist(assay(ntd.d14)) # view histogram 
meanSdPlot(assay(ntd.d14)) # plot using "vsn" bioconductor package - shows the sd y axis (sq root of varaince in all samples) - flat curve may seem like a goals, BUT may be unreasonable in cases with MANY true DEGs from experimental conditions

# Plot heat map ------------------------------------------------------------------------------------------------------ #

##################  heat map of count matrix using "pheatmap" 

save_pheatmap <- function(x, filename, width=1000, height=960) { # template for saving pheatmap outputs
  stopifnot(!missing(x))
  stopifnot(!missing(filename))
  png(filename,width = width, height=height)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}

# settings for the heatmaps
select <- order(rowMeans(counts(dds.d14,normalized=TRUE)), decreasing=TRUE)[1:20] 
df <- as.data.frame(colData(dds.d14)[,c("Primary_Treatment","Second_Treament")])
annotation_colors = list(Treatment = c(A="Blue", M="Orange"))

# rlog heatmap
d14.rlog.heatmap<- pheatmap(assay(rld.d14)[select,], cluster_rows=T, show_rownames=F, # rlog heatmap
                           cluster_cols=T, annotation_col=df, annotation_colors = annotation_colors,
                           main = "rlog_Day14")
save_pheatmap(d14.rlog.heatmap, filename = "RAnalysis/DESeq2/output/Day14/Day14-rlog_heatmap.png")

# ntd heatmap
d14.ntd.heatmap<- pheatmap(assay(ntd.d14)[select,], cluster_rows=T, show_rownames=F, # rlog heatmap
                          cluster_cols=T, annotation_col=df, annotation_colors = annotation_colors,
                          main = "ntd_Day14")
save_pheatmap(d14.ntd.heatmap, filename = "RAnalysis/DESeq2/output/Day14/Day14-ntd_heatmap.png")

##################  heatmap of the sample-to-sample distances

# Colors for plots below; Use RColorBrewer, better and assign mycols variable
mycols.d14 <- brewer.pal(8, "Dark2")[1:length(unique(rld.d14$Primary_Treatment))]
mycols.d14_2 <- brewer.pal(8, "Dark2")[1:length(unique(rld.d14$Second_Treament))]
# Sample distance heatmap
sampleDists.d14 <- as.matrix(dist(t(assay(rld.d14))))
png("RAnalysis/DESeq2/output/Day14/Day14-heatmap.png", 1000, 1000, pointsize=20)
heatmap.2(as.matrix(sampleDists.d14), key=F, trace="none",
          col=colorpanel(100, "black", "white"),
          ColSideColors=mycols.d14_2[rld.d14$Second_Treament], RowSideColors=mycols.d14[rld.d14$Primary_Treatment],
          margin=c(10, 10), main="Day 14 Sample to Sample heatmap")
dev.off()

# Principal components analysis ------------------------------------------------------------------------------------------ #
pcaExplorer(dds = dds.d14)


# ========================================================== 
#
# DAY 21   (3 CPM in 50% samples)                          
# ========================================================== 

# dds.d21

path_out.d21 = 'C:/Users/samjg/Documents/My_Projects/Pgenerosa_TagSeq_Metabolomics/TagSeq/RAnalysis/DESeq2/output/Day21/'

# All treatment (Prim + Second + Third; ref level "A") =============================================================== #

dds.d21 <- DESeq(dds.d21) # wait for this to complete....
resultsNames(dds.d21) # view the names of your results model 

resd21.All.Primary_Treatment_M_vs_A <- results(dds.d21, name="Primary_Treatment_M_vs_A", alpha = 0.05) # FDR 5% 
hist(resd21.All.Primary_Treatment_M_vs_A$pvalue, breaks=20, col="grey") # view histogram,  about 425 genes are likely false positives, alpha = 0.05 FDR is good!
abline(h=(nrow(dds.d21)*0.05),col="red") # add line at expected 5% false positive
table(resd21.All.Primary_Treatment_M_vs_A$padj<0.05) #  221 DEGs
resd21.All.Primary_Treatment_M_vs_A <- resd21.All.Primary_Treatment_M_vs_A[order(resd21.All.Primary_Treatment_M_vs_A$padj), ] ## Order by adjusted p-value
resdata.d21.All.Primary_Treatment_M_vs_A  <- merge(as.data.frame(resd21.All.Primary_Treatment_M_vs_A), as.data.frame(counts(dds.d21, normalized=TRUE)), by="row.names", sort=FALSE) ## Merge with normalized count data
names(resdata.d21.All.Primary_Treatment_M_vs_A)[1] <- "Gene" # assign col 1


resd21.All.Second_Treatment_M_vs_A <- results(dds.d21, name="Second_Treament_M_vs_A", alpha = 0.05) # FDR 5% 
hist(resd21.All.Second_Treatment_M_vs_A$pvalue, breaks=20, col="grey") # view histogram,  about 425 genes are likely false positives, alpha = 0.05 FDR is good!
abline(h=(nrow(dds.d21)*0.05),col="red") # add line at expected 5% false positive
table(resd21.All.Second_Treatment_M_vs_A$padj<0.05) # 4 DEGs
resd21.All.Second_Treatment_M_vs_A <- resd21.All.Second_Treatment_M_vs_A[order(resd21.All.Second_Treatment_M_vs_A$padj), ] ## Order by adjusted p-value
resdata.d21.All.Second_Treatment_M_vs_A  <- merge(as.data.frame(resd21.All.Second_Treatment_M_vs_A), as.data.frame(counts(dds.d21, normalized=TRUE)), by="row.names", sort=FALSE) ## Merge with normalized count data
names(resdata.d21.All.Second_Treatment_M_vs_A)[1] <- "Gene" # assign col 1


resd21.All.Second_Treatment_S_vs_A <- results(dds.d21, name="Second_Treament_S_vs_A", alpha = 0.05) # FDR 5% 
hist(resd21.All.Second_Treatment_S_vs_A$pvalue, breaks=20, col="grey") # view histogram,  about 425 genes are likely false positives, alpha = 0.05 FDR is good!
abline(h=(nrow(dds.d21)*0.05),col="red") # add line at expected 5% false positive
table(resd21.All.Second_Treatment_S_vs_A$padj<0.05) # 0 DEGs
resd21.All.Second_Treatment_S_vs_A <- resd21.All.Second_Treatment_S_vs_A[order(resd21.All.Second_Treatment_S_vs_A$padj), ] ## Order by adjusted p-value
resdata.d21.All.Second_Treatment_S_vs_A <- merge(as.data.frame(resd21.All.Second_Treatment_S_vs_A), as.data.frame(counts(dds.d21, normalized=TRUE)), by="row.names", sort=FALSE) ## Merge with normalized count data
names(resdata.d21.All.Second_Treatment_S_vs_A)[1] <- "Gene" # assign col 1


resd21.All.Third_Treatment_M_vs_A <- results(dds.d21, name="Third_Treatment_M_vs_A", alpha = 0.05) # FDR 5% 
hist(resd21.All.Third_Treatment_M_vs_A$pvalue, breaks=20, col="grey") # view histogram,  about 425 genes are likely false positives, alpha = 0.05 FDR is good!
abline(h=(nrow(dds.d21)*0.05),col="red") # add line at expected 5% false positive
table(resd21.All.Third_Treatment_M_vs_A$padj<0.05) # 0 DEGs
resd21.All.Third_Treatment_M_vs_A <- resd21.All.Third_Treatment_M_vs_A[order(resd21.All.Third_Treatment_M_vs_A$padj), ] ## Order by adjusted p-value
resdata.d21.All.Third_Treatment_M_vs_A  <- merge(as.data.frame(resd21.All.Third_Treatment_M_vs_A), as.data.frame(counts(dds.d21, normalized=TRUE)), by="row.names", sort=FALSE) ## Merge with normalized count data
names(resdata.d21.All.Third_Treatment_M_vs_A)[1] <- "Gene" # assign col 1


## Write results
write.csv(resdata.d21.All.Primary_Treatment_M_vs_A, paste(path_out.d21, "Day21.PrimaryTreament_diffexpr-results.csv"))
write.csv(resdata.d21.All.Second_Treatment_M_vs_A, paste(path_out.d21, "Day21.SecondTreament_MvsA_diffexpr-results.csv"))
write.csv(resdata.d21.All.Second_Treatment_S_vs_A, paste(path_out.d21, "Day21.SecondTreament_SvsA_diffexpr-results.csv"))
write.csv(resdata.d21.All.Third_Treatment_M_vs_A, paste(path_out.d21, "Day21.ThirdTreament_diffexpr-results.csv"))

# volcano plot ------------------------------------------------------------------------------------------------------ #

# Effect of primary treatment 'resd14.all.primary_M_vs_A' - 221 DEGs
png("RAnalysis/DESeq2/output/Day21/Day21.PrimaryTreatment-VolcanoPlot.png", 1000, 1000, pointsize=20)
EnhancedVolcano(resd21.All.Primary_Treatment_M_vs_A,
                lab = rownames(resd21.All.Primary_Treatment_M_vs_A),
                x = 'log2FoldChange',
                y = 'pvalue',
                title = 'Moderate versus Ambient Treatment (Primary)',
                subtitle = "DESeq2 - Differential expression",
                FCcutoff = 2.0,
                pointSize = 4.0,
                labSize = 6.0,
                colAlpha = 1,
                legendPosition = 'right',
                legendLabSize = 12,
                legendIconSize = 4.0,
                widthConnectors = 0.75)
dev.off()


# Effect of second treatment M v. A 'resdata.d14.second_M_vs_A' - 0 DEGs
png("RAnalysis/DESeq2/output/Day21/Day21.SecondTreatment_MvA-VolcanoPlot.png", 1000, 1000, pointsize=20)
EnhancedVolcano(resd21.All.Second_Treatment_M_vs_A,
                lab = rownames(resd21.All.Second_Treatment_M_vs_A),
                x = 'log2FoldChange',
                y = 'pvalue',
                title = 'Moderate versus Ambient Treatment (Second)',
                subtitle = "DESeq2 - Differential expression",
                FCcutoff = 2.0,
                pointSize = 4.0,
                labSize = 6.0,
                colAlpha = 1,
                legendPosition = 'right',
                legendLabSize = 12,
                legendIconSize = 4.0,
                widthConnectors = 0.75)
dev.off()

# Effect of second treatment S v. A 'resd14.all.second_S_vs_A' - 1 DEG
png("RAnalysis/DESeq2/output/Day21/Day21.SecondTreatment_SvA-VolcanoPlot.png", 1000, 1000, pointsize=20)
EnhancedVolcano(resd21.All.Second_Treatment_S_vs_A,
                lab = rownames(resd21.All.Second_Treatment_S_vs_A),
                x = 'log2FoldChange',
                y = 'pvalue',
                title = 'Severe versus Ambient Treatment (Second)',
                subtitle = "DESeq2 - Differential expression",
                FCcutoff = 2.0,
                pointSize = 4.0,
                labSize = 6.0,
                colAlpha = 1,
                legendPosition = 'right',
                legendLabSize = 12,
                legendIconSize = 4.0,
                widthConnectors = 0.75)
dev.off()

# Effect of second treatment S v. A 'resd14.all.second_S_vs_A' - 1 DEG
png("RAnalysis/DESeq2/output/Day21/Day21.ThirdTreatment-VolcanoPlot.png", 1000, 1000, pointsize=20)
EnhancedVolcano(resd21.All.Third_Treatment_M_vs_A,
                lab = rownames(resd21.All.Third_Treatment_M_vs_A),
                x = 'log2FoldChange',
                y = 'pvalue',
                title = 'Moderate versus Ambient Treatment (Third)',
                subtitle = "DESeq2 - Differential expression",
                FCcutoff = 2.0,
                pointSize = 4.0,
                labSize = 6.0,
                colAlpha = 1,
                legendPosition = 'right',
                legendLabSize = 12,
                legendIconSize = 4.0,
                widthConnectors = 0.75)
dev.off()

# Plot dispersions ------------------------------------------------------------------------------------------------------ #

png("RAnalysis/DESeq2/output/Day21/Day21-dispersions.png", 1000, 1000, pointsize=20)
plotDispEsts(dds.d21, main="Day 21 data dispersions")
dev.off()

# Data transformations  ------------------------------------------------------------------------------------------------------ #

# rlog - regularized log transformation of origin count data to log2 scale - fit for each sample and dist. of coefficients in the data
rld.d21<- rlogTransformation(dds.d21) # rlog transform (regularized log)
head(assays(rld.d21)) # view first few rows
hist(assay(rld.d21)) # view histogram 
meanSdPlot(assay(rld.d21)) # plot using "vsn" bioconductor package - shows the sd y axis (sq root of varaince in all samples) - flat curve may seem like a goals, BUT may be unreasonable in cases with MANY true DEGs from experimental conditions
# this gives log2(n + 1)
ntd.d21 <- normTransform(dds.d21)
head(assay(ntd.d21)) # view first few rows
hist(assay(ntd.d21)) # view histogram 
meanSdPlot(assay(ntd.d21)) # plot using "vsn" bioconductor package - shows the sd y axis (sq root of varaince in all samples) - flat curve may seem like a goals, BUT may be unreasonable in cases with MANY true DEGs from experimental conditions

# Plot heat map ------------------------------------------------------------------------------------------------------ #

##################  heat map of count matrix using "pheatmap" 

save_pheatmap <- function(x, filename, width=1000, height=960) { # template for saving pheatmap outputs
  stopifnot(!missing(x))
  stopifnot(!missing(filename))
  png(filename,width = width, height=height)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}

# settings for the heatmaps
select <- order(rowMeans(counts(dds.d21,normalized=TRUE)), decreasing=TRUE)[1:20] 
df <- as.data.frame(colData(dds.d21)[,c("Primary_Treatment","Second_Treament", "Third_Treatment")])
annotation_colors = list(Treatment = c(A="Blue", M="Orange", S= "Red"))

# rlog heatmap
d21.rlog.heatmap<- pheatmap(assay(dds.d21)[select,], cluster_rows=T, show_rownames=F, # rlog heatmap
                            cluster_cols=T, annotation_col=df, annotation_colors = annotation_colors,
                            main = "rlog_Day21")
save_pheatmap(d21.rlog.heatmap, filename = "RAnalysis/DESeq2/output/Day21/Day21-rlog_heatmap.png")

# ntd heatmap
d21.ntd.heatmap<- pheatmap(assay(dds.d21)[select,], cluster_rows=T, show_rownames=F, # rlog heatmap
                           cluster_cols=T, annotation_col=df, annotation_colors = annotation_colors,
                           main = "ntd_Day21")
save_pheatmap(d21.ntd.heatmap, filename = "RAnalysis/DESeq2/output/Day21/Day21-ntd_heatmap.png")

##################  heatmap of the sample-to-sample distances

# Colors for plots below; Use RColorBrewer, better and assign mycols variable
mycols.d21 <- brewer.pal(8, "Dark2")[1:length(unique(rld.d21$Primary_Treatment))]
mycols.d21_2 <- brewer.pal(8, "Dark2")[1:length(unique(rld.d21$Second_Treament))]
# Sample distance heatmap
sampleDists.d21 <- as.matrix(dist(t(assay(rld.d21))))
png("RAnalysis/DESeq2/output/Day21/Day21-heatmap.png", 1000, 1000, pointsize=20)
heatmap.2(as.matrix(sampleDists.d21), key=F, trace="none",
          col=colorpanel(100, "black", "white"),
          ColSideColors=mycols.d21_2[rld.d21$Second_Treament], RowSideColors=mycols.d21[rld.d21$Primary_Treatment],
          margin=c(10, 10), main="Day 21 Sample to Sample heatmap")
dev.off()

# Principal components analysis ------------------------------------------------------------------------------------------ #
pcaExplorer(dds = dds.d21)


#====================================================================================================================== 
#
#                                       GO ANALYSIS (with WEGO)                   
#
#====================================================================================================================== 

# load and call desired columns of the annotation file (Uniprot IDs for each gene) - Roberts Lab (Sam White) compelted this annotation - open online on osf
annotation.txt <- read.delim2(file="Panopea-generosa-genes-annotations.txt", header=F)
annotation.df <- as.data.frame(annotation.txt)
annotation.df <- annotation.df %>% dplyr::select(c('V1','V8')) # call gene name and the GO terms - (Uniprot ID 'V5')
colnames(annotation.df)[1:2] <- c('Gene', 'GO.terms')
annotation.df$GO.terms <- gsub(";", " ", annotation.df$GO.terms) # remove the ; delimiter and replace with nothing - already tabbed 



#====================================================================================================================== 
# Day 0 primary treatment data  ====================================================================================== #
#====================================================================================================================== 

#====================================================================================================================== 
# Day 7 primary treatment data  ====================================================================================== #
#====================================================================================================================== 

resdata.d7.all.prim.M_vs_A # 94 DEGS (log2FoldChange > 1; < -1 & padj < 0.05) on day 14 in response to Primary treatment (initial Mdoerate vs. Ambient conditioning)
DEGS.d7.primary<- resdata.d7.all.prim.M_vs_A %>%  dplyr::filter(padj<0.05)
nrow(DEGS.d7.primary) # 94 - we have all DEGs now

UPREG.d7.primary <- DEGS.d7.primary %>%  dplyr::filter(log2FoldChange > 1) # call upregulated genes
DWNREG.d7.primary <- DEGS.d7.primary %>%  dplyr::filter(log2FoldChange < 1) # call downregulated genes 
nrow(UPREG.d7.primary) + nrow(DWNREG.d7.primary) # should be equal to 94


D7.UPREG_GO <- merge(UPREG.d7.primary, annotation.df, by = "Gene")
D7.UPREG_WEGO <- D7.UPREG_GO %>% dplyr::select(c('Gene', 'GO.terms'))
D7.UPREG_WEGO$GO.terms[is.na(D7.UPREG_WEGO$GO.terms)] <- " " # make NA blank
# View(D7.UPREG_WEGO) # view the data

D7.DWNREG_GO <- merge(DWNREG.d7.primary, annotation.df, by = "Gene")
D7.DWNREG_WEGO <- D7.DWNREG_GO %>% dplyr::select(c('Gene', 'GO.terms'))
D7.DWNREG_WEGO$GO.terms[is.na(D7.DWNREG_WEGO$GO.terms)] <- " " # make NA blank
# View(D7.DWNREG_WEGO) # view the data

# write to GO folder
write.table(D7.UPREG_WEGO, file = "RAnalysis/GO/Day7_Upreg_Primary_WEGO", quote=FALSE, row.names = FALSE, col.names = FALSE) # call at 'Native Format' on  https://wego.genomics.cn/ 
write.table(D7.DWNREG_WEGO, file = "RAnalysis/GO/Day7_Downreg_Primary_WEGO", quote=FALSE, row.names = FALSE, col.names = FALSE) # call at 'Native Format' on  https://wego.genomics.cn/ 



#====================================================================================================================== 
# Day 14  primary treatment data  ====================================================================================== #
#====================================================================================================================== 

resdata.d14.all.prim.M_vs_A # 502 DEGS (log2FoldChange > 1; < -1 & padj < 0.05) on day 14 in response to Primary treatment (initial Mdoerate vs. Ambient conditioning)
DEGS.d14.primary<- resdata.d14.all.prim.M_vs_A %>%  dplyr::filter(padj<0.05)
nrow(DEGS.d14.primary) # 502 - we have all DEGs now

UPREG.d14.primary <- DEGS.d14.primary %>%  dplyr::filter(log2FoldChange > 1) # call upregulated genes
DWNREG.d14.primary <- DEGS.d14.primary %>%  dplyr::filter(log2FoldChange < 1) # call downregulated genes 
nrow(UPREG.d14.primary) + nrow(DWNREG.d14.primary) # should be equal to 502


D14.UPREG_GO <- merge(UPREG.d14.primary, annotation.df, by = "Gene")
D14.UPREG_WEGO <- D14.UPREG_GO %>% dplyr::select(c('Gene', 'GO.terms'))
D14.UPREG_WEGO$GO.terms[is.na(D14.UPREG_WEGO$GO.terms)] <- " " # make NA blank
# View(D14.UPREG_WEGO) # view the data

D14.DWNREG_GO <- merge(DWNREG.d14.primary, annotation.df, by = "Gene")
D14.DWNREG_WEGO <- D14.DWNREG_GO %>% dplyr::select(c('Gene', 'GO.terms'))
D14.DWNREG_WEGO$GO.terms[is.na(D14.DWNREG_WEGO$GO.terms)] <- " " # make NA blank
# View(D14.DWNREG_WEGO) # view the data

# write to GO folder
write.table(D14.UPREG_WEGO, file = "RAnalysis/GO/Day14_Upreg_Primary_WEGO", quote=FALSE, row.names = FALSE, col.names = FALSE) # call at 'Native Format' on  https://wego.genomics.cn/ 
write.table(D14.DWNREG_WEGO, file = "RAnalysis/GO/Day14_Downreg_Primary_WEGO", quote=FALSE, row.names = FALSE, col.names = FALSE) # call at 'Native Format' on  https://wego.genomics.cn/ 


#====================================================================================================================== 
# Day 21 primary treatment data  ====================================================================================== #
#====================================================================================================================== 

resdata.d21.All.Primary_Treatment_M_vs_A # 221 DEGS (log2FoldChange > 1; < -1 & padj < 0.05) on day 21 in response to Primary treatment (initial Mdoerate vs. Ambient conditioning)
DEGS.d21.primary<- resdata.d21.All.Primary_Treatment_M_vs_A %>%  dplyr::filter(padj<0.05)
nrow(DEGS.d21.primary) # 221 - we have all DEGs now

UPREG.d21.primary <- DEGS.d21.primary %>%  dplyr::filter(log2FoldChange > 1) # call upregulated genes
DWNREG.d21.primary <- DEGS.d21.primary %>%  dplyr::filter(log2FoldChange < 1) # call downregulated genes 
nrow(UPREG.d21.primary) + nrow(DWNREG.d21.primary) # should be equal to 221


UPREG_GO <- merge(UPREG.d21.primary, annotation.df, by = "Gene")
UPREG_WEGO <- UPREG_GO %>% dplyr::select(c('Gene', 'GO.terms'))
UPREG_WEGO$GO.terms[is.na(UPREG_WEGO$GO.terms)] <- " " # make NA blank
View(UPREG_WEGO)

DWNREG_GO <- merge(DWNREG.d21.primary, annotation.df, by = "Gene")
DWNREG_WEGO <- DWNREG_GO %>% dplyr::select(c('Gene', 'GO.terms'))
DWNREG_WEGO$GO.terms[is.na(DWNREG_WEGO$GO.terms)] <- " " # make NA blank
View(DWNREG_WEGO)

# write to GO folder
write.table(UPREG_WEGO, file = "RAnalysis/GO/Day21_Upreg_Primary_WEGO", quote=FALSE, row.names = FALSE, col.names = FALSE)
write.table(DWNREG_WEGO, file = "RAnalysis/GO/Day21_Downreg_Primary_WEGO", quote=FALSE, row.names = FALSE, col.names = FALSE)

