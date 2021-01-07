---
# title: "Geoduck_TagSeq_DESeq2"
# author: "Samuel Gurr"
# date: "December 29, 2020"
---
  
# LOAD PACKAGES
library(DESeq2) # note: this was previously installed with the command `BiocManager::install("DESeq2")`
library(ggplot2)
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

# SET WORKING DIRECTORY AND LOAD DATA
setwd("C:/Users/samjg/Documents/My_Projects/Pgenerosa_TagSeq_Metabolomics/TagSeq/")

# LOAD DATA
cts <- read.csv(file="HPC_Bioinf/outputs/transcript_count_matrix.csv", sep=',', header=TRUE) # read the output count matrix from prepDE.py

UT_seq_map <- read.csv(file="20201020_Gurr_TagSeq_UTAustin.csv", sep=',', header=TRUE)
smpl_ref <- read.csv(file="Sample_reference.csv", sep=',', header=TRUE)
treatment_ref <- read.csv(file="Extraction_checklist.csv", sep=',', header=TRUE)

#=====================================================================================
#                             FORMAT EXPERIMENT DESIGN DATAFRAME
#=====================================================================================
# MASTER REFERENCE DATA.FRAME
# format and merge to buld master reference dataframe
smpl_ref$Seq_Pos <- paste(smpl_ref$ï..TagSeq_Plate, smpl_ref$TagSeq_Well, sep="_")
smpl_ref <- smpl_ref %>% select(-c("ï..TagSeq_Plate", "TagSeq_Well"))
UT_seq_map$Seq_Pos  <- paste(UT_seq_map$ï..Plate, UT_seq_map$Well, sep="_")
UT_seq_map <- UT_seq_map %>% select(-c("ï..Plate", "Well"))
Seq.Ref <- merge(smpl_ref, UT_seq_map, by = "Seq_Pos")
treatment_ref <- treatment_ref[, c(3:8)]
Mstr.Ref <- merge(Seq.Ref, treatment_ref, by = "Geoduck_ID")

#======================================== #
# ALL TIMEPOINTS:
# call all experiment design treatments as 'exp.data'
exp.data <- Mstr.Ref[,c("Sample.Name","All_Treatment", "Primary_Treatment", "Second_Treament", "Third_Treatment", "Time")]
#======================================== #
# DAY 0 
exp.data.d0 <- exp.data %>% dplyr::filter(Time %in% 'Day0')
nrow(exp.data.d0) # 8 total samples on Day 0
#======================================== #
# DAY 7 
exp.data.d7 <- exp.data %>% dplyr::filter(Time %in% 'Day7')
nrow(exp.data.d7) # 36 total samples on Day 0
#======================================== #
# DAY 14 
exp.data.d14 <- exp.data %>% dplyr::filter(Time %in% 'DAY14')
nrow(exp.data.d14) # 36 total samples on Day 0
#======================================== #
# DAY 21 
exp.data.d21 <- exp.data %>% dplyr::filter(Time %in% 'DAY21')
nrow(exp.data.d21) # 62 total sampels on day 21



#=====================================================================================  #
#                             FORMAT COUNT MATRIX 
#=====================================================================================
# About: need to call the cts.matrix and call only samples that match the IDs at day 0 and 21 
# Why? DESeq2 can only apply to factors with >=2 levels, so we cannot address 'Day' at just "Day21" for example
# alternatively, we can call subset matrices for the timepoints of interest

# NOTE: this cts matrix was NOT merged by lanes, run'rowsum' on unique delimiters to merge 'cts.merged'
#  Format count matrix (cts has 2x columns per sample on different lanes - sum counts from unique columns + build matrix)
ncol(cts) # 282 samples (not counting gene ID column) - should be 141 samples, need to sum columns by unique ID
cts.merged <- data.frame(cts[,-1], row.names=cts[,1]) # call new dataframe with first column now as row names, now all row values are numeric
names(cts.merged) <- sapply(strsplit(names(cts.merged), "_"), '[', 1) # split the column names by "_" delimiter  and call the first field SG##
cts.merged <- t(rowsum(t(cts.merged), group = colnames(cts.merged), na.rm = T)) # merge all unique columns and sum counts 
ncol(cts.merged)

#======================================== #
# ALL TIMEPOINTS
cts.matrix <-as.matrix(cts.merged, row.names="transcript_id") # call dataframe as matrix
ncol(cts.matrix) # 141 samples
nrow(cts.matrix) # 34947 total genes

cts.merged.as.table <- data.frame(transcript_id = row.names(cts.merged), cts.merged) # add back the rownames 'transcript_ID'
rownames(cts.merged.as.table) <- NULL # ommit the rownames
ncol(cts.merged.as.table) # 142 counting the ID that we want to keep! 

#======================================== #
# DAY 0 
# About: run dyplr 'antijoin' to call cts columns that match 'Sample.Name' in the data frame 'exp.data.d0'
cts.merged.d0 <- cts.merged.as.table[,c(1,na.omit(match(exp.data.d0$Sample.Name, colnames(cts.merged.as.table))))]
cts.merged.d0 <- data.frame(cts.merged.d0[,-1], row.names=cts.merged.d0[,1])
cts.matrix.d0  <-as.matrix(cts.merged.d0, row.names="transcript_id")
ncol(cts.matrix.d0) # 8  samples from just Day 0
colnames(cts.matrix.d0) == exp.data.d0$Sample.Name # check if TRUE, means the same as the exp/design dataframe exp.data.d0
# pre-filtering; genes ommitted if < 3 counts per million reads in 50% of samples
colSums(cts.matrix.d0) # view the colSums of our Day0 samples  - notice the read sums are around 1 million
CPM.d0 <- cpm(cts.matrix.d0) # Obtain CPMs (counts oer million) using egdeR
head(CPM.d0) # Have a look at the output
thresh.d0 <- CPM.d0 > 10 # Which values in myCPM are greater than 3?
head(thresh.d0) # This produces a logical matrix with TRUEs and FALSES
rowSums(head(thresh.d0)) # Summary of how many TRUEs there are in each row
table(rowSums(thresh.d0)) # 9631 genes with TRUE in all 8 samples 
keep.d0 <- rowSums(thresh.d0) >= (ncol(thresh.d0)/2) # we would like to keep genes that have at least 50% TRUES in each row of thresh
summary(keep.d0) # FALSE 19911 & TRUE 15036 -- more than half of the genes did not pass
cts.matrix.d0.filtered <- cts.matrix.d0[keep.d0,] # Subset the rows of countdata to keep the more highly expressed genes
dim(cts.matrix.d0.filtered) # 9091 genes & 8 samples
head(cts.matrix.d0.filtered) # FINAL DAY 0 DATASET
# We will look at a cople samples to check our defined threshold.. 
plot(CPM.d0[,1], cts.matrix.d0[,1], xlab="CPM", ylab="Raw Count", ylim=c(0,50), xlim=c(0,50), main=colnames(CPM.d0)[1]) # check if our threshold of 3 CPM in 50% samples corresponds to a count of about 10-15
abline(v=10) 
abline(h=10)
plot(CPM.d0[,5], cts.matrix.d0[,5], xlab="CPM", ylab="Raw Count", ylim=c(0,50), xlim=c(0,50), main=colnames(CPM.d0)[1]) # check if our threshold of 3 CPM in 50% samples corresponds to a count of about 10-15
abline(v=10) 
abline(h=10)

#======================================== #
# DAY 7
cts.merged.d7 <- cts.merged.as.table[,c(1,na.omit(match(exp.data.d7$Sample.Name, colnames(cts.merged.as.table))))]
cts.merged.d7 <- data.frame(cts.merged.d7[,-1], row.names=cts.merged.d7[,1])
cts.matrix.d7  <-as.matrix(cts.merged.d7, row.names="transcript_id")
ncol(cts.matrix.d7) # 36 samples from just Day 7
colnames(cts.matrix.d7) == exp.data.d7$Sample.Name # check if TRUE, means the same as the exp/design dataframe exp.data.d7
# pre-filtering; genes ommitted if < 3 counts per million reads in 50% of samples
colSums(cts.matrix.d7) # view the colSums of our Day7 samples 
CPM.d7 <- cpm(cts.matrix.d7) # Obtain CPMs (counts oer million) using egdeR
head(CPM.d7) # Have a look at the output
thresh.d7 <- CPM.d7 > 10 # Which values in myCPM are greater than 3?
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


#======================================== #
# DAY 14 
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
thresh.d14 <- CPM.d14 > 10 # Which values in myCPM are greater than 3?
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


#======================================== #
# DAY 21 
cts.merged.d21 <- cts.merged.as.table[,c(1,na.omit(match(exp.data.d21$Sample.Name, colnames(cts.merged.as.table))))]
cts.merged.d21 <- data.frame(cts.merged.d21[,-1], row.names=cts.merged.d21[,1])
cts.matrix.d21  <-as.matrix(cts.merged.d21, row.names="transcript_id")
ncol(cts.matrix.d21) # # 62 total sampels on day 21
colnames(cts.matrix.d21) == exp.data.d21$Sample.Name # check if TRUE, means the same as the exp/design dataframe exp.data.d21
# pre-filtering; genes ommitted if < 3 counts per million reads in 50% of samples
colSums(cts.matrix.d21) # view the colSums of our Day21 samples - notice SOME counts are near 1 million
CPM.d21 <- cpm(cts.matrix.d21) # Obtain CPMs (counts oer million) using egdeR
head(CPM.d21) # Have a look at the output
thresh.d21 <- CPM.d21 > 10 # Which values in myCPM are greater than 3?
head(thresh.d21) # This produces a logical matrix with TRUEs and FALSES
rowSums(head(thresh.d21)) # Summary of how many TRUEs there are in each row
table(rowSums(thresh.d21)) # 3451 genes with TRUE in all 35 samples 
keep.d21 <- rowSums(thresh.d21) >= (ncol(thresh.d21)/2) # we would like to keep genes that have at least 50% TRUES in each row of thresh
summary(keep.d21) # FALSE 26622 & TRUE 8325 -- more than three quarters of the genes did not pass
cts.matrix.d21.filtered <- cts.matrix.d21[keep.d21,] # Subset the rows of countdata to keep the more highly expressed genes
dim(cts.matrix.d21.filtered) # 8541 genes & 35 samples
head(cts.matrix.d21.filtered) # FINAL DAY 21 DATASET
# We will look at a cople samples to check our defined threshold.. 
plot(CPM.d21[,60], cts.matrix.d21[,60], xlab="CPM", ylab="Raw Count", ylim=c(0,50), xlim=c(0,50), main=colnames(CPM.d21)[1]) # check if our threshold of 3 CPM in 50% samples corresponds to a count of about 10-15
abline(v=10) 
abline(h=10)
plot(CPM.d21[,5], cts.matrix.d21[,5], xlab="CPM", ylab="Raw Count", ylim=c(0,50), xlim=c(0,50), main=colnames(CPM.d21)[1]) # check if our threshold of 3 CPM in 50% samples corresponds to a count of about 10-15
abline(v=10) 
abline(h=10)




#=====================================================================================
#  PREP dds files for DESeq2    ~ 'Primary_Treatment' 
#=====================================================================================

#======================================================================================= #
# ALL TIMEPOINTS and ALL TREAMENTS:
#---- INPUT: Experiment/design matrix == 'exp.data' ...formatted to 'exp.data.ALL.TREATMENTS.mtx'
#---- INPUT: TagSeq count matrix == 'cts.matrix'


# format Experiment/design dataframe into matrix
exp.data$All_Treatment <- factor(exp.data$All_Treatment) # change All_Treatment to factor
exp.data$Time <- factor(exp.data$Time) # change Day to a factor
exp.data.ALL.TREATMENTS <- exp.data %>% dplyr::select(c('Sample.Name', 'Time', 'All_Treatment')) # condense dataset to build target matrix
# note: day 7 and day 14 have the same treatment -  line below substr and pastes to 'All_Treatment'
exp.data.ALL.TREATMENTS$All_Treatment <- paste(exp.data.ALL.TREATMENTS$All_Treatment, (substr(exp.data.ALL.TREATMENTS$Time, 4,5)), sep=".d") # substr and add the day to treatment to address separately downstream DESeq2 compairons
exp.data.ALL.TREATMENTS <- data.frame(exp.data.ALL.TREATMENTS[,-1], row.names=exp.data.ALL.TREATMENTS[,1]) # move Sample.Name column as row names  
exp.data.ALL.TREATMENTS.mtx <- as.matrix(exp.data.ALL.TREATMENTS, row.names="Geoduck.ID") # create matrix 
# fix(exp.data.mtx) # view data - remove # to open
# CKECK THE cts.matrix AND THE exp.data.mtx
exp.data.ALL.TREATMENTS.mtx <- exp.data.ALL.TREATMENTS.mtx[match(colnames(cts.merged),rownames(exp.data.ALL.TREATMENTS.mtx)), ]
all(rownames(exp.data.ALL.TREATMENTS.mtx) %in% colnames(cts.matrix)) # should be TRUE
all(rownames(exp.data.ALL.TREATMENTS.mtx) == colnames(cts.matrix)) # should be TRUE
cts.matrix <- cts.matrix[, rownames(exp.data.ALL.TREATMENTS.mtx)]
all(rownames(exp.data.ALL.TREATMENTS.mtx) == colnames(cts.matrix))  # should be TRUE

# build dds
dds.all.treatments <- DESeqDataSetFromMatrix(countData = cts.matrix,
                                      colData = exp.data.ALL.TREATMENTS.mtx,
                                      design = ~ All_Treatment) # DESeq Data Set (dds) - design as~ Primary_Treatment + Day
dds.all.treatments # view dds

# prep for DESeq2
dds.all.treatments$All_Treatment <- relevel(dds.all.treatments$All_Treatment, ref = "A") # specify the reference level for count analysis - A = the control treatment
levels(dds.all.treatments$Time) # NULL 
levels(dds.all.treatments$All_Treatment) #  levels for condition 'Treatment': "A"   "AA"  "AAA" "AAM" "AM"  "AMA" "AMM" "AS"  "ASA" "ASM" "M"   "MA"  "MAA" "MAM" "MM"  "MMA" "MMM" "MS"  "MSA" "MSM"
keep <- rowSums(counts(dds.all.treatments)) >= 5 # PRE-FILTERING... ommits genes where counts between ALL samples (as rows) less than 5 
dds.all.treatments <- dds.all.treatments[keep,] # integrate this criteria in the data 
design(dds.all.treatments) # view the design we have specified 
ncol(dds.all.treatments) # 141 cols - Good
colData(dds.primary)

# DEG ANALYSIS READY!!

#======================================================================================= #
# ALL TIMEPOINTS and PRIMARY TREAMENT:
#---- INPUT: Experiment/design matrix == 'exp.data' ...formatted to 'exp.data.PRIMARY.mtx'
#---- INPUT: TagSeq count matrix == 'cts.matrix'

# format Experiment/design dataframe into matrix
exp.data$Primary_Treatment <- factor(exp.data$Primary_Treatment) # change Primary_Treatment to factor
exp.data$Day <- factor(exp.data$Day) # change Day to a factor
exp.data.PRIMARY <- exp.data %>% dplyr::select(c('Sample.Name', 'Time', 'Primary_Treatment')) # coondense dataset to build target matrix
exp.data.PRIMARY <- data.frame(exp.data.PRIMARY[,-1], row.names=exp.data.PRIMARY[,1]) # move Sample.Name column as row names  
exp.data.PRIMARY.mtx <- as.matrix(exp.data.PRIMARY, row.names="Geoduck.ID") # create matrix 
# fix(exp.data.mtx) # view data - remove # to open
# CKECK THE cts.matrix AND THE exp.data.mtx
exp.data.PRIMARY.mtx <- exp.data.PRIMARY.mtx[match(colnames(cts.merged),rownames(exp.data.PRIMARY.mtx)), ]
all(rownames(exp.data.PRIMARY.mtx) %in% colnames(cts.matrix)) # should be TRUE
all(rownames(exp.data.PRIMARY.mtx) == colnames(cts.matrix)) # should be TRUE
cts.matrix <- cts.matrix[, rownames(exp.data.PRIMARY.mtx)]
all(rownames(exp.data.PRIMARY.mtx) == colnames(cts.matrix))  # should be TRUE

# build dds
dds.primary <- DESeqDataSetFromMatrix(countData = cts.matrix,
                                         colData = exp.data.PRIMARY.mtx,
                                         design = ~ Primary_Treatment + Time) # DESeq Data Set (dds) - design as~ Primary_Treatment + Day
dds.primary # view dds

# prep for DESeq2
dds.primary$Primary_Treatment <- relevel(dds.primary$Primary_Treatment, ref = "A") # specify the reference level for count analysis - A = the control treatment
levels(dds.primary$Time) # levels for condition 'Time': "Day0"  "DAY14" "DAY21" "Day7" 
levels(dds.primary$Primary_Treatment) #  levels for condition 'Treatment': "A" "M"
keep <- rowSums(counts(dds.primary)) >= 5 # PRE-FILTERING... ommits genes where counts between ALL samples (as rows) less than 5 
dds.primary <- dds.primary[keep,] # integrate this criteria in the data 
design(dds.primary) # view the design we have specified 
ncol(dds.primary) # 141 cols - Good
colData(dds.primary)

# DEG ANALYSIS READY!!

#======================================================================================= #
# DAY 0 
#---- INPUT: Experiment/design matrix == 'exp.data.d0' ...formatted to 'exp.data.d0.PRIMARY.mtx'
#---- INPUT: TagSeq count matrix == 'cts.matrix.d0'

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
dds.primary.d0 <- DESeqDataSetFromMatrix(countData = cts.matrix.d0.filtered,
                                      colData = exp.data.d0.PRIMARY.mtx,
                                      design = ~ Primary_Treatment) # DESeq Data Set (dds) - design as ~Primary_Treatment
dds.primary.d0 # view dds

# prep for DESeq2
dds.primary.d0$Primary_Treatment <- relevel(dds.primary.d0$Primary_Treatment, ref = "A") # specify the reference level for count analysis - A = the control treatment
dds.primary.d0$Time <- droplevels(dds.primary.d0$Time)
levels(dds.primary.d0$Time) # NULL
levels(dds.primary.d0$Primary_Treatment) #  levels for condition 'Treatment': "A" "M"

design(dds.primary.d0) # view the design we have specified 
ncol(dds.primary.d0) # 8 cols - Good
nrow(dds.primary.d0) # 9091 cols - Good
colData(dds.primary.d0)

# DEG ANALYSIS READY!!

#======================================================================================= #
# DAY 21 
#---- INPUT: Experiment/design matrix == 'exp.data.d21' ...formatted to 'exp.data.d21.PRIMARY.mtx'
#---- INPUT: TagSeq count matrix == 'cts.matrix.d21'

# format Experiment/design dataframe into matrix
exp.data.d21.PRIMARY <- exp.data.d21 %>% dplyr::select(c('Sample.Name', 'Primary_Treatment', 'Time')) # coondense dataset to build target matrix
exp.data.d21.PRIMARY$Primary_Treatment <- factor(exp.data.d21.PRIMARY$Primary_Treatment) # change Primary_Treatment to factor
exp.data.d21.PRIMARY <- data.frame(exp.data.d21.PRIMARY[,-1], row.names=exp.data.d21.PRIMARY[,1]) # move Sample.Name column as row names  
exp.data.d21.PRIMARY.mtx <- as.matrix(exp.data.d21.PRIMARY, row.names="Geoduck.ID") # create matrix 
# fix(exp.data.mtx) # view data - remove # to open
# CKECK THE cts.matrix AND THE exp.data.mtx
exp.data.d21.PRIMARY.mtx <- exp.data.d21.PRIMARY.mtx[match(colnames(cts.merged.d21),rownames(exp.data.d21.PRIMARY.mtx)), ]
all(rownames(exp.data.d21.PRIMARY.mtx) %in% colnames(cts.matrix.d21)) # should be TRUE
all(rownames(exp.data.d21.PRIMARY.mtx) == colnames(cts.matrix.d21)) # should be TRUE
cts.matrix.d21 <- cts.matrix.d21[, rownames(exp.data.d21.PRIMARY.mtx)]
all(rownames(exp.data.d21.PRIMARY.mtx) == colnames(cts.matrix.d21))  # should be TRUE

# build dds
dds.primary.d21 <- DESeqDataSetFromMatrix(countData = cts.matrix.d21,
                                         colData = exp.data.d21.PRIMARY.mtx,
                                         design = ~ Primary_Treatment) # DESeq Data Set (dds) - design as ~Primary_Treatment
dds.primary.d21 # view dds

# prep for DESeq2
dds.primary.d21$Primary_Treatment <- relevel(dds.primary.d21$Primary_Treatment, ref = "A") # specify the reference level for count analysis - A = the control treatment
dds.primary.d21$Time <- droplevels(dds.primary.d21$Time)
levels(dds.primary.d21$Time) # NULL
levels(dds.primary.d21$Primary_Treatment) #  levels for condition 'Treatment': "A" "M"
keep <- rowSums(counts(dds.primary.d21)) >= 5 # PRE-FILTERING... ommits genes where counts between ALL samples (as rows) less than 5 
colSums(counts(dds.primary.d21))


ncol(dds.primary.d21)
nrow(dds.primary.d21)

dds.primary.d21 <- dds.primary.d21[keep,] # integrate this criteria in the data 
design(dds.primary.d21) # view the design we have specified 
ncol(dds.primary.d21) # 62 cols - Good
colData(dds.primary.d21)

# DEG ANALYSIS READY!!



#======================================================================================= #
# DAY 21 
#---- INPUT: Experiment/design matrix == 'exp.data.d21' ...formatted to 'exp.data.d21.PRIMARY.mtx'
#---- INPUT: TagSeq count matrix == 'cts.matrix.d21'

# format Experiment/design dataframe into matrix
exp.data.d21.ALL <- exp.data.d21 %>% dplyr::select(c('Sample.Name', 'Primary_Treatment', 'Second_Treament', 'Third_Treatment')) # coondense dataset to build target matrix
exp.data.d21.ALL$Primary_Treatment <- factor(exp.data.d21.ALL$Primary_Treatment) # change Primary_Treatment to factor
exp.data.d21.ALL$Second_Treament <- factor(exp.data.d21.ALL$Second_Treament) # change Second_Treament to factor
exp.data.d21.ALL$Third_Treatment <- factor(exp.data.d21.ALL$Third_Treatment) # change Third_Treatment to factor
exp.data.d21.ALL <- data.frame(exp.data.d21.ALL[,-1], row.names=exp.data.d21.ALL[,1]) # move Sample.Name column as row names  
exp.data.d21.ALL.mtx <- as.matrix(exp.data.d21.ALL, row.names="Geoduck.ID") # create matrix 
# fix(exp.data.mtx) # view data - remove # to open
# CKECK THE cts.matrix AND THE exp.data.mtx
exp.data.d21.ALL.mtx <- exp.data.d21.ALL.mtx[match(colnames(cts.merged.d21),rownames(exp.data.d21.ALL.mtx)), ]
all(rownames(exp.data.d21.ALL.mtx) %in% colnames(cts.matrix.d21)) # should be TRUE
all(rownames(exp.data.d21.ALL.mtx) == colnames(cts.matrix.d21)) # should be TRUE
cts.matrix.d21 <- cts.matrix.d21[, rownames(exp.data.d21.ALL.mtx)]
all(rownames(exp.data.d21.ALL.mtx) == colnames(cts.matrix.d21))  # should be TRUE

# build dds
dds.All.d21 <- DESeqDataSetFromMatrix(countData = cts.matrix.d21,
                                          colData = exp.data.d21.ALL.mtx,
                                          design = ~ Primary_Treatment + Second_Treament + Third_Treatment) # DESeq Data Set (dds) - design as ~Primary_Treatment
dds.All.d21 # view dds

# prep for DESeq2
dds.All.d21$Primary_Treatment <- relevel(dds.All.d21$Primary_Treatment, ref = "A") # specify the reference level for count analysis - A = the control treatment
dds.All.d21$Second_Treament <- relevel(dds.All.d21$Second_Treament, ref = "A") # specify the reference level for count analysis - AA = the control treatment
dds.All.d21$Third_Treatment <- relevel(dds.All.d21$Third_Treatment, ref = "A") # specify the reference level for count analysis - AA = the control treatment
keep <- rowSums(counts(dds.All.d21)) >= 5 # PRE-FILTERING... ommits genes where counts between ALL samples (as rows) less than 5 
dds.All.d21 <- dds.All.d21[keep,] # integrate this criteria in the data 
design(dds.All.d21) # view the design we have specified 
ncol(dds.All.d21) # 62 cols - Good
colData(dds.All.d21) # looks good

# DEG ANALYSIS READY!!


#=====================================================================================
# Differential expression analysis
#=====================================================================================
# About: Log fold change shrinkage for visualization and ranking
# Run DESeq : Modeling counts with the specified design(dds) treatment effects 

#======================================================================================= #
# ALL TREATMENTS AND TIMEPOINTS

dds.all.treatments <- DESeq(dds.all.treatments) # wait for this to complete....
resultsNames(dds.all.treatments) # view the names of your results model 

# Get differential expression results
res.ALL <- results(dds.all.treatments)
table(res.ALL$padj<0.05)
## Order by adjusted p-value
res.ALL <- res.ALL[order(res.ALL$padj), ]
## Merge with normalized count data
resdata.ALL <- merge(as.data.frame(res.ALL), as.data.frame(counts(dds.all.treatments, normalized=TRUE)), by="row.names", sort=FALSE)
names(resdata.ALL)[1] <- "Gene"
head(resdata.ALL)
## Write results
write.csv(resdata.d0.primary, file="D0.PrimaryTreatment_diffexpr-results.csv")
## Examine plot of p-values
hist(resd0.primary$pvalue, breaks=50, col="grey") # view histogram


res.test <- results(dds.all.treatments, name="All_Treatment_M.d0_vs_A.d0")
table(res.test.contrasts$padj<0.05)
res.test.contrasts <- results(dds.all.treatments, contrast=c("All_Treatment","M.d0","A.d0")) # alternative way of calling the results (i.e. Day_DAY21_vs_Day0)
res.d0_Ordered <- res.d0[order(res.d0$pvalue),] # how to order based on p-valuel 
res.d0_Ordered # ordered based on pvalue
sum(res.d0_Ordered$padj < 0.1, na.rm=TRUE) # 27 marginal DEGs
sum(res.d0_Ordered$padj < 0.05, na.rm=TRUE) # 14 significant DEGs
summary(res.d0)
# Plot dispersions
png("Day0.PrimaryTreatment-dispersions.png", 1000, 1000, pointsize=20)
plotDispEsts(dds.primary.d0, main="Dispersion plot_M_vs._A_Day0")
dev.off()
# Regularized log transformation for clustering/heatmaps, etc
rld.d0.primary <- rlogTransformation(dds.primary.d0) # rlog transform (regularized log)
head(assay(rld.d0.primary)) # view first few rows
hist(assay(rld.d0.primary)) # view histogram 
# Colors for plots below; Use RColorBrewer, better and assign mycols variable
mycols <- brewer.pal(8, "Dark2")[1:length(unique(rld.d0.primary$Primary_Treatment))]
# Sample distance heatmap
sampleDists <- as.matrix(dist(t(assay(rld.d0.primary))))
png("Day0.PrimaryTreatment-heatmap-samples.png", w=1000, h=1000, pointsize=20)
heatmap.2(as.matrix(sampleDists), key=F, trace="none",
          col=colorpanel(100, "black", "white"),
          ColSideColors=mycols[rld.d0.primary$Primary_Treatment], RowSideColors=mycols[rld.d0.primary$Primary_Treatment],
          margin=c(10, 10), main="Sample Distance Matrix")
dev.off()

# ## Examine independent filtering
# attr(resd0.primary, "filterThreshold")
# plot(attr(resd0.primary,"filterNumRej"), type="b", xlab="quantiles of baseMean", ylab="number of rejections")

# Principal components analysis ------------------------------------------------------------------------------------------------------ #
## Could do with built-in DESeq2 function:
## DESeq2::plotPCA(rld, intgroup="condition")
rld_pca.do.primary <- function (rld.d0.primary, intgroup = "condition", ntop = 500, colors=NULL, legendpos="bottomleft", main="PCA Biplot", textcx=1, ...) {
  require(genefilter)
  require(calibrate)
  require(RColorBrewer)
  rv = rowVars(assay(rld.d0.primary))
  select = order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
  pca = prcomp(t(assay(rld.d0.primary)[select, ]))
  fac = factor(apply(as.data.frame(colData(rld.d0.primary)[, intgroup, drop = FALSE]), 1, paste, collapse = " : "))
  if (is.null(colors)) {
    if (nlevels(fac) >= 3) {
      colors = brewer.pal(nlevels(fac), "Paired")
    }   else {
      colors = c("black", "red")
    }
  }
  pc1var <- round(summary(pca)$importance[2,1]*100, digits=1)
  pc2var <- round(summary(pca)$importance[2,2]*100, digits=1)
  pc1lab <- paste0("PC1 (",as.character(pc1var),"%)")
  pc2lab <- paste0("PC1 (",as.character(pc2var),"%)")
  plot(PC2~PC1, data=as.data.frame(pca$x), bg=colors[fac], pch=21, xlab=pc1lab, ylab=pc2lab, main=main, ...)
  with(as.data.frame(pca$x), textxy(PC1, PC2, labs=rownames(as.data.frame(pca$x)), cex=textcx))
  legend(legendpos, legend=levels(fac), col=colors, pch=20)
  #     rldyplot(PC2 ~ PC1, groups = fac, data = as.data.frame(pca$rld),
  #            pch = 16, cerld = 2, aspect = "iso", col = colours, main = draw.key(key = list(rect = list(col = colours),
  #                                                                                         terldt = list(levels(fac)), rep = FALSE)))
}
png("Day0.Primary.qc-pca.png", 1000, 1000, pointsize=20)
rld_pca.do.primary(rld.d0.primary, colors=mycols, intgroup="Primary_Treatment", xlim=c(-75, 35))
dev.off()

## MA plot ------------------------------------------------------------------------------------------------------ #
## Could do with built-in DESeq2 function:
## DESeq2::plotMA(dds, ylim=c(-1,1), cex=1)
# DO.primary_maplot <- function (resd0.primary, thresh=0.05, labelsig=TRUE, textcx=1, ...) {
#   with(res, plot(baseMean, log2FoldChange, pch=20, cex=.5, log="x", ...))
#   with(subset(resd0.primary, padj<thresh), points(baseMean, log2FoldChange, col="red", pch=20, cex=1.5))
#   if (labelsig) {
#     require(calibrate)
#     with(subset(resd0.primary, padj<thresh), textxy(baseMean, log2FoldChange, labs=Gene, cex=textcx, col=2))
#   }
# }
# png("Day0_primary_diffexpr-maplot.png", 1500, 1000, pointsize=20)
# maplot(resdata.d0.primary)
# dev.off()


## Volcano plot with "significant" genes labeled ------------------------------------------------------------------------------------------------------ #
volcanoplot <- function (resd0.primary, lfcthresh=2, sigthresh=0.05, main="Volcano Plot", legendpos="bottomright", labelsig=TRUE, textcx=1, ...) {
  with(resd0.primary, plot(log2FoldChange, -log10(pvalue), pch=20, main=main, ...))
  with(subset(resd0.primary, padj<sigthresh ), points(log2FoldChange, -log10(pvalue), pch=20, col="red", ...))
  with(subset(resd0.primary, abs(log2FoldChange)>lfcthresh), points(log2FoldChange, -log10(pvalue), pch=20, col="orange", ...))
  with(subset(resd0.primary, padj<sigthresh & abs(log2FoldChange)>lfcthresh), points(log2FoldChange, -log10(pvalue), pch=20, col="green", ...))
  if (labelsig) {
    require(calibrate)
    with(subset(resd0.primary, padj<sigthresh & abs(log2FoldChange)>lfcthresh), textxy(log2FoldChange, -log10(pvalue), labs=Gene, cex=textcx, ...))
  }
  legend(legendpos, xjust=1, yjust=1, legend=c(paste("FDR<",sigthresh,sep=""), paste("|LogFC|>",lfcthresh,sep=""), "both"), pch=20, col=c("red","orange","green"))
}
png("D0.PrimaryTreatment_diffexpr-volcanoplot.png", 1200, 1000, pointsize=20)
volcanoplot(resdata.d0.primary, lfcthresh=1, sigthresh=0.05, textcx=.5, xlim=c(-4, 4))
dev.off()



# convert to dataframe and view the downreg and upregulated genes
dataframe_res.d0.ordered <- as.data.frame(res.d0_Ordered)
dataframe_res.d0.ordered$Gene.ID <- rownames(dataframe_res.d0.ordered)
View(dataframe_res.d0.ordered)

res.d0.DOWNREG <- dataframe_res.d0.ordered %>%  # downregulated genes under moderate stress
  dplyr::filter(log2FoldChange < 0) %>% 
  dplyr::filter(padj < 0.1)
View(res.d0.DOWNREG) # view downreg genes

res.d0.UPREG <- dataframe_res.d0.ordered %>% 
  dplyr::filter(log2FoldChange > 0) %>%   # upregualted genes under moderate stress
  dplyr::filter(padj < 0.1)
View(res.d0.UPREG)  # view upreg genes





# Extracting transformed values 
# vst and rlot provide'blind' argument meaning that the transofmrations can proceed with the 
# sample information WITHOUT having the re-estimate the dispersions 
# with blin d== TRUE this will reestimate the dispersions using only an intercepts - this should 
# be used in order to compare samples in a manner wholly unbiases by the infromation about experimental groups
# HOWEVER the blind dispersion estimation (without sample information) is not appropriate if you expect that many or the majority
# of genes (rows) willhave large differernces in counts explainable by the experimental design and one wishes to transform 
# the data for downstream analysis 
# overall, using blind = False is NOT using the information about which samples were in which experimental group 
# in applying the transformation (an unbiased method) 
vsd <- vst(dds.all.treatments, blind=FALSE) #  transformation functions return an object of class DESeqTransform 
rld <- rlog(dds.all.treatments, blind=FALSE) #  rlog = regularized log - transformed the origin count data to the log2 scale by fitting
# a model with a term for each sample and a prior distrubtion on the coefficients which is esimated from the data 
# captures high dispersions for low counts and therefore these genes exhibit higher shrinkage from the rlog
head(assay(rld), 3) # view the heade of the data
# View the effects of transformations on the variance 
# plots the SD of the transformed data, scross scamples, against the mean, using the log transformation, regualrized log, and 
# the variacle stabalizing transforamtion
# this gives log2(n + 1)
ntd <- normTransform(vsd) # normTransform(object, f = log2, pc = 1) normTransform F = function to apply, pc = a pseudocount 
# without specifying the default if 
library("vsn")
meanSdPlot(assay(ntd)) # log transformation 
meanSdPlot(assay(vsd)) # variance stabalizing transforamtion
meanSdPlot(assay(rld)) # regularized log 
# this gives log2(n + 1)
ntd <- normTransform(dds)
library("vsn")
meanSdPlot(assay(ntd))

plotPCA(vsd, intgroup=c("condition", "type"))

# Note: these are all different ways of looking at the same results... 
res.ALL.TRMTS <- results(dds.all.treatments, name="All_Treatment_MSM.d21_vs_A.d0")
res.ALL.TRMTS.contrasts <- results(dds.all.treatments, contrast=c("All_Treatment","MSM.d21","A.d0")) # alternative way of calling the results (i.e. Day_DAY21_vs_Day0)

res.ALL.TRMTS_Ordered <- res.ALL.TRMTS[order(res.ALL.TRMTS$pvalue),] # how to order based on p-value; chosses EHSM vs A - the last pairwise call 
res.ALL.TRMTS_Ordered # ordered based on pvalue
sum(res.ALL.TRMTS_Ordered$pvalue < 0.1, na.rm=TRUE) # 2165  DEGs between Day0 Amb and Day 21 MSM
sum(res.ALL.TRMTS_Ordered$pvalue < 0.05, na.rm=TRUE) # 1191 DEGs between Day0 Amb and Day 21 MSM

summary(res.ALL.TRMTS)

# CONTROL OVER TIME PAIRWISE COMAPARISONS -------------------- #
# OBJECTIVE: Identify the genes that D/N change over life-stage/development during the experiment 
# these genes can be considered as continuously maintained genes

# here I will call only the 6 pairwise comparisons under the control treatment through time 
# A.d0_AA.d7, A.d0_AA.d14, A.d0_AAa.d21, AA.d7_AA.d14, AA.d7_AAA.d21, AA.d14_AAA.d21

# A.d0_AA.d7 # Ambient Controls Day 0 compared to Day 7
res.ALL.TRMTS <- results(dds.all.treatments, name="All_Treatment_AA.d7_vs_A.d0")
res.A.d0_AA.d7 <- results(dds.all.treatments, contrast=c("All_Treatment","AA.d7","A.d0"))
sum(res.A.d0_AA.d7$padj < 0.1, na.rm=TRUE) # 1 marginal DEGs
sum(res.A.d0_AA.d7$padj < 0.05, na.rm=TRUE) # 1 significant DEGs
res.A.d0_AA.d7.DATAFRAME <- as.data.frame(res.A.d0_AA.d7)
res.A.d0_AA.d7.DEGS <- res.A.d0_AA.d7.DATAFRAME  %>%  # downregulated genes under moderate stress
   dplyr::filter(padj < 0.05)
# View(res.A.d0_AA.d7.DEGS) # 1 total DEG
summary(res.A.d0_AA.d7)

# A.d0_AA.D14 # Ambient Controls Day 0 compared to Day 14
res.A.d0_AA.d14 <- results(dds.all.treatments, contrast=c("All_Treatment","AA.d14","A.d0"))
sum(res.A.d0_AA.d14$padj < 0.1, na.rm=TRUE) # 6 marginal DEGs
sum(res.A.d0_AA.d14$padj < 0.05, na.rm=TRUE) # 5 significant DEGs
res.A.d0_AA.d14.DATAFRAME <- as.data.frame(res.A.d0_AA.d14)
res.A.d0_AA.d14.DEGS <- res.A.d0_AA.d14.DATAFRAME  %>%  # downregulated genes under moderate stress
  dplyr::filter(padj < 0.05)
# View(res.A.d0_AA.d14.DEGS) # 5 total DEGs
# one of these, PGEN_.00g045760, is Fatty acid synthase
summary(res.A.d0_AA.d14)

# # A.d0_AAa.d21 # Ambient Controls Day 0 compared to Day 21
res.A.d0_AAA.d21 <- results(dds.all.treatments, contrast=c("All_Treatment","AAA.d21","A.d0"))
sum(res.A.d0_AAA.d21$padj < 0.1, na.rm=TRUE) # 67 marginal DEGs
sum(res.A.d0_AAA.d21$padj < 0.05, na.rm=TRUE) # 20 significant DEGs
res.A.d0_AAA.d21.DATAFRAME <- as.data.frame(res.A.d0_AAA.d21)
res.A.d0_AAA.d21.DEGS <- res.A.d0_AAA.d21.DATAFRAME  %>%  # downregulated genes under moderate stress
  dplyr::filter(padj < 0.05)
# View(res.A.d0_AAA.d21.DEGS) # 20 total DEGs
summary(res.A.d0_AAA.d21)



# what we find are a total of 26 significantly differentially expressed genes (p < 0.05) 
# lets call these in a dataset for fun using the package 'data.frame'
res.A.d0_AA.d7.NS.GENES <- res.A.d0_AA.d7.DATAFRAME %>%  dplyr::filter(padj > 0.05)
res.A.d0_AA.d7.NS.GENES$Gene.ID <- row.names(res.A.d0_AA.d7.NS.GENES) # change row names to a column'Gene.ID'

res.A.d0_AA.d14.NS.GENES <- res.A.d0_AA.d14.DATAFRAME %>%  dplyr::filter(padj > 0.05)
res.A.d0_AA.d14.NS.GENES$Gene.ID <- row.names(res.A.d0_AA.d14.NS.GENES) # change row names to a column'Gene.ID'

res.A.d0_AAA.d21.NS.GENES <- res.A.d0_AAA.d21.DATAFRAME %>%  dplyr::filter(pvalue > 0.05)
res.A.d0_AAA.d21.NS.GENES$Gene.ID <- row.names(res.A.d0_AAA.d21.NS.GENES) # change row names to a column'Gene.ID'


Amb.DEGs <- lapply(list(res.A.d0_AA.d7.NS.GENES[7], res.A.d0_AA.d14.NS.GENES[7], res.A.d0_AAA.d21.NS.GENES[7]), data.table)
Amb.DEGs.IDs.keep <-rbindlist(lapply(Amb.DEGs, '[', j = 'Gene.ID'))[, .N, by=Gene.ID][N == 3L, 'Gene.ID']
Amb.DEGs.keep <- Reduce(funion, Amb.DEGs)[Gene.ID %in% Amb.DEGs.IDs.keep]


X <- as.data.frame(resultsNames(dds.all.treatments))
X.2 <- X[c(2:13),]
X.2 <-  as.data.frame(X.2)
# for loop
df_total <- data.frame()
for (i in 1:nrow(X.2)) {
  RSLTS <- results(dds.all.treatments, name=X.2[i,])
  sum <- sum(RSLTS$pvalue < 0.05, na.rm=TRUE) 
  Total.Genes <- nrow(RSLTS)
  
  in.loop.table <- data.frame(matrix(nrow = 1, ncol = 3))
  colnames(in.loop.table)<-c('test', 'total.genes', 'DEGs.p0.05')
  in.loop.table$test <- X.2[i,]
  in.loop.table$total.genes <- Total.Genes
  in.loop.table$DEGs.p0.05 <- sum
  
  df <- data.frame(in.loop.table)
  df_total <- rbind(df_total,df)#bind to a cumulative list dataframe
}
df_total


# res.AA.d7_AA.d14 <- results(dds.all.treatments, contrast=c("All_Treatment","AA.d7","A.d0"))
# 
# # AA.d7_AA.d14
# res.A.d0_AA.d7 <- results(dds.all.treatments, contrast=c("All_Treatment","AA.d7","A.d0"))
# 
# # AA.d7_AAA.d21
# res.A.d0_AA.d7 <- results(dds.all.treatments, contrast=c("All_Treatment","AA.d7","A.d0"))
# 
# # AA.d14_AAA.d21
# res.A.d0_AA.d7 <- results(dds.all.treatments, contrast=c("All_Treatment","AA.d7","A.d0"))




#======================================================================================= #
# ALL TIME POINTS - PRIMARY TREATMENT ONLY!!!

dds.primary <- DESeq(dds.primary) # wait for this to complete....
resultsNames(dds.primary) # view the names of your results model 'Primary_Treatment_M_vs_A'

# Note: these are all different ways of looking at the same results... 
res <- results(dds.primary) # view DESeq2 results - calls the last pairwise condition
res <- results(dds.primary, name="Primary_Treatment_M_vs_A")
res.contrasts <- results(dds.primary, contrast=c("Primary_Treatment","M","A")) # alternative way of calling the results (i.e. Day_DAY21_vs_Day0)

res_Ordered <- res[order(res$pvalue),] # how to order based on p-value; chosses EHSM vs A - the last pairwise call 
res_Ordered # ordered based on pvalue
sum(res_Ordered$padj < 0.1, na.rm=TRUE) # 1417 marginal DEGs
sum(res_Ordered$padj < 0.05, na.rm=TRUE) # 1071 significant DEGs

summary(res)



#======================================================================================= #
#======================================================================================= #
# DAY 0
#======================================================================================= #
#======================================================================================= #



dds.primary.d0 <- DESeq(dds.primary.d0) # wait for this to complete....
resultsNames(dds.primary.d0) # view the names of your results model 'Primary_Treatment_M_vs_A'

# Get differential expression results
resd0.primary <- results(dds.primary.d0)
table(resd0.primary$padj<0.05) # 13 DEGs
## Order by adjusted p-value
resd0.primary <- resd0.primary[order(resd0.primary$padj), ]
## Merge with normalized count data
resdata.d0.primary  <- merge(as.data.frame(resd0.primary), as.data.frame(counts(dds.primary.d0, normalized=TRUE)), by="row.names", sort=FALSE)
names(resdata.d0.primary)[1] <- "Gene"
head(resdata.d0.primary)
## Write results
write.csv(resdata.d0.primary, file="Day0.PrimaryTreatment_diffexpr-results.csv")
## Examine plot of p-values
hist(resd0.primary$pvalue, breaks=50, col="grey") # view histogram
res.d0 <- results(dds.primary.d0, name="Primary_Treatment_M_vs_A")
table(res.d0$padj<0.05) # should be the same as 'resd0.primary'
res.d0.contrasts <- results(dds.primary.d0, contrast=c("Primary_Treatment","M","A")) # alternative way of calling the results (i.e. Day_DAY21_vs_Day0)
table(res.d0.contrasts$padj<0.05) # should be the same as 'resd0.primary'
res.d0_Ordered <- res.d0[order(res.d0$pvalue),] # how to order based on p-valuel 
res.d0_Ordered # ordered based on pvalue
sum(res.d0_Ordered$padj < 0.05, na.rm=TRUE) # 14 significant DEGs
summary(res.d0)


# Plot dispersions
png("Day0.PrimaryTreatment-dispersions.png", 1000, 1000, pointsize=20)
plotDispEsts(dds.primary.d0, main="Dispersion plot_M_vs._A_Day0")
dev.off()
# Regularized log transformation for clustering/heatmaps, etc
rld.d0.primary <- rlogTransformation(dds.primary.d0) # rlog transform (regularized log)
head(assay(rld.d0.primary)) # view first few rows
hist(assay(rld.d0.primary)) # view histogram 
# Colors for plots below; Use RColorBrewer, better and assign mycols variable
mycols <- brewer.pal(8, "Dark2")[1:length(unique(rld.d0.primary$Primary_Treatment))]
# Sample distance heatmap
sampleDists <- as.matrix(dist(t(assay(rld.d0.primary))))
png("Day0.PrimaryTreatment-heatmap-samples.png", w=1000, h=1000, pointsize=20)
heatmap.2(as.matrix(sampleDists), key=F, trace="none",
          col=colorpanel(100, "black", "white"),
          ColSideColors=mycols[rld.d0.primary$Primary_Treatment], RowSideColors=mycols[rld.d0.primary$Primary_Treatment],
          margin=c(10, 10), main="Sample Distance Matrix")
dev.off()

# ## Examine independent filtering
# attr(resd0.primary, "filterThreshold")
# plot(attr(resd0.primary,"filterNumRej"), type="b", xlab="quantiles of baseMean", ylab="number of rejections")

# Principal components analysis ------------------------------------------------------------------------------------------------------ #
## Could do with built-in DESeq2 function:
## DESeq2::plotPCA(rld, intgroup="condition")
rld_pca.do.primary <- function (rld.d0.primary, intgroup = "condition", ntop = 500, colors=NULL, legendpos="bottomleft", main="PCA Biplot", textcx=1, ...) {
  require(genefilter)
  require(calibrate)
  require(RColorBrewer)
  rv = rowVars(assay(rld.d0.primary))
  select = order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
  pca = prcomp(t(assay(rld.d0.primary)[select, ]))
  fac = factor(apply(as.data.frame(colData(rld.d0.primary)[, intgroup, drop = FALSE]), 1, paste, collapse = " : "))
  if (is.null(colors)) {
    if (nlevels(fac) >= 3) {
      colors = brewer.pal(nlevels(fac), "Paired")
    }   else {
      colors = c("black", "red")
    }
  }
  pc1var <- round(summary(pca)$importance[2,1]*100, digits=1)
  pc2var <- round(summary(pca)$importance[2,2]*100, digits=1)
  pc1lab <- paste0("PC1 (",as.character(pc1var),"%)")
  pc2lab <- paste0("PC1 (",as.character(pc2var),"%)")
  plot(PC2~PC1, data=as.data.frame(pca$x), bg=colors[fac], pch=21, xlab=pc1lab, ylab=pc2lab, main=main, ...)
  with(as.data.frame(pca$x), textxy(PC1, PC2, labs=rownames(as.data.frame(pca$x)), cex=textcx))
  legend(legendpos, legend=levels(fac), col=colors, pch=20)
  #     rldyplot(PC2 ~ PC1, groups = fac, data = as.data.frame(pca$rld),
  #            pch = 16, cerld = 2, aspect = "iso", col = colours, main = draw.key(key = list(rect = list(col = colours),
  #                                                                                         terldt = list(levels(fac)), rep = FALSE)))
}
png("Day0.Primary.qc-pca.png", 1000, 1000, pointsize=20)
rld_pca.do.primary(rld.d0.primary, colors=mycols, intgroup="Primary_Treatment", xlim=c(-75, 35))
dev.off()

## MA plot ------------------------------------------------------------------------------------------------------ #
## Could do with built-in DESeq2 function:
## DESeq2::plotMA(dds, ylim=c(-1,1), cex=1)
# DO.primary_maplot <- function (resd0.primary, thresh=0.05, labelsig=TRUE, textcx=1, ...) {
#   with(res, plot(baseMean, log2FoldChange, pch=20, cex=.5, log="x", ...))
#   with(subset(resd0.primary, padj<thresh), points(baseMean, log2FoldChange, col="red", pch=20, cex=1.5))
#   if (labelsig) {
#     require(calibrate)
#     with(subset(resd0.primary, padj<thresh), textxy(baseMean, log2FoldChange, labs=Gene, cex=textcx, col=2))
#   }
# }
# png("Day0_primary_diffexpr-maplot.png", 1500, 1000, pointsize=20)
# maplot(resdata.d0.primary)
# dev.off()


## Volcano plot with "significant" genes labeled ------------------------------------------------------------------------------------------------------ #
volcanoplot <- function (resd0.primary, lfcthresh=2, sigthresh=0.05, main="Volcano Plot", legendpos="bottomright", labelsig=TRUE, textcx=1, ...) {
  with(resd0.primary, plot(log2FoldChange, -log10(pvalue), pch=20, main=main, ...))
  with(subset(resd0.primary, padj<sigthresh ), points(log2FoldChange, -log10(pvalue), pch=20, col="red", ...))
  with(subset(resd0.primary, abs(log2FoldChange)>lfcthresh), points(log2FoldChange, -log10(pvalue), pch=20, col="orange", ...))
  with(subset(resd0.primary, padj<sigthresh & abs(log2FoldChange)>lfcthresh), points(log2FoldChange, -log10(pvalue), pch=20, col="green", ...))
  if (labelsig) {
    require(calibrate)
    with(subset(resd0.primary, padj<sigthresh & abs(log2FoldChange)>lfcthresh), textxy(log2FoldChange, -log10(pvalue), labs=Gene, cex=textcx, ...))
  }
  legend(legendpos, xjust=1, yjust=1, legend=c(paste("FDR<",sigthresh,sep=""), paste("|LogFC|>",lfcthresh,sep=""), "both"), pch=20, col=c("red","orange","green"))
}
png("Day0.PrimaryTreatment_diffexpr-volcanoplot.png", 1200, 1000, pointsize=20)
volcanoplot(resdata.d0.primary, lfcthresh=1, sigthresh=0.05, textcx=.5, xlim=c(-4, 4))
dev.off()



# convert to dataframe and view the downreg and upregulated genes
dataframe_res.d0.ordered <- as.data.frame(res.d0_Ordered)
dataframe_res.d0.ordered$Gene.ID <- rownames(dataframe_res.d0.ordered)
View(dataframe_res.d0.ordered)

res.d0.DOWNREG <- dataframe_res.d0.ordered %>%  # downregulated genes under moderate stress
  dplyr::filter(log2FoldChange < 0) %>% 
  dplyr::filter(padj < 0.1)
View(res.d0.DOWNREG) # view downreg genes

res.d0.UPREG <- dataframe_res.d0.ordered %>% 
  dplyr::filter(log2FoldChange > 0) %>%   # upregualted genes under moderate stress
  dplyr::filter(padj < 0.1)
View(res.d0.UPREG)  # view upreg genes




#======================================================================================= #
#======================================================================================= #
#======================================================================================= #
# DAY 21   - Primary Treatment only ==================================================== #
#======================================================================================= #
#======================================================================================= #
#======================================================================================= #


dds.primary.d21 <- DESeq(dds.primary.d21) # wait for this to complete....
resultsNames(dds.primary.d21) # view the names of your results model 'Primary_Treatment_M_vs_A'

# Results ------------------------------------------------------------------------------ #
res.d21 <- results(dds.primary.d21) # view DESeq2 results 
table(res.d21$padj<0.05) 
res.d21 <- results(dds.primary.d21, name="Primary_Treatment_M_vs_A")
table(res.d21$padj<0.05) # should be the same as 'res.d21'
res.D21.contrasts <- results(dds.primary.d21, contrast=c("Primary_Treatment","M","A")) # alternative way of calling the results (i.e. Day_DAY21_vs_Day0)
table(res.D21.contrasts$padj<0.05) # should be the same as 'res.d21'
res.d21_Ordered <- res.d21[order(res.d21$pvalue),] # how to order based on p-value
res.d21_Ordered # ordered based on pvalue
sum(res.d21_Ordered$padj < 0.05, na.rm=TRUE) # 196 significant DEGs
summary(res.d21) # shows pval <0.1 by default 
## Merge with normalized count data 
resdata.d21.primary  <- merge(as.data.frame(res.d21), as.data.frame(counts(dds.primary.d21, normalized=TRUE)), by="row.names", sort=FALSE)
names(resdata.d21.primary)[1] <- "Gene"
head(resdata.d21.primary)
## Write results
write.csv(resdata.d21.primary, file="D21.PrimaryTreatment_diffexpr-results.csv")

# Plot dispersions -------------------------------------------------------------------- #
png("Day21.PrimaryTreatment-dispersions.png", 1000, 1000, pointsize=20)
plotDispEsts(dds.primary.d21, main="Day21.PrimaryTreatment_Dispersion plot")
dev.off()

# Regularized log transformation for clustering/heatmaps, etc ------------------------- #
rld.d21.primary <- rlogTransformation(dds.primary.d21) # rlog transform (regularized log)
head(assay(rld.d21.primary)) # view first few rows
hist(assay(rld.d21.primary)) # view histogram 

# Colors for plots below; Use RColorBrewer, better and assign mycols variable --------- #
mycols <- brewer.pal(8, "Dark2")[1:length(unique(rld.d21.primary$Primary_Treatment))]

# Sample distance heatmap ------------------------------------------------------------- #
sampleDists <- as.matrix(dist(t(assay(rld.d21.primary))))
png("Day21.PrimaryTreatment-heatmap-samples.png", w=1000, h=1000, pointsize=20)
heatmap.2(as.matrix(sampleDists), key=F, trace="none",
          col=colorpanel(100, "black", "white"),
          ColSideColors=mycols[rld.d21.primary$Primary_Treatment], RowSideColors=mycols[rld.d21.primary$Primary_Treatment],
          margin=c(10, 10), main="Sample Distance Matrix")
dev.off()

# Principal components analysis ------------------------------------------------------ #
## Could do with built-in DESeq2 function:
png("Day21.Primary.qc-pca_2.png", 1000, 1000, pointsize=20)
DESeq2::plotPCA(rld.d21.primary, intgroup="Primary_Treatment")
dev.off()

rld_pca.d21.primary <- function (rld.d21.primary, intgroup = "condition", ntop = 500, colors=NULL, legendpos="bottomleft", main="PCA Biplot", textcx=1, ...) {
  require(genefilter)
  require(calibrate)
  require(RColorBrewer)
  rv = rowVars(assay(rld.d21.primary))
  select = order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
  pca = prcomp(t(assay(rld.d21.primary)[select, ]))
  fac = factor(apply(as.data.frame(colData(rld.d21.primary)[, intgroup, drop = FALSE]), 1, paste, collapse = " : "))
  if (is.null(colors)) {
    if (nlevels(fac) >= 3) {
      colors = brewer.pal(nlevels(fac), "Paired")
    }   else {
      colors = c("black", "red")
    }
  }
  pc1var <- round(summary(pca)$importance[2,1]*100, digits=1)
  pc2var <- round(summary(pca)$importance[2,2]*100, digits=1)
  pc1lab <- paste0("PC1 (",as.character(pc1var),"%)")
  pc2lab <- paste0("PC1 (",as.character(pc2var),"%)")
  plot(PC2~PC1, data=as.data.frame(pca$x), bg=colors[fac], pch=21, xlab=pc1lab, ylab=pc2lab, main=main, ...)
  with(as.data.frame(pca$x), textxy(PC1, PC2, labs=rownames(as.data.frame(pca$x)), cex=textcx))
  legend(legendpos, legend=levels(fac), col=colors, pch=20)
  #     rldyplot(PC2 ~ PC1, groups = fac, data = as.data.frame(pca$rld),
  #            pch = 16, cerld = 2, aspect = "iso", col = colours, main = draw.key(key = list(rect = list(col = colours),
  #                                                                                         terldt = list(levels(fac)), rep = FALSE)))
}
png("Day21.Primary.qc-pca.png", 1000, 1000, pointsize=20)
rld_pca.d21.primary(rld.d21.primary, colors=mycols, intgroup="Primary_Treatment", xlim=c(-100, 50))
dev.off()

## Volcano plot with "significant" genes labeled --------------------------------- #
volcanoplot <- function (res.d21, lfcthresh=2, sigthresh=0.05, main="Volcano Plot", legendpos="bottomright", labelsig=TRUE, textcx=1, ...) {
  with(res.d21, plot(log2FoldChange, -log10(pvalue), pch=20, main=main, ...))
  with(subset(res.d21, padj<sigthresh ), points(log2FoldChange, -log10(pvalue), pch=20, col="red", ...))
  with(subset(res.d21, abs(log2FoldChange)>lfcthresh), points(log2FoldChange, -log10(pvalue), pch=20, col="orange", ...))
  with(subset(res.d21, padj<sigthresh & abs(log2FoldChange)>lfcthresh), points(log2FoldChange, -log10(pvalue), pch=20, col="green", ...))
  if (labelsig) {
    require(calibrate)
    with(subset(res.d21, padj<sigthresh & abs(log2FoldChange)>lfcthresh), textxy(log2FoldChange, -log10(pvalue), labs=Gene, cex=textcx, ...))
  }
  legend(legendpos, xjust=1, yjust=1, legend=c(paste("FDR<",sigthresh,sep=""), paste("|LogFC|>",lfcthresh,sep=""), "both"), pch=20, col=c("red","orange","green"))
}
png("Day21.PrimaryTreatment_diffexpr-volcanoplot.png", 1200, 1000, pointsize=20)
volcanoplot(resdata.d21.primary, lfcthresh=1, sigthresh=0.05, textcx=.5, xlim=c(-4, 4))
dev.off()


# convert to dataframe and view the downreg and upregulated genes
dataframe_res.d21.ordered <- as.data.frame(res.d21_Ordered)
dataframe_res.d21.ordered$Gene.ID <- rownames(dataframe_res.d21.ordered)
View(dataframe_res.d21.ordered)

res.d21.DOWNREG <- dataframe_res.d21.ordered %>%  # downregulated genes under moderate stress
  dplyr::filter(log2FoldChange < 0) %>% 
  dplyr::filter(padj < 0.1)
View(res.d21.DOWNREG) # view downreg genes

res.d21.UPREG <- dataframe_res.d21.ordered %>% 
  dplyr::filter(log2FoldChange > 0) %>%   # upregualted genes under moderate stress
  dplyr::filter(padj < 0.1)
View(res.d21.UPREG)  # view upreg genes

#ddsMF <- dds # make a coppy of dds for multi factor designs



#======================================================================================= #
#======================================================================================= #
#======================================================================================= #
# DAY 21   - Primary Treatment only ==================================================== #
#======================================================================================= #
#======================================================================================= #
#======================================================================================= #


dds.All.d21 <- DESeq(dds.All.d21) # wait for this to complete....
resultsNames(dds.All.d21) # view the names of your results models

# Results ------------------------------------------------------------------------------ #
res.d21.All <- results(dds.All.d21) # view DESeq2 results 
table(res.d21.All$padj<0.05) 
res.d21.All <- results(dds.All.d21, name="Primary_Treatment_M_vs_A")
table(res.d21.All$padj<0.05) # should be the same as 'res.d21'
res.D21.contrasts <- results(dds.primary.d21, contrast=c("Primary_Treatment","M","A")) # alternative way of calling the results (i.e. Day_DAY21_vs_Day0)
table(res.D21.contrasts$padj<0.05) # should be the same as 'res.d21'
res.d21_Ordered <- res.d21[order(res.d21$pvalue),] # how to order based on p-value
res.d21_Ordered # ordered based on pvalue
sum(res.d21_Ordered$padj < 0.05, na.rm=TRUE) # 196 significant DEGs
summary(res.d21) # shows pval <0.1 by default 
## Merge with normalized count data 
resdata.d21.primary  <- merge(as.data.frame(res.d21), as.data.frame(counts(dds.primary.d21, normalized=TRUE)), by="row.names", sort=FALSE)
names(resdata.d21.primary)[1] <- "Gene"
head(resdata.d21.primary)
## Write results
write.csv(resdata.d21.primary, file="D21.PrimaryTreatment_diffexpr-results.csv")

# Plot dispersions -------------------------------------------------------------------- #
png("Day21.PrimaryTreatment-dispersions.png", 1000, 1000, pointsize=20)
plotDispEsts(dds.primary.d21, main="Day21.PrimaryTreatment_Dispersion plot")
dev.off()

# Regularized log transformation for clustering/heatmaps, etc ------------------------- #
rld.d21.primary <- rlogTransformation(dds.primary.d21) # rlog transform (regularized log)
head(assay(rld.d21.primary)) # view first few rows
hist(assay(rld.d21.primary)) # view histogram 

# Colors for plots below; Use RColorBrewer, better and assign mycols variable --------- #
mycols <- brewer.pal(8, "Dark2")[1:length(unique(rld.d21.primary$Primary_Treatment))]

# Sample distance heatmap ------------------------------------------------------------- #
sampleDists <- as.matrix(dist(t(assay(rld.d21.primary))))
png("Day21.PrimaryTreatment-heatmap-samples.png", w=1000, h=1000, pointsize=20)
heatmap.2(as.matrix(sampleDists), key=F, trace="none",
          col=colorpanel(100, "black", "white"),
          ColSideColors=mycols[rld.d21.primary$Primary_Treatment], RowSideColors=mycols[rld.d21.primary$Primary_Treatment],
          margin=c(10, 10), main="Sample Distance Matrix")
dev.off()

# Principal components analysis ------------------------------------------------------ #
## Could do with built-in DESeq2 function:
png("Day21.Primary.qc-pca_2.png", 1000, 1000, pointsize=20)
DESeq2::plotPCA(rld.d21.primary, intgroup="Primary_Treatment")
dev.off()

rld_pca.d21.primary <- function (rld.d21.primary, intgroup = "condition", ntop = 500, colors=NULL, legendpos="bottomleft", main="PCA Biplot", textcx=1, ...) {
  require(genefilter)
  require(calibrate)
  require(RColorBrewer)
  rv = rowVars(assay(rld.d21.primary))
  select = order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
  pca = prcomp(t(assay(rld.d21.primary)[select, ]))
  fac = factor(apply(as.data.frame(colData(rld.d21.primary)[, intgroup, drop = FALSE]), 1, paste, collapse = " : "))
  if (is.null(colors)) {
    if (nlevels(fac) >= 3) {
      colors = brewer.pal(nlevels(fac), "Paired")
    }   else {
      colors = c("black", "red")
    }
  }
  pc1var <- round(summary(pca)$importance[2,1]*100, digits=1)
  pc2var <- round(summary(pca)$importance[2,2]*100, digits=1)
  pc1lab <- paste0("PC1 (",as.character(pc1var),"%)")
  pc2lab <- paste0("PC1 (",as.character(pc2var),"%)")
  plot(PC2~PC1, data=as.data.frame(pca$x), bg=colors[fac], pch=21, xlab=pc1lab, ylab=pc2lab, main=main, ...)
  with(as.data.frame(pca$x), textxy(PC1, PC2, labs=rownames(as.data.frame(pca$x)), cex=textcx))
  legend(legendpos, legend=levels(fac), col=colors, pch=20)
  #     rldyplot(PC2 ~ PC1, groups = fac, data = as.data.frame(pca$rld),
  #            pch = 16, cerld = 2, aspect = "iso", col = colours, main = draw.key(key = list(rect = list(col = colours),
  #                                                                                         terldt = list(levels(fac)), rep = FALSE)))
}
png("Day21.Primary.qc-pca.png", 1000, 1000, pointsize=20)
rld_pca.d21.primary(rld.d21.primary, colors=mycols, intgroup="Primary_Treatment", xlim=c(-100, 50))
dev.off()

## Volcano plot with "significant" genes labeled --------------------------------- #
volcanoplot <- function (res.d21, lfcthresh=2, sigthresh=0.05, main="Volcano Plot", legendpos="bottomright", labelsig=TRUE, textcx=1, ...) {
  with(res.d21, plot(log2FoldChange, -log10(pvalue), pch=20, main=main, ...))
  with(subset(res.d21, padj<sigthresh ), points(log2FoldChange, -log10(pvalue), pch=20, col="red", ...))
  with(subset(res.d21, abs(log2FoldChange)>lfcthresh), points(log2FoldChange, -log10(pvalue), pch=20, col="orange", ...))
  with(subset(res.d21, padj<sigthresh & abs(log2FoldChange)>lfcthresh), points(log2FoldChange, -log10(pvalue), pch=20, col="green", ...))
  if (labelsig) {
    require(calibrate)
    with(subset(res.d21, padj<sigthresh & abs(log2FoldChange)>lfcthresh), textxy(log2FoldChange, -log10(pvalue), labs=Gene, cex=textcx, ...))
  }
  legend(legendpos, xjust=1, yjust=1, legend=c(paste("FDR<",sigthresh,sep=""), paste("|LogFC|>",lfcthresh,sep=""), "both"), pch=20, col=c("red","orange","green"))
}
png("Day21.PrimaryTreatment_diffexpr-volcanoplot.png", 1200, 1000, pointsize=20)
volcanoplot(resdata.d21.primary, lfcthresh=1, sigthresh=0.05, textcx=.5, xlim=c(-4, 4))
dev.off()



#=====================================================================================
# EXPLORING AND EXPORTING RESULTS
#=====================================================================================

plotMA(res.d21, ylim=c(-2,2))

# It is more useful visualize the MA-plot for the shrunken log2 fold changes, 
# which remove the noise associated with log2 fold changes from low count genes 
# without requiring arbitrary filtering thresholds.
res.d21_LFC <- lfcShrink(dds.primary.d21, coef="Primary_Treatment_M_vs_A", type="apeglm")
summary(res.d21_LFC)
plotMA(res.d21_LFC, ylim=c(-2,2))
# after calling plotMA you can use function identify to interactively detect the row number of individual genes by clicking on the plot
idx <- identify(res.d21_LFC$baseMean, res.d21_LFC$log2FoldChange)
rownames(res.d21_LFC)[idx]

# plot counts
plotCounts(dds.primary.d21, gene=which.min(res.d21_LFC$padj), intgroup="Primary_Treatment") # this plots the gene which had the smallest p-value
plotCounts(dds.primary.d21, gene='PGEN_.00g108770', intgroup="Primary_Treatment") # putative alternative oxidase!

# sample plot in ggplot 
d <- plotCounts(dds, gene=which.min(res.d21_LFC$padj), intgroup="Treatment", 
                returnData=TRUE)
ggplot(d, aes(x=Treatment, y=count)) + 
  theme_classic() +
  geom_point(position=position_jitter(w=0.1,h=0)) + 
  scale_y_log10(breaks=c(0, 25,100,400,2000))
# heatmap of count matrix
ntd <- normTransform(dds)
select <- order(rowMeans(counts(dds,normalized=TRUE)),
                decreasing=TRUE)[1:20]
df <- as.data.frame(colData(dds)[,c("Treatment","Day")])
pheatmap(assay(ntd)[select,], cluster_rows=T, show_rownames=T,
         cluster_cols=F, annotation_col=df)
