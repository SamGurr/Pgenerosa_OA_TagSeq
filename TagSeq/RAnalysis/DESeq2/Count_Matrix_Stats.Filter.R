---
  # title: "Count_Matrix_Stats.Filter.R"
  # author: "Samuel Gurr"
  # date: "January 24, 2021"
---
  
# LOAD PACKAGES
library(dplyr)

# SET WORKING DIRECTORY AND LOAD DATA
setwd("C:/Users/samjg/Documents/My_Projects/Pgenerosa_TagSeq_Metabolomics/TagSeq/")

# LOAD DATA
# raw_counts <- read.csv(file="HPC_Bioinf/outputs/transcript_count_matrix.csv", sep=',', header=TRUE) # read the output count matrix - NOTE: NOT trimmed at 30 phred threshold!
raw_counts <- read.csv(file="HPC_Bioinf/outputs/transcript_count_matrix_trimpolyx_10.csv", sep=',', header=TRUE) # read the output count matrix - NOTE: TRIMMED at 30 phred threshold!
UT_seq_map <- read.csv(file="20201020_Gurr_TagSeq_UTAustin.csv", sep=',', header=TRUE)
smpl_ref <- read.csv(file="Sample_reference.csv", sep=',', header=TRUE)
treatment_ref <- read.csv(file="Extraction_checklist.csv", sep=',', header=TRUE)

# =====================================================================================================================
#                                                                                                                       
#                             FORMAT EXPERIMENT DESIGN DATAFRAME
# ======================================================================================================================
# path for outputting all .csv filtered count files
path = 'C:/Users/samjg/Documents/My_Projects/Pgenerosa_TagSeq_Metabolomics/TagSeq/RAnalysis/DESeq2/counts_filtered/' # run this for all count matrix outputs!!!
# ========================================================== 
# MASTER REFERENCE DATA.FRAME
# format and merge to buld master reference dataframe
smpl_ref$Seq_Pos <- paste(smpl_ref$ï..TagSeq_Plate, smpl_ref$TagSeq_Well, sep="_")
smpl_ref <- smpl_ref[,-c(1:2)]
UT_seq_map$Seq_Pos  <- paste(UT_seq_map$ï..Plate, UT_seq_map$Well, sep="_")
UT_seq_map <- UT_seq_map[-c(1:2)]
Seq.Ref <- merge(smpl_ref, UT_seq_map, by = "Seq_Pos")
Mstr.Ref <- merge(Seq.Ref, treatment_ref, by = "Geoduck_ID")

#================= #
# ALL TIMEPOINTS:
# call all experiment design treatments as 'exp.data'
exp.data <- Mstr.Ref[,c("Sample.Name","All_Treatment", "Primary_Treatment", "Second_Treament", "Third_Treatment", "Time")]
write.csv(exp.data,paste(path,"all.exp.data.csv"))
#================= #
# DAY 0 
exp.data.d0 <- exp.data %>% dplyr::filter(Time %in% 'Day0') # all data on day 0
nrow(exp.data.d0) # 8 total samples on Day 0
write.csv(exp.data.d0,paste(path,"day0.exp.data.csv"))
#================= #
# DAY 7 
exp.data.d7 <- exp.data %>% dplyr::filter(Time %in% 'Day7') # all data on day 7
nrow(exp.data.d7) # 36 total samples on Day 7
write.csv(exp.data.d7,paste(path,"day7.exp.data.csv"))
#================= #
# DAY 14 
UT_seq_map %>% dplyr::filter(Sample.Name == "SG92") # there was no sample in SG92 for TagSeq; 35 total is correct!

exp.data.d14 <- exp.data %>% dplyr::filter(Time %in% 'DAY14') # all data on day 14
exp.data.d14 <- exp.data.d14 %>% dplyr::filter(!(Sample.Name %in% 'SG92')) # RUN THIS! now 35 rows (samples) ommmits SG92 (NOT sent to for TagSeq!)
nrow(exp.data.d14) #  35 total samples on Day 14
write.csv(exp.data.d14,paste(path," day14.exp.data.csv"))
#================= #
# DAY 21 
exp.data.d21 <- exp.data %>% dplyr::filter(Time %in% 'DAY21') # all data on day 21
nrow(exp.data.d21) # 62 total sampels on day 21
write.csv(exp.data.d21,paste(path," day21.exp.data.csv"))
#======================================================================================================================
#                                                                                                                       
#                             FORMAT COUNT MATRIX
#                                                                                                                       
#====================================================================================================================== 
# path for outputting all .csv filtered count files
path = 'C:/Users/samjg/Documents/My_Projects/Pgenerosa_TagSeq_Metabolomics/TagSeq/RAnalysis/DESeq2/counts_filtered/' # run this for all count matrix outputs!!!

# ========================================================== 
#
# INITIAL FORMATING AND READ COUNTS
# ========================================================== 
# About: need to call the raw_counts.matrix and call only samples that match the IDs at day 0 and 21 
# Why? DESeq2 can only apply to factors with >=2 levels, so we cannot address 'Day' at just "Day21" for example
# alternatively, we can call subset matrices for the timepoints of interest

# NOTE: this raw_counts matrix was NOT merged by lanes, run'rowsum' on unique delimiters to merge 'raw_counts.merged'
#  Format count matrix (raw_counts has 2x columns per sample on different lanes - sum counts from unique columns + build matrix)
ncol(raw_counts) # 282 samples (not counting gene ID column) - should be 141 samples, need to sum columns by unique ID
raw_counts.merged <- data.frame(raw_counts[,-1], row.names=raw_counts[,1]) # call new dataframe with first column now as row names, now all row values are numeric
names(raw_counts.merged) <- sapply(strsplit(names(raw_counts.merged), "_"), '[', 1) # split the column names by "_" delimiter  and call the first field SG##
raw_counts.merged <- t(rowsum(t(raw_counts.merged), group = colnames(raw_counts.merged), na.rm = TRUE)) # merge all unique columns and sum counts 
ncol(raw_counts.merged) # now 141 samples

raw_counts.matrix <-as.matrix(raw_counts.merged, row.names="transcript_id") # call dataframe as matrix
ncol(raw_counts.matrix) # 141 samples
nrow(raw_counts.matrix) # 34947 total genes

raw_counts.merged.as.table <- data.frame(transcript_id = row.names(raw_counts.merged), raw_counts.merged) # add back the rownames 'transcript_ID'
rownames(raw_counts.merged.as.table) <- NULL # ommit the rownames
ncol(raw_counts.merged.as.table) # 142 counting the transcript.ID that we want to keep! 


# READ COUNTS 
dim(raw_counts.matrix) # 34947 total genes 141 samples
sum(transcript_sums) # 114250931 total read counts 

gene_sums <- data.frame(rowSums(raw_counts.matrix))  # all gene.IDs and the sum of unique reads

transcript_sums <- colSums(raw_counts.matrix) # sum of reads for each sample
mean(transcript_sums) # 810290.3 == average raw read counts for each sample

gene_sums_gtr0 <- rowSums(raw_counts.matrix) > 0 # all gene.IDs with at least one unique read
sum(gene_sums_gtr0 == TRUE) # 29335 total genes with unique transcript reads 
( sum(gene_sums_gtr0 == TRUE) / (dim(raw_counts.matrix)[1]) ) *100 # 83.9414 % of genes have a unique mapped read

# ========================================================== 
#
# ALL TIMEPOINTS (3 CPM in 50% samples using edgeR)
# ========================================================== 

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
dim(cts.matrix.all.filtered) # 13870 genes & 141 samples
#head(cts.matrix.all.filtered) # FINAL ALL DATASET
# We will look at a cople samples to check our defined threshold.. 
# plot(CPM.all[,1], cts.matrix.all[,1], xlab="CPM", ylab="Raw Count", ylim=c(0,50), xlim=c(0,50), main=colnames(CPM.all)[1]) # check if our threshold of 3 CPM in 50% samples corresponds to a count of about 10-15
# abline(v=3) 
# abline(h=3)
# plot(CPM.all[,70], cts.matrix.all[,70], xlab="CPM", ylab="Raw Count", ylim=c(0,50), xlim=c(0,50), main=colnames(CPM.all)[1]) # check if our threshold of 3 CPM in 50% samples corresponds to a count of about 10-15
# abline(v=2.37) 
# abline(h=3)

# write csv
write.csv(cts.matrix.all.filtered,paste(path,"all.counts.filtered_matrix.csv")) # 'path' called in previous # write .csv section


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
summary(keep.d0) # FALSE 19798 & TRUE 15149 -- more than half of the genes did not pass
cts.matrix.d0.filtered <- cts.matrix.d0[keep.d0,] # Subset the rows of countdata to keep the more highly expressed genes
dim(cts.matrix.d0.filtered) # 15149 genes & 8 samples
# head(cts.matrix.d0.filtered) # FINAL DAY 0 DATASET
# We will look at a cople samples to check our defined threshold.. 
# plot(CPM.d0[,1], cts.matrix.d0[,1], xlab="CPM", ylab="Raw Count", ylim=c(0,50), xlim=c(0,50), main=colnames(CPM.d0)[1]) # check if our threshold of 3 CPM in 50% samples corresponds to a count of about 10-15
# abline(v=3) 
# abline(h=3)
# plot(CPM.d0[,5], cts.matrix.d0[,5], xlab="CPM", ylab="Raw Count", ylim=c(0,50), xlim=c(0,50), main=colnames(CPM.d0)[1]) # check if our threshold of 3 CPM in 50% samples corresponds to a count of about 10-15
# abline(v=3) 
# abline(h=3)

# write csv
write.csv(cts.matrix.d0.filtered,paste(path,"day0.counts.filtered_matrix.csv")) # 'path' called in previous # write .csv section

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
table(rowSums(thresh.d7)) # 6880 genes with TRUE in all 36 samples 
keep.d7 <- rowSums(thresh.d7) >= (ncol(thresh.d7)/2) # we would like to keep genes that have at least 50% TRUES in each row of thresh
summary(keep.d7) # FALSE 20825 & TRUE 14122 -- more than half of the genes did not pass
cts.matrix.d7.filtered <- cts.matrix.d7[keep.d7,] # Subset the rows of countdata to keep the more highly expressed genes
dim(cts.matrix.d7.filtered) # 14122 genes & 36 samples
# head(cts.matrix.d7.filtered) # FINAL DAY 7 DATASET
# We will look at a cople samples to check our defined threshold.. 
# plot(CPM.d7[,1], cts.matrix.d7[,1], xlab="CPM", ylab="Raw Count", ylim=c(0,50), xlim=c(0,50), main=colnames(CPM.d7)[1]) # check if our threshold of 3 CPM in 50% samples corresponds to a count of about 10-15
# abline(v=3) 
# abline(h=3)
# plot(CPM.d7[,5], cts.matrix.d7[,5], xlab="CPM", ylab="Raw Count", ylim=c(0,50), xlim=c(0,50), main=colnames(CPM.d7)[1]) # check if our threshold of 3 CPM in 50% samples corresponds to a count of about 10-15
# abline(v=3) 
# abline(h=3)

# write csv
write.csv(cts.matrix.d7.filtered,paste(path,"day7.counts.filtered_matrix.csv")) # 'path' called in previous # write .csv section

# ========================================================== 
#
# DAY 14 (3 CPM in 50% samples using edgeR)
# ========================================================== 
cts.merged.d14 <- cts.merged.as.table[,c(1,na.omit(match(exp.data.d14$Sample.Name, colnames(cts.merged.as.table))))]
cts.merged.d14 <- data.frame(cts.merged.d14[,-1], row.names=cts.merged.d14[,1])
cts.matrix.d14  <-as.matrix(cts.merged.d14, row.names="transcript_id")
ncol(cts.matrix.d14) # 35 samples from just Day 14
colnames(cts.matrix.d14) == exp.data.d14$Sample.Name # chec if all TRUE; NOTE: SG92 was ommitted earlier to make sure this reads TRUE
UT_seq_map %>% dplyr::filter(Sample.Name == "SG92") # there was no sample in SG92 for TagSeq; 35 total is correct!
# pre-filtering; genes ommitted if < 3 counts per million reads in 50% of samples
colSums(cts.matrix.d14) # view the colSums of our Day14 samples - notice counts are near 1 million
CPM.d14 <- cpm(cts.matrix.d14) # Obtain CPMs (counts oer million) using egdeR
head(CPM.d14) # Have a look at the output
thresh.d14 <- CPM.d14 > 3 # Which values in myCPM are greater than 3?
head(thresh.d14) # This produces a logical matrix with TRUEs and FALSES
rowSums(head(thresh.d14)) # Summary of how many TRUEs there are in each row
table(rowSums(thresh.d14)) # 6473 genes with TRUE in all 35 samples 
keep.d14 <- rowSums(thresh.d14) >= (ncol(thresh.d14)/2) # we would like to keep genes that have at least 50% TRUES in each row of thresh
summary(keep.d14) # FALSE 20710 & TRUE 14237 -- more than half of the genes did not pass
cts.matrix.d14.filtered <- cts.matrix.d14[keep.d14,] # Subset the rows of countdata to keep the more highly expressed genes
dim(cts.matrix.d14.filtered) # 14437 genes & 35 samples
# head(cts.matrix.d14.filtered) # FINAL DAY 14 DATASET
# We will look at a cople samples to check our defined threshold.. 
# plot(CPM.d14[,1], cts.matrix.d14[,1], xlab="CPM", ylab="Raw Count", ylim=c(0,50), xlim=c(0,50), main=colnames(CPM.d14)[1]) # check if our threshold of 3 CPM in 50% samples corresponds to a count of about 10-15
# abline(v=3) 
# abline(h=3)
# plot(CPM.d14[,5], cts.matrix.d14[,5], xlab="CPM", ylab="Raw Count", ylim=c(0,50), xlim=c(0,50), main=colnames(CPM.d14)[1]) # check if our threshold of 3 CPM in 50% samples corresponds to a count of about 10-15
# abline(v=3) 
# abline(h=3)

# write csv
write.csv(cts.matrix.d14.filtered,paste(path,"day14.counts.filtered_matrix.csv")) # 'path' called in previous # write .csv section

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
table(rowSums(thresh.d21)) # 5219 genes with TRUE in all 62 samples 
keep.d21 <- rowSums(thresh.d21) >= (ncol(thresh.d21)/2) # we would like to keep genes that have at least 50% TRUES in each row of thresh
summary(keep.d21) # FALSE 21104 & TRUE 13843 -- more than three quarters of the genes did not pass
cts.matrix.d21.filtered <- cts.matrix.d21[keep.d21,] # Subset the rows of countdata to keep the more highly expressed genes
dim(cts.matrix.d21.filtered) # 13843 genes &  62 samples
# head(cts.matrix.d21.filtered) # FINAL DAY 21 DATASET
# We will look at a cople samples to check our defined threshold.. 
# plot(CPM.d21[,60], cts.matrix.d21[,60], xlab="CPM", ylab="Raw Count", ylim=c(0,50), xlim=c(0,50), main=colnames(CPM.d21)[1]) # check if our threshold of 3 CPM in 50% samples corresponds to a count of about 10-15
# abline(v=3) 
# abline(h=3)
# plot(CPM.d21[,5], cts.matrix.d21[,5], xlab="CPM", ylab="Raw Count", ylim=c(0,50), xlim=c(0,50), main=colnames(CPM.d21)[1]) # check if our threshold of 3 CPM in 50% samples corresponds to a count of about 10-15
# abline(v=3) 
# abline(h=3)

# write csv
write.csv(cts.matrix.d21.filtered,paste(path,"day21.counts.filtered_matrix.csv")) # 'path' called in previous # write .csv section


