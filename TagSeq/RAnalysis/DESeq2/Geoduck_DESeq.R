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

# SET WORKING DIRECTORY AND LOAD DATA
setwd("C:/Users/samjg/Documents/My_Projects/Pgenerosa_TagSeq_Metabolomics/TagSeq/")

# LOAD DATA
cts <- read.csv(file="HPC_Bioinf/outputs/transcript_count_matrix.csv", sep=',', header=TRUE) # read the output count matrix from prepDE.py

UT_seq_map <- read.csv(file="20201020_Gurr_TagSeq_UTAustin.csv", sep=',', header=TRUE)
smpl_ref <- read.csv(file="Sample_reference.csv", sep=',', header=TRUE)
treatment_ref <- read.csv(file="Extraction_checklist.csv", sep=',', header=TRUE)

# FORMAT COUNT MATRIX (cts has 2x columns per sample on different lanes - sum counts from unique columns + build matrix)
ncol(cts) # 282 samples (not counting gene ID column) - should be 141 samples, need to sum columns by unique ID
cts.merged <- data.frame(cts[,-1], row.names=cts[,1]) # call new dataframe with first column now as row names, now all row values are numeric
names(cts.merged) <- sapply(strsplit(names(cts.merged), "_"), '[', 1) # split the column names by "_" delimiter  and call the first fielf SG##
cts.merged <- t(rowsum(t(cts.merged), group = colnames(cts.merged), na.rm = T)) # merge all unique columns and sum counts 
cts.matrix <-as.matrix(cts.merged, row.names="transcript_id") # call dataframe as matrix
fix(cts.matrix) # view count matrix with abbreviated columns and merged data
ncol(cts.matrix) # 141 samples
nrow(cts.matrix) # 34947 total genes
head(cts.matrix,2)

# MASTER REFERENCE DATA.FRAME
# format and merge to buld master reference dataframe
smpl_ref$Seq_Pos <- paste(smpl_ref$ï..TagSeq_Plate, smpl_ref$TagSeq_Well, sep="_")
smpl_ref <- smpl_ref %>% select(-c("ï..TagSeq_Plate", "TagSeq_Well"))
UT_seq_map$Seq_Pos  <- paste(UT_seq_map$ï..Plate, UT_seq_map$Well, sep="_")
UT_seq_map <- UT_seq_map %>% select(-c("ï..Plate", "Well"))
Seq.Ref <- merge(smpl_ref, UT_seq_map, by = "Seq_Pos")
treatment_ref <- treatment_ref[, c(3:10)]
Mstr.Ref <- merge(Seq.Ref, treatment_ref, by = "Geoduck_ID")
# fix(Mstr.Ref)

# call all experiment design treatments
exp.data <- Mstr.Ref[,c("Sample.Name","ALL_Treatment", "EvA_Treatment", "Day0_Treament", "Day7_Treament", "Day14", "Day21_Treatment", "Day")]


##### TREATMENT = Day0_Treament ##### #
# assign 'Treatment' as EvA_Treatment containing only E and A for the Day 0 
exp.data$ALL_Treatment <- factor(exp.data$ALL_Treatment) 
exp.data$Day <- factor(exp.data$Day)

exp.data <- exp.data %>% dplyr::select(c('Sample.Name', 'Day', 'ALL_Treatment'))
exp.data <- data.frame(exp.data[,-1], row.names=exp.data[,1]) 
exp.data.mtx <- as.matrix(exp.data, row.names="Geoduck.ID")
# fix(exp.data.mtx)

# CKECK THE cts.matrix AND THE exp.data.mtx
exp.data.mtx <- exp.data.mtx[match(colnames(cts.merged),rownames(exp.data.mtx)), ]
all(rownames(exp.data.mtx) %in% colnames(cts.matrix)) # should be TRUE
all(rownames(exp.data.mtx) == colnames(cts.matrix)) # should be TRUE
cts.matrix <- cts.matrix[, rownames(exp.data.mtx)]
all(rownames(exp.data.mtx) == colnames(cts.matrix))  # should be TRUE

# DESeq Data Set (dds) - 'Treatment' as the design 
dds <- DESeqDataSetFromMatrix(countData = cts.matrix,
                              colData = exp.data.mtx,
                              design = ~ ALL_Treatment)
dds

# Differential expression analysis
# Log fold change shrinkage for visualization and ranking
dds <- DESeq(dds)
res <- results(dds)
resultsNames(dds)

res.EvA <- lfcShrink(dds, coef="Treatment_E_vs_A", type="apeglm")
res.EvA.dataframe <- as.data.frame(res.EvA)
sig.genes.EvA <- res.EvA.dataframe %>% dplyr::filter(pvalue < 0.05)


resOrdered <- res[order(res$pvalue),] # how to order based on p-value; chosses EHSM vs A - the last pairwise call 
sum(res$padj < 0.1, na.rm=TRUE) # how many p values (genes) were less than 0.1? (49)
sum(res$padj < 0.05, na.rm=TRUE) # 0.05? (21)
res05 <- results(dds, alpha=0.05)
summary(res05)
sum(res05$padj < 0.05, na.rm=TRUE)

colData(dds)
ddsMF <- dds # make a coppy of dds for multi factor designs


library("pheatmap")



