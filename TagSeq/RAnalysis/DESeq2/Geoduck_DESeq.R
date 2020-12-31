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

#======================================== #
# DAY 0 
# About: run dyplr 'antijoin' to call cts columns that match 'Sample.Name' in the data frame 'exp.data.d0'
cts.merged.as.table <- data.frame(transcript_id = row.names(cts.merged), cts.merged) # add back the rownames 'transcript_ID'
rownames(cts.merged.as.table) <- NULL # ommit the rownames
ncol(cts.merged.as.table) # 142 counting the ID that we want to keep! 
cts.merged.d0 <- cts.merged.as.table[,c(1,na.omit(match(exp.data.d0$Sample.Name, colnames(cts.merged.as.table))))]
cts.merged.d0 <- data.frame(cts.merged.d0[,-1], row.names=cts.merged.d0[,1])
cts.matrix.d0  <-as.matrix(cts.merged.d0, row.names="transcript_id")
ncol(cts.matrix.d0) # 8  samples from just Day 0
colnames(cts.matrix.d0) == exp.data.d0$Sample.Name # check if TRUE, means the same as the exp/design dataframe exp.data.d0

#======================================== #
# DAY 21 
cts.merged.as.table <- data.frame(transcript_id = row.names(cts.merged), cts.merged) # add back the rownames 'transcript_ID'
rownames(cts.merged.as.table) <- NULL # ommit the rownames
ncol(cts.merged.as.table) # 142 counting the ID that we want to keep! 
cts.merged.d21 <- cts.merged.as.table[,c(1,na.omit(match(exp.data.d21$Sample.Name, colnames(cts.merged.as.table))))]
cts.merged.d21 <- data.frame(cts.merged.d21[,-1], row.names=cts.merged.d21[,1])
cts.matrix.d21  <-as.matrix(cts.merged.d21, row.names="transcript_id")
ncol(cts.matrix.d21) # # 62 total sampels on day 21
colnames(cts.matrix.d21) == exp.data.d21$Sample.Name # check if TRUE, means the same as the exp/design dataframe exp.data.d21



#=====================================================================================
#  PREP dds files for DESeq2    ~ 'Primary_Treatment' 
#=====================================================================================

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
exp.data.d0.PRIMARY$Time <- factor(exp.data.d0.PRIMARY$Time) # change Primary_Treatment to factor
exp.data.d0.PRIMARY <- data.frame(exp.data.d0.PRIMARY[,-1], row.names=exp.data.d0.PRIMARY[,1]) # move Sample.Name column as row names  
exp.data.d0.PRIMARY.mtx <- as.matrix(exp.data.d0.PRIMARY, row.names="Geoduck.ID") # create matrix 
# fix(exp.data.mtx) # view data - remove # to open
# CKECK THE cts.matrix AND THE exp.data.mtx
exp.data.d0.PRIMARY.mtx <- exp.data.d0.PRIMARY.mtx[match(colnames(cts.merged.d0),rownames(exp.data.d0.PRIMARY.mtx)), ]
all(rownames(exp.data.d0.PRIMARY.mtx) %in% colnames(cts.matrix.d0)) # should be TRUE
all(rownames(exp.data.d0.PRIMARY.mtx) == colnames(cts.matrix.d0)) # should be TRUE
cts.matrix.d0 <- cts.matrix.d0[, rownames(exp.data.d0.PRIMARY.mtx)]
all(rownames(exp.data.d0.PRIMARY.mtx) == colnames(cts.matrix.d0))  # should be TRUE

# build dds
dds.primary.d0 <- DESeqDataSetFromMatrix(countData = cts.matrix.d0,
                                      colData = exp.data.d0.PRIMARY.mtx,
                                      design = ~ Primary_Treatment) # DESeq Data Set (dds) - design as ~Primary_Treatment
dds.primary.d0 # view dds

# prep for DESeq2
dds.primary.d0$Primary_Treatment <- relevel(dds.primary.d0$Primary_Treatment, ref = "A") # specify the reference level for count analysis - A = the control treatment
dds.primary.d0$Time <- droplevels(dds.primary.d0$Time)
levels(dds.primary.d0$Time) # NULL
levels(dds.primary.d0$Primary_Treatment) #  levels for condition 'Treatment': "A" "M"
keep <- rowSums(counts(dds.primary.d0)) >= 5 # PRE-FILTERING... ommits genes where counts between ALL samples (as rows) less than 5 
dds.primary.d0 <- dds.primary.d0[keep,] # integrate this criteria in the data 
design(dds.primary.d0) # view the design we have specified 
ncol(dds.primary.d0) # 8 cols - Good
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
dds.primary.d21 <- dds.primary.d21[keep,] # integrate this criteria in the data 
design(dds.primary.d21) # view the design we have specified 
ncol(dds.primary.d21) # 62 cols - Good
colData(dds.primary.d21)

# DEG ANALYSIS READY!!


#=====================================================================================
# Differential expression analysis
#=====================================================================================
# About: Log fold change shrinkage for visualization and ranking
# Run DESeq : Modeling counts with the specified design(dds) treatment effects 


#======================================================================================= #
# ALL TIME POINTS

dds.primary <- DESeq(dds.primary) # wait for this to complete....
resultsNames(dds.primary) # view the names of your results model 'Primary_Treatment_M_vs_A'

# Note: these are all different ways of looking at the same results... 
res <- results(dds.primary) # view DESeq2 results 
res <- results(dds.primary, name="Primary_Treatment_M_vs_A")
res.contrasts <- results(dds.primary, contrast=c("Primary_Treatment","M","A")) # alternative way of calling the results (i.e. Day_DAY21_vs_Day0)


res_Ordered <- res[order(res$pvalue),] # how to order based on p-value; chosses EHSM vs A - the last pairwise call 
res_Ordered # ordered based on pvalue
sum(res_Ordered$padj < 0.1, na.rm=TRUE) # 1417 marginal DEGs
sum(res_Ordered$padj < 0.05, na.rm=TRUE) # 1071 significant DEGs

summary(res)

#======================================================================================= #
# DAY 0

dds.primary.d0 <- DESeq(dds.primary.d0) # wait for this to complete....
resultsNames(dds.primary.d0) # view the names of your results model 'Primary_Treatment_M_vs_A'

# Note: these are all different ways of looking at the same results... 
res.d0 <- results(dds.primary.d0) # view DESeq2 results 
res.d0 <- results(dds.primary.d0, name="Primary_Treatment_M_vs_A")
res.d0.contrasts <- results(dds.primary.d0, contrast=c("Primary_Treatment","M","A")) # alternative way of calling the results (i.e. Day_DAY21_vs_Day0)


res.d0_Ordered <- res.d0[order(res.d0$pvalue),] # how to order based on p-value; chosses EHSM vs A - the last pairwise call 
res.d0_Ordered # ordered based on pvalue
sum(res.d0_Ordered$padj < 0.1, na.rm=TRUE) # 27 marginal DEGs
sum(res.d0_Ordered$padj < 0.05, na.rm=TRUE) # 14 significant DEGs

summary(res.d0)

#======================================================================================= #
# DAY 21 

dds.primary.d21 <- DESeq(dds.primary.d21) # wait for this to complete....
resultsNames(dds.primary.d21) # view the names of your results model 'Primary_Treatment_M_vs_A'

# Note: these are all different ways of looking at the same results... 
res.d21 <- results(dds.primary.d21) # view DESeq2 results 
res.d21 <- results(dds, name="Primary_Treatment_M_vs_A")
res.D21.contrasts <- results(dds.primary.d21, contrast=c("Primary_Treatment","M","A")) # alternative way of calling the results (i.e. Day_DAY21_vs_Day0)

res.d21_Ordered <- res.d21[order(res.d21$pvalue),] # how to order based on p-value; chosses EHSM vs A - the last pairwise call 
res.d21_Ordered # ordered based on pvalue
sum(res.d21_Ordered$padj < 0.1, na.rm=TRUE) # 331 marginal DEGs
sum(res.d21_Ordered$padj < 0.05, na.rm=TRUE) # 196 significant DEGs

summary(res.d21)

#ddsMF <- dds # make a coppy of dds for multi factor designs

#=====================================================================================
# EXPLORING AND EXPORTING RESULTS
#=====================================================================================
summary(res.EvA)
plotMA(res.EvA, ylim=c(-2,2))

# It is more useful visualize the MA-plot for the shrunken log2 fold changes, 
# which remove the noise associated with log2 fold changes from low count genes 
# without requiring arbitrary filtering thresholds.
res.EvA_LFC <- lfcShrink(dds, coef="Treatment_E_vs_A", type="apeglm")
summary(res.EvA_LFC)
plotMA(res.EvA_LFC, ylim=c(-2,2))
# after calling plotMA you can use function identify to interactively detect the row number of individual genes by clicking on the plot
idx <- identify(res.EvA_LFC$baseMean, res.EvA_LFC$log2FoldChange)
rownames(res.EvA_LFC)[idx]

# plot counts
plotCounts(dds, gene=which.min(res.EvA_LFC$padj), intgroup="Treatment") # this plots the gene which had the smallest p-value
# sample plot in ggplot 
d <- plotCounts(dds, gene=which.min(res.EvA_LFC$padj), intgroup="Treatment", 
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
