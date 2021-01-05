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
# ALL TREATMENTS AND TIMEPOINTS

dds.all.treatments <- DESeq(dds.all.treatments) # wait for this to complete....
resultsNames(dds.all.treatments) # view the names of your results model 
# Note: these are all different ways of looking at the same results... 
res.ALL.TRMTS <- results(dds.all.treatments) # view DESeq2 result - calls the last pairwise condition 
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

res.A.d0_AAA.d21.NS.GENES <- res.A.d0_AAA.d21.DATAFRAME %>%  dplyr::filter(padj > 0.05)
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
# DAY 0

dds.primary.d0 <- DESeq(dds.primary.d0) # wait for this to complete....
resultsNames(dds.primary.d0) # view the names of your results model 'Primary_Treatment_M_vs_A'

# Note: these are all different ways of looking at the same results... 
res.d0 <- results(dds.primary.d0) # view DESeq2 results 
res.d0 <- results(dds.primary.d0, name="Primary_Treatment_M_vs_A")
res.d0.contrasts <- results(dds.primary.d0, contrast=c("Primary_Treatment","M","A")) # alternative way of calling the results (i.e. Day_DAY21_vs_Day0)

res.d0_Ordered <- res.d0[order(res.d0$pvalue),] # how to order based on p-valuel 
res.d0_Ordered # ordered based on pvalue
sum(res.d0_Ordered$padj < 0.1, na.rm=TRUE) # 27 marginal DEGs
sum(res.d0_Ordered$padj < 0.05, na.rm=TRUE) # 14 significant DEGs

summary(res.d0)

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
# DAY 21 

dds.primary.d21 <- DESeq(dds.primary.d21) # wait for this to complete....
resultsNames(dds.primary.d21) # view the names of your results model 'Primary_Treatment_M_vs_A'

# Note: these are all different ways of looking at the same results... 
res.d21 <- results(dds.primary.d21) # view DESeq2 results 
res.d21 <- results(dds, name="Primary_Treatment_M_vs_A")
res.D21.contrasts <- results(dds.primary.d21, contrast=c("Primary_Treatment","M","A")) # alternative way of calling the results (i.e. Day_DAY21_vs_Day0)

res.d21_Ordered <- res.d21[order(res.d21$pvalue),] # how to order based on p-value
res.d21_Ordered # ordered based on pvalue
sum(res.d21_Ordered$padj < 0.1, na.rm=TRUE) # 331 marginal DEGs
sum(res.d21_Ordered$padj < 0.05, na.rm=TRUE) # 196 significant DEGs

summary(res.d21)

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
