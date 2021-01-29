---
# title: "Geoduck_TagSeq_DESeq2"
# author: "Samuel Gurr"
# date: "January 24, 2021"
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
library(tidybulk)

# SET WORKING DIRECTORY AND LOAD DATA
setwd("C:/Users/samjg/Documents/My_Projects/Pgenerosa_TagSeq_Metabolomics/TagSeq/")

# LOAD DATA
# filtered counts tables - format matrix after upload [from Count_Matrix_Stats.Filter.R] 
all.counts_matrix <- read.csv(file="RAnalysis/DESeq2/counts_filtered/ all.counts.filtered_matrix.csv", sep=',', header=TRUE) 
all.counts_matrix <- data.frame(all.counts_matrix[,-1], row.names=all.counts_matrix[,1]) 
d0.counts_matrix  <- read.csv(file="RAnalysis/DESeq2/counts_filtered/ day0.counts.filtered_matrix.csv", sep=',', header=TRUE) 
d0.counts_matrix <- data.frame(d0.counts_matrix[,-1], row.names=d0.counts_matrix[,1]) 
d7.counts_matrix  <- read.csv(file="RAnalysis/DESeq2/counts_filtered/ day7.counts.filtered_matrix.csv", sep=',', header=TRUE) 
d7.counts_matrix <- data.frame(d7.counts_matrix[,-1], row.names=d7.counts_matrix[,1]) 
d14.counts_matrix <- read.csv(file="RAnalysis/DESeq2/counts_filtered/ day14.counts.filtered_matrix.csv", sep=',', header=TRUE) 
d14.counts_matrix <- data.frame(d14.counts_matrix[,-1], row.names=d14.counts_matrix[,1]) 
d21.counts_matrix <- read.csv(file="RAnalysis/DESeq2/counts_filtered/ day21.counts.filtered_matrix.csv", sep=',', header=TRUE) 
d21.counts_matrix <- data.frame(d21.counts_matrix[,-1], row.names=d21.counts_matrix[,1]) 
# experiment data [from Count_Matrix_Stats.Filter.R] 
all.exp_data  <- read.csv(file="RAnalysis/DESeq2/counts_filtered/ all.exp.data.csv", sep=',', header=TRUE) 
d0.exp_data   <- read.csv(file="RAnalysis/DESeq2/counts_filtered/ day0.exp.data.csv", sep=',', header=TRUE) 
d7.exp_data   <- read.csv(file="RAnalysis/DESeq2/counts_filtered/ day7.exp.data.csv", sep=',', header=TRUE) 
d14.exp_data  <- read.csv(file="RAnalysis/DESeq2/counts_filtered/  day14.exp.data.csv", sep=',', header=TRUE) 
d21.exp_data  <- read.csv(file="RAnalysis/DESeq2/counts_filtered/  day21.exp.data.csv", sep=',', header=TRUE) 

#====================================================================================================================== 
#
#                             PREP dds files for DESeq2    
#
#====================================================================================================================== 



# ========================================================== 
# DAY 0 FULL MODEL ==  design = ~ Primary_Treatment
#---- INPUT: Experiment/design matrix == 'exp.data.d0' 
#---- INPUT: TagSeq count matrix == 'cts.matrix.d0.filtered' (3 CPM in 50% samples)
# ========================================================== 

# format Experiment/design dataframe into matrix
exp.data.d0.PRIMARY <- d0.exp_data %>% dplyr::select(c('Sample.Name', 'Primary_Treatment', 'All_Treatment')) # coondense dataset to build target matrix
exp.data.d0.PRIMARY$Primary_Treatment <- factor(exp.data.d0.PRIMARY$Primary_Treatment) # change Primary_Treatment to factor
exp.data.d0.PRIMARY$All_Treatment <- factor(exp.data.d0.PRIMARY$All_Treatment) # change Time to factor
exp.data.d0.PRIMARY <- data.frame(exp.data.d0.PRIMARY[,-1], row.names=exp.data.d0.PRIMARY[,1]) # move Sample.Name column as row names  
exp.data.d0.PRIMARY.mtx <- as.matrix(exp.data.d0.PRIMARY, row.names="Geoduck.ID") # create matrix 
# fix(exp.data.mtx) # view data - remove # to open
# CKECK THE cts.matrix AND THE exp.data.mtx
exp.data.d0.PRIMARY.mtx <- exp.data.d0.PRIMARY.mtx[match(colnames(d0.counts_matrix),rownames(exp.data.d0.PRIMARY.mtx)), ]
all(rownames(exp.data.d0.PRIMARY.mtx) %in% colnames(d0.counts_matrix)) # should be TRUE
all(rownames(exp.data.d0.PRIMARY.mtx) == colnames(d0.counts_matrix)) # should be TRUE
all(rownames(exp.data.d0.PRIMARY.mtx) == colnames(d0.counts_matrix))  # should be TRUE

# build dds
dds.d0 <- DESeqDataSetFromMatrix(countData = d0.counts_matrix,
                                      colData = exp.data.d0.PRIMARY.mtx,
                                      design = ~ Primary_Treatment) # DESeq Data Set (dds) - design as ~Primary_Treatment
dds.d0 # view dds

# prep for DESeq2
dds.d0$Primary_Treatment <- relevel(dds.d0$Primary_Treatment, ref = "A") # specify the reference level for count analysis - A = the control treatment
levels(dds.d0$Primary_Treatment) #  levels for condition 'Treatment': "A" "M"

design(dds.d0) # view the design we have specified 
ncol(dds.d0) # 8 cols - Good
nrow(dds.d0) # 15149 (cut-off of 3 CPM)
colData(dds.d0) # view the col data

# DEG ANALYSIS READY!!

# ========================================================== 
# DAY 7 - FULL MODEL == design = ~ Primary_Treatment+Second_Treament+Primary_Treatment:Second_Treament
#---- INPUT: Experiment/design matrix == 'exp.data.d7' 
#---- INPUT: TagSeq count matrix == 'cts.matrix.d7.filtered' (3 CPM in 50% samples)
# ==========================================================
d7.exp_data$Primary_Treatment <- factor(d7.exp_data$Primary_Treatment) # change Primary_Treatment to factor
d7.exp_data$Second_Treament <- factor(d7.exp_data$Second_Treament) # change Second_Treament to factor
d7.exp_data$All_Treatment <- factor(d7.exp_data$All_Treatment) # change Time to factor

#==================== #
# Primary*Second (ALL) Treatment
# format Experiment/design dataframe into matrix
exp.data.d7.PRIM_SEC <- d7.exp_data %>% dplyr::select(c('Sample.Name', 'Primary_Treatment', 'Second_Treament', 'All_Treatment')) # coondense dataset to build target matrix
exp.data.d7.PRIM_SEC <- data.frame(exp.data.d7.PRIM_SEC[,-1], row.names=exp.data.d7.PRIM_SEC[,1]) # move Sample.Name column as row names  
exp.data.d7.PRIM_SEC.mtx <- as.matrix(exp.data.d7.PRIM_SEC, row.names="Geoduck.ID") # create matrix 
# CKECK THE cts.matrix AND THE exp.data.mtx
exp.data.d7.PRIM_SEC.mtx <- exp.data.d7.PRIM_SEC.mtx[match(colnames(d7.counts_matrix),rownames(exp.data.d7.PRIM_SEC.mtx)), ]
all(rownames(exp.data.d7.PRIM_SEC.mtx) %in% colnames(d7.counts_matrix)) # should be TRUE
all(rownames(exp.data.d7.PRIM_SEC.mtx) == colnames(d7.counts_matrix)) # should be TRUE
all(rownames(exp.data.d7.PRIM_SEC.mtx) == colnames(d7.counts_matrix))  # should be TRUE
# build dds
# FULL MODEL 
dds.d7 <- DESeqDataSetFromMatrix(countData = d7.counts_matrix,
                                  colData = exp.data.d7.PRIM_SEC,
                                  design = ~ Second_Treament+Primary_Treatment+Second_Treament:Primary_Treatment) # DESeq Data Set (dds) - design as ~Primary_Treatment
design(dds.d7) # ~Second_Treament + Primary_Treatment + Second_Treament:Primary_Treatment
# GROUP
dds.d7.group <- DESeqDataSetFromMatrix(countData = d7.counts_matrix,
                                        colData = exp.data.d7.PRIM_SEC,
                                        design = ~ All_Treatment-1) 
design(dds.d7.group) # ~All_Treatment 
# MAIN EFFECT 
dds.d7.main <- DESeqDataSetFromMatrix(countData = d7.counts_matrix,
                                       colData = exp.data.d7.PRIM_SEC,
                                       design = ~ Second_Treament+Primary_Treatment) 
design(dds.d7.main) # ~Second_Treament + Primary_Treatment

colData(dds.d7) # view the col data

# DEG ANALYSIS READY!!

# ========================================================== 
# DAY 14 -  FULL MODEL == design = ~ Primary_Treatment+Second_Treament+Primary_Treatment:Second_Treament
#---- INPUT: Experiment/design matrix == 'exp.data.d14' 
#---- INPUT: TagSeq count matrix == 'cts.matrix.d14.filtered' (3 CPM in 50% samples)
# ========================================================== 
# format Experiment/design dataframe into matrix
exp.data.d14.PRIM_SEC <- d14.exp_data %>% dplyr::select(c('Sample.Name', 'Primary_Treatment', 'Second_Treament', 'All_Treatment')) # coondense dataset to build target matrix
exp.data.d14.PRIM_SEC <- data.frame(exp.data.d14.PRIM_SEC[,-1], row.names=exp.data.d14.PRIM_SEC[,1]) # move Sample.Name column as row names  
exp.data.d14.PRIM_SEC.mtx <- as.matrix(exp.data.d14.PRIM_SEC, row.names="Geoduck.ID") # create matrix 
# CKECK THE cts.matrix AND THE exp.data.mtx
exp.data.d14.PRIM_SEC.mtx <- exp.data.d14.PRIM_SEC.mtx[match(colnames(d14.counts_matrix),rownames(exp.data.d14.PRIM_SEC.mtx)), ]
all(rownames(exp.data.d14.PRIM_SEC.mtx) %in% colnames(d14.counts_matrix)) # should be TRUE
all(rownames(exp.data.d14.PRIM_SEC.mtx) == colnames(d14.counts_matrix)) # should be TRUE
all(rownames(exp.data.d14.PRIM_SEC.mtx) == colnames(d14.counts_matrix))  # should be TRUE
# build dds
# FULL MODEL 
dds.d14 <- DESeqDataSetFromMatrix(countData = d14.counts_matrix,
                                     colData = exp.data.d14.PRIM_SEC.mtx,
                                     design = ~ Second_Treament+Primary_Treatment+Second_Treament:Primary_Treatment) # DESeq Data Set (dds) - design as ~Primary_Treatment
design(dds.d14) # ~Second_Treament + Primary_Treatment + Second_Treament:Primary_Treatment
# GROUP
dds.d14.group <- DESeqDataSetFromMatrix(countData = d14.counts_matrix,
                                        colData = exp.data.d14.PRIM_SEC.mtx,
                                        design = ~ All_Treatment-1) 
design(dds.d14.group) # ~All_Treatment 
# MAIN EFFECT 
dds.d14.main <- DESeqDataSetFromMatrix(countData = d14.counts_matrix,
                                       colData = exp.data.d14.PRIM_SEC.mtx,
                                       design = ~ Second_Treament+Primary_Treatment) 
design(dds.d14.main) # ~Second_Treament + Primary_Treatment


# DEG ANALYSIS READY!!

# ========================================================== 
# DAY 21 -  FULL MODEL == design = ~Primary_Treatment+Second_Treament+Third_Treatment+Primary_Treatment:Second_Treament+Primary_Treatment:Third_Treatment+Second_Treament:Third_Treatment
#---- INPUT: Experiment/design matrix == 'exp.data.d21'
#---- INPUT: TagSeq count matrix == 'cts.matrix.d21.filtered' (3 CPM in 50% samples)
# ========================================================== 
# format Experiment/design dataframe into matrix
exp.data.d21.PRIM_SEC_THIRD <- d21.exp_data %>% dplyr::select(c('Sample.Name', 'Primary_Treatment','Second_Treament', 'Third_Treatment', 'All_Treatment')) # coondense dataset to build target matrix
exp.data.d21.PRIM_SEC_THIRD <- data.frame(exp.data.d21.PRIM_SEC_THIRD[,-1], row.names=exp.data.d21.PRIM_SEC_THIRD[,1]) # move Sample.Name column as row names  
exp.data.d21.PRIM_SEC_THIRD.mtx <- as.matrix(exp.data.d21.PRIM_SEC_THIRD, row.names="Geoduck.ID") # create matrix 
# CKECK THE cts.matrix AND THE exp.data.mtx
exp.data.d21.PRIM_SEC_THIRD.mtx <- exp.data.d21.PRIM_SEC.mtx[match(colnames(d21.counts_matrix),rownames(exp.data.d21.PRIM_SEC_THIRD.mtx)), ]
all(rownames(exp.data.d21.PRIM_SEC_THIRD.mtx) %in% colnames(d21.counts_matrix)) # should be TRUE
all(rownames(exp.data.d21.PRIM_SEC_THIRD.mtx) == colnames(d21.counts_matrix)) # should be TRUE
all(rownames(exp.data.d21.PRIM_SEC_THIRD.mtx) == colnames(d21.counts_matrix))  # should be TRUE
# build dds
# full model
dds.d21 <- DESeqDataSetFromMatrix(countData = d21.counts_matrix,
                                        colData = exp.data.d21.PRIM_SEC_THIRD.mtx,
                                        design = ~ Primary_Treatment+
                                        Second_Treament+
                                        Third_Treatment+
                                        Primary_Treatment:Second_Treament+
                                        Primary_Treatment:Third_Treatment+
                                        Second_Treament:Third_Treatment) 
# GROUP
dds.d21.group <- DESeqDataSetFromMatrix(countData = d21.counts_matrix,
                                        colData = exp.data.d21.PRIM_SEC_THIRD.mtx,
                                        design = ~ All_Treatment-1) 
# MAIN EFFECT 
dds.d21.main <- DESeqDataSetFromMatrix(countData = d21.counts_matrix,
                                       colData = exp.data.d21.PRIM_SEC_THIRD.mtx,
                                       design = ~ Third_Treatment+Second_Treament+Primary_Treatment) 



# prep for DESeq2
# dds.d21$Primary_Treatment <- relevel(dds.d21$Primary_Treatment, ref = "A") # specify the reference level for count analysis - A = the control treatment
# dds.d21$Second_Treament <- relevel(dds.d21$Second_Treament, ref = "A") 
# dds.d21$Third_Treatment <- relevel(dds.d21$Third_Treatment, ref = "A") 
# levels(dds.d21$Primary_Treatment) #  levels for condition 'Treatment': "A" "M"
# levels(dds.d21$Second_Treament) #  levels for condition 'Treatment': "A" "M" "S"
# levels(dds.d21$Third_Treatment) #  levels for condition 'Treatment': "A" "M"
# design(dds.d21) # view the design we have specified 
# ncol(dds.d21) #  62 cols - Good
# nrow(dds.d21) #13743 (3 CPM)

# DEG ANALYSIS READY!!

#====================================================================================================================== 
#
#                          Differential expression analysis (DESeq2)
#
#====================================================================================================================== 

# ================== #
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

# PRIMARY MODERATE VERSUS AMBEINT HISTORY UNDER THE ENTIRE EXPERIMENT  ========================================================== #
dds.ALL_MvA
#  pre-filtered in edgeR - use independentFiltering=FALSE
dds.ALL_MvA <- DESeq(dds.ALL_MvA) # wait for this to complete....
resultsNames(dds.ALL_MvA) # view the names of your results model "MvA_Primary_History_M_vs_A"
# count  DEGs with pdj threshold < 0.05 and LFC >1 (up) and < -1 (down)
res.ALL.MvA <- results(dds.ALL_MvA, name="MvA_Primary_History_M_vs_A", alpha = 0.05)
hist(res.ALL.MvA$pvalue, breaks=20, col="grey") # view histogram,  about 425 genes are likely false positives, alpha = 0.0
abline(h=(nrow(res.ALL.MvA)*0.05),col="red") # add line at expected 5% false positive
table(res.ALL.MvA$padj<0.05) # 971 DEGs total with padj value < 0.05
res.ALL.MvA <- res.ALL.MvA[order(res.ALL.MvA$padj), ] # Order by adjusted p-value
sum((res.ALL.MvA$log2FoldChange[1:971] >= 1) == TRUE) # 57 DEGs upregulated  (LFC >= 1)
( sum((res.ALL.MvA$log2FoldChange[1:971] >= 1) == TRUE) / (table(res.ALL.MvA$padj<0.05))[2] ) * 100 # 5.870237 % DEGs upregulated (N = 971)
sum((res.ALL.MvA$log2FoldChange[1:971] <= -1) == TRUE) # 255 DEGs downregulated  (LFC <= -1)
( sum((res.ALL.MvA$log2FoldChange[1:971] <= -1) == TRUE) / (table(res.ALL.MvA$padj<0.05))[2] ) * 100 # 26.26159 % DEGs downregulated (N = 971)
# create a dataframe of results 
resdata.ALL.MvA <- merge(as.data.frame(res.ALL.MvA), as.data.frame(counts(dds.ALL_MvA, normalized=TRUE)), by="row.names", sort=FALSE) ## Merge with normalized count data
names(resdata.ALL.MvA)[1] <- "Gene"
# Write results
path_out.d0 = 'C:/Users/samjg/Documents/My_Projects/Pgenerosa_TagSeq_Metabolomics/TagSeq/RAnalysis/DESeq2/output/'
write.csv(resdata.ALL.MvA, paste(path_out.d0, "ALL_Mod_vs_Amb_diffexpr-results.csv"))
# volcano plot ------------------------------------------------------------------------------------------------------ #
png("RAnalysis/DESeq2/output/ALL_ModvsAmb.PrimaryTreatment-VolcanoPlot.png", 1000, 1000, pointsize=20)
EnhancedVolcano(res.ALL.MvA,
                lab = rownames(res.ALL.MvA),
                x = 'log2FoldChange',
                y = 'pvalue',
                title = 'Moderate versus Ambient Treatment (Primary)',
                subtitle = "DESeq2 - Differential expression",
                FCcutoff = 1.0,
                pCutoff = 0.3,
                pointSize = 4.0,
                labSize = 6.0,
                colAlpha = 1,
                legendPosition = 'right',
                legendLabSize = 12,
                legendIconSize = 4.0,
                widthConnectors = 0.75)
dev.off()
# Plot dispersions ------------------------------------------------------------------------------------------------------ #
png("RAnalysis/DESeq2/output/ALL_ModvsAmb.dispersions.png", 1000, 1000, pointsize=20)
plotDispEsts(dds.ALL_MvA, main="All MvA dispersions")
dev.off()
# Heat maps and Principal components analysis ------------------------------------------------------------------------------------------ #
# significantly upregulated genes - padj < 0.05; LFC > 2
res.ALL.MvA_UPREG <- res.ALL.MvA[c(1:971),] 
upreg.T_F <- res.ALL.MvA_UPREG$log2FoldChange >= 1 # we would like to keep genes that have at least 50% TRUES in each row of thresh
res.ALL.MvA_UPREG <- res.ALL.MvA_UPREG[upreg.T_F,]  # call onl upreg genes
dds_UPREG <- dds.ALL_MvA[(rownames(res.ALL.MvA_UPREG))] # call only res.ALL.MvA_UPREG rows (12) 
dim(dds_UPREG) # 57 141 - 57 pass this criteria
rlog.dds_UPREG<- rlogTransformation(dds_UPREG) # Conmplete rlog transformation ....takes a bit of time
sampleDist.rlog.dds_UPREG <- dist(t(assay(rlog.dds_UPREG))) # next assay [access dds data matrix] and t() [transpose rows and columns] in dist() [calculate the Euclidean Distance]
sampleDistMatrix.rlog.dds_UPREG <- as.matrix(sampleDist.rlog.dds_UPREG) # convert to a data matrix for heatmap 
rownames(sampleDistMatrix.rlog.dds_UPREG) <- dds_UPREG$MvA_Primary_History # assign the rownames as the conditions of your samples
pheatmap(sampleDistMatrix.rlog.dds_UPREG, 
         clustering_distance_rows = sampleDist.rlog.dds_UPREG, 
         clustering_distance_cols = sampleDist.rlog.dds_UPREG, 
         col = colors) # create the heatmap 
plotPCA(rlog.dds_UPREG,intgroup="MvA_Primary_History")  # use the DESeq2 call plotPCA() and call the rlog_dds for PCA plotting - Q: difference between replicates or condition explains the most variance in this dataset?


# significantly downregulated genes - padj < 0.05; LFC < 2
res.ALL.MvA_DWNREG <- res.ALL.MvA[c(1:971),] 
dwnreg.T_F <- res.ALL.MvA_DWNREG$log2FoldChange <= -1 # we would like to keep genes that have at least 50% TRUES in each row of thresh
res.ALL.MvA_DWNREG <- res.ALL.MvA_DWNREG[dwnreg.T_F,] # call only downreg genes
dds_DWNREG <- dds.ALL_MvA[(rownames(res.ALL.MvA_DWNREG))] # call only res.ALL.MvA_UPREG rows (12) 
dim(dds_DWNREG) # 255 141 - 255 pass this criteria
rlog.dds_DWNREG<- rlogTransformation(dds_DWNREG) # Conmplete rlog transformation  ....takes a bit of time
sampleDist.rlog.dds_DWNREG <- dist(t(assay(rlog.dds_DWNREG))) # next assay [access dds data matrix] and t() [transpose rows and columns] in dist() [calculate the Euclidean Distance]
sampleDistMatrix.rlog.dds_DWNREG <- as.matrix(sampleDist.rlog.dds_DWNREG) # convert to a data matrix for heatmap 
rownames(sampleDistMatrix.rlog.dds_DWNREG) <- dds_DWNREG$MvA_Primary_History # assign the rownames as the conditions of your samples
pheatmap(sampleDistMatrix.rlog.dds_DWNREG, 
         clustering_distance_rows = sampleDist.rlog.dds_DWNREG, 
         clustering_distance_cols = sampleDist.rlog.dds_DWNREG, 
         col = colors) # create the heatmap 
plotPCA(rlog.dds_DWNREG,intgroup="MvA_Primary_History")  # use the DESeq2 call plotPCA() and call the rlog_dds for PCA plotting - Q: difference between replicates or condition explains the most variance in this dataset?





# AMBIENT RHOGUHTOU ONLY - LOOK FOR TIME-DEPENDENT DEGs UNDER AMBIENT TREATMENT ====================================== #
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
table(resAmb.d0vd21$padj<0.05) # 90 DEGs


# Write results
# path_out.d0 = 'C:/Users/samjg/Documents/My_Projects/Pgenerosa_TagSeq_Metabolomics/TagSeq/RAnalysis/DESeq2/output/Day0/'
# write.csv(resdata.d0.primary, paste(path_out.d0, "Day0.PrimaryTreatment_diffexpr-results.csv"))

# volcano plot ------------------------------------------------------------------------------------------------------ #

# png("RAnalysis/DESeq2/output/Day0/Day0.PrimaryTreatment-VolcanoPlot.png", 1000, 1000, pointsize=20)
EnhancedVolcano(resAmb.d0vd21,
                lab = rownames(resAmb.d0vd21),
                x = 'log2FoldChange',
                y = 'padj',
                title = 'Day 21 v. Day 0 (Ambient Control)',
                subtitle = "DESeq2 - Differential expression",
                FCcutoff = 2.0,
                pCutoff = 0.3,
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
nrow(dds.d0) # 15149 total genes pre filtered
dim(dds.d0) # 15149     8
design(dds.d0) # FULL MODEL ==- ~Primary_Treatment

# RUN DESEQ2 model
dds.d0 <- DESeq(dds.d0) # wait for this to complete....
resultsNames(dds.d0) # view the names of your results model 'Primary_Treatment_M_vs_A'
colData(dds.d0)

# Primary_Treatment_M_vs_A
resd0.primary <- results(dds.d0, alpha = 0.5) # already filtered in edgeR; FDR of 0.05 (~ 651.8 in 15.036 genes)
resd0.primary <- results(dds.d0, name="Primary_Treatment_M_vs_A", alpha = 0.05)

hist(resd0.primary$pvalue, breaks=20, col="grey") # view histogram,  about 425 genes are likely false positives, alpha = 0.05 FDR is good!
abline(h=(nrow(dds.d0)*0.05),col="red") # add line at expected 5% false positive
table(resd0.primary$padj<0.05) #  13  DEGs
resd0.primary <- resd0.primary[order(resd0.primary$padj), ] ## Order by adjusted p-value

# count and % upregulated and downregulated DEGs (first 13 rows of ordered table are the DEGs)
sum((resd0.primary$log2FoldChange[1:11] >= 1) == TRUE) # 3 DEGs upregulated  (LFC >= 1)
( sum((resd0.primary$log2FoldChange[1:11] >= 1) == TRUE) / (table(resd0.primary$padj<0.05))[2] ) * 100 # 27.27273 % DEGs upregulated (N = 971)
sum((resd0.primary$log2FoldChange[1:11] <= -1) == TRUE) # 8 DEGs downregulated  (LFC <= -1)
( sum((resd0.primary$log2FoldChange[1:11] <= -1) == TRUE) / (table(resd0.primary$padj<0.05))[2] ) * 100 # 72.72727 % DEGs downregulated (N = 971)

# Write results - covert to as.data.frame for the ordered results
resdata.d0.primary  <- merge(as.data.frame(resd0.primary), as.data.frame(counts(dds.d0, normalized=TRUE)), by="row.names", sort=FALSE) ## Merge with normalized count data
names(resdata.d0.primary)[1] <- "Gene"
path_out.d0 = 'C:/Users/samjg/Documents/My_Projects/Pgenerosa_TagSeq_Metabolomics/TagSeq/RAnalysis/DESeq2/output/Day0/' # call path
write.csv(resdata.d0.primary, paste(path_out.d0, "Day0.MvA_DESeq2results.csv")) # write

# volcano plot 
png("RAnalysis/DESeq2/output/Day0/Day0.PrimaryTreatment-VolcanoPlot.png", 1000, 1000, pointsize=20)
EnhancedVolcano(resd0.primary,
                lab = rownames(resd0.primary),
                x = 'log2FoldChange',
                y = 'padj',
                title = 'Day 0: Primary Treatment (M v A)',
                subtitle = "DESeq2 - Differential expression",
                FCcutoff = 1,
                pCutoff = 0.05,
                pointSize = 4.0,
                labSize = 6.0,
                colAlpha = 1,
                legendPosition = 'right',
                legendLabSize = 12,
                legendIconSize = 4.0,
                widthConnectors = 0.75)
dev.off()

# Plot dispersions 
png("RAnalysis/DESeq2/output/Day0/Day0.dispersions.png", 1000, 1000, pointsize=20)
plotDispEsts(dds.d0, main="Day 0 dispersions")
dev.off()

# Data transformations for heatmap and PCA visuals ======================================================================================= #
# rlog - regularized log transformation of origin count data to log2 scale - fit for each sample and dist. of coefficients in the data
rlog.d0<- rlogTransformation(dds.d0) # rlog transform (regularized log)

png("RAnalysis/DESeq2/output/Day0/Day0.rlog_histogram.png", 1000, 1000, pointsize=20)# diagnostics of transformation # Histogram and sd plot
hist(assay(rlog.d0)) # view histogram 
dev.off()
png("RAnalysis/DESeq2/output/Day0/Day0.rlog_mean_sd.png", 1000, 1000, pointsize=20)
meanSdPlot(assay(rlog.d0)) # shows the sd y axis (sq root of varaince in all samples) - flat curve may seem like a goals, BUT may be unreasonable in cases with MANY true DEGs from experimental conditions
dev.off()

# PCA plot rlog ------------------------ #
library("ggplot2")
pcaData_d0 <- plotPCA(rlog.d0, intgroup = "Primary_Treatment", returnData = TRUE)
percentVar <- round(100 * attr(pcaData_d0, "percentVar"))
png("RAnalysis/DESeq2/output/Day0/Day0.rlog_PCA.png", 1000, 1000, pointsize=20)
ggplot(pcaData_d0, aes(x = PC1, y = PC2, color = Primary_Treatment, label=name)) +
  scale_shape_manual(values = c(4, 19, 17)) +
  geom_text(aes(label=name),hjust=0.2, vjust=1.4, size=5) +
  geom_point(size =6) +
  theme_classic() +
  theme(text = element_text(size=15)) +
  ggtitle("PCA: Day0 (rlog)") +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  coord_fixed()
dev.off()
# Plot heat map rlog------------------------ #
save_pheatmap <- function(x, filename, width=1000, height=960) { # template for saving pheatmap outputs
  stopifnot(!missing(x))
  stopifnot(!missing(filename))
  png(filename,width = width, height=height)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}
select <- order(rowMeans(counts(dds.d0,normalized=TRUE)), decreasing=TRUE)[1:25] # normalize the counts per row and call the first 25 most expressed genes
df <- as.data.frame(colData(dds.d0)[c("Primary_Treatment")])
annotation_colors = list(Treatment = c(A="Blue", M="Orange"))
d0.rlog.heatmap<- pheatmap(assay(rlog.d0)[select,], cluster_rows=TRUE, show_rownames=TRUE, # rlog heatmap
         cluster_cols=TRUE, annotation_col=df, annotation_colors = annotation_colors,
         main = "Day0.rlog_heatmap")
save_pheatmap(d0.rlog.heatmap, filename = "RAnalysis/DESeq2/output/Day0/Day0.rlog_heatmap.png")

# ========================================================== 
#
# DAY 7    (3 CPM in 50% samples)
# ========================================================== 
# view the design for each of the d7 dds 
design(dds.d7) # ~Second_Treament + Primary_Treatment + Second_Treament:Primary_Treatment
design(dds.d7.group) # ~All_Treatment  (group)
design(dds.d7.main) # ~Second_Treament + Primary_Treatment

# RUN DESEQ2 model - view all the pariwise comparisons
dds.d7 <- DESeq(dds.d7) #  full model                                   wait for this to complete....
dds.d7.group <- DESeq(dds.d7.group) # group model                       wait for this to complete....
dds.d7.main <- DESeq(dds.d7.main) # main effect model                   wait for this to complete....

# model matrix 
# View(model.matrix(~ Second_Treament + Primary_Treatment+Second_Treament*Primary_Treatment, colData(dds.d7)))
# View(model.matrix(~ All_Treatment, colData(dds.d7.group)))
# View(model.matrix(~ Second_Treament + Primary_Treatment, colData(dds.d7.main)))

# ~ FULL MODEL EFFECTS  ======================================================================================================================== #
resultsNames(dds.d7) # view the names of your results model 
dds.d7$Primary_Treatment # ambient is the reference 
dds.d7$Second_Treament # ambient is the reference 

# Q: What is the main effect (Primary treatment) during SECONDARY response to AMBINET? (MA vs AA)
resd7._MAvAA <- results(dds.d7, contrast=c("Primary_Treatment", "M", "A"), alpha = 0.05)
table(resd7._MAvAA$padj<0.05) # 2 DEGs

# Q: What is the main effect (Primary treatment) during SECONDARY response to MODERATE? (MM vs AM)
# This is, by definition, the main effect (Primary treatment) plus the interaction term 
# (the extra condition effect in secondary treatment  'Severe' compared to Secondary Treatment 'Ambient').
resd7._MMvAM <- results(dds.d7, list(c("Primary_Treatment_M_vs_A","Second_TreamentM.Primary_TreatmentM")))
table(resd7._MMvAM$padj<0.05) # 4  DEGs

# Q: What is the main effect (Primary treatment) during SECONDARY response to SEVERE? (MS vs AS)
# This is, by definition, the main effect (Primary treatment) plus the interaction term 
# (the extra condition effect in secondary treatment  'Severe' compared to Secondary Treatment 'Ambient').
resd7._MSvAS <- results(dds.d7, list(c("Primary_Treatment_M_vs_A","Second_TreamentS.Primary_TreatmentM")))
table(resd7._MSvAS$padj<0.05) # 1  DEGs

#Effect of Severe vs. Moderate for the Primary Ambient individuals (AM v AS)
resd7._AMvAS <- results(dds.d7, contrast= c(0,-1,1,0,0,0))
table(resd7._AMvAS$padj<0.05) # 0 DEGs

# Interaction term for the main  effect (Primary Treatment) in Secondary Moderate versus Ambient
resd7._MMvAM_MAvAA<- results(dds.d7, name="Second_TreamentM.Primary_TreatmentM", alpha = 0.05) # (MMvAM × MAvAA)
table(resd7._MMvAM_MAvAA$padj<0.05) # 0 DEGs

# Interaction term for the main  effect (Primary Treatment) in Secondary SEVERE versus Ambient
resd7._MSvAS_MAvAA<- results(dds.d7, name="Second_TreamentS.Primary_TreatmentM", alpha = 0.05) # (MSvAS × MAvAA)
table(resd7._MSvAS_MAvAA$padj<0.05) # 0 DEGs

# Interaction term for the main  effect (Primary Treatment) in Secondary SEVERE  versus MODERATE
resd7._MSvAS_MMvAM= results(dds.d7, contrast=list("Second_TreamentS.Primary_TreatmentM", "Second_TreamentM.Primary_TreatmentM"), alpha = 0.05) # (MSvAS × MMvAM)
table(resd7._MSvAS_MMvAM$padj<0.05) # 0 DEGs


# ~ Main EFFECTS  ======================================================================================================================== #
resultsNames(dds.d7.main) # view the names of your results model 

# Main effect of Primary treatment (WITHOUT ref level in other factor - no interaction term in this model!!!)
resd7._P.MvA<- results(dds.d7.main, name="Primary_Treatment_M_vs_A", alpha = 0.05) 
table(resd7._P.MvA$padj<0.05) # 99 DEGs
resd7._P.MvA <- resd7._P.MvA[order(resd7._P.MvA$padj), ] ## Order by adjusted p-value
sum((resd7._P.MvA$log2FoldChange[1:99] >= 1) == TRUE) # 32 DEGs upregulated  (LFC >= 1)
( sum((resd7._P.MvA$log2FoldChange[1:99] >= 1) == TRUE) / (table(resd7._P.MvA$padj<0.05))[2] ) * 100 # 32.32323 % DEGs upregulated (N = 971)
sum((resd7._P.MvA$log2FoldChange[1:99] <= -1) == TRUE) # 44 DEGs downregulated  (LFC <= -1)
( sum((resd7._P.MvA$log2FoldChange[1:99] <= -1) == TRUE) / (table(resd7._P.MvA$padj<0.05))[2] ) * 100 # 44.44444 % DEGs downregulated (N = 971)

# What is the difference between secondary treatment 'moderate' WITHOUT considering the primary treament history?
resd7._S.MvA<- results(dds.d7.main, name="Second_Treament_M_vs_A", alpha = 0.05)
table(resd7._S.MvA$padj<0.05) # 118 DEGs
resd7._S.MvA <- resd7._S.MvA[order(resd7._S.MvA$padj), ] ## Order by adjusted p-value
sum((resd7._S.MvA$log2FoldChange[1:118] >= 1) == TRUE) # 91 DEGs upregulated  (LFC >= 1)
( sum((resd7._S.MvA$log2FoldChange[1:118] >= 1) == TRUE) / (table(resd7._S.MvA$padj<0.05))[2] ) * 100 # 77.11864 % DEGs upregulated (N = 971)
sum((resd7._S.MvA$log2FoldChange[1:118] <= -1) == TRUE) #  6 DEGs downregulated  (LFC <= -1)
( sum((resd7._S.MvA$log2FoldChange[1:118] <= -1) == TRUE) / (table(resd7._S.MvA$padj<0.05))[2] ) * 100 # 5.084746 % DEGs downregulated (N = 971)

# What is the difference between secondary treatment 'severe' WITHOUT considering the primary treament history?
resd7._S.SvA<- results(dds.d7.main, name="Second_Treament_S_vs_A", alpha = 0.05)
table(resd7._S.SvA$padj<0.05) # 1 DEGs
resd7._S.SvA <- resd7._S.SvA[order(resd7._S.SvA$padj), ] ## Order by adjusted p-value
sum((resd7._S.SvA$log2FoldChange[1] >= 1) == TRUE) # 1 DEGs upregulated  (LFC >= 1)
sum((resd7._S.SvA$log2FoldChange[1] <= -1) == TRUE) #  0 DEGs downregulated  (LFC <= -1)

# What is the difference between secondary treatment 'severe' and second treatment 'moderate' WITHOUT considering the primary treament history?
resd7._S.MvS <- results(dds.d7.main, contrast = list(("Second_Treament_M_vs_A"),("Second_Treament_S_vs_A")),  alpha= 0.05) # Second Severe
table(resd7._S.MvS$padj<0.05) # 14 DEGs
resd7._S.MvS <- resd7._S.MvS[order(resd7._S.MvS$padj), ] ## Order by adjusted p-value
sum((resd7._S.MvS$log2FoldChange[1] >= 1) == TRUE) # 1 DEGs upregulated  (LFC >= 1)
sum((resd7._S.MvS$log2FoldChange[1] <= -1) == TRUE) #  0 DEGs downregulated  (LFC <= -1)


# ~ Group - MAIN EFFECTS(ACCOUNTING FOR SECONDARY-TREATMENT SPECIFIC EFFECTS) ============================================================================================ #
resultsNames(dds.d7.group) # view the names of your results model 

# ALL primary history regarless of subsequent exposure (Ambient versus Moderate)
res.d7_P.MvA_group <- results(dds.d7.group, contrast = list(c("All_TreatmentAA", "All_TreatmentAM", "All_TreatmentAS"),c("All_TreatmentMA", "All_TreatmentMM", "All_TreatmentMS")),  alpha= 0.05)
table(res.d7_P.MvA_group$padj<0.05) # 111 DEGs
res.d7_P.MvA_group <- res.d7_P.MvA_group[order(res.d7_P.MvA_group$padj), ] ## Order by adjusted p-value
sum((res.d7_P.MvA_group$log2FoldChange[1:111] >= 1) == TRUE) # 66 DEGs upregulated  (LFC >= 1)
( sum((res.d7_P.MvA_group$log2FoldChange[1:111] >= 1) == TRUE) / (table(res.d7_P.MvA_group$padj<0.05))[2] ) * 100 # 59.45946 % DEGs upregulated (N = 971)
sum((res.d7_P.MvA_group$log2FoldChange[1:111] <= -1) == TRUE) #  44 DEGs downregulated  (LFC <= -1)
( sum((res.d7_P.MvA_group$log2FoldChange[1:111] <= -1) == TRUE) / (table(res.d7_P.MvA_group$padj<0.05))[2] ) * 100 # 39.63964 % DEGs downregulated (N = 971)

# CONSTANT exposure to SAME treatment
resd7._AAvMM_group <- results(dds.d7.group, contrast = list(("All_TreatmentAA"),("All_TreatmentMM")),  alpha= 0.05) # Constant expore (basically M v A Primaru + 7 more days)
table(resd7._AAvMM_group$padj<0.05) # 28 DEGs
# SECOND exposure by PRIMARY treatment
resd7._AAvMA_group <- results(dds.d7.group, contrast = list(("All_TreatmentAA"),("All_TreatmentMA")),  alpha= 0.05) # Second Ambient
table(resd7._AAvMA_group$padj<0.05) # 2 DEGs
resd7._AMvMM_group <- results(dds.d7.group, contrast = list(("All_TreatmentAM"),("All_TreatmentMM")),  alpha= 0.05) # Second Moderate
table(resd7._AMvMM_group$padj<0.05) # 3 DEGs
resd7._ASvMS_group <- results(dds.d7.group, contrast = list(("All_TreatmentAS"),("All_TreatmentMS")),  alpha= 0.05) # Second Severe
table(resd7._ASvMS_group$padj<0.05) # 0 DEGs
# AMBIENT history and subseqent treatment 
resd7._AAvAM_group <- results(dds.d7.group, contrast = list(("All_TreatmentAA"),("All_TreatmentAM")),  alpha= 0.05) # Amb - response to Moderate
table(resd7._AAvAM_group$padj<0.05) # 1 DEGs
resd7._AAvAS_group <- results(dds.d7.group, contrast = list(("All_TreatmentAA"),("All_TreatmentAS")),  alpha= 0.05) # Amb - response to Severe
table(resd7._AAvAS_group$padj<0.05) # 6 DEGs
resd7._AMvAS_group <- results(dds.d7.group, contrast = list(("All_TreatmentAM"),("All_TreatmentAS")),  alpha= 0.05) # Amb - response to moderate versus Severe
table(resd7._AMvAS_group$padj<0.05) # 0 DEGs
# MODERATE history and subseqent treatment 
resd7._MMvMA_group <- results(dds.d7.group, contrast = list(("All_TreatmentMM"),("All_TreatmentMA")),  alpha= 0.05) # Mod - response to Ambient
table(resd7._MMvMA_group$padj<0.05) # 35 DEGs
resd7._MMvMS_group <- results(dds.d7.group, contrast = list(("All_TreatmentMM"),("All_TreatmentMS")),  alpha= 0.05) # Mod - response to Severe
table(resd7._MMvMS_group$padj<0.05) # 0 DEGs
resd7._MAvMS_group <- results(dds.d7.group, contrast = list(("All_TreatmentMA"),("All_TreatmentMS")),  alpha= 0.05) # Amb - response to ambient versus Severe
table(resd7._MAvMS_group$padj<0.05) # 0 DEGs
# OTHER pairwise interactions
#AA
resd7._AAvMM_group <- results(dds.d7.group, contrast = list(("All_TreatmentAA"),("All_TreatmentMM")),  alpha= 0.05) # Constant expore (basically M v A Primaru + 7 more days)
table(resd7._AAvMM_group$padj<0.05) # 28 DEGs
resd7._AAvMA_group <- results(dds.d7.group, contrast = list(("All_TreatmentAA"),("All_TreatmentMA")),  alpha= 0.05)
table(resd7._AAvMA_group$padj<0.05) # 2 DEGs
resd7._AAvMS_group <- results(dds.d7.group, contrast = list(("All_TreatmentAA"),("All_TreatmentMS")),  alpha= 0.05)
table(resd7._AAvMS_group$padj<0.05) # 8 DEGs
#AM
resd7._AMvMA_group <- results(dds.d7.group, contrast = list(("All_TreatmentAM"),("All_TreatmentMA")),  alpha= 0.05)
table(resd7._AMvMA_group$padj<0.05) # 187 DEGs
resd7._AMvMS_group <- results(dds.d7.group, contrast = list(("All_TreatmentAM"),("All_TreatmentMS")),  alpha= 0.05)
table(resd7._AMvMS_group$padj<0.05) # 51 DEGs
#AS
resd7._ASvMA_group <- results(dds.d7.group, contrast = list(("All_TreatmentAS"),("All_TreatmentMA")),  alpha= 0.05)
table(resd7._ASvMA_group$padj<0.05) # 22 DEGs
resd7._ASvMM_group <- results(dds.d7.group, contrast = list(("All_TreatmentAS"),("All_TreatmentMM")),  alpha= 0.05)
table(resd7._ASvMM_group$padj<0.05) # 25 DEGs
# Second treatment effects (regarless of prior history)
#A v M
resd7._S.AvM_group <- results(dds.d7.group, contrast = list(c("All_TreatmentAM", "All_TreatmentMM"), c("All_TreatmentAA", "All_TreatmentMA")),  alpha= 0.05)
table(resd7._S.AvM_group$padj<0.05) # 114 DEGs
resd7._S.AvM_group <- resd7._S.AvM_group[order(resd7._S.AvM_group$padj), ] ## Order by adjusted p-value
sum((resd7._S.AvM_group$log2FoldChange[1:114] >= 1) == TRUE) # 8 DEGs upregulated  (LFC >= 1)
( sum((resd7._S.AvM_group$log2FoldChange[1:114] >= 1) == TRUE) / (table(resd7._S.AvM_group$padj<0.05))[2] ) * 100 # 7.017544 % DEGs upregulated (N = 971)
sum((resd7._S.AvM_group$log2FoldChange[1:114] <= -1) == TRUE) #  104 DEGs downregulated  (LFC <= -1)
( sum((resd7._S.AvM_group$log2FoldChange[1:114] <= -1) == TRUE) / (table(resd7._S.AvM_group$padj<0.05))[2] ) * 100 # 91.22807 % DEGs downregulated (N = 971)

# A v S
resd7._S.AvS_group <- results(dds.d7.group, contrast = list(c("All_TreatmentAA", "All_TreatmentMA"),c("All_TreatmentAS", "All_TreatmentMS")),  alpha= 0.05)
table(resd7._S.AvS_group$padj<0.05) # 1 DEGs
resd7._S.AvS_group <- resd7._S.AvS_group[order(resd7._S.AvS_group$padj), ] ## Order by adjusted p-value
sum((resd7._S.AvS_group$log2FoldChange[1] >= 1) == TRUE) # 0
sum((resd7._S.AvS_group$log2FoldChange[1] <= -1) == TRUE) # 1

#M v S
resd7._S.MvS_group <- results(dds.d7.group, contrast = list(c("All_TreatmentAM", "All_TreatmentMM"),c("All_TreatmentAS", "All_TreatmentMS")),  alpha= 0.05)
table(resd7._S.MvS_group$padj<0.05) # 13 DEGs
resd7._S.MvS_group <- resd7._S.MvS_group[order(resd7._S.MvS_group$padj), ] ## Order by adjusted p-value
sum((resd7._S.MvS_group$log2FoldChange[1:13] >= 1) == TRUE) # 10 DEGs upregulated  (LFC >= 1)
( sum((resd7._S.MvS_group$log2FoldChange[1:13] >= 1) == TRUE) / (table(resd7._S.MvS_group$padj<0.05))[2] ) * 100 # 76.92308 % DEGs upregulated (N = 971)
sum((resd7._S.MvS_group$log2FoldChange[1:13] <= -1) == TRUE) #  3 DEGs downregulated  (LFC <= -1)
( sum((resd7._S.MvS_group$log2FoldChange[1:13] <= -1) == TRUE) / (table(resd7._S.MvS_group$padj<0.05))[2] ) * 100 # 23.07692 % DEGs downregulated (N = 971)




## Write results
# path_out.d7 = 'C:/Users/samjg/Documents/My_Projects/Pgenerosa_TagSeq_Metabolomics/TagSeq/RAnalysis/DESeq2/output/Day7/'
# write.csv(resdata.d7.all.prim.M_vs_A, paste(path_out.d7, "Day7.PrimaryTreament_diffexpr-results.csv"))
# write.csv(resdata.d7.second_M_vs_A, paste(path_out.d7, "Day7.SecondTreament_MvsA_diffexpr-results.csv"))
# write.csv(resdata.d7.all.second_S_vs_A, paste(path_out.d7, "Day7.SecondTreament_SvsA_diffexpr-results.csv"))


# volcano plot ------------------------------------------------------------------------------------------------------ #
# Effect of primary treatment 'res.d7_P.MvA_group' - 111 DEGs
png("RAnalysis/DESeq2/output/Day7/Day7.Primary_VolcanoPlot.png", 1000, 1000, pointsize=20)
EnhancedVolcano(res.d7_P.MvA_group,
                lab = rownames(res.d7_P.MvA_group),
                x = 'log2FoldChange',
                y = 'padj',
                title = 'Day7 Primary treatment (M v A)',
                subtitle = "DESeq2 - Differential expression",
                FCcutoff = 1,
                pCutoff = 0.05,
                pointSize = 4.0,
                labSize = 6.0,
                colAlpha = 1,
                legendPosition = 'right',
                legendLabSize = 12,
                legendIconSize = 4.0,
                widthConnectors = 0.75)
dev.off()

# Effect of second treatment M v. A 'resd7._S.AvM_group' - 114 DEGs
png("RAnalysis/DESeq2/output/Day7/Day7.SecondAvM_VolcanoPlot.png", 1000, 1000, pointsize=20)
EnhancedVolcano(resd7._S.AvM_group,
                lab = rownames(resd7._S.AvM_group),
                x = 'log2FoldChange',
                y = 'padj',
                title = 'Day7 Second treatment (A v M)',
                subtitle = "DESeq2 - Differential expression",
                FCcutoff = 1,
                pCutoff = 0.05,
                pointSize = 4.0,
                labSize = 6.0,
                colAlpha = 1,
                legendPosition = 'right',
                legendLabSize = 12,
                legendIconSize = 4.0,
                widthConnectors = 0.75)
dev.off()

# Effect of second treatment S v. A 'resd7._S.AvS_group' - 1 DEG
png("RAnalysis/DESeq2/output/Day7/Day7.SecondAvS_VolcanoPlot.png", 1000, 1000, pointsize=20)
EnhancedVolcano(resd7._S.AvS_group,
                lab = rownames(resd7._S.AvS_group),
                x = 'log2FoldChange',
                y = 'padj',
                title = 'Day7 Second treatment (A v S)',
                subtitle = "DESeq2 - Differential expression",
                FCcutoff = 1.0,
                pCutoff = 0.05,
                pointSize = 4.0,
                labSize = 6.0,
                colAlpha = 1,
                legendPosition = 'right',
                legendLabSize = 12,
                legendIconSize = 4.0,
                widthConnectors = 0.75)
dev.off()

# Effect of second treatment M v. S 'resd7._S.MvS_group' - 1 DEG
png("RAnalysis/DESeq2/output/Day7/Day7.SecondMvs_VolcanoPlot.png", 1000, 1000, pointsize=20)
EnhancedVolcano(resd7._S.MvS_group,
                lab = rownames(resd7._S.MvS_group),
                x = 'log2FoldChange',
                y = 'padj',
                title = 'Day7 Second treatment (M v S)',
                subtitle = "DESeq2 - Differential expression",
                FCcutoff = 1.0,
                pCutoff = 0.05,
                pointSize = 4.0,
                labSize = 6.0,
                colAlpha = 1,
                legendPosition = 'right',
                legendLabSize = 12,
                legendIconSize = 4.0,
                widthConnectors = 0.75)
dev.off()

# Effect of second treatment S v. A 'resd7.all.second_S_vs_A' - 1 DEG
png("RAnalysis/DESeq2/output/Day7/Day7.AMvMA_VolcanoPlot.png", 1000, 1000, pointsize=20)
EnhancedVolcano(resd7._AMvMA_group,
                lab = rownames(resd7._AMvMA_group),
                x = 'log2FoldChange',
                y = 'padj',
                title = 'Day 7 Prim.Second (AM v MA)',
                subtitle = "DESeq2 - Differential expression",
                FCcutoff = 2.0,
                pCutoff = 0.05,
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

# Heat maps and Principal components analysis ------   Primary_Treatment_M_vs_A ----------------------------------------- #
# significantly upregulated genes - padj < 0.05; LFC > 2
res7.MvA_UPREG <- res.d7_P.MvA_group[c(1:111),] 
d7.upreg.T_F <- res7.MvA_UPREG$log2FoldChange >= 1 # we would like to keep genes that have at least 50% TRUES in each row of thresh
res7.MvA_UPREG <- res7.MvA_UPREG[d7.upreg.T_F,]  # call onl upreg genes
dds.d7_UPREG <- dds.d7[(rownames(res7.MvA_UPREG))] # call only res.ALL.MvA_UPREG rows (12) 
dim(dds.d7_UPREG) # 66 36 - 66 pass this criteria as upregulated
rlog.dds.d7_UPREG<- rlogTransformation(dds.d7_UPREG) # Conmplete rlog transformation ....takes a bit of time
sampleDist.rlog.dds.d7_UPREG <- dist(t(assay(rlog.dds.d7_UPREG))) # next assay [access dds data matrix] and t() [transpose rows and columns] in dist() [calculate the Euclidean Distance]
sampleDistMatrix.rlog.dds.d7_UPREG <- as.matrix(sampleDist.rlog.dds.d7_UPREG) # convert to a data matrix for heatmap 
rownames(sampleDistMatrix.rlog.dds.d7_UPREG) <- dds.d7_UPREG$Primary_Treatment # assign the rownames as the conditions of your samples
pheatmap(sampleDistMatrix.rlog.dds.d7_UPREG, 
         clustering_distance_rows = sampleDist.rlog.dds.d7_UPREG, 
         clustering_distance_cols = sampleDist.rlog.dds.d7_UPREG) # create the heatmap 
plotPCA(rlog.dds.d7_UPREG,intgroup="Primary_Treatment") # use the DESeq2 call plotPCA() and call the rlog_dds for PCA plotting - Q: difference between replicates or condition explains the most variance in this dataset?


select <- order(rowMeans(counts(dds.d7_UPREG,normalized=TRUE)), decreasing=TRUE)[1:25] # normalize the counts per row and call the first 20 most expressed genes
df <- as.data.frame(colData(dds.d7_UPREG)[c("Primary_Treatment", "Second_Treament")])
annotation_colors = list((Primary_Treatment = c(A="Blue", M="Orange")), Second_Treament = c(A="White", M="Grey", S="Black"))
d7.UPREG.rlog.heatmap<- pheatmap(assay(rlog.dds.d7_UPREG)[select,], cluster_rows=TRUE, show_rownames=TRUE, # rlog heatmap
                           cluster_cols=TRUE, annotation_col=df, annotation_colors = annotation_colors,
                           main = "Day7.UPREG.rlog_heatmap")
save_pheatmap(d7.UPREG.rlog.heatmap, filename = "RAnalysis/DESeq2/output/Day7/Day7.PrimaryMvA_Upregulated_rlog_heatmap.png")


# significantly downregulated genes - padj < 0.05; LFC < 2
res7.MvA_DWNREG <- res.d7_P.MvA_group[c(1:111),] 
d7.dwnreg.T_F <- res7.MvA_DWNREG$log2FoldChange <= -1 # we would like to keep genes that have at least 50% TRUES in each row of thresh
res7.MvA_DWNREG <- res7.MvA_DWNREG[d7.dwnreg.T_F,] # call only downreg genes
dds.d7_DWNREG <- dds.d7[(rownames(res7.MvA_DWNREG))] # call only res.ALL.MvA_UPREG rows (12) 
dim(dds.d7_DWNREG) # 41 36 - 41 pass this criteria
rlog.dds.d7_DWNREG<- rlogTransformation(dds.d7_DWNREG) # Conmplete rlog transformation  ....takes a bit of time
sampleDist.rlog.dds_DWNREG <- dist(t(assay(rlog.dds.d7_DWNREG))) # next assay [access dds data matrix] and t() [transpose rows and columns] in dist() [calculate the Euclidean Distance]
sampleDistMatrix.rlog.dds.d7_DWNREG <- as.matrix(sampleDist.rlog.dds_DWNREG) # convert to a data matrix for heatmap 
rownames(sampleDistMatrix.rlog.dds.d7_DWNREG) <- rlog.dds.d7_DWNREG$Primary_Treatment # assign the rownames as the conditions of your samples
pheatmap(sampleDistMatrix.rlog.dds.d7_DWNREG, 
         clustering_distance_rows = sampleDist.rlog.dds_DWNREG, 
         clustering_distance_cols = sampleDist.rlog.dds_DWNREG) # create the heatmap 
plotPCA(rlog.dds.d7_DWNREG,intgroup="Primary_Treatment")  # use the DESeq2 call plotPCA() and call the rlog_dds for PCA plotting - Q: difference between replicates or condition explains the most variance in this dataset?


select <- order(rowMeans(counts(dds.d7_DWNREG,normalized=TRUE)), decreasing=TRUE)[1:25] # normalize the counts per row and call the first 20 most expressed genes
df <- as.data.frame(colData(dds.d7_DWNREG)[c("Primary_Treatment", "Second_Treament")])
annotation_colors = list((Primary_Treatment = c(A="Blue", M="Orange")), Second_Treament = c(A="White", M="Grey", S="Black"))
d7.DOWNREG.rlog.heatmap<- pheatmap(assay(rlog.dds.d7_DWNREG)[select,], cluster_rows=TRUE, show_rownames=TRUE, # rlog heatmap
                                 cluster_cols=TRUE, annotation_col=df, annotation_colors = annotation_colors,
                                 main = "Day7.DOWNREG.rlog_heatmap")
save_pheatmap(d7.DOWNREG.rlog.heatmap, filename = "RAnalysis/DESeq2/output/Day7/Day7.PrimaryMvA_downregulated_rlog_heatmap.png")


# ======================================== 
# Data transformations for heatmap and PCA visuals 
# ============================================================= 
rlog.d7<- rlogTransformation(dds.d7) #  full model                       wait for this to complete.... 

png("RAnalysis/DESeq2/output/Day7/Day7.rlog_histogram.png", 1000, 1000, pointsize=20)# diagnostics of transformation # Histogram and sd plot
hist(assay(rlog.d7)) # view histogram 
dev.off()
png("RAnalysis/DESeq2/output/Day7/Day7.rlog_mean_sd.png", 1000, 1000, pointsize=20)
meanSdPlot(assay(rlog.d7)) # shows the sd y axis (sq root of varaince in all samples) - flat curve may seem like a goals, BUT may be unreasonable in cases with MANY true DEGs from experimental conditions
dev.off()

# PCA plot rlog ------------------------ #
pcaData_d7 <- plotPCA(rlog.d7, intgroup = c( "Primary_Treatment", "Second_Treament"), returnData = TRUE)
percentVar <- round(100 * attr(pcaData_d7, "percentVar"))
png("RAnalysis/DESeq2/output/Day7/Day7.rlog_PCA.png", 1000, 1000, pointsize=20)
ggplot(pcaData_d7, aes(x = PC1, y = PC2, color = Primary_Treatment, shape = Second_Treament)) +
  scale_shape_manual(values = c(4, 19, 17)) +
  geom_text(aes(label=name),hjust=0.2, vjust=1.4, size=3) +
  geom_point(size =6) +
  theme_classic() +
  theme(text = element_text(size=15)) +
  ggtitle("PCA: Day7 (rlog)") +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  coord_fixed()
dev.off()

# Plot heat map rlog------------------------ #
select <- order(rowMeans(counts(dds.d7,normalized=TRUE)), decreasing=TRUE)[1:20] # normalize the counts per row and call the first 20 most expressed genes
df <- as.data.frame(colData(dds.d7)[c("Primary_Treatment", "Second_Treament")])
annotation_colors = list((Primary_Treatment = c(A="Blue", M="Orange")), Second_Treament = c(A="White", M="Grey", S="Black"))
d7.rlog.heatmap<- pheatmap(assay(rlog.d7)[select,], cluster_rows=TRUE, show_rownames=TRUE, # rlog heatmap
                           cluster_cols=TRUE, annotation_col=df, annotation_colors = annotation_colors,
                           main = "Day7.rlog_heatmap")
save_pheatmap(d7.rlog.heatmap, filename = "RAnalysis/DESeq2/output/Day7/Day7.rlog_heatmap.png")

# ========================================================== 
#
# DAY 14   (3 CPM in 50% samples)
# ========================================================== 

# view the design for each of the d7 dds 
design(dds.d14) # ~Second_Treament + Primary_Treatment + Second_Treament:Primary_Treatment
design(dds.d14.group) # ~All_Treatment - 1 (group)
design(dds.d14.main) # ~Second_Treament + Primary_Treatment

# RUN DESEQ2 model - view all the pariwise comparisons
dds.d14 <- DESeq(dds.d14) #  full model                                   wait for this to complete....
dds.d14.group <- DESeq(dds.d14.group) # group model                       wait for this to complete....
dds.d14.main <- DESeq(dds.d14.main) # main effect model                   wait for this to complete....

# Q: What is the main effect (Primary treatment) during SECONDARY response to AMBINET? (MA vs AA)
resd14._MAvAA <- results(dds.d14, contrast=c("Primary_Treatment", "M", "A"), alpha = 0.05)
table(resd14._MAvAA$padj<0.05) # 4 DEGs

# Q: What is the main effect (Primary treatment) during SECONDARY response to MODERATE? (MM vs AM)
# This is, by definition, the main effect (Primary treatment) plus the interaction term 
# (the extra condition effect in secondary treatment  'Severe' compared to Secondary Treatment 'Ambient').
resd14._MMvAM <- results(dds.d14, list(c("Primary_Treatment_M_vs_A","Second_TreamentM.Primary_TreatmentM")))
table(resd14._MMvAM$padj<0.05) # 7  DEGs

# Q: What is the main effect (Primary treatment) during SECONDARY response to SEVERE? (MS vs AS)
# This is, by definition, the main effect (Primary treatment) plus the interaction term 
# (the extra condition effect in secondary treatment  'Severe' compared to Secondary Treatment 'Ambient').
resd14._MSvAS <- results(dds.d14, list(c("Primary_Treatment_M_vs_A","Second_TreamentS.Primary_TreatmentM")))
table(resd14._MSvAS$padj<0.05) # 2  DEGs

#Effect of Severe vs. Moderate for the Primary Ambient individuals (AM v AS)
resd14._AMvAS <- results(dds.d14, contrast= c(0,-1,1,0,0,0))
table(resd14._AMvAS$padj<0.05) # 0 DEGs

# Interaction term for the main  effect (Primary Treatment) in Secondary Moderate versus Ambient
resd14._MMvAM_MAvAA<- results(dds.d14, name="Second_TreamentM.Primary_TreatmentM", alpha = 0.05) # (MMvAM × MAvAA)
table(resd14._MMvAM_MAvAA$padj<0.05) # 0 DEGs

# Interaction term for the main  effect (Primary Treatment) in Secondary SEVERE versus Ambient
resd14._MSvAS_MAvAA<- results(dds.d14, name="Second_TreamentS.Primary_TreatmentM", alpha = 0.05) # (MSvAS × MAvAA)
table(resd14._MSvAS_MAvAA$padj<0.05) # 0 DEGs

# Interaction term for the main  effect (Primary Treatment) in Secondary SEVERE  versus MODERATE
resd14._MSvAS_MMvAM= results(dds.d14, contrast=list("Second_TreamentS.Primary_TreatmentM", "Second_TreamentM.Primary_TreatmentM"), alpha = 0.05) # (MSvAS × MMvAM)
table(resd14._MSvAS_MMvAM$padj<0.05) # 0 DEGs


# ~ Main EFFECTS  ======================================================================================================================== #
resultsNames(dds.d14.main) # view the names of your results model 

# Main effect of Primary treatment (WITHOUT ref level in other factor - no interaction term in this model!!!)
resd14._P.MvA<- results(dds.d14.main, name="Primary_Treatment_M_vs_A", alpha = 0.05) 
table(resd14._P.MvA$padj<0.05) # 510 DEGs
resd14._P.MvA <- resd14._P.MvA[order(resd14._P.MvA$padj), ] ## Order by adjusted p-value
sum((resd14._P.MvA$log2FoldChange[1:510] >= 1) == TRUE) # 70 DEGs upregulated  (LFC >= 1)
( sum((resd14._P.MvA$log2FoldChange[1:510] >= 1) == TRUE) / (table(resd14._P.MvA$padj<0.05))[2] ) * 100 # 13.72549 % DEGs upregulated (N = 9141)
sum((resd14._P.MvA$log2FoldChange[1:510] <= -1) == TRUE) # 335 DEGs downregulated  (LFC <= -1)
( sum((resd14._P.MvA$log2FoldChange[1:510] <= -1) == TRUE) / (table(resd14._P.MvA$padj<0.05))[2] ) * 100 # 65.68627 % DEGs downregulated (N = 9141)

# What is the difference between secondary treatment 'moderate' WITHOUT considering the primary treament history?
resd14._S.MvA<- results(dds.d14.main, name="Second_Treament_M_vs_A", alpha = 0.05)
table(resd14._S.MvA$padj<0.05) # 0 DEGs

# What is the difference between secondary treatment 'severe' WITHOUT considering the primary treament history?
resd14._S.SvA<- results(dds.d14.main, name="Second_Treament_S_vs_A", alpha = 0.05)
table(resd14._S.SvA$padj<0.05) # 0 DEGs

# ~ Group - MAIN EFFECTS(ACCOUNTING FOR SECONDARY-TREATMENT SPECIFIC EFFECTS) ============================================================================================ #
resultsNames(dds.d14.group) # view the names of your results model 

# ALL primary history regarless of subsequent exposure (Ambient versus Moderate)
res.d14_P.MvA_group <- results(dds.d14.group, contrast = list(c("All_TreatmentMA", "All_TreatmentMM", "All_TreatmentMS"),c("All_TreatmentAA", "All_TreatmentAM", "All_TreatmentAS")),  alpha= 0.05)
table(res.d14_P.MvA_group$padj<0.05) # 511 DEGs
res.d14_P.MvA_group <- res.d14_P.MvA_group[order(res.d14_P.MvA_group$padj), ] ## Order by adjusted p-value
sum((res.d14_P.MvA_group$log2FoldChange[1:511] >= 1) == TRUE) # 111 DEGs upregulated  (LFC >= 1)
( sum((res.d14_P.MvA_group$log2FoldChange[1:511] >= 1) == TRUE) / (table(res.d14_P.MvA_group$padj<0.05))[2] ) * 100 # 21.72211 % DEGs upregulated (N = 9141)
sum((res.d14_P.MvA_group$log2FoldChange[1:511] <= -1) == TRUE) #  396 DEGs downregulated  (LFC <= -1)
( sum((res.d14_P.MvA_group$log2FoldChange[1:511] <= -1) == TRUE) / (table(res.d14_P.MvA_group$padj<0.05))[2] ) * 100 # 77.49511 % DEGs downregulated (N = 9141)

# CONSTANT exposure to SAME treatment
resd14._AAvMM_group <- results(dds.d14.group, contrast = list(("All_TreatmentAA"),("All_TreatmentMM")),  alpha= 0.05) # Constant expore (basically M v A Primaru + 14 more days)
table(resd14._AAvMM_group$padj<0.05) # 11 DEGs
# SECOND exposure by PRIMARY treatment
resd14._AAvMA_group <- results(dds.d14.group, contrast = list(("All_TreatmentAA"),("All_TreatmentMA")),  alpha= 0.05) # Second Ambient
table(resd14._AAvMA_group$padj<0.05) # 4 DEGs
resd14._AMvMM_group <- results(dds.d14.group, contrast = list(("All_TreatmentAM"),("All_TreatmentMM")),  alpha= 0.05) # Second Moderate
table(resd14._AMvMM_group$padj<0.05) # 12 DEGs
resd14._ASvMS_group <- results(dds.d14.group, contrast = list(("All_TreatmentAS"),("All_TreatmentMS")),  alpha= 0.05) # Second Severe
table(resd14._ASvMS_group$padj<0.05) # 2 DEGs
# AMBIENT history and subseqent treatment 
resd14._AAvAM_group <- results(dds.d14.group, contrast = list(("All_TreatmentAA"),("All_TreatmentAM")),  alpha= 0.05) # Amb - response to Moderate
table(resd14._AAvAM_group$padj<0.05) # 0 DEGs
resd14._AAvAS_group <- results(dds.d14.group, contrast = list(("All_TreatmentAA"),("All_TreatmentAS")),  alpha= 0.05) # Amb - response to Severe
table(resd14._AAvAS_group$padj<0.05) # 0 DEGs
resd14._AMvAS_group <- results(dds.d14.group, contrast = list(("All_TreatmentAM"),("All_TreatmentAS")),  alpha= 0.05) # Amb - response to moderate versus Severe
table(resd14._AMvAS_group$padj<0.05) # 0 DEGs
# MODERATE history and subseqent treatment 
resd14._MMvMA_group <- results(dds.d14.group, contrast = list(("All_TreatmentMM"),("All_TreatmentMA")),  alpha= 0.05) # Mod - response to Ambient
table(resd14._MMvMA_group$padj<0.05) # 0 DEGs
resd14._MMvMS_group <- results(dds.d14.group, contrast = list(("All_TreatmentMM"),("All_TreatmentMS")),  alpha= 0.05) # Mod - response to Severe
table(resd14._MMvMS_group$padj<0.05) # 0 DEGs
resd14._MAvMS_group <- results(dds.d14.group, contrast = list(("All_TreatmentMA"),("All_TreatmentMS")),  alpha= 0.05) # Amb - response to ambient versus Severe
table(resd14._MAvMS_group$padj<0.05) # 0 DEGs
# OTHER pairwise interactions
#AA
resd14._AAvMM_group <- results(dds.d14.group, contrast = list(("All_TreatmentAA"),("All_TreatmentMM")),  alpha= 0.05) # Constant expore (basically M v A Primaru + 14 more days)
table(resd14._AAvMM_group$padj<0.05) # 11 DEGs
resd14._AAvMA_group <- results(dds.d14.group, contrast = list(("All_TreatmentAA"),("All_TreatmentMA")),  alpha= 0.05)
table(resd14._AAvMA_group$padj<0.05) # 2 DEGs
resd14._AAvMS_group <- results(dds.d14.group, contrast = list(("All_TreatmentAA"),("All_TreatmentMS")),  alpha= 0.05)
table(resd14._AAvMS_group$padj<0.05) # 39 DEGs
#AM
resd14._AMvMA_group <- results(dds.d14.group, contrast = list(("All_TreatmentAM"),("All_TreatmentMA")),  alpha= 0.05)
table(resd14._AMvMA_group$padj<0.05) # 16 DEGs
resd14._AMvMS_group <- results(dds.d14.group, contrast = list(("All_TreatmentAM"),("All_TreatmentMS")),  alpha= 0.05)
table(resd14._AMvMS_group$padj<0.05) # 0 DEGs
#AS
resd14._ASvMA_group <- results(dds.d14.group, contrast = list(("All_TreatmentAS"),("All_TreatmentMA")),  alpha= 0.05)
table(resd14._ASvMA_group$padj<0.05) # 11 DEGs
resd14._ASvMM_group <- results(dds.d14.group, contrast = list(("All_TreatmentAS"),("All_TreatmentMM")),  alpha= 0.05)
table(resd14._ASvMM_group$padj<0.05) # 1 DEGs
# Second treatment effects (regarless of prior history)
#A v M
resd14._S.AvM_group <- results(dds.d14.group, contrast = list(c("All_TreatmentAA", "All_TreatmentMA"),c("All_TreatmentAM", "All_TreatmentMM")),  alpha= 0.05)
table(resd14._S.AvM_group$padj<0.05) # 0 DEGs

# A v S
resd14._S.AvS_group <- results(dds.d14.group, contrast = list(c("All_TreatmentAA", "All_TreatmentMA"),c("All_TreatmentAS", "All_TreatmentMS")),  alpha= 0.05)
table(resd14._S.AvS_group$padj<0.05) # 0  DEGs

#M v S
resd14._S.MvS_group <- results(dds.d14.group, contrast = list(c("All_TreatmentAM", "All_TreatmentMM"),c("All_TreatmentAS", "All_TreatmentMS")),  alpha= 0.05)
table(resd14._S.MvS_group$padj<0.05) # 0 DEGs



# Write results
# path_out.d14 = 'C:/Users/samjg/Documents/My_Projects/Pgenerosa_TagSeq_Metabolomics/TagSeq/RAnalysis/DESeq2/output/Day14/'
# write.csv(resdata.d14.all.prim.M_vs_A, paste(path_out.d14, "Day14.PrimaryTreament_diffexpr-results.csv"))
# write.csv(resdata.d14.second_M_vs_A, paste(path_out.d14, "Day14.SecondTreament_MvsA_diffexpr-results.csv"))
# write.csv(resdata.d14.all.second_S_vs_A, paste(path_out.d14, "Day14.SecondTreament_SvsA_diffexpr-results.csv"))



# volcano plot ------------------------------------------------------------------------------------------------------ #
# Effect of primary treatment 'res.d14_P.MvA_group' - 111 DEGs
png("RAnalysis/DESeq2/output/Day14/Day14.Primary_VolcanoPlot.png", 1000, 1000, pointsize=20)
EnhancedVolcano(res.d14_P.MvA_group,
                lab = rownames(res.d14_P.MvA_group),
                x = 'log2FoldChange',
                y = 'padj',
                title = 'Day14 Primary treatment (M v A)',
                subtitle = "DESeq2 - Differential expression",
                FCcutoff = 1,
                pCutoff = 0.05,
                pointSize = 4.0,
                labSize = 6.0,
                colAlpha = 1,
                legendPosition = 'right',
                legendLabSize = 12,
                legendIconSize = 4.0,
                widthConnectors = 0.145)
dev.off()

# Effect of second treatment M v. A 'resd14._S.AvM_group' - 114 DEGs
png("RAnalysis/DESeq2/output/Day14/Day14.SecondAvM_VolcanoPlot.png", 1000, 1000, pointsize=20)
EnhancedVolcano(resd14._S.AvM_group,
                lab = rownames(resd14._S.AvM_group),
                x = 'log2FoldChange',
                y = 'padj',
                title = 'Day14 Second treatment (A v M)',
                subtitle = "DESeq2 - Differential expression",
                FCcutoff = 1,
                pCutoff = 0.05,
                pointSize = 4.0,
                labSize = 6.0,
                colAlpha = 1,
                legendPosition = 'right',
                legendLabSize = 12,
                legendIconSize = 4.0,
                widthConnectors = 0.145)
dev.off()

# Effect of second treatment S v. A 'resd14._S.AvS_group' - 1 DEG
png("RAnalysis/DESeq2/output/Day14/Day14.SecondAvS_VolcanoPlot.png", 1000, 1000, pointsize=20)
EnhancedVolcano(resd14._S.AvS_group,
                lab = rownames(resd14._S.AvS_group),
                x = 'log2FoldChange',
                y = 'padj',
                title = 'Day14 Second treatment (A v S)',
                subtitle = "DESeq2 - Differential expression",
                FCcutoff = 1.0,
                pCutoff = 0.05,
                pointSize = 4.0,
                labSize = 6.0,
                colAlpha = 1,
                legendPosition = 'right',
                legendLabSize = 12,
                legendIconSize = 4.0,
                widthConnectors = 0.145)
dev.off()

# Effect of second treatment M v. S 'resd14._S.MvS_group' - 1 DEG
png("RAnalysis/DESeq2/output/Day14/Day14.SecondMvs_VolcanoPlot.png", 1000, 1000, pointsize=20)
EnhancedVolcano(resd14._S.MvS_group,
                lab = rownames(resd14._S.MvS_group),
                x = 'log2FoldChange',
                y = 'padj',
                title = 'Day14 Second treatment (M v S)',
                subtitle = "DESeq2 - Differential expression",
                FCcutoff = 1.0,
                pCutoff = 0.05,
                pointSize = 4.0,
                labSize = 6.0,
                colAlpha = 1,
                legendPosition = 'right',
                legendLabSize = 12,
                legendIconSize = 4.0,
                widthConnectors = 0.145)
dev.off()

# =========================================
# Data transformations for heatmap and PCA visuals 
# ============================================================= 
rlog.d14<- rlogTransformation(dds.d14) #  full model                       wait for this to complete.... 

png("RAnalysis/DESeq2/output/Day14/Day14.rlog_histogram.png", 1000, 1000, pointsize=20)# diagnostics of transformation # Histogram and sd plot
hist(assay(rlog.d14)) # view histogram 
dev.off()
png("RAnalysis/DESeq2/output/Day14/Day14.rlog_mean_sd.png", 1000, 1000, pointsize=20)
meanSdPlot(assay(rlog.d14)) # shows the sd y axis (sq root of varaince in all samples) - flat curve may seem like a goals, BUT may be unreasonable in cases with MANY true DEGs from experimental conditions
dev.off()

# PCA plot rlog ------------------------ #
pcaData_d14 <- plotPCA(rlog.d14, intgroup = c( "Primary_Treatment", "Second_Treament"), returnData = TRUE)
percentVar <- round(100 * attr(pcaData_d14, "percentVar"))
png("RAnalysis/DESeq2/output/Day14/Day14.rlog_PCA.png", 1000, 1000, pointsize=20)
ggplot(pcaData_d14, aes(x = PC1, y = PC2, color = Primary_Treatment, shape = Second_Treament)) +
  scale_shape_manual(values = c(4, 19, 17)) +
  geom_text(aes(label=name),hjust=0.2, vjust=1.4, size=3) +
  geom_point(size =6) +
  theme_classic() +
  theme(text = element_text(size=15)) +
  ggtitle("PCA: Day14 (rlog)") +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  coord_fixed()
dev.off()

# Plot heat map rlog------------------------ #
select <- order(rowMeans(counts(dds.d14,normalized=TRUE)), decreasing=TRUE)[1:20] # normalize the counts per row and call the first 20 most expressed genes
df <- as.data.frame(colData(dds.d14)[c("Primary_Treatment", "Second_Treament")])
annotation_colors = list((Primary_Treatment = c(A="Blue", M="Orange")), Second_Treament = c(A="White", M="Grey", S="Black"))
d14.rlog.heatmap<- pheatmap(assay(rlog.d14)[select,], cluster_rows=TRUE, show_rownames=TRUE, # rlog heatmap
                           cluster_cols=TRUE, annotation_col=df, annotation_colors = annotation_colors,
                           main = "Day7.rlog_heatmap")
save_pheatmap(d14.rlog.heatmap, filename = "RAnalysis/DESeq2/output/Day14/Day14.rlog_heatmap.png")


# Plot dispersions ------------------------------------------------------------------------------------------------------ #
png("RAnalysis/DESeq2/output/Day14/Day14-dispersions.png", 1000, 1000, pointsize=20)
plotDispEsts(dds.d14, main="Day 14 data dispersions")
dev.off()

# Heat maps and Principal components analysis ============================================================================= #

#  Primary_Treatment_M_vs_A 

# ALL 511 DEGS
res14.MvA <- res.d14_P.MvA_group[c(1:511),] # 511 total DEGs - lets call this dataset for a PCA and heat map 
res14.MvA # view the last pdj - should be < 0.05 
dds.d14.MvA<- dds.d14[(rownames(res14.MvA_UPREG))]
rlog.dds.d14.MvA<- rlogTransformation(dds.d14.MvA) # rlog transformation 
# PCA plot rlog 
pcaData_d14 <- plotPCA(rlog.dds.d14.MvA, intgroup = c( "Primary_Treatment", "Second_Treament"), returnData = TRUE)
percentVar <- round(100 * attr(pcaData_d14, "percentVar"))
png("RAnalysis/DESeq2/output/Day14/Day14.DEGs.rlog_PCA.png", 1000, 1000, pointsize=20)
ggplot(pcaData_d14, aes(x = PC1, y = PC2, color = Primary_Treatment, shape = Second_Treament)) +
  scale_shape_manual(values = c(4, 19, 17)) +
  geom_text(aes(label=name),hjust=0.2, vjust=1.4, size=3) +
  geom_point(size =6) +
  theme_classic() +
  theme(text = element_text(size=15)) +
  ggtitle("PCA: Day14 DEGs only (rlog)") +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  coord_fixed()
dev.off()
# Plot heat map rlog
select <- order(rowMeans(counts(dds.d14.MvA,normalized=TRUE)), decreasing=TRUE)[1:25] # normalize the counts per row and call the first 20 most expressed genes
df <- as.data.frame(colData(dds.d14.MvA)[c("Primary_Treatment", "Second_Treament")])
annotation_colors = list((Primary_Treatment = c(A="Blue", M="Orange")), Second_Treament = c(A="White", M="Grey", S="Black"))
d14_DEGs.rlog.heatmap<- pheatmap(assay(rlog.dds.d14.MvA)[select,], cluster_rows=TRUE, show_rownames=TRUE, # rlog heatmap
                           cluster_cols=TRUE, annotation_col=df, annotation_colors = annotation_colors,
                           main = "Day14.DEGsrlog_heatmap")
save_pheatmap(d14_DEGs.rlog.heatmap, filename = "RAnalysis/DESeq2/output/Day14/Day14.DEGsrlog_heatmap.png")

# ========================================================== 
#
# DAY 21   (3 CPM in 50% samples)                          
# ========================================================== 
# view the design for each of the d7 dds 
design(dds.d21) # ~Primary_Treatment + Second_Treament + Third_Treatment + Primary_Treatment:Second_Treament + Primary_Treatment:Third_Treatment + Second_Treament:Third_Treatment
design(dds.d21.group) # ~All_Treatment -1  (group)
design(dds.d21.main) # ~Third_Treatment + Second_Treament + Primary_Treatment

# RUN DESEQ2 model - view all the pariwise comparisons
dds.d21 <- DESeq(dds.d21) #  full model                                   wait for this to complete....
dds.d21.group <- DESeq(dds.d21.group) # group model                       wait for this to complete....
dds.d21.main <- DESeq(dds.d21.main) # main effect model                   wait for this to complete....

# ~ FULL MODEL EFFECTS  ======================================================================================================================== #
resultsNames(dds.d21) # view the names of your results model 

resd21._full_1<- results(dds.d21, name="Primary_Treatment_M_vs_A", alpha = 0.05)
table(resd21._full_1$padj<0.05) # 99 DEGs

resd21._full_2<- results(dds.d21, name="Second_Treament_M_vs_A", alpha = 0.05)
table(resd21._full_2$padj<0.05) # 0 DEGs

resd21._full_3<- results(dds.d21, name="Second_Treament_S_vs_A", alpha = 0.05)
table(resd21._full_3$padj<0.05) # 0 DEGs

resd21._full_4<- results(dds.d21, name="Third_Treatment_M_vs_A", alpha = 0.05)
table(resd21._full_4$padj<0.05) #  0 DEGs

resd21._full_5<- results(dds.d21, name="Primary_TreatmentM.Second_TreamentM", alpha = 0.05)
table(resd21._full_5$padj<0.05) # 89 DEGs





resd21._full_6<- results(dds.d21, list(c("Primary_Treatment_M_vs_A","Primary_TreatmentM.Second_TreamentS")))
table(resd21._full_6$padj<0.05) # 9 DEGs

resd21._full_7<- results(dds.d21, list(c("Primary_Treatment_M_vs_A","Primary_TreatmentM.Third_TreatmentM")))
table(resd21._full_7$padj<0.05) # 10 DEGs

resd21._full_8<- results(dds.d21, list(c("Primary_Treatment_M_vs_A","Second_TreamentM.Third_TreatmentM")))
table(resd21._full_8$padj<0.05) # 34 DEGs

resd21._full_9<- results(dds.d21, list(c("Primary_Treatment_M_vs_A","Second_TreamentS.Third_TreatmentM")))
table(resd21._full_9$padj<0.05) # 178 DEGs




resd21._full_9<- results(dds.d21, list(c("Second_Treament_M_vs_A","Primary_TreatmentM.Second_TreamentS")))
table(resd21._full_9$padj<0.05) # 12 DEGs

resd21._full_9<- results(dds.d21, list(c("Second_Treament_M_vs_A","Primary_TreatmentM.Third_TreatmentM")))
table(resd21._full_9$padj<0.05) # 0 DEGs

resd21._full_9<- results(dds.d21, list(c("Second_Treament_M_vs_A","Second_TreamentM.Third_TreatmentM")))
table(resd21._full_9$padj<0.05) # 15 DEGs

resd21._full_9<- results(dds.d21, list(c("Second_Treament_M_vs_A","Second_TreamentS.Third_TreatmentM"))) 
table(resd21._full_9$padj<0.05) # 182 DEGs



resd21._full_9<- results(dds.d21, list(c("Third_Treatment_M_vs_A","Primary_TreatmentM.Second_TreamentS")))
table(resd21._full_9$padj<0.05) # 57 DEGs

resd21._full_9<- results(dds.d21, list(c("Third_Treatment_M_vs_A","Primary_TreatmentM.Third_TreatmentM")))
table(resd21._full_9$padj<0.05) # 2 DEGs

resd21._full_9<- results(dds.d21, list(c("Primary_Treatment_M_vs_A","Second_TreamentM.Third_TreatmentM")))
table(resd21._full_9$padj<0.05) # 34 DEGs

resd21._full_9<- results(dds.d21, list(c("Primary_Treatment_M_vs_A","Second_TreamentS.Third_TreatmentM")))
table(resd21._full_9$padj<0.05) # 178 DEGs




resd21._full_9<- results(dds.d21, list(c("Third_Treatment_M_vs_A","Second_TreamentS.Third_TreatmentM")))
table(resd21._full_9$padj<0.05) # 0 DEGs




# ~ Main EFFECTS  ======================================================================================================================== #
resultsNames(dds.d21.main) # view the names of your results model 

resd21._main_1<- results(dds.d21.main, name="Primary_Treatment_M_vs_A", alpha = 0.05)
table(resd21._main_1$padj<0.05) # 224 DEGs
resd21._main_1 <- resd21._main_1[order(resd21._main_1$padj), ] ## Order by adjusted p-value
sum((resd21._main_1$log2FoldChange[1:224] >= 1) == TRUE) # 3 DEGs upregulated  (LFC >= 1)
sum((resd21._main_1$log2FoldChange[1:224] <= -1) == TRUE) # 8 DEGs downregulated  (LFC <= -1)


resd21._main_2<- results(dds.d21.main, name="Second_Treament_M_vs_A", alpha = 0.05)
table(resd21._main_2$padj<0.05) # 3 DEGs

resd21._main_3<- results(dds.d21.main, name="Second_Treament_S_vs_A", alpha = 0.05)
table(resd21._main_3$padj<0.05) # 0 DEGs

resd21._main_4<- results(dds.d21.main, name="Third_Treatment_M_vs_A", alpha = 0.05)
table(resd21._main_4$padj<0.05) # 0 DEGs



# ~ Group - \ ============================================================================================ #
resultsNames(dds.d21.group) # view the names of your results model 

# PRIMARY: AMBINET VS MODERATE
res.d21_Prim_MvA_group <- results(dds.d21.group, contrast = list(c("All_TreatmentAAA", "All_TreatmentAAM", "All_TreatmentAMA", "All_TreatmentAMM", "All_TreatmentASA","All_TreatmentASM"), 
                                                                 c("All_TreatmentMAA", "All_TreatmentMAM", "All_TreatmentMMA", "All_TreatmentMMM", "All_TreatmentMSA", "All_TreatmentMSM")),  alpha= 0.05)
table(res.d21_Prim_MvA_group$padj<0.05) # 233 DEGs
res.d21_Prim_MvA_group <- res.d21_Prim_MvA_group[order(res.d21_Prim_MvA_group$padj), ] ## Order by adjusted p-value
sum((res.d21_Prim_MvA_group$log2FoldChange[1:233] >= 1) == TRUE) # 167 DEGs upregulated  (LFC >= 1)
sum((res.d21_Prim_MvA_group$log2FoldChange[1:233] <= -1) == TRUE) # 66 DEGs downregulated  (LFC <= -1)

# SECONDARY: AMBINET VS MODERATE
res.d21_Sec_AvM_group <- results(dds.d21.group, contrast = list(c("All_TreatmentAAA", "All_TreatmentAAM", "All_TreatmentMAA", "All_TreatmentMAM"), 
                                                                 c("All_TreatmentAMA", "All_TreatmentAMM", "All_TreatmentMMA", "All_TreatmentMMM")),  alpha= 0.05)
table(res.d21_Sec_AvM_group$padj<0.05) # 2 DEGs

# SECOND: AMBINET VS SEVERE
res.d21_Sec_AvS_group <- results(dds.d21.group, contrast = list(c("All_TreatmentAAA", "All_TreatmentAAM", "All_TreatmentMAA", "All_TreatmentMAM"), 
                                                                 c("All_TreatmentASA", "All_TreatmentASM", "All_TreatmentMSA", "All_TreatmentMSM")),  alpha= 0.05)
table(res.d21_Sec_AvS_group$padj<0.05) # 2 DEGs


# SECOND: MODERATE VS SEVERE
res.d21_Sec_MvS_group <- results(dds.d21.group, contrast = list(c("All_TreatmentAMA", "All_TreatmentAMM", "All_TreatmentMMA", "All_TreatmentMMM"), 
                                                                c("All_TreatmentASA", "All_TreatmentASM", "All_TreatmentMSA", "All_TreatmentMSM")),  alpha= 0.05)
table(res.d21_Sec_MvS_group$padj<0.05) # 0 DEGs


# THIRD: AMBIENT VS MODERATE
res.d21_Third_AvM_group <- results(dds.d21.group, contrast = list(c("All_TreatmentAAA", "All_TreatmentMAA", "All_TreatmentAMA", "All_TreatmentMMA", "All_TreatmentASA","All_TreatmentMSA"), 
                                                                 c("All_TreatmentAAM", "All_TreatmentMAM", "All_TreatmentAMM", "All_TreatmentMMM", "All_TreatmentASM", "All_TreatmentMSM")),  alpha= 0.05)
table(res.d21_Third_AvM_group$padj<0.05) # 0 DEGs



# PRIMARY A V m - Second Moderate +  Third Moderate
res.d21_Third_AvM_group <- results(dds.d21.group, contrast = list(("All_TreatmentASM"), 
                                                                  ("All_TreatmentMSM")),  alpha= 0.05)
table(res.d21_Third_AvM_group$padj<0.05) # 0 DEGs




# ALL PAIRWISE COMPARISONS
d21.GROUP.DEGs_table = data.frame() # start a data frame
df_d21_group<- data.frame(resultsNames(dds.d21.group))
unique.vars <- df_d21_group[1,]
for(i in 1:nrow(df_d21_group)) {
  
      var <- df_d21_group[i,]
      T_F <- (df_d21_group[,1]!=c(unique.vars,var))
      df_d21_group_2 <- data.frame(df_d21_group[T_F,])

      for(j in 1:nrow(df_d21_group_2)) {
          
          var2 <- df_d21_group_2[j,]
          res <- results(dds.d21.group, contrast = list((var),(var2)),  alpha= 0.05) 
          res.table<-as.data.frame(table(res$padj<0.05))
          DEGs_total <- res.table[2,2]
          DEGs_total[is.na(DEGs_total)] <- 0
          DEGs.table <- data.frame(matrix(nrow = 1, ncol = 7)) # create a new data table
          colnames(DEGs.table)<-c('Var1', 'Var2', 'DEGs', 'num.upreg', 'num.downreg', 'perc.upreg','perc.downreg') 
          res.order <- res[order(res$padj), ] ## Order by adjusted p-value
          DEGs.table$num.upreg <- sum((res.order$log2FoldChange[1:(DEGs_total)] >= 1) == TRUE) 
          DEGs.table$perc.upreg <- (sum((res.order$log2FoldChange[1:(DEGs_total)] >= 1) == TRUE) / (table(res.order$padj<0.05))[2] ) * 100 
          DEGs.table$num.downreg <- sum((res.order$log2FoldChange[1:(DEGs_total)] <= -1) == TRUE) 
          DEGs.table$perc.downreg <- (sum((res.order$log2FoldChange[1:(DEGs_total)] <= -1) == TRUE) / (table(res.order$padj<0.05))[2] ) * 100 
          DEGs.table$Var1 <- substr(var, 14,16) # fill date
          DEGs.table$Var2 <- substr(var2, 14,16) # fill run number 
          DEGs.table$DEGs <- DEGs_total # fill with the chosen alpha value (assigned at start of script)
          df.group <- data.frame(DEGs.table) # name dataframe for this single row
          d21.GROUP.DEGs_table <- rbind(d21.GROUP.DEGs_table, df.group) # bind to a cumulative list dataframe
          print(d21.GROUP.DEGs_table) # show loop progress in the console
          unique.vars <- unique(d21.GROUP.DEGs_table$Var1)
      } 
}
#View(d21.GROUP.DEGs_table)
# knitr::kable(d21.GROUP.DEGs_table, caption = "Day21 All pairwise DEGs in 'group' DESeq2 model")

## Write results
path_out.d21 = 'C:/Users/samjg/Documents/My_Projects/Pgenerosa_TagSeq_Metabolomics/TagSeq/RAnalysis/DESeq2/output/Day21/'# call the path
write.csv(d21.GROUP.DEGs_table, paste(path_out.d21, "Day21_all.pairwise_DEGs.csv"))
# write.csv(resdata.d21.All.Second_Treatment_M_vs_A, paste(path_out.d21, "Day21.SecondTreament_MvsA_diffexpr-results.csv"))
# write.csv(resdata.d21.All.Second_Treatment_S_vs_A, paste(path_out.d21, "Day21.SecondTreament_SvsA_diffexpr-results.csv"))
# write.csv(resdata.d21.All.Third_Treatment_M_vs_A, paste(path_out.d21, "Day21.ThirdTreament_diffexpr-results.csv"))

# volcano plot ------------------------------------------------------------------------------------------------------ #

# Effect of primary treatment 'resd14.all.primary_M_vs_A' - 221 DEGs
png("RAnalysis/DESeq2/output/Day21/Day21.PrimaryTreatment-VolcanoPlot.png", 1000, 1000, pointsize=20)
EnhancedVolcano(resd21.All.Primary_Treatment_M_vs_A,
                lab = rownames(resd21.All.Primary_Treatment_M_vs_A),
                x = 'log2FoldChange',
                y = 'padj',
                title = 'Moderate versus Ambient Treatment (Primary)',
                subtitle = "DESeq2 - Differential expression",
                FCcutoff = 1,
                pCutoff = 0.3,
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


# =========================================
# Data transformations for heatmap and PCA visuals 
# ============================================================= 
rlog.d14<- rlogTransformation(dds.d14) #  full model                       wait for this to complete.... 

png("RAnalysis/DESeq2/output/Day14/Day14.rlog_histogram.png", 1000, 1000, pointsize=20)# diagnostics of transformation # Histogram and sd plot
hist(assay(rlog.d14)) # view histogram 
dev.off()
png("RAnalysis/DESeq2/output/Day14/Day14.rlog_mean_sd.png", 1000, 1000, pointsize=20)
meanSdPlot(assay(rlog.d14)) # shows the sd y axis (sq root of varaince in all samples) - flat curve may seem like a goals, BUT may be unreasonable in cases with MANY true DEGs from experimental conditions
dev.off()

# PCA plot rlog ------------------------ #
pcaData_d14 <- plotPCA(rlog.d14, intgroup = c( "Primary_Treatment", "Second_Treament"), returnData = TRUE)
percentVar <- round(100 * attr(pcaData_d14, "percentVar"))
png("RAnalysis/DESeq2/output/Day14/Day14.rlog_PCA.png", 1000, 1000, pointsize=20)
ggplot(pcaData_d14, aes(x = PC1, y = PC2, color = Primary_Treatment, shape = Second_Treament)) +
  scale_shape_manual(values = c(4, 19, 17)) +
  geom_text(aes(label=name),hjust=0.2, vjust=1.4, size=3) +
  geom_point(size =6) +
  theme_classic() +
  theme(text = element_text(size=15)) +
  ggtitle("PCA: Day14 (rlog)") +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  coord_fixed()
dev.off()

# Plot heat map rlog------------------------ #
select <- order(rowMeans(counts(dds.d14,normalized=TRUE)), decreasing=TRUE)[1:20] # normalize the counts per row and call the first 20 most expressed genes
df <- as.data.frame(colData(dds.d14)[c("Primary_Treatment", "Second_Treament")])
annotation_colors = list((Primary_Treatment = c(A="Blue", M="Orange")), Second_Treament = c(A="White", M="Grey", S="Black"))
d14.rlog.heatmap<- pheatmap(assay(rlog.d14)[select,], cluster_rows=TRUE, show_rownames=TRUE, # rlog heatmap
                            cluster_cols=TRUE, annotation_col=df, annotation_colors = annotation_colors,
                            main = "Day7.rlog_heatmap")
save_pheatmap(d14.rlog.heatmap, filename = "RAnalysis/DESeq2/output/Day14/Day14.rlog_heatmap.png")


# Plot dispersions ------------------------------------------------------------------------------------------------------ #
png("RAnalysis/DESeq2/output/Day14/Day14-dispersions.png", 1000, 1000, pointsize=20)
plotDispEsts(dds.d14, main="Day 14 data dispersions")
dev.off()

# Heat maps and Principal components analysis ============================================================================= #

#  Primary_Treatment_M_vs_A 

# ALL 511 DEGS
res14.MvA <- res.d14_P.MvA_group[c(1:511),] # 511 total DEGs - lets call this dataset for a PCA and heat map 
res14.MvA # view the last pdj - should be < 0.05 
dds.d14.MvA<- dds.d14[(rownames(res14.MvA_UPREG))]
rlog.dds.d14.MvA<- rlogTransformation(dds.d14.MvA) # rlog transformation 
# PCA plot rlog 
pcaData_d14 <- plotPCA(rlog.dds.d14.MvA, intgroup = c( "Primary_Treatment", "Second_Treament"), returnData = TRUE)
percentVar <- round(100 * attr(pcaData_d14, "percentVar"))
png("RAnalysis/DESeq2/output/Day14/Day14.DEGs.rlog_PCA.png", 1000, 1000, pointsize=20)
ggplot(pcaData_d14, aes(x = PC1, y = PC2, color = Primary_Treatment, shape = Second_Treament)) +
  scale_shape_manual(values = c(4, 19, 17)) +
  geom_text(aes(label=name),hjust=0.2, vjust=1.4, size=3) +
  geom_point(size =6) +
  theme_classic() +
  theme(text = element_text(size=15)) +
  ggtitle("PCA: Day14 DEGs only (rlog)") +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  coord_fixed()
dev.off()
# Plot heat map rlog
select <- order(rowMeans(counts(dds.d14.MvA,normalized=TRUE)), decreasing=TRUE)[1:25] # normalize the counts per row and call the first 20 most expressed genes
df <- as.data.frame(colData(dds.d14.MvA)[c("Primary_Treatment", "Second_Treament")])
annotation_colors = list((Primary_Treatment = c(A="Blue", M="Orange")), Second_Treament = c(A="White", M="Grey", S="Black"))
d14_DEGs.rlog.heatmap<- pheatmap(assay(rlog.dds.d14.MvA)[select,], cluster_rows=TRUE, show_rownames=TRUE, # rlog heatmap
                                 cluster_cols=TRUE, annotation_col=df, annotation_colors = annotation_colors,
                                 main = "Day14.DEGsrlog_heatmap")
save_pheatmap(d14_DEGs.rlog.heatmap, filename = "RAnalysis/DESeq2/output/Day14/Day14.DEGsrlog_heatmap.png")


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
# Day ALL Moderate vs. Ambient  ====================================================================================== #
#====================================================================================================================== 

resdata.ALL.MvA

resdata.ALL.MvA # 971 DEGS (log2FoldChange > 1; < -1 & padj < 0.05) on day 14 in response to Primary treatment (initial Mdoerate vs. Ambient conditioning)
DEGS.ALL.MvA<- resdata.ALL.MvA %>%  dplyr::filter(padj<0.05)
nrow(DEGS.ALL.MvA) # 971 - we have all DEGs now

UPREG.ALL.MvA<- DEGS.ALL.MvA %>%  dplyr::filter(log2FoldChange > 1) # call upregulated genes
DWNREG.ALL.MvA <- DEGS.ALL.MvA %>%  dplyr::filter(log2FoldChange < 1) # call downregulated genes 
nrow(UPREG.ALL.MvA) + nrow(DWNREG.ALL.MvA) # should be equal to 971


ALL.MvA.UPREG_GO <- merge(UPREG.ALL.MvA, annotation.df, by = "Gene")
ALL.MvA.UPREG_WEGO <- ALL.MvA.UPREG_GO %>% dplyr::select(c('Gene', 'GO.terms'))
ALL.MvA.UPREG_WEGO$GO.terms[is.na(ALL.MvA.UPREG_WEGO$GO.terms)] <- " " # make NA blank
# View(D7.UPREG_WEGO) # view the data

ALL.MvA.DWNREG_GO <- merge(DWNREG.ALL.MvA, annotation.df, by = "Gene")
ALL.MvA.DWNREG_WEGO <- ALL.MvA.DWNREG_GO %>% dplyr::select(c('Gene', 'GO.terms'))
ALL.MvA.DWNREG_WEGO$GO.terms[is.na(ALL.MvA.DWNREG_WEGO$GO.terms)] <- " " # make NA blank
# View(D7.DWNREG_WEGO) # view the data

# write to GO folder
write.table(ALL.MvA.UPREG_WEGO, file = "RAnalysis/GO/ALL_Upreg_Primary_WEGO", quote=FALSE, row.names = FALSE, col.names = FALSE) # call at 'Native Format' on  https://wego.genomics.cn/ 
write.table(ALL.MvA.DWNREG_WEGO, file = "RAnalysis/GO/ALL_Downreg_Primary_WEGO", quote=FALSE, row.names = FALSE, col.names = FALSE) # call at 'Native Format' on  https://wego.genomics.cn/ 


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

