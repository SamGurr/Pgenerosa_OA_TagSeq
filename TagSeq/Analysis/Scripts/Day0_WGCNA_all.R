---
  # title: "Day0_WGCNA_all"
  # author: "Samuel Gurr"
  # date: "January 8, 2021"
---
  
# LOAD PACKAGES
library(WGCNA) # note: this was previously installed with the command `BiocManager::install("WGCNA")`
library(dplyr)
library(zoo)
library(DESeq2)

# for heatmap 
# library(devtools)
# install_github("jokergoo/ComplexHeatmap") first run these - commented out to avoid running a second time...
library(ComplexHeatmap)
library(circlize)
library(reshape)
library(ggplot2)
library(hrbrthemes)


# SET WORKING DIRECTORY AND LOAD DATA
setwd("C:/Users/samjg/Documents/My_Projects/Pgenerosa_TagSeq_Metabolomics/TagSeq/")
# LOAD DATA
# Tagaseq filtered counts 
day0.counts.matrix <- read.csv(file="Analysis/Data/Filtered_counts/10cpm_50perc/day0.counts.filtered_10cpm50perc.csv", sep=',', header=TRUE)
# Treatment and Phenotype data
Master.Treatment.data <- read.csv(file="Analysis/Data/Experiment_Metadata/Master_Exp.Treatment_Metadata.csv", sep=',', header=TRUE)
d0.Treatment_data     <- Master.Treatment_Phenotype.data %>%  dplyr ::filter(Date %in% 20190724) # split for day 0 data 

# The following setting is important, do not omit.
options(stringsAsFactors = FALSE)

# ===================================================================================
#
#
# Day 0 WGCNA - PREPROCESSING THE DATA INPUT 
#
#  Read here for the pre processing steps using WGCNA!
#  https://bmcbioinformatics.biomedcentral.com/track/pdf/10.1186/1471-2105-9-559.pdf
# ===================================================================================


# filter low values - note: this is pre filteres for < 5 CPM in 50% of samples - less strict than the DESeq2 analysis
dim(day0.counts.matrix) #  9207  rows (genes) -   9  samples  not counting 'x' and 'Gene.ID'
(day0.counts.matrix)[1] # ommit the first line and transpose in the next line 
d0.data = as.data.frame(t(day0.counts.matrix[, -(1)])) # ommit all columns but samples and transpose
dim(d0.data) # 8 9207
# fix(d0.data) # view the data


# trait data ========================================================== #

# Phenotype trait and Treatment data
dim(d0.Treatment_data) #  8 rows and 10 columns
names(d0.Treatment_data) # look at the 10 columns 
d0.Treatment_data_2 = d0.Treatment_data[, -c(1:2,4:5,8:10)]; # remove columns that hold information we do not need. (i.e. tanks, Third treatment, All treatment group, etc.)
dim(d0.Treatment_data_2) # 8 rows and  2 columns
# fix(d0.Treatment_Phenotype.data)  # view the data


# count data  ========================================================== #


# fix(d21.data) # view the data - as you see all columns are now genes but filled as V1, V2, V3...
names(d0.data) = day0.counts.matrix$X # assigns column names (previous jsut numbered) as the gene ID 
rownames(d0.data) = names(day0.counts.matrix)[-(1)]; # assigns the row names as the sample ID
d0.data_matrix <- data.frame(day0.counts.matrix[,-1], row.names=day0.counts.matrix[,1]) 
d0.data_matrix_t <- t(d0.data_matrix)
# fix(d0.data_matrix_t)  # view the data

# create dds objects to transform data 
dds.d0 <- DESeqDataSetFromMatrix(countData = d0.data_matrix,
                                  colData = d0.Treatment_data_2, design = ~ 1) # DESeq Data Set (dds)
# DESeq Data Set (dds)
# NOTE: ~1 stands for no design; user will need to add a design for differential testing
# however for our purpose of just creating an onbject to transform, we do not need a design here...


# transform the data 
dds.d0_vst <- vst(dds.d0) # transform it vst
dds.d0_vst <- assay(dds.d0_vst) # call only the transformed coutns in the dds object
# fix(dds.d0_vst)
dds.d0_vst <- t(dds.d0_vst) # transpose columns to rows and vice versa

# ===================================================================================
#
# Day 0 WGCNA Sample tree 
#
# ===================================================================================

dim(dds.d0_vst) #  9207 genes; 8  samples

gsg = goodSamplesGenes(dds.d0_vst, verbose = 3);  # We first check for genes and samples with too many missing values:
gsg$allOK # If the statement returns TRUE, all genes have passed the cuts. 

# determine outlier samples 
sampleTree = hclust(dist(dds.d0_vst), method = "average") # Next we cluster the samples (in contrast to clustering genes that will come later)  to see if there are any obvious outliers.
sizeGrWindow(12,9) # The user should change the dimensions if the window is too large or too small.
par(cex = 0.6);
par(mar = c(0,4,2,0))

png("Analysis/Output/WGCNA/subseq_treatments_all/Day0/day0_ClusterTree_Precut.png", 1000, 1000, pointsize=20)
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5, 
     cex.axis = 1.5, cex.main = 2) # appears there are two outliers SG55 and SG 105; can remove by hand or an automatic appraoch 
abline(h = 110, col = "red") # add line to plot to show the cut-off od outlier samples (40000) SG105 and SG55
dev.off()

# NOTE: The following script is used to cut outliers from visual observation  of the ClusterTree
# However.. since there are only 8 samples for analysis here (minimum of 15 is recommended by creaters of WGCNA)
# we will proceed WITHOUT omitting samples

# HASHED OUT, NOT USED FOR DAY 0 WGCNA....(copied from day 7 WGCNA)
# # cut ethe tree and ommit from the reads data 
# clust = cutreeStatic(sampleTree, cutHeight = 110, minSize = 10) # Determine cluster under the line
# table(clust) # 0 = cut; 1 = kept; says it will cut 2 and save 60 
# # 'keepsamples' boolean to call the main dataset
# keepSamples = (clust==1) # call the 30 from clust at position 1
# # view the ommited samples from the tree
# omittedsamples = (clust==0)
# omSamples = dds.d0_vst[omittedsamples,] 
# nrow(omSamples) # should be 1
# rownames(omSamples) # should be 6 - "SG90" "SG89" "SG62" "SG21" "SG22" "SG59" What treatments are these??
# omSamplesData <- d0.Treatment_Phenotype.data %>% dplyr::filter(Sample.Name %in% (rownames(omSamples))) %>% dplyr::select(c('Sample.Name', 'Primary_Treatment', 'Second_Treament')) # three Primary Ambient and three Primary Moderate 
# # integrate keepsamples to the main read dataset 
# dds.d0_vst = dds.d0_vst[keepSamples, ] # integreat the boolean 'keepsamples' to ommit oultilers determined in the sample tree above 
# nGenes = ncol(dds.d0_vst) # number of genes == 8548 
# nSamples = nrow(dds.d0_vst) # number of samples == 35  - the cut tree removed 6 samples 

# # plot the tree with the 'keep samples' bollean (T/F) to ommit outlier samples 
#  sampleTree2 = hclust(dist(dds.d0_vst), method = "average") # Next we cluster the samples (in contrast to clustering genes that will come later)  to see if there are any obvious outliers.
#  png("Analysis/Output/WGCNA/subseq_treatments_all/Day0/day0_ClusterTree_Postcut.png", 1000, 1000, pointsize=20)
#  plot(sampleTree2, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5, 
#       cex.axis = 1.5, cex.main = 2)
# # dev.off()


# Form a data frame analogous to expression data that will hold the clinical traits.
d0.Samples = rownames(dds.d0_vst);# start new variable 'd21.Samples' calling the row names of the gene data (sample as 'Mouse' in trait data)
TreatRows = match(d0.Samples, d0.Treatment_data_2$Sample.Name); # match the names 
d0.Traits = d0.Treatment_data_2[TreatRows, -1]; # removes the row numbers
rownames(d0.Traits) = d0.Treatment_data_2[TreatRows, 1]; # inserts the new traitRows - matches sample treatments
dim(d0.Traits) #  8 Samples 2 columns (primary and second treatment) - note: no Third treatment because that occured Days 14-21, after this sampling day
all(rownames(d0.Traits) == rownames(dds.d0_vst))  # should be TRUE
dim(d0.Traits) #  8 2


# ===================================================================================
#
# Prepare Trait data (Treatment groups)
# ===================================================================================

# ALL TRAITS 
d0.Traits

# primary groups 
d0.Traits.Primary <-  d0.Traits %>% dplyr::select('Primary_Treatment')

d0.Traits.Primary$A <- (d0.Traits.Primary$Primary_Treatment == "A")
d0.Traits.Primary$A   <- as.numeric(d0.Traits.Primary$A)
d0.Traits.Primary$A <- as.factor(d0.Traits.Primary$A)

d0.Traits.Primary$M <- (d0.Traits.Primary$Primary_Treatment == "M")
d0.Traits.Primary$M   <- as.numeric(d0.Traits.Primary$M)
d0.Traits.Primary$M <- as.factor(d0.Traits.Primary$M)

d0.Traits.Primary <- d0.Traits.Primary[,c(2:3)] # final dataset of 0,1 for teatment groups - Primary only!

# ===================================================================================
#
# cluster samples by Trait
# ===================================================================================

# Primary treatment Only
png("Analysis/Output/WGCNA/subseq_treatments_all/Day0/day0_ClusterTree_PrimaryTreatment.png", 1000, 1000, pointsize=20)
traitColors_Primary = labels2colors(d0.Traits.Primary); # Convert traits to a color representation: white means low, red means high, grey means missing entry
plotDendroAndColors(sampleTree2, traitColors_Primary, # Plot the sample dendrogram and the colors underneath.
                    groupLabels = names(d0.Traits.Primary), 
                    main = "Sample dendrogram and trait heatmap")
dev.off()

# save data
save(dds.d0_vst, d0.Traits, file = "Analysis/Output/WGCNA/subseq_treatments_all/Day0/d.7-dataInput.RData")

# write the vst transformed data 
write.csv(dds.d0_vst, "Analysis/Output/WGCNA/subseq_treatments_all/Day0/day0_vstTransformed_WGCNAdata.csv") # write

# ===================================================================================
#
# soft threshold
# ===================================================================================

dim(dds.d0_vst) #  8 9207
# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to=20, by=2))
# Call the network topology analysis function
sft = pickSoftThreshold(dds.d0_vst, powerVector = powers, verbose = 5) #...wait for this to finish
# pickSoftThreshold 
#  performs the analysis of network topology and aids the
# user in choosing a proper soft-thresholding power.
# The user chooses a set of candidate powers (the function provides suitable default values)
# function returns a set of network indices that should be inspected
sizeGrWindow(9, 5) # set window size 
# png to output 
png("Analysis/Output/WGCNA/subseq_treatments_all/Day0/day0_ScaleTopology_SoftThresh.png", 1000, 1000, pointsize=20)
par(mfrow = c(1,2));
cex1 = 0.9;
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red") # look at at cut off at power of 6
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
dev.off() # output 


# The left panel... shows the scale-free fit index (y-axis) 
# as a function of the soft-thresholding power (x-axis). 
# We choose the power 9, which is the lowest power for which the scale-free 
# topology index curve attens out upon reaching a high value (in this case, roughly 0.90).
# The right panel.... displays the mean connectivity
# (degree, y-axis) as a function of the soft-thresholding power (x-axis).


#=====================================================================================
#
#  Satrt the step-wise module construction:  
# Step 1 = create adjacency matrix 
# https://peterlangfelder.com/2018/11/25/signed-or-unsigned-which-network-type-is-preferable/
# https://www.rdocumentation.org/packages/WGCNA/10cpm/versions/1.69/topics/adjacency
# https://ramellose.github.io/networktutorials/wgcna.html
# https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/10cpm/TechnicalReports/signedTOM.pdf
#=====================================================================================
softPower = 20 # set your soft threshold based on the plots above 
# NOTE: a softpower of 20 is VERY HIGH relative to the matrices with >15 samples run for Days 7 14 and 21 TagSeq data
# perhaps this is where the 'noise' is shown in the data as claimed by the WGCNA creators. 
# I will continue to run the rest of this script, but keep this in mind moving forward! 

# signed 
adjacency_sign = adjacency(dds.d0_vst, power = softPower, type="signed") # this takes a long time.. just wait...

#=====================================================================================
#
#  Step 2: Turn adjacency into topological overlap
# Calculation of the topological overlap matrix, (TOM)
# and the corresponding dissimilarity, from a given adjacency matrix.
#=====================================================================================
# signed 
TOM_sign = TOMsimilarity(adjacency_sign, TOMType="signed")  # this takes a long time.. just wait...

dissTOM_sign   = 1-TOM_sign

#=====================================================================================
#
#  Step 3:Call the hierarchical clustering function - plot the tree
#
#=====================================================================================
# Call the hierarchical clustering function
geneTree_sign   = hclust(as.dist(dissTOM_sign), method = "average");

# Plot the resulting clustering tree (dendrogram)
sizeGrWindow(12,9)

plot(geneTree_sign, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity - SIGNED",
     labels = FALSE, hang = 0.04);

#=====================================================================================
#
#  Step 4: Set module size and 'cutreeDynamic' to create clusters 
#
#=====================================================================================
# We like large modules, so we set the minimum module size relatively high:
minModuleSize = 100; # set this for the subseqent call...

dynamicMods_sign = cutreeDynamic(dendro = geneTree_sign, distM = dissTOM_sign,
                            deepSplit = 1, pamRespectsDendro = FALSE,
                            minClusterSize = minModuleSize);
table(dynamicMods_sign) # number of genes per module 

#=====================================================================================
#
#  Step 5: convert numeric network to colors and plot the dendrogram
#
#=====================================================================================

# Convert numeric lables into colors
dynamicColors_sign = labels2colors(dynamicMods_sign) # add colors to module labels (previously numbers)
table(dynamicColors_sign) # lets look at this table...
# Plot the dendrogram and colors underneath
sizeGrWindow(8,6)
png("Analysis/Output/WGCNA/subseq_treatments_all/Day0/day0_GeneDendrogram.png", 1000, 1000, pointsize=20)
plotDendroAndColors(geneTree_sign, dynamicColors_sign, "Dynamic Tree Cut - SIGNED",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram and module colors 'SIGNED'")
dev.off()

#=====================================================================================
#
#  Step 6: Calculate Eigengenes - view thier connectivity based on 'MEDiss = 1-cor(MEs)'
#
#=====================================================================================
# Calculate eigengenes
# MEList = moduleEigengenes(dds.d0_vst, colors = dynamicColors)
MEList = moduleEigengenes(dds.d0_vst, colors = dynamicColors_sign)
#MEList = moduleEigengenes(dds.d0_vst, colors = dynamicColors_unsign)
MEs = MEList$eigengenes
# Calculate dissimilarity of module eigengenes
MEDiss = 1-cor(MEs);
# Cluster module eigengenes
METree = hclust(as.dist(MEDiss), method = "average");
# Plot the result
sizeGrWindow(7, 6)
png("Analysis/Output/WGCNA/subseq_treatments_all/Day0/day0_ClusterEigengenes.png", 1000, 1000, pointsize=20)
plot(METree, main = "Clustering of module eigengenes - SIGNED (dissimilarity calc = MEDiss = 1-cor(MEs))",
     xlab = "", sub = "")
MEDissThres = 0.55 
abline(h=MEDissThres, col = "red")
dev.off()

#=====================================================================================
#
#  Step 7: Specify the cut line for the dendrogram (module) - Calc MODULE EIGENGENES (mergeMEs)
#
#=====================================================================================
MEDissThres = 0.55 # outide of the height threshold - will not affect any modules at this point
# Plot the cut line into the dendrogram
abline(h=MEDissThres, col = "red")
# Call an automatic merging function
# merge = mergeCloseModules(dds.d0_vst, dynamicColors, cutHeight = MEDissThres, verbose = 3)
merge = mergeCloseModules(dds.d0_vst, dynamicColors_sign, cutHeight = MEDissThres, verbose = 3)
# The merged module colors
mergedColors = merge$colors;
# Eigengenes of the new merged modules:
mergedMEs = merge$newMEs;

# Cluster module eigengenes
MEDiss2 = 1-cor(mergedMEs);
METree2 = hclust(as.dist(MEDiss2), method = "average");
# Plot the result
sizeGrWindow(7, 6)
png("Analysis/Output/WGCNA/subseq_treatments_all/Day0/day0_ClusterEigengenes_merged.png", 1000, 1000, pointsize=20)
plot(METree2, main = "Clustering of module eigengenes - SIGNED (dissimilarity calc = MEDiss = 1-cor(MEs))",
     xlab = "", sub = "")
dev.off()


#=====================================================================================
#
#  Step 8: Plot dendrogram with the cut line 'MEDissThres' 
#
#=====================================================================================
sizeGrWindow(12, 9)

png("Analysis/Output/WGCNA/subseq_treatments_all/Day0/day0_ClusterDendrogram_signed.png", 1000, 1000, pointsize=20)
plotDendroAndColors(geneTree_sign, cbind(dynamicColors_sign, mergedColors),
                    c("Dynamic Tree Cut", "Merged dynamic"),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()
#=====================================================================================
#
#  Step 9: Commit to mergedcolors as 'MEs' and 'moduleColors'
#
#=====================================================================================
# Rename to moduleColors
moduleColors = mergedColors
# Construct numerical labels corresponding to the colors
colorOrder = c("grey", standardColors(50));
moduleLabels = match(moduleColors, colorOrder)-1;
MEs = mergedMEs;
# Save module colors and labels for use in subsequent parts
save(MEs, moduleLabels, moduleColors, geneTree, file = "Analysis/Output/WGCNA/subseq_treatments_all/Day0/day0-networkConstruction-stepByStep.RData")
# write csv - save the module eigengenes
write.csv(MEs, file = "Analysis/Output/WGCNA/subseq_treatments_all/Day0/d0.WGCNA_ModulEigengenes.csv")
table(mergedColors)

#=====================================================================================
#
#  Constructing the gene network and identifying modules - NOT USED ANYMORE BUT HAS GOOD INFO 
#
#=====================================================================================

# Constructing the gene network and identifying modules is now a simple function call:
# soft thresholding power 4,
# relatively large minimum module size of 30,
# medium sensitivity (deepSplit=2) to cluster splitting.
# mergeCutHeight is the threshold for merging of modules.
# instructed the function to return numeric, rather than color, labels for modules,
# and to save the Topological Overlap Matrix

# NOTE: you will get this message "Error: REAL() can only be applied to a 'numeric', not a 'integer'"
# if your data is not all as numeric - set as numeric with lapply below BEFORE running 'blockwiseModules'
# total_cols <- length(colnames(dds.d21_vst)) # run this if ou data isnt already numeric
# dds.d21_vst[2:total_cols] = lapply(dds.d21_vst[2:total_cols], as.numeric) # run this if ou data isnt already numeric

# https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/10cpm/TechnicalReports/signedTOM.pdf
# The central idea of TOM is to count the direct connection strengths as well as connection strengths
# "mediated" by shared neighbors. The standard, or "unsigned" TOM assumes that neighbor-mediated
# connections can always be considered as "reinforcing" the direct connection. This may not always be the
# case, and the signed TOM is an attempt to take this into account
# In a signed correlation network, nodes with negative correlation are
# considered unconnected (their connection strength is zero or very close to zero). In contrast, in unsigned
# correlation networks, nodes with strong negative correlations have high connection strengths: the unsigned
# network adjacency is based on the absolute value of correlation, so positive and negative correlations are
# treated equally

# Jognson and Kelly "a soft-threshold of 12, a minimum module size of 30, a signed adjacency matrix"

# net = blockwiseModules(dds.d0_vst, #maxBlockSize = 2000,
#                        power = 3, 
#                     #  minModuleSize = 30, 
#                     #  deepSplit=4,
#                        mergeCutHeight = 0.75, #reassignThreshold = 0, 
#                        numericLabels = TRUE,
#                      #  pamRespectsDendro = FALSE,
#                        TOMType = "signed",  
#                        saveTOMs = TRUE,
#                        saveTOMFileBase = "Analysis/Output/WGCNA/subseq_treatments_all/Day0/d0.Geoduck", 
#                        verbose = 3) #... wait for this to finish 
# moduleColors = labels2colors(net$colors) # assign colors to eigengene





# Clustering tree based on eigengene modules - How related are the module eigengenes?
# represent modules by eigengenes and relating them to one another 
# dendrogram showing the 
# datME=moduleEigengenes(dds.d0_vst,moduleColors)$eigengenes
# signif(cor(datME, use="p"), 2)
# dissimME=(1-t(cor(datME, method="p")))/2 # define a dissimilarity measure between module eigengenes that keeps track of the sign of the correlation and uses this to cluster 
# hclustdatME=hclust(as.dist(dissimME), method="average" ) # dendrogram based on the module eigenge connecticity 
# Plot the eigengene dendrogram

# Hub gene of ech module 
# Both functions output a character vector of genes, where the genes are the hub gene picked for each module, 
# \and the names correspond to the module in which each gene is a hub.
# chooseTopHubInEachModule(dds.d0_vst, moduleColors, omitColors = "grey", power = 2, type = "signed")

# 'minKMEtoStay' genes whose eigengene connectivity to their module eigengene is lower than minKMEtoStay are removed from the module.
# 'minGap' minimum cluster gap given as the fraction of the difference between cutHeight and the 5th percentile of joining heights. See cutreeDynamic for more details.


# net$colors # the module assignment
# net$MEs # contains the module eigengenes of the modules.
# table(net$colors)
# there are 19  modules in order of decending size
# note 0 is reserved for genes outside of all modules (2602 genes) 

#=====================================================================================
#
#  Prepare for  module trait associations - Eigengene calc - trait data as factors
#
#=====================================================================================
# identify modules that are signiFcantly
# associated with the measured clinical traits.

# Since we already have a summary profile (eigengene) for each module,
# we simply correlate eigengenes with external traits and look for the most signicant associations:

# Define numbers of genes and samples
nGenes = ncol(dds.d0_vst); # 8548
nSamples = nrow(dds.d0_vst); # 35
# # Recalculate MEs with color labels
MEs0 = moduleEigengenes(dds.d0_vst, moduleColors)$eigengenes
MEs = orderMEs(MEs0) # reorders the columns (colors/modules)
# 
# # write csv - save the module eigengenes
# write.csv(MEs, file = "Analysis/Output/WGCNA/subseq_treatments_all/Day0/d0.WGCNA_ModulEigengenes.csv")
# 
# change chanracter treatments to integers
# ALL TRAIT DATA
d0.Traits$Primary_Treatment <- as.factor(d0.Traits$Primary_Treatment)
d0.Traits$Primary_Treatment <- as.numeric(d0.Traits$Primary_Treatment)

# d0.Traits$Second_Treament <- as.factor(d0.Traits$Second_Treament)
# d0.Traits$Second_Treament <- as.numeric(d0.Traits$Second_Treament)
# 
# # TREATMENT ONLY - break this up into more groups!! 
# d0.Traits.treat$Primary_Treatment <- as.factor(d0.Traits.treat$Primary_Treatment)
# d0.Traits.treat$Primary_Treatment <- as.numeric(d0.Traits.treat$Primary_Treatment)
# 
# d0.Traits.treat$Second_Treament <- as.factor(d0.Traits.treat$Second_Treament)
# d0.Traits.treat$Second_Treament <- as.numeric(d0.Traits.treat$Second_Treament)
# 
# d0.Traits.treat$All <- as.factor(d0.Traits.treat$All)
# d0.Traits.treat$All <- as.numeric(d0.Traits.treat$All)

#=====================================================================================
#
# Module trait correlation
#
#=====================================================================================
# ALL TRAIT DATA

dim(d0.Traits)  # 8 2
dim(MEs)  # 8 6
# moduleTraitCor = cor(MEs, d0.Traits, use = "p");
# moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);

nSamples = nrow(dds.d0_vst)

# primary 
d0.Traits.Primary$A <- as.numeric(d0.Traits.Primary$A)
d0.Traits.Primary$M <- as.numeric(d0.Traits.Primary$M)
moduleTraitCor_Primary = cor(MEs, d0.Traits.Primary, use = "p");
moduleTraitPvalue_Primary = corPvalueStudent(moduleTraitCor_Primary, nSamples);

#=====================================================================================
#
# Heatmaps
#
#=====================================================================================

# PRRIMARY TRETMENT ONLY  ------------------------------------------------------------------ # 

sizeGrWindow(10,10)
# Will display correlations and their p-values
d0.PRIMARYTreatments.matrix <-  paste(signif(moduleTraitCor_Primary, 3), "\n(",
                                       signif(moduleTraitPvalue_Primary, 3), ")", sep = "")
#dim(textMatrix) == dim(moduleTraitCor_treatonly)
par(mar = c(8, 9.5, 5, 3));
# Display the correlation values within a heatmap plot
png("Analysis/Output/WGCNA/subseq_treatments_all/Day0/heatmaps/day0_Treatments_Primary_heatmap2.png", 500, 1000, pointsize=20)
labeledHeatmap(Matrix = moduleTraitCor_Primary,
               xLabels = names(d0.Traits.Primary),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = TRUE,
               colors = blueWhiteRed(50),
               textMatrix = d0.PRIMARYTreatments.matrix,
               setStdMargins = FALSE,
               cex.text = 1,
               zlim = c(-0.6,0.6),
               main = paste("Module-trait relationships - Primary Treatment"))
dev.off()

# this heatmap looks better
d0.PRIMARYtreatment.text <-  as.matrix(signif(moduleTraitPvalue_Primary, 3))
pa = cluster::pam(d0.PRIMARYtreatment.text, k = 3)
col_fun = colorRamp2(c(-0.5, 0, 0.5), c("blue", "white", "red"))
png("Analysis/Output/WGCNA/subseq_treatments_all/Day0/day0_Treatments_Primary_heatmap.png", 500, 1000, pointsize=20)
pdf("Analysis/Output/WGCNA/subseq_treatments_all/Day0/heatmaps/day0_Treatments_Primary_heatmap.pdf", width=5, height=6)
Heatmap(moduleTraitCor_Primary, 
        name = "gene_cor", 
        rect_gp = gpar(col = "grey", lwd = 1),
        column_title = "Day 7 WGCNA - Primary Treatment", 
        column_title_gp = gpar(fontsize = 12, fontface = "bold"),
        # row_title = "WGCNA modules",
        #row_km = 4, 
        column_km = 1,
        row_split = paste0("clstr", pa$clustering),
        row_gap = unit(5, "mm"),
        column_gap = unit(5, "mm"),
        # grid.text(matrix(textMatrix)),
        # border = TRUE,
        border = TRUE,
        col = col_fun,
        cell_fun = function(j, i, x, y, width, height, fill) {
          grid.text(sprintf("%.1f", d0.PRIMARYtreatment.text[i, j]), x, y, gp = gpar(fontsize = 10))
        })
dev.off()





#=====================================================================================
#
# Module eigengene -  MEs boxplots by treatment group
#
#=====================================================================================
MEs_table <- MEs # new table for plotting 
MEs_table$Sample.Name <- rownames(MEs) # call rows as coolumn to merge with treatment data
MEsPlotting <- merge(d0.Treatment_data_2, MEs_table, by = 'Sample.Name') # merge
MEsPlotting <- MEsPlotting[,-c(2)] # ommit the phys data to just plot the module colors 
MEsPlotting_melt <- melt(MEsPlotting, id=c('Sample.Name', 'Primary_Treatment'))
#plot it
png("Analysis/Output/WGCNA/subseq_treatments_all/Day0/day0_ME_Boxplot.png", 600, 1000, pointsize=20)
ggplot(MEsPlotting_melt, aes(x=Primary_Treatment, y=value, fill = factor(Primary_Treatment), shape=Primary_Treatment)) +
  geom_boxplot(aes(middle = mean(value)), position=position_dodge(0.8), outlier.size = 0, alpha = 0.5) + 
  stat_summary(fun.y = mean, color = "black", position = position_dodge(0.75),
               geom = "point", shape = 19, size = 3,
               show.legend = FALSE) +
  ylab("ModuleEigengene") +
  ylim(-0.5,0.5) +
  scale_color_manual(values=c("#56B4E9","#D55E00")) +
  geom_hline(yintercept=0, linetype='dotted', col = 'black', size = 1)+
  theme_bw() +
  theme(legend.position = "none") +
  facet_wrap(~variable)
dev.off()

#=====================================================================================
#
# geneModuleMembership - geneTraitSignificance - GSPvalue
#
#=====================================================================================
# Gene relationship to trait and important modules: 
# Gene Significance and Module Membership

# We quantify associations of individual genes with our trait of interest (TAOC)

# names (colors) of the modules
modNames = substring(names(MEs), 3) # name all the modules, from 3rd character on (first two are ME)

geneModuleMembership = as.data.frame(cor(dds.d0_vst, MEs, use = "p"));
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));

names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");


# A treatment group
A = as.data.frame(d0.Traits.Primary$A); # Define variable containing the desired column 
names(A) = "A"
A_geneTraitSignificance = as.data.frame(cor(dds.d0_vst, A, use = "p"));
A_GSPvalue = as.data.frame(corPvalueStudent(as.matrix(A_geneTraitSignificance), nSamples));
names(A_geneTraitSignificance) = paste("GS.", names(A), sep=""); # MA_geneTraitSignificance - pearsons correlation between reads and the MA grop
names(A_GSPvalue) = paste("p.GS.", names(A), sep=""); # corPvalueStudent

# M treatment group
M = as.data.frame(d0.Traits.Primary$M); # Define variable containing the desired column 
names(M) = "M"
M_geneTraitSignificance = as.data.frame(cor(dds.d0_vst, M, use = "p"));
M_GSPvalue = as.data.frame(corPvalueStudent(as.matrix(M_geneTraitSignificance), nSamples));
names(M_geneTraitSignificance) = paste("GS.", names(M), sep=""); # MA_geneTraitSignificance - pearsons correlation between reads and the MA grop
names(M_GSPvalue) = paste("p.GS.", names(M), sep=""); # corPvalueStudent

#=====================================================================================
#
#   COUNT GENES OF INTEREST IN  MODULES (i.e. violet- refer to heatmap)
#
#=====================================================================================

length(colnames(dds.d0_vst)[moduleColors=="darkgrey"]) # 968 total genes in the violet module
length(colnames(dds.d0_vst)[moduleColors=="pink"]) # 169 total genes in the violet module

#=====================================================================================
#
#  Call annotation data to get module gene data (prep for downstream GO)
#
#=====================================================================================
annot = read.delim2(file = "Seq_details/Panopea-generosa-genes-annotations.txt",header = F);
dim(annot) # 34947     9
names(annot) # V1 are the gene.IDs
probes = names(d0.data)
probes2annot = match(probes, annot$V1)
# The following is the number or probes without annotation:
sum(is.na(probes2annot))
# Should return 0.
#=====================================================================================
#
#  BUILD GENE INFO DATAFRAMES
#
#=====================================================================================
# Create the starting data frame


#   TAOC dataframe - grey60 --------------------------------------------------------------------------- # 
names(A_geneTraitSignificance)
names(A_GSPvalue)
geneInfo_GROUPS = data.frame(substanceBXH = probes,
                                  geneSymbol = annot$V1[probes2annot],
                                  #LocusLinkID = annot$LocusLinkID[probes2annot],
                                  moduleColor = moduleColors,
                                  Uniprot = annot$V5[probes2annot],
                                  HGNC = annot$V6[probes2annot],
                                  GO.terms = annot$V8[probes2annot],
                                  GO.description = annot$V9[probes2annot],
                                  A_geneTraitSignificance, M_geneTraitSignificance, # call this specific to the module and trait of interest
                                  A_GSPvalue,  M_GSPvalue)              # call this specific to the module and trait of interest
View(geneInfo_GROUPS)
modOrder = order(-abs(cor(MEs, A, use = "p"))); # order by the strength of the correlation between module and trait values for each sample

for (mod in 1:ncol(geneModuleMembership)) # Add module membership information in the chosen order
{
  oldNames = names(geneInfo_GROUPS)
  geneInfo_GROUPS = data.frame(geneInfo_GROUPS, geneModuleMembership[, modOrder[mod]], 
                                    MMPvalue[, modOrder[mod]]);
  names(geneInfo_GROUPS) = c(oldNames, paste("A.", modNames[modOrder[mod]], sep=""),
                                  paste("p.A.", modNames[modOrder[mod]], sep=""))
}
# Order the genes in the geneInfo variable first by module color, then by geneTraitSignificance
geneOrder = order(geneInfo_GROUPS$moduleColor, -abs(geneInfo_GROUPS$GS.A));
geneInfo_GROUPS = geneInfo_GROUPS[geneOrder, ]
View(geneInfo_GROUPS)

#=====================================================================================
#
#  Write csv for the modules and corresponding raw read counts
#
#=====================================================================================
# call the module of interest for follow-up GO analysis 

write.csv(geneInfo_GROUPS, file = "Analysis/Output/WGCNA/subseq_treatments_all/Day0/d0.WGCNA_ModulMembership.csv")

#=====================================================================================
#
#  LOAD DATA for the following
#  Plot (1) Linear regressions of eigengenes for each module
#       (2) Eigengen exp by treatment  for modules of interest (review heatmaps)
#       (3) Go terms for modules of interest (review heatmaps)
#=====================================================================================

# SET WORKING DIRECTORY AND LOAD DATA
setwd("C:/Users/samjg/Documents/My_Projects/Pgenerosa_TagSeq_Metabolomics/TagSeq/")

# Load libraries 
library(dplyr)
library(goseq)
library(DESeq2)
library(reshape2)
library(ggplot2)
library(Rmisc)
library(ggpubr)
library(tibble)
library(hrbrthemes)
library(gridExtra)
library(ggplot2)
library(tidyr)
library(forcats) # for plotting later..
# Load data 
d0_ModEigen <- read.csv("Analysis/Output/WGCNA/subseq_treatments_all/Day0/d0.WGCNA_ModulEigengenes.csv")
d0_Annot_ModuleMembership <-  read.csv("Analysis/Output/WGCNA/subseq_treatments_all/Day0/d0.WGCNA_ModulMembership.csv")
d0_vst_data <- read.csv("Analysis/Output/WGCNA/subseq_treatments_all/Day0/day0_vstTransformed_WGCNAdata.csv")
# Master.Treatment_data <- read.csv(file="Analysis/Data/ Master_Phyenotype.and.Exp.Treatment_Metadata.csv", sep=',', header=TRUE)



#=====================================================================================
#
#  Prep and merge datasets for plotting
# 
#=====================================================================================
# Prep data to merge
# module membership
d0_Annot_ModuleMembership <- d0_Annot_ModuleMembership[-c(1:2)] # ommit the redundant gene name columns buut keep 'geneSymbol'
dim(d0_Annot_ModuleMembership) #   9207   22

# module eigengenes
names(d0_ModEigen)[1] <- "Sample.Name"

# normalized count data (same used for the dds object and WGCNA analysis)
d0_vst_data_t <- as.data.frame(t(d0_vst_data[, -(1)])) # trsnpose columns names to rows (genes) 
colnames(d0_vst_data_t) <- d0_vst_data[, 1] # name the columns as the rows in previous dataset (sample IDS)
d0_vst_data_t<- d0_vst_data_t %>% tibble::rownames_to_column("geneSymbol") # "geneSymbol" - create column and remname the gene name column to 'geneSymbol'


# merge Master datasets
# GO terms
d0_GOTermsMaster.Modules<-  merge(d0_Annot_ModuleMembership, d0_vst_data_t, by = "geneSymbol")
dim(d0_GOTermsMaster.Modules) # 9207   30
# Eigengenes and traits
d0_EigenTraitMaster<-  merge(d0_ModEigen, d0.Treatment_data_2, by = "Sample.Name") # d0.Treatment_data_2 created in this script 
dim(d0_EigenTraitMaster) #  8 9



#===================================================================================== 
# 
# EXPLORE THE Expression of each module (for loop plots!) BY TREATMENT
#
# evidence from the heatmaps and the regression shwo strong module membership and eigenenge cor strength with lighcyan 
# in response to treatment (primary specifically) and Total Antioxidant capacity 
#
#=====================================================================================

modcolor <- as.data.frame(unique(d0_Annot_ModuleMembership$moduleColor))
names(modcolor)[1] <- "color"

# experiment treatment and total protein data - narrow the columns 
d0.Treatment_data_2


for(i in 1:nrow(modcolor)) {
  
  # vst read count date - narrow the columns - reshape and rename
  Mod_geneIDs <- d0_Annot_ModuleMembership %>% dplyr::filter(moduleColor %in% modcolor[i,]) %>%  dplyr::select("geneSymbol")
  d0_vst_Mod <- d0_vst_data_t %>% dplyr::filter(geneSymbol %in% Mod_geneIDs[,1])
  d0_vst_Mod_MELT <- melt(d0_vst_Mod, id=("geneSymbol")) # melt using reshape2
  names(d0_vst_Mod_MELT)[(2:3)] <- c('Sample.Name', 'vst_Expression') # change column names
  
  # merge by common row values 'Sample.Name'
  merged_Expdata_Mod <- merge(d0_vst_Mod_MELT, d0.Treatment_data_2, by ='Sample.Name')
  
  # mean Exp response table 
  meanEXp_Mod <- merged_Expdata_Mod %>% 
    select(c('Sample.Name','vst_Expression','Primary_Treatment')) %>% 
    group_by(Sample.Name, Primary_Treatment) %>%
    dplyr::summarize(mean.vstExp = mean(vst_Expression), 
                     sd.vsdtExp = sd(vst_Expression),
                     na.rm=TRUE)
  

  # summarize datasets further by treatment period
  # remember:this is a mean of a mean!! First we complete mean vst exp by sample id (compiling all red module genes) - next all sample IDs by the treatment period (below
  # I will use these for mean SE plots 
  # Primary treatment
  meanEXp_Summary.Prim_Mod <- meanEXp_Mod %>% 
    group_by(Primary_Treatment) %>%
    dplyr::summarize(mean = mean(mean.vstExp), 
                     sd = sd(mean.vstExp),
                     n = n(), 
                     se = sd/sqrt(n))

  # PLOT
  # The errorbars overlapped, so use position_dodge to move them horizontally
  pd <- position_dodge(0.3) # move them .05 to the left and right
  
  # Primary treatment mean sd plot
  min_p1 <- min(meanEXp_Summary.Prim_Mod$mean) - max(meanEXp_Summary.Prim_Mod$se)
  max_p1 <- max(meanEXp_Summary.Prim_Mod$mean) + max(meanEXp_Summary.Prim_Mod$se)
  
  Primary.vst.Mod <- ggplot(meanEXp_Summary.Prim_Mod, aes(x=Primary_Treatment, y=mean, fill=Primary_Treatment)) +  # , colour=supp, group=supp))
    theme_classic() +
    geom_errorbar(aes(ymin=mean-se, ymax=mean+se), colour="black", width=.1, position=pd) +
    geom_line(position=pd) +
    geom_point(position=pd, size = 4, shape=21) +            
    xlab("Primary pCO2 treatment") +
    ylab(NULL) +                 # note the mean was first by sample ID THEN by treatment
    scale_fill_manual(values=c("#56B4E9","#E69F00")) +
    # scale_color_manual(values=c("#56B4E9","#E69F00")) +
    ggtitle(paste("Day 7 WGCNA", modcolor[i,], "Module VST GeneExp", sep =' ')) +
    # expand_limits(y=0) +                                                    # Expand y range
    scale_y_continuous(limits=c((min_p1), (max_p1))) +
    theme(text = element_text(size=15))
  
  # output 
  
  pdf(paste("Analysis/Output/WGCNA/subseq_treatments_all/Day0/ModuleExpression_Treatment/day0_Exp_Module",modcolor[i,],".pdf"), width=6, height=6)
  print(ggarrange(Primary.vst.Mod,         
                  plotlist = NULL,
                  ncol = 1,
                  nrow = 1,
                  labels = NULL))
  dev.off()
  
}

