---
  # title: "Day7_WGCNA"
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
day7.counts.matrix <- read.csv(file="Analysis/Data/filtered_counts/day7.counts.filtered_10cpm50perc.csv", sep=',', header=TRUE)
# Treatment and Phenotype data
Master.Treatment_Phenotype.data <- read.csv(file="Analysis/Data/ Master_Phyenotype.and.Exp.Treatment_Metadata.csv", sep=',', header=TRUE)
d7.Treatment_Phenotype.data     <- Master.Treatment_Phenotype.data %>%  dplyr ::filter(Date %in% 20190731) # split for day 7 data 


# The following setting is important, do not omit.
options(stringsAsFactors = FALSE)

# ===================================================================================
#
#
# Day 7 WGCNA - PREPROCESSING THE DATA INPUT 
#
#  Read here for the pre processing steps using WGCNA!
#  https://bmcbioinformatics.biomedcentral.com/track/pdf/10.1186/1471-2105-9-559.pdf
# ===================================================================================


# filter low values - note: this is pre filteres for < 5 CPM in 50% of samples - less strict than the DESeq2 analysis
dim(day7.counts.matrix) # 12013  rows (genes) -   37  samples  not counting 'x' and 'Gene.ID'
(day7.counts.matrix)[1] # ommit the first line and transpose in the next line 
d7.data = as.data.frame(t(day7.counts.matrix[, -(1)])) # ommit all columns but samples and transpose
dim(d7.data) # 36 8548
# fix(d7.data) # view the data


# trait data ========================================================== #

# Phenotype trait and Treatment data
dim(d7.Treatment_Phenotype.data) #  36 rows and 14 columns
names(d7.Treatment_Phenotype.data) # look at the 14 columns 
d7.Treatment_Phenotype.data = d7.Treatment_Phenotype.data[, -c(1:3,5,8)]; # remove columns that hold information we do not need. (i.e. tankks, Third treatment, All treatment group, etc.)
dim(d7.Treatment_Phenotype.data) # 36 rows and  9 columns
# fix(d7.Treatment_Phenotype.data)  # view the data


# count data  ========================================================== #


# fix(d21.data) # view the data - as you see all columns are now genes but filled as V1, V2, V3...
names(d7.data) = day7.counts.matrix$X # assigns column names (previous jsut numbered) as the gene ID 
rownames(d7.data) = names(day7.counts.matrix)[-(1)]; # assigns the row names as the sample ID
d7.data_matrix <- data.frame(day7.counts.matrix[,-1], row.names=day7.counts.matrix[,1]) 
d7.data_matrix_t <- t(d7.data_matrix)
# fix(d7.data_matrix_t)  # view the data

# create dds objects to transform data 
dds.d7 <- DESeqDataSetFromMatrix(countData = d7.data_matrix,
                                  colData = d7.Treatment_Phenotype.data, design = ~ 1) # DESeq Data Set (dds)
# DESeq Data Set (dds)
# NOTE: ~1 stands for no design; user will need to add a design for differential testing
# however for our purpose of just creating an onbject to transform, we do not need a design here...


# transform the data 
dds.d7_vst <- vst(dds.d7) # transform it vst
dds.d7_vst <- assay(dds.d7_vst) # call only the transformed coutns in the dds object
#fix(dds.d7_vst)
dds.d7_vst <- t(dds.d7_vst) # transpose columns to rows and vice versa

# ===================================================================================
#
# Day 7 WGCNA Sample tree 
#
# ===================================================================================

dim(dds.d7_vst) #  8548 genes; 36  samples

gsg = goodSamplesGenes(dds.d7_vst, verbose = 3);  # We first check for genes and samples with too many missing values:
gsg$allOK # If the statement returns TRUE, all genes have passed the cuts. 

# determine outlier samples 
sampleTree = hclust(dist(dds.d7_vst), method = "average") # Next we cluster the samples (in contrast to clustering genes that will come later)  to see if there are any obvious outliers.
sizeGrWindow(12,9) # The user should change the dimensions if the window is too large or too small.
par(cex = 0.6);
par(mar = c(0,4,2,0))

png("Analysis/Output/WGCNA/Day7/Day7_ClusterTree_Precut.png", 1000, 1000, pointsize=20)
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5, 
     cex.axis = 1.5, cex.main = 2) # appears there are two outliers SG55 and SG 105; can remove by hand or an automatic appraoch 
abline(h = 110, col = "red") # add line to plot to show the cut-off od outlier samples (40000) SG105 and SG55
dev.off()

# cut ethe tree and ommit from the reads data 
clust = cutreeStatic(sampleTree, cutHeight = 110, minSize = 10) # Determine cluster under the line
table(clust) # 0 = cut; 1 = kept; says it will cut 2 and save 60 
# 'keepsamples' boolean to call the main dataset
keepSamples = (clust==1) # call the 30 from clust at position 1
# view the ommited samples from the tree
omittedsamples = (clust==0)
omSamples = dds.d7_vst[omittedsamples,] 
nrow(omSamples) # should be 1
rownames(omSamples) # should be 6 - "SG90" "SG89" "SG62" "SG21" "SG22" "SG59" What treatments are these??
omSamplesData <- d7.Treatment_Phenotype.data %>% dplyr::filter(Sample.Name %in% (rownames(omSamples))) %>% dplyr::select(c('Sample.Name', 'Primary_Treatment', 'Second_Treament')) # three Primary Ambient and three Primary Moderate 
# integrate keepsamples to the main read dataset 
dds.d7_vst = dds.d7_vst[keepSamples, ] # integreat the boolean 'keepsamples' to ommit oultilers determined in the sample tree above 
nGenes = ncol(dds.d7_vst) # number of genes == 8548 
nSamples = nrow(dds.d7_vst) # number of samples == 35  - the cut tree removed 6 samples 



# plot the tree with the 'keep samples' bollean (T/F) to ommit outlier samples 
sampleTree2 = hclust(dist(dds.d7_vst), method = "average") # Next we cluster the samples (in contrast to clustering genes that will come later)  to see if there are any obvious outliers.
png("Analysis/Output/WGCNA/Day7/Day7_ClusterTree_Postcut.png", 1000, 1000, pointsize=20)
plot(sampleTree2, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5, 
     cex.axis = 1.5, cex.main = 2)
dev.off()


# based on outlier removall... call Trait data ===================================================== #
dim(dds.d7_vst) #   35 8548- transformed count data now  has 30 samples - 6 cut in the tree step above 
dim(d7.Treatment_Phenotype.data) # 36  9 - trait data has  36 samples - not yet cut! 

# Form a data frame analogous to expression data that will hold the clinical traits.
d7.Samples = rownames(dds.d7_vst);# start new variable 'd21.Samples' calling the row names of the gene data (sample as 'Mouse' in trait data)
TreatRows = match(d7.Samples, d7.Treatment_Phenotype.data$Sample.Name); # match the names 
d7.Traits = d7.Treatment_Phenotype.data[TreatRows, -1]; # removes the row numbers
rownames(d7.Traits) = d7.Treatment_Phenotype.data[TreatRows, 1]; # inserts the new traitRows - matches sample treatments
dim(d7.Traits) #  35 Samples 8 columns (primary and second treatment) - note: no Third treatment because that occured Days 14-21, after this sampling day
all(rownames(d7.Traits) == rownames(dds.d7_vst))  # should be TRUE
dim(d7.Traits) #  35  8


# ===================================================================================
#
# Prepare Trait data (Phys and Treatment groups)
# ===================================================================================

# ALL TRAITS 
d7.Traits
# ONLY TREATMENTS
d7.Traits.treat<- d7.Traits[,c(1:2)]
d7.Traits.treat$All <- paste (d7.Traits.treat$Primary_Treatment,d7.Traits.treat$Second_Treament, sep ='')
# primary groups 
d7.Traits.Primary <-  d7.Traits.treat %>% dplyr::select('Primary_Treatment')

d7.Traits.Primary$A <- (d7.Traits.Primary$Primary_Treatment == "A")
d7.Traits.Primary$A   <- as.numeric(d7.Traits.Primary$A)
d7.Traits.Primary$A <- as.factor(d7.Traits.Primary$A)

d7.Traits.Primary$M <- (d7.Traits.Primary$Primary_Treatment == "M")
d7.Traits.Primary$M   <- as.numeric(d7.Traits.Primary$M)
d7.Traits.Primary$M <- as.factor(d7.Traits.Primary$M)

d7.Traits.Primary <- d7.Traits.Primary[,c(2:3)] # final dataset of 0,1 for teatment groups - Primary only!

# Secondary groups 
d7.Traits.Secondary <-  d7.Traits.treat %>% dplyr::select('Second_Treament')

d7.Traits.Secondary$A <- (d7.Traits.Secondary$Second_Treament == "A")
d7.Traits.Secondary$A <- as.numeric(d7.Traits.Secondary$A)
d7.Traits.Secondary$A <- as.factor(d7.Traits.Secondary$A)

d7.Traits.Secondary$M <- (d7.Traits.Secondary$Second_Treament == "M")
d7.Traits.Secondary$M <- as.numeric(d7.Traits.Secondary$M)
d7.Traits.Secondary$M <- as.factor(d7.Traits.Secondary$M)

d7.Traits.Secondary$S <- (d7.Traits.Secondary$Second_Treament == "S")
d7.Traits.Secondary$S <- as.numeric(d7.Traits.Secondary$S)
d7.Traits.Secondary$S <- as.factor(d7.Traits.Secondary$S)

d7.Traits.Secondary <- d7.Traits.Secondary[,c(2:4)] # final dataset of 0,1 for teatment groups - Secondary only!

# Group
d7.Traits.Group <-  d7.Traits.treat %>% dplyr::select('All')
d7.Traits.Group$AA <- (d7.Traits.treat$All == "AA")
d7.Traits.Group$AA   <- as.numeric(d7.Traits.Group$AA)
d7.Traits.Group$AA <- as.factor(d7.Traits.Group$AA)

d7.Traits.Group$AM <- (d7.Traits.treat$All == "AM")
d7.Traits.Group$AM   <- as.numeric(d7.Traits.Group$AM)
d7.Traits.Group$AM <- as.factor(d7.Traits.Group$AM)

d7.Traits.Group$AS <- (d7.Traits.treat$All == "AS")
d7.Traits.Group$AS   <- as.numeric(d7.Traits.Group$AS)
d7.Traits.Group$AS <- as.factor(d7.Traits.Group$AS)

d7.Traits.Group$MA <- (d7.Traits.treat$All == "MA")
d7.Traits.Group$MA   <- as.numeric(d7.Traits.Group$MA)
d7.Traits.Group$MA <- as.factor(d7.Traits.Group$MA)

d7.Traits.Group$MM <- (d7.Traits.treat$All == "MM")
d7.Traits.Group$MM   <- as.numeric(d7.Traits.Group$MM)
d7.Traits.Group$MM <- as.factor(d7.Traits.Group$MM)

d7.Traits.Group$MS <- (d7.Traits.treat$All == "MS")
d7.Traits.Group$MS   <- as.numeric(d7.Traits.Group$MS)
d7.Traits.Group$MS <- as.factor(d7.Traits.Group$MS)

d7.Traits.Group <- d7.Traits.Group[,c(2:7)] # final dataset of 0,1 for teatment groups - Secondary groups only!


                         
# ONLY PHYSIOLOGY
d7.Traits.phys<- d7.Traits[,c(3:(ncol(d7.Traits)))] # pphysiologial measuremetns from the same replicate tanks


# ===================================================================================
#
# cluster samples by Trait
# ===================================================================================

# All Traits
png("Analysis/Output/WGCNA/Day7/Day7_ClusterTree_AllTraits.png", 1000, 1000, pointsize=20)
sampleTree2 = hclust(dist(dds.d7_vst), method = "average")
traitColors = labels2colors(d7.Traits); # Convert traits to a color representation: white means low, red means high, grey means missing entry
plotDendroAndColors(sampleTree2, traitColors, # Plot the sample dendrogram and the colors underneath.
                    groupLabels = names(d7.Traits), 
                    main = "Sample dendrogram and trait heatmap")
dev.off()

# Primary treatment Only
png("Analysis/Output/WGCNA/Day7/Day7_ClusterTree_PrimaryTreatment.png", 1000, 1000, pointsize=20)
traitColors_Primary = labels2colors(d7.Traits.Primary); # Convert traits to a color representation: white means low, red means high, grey means missing entry
plotDendroAndColors(sampleTree2, traitColors_Primary, # Plot the sample dendrogram and the colors underneath.
                    groupLabels = names(d7.Traits.Primary), 
                    main = "Sample dendrogram and trait heatmap")
dev.off()

# Second treatment 
png("Analysis/Output/WGCNA/Day7/Day7_ClusterTree_SecondTreatment.png", 1000, 1000, pointsize=20)
traitColors_Second = labels2colors(d7.Traits.Secondary); # Convert traits to a color representation: white means low, red means high, grey means missing entry
plotDendroAndColors(sampleTree2, traitColors_Second, # Plot the sample dendrogram and the colors underneath.
                    groupLabels = names(d7.Traits.Secondary), 
                    main = "Sample dendrogram and trait heatmap")
dev.off()

# gROUP treatment (primary_second)
png("Analysis/Output/WGCNA/Day7/Day7_ClusterTree_GroupTreatment.png", 1000, 1000, pointsize=20)
traitColorsGroup = labels2colors(d7.Traits.Group); # Convert traits to a color representation: white means low, red means high, grey means missing entry
plotDendroAndColors(sampleTree2, traitColorsGroup, # Plot the sample dendrogram and the colors underneath.
                    groupLabels = names(d7.Traits.Group), 
                    main = "Sample dendrogram and trait heatmap")
dev.off()


# Phys Only
png("Analysis/Output/WGCNA/Day7/Day7_ClusterTree_PhysOnly.png", 1000, 1000, pointsize=20)
traitColors_phys = labels2colors(d7.Traits.phys); # Convert traits to a color representation: white means low, red means high, grey means missing entry
plotDendroAndColors(sampleTree2, traitColors_phys, # Plot the sample dendrogram and the colors underneath.
                    groupLabels = names(d7.Traits.phys), 
                    main = "Sample dendrogram and trait heatmap")
dev.off()


# save data
save(dds.d7_vst, d7.Traits, file = "Analysis/Output/WGCNA/Day7/d.7-dataInput.RData")
save(dds.d7_vst, d7.Traits.phys, file = "Analysis/Output/WGCNA/Day7/d.7-dataInput_treatmentonly.RData")
save(dds.d7_vst, d7.Traits.treat, file = "Analysis/Output/WGCNA/Day7/d.7-dataInput_physonly.RData")


# write the vst transformed data 
write.csv(dds.d7_vst, "Analysis/Output/WGCNA/Day7/Day7_vstTransformed_WGCNAdata.csv") # write

#=====================================================================================
#
#  Satrt the step-wise module construction:  
# Step 1 = create adjacency matrix 
#
#=====================================================================================
softPower = 3 # set your soft threshold based on the plots above 
adjacency = adjacency(dds.d7_vst, power = softPower) # this takes a long time.. just wait...
#=====================================================================================
#
#  Step 2: Turn adjacency into topological overlap
#
#=====================================================================================
TOM = TOMsimilarity(adjacency)  # this takes a long time.. just wait...
dissTOM = 1-TOM
#=====================================================================================
#
#  Step 3:Call the hierarchical clustering function - plot the tree
#
#=====================================================================================
# Call the hierarchical clustering function
geneTree = hclust(as.dist(dissTOM), method = "average");
# Plot the resulting clustering tree (dendrogram)
sizeGrWindow(12,9)
plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity",
     labels = FALSE, hang = 0.04);
#=====================================================================================
#
#  Step 4: Set module size and 'cutreeDynamic' to create clusters 
#
#=====================================================================================
# We like large modules, so we set the minimum module size relatively high:
minModuleSize = 100; # set this for the subseqent call...
# Module identification using dynamic tree cut:
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
                            deepSplit = 1, pamRespectsDendro = FALSE,
                            minClusterSize = minModuleSize);
table(dynamicMods) # number of genes per module 
#=====================================================================================
#
#  Step 5: convert numeric network to colors and plot the dendrogram
#
#=====================================================================================
# Convert numeric lables into colors
dynamicColors = labels2colors(dynamicMods) # add colors to module labels (previously numbers)
table(dynamicColors) # lets look at this table...
# Plot the dendrogram and colors underneath
sizeGrWindow(8,6)
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram and module colors")
#=====================================================================================
#
#  Step 6: Calculate Eigengenes - view thier connectivity based on 'MEDiss = 1-cor(MEs)'
#
#=====================================================================================
# Calculate eigengenes
MEList = moduleEigengenes(dds.d7_vst, colors = dynamicColors)
MEs = MEList$eigengenes
# Calculate dissimilarity of module eigengenes
MEDiss = 1-cor(MEs);
# Cluster module eigengenes
METree = hclust(as.dist(MEDiss), method = "average");
# Plot the result
sizeGrWindow(7, 6)
png("Analysis/Output/WGCNA/Day7/Day7_ClusterEigengenes.png", 1000, 1000, pointsize=20)
plot(METree, main = "Clustering of module eigengenes (dissimilarity calc = MEDiss = 1-cor(MEs))",
     xlab = "", sub = "")
dev.off()
#=====================================================================================
#
#  Step 7: Specify the cut line for the dendrogram (module) - Calc MODULE EIGENGENES (mergeMEs)
#
#=====================================================================================
MEDissThres = 0.35
# Plot the cut line into the dendrogram
abline(h=MEDissThres, col = "red")
# Call an automatic merging function
merge = mergeCloseModules(dds.d7_vst, dynamicColors, cutHeight = MEDissThres, verbose = 3)
# The merged module colors
mergedColors = merge$colors;
# Eigengenes of the new merged modules:
mergedMEs = merge$newMEs;
#=====================================================================================
#
#  Step 8: Plot dendrogram with the cut line 'MEDissThres' 
#
#=====================================================================================
sizeGrWindow(12, 9)
png("Analysis/Output/WGCNA/Day7/Day7_ClusterDendrogram.png", 1000, 1000, pointsize=20)
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors),
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
save(MEs, moduleLabels, moduleColors, geneTree, file = "Analysis/Output/WGCNA/Day7/Day7-networkConstruction-stepByStep.RData")
# write csv - save the module eigengenes
write.csv(MEs, file = "Analysis/Output/WGCNA/Day7/d7.WGCNA_ModulEigengenes.csv")
table(mergedColors)

#=====================================================================================
#
#  Constructing the gene network and identifying modules
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

# https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/TechnicalReports/signedTOM.pdf
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

net = blockwiseModules(dds.d7_vst, #maxBlockSize = 2000,
                       power = 3, 
                    #  minModuleSize = 30, 
                    #  deepSplit=4,
                       mergeCutHeight = 0.75, #reassignThreshold = 0, 
                       numericLabels = TRUE,
                     #  pamRespectsDendro = FALSE,
                       TOMType = "signed",  
                       saveTOMs = TRUE,
                       saveTOMFileBase = "Analysis/Output/WGCNA/Day7/d7.Geoduck", 
                       verbose = 3) #... wait for this to finish 
moduleColors = labels2colors(net$colors) # assign colors to eigengene





# Clustering tree based on eigengene modules - How related are the module eigengenes?
# represent modules by eigengenes and relating them to one another 
# dendrogram showing the 
datME=moduleEigengenes(dds.d7_vst,moduleColors)$eigengenes
signif(cor(datME, use="p"), 2)
dissimME=(1-t(cor(datME, method="p")))/2 # define a dissimilarity measure between module eigengenes that keeps track of the sign of the correlation and uses this to cluster 
hclustdatME=hclust(as.dist(dissimME), method="average" ) # dendrogram based on the module eigenge connecticity 
# Plot the eigengene dendrogram

# Hub gene of ech module 
# Both functions output a character vector of genes, where the genes are the hub gene picked for each module, 
# \and the names correspond to the module in which each gene is a hub.
chooseTopHubInEachModule(dds.d7_vst, moduleColors, omitColors = "grey", power = 2, type = "signed")

# 'minKMEtoStay' genes whose eigengene connectivity to their module eigengene is lower than minKMEtoStay are removed from the module.
# 'minGap' minimum cluster gap given as the fraction of the difference between cutHeight and the 5th percentile of joining heights. See cutreeDynamic for more details.


# net$colors # the module assignment
# net$MEs # contains the module eigengenes of the modules.
table(net$colors)
# there are 19  modules in order of decending size
# note 0 is reserved for genes outside of all modules (2602 genes) 

#=====================================================================================
#
#  The hierarchical clustering dendrogram (tree) used for the module identification 
#
#=====================================================================================

# The hierarchical clustering dendrogram (tree) used for the module identification 
# open a graphics window
sizeGrWindow(12, 9)
# Convert labels to colors for plotting
moduleColors = labels2colors(net$colors)
# Plot the dendrogram and the module colors underneath for block 1
png("Analysis/Output/WGCNA/Day7/Day7_ClusterDendrogram.png", 1000, 1000, pointsize=20)
plotDendroAndColors(net$dendrograms[[1]], moduleColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05, abHeight = 0.98)
dev.off()
# Plot the dendrogram and the module colors underneath for block 2
plotDendroAndColors(net$dendrograms[[2]], moduleColors[net$blockGenes[[2]]],
                    "Module colors", main = "Gene dendrogram and module colors in block 2",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05, abHeight = 0.99)
# Plot the dendrogram and the module colors underneath for block 2
plotDendroAndColors(net$dendrograms[[3]], moduleColors[net$blockGenes[[3]]],
                    "Module colors", main = "Gene dendrogram and module colors in block 3",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05, abHeight = 0.99)



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
nGenes = ncol(dds.d7_vst); # 8548
nSamples = nrow(dds.d7_vst); # 35
# Recalculate MEs with color labels
MEs0 = moduleEigengenes(dds.d7_vst, moduleColors)$eigengenes 
MEs = orderMEs(MEs0) # reorders the columns (colors/modules)

# write csv - save the module eigengenes
write.csv(MEs, file = "Analysis/Output/WGCNA/Day7/d7.WGCNA_ModulEigengenes.csv")

# change chanracter treatments to integers
# ALL TRAIT DATA
d7.Traits$Primary_Treatment <- as.factor(d7.Traits$Primary_Treatment)
d7.Traits$Primary_Treatment <- as.numeric(d7.Traits$Primary_Treatment)

d7.Traits$Second_Treament <- as.factor(d7.Traits$Second_Treament)
d7.Traits$Second_Treament <- as.numeric(d7.Traits$Second_Treament)

# TREATMENT ONLY - break this up into more groups!! 
d7.Traits.treat$Primary_Treatment <- as.factor(d7.Traits.treat$Primary_Treatment)
d7.Traits.treat$Primary_Treatment <- as.numeric(d7.Traits.treat$Primary_Treatment)

d7.Traits.treat$Second_Treament <- as.factor(d7.Traits.treat$Second_Treament)
d7.Traits.treat$Second_Treament <- as.numeric(d7.Traits.treat$Second_Treament)

d7.Traits.treat$All <- as.factor(d7.Traits.treat$All)
d7.Traits.treat$All <- as.numeric(d7.Traits.treat$All)

#=====================================================================================
#
# Module trait correlation
#
#=====================================================================================
# ALL TRAIT DATA

dim(d7.Traits)  # 35  8
dim(MEs)  # 35  6
# moduleTraitCor = cor(MEs, d7.Traits, use = "p");
# moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);

# primary 
moduleTraitCor_Primary = cor(MEs, d7.Traits.Primary, use = "p");
moduleTraitPvalue_Primary = corPvalueStudent(moduleTraitCor_Primary, nSamples);
# Second
moduleTraitCor_Secondary = cor(MEs, d7.Traits.Secondary, use = "p");
moduleTraitPvalue_Secondary = corPvalueStudent(moduleTraitCor_Secondary, nSamples);
# Group
moduleTraitCor_Group = cor(MEs, d7.Traits.Group, use = "p");
moduleTraitPvalue_Group = corPvalueStudent(moduleTraitCor_Group, nSamples);

# # PHYSIOLOGY  ONLY 
# moduleTraitCor_physonly = cor(MEs, d7.Traits.phys, use = "p");
# moduleTraitPvalue_physonly = corPvalueStudent(moduleTraitCor_physonly, nSamples);

#=====================================================================================
#
# Heatmaps
#
#=====================================================================================



# PRRIMARY TRETMENT ONLY  ------------------------------------------------------------------ # 




sizeGrWindow(10,10)
# Will display correlations and their p-values
d7.PRIMARYTreatments.matrix <-  paste(signif(moduleTraitCor_Primary, 2), "\n(",
                                       signif(moduleTraitPvalue_Primary, 1), ")", sep = "")
#dim(textMatrix) == dim(moduleTraitCor_treatonly)
par(mar = c(8, 9.5, 5, 3));
# Display the correlation values within a heatmap plot
png("Analysis/Output/WGCNA/Day7/heatmaps/Day7_Treatments_Primary_heatmap2.png", 500, 1000, pointsize=20)
labeledHeatmap(Matrix = moduleTraitCor_Primary,
               xLabels = names(d7.Traits.Primary),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = TRUE,
               colors = blueWhiteRed(50),
               textMatrix = d7.PRIMARYTreatments.matrix,
               setStdMargins = FALSE,
               cex.text = 1,
               zlim = c(-0.6,0.6),
               main = paste("Module-trait relationships - Primary Treatment"))
dev.off()

# this heatmap looks better
d7.PRIMARYtreatment.text <-  as.matrix(signif(moduleTraitPvalue_Primary, 3))
pa = cluster::pam(d7.PRIMARYtreatment.text, k = 3)
col_fun = colorRamp2(c(-0.5, 0, 0.5), c("blue", "white", "red"))
png("Analysis/Output/WGCNA/Day7/Day7_Treatments_Primary_heatmap.png", 500, 1000, pointsize=20)
pdf("Analysis/Output/WGCNA/Day7/heatmaps/Day7_Treatments_Primary_heatmap.pdf", width=5, height=6)
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
          grid.text(sprintf("%.1f", d7.PRIMARYtreatment.text[i, j]), x, y, gp = gpar(fontsize = 10))
        })
dev.off()






# SECOND TREATMENT ONLY ------------------------------------------------------------------------------------- # 






sizeGrWindow(10,10)
# Will display correlations and their p-values
d7.SECONDTreatments.matrix <-  paste(signif(moduleTraitCor_Secondary, 2), "\n(",
                                      signif(moduleTraitPvalue_Secondary, 3), ")", sep = "")
par(mar = c(8, 9.5, 5, 3));
# Display the correlation values within a heatmap plot
png("Analysis/Output/WGCNA/Day7/heatmaps/Day7_Treatments_Second_heatmap2.png", 1000, 500, pointsize=20)
labeledHeatmap(Matrix = moduleTraitCor_Secondary,
               xLabels = names(d7.Traits.Secondary),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = TRUE,
               colors = blueWhiteRed(50),
               textMatrix = d7.SECONDTreatments.matrix,
               setStdMargins = FALSE,
               cex.text = 0.8,
               zlim = c(-0.6,0.6),
               main = paste("Module-trait relationships - Treatment Groups"))
dev.off()

# this heatmap looks better
d7.SECONDtreatment_text<-  as.matrix(signif(moduleTraitPvalue_Secondary, 4))
d7.SECONDtreatment_cor <-  as.matrix(signif(moduleTraitCor_Secondary, 4))
pa = cluster::pam(d7.SECONDtreatment_cor, k = 3)
col_fun = colorRamp2(c(-0.5, 0, 0.5), c("blue", "white", "red"))

png("Analysis/Output/WGCNA/Day7/heatmaps/Day7_Treatments_Second_heatmap.png", 500, 1000, pointsize=20)
pdf("Analysis/Output/WGCNA/Day7/heatmaps/Day7_Treatments_Second_heatmap.pdf", width=5, height=6)
Heatmap(moduleTraitCor_Secondary, 
        name = "gene_cor", 
        rect_gp = gpar(col = "grey", lwd = 1),
        column_title = "Day 7 WGCNA - Treatment Groups", 
        column_title_gp = gpar(fontsize = 12, fontface = "bold"),
        # row_title = "WGCNA modules",
        #row_km = 4, 
        column_km = 2,
        row_split = paste0("clstr", pa$clustering),
        row_gap = unit(5, "mm"),
        column_gap = unit(5, "mm"),
        # grid.text(matrix(textMatrix)),
        # border = TRUE,
        border = TRUE,
        col = col_fun,
        cell_fun = function(j, i, x, y, width, height, fill) {
          grid.text(sprintf("%.1f", d7.SECONDtreatment_text[i, j]), x, y, gp = gpar(fontsize = 10))
        })
dev.off()








# GROUPS ONLY ------------------------------------------------------------------------------------- # 






sizeGrWindow(10,10)
# Will display correlations and their p-values
d7.GROUPTreatments.matrix <-  paste(signif(moduleTraitCor_Group, 2), "\n(",
                                     signif(moduleTraitPvalue_Group, 3), ")", sep = "")
par(mar = c(8, 9.5, 5, 3));
# Display the correlation values within a heatmap plot
png("Analysis/Output/WGCNA/Day7/heatmaps/Day7_Group_heatmap2.png", 1000, 500, pointsize=20)
labeledHeatmap(Matrix = moduleTraitCor_Group,
               xLabels = names(d7.Traits.Group),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = TRUE,
               colors = blueWhiteRed(50),
               textMatrix = d7.GROUPTreatments.matrix,
               setStdMargins = FALSE,
               cex.text = 0.8,
               zlim = c(-0.6,0.6),
               main = paste("Module-trait relationships - Treatment Groups"))
dev.off()

# this heatmap looks better
d7.GROUPtreatment_text<-  as.matrix(signif(moduleTraitPvalue_Group, 4))
d7.GROUPtreatment_cor <-  as.matrix(signif(moduleTraitCor_Group, 4))
pa = cluster::pam(d7.GROUPtreatment_cor, k = 2)
col_fun = colorRamp2(c(-0.5, 0, 0.5), c("blue", "white", "red"))

png("Analysis/Output/WGCNA/Day7/heatmaps/Day7_Group_heatmap.png", 500, 1000, pointsize=20)
pdf("Analysis/Output/WGCNA/Day7/heatmaps/Day7_Group_heatmap.pdf", width=5, height=6)
Heatmap(moduleTraitCor_Group, 
        name = "gene_cor", 
        rect_gp = gpar(col = "grey", lwd = 1),
        column_title = "Day 7 WGCNA - Treatment Groups", 
        column_title_gp = gpar(fontsize = 12, fontface = "bold"),
        # row_title = "WGCNA modules",
        #row_km = 4, 
        column_km = 3,
        row_split = paste0("clstr", pa$clustering),
        row_gap = unit(5, "mm"),
        column_gap = unit(5, "mm"),
        # grid.text(matrix(textMatrix)),
        # border = TRUE,
        border = TRUE,
        col = col_fun,
        cell_fun = function(j, i, x, y, width, height, fill) {
          grid.text(sprintf("%.1f", d7.GROUPtreatment_text[i, j]), x, y, gp = gpar(fontsize = 10))
        })
dev.off()


#=====================================================================================
#
# Module eigengene -  MEs boxplots by treatment group
#
#=====================================================================================
MEs_table <- MEs # new table for plotting 
MEs_table$Sample.Name <- rownames(MEs) # call rows as coolumn to merge with treatment data
MEsPlotting <- merge(d7.Treatment_Phenotype.data, MEs_table, by = 'Sample.Name') # merge
MEsPlotting <- MEsPlotting[,-c(4:9)] # ommit the phys data to just plot the module colors 
MEsPlotting_melt <- melt(MEsPlotting, id=c('Sample.Name', 'Primary_Treatment', 'Second_Treament'))
#plot it
png("Analysis/Output/WGCNA/Day7/Day7_ME_Boxplot.png", 600, 1000, pointsize=20)
ggplot(MEsPlotting_melt, aes(x=Second_Treament, y=value, fill = factor(Primary_Treatment), shape=Primary_Treatment)) +
  geom_boxplot(aes(middle = mean(value)), position=position_dodge(0.8), outlier.size = 0, alpha = 0.5) + 
  stat_summary(fun.y = mean, color = "black", position = position_dodge(0.75),
               geom = "point", shape = 19, size = 3,
               show.legend = FALSE) +
  ylab("ModuleEigengene") +
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

geneModuleMembership = as.data.frame(cor(dds.d7_vst, MEs, use = "p"));
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));

names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");


# AA treatment group
AA = as.data.frame(d7.Traits.Secondary$AA); # Define variable containing the desired column 
names(AA) = "AA"
AA_geneTraitSignificance = as.data.frame(cor(dds.d7_vst, AA, use = "p"));
AA_GSPvalue = as.data.frame(corPvalueStudent(as.matrix(AA_geneTraitSignificance), nSamples));
names(AA_geneTraitSignificance) = paste("GS.", names(AA), sep=""); # MA_geneTraitSignificance - pearsons correlation between reads and the MA grop
names(AA_GSPvalue) = paste("p.GS.", names(AA), sep=""); # corPvalueStudent

# AM treatment group
AM = as.data.frame(d7.Traits.Secondary$AM); # Define variable containing the desired column 
names(AM) = "AM"
AM_geneTraitSignificance = as.data.frame(cor(dds.d7_vst, AM, use = "p"));
AM_GSPvalue = as.data.frame(corPvalueStudent(as.matrix(AM_geneTraitSignificance), nSamples));
names(AM_geneTraitSignificance) = paste("GS.", names(AM), sep=""); # MA_geneTraitSignificance - pearsons correlation between reads and the MA grop
names(AM_GSPvalue) = paste("p.GS.", names(AM), sep=""); # corPvalueStudent

# AS treatment group
AS = as.data.frame(d7.Traits.Secondary$AS); # Define variable containing the desired column 
names(AS) = "AS"
AS_geneTraitSignificance = as.data.frame(cor(dds.d7_vst, AS, use = "p"));
AS_GSPvalue = as.data.frame(corPvalueStudent(as.matrix(AS_geneTraitSignificance), nSamples));
names(AS_geneTraitSignificance) = paste("GS.", names(AS), sep=""); # MA_geneTraitSignificance - pearsons correlation between reads and the MA grop
names(AS_GSPvalue) = paste("p.GS.", names(AS), sep=""); # corPvalueStudent

# MA treatment group
MA = as.data.frame(d7.Traits.Secondary$MA); # Define variable containing the desired column 
names(MA) = "MA"
rownames(dds.d7_vst)
MA_geneTraitSignificance = as.data.frame(cor(dds.d7_vst, MA, use = "p"));
MA_GSPvalue = as.data.frame(corPvalueStudent(as.matrix(MA_geneTraitSignificance), nSamples));
names(MA_geneTraitSignificance) = paste("GS.", names(MA), sep=""); # MA_geneTraitSignificance - pearsons correlation between reads and the MA grop
names(MA_GSPvalue) = paste("p.GS.", names(MA), sep=""); # corPvalueStudent

# MM treatment group
MM = as.data.frame(d7.Traits.Secondary$MM); # Define variable containing the desired column 
names(MM) = "MM"
MM_geneTraitSignificance = as.data.frame(cor(dds.d7_vst, MM, use = "p"));
MM_GSPvalue = as.data.frame(corPvalueStudent(as.matrix(MM_geneTraitSignificance), nSamples));
names(MM_geneTraitSignificance) = paste("GS.", names(MM), sep=""); # MA_geneTraitSignificance - pearsons correlation between reads and the MA grop
names(MM_GSPvalue) = paste("p.GS.", names(MM), sep=""); # corPvalueStudent

# MS treatment group
MS = as.data.frame(d7.Traits.Secondary$MS); # Define variable containing the desired column 
names(MS) = "MS"
MS_geneTraitSignificance = as.data.frame(cor(dds.d7_vst, MS, use = "p"));
MS_GSPvalue = as.data.frame(corPvalueStudent(as.matrix(MS_geneTraitSignificance), nSamples));
names(MS_geneTraitSignificance) = paste("GS.", names(MS), sep=""); # MA_geneTraitSignificance - pearsons correlation between reads and the MA grop
names(MS_GSPvalue) = paste("p.GS.", names(MS), sep=""); # corPvalueStudent

#  PLOT mean.µmol.CRE.g.protein in the MAGENTA module
# unique(moduleColors)
# module = "brown"
# column = match(module, modNames);
# moduleGenes = moduleColors==module;
# 
# sizeGrWindow(7, 7);
# par(mfrow = c(1,1));
# verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
#                    abs(MA_geneTraitSignificance[moduleGenes, 1]),
#                    xlab = paste("Module Membership in", module, "module"),
#                    ylab = "Gene significance for 'MA' treatment group",
#                    main = paste("Day7 'MA' treatment group  'BROWN': Module membership vs. gene significance\n"),
#                    cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = "black")

#=====================================================================================
#
#   COUNT GENES OF INTEREST IN  MODULES (i.e. violet- refer to heatmap)
#
#=====================================================================================

length(colnames(dds.d7_vst)[moduleColors=="brown"]) # 968 total genes in the violet module
length(colnames(dds.d7_vst)[moduleColors=="black"]) # 169 total genes in the violet module

#=====================================================================================
#
#  Call annotation data to get module gene data (prep for downstream GO)
#
#=====================================================================================
annot = read.delim2(file = "Seq_details/Panopea-generosa-genes-annotations.txt",header = F);
dim(annot) # 34947     9
names(annot) # V1 are the gene.IDs
probes = names(d7.data)
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
names(MA_geneTraitSignificance)
names(MA_GSPvalue)
geneInfo_GROUPS = data.frame(substanceBXH = probes,
                                  geneSymbol = annot$V1[probes2annot],
                                  #LocusLinkID = annot$LocusLinkID[probes2annot],
                                  moduleColor = moduleColors,
                                  Uniprot = annot$V5[probes2annot],
                                  HGNC = annot$V6[probes2annot],
                                  GO.terms = annot$V8[probes2annot],
                                  GO.description = annot$V9[probes2annot],
                                  AA_geneTraitSignificance, AM_geneTraitSignificance, AS_geneTraitSignificance,
                                  MA_geneTraitSignificance, MM_geneTraitSignificance, MS_geneTraitSignificance, # call this specific to the module and trait of interest
                                  AA_GSPvalue,  AM_GSPvalue,  AS_GSPvalue,
                                  MA_GSPvalue,  MM_GSPvalue,  MS_GSPvalue)              # call this specific to the module and trait of interest
View(geneInfo_GROUPS)
modOrder = order(-abs(cor(MEs, MA, use = "p"))); # order by the strength of the correlation between module and trait values for each sample

for (mod in 1:ncol(geneModuleMembership)) # Add module membership information in the chosen order
{
  oldNames = names(geneInfo_GROUPS)
  geneInfo_GROUPS = data.frame(geneInfo_GROUPS, geneModuleMembership[, modOrder[mod]], 
                                    MMPvalue[, modOrder[mod]]);
  names(geneInfo_GROUPS) = c(oldNames, paste("MM.", modNames[modOrder[mod]], sep=""),
                                  paste("p.MM.", modNames[modOrder[mod]], sep=""))
}
# Order the genes in the geneInfo variable first by module color, then by geneTraitSignificance
geneOrder = order(geneInfo_GROUPS$moduleColor, -abs(geneInfo_GROUPS$GS.MA));
geneInfo_GROUPS = geneInfo_GROUPS[geneOrder, ]
View(geneInfo_GROUPS)



#=====================================================================================
#
#  PLOTTING LINE GRAPHS BY MODULE COLOR 
#
#=====================================================================================

# view the Treatment groups heatmap to understand rationale for clusters 

# SEPARATE MODULES 

# red module - higher in AS than MM and AM 
length(colnames(dds.d7_vst)[moduleColors=="red"]) # 188 total genes in the brown module
red_mod <- geneInfo_GROUPS %>%  dplyr::filter(moduleColor %in% 'red')
# red <- red_mod[,c(3,4:7)]
# red$Gene.ID  <- rownames(red)
# ExpDataD7 <- as.data.frame(dds.d7_vst)
# ExpDataD7 <- tibble::rownames_to_column(ExpDataD7, "Sample.Name")
# # merge and melt 
# d7.TreatData <-d7.Treatment_Phenotype.data[,c(1:3)] # call treatment and sample ID
# redMAster <- merge(d7.TreatData, ExpDataD7, by ='Sample.Name')  # merge by sample ID
# redMAster_melt <- melt(redMAster,  id=c("Sample.Name", "Primary_Treatment", "Second_Treament"))
# head(redMAster_melt)
# redMAster_melt2 <- redMAster_melt %>% dplyr::filter(variable %in% (rownames(red)))
# 
# # summarize it 
# redMAster_2 <- redMAster_melt2 %>% 
#   group_by(variable, Primary_Treatment, Second_Treament) %>%
#   dplyr::summarize(n = n(), 
#                    mean.vstExp = mean(value), 
#                    sd.vsdtExp = sd(value),
#                    se.vsdtExp = sd.vsdtExp/sqrt(n)) # mean expression for each sample ID 
# unique(redMAster_2$variable)
# # graph it
# head(redMAster_2)
# ggplot(redMAster_2, aes(x=Second_Treament, y=mean.vstExp, group=variable, color = Primary_Treatment)) +
#   theme(panel.grid=element_blank()) +
#   scale_color_manual(values=c("#56B4E9","#D55E00")) +
#   geom_line(size=0.2, alpha=0.8) +
#   theme_bw() +
#   facet_wrap(~Primary_Treatment)
# 
# ggplot(redMAster_2, aes(x=Second_Treament, y=mean.vstExp, group=variable, color = Primary_Treatment)) +
#   theme(panel.grid=element_blank()) +
#   scale_color_manual(values=c("#56B4E9","#D55E00")) +
#   geom_line(size=0.1, alpha=0.08) +
#   theme_bw() +
#   ggtitle("Day 14 WGCNA: Brown module") +
#   ylab("Pearson's correlation of vst Expression and Treatment Group (genes with Module Membership < 0.0001)") +
#   facet_wrap(~Primary_Treatment)

red_mod_05 <- red_mod %>%  dplyr::filter(p.MM.red < 0.001) 
red_mod_05_MM <-red_mod_05[,c(8:13)]
nrow(red_mod_05_MM) # 78 total genes in the red module at 0.0001 module membership 

red_mod_05_MM$Gene.ID <- rownames(red_mod_05_MM)
red_mod_05_MM_melt <- melt(red_mod_05_MM, id=c("Gene.ID"))
red_mod_05_MM_melt$Primary <- substr(red_mod_05_MM_melt$variable, 4,4)
red_mod_05_MM_melt$Second <- substr(red_mod_05_MM_melt$variable, 5,5)

ggplot(red_mod_05_MM_melt, aes(x=Second, y=value, group=Gene.ID, color = Primary)) +
  theme(panel.grid=element_blank()) +
  scale_color_manual(values=c("#56B4E9","#D55E00")) +
  geom_line(size=0.1, alpha=0.2) +
  theme_bw() +
  ggtitle("Day 7 WGCNA: red module") +
  ylab("Pearson's correlation of vst Expression and Treatment Group (genes with Module Membership < 0.001)") +
  facet_wrap(~Primary)


# green module - higher in AS than other modules (oppostie of black module) 
length(colnames(dds.d7_vst)[moduleColors=="green"]) # 595 total genes in the green module
green_mod <- geneInfo_GROUPS %>%  dplyr::filter(moduleColor %in% 'green')
nrow(green_mod) # 595 total genes in the green module
green_mod_05 <- green_mod %>%  dplyr::filter(p.MM.green < 0.001) 
green_mod_05_MM <-green_mod_05[,c(8:13)]
nrow(green_mod_05_MM) # 342 total genes in the green module at 0.0001 module membership 

green_mod_05_MM$Gene.ID <- rownames(green_mod_05_MM)
green_mod_05_MM_melt <- melt(green_mod_05_MM, id=c("Gene.ID"))
green_mod_05_MM_melt$Primary <- substr(green_mod_05_MM_melt$variable, 4,4)
green_mod_05_MM_melt$Second <- substr(green_mod_05_MM_melt$variable, 5,5)

ggplot(green_mod_05_MM_melt, aes(x=Second, y=value, group=Gene.ID, color = Primary)) +
  theme(panel.grid=element_blank()) +
  scale_color_manual(values=c("#56B4E9","#D55E00")) +
  geom_line(size=0.1, alpha=0.2) +
  theme_bw() +
  ggtitle("Day 7 WGCNA: green module") +
  ylab("Pearson's correlation of vst Expression and Treatment Group (genes with Module Membership < 0.001)") +
  facet_wrap(~Primary)


# brown module - higher in AM than MS and MA 
length(colnames(dds.d7_vst)[moduleColors=="brown"]) # 595 total genes in the brown module
brown_mod <- geneInfo_GROUPS %>%  dplyr::filter(moduleColor %in% 'brown')
nrow(brown_mod) # 595 total genes in the brown module
brown_mod_05 <- brown_mod %>%  dplyr::filter(p.MM.brown < 0.001) 
brown_mod_05_MM <-brown_mod_05[,c(8:13)]
nrow(brown_mod_05_MM) # 426 total genes in the brown module at 0.0001 module membership 

brown_mod_05_MM$Gene.ID <- rownames(brown_mod_05_MM)
brown_mod_05_MM_melt <- melt(brown_mod_05_MM, id=c("Gene.ID"))
brown_mod_05_MM_melt$Primary <- substr(brown_mod_05_MM_melt$variable, 4,4)
brown_mod_05_MM_melt$Second <- substr(brown_mod_05_MM_melt$variable, 5,5)

ggplot(brown_mod_05_MM_melt, aes(x=Second, y=value, group=Gene.ID, color = Primary)) +
  theme(panel.grid=element_blank()) +
  scale_color_manual(values=c("#56B4E9","#D55E00")) +
  geom_line(size=0.1, alpha=0.2) +
  theme_bw() +
  ggtitle("Day 7 WGCNA: brown module") +
  ylab("Pearson's correlation of vst Expression and Treatment Group (genes with Module Membership < 0.001)") +
  facet_wrap(~Primary)

# black module - lower in AS than other modules 
length(colnames(dds.d7_vst)[moduleColors=="black"]) # 595 total genes in the black module
black_mod <- geneInfo_GROUPS %>%  dplyr::filter(moduleColor %in% 'black')
nrow(black_mod) # 595 total genes in the black module
black_mod_05 <- black_mod %>%  dplyr::filter(p.MM.black < 0.001) 
black_mod_05_MM <-black_mod_05[,c(8:13)]
nrow(black_mod_05_MM) # 164 total genes in the black module at 0.0001 module membership 

black_mod_05_MM$Gene.ID <- rownames(black_mod_05_MM)
black_mod_05_MM_melt <- melt(black_mod_05_MM, id=c("Gene.ID"))
black_mod_05_MM_melt$Primary <- substr(black_mod_05_MM_melt$variable, 4,4)
black_mod_05_MM_melt$Second <- substr(black_mod_05_MM_melt$variable, 5,5)

ggplot(black_mod_05_MM_melt, aes(x=Second, y=value, group=Gene.ID, color = Primary)) +
  theme(panel.grid=element_blank()) +
  scale_color_manual(values=c("#56B4E9","#D55E00")) +
  geom_line(size=0.1, alpha=0.2) +
  theme_bw() +
  ggtitle("Day 7 WGCNA: black module") +
  ylab("Pearson's correlation of vst Expression and Treatment Group (genes with Module Membership < 0.001)") +
  facet_wrap(~Primary)


#=====================================================================================
#
#  Write csv for the modules and corresponding raw read counts
#
#=====================================================================================
# call the module of interest for follow-up GO analysis 

write.csv(geneInfo_GROUPS, file = "Analysis/Output/WGCNA/Day7/d7.WGCNA_ModulMembership.csv")

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
# Load data 
d7_ModEigen <- read.csv("Analysis/Output/WGCNA/Day7/d7.WGCNA_ModulEigengenes.csv")
d7_Annot_ModuleMembership <-  read.csv("Analysis/Output/WGCNA/Day7/d7.WGCNA_ModulMembership.csv")
d7_vst_data <- read.csv("Analysis/Output/WGCNA/Day7/Day7_vstTransformed_WGCNAdata.csv")
Master.Treatment_Phenotype.data <- read.csv(file="Analysis/Data/ Master_Phyenotype.and.Exp.Treatment_Metadata.csv", sep=',', header=TRUE)



#=====================================================================================
#
#  Prep and merge datasets for plotting
# 
#=====================================================================================
# Prep data to merge
# module membership
d7_Annot_ModuleMembership <- d7_Annot_ModuleMembership[-c(1:2)] # ommit the redundant gene name columns buut keep 'geneSymbol'
dim(d7_Annot_ModuleMembership) #   8548   30

# module eigengenes
names(d7_ModEigen)[1] <- "Sample.Name"

# normalized count data (same used for the dds object and WGCNA analysis)
d7_vst_data_t <- as.data.frame(t(d7_vst_data[, -(1)])) # trsnpose columns names to rows (genes) 
colnames(d7_vst_data_t) <- d7_vst_data[, 1] # name the columns as the rows in previous dataset (sample IDS)
d7_vst_data_t<- d7_vst_data_t %>% tibble::rownames_to_column("geneSymbol") # "geneSymbol" - create column and remname the gene name column to 'geneSymbol'


# merge Master datasets
# GO terms
d7_GOTermsMaster.Modules<-  merge(d7_Annot_ModuleMembership, d7_vst_data_t, by = "geneSymbol")
dim(d7_GOTermsMaster.Modules) # 8548   65
# Eigengenes and traits
d7_EigenTraitMaster<-  merge(d7_ModEigen, Master.Treatment_Phenotype.data, by = "Sample.Name")
dim(d7_EigenTraitMaster) # 35 20




#=====================================================================================
# 
# Linear regression plots Eigenengene expression with Traits for EACH module       
#
#=====================================================================================

# lets loop this through to get a bunch of plots 
phys <- d7_EigenTraitMaster[,c(15:20)]
modules <- d7_EigenTraitMaster[,c(2:7)]


library(ggpmisc)

for(i in 1:ncol(modules)) {
  module_color <- substr(names(modules)[i], 3,12)
  moduledata <- modules[,i]
  
  par(mfrow=c(3,2))
  
  my.formula <- y ~ x
  
  p <- ggplot(d7_EigenTraitMaster, 
              aes(x=moduledata, y=phys[,1])) + # linetype= Primary_Treatment, group = Primary_Treatment
    scale_fill_manual(values=c("#56B4E9", "#D55E00")) +
    geom_smooth(method=lm , color="grey50", alpha = 0.1, se=TRUE, linetype = "dashed", formula = my.formula) + #, linetype = "dashed"
    stat_poly_eq(formula = my.formula,
                 eq.with.lhs = "italic(hat(y))~`=`~",
                 aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
                 parse = TRUE)  + 
    geom_point(aes(colour = factor(Primary_Treatment)), size = 3, shape = 19) +
    scale_color_manual(values=c("#56B4E9", "#D55E00")) +
    labs(x="Eigengene expression by sample", y = paste(names(phys)[1])) +
    theme_bw()
  p1 <- p + ggtitle(paste("WGCNAmodule:", module_color, "&",(names(phys)[1]), sep =' ')) + theme(legend.position="none")
  
  
  p <- ggplot(d7_EigenTraitMaster, 
              aes(x=moduledata, y=phys[,2])) + # linetype= Primary_Treatment, group = Primary_Treatment
    scale_fill_manual(values=c("#56B4E9", "#D55E00")) +
    geom_smooth(method=lm , color="grey50", alpha = 0.1, se=TRUE, linetype = "dashed", formula = my.formula) + #, linetype = "dashed"
    stat_poly_eq(formula = my.formula,
                 eq.with.lhs = "italic(hat(y))~`=`~",
                 aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
                 parse = TRUE)  + 
    geom_point(aes(colour = factor(Primary_Treatment)), size = 3, shape = 19) +
    scale_color_manual(values=c("#56B4E9", "#D55E00")) +
    labs(x="Eigengene expression by sample", y = paste(names(phys)[2])) +
    theme_bw()
  p2 <- p + ggtitle(paste("WGCNAmodule:", module_color, "&",(names(phys)[2]), sep =' ')) + theme(legend.position="none")
  
  
  p <- ggplot(d7_EigenTraitMaster, 
              aes(x=moduledata, y=phys[,3])) + # linetype= Primary_Treatment, group = Primary_Treatment
    scale_fill_manual(values=c("#56B4E9", "#D55E00")) +
    geom_smooth(method=lm , color="grey50", alpha = 0.1, se=TRUE, linetype = "dashed", formula = my.formula) + #, linetype = "dashed"
    stat_poly_eq(formula = my.formula,
                 eq.with.lhs = "italic(hat(y))~`=`~",
                 aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
                 parse = TRUE)  + 
    geom_point(aes(colour = factor(Primary_Treatment)), size = 3, shape = 19) +
    scale_color_manual(values=c("#56B4E9", "#D55E00")) +
    labs(x="Eigengene expression by sample", y = paste(names(phys)[3])) +
    theme_bw()
  p3 <- p + ggtitle(paste("WGCNAmodule:", module_color, "&",(names(phys)[3]), sep =' ')) + theme(legend.position="none")
  
  
  p <- ggplot(d7_EigenTraitMaster, 
              aes(x=moduledata, y=phys[,4])) + # linetype= Primary_Treatment, group = Primary_Treatment
    scale_fill_manual(values=c("#56B4E9", "#D55E00")) +
    geom_smooth(method=lm , color="grey50", alpha = 0.1, se=TRUE, linetype = "dashed", formula = my.formula) + #, linetype = "dashed"
    stat_poly_eq(formula = my.formula,
                 eq.with.lhs = "italic(hat(y))~`=`~",
                 aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
                 parse = TRUE)  + 
    geom_point(aes(colour = factor(Primary_Treatment)), size = 3, shape = 19) +
    scale_color_manual(values=c("#56B4E9", "#D55E00")) +
    labs(x="Eigengene expression by sample", y = paste(names(phys)[4])) +
    theme_bw()
  p4 <- p + ggtitle(paste("WGCNAmodule:", module_color, "&",(names(phys)[4]), sep =' ')) + theme(legend.position="none")
  
  
  p <- ggplot(d7_EigenTraitMaster, 
              aes(x=moduledata, y=phys[,5])) + # linetype= Primary_Treatment, group = Primary_Treatment
    scale_fill_manual(values=c("#56B4E9", "#D55E00")) +
    geom_smooth(method=lm , color="grey50", alpha = 0.1, se=TRUE, linetype = "dashed", formula = my.formula) + #, linetype = "dashed"
    stat_poly_eq(formula = my.formula,
                 eq.with.lhs = "italic(hat(y))~`=`~",
                 aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
                 parse = TRUE)  + 
    geom_point(aes(colour = factor(Primary_Treatment)), size = 3, shape = 19) +
    scale_color_manual(values=c("#56B4E9", "#D55E00")) +
    labs(x="Eigengene expression by sample", y = paste(names(phys)[5])) +
    theme_bw()
  p5 <- p + ggtitle(paste("WGCNAmodule:", module_color, "&",(names(phys)[5]), sep =' ')) + theme(legend.position="none")
  
  p <- ggplot(d7_EigenTraitMaster, 
              aes(x=moduledata, y=phys[,6])) + # linetype= Primary_Treatment, group = Primary_Treatment
    scale_fill_manual(values=c("#56B4E9", "#D55E00")) +
    geom_smooth(method=lm , color="grey50", alpha = 0.1, se=TRUE, linetype = "dashed", formula = my.formula) + #, linetype = "dashed"
    stat_poly_eq(formula = my.formula,
                 eq.with.lhs = "italic(hat(y))~`=`~",
                 aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
                 parse = TRUE)  + 
    geom_point(aes(colour = factor(Primary_Treatment)), size = 3, shape = 19) +
    scale_color_manual(values=c("#56B4E9", "#D55E00")) +
    labs(x="Eigengene expression by sample", y = paste(names(phys)[6])) +
    theme_bw()
  p6 <- p + ggtitle(paste("WGCNAmodule:", module_color, "&",(names(phys)[6]), sep =' ')) + theme(legend.position="none")
  
  # save plot grid 
  #png("Analysis/Output/WGCNA/Day21/Day21_PhysOnly_heatmap2.png", 1000, 1000, pointsize=20)
  ggarrange(p1, p2, p3, p4, p5, p6, ncol = 2, nrow = 3)
  ggsave(paste0("Analysis/Output/WGCNA/Day7/EigengenePlots/Day7_EigengenePlot_", module_color, ".png"))
}



#===================================================================================== 
# 
# EXPLORE THE Expression of each module (for loop plots!) BY TREATMENT
#
# evidence from the heatmaps and the regression shwo strong module membership and eigenenge cor strength with lighcyan 
# in response to treatment (primary specifically) and Total Antioxidant capacity 
#
#=====================================================================================

modcolor <- as.data.frame(unique(d7_Annot_ModuleMembership$moduleColor))
names(modcolor)[1] <- "color"

# experiment treatment and total protein data - narrow the columns 
exp.phys_data <- Master.Treatment_Phenotype.data %>% 
  dplyr::filter(Date %in% 20190731) # filter out Day 7 data only 


for(i in 1:nrow(modcolor)) {
  
  # vst read count date - narrow the columns - reshape and rename
  Mod_geneIDs <- d7_Annot_ModuleMembership %>% dplyr::filter(moduleColor %in% modcolor[i,]) %>%  dplyr::select("geneSymbol")
  d7_vst_Mod <- d7_vst_data_t %>% dplyr::filter(geneSymbol %in% Mod_geneIDs[,1])
  d7_vst_Mod_MELT <- melt(d7_vst_Mod, id=("geneSymbol")) # melt using reshape2
  names(d7_vst_Mod_MELT)[(2:3)] <- c('Sample.Name', 'vst_Expression') # change column names
  
  # merge by common row values 'Sample.Name'
  merged_Expdata_Mod <- merge(d7_vst_Mod_MELT, exp.phys_data, by ='Sample.Name')
  
  # mean Exp response table 
  meanEXp_Mod <- merged_Expdata_Mod %>% 
    select(c('Sample.Name','vst_Expression','Primary_Treatment', 'Second_Treament')) %>% 
    group_by(Sample.Name, Primary_Treatment, Second_Treament) %>%
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
  # Second treatment
  meanEXp_Summary.Sec_Mod <- meanEXp_Mod %>% 
    group_by(Second_Treament,Primary_Treatment) %>%
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
  
  
  # Second treatment mean sd plot
  min_p2 <- min(meanEXp_Summary.Sec_Mod$mean) - max(meanEXp_Summary.Sec_Mod$se)
  max_p2 <- max(meanEXp_Summary.Sec_Mod$mean) + max(meanEXp_Summary.Sec_Mod$se)
  
  Sec.vst.Mod <- ggplot(meanEXp_Summary.Sec_Mod, aes(x=Second_Treament, y=mean, fill=Primary_Treatment, group=Primary_Treatment)) +
    theme_classic() +
    geom_errorbar(aes(ymin=mean-se, ymax=mean+se), colour="black", width=.1, position=pd) +
    geom_line(position=pd) +
    geom_point(position=pd, size = 4, shape=21) +            
    xlab("Second pCO2 treatment") +
    ylab(paste(modcolor[i,]," Module VST Gene Expression (Mean +/- SE)", sep = ' ')) +                 # note the mean was first by sample ID THEN by treatment
    scale_fill_manual(values=c("#56B4E9","#E69F00")) +
    # scale_color_manual(values=c("#56B4E9","#E69F00")) +
    # ggtitle("Day 21 WGCNA red' Module VST GeneExp") +
    # expand_limits(y=0) +                                                    # Expand y range
    scale_y_continuous(limits=c((min_p2), (max_p2))) +
    theme(text = element_text(size=15))

  
  # output 
  
  pdf(paste("Analysis/Output/WGCNA/Day7/ModuleExpression_Treatment/day7_Exp_Module",modcolor[i,],".pdf"), width=6, height=6)
  print(ggarrange(Primary.vst.Mod, Sec.vst.Mod,         
                  plotlist = NULL,
                  ncol = 1,
                  nrow = 2,
                  labels = NULL))
  dev.off()
  
}



#===================================================================================== 
# 
# EXPLORE THE Expression of each module (for loop plots!)  BY PHENOTYPE 
#
# evidence from the heatmaps and the regression shwo strong module membership and eigenenge cor strength with lighcyan 
# in response to treatment (primary specifically) and Total Antioxidant capacity 
#
#=====================================================================================



for(i in 1:nrow(modcolor)) {
  module_color <- modcolor[i,] # call the color character 
  
  # vst read count date - narrow the columns - reshape and rename
  Mod_geneIDs <- d7_Annot_ModuleMembership %>% dplyr::filter(moduleColor %in% module_color) %>%  dplyr::select("geneSymbol") # all geneSymbols (IDs) WITHIN the looped module color 
  d7_vst_Mod <- d7_vst_data_t %>% dplyr::filter(geneSymbol %in% Mod_geneIDs[,1]) # call the transposed vst count data and filer for the module color 
  d7_vst_Mod_MELT <- melt(d7_vst_Mod, id=("geneSymbol")) # melt using reshape2 - ('un-transpose')
  names(d7_vst_Mod_MELT)[(2:3)] <- c('Sample.Name', 'vst_Expression') # change column names to merge next
  
  # merge by common row values 'Sample.Name'
  merged_Expdata_Mod <- merge(d7_vst_Mod_MELT, exp.phys_data, by ='Sample.Name') # merge
  
  # mean Exp response table 
  meanExp_Mod <- merged_Expdata_Mod %>% 
    group_by(Sample.Name) %>%
    dplyr::summarize(n = n(), 
                     mean.vstExp = mean(vst_Expression), 
                     sd.vsdtExp = sd(vst_Expression),
                     se.vsdtExp = sd.vsdtExp/sqrt(n)) # mean expression for each sample ID 
  
  merged_Expdata_phystreat <- merged_Expdata_Mod %>%  # condence down the 'merged_Expdata_Mod' file to merge with 'meanExp_Mod'
    dplyr::select(c("Sample.Name", "Primary_Treatment", "Second_Treament", "mean.µmol.CRE.g.protein",
                    "mean.AFDW", "mean.mgProt.mgAFDW", "mean.Resp.ugLhrindiv", "mean.Biovol", "mean.Shell.Length"))
  
  meanExp_Mod.MASTER <- unique(merge(meanExp_Mod, merged_Expdata_phystreat, by = "Sample.Name")) # call only unqiue values 
 
  phys <- meanExp_Mod.MASTER[, 8:13]  # call only the phys columns for the subseqent for loop - plotting each 
  
  dir.create(paste0("Analysis/Output/WGCNA/Day7/ModuleExpresion_Phenotype/",module_color), showWarnings = FALSE)
  
  for(m in 1:ncol(phys)) {
   # phys[,m] # the column
   # names(phys)[m] # the name of the column 
  
  
  Primary_plot <- ggplot(meanExp_Mod.MASTER, 
          aes(x=phys[,m], y=mean.vstExp, fill = Primary_Treatment, colour = Primary_Treatment, group = Primary_Treatment)) + # linetype= Primary_Treatment, group = Primary_Treatment
     scale_fill_manual(values=c("#56B4E9", "#D55E00")) +
     geom_smooth(method=lm , alpha = 0.1, se=TRUE, linetype = "dashed", formula = my.formula) + #, linetype = "dashed"
     stat_poly_eq(formula = my.formula,
                  eq.with.lhs = "italic(hat(y))~`=`~",
                  aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
                  parse = TRUE)  + 
     geom_point(aes(colour = factor(Primary_Treatment)), size = 3, shape = 19) +
     scale_color_manual(values=c("#56B4E9", "#D55E00")) +
     labs(x=(names(phys)[m]), y = paste("WGCNA module",module_color,'mean.vstExp', sep = ' ')) +
     theme_bw()
 
   PrimSec_facetplot <- ggplot(meanExp_Mod.MASTER, 
                        aes(x=phys[,m], y=mean.vstExp, fill = Primary_Treatment, colour = Primary_Treatment, group = Primary_Treatment)) + # linetype= Primary_Treatment, group = Primary_Treatment
    scale_fill_manual(values=c("#56B4E9", "#D55E00")) +
    geom_smooth(method=lm , alpha = 0.1, se=TRUE, linetype = "dashed", formula = my.formula) + #, linetype = "dashed"
    stat_poly_eq(formula = my.formula,
                 eq.with.lhs = "italic(hat(y))~`=`~",
                 aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
                 parse = TRUE)  + 
    geom_point(aes(colour = factor(Primary_Treatment)), size = 3, shape = 19) +
    scale_color_manual(values=c("#56B4E9", "#D55E00")) +
    labs(x=(names(phys)[m]), y = paste("WGCNA module",module_color,'mean.vstExp', sep = ' ')) +
    theme_bw()+
    facet_wrap(~Second_Treament) 
  PrimSec_facetplot <- PrimSec_facetplot + ggtitle(paste("WGCNAmodule:", module_color, "&",(names(phys)[m]), sep =' ')) + theme(legend.position="none")
  
  png(paste("Analysis/Output/WGCNA/Day7/ModuleExpresion_Phenotype/",module_color,"/Day7_Exp_Module_",module_color, "_",(names(phys)[m]),".png", sep = ''), 600, 1000, pointsize=20)
  print(ggarrange(Primary_plot, PrimSec_facetplot,         
                  plotlist = NULL,
                  ncol = 1,
                  nrow = 2,
                  labels = NULL))
  dev.off()
  }
  
}
   


#===================================================================================== 
# 
# goseq - load the annotation and prepare the four steps for goseq!
#
#===================================================================================== 
library(forcats) # for plotting later..
### Panopea generosa - load .fna ('Geoduck_annotation') and foramt GO terms ('Geoduck_GOterms') and vectors
Geoduck_annotation <- read.delim2(file="C:/Users/samjg/Documents/My_Projects/Pgenerosa_TagSeq_Metabolomics/TagSeq/Seq_details/Panopea-generosa-genes-annotations.txt", header=F)


# build annotation file to merge with the mean LFC tables
annot.condenced <- Geoduck_annotation[,c(1,3:9)]
annot.condenced$gene.length <- annot.condenced$V4 - annot.condenced$V3
annot.condenced <- annot.condenced[,-c(2,3)]
names(annot.condenced) <- c('Gene.ID', 'Uniprot', 'HGNC', 'fxn', 'Go.terms', 'Go.fxns','gene.length')


# Prepare dataframe(s) and vectors for goseq 
### (1) Format 'GO.term' for goseq from the P.generosa annotation .fna file 'Geoduck_annotation'
# GO terms data (ALL)
Geoduck_GOterms <- as.data.frame(Geoduck_annotation) %>% dplyr::select(c('V1','V8'))
colnames(Geoduck_GOterms)[1:2] <- c('transcript.ID', 'GO.terms') # call gene name and the GO terms - (Uniprot ID 'V5')
splitted <- strsplit(as.character(Geoduck_GOterms$GO.terms), ";") #slit into multiple GO ids
GO.terms <- data.frame(v1 = rep.int(Geoduck_GOterms$transcript.ID, sapply(splitted, length)), v2 = unlist(splitted)) #list all genes with each of their GO terms in a single row

### (2) Unique Genes - vector based on all unique mapped reads 
# Construct a named vector of all target genes for goseq
GO_unique.genes.all <- as.vector(unique(Geoduck_annotation$V1)) # call all unique genes for GO analysis (goseq)

### (3) Gene length 
# length vector  
GO_gene.length <- Geoduck_annotation %>% dplyr::mutate(length = V4-V3) %>%  dplyr::select(c("V1","length"))
names(GO_gene.length)[1] <- "Gene.ID"
#GO_gene.length_merge <- merge(GO_gene.length, GO_magenta_genes, by = "Gene.ID")
length_vector <- GO_gene.length$length

### (4) Call modules of interest (i.e. magenta module) - merge to the genes list 'GO_unique.genes.all' to create a binary vector 
# Example: 0 = not in moddulecolor; 1 = in modulecolor

# target modules with significant correlation in the group treatments - modules green, black, brown, and red 
# green module  ------------------------------------------- #

GO_green <- d7_Annot_ModuleMembership %>% dplyr::filter(moduleColor %in% 'green') # %>%  dplyr::select("geneSymbol")
GO_green_genes <- GO_green[1]
names(GO_green_genes)[1] <- "Gene.ID" # 162 genws in the green module 


GO_green_integer <- as.integer(GO_unique.genes.all%in%(GO_green_genes$Gene.ID)) # Day 0 - Primary - Upregulated DEGs
names(GO_green_integer)=GO_unique.genes.all
View(GO_green_integer)

# black module  ------------------------------------------- #

GO_black <- d7_Annot_ModuleMembership %>% dplyr::filter(moduleColor %in% 'black') # %>%  dplyr::select("geneSymbol")
GO_black_genes <- GO_black[1]
names(GO_black_genes)[1] <- "Gene.ID" # 162 genws in the black module 


GO_black_integer <- as.integer(GO_unique.genes.all%in%(GO_black_genes$Gene.ID)) # Day 0 - Primary - Upregulated DEGs
names(GO_black_integer)=GO_unique.genes.all
View(GO_black_integer)

# brown module  ------------------------------------------- #

GO_brown <- d7_Annot_ModuleMembership %>% dplyr::filter(moduleColor %in% 'brown') # %>%  dplyr::select("geneSymbol")
GO_brown_genes <- GO_brown[1]
names(GO_brown_genes)[1] <- "Gene.ID" # 162 genws in the brown module 


GO_brown_integer <- as.integer(GO_unique.genes.all%in%(GO_brown_genes$Gene.ID)) # Day 0 - Primary - Upregulated DEGs
names(GO_brown_integer)=GO_unique.genes.all
View(GO_brown_integer)

# red module  ------------------------------------------- #

GO_red <- d7_Annot_ModuleMembership %>% dplyr::filter(moduleColor %in% 'red') # %>%  dplyr::select("geneSymbol")
GO_red_genes <- GO_red[1]
names(GO_red_genes)[1] <- "Gene.ID" # 162 genws in the red module 


GO_red_integer <- as.integer(GO_unique.genes.all%in%(GO_red_genes$Gene.ID)) # Day 0 - Primary - Upregulated DEGs
names(GO_red_integer)=GO_unique.genes.all
View(GO_red_integer)

# Review what we have for goseq....
# (1) Go terms ---------------- #
#  head(GO.terms)

# (2) ID vector --------------- #
# head(GO_unique.genes.all)

# (3) Length vector ----------- #
# head(length_vector)

# (4) Binary DEG vectors ------ # 
# head(GO_magenta_integer)  

#===================================================================================== .
# It's goseq time!!!
# Calculate Probability Weighting Function
#===================================================================================== 

#Calculate Probability Weighting Function (using 'nullp')
green.pwf<-nullp(GO_green_integer, GO_unique.genes.all, bias.data=length_vector) #weight vector by length of gene
black.pwf<-nullp(GO_black_integer, GO_unique.genes.all, bias.data=length_vector) #weight vector by length of gene
brown.pwf<-nullp(GO_brown_integer, GO_unique.genes.all, bias.data=length_vector) #weight vector by length of gene
red.pwf  <-nullp(GO_red_integer, GO_unique.genes.all, bias.data=length_vector) #weight vector by length of gene

#===================================================================================== 
# It's goseq time!!!
#  Run goseq
#===================================================================================== 

green.goseq<-goseq(green.pwf, ID.vector, gene2cat=GO.terms, test.cats=c("GO:CC", "GO:BP", "GO:MF"), method="Wallenius", use_genes_without_cat=TRUE)
black.goseq<-goseq(black.pwf, ID.vector, gene2cat=GO.terms, test.cats=c("GO:CC", "GO:BP", "GO:MF"), method="Wallenius", use_genes_without_cat=TRUE)
brown.goseq<-goseq(brown.pwf, ID.vector, gene2cat=GO.terms, test.cats=c("GO:CC", "GO:BP", "GO:MF"), method="Wallenius", use_genes_without_cat=TRUE)
red.goseq  <-goseq(red.pwf, ID.vector, gene2cat=GO.terms, test.cats=c("GO:CC", "GO:BP", "GO:MF"), method="Wallenius", use_genes_without_cat=TRUE)

# somegenes <- magenta.goseq %>%  dplyr::filter(ontology %in% 'CC') # xplore why there were no CC terms in the goseq plot
# somegenes0.5CC <- somegenes %>% dplyr::filter(over_represented_pvalue < 0.1)
#===================================================================================== 
# It's goseq time!!!
# PLOTTING
#  GREEN MODULE 
#===================================================================================== 

#Find only enriched GO terms that are statistically significant at cutoff 
green_enriched.GO.05.a<-green.goseq$category[green.goseq$over_represented_pvalue<.05] # change twice here
green_enriched.GO.05<-data.frame(green_enriched.GO.05.a)
colnames(green_enriched.GO.05) <- c("category")
green_enriched.GO.05 <- merge(green_enriched.GO.05, green.goseq, by="category") # change here
green_enriched.GO.05 <- green_enriched.GO.05[order(-green_enriched.GO.05$numDEInCat),]
green_enriched.GO.05$term <- as.factor(green_enriched.GO.05$term)
head(green_enriched.GO.05)

green_MF <- subset(green_enriched.GO.05, ontology=="MF")
green_MF <- green_MF[order(-green_MF$numDEInCat),]
green_CC <- subset(green_enriched.GO.05, ontology=="CC")
green_CC <- green_CC[order(-green_CC$numDEInCat),]
green_BP <- subset(green_enriched.GO.05, ontology=="BP")
green_BP <- green_BP[order(-green_BP$numDEInCat),]

# Molecular Processes
MFplot_green <- green_MF %>% mutate(term = fct_reorder(term, numDEInCat)) %>%
  ggplot( aes(x=term, y=numDEInCat) ) +
  geom_segment( aes(x=term ,xend=term, y=0, yend=numDEInCat), color="grey") +
  geom_point(size=3, color="#69b3a2") +
  coord_flip() +
  theme(
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_blank(),
    legend.position="none"
  ) +
  xlab("") +
  ylab("") +
  ggtitle("Molecular Function") + #add a main title
  theme(plot.title = element_text(face = 'bold', 
                                  size = 12, 
                                  hjust = 0)) +
  theme_bw() + #Set background color 
  theme(panel.border = element_blank(), # Set border
        panel.grid.major = element_blank(), #Set major gridlines
        panel.grid.minor = element_blank(), #Set minor gridlines
        axis.line = element_line(colour = "black"), #Set axes color
        plot.background=element_blank())#Set the plot background

#Cellular Processes
CCplot_green <- green_CC %>% mutate(term = fct_reorder(term, numDEInCat)) %>%
  ggplot( aes(x=term, y=numDEInCat) ) +
  geom_segment( aes(x=term ,xend=term, y=0, yend=numDEInCat), color="grey") +
  geom_point(size=3, color="#69b3a2") +
  coord_flip() +
  theme(
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_blank(),
    legend.position="none"
  ) +
  xlab("") +
  ylab("") +
  ggtitle("Cellular Component") + #add a main title
  theme(plot.title = element_text(face = 'bold', 
                                  size = 12, 
                                  hjust = 0)) +
  theme_bw() + #Set background color 
  theme(panel.border = element_blank(), # Set border
        panel.grid.major = element_blank(), #Set major gridlines
        panel.grid.minor = element_blank(), #Set minor gridlines
        axis.line = element_line(colour = "black"), #Set axes color
        plot.background=element_blank())#Set the plot background

# Biological Processes 
BPplot_green <- green_BP %>% mutate(term = fct_reorder(term, numDEInCat)) %>%
  ggplot( aes(x=term, y=numDEInCat) ) +
  geom_segment( aes(x=term ,xend=term, y=0, yend=numDEInCat), color="grey") +
  geom_point(size=3, color="#69b3a2") +
  coord_flip() +
  theme(
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_blank(),
    legend.position="none") +
  xlab("") +
  ylab("") +
  ggtitle("Biological Process") + #add a main title
  theme(plot.title = element_text(face = 'bold', 
                                  size = 12, 
                                  hjust = 0)) +
  theme_bw() + #Set background color 
  theme(panel.border = element_blank(), # Set border
        panel.grid.major = element_blank(), #Set major gridlines
        panel.grid.minor = element_blank(), #Set minor gridlines
        axis.line = element_line(colour = "black"), #Set axes color
        plot.background=element_blank())#Set the plot background

# merge MF CC BP plots - Light cyan
num_green   <- dim(GO_green_genes)[1] # call num upregulated genes

library(tidyr)
pdf("Analysis/Output/WGCNA/Day7/goseq_modules/Day7_goseq_MEgreen.pdf", width=8, height=5)
green_enriched.GO.05 %>% drop_na(ontology) %>% mutate(term = fct_reorder(term, numDEInCat)) %>%
  mutate(term = fct_reorder(term, ontology)) %>%
  ggplot( aes(x=term, y=numDEInCat) ) +
  geom_segment( aes(x=term ,xend=term, y=0, yend=numDEInCat), color="grey") +
  geom_point(size=3, aes(colour = ontology)) +
  coord_flip() +
  theme(
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_blank(),
    legend.position="bottom"
  ) +
  xlab("") +
  ylab("") +
  ggtitle("GO: Day 7 WGCNA MEgreen") +
  geom_label(aes(x = 2, y = 4, label = paste(num_green, "MEgreen genes"))) +
  theme_bw() + #Set background color 
  theme(panel.border = element_blank(), # Set border
        panel.grid.major = element_blank(), #Set major gridlines
        panel.grid.minor = element_blank(), #Set minor gridlines
        axis.line = element_line(colour = "black"), #Set axes color
        plot.background=element_blank()) #Set the plot background #set title attributes
dev.off()

#===================================================================================== 
# It's goseq time!!!
# PLOTTING
#  black MODULE 
#===================================================================================== 

#Find only enriched GO terms that are statistically significant at cutoff 
black_enriched.GO.05.a<-black.goseq$category[black.goseq$over_represented_pvalue<.05] # change twice here
black_enriched.GO.05<-data.frame(black_enriched.GO.05.a)
colnames(black_enriched.GO.05) <- c("category")
black_enriched.GO.05 <- merge(black_enriched.GO.05, black.goseq, by="category") # change here
black_enriched.GO.05 <- black_enriched.GO.05[order(-black_enriched.GO.05$numDEInCat),]
black_enriched.GO.05$term <- as.factor(black_enriched.GO.05$term)
head(black_enriched.GO.05)

black_MF <- subset(black_enriched.GO.05, ontology=="MF")
black_MF <- black_MF[order(-black_MF$numDEInCat),]
black_CC <- subset(black_enriched.GO.05, ontology=="CC")
black_CC <- black_CC[order(-black_CC$numDEInCat),]
black_BP <- subset(black_enriched.GO.05, ontology=="BP")
black_BP <- black_BP[order(-black_BP$numDEInCat),]

# Molecular Processes
MFplot_black <- black_MF %>% mutate(term = fct_reorder(term, numDEInCat)) %>%
  ggplot( aes(x=term, y=numDEInCat) ) +
  geom_segment( aes(x=term ,xend=term, y=0, yend=numDEInCat), color="grey") +
  geom_point(size=3, color="#69b3a2") +
  coord_flip() +
  theme(
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_blank(),
    legend.position="none"
  ) +
  xlab("") +
  ylab("") +
  ggtitle("Molecular Function") + #add a main title
  theme(plot.title = element_text(face = 'bold', 
                                  size = 12, 
                                  hjust = 0)) +
  theme_bw() + #Set background color 
  theme(panel.border = element_blank(), # Set border
        panel.grid.major = element_blank(), #Set major gridlines
        panel.grid.minor = element_blank(), #Set minor gridlines
        axis.line = element_line(colour = "black"), #Set axes color
        plot.background=element_blank())#Set the plot background

#Cellular Processes
CCplot_black <- black_CC %>% mutate(term = fct_reorder(term, numDEInCat)) %>%
  ggplot( aes(x=term, y=numDEInCat) ) +
  geom_segment( aes(x=term ,xend=term, y=0, yend=numDEInCat), color="grey") +
  geom_point(size=3, color="#69b3a2") +
  coord_flip() +
  theme(
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_blank(),
    legend.position="none"
  ) +
  xlab("") +
  ylab("") +
  ggtitle("Cellular Component") + #add a main title
  theme(plot.title = element_text(face = 'bold', 
                                  size = 12, 
                                  hjust = 0)) +
  theme_bw() + #Set background color 
  theme(panel.border = element_blank(), # Set border
        panel.grid.major = element_blank(), #Set major gridlines
        panel.grid.minor = element_blank(), #Set minor gridlines
        axis.line = element_line(colour = "black"), #Set axes color
        plot.background=element_blank())#Set the plot background

# Biological Processes 
BPplot_black <- black_BP %>% mutate(term = fct_reorder(term, numDEInCat)) %>%
  ggplot( aes(x=term, y=numDEInCat) ) +
  geom_segment( aes(x=term ,xend=term, y=0, yend=numDEInCat), color="grey") +
  geom_point(size=3, color="#69b3a2") +
  coord_flip() +
  theme(
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_blank(),
    legend.position="none") +
  xlab("") +
  ylab("") +
  ggtitle("Biological Process") + #add a main title
  theme(plot.title = element_text(face = 'bold', 
                                  size = 12, 
                                  hjust = 0)) +
  theme_bw() + #Set background color 
  theme(panel.border = element_blank(), # Set border
        panel.grid.major = element_blank(), #Set major gridlines
        panel.grid.minor = element_blank(), #Set minor gridlines
        axis.line = element_line(colour = "black"), #Set axes color
        plot.background=element_blank())#Set the plot background

# merge MF CC BP plots - Light cyan
num_black   <- dim(GO_black_genes)[1] # call num upregulated genes

library(tidyr)
pdf("Analysis/Output/WGCNA/Day7/goseq_modules/Day7_goseq_MEblack.pdf", width=9, height=7)
black_enriched.GO.05 %>% drop_na(ontology) %>% mutate(term = fct_reorder(term, numDEInCat)) %>%
  mutate(term = fct_reorder(term, ontology)) %>%
  ggplot( aes(x=term, y=numDEInCat) ) +
  geom_segment( aes(x=term ,xend=term, y=0, yend=numDEInCat), color="grey") +
  geom_point(size=3, aes(colour = ontology)) +
  coord_flip() +
  theme(
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_blank(),
    legend.position="bottom"
  ) +
  xlab("") +
  ylab("") +
  ggtitle("GO: Day 7 WGCNA MEblack") +
  geom_label(aes(x = 2, y = 15, label = paste(num_black, "MEblack genes"))) +
  theme_bw() + #Set background color 
  theme(panel.border = element_blank(), # Set border
        panel.grid.major = element_blank(), #Set major gridlines
        panel.grid.minor = element_blank(), #Set minor gridlines
        axis.line = element_line(colour = "black"), #Set axes color
        plot.background=element_blank()) #Set the plot background #set title attributes
dev.off()

#===================================================================================== 
# It's goseq time!!!
# PLOTTING
#  BROWN MODULE 
#===================================================================================== 

#Find only enriched GO terms that are statistically significant at cutoff 
brown_enriched.GO.05.a<-brown.goseq$category[brown.goseq$over_represented_pvalue<.05] # change twice here
brown_enriched.GO.05<-data.frame(brown_enriched.GO.05.a)
colnames(brown_enriched.GO.05) <- c("category")
brown_enriched.GO.05 <- merge(brown_enriched.GO.05, brown.goseq, by="category") # change here
brown_enriched.GO.05 <- brown_enriched.GO.05[order(-brown_enriched.GO.05$numDEInCat),]
brown_enriched.GO.05$term <- as.factor(brown_enriched.GO.05$term)
head(brown_enriched.GO.05)

brown_MF <- subset(brown_enriched.GO.05, ontology=="MF")
brown_MF <- brown_MF[order(-brown_MF$numDEInCat),]
brown_CC <- subset(brown_enriched.GO.05, ontology=="CC")
brown_CC <- brown_CC[order(-brown_CC$numDEInCat),]
brown_BP <- subset(brown_enriched.GO.05, ontology=="BP")
brown_BP <- brown_BP[order(-brown_BP$numDEInCat),]

# Molecular Processes
MFplot_brown <- brown_MF %>% mutate(term = fct_reorder(term, numDEInCat)) %>%
  ggplot( aes(x=term, y=numDEInCat) ) +
  geom_segment( aes(x=term ,xend=term, y=0, yend=numDEInCat), color="grey") +
  geom_point(size=3, color="#69b3a2") +
  coord_flip() +
  theme(
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_blank(),
    legend.position="none"
  ) +
  xlab("") +
  ylab("") +
  ggtitle("Molecular Function") + #add a main title
  theme(plot.title = element_text(face = 'bold', 
                                  size = 12, 
                                  hjust = 0)) +
  theme_bw() + #Set background color 
  theme(panel.border = element_blank(), # Set border
        panel.grid.major = element_blank(), #Set major gridlines
        panel.grid.minor = element_blank(), #Set minor gridlines
        axis.line = element_line(colour = "black"), #Set axes color
        plot.background=element_blank())#Set the plot background

#Cellular Processes
CCplot_brown <- brown_CC %>% mutate(term = fct_reorder(term, numDEInCat)) %>%
  ggplot( aes(x=term, y=numDEInCat) ) +
  geom_segment( aes(x=term ,xend=term, y=0, yend=numDEInCat), color="grey") +
  geom_point(size=3, color="#69b3a2") +
  coord_flip() +
  theme(
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_blank(),
    legend.position="none"
  ) +
  xlab("") +
  ylab("") +
  ggtitle("Cellular Component") + #add a main title
  theme(plot.title = element_text(face = 'bold', 
                                  size = 12, 
                                  hjust = 0)) +
  theme_bw() + #Set background color 
  theme(panel.border = element_blank(), # Set border
        panel.grid.major = element_blank(), #Set major gridlines
        panel.grid.minor = element_blank(), #Set minor gridlines
        axis.line = element_line(colour = "black"), #Set axes color
        plot.background=element_blank())#Set the plot background

# Biological Processes 
BPplot_brown <- brown_BP %>% mutate(term = fct_reorder(term, numDEInCat)) %>%
  ggplot( aes(x=term, y=numDEInCat) ) +
  geom_segment( aes(x=term ,xend=term, y=0, yend=numDEInCat), color="grey") +
  geom_point(size=3, color="#69b3a2") +
  coord_flip() +
  theme(
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_blank(),
    legend.position="none") +
  xlab("") +
  ylab("") +
  ggtitle("Biological Process") + #add a main title
  theme(plot.title = element_text(face = 'bold', 
                                  size = 12, 
                                  hjust = 0)) +
  theme_bw() + #Set background color 
  theme(panel.border = element_blank(), # Set border
        panel.grid.major = element_blank(), #Set major gridlines
        panel.grid.minor = element_blank(), #Set minor gridlines
        axis.line = element_line(colour = "black"), #Set axes color
        plot.background=element_blank())#Set the plot background

# merge MF CC BP plots - Light cyan
num_brown   <- dim(GO_brown_genes)[1] # call num upregulated genes

library(tidyr)
pdf("Analysis/Output/WGCNA/Day7/goseq_modules/Day7_goseq_MEbrown.pdf", width=9, height=7)
brown_enriched.GO.05 %>% drop_na(ontology) %>% mutate(term = fct_reorder(term, numDEInCat)) %>%
  mutate(term = fct_reorder(term, ontology)) %>%
  ggplot( aes(x=term, y=numDEInCat) ) +
  geom_segment( aes(x=term ,xend=term, y=0, yend=numDEInCat), color="grey") +
  geom_point(size=3, aes(colour = ontology)) +
  coord_flip() +
  theme(
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_blank(),
    legend.position="bottom"
  ) +
  xlab("") +
  ylab("") +
  ggtitle("GO: Day 7 WGCNA MEbrown") +
  geom_label(aes(x = 2, y = 15, label = paste(num_brown, "MEbrown genes"))) +
  theme_bw() + #Set background color 
  theme(panel.border = element_blank(), # Set border
        panel.grid.major = element_blank(), #Set major gridlines
        panel.grid.minor = element_blank(), #Set minor gridlines
        axis.line = element_line(colour = "black"), #Set axes color
        plot.background=element_blank()) #Set the plot background #set title attributes
dev.off()
#===================================================================================== 
# It's goseq time!!!
# PLOTTING
#  red MODULE 
#===================================================================================== 

#Find only enriched GO terms that are statistically significant at cutoff 
red_enriched.GO.05.a<-red.goseq$category[red.goseq$over_represented_pvalue<.05] # change twice here
red_enriched.GO.05<-data.frame(red_enriched.GO.05.a)
colnames(red_enriched.GO.05) <- c("category")
red_enriched.GO.05 <- merge(red_enriched.GO.05, red.goseq, by="category") # change here
red_enriched.GO.05 <- red_enriched.GO.05[order(-red_enriched.GO.05$numDEInCat),]
red_enriched.GO.05$term <- as.factor(red_enriched.GO.05$term)
head(red_enriched.GO.05)

red_MF <- subset(red_enriched.GO.05, ontology=="MF")
red_MF <- red_MF[order(-red_MF$numDEInCat),]
red_CC <- subset(red_enriched.GO.05, ontology=="CC")
red_CC <- red_CC[order(-red_CC$numDEInCat),]
red_BP <- subset(red_enriched.GO.05, ontology=="BP")
red_BP <- red_BP[order(-red_BP$numDEInCat),]

# Molecular Processes
MFplot_red <- red_MF %>% mutate(term = fct_reorder(term, numDEInCat)) %>%
  ggplot( aes(x=term, y=numDEInCat) ) +
  geom_segment( aes(x=term ,xend=term, y=0, yend=numDEInCat), color="grey") +
  geom_point(size=3, color="#69b3a2") +
  coord_flip() +
  theme(
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_blank(),
    legend.position="none"
  ) +
  xlab("") +
  ylab("") +
  ggtitle("Molecular Function") + #add a main title
  theme(plot.title = element_text(face = 'bold', 
                                  size = 12, 
                                  hjust = 0)) +
  theme_bw() + #Set background color 
  theme(panel.border = element_blank(), # Set border
        panel.grid.major = element_blank(), #Set major gridlines
        panel.grid.minor = element_blank(), #Set minor gridlines
        axis.line = element_line(colour = "red"), #Set axes color
        plot.background=element_blank())#Set the plot background

#Cellular Processes
CCplot_red <- red_CC %>% mutate(term = fct_reorder(term, numDEInCat)) %>%
  ggplot( aes(x=term, y=numDEInCat) ) +
  geom_segment( aes(x=term ,xend=term, y=0, yend=numDEInCat), color="grey") +
  geom_point(size=3, color="#69b3a2") +
  coord_flip() +
  theme(
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_blank(),
    legend.position="none"
  ) +
  xlab("") +
  ylab("") +
  ggtitle("Cellular Component") + #add a main title
  theme(plot.title = element_text(face = 'bold', 
                                  size = 12, 
                                  hjust = 0)) +
  theme_bw() + #Set background color 
  theme(panel.border = element_blank(), # Set border
        panel.grid.major = element_blank(), #Set major gridlines
        panel.grid.minor = element_blank(), #Set minor gridlines
        axis.line = element_line(colour = "red"), #Set axes color
        plot.background=element_blank())#Set the plot background

# Biological Processes 
BPplot_red <- red_BP %>% mutate(term = fct_reorder(term, numDEInCat)) %>%
  ggplot( aes(x=term, y=numDEInCat) ) +
  geom_segment( aes(x=term ,xend=term, y=0, yend=numDEInCat), color="grey") +
  geom_point(size=3, color="#69b3a2") +
  coord_flip() +
  theme(
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_blank(),
    legend.position="none") +
  xlab("") +
  ylab("") +
  ggtitle("Biological Process") + #add a main title
  theme(plot.title = element_text(face = 'bold', 
                                  size = 12, 
                                  hjust = 0)) +
  theme_bw() + #Set background color 
  theme(panel.border = element_blank(), # Set border
        panel.grid.major = element_blank(), #Set major gridlines
        panel.grid.minor = element_blank(), #Set minor gridlines
        axis.line = element_line(colour = "red"), #Set axes color
        plot.background=element_blank())#Set the plot background

# merge MF CC BP plots - Light cyan
num_red   <- dim(GO_red_genes)[1] # call num upregulated genes

library(tidyr)
pdf("Analysis/Output/WGCNA/Day7/goseq_modules/Day7_goseq_MEred.pdf", width=9, height=7)
red_enriched.GO.05 %>% drop_na(ontology) %>% mutate(term = fct_reorder(term, numDEInCat)) %>%
  mutate(term = fct_reorder(term, ontology)) %>%
  ggplot( aes(x=term, y=numDEInCat) ) +
  geom_segment( aes(x=term ,xend=term, y=0, yend=numDEInCat), color="grey") +
  geom_point(size=3, aes(colour = ontology)) +
  coord_flip() +
  theme(
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_blank(),
    legend.position="bottom"
  ) +
  xlab("") +
  ylab("") +
  ggtitle("GO: Day 7 WGCNA MEred") +
  geom_label(aes(x = 2, y = 15, label = paste(num_red, "MEred genes"))) +
  theme_bw() + #Set background color 
  theme(panel.border = element_blank(), # Set border
        panel.grid.major = element_blank(), #Set major gridlines
        panel.grid.minor = element_blank(), #Set minor gridlines
        axis.line = element_line(colour = "black"), #Set axes color
        plot.background=element_blank()) #Set the plot background #set title attributes
dev.off()
