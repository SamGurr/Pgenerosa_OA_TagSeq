---
  # title: "DaysAll_WGCNA_TreatmentTimeseries"
  # author: "Samuel Gurr"
  # date: "January 8, 2021"
---

# LOAD PACKAGES
library(WGCNA) # note: this was previously installed with the command `BiocManager::install("WGCNA")`
library(dplyr)
library(zoo)
library(DESeq2)
library(BiocManager) # if (!requireNamespace("BiocManager", quietly = TRUE)); install.packages("BiocManager"); BiocManager::install(version = "3.13")
library(DESeq2) # BiocManager::install("DESeq2")

# for heatmap 
# library(devtools)
# install_github("jokergoo/ComplexHeatmap") first run these - commented out to avoid running a second time...
library(ComplexHeatmap)
library(circlize)
library(reshape)
library(ggplot2)
library(hrbrthemes)
library(purrr)


# SET WORKING DIRECTORY AND LOAD DATA
setwd("C:/Users/samjg/Documents/Github_repositories/Pgenerosa_TagSeq_Metabolomics/TagSeq/")
# LOAD DATA
# Tagaseq filtered counts - all counts on all days (filtered!)
# NOTE: these count matrices do not have the same genes, as they were filtered separately 
# not avoid omitting under-expressed genes due to lifestage (i.e. gene x on day 7 expressed high but day 21 low - if merged and filtered together it would be omitted; kept separate this gene is present in day 7 matrix!)
day7.counts.matrix  <- read.csv(file="Analysis/Data/Filtered_Counts/10cpm_50perc/day7.counts.filtered_10cpm50perc.csv", sep=',', header=TRUE)   %>%  dplyr::rename(Pgen_ID = X) %>%  mutate(Pgen_ID = as.character(Pgen_ID))
day14.counts.matrix <- read.csv(file="Analysis/Data/Filtered_Counts/10cpm_50perc/day14.counts.filtered_10cpm50perc.csv", sep=',', header=TRUE)  %>%  dplyr::rename(Pgen_ID = X) %>%  mutate(Pgen_ID = as.character(Pgen_ID))
day21.counts.matrix <- read.csv(file="Analysis/Data/Filtered_Counts/10cpm_50perc/day21.counts.filtered_10cpm50perc.csv", sep=',', header=TRUE)  %>%  dplyr::rename(Pgen_ID = X) %>%  mutate(Pgen_ID = as.character(Pgen_ID))

# how many genes (rows) are there in these count matrices??
nrow(day7.counts.matrix)  # 8548 total genes
nrow(day14.counts.matrix) # 8626 total genes
nrow(day21.counts.matrix) # 8421 total genes

# merge together the datasets to one! 
# NOTE: this will omit genes that did not pass the filter in 1 or 2 of the total datasets, meaning that ONLY genes present in ALL three filtered matrices will be present here 
Master.counts.matrix <- list(day7.counts.matrix,
                             day14.counts.matrix,
                             day21.counts.matrix) %>% purrr::reduce(inner_join, by = "Pgen_ID")
nrow(Master.counts.matrix) # 7933 - about 500-600 genes less than the other matrices 

# sanity check that this mege worked correctly....
(day7.counts.matrix %>% dplyr::filter(Pgen_ID %in% Master.counts.matrix$Pgen_ID) %>% dplyr::summarise(n())  ) == 7933  # should output TRUE
(day14.counts.matrix %>% dplyr::filter(Pgen_ID %in% Master.counts.matrix$Pgen_ID) %>% dplyr::summarise(n()) ) == 7933  # should output TRUE
(day21.counts.matrix %>% dplyr::filter(Pgen_ID %in% Master.counts.matrix$Pgen_ID) %>% dplyr::summarise(n()) ) == 7933  # should output TRUE

# now lets calculate the mean sd percent of carryover in these data sets
mean( c( (nrow(Master.counts.matrix) / nrow(day7.counts.matrix) *100),  # day 7  = 92.80533 %
         (nrow(Master.counts.matrix) / nrow(day14.counts.matrix) * 100), # day 14 = 91.96615 %
         (nrow(Master.counts.matrix) / nrow(day21.counts.matrix) * 100) )) # day 21 = 94.20496 %
sd( c(   (nrow(Master.counts.matrix) / nrow(day7.counts.matrix) *100),    # day 7 
         (nrow(Master.counts.matrix) / nrow(day14.counts.matrix) * 100),   # day 14
         (nrow(Master.counts.matrix) / nrow(day21.counts.matrix) * 100) )) # day 21 = 
# mean = 92.99215;sd =  1.131038 



# Treatment and Phenotype data
Master.Exp.Metadata                <- read.csv(file="Analysis/Data/Experiment_Metadata/Master_Phyenotype.and.Exp.Treatment_Metadata.csv", sep=',', header=TRUE) %>% 
                                            dplyr::select(c('Date', 'Sample.Name', 'Primary_Treatment')) %>% 
                                            dplyr::mutate(Date = factor(Date), Primary_Treatment = factor(Primary_Treatment))
levels(Master.Exp.Metadata$Date)   <- c("Day7", "Day14", "Day21") # change the date names to days 
Master.Exp.Metadata$Day_Treatment  <- paste(Master.Exp.Metadata$Date, Master.Exp.Metadata$Primary_Treatment, sep = "_")

# The following setting is important, do not omit.
options(stringsAsFactors = FALSE)

# ===================================================================================
#
#
#  PREPROCESSING THE DATA INPUT 
#
#  Read here for the pre processing steps using WGCNA!
#  https://bmcbioinformatics.biomedcentral.com/track/pdf/10.1186/1471-2105-9-559.pdf
# ===================================================================================


# filter low values - note: this is pre filteres for < 5 CPM in 50% of samples - less strict than the DESeq2 analysis
dim(Master.counts.matrix) # 7933  rows (genes) -   134  samples  not counting 'x' and 'Gene.ID'
(Master.counts.matrix)[1] # ommit the first line and transpose in the next line 
Master.data = as.data.frame(t(Master.counts.matrix[, -(1)])) # omit all columns but samples and transpose
dim(Master.data) # 133 7933
# fix(Master.data)

# trait data ========================================================== #

# Phenotype trait and Treatment data
dim(Master.Exp.Metadata) #  134   rows and 4 columns
names(Master.Exp.Metadata) # look at the 14 columns 
Master.Exp.Metadata_2 = Master.Exp.Metadata[, -1]; # remove columns that hold information we do not need. 
dim(Master.Exp.Metadata_2) # 134   rows and  3 columns
Master.Exp.Metadata_2 <- Master.Exp.Metadata_2 %>% dplyr::filter(!Sample.Name %in% 'SG92') # ommit the extra sample ID
# fix(d14.Treatment_Phenotype.data)

# count data  ========================================================== #

# fix(d21.data) # view the data - as you see all columns are now genes but filled as V1, V2, V3...
names(Master.data) = Master.counts.matrix$Pgen_ID # assigns column names (previous jsut numbered) as the gene ID 
rownames(Master.data) = names(Master.counts.matrix)[-(1)]; # assigns the row names as the sample ID
Master.data_matrix <- data.frame(Master.counts.matrix[,-1], row.names=Master.counts.matrix[,1]) 

Master.data_matrix_t <- t(Master.data_matrix)
fix(Master.data_matrix_t)


# creade dds object  ================================================== # (uses DESeq2 'DESeqDataSetFromMatrix' command)
ncol(Master.data_matrix) == nrow(Master.Exp.Metadata_2) # must be TRUE to proceed

# create dds objects to transform data 
dds.timeseries <- DESeqDataSetFromMatrix(countData = Master.data_matrix,
                                  colData = Master.Exp.Metadata_2, design = ~ 1) # DESeq Data Set (dds)
# DESeq Data Set (dds)
# NOTE: ~1 stands for no design; user will need to add a design for differential testing
# however for our purpose of just creating an onbject to transform, we do not need a design here...


# transform the data 
dds.timeseries_vst <- vst(dds.timeseries) # transform it vst
dds.timeseries_vst <- assay(dds.timeseries_vst) # call only the transformed coutns in the dds object
#fix(dds.d21_vst)
dds.timeseries_vst <- t(dds.timeseries_vst) # transpose columns to rows and vice versa

# ===================================================================================
#
# Day 21 WGCNA Sample tree 
#
# ===================================================================================

dim(dds.timeseries_vst) #  7933 genes;  133   samples

gsg = goodSamplesGenes(dds.timeseries_vst, verbose = 3);  # We first check for genes and samples with too many missing values:
gsg$allOK # If the statement returns TRUE, all genes have passed the cuts. 

# determine outlier samples 
sampleTree = hclust(dist(dds.timeseries_vst), method = "average") # Next we cluster the samples (in contrast to clustering genes that will come later)  to see if there are any obvious outliers.
sizeGrWindow(12,9) # The user should change the dimensions if the window is too large or too small.
par(cex = 0.6);
par(mar = c(0,4,2,0))

png("Analysis/Output/WGCNA/subseq_treatments_all/TimeSeries/Day7.14.21_ClusterTree_Precut.png", 1000, 1000, pointsize=20)
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5, 
     cex.axis = 1.5, cex.main = 2) + # appears there are two outliers SG55 and SG 105; can remove by hand or an automatic appraoch 
     abline(h = 100, col = "red") # add line to plot to show the cut-off od outlier samples (40000) SG105 and SG55
dev.off()

# complete the following if decided to cut ouliers from visualization of the ClusterTree (above)
# cut the tree and ommit from the reads data
clust = cutreeStatic(sampleTree, cutHeight = 100, minSize = 10) # Determine cluster under the line
table(clust) # 0 = cut; 1 = kept; says it will cut 2 and save 60
# 'keepsamples' boolean to call the main dataset
keepSamples = (clust==1) # call the 30 from clust at position 1
# view the ommited samples from the tree
omittedsamples = (clust==0)
omSamples = dds.timeseries_vst[omittedsamples, ]
nrow(omSamples) # should be 12
rownames(omSamples) # should be 12 - "SG141" "SG103" "SG54"  "SG126" "SG133" "SG52"  "SG115" "SG30"  "SG51"  "SG55"  "SG53"  "SG124" "SG29"  "SG110" "SG105" "SG95"  "SG114" What treatments are these??
omSamplesData <- Master.Exp.Metadata_2 %>% dplyr::filter(Sample.Name %in% (rownames(omSamples))) %>% dplyr::select(c('Sample.Name', 'Primary_Treatment', 'Day_Treatment')) # three Primary Ambient and three Primary Moderate
# integrate keepsamples to the main read dataset
dds.timeseries_vst = dds.timeseries_vst[keepSamples, ] # integreat the boolean 'keepsamples' to ommit oultilers determined in the sample tree above
nGenes = ncol(dds.timeseries_vst) # number of genes == 7933
nSamples = nrow(dds.timeseries_vst) # number of samples == 110  - the cut tree removed 12 samples

# plot the tree with the 'keep samples' bollean (T/F) to ommit outlier samples
sampleTree2 = hclust(dist(dds.timeseries_vst), method = "average") # Next we cluster the samples (in contrast to clustering genes that will come later)  to see if there are any obvious outliers.
png("Analysis/Output/WGCNA/subseq_treatments_all/TimeSeries/Day7.14.21_ClusterTree_Postcut.png", 1000, 1000, pointsize=20)
plot(sampleTree2, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2)
dev.off()


# based on outlier removall... call Trait data ===================================================== #
dim(dds.timeseries_vst) #  110 7933
dim(Master.Exp.Metadata_2) # 133   3 - trait data has  133 samples - need to cut to the timeseries_vst

# Form a data frame analogous to expression data that will hold the clinical traits.
TimeSeries.Samples = rownames(dds.timeseries_vst);# start new variable 'd21.Samples' calling the row names of the gene data (sample as 'Mouse' in trait data)
TreatRows = match(TimeSeries.Samples, Master.Exp.Metadata_2$Sample.Name); # match the names
TimeSeries.Traits = Master.Exp.Metadata_2[TreatRows, -1]; # removes the row numbers
rownames(TimeSeries.Traits) = Master.Exp.Metadata_2[TreatRows, 1]; # inserts the new traitRows - matches sample treatments
dim(TimeSeries.Traits) #  110   2 
all(rownames(TimeSeries.Traits) == rownames(dds.timeseries_vst))  # should be TRUE
dim(TimeSeries.Traits) #   110   2


# ===================================================================================
#
# Prepare Trait data  'Master.Exp.Metadata_2'
# ===================================================================================


TimeSeries.Traits # look at our trait data - just treatments

# primary treatment :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
Treatment_Primary <-  TimeSeries.Traits %>% dplyr::select('Primary_Treatment')

Treatment_Primary$A <- (Treatment_Primary$Primary_Treatment == "A")
Treatment_Primary$A   <- as.numeric(Treatment_Primary$A)

Treatment_Primary$M <- (Treatment_Primary$Primary_Treatment == "M")
Treatment_Primary$M   <- as.numeric(Treatment_Primary$M)

Treatment_Primary <- Treatment_Primary[,c(2:3)] # final dataset of 0,1 for teatment groups - Primary only!

# time (days)  :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
Days               <- TimeSeries.Traits %>% dplyr::select('Day_Treatment') 
Days$Day_Treatment <- sub("_.*", "", Days$Day_Treatment)

Days$Day7  <- (Days$Day_Treatment == "Day7")
Days$Day7  <- as.numeric(Days$Day7)

Days$Day14  <- (Days$Day_Treatment == "Day14")
Days$Day14  <- as.numeric(Days$Day14)

Days$Day21  <- (Days$Day_Treatment == "Day21")
Days$Day21  <- as.numeric(Days$Day21)

Days        <- Days[,c(2:4)] # final dataset of 0,1 for teatment groups - Primary only!
 
# group   :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
Day_Treat <-  TimeSeries.Traits %>% dplyr::select('Day_Treatment')


Day_Treat$Day7_A   <- (Day_Treat$Day_Treat == "Day7_A")
Day_Treat$Day7_A   <- as.numeric(Day_Treat$Day7_A)

Day_Treat$Day7_M   <- (Day_Treat$Day_Treat == "Day7_M")
Day_Treat$Day7_M   <- as.numeric(Day_Treat$Day7_M)



Day_Treat$Day14_A   <- (Day_Treat$Day_Treat == "Day14_A")
Day_Treat$Day14_A   <- as.numeric(Day_Treat$Day14_A)

Day_Treat$Day14_M   <- (Day_Treat$Day_Treat == "Day14_M")
Day_Treat$Day14_M   <- as.numeric(Day_Treat$Day14_M)



Day_Treat$Day21_A   <- (Day_Treat$Day_Treat == "Day21_A")
Day_Treat$Day21_A   <- as.numeric(Day_Treat$Day21_A)

Day_Treat$Day21_M   <- (Day_Treat$Day_Treat == "Day21_M")
Day_Treat$Day21_M   <- as.numeric(Day_Treat$Day21_M)

Day_Treat           <- Day_Treat[,c(2:7)] # final dataset of 0,1 for teatment groups - group groups only!


# ===================================================================================
#
# cluster samples by Trait
# ===================================================================================

# Primary treatment  :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
png("Analysis/Output/WGCNA/subseq_treatments_all/TimeSeries/Day7.14.21_ClusterTree_PrimaryTreatment.png", 1000, 1000, pointsize=20)
traitColors_Primary = labels2colors(Treatment_Primary); # Convert traits to a color representation: white means low, red means high, grey means missing entry
plotDendroAndColors(sampleTree2, traitColors_Primary, # Plot the sample dendrogram and the colors underneath.
                    groupLabels = names(traitColors_Primary), 
                    main = "Sample dendrogram and trait heatmap")
dev.off()

# Time/Days  :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
png("Analysis/Output/WGCNA/subseq_treatments_all/TimeSeries/Day7.14.21_ClusterTree_DaysTime.png", 1000, 1000, pointsize=20)
traitColors_Days = labels2colors(Days); # Convert traits to a color representation: white means low, red means high, grey means missing entry
plotDendroAndColors(sampleTree2, traitColors_Days, # Plot the sample dendrogram and the colors underneath.
                    groupLabels = names(traitColors_Days), 
                    main = "Sample dendrogram and trait heatmap")
dev.off()

# Group time and primary treatment (primary x time/days)  :::::::::::::::::::
png("Analysis/Output/WGCNA/subseq_treatments_all/TimeSeries/Day21_ClusterTree_DayTreatment.png", 1000, 1000, pointsize=20)
traitColors_Day_Treat = labels2colors(Day_Treat); # Convert traits to a color representation: white means low, red means high, grey means missing entry
plotDendroAndColors(sampleTree2, traitColors_Day_Treat, # Plot the sample dendrogram and the colors underneath.
                    groupLabels = names(traitColors_Day_Treat), 
                    main = "Sample dendrogram and trait heatmap")
dev.off()

# write the vst transformed data 
write.csv(dds.timeseries_vst, "Analysis/Output/WGCNA/subseq_treatments_all/TimeSeries/Days_7.14.21_vstTransformed_WGCNAdata.csv") # write

# ===================================================================================
#
# soft threshold
# ===================================================================================

dim(dds.timeseries_vst) #  133 7933

# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to=20, by=2))
# Call the network topology analysis function
sft = pickSoftThreshold(dds.timeseries_vst, powerVector = powers, verbose = 5) #...wait for this to finish
# pickSoftThreshold 
#  performs the analysis of network topology and aids the
# user in choosing a proper soft-thresholding power.
# The user chooses a set of candidate powers (the function provides suitable default values)
# function returns a set of network indices that should be inspected
sizeGrWindow(9, 5) # set window size 
# png to output 
png("Analysis/Output/WGCNA/subseq_treatments_all/TimeSeries/Days_7.14.21_ScaleTopology_SoftThresh.png", 1000, 1000, pointsize=20)

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
    
dev.off() # output soft threshold of 9


# The left panel... shows the scale-free fit index (y-axis) 
# as a function of the soft-thresholding power (x-axis). 
# We choose the power 9, which is the lowest power for which the scale-free 
# topology index curve attens out upon reaching a high value (in this case, roughly 0.90).
# The right panel.... displays the mean connectivity
# (degree, y-axis) as a function of the soft-thresholding power (x-axis).

#=====================================================================================
#
#  Start the step-wise module construction:  
# Step 1 = create adjacency matrix 
#
#=====================================================================================
softPower = 9 # set your soft threshold based on the plots above 

# signed 
adjacency_sign = adjacency(dds.timeseries_vst, power = softPower, type="signed") # this takes a long time.. just wait...

#=====================================================================================
#
#  Step 2: Turn adjacency into topological overlap
# Calculation of the topological overlap matrix, (TOM)
# and the corresponding dissimilarity, from a given adjacency matrix.
# Why signed??? Why unsigned???
# https://peterlangfelder.com/2018/11/25/signed-or-unsigned-which-network-type-is-preferable/
# https://www.rdocumentation.org/packages/WGCNA/subseq_treatments_all/versions/1.69/topics/adjacency
# https://ramellose.github.io/networktutorials/wgcna.html
# https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/subseq_treatments_all/TechnicalReports/signedTOM.pdf
#=====================================================================================

# TOM is used to 'reinforce' the direct connection between nodes. 
# in case of 'unsigned' networks, you can have triplet nodes in which the correlation signs are NOT reinforcing and goes UNADDRESSED 
# unisgned == loose information on the sign of underlying correlations - assumes the neigboring cannection can always be 
# considered "reinforcing" the direct connection; considered strong negative correlation as HIGH connection strength
# possible to result in connection between  triplets of nodes with dissimilar signs
# signed == considers negative correaltion between nodes as unconnected (conection strength near 0) 
# an attempt to take into account the cases when the neighboring node is NOT reinforcing 
# signed network the connection strength of nodes with negative correlation is zero so the situation (triplet with anti-reinforcing correlations) cannot arise 

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

plot(geneTree_sign, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity",
     labels = FALSE, hang = 0.04);

#=====================================================================================
#
#  Step 4: Set module size and 'cutreeDynamic' to create clusters 
#
#=====================================================================================
# We like large modules, so we set the minimum module size relatively high:
minModuleSize = 50; # set this for the subseqent call...
# Module identification using dynamic tree cut:

dynamicMods_sign = cutreeDynamic(dendro = geneTree_sign, distM = dissTOM_sign,
                                 deepSplit = 1, pamRespectsDendro = FALSE,
                                 minClusterSize = minModuleSize);
table(dynamicMods_sign) # number of genes per module - signed modules - more of them here all over 300 genes 

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
png("Analysis/Output/WGCNA/subseq_treatments_all/TimeSeries/Day_7.14.21_ClusterDendrogram_signedTOM.png", 1000, 1000, pointsize=20)
plotDendroAndColors(geneTree_sign, dynamicColors_sign, "Dynamic Tree Cut", 
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram and module colors 'SIGNED'") # signed TOM - now segmented well with minimal noise in major modules 
dev.off() # output the signed TOM dendrogram for comparison with the signed TOM (below and after merging - if needed!) 

#=====================================================================================
#
#  Step 6: Calculate Eigengenes - view thier connectivity based on 'MEDiss = 1-cor(MEs)'
#
#=====================================================================================
# Calculate eigengenes
# MEList = moduleEigengenes(dds.d21_vst, colors = dynamicColors)
MEList = moduleEigengenes(dds.timeseries_vst, colors = dynamicColors_sign) # signed TOM 
MEs = MEList$eigengenes # signed TOM  ASSINED EQIGENGENE VARIABE HERE
# Calculate dissimilarity of module eigengenes
MEDiss = 1-cor(MEs); # signed TOM 
# Cluster module eigengenes
METree = hclust(as.dist(MEDiss), method = "average");
# Plot the result
sizeGrWindow(7, 6)
plot(METree, main = "Clustering of module eigengenes (dissimilarity calc = MEDiss = 1-cor(MEs))",
     xlab = "", sub = "")
# change the cut height for the eignengene tree! 
MEDissThres = 0 # call the cut height
# Plot the cut line into the dendrogram
abline(h=MEDissThres, col = "red") # add the line to the tree

# save with the abline added
sizeGrWindow(7, 6)
png("Analysis/Output/WGCNA/subseq_treatments_all/TimeSeries/Day21_ClusterEigengenes_signedTOM.png", 1000, 1000, pointsize=20) # signed TOM 
plot(METree, main = "Clustering of module eigengenes (dissimilarity calc = MEDiss = 1-cor(MEs))",
     xlab = "", sub = "")  +
   abline(h=MEDissThres, col = "red") # add the line to the tree
dev.off()

#=====================================================================================
#
#  Step 7: Specify the cut line for the dendrogram (module) - Calc MODULE EIGENGENES (mergeMEs)
#
#=====================================================================================

# Call an automatic merging function
merge = mergeCloseModules(dds.timeseries_vst, dynamicColors_sign, cutHeight = MEDissThres, verbose = 3) # signed TOM 
# merge = mergeCloseModules(dds.d14_vst, dynamicColors_sign, cutHeight = 0.3, verbose = 3) # signed TOM 
# important! I intentially called cutheaight below the modules (0.3 above) to see the downstream trajetory 
# of these modules before shrunk down. Look below for the dendrogram including the raw and the merged modules 
# appears the module merging is unneccessary considering the low number of moduales already 

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
png("Analysis/Output/WGCNA/subseq_treatments_all/TimeSeries/Day_7.14.21_ClusterDendrogram_signed.png", 1000, 1000, pointsize=20)
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
# moduleColors = mergedColors # if you want to call ther merged data (view dendrogram above)
moduleColors = dynamicColors_sign # if you DO NOT want to merged colors 
# Construct numerical labels corresponding to the colors
colorOrder = c("grey", standardColors(50));
moduleLabels = match(moduleColors, colorOrder)-1;

MEs = mergedMEs # CALL THIS IF YOU WANT TO MERGED MODULE EIGENGENS - IF NOW LEAVE COMMMENTED OUT



# Save module colors and labels for use in subsequent parts
# save(MEs, moduleLabels, moduleColors, geneTree, file = "Analysis/Output/WGCNA/subseq_treatments_all/Day14/Day14-networkConstruction-stepByStep.RData")
# write csv - save the module eigengenes
write.csv(MEs, file = "Analysis/Output/WGCNA/subseq_treatments_all/TimeSeries/d_7.14.21.WGCNA_ModulEigengenes.csv") 

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

# https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/subseq_treatments_all/TechnicalReports/signedTOM.pdf
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


#=====================================================================================
#
#  Prepare for  module trait associations - trait data as factors
#
#=====================================================================================
# identify modules that are signiFcantly
# associated with the measured clinical traits.

# Since we already have a summary profile (eigengene) for each module,
# we simply correlate eigengenes with external traits and look for the most signicant associations:

# Define numbers of genes and samples
nGenes = ncol(dds.timeseries_vst);   # 7933
nSamples = nrow(dds.timeseries_vst); #  133

# change chanracter treatments to integers
# ALL TRAIT DATA
TimeSeries.Traits$Primary_Treatment <- as.factor(TimeSeries.Traits$Primary_Treatment)
TimeSeries.Traits$Primary_Treatment <- as.numeric(TimeSeries.Traits$Primary_Treatment)

TimeSeries.Traits$Day  <- as.factor( sub("_.*","",TimeSeries.Traits$Day_Treatment) )
TimeSeries.Traits$Day  <- as.numeric(TimeSeries.Traits$Day )

TimeSeries.Traits$Day_Treatment <- as.factor(TimeSeries.Traits$Day_Treatment)
TimeSeries.Traits$Day_Treatment <- as.numeric(TimeSeries.Traits$Day_Treatment)

#=====================================================================================
#
# Module trait correlation
#
#=====================================================================================

# Treatment_Primary 
moduleTraitCor_Primary = cor(MEs, Treatment_Primary, use = "p");
moduleTraitPvalue_Primary = corPvalueStudent(moduleTraitCor_Primary, nSamples);
# Day
moduleTraitCor_Days = cor(MEs, Days, use = "p");
moduleTraitPvalue_Days = corPvalueStudent(moduleTraitCor_Days, nSamples);
# Day_Treat
moduleTraitCor_Day_Treat = cor(MEs, Day_Treat, use = "p");
moduleTraitPvalue_Day_Treat = corPvalueStudent(moduleTraitCor_Day_Treat, nSamples);


#=====================================================================================
#
# Module eigengene -  MEs boxplots by treatment group
#
#=====================================================================================

TimeSeries.Traits <- tibble::rownames_to_column(TimeSeries.Traits, "Sample.Name")

library(reshape2)
MEs_table <- MEs # new table for plotting 
MEs_table$Sample.Name <- rownames(MEs) # call rows as coolumn to merge with treatment data
MEsPlotting <- merge(TimeSeries.Traits, MEs_table, by = 'Sample.Name') # merge
#MEsPlotting <- MEsPlotting[,-c(4:7)] # ommit the phys data to just plot the module colors 
MEsPlotting_melt <- melt(MEsPlotting, id=c('Sample.Name', 'Primary_Treatment', 'Day', 'Day_Treatment'))
#plot it
png("Analysis/Output/WGCNA/subseq_treatments_all/TimeSeries/Day_7.14.21_ME_Boxplot.png", 600, 1000, pointsize=20)
ggplot(MEsPlotting_melt, aes(x=Day, y=value, fill = factor(Primary_Treatment))) +
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
# Heatmaps
#
#=====================================================================================



# PRRIMARY TREAMENT ONLY  ------------------------------------------------------------------ # 


sizeGrWindow(10,10)
# Will display correlations and their p-values
PRIMARYTreatments.matrix <-  paste(signif(moduleTraitCor_Primary, 2), "\n(",
                                       signif(moduleTraitPvalue_Primary, 1), ")", sep = "")
#dim(textMatrix) == dim(moduleTraitCor_treatonly)
par(mar = c(8, 9.5, 5, 3));
# Display the correlation values within a heatmap plot
png("Analysis/Output/WGCNA/subseq_treatments_all/TimeSeries/heatmaps/Day2_7.14.21_Treatments_Primary_heatmap2.png", 500, 1000, pointsize=20)
labeledHeatmap(Matrix = moduleTraitCor_Primary,
               xLabels = names(Treatment_Primary),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = TRUE,
               colors = blueWhiteRed(50),
               textMatrix = PRIMARYTreatments.matrix,
               setStdMargins = FALSE,
               cex.text = 1,
               zlim = c(-0.6,0.6),
               main = paste("Module-trait relationships - Primary Treatment"))
dev.off()

# this heatmap looks better
PRIMARYtreatment.text <-  as.matrix(signif(moduleTraitPvalue_Primary, 3))
pa = cluster::pam(PRIMARYtreatment.text, k = 4)
col_fun = colorRamp2(c(-0.5, 0, 0.5), c("blue", "white", "red"))
# png("Analysis/Output/WGCNA/subseq_treatments_all/TimeSeries/heatmaps/Day21_Treatments_Primary_heatmap.png", 500, 1000, pointsize=20)
pdf("Analysis/Output/WGCNA/subseq_treatments_all/TimeSeries/heatmaps/Day21_Treatments_Primary_heatmap.pdf", width=5, height=6)
Heatmap(moduleTraitCor_Primary, 
        name = "gene_cor", 
        rect_gp = gpar(col = "grey", lwd = 1),
        column_title = "Day 21 WGCNA - Primary Treatment", 
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
          grid.text(sprintf("%.1f", PRIMARYtreatment.text[i, j]), x, y, gp = gpar(fontsize = 10))
        })
dev.off()




# Day_Treat group ------------------------------------------------------------------ # 
moduleTraitCor_Day_Treat = cor(MEs, Day_Treat, use = "p");
moduleTraitPvalue_Day_Treat = corPvalueStudent(moduleTraitCor_Day_Treat, nSamples);

sizeGrWindow(10,10)
# Will display correlations and their p-values
Day_Treat.matrix <-  paste(signif(moduleTraitCor_Day_Treat, 2), "\n(",
                                      signif(moduleTraitPvalue_Day_Treat, 3), ")", sep = "")
par(mar = c(8, 9.5, 5, 3));
# Display the correlation values within a heatmap plot
png("Analysis/Output/WGCNA/subseq_treatments_all/TimeSeries/heatmaps/DayTreat_Group_heatmap2.png", 500, 1000, pointsize=20)
labeledHeatmap(Matrix = moduleTraitCor_Day_Treat,
               xLabels = names(Day_Treat),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = TRUE,
               colors = blueWhiteRed(50),
               textMatrix = Day_Treat.matrix,
               setStdMargins = FALSE,
               cex.text = 1,
               zlim = c(-0.6,0.6),
               main = paste("Module-trait relationships - Third Treatment"))
dev.off()

# this heatmap looks better
Day_Treat_text<-  as.matrix(signif(moduleTraitPvalue_Day_Treat, 4))
Day_Treat_cor <-  as.matrix(signif(moduleTraitCor_Day_Treat, 4))
pa = cluster::pam(Day_Treat_cor, k = 3)
col_fun = colorRamp2(c(-0.5, 0, 0.5), c("blue", "white", "red"))
# png("Analysis/Output/WGCNA/subseq_treatments_all/TimeSeries/heatmaps/Day21_Treatment_Third_heatmap.png", 500, 1000, pointsize=20)
pdf("Analysis/Output/WGCNA/subseq_treatments_all/TimeSeries/heatmaps/DayTreat_Group_heatmap.pdf", width=5, height=6)
Heatmap(moduleTraitCor_Day_Treat, 
        name = "gene_cor", 
        rect_gp = gpar(col = "grey", lwd = 1),
        column_title = "Timeseries WGCNA - Day_Treat", 
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
          grid.text(sprintf("%.1f", Day_Treat_text[i, j]), x, y, gp = gpar(fontsize = 10))
        })
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

geneModuleMembership = as.data.frame(cor(dds.timeseries_vst, MEs, use = "p"));
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));

names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");


# AA treatment group
Day7_A = as.data.frame(Day_Treat$Day7_A) # Define variable containing the desired column 
names(Day7_A) = "Day7_A"
Day7_A_geneTraitSignificance = as.data.frame(cor(dds.timeseries_vst, Day7_A, use = "p"));
Day7_A_GSPvalue = as.data.frame(corPvalueStudent(as.matrix(Day7_A_geneTraitSignificance), nSamples));
names(Day7_A_geneTraitSignificance) = paste("GS.", names(Day7_A), sep=""); # MA_geneTraitSignificance - pearsons correlation between reads and the MA grop
names(Day7_A_GSPvalue) = paste("p.GS.", names(Day7_A), sep=""); # corPvalueStudent


#=====================================================================================
#
#   COUNT GENES OF INTEREST IN  MODULES (i.e. violet- refer to heatmap)
#
#=====================================================================================

length(colnames(dds.timeseries_vst)[moduleColors=="brown"]) # 1072 total genes in the brown module
length(colnames(dds.timeseries_vst)[moduleColors=="yellow"]) # 957 total genes in the yellow module

#=====================================================================================
#
#  Call annotation data to get module gene data (prep for downstream GO)
#
#=====================================================================================
annot = read.delim2(file = "Seq_details/Panopea-generosa-genes-annotations.txt",header = F);
dim(annot) # 34947     9
names(annot) # V1 are the gene.IDs
probes = names(Master.data)
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
names(Day7_A_geneTraitSignificance)
names(Day7_A_GSPvalue)
geneInfo_GROUPS = data.frame(substanceBXH = probes,
                             geneSymbol = annot$V1[probes2annot],
                             #LocusLinkID = annot$LocusLinkID[probes2annot],
                             moduleColor = moduleColors,
                             Uniprot = annot$V5[probes2annot],
                             HGNC = annot$V6[probes2annot],
                             GO.terms = annot$V8[probes2annot],
                             GO.description = annot$V9[probes2annot],
                             Day7_A_geneTraitSignificance, Day7_A_GSPvalue)
                             # AA_geneTraitSignificance, AM_geneTraitSignificance, AS_geneTraitSignificance,
                             # MA_geneTraitSignificance, MM_geneTraitSignificance, MS_geneTraitSignificance, # call this specific to the module and trait of interest
                             # AA_GSPvalue,  AM_GSPvalue,  AS_GSPvalue,
                             # MA_GSPvalue,  MM_GSPvalue,  MS_GSPvalue)              # call this specific to the module and trait of interest
View(geneInfo_GROUPS)
modOrder = order(-abs(cor(MEs, Day7_A, use = "p"))); # order by the strength of the correlation between module and trait values for each sample

for (mod in 1:ncol(geneModuleMembership)) # Add module membership information in the chosen order
{
  oldNames = names(geneInfo_GROUPS)
  geneInfo_GROUPS = data.frame(geneInfo_GROUPS, geneModuleMembership[, modOrder[mod]], 
                               MMPvalue[, modOrder[mod]]);
  names(geneInfo_GROUPS) = c(oldNames, paste("MM.", modNames[modOrder[mod]], sep=""),
                             paste("p.MM.", modNames[modOrder[mod]], sep=""))
}
# Order the genes in the geneInfo variable first by module color, then by geneTraitSignificance
geneOrder = order(geneInfo_GROUPS$moduleColor, -abs(geneInfo_GROUPS$GS.Day7_A));
geneInfo_GROUPS = geneInfo_GROUPS[geneOrder, ]
View(geneInfo_GROUPS)



#=====================================================================================
#
#  PLOTTING LINE GRAPHS BY MODULE COLOR 
#
#=====================================================================================

# view the Treatment groups heatmap to understand rationale for clusters 

# green module 
length(colnames(dds.d21_vst)[moduleColors=="green"]) # 595 total genes in the green module
green_mod <- geneInfo_GROUPS %>%  dplyr::filter(moduleColor %in% 'green')
nrow(green_mod) # 630 total genes in the green module
green_mod_05 <- green_mod %>%  dplyr::filter(p.MM.green < 0.05) 
green_mod_05_MM <-green_mod_05[,c(8:13)]
nrow(green_mod_05_MM) # 528 total genes in the green module

green_mod_05_MM$Gene.ID <- rownames(green_mod_05_MM)
green_mod_05_MM_melt <- melt(green_mod_05_MM, id=c("Gene.ID"))
green_mod_05_MM_melt$Primary <- substr(green_mod_05_MM_melt$variable, 4,4)
green_mod_05_MM_melt$Second <- substr(green_mod_05_MM_melt$variable, 5,5)

ggplot(green_mod_05_MM_melt, aes(x=Second, y=value, group=Gene.ID, color = Primary)) +
  theme(panel.grid=element_blank()) +
 # scale_color_manual(values=c("#56B4E9","#D55E00")) +
  geom_line(size=0.2, alpha=0.1) +
  theme_bw() +
  facet_wrap(~Primary)

#=====================================================================================
#
#  Write csv for the modules and corresponding raw read counts
#
#=====================================================================================
# call the module of interest for follow-up GO analysis 

write.csv(geneInfo_GROUPS, file = "Analysis/Output/WGCNA/subseq_treatments_all/TimeSeries/TimeSeries.WGCNA_ModulMembership.csv")

#=====================================================================================
#
#  LOAD DATA for the following
#  Plot (1) Linear regressions of eigengenes for each module
#       (2) Eigengen exp by treatment  for modules of interest (review heatmaps)
#       (3) Go terms for modules of interest (review heatmaps)
#=====================================================================================

# SET WORKING DIRECTORY AND LOAD DATA
setwd("C:/Users/samjg/Documents/My_Projects/Github_repositories/TagSeq/")

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
library(ggpmisc)
# Load data 
ModEigen <- read.csv("Analysis/Output/WGCNA/subseq_treatments_all/TimeSeries/d_7.14.21.WGCNA_ModulEigengenes.csv")
Annot_ModuleMembership <-  read.csv("Analysis/Output/WGCNA/subseq_treatments_all/TimeSeries/TimeSeries.WGCNA_ModulMembership.csv")
vst_data <- read.csv("Analysis/Output/WGCNA/subseq_treatments_all/TimeSeries/Days_7.14.21_vstTransformed_WGCNAdata.csv")
Master.Treatment.data <- read.csv(file="Analysis/Data/Experiment_Metadata/Master_Phyenotype.and.Exp.Treatment_Metadata.csv", sep=',', header=TRUE)



#=====================================================================================
#
#  Prep and merge datasets for plotting
# 
#=====================================================================================
# Prep data to merge
# module membership
Annot_ModuleMembership <- Annot_ModuleMembership[-c(1:2)] # ommit the redundant gene name columns buut keep 'geneSymbol'
dim(Annot_ModuleMembership) # 7933   24

# module eigengenes
names(ModEigen)[1] <- "Sample.Name"

# normalized count data (same used for the dds object and WGCNA analysis)
vst_data_t <- as.data.frame(t(vst_data[, -(1)])) # trsnpose columns names to rows (genes) 
colnames(vst_data_t) <- vst_data[, 1] # name the columns as the rows in previous dataset (sample IDS)
vst_data_t<- vst_data_t %>% tibble::rownames_to_column("geneSymbol") # "geneSymbol" - create column and remname the gene name column to 'geneSymbol'


# merge Master datasets
# GO terms
GOTermsMaster.Modules <-  merge(Annot_ModuleMembership, vst_data_t, by = "geneSymbol")
dim(GOTermsMaster.Modules) #7933  134
# Eigengenes and traits
EigenTraitMaster<-  merge(ModEigen, Master.Treatment.data, by = "Sample.Name")
dim(EigenTraitMaster) # 110  22




#=====================================================================================
# 
# Linear regression plots Eigenengene expression with Traits for EACH module       
#
#=====================================================================================

# lets loop this through to get a bunch of plots 
phys <- EigenTraitMaster[,c(18:23)]
modules <- EigenTraitMaster[,c(2:10)]


for(i in 1:ncol(modules)) {
  module_color <- substr(names(modules)[i], 3,12)
  moduledata <- modules[,i]
  
  par(mfrow=c(3,2))
  
  my.formula <- y ~ x
  
  p <- ggplot(d21_EigenTraitMaster, 
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
  
  
  p <- ggplot(d21_EigenTraitMaster, 
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
  
  
  p <- ggplot(d21_EigenTraitMaster, 
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
  
  
  p <- ggplot(d21_EigenTraitMaster,
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


  p <- ggplot(d21_EigenTraitMaster,
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

  p <- ggplot(d21_EigenTraitMaster,
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
  #png("Analysis/Output/WGCNA/subseq_treatments_all/TimeSeries/Day21_PhysOnly_heatmap2.png", 1000, 1000, pointsize=20)
  ggarrange(p1, p2, p3, 
            p4, p5, p6, 
            ncol = 2, nrow = 3)
  ggsave(paste0("Analysis/Output/WGCNA/subseq_treatments_all/TimeSeries/EigengenePlots/day21_EigengenePlot_", module_color, ".png"))
}



#===================================================================================== 
# 
# EXPLORE THE Expression of each module (for loop plots!) BY TREATMENT
#
# evidence from the heatmaps and the regression shwo strong module membership and eigenenge cor strength with lighcyan 
# in response to treatment (primary specifically) and Total Antioxidant capacity 
#
#=====================================================================================

modcolor <- as.data.frame(unique(d21_Annot_ModuleMembership$moduleColor))
names(modcolor)[1] <- "color"

# experiment treatment and total protein data - narrow the columns 
exp.phys_data <- Master.Treatment_Phenotype.data %>% 
  dplyr::filter(Date %in% 20190814) # filter out Day 7 data only 

for(i in 1:nrow(modcolor)) {
  
  # vst read count date - narrow the columns - reshape and rename
  Mod_geneIDs <- d21_Annot_ModuleMembership %>% dplyr::filter(moduleColor %in% modcolor[i,]) %>%  dplyr::select("geneSymbol")
  d21_vst_Mod <- d21_vst_data_t %>% dplyr::filter(geneSymbol %in% Mod_geneIDs[,1])
  d21_vst_Mod_MELT <- melt(d21_vst_Mod, id=("geneSymbol")) # melt using reshape2
  names(d21_vst_Mod_MELT)[(2:3)] <- c('Sample.Name', 'vst_Expression') # change column names
  
  # merge by common row values 'Sample.Name'
  merged_Expdata_Mod <- merge(d21_vst_Mod_MELT, exp.phys_data, by ='Sample.Name')
  
  # mean Exp response table 
  meanEXp_Mod <- merged_Expdata_Mod %>% 
    select(c('Sample.Name','vst_Expression','Primary_Treatment', 'Second_Treament', 'Third_Treatment', 'All_Treatment')) %>% 
    group_by(Sample.Name, Primary_Treatment, Second_Treament, Third_Treatment, All_Treatment) %>%
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
  # All treatment
  meanEXp_Summary.Group_Mod <- meanEXp_Mod %>% 
    group_by(Third_Treatment,Second_Treament,Primary_Treatment) %>%
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
    ggtitle(paste("Day 14 WGCNA", modcolor[i,], "Module VST GeneExp", sep =' ')) +
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
  
  # Group treatment mean sd plot
  min_p3 <- min(meanEXp_Summary.Group_Mod$mean) - max(meanEXp_Summary.Group_Mod$se)
  max_p3 <- max(meanEXp_Summary.Group_Mod$mean) + max(meanEXp_Summary.Group_Mod$se)
  
  Group.vst.Mod <- ggplot(meanEXp_Summary.Group_Mod, aes(x=Third_Treatment, y=mean, group=Primary_Treatment, fill = Primary_Treatment)) +
    theme_classic() +
    geom_errorbar(aes(ymin=mean-se, ymax=mean+se), colour="black", width=.1, position=pd) +
    geom_line(position=pd) +
    geom_point(position=pd, size = 4, shape=21) +            
    xlab("Third pCO2 treatment") +
   # ylab(paste(modcolor[i,]," Module VST Gene Expression (Mean +/- SE)", sep = ' ')) +                 # note the mean was first by sample ID THEN by treatment
    scale_fill_manual(values=c("#56B4E9","#E69F00")) +
    # scale_color_manual(values=c("#56B4E9","#E69F00")) +
    # ggtitle("Day 21 WGCNA red' Module VST GeneExp") +
    # expand_limits(y=0) +                                                    # Expand y range
    scale_y_continuous(limits=c((min_p3), (max_p3))) +
    theme(text = element_text(size=15)) +
    facet_wrap(~Second_Treament)
  
  # output 
  
  pdf(paste("Analysis/Output/WGCNA/subseq_treatments_all/TimeSeries/ModuleExpression_Treatment/day21_Exp_Module",modcolor[i,],".pdf"), width=6, height=12)
  print(ggarrange(Primary.vst.Mod, Sec.vst.Mod, Group.vst.Mod,        
                  plotlist = NULL,
                  ncol = 1,
                  nrow = 3,
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
  Mod_geneIDs <- d14_Annot_ModuleMembership %>% dplyr::filter(moduleColor %in% module_color) %>%  dplyr::select("geneSymbol") # all geneSymbols (IDs) WITHIN the looped module color 
  d14_vst_Mod <- d14_vst_data_t %>% dplyr::filter(geneSymbol %in% Mod_geneIDs[,1]) # call the transposed vst count data and filer for the module color 
  d14_vst_Mod_MELT <- melt(d14_vst_Mod, id=("geneSymbol")) # melt using reshape2 - ('un-transpose')
  names(d14_vst_Mod_MELT)[(2:3)] <- c('Sample.Name', 'vst_Expression') # change column names to merge next
  
  # merge by common row values 'Sample.Name'
  merged_Expdata_Mod <- merge(d14_vst_Mod_MELT, exp.phys_data, by ='Sample.Name') # merge
  
  # mean Exp response table 
  meanExp_Mod <- merged_Expdata_Mod %>% 
    group_by(Sample.Name) %>%
    dplyr::summarize(n = n(), 
                     mean.vstExp = mean(vst_Expression), 
                     sd.vsdtExp = sd(vst_Expression),
                     se.vsdtExp = sd.vsdtExp/sqrt(n)) # mean expression for each sample ID 
  
  merged_Expdata_phystreat <- merged_Expdata_Mod %>%  # condence down the 'merged_Expdata_Mod' file to merge with 'meanExp_Mod'
    dplyr::select(c("Sample.Name", "Primary_Treatment", "Second_Treament","mean.Resp.ugLhrindiv", "mean.Biovol", "mean.Shell.Length"))
  
  meanExp_Mod.MASTER <- unique(merge(meanExp_Mod, merged_Expdata_phystreat, by = "Sample.Name")) # call only unqiue values 
  
  phys <- meanExp_Mod.MASTER[, 4:6]  # call only the phys columns for the subseqent for loop - plotting each 
  
  dir.create(paste0("Analysis/Output/WGCNA/subseq_treatments_all/day14/ModuleExpresion_Phenotype/",module_color), showWarnings = FALSE)
  
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
    
    png(paste("Analysis/Output/WGCNA/subseq_treatments_all/day14/ModuleExpresion_Phenotype/",module_color,"/day14_Exp_Module_",module_color, "_",(names(phys)[m]),".png", sep = ''), 600, 1000, pointsize=20)
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
Geoduck_annotation <- read.delim2(file="C:/Users/samjg/Documents/My_Projects/Github_repositories/TagSeq/Seq_details/Panopea-generosa-genes-annotations.txt", header=F)


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



#==============================================================================
#
#
#   FOR LOOP THE REST of goseq 
# 
#
#==============================================================================
modcolor <- as.data.frame(unique(d21_Annot_ModuleMembership$moduleColor))
names(modcolor)[1] <- "color"


for(i in 1:nrow(modcolor)) {
  
  # call the genes or goseq
  GO_module <- d21_Annot_ModuleMembership %>% dplyr::filter(moduleColor %in% modcolor[i,]) # %>%  dplyr::select("geneSymbol")
  GO_module_genes <- GO_module[1]
  names(GO_module_genes)[1] <- "Gene.ID" # 162 genws in the green module 
  # convert to integer with all unique genes
  GO_module_integer <- as.integer(GO_unique.genes.all%in%(GO_module_genes$Gene.ID)) # Day 0 - Primary - Upregulated DEGs
  names(GO_module_integer)=GO_unique.genes.all
  
  #Calculate Probability Weighting Function (using 'nullp')
  pwf <-nullp(GO_module_integer, GO_unique.genes.all, bias.data=length_vector) #weight vector by length of gene
  # Run goseq
  mod.goseq  <-goseq(pwf, ID.vector, gene2cat=GO.terms, test.cats=c("GO:CC", "GO:BP", "GO:MF"), method="Wallenius", use_genes_without_cat=TRUE)
  # call enriched GO terms and plot 
  #Find only enriched GO terms that are statistically significant at cutoff 
  mod_enriched.GO.05.a<-mod.goseq$category[mod.goseq$over_represented_pvalue<.05] # change twice here
  mod_enriched.GO.05<-data.frame(mod_enriched.GO.05.a)
  colnames(mod_enriched.GO.05) <- c("category")
  mod_enriched.GO.05 <- merge(mod_enriched.GO.05, mod.goseq, by="category") # change here
  mod_enriched.GO.05 <- mod_enriched.GO.05[order(-mod_enriched.GO.05$numDEInCat),]
  mod_enriched.GO.05$term <- as.factor(mod_enriched.GO.05$term)
  head(mod_enriched.GO.05)
  
  mod_MF <- subset(mod_enriched.GO.05, ontology=="MF")
  mod_MF <- mod_MF[order(-mod_MF$numDEInCat),]
  mod_CC <- subset(mod_enriched.GO.05, ontology=="CC")
  mod_CC <- mod_CC[order(-mod_CC$numDEInCat),]
  mod_BP <- subset(mod_enriched.GO.05, ontology=="BP")
  mod_BP <- mod_BP[order(-mod_BP$numDEInCat),]

  # merge MF CC BP plots - Light cyan
  num_mod   <- dim(GO_module_genes)[1] # call num upregulated genes
  
  # plot the output
  library(tidyr)
  Plot <- mod_enriched.GO.05 %>% drop_na(ontology) %>% mutate(term = fct_reorder(term, (-log10(over_represented_pvalue)) )) %>%
    dplyr::filter(ontology %in% c('MF', 'BP')) %>% 
    mutate(term = fct_reorder(term, ontology)) %>%
    ggplot( aes(x=term, y=(-log10(over_represented_pvalue)) ) ) +
    geom_segment( aes(x=term ,xend=term, y=0, yend=(-log10(over_represented_pvalue)), colour = ontology), size = 2, lineend = "butt", alpha = 0.1) + #
    geom_point(size=2, shape = 15, aes(colour = ontology)) +
    coord_flip() +
    theme(
      panel.grid.minor.y = element_blank(),
      panel.grid.major.y = element_blank(),
      legend.position="bottom"
    ) +
    xlab("") +
    ylab("-log10(pvalue)") +
    ggtitle(paste("GO: Day 21 WGCNA ME", modcolor[i,], sep = "")) +
    geom_label(aes(x = 2, y = 4, label = paste(num_mod, "ME", modcolor[i,], "genes"))) +
    theme_bw() + #Set background color 
    theme(panel.border = element_blank(), # Set border
          panel.grid.major = element_blank(), #Set major gridlines
          panel.grid.minor = element_blank(), #Set minor gridlines
          axis.line = element_line(colour = "black"), #Set axes color
          axis.text.x = element_text(size = 12), # set size of x axis text (log10 vals)
          axis.title.x = element_text(size = 14), # set size of x axis title (log10 vals)
          axis.text.y = element_text(size = 8), # set size of y axis text (terms)
          #legend.position = c(30, 30), 
          plot.background=element_blank()) + #Set the plot background #set title attributes 
    geom_hline(yintercept = 1.3, linetype="dashed", 
               color = "grey", size=1)
  
  pdf(paste("Analysis/Output/WGCNA/subseq_treatments_all/TimeSeries/goseq_modules/Day21_goseq_ME",modcolor[i,],"mod.pdf", sep =''), width=8, height=5)
  print(Plot)
  dev.off()
}


#===================================================================================== 
# It's goseq time!!!
# PLOTTING
#  CLUSTER! BLACK AND PINKK MODULES 
#===================================================================================== 
# Example: 0 = not in moddulecolor; 1 = in modulecolor
magenta.blue_module <- d14_Annot_ModuleMembership %>% dplyr::filter(moduleColor %in% c('magenta', 'blue')) # %>%  dplyr::select("geneSymbol")
magenta.blue_module_genes <- magenta.blue_module[1]
names(magenta.blue_module_genes)[1] <- "Gene.ID" # 162 genws in the green module 

# convert to integer with all unique genes
magenta.blue_module_integer <- as.integer(GO_unique.genes.all%in%(magenta.blue_module_genes$Gene.ID)) # Day 0 - Primary - Upregulated DEGs
names(magenta.blue_module_integer)=GO_unique.genes.all

#Calculate Probability Weighting Function (using 'nullp')
pwf_magenta.blue <-nullp(magenta.blue_module_integer, GO_unique.genes.all, bias.data=length_vector) #weight vector by length of gene


# Run goseq
magenta.blue.goseq  <-goseq(pwf_magenta.blue, ID.vector, gene2cat=GO.terms, test.cats=c("GO:CC", "GO:BP", "GO:MF"), method="Wallenius", use_genes_without_cat=TRUE)

#Find only enriched GO terms that are statistically significant at cutoff 
magenta.blue_enriched.GO.05.a<-magenta.blue.goseq$category[magenta.blue.goseq$over_represented_pvalue<.05] # change twice here
magenta.blue_enriched.GO.05<-data.frame(magenta.blue_enriched.GO.05.a)
colnames(magenta.blue_enriched.GO.05) <- c("category")
magenta.blue_enriched.GO.05 <- merge(magenta.blue_enriched.GO.05, magenta.blue.goseq, by="category") # change here
magenta.blue_enriched.GO.05 <- magenta.blue_enriched.GO.05[order(-magenta.blue_enriched.GO.05$numDEInCat),]
magenta.blue_enriched.GO.05$term <- as.factor(magenta.blue_enriched.GO.05$term)
head(magenta.blue_enriched.GO.05)

magenta.blue_MF <- subset(magenta.blue_enriched.GO.05, ontology=="MF")
magenta.blue_MF <- magenta.blue_MF[order(-magenta.blue_MF$numDEInCat),]
magenta.blue_CC <- subset(magenta.blue_enriched.GO.05, ontology=="CC")
magenta.blue_CC <- magenta.blue_CC[order(-magenta.blue_CC$numDEInCat),]
magenta.blue_BP <- subset(magenta.blue_enriched.GO.05, ontology=="BP")
magenta.blue_BP <- magenta.blue_BP[order(-magenta.blue_BP$numDEInCat),]

# Molecular Processes


MFplot_magenta_blue <- magenta.blue_MF %>% mutate(term = fct_reorder(term, (-log10(over_represented_pvalue)) )) %>%
  ggplot( aes(x=term, y=(-log10(over_represented_pvalue)) ) ) +
  geom_segment( aes(x=term ,xend=term, y=0, yend=(-log10(over_represented_pvalue)) ), color="grey") +
  geom_point(size=3, color="#69b3a2") +
  coord_flip() +
  theme(
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_blank(),
    legend.position="none"
  ) +
  xlab("") +
  ylab("-log10(pvalue) [over represented]") +
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
CCplot_magenta_blue <- magenta.blue_CC %>% mutate(term = fct_reorder(term, (-log10(over_represented_pvalue)) )) %>%
  ggplot( aes(x=term, y=(-log10(over_represented_pvalue)) ) ) +
  geom_segment( aes(x=term ,xend=term, y=0, yend=(-log10(over_represented_pvalue)) ), color="grey") +
  geom_point(size=3, color="#69b3a2") +
  coord_flip() +
  theme(
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_blank(),
    legend.position="none"
  ) +
  xlab("") +
  ylab("-log10(pvalue) [over represented]") +
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
BPplot_magenta_blue <- magenta.blue_BP %>% mutate(term = fct_reorder(term,  (-log10(over_represented_pvalue)) )) %>%
  ggplot( aes(x=term, y= (-log10(over_represented_pvalue)) ) ) +
  geom_segment( aes(x=term ,xend=term, y=0, yend= (-log10(over_represented_pvalue)) ), color="grey") +
  geom_point(size=3, color="#69b3a2") +
  coord_flip() +
  theme(
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_blank(),
    legend.position="none") +
  xlab("") +
  ylab("-log10(pvalue) [over represented]") +
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
num_magenta_blue  <- dim(magenta.blue_module_genes)[1] # call num upregulated genes

library(tidyr)
pdf("Analysis/Output/WGCNA/subseq_treatments_all/TimeSeries/goseq_modules/Day21_goseq_MEmagentablue_cluster.pdf", width=12, height=5)
magenta.blue_enriched.GO.05 %>% drop_na(ontology) %>% mutate(term = fct_reorder(term, (-log10(over_represented_pvalue)) )) %>%
  dplyr::filter(ontology %in% c('MF', 'BP')) %>% 
  mutate(term = fct_reorder(term, ontology)) %>%
  ggplot( aes(x=term, y=(-log10(over_represented_pvalue)) ) ) +
  geom_segment( aes(x=term ,xend=term, y=0, yend=(-log10(over_represented_pvalue)), colour = ontology), size = 3, lineend = "butt", alpha = 0.1) + #
  geom_point(size=3, shape = 15, aes(colour = ontology)) +
  coord_flip() +
  theme(
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_blank(),
    legend.position="bottom"
  ) +
  xlab("") +
  ylab("-log10(pvalue)") +
  ggtitle("GO: Day 21 WGCNA MEmagenta.blue_pink Cluster") +
  geom_label(aes(x = 2, y = 4, label = paste(num_magenta_blue, "MEmagenta.blue_pink genes"))) +
  theme_bw() + #Set background color 
  theme(panel.border = element_blank(), # Set border
        panel.grid.major = element_blank(), #Set major gridlines
        panel.grid.minor = element_blank(), #Set minor gridlines
        axis.line = element_line(colour = "black"), #Set axes color
        axis.text.x = element_text(size = 12), # set size of x axis text (log10 vals)
        axis.title.x = element_text(size = 14), # set size of x axis title (log10 vals)
        axis.text.y = element_text(size = 12), # set size of y axis text (terms)
        #legend.position = c(30, 30), 
        plot.background=element_blank()) + #Set the plot background #set title attributes 
  geom_hline(yintercept = 1.3, linetype="dashed", 
             color = "grey", size=1)
dev.off()
