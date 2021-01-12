---
  # title: "Geoduck_WGCNA"
  # author: "Samuel Gurr"
  # date: "January 8, 2021"
---
  
# LOAD PACKAGES
library(WGCNA) # note: this was previously installed with the command `BiocManager::install("WGCNA")`
library(dplyr)
library(zoo)

# SET WORKING DIRECTORY AND LOAD DATA
setwd("C:/Users/samjg/Documents/My_Projects/Pgenerosa_TagSeq_Metabolomics/TagSeq/")
# LOAD DATA
all.counts.matrix <- read.csv(file=" all.counts.filtered.csv", sep=',', header=TRUE)
day0.counts.matrix <- read.csv(file=" day0.counts.filtered.csv", sep=',', header=TRUE)
day7.counts.matrix <- read.csv(file=" day7.counts.filtered.csv", sep=',', header=TRUE)
day14.counts.matrix <- read.csv(file=" day14.counts.filtered.csv", sep=',', header=TRUE)
day21counts.matrix <- read.csv(file=" day21.counts.filtered.csv", sep=',', header=TRUE)

all.treatment.data <- read.csv(file="RAnalysis/WGCNA/ all.treatment.data.csv", sep=',', header=TRUE)
d0.treatment.data <- read.csv(file="RAnalysis/WGCNA/ d0.treatment.data.csv", sep=',', header=TRUE)
d7.treatment.data <- read.csv(file="RAnalysis/WGCNA/ d7.treatment.data.csv", sep=',', header=TRUE)
d14.treatment.data <- read.csv(file="RAnalysis/WGCNA/ d14.treatment.data.csv", sep=',', header=TRUE)
d21.treatment.data <- read.csv(file="RAnalysis/WGCNA/ d21.treatment.data.csv", sep=',', header=TRUE)

# The following setting is important, do not omit.
options(stringsAsFactors = FALSE)

# =================================================================================== #
#
#
# Day 21 Data 
#
#
# =================================================================================== #
  
dim(day21counts.matrix) # 8230  genes  64 samples; should have 62 columns of just samples 
(day21counts.matrix)[1:2] # ommit this and transpose in the next line 
d21.data = as.data.frame(t(day21counts.matrix[, -c(1:2)])) # ommit all columns but samples and transpose
fix(d21.data) # view the data - as you see all columns are now genes but filled as V1, V2, V3...
names(d21.data) = day21counts.matrix$Gene.ID # assigns column names (previous jsut numbered) as the gene ID 
rownames(d21.data) = names(day21counts.matrix)[-c(1:2)]; # assigns the row names as the sample ID


gsg = goodSamplesGenes(d21.data, verbose = 3);  # We first check for genes and samples with too many missing values:
gsg$allOK # If the statement returns TRUE, all genes have passed the cuts. 


# determine outlier samples 
sampleTree = hclust(dist(d21.data), method = "average") # Next we cluster the samples (in contrast to clustering genes that will come later)  to see if there are any obvious outliers.
sizeGrWindow(12,9) # The user should change the dimensions if the window is too large or too small.
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5, 
     cex.axis = 1.5, cex.main = 2) # appears there are two outliers SG55 and SG 105; can remove by hand or an automatic appraoch 
abline(h = 20000, col = "red") # add line to plot to show the cut-off od outlier samples (20000)
clust = cutreeStatic(sampleTree, cutHeight = 20000, minSize = 10) # Determine cluster under the line
table(clust) # 0 = cut; 1 = kept; says it will cut 2 and save 60 
keepSamples = (clust==1) # call the 60 from clust at position 1
d21.data = d21.data[keepSamples, ]
nGenes = ncol(d21.data) # number of genes == 8230 
nSamples = nrow(d21.data) # number of samples == 60

# Treatment data
dim(d21.treatment.data) # 62 rows and 7 columns
names(d21.treatment.data) # look at the 7 columns 
d21.treatment = d21.treatment.data[, -c(1,3, 7)]; # remove columns that hold information we do not need.
# Form a data frame analogous to expression data that will hold the clinical traits.
d21.Samples = rownames(d21.data);# start new variable 'd21.Samples' calling the row names of the gene data (sample as 'Mouse' in trait data)
TreatRows = match(d21.Samples, d21.treatment$Sample.Name); # match the names 
d21.Traits = d21.treatment[TreatRows, -1]; # removes the row numbers
rownames(d21.Traits) = d21.treatment[TreatRows, 1]; # inserts the new traitRows - matches sample treatments
dim(d21.Traits) # 60 Samples 3 columns (primary, second, and third treatment)
all(rownames(d21.Traits) == rownames(d21.data))  # should be TRUE
fix(d21.data)


d21.treatment = d21.treatment.data[, -c(1,4:7)]; # remove columns that hold information we do not need.
# Form a data frame analogous to expression data that will hold the clinical traits.
d21.Samples = rownames(d21.data);# start new variable 'd21.Samples' calling the row names of the gene data (sample as 'Mouse' in trait data)
TreatRows = match(d21.Samples, d21.treatment$Sample.Name); # match the names 
d21.Traits = as.data.frame(d21.treatment[TreatRows, -1]); # removes the row numbers
d21.Traits
rownames(d21.Traits) = d21.treatment[TreatRows, 1]; # inserts the new traitRows - matches sample treatments
colnames(d21.Traits)[1] <- "Treatment"
d21.Traits$Treatment <- as.factor(d21.Traits$Treatment)
dim(d21.Traits) # 60 Samples 3 columns (primary, second, and third treatment)
all(rownames(d21.Traits) == rownames(d21.data))  # should be TRUE



# Re-cluster samples
sampleTree2 = hclust(dist(d21.data), method = "average")
traitColors = labels2colors(d21.Traits); # Convert traits to a color representation: white means low, red means high, grey means missing entry
plotDendroAndColors(sampleTree2, traitColors, # Plot the sample dendrogram and the colors underneath.
                    groupLabels = names(d21.Traits), 
                    main = "Sample dendrogram and trait heatmap")

save(d21.data, d21.Traits, file = "d.21-dataInput.RData")



# soft threshold ============================================================================== #
# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to=20, by=2))
# Call the network topology analysis function
sft = pickSoftThreshold(d21.data, powerVector = powers, verbose = 5) 
# pickSoftThreshold 
#  performs the analysis of network topology and aids the
# user in choosing a proper soft-thresholding power.
# The user chooses a set of candidate powers (the function provides suitable default values)
# function returns a set of network indices that should be inspected
# Plot the results:
sizeGrWindow(9, 5)
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

# The left panel shows the scale-free fit index (y-axis) 
# as a function of the soft-thresholding power (x-axis). 
# We choose the power 9, which is the lowest power for which the scale-free 
# topology index curve attens out upon reaching a high value (in this case, roughly 0.90).

# The right panel displays the mean connectivity
# (degree, y-axis) as a function of the soft-thresholding power (x-axis).

#=====================================================================================
#
#  CConstructing the gene network and identifying modules
#
#=====================================================================================
# Constructing the gene network and identifying modules is now a simple function call:
# soft thresholding power 7,
# relatively large minimum module size of 30,
# medium sensitivity (deepSplit=2) to cluster splitting.
# mergeCutHeight is the threshold for merging of modules.
# instructed the function to return numeric, rather than color, labels for modules,
# and to save the Topological Overlap Matrix
net = blockwiseModules(d21.data, power = 9,
                       TOMType = "unsigned", minModuleSize = 30,
                       reassignThreshold = 0, mergeCutHeight = 0.25,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs = TRUE,
                       saveTOMFileBase = "d21.Geoduck", 
                       verbose = 3)

net$colors # the module assignmen
net$MEs # contains the module eigengenes of the modules.
table(net$colors)
# there are 11  modules in order of decending size
# note 0 is reserved for genes outside of all modules (2049 genes) 

#=====================================================================================
#
#  The hierarchical clustering dendrogram (tree) used for the module identification 
#
#=====================================================================================
# The hierarchical clustering dendrogram (tree) used for the module identification 
# open a graphics window
sizeGrWindow(12, 9)
# Convert labels to colors for plotting
mergedColors = labels2colors(net$colors)
# Plot the dendrogram and the module colors underneath
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)



#=====================================================================================
#
#  Quantifying module{trait associations
#
#=====================================================================================
# identify modules that are signiFcantly
# associated with the measured clinical traits.

# Since we already have a summary profile (eigengene) for each module,
# we simply correlate eigengenes with external traits and look for the most signicant associations:

# Define numbers of genes and samples
nGenes = ncol(d21.data);
nSamples = nrow(d21.data);
# Recalculate MEs with color labels
MEs0 = moduleEigengenes(d21.data, moduleColors)$eigengenes 
MEs = orderMEs(MEs0) # reorders the columns (colors/modules)
d21.Traits$Treatment <- as.numeric(d21.Traits$Treatment)
moduleTraitCor = cor(MEs, d21.Traits, use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);


#=====================================================================================
#
#  Code chunk 3
#
#=====================================================================================


sizeGrWindow(10,6)
# Will display correlations and their p-values
textMatrix =  paste(signif(moduleTraitCor, 2), "\n(",
                    signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 8.5, 3, 3));
# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(d21.Traits),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = T,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(1,12),
               main = paste("Module-trait relationships"))

# The analysis identifies the several significant 
# module{trait associations. We will concentrate on weight as the trait
# of interest.
