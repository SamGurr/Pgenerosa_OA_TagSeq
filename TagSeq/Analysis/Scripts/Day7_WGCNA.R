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



# SET WORKING DIRECTORY AND LOAD DATA
setwd("C:/Users/samjg/Documents/My_Projects/Pgenerosa_TagSeq_Metabolomics/TagSeq/")
# LOAD DATA
# Tagaseq filtered counts 
day7.counts.matrix <- read.csv(file="Analysis/Data/filtered_counts/day7.counts.filtered_5cpm50perc.csv", sep=',', header=TRUE)
# Treatment and Phenotype data
Master.Treatment_Phenotype.data <- read.csv(file="Analysis/Data/ Master_Phyenotype.and.Exp.Treatment_Metadata.csv", sep=',', header=TRUE)
d7.Treatment_Phenotype.data     <- Master.Treatment_Phenotype.data %>%  dplyr ::filter(Date %in% 20190731) # split for day 7 data 


# The following setting is important, do not omit.
options(stringsAsFactors = FALSE)

# =================================================================================== #
#
#
# Day 7 WGCNA - PREPROCESSING THE DATA INPUT 
#
#  Read here for the pre processing steps using WGCNA!
#  https://bmcbioinformatics.biomedcentral.com/track/pdf/10.1186/1471-2105-9-559.pdf
# =================================================================================== #


# filter low values - note: this is pre filteres for < 5 CPM in 50% of samples - less strict than the DESeq2 analysis
dim(day7.counts.matrix) # 12013  rows (genes) -   37  samples  not counting 'x' and 'Gene.ID'
(day7.counts.matrix)[1] # ommit the first line and transpose in the next line 
d7.data = as.data.frame(t(day7.counts.matrix[, -(1)])) # ommit all columns but samples and transpose
dim(d7.data) # 36 12013



# trait data ========================================================== #

# Phenotype trait and Treatment data
dim(d7.Treatment_Phenotype.data) #  36 rows and 14 columns
names(d7.Treatment_Phenotype.data) # look at the 14 columns 
d7.Treatment_Phenotype.data = d7.Treatment_Phenotype.data[, -c(1:3,5,8)]; # remove columns that hold information we do not need. (i.e. tankks, Third treatment, All treatment group, etc.)
dim(d7.Treatment_Phenotype.data) # 36 rows and  9 columns


# count data  ========================================================== #


# fix(d21.data) # view the data - as you see all columns are now genes but filled as V1, V2, V3...
names(d7.data) = day7.counts.matrix$X # assigns column names (previous jsut numbered) as the gene ID 
rownames(d7.data) = names(day7.counts.matrix)[-(1)]; # assigns the row names as the sample ID
d7.data_matrix <- data.frame(day7.counts.matrix[,-1], row.names=day7.counts.matrix[,1]) 


# create dds objects to transform data 
dds.d7 <- DESeqDataSetFromMatrix(countData = d7.data_matrix,
                                  colData = d7.Treatment_Phenotype.data, design = ~ 1) # DESeq Data Set (dds)
# DESeq Data Set (dds)
# NOTE: ~1 stands for no design; user will need to add a design for differential testing
# however for our purpose of just creating an onbject to transform, we do not need a design here...


# transform the data 
dds.d7_vst <- vst(dds.d7) # transform it vst
dds.d7_vst <- assay(dds.d7_vst) # call only the transformed coutns in the dds object
dds.d7_vst <- t(dds.d7_vst) # transpose columns to rows and vice versa

# =================================================================================== #
#
#
# Day 7 WGCNA Sample tree 
#
#
# =================================================================================== #

dim(dds.d7_vst) #  12013 genes; 36  samples

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
abline(h = 100, col = "red") # add line to plot to show the cut-off od outlier samples (40000) SG105 and SG55
dev.off()

clust = cutreeStatic(sampleTree, cutHeight = 100, minSize = 10) # Determine cluster under the line
table(clust) # 0 = cut; 1 = kept; says it will cut 2 and save 60 
keepSamples = (clust==1) # call the 30 from clust at position 1
dds.d7_vst = dds.d7_vst[keepSamples, ] # integreat the boolean 'keepsamples' to ommit oultilers determined in the sample tree above 
nGenes = ncol(dds.d7_vst) # number of genes == 12013 
nSamples = nrow(dds.d7_vst) # number of samples == 30  - the cut tree removed 6 samples 


sampleTree2 = hclust(dist(dds.d7_vst), method = "average") # Next we cluster the samples (in contrast to clustering genes that will come later)  to see if there are any obvious outliers.
png("Analysis/Output/WGCNA/Day7/Day7_ClusterTree_Postcut.png", 1000, 1000, pointsize=20)
plot(sampleTree2, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5, 
     cex.axis = 1.5, cex.main = 2)
dev.off()


# based on outlier removall... call Trait data ===================================================== #
dim(dds.d7_vst) # transformed count data now  has 30 samples - 6 cut in the tree step above 
dim(d7.Treatment_Phenotype.data) # trait data has  36 samples - not yet cut! 

# Form a data frame analogous to expression data that will hold the clinical traits.
d7.Samples = rownames(dds.d7_vst);# start new variable 'd21.Samples' calling the row names of the gene data (sample as 'Mouse' in trait data)
TreatRows = match(d7.Samples, d7.Treatment_Phenotype.data$Sample.Name); # match the names 
d7.Traits = d7.Treatment_Phenotype.data[TreatRows, -1]; # removes the row numbers
rownames(d7.Traits) = d7.Treatment_Phenotype.data[TreatRows, 1]; # inserts the new traitRows - matches sample treatments
dim(d7.Traits) #  30 Samples 8 columns (primary and second treatment) - note: no Third treatment because that occured Days 14-21, after this sampling day
all(rownames(d7.Traits) == rownames(dds.d7_vst))  # should be TRUE
dim(d7.Traits) #  30  8


# ALL TRAITS 
d7.Traits
# ONLY TREATMENTS
d7.Traits.treat<- d7.Traits[,c(1:2)]
# ONLY PHYSIOLOGY
d7.Traits.phys<- d7.Traits[,c(3:(ncol(d7.Traits)))]




# Re-cluster samples

# All Traits
png("Analysis/Output/WGCNA/Day7/Day7_ClusterTree_AllTraits.png", 1000, 1000, pointsize=20)
sampleTree2 = hclust(dist(dds.d7_vst), method = "average")
traitColors = labels2colors(d7.Traits); # Convert traits to a color representation: white means low, red means high, grey means missing entry
plotDendroAndColors(sampleTree2, traitColors, # Plot the sample dendrogram and the colors underneath.
                    groupLabels = names(d7.Traits), 
                    main = "Sample dendrogram and trait heatmap")
dev.off()

# Treatment Only
png("Analysis/Output/WGCNA/Day7/Day7_ClusterTree_TreatmentOnly.png", 1000, 1000, pointsize=20)
traitColors_treat = labels2colors(d7.Traits.treat); # Convert traits to a color representation: white means low, red means high, grey means missing entry
plotDendroAndColors(sampleTree2, traitColors_treat, # Plot the sample dendrogram and the colors underneath.
                    groupLabels = names(d7.Traits.treat), 
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


# soft threshold ============================================================================== #



# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to=20, by=2))
# Call the network topology analysis function
sft = pickSoftThreshold(dds.d7_vst, powerVector = powers, verbose = 5) #...wait for this to finish
# pickSoftThreshold 
#  performs the analysis of network topology and aids the
# user in choosing a proper soft-thresholding power.
# The user chooses a set of candidate powers (the function provides suitable default values)
# function returns a set of network indices that should be inspected




# pickSoftThreshold ------------------------------------------------------------------------- #



sizeGrWindow(9, 5)
png("Analysis/Output/WGCNA/Day7/Day7_ScaleTopology_SoftThresh.png", 1000, 1000, pointsize=20)
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
dev.off()


# The left panel shows the scale-free fit index (y-axis) 
# as a function of the soft-thresholding power (x-axis). 
# We choose the power 9, which is the lowest power for which the scale-free 
# topology index curve attens out upon reaching a high value (in this case, roughly 0.90).

# The right panel displays the mean connectivity
# (degree, y-axis) as a function of the soft-thresholding power (x-axis).

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

net = blockwiseModules(dds.d7_vst, power = 4,
                       TOMType = "unsigned", minModuleSize = 30,
                       reassignThreshold = 0, mergeCutHeight = 0.25,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs = TRUE,
                       saveTOMFileBase = "Analysis/Output/WGCNA/Day7/d7.Geoduck", 
                       verbose = 3) #... wait for this to finish 

net$colors # the module assignmen
net$MEs # contains the module eigengenes of the modules.
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
# Plot the dendrogram and the module colors underneath
png("Analysis/Output/WGCNA/Day7/Day7_ClusterDendrogram.png", 1000, 1000, pointsize=20)
plotDendroAndColors(net$dendrograms[[1]], moduleColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()

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
nGenes = ncol(dds.d7_vst); # 12013
nSamples = nrow(dds.d7_vst); #  30
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


# TREATMENT ONLY 
d7.Traits.treat$Primary_Treatment <- as.factor(d7.Traits.treat$Primary_Treatment)
d7.Traits.treat$Primary_Treatment <- as.numeric(d7.Traits.treat$Primary_Treatment)

d7.Traits.treat$Second_Treament <- as.factor(d7.Traits.treat$Second_Treament)
d7.Traits.treat$Second_Treament <- as.numeric(d7.Traits.treat$Second_Treament)


# Module trait correlation ------------------------------------------------------- #
# ALL TRAIT DATA
dim(d7.Traits)  # 30  8
dim(MEs)  # 30 38
moduleTraitCor = cor(MEs, d7.Traits, use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);
# TREATMENT ONLY 
moduleTraitCor_treatonly = cor(MEs, d7.Traits.treat, use = "p");
moduleTraitPvalue_treatonly = corPvalueStudent(moduleTraitCor_treatonly, nSamples);
# PHYSIOLOGY  ONLY 
moduleTraitCor_physonly = cor(MEs, d7.Traits.phys, use = "p");
moduleTraitPvalue_physonly = corPvalueStudent(moduleTraitCor_physonly, nSamples);


#=====================================================================================
#
#  labeled heatmap
#
#=====================================================================================

# ALL TRAITS ------------------------------------------------------------------------- # 
sizeGrWindow(10,10)
# Will display correlations and their p-values
d7.AllTraits.matrix <-  paste(signif(moduleTraitCor, 2), "\n(",
                               signif(moduleTraitPvalue, 1), ")", sep = "");

#dim(d21.AllTraits.text) == dim(moduleTraitCor)
par(mar = c(8, 9.5, 5, 3));
# Display the correlation values within a heatmap plot
png("Analysis/Output/WGCNA/Day7/Day7_AllTraits_heatmap.png", 500, 1000, pointsize=20)
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(d7.Traits),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = TRUE,
               colors = blueWhiteRed(50),
               textMatrix = d7.AllTraits.matrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-0.6,0.6),
               main = paste("Module-trait relationships - All traits"))
dev.off()


# this heatmap looks better
d7.AllTraits.text <-  as.matrix(signif(moduleTraitCor, 3))
pa = cluster::pam(moduleTraitCor, k = 3)
col_fun = colorRamp2(c(-0.5, 0, 0.5), c("blue", "white", "red"))
png("Analysis/Output/WGCNA/Day7/Day7_AllTraits_heatmap2.png", 500, 1000, pointsize=20)
Heatmap(moduleTraitCor, 
        name = "gene_cor", 
        rect_gp = gpar(col = "grey", lwd = 1),
        column_title = "Day 7 WGCNA - All traits", 
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
          grid.text(sprintf("%.1f", d7.AllTraits.text[i, j]), x, y, gp = gpar(fontsize = 10))
        })
dev.off()

# The analysis identifies the several significant 
# module{trait associations. We will concentrate on weight as the trait
# of interest.



# TREATMENT ONLY ------------------------------------------------------------------------- # 
sizeGrWindow(10,10)
# Will display correlations and their p-values
d7.Treatments.matrix <-  paste(signif(moduleTraitCor_treatonly, 2), "\n(",
                                signif(moduleTraitCor_treatonly, 1), ")", sep = "")


#dim(textMatrix) == dim(moduleTraitCor_treatonly)
par(mar = c(8, 9.5, 5, 3));
# Display the correlation values within a heatmap plot
png("Analysis/Output/WGCNA/Day7/Day7_Treatments_heatmap.png", 500, 1000, pointsize=20)
labeledHeatmap(Matrix = moduleTraitCor_treatonly,
               xLabels = names(d7.Traits.treat),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = TRUE,
               colors = blueWhiteRed(50),
               textMatrix = d7.Treatments.matrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-0.6,0.6),
               main = paste("Module-trait relationships - Treatments only"))
dev.off()



# this heatmap looks better
d7.Treatments.text <-  as.matrix(signif(moduleTraitCor_treatonly, 3))
pa = cluster::pam(d7.Treatments.text, k = 3)
col_fun = colorRamp2(c(-0.5, 0, 0.5), c("blue", "white", "red"))
png("Analysis/Output/WGCNA/Day7/Day7_Treatments_heatmap2.png", 500, 1000, pointsize=20)
Heatmap(moduleTraitCor_treatonly, 
        name = "gene_cor", 
        rect_gp = gpar(col = "grey", lwd = 1),
        column_title = "Day 7 WGCNA - Treatments only", 
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
          grid.text(sprintf("%.1f", d7.Treatments.text[i, j]), x, y, gp = gpar(fontsize = 10))
        })
dev.off()

# The analysis identifies the several significant 
# module{trait associations. We will concentrate on weight as the trait
# of interest.



# PHYSIOLOGY ONLY  ------------------------------------------------------------------------- # 
sizeGrWindow(10,10)
# Will display correlations and their p-values
d7.Phys.matrix = paste(signif(moduleTraitCor_physonly, 2), "\n(",
                        signif(moduleTraitPvalue_physonly, 1), ")", sep = "")

par(mar = c(8, 9.5, 5, 3));
# Display the correlation values within a heatmap plot
png("Analysis/Output/WGCNA/Day7/Day7_PhysOnly_heatmap.png", 500, 1000, pointsize=20)
labeledHeatmap(Matrix = moduleTraitCor_physonly,
               xLabels = names(d7.Traits.phys),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = TRUE,
               colors = blueWhiteRed(50),
               textMatrix = d7.Phys.matrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-0.5,0.5),
               main = paste("Module-trait relationships - Phys only"))
dev.off()



# this heatmap looks better
d7.Phys.text <-  as.matrix(signif(moduleTraitCor_physonly, 3))
pa = cluster::pam(d7.Phys.text, k = 3)
col_fun = colorRamp2(c(-0.5, 0, 0.5), c("blue", "white", "red"))
png("Analysis/Output/WGCNA/Day7/Day7_PhysOnly_heatmap2.png", 500, 1000, pointsize=20)
Heatmap(moduleTraitCor_physonly, 
        name = "gene_cor", 
        rect_gp = gpar(col = "grey", lwd = 1),
        column_title = "Day 7 WGCNA - Physiology only", 
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
          grid.text(sprintf("%.1f", d7.Phys.text[i, j]), x, y, gp = gpar(fontsize = 10))
        })
dev.off()
# The analysis identifies the several significant 
# module{trait associations. We will concentrate on weight as the trait
# of interest.

#=====================================================================================
#
#  Shell Length in violet module! CREATE DATAFRAMES 
#
#=====================================================================================
# Gene relationship to trait and important modules: 
# Gene Significance and Module Membership

# We quantify associations of individual genes with our trait of interest (TAOC)

# Define variable weight containing the aTAOC column 
ShellLength = as.data.frame(d7.Traits$mean.Shell.Length); # -0.53 in yellow module
names(ShellLength) = "ShellLength"

# names (colors) of the modules
modNames = substring(names(MEs), 3) # name all the modules, from 3rd character on (first two are ME)

geneModuleMembership = as.data.frame(cor(dds.d7_vst, MEs, use = "p"));
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));

names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");

# TAOC
ShellLength_geneTraitSignificance = as.data.frame(cor(dds.d7_vst, ShellLength, use = "p"));
ShellLength_GSPvalue = as.data.frame(corPvalueStudent(as.matrix(ShellLength_geneTraitSignificance), nSamples));

names(ShellLength_geneTraitSignificance) = paste("GS.ShellLength", names(ShellLength), sep="");
names(ShellLength_GSPvalue) = paste("p.GS.ShellLength", names(ShellLength), sep="");



#  PLOT mean.µmol.CRE.g.protein in the MAGENTA module
unique(moduleColors)
module = "violet"
column = match(module, modNames);
moduleGenes = moduleColors==module;

sizeGrWindow(7, 7);
par(mfrow = c(1,1));
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(ShellLength_geneTraitSignificance[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Gene significance for Shell length",
                   main = paste("Day7 shell length 'VIOLET': Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = "black")

#=====================================================================================
#
#   COUNT GENES OF INTEREST IN  MODULES (i.e. violet- refer to heatmap)
#
#=====================================================================================

length(colnames(dds.d7_vst)[moduleColors=="violet"]) # 48 total genes in the violet module

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
names(ShellLength_geneTraitSignificance)
names(ShellLength_GSPvalue)
geneInfo_ShellLength = data.frame(substanceBXH = probes,
                                  geneSymbol = annot$V1[probes2annot],
                                  #LocusLinkID = annot$LocusLinkID[probes2annot],
                                  moduleColor = moduleColors,
                                  Uniprot = annot$V5[probes2annot],
                                  HGNC = annot$V6[probes2annot],
                                  GO.terms = annot$V8[probes2annot],
                                  GO.description = annot$V9[probes2annot],
                                  ShellLength_geneTraitSignificance, # call this specific to the module and trait of interest
                                  ShellLength_GSPvalue)              # call this specific to the module and trait of interest
View(geneInfo_ShellLength)
modOrder = order(-abs(cor(MEs, ShellLength, use = "p"))); # order by the strength of the correlation between module and trait values for each sample

for (mod in 1:ncol(geneModuleMembership)) # Add module membership information in the chosen order
{
  oldNames = names(geneInfo_ShellLength)
  geneInfo_ShellLength = data.frame(geneInfo_ShellLength, geneModuleMembership[, modOrder[mod]], 
                                    MMPvalue[, modOrder[mod]]);
  names(geneInfo_ShellLength) = c(oldNames, paste("MM.", modNames[modOrder[mod]], sep=""),
                                  paste("p.MM.", modNames[modOrder[mod]], sep=""))
}
# Order the genes in the geneInfo variable first by module color, then by geneTraitSignificance
geneOrder = order(geneInfo_ShellLength$moduleColor, -abs(geneInfo_ShellLength$GS.ShellLengthShellLength));
geneInfo_ShellLength = geneInfo_ShellLength[geneOrder, ]
View(geneInfo_ShellLength)



#=====================================================================================
#
#  Write csv for the modules and corresponding raw read counts
#
#=====================================================================================
# call the module of interest for follow-up GO analysis 


write.csv(geneInfo_ShellLength, file = "Analysis/Output/WGCNA/Day7/d7.WGCNA_ModulMembership.csv")



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
dim(d7_Annot_ModuleMembership) #  12013    84

# module eigengenes
names(d7_ModEigen)[1] <- "Sample.Name"

# normalized count data (same used for the dds object and WGCNA analysis)
d7_vst_data_t <- as.data.frame(t(d7_vst_data[, -(1)])) # trsnpose columns names to rows (genes) 
colnames(d7_vst_data_t) <- d7_vst_data[, 1] # name the columns as the rows in previous dataset (sample IDS)
d7_vst_data_t<- d7_vst_data_t %>% tibble::rownames_to_column("geneSymbol") # "geneSymbol" - create column and remname the gene name column to 'geneSymbol'


# merge Master datasets
# GO terms
d7_GOTermsMaster.Modules<-  merge(d7_Annot_ModuleMembership, d7_vst_data_t, by = "geneSymbol")
dim(d7_GOTermsMaster.Modules) # 12013   114
# Eigengenes and traits
d7_EigenTraitMaster<-  merge(d7_ModEigen, Master.Treatment_Phenotype.data, by = "Sample.Name")
dim(d7_EigenTraitMaster) # 30 52




#=====================================================================================
# 
# Linear regression plots Eigenengene expression with Traits for EACH module       
#
#=====================================================================================

# lets loop this through to get a bunch of plots 
phys <- d7_EigenTraitMaster[,c(47:52)]
modules <- d7_EigenTraitMaster[,c(2:39)]


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
    theme_bw() +
    geom_errorbar(aes(ymin=mean-se, ymax=mean+se), colour="black", width=.1, position=pd) +
    geom_line(position=pd) +
    geom_point(position=pd, size = 4, shape=21) +            
    xlab("Primary pCO2 treatment") +
    ylab(NULL) +                 # note the mean was first by sample ID THEN by treatment
    scale_fill_manual(values=c("#56B4E9","#E69F00")) +
    # scale_color_manual(values=c("#56B4E9","#E69F00")) +
    ggtitle(paste("Day 21 WGCNA", modcolor[i,], "Module VST GeneExp", sep =' ')) +
    # expand_limits(y=0) +                                                    # Expand y range
    scale_y_continuous(limits=c((min_p1), (max_p1))) +
    theme(text = element_text(size=15))
  
  
  # Second treatment mean sd plot
  min_p2 <- min(meanEXp_Summary.Sec_Mod$mean) - max(meanEXp_Summary.Sec_Mod$se)
  max_p2 <- max(meanEXp_Summary.Sec_Mod$mean) + max(meanEXp_Summary.Sec_Mod$se)
  
  Sec.vst.Mod <- ggplot(meanEXp_Summary.Sec_Mod, aes(x=Second_Treament, y=mean, fill=Primary_Treatment, group=Primary_Treatment)) +
    theme_bw() +
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
  
  png(paste("Analysis/Output/WGCNA/Day7/ModuleExpression_Treatment/Day7_Exp_Module",modcolor[i,],".png"), 600, 1000, pointsize=20)
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
# magenta MODULE CONTINUED....
# 
# goseq - load the annotation and prepare the four steps for goseq!
#
#===================================================================================== 

### Panopea generosa - load .fna ('Geoduck_annotation') and foramt GO terms ('Geoduck_GOterms') and vectors
Geoduck_annotation <- read.delim2(file="C:/Users/samjg/Documents/My_Projects/Pgenerosa_TagSeq_Metabolomics/TagSeq/Seq_details/Panopea-generosa-genes-annotations.txt", header=F)


# build annotation file to merge with the mean LFC tables
annot.condenced <- Geoduck_annotation[,c(1,3:9)]
annot.condenced$gene.length <- annot.condenced$V4 - annot.condenced$V3
annot.condenced <- annot.condenced[,-c(2,3)]
names(annot.condenced) <- c('Gene.ID', 'Uniprot', 'HGNC', 'fxn', 'Go.terms', 'Go.fxns','gene.length')




# merge with the master gene annotation file 'annot.condenced' to get all desired info of just the magenta genes
GO_magenta_MASTER <- merge(GO_magenta,annot.condenced, by ='Gene.ID')


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
# Example: 0 = not in magenta; 1 = in magenta
# magenta module  ------------------------------------------- #
GO_magenta <- d7_Annot_ModuleMembership %>% dplyr::filter(moduleColor %in% 'magenta') # %>%  dplyr::select("geneSymbol")
GO_magenta_genes <- GO_magenta[1]
names(GO_magenta_genes)[1] <- "Gene.ID" # 162 genws in the magenta module 


GO_magenta_integer <- as.integer(GO_unique.genes.all%in%(GO_magenta_genes$Gene.ID)) # Day 0 - Primary - Upregulated DEGs
names(GO_magenta_integer)=GO_unique.genes.all
View(GO_magenta_integer)


# Review what we have for goseq....

# (1) Go terms ---------------- #
#  head(GO.terms)

# (2) ID vector --------------- #
# head(GO_unique.genes.all)

# (3) Length vector ----------- #
# head(length_vector)

# (4) Binary DEG vectors ------ # 
# head(GO_magenta_integer)  

#===================================================================================== 
# magenta MODULE CONTINUED....
# 
# It's goseq time!!!
# Calculate Probability Weighting Function
#===================================================================================== 

#Calculate Probability Weighting Function (using 'nullp')
magenta.pwf<-nullp(GO_magenta_integer, GO_unique.genes.all, bias.data=length_vector) #weight vector by length of gene

#===================================================================================== 
# magenta MODULE CONTINUED....
# 
# It's goseq time!!!
#  Run goseq
#===================================================================================== 

magenta.goseq<-goseq(magenta.pwf, ID.vector, gene2cat=GO.terms, test.cats=c("GO:CC", "GO:BP", "GO:MF"), method="Wallenius", use_genes_without_cat=TRUE)

# somegenes <- magenta.goseq %>%  dplyr::filter(ontology %in% 'CC') # xplore why there were no CC terms in the goseq plot
# somegenes0.5CC <- somegenes %>% dplyr::filter(over_represented_pvalue < 0.1)
#===================================================================================== 
# magenta MODULE CONTINUED....
# 
# It's goseq time!!!
# PLOTTING
#===================================================================================== 
library(forcats)
#Find only enriched GO terms that are statistically significant at cutoff 
magenta_enriched.GO.05.a<-magenta.goseq$category[magenta.goseq$over_represented_pvalue<.05] # change twice here
magenta_enriched.GO.05<-data.frame(magenta_enriched.GO.05.a)
colnames(magenta_enriched.GO.05) <- c("category")
magenta_enriched.GO.05 <- merge(magenta_enriched.GO.05, magenta.goseq, by="category") # change here
magenta_enriched.GO.05 <- magenta_enriched.GO.05[order(-magenta_enriched.GO.05$numDEInCat),]
magenta_enriched.GO.05$term <- as.factor(magenta_enriched.GO.05$term)
head(magenta_enriched.GO.05)

magenta_MF <- subset(magenta_enriched.GO.05, ontology=="MF")
magenta_MF <- magenta_MF[order(-magenta_MF$numDEInCat),]
magenta_CC <- subset(magenta_enriched.GO.05, ontology=="CC")
magenta_CC <- magenta_CC[order(-magenta_CC$numDEInCat),]
magenta_BP <- subset(magenta_enriched.GO.05, ontology=="BP")
magenta_BP <- magenta_BP[order(-magenta_BP$numDEInCat),]

# Molecular Processes
MFplot_magenta <- magenta_MF %>% mutate(term = fct_reorder(term, numDEInCat)) %>%
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
CCplot_magenta <- magenta_CC %>% mutate(term = fct_reorder(term, numDEInCat)) %>%
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
BPplot_magenta <- magenta_BP %>% mutate(term = fct_reorder(term, numDEInCat)) %>%
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
num_magenta   <- dim(GO_magenta_genes)[1] # call num upregulated genes

library(tidyr)
GOplot.merge_magenta<- magenta_enriched.GO.05 %>% drop_na(ontology) %>% mutate(term = fct_reorder(term, numDEInCat)) %>%
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
  ggtitle("GO: Day 7 WGCNA MEmagenta") +
  geom_label(aes(x = 2, y = 7, label = paste(num_magenta, "MEmagenta genes"))) +
  theme_bw() + #Set background color 
  theme(panel.border = element_blank(), # Set border
        panel.grid.major = element_blank(), #Set major gridlines
        panel.grid.minor = element_blank(), #Set minor gridlines
        axis.line = element_line(colour = "black"), #Set axes color
        plot.background=element_blank()) #Set the plot background #set title attributes
GOplot.merge_magenta
