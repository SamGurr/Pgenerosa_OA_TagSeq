---
  # title: "Geoduck_WGCNA"
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
day21counts.matrix <- read.csv(file="Analysis/Data/filtered_counts/day21.counts.filtered_5cpm50perc.csv", sep=',', header=TRUE)
# Treatment and Phenotype data
Master.Treatment_Phenotype.data <- read.csv(file="Analysis/Data/ Master_Phyenotype.and.Exp.Treatment_Metadata.csv", sep=',', header=TRUE)
d21.Treatment_Phenotype.data     <- Master.Treatment_Phenotype.data %>%  dplyr ::filter(Date %in% 20190814) # split for day 21 data 


# The following setting is important, do not omit.
options(stringsAsFactors = FALSE)

# =================================================================================== #
#
#
# Day 21 WGCNA - PREPROCESSING THE DATA INPUT 
#
#  Read here for the pre processing steps using WGCNA!
#  https://bmcbioinformatics.biomedcentral.com/track/pdf/10.1186/1471-2105-9-559.pdf
# =================================================================================== #

# filter low values - note: this is pre filteres for < 5 CPM in 50% of samples - less strict than the DESeq2 analysis
dim(day21counts.matrix) # 11800  rows (genes) -   63  samples  not counting 'x' and 'Gene.ID'
(day21counts.matrix)[1] # ommit the first line and transpose in the next line 
d21.data = as.data.frame(t(day21counts.matrix[, -(1)])) # ommit all columns but samples and transpose
dim(d21.data)

 
# trait data ========================================================== #

# Phenotype trait and Treatment data
dim(d21.Treatment_Phenotype.data) # 62 rows and 14 columns
names(d21.Treatment_Phenotype.data) # look at the 14 columns 
d21.Treatment_Phenotype.data = d21.Treatment_Phenotype.data[, -c(1:3,5)]; # remove columns that hold information we do not need.
dim(d21.Treatment_Phenotype.data) # 62 rows and 10 columns


# count data  ========================================================== #


# fix(d21.data) # view the data - as you see all columns are now genes but filled as V1, V2, V3...
names(d21.data) = day21counts.matrix$X # assigns column names (previous jsut numbered) as the gene ID 
rownames(d21.data) = names(day21counts.matrix)[-(1)]; # assigns the row names as the sample ID
d21.data_matrix <- data.frame(day21counts.matrix[,-1], row.names=day21counts.matrix[,1]) 


# create dds objects to transform data 
dds.d21 <- DESeqDataSetFromMatrix(countData = d21.data_matrix,
                                            colData = d21.Treatment_Phenotype.data, design = ~ 1) # DESeq Data Set (dds)
 # DESeq Data Set (dds)
# NOTE: ~1 stands for no design; user will need to add a design for differential testing
# however for our purpose of just creating an onbject to transform, we do not need a design here...


# transform the data 
dds.d21_vst <- vst(dds.d21)
dds.d21_vst <- assay(dds.d21_vst)
dds.d21_vst <- t(dds.d21_vst)

# =================================================================================== #
#
#
# Day 21 WGCNA Sample tree 
#
#
# =================================================================================== #
  
dim(dds.d21_vst) #  11800 genes; 62 samples; should have 62 columns of just samples 

gsg = goodSamplesGenes(dds.d21_vst, verbose = 3);  # We first check for genes and samples with too many missing values:
gsg$allOK # If the statement returns TRUE, all genes have passed the cuts. 

# determine outlier samples 
sampleTree = hclust(dist(dds.d21_vst), method = "average") # Next we cluster the samples (in contrast to clustering genes that will come later)  to see if there are any obvious outliers.
sizeGrWindow(12,9) # The user should change the dimensions if the window is too large or too small.
par(cex = 0.6);
par(mar = c(0,4,2,0))

png("Analysis/Output/WGCNA/Day21/Day21_ClusterTree_Precut.png", 1000, 1000, pointsize=20)
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5, 
     cex.axis = 1.5, cex.main = 2) # appears there are two outliers SG55 and SG 105; can remove by hand or an automatic appraoch 
abline(h = 120, col = "red") # add line to plot to show the cut-off od outlier samples (40000) SG105 and SG55
dev.off()

clust = cutreeStatic(sampleTree, cutHeight = 120, minSize = 10) # Determine cluster under the line
table(clust) # 0 = cut; 1 = kept; says it will cut 2 and save 60 
keepSamples = (clust==1) # call the 60 from clust at position 1
dds.d21_vst = dds.d21_vst[keepSamples, ]
nGenes = ncol(dds.d21_vst) # number of genes == 11800 
nSamples = nrow(dds.d21_vst) # number of samples == 45


sampleTree2 = hclust(dist(dds.d21_vst), method = "average") # Next we cluster the samples (in contrast to clustering genes that will come later)  to see if there are any obvious outliers.
png("Analysis/Output/WGCNA/Day21/Day21_ClusterTree_Postcut.png", 1000, 1000, pointsize=20)
plot(sampleTree2, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5, 
     cex.axis = 1.5, cex.main = 2)
dev.off()


# based on outlier removall... call Trait data ===================================================== #
dim(dds.d21_vst) # transformed has 61 samples - one cut in the tree step above 
dim(d21.Treatment_Phenotype.data) # trait data has 62 samples - not yet cut! 

# Form a data frame analogous to expression data that will hold the clinical traits.
d21.Samples = rownames(dds.d21_vst);# start new variable 'd21.Samples' calling the row names of the gene data (sample as 'Mouse' in trait data)
TreatRows = match(d21.Samples, d21.Treatment_Phenotype.data$Sample.Name); # match the names 
d21.Traits = d21.Treatment_Phenotype.data[TreatRows, -1]; # removes the row numbers
rownames(d21.Traits) = d21.Treatment_Phenotype.data[TreatRows, 1]; # inserts the new traitRows - matches sample treatments
dim(d21.Traits) #  61 Samples 9 columns (primary, second, and third treatment)
all(rownames(d21.Traits) == rownames(dds.d21_vst))  # should be TRUE
dim(d21.Traits) #  45  9


# ALL TRAITS 
d21.Traits
# ONLY TREATMENTS
d21.Traits.treat<- d21.Traits[,c(1:3)]
# ONLY PHYSIOLOGY
d21.Traits.phys<- d21.Traits[,c(4:(ncol(d21.Traits)))]




# Re-cluster samples

# All Traits
png("Analysis/Output/WGCNA/Day21/Day21_ClusterTree_AllTraits.png", 1000, 1000, pointsize=20)
sampleTree2 = hclust(dist(dds.d21_vst), method = "average")
traitColors = labels2colors(d21.Traits); # Convert traits to a color representation: white means low, red means high, grey means missing entry
plotDendroAndColors(sampleTree2, traitColors, # Plot the sample dendrogram and the colors underneath.
                    groupLabels = names(d21.Traits), 
                    main = "Sample dendrogram and trait heatmap")
dev.off()

# Treatment Only
png("Analysis/Output/WGCNA/Day21/Day21_ClusterTree_TreatmentOnly.png", 1000, 1000, pointsize=20)
traitColors_treat = labels2colors(d21.Traits.treat); # Convert traits to a color representation: white means low, red means high, grey means missing entry
plotDendroAndColors(sampleTree2, traitColors_treat, # Plot the sample dendrogram and the colors underneath.
                    groupLabels = names(d21.Traits.treat), 
                    main = "Sample dendrogram and trait heatmap")
dev.off()

# Phys Only
png("Analysis/Output/WGCNA/Day21/Day21_ClusterTree_PhysOnly.png", 1000, 1000, pointsize=20)
traitColors_phys = labels2colors(d21.Traits.phys); # Convert traits to a color representation: white means low, red means high, grey means missing entry
plotDendroAndColors(sampleTree2, traitColors_phys, # Plot the sample dendrogram and the colors underneath.
                    groupLabels = names(d21.Traits.phys), 
                    main = "Sample dendrogram and trait heatmap")
dev.off()


# save data
save(dds.d21_vst, d21.Traits, file = "Analysis/Output/WGCNA/Day21/d.21-dataInput.RData")
save(dds.d21_vst, d21.Traits.phys, file = "Analysis/Output/WGCNA/Day21/d.21-dataInput_treatmentonly.RData")
save(dds.d21_vst, d21.Traits.treat, file = "Analysis/Output/WGCNA/Day21/d.21-dataInput_physonly.RData")


# write the vst transformed data 
write.csv(dds.d21_vst, "Analysis/Output/WGCNA/Day21/Day21_vstTransformed_WGCNAdata.csv") # write


# soft threshold ============================================================================== #



# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to=20, by=2))
# Call the network topology analysis function
sft = pickSoftThreshold(dds.d21_vst, powerVector = powers, verbose = 5) #...wait for this to finish
# pickSoftThreshold 
#  performs the analysis of network topology and aids the
# user in choosing a proper soft-thresholding power.
# The user chooses a set of candidate powers (the function provides suitable default values)
# function returns a set of network indices that should be inspected




# pickSoftThreshold ------------------------------------------------------------------------- #



sizeGrWindow(9, 5)
png("Analysis/Output/WGCNA/Day21/Day21_ScaleTopology_SoftThresh.png", 1000, 1000, pointsize=20)
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

net = blockwiseModules(dds.d21_vst, power = 4,
                       TOMType = "unsigned", minModuleSize = 30,
                       reassignThreshold = 0, mergeCutHeight = 0.25,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs = TRUE,
                       saveTOMFileBase = "Analysis/Output/WGCNA/d21.Geoduck", 
                       verbose = 3) #... wait for this to finish 

net$colors # the module assignmen
net$MEs # contains the module eigengenes of the modules.
table(net$colors)
# there are 19  modules in order of decending size
# note 0 is reserved for genes outside of all modules (3675 genes) 

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
plotDendroAndColors(net$dendrograms[[1]], moduleColors[net$blockGenes[[1]]],
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
nGenes = ncol(dds.d21_vst); # 11800
nSamples = nrow(dds.d21_vst); # 45
# Recalculate MEs with color labels
MEs0 = moduleEigengenes(dds.d21_vst, moduleColors)$eigengenes 
MEs = orderMEs(MEs0) # reorders the columns (colors/modules)


# write csv - save the module eigengenes
write.csv(MEs, file = "Analysis/Output/WGCNA/Day21/d21.WGCNA_ModulEigengenes.csv")


# change chanracter treatments to integers
# ALL TRAIT DATA
d21.Traits$Primary_Treatment <- as.factor(d21.Traits$Primary_Treatment)
d21.Traits$Primary_Treatment <- as.numeric(d21.Traits$Primary_Treatment)

d21.Traits$Second_Treament <- as.factor(d21.Traits$Second_Treament)
d21.Traits$Second_Treament <- as.numeric(d21.Traits$Second_Treament)

d21.Traits$Third_Treatment <- as.factor(d21.Traits$Third_Treatment)
d21.Traits$Third_Treatment <- as.numeric(d21.Traits$Third_Treatment)

# TREATMENT ONLY 
d21.Traits.treat$Primary_Treatment <- as.factor(d21.Traits.treat$Primary_Treatment)
d21.Traits.treat$Primary_Treatment <- as.numeric(d21.Traits.treat$Primary_Treatment)

d21.Traits.treat$Second_Treament <- as.factor(d21.Traits.treat$Second_Treament)
d21.Traits.treat$Second_Treament <- as.numeric(d21.Traits.treat$Second_Treament)

d21.Traits.treat$Third_Treatment <- as.factor(d21.Traits.treat$Third_Treatment)
d21.Traits.treat$Third_Treatment <- as.numeric(d21.Traits.treat$Third_Treatment)





# Module trait correlation ------------------------------------------------------- #
# ALL TRAIT DATA
dim(d21.Traits)  # 45  9
dim(MEs)  # 45 27
moduleTraitCor = cor(MEs, d21.Traits, use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);
# TREATMENT ONLY 
moduleTraitCor_treatonly = cor(MEs, d21.Traits.treat, use = "p");
moduleTraitPvalue_treatonly = corPvalueStudent(moduleTraitCor_treatonly, nSamples);
# PHYSIOLOGY  ONLY 
moduleTraitCor_physonly = cor(MEs, d21.Traits.phys, use = "p");
moduleTraitPvalue_physonly = corPvalueStudent(moduleTraitCor_physonly, nSamples);






#=====================================================================================
#
#  labeled heatmap
#
#=====================================================================================

# ALL TRAITS ------------------------------------------------------------------------- # 
sizeGrWindow(10,10)
# Will display correlations and their p-values
d21.AllTraits.matrix <-  paste(signif(moduleTraitCor, 2), "\n(",
                    signif(moduleTraitPvalue, 1), ")", sep = "");

#dim(d21.AllTraits.text) == dim(moduleTraitCor)
par(mar = c(8, 9.5, 5, 3));
# Display the correlation values within a heatmap plot
png("Analysis/Output/WGCNA/Day21/Day21_AllTraits_heatmap.png", 500, 1000, pointsize=20)
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(d21.Traits),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = TRUE,
               colors = blueWhiteRed(50),
               textMatrix = d21.AllTraits.matrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-0.5,0.5),
               main = paste("Module-trait relationships - All traits"))
dev.off()


# this heatmap looks better
d21.AllTraits.text <-  as.matrix(signif(moduleTraitCor, 3))
pa = cluster::pam(moduleTraitCor, k = 3)
col_fun = colorRamp2(c(-0.5, 0, 0.5), c("blue", "white", "red"))
png("Analysis/Output/WGCNA/Day21/Day21_AllTraits_heatmap2.png", 500, 1000, pointsize=20)
Heatmap(moduleTraitCor, 
        name = "gene_cor", 
        rect_gp = gpar(col = "grey", lwd = 1),
        column_title = "Day 21 WGCNA - All traits", 
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
         grid.text(sprintf("%.1f", d21.AllTraits.text[i, j]), x, y, gp = gpar(fontsize = 10))
       })
dev.off()

# The analysis identifies the several significant 
# module{trait associations. We will concentrate on weight as the trait
# of interest.





# TREATMENT ONLY ------------------------------------------------------------------------- # 
sizeGrWindow(10,10)
# Will display correlations and their p-values
d21.Treatments.matrix <-  paste(signif(moduleTraitCor_treatonly, 2), "\n(",
                    signif(moduleTraitCor_treatonly, 1), ")", sep = "")


#dim(textMatrix) == dim(moduleTraitCor_treatonly)
par(mar = c(8, 9.5, 5, 3));
# Display the correlation values within a heatmap plot
png("Analysis/Output/WGCNA/Day21/Day21_Treatments_heatmap.png", 500, 1000, pointsize=20)
labeledHeatmap(Matrix = moduleTraitCor_treatonly,
               xLabels = names(d21.Traits.treat),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = TRUE,
               colors = blueWhiteRed(50),
               textMatrix = d21.Treatments.matrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships - Treatments only"))
dev.off()



# this heatmap looks better
d21.Treatments.text <-  as.matrix(signif(moduleTraitCor_treatonly, 3))
pa = cluster::pam(d21.Treatments.text, k = 3)
col_fun = colorRamp2(c(-0.5, 0, 0.5), c("blue", "white", "red"))
png("Analysis/Output/WGCNA/Day21/Day21_Treatments_heatmap2.png", 500, 1000, pointsize=20)
Heatmap(moduleTraitCor_treatonly, 
        name = "gene_cor", 
        rect_gp = gpar(col = "grey", lwd = 1),
        column_title = "Day 21 WGCNA - Treatments only", 
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
          grid.text(sprintf("%.1f", d21.Treatments.text[i, j]), x, y, gp = gpar(fontsize = 10))
        })
dev.off()

# The analysis identifies the several significant 
# module{trait associations. We will concentrate on weight as the trait
# of interest.





# PHYSIOLOGY ONLY  ------------------------------------------------------------------------- # 
sizeGrWindow(10,10)
# Will display correlations and their p-values
d21.Phys.matrix = paste(signif(moduleTraitCor_physonly, 2), "\n(",
                    signif(moduleTraitPvalue_physonly, 1), ")", sep = "")

par(mar = c(8, 9.5, 5, 3));
# Display the correlation values within a heatmap plot
png("Analysis/Output/WGCNA/Day21/Day21_PhysOnly_heatmap.png", 500, 1000, pointsize=20)
labeledHeatmap(Matrix = moduleTraitCor_physonly,
               xLabels = names(d21.Traits.phys),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = TRUE,
               colors = blueWhiteRed(50),
               textMatrix = d21.Phys.matrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-0.5,0.5),
               main = paste("Module-trait relationships - Phys only"))
dev.off()



# this heatmap looks better
d21.Phys.text <-  as.matrix(signif(moduleTraitCor_physonly, 3))
pa = cluster::pam(d21.Phys.text, k = 3)
col_fun = colorRamp2(c(-0.5, 0, 0.5), c("blue", "white", "red"))
png("Analysis/Output/WGCNA/Day21/Day21_PhysOnly_heatmap2.png", 500, 1000, pointsize=20)
Heatmap(moduleTraitCor_physonly, 
        name = "gene_cor", 
        rect_gp = gpar(col = "grey", lwd = 1),
        column_title = "Day 21 WGCNA - Physiology only", 
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
          grid.text(sprintf("%.1f", d21.Phys.text[i, j]), x, y, gp = gpar(fontsize = 10))
        })
dev.off()
# The analysis identifies the several significant 
# module{trait associations. We will concentrate on weight as the trait
# of interest.




#=====================================================================================
#
#  Total Protein in red module! CREATE DATAFRAMES 
#
#=====================================================================================
# Gene relationship to trait and important modules: 
# Gene Significance and Module Membership

# We quantify associations of individual genes with our trait of interest (TAOC)

# Define variable weight containing the aTAOC column 
TAOC = as.data.frame(d21.Traits$mean.µmol.CRE.g.protein); # -0.53 in yellow module
names(TAOC) = "TAOC"

# names (colors) of the modules
modNames = substring(names(MEs), 3) # name all the modules, from 3rd character on (first two are ME)

geneModuleMembership = as.data.frame(cor(dds.d21_vst, MEs, use = "p"));
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));

names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");

# TAOC
TAOC_geneTraitSignificance = as.data.frame(cor(dds.d21_vst, TAOC, use = "p"));
TAOC_GSPvalue = as.data.frame(corPvalueStudent(as.matrix(TAOC_geneTraitSignificance), nSamples));

names(TAOC_geneTraitSignificance) = paste("GS.TAOC", names(TAOC), sep="");
names(TAOC_GSPvalue) = paste("p.GS.TAOC", names(TAOC), sep="");



#  PLOT mean.µmol.CRE.g.protein in the MAGENTA module
unique(moduleColors)
module = "magenta"
column = match(module, modNames);
moduleGenes = moduleColors==module;

sizeGrWindow(7, 7);
par(mfrow = c(1,1));
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(TAOC_geneTraitSignificance[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Gene significance for TAOC",
                   main = paste("Day21 total protein 'MAGENTA': Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = "black")

#  PLOT mean.µmol.CRE.g.protein in the MAGENTA module
unique(moduleColors)
module = "lightcyan"
column = match(module, modNames);
moduleGenes = moduleColors==module;

sizeGrWindow(7, 7);
par(mfrow = c(1,1));
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(TAOC_geneTraitSignificance[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Gene significance for TAOC",
                   main = paste("Day21 total protein 'LIGHTCYAN': Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = "black")

#=====================================================================================
#
#   COUNT GENES OF INTEREST IN  MODULES (i.e. lighcyan and magenta - refer to heatmap)
#
#=====================================================================================


length(colnames(dds.d21_vst)[moduleColors=="lightcyan"]) # 162 total genes in the lightcyan module
length(colnames(dds.d21_vst)[moduleColors=="magenta"]) # 350 total genes in the magenta module

#=====================================================================================
#
#  Call annotation data to get module gene data (prep for downstream GO)
#
#=====================================================================================


annot = read.delim2(file = "Seq_details/Panopea-generosa-genes-annotations.txt",header = F);
dim(annot) # 34947     9
names(annot) # V1 are the gene.IDs
probes = names(d21.data)
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
names(TAOC_geneTraitSignificance)
names(TAOC_GSPvalue)
geneInfo_TAOC = data.frame(substanceBXH = probes,
                           geneSymbol = annot$V1[probes2annot],
                           #LocusLinkID = annot$LocusLinkID[probes2annot],
                           moduleColor = moduleColors,
                           Uniprot = annot$V5[probes2annot],
                           HGNC = annot$V6[probes2annot],
                           GO.terms = annot$V8[probes2annot],
                           GO.description = annot$V9[probes2annot],
                           TAOC_geneTraitSignificance, # call this specific to the module and trait of interest
                           TAOC_GSPvalue)              # call this specific to the module and trait of interest
View(geneInfo_TAOC)
modOrder = order(-abs(cor(MEs, TAOC, use = "p"))); # order by the strength of the correlation between module and trait values for each sample

for (mod in 1:ncol(geneModuleMembership)) # Add module membership information in the chosen order
{
  oldNames = names(geneInfo_TAOC)
  geneInfo_TAOC = data.frame(geneInfo_TAOC, geneModuleMembership[, modOrder[mod]], 
                             MMPvalue[, modOrder[mod]]);
  names(geneInfo_TAOC) = c(oldNames, paste("MM.", modNames[modOrder[mod]], sep=""),
                           paste("p.MM.", modNames[modOrder[mod]], sep=""))
}
# Order the genes in the geneInfo variable first by module color, then by geneTraitSignificance
geneOrder = order(geneInfo_TAOC$moduleColor, -abs(geneInfo_TAOC$GS.TAOCTAOC));
geneInfo_TAOC = geneInfo_TAOC[geneOrder, ]
View(geneInfo_TAOC)


#=====================================================================================
#
#  Write csv for the modules and corresponding raw read counts
#
#=====================================================================================
# call the module of interest for follow-up GO analysis 


write.csv(geneInfo_TAOC, file = "Analysis/Output/WGCNA/Day21/d21.WGCNA_ModulMembership.csv")



#=====================================================================================
#
#  Plotting expression in the red module
#
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
d21_ModEigen <- read.csv("Analysis/Output/WGCNA/Day21/d21.WGCNA_ModulEigengenes.csv")
d21_Annot_ModuleMembership <-  read.csv("Analysis/Output/WGCNA/Day21/d21.WGCNA_ModulMembership.csv")
d21_vst_data <- read.csv("Analysis/Output/WGCNA/Day21/Day21_vstTransformed_WGCNAdata.csv")
Master.Treatment_Phenotype.data <- read.csv(file="Analysis/Data/ Master_Phyenotype.and.Exp.Treatment_Metadata.csv", sep=',', header=TRUE)

# Prep data to merge
# module membership
d21_Annot_ModuleMembership <- d21_Annot_ModuleMembership[-c(1:2)] # ommit the redundant gene name columns buut keep 'geneSymbol'
dim(d21_Annot_ModuleMembership) #  11800  62

# module eigengenes
names(d21_ModEigen)[1] <- "Sample.Name"

# normalized count data (same used for the dds object and WGCNA analysis)
d21_vst_data_t <- as.data.frame(t(d21_vst_data[, -(1)])) # trsnpose columns names to rows (genes) 
colnames(d21_vst_data_t) <- d21_vst_data[, 1] # name the columns as the rows in previous dataset (sample IDS)
d21_vst_data_t<- d21_vst_data_t %>% tibble::rownames_to_column("geneSymbol") # "geneSymbol" - create column and remname the gene name column to 'geneSymbol'


# merge Master datasets
# GO terms
d21_GOTermsMaster.Modules<-  merge(d21_Annot_ModuleMembership, d21_vst_data_t, by = "geneSymbol")
dim(d21_GOTermsMaster) # 11800   107
# Eigengenes and traits
d21_EigenTraitMaster<-  merge(d21_ModEigen, Master.Treatment_Phenotype.data, by = "Sample.Name")
dim(d21_EigenTraitMaster) # 45 41





# plots 




# exampl plots below
ggplot(d21_EigenTraitMaster, 
       aes(x=MElightcyan, y=mean.µmol.CRE.g.protein, 
           color = Primary_Treatment)) + # linetype= Primary_Treatment, group = Primary_Treatment
  geom_smooth(method=lm , color="black",  fill="grey", se=TRUE, linetype = "dashed") + #, linetype = "dashed"
  geom_point(size = 4, shape = 19) +
  scale_color_manual(values=c("#56B4E9", "#D55E00")) +
  labs(title="WGCNA Module = Lightcyan",x="Eigengene expression by sample", y = "mean TAOC")+
  theme_bw()

ggplot(d21_EigenTraitMaster, 
       aes(x=MEmagenta, y=mean.µmol.CRE.g.protein, 
           color = Primary_Treatment)) + # linetype= Primary_Treatment, group = Primary_Treatment
  geom_smooth(method=lm , color="black",  fill="grey", se=TRUE, linetype = "dashed") + #, linetype = "dashed"
  geom_point(size = 4, shape = 19) +
  scale_color_manual(values=c("#56B4E9", "#D55E00")) +
  labs(title="WGCNA Module = Magenta",x="Eigengene expression by sample", y = "mean TAOC")+
  theme_bw()




# lets loop this through to get a bunch of plots 
phys <- d21_EigenTraitMaster[,c(36:41)]
modules <- d21_EigenTraitMaster[,c(2:28)]

install.packages('ggpmisc')
library(ggpmisc)

for(i in 1:ncol(modules)) {
  module_color <- substr(names(modules)[i], 3,12)
  moduledata <- modules[,i]
  
  par(mfrow=c(3,2))
  
  my.formula <- y ~ x
  
  p <- ggplot(d21_EigenTraitMaster, 
                   aes(x=moduledata, y=phys[,1], 
                       color = Primary_Treatment, fill = Primary_Treatment, group=Primary_Treatment)) + # linetype= Primary_Treatment, group = Primary_Treatment
    scale_fill_manual(values=c("#56B4E9", "#D55E00")) +
    geom_smooth(method=lm , color="grey50", alpha = 0.1, se=TRUE, linetype = "dashed", formula = my.formula) + #, linetype = "dashed"
    stat_poly_eq(formula = my.formula,
                 eq.with.lhs = "italic(hat(y))~`=`~",
                 aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
                 parse = TRUE)  + 
    geom_point(size = 4, shape = 19) +
    scale_color_manual(values=c("#56B4E9", "#D55E00")) +
    labs(x="Eigengene expression by sample", y = paste(names(phys)[1])) +
    theme_bw()
  p1 <- p + ggtitle(paste("WGCNAmodule:", module_color, "&",(names(phys)[1]), sep =' ')) + theme(legend.position="none")
  
  
  p <- ggplot(d21_EigenTraitMaster, 
              aes(x=moduledata, y=phys[,2], 
                  color = Primary_Treatment, fill = Primary_Treatment, group=Primary_Treatment)) + # linetype= Primary_Treatment, group = Primary_Treatment
    scale_fill_manual(values=c("#56B4E9", "#D55E00")) +
    geom_smooth(method=lm , color="grey50", alpha = 0.1, se=TRUE, linetype = "dashed", formula = my.formula) + #, linetype = "dashed"
    stat_poly_eq(formula = my.formula,
                 eq.with.lhs = "italic(hat(y))~`=`~",
                 aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
                 parse = TRUE)  + 
    geom_point(size = 4, shape = 19) +
    scale_color_manual(values=c("#56B4E9", "#D55E00")) +
    labs(x="Eigengene expression by sample", y = paste(names(phys)[2])) +
    theme_bw()
  p2 <- p + ggtitle(paste("WGCNAmodule:", module_color, "&",(names(phys)[2]), sep =' ')) + theme(legend.position="none")
  
  
  p <- ggplot(d21_EigenTraitMaster, 
              aes(x=moduledata, y=phys[,3], 
                  color = Primary_Treatment, fill = Primary_Treatment, group=Primary_Treatment)) + # linetype= Primary_Treatment, group = Primary_Treatment
    scale_fill_manual(values=c("#56B4E9", "#D55E00")) +
    geom_smooth(method=lm , color="grey50", alpha = 0.1, se=TRUE, linetype = "dashed", formula = my.formula) + #, linetype = "dashed"
    stat_poly_eq(formula = my.formula,
                 eq.with.lhs = "italic(hat(y))~`=`~",
                 aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
                 parse = TRUE)  + 
    geom_point(size = 4, shape = 19) +
    scale_color_manual(values=c("#56B4E9", "#D55E00")) +
    labs(x="Eigengene expression by sample", y = paste(names(phys)[3])) +
    theme_bw()
  p3 <- p + ggtitle(paste("WGCNAmodule:", module_color, "&",(names(phys)[3]), sep =' ')) + theme(legend.position="none")
  

  p <- ggplot(d21_EigenTraitMaster, 
              aes(x=moduledata, y=phys[,4], 
                  color = Primary_Treatment, fill = Primary_Treatment, group=Primary_Treatment)) + # linetype= Primary_Treatment, group = Primary_Treatment
    scale_fill_manual(values=c("#56B4E9", "#D55E00")) +
    geom_smooth(method=lm , color="grey50", alpha = 0.1, se=TRUE, linetype = "dashed", formula = my.formula) + #, linetype = "dashed"
    stat_poly_eq(formula = my.formula,
                 eq.with.lhs = "italic(hat(y))~`=`~",
                 aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
                 parse = TRUE)  + 
    geom_point(size = 4, shape = 19) +
    scale_color_manual(values=c("#56B4E9", "#D55E00")) +
    labs(x="Eigengene expression by sample", y = paste(names(phys)[4])) +
    theme_bw()
  p4 <- p + ggtitle(paste("WGCNAmodule:", module_color, "&",(names(phys)[4]), sep =' ')) + theme(legend.position="none")
  
  
  p <- ggplot(d21_EigenTraitMaster, 
              aes(x=moduledata, y=phys[,5], 
                  color = Primary_Treatment, fill = Primary_Treatment, group=Primary_Treatment)) + # linetype= Primary_Treatment, group = Primary_Treatment
    scale_fill_manual(values=c("#56B4E9", "#D55E00")) +
    geom_smooth(method=lm , color="grey50", alpha = 0.1, se=TRUE, linetype = "dashed", formula = my.formula) + #, linetype = "dashed"
    stat_poly_eq(formula = my.formula,
                 eq.with.lhs = "italic(hat(y))~`=`~",
                 aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
                 parse = TRUE)  + 
    geom_point(size = 4, shape = 19) +
    scale_color_manual(values=c("#56B4E9", "#D55E00")) +
    labs(x="Eigengene expression by sample", y = paste(names(phys)[5])) +
    theme_bw()
  p5 <- p + ggtitle(paste("WGCNAmodule:", module_color, "&",(names(phys)[5]), sep =' ')) + theme(legend.position="none")
  
  p <- ggplot(d21_EigenTraitMaster, 
              aes(x=moduledata, y=phys[,6], 
                  color = Primary_Treatment, fill = Primary_Treatment, group=Primary_Treatment)) + # linetype= Primary_Treatment, group = Primary_Treatment
    scale_fill_manual(values=c("#56B4E9", "#D55E00")) +
    geom_smooth(method=lm , color="grey50", alpha = 0.1, se=TRUE, linetype = "dashed", formula = my.formula) + #, linetype = "dashed"
    stat_poly_eq(formula = my.formula,
                 eq.with.lhs = "italic(hat(y))~`=`~",
                 aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
                 parse = TRUE)  + 
    geom_point(size = 4, shape = 19) +
    scale_color_manual(values=c("#56B4E9", "#D55E00")) +
    labs(x="Eigengene expression by sample", y = paste(names(phys)[6])) +
    theme_bw()
  p6 <- p + ggtitle(paste("WGCNAmodule:", module_color, "&",(names(phys)[6]), sep =' ')) + theme(legend.position="none")
  
  # save plot grid 
  #png("Analysis/Output/WGCNA/Day21/Day21_PhysOnly_heatmap2.png", 1000, 1000, pointsize=20)
  ggarrange(p1, p2, p3, p4, p5, p6, ncol = 2, nrow = 3)
  ggsave(paste0("Analysis/Output/WGCNA/Day21/EigengenePlots/Day21_EigengenePlot_", module_color, ".png"))
}
  





# EXPLORE THE LIGHTCYAN MODULE =================================================================================================== #
# evidence from the heatmaps and the regression shwo strong module membership and eigenenge cor strength with lighcyan 
# in response to treatment (primary specifically) and Total Antioxidant capacity 





# Plot the vst transformed by treatment and by Total protein to visualize the effect of the red module genes significant in WGCNA 
# vst read count date - narrow the columns - reshape and rename
lightcyan_geneIDs <- d21_Annot_ModuleMembership %>% dplyr::filter(moduleColor %in% 'lightcyan') %>%  dplyr::select("geneSymbol")
d21_vst_LighCyanModule <- d21_vst_data_t %>% dplyr::filter(geneSymbol %in% lightcyan_geneIDs[,1])
d21_vst_LighCyanModule_MELT <- melt(d21_vst_LighCyanModule, id=("geneSymbol")) # melt using reshape2
names(d21_vst_LighCyanModule_MELT)[(2:3)] <- c('Sample.Name', 'vst_Expression') # change column names

# experiment treatment and total protein data - narrow the columns 
exp.phys_data <- Master.Treatment_Phenotype.data %>% 
  dplyr::filter(Date %in% 20190814) %>%  # filter out Day 21 data only 
  dplyr::select(c('Sample.Name', 'All_Treatment', 'Primary_Treatment', 'Second_Treament', 'Third_Treatment', 'mean.µmol.CRE.g.protein')) # call select columns

# merge by common row values 'Sample.Name'
merged_Expdata_LightCyan <- merge(d21_vst_LighCyanModule_MELT, exp.phys_data, by ='Sample.Name')

# mean Exp response table 
meanEXp_LightCyan <- merged_Expdata_LightCyan %>% 
                        select(c('Sample.Name','vst_Expression','Primary_Treatment', 'Second_Treament','Third_Treatment')) %>% 
                        group_by(Sample.Name, Primary_Treatment, Second_Treament, Third_Treatment) %>%
                        dplyr::summarize(mean.vstExp = mean(vst_Expression), 
                                         sd.vsdtExp = sd(vst_Expression),
                                         na.rm=TRUE)
# summarize datasets further by treatment period
# remember:this is a mean of a mean!! First we complete mean vst exp by sample id (compiling all red module genes) - next all sample IDs by the treatment period (below
# I will use these for mean SE plots 
# Primary treatment
meanEXp_Summary.Prim_LightCyan <- meanEXp_LightCyan %>% 
                                      group_by(Primary_Treatment) %>%
                                      dplyr::summarize(mean = mean(mean.vstExp), 
                                                       sd = sd(mean.vstExp),
                                                       n = n(), 
                                                       se = sd/sqrt(n))
# Second treatment
meanEXp_Summary.Sec_LightCyan <- meanEXp_LightCyan %>% 
                                      group_by(Second_Treament,Primary_Treatment) %>%
                                      dplyr::summarize(mean = mean(mean.vstExp), 
                                                       sd = sd(mean.vstExp),
                                                       n = n(), 
                                                       se = sd/sqrt(n))

# Third treatment
meanEXp_Summary.Third_LightCyan  <- meanEXp_LightCyan %>% 
                                      group_by(Third_Treatment,Second_Treament,Primary_Treatment) %>%
                                      dplyr::summarize(mean = mean(mean.vstExp), 
                                                       sd = sd(mean.vstExp),
                                                       n = n(), 
                                                       se = sd/sqrt(n))


# Primary treatment plot
# PrimaryBoxplot.vst <- ggplot(meanEXp) +
#                           aes(x = Primary_Treatment, y = mean.vstExp, fill=Primary_Treatment) +
#                           theme_classic() +
#                           geom_boxplot(varwidth = FALSE, width = 0.2) + # vary boxes width according to n obs.
#                           geom_jitter(alpha = 0.25, width = 0.1) # adds random noise and limit its width





# Mean SE plots 




# The errorbars overlapped, so use position_dodge to move them horizontally
pd <- position_dodge(0.3) # move them .05 to the left and right


# Primary treatment mean sd plot
Primary.vst.LightCyan <- ggplot(meanEXp_Summary.Prim_LightCyan, aes(x=Primary_Treatment, y=mean, fill=Primary_Treatment)) +  # , colour=supp, group=supp))
                              theme_bw() +
                              geom_errorbar(aes(ymin=mean-se, ymax=mean+se), colour="black", width=.1, position=pd) +
                              geom_line(position=pd) +
                              geom_point(position=pd, size = 4, shape=21) +            
                              xlab("Primary pCO2 treatment") +
                              ylab(NULL) +                 # note the mean was first by sample ID THEN by treatment
                              scale_fill_manual(values=c("#56B4E9","#E69F00")) +
                              # scale_color_manual(values=c("#56B4E9","#E69F00")) +
                              ggtitle("Day 21 WGCNA LightCyan Module VST GeneExp") +
                              # expand_limits(y=0) +                                                    # Expand y range
                              scale_y_continuous(limits=c(3.5, 4.5)) +
                              theme(text = element_text(size=15))


# Second treatment mean sd plot
Sec.vst.LightCyan <- ggplot(meanEXp_Summary.Sec_LightCyan, aes(x=Second_Treament, y=mean, fill=Primary_Treatment, group=Primary_Treatment)) +
                              theme_bw() +
                              geom_errorbar(aes(ymin=mean-se, ymax=mean+se), colour="black", width=.1, position=pd) +
                              geom_line(position=pd) +
                              geom_point(position=pd, size = 4, shape=21) +            
                              xlab("Second pCO2 treatment") +
                              ylab("LightCyan Module VST Gene Expression (Mean +/- SE)") +                 # note the mean was first by sample ID THEN by treatment
                              scale_fill_manual(values=c("#56B4E9","#E69F00")) +
                              # scale_color_manual(values=c("#56B4E9","#E69F00")) +
                              # ggtitle("Day 21 WGCNA red' Module VST GeneExp") +
                              # expand_limits(y=0) +                                                    # Expand y range
                              scale_y_continuous(limits=c(3.5, 4.5)) +
                              theme(text = element_text(size=15))


# Third treatment mean sd plot
meanEXp_Summary.Third_LightCyan$Prim.Sec_Treatment <- paste(meanEXp_Summary.Third_LightCyan$Primary_Treatment, meanEXp_Summary.Third_LightCyan$Second_Treament, sep='')
Third.vst.LightCyan <- ggplot(meanEXp_Summary.Third_LightCyan, aes(x=Third_Treatment, y=mean, fill=Primary_Treatment, group=Prim.Sec_Treatment, alpha = Second_Treament)) +
                              theme_bw() +
                              geom_errorbar(aes(ymin=mean-se, ymax=mean+se), colour="black", width=.1, position=pd,alpha = 1) +
                              geom_line(position=pd, colour="black",alpha = 1) +
                              geom_point(position=pd, size = 4) +            
                              xlab("Third pCO2 treatment") +
                              ylab(NULL) +                 # note the mean was first by sample ID THEN by treatment
                              scale_fill_manual(values=c("#56B4E9","#E69F00")) +
                              # scale_color_manual(values=c("#56B4E9","#E69F00")) +
                              # ggtitle("Day 21 WGCNA red' Module VST GeneExp") +
                              scale_y_continuous(limits=c(3.5, 4.5)) +
                              theme(text = element_text(size=15))

# Total protein by vst expression
meanExp_LightCyan_TAOC <- merged_Expdata_LightCyan %>% 
  group_by(Sample.Name, Primary_Treatment, Second_Treament, Third_Treatment) %>%
  dplyr::summarize(mean.Exp = mean(vst_Expression), 
                   sd.Exp = sd(vst_Expression),
                   mean.TAOC = mean(mean.µmol.CRE.g.protein), 
                   sd.TAOC = sd(mean.µmol.CRE.g.protein))
meanExp_LightCyan_TAOC # notice cd total protein is 0 - this is becasue total protein was the same value regarless of the gene.ID within sample


# two obvious outliers > 50
ggplot(meanExp_total.protein, aes(x=mean.totalprotein, y=mean.Exp)) +
  geom_point(shape=1) +    # Use hollow circles
  geom_smooth()            # Add a loess smoothed fit curve with confidence region

meanExp_total.protein.OM <- meanExp_total.protein %>% filter(mean.totalprotein < 50) # remove the outliers and plot again

# Total protein and expression plot
TAOC_vst.Exp.LightCyan <- ggplot(meanExp_LightCyan_TAOC, aes(x=mean.TAOC, y=mean.Exp,
  color = Primary_Treatment, fill = Primary_Treatment, group=Primary_Treatment)) + # linetype= Primary_Treatment, group = Primary_Treatment
  scale_fill_manual(values=c("#56B4E9", "#D55E00")) +
  geom_smooth(method=lm , color="grey50", alpha = 0.1, se=TRUE, linetype = "dashed", formula = my.formula) + #, linetype = "dashed"
  stat_poly_eq(formula = my.formula,
               eq.with.lhs = "italic(hat(y))~`=`~",
               aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
               parse = TRUE)  + 
  geom_point(size = 4, shape = 19) +
  scale_color_manual(values=c("#56B4E9", "#D55E00")) +
  labs(x="Eigengene expression by sample", y = "mean.TAOC") +
  theme_bw() +
  scale_y_continuous(limits=c(3.5, 4.5)) +
  theme(text = element_text(size=15))

# 
# TAOC_vst.Exp.LightCyan <- ggplot(meanExp_LightCyan_TAOC, aes(x=mean.TAOC, y=mean.Exp,
#  color = Third_Treatment, fill = Third_Treatment, group=Third_Treatment)) + # linetype= Primary_Treatment, group = Primary_Treatment
#   scale_fill_manual(values=c("#56B4E9", "#D55E00")) +
#   geom_smooth(method=lm , color="grey50", alpha = 0.1, se=TRUE, linetype = "dashed", formula = my.formula) + #, linetype = "dashed"
#   stat_poly_eq(formula = my.formula,
#                eq.with.lhs = "italic(hat(y))~`=`~",
#                aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
#                parse = TRUE)  + 
#   geom_point(size = 4, shape = 19) +
#   scale_color_manual(values=c("#56B4E9", "#D55E00")) +
#   labs(x="Eigengene expression by sample", y = "mean.TAOC") +
#   theme_bw() +
#   scale_y_continuous(limits=c(3.5, 4.5)) +
#   theme(text = element_text(size=15))


#  grid arrange the three plots 
# save plot
png("Analysis/Output/WGCNA/Day21/Day21_vst.Exp_LightcyanModule.png", 1000, 1000, pointsize=20)
ggarrange(Primary.vst.LightCyan, Sec.vst.LightCyan, Third.vst.LightCyan, TAOC_vst.Exp.LightCyan, 
          plotlist = NULL,
          ncol = 2,
          nrow = 2,
          labels = NULL)
dev.off()



  
