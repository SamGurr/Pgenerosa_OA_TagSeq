---
  # title: "Day7_WGCNA"
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
# Tagaseq filtered counts 
day7.counts.matrix <- read.csv(file="Analysis/Data/filtered_counts/day7.counts.filtered_10cpm50perc.csv", sep=',', header=TRUE)
# Treatment and Phenotype data
Master.Treatment_Phenotype.data <- read.csv(file="Analysis/Data/ Master_Phyenotype.and.Exp.Treatment_Metadata.csv", sep=',', header=TRUE)
d7.Treatment_Phenotype.data     <- Master.Treatment_Phenotype.data %>%  dplyr ::filter(Date %in% 20190731) # split for day 7 data 

# check the data 

# Day 7 
dim(day7.counts.matrix) # 37 columns - 36 samples not counting 'x' and 'Gene.ID'
names(day7.counts.matrix)[1] <- 'Gene.ID'
dim(d7.Treatment_Phenotype.data) # 36 total rows unique for each sample ID

# The following setting is important, do not omit.
#options(stringsAsFactors = FALSE)

# =================================================================================== #
#
#
# Day 7 WGCNA 
#
#
# =================================================================================== #

dim(day7.counts.matrix) #  8548  genes  37 samples; should have 36 columns of just samples 
(day7.counts.matrix)[1] # ommit this and transpose in the next line 
d7.data = as.data.frame(t(day7.counts.matrix[, -1])) # ommit all columns but samples and transpose
# fix(d7.data) # view the data - as you see all columns are now genes but filled as V1, V2, V3...
names(d7.data) = day7.counts.matrix$Gene.ID # assigns column names (previous jsut numbered) as the gene ID 
rownames(d7.data) = names(day7.counts.matrix)[-1]; # assigns the row names as the sample ID
dim(d7.data) # 36 8548
#fix(d7.data)
gsg = goodSamplesGenes(d7.data, verbose = 3);  # We first check for genes and samples with too many missing values:
gsg$allOK # If the statement returns TRUE, all genes have passed the cuts. 


# determine outlier samples 
sampleTree = hclust(dist(d7.data), method = "average") # Next we cluster the samples (in contrast to clustering genes that will come later)  to see if there are any obvious outliers.
sizeGrWindow(12,9) # The user should change the dimensions if the window is too large or too small.
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5, 
     cex.axis = 1.5, cex.main = 2) # appears there are two outliers SG55 and SG 105; can remove by hand or an automatic appraoch 
abline(h = 30000, col = "red") # add line to plot to show the cut-off od outlier samples (30000) SG90
clust = cutreeStatic(sampleTree, cutHeight = 30000, minSize = 10) # Determine cluster under the line
table(clust) # 0 = cut; 1 = kept; says it will cut 1 and save 35 
keepSamples = (clust==1) # call the 35 from clust at position 1
d7.data = d7.data[keepSamples, ]
nGenes = ncol(d7.data) # number of genes ==  8548 
nSamples = nrow(d7.data) # number of samples ==  35

# Phenotype trait and Treatment data
dim(d7.Treatment_Phenotype.data) # 36 rows and 14 columns
names(d7.Treatment_Phenotype.data) # look at the 14 columns 
d7.Treatment_Phenotype.data = d7.Treatment_Phenotype.data[, -c(1:3, 5, 8)]; # remove columns that hold information we do not need.
dim(d7.Treatment_Phenotype.data) # 36 rows and 9 columns
d7.Treatment_Phenotype.data$Group <- paste(d7.Treatment_Phenotype.data$Primary_Treatment, d7.Treatment_Phenotype.data$Second_Treament,  sep ="")

# Form a data frame analogous to expression data that will hold the clinical traits.
d7.Samples = rownames(d7.data);# start new variable 'd7.Samples' calling the row names of the gene data (sample as 'Mouse' in trait data)
TreatRows = match(d7.Samples, d7.Treatment_Phenotype.data$Sample.Name); # match the names 
d7.Traits = d7.Treatment_Phenotype.data[TreatRows, -1]; # removes the row numbers
rownames(d7.Traits) = d7.Treatment_Phenotype.data[TreatRows, 1]; # inserts the new traitRows - matches sample treatments
dim(d7.Traits) # 35 Samples 9 columns (primary, second, and third treatment)
all(rownames(d7.Traits) == rownames(d7.data))  # should be TRUE
fix(d7.Traits) # view the data


# Form a data frame analogous to expression data that will hold the clinical traits.
# d7.Samples = rownames(d7.data);# start new variable 'd7.Samples' calling the row names of the gene data (sample as 'Mouse' in trait data)
# TreatRows = match(d7.Samples, d7.treatment$Sample.Name); # match the names 
# d7.Traits = as.data.frame(d7.treatment[TreatRows, -1]); # removes the row numbers
# d7.Traits
# rownames(d7.Traits) = d7.treatment[TreatRows, 1]; # inserts the new traitRows - matches sample treatments
# colnames(d7.Traits)[1] <- "Treatment"
# d7.Traits$Treatment <- as.factor(d7.Traits$Treatment)
# dim(d7.Traits) # 60 Samples 3 columns (primary, second, and third treatment)
# all(rownames(d7.Traits) == rownames(d7.data))  # should be TRUE



# Re-cluster samples
sampleTree2 = hclust(dist(d7.data), method = "average")
traitColors = labels2colors(d7.Traits); # Convert traits to a color representation: white means low, red means high, grey means missing entry
plotDendroAndColors(sampleTree2, traitColors, # Plot the sample dendrogram and the colors underneath.
                    groupLabels = names(d7.Traits), 
                    main = "Sample dendrogram and trait heatmap")

save(d7.data, d7.Traits, file = "Analysis/Output/WGCNA/d.7-dataInput.RData")



# soft threshold ============================================================================== #
# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to=20, by=2))
# Call the network topology analysis function
sft = pickSoftThreshold(d7.data, powerVector = powers, verbose = 5) # ...wait for this to finish


# pickSoftThreshold ------------------------------------------------------------------------- #
# performs the analysis of network topology and aids the
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
#  Constructing the gene network and identifying modules
#
#=====================================================================================
# Constructing the gene network and identifying modules is now a simple function call:
# soft thresholding power 12,
# relatively large minimum module size of 30,
# medium sensitivity (deepSplit=2) to cluster splitting.
# mergeCutHeight is the threshold for merging of modules.
# instructed the function to return numeric, rather than color, labels for modules,
# and to save the Topological Overlap Matrix

# NOTE: you will get this message "Error: REAL() can only be applied to a 'numeric', not a 'integer'"
# if your data is not all as numeric - set as numeric with lapply below BEFORE running 'blockwiseModules'
total_cols <- length(colnames(d7.data))
d7.data[2:total_cols] = lapply(d7.data[2:total_cols], as.numeric)

net = blockwiseModules(d7.data, power = 12,
                       TOMType = "unsigned", minModuleSize = 30,
                       reassignThreshold = 0, mergeCutHeight = 0.25,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs = TRUE,
                       saveTOMFileBase = "Analysis/Output/WGCNA/d7.Geoduck", 
                       verbose = 3) #... wait for this to finish 

net$colors # the module assignmen
net$MEs # contains the module eigengenes of the modules.
table(net$colors)
# there are 13  modules in order of decending size
# note 0 is reserved for genes outside of all modules (2308) 

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
#  Quantifying Eigengenes and module trait associations
#
#=====================================================================================
# identify modules that are signiFcantly
# associated with the measured clinical traits.

# Since we already have a summary profile (eigengene) for each module,
# we simply correlate eigengenes with external traits and look for the most signicant associations:

# Define numbers of genes and samples
nGenes = ncol(d7.data); # 8548
nSamples = nrow(d7.data); # 35
# Recalculate MEs with color labels
MEs0 = moduleEigengenes(d7.data, moduleColors)$eigengenes 
MEs = orderMEs(MEs0) # reorders the columns (colors/modules)

# change chanracter treatments to integers
d7.Traits$Primary_Treatment <- as.factor(d7.Traits$Primary_Treatment) # convert to factor before numeric
d7.Traits$Primary_Treatment <- as.numeric(d7.Traits$Primary_Treatment)# now convet to numeric for WGCNA

d7.Traits$Second_Treament <- as.factor(d7.Traits$Second_Treament) # convert to factor before numeric
d7.Traits$Second_Treament <- as.numeric(d7.Traits$Second_Treament)# now convet to numeric for WGCNA

d7.Traits$Group <- as.factor(d7.Traits$Group) # convert to factor before numeric
d7.Traits$Group <- as.numeric(d7.Traits$Group)# now convet to numeric for WGCNA


moduleTraitCor = cor(MEs, d7.Traits, use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);


#=====================================================================================
#
#  CREATE HEATMAP
#
#=====================================================================================


sizeGrWindow(10,10)
# Will display correlations and their p-values
textMatrix =  paste(signif(moduleTraitCor, 2), "\n(",
                    signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(9, 9, 6, 3));
# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(d7.Traits),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = TRUE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.6,
               zlim = c(-.7,.7),
               main = paste("DAY 7: Module-trait relationships"))

# The analysis identifies the several significant 
# module{trait associations. We will concentrate on TAOC as the trait
# of interest.

# concentrate on mean umol CRE g protein as the trait of interest (strongest module yellow)


#=====================================================================================
#
#  TAOC in PINK module! CREATE DATAFRAMES 
#
#=====================================================================================
# Gene relationship to trait and important modules: 
# Gene Significance and Module Membership

# We quantify associations of individual genes with our trait of interest (TAOC)

# Define variable weight containing the weight column of d7.Traits
TAOC = as.data.frame(d7.Traits$mean.µmol.CRE.g.protein); # -0.53 in yellow module
names(TAOC) = "TAOC"

# names (colors) of the modules
modNames = substring(names(MEs), 3) # name all the modules, from 3rd character on (first two are ME)

geneModuleMembership = as.data.frame(cor(d7.data, MEs, use = "p"));
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));

names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");

# TAOC 
TAOC_geneTraitSignificance = as.data.frame(cor(d7.data, TAOC, use = "p"));
TAOC_GSPvalue = as.data.frame(corPvalueStudent(as.matrix(TAOC_geneTraitSignificance), nSamples));

names(TAOC_geneTraitSignificance) = paste("GS.taoc", names(TAOC), sep="");
names(TAOC_GSPvalue) = paste("p.GS.taoc", names(TAOC), sep="");



#  PLOT TAOC IN THE PINK MODULE

module = "pink"
column = match(module, modNames);
moduleGenes = moduleColors==module;

sizeGrWindow(7, 7);
par(mfrow = c(1,1));
pink.TAOC <- verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(TAOC_geneTraitSignificance[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Gene significance for TOAC",
                   main = paste("Day7 TAOC 'pink': Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = "red")

#  PLOT TAOC IN THE PINK MODULE

module = "green"
column = match(module, modNames);
moduleGenes = moduleColors==module;

sizeGrWindow(7, 7);
par(mfrow = c(1,1));
green.TAOC <- verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(TAOC_geneTraitSignificance[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Gene significance for TOAC",
                   main = paste("Day7 TAOC 'green': Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = "green")

#  PLOT TAOC IN THE PINK MODULE

module = "blue"
column = match(module, modNames);
moduleGenes = moduleColors==module;

sizeGrWindow(7, 7);
par(mfrow = c(1,1));
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(TAOC_geneTraitSignificance[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Gene significance for TOAC",
                   main = paste("Day7 TAOC 'blue': Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = "blue")



#=====================================================================================
#
#    CALL THE GENES OF INTEREST IN THE MODULES (i.e. yellow and tan - refer to heatmap)
#
#=====================================================================================


length(names(d7.data)[moduleColors=="pink"]) # 65 total genes in the pink module
length(names(d7.data)[moduleColors=="green"]) # 149 total genes in the green module
length(names(d7.data)[moduleColors=="blue"]) # 510 total genes in the blue module

#=====================================================================================
#
#  Load in the annotation 
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


#   TAOC dataframe --------------------------------------------------------------------------- # 
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
geneOrder = order(geneInfo_TAOC$moduleColor, -abs(geneInfo_TAOC$GS.taocTAOC));
geneInfo_TAOC = geneInfo_TAOC[geneOrder, ]
View(geneInfo_TAOC)


#=====================================================================================
#
#  Write csv for the modules and corresponding raw read counts
#
#=====================================================================================
# call the module of interest for follow-up GO analysis 
# Yellow module and TAOC 
d7.WGCNA_TAOC.Pink.Mod.Gen <- geneInfo_TAOC %>%  dplyr::filter(moduleColor %in% 'pink')
d7.WGCNA_TAOC.Green.Mod.Gen <- geneInfo_TAOC %>%  dplyr::filter(moduleColor %in% 'green')
d7.WGCNA_TAOC.Blue.Mod.Gen <- geneInfo_TAOC %>%  dplyr::filter(moduleColor %in% 'blue')

write.csv(d7.WGCNA_TAOC.Pink.Mod.Gen, file = "Analysis/Output/WGCNA/d7.WGCNA_TAOC.Pink.Mod.Gen.csv")
write.csv(d7.WGCNA_TAOC.Green.Mod.Gen, file = "Analysis/Output/WGCNA/d7.WGCNA_TAOC.Green.Mod.Gen.csv")
write.csv(d7.WGCNA_TAOC.Blue.Mod.Gen, file = "Analysis/Output/WGCNA/d7.WGCNA_TAOC.Blue.Mod.Gen.csv")



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

# Load data 
d7_WGCNA_PinkMod <- read.csv("Analysis/Output/WGCNA/d7.WGCNA_TAOC.Pink.Mod.Gen.csv")
d7_WGCNA_GreenMod <- read.csv("Analysis/Output/WGCNA/d7.WGCNA_TAOC.Green.Mod.Gen.csv")
d7_WGCNA_BlueMod <- read.csv("Analysis/Output/WGCNA/d7.WGCNA_TAOC.Blue.Mod.Gen.csv")

d7_counts_filtered <- read.csv("Analysis/Data/filtered_counts/day7.counts.filtered_10cpm50perc.csv")
Master.Treatment_Phenotype.data <- read.csv(file="Analysis/Data/ Master_Phyenotype.and.Exp.Treatment_Metadata.csv", sep=',', header=TRUE)

# Prep data to merge
d7_WGCNA_PinkMod <- d7_WGCNA_PinkMod[-c(1:2)] # ommit the redundant gene name columns buut keep 'geneSymbol'
dim(d7_WGCNA_PinkMod) #  65 26

names(d7_counts_filtered)[1] <- "geneSymbol" # remname the gene name column to 'geneSymbol'
dim(d7_counts_filtered) #  8548   37
View(d7_counts_filtered)


# transform count data (vst in DESeq2)
d7_counts_filtered.matrix <- data.matrix(d7_counts_filtered[,-1]) # create matrix - remove 'geneSymbol'
d7_counts_matrix.vst <- vst(d7_counts_filtered.matrix)
View(d7_counts_matrix.vst)
d7_vst_counts <- data.frame(d7_counts_matrix.vst)
d7_vst_counts$geneSymbol <- d7_counts_filtered[,1] # insert 'geneSymbol'
View(d7_vst_counts)






###################################### #
# Pink Module 
###################################### #




# merge
d7_pink.module <-  merge(d7_WGCNA_PinkMod, d7_vst_counts, by = "geneSymbol")
dim(d7_pink.module) # 65 62
View(d7_pink.module)

# total dataset ---
colnames(d7_pink.module)

# Plot the vst transformed by treatment and by Total protein to visualize the effect of the pink module genes significant in WGCNA 
# vst read count date - narrow the columns - reshape and rename
d7.pink.module.vst <- d7_pink.module[,c(1, 27:62)]
colnames(d7.pink.module.vst) # check if you have the gene name and all samples 
d7.pink.module.vst_MELT <- melt(d7.pink.module.vst, id=("geneSymbol")) # melt using reshape2

names(d7.pink.module.vst_MELT)[(2:3)] <- c('Sample.Name', 'vst_Expression') # change column names

unique(Master.Treatment_Phenotype.data$Date)
# experiment treatment and total protein data - narrow the columns 
exp.phys_data <- Master.Treatment_Phenotype.data %>% 
  dplyr::filter(Date %in% 20190731) %>%  # filter out Day 21 data only 
  dplyr::select(c('Sample.Name', 'All_Treatment', 'Primary_Treatment', 'Second_Treament', 'mean.µmol.CRE.g.protein')) # call select columns - Day 7 had primary and second treaet - TAOC is the WGCNA of interest

# merge by common row values 'Sample.Name'
merged_pinkModvst_expdata <- merge(d7.pink.module.vst_MELT, exp.phys_data, by ='Sample.Name')

meanEXp <- merged_pinkModvst_expdata %>% 
  select(c('Sample.Name','vst_Expression','Primary_Treatment', 'Second_Treament')) %>% 
  group_by(Sample.Name, Primary_Treatment, Second_Treament) %>%
  dplyr::summarize(mean.vstExp = mean(vst_Expression), 
                   sd.vsdtExp = sd(vst_Expression),
                   na.rm=TRUE)
# summarize datasets further by treatment period
# remember:this is a mean of a mean!! First we complete mean vst exp by sample id (compiling all red module genes) - next all sample IDs by the treatment period (below
# I will use these for mean SE plots 
# Primary treatment
meanEXp_Summary.Prim <- meanEXp %>% 
  group_by(Primary_Treatment) %>%
  dplyr::summarize(mean = mean(mean.vstExp), 
                   sd = sd(mean.vstExp),
                   n = n(), 
                   se = sd/sqrt(n))
# Second treatment
meanEXp_Summary.Sec <- meanEXp %>% 
  group_by(Second_Treament,Primary_Treatment) %>%
  dplyr::summarize(mean = mean(mean.vstExp), 
                   sd = sd(mean.vstExp),
                   n = n(), 
                   se = sd/sqrt(n))



# Mean SE plots --------- #




# The errorbars overlapped, so use position_dodge to move them horizontally
pd <- position_dodge(0.1) # move them .05 to the left and right


# Primary treatment mean sd plot
Primary.vst <- ggplot(meanEXp_Summary.Prim, aes(x=Primary_Treatment, y=mean, fill=Primary_Treatment)) +  # , colour=supp, group=supp))
  theme_bw() +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), colour="black", width=.1, position=pd) +
  geom_line(position=pd) +
  geom_point(position=pd, size = 4, shape=21) +            
  xlab("Primary pCO2 treatment") +
  ylab(NULL) +                 # note the mean was first by sample ID THEN by treatment
  scale_fill_manual(values=c("#56B4E9","#E69F00")) +
  # scale_color_manual(values=c("#56B4E9","#E69F00")) +
  ggtitle("Day 7 WGCNA 'pink' Module VST GeneExp") +
  # expand_limits(y=0) +                                                    # Expand y range
  scale_y_continuous(limits=c(3, 5.5)) +
  theme(text = element_text(size=15))


# Second treatment mean sd plot
Sec.vst <- ggplot(meanEXp_Summary.Sec, aes(x=Second_Treament, y=mean, fill=Primary_Treatment, group=Primary_Treatment)) +
  theme_bw() +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), colour="black", width=.1, position=pd) +
  geom_line(position=pd) +
  geom_point(position=pd, size = 4, shape=21) +            
  xlab("Second pCO2 treatment") +
  ylab("Pink Module VST Gene Expression (Mean +/- SE)") +                 # note the mean was first by sample ID THEN by treatment
  scale_fill_manual(values=c("#56B4E9","#E69F00")) +
  # scale_color_manual(values=c("#56B4E9","#E69F00")) +
  # ggtitle("Day 21 WGCNA red' Module VST GeneExp") +
  # expand_limits(y=0) +                                                    # Expand y range
  scale_y_continuous(limits=c(3, 5.5)) +
  theme(text = element_text(size=15))


# Total protein by vst expression
meanExp_TAOC <- merged_pinkModvst_expdata %>% 
  group_by(Sample.Name) %>%
  dplyr::summarize(mean.Exp = mean(vst_Expression), 
                   sd.Exp = sd(vst_Expression),
                   mean.TAOC = mean(mean.µmol.CRE.g.protein), 
                   sd.TAOC = sd(mean.µmol.CRE.g.protein))
meanExp_TAOC # notice cd total protein is 0 - this is becasue total protein was the same value regarless of the gene.ID within sample

meanExp_TAOC_2 <- merge(meanExp_TAOC,meanEXp, by = "Sample.Name")

# two obvious outliers > 50
ggplot(meanExp_TAOC, aes(x=mean.TAOC, y=mean.Exp)) +
  geom_point(shape=1) +    # Use hollow circles
  geom_smooth()            # Add a loess smoothed fit curve with confidence region

meanExp_total.protein.OM <- meanExp_TAOC %>% filter(mean.TAOC > 0.5) # remove the outliers and plot again

ggplot(meanExp_total.protein.OM, aes(x=mean.TAOC, y=mean.Exp)) +
  geom_point(shape=1) +    # Use hollow circles
  geom_smooth()            # Add a loess smoothed fit curve with confidence region


# Total protein and expression plot
TotalProtein_vst.Exp <- ggplot(meanExp_total.protein.OM, aes(x=mean.totalprotein, y=mean.Exp)) +
  geom_point(shape=1) +    # Use hollow circles
  geom_smooth()   +        # Add a loess smoothed fit curve with confidence region
  theme_bw() +
  scale_y_continuous(limits=c(3, 5.5)) +
  theme(text = element_text(size=15))


#  grid arrange the three plots 
# save plot
png("Analysis/Output/WGCNA/Day7_vst.Exp_PinkModule.png", 500, 1000, pointsize=20)
ggarrange(Primary.vst, Sec.vst,
          plotlist = NULL,
          ncol = 1,
          nrow = 4,
          labels = NULL)
dev.off()




