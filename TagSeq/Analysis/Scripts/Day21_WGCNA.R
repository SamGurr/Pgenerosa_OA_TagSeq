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
# Tagaseq filtered counts 
day21counts.matrix <- read.csv(file="Analysis/Data/filtered_counts/day21.counts.filtered_10cpm50perc.csv", sep=',', header=TRUE)
# Treatment and Phenotype data
Master.Treatment_Phenotype.data <- read.csv(file="Analysis/Data/ Master_Phyenotype.and.Exp.Treatment_Metadata.csv", sep=',', header=TRUE)
d21.Treatment_Phenotype.data     <- Master.Treatment_Phenotype.data %>%  dplyr ::filter(Date %in% 20190814) # split for day 21 data 

# check the data 

# Dau 21 
dim(day21counts.matrix) # 8421 rows - 62 samples not counting 'x' and 'Gene.ID'
dim(d21.Treatment_Phenotype.data) # 62 total rows unique for each sample ID

# The following setting is important, do not omit.
options(stringsAsFactors = FALSE)

# =================================================================================== #
#
#
# Day 21 WGCNA 
#
#
# =================================================================================== #
  
dim(day21counts.matrix) # 8421  genes  63 samples; should have 62 columns of just samples 
(day21counts.matrix)[1] # ommit this and transpose in the next line 
d21.data = as.data.frame(t(day21counts.matrix[, -(1)])) # ommit all columns but samples and transpose
# fix(d21.data) # view the data - as you see all columns are now genes but filled as V1, V2, V3...
names(d21.data) = day21counts.matrix$X # assigns column names (previous jsut numbered) as the gene ID 
rownames(d21.data) = names(day21counts.matrix)[-(1)]; # assigns the row names as the sample ID


gsg = goodSamplesGenes(d21.data, verbose = 3);  # We first check for genes and samples with too many missing values:
gsg$allOK # If the statement returns TRUE, all genes have passed the cuts. 


# determine outlier samples 
sampleTree = hclust(dist(d21.data), method = "average") # Next we cluster the samples (in contrast to clustering genes that will come later)  to see if there are any obvious outliers.
sizeGrWindow(12,9) # The user should change the dimensions if the window is too large or too small.
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5, 
     cex.axis = 1.5, cex.main = 2) # appears there are two outliers SG55 and SG 105; can remove by hand or an automatic appraoch 
abline(h = 40000, col = "red") # add line to plot to show the cut-off od outlier samples (40000) SG105 and SG55
clust = cutreeStatic(sampleTree, cutHeight = 40000, minSize = 10) # Determine cluster under the line
table(clust) # 0 = cut; 1 = kept; says it will cut 2 and save 60 
keepSamples = (clust==1) # call the 60 from clust at position 1
d21.data = d21.data[keepSamples, ]
nGenes = ncol(d21.data) # number of genes == 8421 
nSamples = nrow(d21.data) # number of samples == 60

# Phenotype trait and Treatment data
dim(d21.Treatment_Phenotype.data) # 62 rows and 14 columns
names(d21.Treatment_Phenotype.data) # look at the 14 columns 
d21.Treatment_Phenotype.data = d21.Treatment_Phenotype.data[, -c(1:3, 5)]; # remove columns that hold information we do not need.
dim(d21.Treatment_Phenotype.data) # 62 rows and 10 columns
d21.Treatment_Phenotype.data$Group <- paste(d21.Treatment_Phenotype.data$Primary_Treatment, 
                                            d21.Treatment_Phenotype.data$Second_Treament,  
                                            d21.Treatment_Phenotype.data$Third_Treatment,
                                            sep ="")


# Form a data frame analogous to expression data that will hold the clinical traits.
d21.Samples = rownames(d21.data);# start new variable 'd21.Samples' calling the row names of the gene data (sample as 'Mouse' in trait data)
TreatRows = match(d21.Samples, d21.Treatment_Phenotype.data$Sample.Name); # match the names 
d21.Traits = d21.Treatment_Phenotype.data[TreatRows, -1]; # removes the row numbers
rownames(d21.Traits) = d21.Treatment_Phenotype.data[TreatRows, 1]; # inserts the new traitRows - matches sample treatments
dim(d21.Traits) # 60 Samples 10 columns (primary, second, and third treatment)
all(rownames(d21.Traits) == rownames(d21.data))  # should be TRUE
fix(d21.Traits) # view the data


# Form a data frame analogous to expression data that will hold the clinical traits.
# d21.Samples = rownames(d21.data);# start new variable 'd21.Samples' calling the row names of the gene data (sample as 'Mouse' in trait data)
# TreatRows = match(d21.Samples, d21.treatment$Sample.Name); # match the names 
# d21.Traits = as.data.frame(d21.treatment[TreatRows, -1]); # removes the row numbers
# d21.Traits
# rownames(d21.Traits) = d21.treatment[TreatRows, 1]; # inserts the new traitRows - matches sample treatments
# colnames(d21.Traits)[1] <- "Treatment"
# d21.Traits$Treatment <- as.factor(d21.Traits$Treatment)
# dim(d21.Traits) # 60 Samples 3 columns (primary, second, and third treatment)
# all(rownames(d21.Traits) == rownames(d21.data))  # should be TRUE



# Re-cluster samples
sampleTree2 = hclust(dist(d21.data), method = "average")
traitColors = labels2colors(d21.Traits); # Convert traits to a color representation: white means low, red means high, grey means missing entry
plotDendroAndColors(sampleTree2, traitColors, # Plot the sample dendrogram and the colors underneath.
                    groupLabels = names(d21.Traits), 
                    main = "Sample dendrogram and trait heatmap")

save(d21.data, d21.Traits, file = "Analysis/Output/WGCNA/d.21-dataInput.RData")



# soft threshold ============================================================================== #
# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to=20, by=2))
# Call the network topology analysis function
sft = pickSoftThreshold(d21.data, powerVector = powers, verbose = 5) #...wait for this to finish
# pickSoftThreshold 
#  performs the analysis of network topology and aids the
# user in choosing a proper soft-thresholding power.
# The user chooses a set of candidate powers (the function provides suitable default values)
# function returns a set of network indices that should be inspected


# pickSoftThreshold ------------------------------------------------------------------------- #
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

# NOTE: you will get this message "Error: REAL() can only be applied to a 'numeric', not a 'integer'"
# if your data is not all as numeric - set as numeric with lapply below BEFORE running 'blockwiseModules'
total_cols <- length(colnames(d21.data))
d21.data[2:total_cols] = lapply(d21.data[2:total_cols], as.numeric)

net = blockwiseModules(d21.data, power = 9,
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
# note 0 is reserved for genes outside of all modules (767 genes) 

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
nGenes = ncol(d21.data); # 8421
nSamples = nrow(d21.data); # 60
# Recalculate MEs with color labels
MEs0 = moduleEigengenes(d21.data, moduleColors)$eigengenes 
MEs = orderMEs(MEs0) # reorders the columns (colors/modules)

# change chanracter treatments to integers
d21.Traits$Primary_Treatment <- as.factor(d21.Traits$Primary_Treatment)
d21.Traits$Primary_Treatment <- as.numeric(d21.Traits$Primary_Treatment)

d21.Traits$Second_Treament <- as.factor(d21.Traits$Second_Treament)
d21.Traits$Second_Treament <- as.numeric(d21.Traits$Second_Treament)

d21.Traits$Third_Treatment <- as.factor(d21.Traits$Third_Treatment)
d21.Traits$Third_Treatment <- as.numeric(d21.Traits$Third_Treatment)

d21.Traits$Group <- as.factor(d21.Traits$Group)
d21.Traits$Group <- as.numeric(d21.Traits$Group)

moduleTraitCor = cor(MEs, d21.Traits, use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);


#=====================================================================================
#
#  labeled heatmap
#
#=====================================================================================


sizeGrWindow(10,10)
# Will display correlations and their p-values
textMatrix =  paste(signif(moduleTraitCor, 2), "\n(",
                    signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(8, 9.5, 5, 3));
# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(d21.Traits),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = TRUE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.9,
               zlim = c(-0.5,0.5),
               main = paste("Module-trait relationships"))

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

# Define variable weight containing the weight column of d7.Traits
total_protein = as.data.frame(d21.Traits$mean.mgProt.mgAFDW); # -0.53 in yellow module
names(total_protein) = "total_protein"

# names (colors) of the modules
modNames = substring(names(MEs), 3) # name all the modules, from 3rd character on (first two are ME)

geneModuleMembership = as.data.frame(cor(d21.data, MEs, use = "p"));
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));

names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");

# Total.protein
total.protein_geneTraitSignificance = as.data.frame(cor(d21.data, total_protein, use = "p"));
total.protein_GSPvalue = as.data.frame(corPvalueStudent(as.matrix(total.protein_geneTraitSignificance), nSamples));

names(total.protein_geneTraitSignificance) = paste("GS.protein", names(total_protein), sep="");
names(total.protein_GSPvalue) = paste("p.GS.protein", names(total_protein), sep="");



#  PLOT Total.protein IN THE RED MODULE
unique(moduleColors)
module = "red"
column = match(module, modNames);
moduleGenes = moduleColors==module;

sizeGrWindow(7, 7);
par(mfrow = c(1,1));
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(total.protein_geneTraitSignificance[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Gene significance for TOAC",
                   main = paste("Day21 total protein 'red': Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = "black")


#=====================================================================================
#
#    CALL THE GENES OF INTEREST IN THE MODULES (i.e. red - refer to heatmap)
#
#=====================================================================================


length(names(d21.data)[moduleColors=="red"]) # 104 total genes in the red module


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
names(total.protein_geneTraitSignificance)
names(total.protein_GSPvalue)
geneInfo_total.protein = data.frame(substanceBXH = probes,
                           geneSymbol = annot$V1[probes2annot],
                           #LocusLinkID = annot$LocusLinkID[probes2annot],
                           moduleColor = moduleColors,
                           Uniprot = annot$V5[probes2annot],
                           HGNC = annot$V6[probes2annot],
                           GO.terms = annot$V8[probes2annot],
                           GO.description = annot$V9[probes2annot],
                           total.protein_geneTraitSignificance, # call this specific to the module and trait of interest
                           total.protein_GSPvalue)              # call this specific to the module and trait of interest
View(geneInfo_total.protein)
modOrder = order(-abs(cor(MEs, total_protein, use = "p"))); # order by the strength of the correlation between module and trait values for each sample

for (mod in 1:ncol(geneModuleMembership)) # Add module membership information in the chosen order
{
  oldNames = names(geneInfo_total.protein)
  geneInfo_total.protein = data.frame(geneInfo_total.protein, geneModuleMembership[, modOrder[mod]], 
                             MMPvalue[, modOrder[mod]]);
  names(geneInfo_total.protein) = c(oldNames, paste("MM.", modNames[modOrder[mod]], sep=""),
                           paste("p.MM.", modNames[modOrder[mod]], sep=""))
}
# Order the genes in the geneInfo variable first by module color, then by geneTraitSignificance
geneOrder = order(geneInfo_total.protein$moduleColor, -abs(geneInfo_total.protein$GS.proteintotal_protein));
geneInfo_protein = geneInfo_total.protein[geneOrder, ]
View(geneInfo_protein)


#=====================================================================================
#
#  Write csv for the modules and corresponding raw read counts
#
#=====================================================================================
# call the module of interest for follow-up GO analysis 

# RED module 

# Total protein 
d21.WGCNA_TotalProtein.RED.Mod.Gen <- geneInfo_protein %>%  dplyr::filter(moduleColor %in% 'red')
write.csv(d21.WGCNA_TotalProtein.RED.Mod.Gen, file = "Analysis/Output/WGCNA/d21.WGCNA_TotalProtein.RED.Mod.Gen.csv")


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
d21_WGCNA_redMod <- read.csv("Analysis/Output/WGCNA/d21.WGCNA_TotalProtein.RED.Mod.Gen.csv")
d21_counts_filtered <- read.csv("Analysis/Data/filtered_counts/day21.counts.filtered_10cpm50perc.csv")
Master.Treatment_Phenotype.data <- read.csv(file="Analysis/Data/ Master_Phyenotype.and.Exp.Treatment_Metadata.csv", sep=',', header=TRUE)

# Prep data to merge
d21_WGCNA_redMod <- d21_WGCNA_redMod[-c(1:2)] # ommit the redundant gene name columns buut keep 'geneSymbol'
dim(d21_WGCNA_redMod) #  104  30

names(d21_counts_filtered)[1] <- "geneSymbol" # remname the gene name column to 'geneSymbol'
dim(d21_counts_filtered) #  8421   63
View(d21_counts_filtered)


# transform count data (vst in DESeq2)
d21_counts_filtered.matrix <- data.matrix(d21_counts_filtered[,-1]) # create matrix - remove 'geneSymbol'
d21_counts_matrix.vst <- vst(d21_counts_filtered.matrix)
View(d21_counts_matrix.vst)
d21_vst_counts <- data.frame(d21_counts_matrix.vst)
d21_vst_counts$geneSymbol <- d21_counts_filtered[,1] # insert 'geneSymbol'
View(d21_vst_counts)
# merge
d21_red.module <-  merge(d21_WGCNA_redMod, d21_vst_counts, by = "geneSymbol")
dim(d21_red.module) # 104  92
View(d21_red.module)

# total dataset
d21_red.module

# Plot the vst transformed by treatment and by Total protein to visualize the effect of the red module genes significant in WGCNA 
# vst read count date - narrow the columns - reshape and rename
d21.red.module.vst <- d21_red.module[,c(1, 31:92)]
d21.red.module.vst_MELT <- melt(d21.red.module.vst, id=("geneSymbol")) # melt using reshape2
 'Second_Treament'
names(d21.red.module.vst_MELT)[(2:3)] <- c('Sample.Name', 'vst_Expression') # change column names

# experiment treatment and total protein data - narrow the columns 
exp.phys_data <- Master.Treatment_Phenotype.data %>% 
  dplyr::filter(Date %in% 20190814) %>%  # filter out Day 21 data only 
  dplyr::select(c('Sample.Name', 'All_Treatment', 'Primary_Treatment', 'Second_Treament', 'Third_Treatment', 'mean.mgProt.mgAFDW')) # call select columns

# merge by common row values 'Sample.Name'
merged_redModvst_expdata <- merge(d21.red.module.vst_MELT, exp.phys_data, by ='Sample.Name')

meanEXp <- merged_redModvst_expdata %>% 
                        select(c('Sample.Name','vst_Expression','Primary_Treatment', 'Second_Treament','Third_Treatment')) %>% 
                        group_by(Sample.Name, Primary_Treatment, Second_Treament, Third_Treatment) %>%
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

# Third treatment
meanEXp_Summary.Third <- meanEXp %>% 
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
                  ggtitle("Day 21 WGCNA red' Module VST GeneExp") +
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
                  ylab("Red Module VST Gene Expression (Mean +/- SE)") +                 # note the mean was first by sample ID THEN by treatment
                  scale_fill_manual(values=c("#56B4E9","#E69F00")) +
                  # scale_color_manual(values=c("#56B4E9","#E69F00")) +
                  # ggtitle("Day 21 WGCNA red' Module VST GeneExp") +
                  # expand_limits(y=0) +                                                    # Expand y range
                  scale_y_continuous(limits=c(3, 5.5)) +
                  theme(text = element_text(size=15))


# Third treatment mean sd plot
meanEXp_Summary.Third$Prim.Sec_Treatment <- paste(meanEXp_Summary.Third$Primary_Treatment, meanEXp_Summary.Third$Second_Treament, sep='')
Third.vst <- ggplot(meanEXp_Summary.Third, aes(x=Third_Treatment, y=mean, fill=Primary_Treatment, group=Prim.Sec_Treatment, alpha = Second_Treament)) +
                  theme_bw() +
                  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), colour="black", width=.1, position=pd,alpha = 1) +
                  geom_line(position=pd, colour="black",alpha = 1) +
                  geom_point(position=pd, size = 4, shape=21) +            
                  xlab("Third pCO2 treatment") +
                  ylab(NULL) +                 # note the mean was first by sample ID THEN by treatment
                  scale_fill_manual(values=c("#56B4E9","#E69F00")) +
                  # scale_color_manual(values=c("#56B4E9","#E69F00")) +
                  # ggtitle("Day 21 WGCNA red' Module VST GeneExp") +
                  scale_y_continuous(limits=c(3, 5.5)) +
                  theme(text = element_text(size=15))

# Total protein by vst expression
meanExp_total.protein <- merged_redModvst_expdata %>% 
  group_by(Sample.Name) %>%
  dplyr::summarize(mean.Exp = mean(vst_Expression), 
                   sd.Exp = sd(vst_Expression),
                   mean.totalprotein = mean(mean.mgProt.mgAFDW), 
                   sd.totalprotein = sd(mean.mgProt.mgAFDW))
meanExp_total.protein # notice cd total protein is 0 - this is becasue total protein was the same value regarless of the gene.ID within sample


# two obvious outliers > 50
ggplot(meanExp_total.protein, aes(x=mean.totalprotein, y=mean.Exp)) +
  geom_point(shape=1) +    # Use hollow circles
  geom_smooth()            # Add a loess smoothed fit curve with confidence region

meanExp_total.protein.OM <- meanExp_total.protein %>% filter(mean.totalprotein < 50) # remove the outliers and plot again

# Total protein and expression plot
TotalProtein_vst.Exp <- ggplot(meanExp_total.protein.OM, aes(x=mean.totalprotein, y=mean.Exp)) +
  geom_point(shape=1) +    # Use hollow circles
  geom_smooth()   +        # Add a loess smoothed fit curve with confidence region
  theme_bw() +
  scale_y_continuous(limits=c(3, 5.5)) +
  theme(text = element_text(size=15))


#  grid arrange the three plots 
# save plot
png("Analysis/Output/WGCNA/Day21_vst.Exp_redModule.png", 500, 1000, pointsize=20)
ggarrange(Primary.vst, Sec.vst, Third.vst, TotalProtein_vst.Exp, 
          plotlist = NULL,
          ncol = 1,
          nrow = 4,
          labels = NULL)
dev.off()



  
