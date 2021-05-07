---
# title: "apriori_GOIs"
# author: "Samuel Gurr"
# date: "March 7, 2021"
---
  
# LOAD PACKAGES
library(dplyr)
library(DESeq2)
library(tidyselect)
library(reshape2)
library(ggplot2)
library(gridExtra)
library(ggpubr)
library(stringr)

# SET WORKING DIRECTORY AND LOAD DATA 
setwd("C:/Users/samjg/Documents/My_Projects/Pgenerosa_TagSeq_Metabolomics/TagSeq/")

# LOAD DATA
# count data
as.data.frame(unique(rbind(day7.counts.matrix[,1],day14.counts.matrix[,1]),day21.counts.matrix[,1]))
day14.counts.matrix
day21.counts.matrix


day7.counts.matrix <- read.csv(file="Analysis/Data/Filtered_Counts/10cpm_50perc/day7.counts.filtered_10cpm50perc.csv", sep=',', header=TRUE)
day14.counts.matrix <- read.csv(file="Analysis/Data/Filtered_Counts/10cpm_50perc/day14.counts.filtered_10cpm50perc.csv", sep=',', header=TRUE)
day21.counts.matrix <- read.csv(file="Analysis/Data/Filtered_Counts/10cpm_50perc/day21.counts.filtered_10cpm50perc.csv", sep=',', header=TRUE)
# trait data
Master.Treatment_Phenotype.data <- read.csv(file="Analysis/Data/Experiment_Metadata/Master_Phyenotype.and.Exp.Treatment_Metadata.csv", sep=',', header=TRUE)
# WGCNA results 
D7_WGCNA     <- read.csv(file="Analysis/Output/WGCNA/Day7/d7.WGCNA_ModulMembership.csv", sep=',', header=TRUE)    %>% dplyr::mutate(Day = "Day7") %>%  dplyr::select(c('geneSymbol','moduleColor', 'Day'))
D14_WGCNA    <- read.csv(file="Analysis/Output/WGCNA/Day14/d14.WGCNA_ModulMembership.csv", sep=',', header=TRUE) %>% dplyr::mutate(Day = "Day14") %>%  dplyr::select(c('geneSymbol','moduleColor', 'Day'))
D21_WGCNA    <- read.csv(file="Analysis/Output/WGCNA/Day21/d21.WGCNA_ModulMembership.csv", sep=',', header=TRUE) %>% dplyr::mutate(Day = "Day21") %>%  dplyr::select(c('geneSymbol','moduleColor', 'Day'))
WGCNA_Master <- rbind(D7_WGCNA, D14_WGCNA, D21_WGCNA)

# DESEq2 results 
D7_DEGs_PrimEffect  <- read.csv(file="Analysis/Output/DESeq2/10cpm/Day7/DE_Day7_Primary.csv", sep=',', header=TRUE)  %>% dplyr::select(c('Row.names', 'log2FoldChange', 'pvalue', 'up', 'down'))
D14_DEGs_PrimEffect <- read.csv(file="Analysis/Output/DESeq2/10cpm/Day14/DE_Day14_Primary.csv", sep=',', header=TRUE) %>% dplyr::select(c('Row.names', 'log2FoldChange', 'pvalue', 'up', 'down'))
D21_DEGs_PrimEffect <- read.csv(file="Analysis/Output/DESeq2/10cpm/Day21/DE_Day21_Primary.csv", sep=',', header=TRUE) %>% dplyr::select(c('Row.names', 'log2FoldChange', 'pvalue', 'up', 'down'))
D7_DEGs_PrimEffect$Day  <- "Day7"
D14_DEGs_PrimEffect$Day <- "Day14"
D21_DEGs_PrimEffect$Day <- "Day21"
DEGs_PrimEffect_10cpm <- rbind(D7_DEGs_PrimEffect, D14_DEGs_PrimEffect, D21_DEGs_PrimEffect) # bind dataframes 
names(DEGs_PrimEffect_10cpm)[1] <- "gene.ID" # change name of column 1
# DESEq2 summary table - run the a priori against this table of main effects to look for overable here 
D0_Main.DEGs.Master      <- read.csv(file="Analysis/Output/DESeq2/10cpm/Day0/ Day0.AllPrimaryDEGs_DESeq2results.csv", sep=',', header=TRUE)  
D7_Main.DEGs.Master      <- read.csv(file="Analysis/Output/DESeq2/10cpm/Day7/ Day7.AllMainEffectsDEGs_DESeq2results.csv", sep=',', header=TRUE)  
D14_Main.DEGs.Master     <- read.csv(file="Analysis/Output/DESeq2/10cpm/Day14/ Day14.AllMainEffectsDEGs_DESeq2results.csv", sep=',', header=TRUE)  
D21_Main.DEGs.Master     <- read.csv(file="Analysis/Output/DESeq2/10cpm/Day21/ Day21.AllMainEffectsDEGs_DESeq2results.csv", sep=',', header=TRUE)  
DESeq2_MainEffectsMaster <- rbind(D0_Main.DEGs.Master,D7_Main.DEGs.Master,D14_Main.DEGs.Master,D21_Main.DEGs.Master)
names(DESeq2_MainEffectsMaster)[2] <- 'geneSymbol' # match the WGCNA master file

#Panopea generosa - load .fna ('Geoduck_annotation') and foramt GO terms ('Geoduck_GOterms') and vectors
Geoduck_annotation <- read.delim2(file="C:/Users/samjg/Documents/My_Projects/Pgenerosa_TagSeq_Metabolomics/TagSeq/Seq_details/Panopea-generosa-genes-annotations.txt", header=F)
# build annotation file to merge with the mean LFC tables
annot.condenced <- Geoduck_annotation[,c(1,7)] # load just the PGEN ID and the putative gene terms
annot.condenced$Gene_term    <- sub(" \\(EC.*", "", annot.condenced$V7)  # str_extract(annot.condenced$V7, "[^(]+")
colnames(annot.condenced) <- c("GeneID", "Gene_term_ALL", "Gene_term")
annot.condenced[c(1:2),] # view or condenced table 

# ===================================================================================
#
# Assemble targets GOIs of the a prior hypothesis 
#
# ===================================================================================

# call genes of interact and plot their expression levelsin response to treatment 
# my a priori hypothesis can be found in the end of our bioRxiv in 2020 https://www.biorxiv.org/content/10.1101/2020.08.03.234955v1.full
# 'Stress-induced attentuation of mitchondrial function' == reverse electron transport (RET) overproduces ROS in absence of stress acclimation; I hypothesize
# that NADH dehydrogenase (complex 1) will have greater expression in Ambient vs. Moderate. I further hypothesize that uncoupling proteins and sirtuins will be 
# upregulated in this naive phenotype to alleviate oxidative stress
# 
# 'Adaptive mitchondrial shift under moderate stress'    == lower antioxidant needed in the acclimatized phenotype partially due to upregulated AOX pathway, 
# I hypoothesize that AOX will be upregulated in the moderate conditioned animals 
# plotCounts(dds, gene=which.min(res$padj), intgroup="condition")
# 

# Genes of interest and PGEN ID 

# MITOCHONDRIAL OXPHOS PLAYERS 
# AOX                    == PGEN_.00g108770
# NADH dehydrogenase     == PGEN_.00g299160
# Cytochrome c reductase == PGEN_00g275780
# Uncoupling protein     == PGEN_.00g063670, PGEN_.00g193030, PGEN_.00g230260 

# TRANASCRIPTIONAL REGULATION - proteins involved in methylation and histone modification(s)
# PGEN_.00g074890, PGEN_.00g041430 = 'Uncharacterized_methyltransferase-like_C25B8.10', 'S-adenosylmethionine_synthase_isoform_type-2'
# PGEN_.00g283000 - DNMT_1
# PGEN_.00g029420 - DNMT3A
# PGEN_.00g067800 - DNMT 3b
# PGEN_.00g327700, PGEN_.00g053100, PGEN_.00g053120,  - Histone-lysine_N-methyltransferase_SETD2, histone methyl transferase
# PGEN_.00g064910, PGEN_.00g064920 - positive correlation with histone methylation
# PGEN_.00g066460 - methyltransferase like protein
# PGEN_.00g283340, PGEN_.00g311570, PGEN_.00g320700, PGEN_.00g338440, PGEN_.00g272910 - histone acetyltransferase  [GO:0004402]

# SIRTUINS 
# PGEN_.00g048200 - sirtuin (SIR1)
# PGEN_.00g012340 - sirtuin (SIR2)
# PGEN_.00g149480 - siruitin (SIR4)
# PGEN_.00g153700 - sirtuin (SIR5)
# PGEN_.00g144540 - sirtuin (SIR6)
# PGEN_.00g033970 - sirtuin (SIR7) 

# OXIDATIVE STRESS
# PGEN_.00g010160, PGEN_.00g065700,PGEN_.00g257600- superoxide dismutase
# PGEN_.00g015070 - superoxide dismutase copper chaperone
# PGEN_.00g062450 - superoxide disutase activity Mn
# PGEN_.00g192250,PGEN_.00g180320, PGEN_.00g287800, PGEN_.00g293960 - glutathione peroxidase
# genes <- c('AOX',
#            'NADH_dehydrogenase', 'NADH_dehydrogenase_probable', 'NADH-ubiquinone_oxidoreductase_75kDasubunit', 'NADH_dehydrogenase_alpha_subcomplex_subunit10', 'NADH_dehydrogenase_iron-sulfur_protein8',
#            'Cytochrome_c_reductase','Uncoupling_protein_1','Uncoupling_protein_2','Uncoupling_protein_3',
#            'Uncharacterized_methyltransferase-like_C25B8.10', 'S-adenosylmethionine_synthase_isoform_type-2','DNMT_1','DNMT3A','DNMT_3b', 'Histone-lysine_N-methyltransferase_SETD2', 'HMT_1','HMT_2','positiv_hist_methyl','positiv_hist_methyl_2','methyltransferase','HAT_1','HAT_2','HAT_3','HAT_4','HAT_5',
#            'SIR1','SIR2','SIR4','SIR5','SIR6','SIR7',
#            'SOD_1','SOD_2','SOD_3',
#            'SOD_Cu_chaperone',
#            'SOD_Mn_act',
#            'glutathione_peroxidase_1','glutathione_peroxidase_2','glutathione_peroxidase_3','glutathione_peroxidase_4', 'glutathione_peroxidase_5','glutathione_peroxidase_6',
#            'titin', 'calpain', 'Ecto-NOX_disulfide-thiol_exchanger2')
# 
# 
# geneSymbol <- c('PGEN_.00g108770', 
#                 'PGEN_.00g299160', 'PGEN_.00g220490', 'PGEN_.00g266170', 'PGEN_.00g262440', 'PGEN_.00g133830',
#                 'PGEN_00g275780','PGEN_.00g063670', 'PGEN_.00g193030', 'PGEN_.00g230260', # MITCHONDRIAL PLAYERS 
#                 'PGEN_.00g074890', 'PGEN_.00g041430'  , 'PGEN_.00g283000','PGEN_.00g029420','PGEN_.00g067800', 'PGEN_.00g327700', 'PGEN_.00g053100','PGEN_.00g053120','PGEN_.00g064910', 'PGEN_.00g064920','PGEN_.00g066460', 'PGEN_.00g283340',  'PGEN_.00g311570', 'PGEN_.00g320700', 'PGEN_.00g338440','PGEN_.00g272910', # TRANASCRIPTIONAL REGULATION - proteins involved in methylation and histone modification(s)
#                 'PGEN_.00g048200', 'PGEN_.00g012340', 'PGEN_.00g149480', 'PGEN_.00g153700', 'PGEN_.00g144540', 'PGEN_.00g033970',# SIRTUINS 
#                 'PGEN_.00g010160', 'PGEN_.00g065700', 'PGEN_.00g257600', # OXIDATIVE STRESS
#                 'PGEN_.00g015070', # OXIDATIVE STRESS
#                 'PGEN_.00g062450', # OXIDATIVE STRESS
#                 'PGEN_.00g293960', 'PGEN_.00g287800', 'PGEN_.00g192250','PGEN_.00g180320', 'PGEN_.00g116940','PGEN_.00g049360',# OXIDATIVE STRESS
#                 'PGEN_.00g066340',  'PGEN_.00g014370', 'PGEN_.00g230250') # Growth genes (Hollie asked about these)
# target_GOIs <- data.frame(geneSymbol, genes)

View(Target_GOIs)

Target_GOIs <- dplyr::filter(annot.condenced, grepl('Dnmt3a|Dnmt3b|Dnmt1|S-adenosylmethionine|
                                                     |Jumonji|JmjC|Lysine-specific demethylase| 
                                                     |histone acetyltransferase|histone methyltransferase|Histone-lysine N-methyltransferase|
                                                     |ferric-chelate reductase|
                                                     |carbonic anhydrase|
                                                     |glutathione peroxidase|superoxide dismutase|Glutathione S-transferase|
                                                     |Mitogen-activated protein kinase|
                                                     |Mitochondrial uncoupling protein|ucp|alternative oxidase|NADH dehydrogenase|cytochrome reductase|cytochrome b-c1|
                                                     |ATP synthase|sir1|sir2|sir3|sirt4|sir5|sir6|sir7|NAD-dependent protein deacetylase|sirtuin', Gene_term, ignore.case = TRUE))
Target_GOIs$Translation <- ifelse(grepl("Dnmt3a|Dnmt3b|Dnmt1", Target_GOIs$Gene_term, ignore.case = TRUE), "DNA methyltransferase", 
                           ifelse(grepl("S-adenosylmethionine", Target_GOIs$Gene_term, ignore.case = TRUE), "SAM",
                           ifelse(grepl("acetyltransferase", Target_GOIs$Gene_term, ignore.case = TRUE), "Histone acetyltransferases",
                           ifelse(grepl("histone methyltransferase|Histone-lysine N-methyltransferase", Target_GOIs$Gene_term, ignore.case = TRUE), "Histone methyltransferases",
                           ifelse(grepl("Jumonji|JmjC|Lysine-specific demethylase", Target_GOIs$Gene_term, ignore.case = TRUE), "Histone demethylation",
                           ifelse(grepl("ferric", Target_GOIs$Gene_term, ignore.case = TRUE), "reduction of Fe3.2+ to Fe2+",
                           ifelse(grepl("glutathione|superoxide|Glutathione S-transferase", Target_GOIs$Gene_term, ignore.case = TRUE), "Oxidative stress response",
                           ifelse(grepl("Mitogen-activated protein kinase", Target_GOIs$Gene_term, ignore.case = TRUE), "Stress signaling (MAPK)",
                           ifelse(grepl("NADH dehydrogenase", Target_GOIs$Gene_term, ignore.case = TRUE), "Complex I",
                           ifelse(grepl("cytochrome b-c1", Target_GOIs$Gene_term, ignore.case = TRUE), "Complex III",
                           ifelse(grepl("cytochrome reductase", Target_GOIs$Gene_term, ignore.case = TRUE), "Complex IV",      
                           ifelse(grepl("Mitochondrial uncoupling protein|ucp", Target_GOIs$Gene_term, ignore.case = TRUE), "Mitochondrial uncoupling proteins",
                           ifelse(grepl("alternative oxidase", Target_GOIs$Gene_term, ignore.case = TRUE), "AOX",
                           ifelse(grepl("ATP synthase", Target_GOIs$Gene_term, ignore.case = TRUE), "ATP synthase",
                           ifelse(grepl("carbonic anhydrase", Target_GOIs$Gene_term, ignore.case = TRUE), "acid-base regulation",
                           ifelse(grepl("deacetylase|sirtuin|sirt", Target_GOIs$Gene_term, ignore.case = TRUE), "Sirtuins", "Other"))))))))))))))))
Target_GOIs$Gene_term <- gsub( "/", "", as.character(Target_GOIs$Gene_term)) # removes all occurances of '/' and replaces with a NULL for all gene terms 
Target_GOIs$Gene_term <- gsub( " ", "_", as.character(Target_GOIs$Gene_term)) # removes all occurances of '/' and replaces with a NULL for all gene terms 

# the gsub above allows the for loop of vstExp plots wihout directory error (i.e. thinks that the / in the gene name is a separate folder)
# View(Target_GOIs)

nrow(Target_GOIs %>% dplyr::filter(!Translation %in% 'Stress signaling (MAPK)')) # 105 total genes (in WGCNA (excluding the MAPK proteins))
nrow(Target_GOIs %>% dplyr::filter(Translation %in% 'acid-base regulation')) # 105 total genes (in WGCNA (excluding the MAPK proteins))


# ===================================================================================
# 
#  Loop to make table - look in DESeq2 results and WGCNa results! 
#
# ===================================================================================
# head(WGCNA_Master)
# head(DESeq2_MainEffectsMaster)
# head(target_GOIs)



# summary table of the target GOIs
Target_GOIs_Summary <- Target_GOIs %>% 
  group_by(Translation) %>% 
  tally() 
 
  
  
# WGCNA all target GOIs and the correspoonding WGCNA modules (only significant modules with treatment) in which they are found! 
Sig_WGCNA.modules <- WGCNA_Master %>%  
                    dplyr::mutate(Day_color = paste(Day, moduleColor, sep = '_')) %>% 
                    dplyr::filter(Day_color %in% c('Day7_yellow', 'Day7_green', 'Day7_brown',
                                                   'Day14_brown', 'Day14_black', 'Day14_pink', 'Day14_magenta',
                                                   'Day21_red', 'Day21_blue', 'Day21_yellow', 'Day21_turquoise', 'Day21_black', 'Day21_pink', 'Day21_magenta')) 
colnames(Sig_WGCNA.modules)[1] <- "GeneID"

# How many of the target GOIs are in the WGCNA analysis?
a_priori_WGCNA <- merge(Sig_WGCNA.modules,Target_GOIs[c(1,3:4)], by= "GeneID") # %>% dplyr::select('Sample.Name', 'Primary_Treatment', 'Second_Treament', 'Third_Treatment', any_of(target_GOIs$GeneID))
nrow(a_priori_WGCNA %>% dplyr::filter(!Translation %in% 'Stress signaling (MAPK)')) # 105 total genes (in WGCNA (excluding the MAPK proteins))

# Show only genes associated with Primary Treatment co-expression significance (modules of interest!)
a_priori_WGCNA_TALLY <- a_priori_WGCNA %>%  # summarize occurances of the 'Translation' column for each module_Day
  filter(Day_color %in% c('Day7_brown', 'Day7yellow',                             # primary effect modules day 7
                          'Day14_brown', 'Day14_black',                           # primary effect modules day 14
                          'Day21_blue', 'Day21_magenta', 'Day21_yellow')) %>%     # primary effect modules day 21
  group_by(Day_color,Translation) %>% 
  tally()  # summarize occurances of the 'Translation' column for each module_Day
View(a_priori_WGCNA)


# Tally what functions/moduels are present on which day 
a_priori_WGCNA_TALLY_Day <- a_priori_WGCNA %>%  # summarize occurances of the 'Translation' column for each module_Day
  filter(Day_color %in% c('Day7_brown', 'Day7yellow',                             # primary effect modules day 7
                          'Day14_brown', 'Day14_black',                           # primary effect modules day 14
                          'Day21_blue', 'Day21_magenta', 'Day21_yellow')) %>%     # primary effect modules day 21
  group_by(Day,Translation) %>% 
  tally()  # summarize occurances of the 'Translation' column for each module_Day
View(a_priori_WGCNA_TALLY_Day)

# tally all significant modules for genes with hgihest occurance (N == 3 means this gene occured in a significant module EVERY SAMPLING DAY!)
a_priori_WGCNA_TALLY_gene_Day <- a_priori_WGCNA %>%  # summarize occurances of the 'Translation' column for each module_Day
  dplyr::select(c('GeneID', 'Day', 'Gene_term', 'moduleColor')) %>% 
  dplyr::group_by(GeneID, Gene_term) %>% 
  tally()  # summarize occurances of the 'Translation' column for each module_Day
View(a_priori_WGCNA_TALLY_gene_Day)
# six total genes had 3 occurances (sig moduel assocaiation on days 7 14 and)
# PGEN_.00g041430 S-adenosylmethionine_synthase_isoform_type-2_(AdoMet_synthase_2)
# PGEN_.00g049360 Glutathione_peroxidase_3_(GPx-3)_(GSHPx-3)
# PGEN_.00g116940 Glutathione_peroxidase_(PcGPx)_(Se-PcGPx)
# PGEN_.00g180630 Glutathione_S-transferase_omega-1_(GSTO-1)
# PGEN_.00g219960 Probable_glutathione_S-transferase_7
# PGEN_.00g255890 Carbonic_anhydrase_7


# DESeq 2 main effects 
colnames(DESeq2_MainEffectsMaster)[2] <- "GeneID"
a_priori_DESeq2 <- merge(DESeq2_MainEffectsMaster, Target_GOIs[c(1,3:4)], by = "GeneID")
a_priori_DESeq2 <- a_priori_DESeq2[-2] # ommit the row number column
View(a_priori_DESeq2)


# save these summary data tables 
path_out.apriori = 'C:/Users/samjg/Documents/My_Projects/Pgenerosa_TagSeq_Metabolomics/TagSeq/Analysis/Output/a_priori_hypothesis/' # call path
write.csv(as.data.frame(Target_GOIs_Summary), paste(path_out.apriori,"apriori_TargetGOIs_Table.csv")) # write
write.csv(a_priori_DESeq2, paste(path_out.apriori,"apriori_DESeq2_AllMainEffects.csv")) # write
write.csv(a_priori_WGCNA, paste(path_out.apriori,"apriori_WGCNA_AllSigModules.csv")) # write
write.csv(as.data.frame(a_priori_WGCNA_TALLY), paste(path_out.apriori,"apriori_WGCNA_AllSigModules_Summary.csv")) # write



# ===================================================================================
#
# PLOTTING - DATA PREP FOR PLOTS
#
# ===================================================================================
# treatment data
d7.Treatment_data <- Master.Treatment_Phenotype.data %>%   dplyr::filter(Date %in% 20190731) %>%  dplyr::select('Sample.Name', 'Primary_Treatment', 'Second_Treament') # split for day 7 data 
d14.Treatment_data <- Master.Treatment_Phenotype.data  %>%  dplyr::filter(Date %in% 20190807) %>%  dplyr::select('Sample.Name', 'Primary_Treatment', 'Second_Treament')# split for day 7 data 
d21.Treatment_data <- Master.Treatment_Phenotype.data  %>%  dplyr::filter(Date %in% 20190814) %>%  dplyr::select('Sample.Name', 'Primary_Treatment', 'Second_Treament', 'Third_Treatment')# split for day 7 data 

# ================================================================================== #
# count data - day  7
d7.data = as.data.frame(t(day7.counts.matrix[, -(1)])) # ommit all columns but samples and transpose
names(d7.data) = day7.counts.matrix$X # assigns column names (previous jsut numbered) as the gene ID 
rownames(d7.data) = names(day7.counts.matrix)[-(1)]; # assigns the row names as the sample ID
d7.data_matrix <- data.frame(day7.counts.matrix[,-1], row.names=day7.counts.matrix[,1]) 
d7.data_matrix_t <- t(d7.data_matrix)
dds.d7 <- DESeqDataSetFromMatrix(countData = d7.data_matrix,  colData = d7.Treatment_data, design = ~ 1)
dds.d7_vst <- vst(dds.d7) # transform it vst
dds.d7_vst <- assay(dds.d7_vst) # call only the transformed coutns in the dds object
d7_vst <- t(dds.d7_vst) # transpose columns to rows and vice versa
# fix(d7_vst)

# =================================================================================== #
# count data - day  14
d14.data = as.data.frame(t(day14.counts.matrix[, -(1)])) # ommit all columns but samples and transpose
names(d14.data) = day14.counts.matrix$X # assigns column names (previous jsut numbered) as the gene ID 
rownames(d14.data) = names(day14.counts.matrix)[-(1)]; # assigns the row names as the sample ID
d14.data_matrix <- data.frame(day14.counts.matrix[,-1], row.names=day14.counts.matrix[,1]) 
d14.data_matrix_t <- t(d14.data_matrix)

# match number of samples in treatment data - SG92 was not sequenced  - ommit from the treatment data for the dds object 
d14.Treatment_data_OM <- d14.Treatment_data %>% dplyr::filter(!Sample.Name %in% 'SG92')

dds.d14 <- DESeqDataSetFromMatrix(countData = d14.data_matrix,  colData = d14.Treatment_data_OM, design = ~ 1)
dds.d14_vst <- vst(dds.d14) # transform it vst
dds.d14_vst <- assay(dds.d14_vst) # call only the transformed coutns in the dds object
d14_vst <- t(dds.d14_vst) # transpose columns to rows and vice versa
# fix(d14_vst)

# =================================================================================== #
# count data - day  21
d21.data = as.data.frame(t(day21.counts.matrix[, -(1)])) # ommit all columns but samples and transpose
names(d21.data) = day21.counts.matrix$X # assigns column names (previous jsut numbered) as the gene ID 
rownames(d21.data) = names(day21.counts.matrix)[-(1)]; # assigns the row names as the sample ID
d21.data_matrix <- data.frame(day21.counts.matrix[,-1], row.names=day21.counts.matrix[,1]) 
d21.data_matrix_t <- t(d21.data_matrix)
dds.d21 <- DESeqDataSetFromMatrix(countData = d21.data_matrix,  colData = d21.Treatment_data, design = ~ 1)
dds.d21_vst <- vst(dds.d21) # transform it vst
dds.d21_vst <- assay(dds.d21_vst) # call only the transformed coutns in the dds object
d21_vst <- t(dds.d21_vst) # transpose columns to rows and vice versa
# fix(d21_vst)

# =================================================================================== #
# merge treatment with VST count data 

# create common column to merge treatment by Sample.Name 
d7_vst_counts <- cbind(rownames(d7_vst), data.frame(d7_vst, row.names=NULL))
colnames(d7_vst_counts)[1] <- "Sample.Name"
Day7.ExpVST <- merge(d7_vst_counts, d7.Treatment_data, by = 'Sample.Name') # merge 

d14_vst_counts <- cbind(rownames(d14_vst), data.frame(d14_vst, row.names=NULL))
colnames(d14_vst_counts)[1] <- "Sample.Name"
Day14.ExpVST <- merge(d14_vst_counts, d14.Treatment_data_OM, by = 'Sample.Name') # merge 

d21_vst_counts <- cbind(rownames(d21_vst), data.frame(d21_vst, row.names=NULL))
colnames(d21_vst_counts)[1] <- "Sample.Name"
Day21.ExpVST <- merge(d21_vst_counts, d21.Treatment_data, by = 'Sample.Name') # merge 
# fix(Day21.ExpVST)


# use 'any_of' in tidyselect package to get any of the target GOIs present - may have been cut in the filtering process
Day7.ExpVST_GOIs <- Day7.ExpVST %>% dplyr::select('Sample.Name', 'Primary_Treatment', 'Second_Treament', any_of(Target_GOIs$GeneID)) # integrate the target GOIs here!!!!!
Day7.ExpVST_GOIs$group <- paste(Day7.ExpVST_GOIs$Primary_Treatment , Day7.ExpVST_GOIs$Second_Treament , sep='')

Day14.ExpVST_GOIs <- Day14.ExpVST %>% dplyr::select('Sample.Name', 'Primary_Treatment', 'Second_Treament', any_of(Target_GOIs$GeneID)) # integrate the target GOIs here!!!!!
Day14.ExpVST_GOIs$group <- paste(Day14.ExpVST_GOIs$Primary_Treatment , Day14.ExpVST_GOIs$Second_Treament , sep='')

Day21.ExpVST_GOIs <- Day21.ExpVST %>% dplyr::select('Sample.Name', 'Primary_Treatment', 'Second_Treament', 'Third_Treatment', any_of(Target_GOIs$GeneID)) # integrate the target GOIs here!!!!!
Day21.ExpVST_GOIs$group <- paste(Day21.ExpVST_GOIs$Primary_Treatment , Day21.ExpVST_GOIs$Second_Treament , Day21.ExpVST_GOIs$Third_Treatment,  sep='')






# ===================================================================================
# Day 7 data prep for figures
#
# ===================================================================================
# reshape the data frame and merge geneIDs with tragetGOIs to add back the gene titles - use this to name figures downstream
Day7.ExpVST_GOIs_MELT <- melt(Day7.ExpVST_GOIs, id=(c('Sample.Name', 'Primary_Treatment', 'Second_Treament', 'group'))) # melt using reshape2
names(Day7.ExpVST_GOIs_MELT)[(5:6)] <- c('GeneID', 'vst_Expression') # change column names
Day7_ExpVst_Master <- merge(Target_GOIs, Day7.ExpVST_GOIs_MELT, by = 'GeneID')

# calc the m an expressio by gene ID (add gene title in group but this is the same unique level as GeneID - review target_GOIs)
Day7_meanExpr <- Day7_ExpVst_Master %>% 
  dplyr::select(c('Sample.Name','group', 'vst_Expression', 'GeneID', 'Gene_term')) %>% 
  group_by(GeneID, group, Gene_term) %>%
  dplyr::summarize(mean.vstExp = mean(vst_Expression), 
                   sd.vsdtExp = sd(vst_Expression),
                   n = n(), 
                   se.vsdtExp = sd.vsdtExp/sqrt(n))

# create treatment groups 
Day7_meanExpr$PrimaryTreatment <- substr(Day7_meanExpr$group, 1,1) # primary
Day7_meanExpr$SecondTreatment <- substr(Day7_meanExpr$group, 2,2) # second

# ===================================================================================
# Day 14 data prep for figures
#
# ===================================================================================
# reshape the data frame and merge geneIDs with tragetGOIs to add back the gene titles - use this to name figures downstream
Day14.ExpVST_GOIs_MELT <- melt(Day14.ExpVST_GOIs, id=(c('Sample.Name', 'Primary_Treatment', 'Second_Treament', 'group'))) # melt using reshape2
names(Day14.ExpVST_GOIs_MELT)[(5:6)] <- c('GeneID', 'vst_Expression') # change column names
Day14_ExpVst_Master <- merge(Target_GOIs, Day14.ExpVST_GOIs_MELT, by = 'GeneID')

# calc the m an expressio by gene ID (add gene title in group but this is the same unique level as GeneID - review target_GOIs)
Day14_meanExpr <- Day14_ExpVst_Master %>% 
  dplyr::select(c('Sample.Name','group', 'vst_Expression', 'GeneID', 'Gene_term')) %>% 
  group_by(GeneID, group, Gene_term) %>%
  dplyr::summarize(mean.vstExp = mean(vst_Expression), 
                   sd.vsdtExp = sd(vst_Expression),
                   n = n(), 
                   se.vsdtExp = sd.vsdtExp/sqrt(n))

# create treatment groups 
Day14_meanExpr$PrimaryTreatment <- substr(Day14_meanExpr$group, 1,1) # primary
Day14_meanExpr$SecondTreatment <- substr(Day14_meanExpr$group, 2,2) # second

# ===================================================================================
# Day 21 data prep for figures
#
# ===================================================================================
# reshape the data frame and merge geneIDs with tragetGOIs to add back the gene titles - use this to name figures downstream
Day21.ExpVST_GOIs_MELT <- melt(Day21.ExpVST_GOIs, id=(c('Sample.Name', 'Primary_Treatment', 'Second_Treament', 'Third_Treatment','group'))) # melt using reshape2
names(Day21.ExpVST_GOIs_MELT)[(6:7)] <- c('GeneID', 'vst_Expression') # change column names
Day21_ExpVst_Master <- merge(Target_GOIs, Day21.ExpVST_GOIs_MELT, by = 'GeneID')

# calc the m an expressio by gene ID (add gene title in group but this is the same unique level as GeneID - review target_GOIs)
Day21_meanExpr <- Day21_ExpVst_Master %>% 
  dplyr::select(c('Sample.Name','group', 'vst_Expression', 'GeneID', 'Gene_term')) %>% 
  group_by(GeneID, group, Gene_term) %>%
  dplyr::summarize(mean.vstExp = mean(vst_Expression), 
                   sd.vsdtExp = sd(vst_Expression),
                   n = n(), 
                   se.vsdtExp = sd.vsdtExp/sqrt(n))
# fix(Day21_meanExpr)
# create treatment groups 
Day21_meanExpr$PrimaryTreatment <- substr(Day21_meanExpr$group, 1,1) # primary
Day21_meanExpr$SecondTreatment <- substr(Day21_meanExpr$group, 2,2) # second
Day21_meanExpr$ThirdTreatment <- substr(Day21_meanExpr$group, 3,3) # third

# ===================================================================================
# Now for the for looped plots!
#
# ===================================================================================

# LOAD THE GENES OF INTEREST THAT WERE SUPPORTED IN WGCNA AND DESEQ 2 
apriori_DESeq2PrimEffect <- read.csv(file="Analysis/Output/a_priori_hypothesis/ apriori_DESeq2_AllMainEffects.csv", sep=',', header=TRUE)  
apriori_DESeq2_condenced <- unique(apriori_DESeq2PrimEffect[,c(2,10)])


apriori_WGCNA_AllModules <- read.csv(file="Analysis/Output/a_priori_hypothesis/ apriori_WGCNA_AllSigModules.csv", sep=',', header=TRUE)  


# apriori_DESeq2PrimEffect[1,10] # gene term (i.e. superoxide dismutase)
# apriori_DESeq2PrimEffect[1,2]  # the PGEN ID

# for(i in 1:nrow(target_GOIs)) {

# if ( (nrow(Day7_meanExpr %>%  dplyr::filter(genes %in% target_GOIs[i,1]))) > 0 ) {


    
# ===================================================================================
################################################################################ #       
# CALL A GENE AND PLOT IT ...little plotting machine to view target GOIs
################################################################################ #   
# ===================================================================================



# NOTE here are KEGG targets to plot 
KEGG_targets <- dplyr::filter(annot.condenced, grepl('E3 ubiquitin-protein ligase XIAP|E3 ubiquitin-protein ligase TRIP12|
                                                     |threonine-protein kinase PINK1|threonine-protein kinase TBK1|optineurin|
                                                     |sorting nexin-2|sorting nexin-6|clathrin light chain A|
                                                     |peroxisomal acyl-coenzyme A oxidase 1|peroxisomal acyl-coenzyme A oxidase 3|alcohol dehydrogenase class-3', Gene_term, ignore.case = TRUE))
KEGG_targets$Gene_term <- gsub( "/", "", as.character(KEGG_targets$Gene_term)) # removes all occurances of '/' and replaces with a NULL for all gene terms 
KEGG_targets$Gene_term <- gsub( " ", "_", as.character(KEGG_targets$Gene_term)) # removes all occurances of '/' and replaces with a NULL for all gene terms 
KEGG_priori_WGCNA <- merge(Sig_WGCNA.modules,KEGG_targets[c(1,3)], by= "GeneID") # %>% dplyr::select('Sample.Name', 'Primary_Treatment', 'Second_Treament', 'Third_Treatment', any_of(target_GOIs$GeneID))

# Clathrin_light_chain_A_ == PGEN_.00g313920 
# Sorting_nexin-2 == PGEN_.00g198530
# Peroxisomal_acyl-coenzyme_A_oxidase_3 ==  PGEN_.00g135600
# Peroxisomal_acyl-coenzyme_A_oxidase_1 == PGEN_.00g300030
# Sorting_nexin-6 ==  PGEN_.00g024130
# Alcohol_dehydrogenase_class-3 == PGEN_.00g048930
# E3_ubiquitin-protein_ligase_TRIP12 == PGEN_.00g239600
# Serinethreonine-protein_kinase_PINK1,_mitochondria == PGEN_.00g322450
# Serinethreonine-protein_kinase_TBK1 == PGEN_.00g211820

pd <- position_dodge(0.3)
allmeanExr <- rbind(Day7_meanExpr, Day14_meanExpr, Day21_meanExpr)

GOI_name <- 'PGEN_.00g220600' # call your gene you want to see!

ExpMin <- floor(  min((allmeanExr %>% dplyr::filter(GeneID %in% GOI_name))$mean.vstExp) - max((allmeanExr %>% dplyr::filter(GeneID %in% GOI_name))$sd.vsdtExp)  )
ExpMax <- ceiling(  max((allmeanExr %>% dplyr::filter(GeneID %in% GOI_name))$mean.vstExp) + max((allmeanExr %>% dplyr::filter(GeneID %in% GOI_name))$sd.vsdtExp)  )
if ( (nrow(Day7_meanExpr %>%  dplyr::filter(GeneID %in% GOI_name))) > 0 ) {
d7_GOI_plot <- Day7_meanExpr %>% 
  dplyr::filter(GeneID %in% GOI_name) %>% 
  ggplot(aes(x=SecondTreatment, y=mean.vstExp, fill=PrimaryTreatment)) +  # , colour=supp, group=supp))
  theme_classic() +
  geom_errorbar(aes(ymin=mean.vstExp-se.vsdtExp, ymax=mean.vstExp+se.vsdtExp), colour="black", width=.1, position=pd) +
  geom_point(position=pd, size = 4, shape=21) +            
  xlab("Second pCO2 treatment") +
  ylab('Gene Expression (mean±SE VST transformed)') +                 # note the mean was first by sample ID THEN by treatment
  scale_fill_manual(values=c("white","grey50")) +
  # scale_color_manual(values=c("#56B4E9","#E69F00")) +
  # ggtitle(paste("Day 7:",apriori_DESeq2_condenced[i,1],":", apriori_DESeq2_condenced[i,3],sep='')) +
  ggtitle( (Day7_meanExpr %>% dplyr::filter(GeneID %in% GOI_name))$Gene_term[1]) +
  # expand_limits(y=0) +                                                    # Expand y range
  # scale_y_continuous(limits=c((min_p1), (max_p1))) +
  scale_y_continuous(limits = c(ExpMin+0.5,ExpMax-0.5)) +
  # ylim(ExpMin, ExpMin) +
  theme(axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20),
        axis.ticks.length=unit(.25, "cm"))+
  theme(legend.position = "none")
} else { d7_GOI_plot <- plot.new() }


if ( (nrow(Day14_meanExpr %>%  dplyr::filter(GeneID %in% GOI_name))) > 0 ) {
  d14_GOI_plot <- Day14_meanExpr %>% 
    dplyr::filter(GeneID %in% GOI_name) %>% 
    ggplot(aes(x=SecondTreatment, y=mean.vstExp, fill=PrimaryTreatment)) +  # , colour=supp, group=supp))
    theme_classic() +
    geom_errorbar(aes(ymin=mean.vstExp-se.vsdtExp, ymax=mean.vstExp+se.vsdtExp), colour="black", width=.1, position=pd) +
    geom_point(position=pd, size = 4, shape=21) +            
    xlab("Second pCO2 treatment") +
    ylab("") +                 # note the mean was first by sample ID THEN by treatment
    scale_fill_manual(values=c("white","grey50")) +
    # scale_color_manual(values=c("#56B4E9","#E69F00")) +
    # ggtitle(paste("Day 14:",apriori_DESeq2_condenced[i,1],":", apriori_DESeq2_condenced[i,3], sep='')) +
    # ggtitle(paste("Day 7:",apriori_DESeq2_condenced[i,1],":", apriori_DESeq2_condenced[i,3],sep='')) +
    ggtitle(NULL) +
    # expand_limits(y=0) +                                                    # Expand y range
    # scale_y_continuous(limits=c((min_p1), (max_p1))) +
    scale_y_continuous(limits = c(ExpMin+0.5,ExpMax-0.5)) +
    # ylim(ExpMin, ExpMin) +
    theme(axis.text.x = element_text(size = 20),
          axis.text.y = element_text(size = 20),
          axis.ticks.length=unit(.25, "cm")) +
    theme(legend.position = "none")
} else { d14_GOI_plot <- plot.new() }

if ( (nrow(Day21_meanExpr %>%  dplyr::filter(GeneID %in% GOI_name))) > 0 ) {
  d21_GOI_plot <- Day21_meanExpr %>% 
    dplyr::filter(GeneID %in% GOI_name) %>% 
    ggplot(aes(x=ThirdTreatment, y=mean.vstExp, fill=PrimaryTreatment)) +  # , colour=supp, group=supp))
    theme_classic() +
    geom_errorbar(aes(ymin=mean.vstExp-se.vsdtExp, ymax=mean.vstExp+se.vsdtExp), colour="black", width=.1, position=pd) +
    geom_point(position=pd, size = 4, shape=21) +            
    xlab("Third pCO2 treatment") +
    ylab('') +                 # note the mean was first by sample ID THEN by treatment
    scale_fill_manual(values=c("white","grey50")) +
    # scale_color_manual(values=c("#56B4E9","#E69F00")) +
    # ggtitle(paste("Day 21:",apriori_DESeq2_condenced[i,1], ":", apriori_DESeq2_condenced[i,3], sep='')) +
    # ggtitle(paste("Day 7:",apriori_DESeq2_condenced[i,1],":", apriori_DESeq2_condenced[i,3],sep='')) +
    ggtitle(NULL) +
    # expand_limits(y=0) +                                                    # Expand y range
    # scale_y_continuous(limits=c((min_p1), (max_p1))) +
    scale_y_continuous(limits = c(ExpMin+0.5,ExpMax-0.5)) +
    # ylim(ExpMin, ExpMin) +
    theme(axis.text.x = element_text(size = 20),
          axis.text.y = element_text(size = 20),
          axis.ticks.length=unit(.25, "cm")) +
    theme(legend.position = "none") +
    facet_wrap(~SecondTreatment)
} else { d21_GOI_plot <- plot.new() }

pdf(paste("Analysis/Output/a_priori_hypothesis/other/",GOI_name,"_Optineurin.pdf", sep =''), width=15, height=5)
print(ggarrange(d7_GOI_plot, d14_GOI_plot, d21_GOI_plot,        
                plotlist = NULL,
                ncol = 3,
                nrow = 1,
                labels = NULL))
dev.off()




################################################################################ #       
# For loops plotting DESeq and WGCNA assoaited taret GOIs
#
################################################################################ #    

# set-up for the for loops
pd <- position_dodge(0.3)
allmeanExr <- rbind(Day7_meanExpr, Day14_meanExpr, Day21_meanExpr)

View(apriori_DESeq2_condenced)
# DESeq2  - DEGs overlapped with target GOIs - output all figures 
# run the for loop for plotting by each tagert GOI with WGCNA and DESeq2 results
for(i in 1:nrow(apriori_DESeq2_condenced)) {
  ExpMin <- floor(  min((allmeanExr %>% dplyr::filter(GeneID %in% apriori_DESeq2_condenced[i,1]))$mean.vstExp) - max((allmeanExr %>% dplyr::filter(GeneID %in% apriori_DESeq2_condenced[i,1]))$sd.vsdtExp)  )
  ExpMax <- ceiling(  max((allmeanExr %>% dplyr::filter(GeneID %in% apriori_DESeq2_condenced[i,1]))$mean.vstExp) + max((allmeanExr %>% dplyr::filter(GeneID %in% apriori_DESeq2_condenced[i,1]))$sd.vsdtExp)  )
  if ( (nrow(Day7_meanExpr %>%  dplyr::filter(GeneID %in% apriori_DESeq2_condenced[i,1]))) > 0 ) {
        d7_plots <- Day7_meanExpr %>% 
                      dplyr::filter(GeneID %in% apriori_DESeq2_condenced[i,1]) %>% 
                      ggplot(aes(x=SecondTreatment, y=mean.vstExp, fill=PrimaryTreatment)) +  # , colour=supp, group=supp))
                      theme_classic() +
                      geom_errorbar(aes(ymin=mean.vstExp-se.vsdtExp, ymax=mean.vstExp+se.vsdtExp), colour="black", width=.1, position=pd) +
                      geom_point(position=pd, size = 4, shape=21) +            
                      xlab("Second pCO2 treatment") +
                      ylab('Gene Expression (mean±SE VST transformed)') +                 # note the mean was first by sample ID THEN by treatment
                      scale_fill_manual(values=c("white","grey50")) +
                      # scale_color_manual(values=c("#56B4E9","#E69F00")) +
                      # ggtitle(paste("Day 7:",apriori_DESeq2_condenced[i,1],":", apriori_DESeq2_condenced[i,3],sep='')) +
                      ggtitle(apriori_DESeq2_condenced[i,2]) +
                      # expand_limits(y=0) +                                                    # Expand y range
                      # scale_y_continuous(limits=c((min_p1), (max_p1))) +
                      scale_y_continuous(limits = c(ExpMin+0.5,ExpMax-0.5)) +
                      # ylim(ExpMin, ExpMin) +
                     theme(axis.text.x = element_text(size = 20),
                           axis.text.y = element_text(size = 20),
                           axis.ticks.length=unit(.25, "cm"))+
                      theme(legend.position = "none")
        } else { d7_plots <- plot.new() }
  

if ( (nrow(Day14_meanExpr %>%  dplyr::filter(GeneID %in% apriori_DESeq2_condenced[i,1]))) > 0 ) {
        d14_plots <- Day14_meanExpr %>% 
                      dplyr::filter(GeneID %in% apriori_DESeq2_condenced[i,1]) %>% 
                      ggplot(aes(x=SecondTreatment, y=mean.vstExp, fill=PrimaryTreatment)) +  # , colour=supp, group=supp))
                      theme_classic() +
                      geom_errorbar(aes(ymin=mean.vstExp-se.vsdtExp, ymax=mean.vstExp+se.vsdtExp), colour="black", width=.1, position=pd) +
                      geom_point(position=pd, size = 4, shape=21) +            
                      xlab("Second pCO2 treatment") +
                      ylab("") +                 # note the mean was first by sample ID THEN by treatment
                      scale_fill_manual(values=c("white","grey50")) +
                      # scale_color_manual(values=c("#56B4E9","#E69F00")) +
                      # ggtitle(paste("Day 14:",apriori_DESeq2_condenced[i,1],":", apriori_DESeq2_condenced[i,3], sep='')) +
                      # ggtitle(paste("Day 7:",apriori_DESeq2_condenced[i,1],":", apriori_DESeq2_condenced[i,3],sep='')) +
                      ggtitle(NULL) +
                      # expand_limits(y=0) +                                                    # Expand y range
                      # scale_y_continuous(limits=c((min_p1), (max_p1))) +
                      scale_y_continuous(limits = c(ExpMin+0.5,ExpMax-0.5)) +
                      # ylim(ExpMin, ExpMin) +
                      theme(axis.text.x = element_text(size = 20),
                            axis.text.y = element_text(size = 20),
                            axis.ticks.length=unit(.25, "cm")) +
                      theme(legend.position = "none")
  } else { d14_plots <- plot.new() }
  
  

if ( (nrow(Day21_meanExpr %>%  dplyr::filter(GeneID %in% apriori_DESeq2_condenced[i,1]))) > 0 ) {
          d21_plots <- Day21_meanExpr %>% 
                        dplyr::filter(GeneID %in% apriori_DESeq2_condenced[i,1]) %>% 
                        ggplot(aes(x=ThirdTreatment, y=mean.vstExp, fill=PrimaryTreatment)) +  # , colour=supp, group=supp))
                        theme_classic() +
                        geom_errorbar(aes(ymin=mean.vstExp-se.vsdtExp, ymax=mean.vstExp+se.vsdtExp), colour="black", width=.1, position=pd) +
                        geom_point(position=pd, size = 4, shape=21) +            
                        xlab("Third pCO2 treatment") +
                        ylab('') +                 # note the mean was first by sample ID THEN by treatment
                        scale_fill_manual(values=c("white","grey50")) +
                        # scale_color_manual(values=c("#56B4E9","#E69F00")) +
                        # ggtitle(paste("Day 21:",apriori_DESeq2_condenced[i,1], ":", apriori_DESeq2_condenced[i,3], sep='')) +
                        # ggtitle(paste("Day 7:",apriori_DESeq2_condenced[i,1],":", apriori_DESeq2_condenced[i,3],sep='')) +
                        ggtitle(NULL) +
                        # expand_limits(y=0) +                                                    # Expand y range
                        # scale_y_continuous(limits=c((min_p1), (max_p1))) +
                        scale_y_continuous(limits = c(ExpMin+0.5,ExpMax-0.5)) +
                        # ylim(ExpMin, ExpMin) +
                        theme(axis.text.x = element_text(size = 20),
                              axis.text.y = element_text(size = 20),
                              axis.ticks.length=unit(.25, "cm")) +
                        theme(legend.position = "none") +
                        facet_wrap(~SecondTreatment)
} else { d21_plots <- plot.new() }


pdf(paste("Analysis/Output/a_priori_hypothesis/DESeq2_results/",apriori_DESeq2_condenced[i,2],apriori_DESeq2_condenced[i,1],".pdf", sep =''), width=15, height=5)
print(ggarrange(d7_plots, d14_plots, d21_plots,        
                plotlist = NULL,
                ncol = 3,
                nrow = 1,
                labels = NULL))
dev.off()
}






# WGCNA   - modules overlapped with target GOIs - output all figures 
# run the for loop for plotting by each tagert GOI with WGCNA  results
apriori_WGCNA_AllModules_condenced <- apriori_WGCNA_AllModules[,c(2:7)]
apriori_WGCNA_PrimEffModules <- apriori_WGCNA_AllModules %>% dplyr::filter(Day_color %in% c('Day7_brown', 'Day7_yellow', 
                                                                                           'Day14_brown', 'Day14_black', 
                                                                                           'Day21_blue',  'Day21_magenta', 'Day21_yellow'))
apriori_WGCNA_PrimEffModules_condenced <- apriori_WGCNA_PrimEffModules[,c(2:4,6,7)]

for(i in 1:nrow(apriori_WGCNA_PrimEffModules_condenced)) {
  ExpMin <- floor(  min((allmeanExr %>% dplyr::filter(GeneID %in% apriori_WGCNA_PrimEffModules_condenced[i,1]))$mean.vstExp) - max((allmeanExr %>% dplyr::filter(GeneID %in% apriori_WGCNA_PrimEffModules_condenced[i,1]))$sd.vsdtExp)  )
  ExpMax <- ceiling(  max((allmeanExr %>% dplyr::filter(GeneID %in% apriori_WGCNA_PrimEffModules_condenced[i,1]))$mean.vstExp) + max((allmeanExr %>% dplyr::filter(GeneID %in% apriori_WGCNA_PrimEffModules_condenced[i,1]))$sd.vsdtExp)  )
  if ( (nrow(Day7_meanExpr %>%  dplyr::filter(GeneID %in% apriori_WGCNA_PrimEffModules_condenced[i,1]))) > 0 ) {
    d7_plots <- Day7_meanExpr %>% 
      dplyr::filter(GeneID %in% apriori_WGCNA_PrimEffModules_condenced[i,1]) %>% 
      ggplot(aes(x=SecondTreatment, y=mean.vstExp, fill=PrimaryTreatment)) +  # , colour=supp, group=supp))
      theme_classic() +
      geom_errorbar(aes(ymin=mean.vstExp-se.vsdtExp, ymax=mean.vstExp+se.vsdtExp), colour="black", width=.1, position=pd) +
      geom_point(position=pd, size = 4, shape=21) +            
      xlab("Second pCO2 treatment") +
      ylab('Gene Expression (mean±SE VST transformed)') +                 # note the mean was first by sample ID THEN by treatment
      scale_fill_manual(values=c("white","grey50")) +
      # scale_color_manual(values=c("#56B4E9","#E69F00")) +
      # ggtitle(paste("Day 7:",apriori_DESeq2_condenced[i,1],":", apriori_DESeq2_condenced[i,3],sep='')) +
      ggtitle(paste(apriori_WGCNA_PrimEffModules_condenced[i,3],apriori_WGCNA_PrimEffModules_condenced[i,2], ":", apriori_WGCNA_PrimEffModules_condenced[i,5], sep =' ')) +
      # expand_limits(y=0) +                                                    # Expand y range
      # scale_y_continuous(limits=c((min_p1), (max_p1))) +
      scale_y_continuous(limits = c(ExpMin+0.5,ExpMax-0.5)) +
      # ylim(ExpMin, ExpMin) +
      theme(axis.text.x = element_text(size = 20),
            axis.text.y = element_text(size = 20),
            axis.ticks.length=unit(.25, "cm"))+
      theme(legend.position = "none")
  } else { d7_plots <- plot.new() }
  
  
  if ( (nrow(Day14_meanExpr %>%  dplyr::filter(GeneID %in% apriori_WGCNA_PrimEffModules_condenced[i,1]))) > 0 ) {
    d14_plots <- Day14_meanExpr %>% 
      dplyr::filter(GeneID %in% apriori_WGCNA_PrimEffModules_condenced[i,1]) %>% 
      ggplot(aes(x=SecondTreatment, y=mean.vstExp, fill=PrimaryTreatment)) +  # , colour=supp, group=supp))
      theme_classic() +
      geom_errorbar(aes(ymin=mean.vstExp-se.vsdtExp, ymax=mean.vstExp+se.vsdtExp), colour="black", width=.1, position=pd) +
      geom_point(position=pd, size = 4, shape=21) +            
      xlab("Second pCO2 treatment") +
      ylab("") +                 # note the mean was first by sample ID THEN by treatment
      scale_fill_manual(values=c("white","grey50")) +
      # scale_color_manual(values=c("#56B4E9","#E69F00")) +
      # ggtitle(paste("Day 14:",apriori_DESeq2_condenced[i,1],":", apriori_DESeq2_condenced[i,3], sep='')) +
      # ggtitle(paste("Day 7:",apriori_DESeq2_condenced[i,1],":", apriori_DESeq2_condenced[i,3],sep='')) +
      ggtitle(NULL) +
      # expand_limits(y=0) +                                                    # Expand y range
      # scale_y_continuous(limits=c((min_p1), (max_p1))) +
      scale_y_continuous(limits = c(ExpMin+0.5,ExpMax-0.5)) +
      # ylim(ExpMin, ExpMin) +
      theme(axis.text.x = element_text(size = 20),
            axis.text.y = element_text(size = 20),
            axis.ticks.length=unit(.25, "cm")) +
      theme(legend.position = "none")
  } else { d14_plots <- plot.new() }
  
  
  
  if ( (nrow(Day21_meanExpr %>%  dplyr::filter(GeneID %in% apriori_WGCNA_PrimEffModules_condenced[i,1]))) > 0 ) {
    d21_plots <- Day21_meanExpr %>% 
      dplyr::filter(GeneID %in% apriori_WGCNA_PrimEffModules_condenced[i,1]) %>% 
      ggplot(aes(x=ThirdTreatment, y=mean.vstExp, fill=PrimaryTreatment)) +  # , colour=supp, group=supp))
      theme_classic() +
      geom_errorbar(aes(ymin=mean.vstExp-se.vsdtExp, ymax=mean.vstExp+se.vsdtExp), colour="black", width=.1, position=pd) +
      geom_point(position=pd, size = 4, shape=21) +            
      xlab("Third pCO2 treatment") +
      ylab('') +                 # note the mean was first by sample ID THEN by treatment
      scale_fill_manual(values=c("white","grey50")) +
      # scale_color_manual(values=c("#56B4E9","#E69F00")) +
      # ggtitle(paste("Day 21:",apriori_DESeq2_condenced[i,1], ":", apriori_DESeq2_condenced[i,3], sep='')) +
      # ggtitle(paste("Day 7:",apriori_DESeq2_condenced[i,1],":", apriori_DESeq2_condenced[i,3],sep='')) +
      ggtitle(NULL) +
      # expand_limits(y=0) +                                                    # Expand y range
      # scale_y_continuous(limits=c((min_p1), (max_p1))) +
      scale_y_continuous(limits = c(ExpMin+0.5,ExpMax-0.5)) +
      # ylim(ExpMin, ExpMin) +
      theme(axis.text.x = element_text(size = 20),
            axis.text.y = element_text(size = 20),
            axis.ticks.length=unit(.25, "cm")) +
      theme(legend.position = "none") +
      facet_wrap(~SecondTreatment)
  } else { d21_plots <- plot.new() }
  
  
  pdf(paste("Analysis/Output/a_priori_hypothesis/WGCNA_results/",apriori_WGCNA_PrimEffModules_condenced[i,3],"_", apriori_WGCNA_PrimEffModules_condenced[i,2], "_", apriori_WGCNA_PrimEffModules_condenced[i,1],".pdf", sep =''), width=15, height=5)
  print(ggarrange(d7_plots, d14_plots, d21_plots,        
                  plotlist = NULL,
                  ncol = 3,
                  nrow = 1,
                  labels = NULL))
  dev.off()
}


