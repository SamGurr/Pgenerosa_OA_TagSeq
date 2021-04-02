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
# SET WORKING DIRECTORY AND LOAD DATA
setwd("C:/Users/samjg/Documents/My_Projects/Pgenerosa_TagSeq_Metabolomics/TagSeq/")

# LOAD DATA
# count data
day7.counts.matrix <- read.csv(file="Analysis/Data/Filtered_Counts/5cpm_50perc/day7.counts.filtered_5cpm50perc.csv", sep=',', header=TRUE)
day14.counts.matrix <- read.csv(file="Analysis/Data/Filtered_Counts/5cpm_50perc/day14.counts.filtered_5cpm50perc.csv", sep=',', header=TRUE)
day21.counts.matrix <- read.csv(file="Analysis/Data/Filtered_Counts/5cpm_50perc/day21.counts.filtered_5cpm50perc.csv", sep=',', header=TRUE)
# trait data
Master.Treatment_Phenotype.data <- read.csv(file="Analysis/Data/Experiment_Metadata/Master_Phyenotype.and.Exp.Treatment_Metadata.csv", sep=',', header=TRUE)
# WGCNA results 
D7_WGCNA<- read.csv(file="Analysis/Output/WGCNA/Day7/d7.WGCNA_ModulMembership.csv", sep=',', header=TRUE)
D14_WGCNA<- read.csv(file="Analysis/Output/WGCNA/Day14/d14.WGCNA_ModulMembership.csv", sep=',', header=TRUE)
D21_WGCNA<- read.csv(file="Analysis/Output/WGCNA/Day21/d21.WGCNA_ModulMembership.csv", sep=',', header=TRUE)
D7_WGCNA$Day  <- "Day7"
D14_WGCNA$Day <- "Day14"
D21_WGCNA$Day <- "Day21"
# DESEq2 results 
D7_DEGs_PrimEffect  <- read.csv(file="Analysis/Output/DESeq2/10cpm/Day7/DE_Day7_Primary.csv", sep=',', header=TRUE)  %>% dplyr::select(c('Row.names', 'log2FoldChange', 'pvalue', 'up', 'down'))
D14_DEGs_PrimEffect <- read.csv(file="Analysis/Output/DESeq2/10cpm/Day14/DE_Day14_Primary.csv", sep=',', header=TRUE) %>% dplyr::select(c('Row.names', 'log2FoldChange', 'pvalue', 'up', 'down'))
D21_DEGs_PrimEffect <- read.csv(file="Analysis/Output/DESeq2/10cpm/Day21/DE_Day21_Primary.csv", sep=',', header=TRUE) %>% dplyr::select(c('Row.names', 'log2FoldChange', 'pvalue', 'up', 'down'))
D7_DEGs_PrimEffect$Day  <- "Day7"
D14_DEGs_PrimEffect$Day <- "Day14"
D21_DEGs_PrimEffect$Day <- "Day21"
DEGs_PrimEffect_10cpm <- rbind(D7_DEGs_PrimEffect, D14_DEGs_PrimEffect, D21_DEGs_PrimEffect) # bind dataframes 
names(DEGs_PrimEffect_10cpm)[1] <- "gene.ID" # change name of column 1


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
# PGEN_.00g283000 - DNMT_1
# PGEN_.00g029420 - DNMT3A
# PGEN_.00g067800 - DNMT 3b
# PGEN_.00g053100 + PGEN_.00g053120 - histone methyl transferase
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
# PGEN_.00g192250,PGEN_.00g180320, PGEN_.00g287800, PGEN_.00g293960- glutathione peroxidase

genes <- c('AOX','NADH_dehydrogenase','Cytochrome_c_reductase','Uncoupling_protein_1','Uncoupling_protein_2','Uncoupling_protein_3',
           ' S-adenosylmethionine_synthase_isoform_type-2','DNMT_1','DNMT3A','DNMT_3b','HMT_1','HMT_2','positiv_hist_methyl','positiv_hist_methyl_2','methyltransferase','HAT_1','HAT_2','HAT_3','HAT_4','HAT_5',
           'SIR1','SIR2','SIR4','SIR5','SIR6','SIR7',
           'SOD_1','SOD_2','SOD_3','SOD_Cu_chaperone','SOD_Mn_act','glutathione_peroxidase_1','glutathione_peroxidase_2','glutathione_peroxidase_3','glutathione_peroxidase_4', 'glutathione_peroxidase_5','glutathione_peroxidase_6',
           'titin', 'calpain')


GeneID <- c('PGEN_.00g108770', 'PGEN_.00g299160', 'PGEN_00g275780','PGEN_.00g063670', 'PGEN_.00g193030', 'PGEN_.00g230260', # MITCHONDRIAL PLAYERS 
            'PGEN_.00g041430'  , 'PGEN_.00g283000','PGEN_.00g029420','PGEN_.00g067800','PGEN_.00g053100','PGEN_.00g053120','PGEN_.00g064910', 'PGEN_.00g064920','PGEN_.00g066460', 'PGEN_.00g283340',  'PGEN_.00g311570', 'PGEN_.00g320700', 'PGEN_.00g338440','PGEN_.00g272910', # TRANASCRIPTIONAL REGULATION - proteins involved in methylation and histone modification(s)
            'PGEN_.00g048200', 'PGEN_.00g012340', 'PGEN_.00g149480', 'PGEN_.00g153700', 'PGEN_.00g144540', 'PGEN_.00g033970',# SIRTUINS 
            'PGEN_.00g010160', 'PGEN_.00g065700', 'PGEN_.00g257600', 'PGEN_.00g015070', 'PGEN_.00g062450', 'PGEN_.00g293960', 'PGEN_.00g287800', 'PGEN_.00g192250','PGEN_.00g180320', 'PGEN_.00g116940','PGEN_.00g049360',# OXIDATIVE STRESS
            'PGEN_.00g066340',  'PGEN_.00g014370')
target_GOIs <- data.frame(genes, GeneID)


# ===================================================================================
# 
#  Loop to make table - look in DESeq2 results and WGCNa results! 
#
# ===================================================================================



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
Day7.ExpVST_GOIs <- Day7.ExpVST %>% dplyr::select('Sample.Name', 'Primary_Treatment', 'Second_Treament', any_of(target_GOIs$GeneID))
Day7.ExpVST_GOIs$group <- paste(Day7.ExpVST_GOIs$Primary_Treatment , Day7.ExpVST_GOIs$Second_Treament , sep='')

Day14.ExpVST_GOIs <- Day14.ExpVST %>% dplyr::select('Sample.Name', 'Primary_Treatment', 'Second_Treament', any_of(target_GOIs$GeneID))
Day14.ExpVST_GOIs$group <- paste(Day14.ExpVST_GOIs$Primary_Treatment , Day14.ExpVST_GOIs$Second_Treament , sep='')

Day21.ExpVST_GOIs <- Day21.ExpVST %>% dplyr::select('Sample.Name', 'Primary_Treatment', 'Second_Treament', 'Third_Treatment', any_of(target_GOIs$GeneID))
Day21.ExpVST_GOIs$group <- paste(Day21.ExpVST_GOIs$Primary_Treatment , Day21.ExpVST_GOIs$Second_Treament , Day21.ExpVST_GOIs$Third_Treatment,  sep='')






# ===================================================================================
# Day 7 data prep for figures
#
# ===================================================================================
# reshape the data frame and merge geneIDs with tragetGOIs to add back the gene titles - use this to name figures downstream
Day7.ExpVST_GOIs_MELT <- melt(Day7.ExpVST_GOIs, id=(c('Sample.Name', 'Primary_Treatment', 'Second_Treament', 'group'))) # melt using reshape2
names(Day7.ExpVST_GOIs_MELT)[(5:6)] <- c('GeneID', 'vst_Expression') # change column names
Day7_ExpVst_Master <- merge(target_GOIs, Day7.ExpVST_GOIs_MELT, by = 'GeneID')

# calc the m an expressio by gene ID (add gene title in group but this is the same unique level as GeneID - review target_GOIs)
Day7_meanExpr <- Day7_ExpVst_Master %>% 
  dplyr::select(c('Sample.Name','group', 'vst_Expression', 'GeneID', 'genes')) %>% 
  group_by(GeneID, group, genes) %>%
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
Day14_ExpVst_Master <- merge(target_GOIs, Day14.ExpVST_GOIs_MELT, by = 'GeneID')

# calc the m an expressio by gene ID (add gene title in group but this is the same unique level as GeneID - review target_GOIs)
Day14_meanExpr <- Day14_ExpVst_Master %>% 
  dplyr::select(c('Sample.Name','group', 'vst_Expression', 'GeneID', 'genes')) %>% 
  group_by(GeneID, group, genes) %>%
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
Day21_ExpVst_Master <- merge(target_GOIs, Day21.ExpVST_GOIs_MELT, by = 'GeneID')

# calc the m an expressio by gene ID (add gene title in group but this is the same unique level as GeneID - review target_GOIs)
Day21_meanExpr <- Day21_ExpVst_Master %>% 
  dplyr::select(c('Sample.Name','group', 'vst_Expression', 'GeneID', 'genes')) %>% 
  group_by(GeneID, group, genes) %>%
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
pd <- position_dodge(0.3)
for(i in 1:nrow(target_GOIs)) {

if ( (nrow(Day7_meanExpr %>%  dplyr::filter(genes %in% target_GOIs[i,1]))) > 0 ) {
        d7_plots <- Day7_meanExpr %>% 
                      dplyr::filter(genes %in% target_GOIs[i,1]) %>% 
                      ggplot(aes(x=SecondTreatment, y=mean.vstExp, fill=PrimaryTreatment)) +  # , colour=supp, group=supp))
                      theme_classic() +
                      geom_errorbar(aes(ymin=mean.vstExp-se.vsdtExp, ymax=mean.vstExp+se.vsdtExp), colour="black", width=.1, position=pd) +
                      geom_point(position=pd, size = 4, shape=21) +            
                      xlab("Second pCO2 treatment") +
                      ylab('') +                 # note the mean was first by sample ID THEN by treatment
                      scale_fill_manual(values=c("#56B4E9","#E69F00")) +
                      # scale_color_manual(values=c("#56B4E9","#E69F00")) +
                      ggtitle(paste("Day 7:",target_GOIs[i,1], sep='')) +
                      # expand_limits(y=0) +                                                    # Expand y range
                      #scale_y_continuous(limits=c((min_p1), (max_p1))) +
                      theme(text = element_text(size=15)) +
                      theme(legend.position = "none")
        } else { d7_plots <- plot.new() }
  

if ( (nrow(Day14_meanExpr %>%  dplyr::filter(genes %in% target_GOIs[i,1]))) > 0 ) {
        d14_plots <- Day14_meanExpr %>% 
                      dplyr::filter(genes %in% target_GOIs[i,1]) %>% 
                      ggplot(aes(x=SecondTreatment, y=mean.vstExp, fill=PrimaryTreatment)) +  # , colour=supp, group=supp))
                      theme_classic() +
                      geom_errorbar(aes(ymin=mean.vstExp-se.vsdtExp, ymax=mean.vstExp+se.vsdtExp), colour="black", width=.1, position=pd) +
                      geom_point(position=pd, size = 4, shape=21) +            
                      xlab("Second pCO2 treatment") +
                      ylab("Gene Expression (mean±SE VST transformed)") +                 # note the mean was first by sample ID THEN by treatment
                      scale_fill_manual(values=c("#56B4E9","#E69F00")) +
                      # scale_color_manual(values=c("#56B4E9","#E69F00")) +
                      ggtitle(paste("Day 14:",target_GOIs[i,1], sep='')) +
                      # expand_limits(y=0) +                                                    # Expand y range
                      #scale_y_continuous(limits=c((min_p1), (max_p1))) +
                      theme(text = element_text(size=15)) +
                      theme(legend.position = "none")
  } else { d14_plots <- plot.new() }
  
  

if ( (nrow(Day21_meanExpr %>%  dplyr::filter(genes %in% target_GOIs[i,1]))) > 0 ) {
          d21_plots <- Day21_meanExpr %>% 
                        dplyr::filter(genes %in% target_GOIs[i,1]) %>% 
                        ggplot(aes(x=ThirdTreatment, y=mean.vstExp, fill=PrimaryTreatment)) +  # , colour=supp, group=supp))
                        theme_classic() +
                        geom_errorbar(aes(ymin=mean.vstExp-se.vsdtExp, ymax=mean.vstExp+se.vsdtExp), colour="black", width=.1, position=pd) +
                        geom_point(position=pd, size = 4, shape=21) +            
                        xlab("Third pCO2 treatment") +
                        ylab('') +                 # note the mean was first by sample ID THEN by treatment
                        scale_fill_manual(values=c("#56B4E9","#E69F00")) +
                        # scale_color_manual(values=c("#56B4E9","#E69F00")) +
                        ggtitle(paste("Day 21:",target_GOIs[i,1], sep='')) +
                        # expand_limits(y=0) +                                                    # Expand y range
                        #scale_y_continuous(limits=c((min_p1), (max_p1))) +
                        theme(text = element_text(size=15)) +
                        theme(legend.position = "none") +
                        facet_wrap(~SecondTreatment)
} else { d21_plots <- plot.new() }

getwd()
pdf(paste("Analysis/Output/a_priori_hypothesis/",target_GOIs[i,1],".pdf"), width=10, height=6)
print(ggarrange(d7_plots, d14_plots, d21_plots,        
                plotlist = NULL,
                ncol = 3,
                nrow = 1,
                labels = NULL))
dev.off()
}

