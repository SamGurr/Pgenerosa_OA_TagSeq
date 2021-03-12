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

# SET WORKING DIRECTORY AND LOAD DATA
setwd("C:/Users/samjg/Documents/My_Projects/Pgenerosa_TagSeq_Metabolomics/TagSeq/")

# LOAD DATA
# count data
day7.counts.matrix <- read.csv(file="Analysis/Data/filtered_counts/day7.counts.filtered_5cpm50perc.csv", sep=',', header=TRUE)
day14.counts.matrix <- read.csv(file="Analysis/Data/filtered_counts/day14.counts.filtered_5cpm50perc.csv", sep=',', header=TRUE)
day21.counts.matrix <- read.csv(file="Analysis/Data/filtered_counts/day21.counts.filtered_5cpm50perc.csv", sep=',', header=TRUE)
# trait data
Master.Treatment_Phenotype.data <- read.csv(file="Analysis/Data/ Master_Phyenotype.and.Exp.Treatment_Metadata.csv", sep=',', header=TRUE)

# ===================================================================================
#
# DATA PREP
#
# ===================================================================================

# ===================================================================================
# treatment data
d7.Treatment_data <- Master.Treatment_Phenotype.data %>%   dplyr::filter(Date %in% 20190731) %>%  dplyr::select('Sample.Name', 'Primary_Treatment', 'Second_Treament') # split for day 7 data 
d14.Treatment_data <- Master.Treatment_Phenotype.data  %>%  dplyr::filter(Date %in% 20190807) %>%  dplyr::select('Sample.Name', 'Primary_Treatment', 'Second_Treament')# split for day 7 data 
d21.Treatment_data <- Master.Treatment_Phenotype.data  %>%  dplyr::filter(Date %in% 20190814) %>%  dplyr::select('Sample.Name', 'Primary_Treatment', 'Second_Treament', 'Third_Treatment')# split for day 7 data 

# ==================================================================================
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
fix(d7_vst)

# ===================================================================================
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
fix(d14_vst)

# ===================================================================================
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
fix(d21_vst)

# ===================================================================================
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


# ===================================================================================
#
# a prior hypothesis
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


target_GOIs <- c('PGEN_.00g108770', 'PGEN_.00g299160', 'PGEN_00g275780','PGEN_.00g063670', 'PGEN_.00g193030', 'PGEN_.00g230260', # MITCHONDRIAL PLAYERS 
                 'PGEN_.00g283000','PGEN_.00g029420','PGEN_.00g067800','PGEN_.00g053100','PGEN_.00g064910','PGEN_.00g066460', 'PGEN_.00g283340',  'PGEN_.00g311570', 'PGEN_.00g320700', 'PGEN_.00g338440','PGEN_.00g272910', # TRANASCRIPTIONAL REGULATION - proteins involved in methylation and histone modification(s)
                 'PGEN_.00g048200', 'PGEN_.00g012340', 'PGEN_.00g149480', 'PGEN_.00g153700', 'PGEN_.00g144540', 'PGEN_.00g033970',# SIRTUINS 
                'PGEN_.00g010160', 'PGEN_.00g065700', 'PGEN_.00g257600', 'PGEN_.00g015070', 'PGEN_.00g062450', 'PGEN_.00g293960', 'PGEN_.00g287800', 'PGEN_.00g192250','PGEN_.00g180320')# OXIDATIVE STRESS

# use 'any_of' in tidyselect package to get any of the target GOIs present - may have been cut in the filtering process
Day7.ExpVST_GOIs <- Day7.ExpVST %>% dplyr::select('Sample.Name', 'Primary_Treatment', 'Second_Treament', any_of(target_GOIs))
Day7.ExpVST_GOIs$group <- paste(Day7.ExpVST_GOIs$Primary_Treatment , Day7.ExpVST_GOIs$Second_Treament , sep='')

Day14.ExpVST_GOIs <- Day14.ExpVST %>% dplyr::select('Sample.Name', 'Primary_Treatment', 'Second_Treament', any_of(target_GOIs))
Day14.ExpVST_GOIs$group <- paste(Day14.ExpVST_GOIs$Primary_Treatment , Day14.ExpVST_GOIs$Second_Treament , sep='')

Day21.ExpVST_GOIs <- Day21.ExpVST %>% dplyr::select('Sample.Name', 'Primary_Treatment', 'Second_Treament', 'Third_Treatment', any_of(target_GOIs))
Day21.ExpVST_GOIs$group <- paste(Day21.ExpVST_GOIs$Primary_Treatment , Day21.ExpVST_GOIs$Second_Treament , Day21.ExpVST_GOIs$Third_Treatment,  sep='')



# ===================================================================================
#
# PLOTTING
#
# ===================================================================================

# ===================================================================================
# Day 7
# ===================================================================================
Day7.ExpVST_GOIs_MELT <- melt(Day7.ExpVST_GOIs, id=(c('Sample.Name', 'Primary_Treatment', 'Second_Treament', 'group'))) # melt using reshape2
names(Day7.ExpVST_GOIs_MELT)[(5:6)] <- c('GeneID', 'vst_Expression') # change column names

Day7_meanExpr <- Day7.ExpVST_GOIs_MELT %>% 
  select(c('Sample.Name','group', 'vst_Expression', 'GeneID')) %>% 
  group_by(GeneID, group) %>%
  dplyr::summarize(mean.vstExp = mean(vst_Expression), 
                   sd.vsdtExp = sd(vst_Expression),
                   n = n(), 
                   se.vsdtExp = sd.vsdtExp/sqrt(n))

pd <- position_dodge(0.3) 
Day7_meanExpr %>% 
  dplyr::filter(GeneID %in% 'PGEN_.00g108770') %>% 
  ggplot(aes(x=group, y=mean.vstExp, fill=group)) +  # , colour=supp, group=supp))
  theme_classic() +
  geom_errorbar(aes(ymin=mean.vstExp-se.vsdtExp, ymax=mean.vstExp+se.vsdtExp), colour="black", width=.1, position=pd) +
  geom_line(position=pd) +
  geom_point(position=pd, size = 4, shape=21) +            
  xlab("PrimaryxSecond pCO2 treatment") +
  ylab("Gene Expression (mean±SE VST transformed)") +                 # note the mean was first by sample ID THEN by treatment
  scale_fill_manual(values=c("#56B4E9","#56B4E9","#56B4E9", "#E69F00","#E69F00","#E69F00")) +
  # scale_color_manual(values=c("#56B4E9","#E69F00")) +
  ggtitle("Day 7 vstExpression: Alternative Oxidase ") +
  # expand_limits(y=0) +                                                    # Expand y range
  #scale_y_continuous(limits=c((min_p1), (max_p1))) +
  theme(text = element_text(size=15))

Day7_meanExpr %>% 
  dplyr::filter(GeneID %in% 'PGEN_.00g283000') %>% 
  ggplot(aes(x=group, y=mean.vstExp, fill=group)) +  # , colour=supp, group=supp))
  theme_classic() +
  geom_errorbar(aes(ymin=mean.vstExp-se.vsdtExp, ymax=mean.vstExp+se.vsdtExp), colour="black", width=.1, position=pd) +
  geom_line(position=pd) +
  geom_point(position=pd, size = 4, shape=21) +            
  xlab("PrimaryxSecond pCO2 treatment") +
  ylab("Gene Expression (mean±SE VST transformed)") +                 # note the mean was first by sample ID THEN by treatment
  scale_fill_manual(values=c("#56B4E9","#56B4E9","#56B4E9", "#E69F00","#E69F00","#E69F00")) +
  # scale_color_manual(values=c("#56B4E9","#E69F00")) +
  ggtitle("Day 7 vstExpression: DNMT1") +
  # expand_limits(y=0) +                                                    # Expand y range
  #scale_y_continuous(limits=c((min_p1), (max_p1))) +
  theme(text = element_text(size=15))


Day7_meanExpr %>% 
  dplyr::filter(GeneID %in% 'PGEN_.00g230260') %>% 
  ggplot(aes(x=group, y=mean.vstExp, fill=group)) +  # , colour=supp, group=supp))
  theme_classic() +
  geom_errorbar(aes(ymin=mean.vstExp-se.vsdtExp, ymax=mean.vstExp+se.vsdtExp), colour="black", width=.1, position=pd) +
  geom_line(position=pd) +
  geom_point(position=pd, size = 4, shape=21) +            
  xlab("PrimaryxSecond pCO2 treatment") +
  ylab("Gene Expression (mean±SE VST transformed)") +                 # note the mean was first by sample ID THEN by treatment
  scale_fill_manual(values=c("#56B4E9","#56B4E9","#56B4E9", "#E69F00","#E69F00","#E69F00")) +
  # scale_color_manual(values=c("#56B4E9","#E69F00")) +
  ggtitle("Day 7 vstExpression: SRT2") +
  # expand_limits(y=0) +                                                    # Expand y range
  #scale_y_continuous(limits=c((min_p1), (max_p1))) +
  theme(text = element_text(size=15))


# ===================================================================================
# Day 14
# ===================================================================================
Day14.ExpVST_GOIs_MELT <- melt(Day14.ExpVST_GOIs, id=(c('Sample.Name', 'Primary_Treatment', 'Second_Treament', 'group'))) # melt using reshape2
names(Day14.ExpVST_GOIs_MELT)[(5:6)] <- c('GeneID', 'vst_Expression') # change column names

Day14_meanExpr <- Day14.ExpVST_GOIs_MELT %>% 
  select(c('Sample.Name','group', 'vst_Expression', 'GeneID')) %>% 
  group_by(GeneID, group) %>%
  dplyr::summarize(mean.vstExp = mean(vst_Expression), 
                   sd.vsdtExp = sd(vst_Expression),
                   n = n(), 
                   se.vsdtExp = sd.vsdtExp/sqrt(n))

pd <- position_dodge(0.3) 
Day14_meanExpr %>% 
  dplyr::filter(GeneID %in% 'PGEN_.00g108770') %>% 
  ggplot(aes(x=group, y=mean.vstExp, fill=group)) +  # , colour=supp, group=supp))
  theme_classic() +
  geom_errorbar(aes(ymin=mean.vstExp-se.vsdtExp, ymax=mean.vstExp+se.vsdtExp), colour="black", width=.1, position=pd) +
  geom_line(position=pd) +
  geom_point(position=pd, size = 4, shape=21) +            
  xlab("PrimaryxSecond pCO2 treatment") +
  ylab("Gene Expression (mean±SE VST transformed)") +                 # note the mean was first by sample ID THEN by treatment
  scale_fill_manual(values=c("#56B4E9","#56B4E9","#56B4E9", "#E69F00","#E69F00","#E69F00")) +
  # scale_color_manual(values=c("#56B4E9","#E69F00")) +
  ggtitle("Day 14 vstExpression: Alternative Oxidase ") +
  # expand_limits(y=0) +                                                    # Expand y range
  #scale_y_continuous(limits=c((min_p1), (max_p1))) +
  theme(text = element_text(size=15))

Day14_meanExpr %>% 
  dplyr::filter(GeneID %in% 'PGEN_.00g283000') %>% 
  ggplot(aes(x=group, y=mean.vstExp, fill=group)) +  # , colour=supp, group=supp))
  theme_classic() +
  geom_errorbar(aes(ymin=mean.vstExp-se.vsdtExp, ymax=mean.vstExp+se.vsdtExp), colour="black", width=.1, position=pd) +
  geom_line(position=pd) +
  geom_point(position=pd, size = 4, shape=21) +            
  xlab("PrimaryxSecond pCO2 treatment") +
  ylab("Gene Expression (mean±SE VST transformed)") +                 # note the mean was first by sample ID THEN by treatment
  scale_fill_manual(values=c("#56B4E9","#56B4E9","#56B4E9", "#E69F00","#E69F00","#E69F00")) +
  # scale_color_manual(values=c("#56B4E9","#E69F00")) +
  ggtitle("Day 14 vstExpression: DNMT1") +
  # expand_limits(y=0) +                                                    # Expand y range
  #scale_y_continuous(limits=c((min_p1), (max_p1))) +
  theme(text = element_text(size=15))


Day14_meanExpr %>% 
  dplyr::filter(GeneID %in% 'PGEN_.00g230260') %>% 
  ggplot(aes(x=group, y=mean.vstExp, fill=group)) +  # , colour=supp, group=supp))
  theme_classic() +
  geom_errorbar(aes(ymin=mean.vstExp-se.vsdtExp, ymax=mean.vstExp+se.vsdtExp), colour="black", width=.1, position=pd) +
  geom_line(position=pd) +
  geom_point(position=pd, size = 4, shape=21) +            
  xlab("PrimaryxSecond pCO2 treatment") +
  ylab("Gene Expression (mean±SE VST transformed)") +                 # note the mean was first by sample ID THEN by treatment
  scale_fill_manual(values=c("#56B4E9","#56B4E9","#56B4E9", "#E69F00","#E69F00","#E69F00")) +
  # scale_color_manual(values=c("#56B4E9","#E69F00")) +
  ggtitle("Day 14 vstExpression: SRT2") +
  # expand_limits(y=0) +                                                    # Expand y range
  #scale_y_continuous(limits=c((min_p1), (max_p1))) +
  theme(text = element_text(size=15))


# SIRTUINS 
# PGEN_.00g048200 - sirtuin (SIR1)
# PGEN_.00g012340 - sirtuin (SIR2)
# PGEN_.00g149480 - siruitin (SIR4)
# PGEN_.00g153700 - sirtuin (SIR5)
# PGEN_.00g144540 - sirtuin (SIR6)
# PGEN_.00g033970 - sirtuin (SIR7) 