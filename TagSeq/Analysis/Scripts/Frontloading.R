---
# title: "DESeq2_FrontloadingAnalysis.R"
# author: "Samuel Gurr"
# date: "January 8, 2021"
---
  
  
# LOAD PACKAGES
library(dplyr)
library(reshape2)
# SET WORKING DIRECTORY AND LOAD DATA
setwd("C:/Users/samjg/Documents/Github_repositories/Pgenerosa_TagSeq_Metabolomics/TagSeq/")
  
# Lets explore Frontloading in DESeq2 data!!! 
# inspired by Barshis et al. 2013
# Reference: Barshis, D. J., Ladner, J. T., Oliver, T. A., Seneca, F. O., Traylor-Knowles, N., & Palumbi, S. R. (2013).
# Genomic basis for coral resilience to climate change. Proceedings of the National Academy of Sciences, 110(4), 1387-1392.

# load the P generosa annotated genes file
annot = read.delim2(file = "Seq_details/Panopea-generosa-genes-annotations.txt",header = F);

# Constutive frontloading are the genes showin both...
# (1) higher expression levels by the conditioned-ambient controls relative to naive-ambient controls AND 
# (2) gnees by conditioned animals with REDUCED response to stress (relative to the naive animals)

# Steps 
# (1) ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Need genes of interest, ideally those showing upregualtion ( positive Log fold change) in the naive animals (controls) in response to stress! 
# In our study for example, Day 7 we expsed the conditoined and naive animals to A (ambient), M (moderate), and S (severe) pCO2 conditions 
# thus.. the target upregualted genes are those in the AM vs. AA and AS vs. AA for day 7 (ambient rearing x moderate/severe expoure == AM and AS, respectively)
# however... we only have 12 total upregualted DEGs in this effect! insufficent to tackle this question 
# instead, and considering this repeated stress design, 
# we do have substacial genes upregulated under the moderate rexposure on day 7 as (A x M [second]; 92 downregulated genes - upregualted under moderate stress)
# the one caveat here is that these genes encompass a grouped pairwise effect as [AA + MA] versus [AM + MM] including both histories
# lets try is anyway! 

Day7AvMSecond_totalDEGs <- read.csv(file="Analysis/Output/DESeq2/10cpm/Day7/ Day7.SecondDEGs_AvM_DESeq2results.csv", sep=',', header=TRUE) %>% dplyr::select(!c('X','baseMean'))
Day7AvMSecond_UP        <- Day7AvMSecond_totalDEGs %>% 
  dplyr::filter(padj < 0.05) %>% 
  dplyr::filter(log2FoldChange < 0) %>%  
  mutate(log2FoldChange = abs(log2FoldChange))
  
nrow(Day7AvMSecond_UP) # 92 total genes! 
ncol(Day7AvMSecond_UP[7:ncol(Day7AvMSecond_UP)]) # 36 total samples! 

# Steps 
# (2) ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# we have out 92 genes in question! 'Day7AvMSecond_UP'
# now lets call in the experiment metadata to get the assigned treatment by the sample IDs! 
Master.Exp.Metadata   <- read.csv(file="Analysis/Data/Experiment_Metadata/Master_Phyenotype.and.Exp.Treatment_Metadata.csv", sep=',', header=TRUE) %>% 
  dplyr::filter(Date %in% 20190731) %>% # call only day 7 
  dplyr::select(c('Sample.Name', 'All_Treatment')) %>% # select only the sample names and the treatment groups 
  dplyr::mutate(All_Treatment = factor(All_Treatment)) # make the treatment groups a factor 
nrow(Master.Exp.Metadata) # 36 total samples! 

Day7UP_melted     <- Day7AvMSecond_UP %>% 
                        dplyr::select(!c('log2FoldChange', 'lfcSE','stat','pvalue','padj')) %>% 
                        reshape2::melt(id.var = 'Gene') %>% 
                        dplyr::rename(Sample.Name = variable)

Day7UP_TreatMerge <- merge(Day7UP_melted, Master.Exp.Metadata, by = 'Sample.Name') %>% 
                        dplyr::group_by(Gene, All_Treatment) %>% 
                        dplyr::select(!'Sample.Name') %>% 
                        dplyr::summarise(meanExp = mean(value))
View(Day7UP_TreatMerge)

Day7UP_READY <- dcast(Day7UP_TreatMerge, Gene ~ All_Treatment)


# Steps 
# (3) ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# we have out final dataset ready ro rock! 'Day7UP_READY'
# lets calculate the x and y axis of the frontloading figure following Barshis et al. 2013 criteria 

# X axis - this is the relative fold ratio of the condioned-control to the naive-control as the following: 
# [ (MM/MA) / (AM/AA) ]
# like a ratio of a ratio.. the values <1 will indicate genes that are lower response to stress than the naive animals (oppsite for values >1)

# lets calculate it 
xall_1 <- ( (Day7UP_READY$MM / Day7UP_READY$MA) / (Day7UP_READY$AM / Day7UP_READY$AA) ) # call MA as the control for the MM ratio
xall_2 <- ( (Day7UP_READY$MM / Day7UP_READY$AA) / (Day7UP_READY$AM / Day7UP_READY$AA) ) # call AA as the control for the MM ratio
xall_3 <- (Day7UP_READY$MM / Day7UP_READY$AM)  # call AA as the control for the MM ratio


Day7UP_READY$xall_1 <- xall_1
Day7UP_READY$xall_2 <- xall_2
Day7UP_READY$xall_3 <- xall_3

# Y Axis - this is simply the conditioned control over the naive control 
# ( MA / AA ) 

# lets calculate it 
yall <- (Day7UP_READY$MA / Day7UP_READY$AA)

Day7UP_READY$yall <- yall


# Steps 
# (4) ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# lets plot it! 
library(ggplot2)
P <- Day7UP_READY %>% #dplyr::filter(yall < 6) %>% 
        ggplot(aes(x=xall_1, y=yall)) +
        geom_point() +
        theme_classic() + 
        stat_smooth(method = "lm", 
                    formula = y ~ x + poly(x, 2) - 1) +
        geom_vline(xintercept=1, linetype="dotted") + 
        geom_hline(yintercept=1, linetype="dotted") + 
        labs(y= "Conditioned to naive control ratio", 
             x = "Conditioned to naive foldchange ratio",
             title = "Frontloaded genes") + 
        expand_limits(x = 0, y = 0) + 
  annotate("rect", xmin = 0, xmax = 1, ymin = 1, ymax = 1.5,
           alpha = .2) + 
  annotate("rect", xmin = 0, xmax = 1, ymin = 0, ymax = 1,
           alpha = .5)
P

# call genes with x axis vaulue < 1
frontloaded_genes <- Day7UP_READY %>% 
  dplyr::filter(xall_1 < 1) %>% 
  dplyr::filter(yall > 1) %>% 
  dplyr::select('Gene') # %>% 
# dplyr::filter(yall < 1) 


frontloadprobes = frontloaded_genes$Gene
probes2annot    = match(frontloadprobes, annot$V1)
frontloaded_genesANNOT = data.frame(geneSymbol = annot$V1[probes2annot],
                             Genes = annot$V7[probes2annot])
View(frontloaded_genesANNOT)


# :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: ::::::::::::::::::::::::::::::::::::::::::
# WGCNA data approach to frontloaded genes ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: ::::::::::::::::::::::::::::::::::::::::::

 # LOAD DATA 
day7.ModMem   <- read.csv(file="Analysis/Output/WGCNA/subseq_treatments_all/Day7/d7.WGCNA_ModulMembership.csv", sep=',', header=TRUE)   %>%  dplyr::rename(Pgen_ID = X) %>%  mutate(Pgen_ID = as.character(Pgen_ID)) %>%  dplyr::select(c('Pgen_ID','moduleColor'))
day14.ModMem  <- read.csv(file="Analysis/Output/WGCNA/subseq_treatments_all/Day14/d14.WGCNA_ModulMembership.csv", sep=',', header=TRUE)  %>%  dplyr::rename(Pgen_ID = X) %>%  mutate(Pgen_ID = as.character(Pgen_ID)) %>%  dplyr::select(c('Pgen_ID','moduleColor'))
day21.ModMem  <- read.csv(file="Analysis/Output/WGCNA/subseq_treatments_all/Day21/d21.WGCNA_ModulMembership.csv", sep=',', header=TRUE)  %>%  dplyr::rename(Pgen_ID = X) %>%  mutate(Pgen_ID = as.character(Pgen_ID)) %>%  dplyr::select(c('Pgen_ID','moduleColor'))

# ALL MODULES SHOWING HIGHER EXPRESSION BY PREEXPOSED CLAMS
# NOTE: no genes in these modules thorugh 
# all three days (day 7,14, 21), instead called all genes present on 2/3 days 
# (1) Preexposed_modules # modules with higher expression by preexposed/acclimate/stress-primed animals (relative to naive animals)
day7.yellow   <- day7.ModMem %>% dplyr::filter(moduleColor %in% 'yellow')  %>% dplyr::mutate(day ='Day7')
day14.black   <- day14.ModMem %>% dplyr::filter(moduleColor %in% 'black')  %>% dplyr::mutate(day ='Day14')
day21.yellow  <- day21.ModMem %>% dplyr::filter(moduleColor %in% 'yellow') %>% dplyr::mutate(day ='Day21')

Preexposed_modules   <- rbind(day7.yellow, day14.black, day21.yellow)
PreexpResponse_genes <- Preexposed_modules %>% 
  dplyr::group_by(Pgen_ID) %>% 
  dplyr::summarise(count = n()) %>% 
  dplyr::filter(count == 2)
# View(Frontloaded_genes) 
PreexpResponse_genes       = PreexpResponse_genes$Pgen_ID
PreexpResponse_genes_ANNOT = match(PreexpResponse_genes, annot$V1)
PreexpResponse_genes_data <- data.frame(geneSymbol = annot$V1[PreexpResponse_genes_ANNOT],
                Annotation = annot$V7[PreexpResponse_genes_ANNOT])
nrow(PreexpResponse_genes_data) # 243 genes present in 2/3 days of repeat exposures
# with the same expression pattern (higher expression by preexpoed animals)

# ALL MODULES SHOWING HIGHER EXPRESSION BY NAIVE CLAMS
# (2) Naive_modules # modules with higher expression by naive animals (relaive to preexposed animals)
day7.brown    <- day7.ModMem %>% dplyr::filter(moduleColor %in% 'brown')     %>% dplyr::mutate(day ='Day7')
day14.brown   <- day14.ModMem %>% dplyr::filter(moduleColor %in% 'brown')    %>% dplyr::mutate(day ='Day14')
day21.blue    <- day21.ModMem %>% dplyr::filter(moduleColor %in% 'blue')     %>% dplyr::mutate(day ='Day21')
day21.magenta <- day21.ModMem %>% dplyr::filter(moduleColor %in% 'magenta')  %>% dplyr::mutate(day ='Day21')

Naive_modules <- rbind(day7.brown, day14.brown, day21.blue, day21.magenta)
NaiveResponse_genes <- Naive_modules %>% 
  dplyr::group_by(Pgen_ID) %>% 
  dplyr::summarise(count = n()) %>% 
  dplyr::filter(count == 3)
# View(NaiveResponse_genes) 
NaiveResponse_genes        = NaiveResponse_genes$Pgen_ID
NaiveResponse_genes_ANNOT  = match(NaiveResponse_genes, annot$V1)
NaiveResponse_genes_data   <- data.frame(geneSymbol = annot$V1[NaiveResponse_genes_ANNOT],
                                         Annotation = annot$V7[NaiveResponse_genes_ANNOT])
nrow(NaiveResponse_genes_data) # 315 genes present on all sampling days 
# showing the same expression pattenr (ghigher expression by naive clams) 

# NEXT......
# this set of genes represents those that CONTINUOUSLY present in modules 
# highly correlated with the primary treatment for higher expression by NAIVE animals than preexposed animals 
# Thus, these genes are NOT present  any modules correlated with preexposed animals
# NOTE: are these the same as the TimeSeries WGNCA??
Annot_ModuleMembership <-  read.csv("Analysis/Output/WGCNA/subseq_treatments_all/TimeSeries/TimeSeries.WGCNA_ModulMembership.csv") %>% dplyr::select(c('geneSymbol','moduleColor'))

# Check the TimeSeries module vs. genes in separate modules ... % overlap TimeSeries - these new calls in this script
# (1) Prexposed genes of interest
# filter the 'NaiveResponse_genes_data' for genes in module BROWN in the TimeSeries WGCNA results
# (TimeSeries showing higher expression on all sampling days 7 14 and 21 than the preexposed animals)
Preex_mod_brown <-  Annot_ModuleMembership %>% dplyr::filter(moduleColor %in% 'brown')
nrow(Preex_mod_brown) # 1072 genes in the TimeSeries module brown
PreexpResponse_genes_FILT <- PreexpResponse_genes_data %>% dplyr::filter(geneSymbol %in% Preex_mod_brown$geneSymbol)
( (nrow(PreexpResponse_genes_FILT) / nrow(PreexpResponse_genes_data) )* 100) # 99.04762 only 3 genes not present in module yellow 

# Check the TimeSeries module vs. genes in separate modules ... % overlap TimeSeries - these new calls in this script
# (2) Naive  genes of interest
# filter the 'NaiveResponse_genes_data' for genes in module YELLOW in the TimeSeries WGCNA results
# (TimeSeries showing higher expression on all sampling days 7 14 and 21 than the preexposed animals)
Naive_mod_yellow <-  Annot_ModuleMembership %>% dplyr::filter(moduleColor %in% 'yellow')
nrow(Naive_mod_yellow)
NaiveResponse_genes_FILT <-NaiveResponse_genes_data %>% dplyr::filter(geneSymbol %in% Naive_mod_yellow$geneSymbol)
( (nrow(NaiveResponse_genes_FILT) / nrow(NaiveResponse_genes_data) )* 100) # 99.04762 only 3 genes not present in module yellow 


# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::;
# NOW LET'S EXPLORE FRONTLOADED GENE EXPRESSION!
# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::;

# instead of filtering by genes < padj or threshold fold change....
# filter this dataset of genes on Day 7 by the genes present 
# across ALL DAYS via WGCNA with higher expression by naive clams 
# rationale here is... 
# Genes present in the Day 7 matrix showed a larger quantity of sig DEs showing 
# differential expression under the second treatment A v M, however this treatment encompasses
# BOTH the naive and preexposed animals..
# Thus here we call all of the genes via WGCNA showin ghigher expression by Naive history

# Day 7 gene expression
WGCNA_presistantNaive_genes <- Day7AvMSecond_totalDEGs %>% 
  dplyr::filter(Gene %in% NaiveResponse_genes)
# View(WGCNA_presistantNaive_genes)

WGCNA_naive_melted     <- WGCNA_presistantNaive_genes %>% 
  dplyr::select(!c('log2FoldChange', 'lfcSE','stat','pvalue','padj')) %>% 
  reshape2::melt(id.var = 'Gene') %>% 
  dplyr::rename(Sample.Name = variable)

WGCNA_naive_Merge <- merge(WGCNA_naive_melted, Master.Exp.Metadata, by = 'Sample.Name') %>% 
  dplyr::group_by(Gene, All_Treatment) %>% 
  dplyr::select(!'Sample.Name') %>% 
  dplyr::summarise(meanExp = mean(value))

WGCNA_naive_READY <- dcast(WGCNA_naive_Merge, Gene ~ All_Treatment)

for (i in 1:nrow(WGCNA_naive_READY)) {
  
  # Moderate - higher expression AM > AA
        if (WGCNA_naive_READY$AM[i] > WGCNA_naive_READY$AA[i]) {
          # X axis            
          WGCNA_naive_READY$wgcna.xall_mod[i] <- ( (WGCNA_naive_READY$MM[i] / WGCNA_naive_READY$MA[i]) / (WGCNA_naive_READY$AM[i] / WGCNA_naive_READY$AA[i]) ) # call MA as the control for the MM ratio
          # Y Axis - this is simply the conditioned control over the naive control ( MA / AA ) 
          WGCNA_naive_READY$wgcna.yall_mod[i] <- (WGCNA_naive_READY$MA[i] / WGCNA_naive_READY$AA[i])
              } else {
                # X axis  - call NA          
                WGCNA_naive_READY$wgcna.xall_mod[i] <- NA 
                # Y Axis  - call NA
                WGCNA_naive_READY$wgcna.yall_mod[i] <- NA
              }
  # Severe - higher expression AS > AA
  
          if (WGCNA_naive_READY$AS[i] > WGCNA_naive_READY$AA[i]) {
            # X axis            
            WGCNA_naive_READY$wgcna.xall_sev[i] <- ( (WGCNA_naive_READY$MS[i] / WGCNA_naive_READY$MA[i]) / (WGCNA_naive_READY$AS[i] / WGCNA_naive_READY$AA[i]) ) # call AA as the control for the MM ratio
            # Y Axis - this is simply the conditioned control over the naive control ( MA / AA ) 
            WGCNA_naive_READY$wgcna.yall_sev[i] <- (WGCNA_naive_READY$MA[i] / WGCNA_naive_READY$AA[i])
                } else {
                  # X axis  - call NA          
                  WGCNA_naive_READY$wgcna.xall_sev[i] <- NA
                  # Y Axis  - call NA
                  WGCNA_naive_READY$wgcna.yall_sev[i] <- NA
                }
}
WGCNA_naive_READY



# Steps 
# (4) ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# lets plot it! 
library(ggplot2)
# plot the frontloaded genes formula using x axis for the moderate treatment (second exposure, day 7)
par(mfrow=c(1,2))
P.wgcna.mod <- WGCNA_naive_READY %>% 
  dplyr::select(!c('wgcna.xall_sev', 'wgcna.yall_sev')) %>% 
  na.omit() %>% 
  dplyr::filter(wgcna.yall_mod < 3) %>% 
  ggplot(aes(x=wgcna.xall_mod, y=wgcna.yall_mod)) +
  geom_point() +
  theme_classic() + 
  stat_smooth(method = "lm", 
              formula = y ~ x + poly(x, 2) - 1) +
  geom_vline(xintercept=1, linetype="dotted") + 
  geom_hline(yintercept=1, linetype="dotted") + 
  labs(y= "Conditioned to naive control ratio", 
       x = "Conditioned to naive foldchange ratio",
       title = "Day 7 
       Frontloaded genes; response to moderate pCO2") + 
  expand_limits(x = 0, y = 0) + 
  annotate("rect", xmin = 0, xmax = 1, ymin = 1, ymax = 1.5,
           alpha = .2) + 
  annotate("rect", xmin = 0, xmax = 1, ymin = 0, ymax = 1,
           alpha = .5)
P.wgcna.mod


# plot the frontloaded genes formula using x axis for the severe treatment (second exposure, day 7)
P.wgcna.sev <- WGCNA_naive_READY %>% 
  dplyr::select(!c('wgcna.xall_mod', 'wgcna.yall_mod')) %>% 
  na.omit() %>% 
  dplyr::filter(wgcna.yall_sev < 3) %>% 
  ggplot(aes(x=wgcna.xall_sev, y=wgcna.yall_sev)) +
  geom_point() +
  theme_classic() + 
  stat_smooth(method = "lm", 
              formula = y ~ x + poly(x, 2) - 1) +
  geom_vline(xintercept=1, linetype="dotted") + 
  geom_hline(yintercept=1, linetype="dotted") + 
  labs(y= "Conditioned to naive control ratio", 
       x = "Conditioned to naive foldchange ratio",
       title = "Day 7 Frontloaded genes; response to severe pCO2") + 
  expand_limits(x = 0, y = 0) + 
  annotate("rect", xmin = 0, xmax = 1, ymin = 1, ymax = 1.5,
           alpha = .2) + 
  annotate("rect", xmin = 0, xmax = 1, ymin = 0, ymax = 1,
           alpha = .5)
P.wgcna.sev




# call the genes in the fronloaded quadrant!
annot = read.delim2(file = "Seq_details/Panopea-generosa-genes-annotations.txt",header = F);
# MODERATE second treatment 
wgcna.frontloaded_genes_mod <- WGCNA_naive_READY %>% 
                                dplyr::filter(wgcna.xall_mod < 1) %>% 
                                dplyr::filter(wgcna.yall_mod > 1) %>%
                                dplyr::select('Gene') # %>% 

wgcna.frontloadprobes_mod        = wgcna.frontloaded_genes_mod$Gene
probes2annot_mod                 = match(wgcna.frontloadprobes_mod, annot$V1)
wgcna.frontloaded_mod_ANNOT      = data.frame(geneSymbol = annot$V1[probes2annot_mod],
                                    Genes = annot$V7[probes2annot_mod])
# View(wgcna.frontloaded_mod_ANNOT)

# SEVERE second treatment
wgcna.frontloaded_genes_sev <- WGCNA_naive_READY %>% 
  dplyr::filter(wgcna.xall_sev < 1) %>% 
  dplyr::filter(wgcna.yall_sev > 1) %>%
  dplyr::select('Gene') # %>% 

wgcna.frontloadprobes_sev        = wgcna.frontloaded_genes_sev$Gene
probes2annot_sev                 = match(wgcna.frontloadprobes_sev, annot$V1)
wgcna.frontloaded_sev_ANNOT      = data.frame(geneSymbol = annot$V1[probes2annot_sev],
                                              Genes = annot$V7[probes2annot_sev])
# View(wgcna.frontloaded_sev_ANNOT)


# Q: what is the overlap (f any) between the putative frontloaded genes calles
# as those with moderate-exposed animals and severe-exposed animals
# in other words, do we see frontloading in genes regardless of the subseqent exposure 
# and what are they? 

nrow(wgcna.frontloaded_sev_ANNOT)
ModSev_frontloaded <- wgcna.frontloaded_sev_ANNOT %>% 
  dplyr::filter(geneSymbol %in% wgcna.frontloaded_mod_ANNOT$geneSymbol)
nrow(ModSev_frontloaded) # 35 genes in the sev. dataset with the same gene as the moderate
nrow(wgcna.frontloaded_mod_ANNOT) # 35 total in the moderate; ALL GENES IN RESPONSE TO EITHER MODERATE OR SEVERE (7 extra in the severe)
View(ModSev_frontloaded)

Sec_7extra_frontloaded <- wgcna.frontloaded_sev_ANNOT %>% 
  dplyr::filter(!geneSymbol %in% wgcna.frontloaded_mod_ANNOT$geneSymbol)
nrow(Sec_7extra_frontloaded) # 7 genes in the sev. datasett, NOT shown in the moderate frontloaded
View(Sec_7extra_frontloaded)

# Question - what is the proportion of genes present in the 
# WGCNA method vs. the DESeq2 method 
# NOTE: DESeq2 genes are the A(second) v. M(second) group of sig DEGs < padj. 0.05
# WGCNA gene set for this code was the genes in co expression modules 
# sig associated with primary treatment for HIGHER expression by the naive phenotype 
# AND genes present in this same module pattern on days 7 14 and 21!!!

nrow(frontloaded_genesANNOT) # 13 putative frontloaded genes (DE analysis)
nrow(wgcna.frontloaded_genesANNOT) # 35  putative frontloaded genes (WGCNA modules)

# are the genes from DE analysis also present in the WGCNA set?
answer <- wgcna.frontloaded_genesANNOT %>% 
               dplyr::filter(geneSymbol %in% frontloaded_genesANNOT$geneSymbol)
answer # 5 of the 13 genes fromthe DE frontloaded are in the 70 genes from WGCNA frontloaded
