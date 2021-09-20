---
# title: "WGCNA_frontloaded_genes.R"
# author: "Samuel Gurr"
# date: "January 8, 2021"
---
  
# LOAD PACKAGES
library(WGCNA) # note: this was previously installed with the command `BiocManager::install("WGCNA")`
library(dplyr)
library(tidyr)
library(tidyverse)

# SET WORKING DIRECTORY AND LOAD DATA
setwd("C:/Users/samjg/Documents/Github_repositories/Pgenerosa_TagSeq_Metabolomics/TagSeq/")
# LOAD DATA
annot = read.delim2(file = "Seq_details/Panopea-generosa-genes-annotations.txt",header = F);

day7.ModMem   <- read.csv(file="Analysis/Output/WGCNA/subseq_treatments_all/Day7/d7.WGCNA_ModulMembership.csv", sep=',', header=TRUE)   %>%  dplyr::rename(Pgen_ID = X) %>%  mutate(Pgen_ID = as.character(Pgen_ID)) %>%  dplyr::select(c('Pgen_ID','moduleColor'))
day14.ModMem  <- read.csv(file="Analysis/Output/WGCNA/subseq_treatments_all/Day14/d14.WGCNA_ModulMembership.csv", sep=',', header=TRUE)  %>%  dplyr::rename(Pgen_ID = X) %>%  mutate(Pgen_ID = as.character(Pgen_ID)) %>%  dplyr::select(c('Pgen_ID','moduleColor'))
day21.ModMem  <- read.csv(file="Analysis/Output/WGCNA/subseq_treatments_all/Day21/d21.WGCNA_ModulMembership.csv", sep=',', header=TRUE)  %>%  dplyr::rename(Pgen_ID = X) %>%  mutate(Pgen_ID = as.character(Pgen_ID)) %>%  dplyr::select(c('Pgen_ID','moduleColor'))


# (1) Preexposed_modules # modules with higher expression by preexposed/acclimate/stress-primed animals (relative to naive animals)
day7.yellow   <- day7.ModMem %>% dplyr::filter(moduleColor %in% 'yellow')  %>% dplyr::mutate(day ='Day7')
day14.black   <- day14.ModMem %>% dplyr::filter(moduleColor %in% 'black')  %>% dplyr::mutate(day ='Day14')
day21.yellow  <- day21.ModMem %>% dplyr::filter(moduleColor %in% 'yellow') %>% dplyr::mutate(day ='Day21')

Preexposed_modules <- rbind(day7.yellow, day14.black, day21.yellow)

Frontloaded_genes <- Preexposed_modules %>% 
  dplyr::group_by(Pgen_ID) %>% 
  dplyr::summarise(count = n()) %>% 
  dplyr::filter(count == 2)
View(Frontloaded_genes) 

Frontloaded.probes = Frontloaded_genes$Pgen_ID
Frontloaded.probes_ANNOT = match(Frontloaded.probes, annot$V1)
View(data.frame(geneSymbol = annot$V1[Frontloaded.probes_ANNOT],
                Annotation = annot$V7[Frontloaded.probes_ANNOT]))

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
View(NaiveResponse_genes) 

NaiveResponse_genes        = Persistant_genes$Pgen_ID
NaiveResponse_genes_ANNOT  = match(NaiveResponse_genes, annot$V1)
NaiveResponse_genes_data   <- data.frame(geneSymbol = annot$V1[NaiveResponse_genes_ANNOT],
                                    Annotation = annot$V7[NaiveResponse_genes_ANNOT])
# this set of genes represents those that CONTINUOUSLY present in modules 
# highly correlated with the primary treatment for higher expression by NAIVE animals than preexposed animals 
# Thus, these genes are NOT presentin any modules correlated with preexposed animals
# NOTE: are these the same as the TimeSeries WGNCA??

Annot_ModuleMembership <-  read.csv("Analysis/Output/WGCNA/subseq_treatments_all/TimeSeries/TimeSeries.WGCNA_ModulMembership.csv") %>% dplyr::select(c('geneSymbol','moduleColor'))
Naive_mod_yellow       <-  Annot_ModuleMembership %>% dplyr::filter(moduleColor %in% 'yellow')
nrow(Naive_mod_yellow)
# filter the 'NaiveResponse_genes_data' for genes in module yellow (TimeSeries showing higher expression on all sampling days 7 14 and 21 than the preexposed animals)
NaiveResponse_genes_FILT <-NaiveResponse_genes_data %>% dplyr::filter(geneSymbol %in% Naive_mod_yellow$geneSymbol)
  
( (nrow(NaiveResponse_genes_FILT) / nrow(NaiveResponse_genes_data) )* 100) # 99.04762 only 3 genes not present in module yellow 
