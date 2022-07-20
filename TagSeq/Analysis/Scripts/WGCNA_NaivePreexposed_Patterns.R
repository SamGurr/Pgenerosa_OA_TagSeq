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

PgenID_Gene <- annot %>% dplyr::select(c('V1', 'V7')) %>% dplyr::rename(Pgen_ID = V1) %>% dplyr::rename(gene = V7)

day7.ModMem   <- read.csv(file="Analysis/Output/WGCNA/subseq_treatments_all/Day7/d7.WGCNA_ModulMembership.csv", sep=',', header=TRUE)   %>%  dplyr::rename(Pgen_ID = X) %>%  mutate(Pgen_ID = as.character(Pgen_ID)) %>%  dplyr::select(c('Pgen_ID','moduleColor'))
day14.ModMem  <- read.csv(file="Analysis/Output/WGCNA/subseq_treatments_all/Day14/d14.WGCNA_ModulMembership.csv", sep=',', header=TRUE)  %>%  dplyr::rename(Pgen_ID = X) %>%  mutate(Pgen_ID = as.character(Pgen_ID)) %>%  dplyr::select(c('Pgen_ID','moduleColor'))
day21.ModMem  <- read.csv(file="Analysis/Output/WGCNA/subseq_treatments_all/Day21/d21.WGCNA_ModulMembership.csv", sep=',', header=TRUE)  %>%  dplyr::rename(Pgen_ID = X) %>%  mutate(Pgen_ID = as.character(Pgen_ID)) %>%  dplyr::select(c('Pgen_ID','moduleColor'))









# (1) Preexposed_modules # modules with higher expression by preexposed/acclimate/stress-primed animals (relative to naive animals)
day7.yellow   <- day7.ModMem %>% dplyr::filter(moduleColor %in% 'yellow')  %>% dplyr::mutate(day ='Day7')
day14.black   <- day14.ModMem %>% dplyr::filter(moduleColor %in% 'black')  %>% dplyr::mutate(day ='Day14')
day21.yellow  <- day21.ModMem %>% dplyr::filter(moduleColor %in% 'yellow') %>% dplyr::mutate(day ='Day21')

Preexposed_modules       <- rbind(day7.yellow, day14.black, day21.yellow)
Preexposed_modules_Genes <- merge(PgenID_Gene, Preexposed_modules, by = 'Pgen_ID')


PreexpResponse_Genes <- Preexposed_modules_Genes %>% 
  dplyr::group_by(Pgen_ID, gene) %>% 
  dplyr::summarise(day_count = n(), Days = as.character(list(unique(day)))) %>% 
  dplyr::filter(day_count == 2) 
PreexpResponse_Genes$Days <-  gsub("\"", "", PreexpResponse_Genes$Days)
PreexpResponse_Master <- PreexpResponse_Genes %>% 
  dplyr::mutate(Category = case_when(Days == "c(Day21, Day7)" ~ "rapid_induced_under_exposures_7&21",  
                                     Days == "c(Day7, Day21)" ~ "rapid_induced_under_exposures_7&21", 
                                     Days == "c(Day14, Day7)" ~ "stress_induc_recovery_7&14", 
                                     Days == "c(Day7, Day14)" ~ "stress_induc_recovery_7&14", 
                                     Days == "c(Day21, Day14)" ~ "recovery_prep_reg_14&21",
                                     Days == "c(Day14, Day21)" ~ "recovery_prep_reg_14&21"))
View(PreexpResponse_Master) 
write.csv(PreexpResponse_Master, "C:/Users/samjg/Documents/Github_repositories/Pgenerosa_TagSeq_Metabolomics/TagSeq/Analysis/Output/WGCNA/subseq_treatments_all/PreexposedResponse_WGCNA_Patterns.csv")










# (2) Naive_modules # modules with higher expression by naive animals (relaive to preexposed animals)
day7.brown    <- day7.ModMem %>% dplyr::filter(moduleColor %in% 'brown')     %>% dplyr::mutate(day ='Day7')
day14.brown   <- day14.ModMem %>% dplyr::filter(moduleColor %in% 'brown')    %>% dplyr::mutate(day ='Day14')
day21.blue    <- day21.ModMem %>% dplyr::filter(moduleColor %in% 'blue')     %>% dplyr::mutate(day ='Day21')
day21.magenta <- day21.ModMem %>% dplyr::filter(moduleColor %in% 'magenta')  %>% dplyr::mutate(day ='Day21')

Naive_modules       <- rbind(day7.brown, day14.brown, day21.blue, day21.magenta)
Naive_modules_Genes <- merge(PgenID_Gene, Naive_modules, by = 'Pgen_ID')
NaiveResponse_Genes <- Naive_modules_Genes %>% 
  dplyr::group_by(Pgen_ID, gene) %>% 
  dplyr::summarise(day_count = n(), Days = as.character(list(unique(day)))) %>% 
  dplyr::filter(day_count >1) 

NaiveResponse_Genes$Days <-  gsub("\"", "", NaiveResponse_Genes$Days)
NaiveResponse_Master <- NaiveResponse_Genes %>% 
  dplyr::mutate(Category = case_when(Days == "c(Day21, Day7)" ~ "rapid_induced_under_exposures_7&21",  
                                     Days == "c(Day7, Day21)" ~ "rapid_induced_under_exposures_7&21", 
                                     Days == "c(Day14, Day7)" ~ "stress_induc_recovery_7&14", 
                                     Days == "c(Day7, Day14)" ~ "stress_induc_recovery_7&14", 
                                     Days == "c(Day21, Day14)" ~ "recovery_prep_reg_14&21",
                                     Days == "c(Day14, Day21)" ~ "recovery_prep_reg_14&21",
                                     TRUE                      ~  "all_days"))
View(NaiveResponse_Master) 
write.csv(NaiveResponse_Master, "C:/Users/samjg/Documents/Github_repositories/Pgenerosa_TagSeq_Metabolomics/TagSeq/Analysis/Output/WGCNA/subseq_treatments_all/NaiveResponse_WGCNA_Patterns.csv")

