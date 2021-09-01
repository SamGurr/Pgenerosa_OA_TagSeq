---
  # title: "compileGO_Master"
  # author: "Samuel Gurr"
  # date: "Sept 1, 2021"
---

  
# Purpose: :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

# (1) UNFILTERED GO term data - WGCNA modules only
# The objective of this script is to complile all of the enriched gene ontology terms in significant WGCNa modules into one master
# cumulative csv file. This file will be used as the supplementary material for our manuscript - data in our manuscript contained a 
# filtering step in which only BP terms with >=10 genes and only MF terms with >= 2 genes were used for goslim analysis
  # PRODUCT: All_GOterms_unfiltered_Master.csv

  
# (2) UNFILTERED GO term data - DESEq2 results only
# The objective of this script is to complile all of the enriched gene ontology terms in significant WGCNa modules into one master
# cumulative csv file. This file will be used as the supplementary material for our manuscript - data in our manuscript contained a 
# filtering step in which only BP terms with >=10 genes and only MF terms with >= 2 genes were used for goslim analysis
# PRODUCT: All_GOterms_unfiltered_Master.csv


# (3) GOSlim data - WGCNA modules only
# The next part of this script is to compile all of the GOslim results from the FILTERED GO enrichment analysis (NOT the unfiltered in #1)
  # PRODUCT: All_GOSlim_Master.csv


# SET WORKING DIRECTORY   ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::;: #

setwd("C:/Users/samjg/Documents/Github_repositories/Pgenerosa_TagSeq_Metabolomics/TagSeq/Analysis/Output/GO/WGCNA_goseq/")


# LOAD LIBRARIES  ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::;: #

library(dplyr)


# (1) 
# :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Compile all UNFILTERED GO term data! "...GENE_REFERENCE.csv" files! 
# :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Day 0 

# WGCNA results (all treatments ) ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
d0_midnightblue  <- read.csv("subseq_treatments_all/Day0/GO.05midnightblueModule_unfiltered.csv")
d0_midnightblue <- d0_midnightblue %>% dplyr::filter(ontology %in% c('BP', 'MF')) %>% 
  dplyr::mutate(Day_moduleColor = paste(Day,moduleColor, sep = '_')) %>% 
  dplyr::select(!c('moduleColor' , 'Day')) %>% 
  arrange(ontology, desc(numDEInCat))
# View(d0_midnightblue)


# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Day 7 


d7_brown  <- read.csv("subseq_treatments_all/Day7/GO.05brownModule_unfiltered.csv")
d7_brown <- d7_brown %>% dplyr::filter(ontology %in% c('BP', 'MF')) %>% 
  dplyr::mutate(Day_moduleColor = paste(Day,moduleColor, sep = '_')) %>% 
  dplyr::select(!c('moduleColor' ,'Day')) %>% 
  arrange(ontology, desc(numDEInCat))
# View(d7_brown)

d7_green  <- read.csv("subseq_treatments_all/Day7/GO.05brownModule_unfiltered.csv")
d7_green <- d7_green %>% dplyr::filter(ontology %in% c('BP', 'MF')) %>% 
  dplyr::mutate(Day_moduleColor = paste(Day,moduleColor, sep = '_')) %>% 
  dplyr::select(!c('moduleColor' ,'Day')) %>% 
  arrange(ontology, desc(numDEInCat))
# View(d7_green)

d7_yellow  <- read.csv("subseq_treatments_all/Day7/GO.05brownModule_unfiltered.csv")
d7_yellow <- d7_yellow %>% dplyr::filter(ontology %in% c('BP', 'MF')) %>% 
  dplyr::mutate(Day_moduleColor = paste(Day,moduleColor, sep = '_')) %>% 
  dplyr::select(!c('moduleColor' ,'Day')) %>% 
  arrange(ontology, desc(numDEInCat))
# View(d7_yellow)


# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Day 14 


d14_brown  <- read.csv("subseq_treatments_all/Day14/GO.05brownModule_unfiltered.csv")
d14_brown <- d14_brown %>% dplyr::filter(ontology %in% c('BP', 'MF')) %>% 
  dplyr::mutate(Day_moduleColor = paste(Day,moduleColor, sep = '_')) %>% 
  dplyr::select(!c('moduleColor' ,'Day')) %>% 
  arrange(ontology, desc(numDEInCat))
# View(d14_brown)

d14_black  <- read.csv("subseq_treatments_all/Day14/GO.05blackModule_unfiltered.csv")
d14_black <- d14_black %>% dplyr::filter(ontology %in% c('BP', 'MF')) %>% 
  dplyr::mutate(Day_moduleColor = paste(Day,moduleColor, sep = '_')) %>% 
  dplyr::select(!c('moduleColor' ,'Day')) %>% 
  arrange(ontology, desc(numDEInCat))
# View(d14_black)

d14_magenta  <- read.csv("subseq_treatments_all/Day14/GO.05magentaModule_unfiltered.csv")
d14_magenta <- d14_magenta %>% dplyr::filter(ontology %in% c('BP', 'MF')) %>% 
  dplyr::mutate(Day_moduleColor = paste(Day,moduleColor, sep = '_')) %>% 
  dplyr::select(!c('moduleColor' ,'Day')) %>% 
  arrange(ontology, desc(numDEInCat))
# View(d14_magenta)

d14_pink  <- read.csv("subseq_treatments_all/Day14/GO.05pinkModule_unfiltered.csv")
d14_pink <- d14_pink %>% dplyr::filter(ontology %in% c('BP', 'MF')) %>% 
  dplyr::mutate(Day_moduleColor = paste(Day,moduleColor, sep = '_')) %>% 
  dplyr::select(!c('moduleColor' ,'Day')) %>% 
  arrange(ontology, desc(numDEInCat))
# View(d14_pink)



# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Day 21 


d21_blue  <- read.csv("subseq_treatments_all/Day21/GO.05blueModule_unfiltered.csv")
d21_blue <- d21_blue %>% dplyr::filter(ontology %in% c('BP', 'MF')) %>% 
  dplyr::mutate(Day_moduleColor = paste(Day,moduleColor, sep = '_')) %>% 
  dplyr::select(!c('moduleColor' ,'Day')) %>% 
  arrange(ontology, desc(numDEInCat))
# View(d21_blue)

d21_black  <- read.csv("subseq_treatments_all/Day21/GO.05blackModule_unfiltered.csv")
d21_black <- d21_black %>% dplyr::filter(ontology %in% c('BP', 'MF')) %>% 
  dplyr::mutate(Day_moduleColor = paste(Day,moduleColor, sep = '_')) %>% 
  dplyr::select(!c('moduleColor' ,'Day')) %>% 
  arrange(ontology, desc(numDEInCat))
# View(d21_black)

d21_magenta  <- read.csv("subseq_treatments_all/Day21/GO.05magentaModule_unfiltered.csv")
d21_magenta <- d21_magenta %>% dplyr::filter(ontology %in% c('BP', 'MF')) %>% 
  dplyr::mutate(Day_moduleColor = paste(Day,moduleColor, sep = '_')) %>% 
  dplyr::select(!c('moduleColor' ,'Day')) %>% 
  arrange(ontology, desc(numDEInCat))
# View(d21_magenta)

d21_pink  <- read.csv("subseq_treatments_all/Day21/GO.05pinkModule_unfiltered.csv")
d21_pink <- d21_pink %>% dplyr::filter(ontology %in% c('BP', 'MF')) %>% 
  dplyr::mutate(Day_moduleColor = paste(Day,moduleColor, sep = '_')) %>% 
  dplyr::select(!c('moduleColor' ,'Day')) %>% 
  arrange(ontology, desc(numDEInCat))
# View(d21_pink)

d21_red  <- read.csv("subseq_treatments_all/Day21/GO.05redModule_unfiltered.csv")
d21_red <- d21_red %>% dplyr::filter(ontology %in% c('BP', 'MF')) %>% 
  dplyr::mutate(Day_moduleColor = paste(Day,moduleColor, sep = '_')) %>% 
  dplyr::select(!c('moduleColor' ,'Day')) %>% 
  arrange(ontology, desc(numDEInCat))
# View(d21_red)

d21_turquoise  <- read.csv("subseq_treatments_all/Day21/GO.05turquoiseModule_unfiltered.csv")
d21_turquoise <- d21_turquoise %>% dplyr::filter(ontology %in% c('BP', 'MF')) %>% 
  dplyr::mutate(Day_moduleColor = paste(Day,moduleColor, sep = '_')) %>% 
  dplyr::select(!c('moduleColor' ,'Day')) %>% 
  arrange(ontology, desc(numDEInCat))
# View(d21_turquoise)


# rbind to a cumulative  file
GO_Master_AllTreatments <- rbind(d0_midnightblue,
                                 
                                 d7_brown, 
                                 d7_green, 
                                 d7_yellow, 
                                 
                                 d14_black, 
                                 d14_brown, 
                                 d14_pink,
                                 d14_magenta, 
                                 
                                 d21_black, 
                                 d21_blue,  
                                 d21_pink, 
                                 d21_magenta, 
                                 d21_red, 
                                 d21_turquoise)

GO_Master_AllTreatments <- GO_Master_AllTreatments[,-1] # omit the column 'X' just the row numbers...
write.csv(GO_Master_AllTreatments, file = paste("subseq_treatments_all/All_GOterms_unfiltered_Master.csv", sep ='')) # write the file






# (2) 
# :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Compile all UNFILTERED GO term data for DESeq2 resutls!!!!  "...GENE_REFERENCE.csv" files! 
# :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::


# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Day 0 

d0_PrimEff_downreg  <- read.csv("../DESeq2_goseq/Day0/GO.05.Day0_PrimEffect_Downregulated_unfiltered.csv") %>% dplyr::mutate(Day_effect = 'Day0_Primary_AvM') %>% dplyr::mutate(Dir = 'Downregulated') %>% dplyr::filter(!ontology %in% c('CC', 'NA')) %>% arrange(ontology, desc(numDEInCat))
d0_PrimEff_upreg  <- read.csv("../DESeq2_goseq/Day0/GO.05.Day0_PrimEffect_Upregulated_unfiltered.csv")     %>%  dplyr::mutate(Day_effect = 'Day0_Primary_AvM') %>% dplyr::mutate(Dir = 'Upregulated') %>% dplyr::filter(!ontology %in% c('CC', 'NA'))  %>% arrange(ontology, desc(numDEInCat))

d7_PrimEff_downreg <- read.csv("../DESeq2_goseq/Day7/primary_effect/GO.05.Day7_PrimEffect_Downregulated_unfiltered.csv") %>%  dplyr::mutate(Day_effect = 'Day7_Primary_AvM') %>% dplyr::mutate(Dir = 'Downregulated') %>% dplyr::filter(!ontology %in% c('CC', 'NA')) %>% arrange(ontology, desc(numDEInCat))
d7_PrimEff_upreg   <- read.csv("../DESeq2_goseq/Day7/primary_effect/GO.05.Day7_PrimEffect_Upregulated_unfiltered.csv")   %>%  dplyr::mutate(Day_effect = 'Day7_Primary_AvM') %>% dplyr::mutate(Dir = 'Upregulated')  %>% dplyr::filter(!ontology %in% c('CC', 'NA'))  %>% arrange(ontology, desc(numDEInCat))

d7_SecondEffect_AvM_downreg   <- read.csv("../DESeq2_goseq/Day7/second_effect_AvM/GO.05.Day7_SecondEffect_AvM_Downregulated_unfiltered.csv") %>%  dplyr::mutate(Day_effect = 'Day7_Second_AvM') %>% dplyr::mutate(Dir = 'Downregulated') %>% dplyr::filter(!ontology %in% c('CC', 'NA')) %>% arrange(ontology, desc(numDEInCat))
d7_SecondEffect_AvM_upreg     <- read.csv("../DESeq2_goseq/Day7/second_effect_AvM/GO.05.Day7_SecondEffect_AvM_Upregulated_unfiltered.csv")   %>%  dplyr::mutate(Day_effect = 'Day7_Second_AvM') %>% dplyr::mutate(Dir = 'Upregulated')   %>% dplyr::filter(!ontology %in% c('CC', 'NA')) %>% arrange(ontology, desc(numDEInCat))

d7_GroupEffect_MAvAM_downreg      <- read.csv("../DESeq2_goseq/Day7/group_effect_MAvAM/GO.05.Day7_GroupEffect_MAvAM_Downregulated_unfiltered.csv") %>%  dplyr::mutate(Day_effect = 'Day7_Group_MAvAM') %>% dplyr::mutate(Dir = 'Downregulated') %>% dplyr::filter(!ontology %in% c('CC', 'NA')) %>% arrange(ontology, desc(numDEInCat))
d7_GroupEffect_MAvAM_upreg        <- read.csv("../DESeq2_goseq/Day7/group_effect_MAvAM/GO.05.Day7_GroupEffect_MAvAM_Upregulated_unfiltered.csv")   %>%  dplyr::mutate(Day_effect = 'Day7_Group_MAvAM') %>% dplyr::mutate(Dir = 'Upregulated')   %>% dplyr::filter(!ontology %in% c('CC', 'NA')) %>% arrange(ontology, desc(numDEInCat))

d14_PrimEff_downreg  <- read.csv("../DESeq2_goseq/Day14/GO.05.Day14_PrimEffect_Downregulated_unfiltered.csv") %>%  dplyr::mutate(Day_effect = 'Day14_Primary_AvM') %>% dplyr::mutate(Dir = 'Downregulated') %>% dplyr::filter(!ontology %in% c('CC', 'NA')) %>% arrange(ontology, desc(numDEInCat))
d14_PrimEff_upreg  <- read.csv("../DESeq2_goseq/Day14/GO.05.Day14_PrimEffect_Upregulated_unfiltered.csv")     %>%  dplyr::mutate(Day_effect = 'Day14_Primary_AvM') %>% dplyr::mutate(Dir = 'Upregulated')   %>% dplyr::filter(!ontology %in% c('CC', 'NA')) %>% arrange(ontology, desc(numDEInCat))

d21_PrimEff_downreg  <- read.csv("../DESeq2_goseq/Day21/GO.05.Day21_PrimEffect_Downregulated_unfiltered.csv") %>%  dplyr::mutate(Day_effect = 'Day21_Primary_AvM') %>% dplyr::mutate(Dir = 'Downregulated') %>% dplyr::filter(!ontology %in% c('CC', 'NA')) %>% arrange(ontology, desc(numDEInCat))
d21_PrimEff_upreg  <- read.csv("../DESeq2_goseq/Day21/GO.05.Day21_PrimEffect_Upregulated_unfiltered.csv")     %>%  dplyr::mutate(Day_effect = 'Day21_Primary_AvM') %>% dplyr::mutate(Dir = 'Upregulated')   %>% dplyr::filter(!ontology %in% c('CC', 'NA')) %>% arrange(ontology, desc(numDEInCat))


# compile master file
GO_Master_DESeq2 <- rbind(d0_PrimEff_downreg, d0_PrimEff_upreg,
                      d7_PrimEff_downreg, d7_PrimEff_upreg,
                      d7_SecondEffect_AvM_downreg, d7_SecondEffect_AvM_upreg,
                      d7_GroupEffect_MAvAM_downreg, d7_GroupEffect_MAvAM_upreg,
                      d14_PrimEff_downreg, d14_PrimEff_upreg,
                      d21_PrimEff_downreg, d21_PrimEff_upreg)

View(GO_Master_DESeq2)

write.csv(GO_Master_DESeq2, file = paste("../DESeq2_goseq/DESEq2_GOterms_unfiltered_Master.csv", sep ='')) # write the file


  
  
  
  
  
  

# (3) 
# :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Compile all GOslim data! "...GENE_REFERENCE.csv" files! 
# :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::


# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Day 0 

# WGCNA results (all treatments ) ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
d0_midnightblue_BP  <- read.csv("subseq_treatments_all/Day0/GOterms_and_GOslim_BiolProc_midnightblueModule_GENE_REFERENCE.csv")
d0_midnightblue_BP <- d0_midnightblue_BP %>% 
  select(c('PgenIDs','slim_term','term','Gene_terms',)) %>% 
  dplyr::mutate(day_mod = 'Day0_midnightblue') %>% dplyr::mutate(type = 'biological process') %>% 
  arrange(factor(term))
View(d0_midnightblue_BP)

d0_midnightblue_MF  <- read.csv("subseq_treatments_all/Day0/GOterms_and_GOslim_MolFunction_midnightblueModule_GENE_REFERENCE.csv")
d0_midnightblue_MF <- d0_midnightblue_MF %>% 
  select(c('PgenIDs','slim_term','term','Gene_terms',)) %>% 
  dplyr::mutate(day_mod = 'Day0_midnightblue') %>% dplyr::mutate(type = 'molecular function') %>% 
  arrange(factor(term))
View(d0_midnightblue_MF)


# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Day 7 


d7_brown_BP  <- read.csv("subseq_treatments_all/Day7/GOterms_and_GOslim_BiolProc_brownModule_GENE_REFERENCE.csv")
d7_brown_BP <- d7_brown_BP %>% 
  select(c('PgenIDs','slim_term','term','Gene_terms',)) %>% 
  dplyr::mutate(day_mod = 'Day7_brown') %>% dplyr::mutate(type = 'biological process') %>% 
  arrange(factor(term))
# View(d7_brown_BP)
d7_brown_MF  <- read.csv("subseq_treatments_all/Day7/GOterms_and_GOslim_MolFunction_brownModule_GENE_REFERENCE.csv")
d7_brown_MF <- d7_brown_MF %>% 
  select(c('PgenIDs','slim_term','term','Gene_terms',)) %>% 
  dplyr::mutate(day_mod = 'Day7_brown') %>% dplyr::mutate(type = 'molecular function') %>% 
  arrange(factor(term))
# View(d7_brown_MF)



d7_green_BP  <- read.csv("subseq_treatments_all/Day7/GOterms_and_GOslim_BiolProc_greenModule_GENE_REFERENCE.csv")
d7_green_BP <- d7_green_BP %>% 
  select(c('PgenIDs','slim_term','term','Gene_terms',)) %>% 
  dplyr::mutate(day_mod = 'Day7_green') %>% dplyr::mutate(type = 'biological process') %>% 
  arrange(factor(term))
# View(d7_green_BP)
d7_green_MF  <- read.csv("subseq_treatments_all/Day7/GOterms_and_GOslim_MolFunction_greenModule_GENE_REFERENCE.csv")
d7_green_MF <- d7_green_MF %>% 
  select(c('PgenIDs','slim_term','term','Gene_terms',)) %>% 
  dplyr::mutate(day_mod = 'Day7_green') %>% dplyr::mutate(type = 'molecular function') %>% 
  arrange(factor(term))
# View(d7_green_MF)


d7_yellow_BP  <- read.csv("subseq_treatments_all/Day7/GOterms_and_GOslim_BiolProc_yellowModule_GENE_REFERENCE.csv")
d7_yellow_BP <- d7_yellow_BP %>% 
  select(c('PgenIDs','slim_term','term','Gene_terms',)) %>% 
  dplyr::mutate(day_mod = 'Day7_yellow') %>% dplyr::mutate(type = 'biological process') %>% 
  arrange(factor(term))
# View(d7_yellow_BP)
d7_yellow_MF  <- read.csv("subseq_treatments_all/Day7/GOterms_and_GOslim_MolFunction_yellowModule_GENE_REFERENCE.csv")
d7_yellow_MF <- d7_yellow_MF %>% 
  select(c('PgenIDs','slim_term','term','Gene_terms',)) %>% 
  dplyr::mutate(day_mod = 'Day7_yellow') %>% dplyr::mutate(type = 'molecular function') %>% 
  arrange(factor(term))
# View(d7_yellow_MF)



# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Day 14 


d14_brown_BP  <- read.csv("subseq_treatments_all/Day14/GOterms_and_GOslim_BiolProc_brownModule_GENE_REFERENCE.csv")
d14_brown_BP <- d14_brown_BP %>% 
  select(c('PgenIDs','slim_term','term','Gene_terms',)) %>% 
  dplyr::mutate(day_mod = 'Day14_brown') %>% dplyr::mutate(type = 'biological process') %>% 
  arrange(factor(term))
# View(d14_brown_BP)
d14_brown_MF  <- read.csv("subseq_treatments_all/Day14/GOterms_and_GOslim_MolFunction_brownModule_GENE_REFERENCE.csv")
d14_brown_MF <- d14_brown_MF %>% 
  select(c('PgenIDs','slim_term','term','Gene_terms',)) %>% 
  dplyr::mutate(day_mod = 'Day14_brown') %>% dplyr::mutate(type = 'molecular function') %>% 
  arrange(factor(term))
# View(d14_brown_MF)



d14_black_BP  <- read.csv("subseq_treatments_all/Day14/GOterms_and_GOslim_BiolProc_blackModule_GENE_REFERENCE.csv")
d14_black_BP <- d14_black_BP %>% 
  select(c('PgenIDs','slim_term','term','Gene_terms',)) %>% 
  dplyr::mutate(day_mod = 'Day14_black') %>% dplyr::mutate(type = 'biological process') %>% 
  arrange(factor(term))
# View(d14_black_BP)
d14_black_MF  <- read.csv("subseq_treatments_all/Day14/GOterms_and_GOslim_MolFunction_blackModule_GENE_REFERENCE.csv")
d14_black_MF <- d14_black_MF %>% 
  select(c('PgenIDs','slim_term','term','Gene_terms',)) %>% 
  dplyr::mutate(day_mod = 'Day14_black') %>% dplyr::mutate(type = 'molecular function') %>% 
  arrange(factor(term))
# View(d14_black_MF)


d14_magenta_BP  <- read.csv("subseq_treatments_all/Day14/GOterms_and_GOslim_BiolProc_magentaModule_GENE_REFERENCE.csv")
d14_magenta_BP <- d14_magenta_BP %>% 
  select(c('PgenIDs','slim_term','term','Gene_terms',)) %>% 
  dplyr::mutate(day_mod = 'Day14_magenta') %>% dplyr::mutate(type = 'biological process') %>% 
  arrange(factor(term))
# View(d14_magenta_BP)
d14_magenta_MF  <- read.csv("subseq_treatments_all/Day14/GOterms_and_GOslim_MolFunction_magentaModule_GENE_REFERENCE.csv")
d14_magenta_MF <- d14_magenta_MF %>% 
  select(c('PgenIDs','slim_term','term','Gene_terms',)) %>% 
  dplyr::mutate(day_mod = 'Day14_magenta') %>% dplyr::mutate(type = 'molecular function') %>% 
  arrange(factor(term))
# View(d14_magenta_MF)



d14_pink_BP  <- read.csv("subseq_treatments_all/Day14/GOterms_and_GOslim_BiolProc_pinkModule_GENE_REFERENCE.csv")
d14_pink_BP <- d14_pink_BP %>% 
  select(c('PgenIDs','slim_term','term','Gene_terms',)) %>% 
  dplyr::mutate(day_mod = 'Day14_pink') %>% dplyr::mutate(type = 'biological process') %>% 
  arrange(factor(term))
# View(d14_pink_BP)
d14_pink_MF  <- read.csv("subseq_treatments_all/Day14/GOterms_and_GOslim_MolFunction_pinkModule_GENE_REFERENCE.csv")
d14_pink_MF <- d14_pink_MF %>% 
  select(c('PgenIDs','slim_term','term','Gene_terms',)) %>% 
  dplyr::mutate(day_mod = 'Day14_pink') %>% dplyr::mutate(type = 'molecular function') %>% 
  arrange(factor(term))
# View(d14_pink_MF)



# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Day 21 


d21_blue_BP  <- read.csv("subseq_treatments_all/Day21/GOterms_and_GOslim_BiolProc_blueModule_GENE_REFERENCE.csv")
d21_blue_BP <- d21_blue_BP %>% 
  select(c('PgenIDs','slim_term','term','Gene_terms',)) %>% 
  dplyr::mutate(day_mod = 'Day21_blue') %>% dplyr::mutate(type = 'biological process') %>% 
  arrange(factor(term))
# View(d21_blue_BP)
d21_blue_MF  <- read.csv("subseq_treatments_all/Day21/GOterms_and_GOslim_MolFunction_blueModule_GENE_REFERENCE.csv")
d21_blue_MF <- d21_blue_MF %>% 
  select(c('PgenIDs','slim_term','term','Gene_terms',)) %>% 
  dplyr::mutate(day_mod = 'Day21_blue') %>% dplyr::mutate(type = 'molecular function') %>% 
  arrange(factor(term))
# View(d21_blue_MF)


d21_black_BP  <- read.csv("subseq_treatments_all/Day21/GOterms_and_GOslim_BiolProc_blackModule_GENE_REFERENCE.csv")
d21_black_BP <- d21_black_BP %>% 
  select(c('PgenIDs','slim_term','term','Gene_terms',)) %>% 
  dplyr::mutate(day_mod = 'Day21_black') %>% dplyr::mutate(type = 'biological process') %>% 
  arrange(factor(term))
# View(d21_black_BP)
d21_black_MF  <- read.csv("subseq_treatments_all/Day21/GOterms_and_GOslim_MolFunction_blackModule_GENE_REFERENCE.csv")
d21_black_MF <- d21_black_MF %>% 
  select(c('PgenIDs','slim_term','term','Gene_terms',)) %>% 
  dplyr::mutate(day_mod = 'Day21_black') %>% dplyr::mutate(type = 'molecular function') %>% 
  arrange(factor(term))
# View(d21_black_MF)


d21_magenta_BP  <- read.csv("subseq_treatments_all/Day21/GOterms_and_GOslim_BiolProc_magentaModule_GENE_REFERENCE.csv")
d21_magenta_BP <- d21_magenta_BP %>% 
  select(c('PgenIDs','slim_term','term','Gene_terms',)) %>% 
  dplyr::mutate(day_mod = 'Day21_magenta') %>% dplyr::mutate(type = 'biological process') %>% 
  arrange(factor(term))
# View(d21_magenta_BP)
d21_magenta_MF  <- read.csv("subseq_treatments_all/Day21/GOterms_and_GOslim_MolFunction_magentaModule_GENE_REFERENCE.csv")
d21_magenta_MF <- d21_magenta_MF %>% 
  select(c('PgenIDs','slim_term','term','Gene_terms',)) %>% 
  dplyr::mutate(day_mod = 'Day21_magenta') %>% dplyr::mutate(type = 'molecular function') %>% 
  arrange(factor(term))
# View(d21_magenta_MF)


d21_pink_BP  <- read.csv("subseq_treatments_all/Day21/GOterms_and_GOslim_BiolProc_pinkModule_GENE_REFERENCE.csv")
d21_pink_BP <- d21_pink_BP %>% 
  select(c('PgenIDs','slim_term','term','Gene_terms',)) %>% 
  dplyr::mutate(day_mod = 'Day21_pink') %>% dplyr::mutate(type = 'biological process') %>% 
  arrange(factor(term))
# View(d21_pink_BP)
d21_pink_MF  <- read.csv("subseq_treatments_all/Day21/GOterms_and_GOslim_MolFunction_pinkModule_GENE_REFERENCE.csv")
d21_pink_MF <- d21_pink_MF %>% 
  select(c('PgenIDs','slim_term','term','Gene_terms',)) %>% 
  dplyr::mutate(day_mod = 'Day21_pink') %>% dplyr::mutate(type = 'molecular function') %>% 
  arrange(factor(term))
# View(d21_pink_MF)



d21_red_BP  <- read.csv("subseq_treatments_all/Day21/GOterms_and_GOslim_BiolProc_redModule_GENE_REFERENCE.csv")
d21_red_BP <- d21_red_BP %>% 
  select(c('PgenIDs','slim_term','term','Gene_terms',)) %>% 
  dplyr::mutate(day_mod = 'Day21_red') %>% dplyr::mutate(type = 'biological process') %>% 
  arrange(factor(term))
# View(d21_red_BP)
d21_red_MF  <- read.csv("subseq_treatments_all/Day21/GOterms_and_GOslim_MolFunction_redModule_GENE_REFERENCE.csv")
d21_red_MF <- d21_red_MF %>% 
  select(c('PgenIDs','slim_term','term','Gene_terms',)) %>% 
  dplyr::mutate(day_mod = 'Day21_red') %>% dplyr::mutate(type = 'molecular function') %>% 
  arrange(factor(term))
# View(d21_red_MF)


d21_turquoise_BP  <- read.csv("subseq_treatments_all/Day21/GOterms_and_GOslim_BiolProc_turquoiseModule_GENE_REFERENCE.csv")
d21_turquoise_BP <- d21_turquoise_BP %>% 
  select(c('PgenIDs','slim_term','term','Gene_terms',)) %>% 
  dplyr::mutate(day_mod = 'Day21_turquoise') %>% dplyr::mutate(type = 'biological process') %>% 
  arrange(factor(term))
# View(d21_turquoise_BP)
d21_turquoise_MF  <- read.csv("subseq_treatments_all/Day21/GOterms_and_GOslim_MolFunction_turquoiseModule_GENE_REFERENCE.csv")
d21_turquoise_MF <- d21_turquoise_MF %>% 
  select(c('PgenIDs','slim_term','term','Gene_terms',)) %>% 
  dplyr::mutate(day_mod = 'Day21_turquoise') %>% dplyr::mutate(type = 'molecular function') %>% 
  arrange(factor(term))
# View(d21_turquoise_MF)


GOSlim_Master_AllTreatments <- rbind(d0_midnightblue_BP, d0_midnightblue_MF,
                                 
                                 d7_brown_BP, d7_brown_MF, 
                                 d7_green_BP, d7_green_MF, 
                                 d7_yellow_BP, d7_yellow_MF,
                                 
                                 d14_black_BP, d14_black_MF,
                                 d14_brown_BP, d14_brown_MF, 
                                 d14_pink_BP, d14_pink_MF, 
                                 d14_magenta_BP, d14_magenta_MF,
                                 
                                 d21_black_BP, d21_black_MF,
                                 d21_blue_BP, d21_blue_MF, 
                                 d21_pink_BP, d21_pink_MF, 
                                 d21_magenta_BP, d21_magenta_MF,
                                 d21_red_BP, d21_red_MF,
                                 d21_turquoise_BP, d21_turquoise_MF)
View(GOSlim_Master_AllTreatments)

colnames(GOSlim_Master_AllTreatments) <- c("Pgenerosa_Gene_IDs", "GOslim_term", "GO_term", "Gene_name", "day_module", "GO_type")
write.csv(GOSlim_Master_AllTreatments, file = paste("subseq_treatments_all/All_GOSlim_Master.csv", sep ='')) 
