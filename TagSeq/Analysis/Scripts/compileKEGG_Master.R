# SET WORKING DIRECTORY   ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::;: #
setwd("C:/Users/samjg/Documents/Github_repositories/Pgenerosa_TagSeq_Metabolomics/TagSeq/")

library(dplyr)

# :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# WGCNA results (all treatments ) ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
d0_midnightblue_KEGG   <- read.csv("Analysis/Output/KEGG/subseq_treatments_all/WGCNA/Day0_midnightblue_KEGG_allgenes_unlisted.csv")
d0_midnightblue_KEGG <- d0_midnightblue_KEGG %>% 
  select(c('Cgigas_PathwayID','Pathway_Description','Cgigas_KEGG_IDs','Gene_name',)) %>% 
  dplyr::mutate(day_mod = 'Day0_midnightblue') %>% 
  arrange(factor(Pathway_Description))

d7_brown_KEGG   <- read.csv("Analysis/Output/KEGG/subseq_treatments_all/WGCNA/Day7_brown_KEGG_allgenes_unlisted.csv")
d7_brown_KEGG <- d7_brown_KEGG %>% 
  select(c('Cgigas_PathwayID','Pathway_Description','Cgigas_KEGG_IDs','Gene_name',)) %>% 
  dplyr::mutate(day_mod = 'Day7_brown') %>% 
  arrange(factor(Pathway_Description))

d7_yellow_KEGG   <- read.csv("Analysis/Output/KEGG/subseq_treatments_all/WGCNA/Day7_yellow_KEGG_allgenes_unlisted.csv")
d7_yellow_KEGG <- d7_yellow_KEGG %>% 
  select(c('Cgigas_PathwayID','Pathway_Description','Cgigas_KEGG_IDs','Gene_name',)) %>% 
  dplyr::mutate(day_mod = 'Day7_yellow') %>% 
  arrange(factor(Pathway_Description))


d14_black_KEGG   <- read.csv("Analysis/Output/KEGG/subseq_treatments_all/WGCNA/Day14_black_KEGG_allgenes_unlisted.csv")
d14_black_KEGG <- d14_black_KEGG %>% 
  select(c('Cgigas_PathwayID','Pathway_Description','Cgigas_KEGG_IDs','Gene_name',)) %>% 
  dplyr::mutate(day_mod = 'Day14_black') %>% 
  arrange(factor(Pathway_Description))


d14_brown_KEGG   <- read.csv("Analysis/Output/KEGG/subseq_treatments_all/WGCNA/Day14_brown_KEGG_allgenes_unlisted.csv")
d14_brown_KEGG <- d14_brown_KEGG %>% 
  select(c('Cgigas_PathwayID','Pathway_Description','Cgigas_KEGG_IDs','Gene_name',)) %>% 
  dplyr::mutate(day_mod = 'Day14_brown') %>% 
  arrange(factor(Pathway_Description))


d14_magenta_KEGG   <- read.csv("Analysis/Output/KEGG/subseq_treatments_all/WGCNA/Day14_magenta_KEGG_allgenes_unlisted.csv")
d14_magenta_KEGG <- d14_magenta_KEGG %>% 
  select(c('Cgigas_PathwayID','Pathway_Description','Cgigas_KEGG_IDs','Gene_name',)) %>% 
  dplyr::mutate(day_mod = 'Day14_magenta') %>% 
  arrange(factor(Pathway_Description))


d21_black_KEGG   <- read.csv("Analysis/Output/KEGG/subseq_treatments_all/WGCNA/Day21_black_KEGG_allgenes_unlisted.csv")
d21_black_KEGG <- d21_black_KEGG %>% 
  select(c('Cgigas_PathwayID','Pathway_Description','Cgigas_KEGG_IDs','Gene_name',)) %>% 
  dplyr::mutate(day_mod = 'Day21_black') %>% 
  arrange(factor(Pathway_Description))


d21_blue_KEGG   <- read.csv("Analysis/Output/KEGG/subseq_treatments_all/WGCNA/Day21_blue_KEGG_allgenes_unlisted.csv")
d21_blue_KEGG <- d21_blue_KEGG %>% 
  select(c('Cgigas_PathwayID','Pathway_Description','Cgigas_KEGG_IDs','Gene_name',)) %>% 
  dplyr::mutate(day_mod = 'Day21_blue') %>% 
  arrange(factor(Pathway_Description))


d21_yellow_KEGG   <- read.csv("Analysis/Output/KEGG/subseq_treatments_all/WGCNA/Day21_yellow_KEGG_allgenes_unlisted.csv")
d21_yellow_KEGG <- d21_yellow_KEGG %>% 
  select(c('Cgigas_PathwayID','Pathway_Description','Cgigas_KEGG_IDs','Gene_name',)) %>% 
  dplyr::mutate(day_mod = 'Day21_yellow') %>% 
  arrange(factor(Pathway_Description))


d21_magenta_KEGG   <- read.csv("Analysis/Output/KEGG/subseq_treatments_all/WGCNA/Day21_magenta_KEGG_allgenes_unlisted.csv")
d21_magenta_KEGG <- d21_magenta_KEGG %>% 
  select(c('Cgigas_PathwayID','Pathway_Description','Cgigas_KEGG_IDs','Gene_name',)) %>% 
  dplyr::mutate(day_mod = 'Day21_magenta') %>% 
  arrange(factor(Pathway_Description))


KEGG_Master_AllTreatments <- rbind(d0_midnightblue_KEGG,
                                   d7_brown_KEGG, d7_yellow_KEGG,
                                   d14_black_KEGG, d14_brown_KEGG, d14_magenta_KEGG,
                                   d21_black_KEGG, d21_yellow_KEGG, d21_blue_KEGG, d21_magenta_KEGG)
View(KEGG_Master_AllTreatments)

write.csv(KEGG_Master_AllTreatments, file = paste("Analysis/Output/KEGG/subseq_treatments_all/WGCNA/Master_KEGG_all_genes.csv", sep ='')) 


# :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# WGCNA results (subseqent exposures to low pH ONLY)::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
d7_brown_KEGG   <- read.csv("Analysis/Output/KEGG/subseq_hypercapnia/WGCNA/Day7_brown_KEGG_allgenes_unlisted.csv")
d7_brown_KEGG <- d7_brown_KEGG %>% 
  select(c('Cgigas_PathwayID','Pathway_Description','Cgigas_KEGG_IDs','Gene_name',)) %>% 
  dplyr::mutate(day_mod = 'Day7_brown') %>% 
  arrange(factor(Pathway_Description))

d7_greenyellow_KEGG   <- read.csv("Analysis/Output/KEGG/subseq_hypercapnia/WGCNA/Day7_greenyellow_KEGG_allgenes_unlisted.csv")
d7_greenyellow_KEGG <- d7_greenyellow_KEGG %>% 
  select(c('Cgigas_PathwayID','Pathway_Description','Cgigas_KEGG_IDs','Gene_name',)) %>% 
  dplyr::mutate(day_mod = 'Day7_greenyellow') %>% 
  arrange(factor(Pathway_Description))

d7_turquoise_KEGG   <- read.csv("Analysis/Output/KEGG/subseq_hypercapnia/WGCNA/Day7_turquoise_KEGG_allgenes_unlisted.csv")
d7_turquoise_KEGG <- d7_turquoise_KEGG %>% 
  select(c('Cgigas_PathwayID','Pathway_Description','Cgigas_KEGG_IDs','Gene_name',)) %>% 
  dplyr::mutate(day_mod = 'Day7_turquoise') %>% 
  arrange(factor(Pathway_Description))


d14_blue_KEGG   <- read.csv("Analysis/Output/KEGG/subseq_hypercapnia/WGCNA/Day14_blue_KEGG_allgenes_unlisted.csv")
d14_blue_KEGG <- d14_blue_KEGG %>% 
  select(c('Cgigas_PathwayID','Pathway_Description','Cgigas_KEGG_IDs','Gene_name',)) %>% 
  dplyr::mutate(day_mod = 'Day14_blue') %>% 
  arrange(factor(Pathway_Description))


d14_pink_KEGG   <- read.csv("Analysis/Output/KEGG/subseq_hypercapnia/WGCNA/Day14_pink_KEGG_allgenes_unlisted.csv")
d14_pink_KEGG <- d14_pink_KEGG %>% 
  select(c('Cgigas_PathwayID','Pathway_Description','Cgigas_KEGG_IDs','Gene_name',)) %>% 
  dplyr::mutate(day_mod = 'Day14_pink') %>% 
  arrange(factor(Pathway_Description))


d21_blue_KEGG   <- read.csv("Analysis/Output/KEGG/subseq_hypercapnia/WGCNA/Day21_blue_KEGG_allgenes_unlisted.csv")
d21_blue_KEGG <- d21_blue_KEGG %>% 
  select(c('Cgigas_PathwayID','Pathway_Description','Cgigas_KEGG_IDs','Gene_name',)) %>% 
  dplyr::mutate(day_mod = 'Day21_blue') %>% 
  arrange(factor(Pathway_Description))


d21_magenta_KEGG   <- read.csv("Analysis/Output/KEGG/subseq_hypercapnia/WGCNA/Day21_magenta_KEGG_allgenes_unlisted.csv")
d21_magenta_KEGG <- d21_magenta_KEGG %>% 
  select(c('Cgigas_PathwayID','Pathway_Description','Cgigas_KEGG_IDs','Gene_name',)) %>% 
  dplyr::mutate(day_mod = 'Day21_magenta') %>% 
  arrange(factor(Pathway_Description))


KEGG_Master_SubseqHyp <- rbind(d7_brown_KEGG, d7_greenyellow_KEGG, d7_turquoise_KEGG,
                               d14_blue_KEGG, d14_pink_KEGG,
                               d21_blue_KEGG, d21_magenta_KEGG)
View(KEGG_Master_SubseqHyp)

write.csv(KEGG_Master_SubseqHyp, file = paste("Analysis/Output/KEGG/subseq_hypercapnia/WGCNA/Master_KEGG_all_genes.csv", sep ='')) 

