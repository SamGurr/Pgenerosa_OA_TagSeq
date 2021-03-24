
# SET WORKING DIRECTORY AND LOAD DATA
setwd("C:/Users/samjg/Documents/My_Projects/Pgenerosa_TagSeq_Metabolomics/TagSeq/")

# Load libraries 
library(dplyr)
library(goseq)
library(reshape2)
library(ggplot2)
library(Rmisc)
library(ggpubr)
library(tibble)
library(hrbrthemes)
library(gridExtra)
library(ggplot2)
library(tidyr)
library(forcats) # for plotting later..
library(zoo)
library(ComplexHeatmap)
library(circlize)

# Load data 
d7_Annot_ModuleMembership <-  read.csv("Analysis/Output/WGCNA/Day7/d7.WGCNA_ModulMembership.csv")
d14_Annot_ModuleMembership <-  read.csv("Analysis/Output/WGCNA/Day14/d14.WGCNA_ModulMembership.csv")
d21_Annot_ModuleMembership <-  read.csv("Analysis/Output/WGCNA/Day21/d21.WGCNA_ModulMembership.csv")


#===================================================================================== 
# 
# goseq - load the annotation and prepare the four steps for goseq!
#
#===================================================================================== 

### Panopea generosa - load .fna ('Geoduck_annotation') and foramt GO terms ('Geoduck_GOterms') and vectors
Geoduck_annotation <- read.delim2(file="C:/Users/samjg/Documents/My_Projects/Pgenerosa_TagSeq_Metabolomics/TagSeq/Seq_details/Panopea-generosa-genes-annotations.txt", header=F)

# build annotation file to merge with the mean LFC tables
annot.condenced <- Geoduck_annotation[,c(1,3:9)]
annot.condenced$gene.length <- annot.condenced$V4 - annot.condenced$V3
annot.condenced <- annot.condenced[,-c(2,3)]
names(annot.condenced) <- c('Gene.ID', 'Uniprot', 'HGNC', 'fxn', 'Go.terms', 'Go.fxns','gene.length')


# Prepare dataframe(s) and vectors for goseq 
### (1) Format 'GO.term' for goseq from the P.generosa annotation .fna file 'Geoduck_annotation'
# GO terms data (ALL)
Geoduck_GOterms <- as.data.frame(Geoduck_annotation) %>% dplyr::select(c('V1','V8'))
colnames(Geoduck_GOterms)[1:2] <- c('transcript.ID', 'GO.terms') # call gene name and the GO terms - (Uniprot ID 'V5')
splitted <- strsplit(as.character(Geoduck_GOterms$GO.terms), ";") #slit into multiple GO ids
GO.terms <- data.frame(v1 = rep.int(Geoduck_GOterms$transcript.ID, sapply(splitted, length)), v2 = unlist(splitted)) #list all genes with each of their GO terms in a single row

### (2) Unique Genes - vector based on all unique mapped reads 
# Construct a named vector of all target genes for goseq
GO_unique.genes.all <- as.vector(unique(Geoduck_annotation$V1)) # call all unique genes for GO analysis (goseq)

### (3) Gene length 
# length vector  
GO_gene.length <- Geoduck_annotation %>% dplyr::mutate(length = V4-V3) %>%  dplyr::select(c("V1","length"))
names(GO_gene.length)[1] <- "Gene.ID"
#GO_gene.length_merge <- merge(GO_gene.length, GO_magenta_genes, by = "Gene.ID")
length_vector <- GO_gene.length$
  
#==============================================================================
#
#  PLOTTING GO TERMS  FOR.... PRIMARY TREATMENT EFFECT AMBIENT > MODERATE 
# 
# day 7  : module BROWN
# day 14 : module BROWN
# day 21 : module MAGENTA
# day 21 : module BLUE 
#==============================================================================
Day7_mod_color <- as.data.frame(unique(d7_Annot_ModuleMembership$moduleColor))
names(Day7_mod_color)[1] <- "color"

Day14_mod_color <- as.data.frame(unique(d14_Annot_ModuleMembership$moduleColor))
names(Day14_mod_color)[1] <- "color"

Day21_mod_color <- as.data.frame(unique(d21_Annot_ModuleMembership$moduleColor))
names(Day21_mod_color)[1] <- "color"


#======================================================================= #
# MODULES WITH PRIMARY EFFECT == AMBIENT > MODERATE (vst expression and positive correlation - view heatmap)
# Day 7 brown module 
d7_Mod.Brown <- d7_Annot_ModuleMembership %>% dplyr::filter(moduleColor %in% 'brown') # %>%  dplyr::select("geneSymbol")
d7_Mod.Brown_genes <- d7_Mod.Brown[1]
names(d7_Mod.Brown_genes)[1] <- "Gene.ID" # 162 genws in the green module 
d7_Mod.Brown_integer <- as.integer(GO_unique.genes.all%in%(d7_Mod.Brown_genes$Gene.ID)) # convert to integer with all unique genes
names(d7_Mod.Brown_integer)=GO_unique.genes.all # rename

# Day 14 brown module 
d14_Mod.Brown <- d14_Annot_ModuleMembership %>% dplyr::filter(moduleColor %in% 'brown') # %>%  dplyr::select("geneSymbol")
d14_Mod.Brown_genes <- d14_Mod.Brown[1]
names(d14_Mod.Brown_genes)[1] <- "Gene.ID" # 162 genws in the green module 
d14_Mod.Brown_integer <- as.integer(GO_unique.genes.all%in%(d14_Mod.Brown_genes$Gene.ID)) # convert to integer with all unique genes
names(d14_Mod.Brown_integer)=GO_unique.genes.all # rename

# Day 21 Magenta module 
d21_Mod.Magenta <- d21_Annot_ModuleMembership %>% dplyr::filter(moduleColor %in% 'magenta') # %>%  dplyr::select("geneSymbol")
d21_Mod.Magenta_genes <- d21_Mod.Magenta[1]
names(d21_Mod.Magenta_genes)[1] <- "Gene.ID" # 162 genws in the green module 
d21_Mod.Magenta_integer <- as.integer(GO_unique.genes.all%in%(d21_Mod.Magenta_genes$Gene.ID)) # convert to integer with all unique genes
names(d21_Mod.Magenta_integer)=GO_unique.genes.all # rename

# Day 21 blue module 
d21_Mod.blue <- d21_Annot_ModuleMembership %>% dplyr::filter(moduleColor %in% 'blue') # %>%  dplyr::select("geneSymbol")
d21_Mod.blue_genes <- d21_Mod.blue[1]
names(d21_Mod.blue_genes)[1] <- "Gene.ID" # 162 genws in the green module 
d21_Mod.blue_integer <- as.integer(GO_unique.genes.all%in%(d21_Mod.blue_genes$Gene.ID)) # convert to integer with all unique genes
names(d21_Mod.blue_integer)=GO_unique.genes.all # rename

#======================================================================= #
#Calculate Probability Weighting Function (using 'nullp')
d7_Mod.Brown_pwf <-nullp(d7_Mod.Brown_integer, GO_unique.genes.all, bias.data=length_vector) #weight vector by length of gene
d14_Mod.Brown_pwf <-nullp(d14_Mod.Brown_integer, GO_unique.genes.all, bias.data=length_vector) #weight vector by length of gene
d21_Mod.Magenta_pwf <-nullp(d21_Mod.Magenta_integer, GO_unique.genes.all, bias.data=length_vector) #weight vector by length of gene
d21_Mod.Blue_pwf <-nullp(d21_Mod.blue_integer, GO_unique.genes.all, bias.data=length_vector) #weight vector by length of gene

#======================================================================= #
# Run goseq
d7_Mod.Brown.goseq  <-goseq(d7_Mod.Brown_pwf, ID.vector, gene2cat=GO.terms, test.cats=c("GO:CC", "GO:BP", "GO:MF"), method="Wallenius", use_genes_without_cat=TRUE)
d14_Mod.Brown.goseq  <-goseq(d14_Mod.Brown_pwf, ID.vector, gene2cat=GO.terms, test.cats=c("GO:CC", "GO:BP", "GO:MF"), method="Wallenius", use_genes_without_cat=TRUE)
d21_Mod.Magenta.goseq  <-goseq(d21_Mod.Magenta_pwf, ID.vector, gene2cat=GO.terms, test.cats=c("GO:CC", "GO:BP", "GO:MF"), method="Wallenius", use_genes_without_cat=TRUE)
d21_Mod.Blue.goseq  <-goseq(d21_Mod.Blue_pwf, ID.vector, gene2cat=GO.terms, test.cats=c("GO:CC", "GO:BP", "GO:MF"), method="Wallenius", use_genes_without_cat=TRUE)

#======================================================================= #
# call enriched GO terms and plot 
# day 7 Module Brown
d7_Mod.Brown.GO.05.a<-d7_Mod.Brown.goseq$category[d7_Mod.Brown.goseq$over_represented_pvalue<.05] # change twice here
d7_Mod.Brown.GO.05<-data.frame(d7_Mod.Brown.GO.05.a)
colnames(d7_Mod.Brown.GO.05) <- c("category")
d7_Mod.Brown.GO.05 <- merge(d7_Mod.Brown.GO.05, d7_Mod.Brown.goseq, by="category") # change here
d7_Mod.Brown.GO.05 <- d7_Mod.Brown.GO.05[order(-d7_Mod.Brown.GO.05$numDEInCat),]
d7_Mod.Brown.GO.05$term <- as.factor(d7_Mod.Brown.GO.05$term)
head(d7_Mod.Brown.GO.05)

d7_Mod.Brown_MF <- d7_Mod.Brown.GO.05 %>% 
                        drop_na(ontology) %>% 
                        mutate(term = fct_reorder(term, (-log10(over_represented_pvalue)) ),
                               ontology = paste(ontology, 'd7_brown', sep = "_")) %>%
                        dplyr::filter(ontology %in% ('MF_d7_brown'))
d7_Mod.Brown_BP <- d7_Mod.Brown.GO.05 %>% 
                        drop_na(ontology) %>% 
                        mutate(term = fct_reorder(term, (-log10(over_represented_pvalue)) ),
                               ontology = paste(ontology, 'd7_brown', sep = "_")) %>%
                        dplyr::filter(ontology %in% ('BP_d7_brown'))


# day 14 Module Brown
d14_Mod.Brown.GO.05.a<-d14_Mod.Brown.goseq$category[d14_Mod.Brown.goseq$over_represented_pvalue<.05] # change twice here
d14_Mod.Brown.GO.05<-data.frame(d14_Mod.Brown.GO.05.a)
colnames(d14_Mod.Brown.GO.05) <- c("category")
d14_Mod.Brown.GO.05 <- merge(d14_Mod.Brown.GO.05, d14_Mod.Brown.goseq, by="category") # change here
d14_Mod.Brown.GO.05 <- d14_Mod.Brown.GO.05[order(-d14_Mod.Brown.GO.05$numDEInCat),]
d14_Mod.Brown.GO.05$term <- as.factor(d14_Mod.Brown.GO.05$term)
head(d14_Mod.Brown.GO.05)

d14_Mod.Brown_MF <- d14_Mod.Brown.GO.05 %>% 
                          drop_na(ontology) %>% 
                          mutate(term = fct_reorder(term, (-log10(over_represented_pvalue)) ),
                                 ontology = paste(ontology, 'd14_brown', sep = "_")) %>%
                          dplyr::filter(ontology %in% ('MF_d14_brown'))
d14_Mod.Brown_BP <- d14_Mod.Brown.GO.05 %>% 
                          drop_na(ontology) %>% 
                          mutate(term = fct_reorder(term, (-log10(over_represented_pvalue)) ),
                                 ontology = paste(ontology, 'd14_brown', sep = "_")) %>%
                          dplyr::filter(ontology %in% ('BP_d14_brown'))

# day 21 module magenta
d21_Mod.Magenta.GO.05.a<-d21_Mod.Magenta.goseq$category[d21_Mod.Magenta.goseq$over_represented_pvalue<.05] # change twice here
d21_Mod.Magenta.GO.05<-data.frame(d21_Mod.Magenta.GO.05.a)
colnames(d21_Mod.Magenta.GO.05) <- c("category")
d21_Mod.Magenta.GO.05 <- merge(d21_Mod.Magenta.GO.05, d21_Mod.Magenta.goseq, by="category") # change here
d21_Mod.Magenta.GO.05 <- d21_Mod.Magenta.GO.05[order(-d21_Mod.Magenta.GO.05$numDEInCat),]
d21_Mod.Magenta.GO.05$term <- as.factor(d21_Mod.Magenta.GO.05$term)
head(d21_Mod.Magenta.GO.05)

d21_Mod.Magenta_MF <- d21_Mod.Magenta.GO.05 %>% 
                          drop_na(ontology) %>% 
                          mutate(term = fct_reorder(term, (-log10(over_represented_pvalue)) ),
                                 ontology = paste(ontology, 'd21_blue.magenta', sep = "_")) %>%
                          dplyr::filter(ontology %in% ('MF_d21_blue.magenta'))
d21_Mod.Magenta_BP <- d21_Mod.Magenta.GO.05 %>% 
                          drop_na(ontology) %>% 
                          mutate(term = fct_reorder(term, (-log10(over_represented_pvalue)) ),
                                # ontology = paste(ontology, 'd21_magenta', sep = "_")) %>%
                          # dplyr::filter(ontology %in% ('BP_d21_magenta'))
                                ontology = paste(ontology, 'd21_blue.magenta', sep = "_")) %>%
                          dplyr::filter(ontology %in% ('BP_d21_blue.magenta'))

# day 21 module Blue
d21_Mod.Blue.GO.05.a<-d21_Mod.Blue.goseq$category[d21_Mod.Blue.goseq$over_represented_pvalue<.05] # change twice here
d21_Mod.Blue.GO.05<-data.frame(d21_Mod.Blue.GO.05.a)
colnames(d21_Mod.Blue.GO.05) <- c("category")
d21_Mod.Blue.GO.05 <- merge(d21_Mod.Blue.GO.05, d21_Mod.Blue.goseq, by="category") # change here
d21_Mod.Blue.GO.05 <- d21_Mod.Blue.GO.05[order(-d21_Mod.Blue.GO.05$numDEInCat),]
d21_Mod.Blue.GO.05$term <- as.factor(d21_Mod.Blue.GO.05$term)
head(d21_Mod.Blue.GO.05)

d21_Mod.Blue_MF <- d21_Mod.Blue.GO.05 %>% 
                          drop_na(ontology) %>% 
                          mutate(term = fct_reorder(term, (-log10(over_represented_pvalue)) ),
                                 ontology = paste(ontology, 'd21_blue.magenta', sep = "_")) %>%
                          dplyr::filter(ontology %in% ('MF_d21_blue.magenta'))
d21_Mod.Blue_BP <- d21_Mod.Blue.GO.05 %>% 
                          drop_na(ontology) %>% 
                          mutate(term = fct_reorder(term, (-log10(over_represented_pvalue)) ),
                                 ontology = paste(ontology, 'd21_blue.magenta', sep = "_")) %>%
                          dplyr::filter(ontology %in% ('BP_d21_blue.magenta'))

#======================================================================= #
# bind the rows of MF and BP in the clustered modules (those twith primary treatment effect in same pattern/directionality)

# Primary effect: WGCNA significant module showing Ambient > Moderate - bind together day 7 brown, day 14 brown 
# MF
MF_Amb_effect <- rbind(d7_Mod.Brown_MF, d14_Mod.Brown_MF, d21_Mod.Magenta_MF, d21_Mod.Blue_MF) # d21_Mod.Magenta_MF,
# BP
BP_Amb_effect <- rbind(d7_Mod.Brown_BP, d14_Mod.Brown_BP, d21_Mod.Magenta_BP, d21_Mod.Blue_BP)


#======================================================================= #
# PLOT

MF_tile_Amb_effect  <- MF_Amb_effect %>% 
                          mutate(term = fct_reorder(term, ontology)) %>%
                          ggplot(aes(factor(ontology,level = c('MF_d7_brown', 'MF_d14_brown', 'MF_d21_blue.magenta')), term, fill= factor(ontology), alpha = over_represented_pvalue)) + 
                          geom_tile(fill = '#00BFC4') +
                          scale_alpha_continuous(range = c(1, 0.3)) +
                          theme_classic2() + 
                          ggtitle("GO Mol Fxn: Priming effect (> Exp in Amb history)") +
                          xlab("GO.ontology_Sampling.Day") +
                          ylab('') +
                          theme(legend.position="none") 
BP_tile_Amb_effect <- BP_Amb_effect %>% 
                          mutate(term = fct_reorder(term, ontology)) %>%
                          ggplot(aes(factor(ontology,level = c('BP_d7_brown', 'BP_d14_brown', 'BP_d21_blue.magenta')), term, fill= factor(ontology), alpha = over_represented_pvalue)) + 
                          geom_tile(fill = '#F8766D') +
                          scale_alpha_continuous(range = c(1, 0.3)) +
                          theme_classic2() + 
                          ggtitle("GO Biol Proc: Priming effect (> Exp in Amb history)") +
                          xlab("GO.ontology_Sampling.Day") +
                          ylab('') +
                          theme(legend.position="none") 
# SAVE PLOT
pdf(paste("Analysis/Output/WGCNA/GO_SigEnrich_PrimaryAmbient.pdf", sep =''), width=20, height=15)
print(ggarrange(MF_tile_Amb_effect, BP_tile_Amb_effect,         
                plotlist = NULL,
                ncol = 2,
                nrow = 1,
                labels = NULL))
dev.off()  

#==============================================================================
#
#  PLOTTING GO TERMS  FOR.... PRIMARY TREATMENT EFFECT MODERATE > AMBIENT 
# day 7  : module YELLOW
# day 14 : module BLACK
# day 21 : module YELLOW
#==============================================================================

#======================================================================= #
# MODULES WITH PRIMARY EFFECT == AMBIENT > MODERATE (vst expression and positive correlation - view heatmap)
# Day 7 Yellow module 
d7_Mod.Yellow <- d7_Annot_ModuleMembership %>% dplyr::filter(moduleColor %in% 'yellow') # %>%  dplyr::select("geneSymbol")
d7_Mod.Yellow_genes <- d7_Mod.Yellow[1]
names(d7_Mod.Yellow_genes)[1] <- "Gene.ID" # 162 genws in the green module 
d7_Mod.Yellow_integer <- as.integer(GO_unique.genes.all%in%(d7_Mod.Yellow_genes$Gene.ID)) # convert to integer with all unique genes
names(d7_Mod.Yellow_integer)=GO_unique.genes.all # rename

# Day 14 Black module 
d14_Mod.Black <- d14_Annot_ModuleMembership %>% dplyr::filter(moduleColor %in% 'black') # %>%  dplyr::select("geneSymbol")
d14_Mod.Black_genes <- d14_Mod.Black[1]
names(d14_Mod.Black_genes)[1] <- "Gene.ID" # 162 genws in the green module 
d14_Mod.Black_integer <- as.integer(GO_unique.genes.all%in%(d14_Mod.Black_genes$Gene.ID)) # convert to integer with all unique genes
names(d14_Mod.Black_integer)=GO_unique.genes.all # rename

# Day 21 Yellow module 
d21_Mod.Yellow <- d21_Annot_ModuleMembership %>% dplyr::filter(moduleColor %in% 'yellow') # %>%  dplyr::select("geneSymbol")
d21_Mod.Yellow_genes <- d21_Mod.Yellow[1]
names(d21_Mod.Yellow_genes)[1] <- "Gene.ID" # 162 genws in the green module 
d21_Mod.Yellow_integer <- as.integer(GO_unique.genes.all%in%(d21_Mod.Yellow_genes$Gene.ID)) # convert to integer with all unique genes
names(d21_Mod.Yellow_integer)=GO_unique.genes.all # rename

#======================================================================= #
#Calculate Probability Weighting Function (using 'nullp')
d7_Mod.Yellow_pwf <-nullp(d7_Mod.Yellow_integer, GO_unique.genes.all, bias.data=length_vector) #weight vector by length of gene
d14_Mod.Black_pwf <-nullp(d14_Mod.Black_integer, GO_unique.genes.all, bias.data=length_vector) #weight vector by length of gene
d21_Mod.Yellow_pwf <-nullp(d21_Mod.Yellow_integer, GO_unique.genes.all, bias.data=length_vector) #weight vector by length of gene

#======================================================================= #
# Run goseq
d7_Mod.Yellow.goseq  <-goseq(d7_Mod.Yellow_pwf, ID.vector, gene2cat=GO.terms, test.cats=c("GO:CC", "GO:BP", "GO:MF"), method="Wallenius", use_genes_without_cat=TRUE)
d14_Mod.Black.goseq  <-goseq(d14_Mod.Black_pwf, ID.vector, gene2cat=GO.terms, test.cats=c("GO:CC", "GO:BP", "GO:MF"), method="Wallenius", use_genes_without_cat=TRUE)
d21_Mod.Yellow.goseq  <-goseq(d21_Mod.Yellow_pwf, ID.vector, gene2cat=GO.terms, test.cats=c("GO:CC", "GO:BP", "GO:MF"), method="Wallenius", use_genes_without_cat=TRUE)

#======================================================================= #
# call enriched GO terms and plot 
# day 7 Module Yellow
d7_Mod.Yellow.GO.05.a<-d7_Mod.Yellow.goseq$category[d7_Mod.Yellow.goseq$over_represented_pvalue<.05] # change twice here
d7_Mod.Yellow.GO.05<-data.frame(d7_Mod.Yellow.GO.05.a)
colnames(d7_Mod.Yellow.GO.05) <- c("category")
d7_Mod.Yellow.GO.05 <- merge(d7_Mod.Yellow.GO.05, d7_Mod.Yellow.goseq, by="category") # change here
d7_Mod.Yellow.GO.05 <- d7_Mod.Yellow.GO.05[order(-d7_Mod.Yellow.GO.05$numDEInCat),]
d7_Mod.Yellow.GO.05$term <- as.factor(d7_Mod.Yellow.GO.05$term)
head(d7_Mod.Yellow.GO.05)

d7_Mod.Yellow_MF <- d7_Mod.Yellow.GO.05 %>% 
  drop_na(ontology) %>% 
  mutate(term = fct_reorder(term, (-log10(over_represented_pvalue)) ),
         ontology = paste(ontology, 'd7_Yellow', sep = "_")) %>%
  dplyr::filter(ontology %in% ('MF_d7_Yellow'))
d7_Mod.Yellow_BP <- d7_Mod.Yellow.GO.05 %>% 
  drop_na(ontology) %>% 
  mutate(term = fct_reorder(term, (-log10(over_represented_pvalue)) ),
         ontology = paste(ontology, 'd7_Yellow', sep = "_")) %>%
  dplyr::filter(ontology %in% ('BP_d7_Yellow'))


# day 14 Module Black
d14_Mod.Black.GO.05.a<-d14_Mod.Black.goseq$category[d14_Mod.Black.goseq$over_represented_pvalue<.05] # change twice here
d14_Mod.Black.GO.05<-data.frame(d14_Mod.Black.GO.05.a)
colnames(d14_Mod.Black.GO.05) <- c("category")
d14_Mod.Black.GO.05 <- merge(d14_Mod.Black.GO.05, d14_Mod.Black.goseq, by="category") # change here
d14_Mod.Black.GO.05 <- d14_Mod.Black.GO.05[order(-d14_Mod.Black.GO.05$numDEInCat),]
d14_Mod.Black.GO.05$term <- as.factor(d14_Mod.Black.GO.05$term)
head(d14_Mod.Black.GO.05)

d14_Mod.Black_MF <- d14_Mod.Black.GO.05 %>% 
  drop_na(ontology) %>% 
  mutate(term = fct_reorder(term, (-log10(over_represented_pvalue)) ),
         ontology = paste(ontology, 'd14_Black', sep = "_")) %>%
  dplyr::filter(ontology %in% ('MF_d14_Black'))
d14_Mod.Black_BP <- d14_Mod.Black.GO.05 %>% 
  drop_na(ontology) %>% 
  mutate(term = fct_reorder(term, (-log10(over_represented_pvalue)) ),
         ontology = paste(ontology, 'd14_Black', sep = "_")) %>%
  dplyr::filter(ontology %in% ('BP_d14_Black'))

# day 21 module Yellow
d21_Mod.Yellow.GO.05.a<-d21_Mod.Yellow.goseq$category[d21_Mod.Yellow.goseq$over_represented_pvalue<.05] # change twice here
d21_Mod.Yellow.GO.05<-data.frame(d21_Mod.Yellow.GO.05.a)
colnames(d21_Mod.Yellow.GO.05) <- c("category")
d21_Mod.Yellow.GO.05 <- merge(d21_Mod.Yellow.GO.05, d21_Mod.Yellow.goseq, by="category") # change here
d21_Mod.Yellow.GO.05 <- d21_Mod.Yellow.GO.05[order(-d21_Mod.Yellow.GO.05$numDEInCat),]
d21_Mod.Yellow.GO.05$term <- as.factor(d21_Mod.Yellow.GO.05$term)
head(d21_Mod.Yellow.GO.05)

d21_Mod.Yellow_MF <- d21_Mod.Yellow.GO.05 %>% 
  drop_na(ontology) %>% 
  mutate(term = fct_reorder(term, (-log10(over_represented_pvalue)) ),
         ontology = paste(ontology, 'd21_Yellow', sep = "_")) %>%
  dplyr::filter(ontology %in% ('MF_d21_Yellow'))
d21_Mod.Yellow_BP <- d21_Mod.Yellow.GO.05 %>% 
  drop_na(ontology) %>% 
  mutate(term = fct_reorder(term, (-log10(over_represented_pvalue)) ),
         ontology = paste(ontology, 'd21_Yellow', sep = "_")) %>%
  dplyr::filter(ontology %in% ('BP_d21_Yellow'))


#======================================================================= #
# bind the rows of MF and BP in the clustered modules (those twith primary treatment effect in same pattern/directionality)

# Primary effect: WGCNA significant module showing Ambient > Moderate - bind together day 7 brown, day 14 brown 
# MF
MF_Mod_effect <- rbind(d7_Mod.Yellow_MF, d14_Mod.Black_MF, d21_Mod.Yellow_MF) # d21_Mod.Magenta_MF,
# BP
BP_Mod_effect <- rbind(d7_Mod.Yellow_BP, d14_Mod.Black_BP, d21_Mod.Yellow_BP)


#======================================================================= #
# PLOT

MF_tile_Mod_effect  <- MF_Mod_effect %>% 
                          mutate(term = fct_reorder(term, ontology)) %>%
                          ggplot(aes(factor(ontology,level = c('MF_d7_Yellow', 'MF_d14_Black', 'MF_d21_Yellow')), term, fill= factor(ontology), alpha = over_represented_pvalue)) + 
                          geom_tile(fill = '#00BFC4') +
                          scale_alpha_continuous(range = c(1, 0.3)) +
                          theme_classic2() + 
                              ggtitle("GO Mol Fxn: Priming effect (> Exp in Moderate history)") +
                              xlab("GO.ontology_Sampling.Day") +
                              ylab('') +
                              theme(legend.position="none") 
BP_tile_Mod_effect <- BP_Mod_effect %>% 
                          mutate(term = fct_reorder(term, ontology)) %>%
                          ggplot(aes(factor(ontology,level = c('BP_d7_Yellow', 'BP_d14_Black', 'BP_d21_Yellow')), term, fill= factor(ontology), alpha = over_represented_pvalue)) + 
                          geom_tile(fill = '#F8766D') +
                          scale_alpha_continuous(range = c(1, 0.3)) +
                          theme_classic2() + 
                              ggtitle("GO Biol Proc: Priming effect (> Exp in Moderate history)") +
                              xlab("GO.ontology_Sampling.Day") +
                              ylab('') +
                              theme(legend.position="none") 
# SAVE PLOT
pdf(paste("Analysis/Output/WGCNA/GO_SigEnrich_PrimaryModerate.pdf", sep =''), width=20, height=15)
print(ggarrange(MF_tile_Mod_effect, BP_tile_Mod_effect,         
                plotlist = NULL,
                ncol = 2,
                nrow = 1,
                labels = NULL))
dev.off() 
