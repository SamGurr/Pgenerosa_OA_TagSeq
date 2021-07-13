---
  # title: "GO_Analysis_WGCNA_all"
  # author: "Samuel Gurr"
  # date: "April 12, 2021"
---
  

# SET WORKING DIRECTORY AND LOAD DATA
setwd("C:/Users/samjg/Documents/My_Projects/Pgenerosa_TagSeq_Metabolomics/TagSeq/")

# Load libraries 
library(dplyr)
library(tidyr)
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
library(GO.db)
library(GSEABase)
library(data.table) 
library(stringr)

# =================================================================================================
# LOAD DATA -  WGCNA Module membership +  Module colors + goslim_generic.obo 
#
# =================================================================================================
# Build a master list of all genes and GO terms 'Pgen_GOterms2'
Geoduck_annotation <- read.delim2(file="C:/Users/samjg/Documents/My_Projects/Pgenerosa_TagSeq_Metabolomics/TagSeq/Seq_details/Panopea-generosa-genes-annotations.txt", header=F) # Load themaster Pgenerosa gene list 
Pgen_GOterms <- Geoduck_annotation %>% dplyr::select(c('V1','V8')) # select only two columns - those with the gene IDs and those with the GO terms
Pgen_GOterms2 <- strsplit(Pgen_GOterms$V8, split = "; ") # create a string splitting by delimiter '; ' - view the data to see that this separates each GO term entry in the string
Pgen_GOterms2 <- data.frame(gene.ID = rep(Pgen_GOterms$V1, sapply(Pgen_GOterms2, length)), Go.terms = unlist(Pgen_GOterms2)) # create new dataframe 'Pgen_GOterms2' listing genes for each GO term (MUCH longer!)
Pgen_GOterms2 <- na.omit(Pgen_GOterms2) # ommit the NAs  - genes without GO annotation


# WGCNa data results 
d0_Annot_ModuleMembership      <-  read.csv("Analysis/Output/WGCNA/subseq_treatments_all/Day0/d0.WGCNA_ModulMembership.csv")   # WGCNA results day 7  - Module membership 
d0_Annot_ModuleMembership      <- d0_Annot_ModuleMembership[,c(3:4,7:8)] # call gene ID, module color, and GO annotation
d0_Annot_ModuleMembership$Day  <- "Day0"  # common column to divide master dataset
d0ModCols                      <- data.frame(moduleColor = unique(d0_Annot_ModuleMembership$moduleColor)) # call all unique module colors 
d0ModCols                      <- d0ModCols %>% filter(moduleColor %in% 'midnightblue') # MODULES WITH SIG CORR WITH TREATMENT
d0ModCols$Day                  <- "Day0" # common column for the for loop

d7_Annot_ModuleMembership      <-  read.csv("Analysis/Output/WGCNA/subseq_treatments_all/Day7/d7.WGCNA_ModulMembership.csv")   # WGCNA results day 7  - Module membership 
d7_Annot_ModuleMembership      <- d7_Annot_ModuleMembership[,c(3:4,7:8)] # call gene ID, module color, and GO annotation
d7_Annot_ModuleMembership$Day  <- "Day7"  # common column to divide master dataset
d7ModCols                      <- data.frame(moduleColor = unique(d7_Annot_ModuleMembership$moduleColor)) # call all unique module colors 
d7ModCols                      <- d7ModCols %>% filter(moduleColor %in% c('brown', 'yellow', 'green')) # MODULES WITH SIG CORR WITH TREATMENT
d7ModCols$Day                  <- "Day7" # common column for the for loop

d14_Annot_ModuleMembership     <-  read.csv("Analysis/Output/WGCNA/subseq_treatments_all/Day14/d14.WGCNA_ModulMembership.csv") # WGCNA results day 14 - Module membership 
d14_Annot_ModuleMembership     <- d14_Annot_ModuleMembership[,c(3:4,7:8)] # call gene ID, module color, and GO annotation
d14_Annot_ModuleMembership$Day <- "Day14"  # common column to divide master dataset
d14ModCols                     <- data.frame(moduleColor = unique(d14_Annot_ModuleMembership$moduleColor)) # call all unique module colors 
d14ModCols                     <- d14ModCols %>% filter(moduleColor %in% c('brown', 'black', 'pink', 'magenta')) # MODULES WITH SIG CORR WITH TREATMENT
d14ModCols$Day                 <- "Day14" # common column for the for loop

d21_Annot_ModuleMembership     <-  read.csv("Analysis/Output/WGCNA/subseq_treatments_all/Day21/d21.WGCNA_ModulMembership.csv") # WGCNA results day 21 - Module membership 
d21_Annot_ModuleMembership     <- d21_Annot_ModuleMembership[,c(3:4,7:8)] # call gene ID, module color, and GO annotation
d21_Annot_ModuleMembership$Day <- "Day21" # common column to divide master dataset
d21ModCols                     <- data.frame(moduleColor = unique(d21_Annot_ModuleMembership$moduleColor)) # call all unique module colors 
d21ModCols                     <- d21ModCols %>% filter(moduleColor %in% c('magenta', 'blue', 'yellow', 'red', 'black', 'pink', 'turquoise')) # MODULES WITH SIG CORR WITH TREATMENT
d21ModCols$Day                 <- "Day21" # common column for the for loop
 
WGCNA_MasterModData   <-  rbind(d0_Annot_ModuleMembership, d7_Annot_ModuleMembership, d14_Annot_ModuleMembership, d21_Annot_ModuleMembership) # master WGCNA data table 
WGCNA_ColorList       <-  rbind(d0ModCols, d7ModCols, d14ModCols, d21ModCols) # master WGCNA color list - use this to loop all the analysis 

slim <- getOBOCollection("http://current.geneontology.org/ontology/subsets/goslim_generic.obo") #get GO database - # call goslim_generic.obo terms as 'slim'

#===================================================================================== # 
# LOAD DATA - Raw Count Data - FILTERED COUNT MATRICES USED IN WGCNA- 10CPM in 50% of samples
#
#===================================================================================== # 
# day7 filtered 10cpm in 50% samples ----------------------------- # 
Day0_all.counts <- read.csv(file="C:/Users/samjg/Documents/My_Projects/Pgenerosa_TagSeq_Metabolomics/TagSeq/Analysis/Data/Filtered_Counts/10cpm_50perc/day0.counts.filtered_10cpm50perc.csv", sep=',', header=TRUE) 
colnames(Day0_all.counts)[1] <- "gene.ID"# rename Pgen gene ID column

# day7 filtered 10cpm in 50% samples ----------------------------- # 
Day7_all.counts <- read.csv(file="C:/Users/samjg/Documents/My_Projects/Pgenerosa_TagSeq_Metabolomics/TagSeq/Analysis/Data/Filtered_Counts/10cpm_50perc/day7.counts.filtered_10cpm50perc.csv", sep=',', header=TRUE) 
colnames(Day7_all.counts)[1] <- "gene.ID"# rename Pgen gene ID column

# day14 filtered 10cpm in 50% samples ----------------------------- # 
Day14_all.counts <- read.csv(file="C:/Users/samjg/Documents/My_Projects/Pgenerosa_TagSeq_Metabolomics/TagSeq/Analysis/Data/Filtered_Counts/10cpm_50perc/day14.counts.filtered_10cpm50perc.csv", sep=',', header=TRUE) 
colnames(Day14_all.counts)[1] <- "gene.ID"# rename Pgen gene ID column

# day21 filtered 10cpm in 50% samples ----------------------------- # 
Day21_all.counts <- read.csv(file="C:/Users/samjg/Documents/My_Projects/Pgenerosa_TagSeq_Metabolomics/TagSeq/Analysis/Data/Filtered_Counts/10cpm_50perc/day21.counts.filtered_10cpm50perc.csv", sep=',', header=TRUE) 
colnames(Day21_all.counts)[1] <- "gene.ID"# rename Pgen gene ID column


#===================================================================================== # 
# LOAD DATA - goseq; load the annotation and prepare the fouressentail steps for goseq
#
#===================================================================================== #
#Panopea generosa - load .fna ('Geoduck_annotation') and foramt GO terms ('Geoduck_GOterms') and vectors
Geoduck_annotation <- read.delim2(file="C:/Users/samjg/Documents/My_Projects/Pgenerosa_TagSeq_Metabolomics/TagSeq/Seq_details/Panopea-generosa-genes-annotations.txt", header=F)

# build annotation file to merge with the mean LFC tables
annot.condenced <- Geoduck_annotation[,c(1,3:9)]
annot.condenced$gene.length <- annot.condenced$V4 - annot.condenced$V3
annot.condenced <- annot.condenced[,-c(2,3)]
names(annot.condenced) <- c('Gene.ID', 'Uniprot', 'HGNC', 'fxn', 'Go.terms', 'Go.fxns','gene.length')

# Prepare dataframe(s) and vectors for goseq 
# (1) Format 'GO.term' for goseq from the P.generosa annotation .fna file 'Geoduck_annotation'
Geoduck_GOterms <- as.data.frame(Geoduck_annotation) %>% dplyr::select(c('V1','V8'))
colnames(Geoduck_GOterms)[1:2] <- c('transcript.ID', 'GO.terms') # call gene name and the GO terms - (Uniprot ID 'V5')
splitted <- strsplit(as.character(Geoduck_GOterms$GO.terms), "; ") #slit into multiple GO ids
GO.terms <- data.frame(v1 = rep.int(Geoduck_GOterms$transcript.ID, sapply(splitted, length)), v2 = unlist(splitted)) #list all genes with each of their GO terms in a single row
GO.terms$Term <- Term(GO.terms$v2)
GO.terms$Ontology  <- Ontology(GO.terms$v2)

# (2) Unique Genes - vector based on all unique mapped reads 
# Construct a named vector of all target genes for goseq
GO_unique.genes.all <- as.vector(unique(Geoduck_annotation$V1)) # call all unique genes for GO analysis (goseq)
IDvector.d0         <- as.vector(unique(Day0_all.counts$gene.ID))  # call unique genes (those filtered and used in DESEq2) on day0 - 'IDvector'
IDvector.d7         <- as.vector(unique(Day7_all.counts$gene.ID))  # call unique genes (those filtered and used in DESEq2) on day7 - 'IDvector'
IDvector.d14        <- as.vector(unique(Day14_all.counts$gene.ID)) # call unique genes (those filtered and used in DESEq2) on day14 - 'IDvector'
IDvector.d21        <- as.vector(unique(Day21_all.counts$gene.ID)) # call unique genes (those filtered and used in DESEq2) on day21 - 'IDvector'

# (3) Gene length 
# length vector  
GO_gene.length <- Geoduck_annotation %>% dplyr::mutate(length = V4-V3) %>%  dplyr::select(c("V1","length"))
names(GO_gene.length)[1] <- "gene.ID"
# merge length with counts data
length_vector   <- GO_gene.length$length
GeneLength.d0   <- merge(GO_gene.length, Day0_all.counts, by = "gene.ID")  # merge day0 counts with 'GO_gene.length' 
GeneLength.d7   <- merge(GO_gene.length, Day7_all.counts, by = "gene.ID")  # merge day7 counts with 'GO_gene.length' 
GeneLength.d14  <- merge(GO_gene.length, Day14_all.counts, by = "gene.ID")  # merge day14 counts with 'GO_gene.length' 
GeneLength.d21  <- merge(GO_gene.length, Day21_all.counts, by = "gene.ID")  # merge day21 counts with 'GO_gene.length'
# call length values for goseq - confirms that the IDvector and length_vector are the same!!!
length_vector.d0 <- GeneLength.d0$length    # length vector for all unique reads address in WGCNA on day 0
sum(sapply(length_vector.d0,length)) == dim(Day0_all.counts)[1] #should be TRUE
length_vector.d7 <- GeneLength.d7$length    # length vector for all unique reads address in WGCNA on day 7
sum(sapply(length_vector.d7,length)) == dim(Day7_all.counts)[1] #should be TRUE
length_vector.d14 <- GeneLength.d14$length  # length vector for all unique reads address in WGCNA on day 14
sum(sapply(length_vector.d14,length)) == dim(Day14_all.counts)[1] #should be TRUE
length_vector.d21 <- GeneLength.d21$length  # length vector for all unique reads address in WGCNA on day 21
sum(sapply(length_vector.d21,length)) == dim(Day21_all.counts)[1] #should be TRUE 





# =================================================================================================
# FOR LOOP goseq
# for loop for all goseq analysis - if/else to loop through the day 7, 14, and 21 WGCNA modules separately 
# Objectives: 
# - run goseq - go enrichment analysis calling just those with adj P value < 0.05
# - run GOslim
# - plot all modules (within each sampling day) with heatmaos fshowing the number of genes within each GOslim
#
# =================================================================================================




for (i in 1:nrow(WGCNA_ColorList)) {
        if (WGCNA_ColorList[i,2] == "Day0") {
        Mod <- d0_Annot_ModuleMembership %>% dplyr::filter(moduleColor %in% WGCNA_ColorList[i,1]) # call the WGCNA Day - essential here!
        Modgenes <- Mod[1]
        names(Modgenes)[1] <- "Gene.ID" # 162 genws in the green module 
        # Mod_integer <- as.integer(IDvector.d7 %in% (Modgenes$Gene.ID)) # call the day-specific ID vector 
        # names(Mod_integer)=IDvector.d7 # rename
        Mod_integer <- as.integer(GO_unique.genes.all %in% (Modgenes$Gene.ID)) # w/o day-specific ID vector
        names(Mod_integer)=GO_unique.genes.all # rename
        
        #pwf       <- nullp(Mod_integer,    id=IDvector.d7, bias.data=length_vector.d7) # make figure margins large enough for this to run...
        pwf       <- nullp(Mod_integer,    id=GO_unique.genes.all, bias.data=length_vector) # make figure margins large enough for this to run...

        goseq     <- goseq(pwf, gene2cat=GO.terms, test.cats=c("GO:CC", "GO:BP", "GO:MF"), method="Wallenius", use_genes_without_cat=TRUE)
        
        GO.05.a         <- goseq$category[goseq$over_represented_pvalue<.05] # change twice here
        GO.05           <- data.frame(GO.05.a)
        colnames(GO.05) <- c("category")
        GO.05           <- merge(GO.05, goseq, by="category") # change here
        GO.05           <- GO.05[order(GO.05$ontology, GO.05$over_represented_pvalue,-GO.05$numDEInCat),]
        GO.05$term      <- as.factor(GO.05$term)
        GO.05$moduleColor <- WGCNA_ColorList[i,1]
        GO.05$Day       <- "Day0"
        
        # remove Biological Process GO terms with < 10 genes in the module  (with that term) and ommit Molecular Function terms with < 3 genes in the module (with that term)
        GO.05_filtered <- GO.05 %>% filter(!(numDEInCat<10 & ontology == "BP"), !(numDEInCat<2 & ontology == "MF"))
        
        write.csv(GO.05_filtered, file = paste("Analysis/Output/GO/WGCNA_goseq/subseq_treatments_all/Day0/GO.05",WGCNA_ColorList[i,1], "Module.csv", sep ='')) # save csv file
        
        print(paste(WGCNA_ColorList[i,2], WGCNA_ColorList[i,1], "done", sep = ' '))
        
        
        } else if (WGCNA_ColorList[i,2] == "Day7") {
          Mod <- d7_Annot_ModuleMembership %>% dplyr::filter(moduleColor %in% WGCNA_ColorList[i,1]) # call the WGCNA Day - essential here!
          Modgenes <- Mod[1]
          names(Modgenes)[1] <- "Gene.ID" # 162 genws in the green module 
          # Mod_integer <- as.integer(IDvector.d14 %in% (Modgenes$Gene.ID)) # call the day-specific ID vector 
          # names(Mod_integer)=IDvector.d14 # rename 
          Mod_integer <- as.integer(GO_unique.genes.all %in% (Modgenes$Gene.ID)) # w/o day-specific ID vector
          names(Mod_integer)=GO_unique.genes.all # rename
          
          #pwf       <- nullp(Mod_integer,    id=IDvector.d14, bias.data=length_vector.d14) # make figure margins large enough for this to run...
          pwf       <- nullp(Mod_integer,    id=GO_unique.genes.all, bias.data=length_vector) # make figure margins large enough for this to run...
          
          goseq     <- goseq(pwf, gene2cat=GO.terms, test.cats=c("GO:CC", "GO:BP", "GO:MF"), method="Wallenius", use_genes_without_cat=TRUE)
          
          GO.05.a         <- goseq$category[goseq$over_represented_pvalue<.05] # change twice here
          GO.05           <- data.frame(GO.05.a)
          colnames(GO.05) <- c("category")
          GO.05           <- merge(GO.05, goseq, by="category") # change here
          GO.05           <- GO.05[order(GO.05$ontology, GO.05$over_represented_pvalue,-GO.05$numDEInCat),]
          GO.05$term      <- as.factor(GO.05$term)
          GO.05$moduleColor <- WGCNA_ColorList[i,1]
          GO.05$Day       <- "Day7"
          
          # remove Biological Process GO terms with < 10 genes in the module  (with that term) and ommit Molecular Function terms with < 3 genes in the module (with that term)
          GO.05_filtered <- GO.05 %>% filter(!(numDEInCat<10 & ontology == "BP"), !(numDEInCat<2 & ontology == "MF"))
          
          write.csv(GO.05_filtered, file = paste("Analysis/Output/GO/WGCNA_goseq/subseq_treatments_all/Day7/GO.05",WGCNA_ColorList[i,1], "Module.csv", sep ='')) # save csv file
          
          print(paste(WGCNA_ColorList[i,2], WGCNA_ColorList[i,1], "done", sep = ' '))
          
          
             } else if (WGCNA_ColorList[i,2] == "Day14") {
                Mod <- d14_Annot_ModuleMembership %>% dplyr::filter(moduleColor %in% WGCNA_ColorList[i,1]) # call the WGCNA Day - essential here!
                Modgenes <- Mod[1]
                names(Modgenes)[1] <- "Gene.ID" # 162 genws in the green module 
                # Mod_integer <- as.integer(IDvector.d14 %in% (Modgenes$Gene.ID)) # call the day-specific ID vector 
                # names(Mod_integer)=IDvector.d14 # rename 
                Mod_integer <- as.integer(GO_unique.genes.all %in% (Modgenes$Gene.ID)) # w/o day-specific ID vector
                names(Mod_integer)=GO_unique.genes.all # rename
                
                #pwf       <- nullp(Mod_integer,    id=IDvector.d14, bias.data=length_vector.d14) # make figure margins large enough for this to run...
                pwf       <- nullp(Mod_integer,    id=GO_unique.genes.all, bias.data=length_vector) # make figure margins large enough for this to run...
                
                goseq     <- goseq(pwf, gene2cat=GO.terms, test.cats=c("GO:CC", "GO:BP", "GO:MF"), method="Wallenius", use_genes_without_cat=TRUE)
                
                GO.05.a         <- goseq$category[goseq$over_represented_pvalue<.05] # change twice here
                GO.05           <- data.frame(GO.05.a)
                colnames(GO.05) <- c("category")
                GO.05           <- merge(GO.05, goseq, by="category") # change here
                GO.05           <- GO.05[order(GO.05$ontology, GO.05$over_represented_pvalue,-GO.05$numDEInCat),]
                GO.05$term      <- as.factor(GO.05$term)
                GO.05$moduleColor <- WGCNA_ColorList[i,1]
                GO.05$Day       <- "Day14"
                
                # remove Biological Process GO terms with < 10 genes in the module  (with that term) and ommit Molecular Function terms with < 3 genes in the module (with that term)
                GO.05_filtered <- GO.05 %>% filter(!(numDEInCat<10 & ontology == "BP"), !(numDEInCat<2 & ontology == "MF"))
                
                write.csv(GO.05_filtered, file = paste("Analysis/Output/GO/WGCNA_goseq/subseq_treatments_all/Day14/GO.05",WGCNA_ColorList[i,1], "Module.csv", sep ='')) # save csv file
                
                print(paste(WGCNA_ColorList[i,2], WGCNA_ColorList[i,1], "done", sep = ' '))
                
                
                } else {
                  Mod <- d21_Annot_ModuleMembership %>% dplyr::filter(moduleColor %in% WGCNA_ColorList[i,1]) # call the WGCNA Day - essential here!
                  Modgenes <- Mod[1]
                  names(Modgenes)[1] <- "Gene.ID" # 162 genws in the green module 
                  # Mod_integer <- as.integer(IDvector.d21 %in% (Modgenes$Gene.ID)) # call the day-specific ID vector 
                  # names(Mod_integer)=IDvector.d21 # rename
                  Mod_integer <- as.integer(GO_unique.genes.all %in% (Modgenes$Gene.ID)) # w/o day-specific ID vector
                  names(Mod_integer)=GO_unique.genes.all # rename
                  
                  #pwf       <- nullp(Mod_integer,    id=IDvector.d21, bias.data=length_vector.d21) # make figure margins large enough for this to run...
                  pwf       <- nullp(Mod_integer, id=GO_unique.genes.all, bias.data=length_vector) # make figure margins large enough for this to run...
                  
                  goseq     <- goseq(pwf, gene2cat=GO.terms, test.cats=c("GO:CC", "GO:BP", "GO:MF"), method="Wallenius", use_genes_without_cat=TRUE)
                  
                  GO.05.a         <- goseq$category[goseq$over_represented_pvalue<.05] # change twice here
                  GO.05           <- data.frame(GO.05.a)
                  colnames(GO.05) <- c("category")
                  GO.05           <- merge(GO.05, goseq, by="category") # change here
                  GO.05           <- GO.05[order(GO.05$ontology, GO.05$over_represented_pvalue,-GO.05$numDEInCat),]
                  GO.05$term      <- as.factor(GO.05$term)
                  GO.05$moduleColor <- WGCNA_ColorList[i,1]
                  GO.05$Day       <- "Day21"
                  
                  # remove Biological Process GO terms with < 10 genes in the module  (with that term) and ommit Molecular Function terms with < 3 genes in the module (with that term)
                  GO.05_filtered <- GO.05 %>% filter(!(numDEInCat<10 & ontology == "BP"), !(numDEInCat<2 & ontology == "MF"))
    
                  write.csv(GO.05_filtered, file = paste("Analysis/Output/GO/WGCNA_goseq/subseq_treatments_all/Day21/GO.05",WGCNA_ColorList[i,1], "Module.csv", sep ='')) # save csv file              
                  
                  print(paste(WGCNA_ColorList[i,2], WGCNA_ColorList[i,1], "done", sep = ' '))
                  
                } #of if statement
    
} # end of for loop

#===================================================================================================
#
#
#  GOslim analysis 
# - data reduction from the many GO terms to more broad functions/processes
#
#
#===================================================================================================
# load the output from the previous for loop
d0_GO.05midnightbluenModule   <- read.csv("Analysis/Output/GO/WGCNA_goseq/subseq_treatments_all/Day0/GO.05midnightblueModule.csv")

d7_GO.05brownModule   <- read.csv("Analysis/Output/GO/WGCNA_goseq/subseq_treatments_all/Day7/GO.05brownModule.csv")
d7_GO.05greenModule   <- read.csv("Analysis/Output/GO/WGCNA_goseq/subseq_treatments_all/Day7/GO.05greenModule.csv")
d7_GO.05yellowModule  <- read.csv("Analysis/Output/GO/WGCNA_goseq/subseq_treatments_all/Day7/GO.05yellowModule.csv")

d14_GO.05blackModule   <- read.csv("Analysis/Output/GO/WGCNA_goseq/subseq_treatments_all/Day14/GO.05blackModule.csv")
d14_GO.05brownModule   <- read.csv("Analysis/Output/GO/WGCNA_goseq/subseq_treatments_all/Day14/GO.05brownModule.csv")
d14_GO.05magentaModule <- read.csv("Analysis/Output/GO/WGCNA_goseq/subseq_treatments_all/Day14/GO.05magentaModule.csv")
d14_GO.05pinkModule    <- read.csv("Analysis/Output/GO/WGCNA_goseq/subseq_treatments_all/Day14/GO.05pinkModule.csv")

d21_GO.05blackModule    <- read.csv("Analysis/Output/GO/WGCNA_goseq/subseq_treatments_all/Day21/GO.05blackModule.csv")
d21_GO.05blueModule     <- read.csv("Analysis/Output/GO/WGCNA_goseq/subseq_treatments_all/Day21/GO.05blueModule.csv")
d21_GO.05magentaModule  <- read.csv("Analysis/Output/GO/WGCNA_goseq/subseq_treatments_all/Day21/GO.05magentaModule.csv")
d21_GO.05pinkModule     <- read.csv("Analysis/Output/GO/WGCNA_goseq/subseq_treatments_all/Day21/GO.05pinkModule.csv")
d21_GO.05redModule      <- read.csv("Analysis/Output/GO/WGCNA_goseq/subseq_treatments_all/Day21/GO.05redModule.csv")
d21_GO.05yellowModule   <- read.csv("Analysis/Output/GO/WGCNA_goseq/subseq_treatments_all/Day21/GO.05yellowModule.csv")
d21_GO.05turquoiseModule <- read.csv("Analysis/Output/GO/WGCNA_goseq/subseq_treatments_all/Day21/GO.05turquoiseModule.csv")
# View(d21_GO.05yellowModule)
Geoduck_annotation <- read.delim2(file="C:/Users/samjg/Documents/My_Projects/Pgenerosa_TagSeq_Metabolomics/TagSeq/Seq_details/Panopea-generosa-genes-annotations.txt", header=F)
# build annotation file to merge with the mean LFC tables
Pgen_condenced <- Geoduck_annotation[,c(1,7)] # load just the PGEN ID and the putative gene terms
Pgen_condenced$Gene_term    <- sub(" \\(EC.*", "", Pgen_condenced$V7)  # str_extract(annot.condenced$V7, "[^(]+")
Pgen_reference <- na.omit(Pgen_condenced[c(1,3)])
names(Pgen_reference)  <- c('PgenIDs', 'Gene_terms')


# Master ALL WGCNA significant modules with treatment
Master_goseq_results      <- rbind(d0_GO.05midnightbluenModule,
                              d7_GO.05brownModule, d7_GO.05greenModule, d7_GO.05yellowModule,
                              d14_GO.05blackModule, d14_GO.05brownModule, d14_GO.05magentaModule, d14_GO.05pinkModule,
                              d21_GO.05blackModule, d21_GO.05blueModule, d21_GO.05magentaModule, d21_GO.05pinkModule, d21_GO.05redModule, d21_GO.05yellowModule, d21_GO.05turquoiseModule)
# View(Master_goseq_results)


# call all the module colors and days to loop GOslim analysis 
GOslimLoop_vars <- unique(Master_goseq_results[c(9,10)]) 



for (i in 1:nrow(GOslimLoop_vars)) {
  # call the target dataset
  goseq_res       <- Master_goseq_results %>%  dplyr::filter(Day %in% GOslimLoop_vars[i,1], moduleColor %in% GOslimLoop_vars[i,2])
  WGCNA_res       <- WGCNA_MasterModData  %>%  dplyr::filter(Day %in% GOslimLoop_vars[i,1], moduleColor %in% GOslimLoop_vars[i,2])
  gene_names      <- WGCNA_res$geneSymbol # all gene IDs in the particular WGCNA module 
  
  # Biological Function - run GOslim
  goseq_res_BP        <- goseq_res %>%  filter(ontology=="BP") # BP - all GO terms upregulated
  BP_GOcollection     <- GOCollection(goseq_res_BP$category)
  GOslims_BP          <- data.frame(goSlim(BP_GOcollection, slim, "BP")) #Find common parent terms to slim down our list
  GOslims_BP$category <- row.names(GOslims_BP) #save rownames as category
  
  # Molecular Function - run GOslim
  goseq_res_MF        <- goseq_res %>%  filter(ontology=="MF") # BP - all GO terms upregulated
  MF_GOcollection     <- GOCollection(goseq_res_MF$category)
  GOslims_MF          <- data.frame(goSlim(MF_GOcollection, slim, "MF")) #Find common parent terms to slim down our list
  GOslims_MF$category <- row.names(GOslims_MF) #save rownames as category
  
  # ====================================================================================
  # Get mapped terms - add to the GOslims datatable 
  # from Sam White's Biostars [post](https://support.bioconductor.org/p/128407/#128409).
  # ====================================================================================
  # Write function mappedIds to get the query terms that mapped to the slim categories 
  # ...in other words, add a column to your slim dataframe with all the GO terms from goseq
  mappedIds <-  function(df, collection, OFFSPRING) {  #the command to run requires a dataframe of slim terms, like slims_MF above, your list of query terms, and the offspring from the GOCollection by goSlim
    map <- as.list(OFFSPRING[rownames(df)]) # Subset GOcollection offspring by the rownames of your dataframe
    mapped <- lapply(map, intersect, ids(collection)) #Find the terms that intersect between the subset made above of your query terms and the GOids from the GO collection
    df[["go_terms"]] <- vapply(unname(mapped), paste, collapse = ";", character(1L)) #Add column "go_terms" with matching terms 
    df #show resulting dataframe
  }
   
  BPslim_Mapped <- mappedIds(GOslims_BP, BP_GOcollection, GOBPOFFSPRING)
  MFslim_Mapped <- mappedIds(GOslims_MF, MF_GOcollection, GOMFOFFSPRING)
  
  
  # BIOLOGICAL PROCESS
  BPslim             <- filter(BPslim_Mapped, Term!="biological_process") #filter out empty slims and term "biological process" and slims with < 2 GO terms (omitted the Count>=2)
  BPsplitted         <- strsplit(as.character(BPslim$go_terms), ";") #split into multiple GO ids
  BPslim$BPsplitted  <- BPsplitted
  for (n in 1:nrow(BPslim)) {
    table       <- data.frame(GOlist = unlist(BPslim[n,6])) # call the BPsplitted column of characters and create a small table to filter
    table       <- unique(table)
    
    Pgen_module  <- Pgen_GOterms2 %>% dplyr::filter(gene.ID %in% gene_names) # d7_Mod.Brown$X calls the gene names in the origin Module membership dataframe at the start of this script
    Pgen_loop    <- Pgen_module %>% dplyr::filter(Go.terms %in% table$GOlist) # filter Gene IDs with the GO term
    Pgen_geneIDs <- Pgen_loop[-2] # ommit the GO.terms to call unique gene calls 
    Pgen_geneIDs <- unique(Pgen_geneIDs) # call unique Gene calls (unique genes that had the GO term within each of the GOslim bins)
   
    BPslim$Gene.Count[n] <- nrow(Pgen_geneIDs)  # count of unique GeneIDs in each GOslim bin
    BPslim$Gene.IDs[[n]] <- vapply((Pgen_geneIDs$gene.ID), paste, collapse = ";", character(1L))
    } # end n in 1:nrow
  
  BPslim_A <- data.frame(Term = rep.int(BPslim$Term, sapply(BPsplitted, length)), go_term = unlist(BPsplitted)) #list all
  BPslim_B <- merge(BPslim_A, BPslim, by="Term") #Add back counts, term, and category info
  BPslim_C <- unique(setDT(BPslim_B)[order(go_term, -Gene.Count)], by = "category") #remove duplicate offspring terms, keeping only those in the larger umbrella term (Count number)
  BPslim_final <- BPslim_C[,c(1,5,3,2,6,8:9)]
  colnames(BPslim_final) <- c("slim_term", "slim_cat", "GO_count", "GO_terms", "GO_list", "Gene_count", "Gene_IDs")
  BPslim_final[["Gene_IDs"]] <- vapply(unname(BPslim_final$Gene_IDs), paste, collapse = ";", character(1L)) # convert from a list to simply a charaacter string with ; delimiter
  #BPslim_final <- data.frame(slim_term=BPslim_C$Term, slim_cat=BPslim_C$category, category=BPslim_C$go_term, Gene.Count=BPslim_C$Gene.Count, GO.Count=BPslim_C$Count, Gene.IDs=BPslim_C$Gene.IDs) #rename columns) #rename columns
  BPslim_final$module_day <- paste(GOslimLoop_vars[i,2], GOslimLoop_vars[i,1], sep = '_')
  
            if (nrow(goseq_res_BP) >0) {
            BPslim_GOterm_summary       <- BPslim_final[,c(1,2,5,7)]
            s <- strsplit(BPslim_GOterm_summary$GO_list, split = ";")
            BPslim_GOterm_summary_2     <- data.frame(slim_term = rep(BPslim_GOterm_summary$slim_term, sapply(s, length)),
                                                      Gene_IDs  = rep(BPslim_GOterm_summary$Gene_IDs, sapply(s, length)),
                                                       GO_terms = unlist(s))
            colnames(goseq_res_BP)[2]   <- 'GO_terms'
            BPslim_GOterm_summary_final <- merge(goseq_res_BP, BPslim_GOterm_summary_2, by = 'GO_terms')
            
            
            s_2 <- strsplit(BPslim_GOterm_summary_final$Gene_IDs, split = ";")
            BPslim_GOterm_gene_annotation <- data.frame(slim_term      = rep(BPslim_GOterm_summary_final$slim_term, sapply(s_2, length)),
                                                        over_represented_pvalue = rep(BPslim_GOterm_summary_final$over_represented_pvalue, sapply(s_2, length)),
                                                        GO_terms       = rep(BPslim_GOterm_summary_final$GO_terms, sapply(s_2, length)),
                                                        term           = rep(BPslim_GOterm_summary_final$term, sapply(s_2, length)),
                                                        ontology       = rep(BPslim_GOterm_summary_final$ontology, sapply(s_2, length)),
                                                        moduleColor    = rep(BPslim_GOterm_summary_final$moduleColor, sapply(s_2, length)),
                                                        PgenIDs        = unlist(s_2))
            BP_master_gene_reference <- merge(BPslim_GOterm_gene_annotation, Pgen_reference, by = "PgenIDs")
            
            } else (c(BPslim_GOterm_summary_final = NULL, BP_master_gene_reference = NULL))
  
  # MOLECULAR FUNCTION
  MFslim             <- filter(MFslim_Mapped, Term!="molecular_function") #filter out empty slims and term "biological process" (omitted the Count>=2)
  MFsplitted         <- strsplit(as.character(MFslim$go_terms), ";") #split into multiple GO ids
  MFslim$MFsplitted  <- MFsplitted
  for (m in 1:nrow(MFslim)) {
    table       <- data.frame(GOlist = unlist(MFslim[m,6])) # call the MFsplitted column of characters and create a small table to filter
    table       <- unique(table)
    
    Pgen_module  <- Pgen_GOterms2 %>% dplyr::filter(gene.ID %in% gene_names) # d7_Mod.Brown$X calls the gene names in the origin Module membership dataframe at the start of this script
    Pgen_loop    <- Pgen_module %>% dplyr::filter(Go.terms %in% table$GOlist) # filter Gene IDs with the GO term
    Pgen_geneIDs <- Pgen_loop[-2] # ommit the GO.terms to call unique gene calls 
    Pgen_geneIDs <- unique(Pgen_geneIDs) # call unique Gene calls (unique genes that had the GO term within each of the GOslim bins)
    
    MFslim$Gene.Count[m] <- nrow(Pgen_geneIDs)  # count of unique GeneIDs in each GOslim bin
    MFslim$Gene.IDs[[m]] <- vapply((Pgen_geneIDs$gene.ID), paste, collapse = ";", character(1L))
    } # end m in 1:nrow
  
    MFslim_A <- data.frame(Term = rep.int(MFslim$Term, sapply(MFsplitted, length)), go_term = unlist(MFsplitted)) #list all
    MFslim_B <- merge(MFslim_A, MFslim, by="Term") #Add back counts, term, and category info
    MFslim_C <- unique(setDT(MFslim_B)[order(go_term, -Gene.Count)], by = "category") #remove duplicate offspring terms, keeping only those in the larger umbrella term (Count number)
    MFslim_final <- MFslim_C[,c(1,5,3,2,6,8:9)]
    colnames(MFslim_final) <- c("slim_term", "slim_cat", "GO_count", "GO_terms", "GO_list", "Gene_count", "Gene_IDs")
    MFslim_final[["Gene_IDs"]] <- vapply(unname(MFslim_final$Gene_IDs), paste, collapse = ";", character(1L)) # convert from a list to simply a charaacter string with ; delimiter
    # MFslim_final <- data.frame(slim_term=MFslim_C$Term, slim_cat=MFslim_C$category, category=MFslim_C$go_term, Gene.Count=MFslim_C$Gene.Count, GO.Count=MFslim_C$Count) #rename columns) #rename columns
    MFslim_final$module_day <- paste(GOslimLoop_vars[i,2], GOslimLoop_vars[i,1], sep = '_')
  
            if (nrow(goseq_res_MF) >0) {
            MFslim_GOterm_summary       <- MFslim_final[,c(1,2,5,7)]
            s <- strsplit(MFslim_GOterm_summary$GO_list, split = ";")
            MFslim_GOterm_summary_2     <- data.frame(slim_term = rep(MFslim_GOterm_summary$slim_term, sapply(s, length)),
                                                      Gene_IDs  = rep(MFslim_GOterm_summary$Gene_IDs, sapply(s, length)),
                                                      GO_terms = unlist(s))
            colnames(goseq_res_MF)[2]   <- 'GO_terms'
            MFslim_GOterm_summary_final <- merge(goseq_res_MF, MFslim_GOterm_summary_2, by = 'GO_terms')
            
            s_2 <- strsplit(MFslim_GOterm_summary_final$Gene_IDs, split = ";")
            MFslim_GOterm_gene_annotation <- data.frame(slim_term      = rep(MFslim_GOterm_summary_final$slim_term, sapply(s_2, length)),
                                                        over_represented_pvalue = rep(MFslim_GOterm_summary_final$over_represented_pvalue, sapply(s_2, length)),
                                                        GO_terms       = rep(MFslim_GOterm_summary_final$GO_terms, sapply(s_2, length)),
                                                        term           = rep(MFslim_GOterm_summary_final$term, sapply(s_2, length)),
                                                        ontology       = rep(MFslim_GOterm_summary_final$ontology, sapply(s_2, length)),
                                                        moduleColor    = rep(MFslim_GOterm_summary_final$moduleColor, sapply(s_2, length)),
                                                        PgenIDs        = unlist(s_2))
            MF_master_gene_reference <- merge(MFslim_GOterm_gene_annotation, Pgen_reference, by = "PgenIDs")
    
    } else (c(MFslim_GOterm_summary_final = NULL, MF_master_gene_reference = NULL))
    
    # save GOslim final datasets for BP and MF of each  module in their respective folder(s) by Day
    write.csv(MFslim_final, file = paste("Analysis/Output/GO/WGCNA_goseq/subseq_treatments_all/", GOslimLoop_vars[i,1], "/GOslim_MolFunction_",GOslimLoop_vars[i,2], "Module.csv", sep ='')) # save csv file              
    write.csv(BPslim_final, file = paste("Analysis/Output/GO/WGCNA_goseq/subseq_treatments_all/", GOslimLoop_vars[i,1], "/GOslim_BiolProc_",GOslimLoop_vars[i,2], "Module.csv", sep ='')) # save csv file       

    write.csv(MFslim_GOterm_summary_final, file = paste("Analysis/Output/GO/WGCNA_goseq/subseq_treatments_all/", GOslimLoop_vars[i,1], "/GOterms_and_GOslim_MolFunction_",GOslimLoop_vars[i,2], "Module.csv", sep ='')) # save csv file              
    write.csv(BPslim_GOterm_summary_final, file = paste("Analysis/Output/GO/WGCNA_goseq/subseq_treatments_all/", GOslimLoop_vars[i,1], "/GOterms_and_GOslim_BiolProc_",GOslimLoop_vars[i,2], "Module.csv", sep ='')) # save csv file       
    
    write.csv(MF_master_gene_reference, file = paste("Analysis/Output/GO/WGCNA_goseq/subseq_treatments_all/", GOslimLoop_vars[i,1], "/GOterms_and_GOslim_MolFunction_",GOslimLoop_vars[i,2], "Module_GENE_REFERENCE.csv", sep ='')) # save csv file              
    write.csv(BP_master_gene_reference, file = paste("Analysis/Output/GO/WGCNA_goseq/subseq_treatments_all/", GOslimLoop_vars[i,1], "/GOterms_and_GOslim_BiolProc_",GOslimLoop_vars[i,2], "Module_GENE_REFERENCE.csv", sep ='')) # save csv file       
  
    
    print(paste(GOslimLoop_vars[i,2], GOslimLoop_vars[i,1], "done", sep = ' '))
} # end of GOslim for loop! (i in 1:nrow)

#===================================================================================================
#
#
#  Tally up persistant and transient gene groups 
# 
#  Why? This will allow interpretation of the resutls downstream
#
#===================================================================================================
# load the output from the previous for loop
# BIOLOGICAL PROCESS  GO  SLIMS
# primary effect modules; ambient > moderate
d7_slimBP_brownModule   <- read.csv("Analysis/Output/GO/WGCNA_goseq/subseq_treatments_all/Day7/GOterms_and_GOslim_BiolProc_brownModule_GENE_REFERENCE.csv")
d7_slimBP_brownModule$Day <- "Day7"

d14_slimBP_brownModule  <- read.csv("Analysis/Output/GO/WGCNA_goseq/subseq_treatments_all/Day14/GOterms_and_GOslim_BiolProc_brownModule_GENE_REFERENCE.csv")
d14_slimBP_brownModule$Day <- "Day14"

d21_slimBP_blueModule   <- read.csv("Analysis/Output/GO/WGCNA_goseq/subseq_treatments_all/Day21/GOterms_and_GOslim_BiolProc_blueModule_GENE_REFERENCE.csv")
d21_slimBP_blueModule$Day <- "Day21"
# d21_slimBP_magentaModule     <- read.csv("Analysis/Output/GO/WGCNA_goseq/subseq_treatments_all/Day21/GOterms_and_GOslim_BiolProc_magentaModule.csv") # no BP terms!

# primary effect modules; moderate > ambient 
d7_slimBP_yellowModule  <- read.csv("Analysis/Output/GO/WGCNA_goseq/subseq_treatments_all/Day7/GOterms_and_GOslim_BiolProc_yellowModule_GENE_REFERENCE.csv")
d7_slimBP_yellowModule$Day <- "Day7"

d14_slimBP_blackModule  <- read.csv("Analysis/Output/GO/WGCNA_goseq/subseq_treatments_all/Day14/GOterms_and_GOslim_BiolProc_blackModule_GENE_REFERENCE.csv")
d14_slimBP_blackModule$Day <- "Day14"

d21_slimBP_yellowModule <- read.csv("Analysis/Output/GO/WGCNA_goseq/subseq_treatments_all/Day21/GOterms_and_GOslim_BiolProc_yellowModule_GENE_REFERENCE.csv")
d21_slimBP_yellowModule$Day <- "Day21"


# MOLECULAR FUNCTION GO  SLIMS
# primary effect modules; ambient > moderate
d7_slimMF_brownModule   <- read.csv("Analysis/Output/GO/WGCNA_goseq/subseq_treatments_all/Day7/GOterms_and_GOslim_MolFunction_brownModule_GENE_REFERENCE.csv")
d7_slimMF_brownModule$Day <- "Day7"

d14_slimMF_brownModule  <- read.csv("Analysis/Output/GO/WGCNA_goseq/subseq_treatments_all/Day14/GOterms_and_GOslim_MolFunction_brownModule_GENE_REFERENCE.csv")
d14_slimMF_brownModule$Day <- "Day14"

d21_slimMF_blueModule   <- read.csv("Analysis/Output/GO/WGCNA_goseq/subseq_treatments_all/Day21/GOterms_and_GOslim_MolFunction_blueModule_GENE_REFERENCE.csv")
d21_slimMF_blueModule$Day <- "Day21"

d21_slimMF_magentaModule<- read.csv("Analysis/Output/GO/WGCNA_goseq/subseq_treatments_all/Day21/GOterms_and_GOslim_MolFunction_magentaModule_GENE_REFERENCE.csv")
d21_slimMF_magentaModule$Day <- "Day21"

# primary effect modules; moderate > ambient 
d7_slimMF_yellowModule  <- read.csv("Analysis/Output/GO/WGCNA_goseq/subseq_treatments_all/Day7/GOterms_and_GOslim_MolFunction_yellowModule_GENE_REFERENCE.csv")
d7_slimMF_yellowModule$Day <- "Day7"

d14_slimMF_blackModule  <- read.csv("Analysis/Output/GO/WGCNA_goseq/subseq_treatments_all/Day14/GOterms_and_GOslim_MolFunction_blackModule_GENE_REFERENCE.csv")
d14_slimMF_blackModule$Day <- "Day14"

d21_slimMF_yellowModule <- read.csv("Analysis/Output/GO/WGCNA_goseq/subseq_treatments_all/Day21/GOterms_and_GOslim_MolFunction_yellowModule_GENE_REFERENCE.csv")
d21_slimMF_yellowModule$Day <- "Day21"

# BIND DATA TO TALLY BY GROUP USING 'dplyr'  ---------------------------------------------------------- #
BP_Amb <- rbind(d7_slimBP_brownModule, d14_slimBP_brownModule,  d21_slimBP_blueModule) # merge the data by binding rows 
MF_Amb <- rbind(d7_slimMF_brownModule, d14_slimMF_brownModule, d21_slimMF_magentaModule, d21_slimMF_blueModule) # merge the data by binding rows 

# (2.B) Priming Effect MODERATE > AMBIENT 
BP_Mod <- rbind(d7_slimBP_yellowModule, d14_slimBP_blackModule, d21_slimBP_yellowModule) # merge the data by binding rows 
MF_Mod <- rbind(d7_slimMF_yellowModule, d14_slimMF_blackModule, d21_slimMF_yellowModule) # merge the data by binding rows 

# Create datasets for interpreting patterns in response to priamry treatment modules in WGCNA
library(stringr)




# Biological Process - Ambient > Moderate 'BP_Amb_AGGREG_2' == end product for output 
BP_Amb_2 <- BP_Amb %>%  dplyr::select(c('slim_term', 'Day', 'Gene_terms')) # selct columns of interest
BP_Amb_uniq <- unique(BP_Amb_2) # call all unique terms ()
BP_Amb_uniq$Gene_terms    <- sub(" \\(.*", "", BP_Amb_uniq$Gene_terms) 

BP_Amb_TALLY <- BP_Amb_uniq %>% dplyr::group_by(slim_term,  Gene_terms) %>% dplyr::mutate(count = n(), Day = Day) # tally the number of occurnace (count = repeated Days with the gene/term occured in the module)
BP_Amb_AGGREG <- aggregate(Day ~ Gene_terms*slim_term*count, BP_Amb_TALLY,function(x) toString(unique(x))) # integrate the DAYS corresonding with this count (i.e. count = 2; Days 7 and day 14)

BP_Amb_TALLY_2 <- BP_Amb_AGGREG %>% dplyr::group_by(slim_term, Day) %>% dplyr::mutate(count = n()) # now tally the number of unique gene terms corresponding with the unique or persistant occurances by day (count now == the number of unique genes per go slim term)
BP_Amb_AGGREG_2 <- aggregate(Gene_terms ~ Day*slim_term*count, BP_Amb_TALLY_2,function(x) toString(unique(x))) # integrate the GENE TERMS corresonding with this count

BP_Amb_AGGREG_2$Pattern <- ifelse(grepl("\\<Day7, Day21\\>", BP_Amb_AGGREG_2$Day, ignore.case = TRUE), "Transient genes: consistantly induced by stress", 
                                  ifelse(grepl("\\<Day7, Day14, Day21\\>", BP_Amb_AGGREG_2$Day, ignore.case = TRUE), "Persistant genes",  
                                  ifelse(grepl("\\<Day14, Day21\\>", BP_Amb_AGGREG_2$Day, ignore.case = TRUE), "Transient genes: ambient recovery + stress_3 = preparatory regulation?", 
                                  ifelse(grepl("\\<Day7, Day14\\>", BP_Amb_AGGREG_2$Day, ignore.case = TRUE), "Transient genes: stress_2 + ambient recovery - frontload?",
                                  ifelse(grepl("\\<Day14\\>", BP_Amb_AGGREG_2$Day, ignore.case = TRUE), "Transient genes: ambient recovery = preparatory regulation?", 
                                  ifelse(grepl("\\<Day7\\>|\\<Day21\\>", BP_Amb_AGGREG_2$Day, ignore.case = TRUE), "Transient genes: unique stress response","Other"))))))





# Biological Process - Moderate > Ambient 'BP_Mod_AGGREG_2' == end product for output 
BP_Mod_2 <- BP_Mod %>%  dplyr::select(c('slim_term', 'Day', 'Gene_terms')) # selct columns of interest
BP_Mod_uniq <- unique(BP_Mod_2) # call all unique terms ()
BP_Mod_uniq$Gene_terms    <- sub(" \\(.*", "", BP_Mod_uniq$Gene_terms) 

BP_Mod_TALLY <- BP_Mod_uniq %>% dplyr::group_by(slim_term,  Gene_terms) %>% dplyr::mutate(count = n(), Day = Day) # tally the number of occurnace (count = repeated Days with the gene/term occured in the module)
BP_Mod_AGGREG <- aggregate(Day ~ Gene_terms*slim_term*count, BP_Mod_TALLY,function(x) toString(unique(x))) # integrate the DAYS corresonding with this count (i.e. count = 2; Days 7 and day 14)

BP_Mod_TALLY_2 <- BP_Mod_AGGREG %>% dplyr::group_by(slim_term, Day) %>% dplyr::mutate(count = n()) # now tally the number of unique gene terms corresponding with the unique or persistant occurances by day (count now == the number of unique genes per go slim term)
BP_Mod_AGGREG_2 <- aggregate(Gene_terms ~ Day*slim_term*count, BP_Mod_TALLY_2,function(x) toString(unique(x))) # integrate the GENE TERMS corresonding with this count

BP_Mod_AGGREG_2$Pattern <- ifelse(grepl("\\<Day7, Day21\\>", BP_Mod_AGGREG_2$Day, ignore.case = TRUE), "Transient genes: consistantly induced by stress", 
                                  ifelse(grepl("\\<Day7, Day14, Day21\\>", BP_Mod_AGGREG_2$Day, ignore.case = TRUE), "Persistant genes",  
                                         ifelse(grepl("\\<Day14, Day21\\>", BP_Mod_AGGREG_2$Day, ignore.case = TRUE), "Transient genes: Modient recovery + stress_3 = preparatory regulation?", 
                                                ifelse(grepl("\\<Day7, Day14\\>", BP_Mod_AGGREG_2$Day, ignore.case = TRUE), "Transient genes: stress_2 + Modient recovery - frontload?",
                                                       ifelse(grepl("\\<Day14\\>", BP_Mod_AGGREG_2$Day, ignore.case = TRUE), "Transient genes: Ambient recovery = preparatory regulation?", 
                                                              ifelse(grepl("\\<Day7\\>|\\<Day21\\>", BP_Mod_AGGREG_2$Day, ignore.case = TRUE), "Transient genes: unique stress response","Other"))))))







# Biological Process - Ambient > Moderate 'MF_Amb_AGGREG_2' == end product for output 
MF_Amb_2 <- MF_Amb %>%  dplyr::select(c('slim_term', 'Day', 'Gene_terms')) # selct columns of interest
MF_Amb_uniq <- unique(MF_Amb_2) # call all unique terms ()
MF_Amb_uniq$Gene_terms    <- sub(" \\(.*", "", MF_Amb_uniq$Gene_terms) 

MF_Amb_TALLY <- MF_Amb_uniq %>% dplyr::group_by(slim_term,  Gene_terms) %>% dplyr::mutate(count = n(), Day = Day) # tally the number of occurnace (count = repeated Days with the gene/term occured in the Ambule)
MF_Amb_AGGREG <- aggregate(Day ~ Gene_terms*slim_term*count, MF_Amb_TALLY,function(x) toString(unique(x))) # integrate the DAYS corresonding with this count (i.e. count = 2; Days 7 and day 14)

MF_Amb_TALLY_2 <- MF_Amb_AGGREG %>% dplyr::group_by(slim_term, Day) %>% dplyr::mutate(count = n()) # now tally the number of unique gene terms corresponding with the unique or persistant occurances by day (count now == the number of unique genes per go slim term)
MF_Amb_AGGREG_2 <- aggregate(Gene_terms ~ Day*slim_term*count, MF_Amb_TALLY_2,function(x) toString(unique(x))) # integrate the GENE TERMS corresonding with this count

MF_Amb_AGGREG_2$Pattern <- ifelse(grepl("\\<Day7, Day21\\>", MF_Amb_AGGREG_2$Day, ignore.case = TRUE), "Transient genes: consistantly induced by stress", 
                                  ifelse(grepl("\\<Day7, Day14, Day21\\>", MF_Amb_AGGREG_2$Day, ignore.case = TRUE), "Persistant genes",  
                                         ifelse(grepl("\\<Day14, Day21\\>", MF_Amb_AGGREG_2$Day, ignore.case = TRUE), "Transient genes: Ambient recovery + stress_3 = preparatory regulation?", 
                                                ifelse(grepl("\\<Day7, Day14\\>", MF_Amb_AGGREG_2$Day, ignore.case = TRUE), "Transient genes: stress_2 + Ambient recovery - frontload?",
                                                       ifelse(grepl("\\<Day14\\>", MF_Amb_AGGREG_2$Day, ignore.case = TRUE), "Transient genes: Ambient recovery = preparatory regulation?", 
                                                              ifelse(grepl("\\<Day7\\>|\\<Day21\\>", MF_Amb_AGGREG_2$Day, ignore.case = TRUE), "Transient genes: unique stress response","Other"))))))






# Biological Process - Moderate > Ambient 'MF_Mod_AGGREG_2' == end product for output 
MF_Mod_2 <- MF_Mod %>%  dplyr::select(c('slim_term', 'Day', 'Gene_terms')) # selct columns of interest
MF_Mod_uniq <- unique(MF_Mod_2) # call all unique terms ()
MF_Mod_uniq$Gene_terms    <- sub(" \\(.*", "", MF_Mod_uniq$Gene_terms) 

MF_Mod_TALLY <- MF_Mod_uniq %>% dplyr::group_by(slim_term,  Gene_terms) %>% dplyr::mutate(count = n(), Day = Day) # tally the number of occurnace (count = repeated Days with the gene/term occured in the module)
MF_Mod_AGGREG <- aggregate(Day ~ Gene_terms*slim_term*count, MF_Mod_TALLY,function(x) toString(unique(x))) # integrate the DAYS corresonding with this count (i.e. count = 2; Days 7 and day 14)

MF_Mod_TALLY_2 <- MF_Mod_AGGREG %>% dplyr::group_by(slim_term, Day) %>% dplyr::mutate(count = n()) # now tally the number of unique gene terms corresponding with the unique or persistant occurances by day (count now == the number of unique genes per go slim term)
MF_Mod_AGGREG_2 <- aggregate(Gene_terms ~ Day*slim_term*count, MF_Mod_TALLY_2,function(x) toString(unique(x))) # integrate the GENE TERMS corresonding with this count

MF_Mod_AGGREG_2$Pattern <- ifelse(grepl("\\<Day7, Day21\\>", MF_Mod_AGGREG_2$Day, ignore.case = TRUE), "Transient genes: consistantly induced by stress", 
                                  ifelse(grepl("\\<Day7, Day14, Day21\\>", MF_Mod_AGGREG_2$Day, ignore.case = TRUE), "Persistant genes",  
                                         ifelse(grepl("\\<Day14, Day21\\>", MF_Mod_AGGREG_2$Day, ignore.case = TRUE), "Transient genes: Modient recovery + stress_3 = preparatory regulation?", 
                                                ifelse(grepl("\\<Day7, Day14\\>", MF_Mod_AGGREG_2$Day, ignore.case = TRUE), "Transient genes: stress_2 + Ambient recovery - frontload?",
                                                       ifelse(grepl("\\<Day14\\>", MF_Mod_AGGREG_2$Day, ignore.case = TRUE), "Transient genes: Modient recovery = preparatory regulation?", 
                                                              ifelse(grepl("\\<Day7\\>|\\<Day21\\>", MF_Mod_AGGREG_2$Day, ignore.case = TRUE), "Transient genes: unique stress response","Other"))))))





# write csv files 
write.csv(BP_Amb_AGGREG_2, file = paste("Analysis/Output/GO/WGCNA_goseq/BiologicalProcess_PrimEff_Ambient_Master.csv", sep ='')) # save csv file       
write.csv(MF_Amb_AGGREG_2, file = paste("Analysis/Output/GO/WGCNA_goseq/MolecularFunction_PrimEff_Ambient_Master.csv", sep ='')) # save csv file   

write.csv(BP_Mod_AGGREG_2, file = paste("Analysis/Output/GO/WGCNA_goseq/BiologicalProcess_PrimEff_Moderate_Master.csv", sep ='')) # save csv file       
write.csv(MF_Mod_AGGREG_2, file = paste("Analysis/Output/GO/WGCNA_goseq/MolecularFunction_PrimEff_Moderate_Master.csv", sep ='')) # save csv file   





#===================================================================================================
#
#
#  Youre done! Now time to load in those final GOslim csvs and plot!!!!
# 
#   PLOTTING
#
#===================================================================================================
# load the output from the previous for loop
# BIOLOGICAL PROCESS  GO  SLIMS
# View(d7_slimBP_brownModule)
d7_slimBP_brownModule   <- read.csv("Analysis/Output/GO/WGCNA_goseq/subseq_treatments_all/Day7/GOterms_and_GOslim_BiolProc_brownModule.csv")
d7_slimBP_greenModule   <- read.csv("Analysis/Output/GO/WGCNA_goseq/subseq_treatments_all/Day7/GOterms_and_GOslim_BiolProc_greenModule.csv")
d7_slimBP_yellowModule  <- read.csv("Analysis/Output/GO/WGCNA_goseq/subseq_treatments_all/Day7/GOterms_and_GOslim_BiolProc_yellowModule.csv")

d14_slimBP_blackModule   <- read.csv("Analysis/Output/GO/WGCNA_goseq/subseq_treatments_all/Day14/GOterms_and_GOslim_BiolProc_blackModule.csv")
d14_slimBP_brownModule   <- read.csv("Analysis/Output/GO/WGCNA_goseq/subseq_treatments_all/Day14/GOterms_and_GOslim_BiolProc_brownModule.csv")
d14_slimBP_magentaModule <- read.csv("Analysis/Output/GO/WGCNA_goseq/subseq_treatments_all/Day14/GOterms_and_GOslim_BiolProc_magentaModule.csv")
d14_slimBP_pinkModule    <- read.csv("Analysis/Output/GO/WGCNA_goseq/subseq_treatments_all/Day14/GOterms_and_GOslim_BiolProc_pinkModule.csv")

d21_slimBP_blackModule       <- read.csv("Analysis/Output/GO/WGCNA_goseq/subseq_treatments_all/Day21/GOterms_and_GOslim_BiolProc_blackModule.csv")
d21_slimBP_blueModule        <- read.csv("Analysis/Output/GO/WGCNA_goseq/subseq_treatments_all/Day21/GOterms_and_GOslim_BiolProc_blueModule.csv")
# d21_slimBP_magentaModule     <- read.csv("Analysis/Output/GO/WGCNA_goseq/subseq_treatments_all/Day21/GOterms_and_GOslim_BiolProc_magentaModule.csv") # no BP terms!
d21_slimBP_pinkModule        <- read.csv("Analysis/Output/GO/WGCNA_goseq/subseq_treatments_all/Day21/GOterms_and_GOslim_BiolProc_pinkModule.csv")
d21_slimBP_redModule         <- read.csv("Analysis/Output/GO/WGCNA_goseq/subseq_treatments_all/Day21/GOterms_and_GOslim_BiolProc_redModule.csv")
d21_slimBP_yellowModule      <- read.csv("Analysis/Output/GO/WGCNA_goseq/subseq_treatments_all/Day21/GOterms_and_GOslim_BiolProc_yellowModule.csv")
d21_slimBP_turquoiseModule   <- read.csv("Analysis/Output/GO/WGCNA_goseq/subseq_treatments_all/Day21/GOterms_and_GOslim_BiolProc_turquoiseModule.csv")

# MOLECULAR FUNCTION GO  SLIMS
d7_slimMF_brownModule   <- read.csv("Analysis/Output/GO/WGCNA_goseq/subseq_treatments_all/Day7/GOterms_and_GOslim_MolFunction_brownModule.csv")
d7_slimMF_greenModule   <- read.csv("Analysis/Output/GO/WGCNA_goseq/subseq_treatments_all/Day7/GOterms_and_GOslim_MolFunction_greenModule.csv")
d7_slimMF_yellowModule  <- read.csv("Analysis/Output/GO/WGCNA_goseq/subseq_treatments_all/Day7/GOterms_and_GOslim_MolFunction_yellowModule.csv")

d14_slimMF_blackModule   <- read.csv("Analysis/Output/GO/WGCNA_goseq/subseq_treatments_all/Day14/GOterms_and_GOslim_MolFunction_blackModule.csv")
d14_slimMF_brownModule   <- read.csv("Analysis/Output/GO/WGCNA_goseq/subseq_treatments_all/Day14/GOterms_and_GOslim_MolFunction_brownModule.csv")
d14_slimMF_magentaModule <- read.csv("Analysis/Output/GO/WGCNA_goseq/subseq_treatments_all/Day14/GOterms_and_GOslim_MolFunction_magentaModule.csv")
d14_slimMF_pinkModule    <- read.csv("Analysis/Output/GO/WGCNA_goseq/subseq_treatments_all/Day14/GOterms_and_GOslim_MolFunction_pinkModule.csv")

d21_slimMF_blackModule       <- read.csv("Analysis/Output/GO/WGCNA_goseq/subseq_treatments_all/Day21/GOterms_and_GOslim_MolFunction_blackModule.csv")
d21_slimMF_blueModule        <- read.csv("Analysis/Output/GO/WGCNA_goseq/subseq_treatments_all/Day21/GOterms_and_GOslim_MolFunction_blueModule.csv")
d21_slimMF_magentaModule     <- read.csv("Analysis/Output/GO/WGCNA_goseq/subseq_treatments_all/Day21/GOterms_and_GOslim_MolFunction_magentaModule.csv")
d21_slimMF_pinkModule        <- read.csv("Analysis/Output/GO/WGCNA_goseq/subseq_treatments_all/Day21/GOterms_and_GOslim_MolFunction_pinkModule.csv")
d21_slimMF_redModule         <- read.csv("Analysis/Output/GO/WGCNA_goseq/subseq_treatments_all/Day21/GOterms_and_GOslim_MolFunction_redModule.csv")
d21_slimMF_yellowModule      <- read.csv("Analysis/Output/GO/WGCNA_goseq/subseq_treatments_all/Day21/GOterms_and_GOslim_MolFunction_yellowModule.csv")
d21_slimMF_turquoiseModule   <- read.csv("Analysis/Output/GO/WGCNA_goseq/subseq_treatments_all/Day21/GOterms_and_GOslim_MolFunction_turquoiseModule.csv")

# BIND DATA TO FACET THE PLOTS    ---------------------------------------------------------- #
# (1) By sampling date 0 Day 7 14 and 21 separately
BP_D7  <- rbind(d7_slimBP_brownModule, d7_slimBP_greenModule, d7_slimBP_yellowModule)
BP_D14 <- rbind(d14_slimBP_blackModule, d14_slimBP_brownModule, d14_slimBP_magentaModule, d14_slimBP_pinkModule)
BP_D21 <- rbind(d21_slimBP_blackModule, d21_slimBP_blueModule, d21_slimBP_pinkModule, d21_slimBP_redModule, d21_slimBP_yellowModule,d21_slimBP_turquoiseModule)

MF_D7  <- rbind(d7_slimMF_brownModule, d7_slimMF_greenModule, d7_slimMF_yellowModule)
MF_D14 <- rbind(d14_slimMF_blackModule, d14_slimMF_brownModule, d14_slimMF_magentaModule, d14_slimMF_pinkModule)
MF_D21 <- rbind(d21_slimMF_blackModule, d21_slimMF_blueModule, d21_slimMF_magentaModule, d21_slimMF_pinkModule, d21_slimMF_redModule, d21_slimMF_yellowModule, d21_slimMF_turquoiseModule)

# (2) By trend in vst expression pattern ACROSS all days 
# (2.A) Priming Effect AMBIENT > MODERATE 
BP_Amb <- rbind(d7_slimBP_brownModule, d14_slimBP_brownModule,  d21_slimBP_blueModule) # merge the data by binding rows 
MF_Amb <- rbind(d7_slimMF_brownModule, d14_slimMF_brownModule, d21_slimMF_magentaModule, d21_slimMF_blueModule) # merge the data by binding rows 

# (2.B) Priming Effect MODERATE > AMBIENT 
BP_Mod <- rbind(d7_slimBP_yellowModule, d14_slimBP_blackModule, d21_slimBP_yellowModule) # merge the data by binding rows 
MF_Mod <- rbind(d7_slimMF_yellowModule, d14_slimMF_blackModule, d21_slimMF_yellowModule) # merge the data by binding rows 

BP_D7$Ont        <- "BP" # create a new common clumn to call in the plot 
BP_D7$module_day <- factor(BP_D7$module_day, levels = c("Day7_brown", "Day7_yellow", "Day7_green"))# reorder the facotr level for the facet wrap plot 
# BP_D7_filtered <- BP_D7 %>%  dplyr::filter(Gene_count >= 10) # ommit all terms gene counts 51
# BP_D7_filtered$slim_term <- factor(BP_D7_filtered$slim_term) # make slim term alphabetical for plotting
BP_D7$slim_term  <- factor(BP_D7$slim_term) # make slim term alphabetical for plotting





############################################################################################################## 3
#  WCNA Day 7
BP_D7$Day_moduleColor <- paste(BP_D7$Day, BP_D7$moduleColor, sep = '_')
BP_D7$Day_moduleColor <- factor(BP_D7$Day_moduleColor, levels = c("Day7_brown", "Day7_yellow", "Day7_green"))
BP_D7_Plot <-ggplot(data = BP_D7, aes(x = Day_moduleColor, y = forcats::fct_rev(term))) + 
  geom_tile(aes(fill= -log10(over_represented_pvalue), width = 1)) + 
  scale_fill_gradient(low = "skyblue", high = "navyblue",limits=c(0, (ceiling(max(-log10(BP_D7$over_represented_pvalue))))), breaks=seq(0, (ceiling(max(-log10(BP_D7$over_represented_pvalue)))),by=2))  +
  facet_grid(~module_day, labeller = label_wrap_gen(width = 10, multi_line = TRUE))+
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                     strip.text.y = element_text(angle=0, size = 8, face = "bold"),
                     strip.text.x = element_text(size = 8, face = "bold"),
                     axis.title.x = element_blank(),
                     axis.title.y = element_text(size=8),
                     axis.text = element_text(size = 8), legend.position = "right",
                     plot.margin = unit(c(0,1,0,0.25), "cm"))+
  ggtitle('Biological Process: Day 7 WGCNA') +
  # facet_wrap(~ (substr(slim_term, 1,30)))
  facet_wrap((substr(slim_term, 1,30)) ~., scales="free_y", ncol= 1, strip.position="right")

pdf(paste("Analysis/Output/GO/WGCNA_goseq/GOterm_enrichment_plots/Biol_Process_Day7_heatmap.pdf", sep =''), width=15, height=25)
print(BP_D7_Plot)
dev.off()


MF_D7$Day_moduleColor <- paste(MF_D7$Day, MF_D7$moduleColor, sep = '_')
MF_D7$Day_moduleColor <- factor(MF_D7$Day_moduleColor, levels = c("Day7_brown", "Day7_yellow", "Day7_green"))
MF_D7_Plot <-ggplot(data = MF_D7, aes(x = Day_moduleColor, y = forcats::fct_rev(term))) + 
  geom_tile(aes(fill= -log10(over_represented_pvalue), width = 1)) + 
  scale_fill_gradient(low = "skyblue", high = "navyblue",limits=c(0, (ceiling(max(-log10(MF_D7$over_represented_pvalue))))), breaks=seq(0, (ceiling(max(-log10(MF_D7$over_represented_pvalue)))),by=2))  +
  facet_grid(~module_day, labeller = label_wrap_gen(width = 10, multi_line = TRUE))+
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                     strip.text.y = element_text(angle=0, size = 8, face = "bold"),
                     strip.text.x = element_text(size = 8, face = "bold"),
                     axis.title.x = element_blank(),
                     axis.title.y = element_text(size=8),
                     axis.text = element_text(size = 8), legend.position = "right",
                     plot.margin = unit(c(0,1,0,0.25), "cm"))+
  ggtitle('Molecular Function: Primary Moderate Effect Day 7') +
  # facet_wrap(~ (substr(slim_term, 1,30)))
  facet_wrap((substr(slim_term, 1,30)) ~., scales="free_y", ncol= 1, strip.position="right")


pdf(paste("Analysis/Output/GO/WGCNA_goseq/GOterm_enrichment_plots/Mol_Function_Day7_heatmap.pdf", sep =''), width=28, height=45)
print(MF_D7_Plot)
dev.off()






############################################################################################################## 3
#  WCNA Day 14
BP_D14$Day_moduleColor <- paste(BP_D14$Day, BP_D14$moduleColor, sep = '_')
BP_D14$Day_moduleColor <- factor(BP_D14$Day_moduleColor, levels = c("Day14_brown", "Day14_black", "Day14_pink", "Day14_magenta"))
BP_D14_Plot <-ggplot(data = BP_D14, aes(x = Day_moduleColor, y = forcats::fct_rev(term))) + 
  geom_tile(aes(fill= -log10(over_represented_pvalue), width = 1)) + 
  scale_fill_gradient(low = "skyblue", high = "navyblue",limits=c(0, (ceiling(max(-log10(BP_D14$over_represented_pvalue))))), breaks=seq(0, (ceiling(max(-log10(BP_D14$over_represented_pvalue)))),by=2))  +
  facet_grid(~module_day, labeller = label_wrap_gen(width = 10, multi_line = TRUE))+
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                     strip.text.y = element_text(angle=0, size = 8, face = "bold"),
                     strip.text.x = element_text(size = 8, face = "bold"),
                     axis.title.x = element_blank(),
                     axis.title.y = element_text(size=8),
                     axis.text = element_text(size = 8), legend.position = "right",
                     plot.margin = unit(c(0,1,0,0.25), "cm"))+
  ggtitle('Biological Process: Day 14 WGCNA') +
  # facet_wrap(~ (substr(slim_term, 1,30)))
  facet_wrap((substr(slim_term, 1,30)) ~., scales="free_y", ncol= 1, strip.position="right")

pdf(paste("Analysis/Output/GO/WGCNA_goseq/GOterm_enrichment_plots/Biol_Process_Day14_heatmap.pdf", sep =''), width=15, height=25)
print(BP_D14_Plot)
dev.off()

MF_D14$Day_moduleColor <- paste(MF_D14$Day, MF_D14$moduleColor, sep = '_')
MF_D14$Day_moduleColor <- factor(MF_D14$Day_moduleColor, levels = c("Day14_brown", "Day14_black", "Day14_pink", "Day14_magenta"))
MF_D14_Plot <-ggplot(data = MF_D14, aes(x = Day_moduleColor, y = forcats::fct_rev(term))) + 
  geom_tile(aes(fill= -log10(over_represented_pvalue), width = 1)) + 
  scale_fill_gradient(low = "skyblue", high = "navyblue",limits=c(0, (ceiling(max(-log10(MF_D14$over_represented_pvalue))))), breaks=seq(0, (ceiling(max(-log10(MF_D14$over_represented_pvalue)))),by=2))  +
  facet_grid(~module_day, labeller = label_wrap_gen(width = 10, multi_line = TRUE))+
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                     strip.text.y = element_text(angle=0, size = 8, face = "bold"),
                     strip.text.x = element_text(size = 8, face = "bold"),
                     axis.title.x = element_blank(),
                     axis.title.y = element_text(size=8),
                     axis.text = element_text(size = 8), legend.position = "right",
                     plot.margin = unit(c(0,1,0,0.25), "cm"))+
  ggtitle('Molecular Function: Primary Moderate Effect Day 14') +
  # facet_wrap(~ (substr(slim_term, 1,30)))
  facet_wrap((substr(slim_term, 1,30)) ~., scales="free_y", ncol= 1, strip.position="right")


pdf(paste("Analysis/Output/GO/WGCNA_goseq/GOterm_enrichment_plots/Mol_Function_Day14_heatmap.pdf", sep =''), width=28, height=45)
print(MF_D14_Plot)
dev.off()




############################################################################################################## 3
#  WCNA Day 21
BP_D21$Day_moduleColor <- paste(BP_D21$Day, BP_D21$moduleColor, sep = '_')
BP_D21$Day_moduleColor <- factor(BP_D21$Day_moduleColor, levels = c("Day21_blue", "Day21_yellow", "Day21_red", "Day21_black", "Day21_pink", "Day21_turquoise"))
BP_D21_Plot <-ggplot(data = BP_D21, aes(x = Day_moduleColor, y = forcats::fct_rev(term))) + 
  geom_tile(aes(fill= -log10(over_represented_pvalue), width = 1)) + 
  scale_fill_gradient(low = "skyblue", high = "navyblue",limits=c(0, (ceiling(max(-log10(BP_D21$over_represented_pvalue))))), breaks=seq(0, (ceiling(max(-log10(BP_D21$over_represented_pvalue)))),by=2))  +
  facet_grid(~module_day, labeller = label_wrap_gen(width = 10, multi_line = TRUE))+
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                     strip.text.y = element_text(angle=0, size = 8, face = "bold"),
                     strip.text.x = element_text(size = 8, face = "bold"),
                     axis.title.x = element_blank(),
                     axis.title.y = element_text(size=8),
                     axis.text = element_text(size = 8), legend.position = "right",
                     plot.margin = unit(c(0,1,0,0.25), "cm"))+
  ggtitle('Biological Process: Day 21 WGCNA') +
  # facet_wrap(~ (substr(slim_term, 1,30)))
  facet_wrap((substr(slim_term, 1,30)) ~., scales="free_y", ncol= 1, strip.position="right")


pdf(paste("Analysis/Output/GO/WGCNA_goseq/GOterm_enrichment_plots/Biol_Process_Day21_heatmap.pdf", sep =''), width=15, height=25)
print(BP_D21_Plot)
dev.off()


MF_D21$Day_moduleColor <- paste(MF_D21$Day, MF_D21$moduleColor, sep = '_')
MF_D21$Day_moduleColor <- factor(MF_D21$Day_moduleColor, levels = c("Day21_magenta","Day21_blue", "Day21_yellow", "Day21_red", "Day21_black", "Day21_pink", "Day21_turquoise"))
MF_D21_Plot <-ggplot(data = MF_D21, aes(x = Day_moduleColor, y = forcats::fct_rev(term))) + 
  geom_tile(aes(fill= -log10(over_represented_pvalue), width = 1)) + 
  scale_fill_gradient(low = "skyblue", high = "navyblue",limits=c(0, (ceiling(max(-log10(MF_D21$over_represented_pvalue))))), breaks=seq(0, (ceiling(max(-log10(MF_D21$over_represented_pvalue)))),by=2))  +
  facet_grid(~module_day, labeller = label_wrap_gen(width = 10, multi_line = TRUE))+
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                     strip.text.y = element_text(angle=0, size = 8, face = "bold"),
                     strip.text.x = element_text(size = 8, face = "bold"),
                     axis.title.x = element_blank(),
                     axis.title.y = element_text(size=8),
                     axis.text = element_text(size = 8), legend.position = "right",
                     plot.margin = unit(c(0,1,0,0.25), "cm"))+
  ggtitle('Molecular Function: Primary Moderate Effect Day 21') +
  # facet_wrap(~ (substr(slim_term, 1,30)))
  facet_wrap((substr(slim_term, 1,30)) ~., scales="free_y", ncol= 1, strip.position="right")


pdf(paste("Analysis/Output/GO/WGCNA_goseq/GOterm_enrichment_plots/Mol_Function_Day21_heatmap.pdf", sep =''), width=15, height=40)
print(MF_D21_Plot)
dev.off()


############################################################################################################## 3
# moderate WCNA module effects day 7 yellow, day 14 black, day 21 yellow
BP_Mod$Day_moduleColor <- paste(BP_Mod$Day, BP_Mod$moduleColor, sep = '_')
BP_Mod$Day_moduleColor <- factor(BP_Mod$Day_moduleColor, levels = c("Day7_yellow", "Day14_black", "Day21_yellow"))
BP_Mod_Plot <-ggplot(data = BP_Mod, aes(x = Day_moduleColor, y = forcats::fct_rev(term))) + 
  geom_tile(aes(fill= -log10(over_represented_pvalue), width = 1)) + 
  scale_fill_gradient(low = "skyblue", high = "navyblue",limits=c(0, (ceiling(max(-log10(BP_Mod$over_represented_pvalue))))), breaks=seq(0, (ceiling(max(-log10(BP_Mod$over_represented_pvalue)))),by=2))  +
  facet_grid(~module_day, labeller = label_wrap_gen(width = 10, multi_line = TRUE))+
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                     strip.text.y = element_text(angle=0, size = 8, face = "bold"),
                     strip.text.x = element_text(size = 8, face = "bold"),
                     axis.title.x = element_blank(),
                     axis.title.y = element_text(size=8),
                     axis.text = element_text(size = 8), legend.position = "right",
                     plot.margin = unit(c(0,1,0,0.25), "cm"))+
  ggtitle('Biological Process: Primary Moderate Effect WGCNA') +
  # facet_wrap(~ (substr(slim_term, 1,30)))
  facet_wrap((substr(slim_term, 1,30)) ~., scales="free_y", ncol= 1, strip.position="right")





BP_Mod_Plot_condenced <- BP_Mod %>%  dplyr::filter(slim_term %in% c('response to stress', 'cell motility', 
                                                                    'signal transduction' ,'immune system process' ,
                                                                    'cellular nitrogen compound metaboli...', 'cell death',
                                                                    'cellular protein modification proce...')) %>% 
  ggplot(aes(x = Day_moduleColor, y = forcats::fct_rev(term))) + 
  geom_tile(aes(fill= -log10(over_represented_pvalue), width = 1)) + 
  scale_fill_gradient(low = "skyblue", high = "navyblue",limits=c(0, (ceiling(max(-log10(BP_Mod$over_represented_pvalue))))), breaks=seq(0, (ceiling(max(-log10(BP_Mod$over_represented_pvalue)))),by=2))  +
  facet_grid(~module_day, labeller = label_wrap_gen(width = 10, multi_line = TRUE))+
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                     strip.text.y = element_text(angle=0, size = 8, face = "bold"),
                     strip.text.x = element_text(size = 8, face = "bold"),
                     axis.title.x = element_blank(),
                     axis.title.y = element_text(size=8),
                     axis.text = element_text(size = 8), legend.position = "right",
                     plot.margin = unit(c(0,1,0,0.25), "cm"))+
  ggtitle('Biological Process: Primary Moderate Effect WGCNA') +
  # facet_wrap(~ (substr(slim_term, 1,30)))
  facet_wrap((substr(slim_term, 1,30)) ~., scales="free_y", ncol= 1, strip.position="right")







pdf(paste("Analysis/Output/GO/WGCNA_goseq/GOterm_enrichment_plots/Biol_Process_Moderate_heatmap.pdf", sep =''), width=10, height=18)
print(BP_Mod_Plot)
dev.off()

pdf(paste("Analysis/Output/GO/WGCNA_goseq/GOterm_enrichment_plots/Biol_Process_Moderate_heatmap_CONDENCED.pdf", sep =''), width=10, height=5)
print(BP_Mod_Plot_condenced)
dev.off()


MF_Mod$term    <- sub(" \\with incorporation.*", "", MF_Mod$term) # ommit the really long term to assist plotting

MF_Mod$Day_moduleColor <- paste(MF_Mod$Day, MF_Mod$moduleColor, sep = '_')
MF_Mod$Day_moduleColor <- factor(MF_Mod$Day_moduleColor, levels = c("Day7_yellow", "Day14_black", "Day21_yellow"))
MF_Mod_Plot <-ggplot(data = MF_Mod, aes(x = Day_moduleColor, y = forcats::fct_rev(term))) + 
  geom_tile(aes(fill= -log10(over_represented_pvalue), width = 1)) + 
  scale_fill_gradient(low = "skyblue", high = "navyblue",limits=c(0, (ceiling(max(-log10(MF_Mod$over_represented_pvalue))))), breaks=seq(0, (ceiling(max(-log10(MF_Mod$over_represented_pvalue)))),by=2))  +
  facet_grid(~module_day, labeller = label_wrap_gen(width = 10, multi_line = TRUE))+
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                     strip.text.y = element_text(angle=0, size = 8, face = "bold"),
                     strip.text.x = element_text(size = 8, face = "bold"),
                     axis.title.x = element_blank(),
                     axis.title.y = element_text(size=8),
                     axis.text = element_text(size = 8), legend.position = "right",
                     plot.margin = unit(c(0,1,0,0.25), "cm"))+
  ggtitle('Molecular Function: Primary Moderate Effect WGCNA') +
  # facet_wrap(~ (substr(slim_term, 1,30)))
  facet_wrap((substr(slim_term, 1,30)) ~., scales="free_y", ncol= 1, strip.position="right")


MF_Mod_Plot_condenced <- MF_Mod %>%  dplyr::filter(slim_term %in% c('oxidoreductase activity', 'transmembrane transporter activity', 
                                                                    'methyltransferase activity' ,'ion binding' ,
                                                                    'enzyme binding', 'kinase activity', 'transcription factor binding', 
                                                                    'cytoskeletal protein binding', 'RNA binding')) %>% 
  ggplot(aes(x = Day_moduleColor, y = forcats::fct_rev(term))) + 
  geom_tile(aes(fill= -log10(over_represented_pvalue), width = 1)) + 
  scale_fill_gradient(low = "skyblue", high = "navyblue",limits=c(0, (ceiling(max(-log10(MF_Mod$over_represented_pvalue))))), breaks=seq(0, (ceiling(max(-log10(MF_Mod$over_represented_pvalue)))),by=2))  +
  facet_grid(~module_day, labeller = label_wrap_gen(width = 10, multi_line = TRUE))+
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                     strip.text.y = element_text(angle=0, size = 8, face = "bold"),
                     strip.text.x = element_text(size = 8, face = "bold"),
                     axis.title.x = element_blank(),
                     axis.title.y = element_text(size=8),
                     axis.text = element_text(size = 8), legend.position = "right",
                     plot.margin = unit(c(0,1,0,0.25), "cm"))+
  ggtitle('Biological Process: Primary Moderate Effect WGCNA') +
  # facet_wrap(~ (substr(slim_term, 1,30)))
  facet_wrap((substr(slim_term, 1,30)) ~., scales="free_y", ncol= 1, strip.position="right")
MF_Mod_Plot_condenced


pdf(paste("Analysis/Output/GO/WGCNA_goseq/GOterm_enrichment_plots/Mol_Function_Moderate_heatmap.pdf", sep =''), width=15, height=35)
print(MF_Mod_Plot)
dev.off()

pdf(paste("Analysis/Output/GO/WGCNA_goseq/GOterm_enrichment_plots/Mol_Function_Moderate_heatmap_CONDENCED.pdf", sep =''), width=12, height=12)
print(MF_Mod_Plot_condenced)
dev.off()












############################################################################################################## 3
# ambient WCNA module effects day 7 brown, day 14 brown, day 21 blue and magenta
BP_Amb$Day_moduleColor <- paste(BP_Amb$Day, BP_Amb$moduleColor, sep = '_')
BP_Amb$Day_moduleColor <- factor(BP_Amb$Day_moduleColor, levels = c("Day7_brown", "Day14_brown", "Day21_blue"))
BP_Amb_Plot <-ggplot(data = BP_Amb, aes(x = Day_moduleColor, y = forcats::fct_rev(term))) + 
  geom_tile(aes(fill= -log10(over_represented_pvalue), width = 1)) + 
  scale_fill_gradient(low = "skyblue", high = "navyblue",limits=c(0, (ceiling(max(-log10(BP_Amb$over_represented_pvalue))))), breaks=seq(0, (ceiling(max(-log10(BP_Amb$over_represented_pvalue)))),by=2))  +
  facet_grid(~module_day, labeller = label_wrap_gen(width = 10, multi_line = TRUE))+
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                     strip.text.y = element_text(angle=0, size = 8, face = "bold"),
                     strip.text.x = element_text(size = 8, face = "bold"),
                     axis.title.x = element_blank(),
                     axis.title.y = element_text(size=8),
                     axis.text = element_text(size = 8), legend.position = "right",
                     plot.margin = unit(c(0,1,0,0.25), "cm"))+
  ggtitle('Biological Process: Primary Ambient Effect WGCNA') +
 # facet_wrap(~ (substr(slim_term, 1,30)))
  facet_wrap((substr(slim_term, 1,30)) ~., scales="free_y", ncol= 1, strip.position="right")


BP_Amb_Plot_2 <- BP_Amb %>%  dplyr::filter(slim_term %in% c('protein transport', 'transport', 
                                                            'vesicle-mediated transport','catabolic process',
                                                            'immune system process','response to stress')) %>% 
  ggplot( aes(x = Day_moduleColor, y = forcats::fct_rev(term))) + 
  geom_tile(aes(fill= -log10(over_represented_pvalue), width = 1)) + 
  scale_fill_gradient(low = "skyblue", high = "navyblue",limits=c(0, (ceiling(max(-log10(BP_Amb$over_represented_pvalue))))), breaks=seq(0, (ceiling(max(-log10(BP_Amb$over_represented_pvalue)))),by=2))  +
  facet_grid(~module_day, labeller = label_wrap_gen(width = 10, multi_line = TRUE))+
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                     strip.text.y = element_text(angle=0, size = 8, face = "bold"),
                     strip.text.x = element_text(size = 8, face = "bold"),
                     axis.title.x = element_blank(),
                     axis.title.y = element_text(size=8),
                     axis.text = element_text(size = 8), legend.position = "right",
                     plot.margin = unit(c(0,1,0,0.25), "cm"))+
  ggtitle('Biological Process: Primary Ambient Effect WGCNA') +
  # facet_wrap(~ (substr(slim_term, 1,30)))
  facet_wrap((substr(slim_term, 1,30)) ~., scales="free_y", ncol= 1, strip.position="right")

pdf(paste("Analysis/Output/GO/WGCNA_goseq/GOterm_enrichment_plots/Biol_Process_Ambient_heatmap.pdf", sep =''), width=15, height=35)
print(BP_Amb_Plot)
dev.off()


pdf(paste("Analysis/Output/GO/WGCNA_goseq/GOterm_enrichment_plots/Biol_Process_Ambient_heatmap_CONDENCED.pdf", sep =''), width=10, height=7)
print(BP_Amb_Plot_2)
dev.off()




MF_Amb$Day_moduleColor <- paste(MF_Amb$Day, MF_Amb$moduleColor, sep = '_')
MF_Amb$Day_moduleColor <- factor(MF_Amb$Day_moduleColor, levels = c("Day7_brown", "Day14_brown", "Day21_blue", "Day21_magenta"))
MF_Amb_Plot <-ggplot(data = MF_Amb, aes(x = Day_moduleColor, y = forcats::fct_rev(term))) + 
  geom_tile(aes(fill= -log10(over_represented_pvalue), width = 1)) + 
  scale_fill_gradient(low = "skyblue", high = "navyblue",limits=c(0, (ceiling(max(-log10(MF_Amb$over_represented_pvalue))))), breaks=seq(0, (ceiling(max(-log10(MF_Amb$over_represented_pvalue)))),by=2))  +
  facet_grid(~module_day, labeller = label_wrap_gen(width = 10, multi_line = TRUE))+
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                     strip.text.y = element_text(angle=0, size = 8, face = "bold"),
                     strip.text.x = element_text(size = 8, face = "bold"),
                     axis.title.x = element_blank(),
                     axis.title.y = element_text(size=8),
                     axis.text = element_text(size = 8), legend.position = "right",
                     plot.margin = unit(c(0,1,0,0.25), "cm"))+
  ggtitle('Molecular Function: Primary Ambient Effect WGCNA') +
  # facet_wrap(~ (substr(slim_term, 1,30)))
  facet_wrap((substr(slim_term, 1,30)) ~., scales="free_y", ncol= 1, strip.position="right")


MF_Amb_Plot_2 <- MF_Amb %>%  dplyr::filter(slim_term %in% c('enzyme building', 'ligase activity', 'enzyme regulator activity','enzyme binding',
                                                            'ion binding','lipid binding', 'oxidoreductase activity', 
                                                            'peptidase activity', 'transmembrane transporter activity')) %>% 
  ggplot( aes(x = Day_moduleColor, y = forcats::fct_rev(term))) + 
  geom_tile(aes(fill= -log10(over_represented_pvalue), width = 1)) + 
  scale_fill_gradient(low = "skyblue", high = "navyblue",limits=c(0, (ceiling(max(-log10(MF_Amb$over_represented_pvalue))))), breaks=seq(0, (ceiling(max(-log10(MF_Amb$over_represented_pvalue)))),by=2))  +
  facet_grid(~module_day, labeller = label_wrap_gen(width = 10, multi_line = TRUE))+
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                     strip.text.y = element_text(angle=0, size = 8, face = "bold"),
                     strip.text.x = element_text(size = 8, face = "bold"),
                     axis.title.x = element_blank(),
                     axis.title.y = element_text(size=8),
                     axis.text = element_text(size = 8), legend.position = "right",
                     plot.margin = unit(c(0,1,0,0.25), "cm"))+
  ggtitle('Molecular Function: Primary Ambient Effect WGCNA') +
  # facet_wrap(~ (substr(slim_term, 1,30)))
  facet_wrap((substr(slim_term, 1,30)) ~., scales="free_y", ncol= 1, strip.position="right")


pdf(paste("Analysis/Output/GO/WGCNA_goseq/GOterm_enrichment_plots/Mol_Function_Ambient_heatmap.pdf", sep =''), width=15, height=35)
print(MF_Amb_Plot)
dev.off()

pdf(paste("Analysis/Output/GO/WGCNA_goseq/GOterm_enrichment_plots/Mol_Function_Ambient_heatmap_CONDENCED.pdf", sep =''), width=13, height=16)
print(MF_Amb_Plot_2)
dev.off()






#===================================================================================================
#
#
#  Youre done! Now time to load in those final GOslim csvs and plot!!!!
# 
#   PLOTTING
#
#===================================================================================================
# load the output from the previous for loop
# BIOLOGICAL PROCESS  GO  SLIMS
View(d7_slimBP_brownModule)
d7_slimBP_brownModule   <- read.csv("Analysis/Output/GO/WGCNA_goseq/subseq_treatments_all/Day7/GOslim_BiolProc_brownModule.csv")
d7_slimBP_greenModule   <- read.csv("Analysis/Output/GO/WGCNA_goseq/subseq_treatments_all/Day7/GOslim_BiolProc_greenModule.csv")
d7_slimBP_yellowModule  <- read.csv("Analysis/Output/GO/WGCNA_goseq/subseq_treatments_all/Day7/GOslim_BiolProc_yellowModule.csv")

d14_slimBP_blackModule   <- read.csv("Analysis/Output/GO/WGCNA_goseq/subseq_treatments_all/Day14/GOslim_BiolProc_blackModule.csv")
d14_slimBP_brownModule   <- read.csv("Analysis/Output/GO/WGCNA_goseq/subseq_treatments_all/Day14/GOslim_BiolProc_brownModule.csv")
d14_slimBP_magentaModule <- read.csv("Analysis/Output/GO/WGCNA_goseq/subseq_treatments_all/Day14/GOslim_BiolProc_magentaModule.csv")
d14_slimBP_pinkModule    <- read.csv("Analysis/Output/GO/WGCNA_goseq/subseq_treatments_all/Day14/GOslim_BiolProc_pinkModule.csv")

d21_slimBP_blackModule       <- read.csv("Analysis/Output/GO/WGCNA_goseq/subseq_treatments_all/Day21/GOslim_BiolProc_blackModule.csv")
d21_slimBP_blueModule        <- read.csv("Analysis/Output/GO/WGCNA_goseq/subseq_treatments_all/Day21/GOslim_BiolProc_blueModule.csv")
d21_slimBP_magentaModule     <- read.csv("Analysis/Output/GO/WGCNA_goseq/subseq_treatments_all/Day21/GOslim_BiolProc_magentaModule.csv")
d21_slimBP_pinkModule        <- read.csv("Analysis/Output/GO/WGCNA_goseq/subseq_treatments_all/Day21/GOslim_BiolProc_pinkModule.csv")
d21_slimBP_redModule         <- read.csv("Analysis/Output/GO/WGCNA_goseq/subseq_treatments_all/Day21/GOslim_BiolProc_redModule.csv")
d21_slimBP_yellowModule      <- read.csv("Analysis/Output/GO/WGCNA_goseq/subseq_treatments_all/Day21/GOslim_BiolProc_yellowModule.csv")
d21_slimBP_turquoiseModule   <- read.csv("Analysis/Output/GO/WGCNA_goseq/subseq_treatments_all/Day21/GOslim_BiolProc_turquoiseModule.csv")

# MOLECULAR FUNCTION GO  SLIMS
d7_slimMF_brownModule   <- read.csv("Analysis/Output/GO/WGCNA_goseq/subseq_treatments_all/Day7/GOslim_MolFunction_brownModule.csv")
d7_slimMF_greenModule   <- read.csv("Analysis/Output/GO/WGCNA_goseq/subseq_treatments_all/Day7/GOslim_MolFunction_greenModule.csv")
d7_slimMF_yellowModule  <- read.csv("Analysis/Output/GO/WGCNA_goseq/subseq_treatments_all/Day7/GOslim_MolFunction_yellowModule.csv")

d14_slimMF_blackModule   <- read.csv("Analysis/Output/GO/WGCNA_goseq/subseq_treatments_all/Day14/GOslim_MolFunction_blackModule.csv")
d14_slimMF_brownModule   <- read.csv("Analysis/Output/GO/WGCNA_goseq/subseq_treatments_all/Day14/GOslim_MolFunction_brownModule.csv")
d14_slimMF_magentaModule <- read.csv("Analysis/Output/GO/WGCNA_goseq/subseq_treatments_all/Day14/GOslim_MolFunction_magentaModule.csv")
d14_slimMF_pinkModule    <- read.csv("Analysis/Output/GO/WGCNA_goseq/subseq_treatments_all/Day14/GOslim_MolFunction_pinkModule.csv")

d21_slimMF_blackModule       <- read.csv("Analysis/Output/GO/WGCNA_goseq/subseq_treatments_all/Day21/GOslim_MolFunction_blackModule.csv")
d21_slimMF_blueModule        <- read.csv("Analysis/Output/GO/WGCNA_goseq/subseq_treatments_all/Day21/GOslim_MolFunction_blueModule.csv")
d21_slimMF_magentaModule     <- read.csv("Analysis/Output/GO/WGCNA_goseq/subseq_treatments_all/Day21/GOslim_MolFunction_magentaModule.csv")
d21_slimMF_pinkModule        <- read.csv("Analysis/Output/GO/WGCNA_goseq/subseq_treatments_all/Day21/GOslim_MolFunction_pinkModule.csv")
d21_slimMF_redModule         <- read.csv("Analysis/Output/GO/WGCNA_goseq/subseq_treatments_all/Day21/GOslim_MolFunction_redModule.csv")
d21_slimMF_yellowModule      <- read.csv("Analysis/Output/GO/WGCNA_goseq/subseq_treatments_all/Day21/GOslim_MolFunction_yellowModule.csv")
d21_slimMF_turquoiseModule   <- read.csv("Analysis/Output/GO/WGCNA_goseq/subseq_treatments_all/Day21/GOslim_MolFunction_turquoiseModule.csv")

# BIND DATA TO FACET THE PLOTS    ---------------------------------------------------------- #
# (1) By sampling date 0 Day 7 14 and 21 separately
BP_D7  <- rbind(d7_slimBP_brownModule, d7_slimBP_greenModule, d7_slimBP_yellowModule)
BP_D14 <- rbind(d14_slimBP_blackModule, d14_slimBP_brownModule, d14_slimBP_magentaModule, d14_slimBP_pinkModule)
BP_D21 <- rbind(d21_slimBP_blackModule, d21_slimBP_blueModule, d21_slimBP_magentaModule, d21_slimBP_pinkModule, d21_slimBP_redModule, d21_slimBP_yellowModule,d21_slimBP_turquoiseModule)

MF_D7  <- rbind(d7_slimMF_brownModule, d7_slimMF_greenModule, d7_slimMF_yellowModule)
MF_D14 <- rbind(d14_slimMF_blackModule, d14_slimMF_brownModule, d14_slimMF_magentaModule, d14_slimMF_pinkModule)
MF_D21 <- rbind(d21_slimMF_blackModule, d21_slimMF_blueModule, d21_slimMF_magentaModule, d21_slimMF_pinkModule, d21_slimMF_redModule, d21_slimMF_yellowModule, d21_slimMF_turquoiseModule)

# (2) By trend in vst expression pattern ACROSS all days 
# (2.A) Priming Effect AMBIENT > MODERATE 
BP_Amb <- rbind(d7_slimBP_brownModule, d14_slimBP_brownModule, d21_slimBP_magentaModule, d21_slimBP_blueModule) # merge the data by binding rows 
MF_Amb <- rbind(d7_slimMF_brownModule, d14_slimMF_brownModule, d21_slimMF_magentaModule, d21_slimMF_blueModule) # merge the data by binding rows 

# (2.B) Priming Effect MODERATE > AMBIENT 
BP_Mod <- rbind(d7_slimBP_yellowModule, d14_slimBP_blackModule, d21_slimBP_yellowModule) # merge the data by binding rows 
MF_Mod <- rbind(d7_slimMF_yellowModule, d14_slimMF_blackModule, d21_slimMF_yellowModule) # merge the data by binding rows 





# PLOTTING --------------------------------------------------------------------------------------------------- #
# (1) By sampling date 0 Day 7 14 and 21 separately




#  DAY 7 



BP_D7$Ont        <- "BP" # create a new common clumn to call in the plot 
BP_D7$module_day <- factor(BP_D7$module_day, levels = c("Day7_brown", "Day7_yellow", "Day7_green"))# reorder the facotr level for the facet wrap plot 
# BP_D7_filtered <- BP_D7 %>%  dplyr::filter(Gene_count >= 10) # ommit all terms gene counts 51
# BP_D7_filtered$slim_term <- factor(BP_D7_filtered$slim_term) # make slim term alphabetical for plotting
BP_D7$slim_term  <- factor(BP_D7$slim_term) # make slim term alphabetical for plotting

BP_D7_Plot <-ggplot(data = BP_D7, aes(x = Ont, y = forcats::fct_rev(slim_term))) + 
              geom_tile(aes(fill=Gene_count, width = 1)) + 
              scale_fill_gradient(low = "grey95", high = "grey10",limits=c(0, 100), breaks=seq(0,100,by=25))  +
              facet_grid(~module_day, labeller = label_wrap_gen(width = 10, multi_line = TRUE))+
              theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                                 panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                                 strip.text.y = element_text(angle=0, size = 8, face = "bold"),
                                 strip.text.x = element_text(size = 8, face = "bold"),
                                 axis.title.x = element_blank(),
                                 axis.title.y = element_text(size=8),
                                 axis.text = element_text(size = 8), legend.position = "right",
                                 plot.margin = unit(c(0,1,0,0.25), "cm"))+
              ggtitle('GOslim Biological Process: WGCNA Day 7')

MF_D7$Ont <- "MF" # create a new common clumn to call in the plot 
MF_D7$module_day <- factor(MF_D7$module_day, levels = c("Day7_brown", "Day7_yellow", "Day7_green"))# reorder the facotr level for the facet wrap plot 
# MF_D7_filtered <- MF_D7 %>%  dplyr::filter(Gene_count >= 3) # ommit all with gene counts <1
# MF_D7_filtered$slim_term <- factor(MF_D7_filtered$slim_term) # make slim term alphabetical for plotting
MF_D7$slim_term <- factor(MF_D7$slim_term) # make slim term alphabetical for plotting

MF_D7_Plot <-ggplot(data = MF_D7, aes(x = Ont, y = forcats::fct_rev(slim_term))) + 
              geom_tile(aes(fill=Gene_count, width = 1)) + 
              scale_fill_gradient(low = "grey95", high = "grey10",limits=c(0, 100), breaks=seq(0,100,by=25))  +
              facet_grid(~module_day, labeller = label_wrap_gen(width = 10, multi_line = TRUE))+
              theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                                 panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                                 strip.text.y = element_text(angle=0, size = 8, face = "bold"),
                                 strip.text.x = element_text(size = 8, face = "bold"),
                                 axis.title.x = element_blank(),
                                 axis.title.y = element_text(size=8),
                                 axis.text = element_text(size = 8), legend.position = "right",
                                 plot.margin = unit(c(0,1,0,0.25), "cm"))+
             ggtitle('GOslim Molecular Function: WGCNA Day 7')




# DAY 14




BP_D14$Ont <- "BP" # create a new common clumn to call in the plot 
BP_D14$module_day <- factor(BP_D14$module_day, levels = c("Day14_brown", "Day14_black", "Day14_pink", "Day14_magenta"))# reorder the facotr level for the facet wrap plot 
# BP_D14_filtered <- BP_D14 %>%  dplyr::filter(Gene_count >= 10) # ommit all with gene counts <1
# BP_D14_filtered$slim_term <- factor(BP_D14_filtered$slim_term) # make slim term alphabetical for plotting
BP_D14$slim_term <- factor(BP_D14$slim_term) # make slim term alphabetical for plotting

BP_D14_Plot <-ggplot(data = BP_D14, aes(x = Ont, y = forcats::fct_rev(slim_term))) + 
                geom_tile(aes(fill=Gene_count, width = 1)) + 
                scale_fill_gradient(low = "grey95", high = "grey10",limits=c(0, 100), breaks=seq(0,100,by=25)) +
                facet_grid(~module_day, labeller = label_wrap_gen(width = 10, multi_line = TRUE))+
                theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                                   panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                                   strip.text.y = element_text(angle=0, size = 8, face = "bold"),
                                   strip.text.x = element_text(size = 8, face = "bold"),
                                   axis.title.x = element_blank(),
                                   axis.title.y = element_text(size=8),
                                   axis.text = element_text(size = 8), legend.position = "right",
                                   plot.margin = unit(c(0,1,0,0.25), "cm"))+
                ggtitle('GOslim Biological Process: WGCNA Day 14')

MF_D14$Ont <- "MF" # create a new common clumn to call in the plot 
MF_D14$module_day <- factor(MF_D14$module_day, levels = c("Day14_brown", "Day14_black", "Day14_pink", "Day14_magenta"))# reorder the facotr level for the facet wrap plot 
# MF_D14_filtered <- MF_D14 %>%  dplyr::filter(Gene_count >= 3) # ommit all with gene counts <1
# MF_D14_filtered$slim_term <- factor(MF_D14_filtered$slim_term) # make slim term alphabetical for plotting
MF_D14$slim_term <- factor(MF_D14$slim_term) # make slim term alphabetical for plotting

MF_D14_Plot <-ggplot(data = MF_D14, aes(x = Ont, y = forcats::fct_rev(slim_term))) + 
              geom_tile(aes(fill=Gene_count, width = 1)) + 
              scale_fill_gradient(low = "grey95", high = "grey10",limits=c(0, 100), breaks=seq(0,100,by=25)) +
              facet_grid(~module_day, labeller = label_wrap_gen(width = 10, multi_line = TRUE))+
              theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                                 panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                                 strip.text.y = element_text(angle=0, size = 8, face = "bold"),
                                 strip.text.x = element_text(size = 8, face = "bold"),
                                 axis.title.x = element_blank(),
                                 axis.title.y = element_text(size=8),
                                 axis.text = element_text(size = 8), legend.position = "right",
                                 plot.margin = unit(c(0,1,0,0.25), "cm"))+
              ggtitle('GOslim Molecular Function: WGCNA Day 14')




# DAY 21




BP_D21$Ont <- "BP" # create a new common clumn to call in the plot 
BP_D21$module_day <- factor(BP_D21$module_day, levels = c("Day21_magenta",  "Day21_blue", "Day21_yellow",  "Day21_red", "Day21_black", "Day21_pink", "Day21_turquoise"))# reorder the facotr level for the facet wrap plot 
# BP_D21_filtered <- BP_D21 %>%  dplyr::filter(Gene_count >= 10) # ommit all with gene counts <1
# BP_D21_filtered$slim_term <- factor(BP_D21_filtered$slim_term) # make slim term alphabetical for plotting
BP_D21$slim_term <- factor(BP_D21$slim_term) # make slim term alphabetical for plotting

BP_D21_Plot <-ggplot(data = BP_D21, aes(x = Ont, y = forcats::fct_rev(slim_term))) + 
              geom_tile(aes(fill=Gene_count, width = 1)) + 
              scale_fill_gradient(low = "grey95", high = "grey10",limits=c(0, 250), breaks=seq(0,250,by=50)) +
              facet_grid(~module_day, labeller = label_wrap_gen(width = 10, multi_line = TRUE))+
              theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                                 panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                                 strip.text.y = element_text(angle=0, size = 8, face = "bold"),
                                 strip.text.x = element_text(size = 8, face = "bold"),
                                 axis.title.x = element_blank(),
                                 axis.title.y = element_text(size=8),
                                 axis.text = element_text(size = 8), legend.position = "right",
                                 plot.margin = unit(c(0,1,0,0.25), "cm"))+
              ggtitle('GOslim Biological Process: WGCNA Day 21')

MF_D21$Ont <- "MF" # create a new common clumn to call in the plot 
MF_D21$module_day <- factor(MF_D21$module_day, levels = c("Day21_magenta",  "Day21_blue", "Day21_yellow",  "Day21_red", "Day21_black", "Day21_pink", "Day21_turquoise"))# reorder the facotr level for the facet wrap plot 
# MF_D21_filtered <- MF_D21 %>%  dplyr::filter(Gene_count >= 3) # ommit all with gene counts <1
# MF_D21_filtered$slim_term <- factor(MF_D21_filtered$slim_term)  # make slim term alphabetical for plotting
MF_D21$slim_term <- factor(MF_D21$slim_term)  # make slim term alphabetical for plotting

MF_D21_Plot <-ggplot(data = MF_D21, aes(x = Ont, y = forcats::fct_rev(slim_term))) + 
              geom_tile(aes(fill=Gene_count, width = 1)) + 
              scale_fill_gradient(low = "grey95", high = "grey10",limits=c(0, 250), breaks=seq(0,250,by=50)) +
              facet_grid(~module_day, labeller = label_wrap_gen(width = 10, multi_line = TRUE))+
              theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                                 panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                                 strip.text.y = element_text(angle=0, size = 8, face = "bold"),
                                 strip.text.x = element_text(size = 8, face = "bold"),
                                 axis.title.x = element_blank(),
                                 axis.title.y = element_text(size=8),
                                 axis.text = element_text(size = 8), legend.position = "right",
                                 plot.margin = unit(c(0,1,0,0.25), "cm"))+
              #scale_y_discrete("slim_term", trans = 'reverse') +
              ggtitle('GOslim Molecular Function: WGCNA Day 21')



# SAVE PLOTS  ------------------------------------------------------------------------------------------------- #
# WGCNA Day 7 - all significant modules (correlated with treatment(s))
pdf(paste("Analysis/Output/GO/WGCNA_goseq/subseq_treatments_all/Day7/GOslim_day7.pdf", sep =''), width=18, height=8)
print(ggarrange(BP_D7_Plot, MF_D7_Plot,         
                 plotlist = NULL,
                 ncol = 2,
                 nrow = 1,
                 labels = NULL))
dev.off()     
# WGCNA Day 14 - all significant modules (correlated with treatment(s))
pdf(paste("Analysis/Output/GO/WGCNA_goseq/subseq_treatments_all/Day14/GOslim_day14.pdf", sep =''), width=18, height=8)
print(ggarrange(BP_D14_Plot, MF_D14_Plot,         
                plotlist = NULL,
                ncol = 2,
                nrow = 1,
                labels = NULL))
dev.off()     
# WGCNA Day 21 - all significant modules (correlated with treatment(s))
pdf(paste("Analysis/Output/GO/WGCNA_goseq/subseq_treatments_all/Day21/GOslim_day21.pdf", sep =''), width=20, height=8)
print(ggarrange(BP_D21_Plot, MF_D21_Plot,         
                plotlist = NULL,
                ncol = 2,
                nrow = 1,
                labels = NULL))
dev.off() 








# (2.A) Priming Effect AMBIENT > MODERATE 


# BP  - call the sig moduels with A > M main effects 
BP_Amb$Ont <- "BP" # create a new common clumn to call in the plot 
BP_Amb$module_day <- factor(BP_Amb$module_day, levels = c("Day7_brown", "Day14_brown", "Day21_blue", "Day21_magenta"))# reorder the facotr level for the facet wrap plot 
BP_Amb_filtered <- BP_Amb %>%  dplyr::filter(Gene_count > 1) # ommit all with gene counts <1
# MF   - call the sig moduels with A > M main effects 
MF_Amb$Ont <- "MF" # create a new common clumn to call in the plot 
MF_Amb$module_day <- factor(MF_Amb$module_day, levels = c("Day7_brown", "Day14_brown", "Day21_blue", "Day21_magenta")) # reorder the facotr level for the facet wrap plot 
MF_Amb_filtered <- MF_Amb %>%  dplyr::filter(Gene_count > 1) # ommit all with gene counts <1



# plots
BP_Amb_Plot <-ggplot(data = BP_Amb_filtered, aes(x = Ont, y = forcats::fct_rev(slim_term))) + 
              geom_tile(aes(fill=Gene_count, width = 1)) + 
              scale_fill_gradient(low = "grey95", high = "grey10",limits=c(0, 200), breaks=seq(0,200,by=50)) +
              facet_grid(~module_day, labeller = label_wrap_gen(width = 10, multi_line = TRUE))+
              theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                                 panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                                 strip.text.y = element_text(angle=0, size = 8, face = "bold"),
                                 strip.text.x = element_text(size = 8, face = "bold"),
                                 axis.title.x = element_blank(),
                                 axis.title.y = element_text(size=8),
                                 axis.text = element_text(size = 8), legend.position = "right",
                                 plot.margin = unit(c(0,1,0,0.25), "cm"))+
              ggtitle('GOslim Biological Process: WGCNA all modules Ambient > Moderate')

# plot MF
MF_Amb_Plot <-ggplot(data = MF_Amb_filtered, aes(x = Ont, y = forcats::fct_rev(slim_term))) + 
              geom_tile(aes(fill=Gene_count, width = 1)) + 
              scale_fill_gradient(low = "grey95", high = "grey10",limits=c(0, 200), breaks=seq(0,200,by=50)) +
              facet_grid(~module_day, labeller = label_wrap_gen(width = 10, multi_line = TRUE))+
              theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                                 panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                                 strip.text.y = element_text(angle=0, size = 8, face = "bold"),
                                 strip.text.x = element_text(size = 8, face = "bold"),
                                 axis.title.x = element_blank(),
                                 axis.title.y = element_text(size=8),
                                 axis.text = element_text(size = 8), legend.position = "right",
                                 plot.margin = unit(c(0,1,0,0.25), "cm"))+
              ggtitle('GOslim Molecular Function: WGCNA all modules Ambient > Moderate')



# SAVE PLOTS  ------------------------------------------------------------------------------------------------- #
pdf(paste("Analysis/Output/GO/WGCNA_goseq/GOslim_Modules_Ambient.pdf", sep =''), width=20, height=8)
print(ggarrange(BP_Amb_Plot, MF_Amb_Plot,         
                plotlist = NULL,
                ncol = 2,
                nrow = 1,
                labels = NULL))
dev.off() 









# (2.B) Priming Effect MODERATE > AMBIENT 








# BP  - call the sig moduels with A > M main effects 
BP_Mod$Ont <- "BP" # create a new common clumn to call in the plot 
BP_Mod$module_day <- factor(BP_Mod$module_day, levels = c("Day7_yellow", "Day14_black", "Day21_yellow"))# reorder the facotr level for the facet wrap plot 
BP_Mod_filtered <- BP_Mod %>%  dplyr::filter(Gene_count > 1) # ommit all with gene counts <1
# MF   - call the sig moduels with A > M main effects 
MF_Mod$Ont <- "MF" # create a new common clumn to call in the plot 
MF_Mod$module_day <- factor(MF_Mod$module_day, levels = c("Day7_yellow", "Day14_black", "Day21_yellow")) # reorder the facotr level for the facet wrap plot 
MF_Mod_filtered <- MF_Mod %>%  dplyr::filter(Gene_count > 1) # ommit all with gene counts <1



# plots
BP_Mod_Plot <-ggplot(data = BP_Mod, aes(x = Ont, y = forcats::fct_rev(slim_term))) + 
                geom_tile(aes(fill=Gene_count, width = 1)) + 
                scale_fill_gradient(low = "grey95", high = "grey10",limits=c(0, 150), breaks=seq(0,150,by=25)) +
                facet_grid(~module_day, labeller = label_wrap_gen(width = 10, multi_line = TRUE))+
                theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                                   panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                                   strip.text.y = element_text(angle=0, size = 8, face = "bold"),
                                   strip.text.x = element_text(size = 8, face = "bold"),
                                   axis.title.x = element_blank(),
                                   axis.title.y = element_text(size=8),
                                   axis.text = element_text(size = 8), legend.position = "right",
                                   plot.margin = unit(c(0,1,0,0.25), "cm"))+
                ggtitle('GOslim Biological Process: WGCNA all modules Modient > Ambient')

# plot MF
MF_Mod_Plot <-ggplot(data = MF_Mod, aes(x = Ont, y = forcats::fct_rev(slim_term))) + 
                geom_tile(aes(fill=Gene_count, width = 1)) + 
                #scale_fill_continuous(limits=c(0, 200), breaks=seq(0,200,by=25)) +
                scale_fill_gradient(low = "grey95", high = "grey10",limits=c(0, 150), breaks=seq(0,150,by=25)) +
                facet_grid(~module_day, labeller = label_wrap_gen(width = 10, multi_line = TRUE))+
                theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                                   panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                                   strip.text.y = element_text(angle=0, size = 8, face = "bold"),
                                   strip.text.x = element_text(size = 8, face = "bold"),
                                   axis.title.x = element_blank(),
                                   axis.title.y = element_text(size=8),
                                   axis.text = element_text(size = 8), legend.position = "right",
                                   plot.margin = unit(c(0,1,0,0.25), "cm"))+
                ggtitle('GOslim Molecular Function: WGCNA all modules Modient > Ambient')



# SAVE PLOTS  ------------------------------------------------------------------------------------------------- #
pdf(paste("Analysis/Output/GO/WGCNA_goseq/GOslim_Modules_Moderate.pdf", sep =''), width=20, height=8)
print(ggarrange(BP_Mod_Plot, MF_Mod_Plot,         
                plotlist = NULL,
                ncol = 2,
                nrow = 1,
                labels = NULL))
dev.off() 






#===================================================================================================
#
#
#  After lookinf over the visuals - I have dtermined several lines of interest to explore further...
# 
#   MORE PLOTTING 
#
#===================================================================================================
library(dplyr)
library(VennDiagram)
library(ggVennDiagram)
library(ggvenn)
library(gridExtra)
library(ggpubr)
library(reshape2)
library(DESeq2)
library(stringr)
# =================================================================================== #
#
#  LOAD DATA AND  DATA PREP FOR PLOTS
#
# =================================================================================== # 
# LOAD DATA
# BIOLOGICAL PROCESS  GO  SLIMS
d7_slimBP_brownModule   <- read.csv("Analysis/Output/GO/WGCNA_goseq/subseq_treatments_all/Day7/GOslim_BiolProc_brownModule.csv")
d7_slimBP_greenModule   <- read.csv("Analysis/Output/GO/WGCNA_goseq/subseq_treatments_all/Day7/GOslim_BiolProc_greenModule.csv")
d7_slimBP_yellowModule  <- read.csv("Analysis/Output/GO/WGCNA_goseq/subseq_treatments_all/Day7/GOslim_BiolProc_yellowModule.csv")

d14_slimBP_blackModule   <- read.csv("Analysis/Output/GO/WGCNA_goseq/subseq_treatments_all/Day14/GOslim_BiolProc_blackModule.csv")
d14_slimBP_brownModule   <- read.csv("Analysis/Output/GO/WGCNA_goseq/subseq_treatments_all/Day14/GOslim_BiolProc_brownModule.csv")
d14_slimBP_magentaModule <- read.csv("Analysis/Output/GO/WGCNA_goseq/subseq_treatments_all/Day14/GOslim_BiolProc_magentaModule.csv")
d14_slimBP_pinkModule    <- read.csv("Analysis/Output/GO/WGCNA_goseq/subseq_treatments_all/Day14/GOslim_BiolProc_pinkModule.csv")

d21_slimBP_blackModule       <- read.csv("Analysis/Output/GO/WGCNA_goseq/subseq_treatments_all/Day21/GOslim_BiolProc_blackModule.csv")
d21_slimBP_blueModule        <- read.csv("Analysis/Output/GO/WGCNA_goseq/subseq_treatments_all/Day21/GOslim_BiolProc_blueModule.csv")
d21_slimBP_magentaModule     <- read.csv("Analysis/Output/GO/WGCNA_goseq/subseq_treatments_all/Day21/GOslim_BiolProc_magentaModule.csv")
d21_slimBP_pinkModule        <- read.csv("Analysis/Output/GO/WGCNA_goseq/subseq_treatments_all/Day21/GOslim_BiolProc_pinkModule.csv")
d21_slimBP_redModule         <- read.csv("Analysis/Output/GO/WGCNA_goseq/subseq_treatments_all/Day21/GOslim_BiolProc_redModule.csv")
d21_slimBP_yellowModule      <- read.csv("Analysis/Output/GO/WGCNA_goseq/subseq_treatments_all/Day21/GOslim_BiolProc_yellowModule.csv")
d21_slimBP_turquoiseModule   <- read.csv("Analysis/Output/GO/WGCNA_goseq/subseq_treatments_all/Day21/GOslim_BiolProc_turquoiseModule.csv")

# MOLECULAR FUNCTION GO  SLIMS
d7_slimMF_brownModule   <- read.csv("Analysis/Output/GO/WGCNA_goseq/subseq_treatments_all/Day7/GOslim_MolFunction_brownModule.csv")
d7_slimMF_greenModule   <- read.csv("Analysis/Output/GO/WGCNA_goseq/subseq_treatments_all/Day7/GOslim_MolFunction_greenModule.csv")
d7_slimMF_yellowModule  <- read.csv("Analysis/Output/GO/WGCNA_goseq/subseq_treatments_all/Day7/GOslim_MolFunction_yellowModule.csv")

d14_slimMF_blackModule   <- read.csv("Analysis/Output/GO/WGCNA_goseq/subseq_treatments_all/Day14/GOslim_MolFunction_blackModule.csv")
d14_slimMF_brownModule   <- read.csv("Analysis/Output/GO/WGCNA_goseq/subseq_treatments_all/Day14/GOslim_MolFunction_brownModule.csv")
d14_slimMF_magentaModule <- read.csv("Analysis/Output/GO/WGCNA_goseq/subseq_treatments_all/Day14/GOslim_MolFunction_magentaModule.csv")
d14_slimMF_pinkModule    <- read.csv("Analysis/Output/GO/WGCNA_goseq/subseq_treatments_all/Day14/GOslim_MolFunction_pinkModule.csv")

d21_slimMF_blackModule       <- read.csv("Analysis/Output/GO/WGCNA_goseq/subseq_treatments_all/Day21/GOslim_MolFunction_blackModule.csv")
d21_slimMF_blueModule        <- read.csv("Analysis/Output/GO/WGCNA_goseq/subseq_treatments_all/Day21/GOslim_MolFunction_blueModule.csv")
d21_slimMF_magentaModule     <- read.csv("Analysis/Output/GO/WGCNA_goseq/subseq_treatments_all/Day21/GOslim_MolFunction_magentaModule.csv")
d21_slimMF_pinkModule        <- read.csv("Analysis/Output/GO/WGCNA_goseq/subseq_treatments_all/Day21/GOslim_MolFunction_pinkModule.csv")
d21_slimMF_redModule         <- read.csv("Analysis/Output/GO/WGCNA_goseq/subseq_treatments_all/Day21/GOslim_MolFunction_redModule.csv")
d21_slimMF_yellowModule      <- read.csv("Analysis/Output/GO/WGCNA_goseq/subseq_treatments_all/Day21/GOslim_MolFunction_yellowModule.csv")
d21_slimMF_turquoiseModule   <- read.csv("Analysis/Output/GO/WGCNA_goseq/subseq_treatments_all/Day21/GOslim_MolFunction_turquoiseModule.csv")

# Ambient and Moderate effect modules 
# Ambient 
BP_Amb <- rbind(d7_slimBP_brownModule, d14_slimBP_brownModule, d21_slimBP_magentaModule, d21_slimBP_blueModule) # merge the data by binding rows 
MF_Amb <- rbind(d7_slimMF_brownModule, d14_slimMF_brownModule, d21_slimMF_magentaModule, d21_slimMF_blueModule) # merge the data by binding rows 
# Moderate
BP_Mod <- rbind(d7_slimBP_yellowModule, d14_slimBP_blackModule, d21_slimBP_yellowModule) # merge the data by binding rows 
MF_Mod <- rbind(d7_slimMF_yellowModule, d14_slimMF_blackModule, d21_slimMF_yellowModule) # merge the data by binding rows 

# Load the annotation file 
Geoduck_annotation      <- read.delim2(file="C:/Users/samjg/Documents/My_Projects/Pgenerosa_TagSeq_Metabolomics/TagSeq/Seq_details/Panopea-generosa-genes-annotations.txt", header=F)
annot.condenced         <- Geoduck_annotation[,c(1,5,7)] # build annotation file to merge with master: gene ID, uniprot ID, gene terms with EC (for KEGG)
names(annot.condenced)  <- c('Gene_IDs', 'Uniprot', 'Gene_term_EC') # RENAME THE COLUMNS 
annot.condenced$gene    <- str_extract(annot.condenced$Gene_term_EC, "[^(]+")

# count data
day7.counts.matrix  <- read.csv(file="Analysis/Data/Filtered_Counts/10cpm_50perc/day7.counts.filtered_10cpm50perc.csv", sep=',', header=TRUE)
day14.counts.matrix <- read.csv(file="Analysis/Data/Filtered_Counts/10cpm_50perc/day14.counts.filtered_10cpm50perc.csv", sep=',', header=TRUE)
day21.counts.matrix <- read.csv(file="Analysis/Data/Filtered_Counts/10cpm_50perc/day21.counts.filtered_10cpm50perc.csv", sep=',', header=TRUE)

# trait data
Master.Treatment_Phenotype.data <- read.csv(file="Analysis/Data/Experiment_Metadata/Master_Phyenotype.and.Exp.Treatment_Metadata.csv", sep=',', header=TRUE)

# treatment data
d7.Treatment_data <- Master.Treatment_Phenotype.data %>%   dplyr::filter(Date %in% 20190731) %>%  dplyr::select('Sample.Name', 'Primary_Treatment', 'Second_Treament') # split for day 7 data 
d14.Treatment_data <- Master.Treatment_Phenotype.data  %>%  dplyr::filter(Date %in% 20190807) %>%  dplyr::select('Sample.Name', 'Primary_Treatment', 'Second_Treament')# split for day 7 data 
d21.Treatment_data <- Master.Treatment_Phenotype.data  %>%  dplyr::filter(Date %in% 20190814) %>%  dplyr::select('Sample.Name', 'Primary_Treatment', 'Second_Treament', 'Third_Treatment')# split for day 7 data 

# PREPARE THE DATA FOR THE 3-PANEL TREATMENT PLOTS 
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
# ================================================================================== #
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


















# BIOLOGICAL PROCESS:  HIGHER expression under AMBIENT CONDITIONING ========================================================== #
BP_Amb_GOslim <- BP_Amb %>% 
                  group_by(slim_term) %>% 
                  tally() %>% 
                  dplyr::filter(n == 3) # create the dataset to loop Venn diagrams 

for (i in 1:nrow(BP_Amb_GOslim)) {
  
  # PREP THE LOOP DATA CALLS
  slim <- BP_Amb_GOslim[i,1] # call the term for the loop - creates Venn and Figures 
  # prep data
  slim_BP <- BP_Amb %>% dplyr::filter(slim_term %in% slim)
  slim_BP$slimSplitted  <- strsplit(as.character(slim_BP$Gene_IDs), ";") # split the Gene IDs for Venn 
  slim_d7  <- data.frame(Day = 'Day7', Gene.ID = unlist(slim_BP[1,10]))  # name the first row
  slim_d14 <- data.frame(Day = 'Day14', Gene.ID = unlist(slim_BP[2,10])) # name the second row
  slim_d21 <- data.frame(Day = 'Day21', Gene.ID = unlist(slim_BP[3,10])) # name the third row
  slim_df  <- list(
    Day7_BrownMod           = slim_d7$Gene.ID, 
    Day14_BrownMod          = slim_d14$Gene.ID, 
    Day21_BlueMod           = slim_d21$Gene.ID) # create a list to call in using ggvenn 
  
  # Venn diagram
  VennDiag <- ggvenn(slim_df,  # Venn diagram
                    fill_color = c("white", "white", "white"),
                    stroke_size = 0.5, set_name_size = 4) +
                    ggtitle(paste("WGCNA: ",  slim, sep = ''))
  pdf(paste("Analysis/Output/GO/WGCNA_goseq/PrimaryEffect_Figures/GOslim_venn_tables/Modules_Ambient/all_shared/Ambient_BP_",slim,".pdf", sep =''))
  print(VennDiag)
  # grid.arrange(Day7.volcano.primary, Day14.volcano.primary,ncol=1, nrow=2, clip="off")
  dev.off()
  
  # write a csv for the corresponding Venn
  bind_data <- rbind(slim_d7, slim_d14, slim_d21)
  colnames(bind_data)[2] <- 'Gene_IDs'
  merged_with_annot <- merge(bind_data, annot.condenced[c(1,2,4)], by = 'Gene_IDs') %>% 
                        group_by(Gene_IDs, gene) %>% 
                        mutate(count = n()) %>% 
                        dplyr::arrange(desc(count))
  # write csv
  write.csv(merged_with_annot, paste("Analysis/Output/GO/WGCNA_goseq/PrimaryEffect_Figures/GOslim_venn_tables/Modules_Ambient/all_shared/Ambient_BP_",slim,".csv", sep =''))

  
  panelfig_data <- merged_with_annot %>% dplyr::filter(count >=2)
  
  
  if (nrow(panelfig_data) > 0) {
    Day7.ExpVST_GOIs <- Day7.ExpVST %>% dplyr::select('Sample.Name', 'Primary_Treatment', 'Second_Treament', any_of(panelfig_data$Gene_IDs))
    Day7.ExpVST_GOIs$group <- paste(Day7.ExpVST_GOIs$Primary_Treatment , Day7.ExpVST_GOIs$Second_Treament , sep='')
    
    Day14.ExpVST_GOIs <- Day14.ExpVST %>% dplyr::select('Sample.Name', 'Primary_Treatment', 'Second_Treament', any_of(panelfig_data$Gene_IDs))
    Day14.ExpVST_GOIs$group <- paste(Day14.ExpVST_GOIs$Primary_Treatment , Day14.ExpVST_GOIs$Second_Treament , sep='')
    
    Day21.ExpVST_GOIs <- Day21.ExpVST %>% dplyr::select('Sample.Name', 'Primary_Treatment', 'Second_Treament', 'Third_Treatment', any_of(panelfig_data$Gene_IDs))
    Day21.ExpVST_GOIs$group <- paste(Day21.ExpVST_GOIs$Primary_Treatment , Day21.ExpVST_GOIs$Second_Treament , Day21.ExpVST_GOIs$Third_Treatment,  sep='')
    
    
    # ===================================================================================
    # Day 7 data prep for figures
    #
    # ===================================================================================
    # reshape the data frame and merge geneIDs with tragetGOIs to add back the gene titles - use this to name figures downstream
    Day7.ExpVST_GOIs_MELT <- reshape2::melt(Day7.ExpVST_GOIs, id=(c('Sample.Name', 'Primary_Treatment', 'Second_Treament', 'group'))) # melt using reshape2
    names(Day7.ExpVST_GOIs_MELT)[(5:6)] <- c('Gene_IDs', 'vst_Expression') # change column names
    Day7_ExpVst_Master <- merge(panelfig_data, Day7.ExpVST_GOIs_MELT, by = 'Gene_IDs')
    
    # calc the m an expressio by gene ID (add gene title in group but this is the same unique level as GeneID - review target_GOIs)
    if(nrow(as.data.frame(unique(Day7_ExpVst_Master$Gene_IDs))) == 1) {
      Day7_meanExpr <- Day7_ExpVst_Master %>% 
        dplyr::filter(Day %in% 'Day7') %>% # if panelfig_data is set to >=2, this will plot the mean SE of the total sum of shared genes (as opposed to the exact cetner of the Venn)
        dplyr::select(c('Sample.Name','group', 'vst_Expression', 'Gene_IDs', 'gene')) %>% 
        group_by(group) %>%
        dplyr::summarize(mean.vstExp = mean(vst_Expression), 
                         sd.vsdtExp = sd(vst_Expression),
                         n = n(), 
                         se.vsdtExp = sd.vsdtExp/sqrt(n))
    } else if ( nrow(as.data.frame(unique(Day7_ExpVst_Master$Gene_IDs))) >= 2) {
      Day7_meanExpr <- Day7_ExpVst_Master %>% 
        dplyr::filter(Day %in% 'Day7') %>% # if panelfig_data is set to >=2, this will plot the mean SE of the total sum of shared genes (as opposed to the exact cetner of the Venn)
        dplyr::select(c('Sample.Name','group', 'vst_Expression', 'Gene_IDs', 'gene')) %>% 
        group_by(Sample.Name,group) %>%
        dplyr::summarize(mean.vstbygene = mean(vst_Expression),number_genes = n()) %>% 
        group_by(group) %>% 
        dplyr::summarize(mean.vstExp = mean(mean.vstbygene),
                         sd.vsdtExp = sd(mean.vstbygene),
                         n = n(), 
                         se.vsdtExp = sd.vsdtExp/sqrt(n))
    } else { Day7_meanExpr <- as.data.frame('NULL') }
    
    # create treatment groups 
    Day7_meanExpr$PrimaryTreatment <- substr(Day7_meanExpr$group, 1,1) # primary
    Day7_meanExpr$SecondTreatment <- substr(Day7_meanExpr$group, 2,2) # second
    
    # ===================================================================================
    # Day 14 data prep for figures
    #
    # ===================================================================================
    # reshape the data frame and merge geneIDs with tragetGOIs to add back the gene titles - use this to name figures downstream
    Day14.ExpVST_GOIs_MELT <- reshape2::melt(Day14.ExpVST_GOIs, id=(c('Sample.Name', 'Primary_Treatment', 'Second_Treament', 'group'))) # melt using reshape2
    names(Day14.ExpVST_GOIs_MELT)[(5:6)] <- c('Gene_IDs', 'vst_Expression') # change column names
    Day14_ExpVst_Master <- merge(panelfig_data, Day14.ExpVST_GOIs_MELT, by = 'Gene_IDs')
    
    # calc the m an expressio by gene ID (add gene title in group but this is the same unique level as GeneID - review target_GOIs)
    if(nrow(as.data.frame(unique(Day14_ExpVst_Master$Gene_IDs))) == 1) {
      Day14_meanExpr <- Day14_ExpVst_Master %>% 
        dplyr::filter(Day %in% 'Day14') %>% # if panelfig_data is set to >=2, this will plot the mean SE of the total sum of shared genes (as opposed to the exact cetner of the Venn)
        dplyr::select(c('Sample.Name','group', 'vst_Expression', 'Gene_IDs', 'gene')) %>% 
        group_by(group) %>%
        dplyr::summarize(mean.vstExp = mean(vst_Expression), 
                         sd.vsdtExp = sd(vst_Expression),
                         n = n(), 
                         se.vsdtExp = sd.vsdtExp/sqrt(n))
    } else if ( nrow(as.data.frame(unique(Day14_ExpVst_Master$Gene_IDs))) >= 2) {
      Day14_meanExpr <- Day14_ExpVst_Master %>% 
        dplyr::filter(Day %in% 'Day14') %>% # if panelfig_data is set to >=2, this will plot the mean SE of the total sum of shared genes (as opposed to the exact cetner of the Venn)
        dplyr::select(c('Sample.Name','group', 'vst_Expression', 'Gene_IDs', 'gene')) %>% 
        group_by(Sample.Name,group) %>%
        dplyr::summarize(mean.vstbygene = mean(vst_Expression),number_genes = n()) %>% 
        group_by(group) %>% 
        dplyr::summarize(mean.vstExp = mean(mean.vstbygene),
                         sd.vsdtExp = sd(mean.vstbygene),
                         n = n(), 
                         se.vsdtExp = sd.vsdtExp/sqrt(n))
    } else { Day14_meanExpr <- as.data.frame('NULL') }
    
    # create treatment groups 
    Day14_meanExpr$PrimaryTreatment <- substr(Day14_meanExpr$group, 1,1) # primary
    Day14_meanExpr$SecondTreatment <- substr(Day14_meanExpr$group, 2,2) # second
    
    # ===================================================================================
    # Day 21 data prep for figures
    #
    # ===================================================================================
    # reshape the data frame and merge geneIDs with tragetGOIs to add back the gene titles - use this to name figures downstream
    Day21.ExpVST_GOIs_MELT <- reshape2::melt(Day21.ExpVST_GOIs, id=(c('Sample.Name', 'Primary_Treatment', 'Second_Treament', 'Third_Treatment','group'))) # melt using reshape2
    names(Day21.ExpVST_GOIs_MELT)[(6:7)] <- c('Gene_IDs', 'vst_Expression') # change column names
    Day21_ExpVst_Master <- merge(panelfig_data, Day21.ExpVST_GOIs_MELT, by = 'Gene_IDs')
    
    
    # calc the m an expressio by gene ID (add gene title in group but this is the same unique level as GeneID - review target_GOIs)
    if(nrow(as.data.frame(unique(Day21_ExpVst_Master$Gene_IDs))) == 1) {
      Day21_meanExpr <- Day21_ExpVst_Master %>% 
        dplyr::filter(Day %in% 'Day21') %>% # if panelfig_data is set to >=2, this will plot the mean SE of the total sum of shared genes (as opposed to the exact cetner of the Venn)
        dplyr::select(c('Sample.Name','group', 'vst_Expression', 'Gene_IDs', 'gene')) %>% 
        group_by(group) %>%
        dplyr::summarize(mean.vstExp = mean(vst_Expression), 
                         sd.vsdtExp = sd(vst_Expression),
                         n = n(), 
                         se.vsdtExp = sd.vsdtExp/sqrt(n))
    } else if ( nrow(as.data.frame(unique(Day21_ExpVst_Master$Gene_IDs))) >= 2) {
      Day21_meanExpr <- Day21_ExpVst_Master %>% 
        dplyr::filter(Day %in% 'Day21') %>% # if panelfig_data is set to >=2, this will plot the mean SE of the total sum of shared genes (as opposed to the exact cetner of the Venn)
        dplyr::select(c('Sample.Name','group', 'vst_Expression', 'Gene_IDs', 'gene')) %>% 
        group_by(Sample.Name,group) %>%
        dplyr::summarize(mean.vstbygene = mean(vst_Expression),number_genes = n()) %>% 
        group_by(group) %>% 
        dplyr::summarize(mean.vstExp = mean(mean.vstbygene),
                         sd.vsdtExp = sd(mean.vstbygene),
                         n = n(), 
                         se.vsdtExp = sd.vsdtExp/sqrt(n))
    } else { Day21_meanExpr <- as.data.frame('NULL') }
    
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
    # for(i in 1:nrow(target_GOIs)) {
    #   
    all <- rbind(Day7_meanExpr[,c(1:5)], Day14_meanExpr[,c(1:5)], Day21_meanExpr[,c(1:5)])
    ExpMin <- round( ( min(all$mean.vstExp) - max(all$se.vsdtExp)  ), digits =1)
    ExpMax <- round( ( max(all$mean.vstExp) + max(all$se.vsdtExp)   ), digits =1)
    
    
    if(nrow(Day7_meanExpr) > 1 ) {
      d7_plots <- Day7_meanExpr %>% 
        ggplot(aes(x=SecondTreatment, y=mean.vstExp, fill=PrimaryTreatment)) +  # , colour=supp, group=supp))
        theme_classic() +
        geom_errorbar(aes(ymin=mean.vstExp-se.vsdtExp, ymax=mean.vstExp+se.vsdtExp), colour="black", width=.1, position=pd) +
        geom_point(position=pd, size = 4, shape=21) +            
        xlab("Second pCO2 treatment") +
        ylab('') +                 # note the mean was first by sample ID THEN by treatment
        scale_fill_manual(values=c("white","grey50")) +
        # scale_color_manual(values=c("white","grey50")) +
        ggtitle(paste("Day 7:",substr(slim,1,13), " (N = ", nrow(unique(as.data.frame((Day7_ExpVst_Master %>% dplyr::filter(Day %in% 'Day7'))$Gene_IDs))), ")", sep='')) +
        # expand_limits(y=0) +                                                    # Expand y range
        #scale_y_continuous(limits=c((min_p1), (max_p1))) +
        theme(text = element_text(size=15)) +
        scale_y_continuous(limits = c(ExpMin,ExpMax)) +
        theme(legend.position = "none")
    } else { d7_plots <- plot.new() }
    # 
    #   
    if(nrow(Day14_meanExpr) > 1 ) {
      d14_plots <- Day14_meanExpr %>% 
        # dplyr::filter(genes %in% target_GOIs[i,1]) %>% 
        ggplot(aes(x=SecondTreatment, y=mean.vstExp, fill=PrimaryTreatment)) +  # , colour=supp, group=supp))
        theme_classic() +
        geom_errorbar(aes(ymin=mean.vstExp-se.vsdtExp, ymax=mean.vstExp+se.vsdtExp), colour="black", width=.1, position=pd) +
        geom_point(position=pd, size = 4, shape=21) +            
        xlab("Second pCO2 treatment") +
        ylab("Gene Expression (mean±SE VST transformed)") +                 # note the mean was first by sample ID THEN by treatment
        scale_fill_manual(values=c("white","grey50")) +
        # scale_color_manual(values=c("white","grey50")) +
        ggtitle(paste("Day 14:",substr(slim,1,13), " (N = ", nrow(unique(as.data.frame((Day14_ExpVst_Master %>% dplyr::filter(Day %in% 'Day14'))$Gene_IDs))), ")", sep='')) +
        # expand_limits(y=0) +                                                    # Expand y range
        #scale_y_continuous(limits=c((min_p1), (max_p1))) +
        theme(text = element_text(size=15)) +
        scale_y_continuous(limits = c(ExpMin,ExpMax)) +
        theme(legend.position = "none")
    } else { d14_plots <- plot.new() }
    # 
    # 
    #       
    if(nrow(Day21_meanExpr) > 1 ) {
      d21_plots <- Day21_meanExpr %>% 
        # dplyr::filter(genes %in% target_GOIs[i,1]) %>% 
        ggplot(aes(x=ThirdTreatment, y=mean.vstExp, fill=PrimaryTreatment)) +  # , colour=supp, group=supp))
        theme_classic() +
        geom_errorbar(aes(ymin=mean.vstExp-se.vsdtExp, ymax=mean.vstExp+se.vsdtExp), colour="black", width=.1, position=pd) +
        geom_point(position=pd, size = 4, shape=21) +            
        xlab("Third pCO2 treatment") +
        ylab('') +                 # note the mean was first by sample ID THEN by treatment
        scale_fill_manual(values=c("white","grey50")) +
        # scale_color_manual(values=c("white","grey50")) +
        ggtitle(NULL) +
        # expand_limits(y=0) +                                                    # Expand y range
        #scale_y_continuous(limits=c((min_p1), (max_p1))) +
        theme(text = element_text(size=15)) +
        theme(legend.position = "none") +
        scale_y_continuous(limits = c(ExpMin,ExpMax)) +
        facet_wrap(~SecondTreatment)
    } else { d21_plots <- plot.new() }
    
    pdf(paste("Analysis/Output/GO/WGCNA_goseq/PrimaryEffect_Figures/GOslim_venn_tables/Modules_Ambient/all_shared/Ambient_BP_",substr(slim,1,13),"_vstExpression.pdf", sep =''), width=15, height=6)
    print(ggarrange(d7_plots, d14_plots, d21_plots,        
                    plotlist = NULL,
                    ncol = 3,
                    nrow = 1,
                    labels = NULL))
    dev.off()
  } # plotting for loop
  
  else {}
  
}
































# MOLECULAR FUNCTION:  HIGHER expression under AMBIENT CONDITIONING ========================================================== #
MF_Amb_GOslim <- MF_Amb %>% 
  filter(!module_day %in% "Day21_magenta", Gene_count >0) %>% 
  # filter(Gene_count >0) %>% 
  group_by(slim_term) %>% 
  tally() %>% 
  dplyr::filter(n >= 3) # create the dataset to loop Venn diagrams 

for (i in 1:nrow(MF_Amb_GOslim)) {
  slim <- MF_Amb_GOslim[i,1] # call the term for the loop - creates Venn and Figures 
  
  # prep data
  slim_MF <- MF_Amb  %>% dplyr::filter(slim_term %in% slim, Gene_count >0)
  slim_MF$slimSplitted  <- strsplit(as.character(slim_MF$Gene_IDs), ";") # split the Gene IDs for Venn 
  slim_d7  <- data.frame(Day = 'Day7',  Gene.ID = unlist(slim_MF[1,10]))  # name the first row
  slim_d14 <- data.frame(Day = 'Day14', Gene.ID = unlist(slim_MF[2,10])) # name the second row
  slim_d21 <- data.frame(Day = 'Day21', Gene.ID = unlist(slim_MF[3,10])) # name the third row
  slim_df  <- list(
    Day7_BrownMod           = slim_d7$Gene.ID, 
    Day14_BrownMod          = slim_d14$Gene.ID, 
    Day21_BlueMagentaMod    = slim_d21$Gene.ID) # create a list to call in using ggvenn 
  
  # Venn diagram
  VennDiag <- ggvenn(slim_df,  # Venn diagram
                     fill_color = c("white", "white", "white"),
                     stroke_size = 0.5, set_name_size = 4) +
    ggtitle(paste("WGCNA: ",  slim, sep = ''))
  pdf(paste("Analysis/Output/GO/WGCNA_goseq/PrimaryEffect_Figures/GOslim_venn_tables/Modules_Ambient/all_shared/Ambient_MF_",substr(slim,1,13),".pdf", sep =''))
  print(VennDiag)
  # grid.arrange(Day7.volcano.primary, Day14.volcano.primary,ncol=1, nrow=2, clip="off")
  dev.off()
  
  # write a csv for the corresponding Venn
  bind_data <- rbind(slim_d7, slim_d14, slim_d21)
  colnames(bind_data)[2] <- 'Gene_IDs'
  merged_with_annot <- merge(bind_data, annot.condenced[c(1,2,4)], by = 'Gene_IDs') %>% 
    group_by(Gene_IDs, gene) %>% 
    mutate(count = n()) %>% 
    dplyr::arrange(desc(count))
  # write csv
  write.csv(merged_with_annot, paste("Analysis/Output/GO/WGCNA_goseq/PrimaryEffect_Figures/GOslim_venn_tables/Modules_Ambient/all_shared/Ambient_MF_",substr(slim,1,13),".csv", sep =''))
  
  panelfig_data <- merged_with_annot %>% dplyr::filter(count >=2)
  
  
  if (nrow(panelfig_data) > 0) {
    Day7.ExpVST_GOIs <- Day7.ExpVST %>% dplyr::select('Sample.Name', 'Primary_Treatment', 'Second_Treament', any_of(panelfig_data$Gene_IDs))
    Day7.ExpVST_GOIs$group <- paste(Day7.ExpVST_GOIs$Primary_Treatment , Day7.ExpVST_GOIs$Second_Treament , sep='')
    
    Day14.ExpVST_GOIs <- Day14.ExpVST %>% dplyr::select('Sample.Name', 'Primary_Treatment', 'Second_Treament', any_of(panelfig_data$Gene_IDs))
    Day14.ExpVST_GOIs$group <- paste(Day14.ExpVST_GOIs$Primary_Treatment , Day14.ExpVST_GOIs$Second_Treament , sep='')
    
    Day21.ExpVST_GOIs <- Day21.ExpVST %>% dplyr::select('Sample.Name', 'Primary_Treatment', 'Second_Treament', 'Third_Treatment', any_of(panelfig_data$Gene_IDs))
    Day21.ExpVST_GOIs$group <- paste(Day21.ExpVST_GOIs$Primary_Treatment , Day21.ExpVST_GOIs$Second_Treament , Day21.ExpVST_GOIs$Third_Treatment,  sep='')
    
    
    # ===================================================================================
    # Day 7 data prep for figures
    #
    # ===================================================================================
    # reshape the data frame and merge geneIDs with tragetGOIs to add back the gene titles - use this to name figures downstream
    Day7.ExpVST_GOIs_MELT <- reshape2::melt(Day7.ExpVST_GOIs, id=(c('Sample.Name', 'Primary_Treatment', 'Second_Treament', 'group'))) # melt using reshape2
    names(Day7.ExpVST_GOIs_MELT)[(5:6)] <- c('Gene_IDs', 'vst_Expression') # change column names
    Day7_ExpVst_Master <- merge(panelfig_data, Day7.ExpVST_GOIs_MELT, by = 'Gene_IDs')
    
    # calc the m an expressio by gene ID (add gene title in group but this is the same unique level as GeneID - review target_GOIs)
    if(nrow(as.data.frame(unique(Day7_ExpVst_Master$Gene_IDs))) == 1) {
      Day7_meanExpr <- Day7_ExpVst_Master %>% 
        dplyr::filter(Day %in% 'Day7') %>% # if panelfig_data is set to >=2, this will plot the mean SE of the total sum of shared genes (as opposed to the exact cetner of the Venn)
        dplyr::select(c('Sample.Name','group', 'vst_Expression', 'Gene_IDs', 'gene')) %>% 
        group_by(group) %>%
        dplyr::summarize(mean.vstExp = mean(vst_Expression), 
                         sd.vsdtExp = sd(vst_Expression),
                         n = n(), 
                         se.vsdtExp = sd.vsdtExp/sqrt(n))
    } else if ( nrow(as.data.frame(unique(Day7_ExpVst_Master$Gene_IDs))) >= 2) {
      Day7_meanExpr <- Day7_ExpVst_Master %>% 
        dplyr::filter(Day %in% 'Day7') %>% # if panelfig_data is set to >=2, this will plot the mean SE of the total sum of shared genes (as opposed to the exact cetner of the Venn)
        dplyr::select(c('Sample.Name','group', 'vst_Expression', 'Gene_IDs', 'gene')) %>% 
        group_by(Sample.Name,group) %>%
        dplyr::summarize(mean.vstbygene = mean(vst_Expression),number_genes = n()) %>% 
        group_by(group) %>% 
        dplyr::summarize(mean.vstExp = mean(mean.vstbygene),
                         sd.vsdtExp = sd(mean.vstbygene),
                         n = n(), 
                         se.vsdtExp = sd.vsdtExp/sqrt(n))
    } else { Day7_meanExpr <- as.data.frame('NULL') }
    
    # create treatment groups 
    Day7_meanExpr$PrimaryTreatment <- substr(Day7_meanExpr$group, 1,1) # primary
    Day7_meanExpr$SecondTreatment <- substr(Day7_meanExpr$group, 2,2) # second
    
    # ===================================================================================
    # Day 14 data prep for figures
    #
    # ===================================================================================
    # reshape the data frame and merge geneIDs with tragetGOIs to add back the gene titles - use this to name figures downstream
    Day14.ExpVST_GOIs_MELT <- reshape2::melt(Day14.ExpVST_GOIs, id=(c('Sample.Name', 'Primary_Treatment', 'Second_Treament', 'group'))) # melt using reshape2
    names(Day14.ExpVST_GOIs_MELT)[(5:6)] <- c('Gene_IDs', 'vst_Expression') # change column names
    Day14_ExpVst_Master <- merge(panelfig_data, Day14.ExpVST_GOIs_MELT, by = 'Gene_IDs')
    
    # calc the m an expressio by gene ID (add gene title in group but this is the same unique level as GeneID - review target_GOIs)
    if(nrow(as.data.frame(unique(Day14_ExpVst_Master$Gene_IDs))) == 1) {
      Day14_meanExpr <- Day14_ExpVst_Master %>% 
        dplyr::filter(Day %in% 'Day14') %>% # if panelfig_data is set to >=2, this will plot the mean SE of the total sum of shared genes (as opposed to the exact cetner of the Venn)
        dplyr::select(c('Sample.Name','group', 'vst_Expression', 'Gene_IDs', 'gene')) %>% 
        group_by(group) %>%
        dplyr::summarize(mean.vstExp = mean(vst_Expression), 
                         sd.vsdtExp = sd(vst_Expression),
                         n = n(), 
                         se.vsdtExp = sd.vsdtExp/sqrt(n))
    } else if ( nrow(as.data.frame(unique(Day14_ExpVst_Master$Gene_IDs))) >= 2) {
      Day14_meanExpr <- Day14_ExpVst_Master %>% 
        dplyr::filter(Day %in% 'Day14') %>% # if panelfig_data is set to >=2, this will plot the mean SE of the total sum of shared genes (as opposed to the exact cetner of the Venn)
        dplyr::select(c('Sample.Name','group', 'vst_Expression', 'Gene_IDs', 'gene')) %>% 
        group_by(Sample.Name,group) %>%
        dplyr::summarize(mean.vstbygene = mean(vst_Expression),number_genes = n()) %>% 
        group_by(group) %>% 
        dplyr::summarize(mean.vstExp = mean(mean.vstbygene),
                         sd.vsdtExp = sd(mean.vstbygene),
                         n = n(), 
                         se.vsdtExp = sd.vsdtExp/sqrt(n))
    } else { Day14_meanExpr <- as.data.frame('NULL') }
    
    # create treatment groups 
    Day14_meanExpr$PrimaryTreatment <- substr(Day14_meanExpr$group, 1,1) # primary
    Day14_meanExpr$SecondTreatment <- substr(Day14_meanExpr$group, 2,2) # second
    
    # ===================================================================================
    # Day 21 data prep for figures
    #
    # ===================================================================================
    # reshape the data frame and merge geneIDs with tragetGOIs to add back the gene titles - use this to name figures downstream
    Day21.ExpVST_GOIs_MELT <- reshape2::melt(Day21.ExpVST_GOIs, id=(c('Sample.Name', 'Primary_Treatment', 'Second_Treament', 'Third_Treatment','group'))) # melt using reshape2
    names(Day21.ExpVST_GOIs_MELT)[(6:7)] <- c('Gene_IDs', 'vst_Expression') # change column names
    Day21_ExpVst_Master <- merge(panelfig_data, Day21.ExpVST_GOIs_MELT, by = 'Gene_IDs')
   
    
    # calc the m an expressio by gene ID (add gene title in group but this is the same unique level as GeneID - review target_GOIs)
    if(nrow(as.data.frame(unique(Day21_ExpVst_Master$Gene_IDs))) == 1) {
      Day21_meanExpr <- Day21_ExpVst_Master %>% 
        dplyr::filter(Day %in% 'Day21') %>% # if panelfig_data is set to >=2, this will plot the mean SE of the total sum of shared genes (as opposed to the exact cetner of the Venn)
        dplyr::select(c('Sample.Name','group', 'vst_Expression', 'Gene_IDs', 'gene')) %>% 
        group_by(group) %>%
        dplyr::summarize(mean.vstExp = mean(vst_Expression), 
                         sd.vsdtExp = sd(vst_Expression),
                         n = n(), 
                         se.vsdtExp = sd.vsdtExp/sqrt(n))
    } else if ( nrow(as.data.frame(unique(Day21_ExpVst_Master$Gene_IDs))) >= 2) {
      Day21_meanExpr <- Day21_ExpVst_Master %>% 
        dplyr::filter(Day %in% 'Day21') %>% # if panelfig_data is set to >=2, this will plot the mean SE of the total sum of shared genes (as opposed to the exact cetner of the Venn)
        dplyr::select(c('Sample.Name','group', 'vst_Expression', 'Gene_IDs', 'gene')) %>% 
        group_by(Sample.Name,group) %>%
        dplyr::summarize(mean.vstbygene = mean(vst_Expression),number_genes = n()) %>% 
        group_by(group) %>% 
        dplyr::summarize(mean.vstExp = mean(mean.vstbygene),
                         sd.vsdtExp = sd(mean.vstbygene),
                         n = n(), 
                         se.vsdtExp = sd.vsdtExp/sqrt(n))
    } else { Day21_meanExpr <- as.data.frame(NULL) }

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
    # for(i in 1:nrow(target_GOIs)) {
    #   
    all <- rbind(Day7_meanExpr[,c(1:5)], Day14_meanExpr[,c(1:5)], Day21_meanExpr[,c(1:5)])
    ExpMin <- round( ( min(all$mean.vstExp) - max(all$se.vsdtExp)  ), digits =1)
    ExpMax <- round( ( max(all$mean.vstExp) + max(all$se.vsdtExp)   ), digits =1)
    
    
    if(nrow(Day7_meanExpr) > 1 ) {
    d7_plots <- Day7_meanExpr %>% 
      ggplot(aes(x=SecondTreatment, y=mean.vstExp, fill=PrimaryTreatment)) +  # , colour=supp, group=supp))
      theme_classic() +
      geom_errorbar(aes(ymin=mean.vstExp-se.vsdtExp, ymax=mean.vstExp+se.vsdtExp), colour="black", width=.1, position=pd) +
      geom_point(position=pd, size = 4, shape=21) +            
      xlab("Second pCO2 treatment") +
      ylab('') +                 # note the mean was first by sample ID THEN by treatment
      scale_fill_manual(values=c("white","grey50")) +
      # scale_color_manual(values=c("white","grey50")) +
      ggtitle(paste("Day 7:",substr(slim,1,13), " (N = ", nrow(unique(as.data.frame((Day7_ExpVst_Master %>% dplyr::filter(Day %in% 'Day7'))$Gene_IDs))), ")", sep='')) +
      # expand_limits(y=0) +                                                    # Expand y range
      #scale_y_continuous(limits=c((min_p1), (max_p1))) +
      theme(text = element_text(size=15)) +
      scale_y_continuous(limits = c(ExpMin,ExpMax)) +
      theme(legend.position = "none")
    } else { d7_plots <- plot.new() }
    # 
    #   
    if(nrow(Day14_meanExpr) > 1 ) {
    d14_plots <- Day14_meanExpr %>% 
      # dplyr::filter(genes %in% target_GOIs[i,1]) %>% 
      ggplot(aes(x=SecondTreatment, y=mean.vstExp, fill=PrimaryTreatment)) +  # , colour=supp, group=supp))
      theme_classic() +
      geom_errorbar(aes(ymin=mean.vstExp-se.vsdtExp, ymax=mean.vstExp+se.vsdtExp), colour="black", width=.1, position=pd) +
      geom_point(position=pd, size = 4, shape=21) +            
      xlab("Second pCO2 treatment") +
      ylab("Gene Expression (mean±SE VST transformed)") +                 # note the mean was first by sample ID THEN by treatment
      scale_fill_manual(values=c("white","grey50")) +
      # scale_color_manual(values=c("white","grey50")) +
      ggtitle(paste("Day 14:",substr(slim,1,13), " (N = ", nrow(unique(as.data.frame((Day14_ExpVst_Master %>% dplyr::filter(Day %in% 'Day14'))$Gene_IDs))), ")", sep='')) +
      # expand_limits(y=0) +                                                    # Expand y range
      #scale_y_continuous(limits=c((min_p1), (max_p1))) +
      theme(text = element_text(size=15)) +
      scale_y_continuous(limits = c(ExpMin,ExpMax)) +
      theme(legend.position = "none")
    } else { d14_plots <- plot.new() }
    # 
    # 
    #       
    if(nrow(Day21_meanExpr) > 1 ) {
    d21_plots <- Day21_meanExpr %>% 
      # dplyr::filter(genes %in% target_GOIs[i,1]) %>% 
      ggplot(aes(x=ThirdTreatment, y=mean.vstExp, fill=PrimaryTreatment)) +  # , colour=supp, group=supp))
      theme_classic() +
      geom_errorbar(aes(ymin=mean.vstExp-se.vsdtExp, ymax=mean.vstExp+se.vsdtExp), colour="black", width=.1, position=pd) +
      geom_point(position=pd, size = 4, shape=21) +            
      xlab("Third pCO2 treatment") +
      ylab('') +                 # note the mean was first by sample ID THEN by treatment
      scale_fill_manual(values=c("white","grey50")) +
      # scale_color_manual(values=c("white","grey50")) +
      ggtitle(NULL) +
      # expand_limits(y=0) +                                                    # Expand y range
      #scale_y_continuous(limits=c((min_p1), (max_p1))) +
      theme(text = element_text(size=15)) +
      theme(legend.position = "none") +
      scale_y_continuous(limits = c(ExpMin,ExpMax)) +
      facet_wrap(~SecondTreatment)
    } else { d21_plots <- plot.new() }
    
    pdf(paste("Analysis/Output/GO/WGCNA_goseq/PrimaryEffect_Figures/GOslim_venn_tables/Modules_Ambient/all_shared/Ambient_MF_",substr(slim,1,13),"_vstExpression.pdf", sep =''), width=15, height=6)
    print(ggarrange(d7_plots, d14_plots, d21_plots,        
                    plotlist = NULL,
                    ncol = 3,
                    nrow = 1,
                    labels = NULL))
    dev.off()
  } # plotting for loop
  
  else {}
  
}

































# Biological Process:  HIGHER expression under Moderate CONDITIONING ========================================================== #
BP_Mod_GOslim <- BP_Mod %>% 
  group_by(slim_term) %>% 
  tally() %>% 
  dplyr::filter(n == 3) # create the dataset to loop Venn diagrams 

for (i in 1:nrow(BP_Mod_GOslim)) {
  slim <- BP_Mod_GOslim[i,1] # call the term for the loop - creates Venn and Figures 
  
  # prep data
  slim_BP <- BP_Mod %>% dplyr::filter(slim_term %in% slim)
  slim_BP$slimSplitted  <- strsplit(as.character(slim_BP$Gene_IDs), ";") # split the Gene IDs for Venn 
  slim_d7  <- data.frame(Day = 'Day7', Gene.ID = unlist(slim_BP[1,10]))  # name the first row
  slim_d14 <- data.frame(Day = 'Day14', Gene.ID = unlist(slim_BP[2,10])) # name the second row
  slim_d21 <- data.frame(Day = 'Day21', Gene.ID = unlist(slim_BP[3,10])) # name the third row
  slim_df  <- list(
    Day7_YellowMod          = slim_d7$Gene.ID, 
    Day14_BlackMod          = slim_d14$Gene.ID, 
    Day21_YellowMod        = slim_d21$Gene.ID) # create a list to call in using ggvenn 
  
  # Venn diagram
  VennDiag <- ggvenn(slim_df,  # Venn diagram
                     fill_color = c("white", "white", "white"),
                     stroke_size = 0.5, set_name_size = 4) +
    ggtitle(paste("WGCNA: ",  slim, sep = ''))
  pdf(paste("Analysis/Output/GO/WGCNA_goseq/PrimaryEffect_Figures/GOslim_venn_tables/Modules_Moderate/Moderate_BP_",substr(slim,1,13),".pdf", sep =''))
  print(VennDiag)
  # grid.arrange(Day7.volcano.primary, Day14.volcano.primary,ncol=1, nrow=2, clip="off")
  dev.off()
  
  # write a csv for the corresponding Venn
  bind_data <- rbind(slim_d7, slim_d14, slim_d21)
  colnames(bind_data)[2] <- 'Gene_IDs'
  merged_with_annot <- merge(bind_data, annot.condenced[c(1,2,4)], by = 'Gene_IDs') %>% 
    group_by(Gene_IDs, gene) %>% 
    mutate(count = n()) %>% 
    dplyr::arrange(desc(count))
  merged_with_annot <- as.data.frame(merged_with_annot)
  # write csv
  write.csv(merged_with_annot, paste("Analysis/Output/GO/WGCNA_goseq/PrimaryEffect_Figures/GOslim_venn_tables/Modules_Moderate/Moderate_BP_",substr(slim,1,13),".csv", sep =''))
  
  panelfig_data <- merged_with_annot %>% dplyr::filter(count >=2)
  
  
  if (nrow(panelfig_data) > 0) {
    Day7.ExpVST_GOIs <- Day7.ExpVST %>% dplyr::select('Sample.Name', 'Primary_Treatment', 'Second_Treament', any_of(panelfig_data$Gene_IDs))
    Day7.ExpVST_GOIs$group <- paste(Day7.ExpVST_GOIs$Primary_Treatment , Day7.ExpVST_GOIs$Second_Treament , sep='')
    
    Day14.ExpVST_GOIs <- Day14.ExpVST %>% dplyr::select('Sample.Name', 'Primary_Treatment', 'Second_Treament', any_of(panelfig_data$Gene_IDs))
    Day14.ExpVST_GOIs$group <- paste(Day14.ExpVST_GOIs$Primary_Treatment , Day14.ExpVST_GOIs$Second_Treament , sep='')
    
    Day21.ExpVST_GOIs <- Day21.ExpVST %>% dplyr::select('Sample.Name', 'Primary_Treatment', 'Second_Treament', 'Third_Treatment', any_of(panelfig_data$Gene_IDs))
    Day21.ExpVST_GOIs$group <- paste(Day21.ExpVST_GOIs$Primary_Treatment , Day21.ExpVST_GOIs$Second_Treament , Day21.ExpVST_GOIs$Third_Treatment,  sep='')
    
    
    # ===================================================================================
    # Day 7 data prep for figures
    #
    # ===================================================================================
    # reshape the data frame and merge geneIDs with tragetGOIs to add back the gene titles - use this to name figures downstream
    Day7.ExpVST_GOIs_MELT <- reshape2::melt(Day7.ExpVST_GOIs, id=(c('Sample.Name', 'Primary_Treatment', 'Second_Treament', 'group'))) # melt using reshape2
    names(Day7.ExpVST_GOIs_MELT)[(5:6)] <- c('Gene_IDs', 'vst_Expression') # change column names
    Day7_ExpVst_Master <- merge(panelfig_data, Day7.ExpVST_GOIs_MELT, by = 'Gene_IDs')
    
    # calc the m an expressio by gene ID (add gene title in group but this is the same unique level as GeneID - review target_GOIs)
    if(nrow(as.data.frame(unique(Day7_ExpVst_Master$Gene_IDs))) == 1) {
      Day7_meanExpr <- Day7_ExpVst_Master %>% 
        dplyr::filter(Day %in% 'Day7') %>% # if panelfig_data is set to >=2, this will plot the mean SE of the total sum of shared genes (as opposed to the exact cetner of the Venn)
        dplyr::select(c('Sample.Name','group', 'vst_Expression', 'Gene_IDs', 'gene')) %>% 
        group_by(group) %>%
        dplyr::summarize(mean.vstExp = mean(vst_Expression), 
                         sd.vsdtExp = sd(vst_Expression),
                         n = n(), 
                         se.vsdtExp = sd.vsdtExp/sqrt(n))
    } else if ( nrow(as.data.frame(unique(Day7_ExpVst_Master$Gene_IDs))) >= 2) {
      Day7_meanExpr <- Day7_ExpVst_Master %>% 
        dplyr::filter(Day %in% 'Day7') %>% # if panelfig_data is set to >=2, this will plot the mean SE of the total sum of shared genes (as opposed to the exact cetner of the Venn)
        dplyr::select(c('Sample.Name','group', 'vst_Expression', 'Gene_IDs', 'gene')) %>% 
        group_by(Sample.Name,group) %>%
        dplyr::summarize(mean.vstbygene = mean(vst_Expression),number_genes = n()) %>% 
        group_by(group) %>% 
        dplyr::summarize(mean.vstExp = mean(mean.vstbygene),
                         sd.vsdtExp = sd(mean.vstbygene),
                         n = n(), 
                         se.vsdtExp = sd.vsdtExp/sqrt(n))
    } else { Day7_meanExpr <- as.data.frame('NULL') }
    
    # create treatment groups 
    Day7_meanExpr$PrimaryTreatment <- substr(Day7_meanExpr$group, 1,1) # primary
    Day7_meanExpr$SecondTreatment <- substr(Day7_meanExpr$group, 2,2) # second
    
    # ===================================================================================
    # Day 14 data prep for figures
    #
    # ===================================================================================
    # reshape the data frame and merge geneIDs with tragetGOIs to add back the gene titles - use this to name figures downstream
    Day14.ExpVST_GOIs_MELT <- reshape2::melt(Day14.ExpVST_GOIs, id=(c('Sample.Name', 'Primary_Treatment', 'Second_Treament', 'group'))) # melt using reshape2
    names(Day14.ExpVST_GOIs_MELT)[(5:6)] <- c('Gene_IDs', 'vst_Expression') # change column names
    Day14_ExpVst_Master <- merge(panelfig_data, Day14.ExpVST_GOIs_MELT, by = 'Gene_IDs')
    
    # calc the m an expressio by gene ID (add gene title in group but this is the same unique level as GeneID - review target_GOIs)
    if(nrow(as.data.frame(unique(Day14_ExpVst_Master$Gene_IDs))) == 1) {
      Day14_meanExpr <- Day14_ExpVst_Master %>% 
        dplyr::filter(Day %in% 'Day14') %>% # if panelfig_data is set to >=2, this will plot the mean SE of the total sum of shared genes (as opposed to the exact cetner of the Venn)
        dplyr::select(c('Sample.Name','group', 'vst_Expression', 'Gene_IDs', 'gene')) %>% 
        group_by(group) %>%
        dplyr::summarize(mean.vstExp = mean(vst_Expression), 
                         sd.vsdtExp = sd(vst_Expression),
                         n = n(), 
                         se.vsdtExp = sd.vsdtExp/sqrt(n))
    } else if ( nrow(as.data.frame(unique(Day14_ExpVst_Master$Gene_IDs))) >= 2) {
      Day14_meanExpr <- Day14_ExpVst_Master %>% 
        dplyr::filter(Day %in% 'Day14') %>% # if panelfig_data is set to >=2, this will plot the mean SE of the total sum of shared genes (as opposed to the exact cetner of the Venn)
        dplyr::select(c('Sample.Name','group', 'vst_Expression', 'Gene_IDs', 'gene')) %>% 
        group_by(Sample.Name,group) %>%
        dplyr::summarize(mean.vstbygene = mean(vst_Expression),number_genes = n()) %>% 
        group_by(group) %>% 
        dplyr::summarize(mean.vstExp = mean(mean.vstbygene),
                         sd.vsdtExp = sd(mean.vstbygene),
                         n = n(), 
                         se.vsdtExp = sd.vsdtExp/sqrt(n))
    } else { Day14_meanExpr <- as.data.frame('NULL') }
    
    # create treatment groups 
    Day14_meanExpr$PrimaryTreatment <- substr(Day14_meanExpr$group, 1,1) # primary
    Day14_meanExpr$SecondTreatment <- substr(Day14_meanExpr$group, 2,2) # second
    
    # ===================================================================================
    # Day 21 data prep for figures
    #
    # ===================================================================================
    # reshape the data frame and merge geneIDs with tragetGOIs to add back the gene titles - use this to name figures downstream
    Day21.ExpVST_GOIs_MELT <- reshape2::melt(Day21.ExpVST_GOIs, id=(c('Sample.Name', 'Primary_Treatment', 'Second_Treament', 'Third_Treatment','group'))) # melt using reshape2
    names(Day21.ExpVST_GOIs_MELT)[(6:7)] <- c('Gene_IDs', 'vst_Expression') # change column names
    Day21_ExpVst_Master <- merge(panelfig_data, Day21.ExpVST_GOIs_MELT, by = 'Gene_IDs')
    
    
    # calc the m an expressio by gene ID (add gene title in group but this is the same unique level as GeneID - review target_GOIs)
    if(nrow(as.data.frame(unique(Day21_ExpVst_Master$Gene_IDs))) == 1) {
      Day21_meanExpr <- Day21_ExpVst_Master %>% 
        dplyr::filter(Day %in% 'Day21') %>% # if panelfig_data is set to >=2, this will plot the mean SE of the total sum of shared genes (as opposed to the exact cetner of the Venn)
        dplyr::select(c('Sample.Name','group', 'vst_Expression', 'Gene_IDs', 'gene')) %>% 
        group_by(group) %>%
        dplyr::summarize(mean.vstExp = mean(vst_Expression), 
                         sd.vsdtExp = sd(vst_Expression),
                         n = n(), 
                         se.vsdtExp = sd.vsdtExp/sqrt(n))
    } else if ( nrow(as.data.frame(unique(Day21_ExpVst_Master$Gene_IDs))) >= 2) {
      Day21_meanExpr <- Day21_ExpVst_Master %>% 
        dplyr::filter(Day %in% 'Day21') %>% # if panelfig_data is set to >=2, this will plot the mean SE of the total sum of shared genes (as opposed to the exact cetner of the Venn)
        dplyr::select(c('Sample.Name','group', 'vst_Expression', 'Gene_IDs', 'gene')) %>% 
        group_by(Sample.Name,group) %>%
        dplyr::summarize(mean.vstbygene = mean(vst_Expression),number_genes = n()) %>% 
        group_by(group) %>% 
        dplyr::summarize(mean.vstExp = mean(mean.vstbygene),
                         sd.vsdtExp = sd(mean.vstbygene),
                         n = n(), 
                         se.vsdtExp = sd.vsdtExp/sqrt(n))
    } else { Day21_meanExpr <- as.data.frame(NULL) }
    
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
    # for(i in 1:nrow(target_GOIs)) {
    #   
    all <- rbind(Day7_meanExpr[,c(1:5)], Day14_meanExpr[,c(1:5)], Day21_meanExpr[,c(1:5)])
    ExpMin <- round( ( min(all$mean.vstExp) - max(all$se.vsdtExp)  ), digits =1)
    ExpMax <- round( ( max(all$mean.vstExp) + max(all$se.vsdtExp)   ), digits =1)
    
    
    if(nrow(Day7_meanExpr) > 1 ) {
      d7_plots <- Day7_meanExpr %>% 
        ggplot(aes(x=SecondTreatment, y=mean.vstExp, fill=PrimaryTreatment)) +  # , colour=supp, group=supp))
        theme_classic() +
        geom_errorbar(aes(ymin=mean.vstExp-se.vsdtExp, ymax=mean.vstExp+se.vsdtExp), colour="black", width=.1, position=pd) +
        geom_point(position=pd, size = 4, shape=21) +            
        xlab("Second pCO2 treatment") +
        ylab('') +                 # note the mean was first by sample ID THEN by treatment
        scale_fill_manual(values=c("white","grey50")) +
        # scale_color_manual(values=c("white","grey50")) +
        ggtitle(paste("Day 7:",substr(slim,1,13), " (N = ", nrow(unique(as.data.frame((Day7_ExpVst_Master %>% dplyr::filter(Day %in% 'Day7'))$Gene_IDs))), ")", sep='')) +
        # expand_limits(y=0) +                                                    # Expand y range
        #scale_y_continuous(limits=c((min_p1), (max_p1))) +
        theme(text = element_text(size=15)) +
        scale_y_continuous(limits = c(ExpMin,ExpMax)) +
        theme(legend.position = "none")
    } else { d7_plots <- plot.new() }
    # 
    #   
    if(nrow(Day14_meanExpr) > 1 ) {
      d14_plots <- Day14_meanExpr %>% 
        # dplyr::filter(genes %in% target_GOIs[i,1]) %>% 
        ggplot(aes(x=SecondTreatment, y=mean.vstExp, fill=PrimaryTreatment)) +  # , colour=supp, group=supp))
        theme_classic() +
        geom_errorbar(aes(ymin=mean.vstExp-se.vsdtExp, ymax=mean.vstExp+se.vsdtExp), colour="black", width=.1, position=pd) +
        geom_point(position=pd, size = 4, shape=21) +            
        xlab("Second pCO2 treatment") +
        ylab("Gene Expression (mean±SE VST transformed)") +                 # note the mean was first by sample ID THEN by treatment
        scale_fill_manual(values=c("white","grey50")) +
        # scale_color_manual(values=c("white","grey50")) +
        ggtitle(paste("Day 14:",substr(slim,1,13), " (N = ", nrow(unique(as.data.frame((Day14_ExpVst_Master %>% dplyr::filter(Day %in% 'Day14'))$Gene_IDs))), ")", sep='')) +
        # expand_limits(y=0) +                                                    # Expand y range
        #scale_y_continuous(limits=c((min_p1), (max_p1))) +
        theme(text = element_text(size=15)) +
        scale_y_continuous(limits = c(ExpMin,ExpMax)) +
        theme(legend.position = "none")
    } else { d14_plots <- plot.new() }
    # 
    # 
    #       
    if(nrow(Day21_meanExpr) > 1 ) {
      d21_plots <- Day21_meanExpr %>% 
        # dplyr::filter(genes %in% target_GOIs[i,1]) %>% 
        ggplot(aes(x=ThirdTreatment, y=mean.vstExp, fill=PrimaryTreatment)) +  # , colour=supp, group=supp))
        theme_classic() +
        geom_errorbar(aes(ymin=mean.vstExp-se.vsdtExp, ymax=mean.vstExp+se.vsdtExp), colour="black", width=.1, position=pd) +
        geom_point(position=pd, size = 4, shape=21) +            
        xlab("Third pCO2 treatment") +
        ylab('') +                 # note the mean was first by sample ID THEN by treatment
        scale_fill_manual(values=c("white","grey50")) +
        # scale_color_manual(values=c("white","grey50")) +
        ggtitle(NULL) +
        # expand_limits(y=0) +                                                    # Expand y range
        #scale_y_continuous(limits=c((min_p1), (max_p1))) +
        theme(text = element_text(size=15)) +
        theme(legend.position = "none") +
        scale_y_continuous(limits = c(ExpMin,ExpMax)) +
        facet_wrap(~SecondTreatment)
    } else { d21_plots <- plot.new() }
    
    pdf(paste("Analysis/Output/GO/WGCNA_goseq/PrimaryEffect_Figures/GOslim_venn_tables/Modules_Moderate/Moderate_BP_",substr(slim,1,13),"_vstExpression.pdf", sep =''), width=15, height=6)
    print(ggarrange(d7_plots, d14_plots, d21_plots,        
                    plotlist = NULL,
                    ncol = 3,
                    nrow = 1,
                    labels = NULL))
    dev.off()
  } # plotting for loop
  
  else {}
  
}















































# MOLECULAR FUNCTION:  HIGHER expression under Moderate CONDITIONING ========================================================== #
MF_Mod_GOslim <- MF_Mod %>% 
  group_by(slim_term) %>% 
  tally() %>% 
  dplyr::filter(n == 3) # create the dataset to loop Venn diagrams 

for (i in 1:nrow(MF_Mod_GOslim)) {
  slim <- MF_Mod_GOslim[i,1] # call the term for the loop - creates Venn and Figures 
  
  # prep data
  slim_MF <- MF_Mod %>% dplyr::filter(slim_term %in% slim)
  slim_MF$slimSplitted  <- strsplit(as.character(slim_MF$Gene_IDs), ";") # split the Gene IDs for Venn 
  slim_d7  <- data.frame(Day = 'Day7', Gene.ID = unlist(slim_MF[1,10]))  # name the first row
  slim_d14 <- data.frame(Day = 'Day14', Gene.ID = unlist(slim_MF[2,10])) # name the second row
  slim_d21 <- data.frame(Day = 'Day21', Gene.ID = unlist(slim_MF[3,10])) # name the third row
  slim_df  <- list(
    Day7_YellowMod          = slim_d7$Gene.ID, 
    Day14_BlackMod          = slim_d14$Gene.ID, 
    Day21_YellowMod        = slim_d21$Gene.ID) # create a list to call in using ggvenn 
  
  # Venn diagram
  VennDiag <- ggvenn(slim_df,  # Venn diagram
                     fill_color = c("white", "white", "white"),
                     stroke_size = 0.5, set_name_size = 4) +
    ggtitle(paste("WGCNA: ",  slim, sep = ''))
  pdf(paste("Analysis/Output/GO/WGCNA_goseq/PrimaryEffect_Figures/GOslim_venn_tables/Modules_Moderate/Moderate_MF_",substr(slim,1,13),".pdf", sep =''))
  print(VennDiag)
  # grid.arrange(Day7.volcano.primary, Day14.volcano.primary,ncol=1, nrow=2, clip="off")
  dev.off()
  
  # write a csv for the corresponding Venn
  bind_data <- rbind(slim_d7, slim_d14, slim_d21)
  colnames(bind_data)[2] <- 'Gene_IDs'
  merged_with_annot <- merge(bind_data, annot.condenced[c(1,2,4)], by = 'Gene_IDs') %>% 
    group_by(Gene_IDs, gene) %>% 
    mutate(count = n()) %>% 
    dplyr::arrange(desc(count))
  merged_with_annot <- as.data.frame(merged_with_annot)
  # write csv
  write.csv(merged_with_annot, paste("Analysis/Output/GO/WGCNA_goseq/PrimaryEffect_Figures/GOslim_venn_tables/Modules_Moderate/Moderate_MF_",substr(slim,1,13),".csv", sep =''))
  
  panelfig_data <- merged_with_annot %>% dplyr::filter(count >=2)
  
  
  if (nrow(panelfig_data) > 0) {
    Day7.ExpVST_GOIs <- Day7.ExpVST %>% dplyr::select('Sample.Name', 'Primary_Treatment', 'Second_Treament', any_of(panelfig_data$Gene_IDs))
    Day7.ExpVST_GOIs$group <- paste(Day7.ExpVST_GOIs$Primary_Treatment , Day7.ExpVST_GOIs$Second_Treament , sep='')
    
    Day14.ExpVST_GOIs <- Day14.ExpVST %>% dplyr::select('Sample.Name', 'Primary_Treatment', 'Second_Treament', any_of(panelfig_data$Gene_IDs))
    Day14.ExpVST_GOIs$group <- paste(Day14.ExpVST_GOIs$Primary_Treatment , Day14.ExpVST_GOIs$Second_Treament , sep='')
    
    Day21.ExpVST_GOIs <- Day21.ExpVST %>% dplyr::select('Sample.Name', 'Primary_Treatment', 'Second_Treament', 'Third_Treatment', any_of(panelfig_data$Gene_IDs))
    Day21.ExpVST_GOIs$group <- paste(Day21.ExpVST_GOIs$Primary_Treatment , Day21.ExpVST_GOIs$Second_Treament , Day21.ExpVST_GOIs$Third_Treatment,  sep='')
    
    
    # ===================================================================================
    # Day 7 data prep for figures
    #
    # ===================================================================================
    # reshape the data frame and merge geneIDs with tragetGOIs to add back the gene titles - use this to name figures downstream
    Day7.ExpVST_GOIs_MELT <- reshape2::melt(Day7.ExpVST_GOIs, id=(c('Sample.Name', 'Primary_Treatment', 'Second_Treament', 'group'))) # melt using reshape2
    names(Day7.ExpVST_GOIs_MELT)[(5:6)] <- c('Gene_IDs', 'vst_Expression') # change column names
    Day7_ExpVst_Master <- merge(panelfig_data, Day7.ExpVST_GOIs_MELT, by = 'Gene_IDs')
    
    # calc the m an expressio by gene ID (add gene title in group but this is the same unique level as GeneID - review target_GOIs)
    if(nrow(as.data.frame(unique(Day7_ExpVst_Master$Gene_IDs))) == 1) {
      Day7_meanExpr <- Day7_ExpVst_Master %>% 
        dplyr::filter(Day %in% 'Day7') %>% # if panelfig_data is set to >=2, this will plot the mean SE of the total sum of shared genes (as opposed to the exact cetner of the Venn)
        dplyr::select(c('Sample.Name','group', 'vst_Expression', 'Gene_IDs', 'gene')) %>% 
        group_by(group) %>%
        dplyr::summarize(mean.vstExp = mean(vst_Expression), 
                         sd.vsdtExp = sd(vst_Expression),
                         n = n(), 
                         se.vsdtExp = sd.vsdtExp/sqrt(n))
    } else if ( nrow(as.data.frame(unique(Day7_ExpVst_Master$Gene_IDs))) >= 2) {
      Day7_meanExpr <- Day7_ExpVst_Master %>% 
        dplyr::filter(Day %in% 'Day7') %>% # if panelfig_data is set to >=2, this will plot the mean SE of the total sum of shared genes (as opposed to the exact cetner of the Venn)
        dplyr::select(c('Sample.Name','group', 'vst_Expression', 'Gene_IDs', 'gene')) %>% 
        group_by(Sample.Name,group) %>%
        dplyr::summarize(mean.vstbygene = mean(vst_Expression),number_genes = n()) %>% 
        group_by(group) %>% 
        dplyr::summarize(mean.vstExp = mean(mean.vstbygene),
                         sd.vsdtExp = sd(mean.vstbygene),
                         n = n(), 
                         se.vsdtExp = sd.vsdtExp/sqrt(n))
    } else { Day7_meanExpr <- as.data.frame('NULL') }
    
    # create treatment groups 
    Day7_meanExpr$PrimaryTreatment <- substr(Day7_meanExpr$group, 1,1) # primary
    Day7_meanExpr$SecondTreatment <- substr(Day7_meanExpr$group, 2,2) # second
    
    # ===================================================================================
    # Day 14 data prep for figures
    #
    # ===================================================================================
    # reshape the data frame and merge geneIDs with tragetGOIs to add back the gene titles - use this to name figures downstream
    Day14.ExpVST_GOIs_MELT <- reshape2::melt(Day14.ExpVST_GOIs, id=(c('Sample.Name', 'Primary_Treatment', 'Second_Treament', 'group'))) # melt using reshape2
    names(Day14.ExpVST_GOIs_MELT)[(5:6)] <- c('Gene_IDs', 'vst_Expression') # change column names
    Day14_ExpVst_Master <- merge(panelfig_data, Day14.ExpVST_GOIs_MELT, by = 'Gene_IDs')
    
    # calc the m an expressio by gene ID (add gene title in group but this is the same unique level as GeneID - review target_GOIs)
    if(nrow(as.data.frame(unique(Day14_ExpVst_Master$Gene_IDs))) == 1) {
      Day14_meanExpr <- Day14_ExpVst_Master %>% 
        dplyr::filter(Day %in% 'Day14') %>% # if panelfig_data is set to >=2, this will plot the mean SE of the total sum of shared genes (as opposed to the exact cetner of the Venn)
        dplyr::select(c('Sample.Name','group', 'vst_Expression', 'Gene_IDs', 'gene')) %>% 
        group_by(group) %>%
        dplyr::summarize(mean.vstExp = mean(vst_Expression), 
                         sd.vsdtExp = sd(vst_Expression),
                         n = n(), 
                         se.vsdtExp = sd.vsdtExp/sqrt(n))
    } else if ( nrow(as.data.frame(unique(Day14_ExpVst_Master$Gene_IDs))) >= 2) {
      Day14_meanExpr <- Day14_ExpVst_Master %>% 
        dplyr::filter(Day %in% 'Day14') %>% # if panelfig_data is set to >=2, this will plot the mean SE of the total sum of shared genes (as opposed to the exact cetner of the Venn)
        dplyr::select(c('Sample.Name','group', 'vst_Expression', 'Gene_IDs', 'gene')) %>% 
        group_by(Sample.Name,group) %>%
        dplyr::summarize(mean.vstbygene = mean(vst_Expression),number_genes = n()) %>% 
        group_by(group) %>% 
        dplyr::summarize(mean.vstExp = mean(mean.vstbygene),
                         sd.vsdtExp = sd(mean.vstbygene),
                         n = n(), 
                         se.vsdtExp = sd.vsdtExp/sqrt(n))
    } else { Day14_meanExpr <- as.data.frame('NULL') }
    
    # create treatment groups 
    Day14_meanExpr$PrimaryTreatment <- substr(Day14_meanExpr$group, 1,1) # primary
    Day14_meanExpr$SecondTreatment <- substr(Day14_meanExpr$group, 2,2) # second
    
    # ===================================================================================
    # Day 21 data prep for figures
    #
    # ===================================================================================
    # reshape the data frame and merge geneIDs with tragetGOIs to add back the gene titles - use this to name figures downstream
    Day21.ExpVST_GOIs_MELT <- reshape2::melt(Day21.ExpVST_GOIs, id=(c('Sample.Name', 'Primary_Treatment', 'Second_Treament', 'Third_Treatment','group'))) # melt using reshape2
    names(Day21.ExpVST_GOIs_MELT)[(6:7)] <- c('Gene_IDs', 'vst_Expression') # change column names
    Day21_ExpVst_Master <- merge(panelfig_data, Day21.ExpVST_GOIs_MELT, by = 'Gene_IDs')
    
    
    # calc the m an expressio by gene ID (add gene title in group but this is the same unique level as GeneID - review target_GOIs)
    if(nrow(as.data.frame(unique(Day21_ExpVst_Master$Gene_IDs))) == 1) {
      Day21_meanExpr <- Day21_ExpVst_Master %>% 
        dplyr::filter(Day %in% 'Day21') %>% # if panelfig_data is set to >=2, this will plot the mean SE of the total sum of shared genes (as opposed to the exact cetner of the Venn)
        dplyr::select(c('Sample.Name','group', 'vst_Expression', 'Gene_IDs', 'gene')) %>% 
        group_by(group) %>%
        dplyr::summarize(mean.vstExp = mean(vst_Expression), 
                         sd.vsdtExp = sd(vst_Expression),
                         n = n(), 
                         se.vsdtExp = sd.vsdtExp/sqrt(n))
    } else if ( nrow(as.data.frame(unique(Day21_ExpVst_Master$Gene_IDs))) >= 2) {
      Day21_meanExpr <- Day21_ExpVst_Master %>% 
        dplyr::filter(Day %in% 'Day21') %>% # if panelfig_data is set to >=2, this will plot the mean SE of the total sum of shared genes (as opposed to the exact cetner of the Venn)
        dplyr::select(c('Sample.Name','group', 'vst_Expression', 'Gene_IDs', 'gene')) %>% 
        group_by(Sample.Name,group) %>%
        dplyr::summarize(mean.vstbygene = mean(vst_Expression),number_genes = n()) %>% 
        group_by(group) %>% 
        dplyr::summarize(mean.vstExp = mean(mean.vstbygene),
                         sd.vsdtExp = sd(mean.vstbygene),
                         n = n(), 
                         se.vsdtExp = sd.vsdtExp/sqrt(n))
    } else { Day21_meanExpr <- as.data.frame(NULL) }
    
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
    # for(i in 1:nrow(target_GOIs)) {
    #   
    all <- rbind(Day7_meanExpr[,c(1:5)], Day14_meanExpr[,c(1:5)], Day21_meanExpr[,c(1:5)])
    ExpMin <- round( ( min(all$mean.vstExp) - max(all$se.vsdtExp)  ), digits =1)
    ExpMax <- round( ( max(all$mean.vstExp) + max(all$se.vsdtExp)   ), digits =1)

    if(nrow(Day7_meanExpr) > 1 ) {
      d7_plots <- Day7_meanExpr %>% 
        ggplot(aes(x=SecondTreatment, y=mean.vstExp, fill=PrimaryTreatment)) +  # , colour=supp, group=supp))
        theme_classic() +
        geom_errorbar(aes(ymin=mean.vstExp-se.vsdtExp, ymax=mean.vstExp+se.vsdtExp), colour="black", width=.1, position=pd) +
        geom_point(position=pd, size = 4, shape=21) +            
        xlab("Second pCO2 treatment") +
        ylab('') +                 # note the mean was first by sample ID THEN by treatment
        scale_fill_manual(values=c("white","grey50")) +
        # scale_color_manual(values=c("white","grey50")) +
        ggtitle(paste("Day 7:",substr(slim,1,13), " (N = ", nrow(unique(as.data.frame((Day7_ExpVst_Master %>% dplyr::filter(Day %in% 'Day7'))$Gene_IDs))), ")", sep='')) +
        # expand_limits(y=0) +                                                    # Expand y range
        #scale_y_continuous(limits=c((min_p1), (max_p1))) +
        theme(text = element_text(size=15)) +
        scale_y_continuous(limits = c(ExpMin,ExpMax)) +
        theme(legend.position = "none")
    } else { d7_plots <- plot.new() }
    # 
    #   
    if(nrow(Day14_meanExpr) > 1 ) {
      d14_plots <- Day14_meanExpr %>% 
        # dplyr::filter(genes %in% target_GOIs[i,1]) %>% 
        ggplot(aes(x=SecondTreatment, y=mean.vstExp, fill=PrimaryTreatment)) +  # , colour=supp, group=supp))
        theme_classic() +
        geom_errorbar(aes(ymin=mean.vstExp-se.vsdtExp, ymax=mean.vstExp+se.vsdtExp), colour="black", width=.1, position=pd) +
        geom_point(position=pd, size = 4, shape=21) +            
        xlab("Second pCO2 treatment") +
        ylab("Gene Expression (mean±SE VST transformed)") +                 # note the mean was first by sample ID THEN by treatment
        scale_fill_manual(values=c("white","grey50")) +
        # scale_color_manual(values=c("white","grey50")) +
        ggtitle(paste("Day 14:",substr(slim,1,13), " (N = ", nrow(unique(as.data.frame((Day14_ExpVst_Master %>% dplyr::filter(Day %in% 'Day14'))$Gene_IDs))), ")", sep='')) +
        # expand_limits(y=0) +                                                    # Expand y range
        #scale_y_continuous(limits=c((min_p1), (max_p1))) +
        theme(text = element_text(size=15)) +
        scale_y_continuous(limits = c(ExpMin,ExpMax)) +
        theme(legend.position = "none")
    } else { d14_plots <- plot.new() }
    # 
    # 
    #       
    if(nrow(Day21_meanExpr) > 1 ) {
      d21_plots <- Day21_meanExpr %>% 
        # dplyr::filter(genes %in% target_GOIs[i,1]) %>% 
        ggplot(aes(x=ThirdTreatment, y=mean.vstExp, fill=PrimaryTreatment)) +  # , colour=supp, group=supp))
        theme_classic() +
        geom_errorbar(aes(ymin=mean.vstExp-se.vsdtExp, ymax=mean.vstExp+se.vsdtExp), colour="black", width=.1, position=pd) +
        geom_point(position=pd, size = 4, shape=21) +            
        xlab("Third pCO2 treatment") +
        ylab('') +                 # note the mean was first by sample ID THEN by treatment
        scale_fill_manual(values=c("white","grey50")) +
        # scale_color_manual(values=c("white","grey50")) +
        ggtitle(NULL) +
        # expand_limits(y=0) +                                                    # Expand y range
        #scale_y_continuous(limits=c((min_p1), (max_p1))) +
        theme(text = element_text(size=15)) +
        theme(legend.position = "none") +
        scale_y_continuous(limits = c(ExpMin,ExpMax)) +
        facet_wrap(~SecondTreatment)
    } else { d21_plots <- plot.new() }
    
    pdf(paste("Analysis/Output/GO/WGCNA_goseq/PrimaryEffect_Figures/GOslim_venn_tables/Modules_Moderate/Moderate_MF_",substr(slim,1,13),"_vstExpression.pdf", sep =''), width=15, height=6)
    print(ggarrange(d7_plots, d14_plots, d21_plots,        
                    plotlist = NULL,
                    ncol = 3,
                    nrow = 1,
                    labels = NULL))
    dev.off()
  } # plotting for loop
  
  else {}
  
}

  
  
















#===================================================================================================
#
#   Unnest the genes in each GOslim to explore Functiona envrichment of interest! 
#  load the GOslim csvs, unnest (tidyr) the gene ID list and merge with Gene annnotation
# 
#   
#
#===================================================================================================
library(stringr)
library(DESeq2)
# BIOLOGICAL PROCESS  GO  SLIMS
d7_slimBP_brownModule   <- read.csv("Analysis/Output/GO/WGCNA_goseq/subseq_treatments_all/Day7/GOslim_BiolProc_brownModule.csv")
d7_slimBP_greenModule   <- read.csv("Analysis/Output/GO/WGCNA_goseq/subseq_treatments_all/Day7/GOslim_BiolProc_greenModule.csv")
d7_slimBP_yellowModule  <- read.csv("Analysis/Output/GO/WGCNA_goseq/subseq_treatments_all/Day7/GOslim_BiolProc_yellowModule.csv")

d14_slimBP_blackModule   <- read.csv("Analysis/Output/GO/WGCNA_goseq/subseq_treatments_all/Day14/GOslim_BiolProc_blackModule.csv")
d14_slimBP_brownModule   <- read.csv("Analysis/Output/GO/WGCNA_goseq/subseq_treatments_all/Day14/GOslim_BiolProc_brownModule.csv")
d14_slimBP_magentaModule <- read.csv("Analysis/Output/GO/WGCNA_goseq/subseq_treatments_all/Day14/GOslim_BiolProc_magentaModule.csv")
d14_slimBP_pinkModule    <- read.csv("Analysis/Output/GO/WGCNA_goseq/subseq_treatments_all/Day14/GOslim_BiolProc_pinkModule.csv")

d21_slimBP_blackModule    <- read.csv("Analysis/Output/GO/WGCNA_goseq/subseq_treatments_all/Day21/GOslim_BiolProc_blackModule.csv")
d21_slimBP_blueModule     <- read.csv("Analysis/Output/GO/WGCNA_goseq/subseq_treatments_all/Day21/GOslim_BiolProc_blueModule.csv")
d21_slimBP_magentaModule  <- read.csv("Analysis/Output/GO/WGCNA_goseq/subseq_treatments_all/Day21/GOslim_BiolProc_magentaModule.csv")
d21_slimBP_pinkModule     <- read.csv("Analysis/Output/GO/WGCNA_goseq/subseq_treatments_all/Day21/GOslim_BiolProc_pinkModule.csv")
d21_slimBP_redModule      <- read.csv("Analysis/Output/GO/WGCNA_goseq/subseq_treatments_all/Day21/GOslim_BiolProc_redModule.csv")
d21_slimBP_yellowModule   <- read.csv("Analysis/Output/GO/WGCNA_goseq/subseq_treatments_all/Day21/GOslim_BiolProc_yellowModule.csv")

# MOLECULAR FUNCTION GO  SLIMS
d7_slimMF_brownModule   <- read.csv("Analysis/Output/GO/WGCNA_goseq/subseq_treatments_all/Day7/GOslim_MolFunction_brownModule.csv")
d7_slimMF_greenModule   <- read.csv("Analysis/Output/GO/WGCNA_goseq/subseq_treatments_all/Day7/GOslim_MolFunction_greenModule.csv")
d7_slimMF_yellowModule  <- read.csv("Analysis/Output/GO/WGCNA_goseq/subseq_treatments_all/Day7/GOslim_MolFunction_yellowModule.csv")

d14_slimMF_blackModule   <- read.csv("Analysis/Output/GO/WGCNA_goseq/subseq_treatments_all/Day14/GOslim_MolFunction_blackModule.csv")
d14_slimMF_brownModule   <- read.csv("Analysis/Output/GO/WGCNA_goseq/subseq_treatments_all/Day14/GOslim_MolFunction_brownModule.csv")
d14_slimMF_magentaModule <- read.csv("Analysis/Output/GO/WGCNA_goseq/subseq_treatments_all/Day14/GOslim_MolFunction_magentaModule.csv")
d14_slimMF_pinkModule    <- read.csv("Analysis/Output/GO/WGCNA_goseq/subseq_treatments_all/Day14/GOslim_MolFunction_pinkModule.csv")

d21_slimMF_blackModule    <- read.csv("Analysis/Output/GO/WGCNA_goseq/subseq_treatments_all/Day21/GOslim_MolFunction_blackModule.csv")
d21_slimMF_blueModule     <- read.csv("Analysis/Output/GO/WGCNA_goseq/subseq_treatments_all/Day21/GOslim_MolFunction_blueModule.csv")
d21_slimMF_magentaModule  <- read.csv("Analysis/Output/GO/WGCNA_goseq/subseq_treatments_all/Day21/GOslim_MolFunction_magentaModule.csv")
d21_slimMF_pinkModule     <- read.csv("Analysis/Output/GO/WGCNA_goseq/subseq_treatments_all/Day21/GOslim_MolFunction_pinkModule.csv")
d21_slimMF_redModule      <- read.csv("Analysis/Output/GO/WGCNA_goseq/subseq_treatments_all/Day21/GOslim_MolFunction_redModule.csv")
d21_slimMF_yellowModule   <- read.csv("Analysis/Output/GO/WGCNA_goseq/subseq_treatments_all/Day21/GOslim_MolFunction_yellowModule.csv")

# Load the annotation file 
Geoduck_annotation      <- read.delim2(file="C:/Users/samjg/Documents/My_Projects/Pgenerosa_TagSeq_Metabolomics/TagSeq/Seq_details/Panopea-generosa-genes-annotations.txt", header=F)
annot.condenced         <- Geoduck_annotation[,c(1,5,7)] # build annotation file to merge with master: gene ID, uniprot ID, gene terms with EC (for KEGG)
names(annot.condenced)  <- c('Gene_IDs', 'Uniprot', 'Gene_term_EC') # RENAME THE COLUMNS 
annot.condenced$gene    <- str_extract(annot.condenced$Gene_term_EC, "[^(]+") # condence the gene ID to all characters before the first occurance of a open parenthese '('
# View(annot.condenced) # view the data 


# MASTER FILES WITH AND WITHOUT THE GENE ANNOTATION FOR EACH GENE AND SLIM TERM
# Biological Process ====================================== #
MasterBP_GOslim <- rbind(d7_slimBP_brownModule, d7_slimBP_greenModule, d7_slimBP_yellowModule, # rbind all sig modules - biological process
                              d14_slimBP_blackModule, d14_slimBP_brownModule, d14_slimBP_magentaModule, d14_slimBP_pinkModule,
                              d21_slimBP_blackModule, d21_slimBP_blueModule, d21_slimBP_magentaModule, d21_slimBP_pinkModule, d21_slimBP_redModule, d21_slimBP_yellowModule)
MasterBP_GOslim_OM.GOterms          <- MasterBP_GOslim[-c(1,4:6)] # ommit undesired columsn: rows, GO count, GO term, GO ID list
MasterBP_GOslim_OM.GOterms$Gene_IDs <- strsplit(MasterBP_GOslim_OM.GOterms$Gene_IDs, ";") # convert the ';' delimited rows of gene IDs into a list
MasterBP_unnest                     <- tidyr::unnest(MasterBP_GOslim_OM.GOterms, Gene_IDs) # here we go! unnest the list of genes - this will duplicate the other columns for each gene ID in common cell
MasterBP_unnest$Day                 <- str_extract(MasterBP_unnest$module_day, "[^_]+") # use stringr to extract dat from module_day
MasterBP_unnest$moduleColor         <- word(MasterBP_unnest$module_day, 2, sep = "_")
# highly recommend using View() for the unnested v. the previous master files to make sure this worked correctly...
# View(MasterBP_GOslim_OM.GOterms)
# View(MasterBP_unnest)
MasterBP_FINAL <- merge(MasterBP_unnest, annot.condenced, by = 'Gene_IDs') # save this!! 
View(MasterBP_FINAL)
# Molecular Function  ====================================== #
MasterMF_GOslim <- rbind(d7_slimMF_brownModule, d7_slimMF_greenModule, d7_slimMF_yellowModule, # rbind all sig modules - biological process
                         d14_slimMF_blackModule, d14_slimMF_brownModule, d14_slimMF_magentaModule, d14_slimMF_pinkModule,
                         d21_slimMF_blackModule, d21_slimMF_blueModule, d21_slimMF_magentaModule, d21_slimMF_pinkModule, d21_slimMF_redModule, d21_slimMF_yellowModule)
MasterMF_GOslim_OM.GOterms          <- MasterMF_GOslim[-c(1,4:6)] # ommit undesired columsn: rows, GO count, GO term, GO ID list
MasterMF_GOslim_OM.GOterms$Gene_IDs <- strsplit(MasterMF_GOslim_OM.GOterms$Gene_IDs, ";") # convert the ';' delimited rows of gene IDs into a list
MasterMF_unnest                     <- tidyr::unnest(MasterMF_GOslim_OM.GOterms, Gene_IDs) # here we go! unnest the list of genes - this will duplicate the other columns for each gene ID in common cell
MasterMF_unnest$Day                 <- str_extract(MasterMF_unnest$module_day, "[^_]+") # use stringr to extract dat from module_day
MasterMF_unnest$moduleColor         <- word(MasterMF_unnest$module_day, 2, sep = "_")
# highly recommend using View() for the unnested v. the previous master files to make sure this worked correctly...
# View(MasterMF_GOslim_OM.GOterms)
# View(MasterMF_unnest)
MasterMF_FINAL <- merge(MasterMF_unnest, annot.condenced, by = 'Gene_IDs')
# save the master files 
write.csv(MasterBP_FINAL, file = paste("Analysis/Output/GO/WGCNA_goseq/BiologicalProcess_GOslim_Master.csv", sep ='')) # save csv file       
write.csv(MasterMF_FINAL, file = paste("Analysis/Output/GO/WGCNA_goseq/MolecularFunction_GOslim_Master.csv", sep ='')) # save csv file              








