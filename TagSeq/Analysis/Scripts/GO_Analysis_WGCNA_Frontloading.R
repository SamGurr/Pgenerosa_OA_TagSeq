---
  # title: "GO_Analysis_WGCNA_all"
  # author: "Samuel Gurr"
  # date: "April 12, 2021"
---
  
  
# SET WORKING DIRECTORY AND LOAD DATA
setwd("C:/Users/samjg/Documents/Github_repositories/Pgenerosa_TagSeq_Metabolomics/TagSeq/")

# Load libraries 
library(dplyr)
library(tidyr)
library(goseq) # BiocManager::install('goseq')
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
Geoduck_annotation <- read.delim2(file="C:/Users/samjg/Documents/Github_repositories/Pgenerosa_TagSeq_Metabolomics/TagSeq/Seq_details/Panopea-generosa-genes-annotations.txt", header=F) # Load themaster Pgenerosa gene list 
Pgen_GOterms <- Geoduck_annotation %>% dplyr::select(c('V1','V8')) # select only two columns - those with the gene IDs and those with the GO terms
Pgen_GOterms2 <- strsplit(Pgen_GOterms$V8, split = "; ") # create a string splitting by delimiter '; ' - view the data to see that this separates each GO term entry in the string
Pgen_GOterms2 <- data.frame(gene.ID = rep(Pgen_GOterms$V1, sapply(Pgen_GOterms2, length)), Go.terms = unlist(Pgen_GOterms2)) # create new dataframe 'Pgen_GOterms2' listing genes for each GO term (MUCH longer!)
Pgen_GOterms2 <- na.omit(Pgen_GOterms2) # ommit the NAs  - genes without GO annotation


# FRrntloaded gene sets
d7_frontloaded_moderate  <- read.csv("Analysis/Output/WGCNA/subseq_treatments_all/Frontloading/Day7_ModerateSevere_FrontloadedGenes.csv") %>% dplyr::filter(Frontloaded_Moderate %in% 'frontloaded') %>% dplyr::rename(geneSymbol = Gene)
d7_frontloaded_severe    <- read.csv("Analysis/Output/WGCNA/subseq_treatments_all/Frontloading/Day7_ModerateSevere_FrontloadedGenes.csv") %>% dplyr::filter(Frontloaded_Severe %in% 'frontloaded') %>% dplyr::rename(geneSymbol = Gene)
d21_frontloaded          <- read.csv("Analysis/Output/WGCNA/subseq_treatments_all/Frontloading/Day21_Moderate_FrontloadedGenes.csv") %>% dplyr::filter(Frontloaded %in% 'frontloaded') %>% dplyr::rename(geneSymbol = Gene)
d7.d21_frontloaded_moderate          <- read.csv("Analysis/Output/WGCNA/subseq_treatments_all/Frontloading/Day7_Day21_Shared_Moderate_FrontloadedGenes.csv") 



slim <- getOBOCollection("http://current.geneontology.org/ontology/subsets/goslim_generic.obo") #get GO database - # call goslim_generic.obo terms as 'slim'

#===================================================================================== # 
# LOAD DATA - Raw Count Data - FILTERED COUNT MATRICES USED IN WGCNA- 10CPM in 50% of samples
#
#===================================================================================== # 

# day7 filtered 10cpm in 50% samples ----------------------------- # 
Day7_all.counts <- read.csv(file="C:/Users/samjg/Documents/Github_repositories/Pgenerosa_TagSeq_Metabolomics/TagSeq/Analysis/Data/Filtered_Counts/10cpm_50perc/day7.counts.filtered_10cpm50perc.csv", sep=',', header=TRUE) 
colnames(Day7_all.counts)[1] <- "gene.ID"# rename Pgen gene ID column


# day21 filtered 10cpm in 50% samples ----------------------------- # 
Day21_all.counts <- read.csv(file="C:/Users/samjg/Documents/Github_repositories/Pgenerosa_TagSeq_Metabolomics/TagSeq/Analysis/Data/Filtered_Counts/10cpm_50perc/day21.counts.filtered_10cpm50perc.csv", sep=',', header=TRUE) 
colnames(Day21_all.counts)[1] <- "gene.ID"# rename Pgen gene ID column



#===================================================================================== # 
# LOAD DATA - goseq; load the annotation and prepare the fouressentail steps for goseq
#
#===================================================================================== #
#Panopea generosa - load .fna ('Geoduck_annotation') and foramt GO terms ('Geoduck_GOterms') and vectors
Geoduck_annotation <- read.delim2(file="C:/Users/samjg/Documents/Github_repositories/Pgenerosa_TagSeq_Metabolomics/TagSeq/Seq_details/Panopea-generosa-genes-annotations.txt", header=F)

# build annotation file to merge with the mean LFC tables
annot.condenced <- Geoduck_annotation[,c(1,3:9)]
annot.condenced$gene.length <- annot.condenced$V4 - annot.condenced$V3
annot.condenced <- annot.condenced[,-c(2,3)]
names(annot.condenced) <- c('Gene.ID', 'Uniprot', 'HGNC', 'fxn', 'Go.terms', 'Go.fxns','gene.length')


Geoduck_annotation <- read.delim2(file="C:/Users/samjg/Documents/Github_repositories/Pgenerosa_TagSeq_Metabolomics/TagSeq/Seq_details/Panopea-generosa-genes-annotations.txt", header=F)
# build annotation file to merge with the mean LFC tables
Pgen_condenced <- Geoduck_annotation[,c(1,7)] # load just the PGEN ID and the putative gene terms
Pgen_condenced$Gene_term    <- sub(" \\(EC.*", "", Pgen_condenced$V7)  # str_extract(annot.condenced$V7, "[^(]+")
Pgen_reference <- na.omit(Pgen_condenced[c(1,3)])
names(Pgen_reference)  <- c('PgenIDs', 'Gene_terms')


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

# (3) Gene length 
# length vector  
GO_gene.length <- Geoduck_annotation %>% dplyr::mutate(length = V4-V3) %>%  dplyr::select(c("V1","length"))
names(GO_gene.length)[1] <- "gene.ID"
# merge length with counts data
length_vector   <- GO_gene.length$length


#===================================================================================== # 
# Frontloaded genes from Day7 module yellow (higher by preexposded clams) and in response to MODERATELY-elevated pCO2
#===================================================================================== # 
FrontloadedDay7_moderate                       <- d7_frontloaded_moderate[2]
names(FrontloadedDay7_moderate)[1]             <- "Gene.ID" # 162 genws in the green module 
FrontloadModD7_integer                         <- as.integer(GO_unique.genes.all %in% (FrontloadedDay7_moderate$Gene.ID)) # w/o day-specific ID vector
names(FrontloadModD7_integer)                  =  GO_unique.genes.all # rename
pwf_FrontloadModD7                             <- nullp(FrontloadModD7_integer,    id=GO_unique.genes.all, bias.data=length_vector) # make figure margins large enough for this to run...
goseq_FrontloadModD7                           <- goseq(pwf_FrontloadModD7, gene2cat=GO.terms, test.cats=c("GO:CC", "GO:BP", "GO:MF"), method="Wallenius", use_genes_without_cat=TRUE)
GO.05.a_FrontloadModD7                         <- goseq_FrontloadModD7$category[goseq_FrontloadModD7$over_represented_pvalue<.05] # change twice here
GO.05_FrontloadModD7                           <- data.frame(GO.05.a_FrontloadModD7 )
colnames(GO.05_FrontloadModD7 )                <- c("category")
GO.05_FrontloadModD7                           <- merge(GO.05_FrontloadModD7 , goseq_FrontloadModD7, by="category") # change here
GO.05_FrontloadModD7                           <- GO.05_FrontloadModD7[order(GO.05_FrontloadModD7$ontology, GO.05_FrontloadModD7$over_represented_pvalue,-GO.05_FrontloadModD7$numDEInCat),]
GO.05_FrontloadModD7$term                     <- as.factor(GO.05_FrontloadModD7$term)
GO.05_FrontloadModD7$Day                      <- "Day7"
# remove Biological Process GO terms with < 10 genes in the module  (with that term) and ommit Molecular Function terms with < 3 genes in the module (with that term)
GO.05_filtered_FrontloadedDay7_moderate <- GO.05_FrontloadModD7  %>% filter(!(ontology %in% "CC"), !(numDEInCat<=2 & ontology == "BP"), !(numDEInCat<=2 & ontology == "MF"))
View(GO.05_filtered_FrontloadedDay7_moderate)

#===================================================================================== # 
# Frontloaded genes from Day7 module yellow (higher by preexposded clams) and in response to SEVERELY-elevated pCO2
#===================================================================================== # 
FrontloadedDay7_severe                       <- d7_frontloaded_severe[2]
names(FrontloadedDay7_severe)[1]             <- "Gene.ID" # 162 genws in the green module 
FrontloadSevD7_integer                         <- as.integer(GO_unique.genes.all %in% (FrontloadedDay7_severe$Gene.ID)) # w/o day-specific ID vector
names(FrontloadSevD7_integer)                  =  GO_unique.genes.all # rename
pwf_FrontloadSevD7                             <- nullp(FrontloadSevD7_integer,    id=GO_unique.genes.all, bias.data=length_vector) # make figure margins large enough for this to run...
goseq_FrontloadSevD7                           <- goseq(pwf_FrontloadSevD7, gene2cat=GO.terms, test.cats=c("GO:CC", "GO:BP", "GO:MF"), method="Wallenius", use_genes_without_cat=TRUE)
GO.05.a_FrontloadSevD7                         <- goseq_FrontloadSevD7$category[goseq_FrontloadSevD7$over_represented_pvalue<.05] # change twice here
GO.05_FrontloadSevD7                           <- data.frame(GO.05.a_FrontloadSevD7 )
colnames(GO.05_FrontloadSevD7 )                <- c("category")
GO.05_FrontloadSevD7                           <- merge(GO.05_FrontloadSevD7 , goseq_FrontloadSevD7, by="category") # change here
GO.05_FrontloadSevD7                           <- GO.05_FrontloadSevD7[order(GO.05_FrontloadSevD7$ontology, GO.05_FrontloadSevD7$over_represented_pvalue,-GO.05_FrontloadSevD7$numDEInCat),]
GO.05_FrontloadSevD7$term                     <- as.factor(GO.05_FrontloadSevD7$term)
GO.05_FrontloadSevD7$Day                      <- "Day7"
# remove Biological Process GO terms with < 10 genes in the module  (with that term) and ommit Molecular Function terms with < 3 genes in the module (with that term)
GO.05_filtered_FrontloadedDay7_severe <- GO.05_FrontloadSevD7  %>% filter(!(ontology %in% "CC"), !(numDEInCat<=2 & ontology == "BP"), !(numDEInCat<=2 & ontology == "MF"))
View(GO.05_filtered_FrontloadedDay7_severe)


#===================================================================================== # 
# Frontloaded genes from Day21 module yellow (higher by preexposded clams) and in response to SEVERELY-elevated pCO2
#==================================================================================== # 
FrontloadedDay21                             <- d21_frontloaded[2]
names(FrontloadedDay21)[1]                   <- "Gene.ID" # 162 genws in the green module 
FrontloadD21_integer                         <- as.integer(GO_unique.genes.all %in% (FrontloadedDay21$Gene.ID)) # w/o day-specific ID vector
names(FrontloadD21_integer)                  =  GO_unique.genes.all # rename
pwf_FrontloadD21                             <- nullp(FrontloadD21_integer,    id=GO_unique.genes.all, bias.data=length_vector) # make figure margins large enough for this to run...
goseq_FrontloadD21                           <- goseq(pwf_FrontloadD21, gene2cat=GO.terms, test.cats=c("GO:CC", "GO:BP", "GO:MF"), method="Wallenius", use_genes_without_cat=TRUE)
GO.05.a_FrontloadD21                         <- goseq_FrontloadD21$category[goseq_FrontloadD21$over_represented_pvalue<.05] # change twice here
GO.05_FrontloadD21                           <- data.frame(GO.05.a_FrontloadD21 )
colnames(GO.05_FrontloadD21 )                <- c("category")
GO.05_FrontloadD21                           <- merge(GO.05_FrontloadD21 , goseq_FrontloadD21, by="category") # change here
GO.05_FrontloadD21                           <- GO.05_FrontloadD21[order(GO.05_FrontloadD21$ontology, GO.05_FrontloadD21$over_represented_pvalue,-GO.05_FrontloadD21$numDEInCat),]
GO.05_FrontloadD21$term                     <- as.factor(GO.05_FrontloadD21$term)
GO.05_FrontloadD21$Day                      <- "Day21"
# remove Biological Process GO terms with < 10 genes in the module  (with that term) and ommit Molecular Function terms with < 3 genes in the module (with that term)
GO.05_filtered_FrontloadedDay21 <- GO.05_FrontloadD21  %>% filter(!(ontology %in% "CC"), !(numDEInCat<=2 & ontology == "BP"), !(numDEInCat<=2 & ontology == "MF"))
View(GO.05_filtered_FrontloadedDay21)


#===================================================================================== # 
# Frontloaded genes from Day7 module yellow (higher by preexposded clams) SHARED BETWEEN MODERATE AND SEVERE FRONTLOADED GENES
#===================================================================================== # 

d7_frontloaded_shared <-  as.data.frame(rbind(d7_frontloaded_moderate[,c(2,3)], d7_frontloaded_severe[,c(2,3)]) %>% dplyr::group_by(geneSymbol) %>%  dplyr::summarise(count = n()) %>% dplyr::filter(count == 2))

FrontloadedDay7_shared                       <- d7_frontloaded_shared[1]
names(FrontloadedDay7_shared)[1]             <- "Gene.ID" # 162 genws in the green module 
FrontloadsharedD7_integer                         <- as.integer(GO_unique.genes.all %in% (FrontloadedDay7_shared$Gene.ID)) # w/o day-specific ID vector
names(FrontloadsharedD7_integer)                  =  GO_unique.genes.all # rename
pwf_FrontloadsharedD7                             <- nullp(FrontloadsharedD7_integer,    id=GO_unique.genes.all, bias.data=length_vector) # make figure margins large enough for this to run...
goseq_FrontloadsharedD7                           <- goseq(pwf_FrontloadsharedD7, gene2cat=GO.terms, test.cats=c("GO:CC", "GO:BP", "GO:MF"), method="Wallenius", use_genes_without_cat=TRUE)
GO.05.a_FrontloadsharedD7                         <- goseq_FrontloadsharedD7$category[goseq_FrontloadsharedD7$over_represented_pvalue<.05] # change twice here
GO.05_FrontloadsharedD7                           <- data.frame(GO.05.a_FrontloadsharedD7 )
colnames(GO.05_FrontloadsharedD7 )                <- c("category")
GO.05_FrontloadsharedD7                           <- merge(GO.05_FrontloadsharedD7 , goseq_FrontloadsharedD7, by="category") # change here
GO.05_FrontloadsharedD7                           <- GO.05_FrontloadsharedD7[order(GO.05_FrontloadsharedD7$ontology, GO.05_FrontloadsharedD7$over_represented_pvalue,-GO.05_FrontloadsharedD7$numDEInCat),]
GO.05_FrontloadsharedD7$term                     <- as.factor(GO.05_FrontloadsharedD7$term)
GO.05_FrontloadsharedD7$Day                      <- "Day7"
# remove Biological Process GO terms with < 10 genes in the module  (with that term) and ommit Molecular Function terms with < 3 genes in the module (with that term)
GO.05_filtered_FrontloadedDay7_shared <- GO.05_FrontloadsharedD7  %>% filter(!(ontology %in% "CC"), !(numDEInCat<=2 & ontology == "BP"), !(numDEInCat<=2 & ontology == "MF"))
View(GO.05_filtered_FrontloadedDay7_shared)




#===================================================================================== # 
# Frontloaded genes from Day7 module yellow (higher by preexposded clams) SHARED BETWEEN MODERATE AND SEVERE FRONTLOADED GENES
#===================================================================================== # 
View(FrontloadsharedD7_D21_mod_integer)
d7.d21_frontloaded_moderate                               <- d7.d21_frontloaded_moderate[2]
names(d7.d21_frontloaded_moderate)[1]                     <- "Gene.ID" # 162 genws in the green module 
FrontloadsharedD7_D21_mod_integer                         <- as.integer(GO_unique.genes.all %in% (d7.d21_frontloaded_moderate$Gene.ID)) # w/o day-specific ID vector
names(FrontloadsharedD7_D21_mod_integer)                  =  GO_unique.genes.all # rename
pwf_FrontloadsharedD7_D21_mod                             <- nullp(FrontloadsharedD7_D21_mod_integer,    id=GO_unique.genes.all, bias.data=length_vector) # make figure margins large enough for this to run...
goseq_FrontloadsharedD7_D21_mod                           <- goseq(pwf_FrontloadsharedD7_D21_mod, gene2cat=GO.terms, test.cats=c("GO:CC", "GO:BP", "GO:MF"), method="Wallenius", use_genes_without_cat=TRUE)
GO.05.a_FrontloadsharedD7_D21_mod                         <- goseq_FrontloadsharedD7_D21_mod$category[goseq_FrontloadsharedD7_D21_mod$over_represented_pvalue<.05] # change twice here
GO.05_FrontloadsharedD7_D21_mod                           <- data.frame(GO.05.a_FrontloadsharedD7_D21_mod )
colnames(GO.05_FrontloadsharedD7_D21_mod )                <- c("category")
GO.05_FrontloadsharedD7_D21_mod                           <- merge(GO.05_FrontloadsharedD7_D21_mod , goseq_FrontloadsharedD7_D21_mod, by="category") # change here
GO.05_FrontloadsharedD7_D21_mod                           <- GO.05_FrontloadsharedD7_D21_mod[order(GO.05_FrontloadsharedD7_D21_mod$ontology, GO.05_FrontloadsharedD7_D21_mod$over_represented_pvalue,-GO.05_FrontloadsharedD7_D21_mod$numDEInCat),]
GO.05_FrontloadsharedD7_D21_mod$term                     <- as.factor(GO.05_FrontloadsharedD7_D21_mod$term)
GO.05_FrontloadsharedD7_D21_mod$Day                      <- "D7_D21_mod"
# remove Biological Process GO terms with < 10 genes in the module  (with that term) and ommit Molecular Function terms with < 3 genes in the module (with that term)
GO.05_filtered_FrontloadedD7_D21_mod <- GO.05_FrontloadsharedD7_D21_mod  %>% filter(!(ontology %in% "CC"), !(numDEInCat<=2 & ontology == "BP"), !(numDEInCat<=2 & ontology == "MF"))
View(GO.05_filtered_FrontloadedD7_D21_mod)







#===================================================================================== # 
# Frontloaded response to MODERATE  DAY7  :::::::::::::::::::::::::::::::::::::::::::::::::::
#===================================================================================== # 

# Biological Function - run GOslim 
goseq_res_BP        <- GO.05_filtered_FrontloadedDay7_moderate %>%  filter(ontology=="BP") # BP - all GO terms upregulated   ONLY LINE TO CHANGE!!!
BP_GOcollection     <- GOCollection(goseq_res_BP$category)
GOslims_BP          <- data.frame(goSlim(BP_GOcollection, slim, "BP")) #Find common parent terms to slim down our list
GOslims_BP$category <- row.names(GOslims_BP) #save rownames as category

# Molecular Function - run GOslim
goseq_res_MF        <- GO.05_filtered_FrontloadedDay7_moderate %>%  filter(ontology=="MF") # BP - all GO terms upregulated   ONLY LINE TO CHANGE!!!
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
gene_names    <- FrontloadedDay7_moderate$Gene.ID # all gene IDs in the particular WGCNA module 


# BIOLOGICAL PROCESS  :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
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

BPslim_GOterm_summary       <- BPslim_final[,c(1,2,5,7)]
s <- strsplit(BPslim_GOterm_summary$GO_list, split = ";")
BPslim_GOterm_summary_2     <- data.frame(slim_term = rep(BPslim_GOterm_summary$slim_term, sapply(s, length)),
                                          Gene_IDs  = rep(BPslim_GOterm_summary$Gene_IDs, sapply(s, length)),
                                          GO_terms = unlist(s))
colnames(goseq_res_BP)[1]   <- 'GO_terms'
BPslim_GOterm_summary_final <- merge(goseq_res_BP, BPslim_GOterm_summary_2, by = 'GO_terms')


s_2 <- strsplit(BPslim_GOterm_summary_final$Gene_IDs, split = ";")
BPslim_GOterm_gene_annotation <- data.frame(slim_term      = rep(BPslim_GOterm_summary_final$slim_term, sapply(s_2, length)),
                                            over_represented_pvalue = rep(BPslim_GOterm_summary_final$over_represented_pvalue, sapply(s_2, length)),
                                            GO_terms       = rep(BPslim_GOterm_summary_final$GO_terms, sapply(s_2, length)),
                                            term           = rep(BPslim_GOterm_summary_final$term, sapply(s_2, length)),
                                            ontology       = rep(BPslim_GOterm_summary_final$ontology, sapply(s_2, length)),
                                            PgenIDs        = unlist(s_2))
BP_master_gene_reference <- merge(BPslim_GOterm_gene_annotation, Pgen_reference, by = "PgenIDs")


# MOLECULAR FUNCTION :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
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

MFslim_GOterm_summary       <- MFslim_final[,c(1,2,5,7)]
s <- strsplit(MFslim_GOterm_summary$GO_list, split = ";")
MFslim_GOterm_summary_2     <- data.frame(slim_term = rep(MFslim_GOterm_summary$slim_term, sapply(s, length)),
                                          Gene_IDs  = rep(MFslim_GOterm_summary$Gene_IDs, sapply(s, length)),
                                          GO_terms = unlist(s))
colnames(goseq_res_MF)[1]   <- 'GO_terms'
MFslim_GOterm_summary_final <- merge(goseq_res_MF, MFslim_GOterm_summary_2, by = 'GO_terms')

s_2 <- strsplit(MFslim_GOterm_summary_final$Gene_IDs, split = ";")
MFslim_GOterm_gene_annotation <- data.frame(slim_term      = rep(MFslim_GOterm_summary_final$slim_term, sapply(s_2, length)),
                                            over_represented_pvalue = rep(MFslim_GOterm_summary_final$over_represented_pvalue, sapply(s_2, length)),
                                            GO_terms       = rep(MFslim_GOterm_summary_final$GO_terms, sapply(s_2, length)),
                                            term           = rep(MFslim_GOterm_summary_final$term, sapply(s_2, length)),
                                            ontology       = rep(MFslim_GOterm_summary_final$ontology, sapply(s_2, length)),
                                            PgenIDs        = unlist(s_2))
MF_master_gene_reference <- merge(MFslim_GOterm_gene_annotation, Pgen_reference, by = "PgenIDs")

write.csv(BPslim_final, file = paste("Analysis/Output/GO/WGCNA_goseq/subseq_treatments_all/Frontloaded_genes/Day7_Frontload_MODERATE_TableGOslim_BiolProc.csv", sep ='')) # save csv file       write.csv(BPslim_final, file = paste("Analysis/Output/GO/WGCNA_goseq/subseq_treatments_all/Frontloaded_genes/Day7_GOterms_and_BOslim_BiolProc.csv", sep ='')) # save csv file       
write.csv(MFslim_final, file = paste("Analysis/Output/GO/WGCNA_goseq/subseq_treatments_all/Frontloaded_genes/Day7_Frontload_MODERATE_TableGOslim_MolFunction.csv", sep ='')) # save csv file       write.csv(BPslim_final, file = paste("Analysis/Output/GO/WGCNA_goseq/subseq_treatments_all/Frontloaded_genes/Day7_GOterms_and_BOslim_BiolProc.csv", sep ='')) # save csv file       


write.csv(BP_master_gene_reference, file = paste("Analysis/Output/GO/WGCNA_goseq/subseq_treatments_all/Frontloaded_genes/Day7_Frontload_MODERATE_GOterms_and_GOslim_BiolProc.csv", sep ='')) # save csv file       write.csv(BPslim_final, file = paste("Analysis/Output/GO/WGCNA_goseq/subseq_treatments_all/Frontloaded_genes/Day7_GOterms_and_BOslim_BiolProc.csv", sep ='')) # save csv file       
write.csv(MF_master_gene_reference, file = paste("Analysis/Output/GO/WGCNA_goseq/subseq_treatments_all/Frontloaded_genes/Day7_Frontload_MODERATE_GOterms_and_GOslim_MolFunction.csv", sep ='')) # save csv file       write.csv(BPslim_final, file = paste("Analysis/Output/GO/WGCNA_goseq/subseq_treatments_all/Frontloaded_genes/Day7_GOterms_and_BOslim_BiolProc.csv", sep ='')) # save csv file       












#===================================================================================== # 
# Frontloaded response to SEVERE  DAY7 :::::::::::::::::::::::::::::::::::::::::::::::::::
#===================================================================================== # 

# Biological Function - run GOslim 
goseq_res_BP        <- GO.05_filtered_FrontloadedDay7_severe %>%  filter(ontology=="BP") # BP - all GO terms upregulated   ONLY LINE TO CHANGE!!!
BP_GOcollection     <- GOCollection(goseq_res_BP$category)
GOslims_BP          <- data.frame(goSlim(BP_GOcollection, slim, "BP")) #Find common parent terms to slim down our list
GOslims_BP$category <- row.names(GOslims_BP) #save rownames as category

# Molecular Function - run GOslim
goseq_res_MF        <- GO.05_filtered_FrontloadedDay7_severe %>%  filter(ontology=="MF") # BP - all GO terms upregulated   ONLY LINE TO CHANGE!!!
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
gene_names    <- FrontloadedDay7_moderate$Gene.ID # all gene IDs in the particular WGCNA module 


# BIOLOGICAL PROCESS  :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
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

BPslim_GOterm_summary       <- BPslim_final[,c(1,2,5,7)]
s <- strsplit(BPslim_GOterm_summary$GO_list, split = ";")
BPslim_GOterm_summary_2     <- data.frame(slim_term = rep(BPslim_GOterm_summary$slim_term, sapply(s, length)),
                                          Gene_IDs  = rep(BPslim_GOterm_summary$Gene_IDs, sapply(s, length)),
                                          GO_terms = unlist(s))
colnames(goseq_res_BP)[1]   <- 'GO_terms'
BPslim_GOterm_summary_final <- merge(goseq_res_BP, BPslim_GOterm_summary_2, by = 'GO_terms')


s_2 <- strsplit(BPslim_GOterm_summary_final$Gene_IDs, split = ";")
BPslim_GOterm_gene_annotation <- data.frame(slim_term      = rep(BPslim_GOterm_summary_final$slim_term, sapply(s_2, length)),
                                            over_represented_pvalue = rep(BPslim_GOterm_summary_final$over_represented_pvalue, sapply(s_2, length)),
                                            GO_terms       = rep(BPslim_GOterm_summary_final$GO_terms, sapply(s_2, length)),
                                            term           = rep(BPslim_GOterm_summary_final$term, sapply(s_2, length)),
                                            ontology       = rep(BPslim_GOterm_summary_final$ontology, sapply(s_2, length)),
                                            PgenIDs        = unlist(s_2))
BP_master_gene_reference <- merge(BPslim_GOterm_gene_annotation, Pgen_reference, by = "PgenIDs")


# MOLECULAR FUNCTION :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
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

MFslim_GOterm_summary       <- MFslim_final[,c(1,2,5,7)]
s <- strsplit(MFslim_GOterm_summary$GO_list, split = ";")
MFslim_GOterm_summary_2     <- data.frame(slim_term = rep(MFslim_GOterm_summary$slim_term, sapply(s, length)),
                                          Gene_IDs  = rep(MFslim_GOterm_summary$Gene_IDs, sapply(s, length)),
                                          GO_terms = unlist(s))
colnames(goseq_res_MF)[1]   <- 'GO_terms'
MFslim_GOterm_summary_final <- merge(goseq_res_MF, MFslim_GOterm_summary_2, by = 'GO_terms')

s_2 <- strsplit(MFslim_GOterm_summary_final$Gene_IDs, split = ";")
MFslim_GOterm_gene_annotation <- data.frame(slim_term      = rep(MFslim_GOterm_summary_final$slim_term, sapply(s_2, length)),
                                            over_represented_pvalue = rep(MFslim_GOterm_summary_final$over_represented_pvalue, sapply(s_2, length)),
                                            GO_terms       = rep(MFslim_GOterm_summary_final$GO_terms, sapply(s_2, length)),
                                            term           = rep(MFslim_GOterm_summary_final$term, sapply(s_2, length)),
                                            ontology       = rep(MFslim_GOterm_summary_final$ontology, sapply(s_2, length)),
                                            PgenIDs        = unlist(s_2))
MF_master_gene_reference <- merge(MFslim_GOterm_gene_annotation, Pgen_reference, by = "PgenIDs")


write.csv(BPslim_final, file = paste("Analysis/Output/GO/WGCNA_goseq/subseq_treatments_all/Frontloaded_genes/Day7_Frontload_SEVERE_TableGOslim_BiolProc.csv", sep ='')) # save csv file       write.csv(BPslim_final, file = paste("Analysis/Output/GO/WGCNA_goseq/subseq_treatments_all/Frontloaded_genes/Day7_GOterms_and_BOslim_BiolProc.csv", sep ='')) # save csv file       
write.csv(MFslim_final, file = paste("Analysis/Output/GO/WGCNA_goseq/subseq_treatments_all/Frontloaded_genes/Day7_Frontload_SEVERE_TableGOslim_MolFunction.csv", sep ='')) # save csv file       write.csv(BPslim_final, file = paste("Analysis/Output/GO/WGCNA_goseq/subseq_treatments_all/Frontloaded_genes/Day7_GOterms_and_BOslim_BiolProc.csv", sep ='')) # save csv file       


write.csv(BP_master_gene_reference, file = paste("Analysis/Output/GO/WGCNA_goseq/subseq_treatments_all/Frontloaded_genes/Day7_Frontload_SEVERE_GOterms_and_GOslim_BiolProc.csv", sep ='')) # save csv file       write.csv(BPslim_final, file = paste("Analysis/Output/GO/WGCNA_goseq/subseq_treatments_all/Frontloaded_genes/Day7_GOterms_and_BOslim_BiolProc.csv", sep ='')) # save csv file       
write.csv(MF_master_gene_reference, file = paste("Analysis/Output/GO/WGCNA_goseq/subseq_treatments_all/Frontloaded_genes/Day7_Frontload_SEVERE_GOterms_and_GOslim_MolFunction.csv", sep ='')) # save csv file       write.csv(BPslim_final, file = paste("Analysis/Output/GO/WGCNA_goseq/subseq_treatments_all/Frontloaded_genes/Day7_GOterms_and_BOslim_BiolProc.csv", sep ='')) # save csv file       













#===================================================================================== # 
# Frontloaded SHARED BETWEEN MODERATE AND SEVERE DAY7 !!!  :::::::::::::::::::::::::::::::::::::::::::::::::::
#===================================================================================== # 

# Biological Function - run GOslim 
goseq_res_BP        <- GO.05_filtered_FrontloadedDay7_shared %>%  filter(ontology=="BP") # BP - all GO terms upregulated   ONLY LINE TO CHANGE!!!
BP_GOcollection     <- GOCollection(goseq_res_BP$category)
GOslims_BP          <- data.frame(goSlim(BP_GOcollection, slim, "BP")) #Find common parent terms to slim down our list
GOslims_BP$category <- row.names(GOslims_BP) #save rownames as category

# Molecular Function - run GOslim
goseq_res_MF        <- GO.05_filtered_FrontloadedDay7_shared %>%  filter(ontology=="MF") # BP - all GO terms upregulated   ONLY LINE TO CHANGE!!!
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
gene_names    <- FrontloadedDay7_moderate$Gene.ID # all gene IDs in the particular WGCNA module 


# BIOLOGICAL PROCESS  :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
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

BPslim_GOterm_summary       <- BPslim_final[,c(1,2,5,7)]
s <- strsplit(BPslim_GOterm_summary$GO_list, split = ";")
BPslim_GOterm_summary_2     <- data.frame(slim_term = rep(BPslim_GOterm_summary$slim_term, sapply(s, length)),
                                          Gene_IDs  = rep(BPslim_GOterm_summary$Gene_IDs, sapply(s, length)),
                                          GO_terms = unlist(s))
colnames(goseq_res_BP)[1]   <- 'GO_terms'
BPslim_GOterm_summary_final <- merge(goseq_res_BP, BPslim_GOterm_summary_2, by = 'GO_terms')


s_2 <- strsplit(BPslim_GOterm_summary_final$Gene_IDs, split = ";")
BPslim_GOterm_gene_annotation <- data.frame(slim_term      = rep(BPslim_GOterm_summary_final$slim_term, sapply(s_2, length)),
                                            over_represented_pvalue = rep(BPslim_GOterm_summary_final$over_represented_pvalue, sapply(s_2, length)),
                                            GO_terms       = rep(BPslim_GOterm_summary_final$GO_terms, sapply(s_2, length)),
                                            term           = rep(BPslim_GOterm_summary_final$term, sapply(s_2, length)),
                                            ontology       = rep(BPslim_GOterm_summary_final$ontology, sapply(s_2, length)),
                                            PgenIDs        = unlist(s_2))
BP_master_gene_reference <- merge(BPslim_GOterm_gene_annotation, Pgen_reference, by = "PgenIDs")


# MOLECULAR FUNCTION :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
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

MFslim_GOterm_summary       <- MFslim_final[,c(1,2,5,7)]
s <- strsplit(MFslim_GOterm_summary$GO_list, split = ";")
MFslim_GOterm_summary_2     <- data.frame(slim_term = rep(MFslim_GOterm_summary$slim_term, sapply(s, length)),
                                          Gene_IDs  = rep(MFslim_GOterm_summary$Gene_IDs, sapply(s, length)),
                                          GO_terms = unlist(s))
colnames(goseq_res_MF)[1]   <- 'GO_terms'
MFslim_GOterm_summary_final <- merge(goseq_res_MF, MFslim_GOterm_summary_2, by = 'GO_terms')

s_2 <- strsplit(MFslim_GOterm_summary_final$Gene_IDs, split = ";")
MFslim_GOterm_gene_annotation <- data.frame(slim_term      = rep(MFslim_GOterm_summary_final$slim_term, sapply(s_2, length)),
                                            over_represented_pvalue = rep(MFslim_GOterm_summary_final$over_represented_pvalue, sapply(s_2, length)),
                                            GO_terms       = rep(MFslim_GOterm_summary_final$GO_terms, sapply(s_2, length)),
                                            term           = rep(MFslim_GOterm_summary_final$term, sapply(s_2, length)),
                                            ontology       = rep(MFslim_GOterm_summary_final$ontology, sapply(s_2, length)),
                                            PgenIDs        = unlist(s_2))
MF_master_gene_reference <- merge(MFslim_GOterm_gene_annotation, Pgen_reference, by = "PgenIDs")


write.csv(BPslim_final, file = paste("Analysis/Output/GO/WGCNA_goseq/subseq_treatments_all/Frontloaded_genes/Day7_Frontload_SHARED_TableGOslim_BiolProc.csv", sep ='')) # save csv file       write.csv(BPslim_final, file = paste("Analysis/Output/GO/WGCNA_goseq/subseq_treatments_all/Frontloaded_genes/Day7_GOterms_and_BOslim_BiolProc.csv", sep ='')) # save csv file       
write.csv(MFslim_final, file = paste("Analysis/Output/GO/WGCNA_goseq/subseq_treatments_all/Frontloaded_genes/Day7_Frontload_SHARED_TableGOslim_MolFunction.csv", sep ='')) # save csv file       write.csv(BPslim_final, file = paste("Analysis/Output/GO/WGCNA_goseq/subseq_treatments_all/Frontloaded_genes/Day7_GOterms_and_BOslim_BiolProc.csv", sep ='')) # save csv file       


write.csv(BP_master_gene_reference, file = paste("Analysis/Output/GO/WGCNA_goseq/subseq_treatments_all/Frontloaded_genes/Day7_Frontload_SHARED_GOterms_and_GOslim_BiolProc.csv", sep ='')) # save csv file       write.csv(BPslim_final, file = paste("Analysis/Output/GO/WGCNA_goseq/subseq_treatments_all/Frontloaded_genes/Day7_GOterms_and_BOslim_BiolProc.csv", sep ='')) # save csv file       
write.csv(MF_master_gene_reference, file = paste("Analysis/Output/GO/WGCNA_goseq/subseq_treatments_all/Frontloaded_genes/Day7_Frontload_SHARED_GOterms_and_GOslim_MolFunction.csv", sep ='')) # save csv file       write.csv(BPslim_final, file = paste("Analysis/Output/GO/WGCNA_goseq/subseq_treatments_all/Frontloaded_genes/Day7_GOterms_and_BOslim_BiolProc.csv", sep ='')) # save csv file     













#===================================================================================== # 
# Frontloaded response to MOERATE Day 21  :::::::::::::::::::::::::::::::::::::::::::::::::::
#===================================================================================== # 

# Biological Function - run GOslim 
goseq_res_BP        <- GO.05_filtered_FrontloadedDay21 %>%  filter(ontology=="BP") # BP - all GO terms upregulated   ONLY LINE TO CHANGE!!!
BP_GOcollection     <- GOCollection(goseq_res_BP$category)
GOslims_BP          <- data.frame(goSlim(BP_GOcollection, slim, "BP")) #Find common parent terms to slim down our list
GOslims_BP$category <- row.names(GOslims_BP) #save rownames as category

# Molecular Function - run GOslim
goseq_res_MF        <- GO.05_filtered_FrontloadedDay21 %>%  filter(ontology=="MF") # BP - all GO terms upregulated   ONLY LINE TO CHANGE!!!
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
gene_names    <- FrontloadedDay21$Gene.ID # all gene IDs in the particular WGCNA module 


# BIOLOGICAL PROCESS  :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
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

BPslim_GOterm_summary       <- BPslim_final[,c(1,2,5,7)]
s <- strsplit(BPslim_GOterm_summary$GO_list, split = ";")
BPslim_GOterm_summary_2     <- data.frame(slim_term = rep(BPslim_GOterm_summary$slim_term, sapply(s, length)),
                                          Gene_IDs  = rep(BPslim_GOterm_summary$Gene_IDs, sapply(s, length)),
                                          GO_terms = unlist(s))
colnames(goseq_res_BP)[1]   <- 'GO_terms'
BPslim_GOterm_summary_final <- merge(goseq_res_BP, BPslim_GOterm_summary_2, by = 'GO_terms')


s_2 <- strsplit(BPslim_GOterm_summary_final$Gene_IDs, split = ";")
BPslim_GOterm_gene_annotation <- data.frame(slim_term      = rep(BPslim_GOterm_summary_final$slim_term, sapply(s_2, length)),
                                            over_represented_pvalue = rep(BPslim_GOterm_summary_final$over_represented_pvalue, sapply(s_2, length)),
                                            GO_terms       = rep(BPslim_GOterm_summary_final$GO_terms, sapply(s_2, length)),
                                            term           = rep(BPslim_GOterm_summary_final$term, sapply(s_2, length)),
                                            ontology       = rep(BPslim_GOterm_summary_final$ontology, sapply(s_2, length)),
                                            PgenIDs        = unlist(s_2))
BP_master_gene_reference <- merge(BPslim_GOterm_gene_annotation, Pgen_reference, by = "PgenIDs")


# MOLECULAR FUNCTION :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
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

MFslim_GOterm_summary       <- MFslim_final[,c(1,2,5,7)]
s <- strsplit(MFslim_GOterm_summary$GO_list, split = ";")
MFslim_GOterm_summary_2     <- data.frame(slim_term = rep(MFslim_GOterm_summary$slim_term, sapply(s, length)),
                                          Gene_IDs  = rep(MFslim_GOterm_summary$Gene_IDs, sapply(s, length)),
                                          GO_terms = unlist(s))
colnames(goseq_res_MF)[1]   <- 'GO_terms'
MFslim_GOterm_summary_final <- merge(goseq_res_MF, MFslim_GOterm_summary_2, by = 'GO_terms')

s_2 <- strsplit(MFslim_GOterm_summary_final$Gene_IDs, split = ";")
MFslim_GOterm_gene_annotation <- data.frame(slim_term      = rep(MFslim_GOterm_summary_final$slim_term, sapply(s_2, length)),
                                            over_represented_pvalue = rep(MFslim_GOterm_summary_final$over_represented_pvalue, sapply(s_2, length)),
                                            GO_terms       = rep(MFslim_GOterm_summary_final$GO_terms, sapply(s_2, length)),
                                            term           = rep(MFslim_GOterm_summary_final$term, sapply(s_2, length)),
                                            ontology       = rep(MFslim_GOterm_summary_final$ontology, sapply(s_2, length)),
                                            PgenIDs        = unlist(s_2))
MF_master_gene_reference <- merge(MFslim_GOterm_gene_annotation, Pgen_reference, by = "PgenIDs")


write.csv(BPslim_final, file = paste("Analysis/Output/GO/WGCNA_goseq/subseq_treatments_all/Frontloaded_genes/Day21_Frontload_MODERATEE_TableGOslim_BiolProc.csv", sep ='')) # save csv file       write.csv(BPslim_final, file = paste("Analysis/Output/GO/WGCNA_goseq/subseq_treatments_all/Frontloaded_genes/Day7_GOterms_and_BOslim_BiolProc.csv", sep ='')) # save csv file       
write.csv(MFslim_final, file = paste("Analysis/Output/GO/WGCNA_goseq/subseq_treatments_all/Frontloaded_genes/Day21_Frontload_MODERATE_TableGOslim_MolFunction.csv", sep ='')) # save csv file       write.csv(BPslim_final, file = paste("Analysis/Output/GO/WGCNA_goseq/subseq_treatments_all/Frontloaded_genes/Day7_GOterms_and_BOslim_BiolProc.csv", sep ='')) # save csv file       


write.csv(BP_master_gene_reference, file = paste("Analysis/Output/GO/WGCNA_goseq/subseq_treatments_all/Frontloaded_genes/Day21_Frontload_MODERATE_GOterms_and_GOslim_BiolProc.csv", sep ='')) # save csv file       write.csv(BPslim_final, file = paste("Analysis/Output/GO/WGCNA_goseq/subseq_treatments_all/Frontloaded_genes/Day7_GOterms_and_BOslim_BiolProc.csv", sep ='')) # save csv file       
write.csv(MF_master_gene_reference, file = paste("Analysis/Output/GO/WGCNA_goseq/subseq_treatments_all/Frontloaded_genes/Day21_Frontload_MODERATE_GOterms_and_GOslim_MolFunction.csv", sep ='')) # save csv file       write.csv(BPslim_final, file = paste("Analysis/Output/GO/WGCNA_goseq/subseq_treatments_all/Frontloaded_genes/Day7_GOterms_and_BOslim_BiolProc.csv", sep ='')) # save csv file       












#===================================================================================== # 
# Frontloaded response to MOERATE Day 21  :::::::::::::::::::::::::::::::::::::::::::::::::::
#===================================================================================== # 

# Biological Function - run GOslim 
goseq_res_BP        <- GO.05_filtered_FrontloadedD7_D21_mod %>%  filter(ontology=="BP") # BP - all GO terms upregulated   ONLY LINE TO CHANGE!!!
BP_GOcollection     <- GOCollection(goseq_res_BP$category)
GOslims_BP          <- data.frame(goSlim(BP_GOcollection, slim, "BP")) #Find common parent terms to slim down our list
GOslims_BP$category <- row.names(GOslims_BP) #save rownames as category

# Molecular Function - run GOslim
goseq_res_MF        <- GO.05_filtered_FrontloadedD7_D21_mod %>%  filter(ontology=="MF") # BP - all GO terms upregulated   ONLY LINE TO CHANGE!!!
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
gene_names    <- d7.d21_frontloaded_moderate$Gene.ID # all gene IDs in the particular WGCNA module 


# BIOLOGICAL PROCESS  :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
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

BPslim_GOterm_summary       <- BPslim_final[,c(1,2,5,7)]
s <- strsplit(BPslim_GOterm_summary$GO_list, split = ";")
BPslim_GOterm_summary_2     <- data.frame(slim_term = rep(BPslim_GOterm_summary$slim_term, sapply(s, length)),
                                          Gene_IDs  = rep(BPslim_GOterm_summary$Gene_IDs, sapply(s, length)),
                                          GO_terms = unlist(s))
colnames(goseq_res_BP)[1]   <- 'GO_terms'
BPslim_GOterm_summary_final <- merge(goseq_res_BP, BPslim_GOterm_summary_2, by = 'GO_terms')


s_2 <- strsplit(BPslim_GOterm_summary_final$Gene_IDs, split = ";")
BPslim_GOterm_gene_annotation <- data.frame(slim_term      = rep(BPslim_GOterm_summary_final$slim_term, sapply(s_2, length)),
                                            over_represented_pvalue = rep(BPslim_GOterm_summary_final$over_represented_pvalue, sapply(s_2, length)),
                                            GO_terms       = rep(BPslim_GOterm_summary_final$GO_terms, sapply(s_2, length)),
                                            term           = rep(BPslim_GOterm_summary_final$term, sapply(s_2, length)),
                                            ontology       = rep(BPslim_GOterm_summary_final$ontology, sapply(s_2, length)),
                                            PgenIDs        = unlist(s_2))
BP_master_gene_reference <- merge(BPslim_GOterm_gene_annotation, Pgen_reference, by = "PgenIDs")
View(BP_master_gene_reference)

# MOLECULAR FUNCTION :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
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

MFslim_GOterm_summary       <- MFslim_final[,c(1,2,5,7)]
s <- strsplit(MFslim_GOterm_summary$GO_list, split = ";")
MFslim_GOterm_summary_2     <- data.frame(slim_term = rep(MFslim_GOterm_summary$slim_term, sapply(s, length)),
                                          Gene_IDs  = rep(MFslim_GOterm_summary$Gene_IDs, sapply(s, length)),
                                          GO_terms = unlist(s))
colnames(goseq_res_MF)[1]   <- 'GO_terms'
MFslim_GOterm_summary_final <- merge(goseq_res_MF, MFslim_GOterm_summary_2, by = 'GO_terms')

s_2 <- strsplit(MFslim_GOterm_summary_final$Gene_IDs, split = ";")
MFslim_GOterm_gene_annotation <- data.frame(slim_term      = rep(MFslim_GOterm_summary_final$slim_term, sapply(s_2, length)),
                                            over_represented_pvalue = rep(MFslim_GOterm_summary_final$over_represented_pvalue, sapply(s_2, length)),
                                            GO_terms       = rep(MFslim_GOterm_summary_final$GO_terms, sapply(s_2, length)),
                                            term           = rep(MFslim_GOterm_summary_final$term, sapply(s_2, length)),
                                            ontology       = rep(MFslim_GOterm_summary_final$ontology, sapply(s_2, length)),
                                            PgenIDs        = unlist(s_2))
MF_master_gene_reference <- merge(MFslim_GOterm_gene_annotation, Pgen_reference, by = "PgenIDs")


write.csv(BPslim_final, file = paste("Analysis/Output/GO/WGCNA_goseq/subseq_treatments_all/Frontloaded_genes/Day7.Day21_Frontload_Shared_MODERATE_TableGOslim_BiolProc.csv", sep ='')) # save csv file       write.csv(BPslim_final, file = paste("Analysis/Output/GO/WGCNA_goseq/subseq_treatments_all/Frontloaded_genes/Day7_GOterms_and_BOslim_BiolProc.csv", sep ='')) # save csv file       

write.csv(BP_master_gene_reference, file = paste("Analysis/Output/GO/WGCNA_goseq/subseq_treatments_all/Frontloaded_genes/Day7.Day21_Frontload_Shared_MODERATEE_GOterms_and_GOslim_BiolProc.csv", sep ='')) # save csv file       write.csv(BPslim_final, file = paste("Analysis/Output/GO/WGCNA_goseq/subseq_treatments_all/Frontloaded_genes/Day7_GOterms_and_BOslim_BiolProc.csv", sep ='')) # save csv file       



