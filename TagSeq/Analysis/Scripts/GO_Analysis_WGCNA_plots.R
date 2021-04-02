
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
library(GO.db)
library(GSEABase)
library(data.table) 
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
d7_Annot_ModuleMembership      <-  read.csv("Analysis/Output/WGCNA/Day7/d7.WGCNA_ModulMembership.csv")   # WGCNA results day 7  - Module membership 
d7_Annot_ModuleMembership      <- d7_Annot_ModuleMembership[,c(3:4,7:8)] # call gene ID, module color, and GO annotation
d7_Annot_ModuleMembership$Day  <- "Day7"  # common column to divide master dataset
d7ModCols                      <- data.frame(moduleColor = unique(d7_Annot_ModuleMembership$moduleColor)) # call all unique module colors 
d7ModCols                      <- d7ModCols %>% filter(moduleColor %in% c('brown', 'yellow', 'green')) # MODULES WITH SIG CORR WITH TREATMENT
d7ModCols$Day                  <- "Day7" # common column for the for loop

d14_Annot_ModuleMembership     <-  read.csv("Analysis/Output/WGCNA/Day14/d14.WGCNA_ModulMembership.csv") # WGCNA results day 14 - Module membership 
d14_Annot_ModuleMembership     <- d14_Annot_ModuleMembership[,c(3:4,7:8)] # call gene ID, module color, and GO annotation
d14_Annot_ModuleMembership$Day <- "Day14"  # common column to divide master dataset
d14ModCols                     <- data.frame(moduleColor = unique(d14_Annot_ModuleMembership$moduleColor)) # call all unique module colors 
d14ModCols                     <- d14ModCols %>% filter(moduleColor %in% c('brown', 'black', 'pink', 'magenta')) # MODULES WITH SIG CORR WITH TREATMENT
d14ModCols$Day                 <- "Day14" # common column for the for loop

d21_Annot_ModuleMembership     <-  read.csv("Analysis/Output/WGCNA/Day21/d21.WGCNA_ModulMembership.csv") # WGCNA results day 21 - Module membership 
d21_Annot_ModuleMembership     <- d21_Annot_ModuleMembership[,c(3:4,7:8)] # call gene ID, module color, and GO annotation
d21_Annot_ModuleMembership$Day <- "Day21" # common column to divide master dataset
d21ModCols                     <- data.frame(moduleColor = unique(d21_Annot_ModuleMembership$moduleColor)) # call all unique module colors 
d21ModCols                     <- d21ModCols %>% filter(moduleColor %in% c('magenta', 'blue', 'yellow', 'red', 'black', 'pink')) # MODULES WITH SIG CORR WITH TREATMENT
d21ModCols$Day                 <- "Day21" # common column for the for loop
 
WGCNA_MasterModData   <-  rbind(d7_Annot_ModuleMembership, d14_Annot_ModuleMembership, d21_Annot_ModuleMembership) # master WGCNA data table 
WGCNA_ColorList       <-  rbind(d7ModCols, d14ModCols, d21ModCols) # master WGCNA color list - use this to loop all the analysis 

slim <- getOBOCollection("http://current.geneontology.org/ontology/subsets/goslim_generic.obo") #get GO database - # call goslim_generic.obo terms as 'slim'

#===================================================================================== # 
# LOAD DATA - Raw Count Data - FILTERED COUNT MATRICES USED IN WGCNA- 10CPM in 50% of samples
#
#===================================================================================== # 
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
IDvector.d7         <- as.vector(unique(Day7_all.counts$gene.ID))  # call unique genes (those filtered and used in DESEq2) on day7 - 'IDvector'
IDvector.d14        <- as.vector(unique(Day14_all.counts$gene.ID)) # call unique genes (those filtered and used in DESEq2) on day14 - 'IDvector'
IDvector.d21        <- as.vector(unique(Day21_all.counts$gene.ID)) # call unique genes (those filtered and used in DESEq2) on day21 - 'IDvector'

# (3) Gene length 
# length vector  
GO_gene.length <- Geoduck_annotation %>% dplyr::mutate(length = V4-V3) %>%  dplyr::select(c("V1","length"))
names(GO_gene.length)[1] <- "gene.ID"
# merge length with counts data
length_vector   <- GO_gene.length$length
GeneLength.d7   <- merge(GO_gene.length, Day7_all.counts, by = "gene.ID")  # merge day7 counts with 'GO_gene.length' 
GeneLength.d14  <- merge(GO_gene.length, Day14_all.counts, by = "gene.ID")  # merge day14 counts with 'GO_gene.length' 
GeneLength.d21  <- merge(GO_gene.length, Day21_all.counts, by = "gene.ID")  # merge day21 counts with 'GO_gene.length'
# call length values for goseq - confirms that the IDvector and length_vector are the same!!!
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
        if (WGCNA_ColorList[i,2] == "Day7") {
        Mod <- d7_Annot_ModuleMembership %>% dplyr::filter(moduleColor %in% WGCNA_ColorList[i,1]) # call the WGCNA Day - essential here!
        Modgenes <- Mod[1]
        names(Modgenes)[1] <- "Gene.ID" # 162 genws in the green module 
        # Mod_integer <- as.integer(IDvector.d7 %in% (Modgenes$Gene.ID)) # call the day-specific ID vector 
        # names(Mod_integer)=IDvector.d7 # rename
        Mod_integer <- as.integer(GO_unique.genes.all %in% (Modgenes$Gene.ID)) # call the day-specific ID vector 
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
        GO.05$Day       <- "Day7"
        
        write.csv(GO.05, file = paste("Analysis/Output/GO/WGCNA_goseq/Day7/GO.05",WGCNA_ColorList[i,1], "Module.csv", sep ='')) # save csv file
        
        
         } else if (WGCNA_ColorList[i,2] == "Day14") {
            Mod <- d14_Annot_ModuleMembership %>% dplyr::filter(moduleColor %in% WGCNA_ColorList[i,1]) # call the WGCNA Day - essential here!
            Modgenes <- Mod[1]
            names(Modgenes)[1] <- "Gene.ID" # 162 genws in the green module 
            # Mod_integer <- as.integer(IDvector.d14 %in% (Modgenes$Gene.ID)) # call the day-specific ID vector 
            # names(Mod_integer)=IDvector.d14 # rename 
            Mod_integer <- as.integer(GO_unique.genes.all %in% (Modgenes$Gene.ID)) # call the day-specific ID vector 
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
            
            write.csv(GO.05, file = paste("Analysis/Output/GO/WGCNA_goseq/Day14/GO.05",WGCNA_ColorList[i,1], "Module.csv", sep ='')) # save csv file
            
            } else {
              Mod <- d21_Annot_ModuleMembership %>% dplyr::filter(moduleColor %in% WGCNA_ColorList[i,1]) # call the WGCNA Day - essential here!
              Modgenes <- Mod[1]
              names(Modgenes)[1] <- "Gene.ID" # 162 genws in the green module 
              # Mod_integer <- as.integer(IDvector.d21 %in% (Modgenes$Gene.ID)) # call the day-specific ID vector 
              # names(Mod_integer)=IDvector.d21 # rename
              Mod_integer <- as.integer(GO_unique.genes.all %in% (Modgenes$Gene.ID)) # call the day-specific ID vector 
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
              
              write.csv(GO.05, file = paste("Analysis/Output/GO/WGCNA_goseq/Day21/GO.05",WGCNA_ColorList[i,1], "Module.csv", sep ='')) # save csv file              
              
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
d7_GO.05brownModule   <- read.csv("Analysis/Output/GO/WGCNA_goseq/Day7/GO.05brownModule.csv")
d7_GO.05greenModule   <- read.csv("Analysis/Output/GO/WGCNA_goseq/Day7/GO.05greenModule.csv")
d7_GO.05yellowModule  <- read.csv("Analysis/Output/GO/WGCNA_goseq/Day7/GO.05yellowModule.csv")

d14_GO.05blackModule   <- read.csv("Analysis/Output/GO/WGCNA_goseq/Day14/GO.05blackModule.csv")
d14_GO.05brownModule   <- read.csv("Analysis/Output/GO/WGCNA_goseq/Day14/GO.05brownModule.csv")
d14_GO.05magentaModule <- read.csv("Analysis/Output/GO/WGCNA_goseq/Day14/GO.05magentaModule.csv")
d14_GO.05pinkModule    <- read.csv("Analysis/Output/GO/WGCNA_goseq/Day14/GO.05pinkModule.csv")

d21_GO.05blackModule    <- read.csv("Analysis/Output/GO/WGCNA_goseq/Day21/GO.05blackModule.csv")
d21_GO.05blueModule     <- read.csv("Analysis/Output/GO/WGCNA_goseq/Day21/GO.05blueModule.csv")
d21_GO.05magentaModule  <- read.csv("Analysis/Output/GO/WGCNA_goseq/Day21/GO.05magentaModule.csv")
d21_GO.05pinkModule     <- read.csv("Analysis/Output/GO/WGCNA_goseq/Day21/GO.05pinkModule.csv")
d21_GO.05redModule      <- read.csv("Analysis/Output/GO/WGCNA_goseq/Day21/GO.05redModule.csv")
d21_GO.05yellowModule   <- read.csv("Analysis/Output/GO/WGCNA_goseq/Day21/GO.05yellowModule.csv")

Master_goseq_results <- rbind(d7_GO.05brownModule, d7_GO.05greenModule, d7_GO.05yellowModule,
                              d14_GO.05blackModule, d14_GO.05brownModule, d14_GO.05magentaModule, d14_GO.05pinkModule,
                              d21_GO.05blackModule, d21_GO.05blueModule, d21_GO.05magentaModule, d21_GO.05pinkModule, d21_GO.05redModule, d21_GO.05yellowModule)

GOslimLoop_vars <- unique(Master_goseq_results[c(9,10)])



for (i in 1:nrow(GOslimLoop_vars)) {
  # call the target dataset
  goseq_res       <- Master_goseq_results %>%  dplyr::filter(Day %in% GOslimLoop_vars[i,2], moduleColor %in% GOslimLoop_vars[i,1])
  WGCNA_res       <- WGCNA_MasterModData  %>%  dplyr::filter(Day %in% GOslimLoop_vars[i,2], moduleColor %in% GOslimLoop_vars[i,1])
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
  BPslim             <- filter(BPslim_Mapped, Count>5 & Term!="biological_process") #filter out empty slims and term "biological process"
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
  BPslim_final <- data.frame(slim_term=BPslim_C$Term, slim_cat=BPslim_C$category, category=BPslim_C$go_term, Gene.Count=BPslim_C$Gene.Count, GO.Count=BPslim_C$Count) #rename columns) #rename columns
  BPslim_final$module_day <- paste(GOslimLoop_vars[i,2], GOslimLoop_vars[i,1], sep = '_')
  
  # MOLECULAR FUNCTION
  MFslim             <- filter(MFslim_Mapped, Count>2 & Term!="molecular_function") #filter out empty slims and term "biological process"
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
    MFslim_final <- data.frame(slim_term=MFslim_C$Term, slim_cat=MFslim_C$category, category=MFslim_C$go_term, Gene.Count=MFslim_C$Gene.Count, GO.Count=MFslim_C$Count) #rename columns) #rename columns
    MFslim_final$module_day <- paste(GOslimLoop_vars[i,2], GOslimLoop_vars[i,1], sep = '_')
    
    # save GOslim final datasets for BP and MF of each  module in their respective folder(s) by Day
    write.csv(MFslim_final, file = paste("Analysis/Output/GO/WGCNA_goseq/", GOslimLoop_vars[i,2], "/GOslim_MolFunction_",GOslimLoop_vars[i,1], "Module.csv", sep ='')) # save csv file              
    write.csv(BPslim_final, file = paste("Analysis/Output/GO/WGCNA_goseq/", GOslimLoop_vars[i,2], "/GOslim_BiolProc_",GOslimLoop_vars[i,1], "Module.csv", sep ='')) # save csv file       

} # end of GOslim for loop! (i in 1:nrow)


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
d7_slimBP_brownModule   <- read.csv("Analysis/Output/GO/WGCNA_goseq/Day7/GOslim_BiolProc_brownModule.csv")
d7_slimBP_greenModule   <- read.csv("Analysis/Output/GO/WGCNA_goseq/Day7/GOslim_BiolProc_greenModule.csv")
d7_slimBP_yellowModule  <- read.csv("Analysis/Output/GO/WGCNA_goseq/Day7/GOslim_BiolProc_yellowModule.csv")

d14_slimBP_blackModule   <- read.csv("Analysis/Output/GO/WGCNA_goseq/Day14/GOslim_BiolProc_blackModule.csv")
d14_slimBP_brownModule   <- read.csv("Analysis/Output/GO/WGCNA_goseq/Day14/GOslim_BiolProc_brownModule.csv")
d14_slimBP_magentaModule <- read.csv("Analysis/Output/GO/WGCNA_goseq/Day14/GOslim_BiolProc_magentaModule.csv")
d14_slimBP_pinkModule    <- read.csv("Analysis/Output/GO/WGCNA_goseq/Day14/GOslim_BiolProc_pinkModule.csv")

d21_slimBP_blackModule    <- read.csv("Analysis/Output/GO/WGCNA_goseq/Day21/GOslim_BiolProc_blackModule.csv")
d21_slimBP_blueModule     <- read.csv("Analysis/Output/GO/WGCNA_goseq/Day21/GOslim_BiolProc_blueModule.csv")
d21_slimBP_magentaModule  <- read.csv("Analysis/Output/GO/WGCNA_goseq/Day21/GOslim_BiolProc_magentaModule.csv")
d21_slimBP_pinkModule     <- read.csv("Analysis/Output/GO/WGCNA_goseq/Day21/GOslim_BiolProc_pinkModule.csv")
d21_slimBP_redModule      <- read.csv("Analysis/Output/GO/WGCNA_goseq/Day21/GOslim_BiolProc_redModule.csv")
d21_slimBP_yellowModule   <- read.csv("Analysis/Output/GO/WGCNA_goseq/Day21/GOslim_BiolProc_yellowModule.csv")
# MOLECULAR FUNCTION GO  SLIMS
d7_slimMF_brownModule   <- read.csv("Analysis/Output/GO/WGCNA_goseq/Day7/GOslim_MolFunction_brownModule.csv")
d7_slimMF_greenModule   <- read.csv("Analysis/Output/GO/WGCNA_goseq/Day7/GOslim_MolFunction_greenModule.csv")
d7_slimMF_yellowModule  <- read.csv("Analysis/Output/GO/WGCNA_goseq/Day7/GOslim_MolFunction_yellowModule.csv")

d14_slimMF_blackModule   <- read.csv("Analysis/Output/GO/WGCNA_goseq/Day14/GOslim_MolFunction_blackModule.csv")
d14_slimMF_brownModule   <- read.csv("Analysis/Output/GO/WGCNA_goseq/Day14/GOslim_MolFunction_brownModule.csv")
d14_slimMF_magentaModule <- read.csv("Analysis/Output/GO/WGCNA_goseq/Day14/GOslim_MolFunction_magentaModule.csv")
d14_slimMF_pinkModule    <- read.csv("Analysis/Output/GO/WGCNA_goseq/Day14/GOslim_MolFunction_pinkModule.csv")

d21_slimMF_blackModule    <- read.csv("Analysis/Output/GO/WGCNA_goseq/Day21/GOslim_MolFunction_blackModule.csv")
d21_slimMF_blueModule     <- read.csv("Analysis/Output/GO/WGCNA_goseq/Day21/GOslim_MolFunction_blueModule.csv")
d21_slimMF_magentaModule  <- read.csv("Analysis/Output/GO/WGCNA_goseq/Day21/GOslim_MolFunction_magentaModule.csv")
d21_slimMF_pinkModule     <- read.csv("Analysis/Output/GO/WGCNA_goseq/Day21/GOslim_MolFunction_pinkModule.csv")
d21_slimMF_redModule      <- read.csv("Analysis/Output/GO/WGCNA_goseq/Day21/GOslim_MolFunction_redModule.csv")
d21_slimMF_yellowModule   <- read.csv("Analysis/Output/GO/WGCNA_goseq/Day21/GOslim_MolFunction_yellowModule.csv")


# BIND DATA TO FACET THE PLOTS    ---------------------------------------------------------- #
# (1) By sampling date 0 Day 7 14 and 21 separately
BP_D7  <- rbind(d7_slimBP_brownModule, d7_slimBP_greenModule, d7_slimBP_yellowModule)
BP_D14 <- rbind(d14_slimBP_blackModule, d14_slimBP_brownModule, d14_slimBP_magentaModule, d14_slimBP_pinkModule)
BP_D21 <- rbind(d21_slimBP_blackModule, d21_slimBP_blueModule, d21_slimBP_magentaModule, d21_slimBP_pinkModule, d21_slimBP_redModule, d21_slimBP_yellowModule)

MF_D7  <- rbind(d7_slimMF_brownModule, d7_slimMF_greenModule, d7_slimMF_yellowModule)
MF_D14 <- rbind(d14_slimMF_blackModule, d14_slimMF_brownModule, d14_slimMF_magentaModule, d14_slimMF_pinkModule)
MF_D21 <- rbind(d21_slimMF_blackModule, d21_slimMF_blueModule, d21_slimMF_magentaModule, d21_slimMF_pinkModule, d21_slimMF_redModule, d21_slimMF_yellowModule)



# (2) By trend in vst expression pattern ACROSS all days 
# (2.A) Priming Effect AMBIENT > MODERATE 
# (2.B) Priming Effect MODERATE > AMBIENT 


# PLOTTING --------------------------------------------------------------------------------------------------- #
# (1) By sampling date 0 Day 7 14 and 21 separately



#  DAY 7 



BP_D7$Ont <- "BP" # create a new common clumn to call in the plot 
BP_D7$module_day <- factor(BP_D7$module_day, levels = c("Day7_brown", "Day7_yellow", "Day7_green"))# reorder the facotr level for the facet wrap plot 
BP_D7_filtered <- BP_D7 %>%  dplyr::filter(Gene.Count > 1) # ommit all with gene counts <1
BP_D7_filtered$slim_term <- factor(BP_D7_filtered$slim_term ,levels=rev(unique(BP_D7_filtered$slim_term))) # make slim term alphabetical for plotting

BP_D7_Plot <-ggplot(data = BP_D7_filtered, aes(x = Ont, y = slim_term)) + 
              geom_tile(aes(fill=Gene.Count, width = 1)) + 
              scale_fill_gradient(low = "thistle1", high = "steelblue4") +
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
MF_D7_filtered <- MF_D7 %>%  dplyr::filter(Gene.Count > 1) # ommit all with gene counts <1
MF_D7_filtered$slim_term <- factor(MF_D7_filtered$slim_term ,levels=rev(unique(MF_D7_filtered$slim_term))) # make slim term alphabetical for plotting

MF_D7_Plot <-ggplot(data = MF_D7_filtered, aes(x = Ont, y = slim_term)) + 
              geom_tile(aes(fill=Gene.Count, width = 1)) + 
              scale_fill_gradient(low = "azure2", high = "springgreen4") +
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
BP_D14_filtered <- BP_D14 %>%  dplyr::filter(Gene.Count > 1) # ommit all with gene counts <1
BP_D14_filtered$slim_term <- factor(BP_D14_filtered$slim_term ,levels=rev(unique(BP_D14_filtered$slim_term))) # make slim term alphabetical for plotting

BP_D14_Plot <-ggplot(data = BP_D14_filtered, aes(x = Ont, y = slim_term)) + 
                geom_tile(aes(fill=Gene.Count, width = 1)) + 
                scale_fill_gradient(low = "thistle1", high = "steelblue4") +
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
MF_D14_filtered <- MF_D14 %>%  dplyr::filter(Gene.Count > 1) # ommit all with gene counts <1
MF_D14_filtered$slim_term <- factor(MF_D14_filtered$slim_term ,levels=rev(unique(MF_D14_filtered$slim_term))) # make slim term alphabetical for plotting

MF_D14_Plot <-ggplot(data = MF_D14_filtered, aes(x = Ont, y = slim_term)) + 
              geom_tile(aes(fill=Gene.Count, width = 1)) + 
              scale_fill_gradient(low = "azure2", high = "springgreen4") +
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
BP_D21$module_day <- factor(BP_D21$module_day, levels = c("Day21_magenta",  "Day21_blue", "Day21_yellow",  "Day21_red", "Day21_black", "Day21_pink"))# reorder the facotr level for the facet wrap plot 
BP_D21_filtered <- BP_D21 %>%  dplyr::filter(Gene.Count > 1) # ommit all with gene counts <1
BP_D21_filtered$slim_term <- factor(BP_D21_filtered$slim_term ,levels=rev(unique(BP_D21_filtered$slim_term))) # make slim term alphabetical for plotting

BP_D21_Plot <-ggplot(data = BP_D21_filtered, aes(x = Ont, y = slim_term)) + 
              geom_tile(aes(fill=Gene.Count, width = 1)) + 
              scale_fill_gradient(low = "thistle1", high = "steelblue4") +
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
MF_D21$module_day <- factor(MF_D21$module_day, levels = c("Day21_magenta",  "Day21_blue", "Day21_yellow",  "Day21_red", "Day21_black", "Day21_pink"))# reorder the facotr level for the facet wrap plot 
MF_D21_filtered <- MF_D21 %>%  dplyr::filter(Gene.Count > 1) # ommit all with gene counts <1
MF_D21_filtered$slim_term <- factor(MF_D21_filtered$slim_term ,levels=rev(unique(MF_D21_filtered$slim_term))) # make slim term alphabetical for plotting

MF_D21_Plot <-ggplot(data = MF_D21_filtered, aes(x = Ont, y = slim_term)) + 
              geom_tile(aes(fill=Gene.Count, width = 1)) + 
              scale_fill_gradient(low = "azure2", high = "springgreen4") +
              facet_grid(~module_day, labeller = label_wrap_gen(width = 10, multi_line = TRUE))+
              theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                                 panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                                 strip.text.y = element_text(angle=0, size = 8, face = "bold"),
                                 strip.text.x = element_text(size = 8, face = "bold"),
                                 axis.title.x = element_blank(),
                                 axis.title.y = element_text(size=8),
                                 axis.text = element_text(size = 8), legend.position = "right",
                                 plot.margin = unit(c(0,1,0,0.25), "cm"))+
              ggtitle('GOslim Molecular Function: WGCNA Day 21')

# SAVE PLOTS  ------------------------------------------------------------------------------------------------- #
# WGCNA Day 7 - all significant modules (correlated with treatment(s))
pdf(paste("Analysis//Output/GO/WGCNA_goseq/Day7/GOslim_day7.pdf", sep =''), width=18, height=8)
print(ggarrange(BP_D7_Plot, MF_D7_Plot,         
                 plotlist = NULL,
                 ncol = 2,
                 nrow = 1,
                 labels = NULL))
dev.off()     
# WGCNA Day 14 - all significant modules (correlated with treatment(s))
pdf(paste("Analysis//Output/GO/WGCNA_goseq/Day14/GOslim_day14.pdf", sep =''), width=18, height=8)
print(ggarrange(BP_D14_Plot, MF_D14_Plot,         
                plotlist = NULL,
                ncol = 2,
                nrow = 1,
                labels = NULL))
dev.off()     
# WGCNA Day 21 - all significant modules (correlated with treatment(s))
pdf(paste("Analysis//Output/GO/WGCNA_goseq/Day21/GOslim_day21.pdf", sep =''), width=20, height=8)
print(ggarrange(BP_D21_Plot, MF_D21_Plot,         
                plotlist = NULL,
                ncol = 2,
                nrow = 1,
                labels = NULL))
dev.off() 












# Day 7 sig WGCNA modules: brown, yellow, green
d7_Mod.Brown <- d7_Annot_ModuleMembership %>% dplyr::filter(moduleColor %in% 'brown') # %>%  dplyr::select("geneSymbol")


d7_Mod.Brown_genes <- d7_Mod.Brown[1]

d7_Mod.Yellow <- d7_Annot_ModuleMembership %>% dplyr::filter(moduleColor %in% 'yellow') # %>%  dplyr::select("geneSymbol")
d7_Mod.Yellow_genes <- d7_Mod.Yellow[1]

d7_Mod.Green <- d7_Annot_ModuleMembership %>% dplyr::filter(moduleColor %in% 'green') # %>%  dplyr::select("geneSymbol")
d7_Mod.Green_genes <- d7_Mod.Green[1]

# Day 14 sig WGCNA modules: brown, black, pink, magenta
d14_Mod.Brown   <- d14_Annot_ModuleMembership %>% dplyr::filter(moduleColor %in% 'brown') # %>%  dplyr::select("geneSymbol")
d14_Mod.Brown_genes   <- d14_Mod.Brown[1]

d14_Mod.Black   <- d14_Annot_ModuleMembership %>% dplyr::filter(moduleColor %in% 'black') # %>%  dplyr::select("geneSymbol")
d14_Mod.Black_genes   <- d14_Mod.Black[1]

d14_Mod.Pink    <- d14_Annot_ModuleMembership %>% dplyr::filter(moduleColor %in% 'pink') # %>%  dplyr::select("geneSymbol")
d14_Mod.Pink_genes    <- d14_Mod.Pink[1]

d14_Mod.Magenta <- d14_Annot_ModuleMembership %>% dplyr::filter(moduleColor %in% 'magenta') # %>%  dplyr::select("geneSymbol")
d14_Mod.Magenta_genes <- d14_Mod.Magenta[1]

# Day 21 sig WGCNA modules: magenta, blue, yellow, red, black, pink
d21_Mod.Magenta   <- d21_Annot_ModuleMembership %>% dplyr::filter(moduleColor %in% 'magenta') # %>%  dplyr::select("geneSymbol")
d21_Mod.Magenta_genes   <- d21_Mod.Magenta[1]

d21_Mod.Blue      <- d21_Annot_ModuleMembership %>% dplyr::filter(moduleColor %in% 'blue') # %>%  dplyr::select("geneSymbol")
d21_Mod.Blue_genes      <- d21_Mod.Blue[1]

d21_Mod.Yellow    <- d21_Annot_ModuleMembership %>% dplyr::filter(moduleColor %in% 'yellow') # %>%  dplyr::select("geneSymbol")
d21_Mod.Yellow_genes    <- d21_Mod.Yellow[1]

d21_Mod.Red       <- d21_Annot_ModuleMembership %>% dplyr::filter(moduleColor %in% 'red') # %>%  dplyr::select("geneSymbol")
d21_Mod.Red_genes       <- d21_Mod.Red[1]

d21_Mod.Black     <- d21_Annot_ModuleMembership %>% dplyr::filter(moduleColor %in% 'black') # %>%  dplyr::select("geneSymbol")
d21_Mod.Black_genes     <- d21_Mod.Black[1]

d21_Mod.Pink      <- d21_Annot_ModuleMembership %>% dplyr::filter(moduleColor %in% 'pink') # %>%  dplyr::select("geneSymbol")
d21_Mod.Pink_genes     <- d21_Mod.Pink[1]

# filtered raw count matrices used in WGCNA - 10CPM in 50% of samples
# day7 filtered 10cpm in 50% samples ----------------------------- # 
Day7_all.counts <- read.csv(file="C:/Users/samjg/Documents/My_Projects/Pgenerosa_TagSeq_Metabolomics/TagSeq/Analysis/Data/Filtered_Counts/10cpm_50perc/day7.counts.filtered_10cpm50perc.csv", sep=',', header=TRUE) 
colnames(Day7_all.counts)[1] <- "gene.ID"# rename Pgen gene ID column

# day14 filtered 10cpm in 50% samples ----------------------------- # 
Day14_all.counts <- read.csv(file="C:/Users/samjg/Documents/My_Projects/Pgenerosa_TagSeq_Metabolomics/TagSeq/Analysis/Data/Filtered_Counts/10cpm_50perc/day14.counts.filtered_10cpm50perc.csv", sep=',', header=TRUE) 
colnames(Day14_all.counts)[1] <- "gene.ID"# rename Pgen gene ID column

# day21 filtered 10cpm in 50% samples ----------------------------- # 
Day21_all.counts <- read.csv(file="C:/Users/samjg/Documents/My_Projects/Pgenerosa_TagSeq_Metabolomics/TagSeq/Analysis/Data/Filtered_Counts/10cpm_50perc/day21.counts.filtered_10cpm50perc.csv", sep=',', header=TRUE) 
colnames(Day21_all.counts)[1] <- "gene.ID"# rename Pgen gene ID column




#==============================================================================
#
#  PLOTTING GO TERMS  FOR.... PRIMARY TREATMENT EFFECT AMBIENT > MODERATE 
# 
# day 7  : module BROWN
# day 14 : module BROWN
# day 21 : module MAGENTA
# day 21 : module BLUE 
#==============================================================================
#======================================================================= #
# MODULES WITH PRIMARY EFFECT == AMBIENT > MODERATE (vst expression and positive correlation - view heatmap)
# Day 7 brown module 
names(d7_Mod.Brown_genes)[1] <- "Gene.ID" # 162 genws in the green module 
d7_Mod.Brown_integer <- as.integer(IDvector.d7%in%(d7_Mod.Brown_genes$Gene.ID)) # convert to integer with all unique genes
names(d7_Mod.Brown_integer)=IDvector.d7 # rename

# Day 14 brown module 
names(d14_Mod.Brown_genes)[1] <- "Gene.ID" # 162 genws in the green module 
d14_Mod.Brown_integer <- as.integer(IDvector.d14%in%(d14_Mod.Brown_genes$Gene.ID)) # convert to integer with all unique genes
names(d14_Mod.Brown_integer)=IDvector.d14 # rename

# Day 21 Magenta module 
d21_Mod.Magenta <- d21_Annot_ModuleMembership %>% dplyr::filter(moduleColor %in% 'magenta') # %>%  dplyr::select("geneSymbol")
d21_Mod.Magenta_genes <- d21_Mod.Magenta[1]
names(d21_Mod.Magenta_genes)[1] <- "Gene.ID" # 162 genws in the green module 
d21_Mod.Magenta_integer <- as.integer(IDvector.d21%in%(d21_Mod.Magenta_genes$Gene.ID)) # convert to integer with all unique genes
names(d21_Mod.Magenta_integer)=IDvector.d21 # rename

# Day 21 blue module 
d21_Mod.blue <- d21_Annot_ModuleMembership %>% dplyr::filter(moduleColor %in% 'blue') # %>%  dplyr::select("geneSymbol")
d21_Mod.blue_genes <- d21_Mod.blue[1]
names(d21_Mod.blue_genes)[1] <- "Gene.ID" # 162 genws in the green module 
d21_Mod.blue_integer <- as.integer(IDvector.d21%in%(d21_Mod.blue_genes$Gene.ID)) # convert to integer with all unique genes
names(d21_Mod.blue_integer)=IDvector.d21 # rename

#======================================================================= #
#Calculate Probability Weighting Function (using 'nullp')
d7_Mod.Brown_pwf    <-nullp(d7_Mod.Brown_integer,    id=IDvector.d7, bias.data=length_vector.d7) #weight vector by length of gene - Day 7 module 
d14_Mod.Brown_pwf   <-nullp(d14_Mod.Brown_integer,   id=IDvector.d14, bias.data=length_vector.d14) #weight vector by length of gene - Day 14 module 
d21_Mod.Magenta_pwf <-nullp(d21_Mod.Magenta_integer, id=IDvector.d21, bias.data=length_vector.d21) #weight vector by length of gene - Day 21 module 
d21_Mod.Blue_pwf    <-nullp(d21_Mod.blue_integer,    id=IDvector.d21, bias.data=length_vector.d21) #weight vector by length of gene - Day 21 module 

#======================================================================= #
# Run goseq
d7_Mod.Brown.goseq     <-goseq(d7_Mod.Brown_pwf, gene2cat=GO.terms, test.cats=c("GO:CC", "GO:BP", "GO:MF"), method="Wallenius", use_genes_without_cat=TRUE)
d14_Mod.Brown.goseq    <-goseq(d14_Mod.Brown_pwf,gene2cat=GO.terms, test.cats=c("GO:CC", "GO:BP", "GO:MF"), method="Wallenius", use_genes_without_cat=TRUE)
d21_Mod.Magenta.goseq  <-goseq(d21_Mod.Magenta_pwf, ID.vector, gene2cat=GO.terms, test.cats=c("GO:CC", "GO:BP", "GO:MF"), method="Wallenius", use_genes_without_cat=TRUE)
d21_Mod.Blue.goseq     <-goseq(d21_Mod.Blue_pwf, ID.vector, gene2cat=GO.terms, test.cats=c("GO:CC", "GO:BP", "GO:MF"), method="Wallenius", use_genes_without_cat=TRUE)

#======================================================================= #
# call enriched GO terms and plot 
# day 7 Module Brown
d7_Mod.Brown.GO.05.a<-d7_Mod.Brown.goseq$category[d7_Mod.Brown.goseq$over_represented_pvalue<.05] # change twice here
d7_Mod.Brown.GO.05<-data.frame(d7_Mod.Brown.GO.05.a)
colnames(d7_Mod.Brown.GO.05) <- c("category")
d7_Mod.Brown.GO.05 <- merge(d7_Mod.Brown.GO.05, d7_Mod.Brown.goseq, by="category") # change here
d7_Mod.Brown.GO.05 <- d7_Mod.Brown.GO.05[order(d7_Mod.Brown.GO.05$ontology, d7_Mod.Brown.GO.05$over_represented_pvalue,-d7_Mod.Brown.GO.05$numDEInCat),]
d7_Mod.Brown.GO.05$term <- as.factor(d7_Mod.Brown.GO.05$term)
head(d7_Mod.Brown.GO.05)

# divide into MF and BP datasets
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
#=======================================================================
# Save all the significant go terms in these modules - AMBIENT > MODERATE 
#=======================================================================
write.csv(d7_Mod.Brown.GO.05, file     = "Analysis/Output/GO/WGCNA_goseq/Day7/GO.05.BrownMod.csv", row.names = FALSE)# save significant terms
write.csv(d14_Mod.Brown.GO.05, file    = "Analysis/Output/GO/WGCNA_goseq/Day14/GO.05.BrownMod.csv", row.names = FALSE)# save significant terms
write.csv(d21_Mod.Magenta.GO.05, file  = "Analysis/Output/GO/WGCNA_goseq/Day21/GO.05.MagentaMod.csv", row.names = FALSE)# save significant terms
write.csv(d21_Mod.Blue.GO.05, file     = "Analysis/Output/GO/WGCNA_goseq/Day21/GO.05.BlueMod.csv", row.names = FALSE)# save significant terms

#======================================================================= #
# bind the rows of MF and BP in the clustered modules (those twith primary treatment effect in same pattern/directionality)
# Primary effect: WGCNA significant module showing Ambient > Moderate - bind together day 7 brown, day 14 brown 
# MF
MF_Amb_effect <- rbind(d7_Mod.Brown_MF, d14_Mod.Brown_MF, d21_Mod.Magenta_MF, d21_Mod.Blue_MF) # d21_Mod.Magenta_MF,
# BP
BP_Amb_effect <- rbind(d7_Mod.Brown_BP, d14_Mod.Brown_BP, d21_Mod.Magenta_BP, d21_Mod.Blue_BP)

#=======================================================================
# PLOT and save plots - AMBIENT > MODERATE 
# (heatmap of GO significance by Goterm & sampling day)
#=======================================================================
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
pdf(paste("Analysis/Output/WGCNA/GO_SigEnrich_PrimaryAmbient.pdf", sep =''), width=30, height=60)
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
d7_Mod.Yellow_pwf  <-nullp(d7_Mod.Yellow_integer, GO_unique.genes.all, bias.data=length_vector) #weight vector by length of gene
d14_Mod.Black_pwf  <-nullp(d14_Mod.Black_integer, GO_unique.genes.all, bias.data=length_vector) #weight vector by length of gene
d21_Mod.Yellow_pwf <-nullp(d21_Mod.Yellow_integer, GO_unique.genes.all, bias.data=length_vector) #weight vector by length of gene

#======================================================================= #
# Run goseq
d7_Mod.Yellow.goseq   <-goseq(d7_Mod.Yellow_pwf, ID.vector, gene2cat=GO.terms, test.cats=c("GO:CC", "GO:BP", "GO:MF"), method="Wallenius", use_genes_without_cat=TRUE)
d14_Mod.Black.goseq   <-goseq(d14_Mod.Black_pwf, ID.vector, gene2cat=GO.terms, test.cats=c("GO:CC", "GO:BP", "GO:MF"), method="Wallenius", use_genes_without_cat=TRUE)
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

#=======================================================================
# Save all the significant go terms in these modules -  MODERATE > AMBIENT 
#=======================================================================
write.csv(d7_Mod.Yellow.GO.05, file    = "Analysis/Output/GO/goseq/Day7/GO.05.YellowMod.csv", row.names = FALSE)# save significant terms
write.csv(d14_Mod.Black.GO.05, file    = "Analysis/Output/GO/goseq/Day14/GO.05.BlackMod.csv", row.names = FALSE)# save significant terms
write.csv(d21_Mod.Yellow.GO.05, file   = "Analysis/Output/GO/goseq/Day21/GO.05.YellowMod.csv", row.names = FALSE)# save significant terms

#======================================================================= #
# bind the rows of MF and BP in the clustered modules (those twith primary treatment effect in same pattern/directionality)
# Primary effect: WGCNA significant module showing Ambient > Moderate - bind together day 7 brown, day 14 brown 
# MF
MF_Mod_effect <- rbind(d7_Mod.Yellow_MF, d14_Mod.Black_MF, d21_Mod.Yellow_MF) # d21_Mod.Magenta_MF,
# BP
BP_Mod_effect <- rbind(d7_Mod.Yellow_BP, d14_Mod.Black_BP, d21_Mod.Yellow_BP)


#=======================================================================
# PLOT and save plots - MODERATE > AMBIENT 
# (heatmap of GO significance by Goterm & sampling day)
#=======================================================================
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
pdf(paste("Analysis/Output/WGCNA/GO_SigEnrich_PrimaryModerate.pdf", sep =''), width=30, height=60)
print(ggarrange(MF_tile_Mod_effect, BP_tile_Mod_effect,         
                plotlist = NULL,
                ncol = 2,
                nrow = 1,
                labels = NULL))
dev.off() 


#=======================================================================
#
#
# GOslim analysis 
# - data reduction from the many GO terms to more broad functions/processes
#
#=======================================================================
# Read in sig GO terms from modules above; review 'Save all the significant go terms in these modules'
# Modules suggesing Exp   Ambient  > Modearte  (Review GO heatmaps in the ouput folder)
D7BrownAmb       <- read.csv("Analysis/Output/GO/WGCNA_goseq/Day7/GO.05.BrownMod.csv")
D14BrownAmb      <- read.csv("Analysis/Output/GO/WGCNA_goseq/Day14/GO.05.BrownMod.csv")
D21MagentaAmb    <- read.csv("Analysis/Output/GO/WGCNA_goseq/Day21/GO.05.MagentaMod.csv")
D21BlueAmb       <- read.csv("Analysis/Output/GO/WGCNA_goseq/Day21/GO.05.BlueMod.csv")

# Modules suggesing Exp   Moderate  > Ambient  (Review GO heatmaps in the ouput folder)
D7YellowMod    <- read.csv("Analysis/Output/GO/WGCNA_goseq/Day7/GO.05.YellowMod.csv")
D14BlackMod    <- read.csv("Analysis/Output/GO/WGCNA_goseq/Day14/GO.05.BlackMod.csv")
D21YellowMod   <- read.csv("Analysis/Output/GO/WGCNA_goseq/Day21/GO.05.YellowMod.csv")

# Load GSEABase to call remote GOslim and conduct GO slim analysis 
library(GSEABase)
# call goslim_generic.obo terms as 'slim'
slim <- getOBOCollection("http://current.geneontology.org/ontology/subsets/goslim_generic.obo") #get GO database


# ====================================================================================
# Prep module GO terms for analysis - Separate MF from BP 
#
# ====================================================================================
# ALL MODULES AMBIENT > MODERATE

# D7BrownAmb
D7BrownAmb_BP <- D7BrownAmb %>% # BP - all GO terms upregulated
  filter(ontology=="BP")
D7BrownAmb_BP_GO_collection <- GOCollection(D7BrownAmb_BP$category) #Make library of query terms
D7BrownAmb_GOslims_BP <- data.frame(goSlim(D7BrownAmb_BP_GO_collection, slim, "BP")) #Find common parent terms to slim down our list
D7BrownAmb_GOslims_BP$category <- row.names(D7BrownAmb_GOslims_BP) #save rownames as category

D7BrownAmb_MF <- D7BrownAmb %>% # MF - all GO terms upregulated
  filter(ontology=="MF")
D7BrownAmb_MF_GO_collection <- GOCollection(D7BrownAmb_MF$category) #Make library of query terms
D7BrownAmb_GOslims_MF <- data.frame(goSlim(D7BrownAmb_MF_GO_collection, slim, "MF")) #Find common parent terms to slim down our list
D7BrownAmb_GOslims_MF$category <- row.names(D7BrownAmb_GOslims_MF) #save rownames as category

# D14BrownAmb
D14BrownAmb_BP <- D14BrownAmb %>% # BP - all GO terms upregulated
  filter(ontology=="BP")
D14BrownAmb_BP_GO_collection <- GOCollection(D14BrownAmb_BP$category) #Make library of query terms
D14BrownAmb_GOslims_BP <- data.frame(goSlim(D14BrownAmb_BP_GO_collection, slim, "BP")) #Find common parent terms to slim down our list
D14BrownAmb_GOslims_BP$category <- row.names(D14BrownAmb_GOslims_BP) #save rownames as category

D14BrownAmb_MF <- D14BrownAmb %>% # MF - all GO terms upregulated
  filter(ontology=="MF")
D14BrownAmb_MF_GO_collection <- GOCollection(D14BrownAmb_MF$category) #Make library of query terms
D14BrownAmb_GOslims_MF <- data.frame(goSlim(D14BrownAmb_MF_GO_collection, slim, "MF")) #Find common parent terms to slim down our list
D14BrownAmb_GOslims_MF$category <- row.names(D14BrownAmb_GOslims_MF) #save rownames as category


# D21MagentaAmb
D21MagentaAmb_BP <- D21MagentaAmb %>% # BP - all GO terms upregulated
  filter(ontology=="BP")
D21MagentaAmb_BP_GO_collection <- GOCollection(D21MagentaAmb_BP$category) #Make library of query terms
D21MagentaAmb_GOslims_BP <- data.frame(goSlim(D21MagentaAmb_BP_GO_collection, slim, "BP")) #Find common parent terms to slim down our list
D21MagentaAmb_GOslims_BP$category <- row.names(D21MagentaAmb_GOslims_BP) #save rownames as category

D21MagentaAmb_MF <- D21MagentaAmb %>% # MF - all GO terms upregulated
  filter(ontology=="MF")
D21MagentaAmb_MF_GO_collection <- GOCollection(D21MagentaAmb_MF$category) #Make library of query terms
D21MagentaAmb_GOslims_MF <- data.frame(goSlim(D21MagentaAmb_MF_GO_collection, slim, "MF")) #Find common parent terms to slim down our list
D21MagentaAmb_GOslims_MF$category <- row.names(D21MagentaAmb_GOslims_MF) #save rownames as category


# D21BlueAmb
D21BlueAmb_BP <- D21BlueAmb %>% # BP - all GO terms upregulated
  filter(ontology=="BP")
D21BlueAmb_BP_GO_collection <- GOCollection(D21BlueAmb_BP$category) #Make library of query terms
D21BlueAmb_GOslims_BP <- data.frame(goSlim(D21BlueAmb_BP_GO_collection, slim, "BP")) #Find common parent terms to slim down our list
D21BlueAmb_GOslims_BP$category <- row.names(D21BlueAmb_GOslims_BP) #save rownames as category

D21BlueAmb_MF <- D21BlueAmb %>% # MF - all GO terms upregulated
  filter(ontology=="MF")
D21BlueAmb_MF_GO_collection <- GOCollection(D21BlueAmb_MF$category) #Make library of query terms
D21BlueAmb_GOslims_MF <- data.frame(goSlim(D21BlueAmb_MF_GO_collection, slim, "MF")) #Find common parent terms to slim down our list
D21BlueAmb_GOslims_MF$category <- row.names(D21BlueAmb_GOslims_MF) #save rownames as category


# ALL MODULES MODERATE > AMBEINT


# D7YellowMod
D7YellowMod_BP <- D7YellowMod %>% # BP - all GO terms upregulated
  filter(ontology=="BP")
D7YellowMod_BP_GO_collection <- GOCollection(D7YellowMod_BP$category) #Make library of query terms
D7YellowMod_GOslims_BP <- data.frame(goSlim(D7YellowMod_BP_GO_collection, slim, "BP")) #Find common parent terms to slim down our list
D7YellowMod_GOslims_BP$category <- row.names(D7YellowMod_GOslims_BP) #save rownames as category

D7YellowMod_MF <- D7YellowMod %>% # MF - all GO terms upregulated
  filter(ontology=="MF")
D7YellowMod_MF_GO_collection <- GOCollection(D7YellowMod_MF$category) #Make library of query terms
D7YellowMod_GOslims_MF <- data.frame(goSlim(D7YellowMod_MF_GO_collection, slim, "MF")) #Find common parent terms to slim down our list
D7YellowMod_GOslims_MF$category <- row.names(D7YellowMod_GOslims_MF) #save rownames as category


# D14BlackMod
D14BlackMod_BP <- D14BlackMod %>% # BP - all GO terms upregulated
  filter(ontology=="BP")
D14BlackMod_BP_GO_collection <- GOCollection(D14BlackMod_BP$category) #Make library of query terms
D14BlackMod_GOslims_BP <- data.frame(goSlim(D14BlackMod_BP_GO_collection, slim, "BP")) #Find common parent terms to slim down our list
D14BlackMod_GOslims_BP$category <- row.names(D14BlackMod_GOslims_BP) #save rownames as category

D14BlackMod_MF <- D14BlackMod %>% # MF - all GO terms upregulated
  filter(ontology=="MF")
D14BlackMod_MF_GO_collection <- GOCollection(D14BlackMod_MF$category) #Make library of query terms
D14BlackMod_GOslims_MF <- data.frame(goSlim(D14BlackMod_MF_GO_collection, slim, "MF")) #Find common parent terms to slim down our list
D14BlackMod_GOslims_MF$category <- row.names(D14BlackMod_GOslims_MF) #save rownames as category


# D21YellowMod
D21YellowMod_BP <- D21YellowMod %>% # BP - all GO terms upregulated
  filter(ontology=="BP")
D21YellowMod_BP_GO_collection <- GOCollection(D21YellowMod_BP$category) #Make library of query terms
D21YellowMod_GOslims_BP <- data.frame(goSlim(D21YellowMod_BP_GO_collection, slim, "BP")) #Find common parent terms to slim down our list
D21YellowMod_GOslims_BP$category <- row.names(D21YellowMod_GOslims_BP) #save rownames as category

D21YellowMod_MF <- D21YellowMod %>% # MF - all GO terms upregulated
  filter(ontology=="MF")
D21YellowMod_MF_GO_collection <- GOCollection(D21YellowMod_MF$category) #Make library of query terms
D21YellowMod_GOslims_MF <- data.frame(goSlim(D21YellowMod_MF_GO_collection, slim, "MF")) #Find common parent terms to slim down our list
D21YellowMod_GOslims_MF$category <- row.names(D21YellowMod_GOslims_MF) #save rownames as category

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
#Run function for MF and BP terms
# D7BrownAmb
D7BrownAmb_BPslim <- mappedIds(D7BrownAmb_GOslims_BP, D7BrownAmb_BP_GO_collection, GOBPOFFSPRING)
D7BrownAmb_MFslim <- mappedIds(D7BrownAmb_GOslims_MF, D7BrownAmb_MF_GO_collection, GOMFOFFSPRING)
# D14BrownAmb
D14BrownAmb_BPslim <- mappedIds(D14BrownAmb_GOslims_BP, D14BrownAmb_BP_GO_collection, GOBPOFFSPRING)
D14BrownAmb_MFslim <- mappedIds(D14BrownAmb_GOslims_MF, D14BrownAmb_MF_GO_collection, GOMFOFFSPRING)
# D21MagentaAmb
D21MagentaAmb_BPslim <- mappedIds(D21MagentaAmb_GOslims_BP, D21MagentaAmb_BP_GO_collection, GOBPOFFSPRING)
D21MagentaAmb_MFslim <- mappedIds(D21MagentaAmb_GOslims_MF, D21MagentaAmb_MF_GO_collection, GOMFOFFSPRING)
# D21BlueAmb
D21BlueAmb_BPslim <- mappedIds(D21BlueAmb_GOslims_BP, D21BlueAmb_BP_GO_collection, GOBPOFFSPRING)
D21BlueAmb_MFslim <- mappedIds(D21BlueAmb_GOslims_MF, D21BlueAmb_MF_GO_collection, GOMFOFFSPRING)


# D7YellowMod
D7YellowMod_BPslim <- mappedIds(D7YellowMod_GOslims_BP, D7YellowMod_BP_GO_collection, GOBPOFFSPRING)
D7YellowMod_MFslim <- mappedIds(D7YellowMod_GOslims_MF, D7YellowMod_MF_GO_collection, GOMFOFFSPRING)
# D14BlackMod
D14BlackMod_BPslim <- mappedIds(D14BlackMod_GOslims_BP, D14BlackMod_BP_GO_collection, GOBPOFFSPRING)
D14BlackMod_MFslim <- mappedIds(D14BlackMod_GOslims_MF, D14BlackMod_MF_GO_collection, GOMFOFFSPRING)
# D21YellowMod
D21YellowMod_BPslim <- mappedIds(D21YellowMod_GOslims_BP, D21YellowMod_BP_GO_collection, GOBPOFFSPRING)
D21YellowMod_MFslim <- mappedIds(D21YellowMod_GOslims_MF, D21YellowMod_MF_GO_collection, GOMFOFFSPRING)

# ====================================================================================
# Build final GO slim data set for plotting
# - call the annotation file - build a master sheet of unique rows for genes and GO terms (multiple gene IDs for each GO term annotated)
# - use this datasetto filter by upreg or down reg DE genes (DESEq2 directionality) - then for loop into the slim data to call unique genes in each GOslim bin
# - output column 'Gene.Count' to the final GOslim data table (i.e. All.UP_BPslim_final == all DESeq2 primary effect genes table, upregualted DEGs only, Biological Process GOslim/GOterms only)
# ====================================================================================
library(data.table) # for the setDT fxn in this cluster... 

# Build a master list of all genes and GO terms
Geoduck_annotation <- read.delim2(file="C:/Users/samjg/Documents/My_Projects/Pgenerosa_TagSeq_Metabolomics/TagSeq/Seq_details/Panopea-generosa-genes-annotations.txt", header=F) # Load themaster Pgenerosa gene list 
Pgen_GOterms <- Geoduck_annotation %>% dplyr::select(c('V1','V8')) # select only two columns - those with the gene IDs and those with the GO terms
Pgen_GOterms2 <- strsplit(Pgen_GOterms$V8, split = "; ") # create a string splitting by delimiter '; ' - view the data to see that this separates each GO term entry in the string
Pgen_GOterms2 <- data.frame(gene.ID = rep(Pgen_GOterms$V1, sapply(Pgen_GOterms2, length)), Go.terms = unlist(Pgen_GOterms2)) # create new dataframe 'Pgen_GOterms2' listing genes for each GO term (MUCH longer!)
Pgen_GOterms2 <- na.omit(Pgen_GOterms2) # ommit the NAs  - genes without GO annotation


# ====================================================================================
# Ambient > Moderate Modules - prep GO slim final table (number of genes in each GO slim bin)
#
# ...to do this I call all unique genes in the module that contain the go term before merging by GOslim
# ====================================================================================

# D7BrownAmb  ================================================================= #
head(d7_Mod.Brown) # ESSENTIAL THAT YOU HAVE THIS LOADED! - used to filter 'Pgen_GOterms2' in the for loop
#BP
D7BrownAmb_BPslim <- filter(D7BrownAmb_BPslim, Count>5 & Term!="biological_process") #filter out empty slims and term "biological process"
BPsplitted <- strsplit(as.character(D7BrownAmb_BPslim$go_terms), ";") #split into multiple GO ids
D7BrownAmb_BPslim$BPsplitted <- BPsplitted
for (i in 1:nrow(D7BrownAmb_BPslim)) {
  table <- data.frame(GOlist = unlist(D7BrownAmb_BPslim[i,6])) # call the BPsplitted column of characters and create a small table to filter
  table <- unique(table)
  Pgen_module <- Pgen_GOterms2 %>% dplyr::filter(gene.ID %in% d7_Mod.Brown$X) # d7_Mod.Brown$X calls the gene names in the origin Module membership dataframe at the start of this script
  Pgen_loop <- Pgen_module %>% dplyr::filter(Go.terms %in% table$GOlist) # filter Gene IDs with the GO term
  Pgen_geneIDs <- Pgen_loop[-2] # ommit the GO.terms to call unique gene calls 
  Pgen_geneIDs <- unique(Pgen_geneIDs) # call unique Gene calls (unique genes that had the GO term within each of the GOslim bins)
  D7BrownAmb_BPslim$Gene.Count[i] <- nrow(Pgen_geneIDs)  # count of unique GeneIDs in each GOslim bin
  D7BrownAmb_BPslim$Gene.IDs[[i]] <- vapply((Pgen_geneIDs$gene.ID), paste, collapse = ";", character(1L))} # name of each unique gene.id in each GOslim bin
D7BrownAmb_BPslim_A <- data.frame(Term = rep.int(D7BrownAmb_BPslim$Term, sapply(BPsplitted, length)), go_term = unlist(BPsplitted)) #list all
D7BrownAmb_BPslim_B <- merge(D7BrownAmb_BPslim_A, D7BrownAmb_BPslim, by="Term") #Add back counts, term, and category info
D7BrownAmb_BPslim_C <- unique(setDT(D7BrownAmb_BPslim_B)[order(go_term, -Gene.Count)], by = "category") #remove duplicate offspring terms, keeping only those in the larger umbrella term (Count number)
D7BrownAmb_BPslim_final <- data.frame(slim_term=D7BrownAmb_BPslim_C$Term, slim_cat=D7BrownAmb_BPslim_C$category, category=D7BrownAmb_BPslim_C$go_term, Gene.Count=D7BrownAmb_BPslim_C$Gene.Count, GO.Count=D7BrownAmb_BPslim_C$Count) #rename columns) #rename columns
#MF
D7BrownAmb_MFslim <- filter(D7BrownAmb_MFslim, Count>2 & Term!="molecular_function") #filter out empty slims and term "biological process"
MFsplitted <- strsplit(as.character(D7BrownAmb_MFslim$go_terms), ";") #split into multiple GO ids
D7BrownAmb_MFslim$MFsplitted <- MFsplitted
for (i in 1:nrow(D7BrownAmb_MFslim)) {
  table <- data.frame(GOlist = unlist(D7BrownAmb_MFslim[i,6])) # call the MFsplitted column of characters and create a small table to filter
  table <- unique(table)
  Pgen_module <- Pgen_GOterms2 %>% dplyr::filter(gene.ID %in% d7_Mod.Brown$X) # d7_Mod.Brown$X calls the gene names in the origin Module membership dataframe at the start of this script
  Pgen_loop <- Pgen_module %>% dplyr::filter(Go.terms %in% table$GOlist) # filter Gene IDs with the GO term
  Pgen_geneIDs <- Pgen_loop[-2] # ommit the GO.terms to call unique gene calls 
  Pgen_geneIDs <- unique(Pgen_geneIDs) # call unique Gene calls (unique genes that had the GO term within each of the GOslim bins)
  D7BrownAmb_MFslim$Gene.Count[i] <- nrow(Pgen_geneIDs)  # count of unique GeneIDs in each GOslim bin
  D7BrownAmb_MFslim$Gene.IDs[[i]] <- vapply((Pgen_geneIDs$gene.ID), paste, collapse = ";", character(1L))} # name of each unique gene.id in each GOslim bin
D7BrownAmb_MFslim_A <- data.frame(Term = rep.int(D7BrownAmb_MFslim$Term, sapply(MFsplitted, length)), go_term = unlist(MFsplitted)) #list all
D7BrownAmb_MFslim_B <- merge(D7BrownAmb_MFslim_A, D7BrownAmb_MFslim, by="Term") #Add back counts, term, and category info
D7BrownAmb_MFslim_C <- unique(setDT(D7BrownAmb_MFslim_B)[order(go_term, -Gene.Count)], by = "category") #remove duplicate offspring terms, keeping only those in the larger umbrella term (Count number)
D7BrownAmb_MFslim_final <- data.frame(slim_term=D7BrownAmb_MFslim_C$Term, slim_cat=D7BrownAmb_MFslim_C$category, category=D7BrownAmb_MFslim_C$go_term, Gene.Count=D7BrownAmb_MFslim_C$Gene.Count, GO.Count=D7BrownAmb_MFslim_C$Count) #rename columns) #rename columns


# D14BrownAmb  ================================================================= #
head(d14_Mod.Brown) # ESSENTIAL THAT YOU HAVE THIS LOADED! - used to filter 'Pgen_GOterms2' in the for loop
#BP
D14BrownAmb_BPslim <- filter(D14BrownAmb_BPslim, Count>5 & Term!="biological_process") #filter out empty slims and term "biological process"
BPsplitted <- strsplit(as.character(D14BrownAmb_BPslim$go_terms), ";") #split into multiple GO ids
D14BrownAmb_BPslim$BPsplitted <- BPsplitted
for (i in 1:nrow(D14BrownAmb_BPslim)) {
  table <- data.frame(GOlist = unlist(D14BrownAmb_BPslim[i,6])) # call the BPsplitted column of characters and create a small table to filter
  table <- unique(table)
  Pgen_module <- Pgen_GOterms2 %>% dplyr::filter(gene.ID %in% d14_Mod.Brown$X) # d14_Mod.Brown$X calls the gene names in the origin Module membership dataframe at the start of this script
  Pgen_loop <- Pgen_module %>% dplyr::filter(Go.terms %in% table$GOlist) # filter Gene IDs with the GO term
  Pgen_geneIDs <- Pgen_loop[-2] # ommit the GO.terms to call unique gene calls 
  Pgen_geneIDs <- unique(Pgen_geneIDs) # call unique Gene calls (unique genes that had the GO term within each of the GOslim bins)
  D14BrownAmb_BPslim$Gene.Count[i] <- nrow(Pgen_geneIDs)  # count of unique GeneIDs in each GOslim bin
  D14BrownAmb_BPslim$Gene.IDs[[i]] <- vapply((Pgen_geneIDs$gene.ID), paste, collapse = ";", character(1L))} # name of each unique gene.id in each GOslim bin
D14BrownAmb_BPslim_A <- data.frame(Term = rep.int(D14BrownAmb_BPslim$Term, sapply(BPsplitted, length)), go_term = unlist(BPsplitted)) #list all
D14BrownAmb_BPslim_B <- merge(D14BrownAmb_BPslim_A, D14BrownAmb_BPslim, by="Term") #Add back counts, term, and category info
D14BrownAmb_BPslim_C <- unique(setDT(D14BrownAmb_BPslim_B)[order(go_term, -Gene.Count)], by = "category") #remove duplicate offspring terms, keeping only those in the larger umbrella term (Count number)
D14BrownAmb_BPslim_final <- data.frame(slim_term=D14BrownAmb_BPslim_C$Term, slim_cat=D14BrownAmb_BPslim_C$category, category=D14BrownAmb_BPslim_C$go_term, Gene.Count=D14BrownAmb_BPslim_C$Gene.Count, GO.Count=D14BrownAmb_BPslim_C$Count) #rename columns) #rename columns
#MF
D14BrownAmb_MFslim <- filter(D14BrownAmb_MFslim, Count>2 & Term!="molecular_function") #filter out empty slims and term "biological process"
MFsplitted <- strsplit(as.character(D14BrownAmb_MFslim$go_terms), ";") #split into multiple GO ids
D14BrownAmb_MFslim$MFsplitted <- MFsplitted
for (i in 1:nrow(D14BrownAmb_MFslim)) {
  table <- data.frame(GOlist = unlist(D14BrownAmb_MFslim[i,6])) # call the MFsplitted column of characters and create a small table to filter
  table <- unique(table)
  Pgen_module <- Pgen_GOterms2 %>% dplyr::filter(gene.ID %in% d14_Mod.Brown$X) # d14_Mod.Brown$X calls the gene names in the origin Module membership dataframe at the start of this script
  Pgen_loop <- Pgen_module %>% dplyr::filter(Go.terms %in% table$GOlist) # filter Gene IDs with the GO term
  Pgen_geneIDs <- Pgen_loop[-2] # ommit the GO.terms to call unique gene calls 
  Pgen_geneIDs <- unique(Pgen_geneIDs) # call unique Gene calls (unique genes that had the GO term within each of the GOslim bins)
  D14BrownAmb_MFslim$Gene.Count[i] <- nrow(Pgen_geneIDs)  # count of unique GeneIDs in each GOslim bin
  D14BrownAmb_MFslim$Gene.IDs[[i]] <- vapply((Pgen_geneIDs$gene.ID), paste, collapse = ";", character(1L))} # name of each unique gene.id in each GOslim bin
D14BrownAmb_MFslim_A <- data.frame(Term = rep.int(D14BrownAmb_MFslim$Term, sapply(MFsplitted, length)), go_term = unlist(MFsplitted)) #list all
D14BrownAmb_MFslim_B <- merge(D14BrownAmb_MFslim_A, D14BrownAmb_MFslim, by="Term") #Add back counts, term, and category info
D14BrownAmb_MFslim_C <- unique(setDT(D14BrownAmb_MFslim_B)[order(go_term, -Gene.Count)], by = "category") #remove duplicate offspring terms, keeping only those in the larger umbrella term (Count number)
D14BrownAmb_MFslim_final <- data.frame(slim_term=D14BrownAmb_MFslim_C$Term, slim_cat=D14BrownAmb_MFslim_C$category, category=D14BrownAmb_MFslim_C$go_term, Gene.Count=D14BrownAmb_MFslim_C$Gene.Count, GO.Count=D14BrownAmb_MFslim_C$Count) #rename columns) #rename columns


# D21MagentaAmb  ================================================================= #
head(d21_Mod.Magenta) # ESSENTIAL THAT YOU HAVE THIS LOADED! - used to filter 'Pgen_GOterms2' in the for loop
#BP
D21MagentaAmb_BPslim <- filter(D21MagentaAmb_BPslim, Count>5 & Term!="biological_process") #filter out empty slims and term "biological process"
BPsplitted <- strsplit(as.character(D21MagentaAmb_BPslim$go_terms), ";") #split into multiple GO ids
D21MagentaAmb_BPslim$BPsplitted <- BPsplitted
for (i in 1:nrow(D21MagentaAmb_BPslim)) {
  table <- data.frame(GOlist = unlist(D21MagentaAmb_BPslim[i,6])) # call the BPsplitted column of characters and create a small table to filter
  table <- unique(table)
  Pgen_module <- Pgen_GOterms2 %>% dplyr::filter(gene.ID %in% d21_Mod.Magenta$X) # d21_Mod.Magenta$X calls the gene names in the origin Module membership dataframe at the start of this script
  Pgen_loop <- Pgen_module %>% dplyr::filter(Go.terms %in% table$GOlist) # filter Gene IDs with the GO term
  Pgen_geneIDs <- Pgen_loop[-2] # ommit the GO.terms to call unique gene calls 
  Pgen_geneIDs <- unique(Pgen_geneIDs) # call unique Gene calls (unique genes that had the GO term within each of the GOslim bins)
  D21MagentaAmb_BPslim$Gene.Count[i] <- nrow(Pgen_geneIDs)  # count of unique GeneIDs in each GOslim bin
  D21MagentaAmb_BPslim$Gene.IDs[[i]] <- vapply((Pgen_geneIDs$gene.ID), paste, collapse = ";", character(1L))} # name of each unique gene.id in each GOslim bin
D21MagentaAmb_BPslim_A <- data.frame(Term = rep.int(D21MagentaAmb_BPslim$Term, sapply(BPsplitted, length)), go_term = unlist(BPsplitted)) #list all
D21MagentaAmb_BPslim_B <- merge(D21MagentaAmb_BPslim_A, D21MagentaAmb_BPslim, by="Term") #Add back counts, term, and category info
D21MagentaAmb_BPslim_C <- unique(setDT(D21MagentaAmb_BPslim_B)[order(go_term, -Gene.Count)], by = "category") #remove duplicate offspring terms, keeping only those in the larger umbrella term (Count number)
D21MagentaAmb_BPslim_final <- data.frame(slim_term=D21MagentaAmb_BPslim_C$Term, slim_cat=D21MagentaAmb_BPslim_C$category, category=D21MagentaAmb_BPslim_C$go_term, Gene.Count=D21MagentaAmb_BPslim_C$Gene.Count, GO.Count=D21MagentaAmb_BPslim_C$Count) #rename columns) #rename columns
#MF
D21MagentaAmb_MFslim <- filter(D21MagentaAmb_MFslim, Count>2 & Term!="molecular_function") #filter out empty slims and term "biological process"
MFsplitted <- strsplit(as.character(D21MagentaAmb_MFslim$go_terms), ";") #split into multiple GO ids
D21MagentaAmb_MFslim$MFsplitted <- MFsplitted
for (i in 1:nrow(D21MagentaAmb_MFslim)) {
  table <- data.frame(GOlist = unlist(D21MagentaAmb_MFslim[i,6])) # call the MFsplitted column of characters and create a small table to filter
  table <- unique(table)
  Pgen_module <- Pgen_GOterms2 %>% dplyr::filter(gene.ID %in% d21_Mod.Magenta$X) # d21_Mod.Magenta$X calls the gene names in the origin Module membership dataframe at the start of this script
  Pgen_loop <- Pgen_module %>% dplyr::filter(Go.terms %in% table$GOlist) # filter Gene IDs with the GO term
  Pgen_geneIDs <- Pgen_loop[-2] # ommit the GO.terms to call unique gene calls 
  Pgen_geneIDs <- unique(Pgen_geneIDs) # call unique Gene calls (unique genes that had the GO term within each of the GOslim bins)
  D21MagentaAmb_MFslim$Gene.Count[i] <- nrow(Pgen_geneIDs)  # count of unique GeneIDs in each GOslim bin
  D21MagentaAmb_MFslim$Gene.IDs[[i]] <- vapply((Pgen_geneIDs$gene.ID), paste, collapse = ";", character(1L))} # name of each unique gene.id in each GOslim bin
D21MagentaAmb_MFslim_A <- data.frame(Term = rep.int(D21MagentaAmb_MFslim$Term, sapply(MFsplitted, length)), go_term = unlist(MFsplitted)) #list all
D21MagentaAmb_MFslim_B <- merge(D21MagentaAmb_MFslim_A, D21MagentaAmb_MFslim, by="Term") #Add back counts, term, and category info
D21MagentaAmb_MFslim_C <- unique(setDT(D21MagentaAmb_MFslim_B)[order(go_term, -Gene.Count)], by = "category") #remove duplicate offspring terms, keeping only those in the larger umbrella term (Count number)
D21MagentaAmb_MFslim_final <- data.frame(slim_term=D21MagentaAmb_MFslim_C$Term, slim_cat=D21MagentaAmb_MFslim_C$category, category=D21MagentaAmb_MFslim_C$go_term, Gene.Count=D21MagentaAmb_MFslim_C$Gene.Count, GO.Count=D21MagentaAmb_MFslim_C$Count) #rename columns) #rename columns



# D21BlueAmb  ================================================================= #
head(d21_Mod.blue) # ESSENTIAL THAT YOU HAVE THIS LOADED! - used to filter 'Pgen_GOterms2' in the for loop
#BP
D21BlueAmb_BPslim <- filter(D21BlueAmb_BPslim, Count>5 & Term!="biological_process") #filter out empty slims and term "biological process"
BPsplitted <- strsplit(as.character(D21BlueAmb_BPslim$go_terms), ";") #split into multiple GO ids
D21BlueAmb_BPslim$BPsplitted <- BPsplitted
for (i in 1:nrow(D21BlueAmb_BPslim)) {
  table <- data.frame(GOlist = unlist(D21BlueAmb_BPslim[i,6])) # call the BPsplitted column of characters and create a small table to filter
  table <- unique(table)
  Pgen_module <- Pgen_GOterms2 %>% dplyr::filter(gene.ID %in% d21_Mod.blue$X) # d21_Mod.Blue$X calls the gene names in the origin Module membership dataframe at the start of this script
  Pgen_loop <- Pgen_module %>% dplyr::filter(Go.terms %in% table$GOlist) # filter Gene IDs with the GO term
  Pgen_geneIDs <- Pgen_loop[-2] # ommit the GO.terms to call unique gene calls 
  Pgen_geneIDs <- unique(Pgen_geneIDs) # call unique Gene calls (unique genes that had the GO term within each of the GOslim bins)
  D21BlueAmb_BPslim$Gene.Count[i] <- nrow(Pgen_geneIDs)  # count of unique GeneIDs in each GOslim bin
  D21BlueAmb_BPslim$Gene.IDs[[i]] <- vapply((Pgen_geneIDs$gene.ID), paste, collapse = ";", character(1L))} # name of each unique gene.id in each GOslim bin
D21BlueAmb_BPslim_A <- data.frame(Term = rep.int(D21BlueAmb_BPslim$Term, sapply(BPsplitted, length)), go_term = unlist(BPsplitted)) #list all
D21BlueAmb_BPslim_B <- merge(D21BlueAmb_BPslim_A, D21BlueAmb_BPslim, by="Term") #Add back counts, term, and category info
D21BlueAmb_BPslim_C <- unique(setDT(D21BlueAmb_BPslim_B)[order(go_term, -Gene.Count)], by = "category") #remove duplicate offspring terms, keeping only those in the larger umbrella term (Count number)
D21BlueAmb_BPslim_final <- data.frame(slim_term=D21BlueAmb_BPslim_C$Term, slim_cat=D21BlueAmb_BPslim_C$category, category=D21BlueAmb_BPslim_C$go_term, Gene.Count=D21BlueAmb_BPslim_C$Gene.Count, GO.Count=D21BlueAmb_BPslim_C$Count) #rename columns) #rename columns
#MF
D21BlueAmb_MFslim <- filter(D21BlueAmb_MFslim, Count>2 & Term!="molecular_function") #filter out empty slims and term "biological process"
MFsplitted <- strsplit(as.character(D21BlueAmb_MFslim$go_terms), ";") #split into multiple GO ids
D21BlueAmb_MFslim$MFsplitted <- MFsplitted
for (i in 1:nrow(D21BlueAmb_MFslim)) {
  table <- data.frame(GOlist = unlist(D21BlueAmb_MFslim[i,6])) # call the MFsplitted column of characters and create a small table to filter
  table <- unique(table)
  Pgen_module <- Pgen_GOterms2 %>% dplyr::filter(gene.ID %in% d21_Mod.blue$X) # d21_Mod.Blue$X calls the gene names in the origin Module membership dataframe at the start of this script
  Pgen_loop <- Pgen_module %>% dplyr::filter(Go.terms %in% table$GOlist) # filter Gene IDs with the GO term
  Pgen_geneIDs <- Pgen_loop[-2] # ommit the GO.terms to call unique gene calls 
  Pgen_geneIDs <- unique(Pgen_geneIDs) # call unique Gene calls (unique genes that had the GO term within each of the GOslim bins)
  D21BlueAmb_MFslim$Gene.Count[i] <- nrow(Pgen_geneIDs)  # count of unique GeneIDs in each GOslim bin
  D21BlueAmb_MFslim$Gene.IDs[[i]] <- vapply((Pgen_geneIDs$gene.ID), paste, collapse = ";", character(1L))} # name of each unique gene.id in each GOslim bin
D21BlueAmb_MFslim_A <- data.frame(Term = rep.int(D21BlueAmb_MFslim$Term, sapply(MFsplitted, length)), go_term = unlist(MFsplitted)) #list all
D21BlueAmb_MFslim_B <- merge(D21BlueAmb_MFslim_A, D21BlueAmb_MFslim, by="Term") #Add back counts, term, and category info
D21BlueAmb_MFslim_C <- unique(setDT(D21BlueAmb_MFslim_B)[order(go_term, -Gene.Count)], by = "category") #remove duplicate offspring terms, keeping only those in the larger umbrella term (Count number)
D21BlueAmb_MFslim_final <- data.frame(slim_term=D21BlueAmb_MFslim_C$Term, slim_cat=D21BlueAmb_MFslim_C$category, category=D21BlueAmb_MFslim_C$go_term, Gene.Count=D21BlueAmb_MFslim_C$Gene.Count, GO.Count=D21BlueAmb_MFslim_C$Count) #rename columns) #rename columns


# ====================================================================================
# Ambient > Moderate Modules - STACKED BAR PLOT 
#
# 
# ====================================================================================
# create color palette  --------------------------------------------------------------------------------------- #
library(RColorBrewer)
# cbp1 <- c("#999999", "#E69F00", "#56B4E9", "#009E73",
#           "#F0E442", "#0072B2", "#D55E00", "#CC79A7") # color blind pallette
OrangeRed <- brewer.pal(3, "OrRd") 
GreenBlue <- brewer.pal(3, "GnBu") 

# add ontology column  ---------------------------------------------------------------------------------------- #
D7BrownAmb_MFslim_final$Ontolgy     <- "D7Brown_MF"
D14BrownAmb_MFslim_final$Ontolgy    <- "D14Brown_MF"
D21MagentaAmb_MFslim_final$Ontolgy  <- "D21Magenta_MF"
D21BlueAmb_MFslim_final$Ontolgy     <- "D21Blue_MF"

D7BrownAmb_BPslim_final$Ontolgy     <- "D7Brown_BP"
D14BrownAmb_BPslim_final$Ontolgy    <- "D14Brown_BP"
D21MagentaAmb_BPslim_final$Ontolgy  <- "D21Magenta_BP"
D21BlueAmb_BPslim_final$Ontolgy     <- "D21Blue_BP"


# BIND DATA TO FACET THE PLOTS BY UP AND DOWN REG   ---------------------------------------------------------- #
# bp
BP_Amb <- rbind(D7BrownAmb_BPslim_final, D14BrownAmb_BPslim_final, D21MagentaAmb_BPslim_final, D21BlueAmb_BPslim_final) # merge the data by binding rows 
BP_Amb$Ont <- "BP" # create a new common clumn to call in the plot 
BP_Amb$Ontolgy <- factor(BP_Amb$Ontolgy, levels = c("D7Brown_BP", "D14Brown_BP", "D21Magenta_BP", "D21Blue_BP"))# reorder the facotr level for the facet wrap plot 
BP_Amb_filtered <- BP_Amb %>%  dplyr::filter(Gene.Count > 1) # ommit all with gene counts <1
# mf
MF_Amb <- rbind(D7BrownAmb_MFslim_final, D14BrownAmb_MFslim_final, D21MagentaAmb_MFslim_final, D21BlueAmb_MFslim_final) # merge the data by binding rows 
MF_Amb$Ont <- "MF" # create a new common clumn to call in the plot 
MF_Amb$Ontolgy <- factor(MF_Amb$Ontolgy, levels = c("D7Brown_MF", "D14Brown_MF", "D21Magenta_MF", "D21Blue_MF")) # reorder the facotr level for the facet wrap plot 
MF_Amb_filtered <- MF_Amb %>%  dplyr::filter(Gene.Count > 1) # ommit all with gene counts <1

# PLOTTING --------------------------------------------------------------------------------------------------- #
BP_Amb_filtered$slim_term <- factor(BP_Amb_filtered$slim_term ,levels=rev(unique(BP_Amb_filtered$slim_term))) # make slim term alphabetical for plotting
BP_Plot_PrimEffect_AMB <-ggplot(data = BP_Amb_filtered, aes(x = Ont, y = slim_term)) + 
                          geom_tile(aes(fill=Gene.Count, width = 1)) + 
                          scale_fill_gradient(low = "thistle1", high = "steelblue4") +
                          facet_grid(~Ontolgy, labeller = label_wrap_gen(width = 10, multi_line = TRUE))+
                          theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                                             panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                                             strip.text.y = element_text(angle=0, size = 11, face = "bold"),
                                             strip.text.x = element_text(size = 12, face = "bold"),
                                             axis.title.x = element_blank(),
                                             axis.title.y = element_text(size=15),
                                             axis.text = element_text(size = 12), legend.position = "right",
                                             plot.margin = unit(c(0,1,0,0.25), "cm"))+
                          ggtitle('GOslim Biological Process: WGCNA Ambient > Moderate Primary Effect')

MF_Amb_filtered$slim_term <- factor(MF_Amb_filtered$slim_term ,levels=rev(unique(MF_Amb_filtered$slim_term))) # make slim term alphabetical for plotting
MF_Plot_PrimEffect_AMB <- ggplot(data = MF_Amb_filtered, aes(x = Ont, y = slim_term)) + 
                          geom_tile(aes(fill=Gene.Count, width = 1)) + 
                          scale_fill_gradient(low = "thistle1", high = "tomato4") +
                          facet_grid(~Ontolgy, scales = "free_y", labeller = label_wrap_gen(width = 10, multi_line = TRUE))+
                          theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                                             panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                                             strip.text.y = element_text(angle=0, size = 11, face = "bold"),
                                             strip.text.x = element_text(size = 12, face = "bold"),
                                             axis.title.x = element_blank(),
                                             axis.title.y = element_text(size=15),
                                             axis.text = element_text(size = 12), legend.position = "right",
                                             plot.margin = unit(c(0,1,0,0.25), "cm"))+
                          ggtitle('GOslim Molecular Function: WGCNA Ambient > Moderate Primary Effect')

# SAVE PLOTS  ------------------------------------------------------------------------------------------------- #
pdf(paste("Analysis//Output/GO/WGCNA_goseq/GOslim_PrimEff_Ambient_BP.pdf", sep =''), width=10, height=5)
print(ggarrange(BP_Plot_PrimEffect_AMB,         
                plotlist = NULL,
                ncol = 1,
                nrow = 1,
                labels = NULL))
dev.off()     

pdf(paste("Analysis//Output/GO/WGCNA_goseq/GOslim_PrimEff_Ambient_MF.pdf", sep =''), width=10, height=5)
print(ggarrange(MF_Plot_PrimEffect_AMB,         
                plotlist = NULL,
                ncol = 1,
                nrow = 1,
                labels = NULL))
dev.off()     






# ====================================================================================
# Moderate > Ambient Modules - prep GO slim final table (number of genes in each GO slim bin)
#
# ...to do this I call all unique genes in the module that contain the go term before merging by GOslim
# ====================================================================================

# D7YellowAmb  ================================================================= #
head(d7_Mod.Yellow) # ESSENTIAL THAT YOU HAVE THIS LOADED! - used to filter 'Pgen_GOterms2' in the for loop
#BP
D7YellowMod_BPslim <- filter(D7YellowMod_BPslim, Count>5 & Term!="biological_process") #filter out empty slims and term "biological process"
BPsplitted <- strsplit(as.character(D7YellowMod_BPslim$go_terms), ";") #split into multiple GO ids
D7YellowMod_BPslim$BPsplitted <- BPsplitted
for (i in 1:nrow(D7YellowMod_BPslim)) {
  table <- data.frame(GOlist = unlist(D7YellowMod_BPslim[i,6])) # call the BPsplitted column of characters and create a small table to filter
  table <- unique(table)
  Pgen_module <- Pgen_GOterms2 %>% dplyr::filter(gene.ID %in% d7_Mod.Yellow$X) # d7_Mod.Yellow$X calls the gene names in the origin Module membership dataframe at the start of this script
  Pgen_loop <- Pgen_module %>% dplyr::filter(Go.terms %in% table$GOlist) # filter Gene IDs with the GO term
  Pgen_geneIDs <- Pgen_loop[-2] # ommit the GO.terms to call unique gene calls 
  Pgen_geneIDs <- unique(Pgen_geneIDs) # call unique Gene calls (unique genes that had the GO term within each of the GOslim bins)
  D7YellowMod_BPslim$Gene.Count[i] <- nrow(Pgen_geneIDs)  # count of unique GeneIDs in each GOslim bin
  D7YellowMod_BPslim$Gene.IDs[[i]] <- vapply((Pgen_geneIDs$gene.ID), paste, collapse = ";", character(1L))} # name of each unique gene.id in each GOslim bin
D7YellowMod_BPslim_A <- data.frame(Term = rep.int(D7YellowMod_BPslim$Term, sapply(BPsplitted, length)), go_term = unlist(BPsplitted)) #list all
D7YellowMod_BPslim_B <- merge(D7YellowMod_BPslim_A, D7YellowMod_BPslim, by="Term") #Add back counts, term, and category info
D7YellowMod_BPslim_C <- unique(setDT(D7YellowMod_BPslim_B)[order(go_term, -Gene.Count)], by = "category") #remove duplicate offspring terms, keeping only those in the larger umbrella term (Count number)
D7YellowMod_BPslim_final <- data.frame(slim_term=D7YellowMod_BPslim_C$Term, slim_cat=D7YellowMod_BPslim_C$category, category=D7YellowMod_BPslim_C$go_term, Gene.Count=D7YellowMod_BPslim_C$Gene.Count, GO.Count=D7YellowMod_BPslim_C$Count) #rename columns) #rename columns
#MF
D7YellowMod_MFslim <- filter(D7YellowMod_MFslim, Count>2 & Term!="molecular_function") #filter out empty slims and term "biological process"
MFsplitted <- strsplit(as.character(D7YellowMod_MFslim$go_terms), ";") #split into multiple GO ids
D7YellowMod_MFslim$MFsplitted <- MFsplitted
for (i in 1:nrow(D7YellowMod_MFslim)) {
  table <- data.frame(GOlist = unlist(D7YellowMod_MFslim[i,6])) # call the MFsplitted column of characters and create a small table to filter
  table <- unique(table)
  Pgen_module <- Pgen_GOterms2 %>% dplyr::filter(gene.ID %in% d7_Mod.Yellow$X) # d7_Mod.Yellow$X calls the gene names in the origin Module membership dataframe at the start of this script
  Pgen_loop <- Pgen_module %>% dplyr::filter(Go.terms %in% table$GOlist) # filter Gene IDs with the GO term
  Pgen_geneIDs <- Pgen_loop[-2] # ommit the GO.terms to call unique gene calls 
  Pgen_geneIDs <- unique(Pgen_geneIDs) # call unique Gene calls (unique genes that had the GO term within each of the GOslim bins)
  D7YellowMod_MFslim$Gene.Count[i] <- nrow(Pgen_geneIDs)  # count of unique GeneIDs in each GOslim bin
  D7YellowMod_MFslim$Gene.IDs[[i]] <- vapply((Pgen_geneIDs$gene.ID), paste, collapse = ";", character(1L))} # name of each unique gene.id in each GOslim bin
D7YellowMod_MFslim_A <- data.frame(Term = rep.int(D7YellowMod_MFslim$Term, sapply(MFsplitted, length)), go_term = unlist(MFsplitted)) #list all
D7YellowMod_MFslim_B <- merge(D7YellowMod_MFslim_A, D7YellowMod_MFslim, by="Term") #Add back counts, term, and category info
D7YellowMod_MFslim_C <- unique(setDT(D7YellowMod_MFslim_B)[order(go_term, -Gene.Count)], by = "category") #remove duplicate offspring terms, keeping only those in the larger umbrella term (Count number)
D7YellowMod_MFslim_final <- data.frame(slim_term=D7YellowMod_MFslim_C$Term, slim_cat=D7YellowMod_MFslim_C$category, category=D7YellowMod_MFslim_C$go_term, Gene.Count=D7YellowMod_MFslim_C$Gene.Count, GO.Count=D7YellowMod_MFslim_C$Count) #rename columns) #rename columns


# D14BlackAmb  ================================================================= #
head(d14_Mod.Black) # ESSENTIAL THAT YOU HAVE THIS LOADED! - used to filter 'Pgen_GOterms2' in the for loop
#BP
D14BlackMod_BPslim <- filter(D14BlackMod_BPslim, Count>5 & Term!="biological_process") #filter out empty slims and term "biological process"
BPsplitted <- strsplit(as.character(D14BlackMod_BPslim$go_terms), ";") #split into multiple GO ids
D14BlackMod_BPslim$BPsplitted <- BPsplitted
for (i in 1:nrow(D14BlackMod_BPslim)) {
  table <- data.frame(GOlist = unlist(D14BlackMod_BPslim[i,6])) # call the BPsplitted column of characters and create a small table to filter
  table <- unique(table)
  Pgen_module <- Pgen_GOterms2 %>% dplyr::filter(gene.ID %in% d14_Mod.Black$X) # d14_Mod.Black$X calls the gene names in the origin Module membership dataframe at the start of this script
  Pgen_loop <- Pgen_module %>% dplyr::filter(Go.terms %in% table$GOlist) # filter Gene IDs with the GO term
  Pgen_geneIDs <- Pgen_loop[-2] # ommit the GO.terms to call unique gene calls 
  Pgen_geneIDs <- unique(Pgen_geneIDs) # call unique Gene calls (unique genes that had the GO term within each of the GOslim bins)
  D14BlackMod_BPslim$Gene.Count[i] <- nrow(Pgen_geneIDs)  # count of unique GeneIDs in each GOslim bin
  D14BlackMod_BPslim$Gene.IDs[[i]] <- vapply((Pgen_geneIDs$gene.ID), paste, collapse = ";", character(1L))} # name of each unique gene.id in each GOslim bin
D14BlackMod_BPslim_A <- data.frame(Term = rep.int(D14BlackMod_BPslim$Term, sapply(BPsplitted, length)), go_term = unlist(BPsplitted)) #list all
D14BlackMod_BPslim_B <- merge(D14BlackMod_BPslim_A, D14BlackMod_BPslim, by="Term") #Add back counts, term, and category info
D14BlackMod_BPslim_C <- unique(setDT(D14BlackMod_BPslim_B)[order(go_term, -Gene.Count)], by = "category") #remove duplicate offspring terms, keeping only those in the larger umbrella term (Count number)
D14BlackMod_BPslim_final <- data.frame(slim_term=D14BlackMod_BPslim_C$Term, slim_cat=D14BlackMod_BPslim_C$category, category=D14BlackMod_BPslim_C$go_term, Gene.Count=D14BlackMod_BPslim_C$Gene.Count, GO.Count=D14BlackMod_BPslim_C$Count) #rename columns) #rename columns
#MF
D14BlackMod_MFslim <- filter(D14BlackMod_MFslim, Count>2 & Term!="molecular_function") #filter out empty slims and term "biological process"
MFsplitted <- strsplit(as.character(D14BlackMod_MFslim$go_terms), ";") #split into multiple GO ids
D14BlackMod_MFslim$MFsplitted <- MFsplitted
for (i in 1:nrow(D14BlackMod_MFslim)) {
  table <- data.frame(GOlist = unlist(D14BlackMod_MFslim[i,6])) # call the MFsplitted column of characters and create a small table to filter
  table <- unique(table)
  Pgen_module <- Pgen_GOterms2 %>% dplyr::filter(gene.ID %in% d14_Mod.Black$X) # d14_Mod.Black$X calls the gene names in the origin Module membership dataframe at the start of this script
  Pgen_loop <- Pgen_module %>% dplyr::filter(Go.terms %in% table$GOlist) # filter Gene IDs with the GO term
  Pgen_geneIDs <- Pgen_loop[-2] # ommit the GO.terms to call unique gene calls 
  Pgen_geneIDs <- unique(Pgen_geneIDs) # call unique Gene calls (unique genes that had the GO term within each of the GOslim bins)
  D14BlackMod_MFslim$Gene.Count[i] <- nrow(Pgen_geneIDs)  # count of unique GeneIDs in each GOslim bin
  D14BlackMod_MFslim$Gene.IDs[[i]] <- vapply((Pgen_geneIDs$gene.ID), paste, collapse = ";", character(1L))} # name of each unique gene.id in each GOslim bin
D14BlackMod_MFslim_A <- data.frame(Term = rep.int(D14BlackMod_MFslim$Term, sapply(MFsplitted, length)), go_term = unlist(MFsplitted)) #list all
D14BlackMod_MFslim_B <- merge(D14BlackMod_MFslim_A, D14BlackMod_MFslim, by="Term") #Add back counts, term, and category info
D14BlackMod_MFslim_C <- unique(setDT(D14BlackMod_MFslim_B)[order(go_term, -Gene.Count)], by = "category") #remove duplicate offspring terms, keeping only those in the larger umbrella term (Count number)
D14BlackMod_MFslim_final <- data.frame(slim_term=D14BlackMod_MFslim_C$Term, slim_cat=D14BlackMod_MFslim_C$category, category=D14BlackMod_MFslim_C$go_term, Gene.Count=D14BlackMod_MFslim_C$Gene.Count, GO.Count=D14BlackMod_MFslim_C$Count) #rename columns) #rename columns


# D21YellowAmb  ================================================================= #
head(d21_Mod.Yellow) # ESSENTIAL THAT YOU HAVE THIS LOADED! - used to filter 'Pgen_GOterms2' in the for loop
#BP
D21YellowMod_BPslim <- filter(D21YellowMod_BPslim, Count>5 & Term!="biological_process") #filter out empty slims and term "biological process"
BPsplitted <- strsplit(as.character(D21YellowMod_BPslim$go_terms), ";") #split into multiple GO ids
D21YellowMod_BPslim$BPsplitted <- BPsplitted
for (i in 1:nrow(D21YellowMod_BPslim)) {
  table <- data.frame(GOlist = unlist(D21YellowMod_BPslim[i,6])) # call the BPsplitted column of characters and create a small table to filter
  table <- unique(table)
  Pgen_module <- Pgen_GOterms2 %>% dplyr::filter(gene.ID %in% d21_Mod.Yellow$X) # d21_Mod.Yellow$X calls the gene names in the origin Module membership dataframe at the start of this script
  Pgen_loop <- Pgen_module %>% dplyr::filter(Go.terms %in% table$GOlist) # filter Gene IDs with the GO term
  Pgen_geneIDs <- Pgen_loop[-2] # ommit the GO.terms to call unique gene calls 
  Pgen_geneIDs <- unique(Pgen_geneIDs) # call unique Gene calls (unique genes that had the GO term within each of the GOslim bins)
  D21YellowMod_BPslim$Gene.Count[i] <- nrow(Pgen_geneIDs)  # count of unique GeneIDs in each GOslim bin
  D21YellowMod_BPslim$Gene.IDs[[i]] <- vapply((Pgen_geneIDs$gene.ID), paste, collapse = ";", character(1L))} # name of each unique gene.id in each GOslim bin
D21YellowMod_BPslim_A <- data.frame(Term = rep.int(D21YellowMod_BPslim$Term, sapply(BPsplitted, length)), go_term = unlist(BPsplitted)) #list all
D21YellowMod_BPslim_B <- merge(D21YellowMod_BPslim_A, D21YellowMod_BPslim, by="Term") #Add back counts, term, and category info
D21YellowMod_BPslim_C <- unique(setDT(D21YellowMod_BPslim_B)[order(go_term, -Gene.Count)], by = "category") #remove duplicate offspring terms, keeping only those in the larger umbrella term (Count number)
D21YellowMod_BPslim_final <- data.frame(slim_term=D21YellowMod_BPslim_C$Term, slim_cat=D21YellowMod_BPslim_C$category, category=D21YellowMod_BPslim_C$go_term, Gene.Count=D21YellowMod_BPslim_C$Gene.Count, GO.Count=D21YellowMod_BPslim_C$Count) #rename columns) #rename columns
#MF
D21YellowMod_MFslim <- filter(D21YellowMod_MFslim, Count>2 & Term!="molecular_function") #filter out empty slims and term "biological process"
MFsplitted <- strsplit(as.character(D21YellowMod_MFslim$go_terms), ";") #split into multiple GO ids
D21YellowMod_MFslim$MFsplitted <- MFsplitted
for (i in 1:nrow(D21YellowMod_MFslim)) {
  table <- data.frame(GOlist = unlist(D21YellowMod_MFslim[i,6])) # call the MFsplitted column of characters and create a small table to filter
  table <- unique(table)
  Pgen_module <- Pgen_GOterms2 %>% dplyr::filter(gene.ID %in% d21_Mod.Yellow$X) # d21_Mod.Yellow$X calls the gene names in the origin Module membership dataframe at the start of this script
  Pgen_loop <- Pgen_module %>% dplyr::filter(Go.terms %in% table$GOlist) # filter Gene IDs with the GO term
  Pgen_geneIDs <- Pgen_loop[-2] # ommit the GO.terms to call unique gene calls 
  Pgen_geneIDs <- unique(Pgen_geneIDs) # call unique Gene calls (unique genes that had the GO term within each of the GOslim bins)
  D21YellowMod_MFslim$Gene.Count[i] <- nrow(Pgen_geneIDs)  # count of unique GeneIDs in each GOslim bin
  D21YellowMod_MFslim$Gene.IDs[[i]] <- vapply((Pgen_geneIDs$gene.ID), paste, collapse = ";", character(1L))} # name of each unique gene.id in each GOslim bin
D21YellowMod_MFslim_A <- data.frame(Term = rep.int(D21YellowMod_MFslim$Term, sapply(MFsplitted, length)), go_term = unlist(MFsplitted)) #list all
D21YellowMod_MFslim_B <- merge(D21YellowMod_MFslim_A, D21YellowMod_MFslim, by="Term") #Add back counts, term, and category info
D21YellowMod_MFslim_C <- unique(setDT(D21YellowMod_MFslim_B)[order(go_term, -Gene.Count)], by = "category") #remove duplicate offspring terms, keeping only those in the larger umbrella term (Count number)
D21YellowMod_MFslim_final <- data.frame(slim_term=D21YellowMod_MFslim_C$Term, slim_cat=D21YellowMod_MFslim_C$category, category=D21YellowMod_MFslim_C$go_term, Gene.Count=D21YellowMod_MFslim_C$Gene.Count, GO.Count=D21YellowMod_MFslim_C$Count) #rename columns) #rename columns

# ====================================================================================
# Moderate > Ambient Modules - STACKED BAR PLOT 
#
# 
# ====================================================================================
# create color palette  --------------------------------------------------------------------------------------- #
library(RColorBrewer)
# cbp1 <- c("#999999", "#E69F00", "#56B4E9", "#009E73",
#           "#F0E442", "#0072B2", "#D55E00", "#CC79A7") # color blind pallette
OrangeRed <- brewer.pal(3, "OrRd") 
GreenBlue <- brewer.pal(3, "GnBu") 

# add ontology column  ---------------------------------------------------------------------------------------- #
D7YellowMod_MFslim_final$Ontolgy     <- "D7Yellow_MF"
D14BlackMod_MFslim_final$Ontolgy     <- "D14Black_MF"
D21YellowMod_MFslim_final$Ontolgy    <- "D21Yellow_MF"

D7YellowMod_BPslim_final$Ontolgy     <- "D7Yellow_BP"
D14BlackMod_BPslim_final$Ontolgy     <- "D14Black_BP"
D21YellowMod_BPslim_final$Ontolgy    <- "D21Yellow_BP"


# BIND DATA TO FACET THE PLOTS BY UP AND DOWN REG   ---------------------------------------------------------- #
# bp
BP_Mod <- rbind(D7YellowMod_BPslim_final, D14BlackMod_BPslim_final, D21YellowMod_BPslim_final) # merge the data by binding rows 
BP_Mod$Ont <- "BP" # create a new common clumn to call in the plot 
BP_Mod$Ontolgy <- factor(BP_Mod$Ontolgy, levels = c("D7Yellow_BP", "D14Black_BP", "D21Yellow_BP"))# reorder the facotr level for the facet wrap plot 
BP_Mod_filtered <- BP_Mod %>%  dplyr::filter(Gene.Count > 1) # ommit all with gene counts <1
# mf
MF_Mod <- rbind(D7YellowMod_MFslim_final,D14BlackMod_MFslim_final, D21YellowMod_MFslim_final) # merge the data by binding rows 
MF_Mod$Ont <- "MF" # create a new common clumn to call in the plot 
MF_Mod$Ontolgy <- factor(MF_Mod$Ontolgy, levels = c("D7Yellow_MF", "D14Black_MF", "D21Yellow_MF")) # reorder the facotr level for the facet wrap plot 
MF_Mod_filtered <- MF_Mod %>%  dplyr::filter(Gene.Count > 1) # ommit all with gene counts <1

# PLOTTING --------------------------------------------------------------------------------------------------- #
BP_Plot_PrimEffect_MOD <- ggplot(data = BP_Mod_filtered, aes(x = Ont, y = Gene.Count, fill=forcats::fct_reorder(slim_term,Gene.Count,.desc = TRUE))) + # BP plot
  geom_bar(color = "white",size=1.5,stat="identity") + # call bar and create lines between each slim term
  scale_fill_manual(values =  rev( colorRampPalette(OrangeRed)( (nrow(BP_Mod_filtered)) ) ) )+ # fill with pallette accomodating the max nmber of slim terms
  geom_text(aes(label = forcats::fct_reorder(slim_term,Gene.Count,.desc = TRUE)), colour = "black", position="stack", size = 5, vjust = 1.25) + # add labels
  theme_classic() + # classic theme
  ylim(0, 1000) + 
  theme(legend.position = "none",text = element_text(size=25)) +
  facet_wrap(~Ontolgy) # facet by the upregulated and downregulated datasets

MF_Plot_PrimEffect_MOD <- ggplot(data = MF_Mod_filtered, aes(x = Ont, y = Gene.Count, fill=forcats::fct_reorder(slim_term,Gene.Count,.desc = TRUE))) + # MF plot
  geom_bar(color = "white",size=1.5,stat="identity") + # call bar and create lines between each slim term
  scale_fill_manual(values =  rev( colorRampPalette(GreenBlue)( (nrow(MF_Mod_filtered)) ) ) )+ # fill with pallette accomodating the max nmber of slim terms
  geom_text(aes(label = forcats::fct_reorder(slim_term,Gene.Count,.desc = TRUE)), colour = "black", position="stack", size = 5, vjust = 1.25) + # add labels
  theme_classic() + # classic theme
  ylim(0, 300) +
  theme(legend.position = "none",text = element_text(size=25)) +
  facet_wrap(~Ontolgy) # facet by the upregulated and downregulated datasets

# SAVE PLOTS  ------------------------------------------------------------------------------------------------- #
pdf(paste("Analysis//Output/GO/WGCNA_goseq/GOslim_PrimEff_Moderate_BP.pdf", sep =''), width=15, height=20)
print(ggarrange(BP_Plot_PrimEffect_MOD,         
                plotlist = NULL,
                ncol = 1,
                nrow = 1,
                labels = NULL))
dev.off()     

pdf(paste("Analysis//Output/GO/WGCNA_goseq/GOslim_PrimEff_Moderate_MF.pdf", sep =''), width=15, height=20)
print(ggarrange(MF_Plot_PrimEffect_MOD,         
                plotlist = NULL,
                ncol = 1,
                nrow = 1,
                labels = NULL))
dev.off()     
