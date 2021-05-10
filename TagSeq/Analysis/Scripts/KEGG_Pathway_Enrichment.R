---
  # title: "KEGG analysis"
  # author: "Samuel Gurr"
  # date: "April 23, 2021"
---
# INFORMATION FOR KEGG IN R FOUND HERE: (http://yulab-smu.top/clusterProfiler-book/chapter6.html#kegg-over-representation-test)

# LOAD PACKAGES
library(KEGGprofile) # BiocManager::install("KEGGprofile")
library(clusterProfiler)
library(KEGGREST)
library(tidyr)
library(stringr)
library(tidyr)

# SET WORKING DIRECTORY AND LOAD DATA   ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::;: #
setwd("C:/Users/samjg/Documents/My_Projects/Pgenerosa_TagSeq_Metabolomics/TagSeq/")




# Calland prep the reference annotation file...    :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::;: #


Geoduck_annotation <- read.delim2(file="C:/Users/samjg/Documents/My_Projects/Pgenerosa_TagSeq_Metabolomics/TagSeq/Seq_details/Panopea-generosa-genes-annotations.txt", header=F)
annot.condenced <- Geoduck_annotation[,c(1,7)] # load just the PGEN ID and the putative gene terms
annot.condenced$Gene_term    <- sub(" \\(EC.*", "", annot.condenced$V7)  # call the gene term BEFORe the EC is listed
Pgen_reference <- na.omit(annot.condenced[c(1,3)]) # ommit unannotated genes
names(Pgen_reference)  <- c('Pgenerosa_Gene_IDs', 'Gene_terms') # rename the columns to merge with the modules below...



# PREPARE THE KEGG  Crassostra gigas (Pacific oyster) genome 'Crass_gigas_genome_dataframe' will be called in the KEGG analysis for loops::::::::::::::::::;: #



# Crassostrea gigas - prep for merge with genes of interest (i.e. genes in modules with significant GO enrichment)
#NOTE: view the Pgenerasa genome - there are no instance of '-like', '-like isoform', '-like protein precursor', 'isoform X', 'precursor' at the end of terms,
Crass_gigas_genome <- keggList("crg") # call the C. gigas genome! - notice the csa terms are rownames!
Crass_gigas_genome_dataframe <- as.data.frame(Crass_gigas_genome) # with will allow us to merge 
colnames(Crass_gigas_genome_dataframe)[1] <- 'Gene_terms' # rename the terms to bind with the module of interest
Crass_gigas_genome_dataframe <- tibble::rownames_to_column(Crass_gigas_genome_dataframe, "Cgigas_KEGG_IDs") # tubbel to make th rows a column
Crass_gigas_genome_dataframe$Gene_terms  <- sub(" \\-like isoform|-like precursor|-like protein precursor| precursor| isoform X.*", "", Crass_gigas_genome_dataframe$Gene_terms) # none of the Pgenerosa genome has "isoform X1" as the putative term - remove this and merge to see if this increases mapping efficieny
Crass_gigas_genome_dataframe$Gene_terms  <- sub("-like$","",Crass_gigas_genome_dataframe$Gene_terms) # remove all occurance of "-like" only when they occur at the end of the string using '$' to call 'only at the end of the string
head(Crass_gigas_genome_dataframe) # look at the dataframe, ready to try a merge! 


###############################################################################################################################################
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::;: #
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::;: #
#   DESeq2 DATA: KEGG ANALYSIS OF ALL GENES (upregualted adn downregulated) FROM PAIRWISE DIFFERENTIAL GENE EXPRESSION ANALYSIS
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::;: #
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::;: #
###############################################################################################################################################

# DESEq2 summary table - run the a priori against this table of main effects to look for overable here 
DESeq2_PrimaryEffects    <- read.csv(file="Analysis/Output/DESeq2/10cpm/DE_PrimaryTreatment.All.csv", sep=',', header=TRUE)  

d7_PrimaryEffects_up   <- DESeq2_PrimaryEffects %>% dplyr::filter(Day_Trmt %in% 'Day7_Primary_AvM',  up == "TRUE")   %>% dplyr::arrange(desc(log2FoldChange)) %>% dplyr::mutate(DE_dir = 'up')
d7_PrimaryEffects_down <- DESeq2_PrimaryEffects %>% dplyr::filter(Day_Trmt %in% 'Day7_Primary_AvM',  down == "TRUE") %>% dplyr::arrange(desc(abs(log2FoldChange))) %>% dplyr::mutate(DE_dir = 'down')

d14_PrimaryEffects_up   <- DESeq2_PrimaryEffects %>% dplyr::filter(Day_Trmt %in% 'Day14_Primary_AvM',  up == "TRUE")   %>% dplyr::arrange(desc(log2FoldChange)) %>% dplyr::mutate(DE_dir = 'up')
d14_PrimaryEffects_down <- DESeq2_PrimaryEffects %>% dplyr::filter(Day_Trmt %in% 'Day14_Primary_AvM',  down == "TRUE") %>% dplyr::arrange(desc(abs(log2FoldChange))) %>% dplyr::mutate(DE_dir = 'down')

d21_PrimaryEffects_up   <- DESeq2_PrimaryEffects %>% dplyr::filter(Day_Trmt %in% 'Day21_Primary_AvM',  up == "TRUE")   %>% dplyr::arrange(desc(log2FoldChange)) %>% dplyr::mutate(DE_dir = 'up')
d21_PrimaryEffects_down <- DESeq2_PrimaryEffects %>% dplyr::filter(Day_Trmt %in% 'Day21_Primary_AvM',  down == "TRUE") %>% dplyr::arrange(desc(abs(log2FoldChange))) %>% dplyr::mutate(DE_dir = 'down')

DESeq2_PrimaryEffects_master <- rbind(d7_PrimaryEffects_up, d7_PrimaryEffects_down, 
                                       d14_PrimaryEffects_up, d14_PrimaryEffects_down,
                                        d21_PrimaryEffects_up, d21_PrimaryEffects_down)


# For loop ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::;;; #
days <- c('Day7_Primary_AvM','Day7_Primary_AvM', 'Day14_Primary_AvM', 'Day14_Primary_AvM', 'Day21_Primary_AvM', 'Day21_Primary_AvM')
DE_dir <- c('up','down', 'up','down', 'up','down')
PrimEff_loop <- data.frame(days, DE_dir)
PrimEff_loop

for (i in 1:nrow(PrimEff_loop)) {
  day_trmt <- PrimEff_loop[i,1]
  dir      <- PrimEff_loop[i,2]
  
  DEGs_loop <- DESeq2_PrimaryEffects_master %>%  dplyr::filter(Day_Trmt == day_trmt, DE_dir ==dir) # note this dataset was already sorted by descending LFC
  total_DEGs <- nrow(DEGs_loop)
  colnames(DEGs_loop)[2] <- 'Pgenerosa_Gene_IDs'
  # merge with the Pgen reference annotation and custom the gene terms to later merge with Cgigas KEGG database
  DEGs_Pgen_ref <- merge(Pgen_reference, DEGs_loop, by ='Pgenerosa_Gene_IDs') # merge with the Pgen reference (Gene terms before the EC is stated!)
  DEGs_Pgen_ref$Gene_terms    <- sub(" \\(.*", "", DEGs_Pgen_ref$Gene_terms) # removes all of gene term before the parenthesis 
  DEGs_Pgen_ref$Gene_terms    <- sub(" \\[.*", "", DEGs_Pgen_ref$Gene_terms) # calls al of the string BEFORE an open bracket (in some cases the sttring has just an open bracket without the closed 
  DEGs_Pgen_ref$Gene_terms    <- gsub("\\[.*?\\]", "", DEGs_Pgen_ref$Gene_terms) # now removes all bracketted terms 
  DEGs_Pgen_ref$Gene_terms    <- tolower(DEGs_Pgen_ref$Gene_terms) # convert all to lower case to merge with Cgigas
  DEGs_Pgen_ref_unique <- unique(DEGs_Pgen_ref) # call only unique occurances 

  Pgenerosa_gene_calls <- nrow(DEGs_Pgen_ref) # total unique DEGs (either up or down) for the print out 

  # merge and calculate the percent mapping
  # Lets merge with the Cgigas genome! 
  vector                 <- as.vector(DEGs_Pgen_ref$Gene_terms) # call the gene term as a vector 
  vector_with_asserts    <- paste("^", vector, "$", sep='') # assets help to call the exact beginning and end of the gene term in grep1 (exact mtatch!) - otheriwse grep will use an 'in it' like funciton and call any genes that contain the whole term in addition to extra string chracters
  DEGs_Cgigas <- subset(Crass_gigas_genome_dataframe, grepl(paste("^",vector_with_asserts,"$", sep="", collapse= "|"), Gene_terms, ignore.case = TRUE)) # subset the C gigas genome by the common and exact match to gene terms (with sig GO enrichmetn from 'goseq'!) in the module of interest 
  length(unique(DEGs_Cgigas$Gene_terms)) # number of genes that aligned with Cgigas
  DEGs_Cgigass_removeduplicates = DEGs_Cgigas[!duplicated(DEGs_Cgigas$Gene_terms),]
  
  Cgigas_gene_calls <- nrow(DEGs_Cgigass_removeduplicates) # total Pgen gens mapped to the Cgigas KEGG genome  - use this to learn the percent merged with the C.gigas KEGG dataset
  
  # calaculate the percent mapped and print this...
  Percent_mapped <- (Cgigas_gene_calls / Pgenerosa_gene_calls) * 100
  print(paste(day_trmt, " ", dir, ":", Percent_mapped, "percent mapped to Cgigas KEGG of N =", Pgenerosa_gene_calls, "genes (of ",total_DEGs, ")",  sep = ' '))
  
  # Run KEGG analysis
  KEGG_vector_Pgen_Cgigas   <- as.vector(gsub(".*:","",DEGs_Cgigass_removeduplicates$Cgigas_KEGG_IDs)) # ommit the 'crg:' before the actual terms
  
  KEGG_cgigas <- enrichKEGG(gene = KEGG_vector_Pgen_Cgigas, 
                            organism  = 'crg', # 'hsa' is human 'crg' is pacific oyster 
                            pvalueCutoff = 0.05) 
  # if lloop to output the KEGG enrichment analysis ONLY if genes were successfully mapped...
  if (  nrow(as.data.frame(head(KEGG_cgigas))) > 0 ) {
    # creat dateframe and write the csv file out 
    df <- as.data.frame(head(KEGG_cgigas))
    rownames(df) <- c()
    KEGGoutput <- as.data.frame(do.call(cbind.data.frame, df))
    KEGGoutput$GeneRatio <- gsub("/"," of ", KEGGoutput$GeneRatio)
    write.csv(KEGGoutput, file = paste("Analysis/Output/KEGG/DESeq2/",unlist(strsplit(day_trmt, "_"))[1],"_", dir,"regulated_KEGG.csv", sep ='')) 
    
    # stringsplit and unnest for a data set of genes and IDs associated with each pathway 
    df_2 <- as.data.frame(KEGG_cgigas)[c(1:2,8:9)]
    df_2$gene_IDs <- as.vector(strsplit(as.character(df_2$geneID), "/"))
    colnames(df_2) <- c("Cgigas_PathwayID", "Pathway_Description", "Pathway_gene_list", "Pathway_gene_count", "gene_IDs")
    df_3 <- unnest(df_2, gene_IDs)
    df_3$Cgigas_KEGG_IDs <- paste("crg:", df_3$gene_IDs, sep='')
    df_final <- merge(df_3, Crass_gigas_genome_dataframe, by='Cgigas_KEGG_IDs')
    write.csv(df_final, file = paste("Analysis/Output/KEGG/DESeq2/",unlist(strsplit(day_trmt, "_"))[1],"_", dir,"regulated_KEGG_unlisted.csv", sep ='')) 
    
  } else {}
  
  print(paste(unlist(strsplit(day_trmt, "_"))[1], "_",dir,"regulated - COMPLETED!", sep = ""))
}



###############################################################################################################################################
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::;: #
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::;: #
#   WGCNA DATA: KEGG ANALYSIS OF ALL GENES IN SIGNIFICATN MODULES 
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::;: #
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::;: #
###############################################################################################################################################

# Notes on this pseudo-cluster....
# this cluster was writteni in May  2021 to tackle the KEGG analysis of ALL genes in significant WGCNA modules


# LOAD DATA    :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::;: #
# BIOLOGICAL PROCESS  GO  SLIMS
d7_WGCNA_all   <- read.csv("Analysis/Output/WGCNA/Day7/d7.WGCNA_ModulMembership.csv")
d14_WGCNA_all   <- read.csv("Analysis/Output/WGCNA/Day14/d14.WGCNA_ModulMembership.csv")
d21_WGCNA_all  <- read.csv("Analysis/Output/WGCNA/Day21/d21.WGCNA_ModulMembership.csv")


# Day 7 for loop ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::;;; #
Day7_WGCNA_sigmodules <- as.data.frame(c('brown','yellow', 'green'))
for (i in 1:nrow(Day7_WGCNA_sigmodules)) {
  # start with loop by calling the row value common with the 'Master_KEGG_BPTerms' data frind from rbind above 
  modColor <- Day7_WGCNA_sigmodules[i,1]
  
  # call the module color in the Day 7 data
  ModuleLoop <- d7_WGCNA_all %>% dplyr::filter(moduleColor %in% modColor)
  genes_per_module <- nrow(ModuleLoop) # use this for the looped print out 
  # Prep the enriched genes of interest wihtin each module...
  #NOTE: the master file contains PGEN gene IDs separated by ';' - each row of PGEN IDs contains all genes associated with signifcant GO terms following  GO enrichment (in 'goseq')
  # thus, these genes are called as important due to their GO annotation and significant calls  from functional enrichment analysis 
  EnrichedGenes <- as.data.frame(ModuleLoop$geneSymbol) # split the genes for the significant goseq enrichment 
  colnames(EnrichedGenes)[1] <- 'Pgenerosa_Gene_IDs' # rename the single column in this call
  EnrichedGenes_2 <- merge(Pgen_reference, EnrichedGenes, by ='Pgenerosa_Gene_IDs') # merge with the Pgen reference (Gene terms before the EC is stated!)
  EnrichedGenes_2$Gene_terms    <- sub(" \\(.*", "", EnrichedGenes_2$Gene_terms) # removes all of gene term before the parenthesis 
  EnrichedGenes_2$Gene_terms    <- sub(" \\[.*", "", EnrichedGenes_2$Gene_terms) # calls al of the string BEFORE an open bracket (in some cases the sttring has just an open bracket without the closed 
  EnrichedGenes_2$Gene_terms    <- gsub("\\[.*?\\]", "", EnrichedGenes_2$Gene_terms) # now removes all bracketted terms 
  EnrichedGenes_2$Gene_terms    <- tolower(EnrichedGenes_2$Gene_terms) # convert all to lower case to merge with Cgigas
  EnrichedGenes_2 <- unique(EnrichedGenes_2) # call only unique occurances 
  
  Pgenerosa_gene_calls <- nrow(EnrichedGenes_2) # total genes called - use this to learn the percent merged with the C.gigas KEGG dataset
  
  # merge and calculate the percent mapping
  # Lets merge with the Cgigas genome! 
  vector                 <- as.vector(EnrichedGenes_2$Gene_terms) # call the gene term as a vector 
  vector_with_asserts    <- paste("^", vector, "$", sep='') # assets help to call the exact beginning and end of the gene term in grep1 (exact mtatch!) - otheriwse grep will use an 'in it' like funciton and call any genes that contain the whole term in addition to extra string chracters
  EnrichedGenes_Cgigas <- subset(Crass_gigas_genome_dataframe, grepl(paste("^",vector_with_asserts,"$", sep="", collapse= "|"), Gene_terms, ignore.case = TRUE)) # subset the C gigas genome by the common and exact match to gene terms (with sig GO enrichmetn from 'goseq'!) in the module of interest 
  length(unique(EnrichedGenes_Cgigas$Gene_terms)) # 243 calls
  EnrichedGenes_Cgigas_removeduplicates = EnrichedGenes_Cgigas[!duplicated(EnrichedGenes_Cgigas$Gene_terms),]
  
  Cgigas_gene_calls <- nrow(EnrichedGenes_Cgigas_removeduplicates) # total Pgen gens mapped to the Cgigas KEGG genome  - use this to learn the percent merged with the C.gigas KEGG dataset
  
  # calaculate the percent mapped and print this...
  Percent_mapped <- (Cgigas_gene_calls / Pgenerosa_gene_calls) * 100
  print(paste("Day7", modColor, " ", genes_per_module, " genes per module -->", Percent_mapped, "percent mapped to Cgigas KEGG of N =", Pgenerosa_gene_calls, "total genes", sep = ' '))
  
  # Run KEGG analysis
  KEGG_vector_Pgen_Cgigas   <- as.vector(gsub(".*:","",EnrichedGenes_Cgigas_removeduplicates$Cgigas_KEGG_IDs)) # ommit the 'crg:' before the actual terms
  
  KEGG_cgigas <- enrichKEGG(gene = KEGG_vector_Pgen_Cgigas, 
                            organism  = 'crg', # 'hsa' is human 'crg' is pacific oyster 
                            pvalueCutoff = 0.05) 
  # if lloop to output the KEGG enrichment analysis ONLY if genes were successfully mapped...
  if (  nrow(as.data.frame(head(KEGG_cgigas))) > 0 ) {
    # creat dateframe and write the csv file out 
    df <- as.data.frame(head(KEGG_cgigas))
    rownames(df) <- c()
    KEGGoutput <- as.data.frame(do.call(cbind.data.frame, df))
    KEGGoutput$GeneRatio <- gsub("/"," of ", KEGGoutput$GeneRatio)
    write.csv(KEGGoutput, file = paste("Analysis/Output/KEGG/Day7_",modColor,"_KEGG_allgenes.csv", sep ='')) 
    
    # stringsplit and unnest for a data set of genes and IDs associated with each pathway 
    df_2 <- as.data.frame(KEGG_cgigas)[c(1:2,8:9)]
    df_2$gene_IDs <- as.vector(strsplit(as.character(df_2$geneID), "/"))
    colnames(df_2) <- c("Cgigas_PathwayID", "Pathway_Description", "Pathway_gene_list", "Pathway_gene_count", "gene_IDs")
    df_3 <- unnest(df_2, gene_IDs)
    df_3$Cgigas_KEGG_IDs <- paste("crg:", df_3$gene_IDs, sep='')
    df_final <- merge(df_3, Crass_gigas_genome_dataframe, by='Cgigas_KEGG_IDs')
    write.csv(df_final, file = paste("Analysis/Output/KEGG/Day7_",modColor,"_KEGG_allgenes_unlisted.csv", sep ='')) 
    
  } else {}
  
  print(paste("Finished! Day7 module = ", modColor, sep = " "))
}





# Day 14 for loop ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::;;; #
Day14_WGCNA_sigmodules <- as.data.frame(c('brown','black', 'magenta', 'pink'))
for (i in 1:nrow(Day14_WGCNA_sigmodules)) {
  # start with loop by calling the row value common with the 'Master_KEGG_BPTerms' data frind from rbind above 
  modColor <- Day14_WGCNA_sigmodules[i,1]
  
  # call the module color in the Day 7 data
  ModuleLoop <- d14_WGCNA_all %>% dplyr::filter(moduleColor %in% modColor)
  colnames(ModuleLoop)
  genes_per_module <- nrow(ModuleLoop) # use this for the looped print out 
  # Prep the enriched genes of interest wihtin each module...
  #NOTE: the master file contains PGEN gene IDs separated by ';' - each row of PGEN IDs contains all genes associated with signifcant GO terms following  GO enrichment (in 'goseq')
  # thus, these genes are called as important due to their GO annotation and significant calls  from functional enrichment analysis 
  EnrichedGenes <- as.data.frame(ModuleLoop$geneSymbol) # split the genes for the significant goseq enrichment 
  colnames(EnrichedGenes)[1] <- 'Pgenerosa_Gene_IDs' # rename the single column in this call
  EnrichedGenes_2 <- merge(Pgen_reference, EnrichedGenes, by ='Pgenerosa_Gene_IDs') # merge with the Pgen reference (Gene terms before the EC is stated!)
  EnrichedGenes_2$Gene_terms    <- sub(" \\(.*", "", EnrichedGenes_2$Gene_terms) # removes all of gene term before the parenthesis 
  EnrichedGenes_2$Gene_terms    <- sub(" \\[.*", "", EnrichedGenes_2$Gene_terms) # calls al of the string BEFORE an open bracket (in some cases the sttring has just an open bracket without the closed 
  EnrichedGenes_2$Gene_terms    <- gsub("\\[.*?\\]", "", EnrichedGenes_2$Gene_terms) # now removes all bracketted terms 
  EnrichedGenes_2$Gene_terms    <- tolower(EnrichedGenes_2$Gene_terms) # convert all to lower case to merge with Cgigas
  EnrichedGenes_2 <- unique(EnrichedGenes_2) # call only unique occurances 
  
  Pgenerosa_gene_calls <- nrow(EnrichedGenes_2) # total genes called - use this to learn the percent merged with the C.gigas KEGG dataset
  
  # merge and calculate the percent mapping
  # Lets merge with the Cgigas genome! 
  vector                 <- as.vector(EnrichedGenes_2$Gene_terms) # call the gene term as a vector 
  vector_with_asserts    <- paste("^", vector, "$", sep='') # assets help to call the exact beginning and end of the gene term in grep1 (exact mtatch!) - otheriwse grep will use an 'in it' like funciton and call any genes that contain the whole term in addition to extra string chracters
  EnrichedGenes_Cgigas <- subset(Crass_gigas_genome_dataframe, grepl(paste("^",vector_with_asserts,"$", sep="", collapse= "|"), Gene_terms, ignore.case = TRUE)) # subset the C gigas genome by the common and exact match to gene terms (with sig GO enrichmetn from 'goseq'!) in the module of interest 
  length(unique(EnrichedGenes_Cgigas$Gene_terms)) # 243 calls
  EnrichedGenes_Cgigas_removeduplicates = EnrichedGenes_Cgigas[!duplicated(EnrichedGenes_Cgigas$Gene_terms),]
  
  Cgigas_gene_calls <- nrow(EnrichedGenes_Cgigas_removeduplicates) # total Pgen gens mapped to the Cgigas KEGG genome  - use this to learn the percent merged with the C.gigas KEGG dataset
  
  # calculate the percent mapped and print this...
  Percent_mapped <- (Cgigas_gene_calls / Pgenerosa_gene_calls) * 100
  print(paste("Day14", modColor, " ", genes_per_module, " genes per module -->", Percent_mapped, "percent mapped to Cgigas KEGG of N =", Pgenerosa_gene_calls, "total genes", sep = ' '))
  
  # Run KEGG analysis
  KEGG_vector_Pgen_Cgigas   <- as.vector(gsub(".*:","",EnrichedGenes_Cgigas_removeduplicates$Cgigas_KEGG_IDs)) # ommit the 'crg:' before the actual terms
  
  KEGG_cgigas <- enrichKEGG(gene = KEGG_vector_Pgen_Cgigas, 
                            organism  = 'crg', # 'hsa' is human 'crg' is pacific oyster 
                            pvalueCutoff = 0.05) 
  
  
  # if lloop to output the KEGG enrichment analysis ONLY if genes were successfully mapped...
  if (  nrow(as.data.frame(head(KEGG_cgigas))) > 0 ) {
    # creat dateframe and write the csv file out 
    df <- as.data.frame(head(KEGG_cgigas))
    rownames(df) <- c()
    KEGGoutput <- as.data.frame(do.call(cbind.data.frame, df))
    KEGGoutput$GeneRatio <- gsub("/"," of ", KEGGoutput$GeneRatio)
    write.csv(KEGGoutput, file = paste("Analysis/Output/KEGG/Day14_",modColor,"_KEGG_allgenes.csv", sep ='')) 
    
    # stringsplit and unnest for a data set of genes and IDs associated with each pathway 
    df_2 <- as.data.frame(KEGG_cgigas)[c(1:2,8:9)]
    df_2$gene_IDs <- as.vector(strsplit(as.character(df_2$geneID), "/"))
    colnames(df_2) <- c("Cgigas_PathwayID", "Pathway_Description", "Pathway_gene_list", "Pathway_gene_count", "gene_IDs")
    df_3 <- unnest(df_2, gene_IDs)
    df_3$Cgigas_KEGG_IDs <- paste("crg:", df_3$gene_IDs, sep='')
    df_final <- merge(df_3, Crass_gigas_genome_dataframe, by='Cgigas_KEGG_IDs')
    write.csv(df_final, file = paste("Analysis/Output/KEGG/Day14_",modColor,"_KEGG_allgenes_unlisted.csv", sep ='')) 
    
  } else {}
  
  print(paste("Finished! Day14 module = ", modColor, sep = " "))
}




# Day 21 for loop ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::;;; #
Day21_WGCNA_sigmodules <- as.data.frame(c('blue','magenta', 'yellow', 'black', 'pink', 'red', 'turquoise'))
for (i in 1:nrow(Day21_WGCNA_sigmodules)) {
  # start with loop by calling the row value common with the 'Master_KEGG_BPTerms' data frind from rbind above 
  modColor <- Day21_WGCNA_sigmodules[7,1]
  
  # call the module color in the Day 7 data
  ModuleLoop <- d21_WGCNA_all %>% dplyr::filter(moduleColor %in% modColor)
  genes_per_module <- nrow(ModuleLoop) # use this for the looped print out 
  
  # Prep the enriched genes of interest wihtin each module...
  #NOTE: the master file contains PGEN gene IDs separated by ';' - each row of PGEN IDs contains all genes associated with signifcant GO terms following  GO enrichment (in 'goseq')
  # thus, these genes are called as important due to their GO annotation and significant calls  from functional enrichment analysis 
  EnrichedGenes <- as.data.frame(ModuleLoop$geneSymbol) # split the genes for the significant goseq enrichment 
  colnames(EnrichedGenes)[1] <- 'Pgenerosa_Gene_IDs' # rename the single column in this call
  EnrichedGenes_2 <- merge(Pgen_reference, EnrichedGenes, by ='Pgenerosa_Gene_IDs') # merge with the Pgen reference (Gene terms before the EC is stated!)
  EnrichedGenes_2$Gene_terms    <- sub(" \\(.*", "", EnrichedGenes_2$Gene_terms) # removes all of gene term before the parenthesis 
  EnrichedGenes_2$Gene_terms    <- sub(" \\[.*", "", EnrichedGenes_2$Gene_terms) # calls al of the string BEFORE an open bracket (in some cases the sttring has just an open bracket without the closed 
  EnrichedGenes_2$Gene_terms    <- gsub("\\[.*?\\]", "", EnrichedGenes_2$Gene_terms) # now removes all bracketted terms 
  EnrichedGenes_2$Gene_terms    <- tolower(EnrichedGenes_2$Gene_terms) # convert all to lower case to merge with Cgigas
  EnrichedGenes_2 <- unique(EnrichedGenes_2) # call only unique occurances 
  Pgenerosa_gene_calls <- nrow(EnrichedGenes_2) # total genes called - use this to learn the percent merged with the C.gigas KEGG dataset


  # merge and calculate the percent mapping
  # Lets merge with the Cgigas genome! 
  # EnrichedGenes_2$Gene_terms <- gsub("[[]", "", EnrichedGenes_2$Gene_terms) - need to run this for turquoise - an open bracket [ in one of the genes
  vector                 <- as.vector(EnrichedGenes_2$Gene_terms) # call the gene term as a vector 
  vector_with_asserts    <- paste("^", vector, "$", sep='') # assets help to call the exact beginning and end of the gene term in grep1 (exact mtatch!) - otheriwse grep will use an 'in it' like funciton and call any genes that contain the whole term in addition to extra string chracters

  EnrichedGenes_Cgigas <- subset(Crass_gigas_genome_dataframe, grepl(paste("^",vector_with_asserts,"$", sep="", collapse= "|"), Gene_terms, ignore.case = TRUE)) # subset the C gigas genome by the common and exact match to gene terms (with sig GO enrichmetn from 'goseq'!) in the module of interest 
  length(unique(EnrichedGenes_Cgigas$Gene_terms)) # 243 calls
  EnrichedGenes_Cgigas_removeduplicates = EnrichedGenes_Cgigas[!duplicated(EnrichedGenes_Cgigas$Gene_terms),]
  
  Cgigas_gene_calls <- nrow(EnrichedGenes_Cgigas_removeduplicates) # total Pgen gens mapped to the Cgigas KEGG genome  - use this to learn the percent merged with the C.gigas KEGG dataset

  # calculate the percent mapped and print this...
  Percent_mapped <- (Cgigas_gene_calls / Pgenerosa_gene_calls) * 100
  print(paste("Day21", modColor, " ", genes_per_module, " genes per module -->", Percent_mapped, "percent mapped to Cgigas KEGG of N =", Pgenerosa_gene_calls, "total genes", sep = ' '))
  
  # Run KEGG analysis
  KEGG_vector_Pgen_Cgigas   <- as.vector(gsub(".*:","",EnrichedGenes_Cgigas_removeduplicates$Cgigas_KEGG_IDs)) # ommit the 'crg:' before the actual terms
  
  KEGG_cgigas <- enrichKEGG(gene = KEGG_vector_Pgen_Cgigas, 
                            organism  = 'crg', # 'hsa' is human 'crg' is pacific oyster 
                            pvalueCutoff = 0.05) 
  # if lloop to output the KEGG enrichment analysis ONLY if genes were successfully mapped...
  if (  nrow(as.data.frame(head(KEGG_cgigas))) > 0 ) {
    # creat dateframe and write the csv file out 
    df <- as.data.frame(head(KEGG_cgigas))
    rownames(df) <- c()
    KEGGoutput <- as.data.frame(do.call(cbind.data.frame, df))
    KEGGoutput$GeneRatio <- gsub("/"," of ", KEGGoutput$GeneRatio)
    write.csv(KEGGoutput, file = paste("Analysis/Output/KEGG/Day21_",modColor,"_KEGG_allgenes.csv", sep ='')) 
    
    # stringsplit and unnest for a data set of genes and IDs associated with each pathway 
    df_2 <- as.data.frame(KEGG_cgigas)[c(1:2,8:9)]
    df_2$gene_IDs <- as.vector(strsplit(as.character(df_2$geneID), "/"))
    colnames(df_2) <- c("Cgigas_PathwayID", "Pathway_Description", "Pathway_gene_list", "Pathway_gene_count", "gene_IDs")
    df_3 <- unnest(df_2, gene_IDs)
    df_3$Cgigas_KEGG_IDs <- paste("crg:", df_3$gene_IDs, sep='')
    df_final <- merge(df_3, Crass_gigas_genome_dataframe, by='Cgigas_KEGG_IDs')
    write.csv(df_final, file = paste("Analysis/Output/KEGG/Day21_",modColor,"_KEGG_allgenes_unlisted.csv", sep ='')) 
    
  } else {}
  
  print(paste("Finished! Day21 module = ", modColor, sep = " "))
}



###############################################################################################################################################
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::;: #
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::;: #
#   WGCNA DATA: KEGG ANALYSIS OF FILTERED GENES FROM SIGNIFICANT MODULES  (only those assocaiated with significant GO enrichment - decided against this analysis on  05/2021) 
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::;: #
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::;: #
###############################################################################################################################################

# Notes on this pseudo-cluster....
# this cluster was writteni n April 2021 to tackle the KEGG analysis of genes assocaited with significatn GO enrichement in each module 
# Initial rationale was a filtering system to narrow the genes down to those that are already showing a signiiantly functional enrichment 
# however, as a team we decided against this route to have a less bias view of our data and simple run KEGG against ALL module genes (above) 
# feel free to run this if you are interested.. this script is useful in unlisting a row of delimited gene IDs into a column to merge with KEGG genome data (e.g. oyster) 
# it also writes a for loop that prints the genes called and enriched pathway statistics for each significant module 

# because this is writetn following goseq and filtering only genes of enriched GO terms... the loop is written fro both biological process and molecular function terms!

# LOAD DATA    :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::;: #
# BIOLOGICAL PROCESS  GO  SLIMS
d7_slimBP_brownModule   <- read.csv("Analysis/Output/GO/WGCNA_goseq/Day7/GOslim_BiolProc_brownModule.csv")
d7_slimBP_greenModule   <- read.csv("Analysis/Output/GO/WGCNA_goseq/Day7/GOslim_BiolProc_greenModule.csv")
d7_slimBP_yellowModule  <- read.csv("Analysis/Output/GO/WGCNA_goseq/Day7/GOslim_BiolProc_yellowModule.csv")
  
d14_slimBP_blackModule   <- read.csv("Analysis/Output/GO/WGCNA_goseq/Day14/GOslim_BiolProc_blackModule.csv")
d14_slimBP_brownModule   <- read.csv("Analysis/Output/GO/WGCNA_goseq/Day14/GOslim_BiolProc_brownModule.csv")
d14_slimBP_magentaModule <- read.csv("Analysis/Output/GO/WGCNA_goseq/Day14/GOslim_BiolProc_magentaModule.csv")
d14_slimBP_pinkModule    <- read.csv("Analysis/Output/GO/WGCNA_goseq/Day14/GOslim_BiolProc_pinkModule.csv")
  
d21_slimBP_blackModule       <- read.csv("Analysis/Output/GO/WGCNA_goseq/Day21/GOslim_BiolProc_blackModule.csv")
d21_slimBP_blueModule        <- read.csv("Analysis/Output/GO/WGCNA_goseq/Day21/GOslim_BiolProc_blueModule.csv")
d21_slimBP_magentaModule     <- read.csv("Analysis/Output/GO/WGCNA_goseq/Day21/GOslim_BiolProc_magentaModule.csv")
d21_slimBP_pinkModule        <- read.csv("Analysis/Output/GO/WGCNA_goseq/Day21/GOslim_BiolProc_pinkModule.csv")
d21_slimBP_redModule         <- read.csv("Analysis/Output/GO/WGCNA_goseq/Day21/GOslim_BiolProc_redModule.csv")
d21_slimBP_yellowModule      <- read.csv("Analysis/Output/GO/WGCNA_goseq/Day21/GOslim_BiolProc_yellowModule.csv")
d21_slimBP_turquoiseModule   <- read.csv("Analysis/Output/GO/WGCNA_goseq/Day21/GOslim_BiolProc_turquoiseModule.csv")
  
# MOLECULAR FUNCTION GO  SLIMS
d7_slimMF_brownModule   <- read.csv("Analysis/Output/GO/WGCNA_goseq/Day7/GOslim_MolFunction_brownModule.csv")
d7_slimMF_greenModule   <- read.csv("Analysis/Output/GO/WGCNA_goseq/Day7/GOslim_MolFunction_greenModule.csv")
d7_slimMF_yellowModule  <- read.csv("Analysis/Output/GO/WGCNA_goseq/Day7/GOslim_MolFunction_yellowModule.csv")
  
d14_slimMF_blackModule   <- read.csv("Analysis/Output/GO/WGCNA_goseq/Day14/GOslim_MolFunction_blackModule.csv")
d14_slimMF_brownModule   <- read.csv("Analysis/Output/GO/WGCNA_goseq/Day14/GOslim_MolFunction_brownModule.csv")
d14_slimMF_magentaModule <- read.csv("Analysis/Output/GO/WGCNA_goseq/Day14/GOslim_MolFunction_magentaModule.csv")
d14_slimMF_pinkModule    <- read.csv("Analysis/Output/GO/WGCNA_goseq/Day14/GOslim_MolFunction_pinkModule.csv")
  
d21_slimMF_blackModule       <- read.csv("Analysis/Output/GO/WGCNA_goseq/Day21/GOslim_MolFunction_blackModule.csv")
d21_slimMF_blueModule        <- read.csv("Analysis/Output/GO/WGCNA_goseq/Day21/GOslim_MolFunction_blueModule.csv")
d21_slimMF_magentaModule     <- read.csv("Analysis/Output/GO/WGCNA_goseq/Day21/GOslim_MolFunction_magentaModule.csv")
d21_slimMF_pinkModule        <- read.csv("Analysis/Output/GO/WGCNA_goseq/Day21/GOslim_MolFunction_pinkModule.csv")
d21_slimMF_redModule         <- read.csv("Analysis/Output/GO/WGCNA_goseq/Day21/GOslim_MolFunction_redModule.csv")
d21_slimMF_yellowModule      <- read.csv("Analysis/Output/GO/WGCNA_goseq/Day21/GOslim_MolFunction_yellowModule.csv")
d21_slimMF_turquoiseModule   <- read.csv("Analysis/Output/GO/WGCNA_goseq/Day21/GOslim_MolFunction_turquoiseModule.csv")



# BIOLOGICAL PPROCESS FOR LOOP FOR KEGG ENRICHMENT :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::;: #



# ESSENTIAL PREP for the for KEGG loop below! 
# Biolgical process terms master file 
Master_KEGG_BPTerms <- rbind(d7_slimBP_brownModule, d7_slimBP_yellowModule, d7_slimBP_greenModule, # Day 7 modules 
                            d14_slimBP_blackModule, d14_slimBP_brownModule, d14_slimBP_magentaModule, d14_slimBP_pinkModule, # Day 14 modules 
                            d21_slimBP_blackModule, d21_slimBP_blueModule, d21_slimBP_pinkModule,d21_slimBP_redModule, d21_slimBP_yellowModule, d21_slimBP_turquoiseModule)  # Day 21 primary effect (note: d21_slimBP_magentaModule == NULL)
BP_daymodule_forloop <- as.data.frame(unique(Master_KEGG_BPTerms$module_day))


# Biological process for loop
for (i in 1:nrow(BP_daymodule_forloop)) {
  # start with loop by calling the row value common with the 'Master_KEGG_BPTerms' data frind from rbind above 
  day_module <- BP_daymodule_forloop[i,1]
  day <- sub("\\_.*", "", day_module) # call the day - help for name of csv in output file
  module <- sub(".*\\_", "", day_module) # call the module - help for name of csv in output file
  
  # call the module and day
  ModuleLoop <- Master_KEGG_BPTerms %>% dplyr::filter(module_day %in% day_module)
  
  # Prep the enriched genes of interest wihtin each module...
  #NOTE: the master file contains PGEN gene IDs separated by ';' - each row of PGEN IDs contains all genes associated with signifcant GO terms following  GO enrichment (in 'goseq')
  # thus, these genes are called as important due to their GO annotation and significant calls  from functional enrichment analysis 
  EnrichedGenes <- as.data.frame(unlist( (strsplit(as.character(ModuleLoop$Gene_IDs), ";")) )) # split the genes for the significant goseq enrichment 
  colnames(EnrichedGenes)[1] <- 'Pgenerosa_Gene_IDs' # rename the single column in this call
  EnrichedGenes_2 <- merge(Pgen_reference, EnrichedGenes, by ='Pgenerosa_Gene_IDs') # merge with the Pgen reference (Gene terms before the EC is stated!)
  EnrichedGenes_2$Gene_terms    <- sub(" \\(.*", "", EnrichedGenes_2$Gene_terms) # removes all of gene term before the parenthesis 
  EnrichedGenes_2$Gene_terms    <- sub(" \\[.*", "", EnrichedGenes_2$Gene_terms) # calls al of the string BEFORE an open bracket (in some cases the sttring has just an open bracket without the closed 
  EnrichedGenes_2$Gene_terms    <- gsub("\\[.*?\\]", "", EnrichedGenes_2$Gene_terms) # now removes all bracketted terms 
  EnrichedGenes_2$Gene_terms    <- tolower(EnrichedGenes_2$Gene_terms) # convert all to lower case to merge with Cgigas
  EnrichedGenes_2 <- unique(EnrichedGenes_2) # call only unique occurances 
  
  Pgenerosa_gene_calls <- length(unique(EnrichedGenes_2$Gene_terms)) # total genes called - use this to learn the percent merged with the C.gigas KEGG dataset
  
  # merge and calculate the percent mapping
  # Lets merge with the Cgigas genome! 
  vector                 <- as.vector(unique(EnrichedGenes_2$Gene_terms)) # call the gene term as a vector 
  vector_with_asserts    <- paste("^", vector, "$", sep='') # assets help to call the exact beginning and end of the gene term in grep1 (exact mtatch!) - otheriwse grep will use an 'in it' like funciton and call any genes that contain the whole term in addition to extra string chracters
  EnrichedGenes_Cgigas <- subset(Crass_gigas_genome_dataframe, grepl(paste("^",vector_with_asserts,"$", sep="", collapse= "|"), Gene_terms, ignore.case = TRUE)) # subset the C gigas genome by the common and exact match to gene terms (with sig GO enrichmetn from 'goseq'!) in the module of interest 
  length(unique(EnrichedGenes_Cgigas$Gene_terms)) # 
  EnrichedGenes_Cgigas_removeduplicates = EnrichedGenes_Cgigas[!duplicated(EnrichedGenes_Cgigas$Gene_terms),]
  
  Cgigas_gene_calls <- nrow(EnrichedGenes_Cgigas_removeduplicates) # total Pgen gens mapped to the Cgigas KEGG genome  - use this to learn the percent merged with the C.gigas KEGG dataset
  
  # Run KEGG analysis
  KEGG_vector_Pgen_Cgigas   <- as.vector(gsub(".*:","",EnrichedGenes_Cgigas_removeduplicates$Cgigas_KEGG_IDs)) # ommit the 'crg:' before the actual terms
  
  KEGG_cgigas <- enrichKEGG(gene = KEGG_vector_Pgen_Cgigas, 
                            organism  = 'crg', # 'hsa' is human 'crg' is pacific oyster 
                            pvalueCutoff = 0.05) 
  
  # calaculate the percent mapped and print this...
  Percent_mapped <- (Cgigas_gene_calls / Pgenerosa_gene_calls) * 100
  print(paste(day_module, "-->", Percent_mapped, "percent mapped to Cgigas KEGG of N =", Pgenerosa_gene_calls, "total genes", sep = ' '))
  print(head(KEGG_cgigas)[c(2:4,6)])
  
  # if lloop to output the KEGG enrichment analysis ONLY if genes were successfully mapped...
  if (  nrow(as.data.frame(head(KEGG_cgigas))) > 0 ) {
    # creat dateframe and write the csv file out 
    df <- as.data.frame(head(KEGG_cgigas))
    rownames(df) <- c()
    KEGGoutput <- as.data.frame(do.call(cbind.data.frame, df))
    KEGGoutput$GeneRatio <- gsub("/"," of ", KEGGoutput$GeneRatio)
    write.csv(KEGGoutput, file = paste("Analysis/Output/GO/WGCNA_goseq/",day,"/KEGG/",day,"_module_",module,"_KEGG_BiologicalProcess.csv", sep ='')) 
  } else {}

  print(paste("Finished!", day_module, sep = " "))
}



# Molecular Function FOR LOOP FOR KEGG ENRICHMENT :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::;: #



# ESSENTIAL PREP for the for KEGG loop below! 
# Biolgical process terms master file 
Master_KEGG_MFTerms <- rbind(d7_slimMF_brownModule, d7_slimMF_yellowModule, d7_slimMF_greenModule, # Day 7 modules 
                             d14_slimMF_blackModule, d14_slimMF_brownModule, d14_slimMF_magentaModule, d14_slimMF_pinkModule, # Day 14 modules 
                             d21_slimMF_blackModule, d21_slimMF_blueModule, d21_slimMF_magentaModule, d21_slimMF_pinkModule,d21_slimMF_redModule, d21_slimMF_yellowModule, d21_slimMF_turquoiseModule)  # Day 21 primary effect 
MF_daymodule_forloop <- as.data.frame(unique(Master_KEGG_MFTerms$module_day))


# Biological process for loop
for (i in 1:nrow(MF_daymodule_forloop)) {
  # start with loop by calling the row value common with the 'Master_KEGG_BPTerms' data frind from rbind above 
  day_module <- MF_daymodule_forloop[i,1]
  day <- sub("\\_.*", "", day_module) # call the day - help for name of csv in output file
  module <- sub(".*\\_", "", day_module) # call the module - help for name of csv in output file
  
  # call the module and day
  ModuleLoop <- Master_KEGG_MFTerms %>% dplyr::filter(module_day %in% day_module)
  
  # Prep the enriched genes of interest wihtin each module...
  #NOTE: the master file contains PGEN gene IDs separated by ';' - each row of PGEN IDs contains all genes associated with signifcant GO terms following  GO enrichment (in 'goseq')
  # thus, these genes are called as important due to their GO annotation and significant calls  from functional enrichment analysis 
  EnrichedGenes <- as.data.frame(unlist( (strsplit(as.character(ModuleLoop$Gene_IDs), ";")) )) # split the genes for the significant goseq enrichment 
  colnames(EnrichedGenes)[1] <- 'Pgenerosa_Gene_IDs' # rename the single column in this call
  EnrichedGenes_2 <- merge(Pgen_reference, EnrichedGenes, by ='Pgenerosa_Gene_IDs') # merge with the Pgen reference (Gene terms before the EC is stated!)
  EnrichedGenes_2$Gene_terms    <- sub(" \\(.*", "", EnrichedGenes_2$Gene_terms) # removes all of gene term before the parenthesis 
  EnrichedGenes_2$Gene_terms    <- sub(" \\[.*", "", EnrichedGenes_2$Gene_terms) # calls al of the string BEFORE an open bracket (in some cases the sttring has just an open bracket without the closed 
  EnrichedGenes_2$Gene_terms    <- gsub("\\[.*?\\]", "", EnrichedGenes_2$Gene_terms) # now removes all bracketted terms 
  EnrichedGenes_2$Gene_terms    <- tolower(EnrichedGenes_2$Gene_terms) # convert all to lower case to merge with Cgigas
  EnrichedGenes_2 <- unique(EnrichedGenes_2) # call only unique occurances 
  
  Pgenerosa_gene_calls <- length(unique(EnrichedGenes_2$Gene_terms)) # total genes called - use this to learn the percent merged with the C.gigas KEGG dataset
  
  # merge and calculate the percent mapping
  # Lets merge with the Cgigas genome! 
  vector                 <- as.vector(unique(EnrichedGenes_2$Gene_terms)) # call the gene term as a vector 
  vector_with_asserts    <- paste("^", vector, "$", sep='') # assets help to call the exact beginning and end of the gene term in grep1 (exact mtatch!) - otheriwse grep will use an 'in it' like funciton and call any genes that contain the whole term in addition to extra string chracters
  EnrichedGenes_Cgigas <- subset(Crass_gigas_genome_dataframe, grepl(paste("^",vector_with_asserts,"$", sep="", collapse= "|"), Gene_terms, ignore.case = TRUE)) # subset the C gigas genome by the common and exact match to gene terms (with sig GO enrichmetn from 'goseq'!) in the module of interest 
  length(unique(EnrichedGenes_Cgigas$Gene_terms)) # 
  EnrichedGenes_Cgigas_removeduplicates = EnrichedGenes_Cgigas[!duplicated(EnrichedGenes_Cgigas$Gene_terms),]
  
  Cgigas_gene_calls <- nrow(EnrichedGenes_Cgigas_removeduplicates) # total Pgen gens mapped to the Cgigas KEGG genome  - use this to learn the percent merged with the C.gigas KEGG dataset
  
  # Run KEGG analysis
  KEGG_vector_Pgen_Cgigas   <- as.vector(gsub(".*:","",EnrichedGenes_Cgigas_removeduplicates$Cgigas_KEGG_IDs)) # ommit the 'crg:' before the actual terms
  
  KEGG_cgigas <- enrichKEGG(gene = KEGG_vector_Pgen_Cgigas, 
                            organism  = 'crg', # 'hsa' is human 'crg' is pacific oyster 
                            pvalueCutoff = 0.05) 
  
  # calaculate the percent mapped and print this...
  Percent_mapped <- (Cgigas_gene_calls / Pgenerosa_gene_calls) * 100
  print(paste(day_module, "-->", Percent_mapped, "percent mapped to Cgigas KEGG of N =", Pgenerosa_gene_calls, "total genes", sep = ' '))
  print(head(KEGG_cgigas)[c(2:4,6)])
  
  # if lloop to output the KEGG enrichment analysis ONLY if genes were successfully mapped...
  if (  nrow(as.data.frame(head(KEGG_cgigas))) > 0 ) {
    # creat dateframe and write the csv file out 
    df <- as.data.frame(head(KEGG_cgigas))
    rownames(df) <- c()
    KEGGoutput <- as.data.frame(do.call(cbind.data.frame, df))
    KEGGoutput$GeneRatio <- gsub("/"," of ", KEGGoutput$GeneRatio)
    write.csv(KEGGoutput, file = paste("Analysis/Output/GO/WGCNA_goseq/",day,"/KEGG/",day,"_module_",module,"_KEGG_MolecularFunction.csv", sep ='')) 
  } else {}
  
  print(paste("Finished!", day_module, sep = " "))
}
























































# 
# # run with out data! 
# KEGG_vector_Pgen_Cgigas   <- as.vector(gsub(".*:","",d7_brown_Cgigas_removeduplicates$Cgigas_KEGG_IDs)) # ommit the 'crg:' before the actual terms
# 
# kk_cgigas <- enrichKEGG(gene = KEGG_vector_Pgen_Cgigas, 
#                  organism  = 'crg', # 'hsa' is human 'crg' is pacific oyster 
#                  pvalueCutoff = 0.05) 
# as.data.frame(head(kk_cgigas))
# 
# 
# # KEGG Module over-representation test
# # mkk <- enrichMKEGG(gene = gene,
# #                    organism = 'hsa') # 'hsa' is human 'crg' is pacific oyster 
# 
# mkk_cgigas <- enrichMKEGG(gene = KEGG_vector_Pgen_Cgigas,
#                    organism = 'crg') # 'hsa' is human 'crg' is pacific oyster 
# 
# head(mkk_cgigas)





















































# KEGG gene set enrichemnt analysis 
kk2 <- gseKEGG(geneList     = geneList,
               organism     = 'hsa', # 'hsa' is human 'crg' is pacific oyster 
               nPerm        = 1000,
               minGSSize    = 120,
               pvalueCutoff = 0.05,
               verbose      = FALSE)
head(kk2)


kk2_cgigas <- gseKEGG(geneList = KEGG_vector_Pgen_Cgigas,
                      organism     = 'crg', # 'hsa' is human 'crg' is pacific oyster 
                      nPerm        = 1000,
                      minGSSize    = 120,
                      pvalueCutoff = 0.05,
                      verbose      = FALSE)
head(kk2_cgigas)



# KEGG overrepresentation test
data(geneList, package="DOSE")
View(geneList)  # list of gene names (row names) and their log fold change (column 1)
gene <- names(geneList)[abs(geneList) > 2] # call only genes with a LFC > 2 
kk <- enrichKEGG(gene = gene, 
                 organism  = 'hsa', # 'hsa' is human 'crg' is pacific oyster 
                 pvalueCutoff = 0.05) 
head(kk)


# KEGG Module Gene Set Enrichment Analysis
mkk2 <- gseMKEGG(geneList = geneList,
                 organism = 'hsa') # 'hsa' is human 'crg' is pacific oyster 

head(mkk2)
