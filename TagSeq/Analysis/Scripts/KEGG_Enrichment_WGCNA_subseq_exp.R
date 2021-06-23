---
  # title: "KEGG analysis" - for subseqent exposures WGGCNA results (without ambient)
  # author: "Samuel Gurr"
  # date: "April 23, 2021"
---
# INFORMATION FOR KEGG IN R FOUND HERE: (http://yulab-smu.top/clusterProfiler-book/chapter6.html#kegg-over-representation-test)

  
  # NOTE: what is the purpose of this script and how does it differ from the full dataset ('_all' version)? 
  # One can think of this analysis as a sanity test relative to the full analysis of controls (ambient exposures)
  #  to hone in on the our biological hypothesis; the goal of this experimental design is to determine whether 
  #  stress history (acclimation) affects gene expression under subseqent stress encounters. 
  #  Thus, to proposerly address this it is critical to let the data do the work and see the broad patterns (completed in other script)
  #  and run a sanity check against target samples core to this question 
  
  # To do this....
  # We will focus on the WGCNA results from the scripts '..._subseq_exp' in which samples from ambient subseqwent exposures were ommitted 
  # WHY? we are interest in the how the initial acclimation ( M and A) affected gene expression patterns 
  # during the subseqent stress. Ambient concurrent conditions (i.e. Second A and Third A) to assess the response specifically under  elevated pCO2



# LOAD PACKAGES
library(KEGGprofile) # BiocManager::install("KEGGprofile")
library(clusterProfiler)
library(KEGGREST)
library(tidyr)
library(stringr)
library(forcats)
library(ggplot2)
library(scales)
library(ape)
library(data.table)
library(tidyverse)
library(fBasics)

# SET WORKING DIRECTORY   ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::;: #
setwd("C:/Users/samjg/Documents/My_Projects/Pgenerosa_TagSeq_Metabolomics/TagSeq/")


# LOAD DATA  :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::;: #


# Crass_gigas_genome_dataframe - what is this? after KEGG enrichment analysis this script with unlist genes ithin enriched pathways and merge IDs to obtain the gene names 
#NOTE: view the Pgenerasa genome - there are no instance of '-like', '-like isoform', '-like protein precursor', 'isoform X', 'precursor' at the end of terms,
Crass_gigas_genome <- keggList("crg") # call the C. gigas genome! - notice the csa terms are rownames!
Crass_gigas_genome_dataframe <- as.data.frame(Crass_gigas_genome) %>%  rownames_to_column() # with will allow us to merge 
colnames(Crass_gigas_genome_dataframe) <- c('sseqid', 'Gene_name') # rename the columns - call the KP term 'sseqid' to merge with the blast data - youll see why donwstream below...
# Crass_gigas_genome_dataframe2 <- tibble::rownames_to_column(Crass_gigas_genome_dataframe, "Cgigas_KEGG_IDs") # tubbel to make th rows a column
# Crass_gigas_genome_dataframe2$Gene_terms  <- sub(" \\-like isoform|-like precursor|-like protein precursor| precursor| isoform X.*", "", Crass_gigas_genome_dataframe$Gene_terms) # none of the Pgenerosa genome has "isoform X1" as the putative term - remove this and merge to see if this increases mapping efficieny
# Crass_gigas_genome_dataframe2$Gene_terms  <- sub("-like$","",Crass_gigas_genome_dataframe$Gene_terms) # remove all occurance of "-like" only when they occur at the end of the string using '$' to call 'only at the end of the string
# hashed out lines are for reducing the name down to core ID (wihtout iosofrm, precorcsors -like, etc. in putative gene names)



# Pgen annotation files and the blastx (using DIAMONd) hits to Cgigas to use downstream in KEGG pathway analysis 
Geoduck_annotation      <- read.delim2(file="C:/Users/samjg/Documents/My_Projects/Pgenerosa_TagSeq_Metabolomics/TagSeq/Seq_details/Panopea-generosa-genes-annotations.txt", header=F)
ref_update_20210602     <- read.gff("C:/Users/samjg/Documents/My_Projects/Pgenerosa_TagSeq_Metabolomics/TagSeq/Seq_details/Panopea-generosa-v1.0.a4.gene.gff3", GFF3 = TRUE) # use library(ape) to import a gff3 as a datatable


crgKEGG_Pgenref_DIAMOND <- read.table(file ="C:/Users/samjg/Documents/My_Projects/Pgenerosa_TagSeq_Metabolomics/TagSeq/HPC_work/Output/crgKEGG_diamond_out.txt", sep = '\t', header = F)
colnames(crgKEGG_Pgenref_DIAMOND) <- c('qseqid','sseqid','pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore') # change the names of the columns if the best blast hits


# WGCNA resilts (all treatments)
d7_WGCNA_all    <- read.csv("Analysis/Output/WGCNA/subseq_hypercapnia/Day7/d7.WGCNA_ModulMembership.csv")
d14_WGCNA_all   <- read.csv("Analysis/Output/WGCNA/subseq_hypercapnia/Day14/d14.WGCNA_ModulMembership.csv")
d21_WGCNA_all   <- read.csv("Analysis/Output/WGCNA/subseq_hypercapnia/Day21/d21.WGCNA_ModulMembership.csv")


table(d21_WGCNA_all$moduleColor)
# DESEq2 results - load data and assign up and downregualted genes!
DESeq2_PrimaryEffects    <- read.csv(file="Analysis/Output/DESeq2/10cpm/DE_PrimaryTreatment.All.csv", sep=',', header=TRUE)  
      # Day 7
d7_PrimaryEffects_up   <- DESeq2_PrimaryEffects %>% dplyr::filter(Day_Trmt %in% 'Day7_Primary_AvM',  up == "TRUE")   %>% dplyr::arrange(desc(log2FoldChange)) %>% dplyr::mutate(DE_dir = 'up')
d7_PrimaryEffects_down <- DESeq2_PrimaryEffects %>% dplyr::filter(Day_Trmt %in% 'Day7_Primary_AvM',  down == "TRUE") %>% dplyr::arrange(desc(abs(log2FoldChange))) %>% dplyr::mutate(DE_dir = 'down')
      # Day 14
d14_PrimaryEffects_up   <- DESeq2_PrimaryEffects %>% dplyr::filter(Day_Trmt %in% 'Day14_Primary_AvM',  up == "TRUE")   %>% dplyr::arrange(desc(log2FoldChange)) %>% dplyr::mutate(DE_dir = 'up')
d14_PrimaryEffects_down <- DESeq2_PrimaryEffects %>% dplyr::filter(Day_Trmt %in% 'Day14_Primary_AvM',  down == "TRUE") %>% dplyr::arrange(desc(abs(log2FoldChange))) %>% dplyr::mutate(DE_dir = 'down')
      # Day 21
d21_PrimaryEffects_up   <- DESeq2_PrimaryEffects %>% dplyr::filter(Day_Trmt %in% 'Day21_Primary_AvM',  up == "TRUE")   %>% dplyr::arrange(desc(log2FoldChange)) %>% dplyr::mutate(DE_dir = 'up')
d21_PrimaryEffects_down <- DESeq2_PrimaryEffects %>% dplyr::filter(Day_Trmt %in% 'Day21_Primary_AvM',  down == "TRUE") %>% dplyr::arrange(desc(abs(log2FoldChange))) %>% dplyr::mutate(DE_dir = 'down')
      # merge for a master file of all primary effects (acclimation 110 days under ambient pCO2 vs. moderate pCO2)  on DE
DESeq2_PrimaryEffects_master <- rbind(d7_PrimaryEffects_up, d7_PrimaryEffects_down, 
                                      d14_PrimaryEffects_up, d14_PrimaryEffects_down,
                                      d21_PrimaryEffects_up, d21_PrimaryEffects_down)







# INITAL SANITY CHECKS AND FILTER FOR BEST BLAST HITS  :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::;: #



#  INITIAL ASSESSMENTS
#  in the annotated geneome relative to the number of genes that hat a blast hit to Cgigas proteins (using blastx in DIAMOND - review the HPC script and output folder....)
# lets see the coverage of the blast hist to the crg KEGG genome
length(na.omit(Geoduck_annotation$V7)) # 14671 total genes with a gene name
length(unique(crgKEGG_Pgenref_DIAMOND$qseqid)) # 18027 of 14671
(( length(unique(crgKEGG_Pgenref_DIAMOND$qseqid)) ) / ( length(na.omit(Geoduck_annotation$V7)) ) ) * 100 # 122.8751 % of annotated unigenes genes had a blast hit on default settings!




###############################################################################################################################################
#::::::::::::::::::::::::::::::::::: MERGE PGENEROSA WITH CGIGAS BY GENE SEQ RELATEDNESS :::::::::::::::::::::::::::::::::::::::::::::::::::;: #
#:::::::::::::::::::::::::::::::::::::: BEST OPTION!!! (next cluster is by gene name) ::::::::::::::::::::::::::::::::::::::::::::::::::::::;: #
# FILTER BEST BLAST HITS - Clean the blast (using DIAMOND) hits ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::;: #
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::;: #
#::::::::::::::::::::::::::::::::::::::::::::::::::::::  crgKEGG_PgenREF_besthits  :::::::::::::::::::::::::::::::::::::::::::::::::::::::::;: #
###############################################################################################################################################



# merge with the cgigas KEGG genome IDs 
crgKEGG_Pgenref_DIAMOND_merge  <- merge(crgKEGG_Pgenref_DIAMOND, Crass_gigas_genome_dataframe, by = 'sseqid')

# PERCENT OF ANNOTATED (GENE NAME AND GO TERM) UNIGENS OF PGENEROSA WITH A BLAST HIT
# lets see how many of the gens in the Pgneorsa geome with putative gene name (unigenes annotated) have a successful blast hit! 
Geoduck_annotation_OM <- na.omit(Geoduck_annotation) # remove all instances of NA (genes without gene name!) 
test <- Geoduck_annotation_OM %>% dplyr::filter(V1 %in% crgKEGG_Pgenref_DIAMOND_merge$qseqid) # filter this dataset by the names in the crg KEGG blast result set
length(test$V1) # 12811 total unigenes from the blast hit in the pool of unigenes iwth gene name 
( length(test$V1) / length(Geoduck_annotation_OM$V1) ) * 100 # 87.32193 % of the Pgen genome (unigenes with gene name) had a successful blast hit to the Cgigas protein database!



# NOTE: we see that there are many genes labeled as "uncharacterized protein ..." 
# In some cases this gene may have the best bti sccore (lowest evalue) and not contribute to the pathway enrichment, whereas another significant gene 
# with putative gene name/annotation can be overlooked! Thus, we will remoe all instanced in  which a blast hit contained an unchracterized protein occurance
# before we call for the blast hit with greatest stttrongt (lowest evalue,highest bitc score)

# ommit all rows that contain the specifi string"uncharacterized protein" in the 'Gene_name' column
crgKEGG_Pgenref_DIAMOND_OM <- crgKEGG_Pgenref_DIAMOND_merge[- grep("uncharacterized", crgKEGG_Pgenref_DIAMOND_merge$Gene_name),]
length(unique(crgKEGG_Pgenref_DIAMOND_OM$qseqid)) # 13726 of 14671 - narrowed down about 5000 hits!

# FILTER BEST BLAST HITS - Clean the blast (using DIAMOND) hits 
# call unique rows (containg unique Pgenerosa IDs and the corresponsding 'crg' KEGG number or 'KO' ID)
# w/ the hightest bit score 
crgKEGG_PgenREF           <- as.data.table(crgKEGG_Pgenref_DIAMOND_OM)
crgKEGG_PgenREF_besthits  <- crgKEGG_PgenREF[,.SD[which.max(bitscore)],by=qseqid] # best hits by highest bit score - do sanity check below..
# sanity check - run by max bitcore and by min evalue - should all say TRUE - use View(booleans)
bybitscore  <- crgKEGG_PgenREF[,.SD[which.max(bitscore)],by=qseqid] # max bitscore
byevalue    <- crgKEGG_PgenREF[,.SD[which.min(evalue)],by=qseqid] # min evalue
# View(bybitscore == byevalue) # use View(booleans). should say TRUE everywhere - this confims the assumption that whether btcore or evalue, you get the same strongest hit

View(crgKEGG_PgenREF_besthits) # file does not contain any ocurance s of 'uncharacterized..." that will not provide insight on pathway enrichment

# WHAT IS THE MEANSD BIT SCORE, EVALUE AND PERCENTID OF THE best hits (note: also NON 'uncharacteried protein or uncharacterized family)
mean(crgKEGG_PgenREF_besthits$pident) # 55.21958
colnames(crgKEGG_PgenREF_besthits) #  "qseqid"    "sseqid"    "pident"    "length"    "mismatch"  "gapopen"   "qstart"    "qend"      "sstart"    "send"      "evalue"    "bitscore"  "Gene_name"
colMeans(crgKEGG_PgenREF_besthits[,c('pident', 'length', 'evalue', 'bitscore')], na.rm=TRUE)  # pident       length       evalue     bitscore   # 5.521958e+01 1.505178e+02 1.126924e-05 1.290717e+02
colStdevs(crgKEGG_PgenREF_besthits[,c('pident', 'length', 'evalue', 'bitscore')], na.rm=TRUE) # pident       length       evalue     bitscore  # 1.776542e+01 1.842089e+02 7.259175e-05 1.105934e+02 

###############################################################################################################################################
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::;: #
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::;: #
#   WGCNA DATA: KEGG ANALYSIS OF ONLY GENES THAT ALIGNED WITH OYSTER KEGG DATABASE (BLASTED AGINST THE OYSTER)
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::;: #
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::;: #
###############################################################################################################################################

# Notes on this pseudo-cluster....
# this cluster was writteni in May  2021 to tackle the KEGG analysis of ALL genes in significant WGCNA modules


# LOAD DATA    :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::;: #

crgKEGG_PgenREF_besthits$qseqid
colnames(crgKEGG_PgenREF_besthits)[1:2] <- c('geneSymbol', 'crg_KO') # reanme the first two columns to merge witht he WGCNA results!

d7_WGCNA_crgKEGGhits  <- merge(d7_WGCNA_all, crgKEGG_PgenREF_besthits[,c(1:2)], by ='geneSymbol')
d14_WGCNA_crgKEGGhits <- merge(d14_WGCNA_all, crgKEGG_PgenREF_besthits[,c(1:2)], by ='geneSymbol')
d21_WGCNA_crgKEGGhits <- merge(d21_WGCNA_all, crgKEGG_PgenREF_besthits[,c(1:2)], by ='geneSymbol')


# Day 7 for loop ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::;;; #
Day7_WGCNA_sigmodules <- as.data.frame(c('blue','brown', 'greenyellow', 'turquoise'))
for (i in 1:nrow(Day7_WGCNA_sigmodules)) {
  # start with loop by calling the row value common with the 'Master_KEGG_BPTerms' data frind from rbind above 
  modColor <- Day7_WGCNA_sigmodules[i,1]
  
  # call the module color in the Day 7 data
  module_without_filter   <- d7_WGCNA_all %>% dplyr::filter(moduleColor %in% modColor)
  genes_per_module        <- length(unique(module_without_filter$geneSymbol)) # nrow(ModuleLoop) # use this for the looped print out 
  annotgenes_per_module   <- length(na.omit(module_without_filter$HGNC)) # nrow(ModuleLoop) # use this for the looped print out 
  
  
  ModuleLoop_blasthit            <- d7_WGCNA_crgKEGGhits %>% dplyr::filter(moduleColor %in% modColor)
  genes_per_module_blasthit      <- na.omit(module_without_filter)  %>% # ommit  genes without gene name annotation 
    dplyr::filter(moduleColor %in% modColor) %>%   # filter for the module loop
    dplyr::filter(geneSymbol %in% d7_WGCNA_crgKEGGhits$geneSymbol) %>%  # call only genes that have a blast hit      
    nrow() # 
  perc_annot_genes_with_blasthit <- ( genes_per_module_blasthit / annotgenes_per_module ) *100
  
  # calculate the percent mapped and print this...
  print(paste("Day7", modColor, " ", 
              genes_per_module, "genes per module", 
              annotgenes_per_module, " annotated; ", 
              genes_per_module_blasthit, " or", perc_annot_genes_with_blasthit,"% annotated genes with blasthit", sep = ' '))
  
  
  # Run KEGG analysis
  KEGG_vector_Pgen_Cgigas   <- as.vector(gsub(".*:","",ModuleLoop_blasthit$crg_KO)) # ommit the 'crg:' before the actual terms
  KEGG_cgigas <- enrichKEGG(gene = KEGG_vector_Pgen_Cgigas, 
                            organism  = 'crg', # 'hsa' is human 'crg' is pacific oyster 
                            pvalueCutoff = 0.05) 
  # if loop to output the KEGG enrichment analysis ONLY if genes were successfully mapped...
  if (  nrow(as.data.frame(head(KEGG_cgigas))) > 0 ) {
    # creat dateframe and write the csv file out 
    df <- as.data.frame(head(KEGG_cgigas))
    rownames(df) <- c()
    KEGGoutput <- as.data.frame(do.call(cbind.data.frame, df))
    KEGGoutput$GeneRatio_2 <- gsub("/"," of ", KEGGoutput$GeneRatio)
    KEGGoutput$Rich_Factor <- (  (as.numeric(sub("/.*", "", KEGGoutput$GeneRatio))) / (as.numeric(sub("/.*", "", KEGGoutput$BgRatio)))  ) 
    
    write.csv(KEGGoutput, file = paste("Analysis/Output/KEGG/subseq_hypercapnia/WGCNA/Day7_",modColor,"_KEGG_allgenes.csv", sep ='')) 
    
    # Plot
    theme_set(theme_classic())
    plot<- KEGGoutput %>%  
      ggplot(aes(x=reorder(Description, Rich_Factor), y= Rich_Factor)) + 
      geom_point( aes(col=qvalue, size=Count)) +   # Draw points
      geom_segment(aes(x=Description, 
                       xend=Description, 
                       y=min(Rich_Factor), 
                       yend=max(Rich_Factor)),  
                   linetype=NA, 
                   size=0) +   # Draw dashed lines
      labs(title="Day 7", 
           x = "Pathway",
           y = "Rich Factor",
           subtitle=paste("WGCNA Module:", modColor, sep =' ')) +
      coord_flip()
    pdf(paste("Analysis/Output/KEGG/subseq_hypercapnia/WGCNA/Day7_",modColor,"_RichFactorPlot.pdf", sep =''), width=5, height=6)
    print(plot)
    dev.off()
    
    
    # stringsplit and unnest for a data set of genes and IDs associated with each pathway 
    df_2 <- as.data.frame(KEGG_cgigas)[c(1:2,8:9)]
    df_2$gene_IDs <- as.vector(strsplit(as.character(df_2$geneID), "/"))
    colnames(df_2) <- c("Cgigas_PathwayID", "Pathway_Description", "Pathway_gene_list", "Pathway_gene_count", "gene_IDs")
    df_3 <- unnest(df_2, gene_IDs)
    df_3$Cgigas_KEGG_IDs <- paste("crg:", df_3$gene_IDs, sep='')
    Crass_gigas_ref <- Crass_gigas_genome_dataframe %>% mutate(Cgigas_KEGG_IDs = Crass_gigas_genome_dataframe$sseqid) %>% select(c('Cgigas_KEGG_IDs','Gene_name'))
    df_final <- merge(df_3, Crass_gigas_ref, by='Cgigas_KEGG_IDs')
    write.csv(df_final, file = paste("Analysis/Output/KEGG/subseq_hypercapnia/WGCNA/Day7_",modColor,"_KEGG_allgenes_unlisted.csv", sep ='')) 
    
  } else {}
  
  print(paste("Finished! Day7 module = ", modColor, sep = " "))
}





# Day 14 for loop ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::;;; #
Day14_WGCNA_sigmodules <- as.data.frame(c('pink','blue'))
for (i in 1:nrow(Day14_WGCNA_sigmodules)) {
  # start with loop by calling the row value common with the 'Master_KEGG_BPTerms' data frind from rbind above 
  modColor <- Day14_WGCNA_sigmodules[i,1]
  
  # call the module color in the Day 7 data
  module_without_filter   <- d14_WGCNA_all %>% dplyr::filter(moduleColor %in% modColor)
  genes_per_module        <- length(unique(module_without_filter$geneSymbol)) # nrow(ModuleLoop) # use this for the looped print out 
  annotgenes_per_module   <- length(na.omit(module_without_filter$HGNC)) # nrow(ModuleLoop) # use this for the looped print out 
  
  
  ModuleLoop_blasthit            <- d14_WGCNA_crgKEGGhits %>% dplyr::filter(moduleColor %in% modColor)
  genes_per_module_blasthit      <- na.omit(module_without_filter)  %>% # ommit  genes without gene name annotation 
    dplyr::filter(moduleColor %in% modColor) %>%   # filter for the module loop
    dplyr::filter(geneSymbol %in% d14_WGCNA_crgKEGGhits$geneSymbol) %>%  # call only genes that have a blast hit      
    nrow() # 
  perc_annot_genes_with_blasthit <- ( genes_per_module_blasthit / annotgenes_per_module ) *100
  
  # calculate the percent mapped and print this...
  print(paste("Day14", modColor, " ", 
              genes_per_module, "genes per module", 
              annotgenes_per_module, " annotated; ", 
              genes_per_module_blasthit, " or", perc_annot_genes_with_blasthit,"% annotated genes with blasthit", sep = ' '))
  
  
  # Run KEGG analysis
  KEGG_vector_Pgen_Cgigas   <- as.vector(gsub(".*:","",ModuleLoop_blasthit$crg_KO)) # ommit the 'crg:' before the actual terms
  KEGG_cgigas <- enrichKEGG(gene = KEGG_vector_Pgen_Cgigas, 
                            organism  = 'crg', # 'hsa' is human 'crg' is pacific oyster 
                            pvalueCutoff = 0.05)
  
  # if lloop to output the KEGG enrichment analysis ONLY if genes were successfully mapped...
  if (  nrow(as.data.frame(head(KEGG_cgigas))) > 0 ) {
    # creat dateframe and write the csv file out 
    df <- as.data.frame(head(KEGG_cgigas))
    rownames(df) <- c()
    KEGGoutput <- as.data.frame(do.call(cbind.data.frame, df))
    KEGGoutput$GeneRatio_2 <- gsub("/"," of ", KEGGoutput$GeneRatio)
    KEGGoutput$Rich_Factor <- (  (as.numeric(sub("/.*", "", KEGGoutput$GeneRatio))) / (as.numeric(sub("/.*", "", KEGGoutput$BgRatio)))  ) 
    
    write.csv(KEGGoutput, file = paste("Analysis/Output/KEGG/subseq_hypercapnia/WGCNA/Day14_",modColor,"_KEGG_allgenes.csv", sep ='')) 
    
    # Plot
    theme_set(theme_classic())
    plot<- KEGGoutput %>%  
      ggplot(aes(x=reorder(Description, Rich_Factor), y= Rich_Factor)) + 
      geom_point( aes(col=qvalue, size=Count)) +   # Draw points
      geom_segment(aes(x=Description, 
                       xend=Description, 
                       y=min(Rich_Factor), 
                       yend=max(Rich_Factor)),  
                   linetype=NA, 
                   size=0) +   # Draw dashed lines
      labs(title="Day 14", 
           x = "Pathway",
           y = "Rich Factor",
           subtitle=paste("WGCNA Module:", modColor, sep =' ')) +
      coord_flip()
    pdf(paste("Analysis/Output/KEGG/subseq_hypercapnia/WGCNA/Day14_",modColor,"_RichFactorPlot.pdf", sep =''), width=5, height=6)
    print(plot)
    dev.off()
    
    
    # stringsplit and unnest for a data set of genes and IDs associated with each pathway 
    df_2 <- as.data.frame(KEGG_cgigas)[c(1:2,8:9)]
    df_2$gene_IDs <- as.vector(strsplit(as.character(df_2$geneID), "/"))
    colnames(df_2) <- c("Cgigas_PathwayID", "Pathway_Description", "Pathway_gene_list", "Pathway_gene_count", "gene_IDs")
    df_3 <- unnest(df_2, gene_IDs)
    df_3$Cgigas_KEGG_IDs <- paste("crg:", df_3$gene_IDs, sep='')
    Crass_gigas_ref <- Crass_gigas_genome_dataframe %>% mutate(Cgigas_KEGG_IDs = Crass_gigas_genome_dataframe$sseqid) %>% select(c('Cgigas_KEGG_IDs','Gene_name'))
    df_final <- merge(df_3, Crass_gigas_ref, by='Cgigas_KEGG_IDs')
    write.csv(df_final, file = paste("Analysis/Output/KEGG/subseq_hypercapnia/WGCNA/Day14_",modColor,"_KEGG_allgenes_unlisted.csv", sep ='')) 
    
  } else {}
  
  print(paste("Finished! Day14 module = ", modColor, sep = " "))
}




# Day 21 for loop ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::;;; #
Day21_WGCNA_sigmodules <- as.data.frame(c('yellow','magenta', 'blue', 'pink'))
for (i in 1:nrow(Day21_WGCNA_sigmodules)) {
  # start with loop by calling the row value common with the 'Master_KEGG_BPTerms' data frind from rbind above 
  modColor <- Day21_WGCNA_sigmodules[i,1]
  
  # call the module color in the Day 7 data
  module_without_filter   <- d21_WGCNA_all %>% dplyr::filter(moduleColor %in% modColor)
  genes_per_module        <- length(unique(module_without_filter$geneSymbol)) # nrow(ModuleLoop) # use this for the looped print out 
  annotgenes_per_module   <- length(na.omit(module_without_filter$HGNC)) # nrow(ModuleLoop) # use this for the looped print out 
  
  
  ModuleLoop_blasthit            <- d21_WGCNA_crgKEGGhits %>% dplyr::filter(moduleColor %in% modColor)
  genes_per_module_blasthit      <- na.omit(module_without_filter)  %>% # ommit  genes without gene name annotation 
    dplyr::filter(moduleColor %in% modColor) %>%   # filter for the module loop
    dplyr::filter(geneSymbol %in% d21_WGCNA_crgKEGGhits$geneSymbol) %>%  # call only genes that have a blast hit      
    nrow() # 
  perc_annot_genes_with_blasthit <- ( genes_per_module_blasthit / annotgenes_per_module ) *100
  
  # calculate the percent mapped and print this...
  print(paste("Day21", modColor, " ", 
              genes_per_module, "genes per module", 
              annotgenes_per_module, " annotated; ", 
              genes_per_module_blasthit, " or", perc_annot_genes_with_blasthit,"% annotated genes with blasthit", sep = ' '))
  
  
  # Run KEGG analysis
  KEGG_vector_Pgen_Cgigas   <- as.vector(gsub(".*:","",ModuleLoop_blasthit$crg_KO)) # ommit the 'crg:' before the actual terms
  KEGG_cgigas <- enrichKEGG(gene = KEGG_vector_Pgen_Cgigas, 
                            organism  = 'crg', # 'hsa' is human 'crg' is pacific oyster 
                            pvalueCutoff = 0.05)
  
  # if lloop to output the KEGG enrichment analysis ONLY if genes were successfully mapped...
  if (  nrow(as.data.frame(head(KEGG_cgigas))) > 0 ) {
    # creat dateframe and write the csv file out 
    df <- as.data.frame(head(KEGG_cgigas))
    rownames(df) <- c()
    KEGGoutput <- as.data.frame(do.call(cbind.data.frame, df))
    KEGGoutput$GeneRatio_2 <- gsub("/"," of ", KEGGoutput$GeneRatio)
    KEGGoutput$Rich_Factor <- (  (as.numeric(sub("/.*", "", KEGGoutput$GeneRatio))) / (as.numeric(sub("/.*", "", KEGGoutput$BgRatio)))  ) 
    
    
    
    write.csv(KEGGoutput, file = paste("Analysis/Output/KEGG/subseq_hypercapnia/WGCNA/Day21_",modColor,"_KEGG_allgenes.csv", sep ='')) 
    
    
    
    # Plot
    theme_set(theme_classic())
    
    plot<- KEGGoutput %>%  
      ggplot(aes(x=reorder(Description, Rich_Factor), y= Rich_Factor)) + 
      geom_point( aes(col=qvalue, size=Count)) +   # Draw points
      geom_segment(aes(x=Description, 
                       xend=Description, 
                       y=min(Rich_Factor), 
                       yend=max(Rich_Factor)),  
                   linetype=NA, 
                   size=0) +   # Draw dashed lines
      labs(title="Day 21", 
           x = "Pathway",
           y = "Rich Factor",
           subtitle=paste("WGCNA Module:", modColor, sep =' ')) +
      coord_flip()
    
    
    pdf(paste("Analysis/Output/KEGG/subseq_hypercapnia/WGCNA/Day21_",modColor,"_RichFactorPlot.pdf", sep =''), width=5, height=6)
    print(plot)
    dev.off()
    
    
    # stringsplit and unnest for a data set of genes and IDs associated with each pathway 
    df_2 <- as.data.frame(KEGG_cgigas)[c(1:2,8:9)]
    df_2$gene_IDs <- as.vector(strsplit(as.character(df_2$geneID), "/"))
    colnames(df_2) <- c("Cgigas_PathwayID", "Pathway_Description", "Pathway_gene_list", "Pathway_gene_count", "gene_IDs")
    df_3 <- unnest(df_2, gene_IDs)
    df_3$Cgigas_KEGG_IDs <- paste("crg:", df_3$gene_IDs, sep='')
    Crass_gigas_ref <- Crass_gigas_genome_dataframe %>% mutate(Cgigas_KEGG_IDs = Crass_gigas_genome_dataframe$sseqid) %>% select(c('Cgigas_KEGG_IDs','Gene_name'))
    df_final <- merge(df_3, Crass_gigas_ref, by='Cgigas_KEGG_IDs')
    write.csv(df_final, file = paste("Analysis/Output/KEGG/subseq_hypercapnia/WGCNA/Day21_",modColor,"_KEGG_allgenes_unlisted.csv", sep ='')) 
    
  } else {}
  
  print(paste("Finished! Day21 module = ", modColor, sep = " "))
}
























###############################################################################################################################################
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::;: #
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::;: #
#   BETTER PLOTS FOR THE KEGG ENRICHMENT ANALYSIS (rich factor plots) 
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::;: #
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::;: #
###############################################################################################################################################


# LOAD DATA    :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::;: #
# Day7 modules
d7_KEGG_brown               <- read.csv("Analysis/Output/KEGG/subseq_hypercapnia/WGCNA/Day7_brown_KEGG_allgenes.csv")
d7_KEGG_brown$ModDay        <- "Day7_brown"
d7_KEGG_greenyellow         <- read.csv("Analysis/Output/KEGG/subseq_hypercapnia/WGCNA/Day7_greenyellow_KEGG_allgenes.csv")
d7_KEGG_greenyellow$ModDay  <- "Day7_greenyellow"
d7_KEGG_turquoise           <- read.csv("Analysis/Output/KEGG/subseq_hypercapnia/WGCNA/Day7_turquoise_KEGG_allgenes.csv")
d7_KEGG_turquoise$ModDay    <- "Day7_turquoise"


# SET-UP AND PLOT  :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::;: #
KEGG_Day7_mods        <- rbind(d7_KEGG_brown, d7_KEGG_greenyellow, d7_KEGG_turquoise)
KEGG_Day7_mods$ModDay <- factor(KEGG_Day7_mods$ModDay , levels = c("Day7_turquoise", "Day7_brown", "Day7_greenyellow")) # for the correct order of facets in the plot below
KEGG_Day7_mods_WGCNA  <- KEGG_Day7_mods  %>%  
  ggplot(aes(x=reorder(Description, Rich_Factor), y= Rich_Factor)) + 
  geom_point( aes(col=qvalue, size=Count)) +   # Draw points
  scale_colour_gradient(low = "#D55E00", high = "#56B4E9", limits = c(0,0.05)) +
  geom_segment(aes(x=Description, 
                   xend=Description, 
                   y=min(Rich_Factor), 
                   yend=max(Rich_Factor)),  
               linetype=NA, 
               size=0) +   # Draw dashed lines
  labs(title="KEGG pathway enrichment analysis", 
       x = "Pathway",
       y = "Rich Factor",
       subtitle="Day 7 subseqent Exposure Modules (WGCNA)") +
  theme_bw() +
  coord_flip() +
  facet_wrap(~ModDay,scales="free_y")
KEGG_Day7_mods_WGCNA


# Day14 modules
d14_KEGG_pink            <- read.csv("Analysis/Output/KEGG/subseq_hypercapnia/WGCNA/Day14_pink_KEGG_allgenes.csv")
d14_KEGG_pink$ModDay     <- "Day14_pink"
d14_KEGG_blue            <- read.csv("Analysis/Output/KEGG/subseq_hypercapnia/WGCNA/Day14_blue_KEGG_allgenes.csv")
d14_KEGG_blue$ModDay     <- "Day14_blue"

# SET-UP AND PLOT  :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::;: #
KEGG_Day14_mods        <- rbind(d14_KEGG_pink, d14_KEGG_blue)
KEGG_Day14_mods$ModDay <- factor(KEGG_Day14_mods$ModDay , levels = c("Day14_blue", "Day14_pink")) # for the correct order of facets in the plot below
KEGG_Day14_mods_WGCNA <- KEGG_Day14_mods  %>%  
  ggplot(aes(x=reorder(Description, Rich_Factor), y= Rich_Factor)) + 
  geom_point( aes(col=qvalue, size=Count)) +   # Draw points
  scale_colour_gradient(low = "#D55E00", high = "#56B4E9", limits = c(0,0.05)) +
  geom_segment(aes(x=Description, 
                   xend=Description, 
                   y=min(Rich_Factor), 
                   yend=max(Rich_Factor)),  
               linetype=NA, 
               size=0) +   # Draw dashed lines
  labs(title="KEGG pathway enrichment analysis", 
       x = "Pathway",
       y = "Rich Factor",
       subtitle="Moderate-Effect Modules (WGCNA)") +
  theme_bw() +
  coord_flip() +
  facet_wrap(~ModDay,scales="free_y")
KEGG_Day14_mods_WGCNA




# Day21 modules
d21_KEGG_magenta         <- read.csv("Analysis/Output/KEGG/subseq_hypercapnia/WGCNA/Day21_magenta_KEGG_allgenes.csv")
d21_KEGG_magenta$ModDay  <- "Day21_magenta"
d21_KEGG_blue            <- read.csv("Analysis/Output/KEGG/subseq_hypercapnia/WGCNA/Day21_blue_KEGG_allgenes.csv")
d21_KEGG_blue$ModDay     <- "Day21_blue"

# SET-UP AND PLOT  :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::;: #
KEGG_Day21_mods        <- rbind(d21_KEGG_magenta, d21_KEGG_blue)
KEGG_Day21_mods$ModDay <- factor(KEGG_Day21_mods$ModDay , levels = c("Day21_blue", "Day21_magenta")) # for the correct order of facets in the plot below
KEGG_Day21_mods_WGCNA <- KEGG_Day21_mods  %>%  
  ggplot(aes(x=reorder(Description, Rich_Factor), y= Rich_Factor)) + 
  geom_point( aes(col=qvalue, size=Count)) +   # Draw points
  scale_colour_gradient(low = "#D55E00", high = "#56B4E9", limits = c(0,0.05)) +
  geom_segment(aes(x=Description, 
                   xend=Description, 
                   y=min(Rich_Factor), 
                   yend=max(Rich_Factor)),  
               linetype=NA, 
               size=0) +   # Draw dashed lines
  labs(title="KEGG pathway enrichment analysis", 
       x = "Pathway",
       y = "Rich Factor",
       subtitle="Moderate-Effect Modules (WGCNA)") +
  theme_bw() +
  coord_flip() +
  facet_wrap(~ModDay,scales="free_y")
KEGG_Day21_mods_WGCNA


# ALL MODULES for days 7 14 and 21!!!   :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::;: #
KEGG_subseq_eff_mods      <- rbind(d7_KEGG_brown, d7_KEGG_greenyellow, d7_KEGG_turquoise,
                                   d14_KEGG_pink, d14_KEGG_blue, 
                                   d21_KEGG_magenta, d21_KEGG_blue)
KEGG_subseq_eff_mods$ModDay <- factor(KEGG_subseq_eff_mods$ModDay , levels = c("Day7_turquoise", "Day7_brown", "Day7_greenyellow",
                                                                                    "Day14_blue", "Day14_pink",
                                                                                    "Day21_blue", "Day21_magenta")) # for the correct order of facets in the plot below
KEGG_subseq_eff_mods_WGCNA <- KEGG_subseq_eff_mods  %>%  
  ggplot(aes(x=reorder(Description, Rich_Factor), y= Rich_Factor)) + 
  geom_point( aes(col=qvalue, size=Count)) +   # Draw points
  scale_colour_gradient(low = "#D55E00", high = "#56B4E9", limits = c(0,0.05)) +
  geom_segment(aes(x=Description, 
                   xend=Description, 
                   y=min(Rich_Factor), 
                   yend=max(Rich_Factor)),  
               linetype=NA, 
               size=0) +   # Draw dashed lines
  labs(title="KEGG pathway enrichment analysis", 
       x = "Pathway",
       y = "Rich Factor",
       subtitle="Primaryt-Effect Modules (WGCNA)") +
  theme_bw() +
  coord_flip() +
  facet_wrap(~ModDay,scales="free_y")
KEGG_subseq_eff_mods_WGCNA


pdf(paste("Analysis/Output/KEGG/subseq_hypercapnia/WGCNA/AllModules_RichFactor.pdf", sep =''), width=12, height=10)
print(KEGG_subseq_eff_mods_WGCNA)
dev.off()




