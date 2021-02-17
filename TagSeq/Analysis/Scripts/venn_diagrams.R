#Author: Sam Gurr
#Edit by: Sam Gurr
#Date Last Modified: 20210216

rm(list=ls())

# Load packages and pacage version/date/import/depends info
library(dplyr)
library(VennDiagram)
library("ggVennDiagram")
library(ggvenn)
library(gridExtra)

#set working directory--------------------------------------------------------------------------------------------------
# SET WORKING DIRECTORY AND LOAD DATA
setwd("C:/Users/samjg/Documents/My_Projects/Pgenerosa_TagSeq_Metabolomics/TagSeq/")

# upload data
day0.primaryDE<-read.csv("Analysis/Output/DESeq2/Day0/ Day0.Primary_DESeq2results.csv", header=TRUE, sep=",", na.string="NA", as.is=TRUE) 
day7.primaryDE<-read.csv("Analysis/Output/DESeq2/Day7/ Day7.Primary_DESeq2results.csv", header=TRUE, sep=",", na.string="NA", as.is=TRUE) 
day14.primaryDE<-read.csv("Analysis/Output/DESeq2/Day14/ Day14.Primary_DESeq2results.csv", header=TRUE, sep=",", na.string="NA", as.is=TRUE) 
day21.primaryDE<-read.csv("Analysis/Output/DESeq2/Day21/ Day21.Primary_DESeq2results.csv", header=TRUE, sep=",", na.string="NA", as.is=TRUE) 

# call only sig DEGs - osolate just the gene names - rename column to day
day0.primaryDE_UP <- day0.primaryDE %>% dplyr::filter(padj < 0.05) %>% dplyr::filter(log2FoldChange > 1) %>%  dplyr::select('Gene')
day0.primaryDE_UP$Dir <- "upregulated"
day0.primaryDE_UP$Day <- "Day0"
day0.primaryDE_DOWN <- day0.primaryDE %>% dplyr::filter(padj < 0.05) %>% dplyr::filter(log2FoldChange < -1)  %>%  dplyr::select('Gene')
day0.primaryDE_DOWN$Dir <- "downregulated"
day0.primaryDE_DOWN$Day <- "Day0"

day7.primaryDE_UP <- day7.primaryDE %>% dplyr::filter(padj < 0.05) %>% dplyr::filter(log2FoldChange > 1) %>%  dplyr::select('Gene')
day7.primaryDE_UP$Dir <- "upregulated"
day7.primaryDE_UP$Day <- "Day7"
day7.primaryDE_DOWN <- day7.primaryDE %>% dplyr::filter(padj < 0.05) %>% dplyr::filter(log2FoldChange < -1)  %>%  dplyr::select('Gene')
day7.primaryDE_DOWN$Dir <- "downregulated"
day7.primaryDE_DOWN$Day <- "Day7"

day14.primaryDE_UP <- day14.primaryDE %>% dplyr::filter(padj < 0.05) %>% dplyr::filter(log2FoldChange > 1) %>%  dplyr::select('Gene')
day14.primaryDE_UP$Dir <- "upregulated"
day14.primaryDE_UP$Day <- "Day14"
day14.primaryDE_DOWN <- day14.primaryDE %>% dplyr::filter(padj < 0.05) %>% dplyr::filter(log2FoldChange < -1) %>%  dplyr::select('Gene')
day14.primaryDE_DOWN$Dir <- "downregulated"
day14.primaryDE_DOWN$Day <- "Day14"

day21.primaryDE_UP <- day21.primaryDE %>% dplyr::filter(padj < 0.05) %>% dplyr::filter(log2FoldChange > 1)  %>%  dplyr::select('Gene')
day21.primaryDE_UP$Dir <- "upregulated"
day21.primaryDE_UP$Day <- "Day21"
day21.primaryDE_DOWN <- day21.primaryDE %>% dplyr::filter(padj < 0.05) %>% dplyr::filter(log2FoldChange < -1)  %>%  dplyr::select('Gene')
day21.primaryDE_DOWN$Dir <- "downregulated"
day21.primaryDE_DOWN$Day <- "Day21"

# without Day 0 (just days 7 15 and 21)
UPREG.primary <- list(
  Day7    = day7.primaryDE_UP$Gene, 
  Day14   = day14.primaryDE_UP$Gene, 
  Day21   = day21.primaryDE_UP$Gene
)
DWONREG.primary <- list(
  Day7  = day7.primaryDE_DOWN$Gene, 
  Day14 = day14.primaryDE_DOWN$Gene, 
  Day21 = day21.primaryDE_DOWN$Gene
)
# with Day 0
UPREG.primary_withday0 <- list(
  Day0    = day0.primaryDE_UP$Gene, 
  Day7    = day7.primaryDE_UP$Gene, 
  Day14   = day14.primaryDE_UP$Gene, 
  Day21   = day21.primaryDE_UP$Gene
)
DWONREG.primary_withday0 <- list(
  Day0  = day0.primaryDE_DOWN$Gene, 
  Day7  = day7.primaryDE_DOWN$Gene, 
  Day14 = day14.primaryDE_DOWN$Gene, 
  Day21 = day21.primaryDE_DOWN$Gene
)



# Venn diagram - ggVennDiagram
up.venn1 <- ggVennDiagram(UPREG.primary, 
              label_alpha = 0)
up.venn1 <- up.venn1 + ggtitle("UPREGULATED: Primary treatment")
down.venn1 <- ggVennDiagram(DWONREG.primary, 
              label_alpha = 0)
down.venn1 <- down.venn1 + ggtitle("DOWNREGULATED: Primary treatment")
Venn_primary <- grid.arrange(up.venn1, down.venn1, ncol=2, nrow=1, clip="off")

# with Day0
up.venn1_w0 <- ggVennDiagram(UPREG.primary_withday0, 
                          label_alpha = 0)
up.venn1_w0 <- up.venn1_w0 + ggtitle("UPREGULATED: Primary treatment")
down.venn1_w0 <- ggVennDiagram(DWONREG.primary_withday0, 
                            label_alpha = 0)
down.venn1_w0 <- down.venn1_w0 + ggtitle("DOWNREGULATED: Primary treatment")
grid.arrange(up.venn1_w0, down.venn1_w0, ncol=2, nrow=1, clip="off")





# Venn diagram - ggvenn
up.venn2 <- ggvenn(UPREG.primary, 
                   fill_color = c("#999999", "#E69F00", "#56B4E9", "#009E73"),
                   stroke_size = 0.5, set_name_size = 4)
up.venn2 <- up.venn2 + ggtitle("UPREGULATED: Primary treatment")
down.venn2 <- ggvenn(DWONREG.primary, 
                     fill_color = c("#999999", "#E69F00", "#56B4E9", "#009E73"),
                     stroke_size = 0.5, set_name_size = 4)
down.venn2 <- down.venn2 + ggtitle("DOWNREGULATED: Primary treatment")
grid.arrange(up.venn2, down.venn2, ncol=2, nrow=1, clip="off")

# with day 0
up.venn2_w0 <- ggvenn(UPREG.primary_withday0, 
                   fill_color = c("white", "#E69F00", "#56B4E9", "#009E73"),
                   stroke_size = 0.5, set_name_size = 4)
up.venn2_w0 <- up.venn2_w0 + ggtitle("UPREGULATED: Primary treatment")
down.venn2_w0 <- ggvenn(DWONREG.primary_withday0, 
                     fill_color = c("white", "#E69F00", "#56B4E9", "#009E73"),
                     stroke_size = 0.5, set_name_size = 4)
down.venn2_w0 <- down.venn2_w0 + ggtitle("DOWNREGULATED: Primary treatment")
Venn_primary <- grid.arrange(up.venn2_w0, down.venn2_w0, ncol=2, nrow=1, clip="off")


# output
pdf("Analysis/Output/DESeq2/Venn_DE.pdf")
grid.arrange(up.venn2_w0, down.venn2_w0, ncol=2, nrow=1, clip="off")
dev.off()
