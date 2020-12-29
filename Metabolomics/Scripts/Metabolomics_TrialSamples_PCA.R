# Title: Metabolomics_TrialSamples_PCA
# Written by Sam JGurr
# Data updated: 20201211

# load packages:
library("factoextra")
library("ggfortify")
library("dplyr")
library("tidyr")
library("tidyverse")
# set wd
setwd("C:/Users/samjg/Documents/My_Projects/Pgenerosa_RNAseq_Metabolomics/")
# load data
TrialNegPeaks <- read.csv(file="Data/Metabolomics/20200911_TrialSamples_Peaks_Neg.csv", sep =",") 
TrialPosPeaks <- read.csv(file="Data/Metabolomics/20200911_TrialSamples_Peaks_Pos.csv", sep =",") 
SampleReference <- read.csv(file="Data/Metabolomics/References.csv", sep =",") 
# choose columns
list(colnames(TrialNegPeaks))
TrialNegPeaks_2 <- TrialNegPeaks[,c(9,17:25)]
TrialPosPeaks_2 <- TrialPosPeaks[,c(9,17:25)]
# merge with reference
TrialNegPeaks_long <- TrialNegPeaks_2 %>% tidyr::pivot_longer(!compound, names_to = "Metabolomics_ID", values_to = "values")
TrialNegPeaks_long_refs <- merge(TrialNegPeaks_long, SampleReference, by='Metabolomics_ID')
TrialNegPeaks_wide_refs <- TrialNegPeaks_long_refs %>% pivot_wider(names_from = compound, values_from = values)


TrialPosPeaks_long <- TrialPosPeaks_2 %>% tidyr::pivot_longer(!compound, names_to = "Metabolomics_ID", values_to = "values")
TrialPosPeaks_long_refs <- merge(TrialPosPeaks_long, SampleReference, by='Metabolomics_ID')
TrialPosPeaks_wide_refs <- TrialPosPeaks_long_refs %>% pivot_wider(names_from = compound, values_from = values)
# PCA
list(colnames(TrialNegPeaks_wide_refs)) # negative peaks
list(colnames(TrialPosPeaks_wide_refs)) # positive peaks
# NEGATIVE PEAKS
TrialNegPeaks_treatprimary <- TrialNegPeaks_wide_refs[,c(1,8:106)]
colnames(TrialNegPeaks_treatprimary)
TrialNegPeaks_treatprimary2 <- TrialNegPeaks_treatprimary %>% remove_rownames %>% column_to_rownames(var="Metabolomics_ID")
res.pca <- prcomp(TrialNegPeaks_treatprimary2, scale = TRUE)
fviz_eig(res.pca)
fviz_pca_ind(res.pca,
             col.ind = "cos2", # Color by the quality of representation
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
)
fviz_pca_var(res.pca,
             col.var = "contrib", # Color by contributions to the PC
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
)
fviz_pca_biplot(res.pca, repel = TRUE,
                col.var = "#2E9FDF", # Variables color
                col.ind = "#696969"  # Individuals color
)

test <- TrialNegPeaks_wide_refs[c(4:9),c(6,8:106)] # only days 4 and 11 and calling treatment
tapply(df$speed, df$dive, mean)

  
test2 <- test[, lapply(mean), by=Treatment_primary]


S_negpeaks <- cov(TrialNegPeaks_wide_refs[,8:106])  # Find the covariance matrix S of the data. The grouping column is not included.
sum(diag(S_negpeaks)) # The total variance is defined as: ???j=1ksjj Which is also equal to the sum of the eigenvalues of S
s.eigen_negpeaks <- eigen(S_negpeaks) # Compute the eigenvalues and corresponding eigenvectors of S.
s.eigen_negpeaks
# The eigenvectors represent the principal components of S. 
# The eigenvalues of S are used to find the proportion of the total variance explained by the components.
for (s in s.eigen_negpeaks$values) {
  print(s / sum(s.eigen_negpeaks$values))
}
# 0.8941043 PCA1 
# 0.09588464 PCA 2 
# The first two principal components account for 90% of the total variance. 
# A scree graph of the eigenvalues can be plotted to visualize the proportion of 
# variance explained by each subsequential eigenvalue.
# plotting PCA and summary status
plot(s.eigen_negpeaks$values, xlab = 'Eigenvalue Number', ylab = 'Eigenvalue Size', main = 'Scree Graph')
lines(s.eigen_negpeaks$values) # you can see a high drop off from PCA 1 to PCA 2
s.eigen_negpeaks$vectors
negpeaks.pca <- prcomp(TrialNegPeaks_wide_refs[,8:106])
negpeaks.pca
summary(negpeaks.pca)
negpeaks.plot <- autoplot(negpeaks.pca, data = TrialNegPeaks_wide_refs, colour = 'compound')
negpeaks.plot
#Interpretation of principal components is still a heavily researched topic in statistics, and although the components
# may be readily interpreted in most settings, this is not always the case (Joliffe, 2002).
#One method of interpretation of the principal components is to calculate the correlation 
# between the original data and the component. The autoplot() function also generates a nice 
# data table with the original variables and the calculated PCs, which we will use here to find the correlations.
# First, compute the correlations between the data and the calculated components of the covariance matrix S.
comps <- negpeaks.plot$data[,c(1:2)]
cor(comps[,3:8], comps[,c(1:2,10)])