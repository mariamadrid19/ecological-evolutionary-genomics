#Script modified from Swaegers, J 2024

###############################################
###############################################
##### Analyzing RNA-seq data with DESeq2 ######
###############################################
###############################################

# script partly based on http://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html and 'Beginner's guide to using the DESeq2 package', Love et al. 2014

# data from:
# Swaegers, J., Spanier, K. I., & Stoks, R. (2020). Genetic compensation rather than 
# genetic assimilation drives the evolution of plasticity in response to mild warming 
# across latitudes in a damselfly. Molecular Ecology. 2020 doi: 10.1111/mec.15676

setwd("/Users/mariamadrid/Documents/Eco-Evo-Genomics/lab9")

# In this script we explore the expression variation between the samples, 
# infer differentially expressed genes (DEGs) from several contrasts
# and visualise the DEG analysis

# ! Attention: in the paper sample replicates were collapsed for analysis purposes further on.
# Here we keep them as seperate samples. Differences in results hence occur.

# install DESeq2
#if (!requireNamespace("BiocManager", quietly = TRUE))
  #install.packages("BiocManager")
#BiocManager::install("DESeq2")
library(DESeq2)
library(ggplot2) 
library(egg)


# read count matrix and metadata
count_data <- read.table("count_data_exercise.txt", row.names = 1, header = T)
str(count_data)
count_data_matrix <- as.matrix(count_data)                                                        
mode(count_data_matrix) <- "integer"    

design <- read.table("design.txt", row.names = 1, header = T)
design$treat <- factor(design$treat)
design$lat <- factor(design$lat)
group2<-factor(design$group2)

# run DESeq2

# this pipleline consists of the estimation of size factors 
# (which control for differences in the library size of the sequencing experiments),
# the estimation of dispersion for each gene (a measure of spread or variability in the data), 
# and fitting a generalized linear model.

dds <- DESeqDataSetFromMatrix(countData = count_data,
                              colData = design,
                              design = ~ group2)
dds <- DESeq(dds)

# Principal component analysis

# First apply a variance stabilizing transformation on dds
vsd=vst(dds)

# PCA
pcaDataDEG <- plotPCA(vsd,intgroup=c("treat","lat"),returnData=TRUE)
percentVar <- round(100 * attr(pcaDataDEG, "percentVar"))

ggplot(pcaDataDEG, aes(x = PC1, y = PC2, color = treat, shape = lat)) +
  geom_point(size =3) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  coord_fixed()+theme_classic()

# Detection of DEGs per contrast

# N24 vs N20
res_N24_vsN20 <- results(dds,contrast = c("group2","N24","N20"))
sig_N24_vsN20 <- res_N24_vsN20[ which(res_N24_vsN20$padj < 0.05 ), ] 

# S24 vs S20
res_S24_vsS20 <- results(dds, contrast = c("group2", "S24", "S20"))
sig_S24_vsS20 <- res_S24_vsS20[ which(res_S24_vsS20$padj < 0.05 ), ] 

# S24 vs N24
res_S24_vsN24 <- results(dds, contrast = c("group2", "S24", "N24"))
sig_S24_vsN24 <- res_S24_vsN24[ which(res_S24_vsN24$padj < 0.05 ), ]

# S20 vs N20
res_S20_vsN20 <- results(dds, contrast = c("group2", "S20", "N20"))
sig_S20_vsN20 <- res_S20_vsN20[ which(res_S20_vsN20$padj < 0.05 ), ]

# explore the result table:
# baseMean is the average of the normalized count values, dividing by size factors, 
# taken over all samples. 
# The column log2FoldChange is the effect size estimate. It tells us how much the contig's expression
# seems to have changed due to the temperature change or latitude. This value is reported on
# a logarithmic scale to base 2: for example, a log2 fold change of 1.5 means that the contig's expression
# is increased by a multiplicative factor of 2^1.5, which is 2.82. 
# upregulation: log.fc > 0
# downregulation: log.fc < 0

# Make a MA-plot

# In DESeq2, the function plotMA shows the log2 fold changes attributable to a given 
# variable over the mean of normalized counts for all the samples in the DESeqDataSet. 
# Points which fall out of the window are plotted as open triangles pointing either up or down.

plot.new()
# Set up the plotting area to show 4 plots in one figure
par(mfrow=c(4, 1), mar=c(4, 4, 2, 1))

# Plot each of the MA plots in sequence
plotMA(res_N24_vsN20, ylim=c(-2,2), xlab="", main="N24 vs N20")
plotMA(res_S24_vsS20, ylim=c(-2,2), xlab="", main="S24 vs S20")
plotMA(res_S24_vsN24, ylim=c(-2,2), xlab="", main="S24 vs N24")
plotMA(res_S20_vsN20, ylim=c(-2,2), xlab="Mean Expression", main="S20 vs N20")

# Reset plotting parameters to default
par(mfrow=c(1, 1))



# Make a volcano plot
#BiocManager::install('EnhancedVolcano')
library(EnhancedVolcano)

#North 24 vs North 20
plot.new()
par(mfrow=c(2, 2), mar=c(4, 4, 2, 1))
EnhancedVolcano(res_N24_vsN20,
                lab = rownames(res_N24_vsN20),
                x = 'log2FoldChange',
                y = 'padj',
                labSize = 1.5)


#South 24 vs South 20
EnhancedVolcano(res_S24_vsS20,
                lab = rownames(res_S24_vsS20),
                x = 'log2FoldChange',
                y = 'padj',
                labSize = 1.5)

#South 24 vs North 24
EnhancedVolcano(res_S24_vsN24,
                lab = rownames(res_S24_vsN24),
                x = 'log2FoldChange',
                y = 'padj',
                labSize = 1.5)

#South 20 vs North 20
EnhancedVolcano(res_S20_vsN20,
                lab = rownames(res_S20_vsN20),
                x = 'log2FoldChange',
                y = 'padj',
                labSize = 1.5)


# Make a heatmap with all DEGs

#install.packages("pheatmap")
library("pheatmap")

# Combine all unique differentially expressed genes (DEGs) for all contrasts 
# (with unique = each DEG transcript ID occuring only once)

DEGs_N24_vsN20=rownames(sig_N24_vsN20)
DEGs_S24_vsS20=rownames(sig_S24_vsS20)
DEGs_S24_vsN24=rownames(sig_S24_vsN24)
DEGs_S20_vsN20=rownames(sig_S20_vsN20)
allDEGs<-c(DEGs_N24_vsN20,DEGs_S24_vsS20,DEGs_S24_vsN24,DEGs_S20_vsN20)
allDEGsUnique=unique(allDEGs)
select=allDEGsUnique

nt <- normTransform(dds) # defaults to log2(x+1)
log2.norm.counts <- assay(nt)[select,]
df <- as.data.frame(colData(dds)[,c("treat","lat","pop")])

log2.norm.counts.scaled <- t(scale(t(log2.norm.counts))) #scale your input data to standardised scores
 
pheatmap(log2.norm.counts.scaled, cluster_rows=FALSE, show_rownames=FALSE,
         cluster_cols=FALSE, annotation_col=df)

# Venn diagrams
# make a Venn diagram where you combine the DEGs for both temperature contrasts
# and one where you combine the DEGs for both latitude contrasts

# extract the names of the DEGs like this (do also for other contrasts)
write.table(DEGs_N24_vsN20,file="DEGs_N24_vsN20.txt",row.names=FALSE,col.names=FALSE, sep="\t", quote = FALSE)
write.table(DEGs_S24_vsS20,file="DEGs_S24_vsS20.txt",row.names=FALSE,col.names=FALSE, sep="\t", quote = FALSE)
write.table(DEGs_S24_vsN24,file="DEGs_S24_vsN24.txt",row.names=FALSE,col.names=FALSE, sep="\t", quote = FALSE)
write.table(DEGs_S20_vsN20,file="DEGs_S20_vsN20.txt",row.names=FALSE,col.names=FALSE, sep="\t", quote = FALSE)

#https://bioinfogp.cnb.csic.es/tools/venny/


# It can also be useful to examine the counts of reads for a single contig across the groups.
# A simple function for making this plot is plotCounts, which normalizes counts by sequencing depth 
# and adds a pseudocount of 1/2 to allow for log scale plotting
# The counts are grouped by the variables in intgroup, where more than one variable can be specified. 
# Here we specify the gene which had the smallest p value from the results table created above. 
# You can select the gene to plot by rowname or by numeric index.

# example
d1 <- plotCounts(dds, gene=which.min(res_N24_vsN20$padj), intgroup="group2", 
                returnData=TRUE)
plotd1 <- ggplot(d1, aes(x=group2, y=count)) + 
  geom_point(position=position_jitter(w=0.1,h=0)) + 
  scale_y_log10(breaks=c(25,2000,3000))+theme_minimal()+ylab("log10 counts")


d2 <- plotCounts(dds, gene=which.min(res_S24_vsS20$padj), intgroup="group2", 
                 returnData=TRUE)
plotd2 <- ggplot(d2, aes(x=group2, y=count)) + 
  geom_point(position=position_jitter(w=0.1,h=0)) + 
  scale_y_log10(breaks=c(0,400,800))+theme_minimal()+ylab("log10 counts")


d3 <- plotCounts(dds, gene=which.min(res_S24_vsN24$padj), intgroup="group2", 
                 returnData=TRUE)
plotd3 <- ggplot(d3, aes(x=group2, y=count)) + 
  geom_point(position=position_jitter(w=0.1,h=0)) + 
  scale_y_log10(breaks=c(25,100,400))+theme_minimal()+ylab("log10 counts")


d4 <- plotCounts(dds, gene=which.min(res_S20_vsN20$padj), intgroup="group2", 
                 returnData=TRUE)
plotd4 <- ggplot(d4, aes(x=group2, y=count)) + 
  geom_point(position=position_jitter(w=0.1,h=0)) + 
  scale_y_log10(breaks=c(500,2000,5000))+theme_minimal()+ylab("log10 counts")


ggarrange(plotd1+rremove("xlab"), plotd2+rremove("xlab")+rremove("ylab"), plotd3+rremove("xlab"), plotd4+rremove("xlab")+rremove("ylab"), ncol = 2, nrow = 2)
