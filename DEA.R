library(tidyverse)

#Load in raw counts generated from raw sequences using mRatBN7.2 reference genome and Rsubread package from Bioconductor.
#Combine to make one dataset with rownames as gene names 
Counts<-cbind(VEH,THC)

#create colData table; column names are sample number
#Samples 1-4 = "VEH", Samples 11-14 = "THC"
condition<-factor(c(rep("VEH",4), rep("THC", 4)))
coldata<-data.frame(row.names=colnames(Counts), condition)
#Remember that DESeq dataset can only be made from numeric matrix 
Counts<-as.matrix(Counts)
mode(Counts)<-"integer"

#Set up DESeq2 objeCounts#Set up DESeq2 object
library(DESeq2)
dds <-DESeqDataSetFromMatrix(countData=Counts, colData=coldata, design=~condition)

#Filter to remove low count genes 
smallestGroupSize<-4
keep<-rowSums(counts(dds) >=10) >= smallestGroupSize 
dds<-dds[keep,]
summary(keep)
#19217 genes left after filtering 


#Run DESeq with design considering 2 groups (use Wald test first for pair-wise comparisons)
dds<-DESeq(dds, test = "Wald")
res_Wald<-results(dds, independentFiltering=TRUE, alpha=0.05)
#Any significant genes with FDR of 5%?
sig<-subset(res_Wald, padj<0.05)
summary(sig)
#49 significant genes

#Apply variance-stabilizing transformation to account for highly expressed genes being more variable, masking patterns in lowly expressed genes
#Use blind = False modifier to allow vst to account for groups (i.e., not blinded)
vsd<-vst(dds, blind = FALSE)

#PCA plot
plotPCA(vsd, intgroup = "condition")
#Note that this function is using ntop=500 top features by variance 
#Sample distances (Euclidean) and plot as heatmap
distances<-dist(t(assay(vsd)))
heatmap(as.matrix(distances))
#Interpreting: darker = samples more different in expression profile, lighter = more similar

#Use LRT instead of Wald test to run omnibus test ("Does expression differ between the treatment groups?")
dds_lrt<-DESeq(dds, test="LRT", reduced = ~1)
#Results 
res_lrt<-results(dds_lrt, independentFiltering=TRUE, alpha=0.05)
summary(res_lrt)
#Any significant genes with FDR of 5%?
sig<-subset(res_lrt, padj<0.05)
summary(sig)
#46 genes with adjusted p-value of < 0.05 from DEA between THC and VEH groups


#Create lists of padj<0.05 genes for both Wald and LRT DESeq tests, with log fold change values (for Github repo)
sigWald<-subset(res_Wald, padj<0.05)
sigLRT<-subset(res_lrt, padj<0.05)
#Export these results 
sigWald<-as.data.frame(sigWald)
write.csv(sigWald, file="DESeq_sig_genes_Wald_test.csv")
sigLRT<-as.data.frame(sigLRT)
write.csv(sigLRT, file="DESeq_sig_genes_LRT_test.csv")

#For graphing/publication: although Wald test is default for 2 groups, however LRT is more robust
#Use LRT results for publication 
library(ggplot2)
#MA plot using significant genes from LRT results
ggplot(sigLRT, aes(baseMean, log2FoldChange, color = padj < 0.05)) +
  geom_point(alpha = 0.5, size = 1.5) +
  scale_x_log10() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(x = "Mean expression", y = "Log2 fold change") +
  scale_color_manual(values = c("grey6", "red")) +
  theme_classic(base_size = 14)

#MA plot with all data, highlighting significant genes
LRT<-as.data.frame(res_lrt)
ggplot(LRT, aes(baseMean, log2FoldChange, color = padj < 0.05)) +
  geom_point(alpha = 0.5, size = 1.5) +
  scale_x_log10() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(x = "Mean expression", y = "Log2 fold change") +
  scale_color_manual(values = c("grey6", "red")) +
  theme_classic(base_size = 14)

#Heatmap of significant DEGs
#Use normalized counts 
vsd_LRT<-vst(dds_lrt, blind=FALSE)
#Plot top 50 genes by variance
topLRTgenes<-head(order(rowVars(assay(vsd)),decreasing=TRUE), 50)
library(pheatmap)
pheatmap(assay(vsd_LRT)[topLRTgenes, ],
         cluster_rows = TRUE,
         cluster_cols = FALSE,
         show_rownames = TRUE,
         color = colorRampPalette(c("navy","white","firebrick3"))(50))

#Plot top significant genes from LRT analysis 
sigLRTgenes<-rownames(subset(res_lrt, padj<0.05 & !is.na(padj)))
LRT <- assay(vsd_LRT)[sigLRTgenes, ]
pheatmap(LRT,
         cluster_rows = TRUE,
         cluster_cols = FALSE,
         show_rownames = TRUE,
         color = colorRampPalette(c("navy","white","firebrick3"))(50))

#Colours are very similar -> use z-score to scale each gene 
LRT<-t(scale(t(LRT)))
pheatmap(LRT,
         cluster_rows = TRUE,
         cluster_cols = FALSE,
         show_rownames = TRUE,
         color = colorRampPalette(c("navy","white","firebrick3"))(50))
#Much better colour differentiation -> use for publication