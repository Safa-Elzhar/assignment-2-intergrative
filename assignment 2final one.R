
library(BiocManager)
library(DESeq2)
library(ggplot2)
library(RColorBrewer)
library(NMF)
•	to read the file 
data = as.matrix(read.csv("~/Desktop/Intergrative/lec 2/assignment 1/bladder_counts.csv", row.names=1, header = T))
View(data)
pheno = read.csv("~/Desktop/Intergrative/lec 2/assignment 1/bladder_sample_sheet.csv", row.names=1)
View(pheno)
•	To sort column
table(pheno$Sample.Type)
•	to have the dimension of data
dim(data)
dim(pheno)
•	to check the normality of data
hist(data, col = "blue", main="Histogram")
hist(log2(data+1), col = "blue", main="Histogram")
boxplot(log2(data[1:5,]+1))
qqnorm(data[1,])
qqline(data[1,])
shapiro.test(data[100,])
•	to check if there is any missing data
is.na(data)
is.null(data)
is.nan(data)
•	to convert data to integers
data=apply(round(data),2,as.integer)
•	to create a deseq dataset object
dds= DESeqDataSetFromMatrix( countData = data , colData = pheno, design = ~ Sample.Type)
View(dds)
dds.run = DESeq(dds)
•	to specify conditions for comparsion
cond1="Primary Tumor"
cond2="Solid Tissue Normal"
res=results(dds.run, contrast = c("Sample.Type",cond1 ,cond2))
•	to remove nulls
res=as.data.frame(res[complete.cases(res), ])
View(res)
•	To choose the statistically significant in differentially expressed genes (DEGs) based
on the p adjusted value less than 0.05 and biological significance base
deseq.deg=res[res$padj < 0.05 & abs(res$log2FoldChange)>1.5,]
•	to save the results as a file
write.csv(as.matrix(deseq.deg),file="deseq.deg.csv", quote=F,row.names=T)

•	to draw DEGs volcano plot
with(res, plot(log2FoldChange, -log10(padj), pch=20, main="normal vs tumor DEGs"))
with(subset(res, padj<.05 & (log2FoldChange)>1.5), points(log2FoldChange, -log10(padj), pch=20, col="green"))
with(subset(res, padj<.05 & (log2FoldChange)< -1.1), points(log2FoldChange, -log10(padj), pch=20, col="purple"))
•	to check and ensure the normalization before heatmap step
dds2 <- estimateSizeFactors(dds)
normalized_counts <- as.data.frame(counts(dds2, normalized=TRUE))
View(normalized_counts)

To extract counts values of DEGs
exp.degs=as.matrix(normalized_counts[rownames(normalized_counts) %in% rownames(deseq.deg), ])
View(exp.degs)
•	Draw heat map
aheatmap(log2(exp.degs+1), annCol =pheno$Sample.Type, col = rev(brewer.pal(9,"RdBu")), main="mRNA Control vs tumor")

heatmap(log2(exp.degs+1))

