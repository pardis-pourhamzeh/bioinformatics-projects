
##This file contains exploratory analyses on Gtex tissue data, as processed by Recount3.



if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

BiocManager::install("recount3")


if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
BiocManager::install(c(
  "recount", "GenomicRanges", "limma", "edgeR", "DESeq2",
  "regionReport", "clusterProfiler", "org.Hs.eg.db", "gplots",
  "derfinder", "GenomicState", "bumphunter", "derfinderPlot", "sessioninfo"
))


install.packages("tidyverse")

install.packages('Seurat')
library(Seurat)

devtools::install_github("thomasp85/patchwork")

library(recount3)

#load RDS files in R
rse_brain <- readRDS("rse_brain.RDS")
rse_brain

#but the counts are in “coverage” format, not count of reads.
##To transform them in counts as we need for DE analysis:
assays(rse_brain)$counts <- transform_counts(rse_brain)


#Thus, we have now two tables: “raw_counts” (coverage) and “counts” 
##- the actual read counts as we know them. 
rse_brain

#We might want to add another table containing the counts converted into TPM values
names(colData(rse_brain))

#Some of the column info are derived from the Gtex repository itself (they begin with “gtex”).
##For example, we can see the RIN of each sample:
head(colData(rse_brain)$gtex.smrin)

#And some info about the source of each sample (e.g. which part of the organ, sex and age range of the donor) 
## we use “table” to summarize the annotations present in each of the following fields in each sample.
table(colData(rse_brain)$gtex.smtsd)

table(colData(rse_brain)$gtex.sex)

table(colData(rse_brain)$gtex.age)

# the number of reads in each sample

head(colData(rse_brain)$"recount_qc.star.number_of_input_reads_both")

#and the percentage of uniquely mapped reads 
## the ones employed for computing the read counts:
head(colData(rse_brain)$'recount_qc.star.uniquely_mapped_reads_%_both')

#In both cases, the field name ends with “both” because the sequencing was paired end.
## So reads and % of mapped reads are counted as pairs.
###Plotting the distribution of reads per sample (all reads available, not only the mapped ones):

inputreads <- colData(rse_brain)$"recount_qc.star.number_of_input_reads_both"
boxplot(inputreads)

# Most of the samples are of 50M paired end reads, but there are some strange outliers.
## Let’s see the minimum and maximum:
minlib <- min(inputreads)
maxlib <- max(inputreads)
minlib

#plotting the % of mapped reads
mapped <- colData(rse_brain)$"recount_qc.star.uniquely_mapped_reads_%_both"
boxplot(mapped)

#As we can see, there are a lot of problematic samples with a very low mapping percentage.
##The minumum percentage is obtained for this sample:
min(mapped)
min_map <- which.min(mapped)
colData(rse_brain)[min_map,]
#NOTE:Among the info there are also the % of the reads unmapped, split according to the cause.
##With a result like this, it would be advisable to go back to the original fastq files and perform quality controls on the reads.

#Also the quality of RNA (the RIN integrity number) is overall just above the passing grade (6) for most of the samples:
rinplot <- colData(rse_brain)$gtex.smrin
boxplot(rinplot)


#We can also retrieve the estimated % of reads coming from rRNA genes
head(colData(rse_brain)$gtex.smrrnart)
boxplot(colData(rse_brain)$gtex.smrrnart)

#That is surprisingly correlated with the % of reads coming from mitochondrial genes (another info available in the rse object):
head(colData(rse_brain)$"recount_qc.aligned_reads%.chrm") 
boxplot(colData(rse_brain)$"recount_qc.aligned_reads%.chrm")

plot(colData(rse_brain)$gtex.smrrnart, colData(rse_brain)$"recount_qc.aligned_reads%.chrm")


#The rows contain info related to genes instead:
names(rowData(rse_brain))


#There are different fields all related to the gene name.
##Since we are using here the RefSeq annotation, the name given to the gene by RefSeq annotators is in “Name”,
###while the “official gene symbol” is in “gene_name”. “gbkey” tells us what type of gene each one is supposed to be:
table(rowData(rse_brain)$gbkey)


#Since it is a “ranged” summarized experiment, we have an additional field, the “rowRanges”:
##NOTE:telling us the genomic coordinates of each gene, as well as its length (sum of the lengths of its annotated exons).
rowRanges(rse_brain)


# As we can see, also non-canonical chromosomes were unfortunately included in the annotation. 
##It means that some genes might appear both on the canonical and on the alternative chromosomes.
table(rowRanges(rse_brain)@seqnames)



#To see if genes on the chrM (mitochondrial RNA) are included in the annotation:
mito <- rse_brain[rowRanges(rse_brain)@seqnames == 'chrM']
dim(mito)
mito
# Note:We can see that 24 genes are annotated on the mitochondrial chromosome, with a lot of reads mapping on them.
##They are not 37 because, strangely, the genes included are the tRNA and rRNA genes (24), not the protein coding ones!
rowRanges(mito)


#For choosing the three replicates we have to work on, we have to check the RIN,
##the overall number of reads of each replicate, % of rRNA and % of mapped reads of the corresponding samples.
###My replicates would be in columns 40,41,42.
colData(rse_brain)$gtex.smrin[40]
colData(rse_brain)$gtex.smrin[41]
colData(rse_brain)$gtex.smrin[42]
#NOTE:What is the minimum RIN? we usually go with at least 7, but 6 or higher is usually considered to be “acceptable”.
##So all the three samples pass the threshold.

#How many reads in each replicate (before mapping):
colData(rse_brain)$"recount_qc.star.number_of_input_reads_both"[40] 
colData(rse_brain)$"recount_qc.star.number_of_input_reads_both"[41]
colData(rse_brain)$"recount_qc.star.number_of_input_reads_both"[42]

#NOTE: To be on the safe side, we better discard any sample with less than 30M reads.

#looking at the estimated fraction of rRNA in the three replicates:
colData(rse_brain)$gtex.smrrnart[40]
colData(rse_brain)$gtex.smrrnart[41]
colData(rse_brain)$gtex.smrrnart[42]
#NOTE:This should be very low, never anyway higher than 10% (or 0.1 since here is the fraction to be reported).
##Once again, all three samples pass the threshold


# The percentage of mapped reads:
colData(rse_brain)$"recount_qc.star.uniquely_mapped_reads_%_both"[40]
colData(rse_brain)$"recount_qc.star.uniquely_mapped_reads_%_both"[41]
colData(rse_brain)$"recount_qc.star.uniquely_mapped_reads_%_both"[42]
#NOTE: We want here at least 80% of the reads uniquely mapped, since it is a human sample.
##Once again, all the three samples pass this threshold.



#Once we have selected the three samples/replicates for the DE analysis, we can build a rse and the count table containing only them.
##The count table will be fed to edgeR.

rse_brain_selected <- rse_brain[,c(40,41,42)]

counts_brain_selected <- assays(rse_brain_selected)$counts

#It is better to give to the genes in the count table their “official gene symbol” as name:
rownames(counts_brain_selected) <- rowData(rse_brain)$gene_name

#Finding DE genes in data from Gtex 
library(recount3)
library(recount)
library(edgeR)
rse_brain <- readRDS("rse_brain.RDS")
rse_heart <- readRDS("rse_heart.RDS")
rse_pancreas <- readRDS("rse_pancreas.RDS")
assays(rse_brain)$counts <- transform_counts(rse_brain)
assays(rse_heart)$counts <- transform_counts(rse_heart)
assays(rse_pancreas)$counts <- transform_counts(rse_pancreas)


#We now select the samples (three per tissue)
rse_brain_selected <- rse_brain[,c(75,76,77)]
rse_heart_selected <- rse_heart[,c(75,76,77)]
rse_pancreas_selected <- rse_pancreas[,c(75,76,77)]
counts_brain_selected <- assays(rse_brain_selected)$counts
counts_heart_selected <- assays(rse_heart_selected)$counts
counts_pancreas_selected <- assays(rse_pancreas_selected)$counts


#We build the count table, give more meaningful names to the columns,
##use for row names to the gene symbols, and finally build the DGE object for edgeR.
x <- cbind(counts_brain_selected,counts_heart_selected,counts_pancreas_selected)

colnames(x) <- c("Brain75", "Brain76","Brain77","Heart75", "Heart76","Heart77","Pancreas75","Pancreas76","Pancreas77")

rownames(x) <- rowData(rse_brain_selected)$gene_name

y <- DGEList(counts=x)

#We define how replicates are grouped:
group <- as.factor(c("Brain","Brain","Brain","Heart","Heart","Heart","Pancreas","Pancreas","Pancreas"))

y$samples$group <- group

#We add to the samples info the “quality” information that we employed to select them.
#For this we can add new labels to the “samples” field of the DGE object, like we just did for the groups.
y$samples$rin <- as.factor(c(colData(rse_brain_selected)$gtex.smrin,colData(rse_heart_selected)$gtex.smrin,colData(rse_pancreas_selected)$gtex.smrin))

y$samples$slice <- as.factor(c(colData(rse_brain_selected)$gtex.smtsd,colData(rse_heart_selected)$gtex.smtsd,colData(rse_pancreas_selected)$gtex.smtsd))

y$samples$sex <- as.factor(c(colData(rse_brain_selected)$gtex.sex,colData(rse_heart_selected)$gtex.sex,colData(rse_pancreas_selected)$gtex.sex))

y$samples$age <- as.factor(c(colData(rse_brain_selected)$gtex.age,colData(rse_heart_selected)$gtex.age,colData(rse_pancreas_selected)$gtex.age))

y$samples$rRNA <- as.factor(c(colData(rse_brain_selected)$gtex.smrrnart,colData(rse_heart_selected)$gtex.smrrnart,colData(rse_pancreas_selected)$gtex.smrrnart))

y$samples$mapped <- as.factor(c(colData(rse_brain_selected)$"recount_qc.star.uniquely_mapped_reads_%_both", colData(rse_heart_selected)$"recount_qc.star.uniquely_mapped_reads_%_both",colData(rse_pancreas_selected)$"recount_qc.star.uniquely_mapped_reads_%_both"))

y$samples$chrm <- as.factor(c(colData(rse_brain_selected)$"recount_qc.aligned_reads%.chrm", colData(rse_heart_selected)$"recount_qc.aligned_reads%.chrm",colData(rse_pancreas_selected)$"recount_qc.aligned_reads%.chrm"))
y

#NOTE:If we look at the library sizes, there are some samples (see heart) with a very low number of reads reported (10 million!).
## I’ve just discovered a strange behavior of the “transform_counts” function of Recount3.
##However, let us assume that if the quality parameters are ok, the counts are ok as well.
table(rowSums(y$counts==0)==9)

#A lot of genes with zero counts in all samples.
keep.exprs <- filterByExpr(y, group=group)
y <- y[keep.exprs,, keep.lib.sizes=FALSE]

#NOTE:A lot of genes removed because of zero or low expression in the samples, with slightly more than 20,000 remaining:
dim(y)

#edgeR includes several useful functions for transforming the counts into counts per million, FPKM, TPM etc.
##here we extract and store in a vector the log of the counts per million before normalization with the “cpm” function, and then normalize them:
logcpm_before <- cpm(y, log=TRUE)
y <- calcNormFactors(y, method = "TMM")
y

#Notice that one of the pancreas samples has a very small normalization factor.
##to compare the distribution of count values before and after normalization:
logcpm_after <- cpm(y, log=TRUE)
boxplot(logcpm_before, notch=T, las=2)

boxplot(logcpm_after, notch=T, las=2)
#NOTE:the shift of some of the boxes is quite evident


#Now we must design the linear model.
design <- model.matrix(~0+group, data=y$samples)
colnames(design) <- levels(y$samples$group)
design

#To see if/how the samples/replicates cluster together.
logcpm <- cpm(y, log=TRUE)
plotMDS(logcpm, labels=group)

#Replicates cluster well. One of the heart samples is a bit distant from the other two.
##By trying labelling the points in the plot with different info, I could single out what seem to be the most relevant factors:
plotMDS(logcpm, labels=y$samples$rRNA)

plotMDS(logcpm, labels=y$samples$chrm)

plotMDS(logcpm, labels=y$samples$age)

#That is - the more distant heart replicate differs from the other two for the % of rRNA and of mitochondrial RNA,
##and the age of the donor.
y <- estimateDisp(y, design)
plotBCV(y)


fit <- glmQLFit(y, design)

#heart (top) vs brain (bottom)
qlfHB <- glmQLFTest(fit, contrast=c(-1,1,0))
#pancreas (top) vs brain (bottom)
qlfPB <- glmQLFTest(fit, contrast=c(-1,0,1))
#pancreas (top) vs heart (bottom)
qlfPH <- glmQLFTest(fit, contrast=c(0,-1,1))

#The three comparisons return a data structure with all the info regarding the comparisons themselves,
##including a table with the result of the test for each gene.
##Looking at the numbers (-1/+1), tells which sample will correspond to a negative and a positive variation of expression.
qlfHB

head(qlfHB$table)

#“topTags” extracts the table, sorted by p-value, and adds the FDR:
topTags(qlfHB, n=10,adjust.method = "BH", sort.by = "PValue")

#By default it outputs just the top 10 genes (the n parameter).
resultsHB <- topTags(qlfHB, n = 10000000, adjust.method = "BH", sort.by = "PValue", p.value = 1)

#That can be saved in a file, e.g. for importing it in Excel:
write.table(resultsHB, "resultsHB.txt")

#Finally, a quick way to have an idea on how many DEG genes there are,
##according to different FDR/FC thresholds:
summary(decideTests(qlfHB, p.value=0.05, adjust.method = "BH", lfc=0))

summary(decideTests(qlfHB, p.value=0.05, adjust.method = "BH", lfc=1))

summary(decideTests(qlfHB, p.value=0.01, adjust.method = "BH", lfc=0))

summary(decideTests(qlfHB, p.value=0.01, adjust.method = "BH",lfc=1))

#We can save also the results of the other two comparisons:
resultsPB <- topTags(qlfPB, n = 10000000, adjust.method = "BH", sort.by = "PValue", p.value = 1)

resultsPH <- topTags(qlfPH, n = 10000000, adjust.method = "BH", sort.by = "PValue", p.value = 1)

write.table(resultsPH, "resultsPH.txt")

write.table(resultsPB, "resultsPB.txt")

#Finally, to check whether one of the genes that were DE among these samples are still DE if we consider the complete tissue datasets.
##For example, gene “PCDH10” is on the top of the DE lists, over-expressed in brain with respect to the other two tissues.
## to plot the distribution of expression across the three complete datasets, as TPM:

assays(rse_brain)$TPM <- recount::getTPM(rse_brain)

assays(rse_heart)$TPM <- recount::getTPM(rse_heart)

assays(rse_pancreas)$TPM <- recount::getTPM(rse_pancreas)

which(rowData(rse_brain)$gene_name == "PCDH10")

#In this way we got the row number corresponding to gene “PCDH10”.
boxplot(assays(rse_brain)$TPM[39995,],assays(rse_heart)$TPM[39995,], assays(rse_pancreas)$TPM[39995,], outline=F , names=c("Brain","Heart","Pancreas"))

#Now, the fact that the gene is over-expressed in brain across all samples is obvious.
##to confirm it anyway with a non-parametric test, using the TPM values of the gene in all the samples of the three tissues:
wilcox.test(assays(rse_brain)$TPM[39995,], assays(rse_heart)$TPM[39995,])

#And
wilcox.test(assays(rse_brain)$TPM[39995,], assays(rse_pancreas)$TPM[39995,])

