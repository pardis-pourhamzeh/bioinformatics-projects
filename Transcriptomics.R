
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
