library(tidyverse)
library(DESeq2)

#prepare the count data
count <- read.delim("https://www.ebi.ac.uk/gxa/experiments-content/E-MTAB-5244/resources/DifferentialSecondaryDataFiles.RnaSeq/raw-counts")
count <- as.data.frame(count)
view(count)
rownames(count) <- count$Gene.ID
count<- count[,-c(1,2)]

#prepare the sample data
sample <- read.delim("https://www.ebi.ac.uk/gxa/experiments-content/E-MTAB-5244/resources/ExperimentDesignFile.RnaSeq/experiment-design")
sample <- as.data.frame(sample)
rownames(sample) <- sample$Run
sample <- sample[,c("Run","Sample.Characteristic.genotype.")]

colnames(sample)[2] <- "genotype"
sample$genotype <- factor(gsub(" ", "_", sample$genotype))

#check the col and rows of the sample and count data
all(colnames(count) == rownames(sample))

#make the DEseq2 matrix
dd <- DESeqDataSetFromMatrix(countData = count, 
                             colData = sample,
                             design = ~genotype)
#pre filtering
ke <- rowSums(counts(dd)) >= 10
dd <- dd[ke,]

#run the deseq
dds <- DESeq(dd)
#result 

result <- results(dds)