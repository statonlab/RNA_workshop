library(tidyverse)
library(tximport)
library(GenomicFeatures)
library(DESeq2)
library(pheatmap)

# create a workding directory and choose a place to write your results to
setwd("")
dir <- ""

# load the gff3 file, then create a transcript database/dataframe for use with deseq
txdb <- makeTxDbFromGFF("Athaliana_447_Araport11.gene_exons.gff3")
keytypes(txdb)
k <- keys(txdb, keytype = "TXNAME")
txdf = AnnotationDbi::select(txdb, k, "GENEID", "TXNAME")

# load in the metadata
samples <- read_csv("samples_file.csv")
Qfiles <- file.path(dir, samples$quant_file)

# this step imports the count data from salmon
txi <- tximport(files = Qfiles, type = "salmon", tx2gene = txdf)
colnames(txi$counts) <- samples$sample_id
names(txi)
head(txi$counts)
summary(txi)

# now we convert the txi object into a deseq-formatted object
dds_data <- DESeqDataSetFromTximport(txi = txi, colData = samples, design = ~condition)
dds <- DESeq(dds_data)

# plot dispersion
plotDispEsts(dds)

# summarize results
res <- results(dds)
head(res)

# create a contrast with an alpha cutoff (first list item is condition from samples object)
res_sig <- results(dds, alpha = 0.05, contrast = c("condition", "max2", "control"))
max_res_sig <- as.data.frame(res_sig[ which(res_sig$padj < 0.05),])
write_csv(max_res_sig, file=paste0(dir,"salmon_sig_results.csv"))
summary(res_sig)
plotMA(res_sig, ylim=c(-12,12))

# these steps can be used to find individual points on the graph
# after running "identify", click on the plot, then hit "finish" button in top right of plot
idx <- identify(res$baseMean, res$log2FoldChange)
rownames(res)[idx]

# create a plot for a single gene
plotCounts(dds, gene="AT1G53480", intgroup="condition")

# create a PCA plot
rld <- rlog(dds_data, blind = FALSE)
plotPCA(rld, intgroup = c("condition"))
res_lfc <- subset(res_sig, abs(log2FoldChange) > 1) 

# create a heatmap
# first run vst (quickly estimate dispersion trend and apply a variance stabilizing transformation)
vsd <- vst(dds)
genes <- order(res_lfc$log2FoldChange, decreasing=TRUE)[1:50]
pheatmap(assay(vsd)[genes, ], cluster_rows=TRUE, show_rownames=TRUE,
         cluster_cols=TRUE)



