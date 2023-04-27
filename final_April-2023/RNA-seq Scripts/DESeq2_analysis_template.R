getwd()

# follow this tutorial:
# http://bioconductor.org/help/course-materials/2016/CSAMA/lab-3-rnaseq/rnaseq_gene_CSAMA2016.pdf

# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("GenomicFeatures")
# code to import all written RNA-Seq functions to use 
source("F:/Connor/20220505_NextSeq_RNA-seq/scripts/RNA-seq_functions.R")

library('Rsamtools')
library('GenomicFeatures')
library('GenomicAlignments')
library('BiocParallel')
library('DESeq2')
library('PCAtools')
library('ggplot2')
library('pheatmap')
library('AnnotationDbi')
library('Homo.sapiens')
library('topGO')
library('RColorBrewer')
library('ggrepel')
library('readr')
library('tibble')
library('grid')
library('lattice')
library('gridExtra')
library('patchwork')
library('dplyr')
library('tximport')

# read in sample table that contains metadata (e.g. cell type, treatment, etc.)
dir <- 'F:/Connor/20220505_NextSeq_RNA-seq/'
csvfile <- file.path(dir, 'sample_table.csv')
sampleTable <- read.csv(csvfile,row.names=1)

# get file names of bam files - the roots should match what is in the sample table above
dir2 <- 'F:/Connor/20220505_NextSeq_RNA-seq/bam'
filenames1 <- file.path(dir2, paste0(sampleTable$Run, '.sorted.bam'))
file.exists(filenames1)

bamfiles <- BamFileList(filenames1)
seqinfo(bamfiles[1])

# get reference genome annotation file that contains gene names/positions for counting
dir3 = 'F:/Connor/transcriptomes/' 
gtffile <- file.path(dir3, 'GRCh38_20220506_HHVs.gtf')
txdb <- makeTxDbFromGFF(gtffile, format='gtf')

# https://sbc.shef.ac.uk/workshops/2019-01-14-rna-seq-r/rna-seq-preprocessing.nb.html#fixing_the_problem_with_transcript_names_-_the_hard_way
# k <- keys(txdb, keytype="TXNAME")
# tx_map <- AnnotationDbi::select(txdb, keys = k, columns="GENEID", keytype = "TXNAME")
# head(tx_map)
# tx2gene <- tx_map
# write.csv(tx2gene, file="F:/Connor/20220505_NextSeq_RNA-seq/tx2gene.csv", row.names=FALSE, quote=FALSE)
# txi <- tximport(quant_files, type="salmon",tx2gene = tx2gene)


ebg <- exonsBy(txdb, by='gene')
seqnames(ebg)

#use 8 cores - 16 crashes machine
register(SnowParam(6))

# create a SummarizedExperiment object with summarizeOverlaps
# note that we use Union mode, but other options include IntersectionStrict and IntersectionNotEmpty
se <- summarizeOverlaps(features=ebg,
                        reads=bamfiles,
                        mode='Union',
                        singleEnd=FALSE,
                        ignore.strand=TRUE,
                        fragments=FALSE)

# save counts matrix as csv as a record
colData(se) <- DataFrame(sampleTable)
write.table(assay(se), file='F:/Connor/20220505_NextSeq_RNA-seq/se_counts.csv', sep=',', quote=FALSE, row.names=TRUE, col.names=TRUE)

# read in counts as se so not to re-run the code (general purposes)
countsFile <- read.csv(file="F:/Connor/20220505_NextSeq_RNA-seq/se_counts.csv", stringsAsFactors=FALSE)
se <- SummarizedExperiment(data.matrix(countsFile))
colData(se) <- DataFrame(sampleTable)
se
assay(se)
# se2 <- se
# dds <- get_dds(1:24, 'plasmid', 'Rep')
# fpkm <- fpkm(dds, onlyZeroes=TRUE)
# fpkm
sampleTable


se2 <- se[rowSums(assay(se)) >=5, ]

# read in counts as se so not to re-run the code, and subset to remove transgenes so that these don't drive the variance during PCA/DESeq
countsFile <- read.csv(file="F:/Connor/20220505_NextSeq_RNA-seq/se_counts.csv", stringsAsFactors=FALSE)
countsFile <- countsFile[-c(1:4), ] # add this to remove EBNA counts
countsFile <- head(countsFile, -7)
se <- SummarizedExperiment(data.matrix(countsFile))
colData(se) <- DataFrame(sampleTable)
se
assay(se)
se2 <- se[rowSums(assay(se)) >=5, ]

# the following allows us to save the actual PCA values and percent variance for plotting in python
se1 <- se2[,15:22]
dds <- get_dds(15:22, 'plasmid', 'Rep')
vsd <- vst(dds)
p <- pca(assay(vsd), metadata = colData(se1))
write.csv(p$rotated, file=file.path(dir, 'EBNA2_PCA_table_vals-to-plot.csv'))
write.csv(p$variance, file=file.path(dir, 'EBNA2_PCA_table_perc-variance.csv'))



# general
make_PCA_plots(1:24, 'plasmid', 'Rep')
make_eigencorplot(1:24, 'plasmid', 'Rep')

# mCitrine +/-dox <-- this is a control to see how expressing a theoretically innocuous protein affects host gene expression
get_DESeq(1:4, 'dox', 'Rep', 'dox', 'nodox')
get_ordered_results(1:4, 'dox', 'Rep', 'dox', 'nodox', 'CL048_dox-nodox.csv')
make_volcano_plot(1:4, 'dox', 'Rep', 'dox', 'nodox')
get_GO_up(1:4, 'dox', 'Rep', 'dox', 'nodox', 'CL048_dox-nodox_GO_up.csv')
get_GO_down(1:4, 'dox', 'Rep', 'dox', 'nodox', 'CL048_dox-nodox_GO_down.csv')

# VIRF3 vs Citrine (dox)
get_DESeq(c(1:2, 5:6), 'plasmid', 'Rep', 'CL194', 'CL048')
get_ordered_results(c(1:2, 5:6), 'plasmid', 'Rep', 'CL194', 'CL048', 'CL194-CL048_dox.csv')
make_volcano_plot(c(1:2, 5:6), 'plasmid', 'Rep', 'CL194', 'CL048')
get_GO_up(c(1:2, 5:6), 'plasmid', 'Rep', 'CL194', 'CL048', 'CL194-CL048_dox_GO_up.csv')
get_GO_down(c(1:2, 5:6), 'plasmid', 'Rep', 'CL194', 'CL048', 'CL194-CL048_dox_GO_down.csv')

# VIRF4 vs Citrine (dox)
get_DESeq(c(1:2, 7:8), 'plasmid', 'Rep', 'CL195', 'CL048')
get_ordered_results(c(1:2, 7:8), 'plasmid', 'Rep', 'CL195', 'CL048', 'CL195-CL048_dox.csv')
make_volcano_plot(c(1:2, 7:8), 'plasmid', 'Rep', 'CL195', 'CL048')
get_GO_up(c(1:2, 7:8), 'plasmid', 'Rep', 'CL195', 'CL048', 'CL195-CL048_dox_GO_up.csv')
get_GO_down(c(1:2, 7:8), 'plasmid', 'Rep', 'CL195', 'CL048', 'CL195-CL048_dox_GO_down.csv')

# ORF10 vs Citrine (dox)
get_DESeq(c(1:2, 9:10), 'plasmid', 'Rep', 'CL246', 'CL048')
get_ordered_results(c(1:2, 9:10), 'plasmid', 'Rep', 'CL246', 'CL048', 'CL246-CL048_dox.csv')
make_volcano_plot(c(1:2, 9:10), 'plasmid', 'Rep', 'CL246', 'CL048')

# VIRF2 vs Citrine (dox)
get_DESeq(c(1:2, 11:12), 'plasmid', 'Rep', 'CL290', 'CL048')
get_ordered_results(c(1:2, 11:12), 'plasmid', 'Rep', 'CL290', 'CL048', 'CL290-CL048_dox.csv')
make_volcano_plot(c(1:2, 11:12), 'plasmid', 'Rep', 'CL290', 'CL048')
get_GO_up(c(1:2, 11:12), 'plasmid', 'Rep', 'CL290', 'CL048', 'CL290-CL048_dox_GO_up.csv')
get_GO_down(c(1:2, 11:12), 'plasmid', 'Rep', 'CL290', 'CL048', 'CL290-CL048_dox_GO_down.csv')

# IE1 vs Citrine (dox)
get_DESeq(c(1:2, 13:14), 'plasmid', 'Rep', 'CL402', 'CL048')
get_ordered_results(c(1:2, 13:14), 'plasmid', 'Rep', 'CL402', 'CL048', 'CL402-CL048_dox.csv')
make_volcano_plot(c(1:2, 13:14), 'plasmid', 'Rep', 'CL402', 'CL048')

# EBNA2_B-B-B vs Citrine (dox)
get_DESeq(c(1:2, 15:16), 'plasmid', 'Rep', 'CL403', 'CL048')
get_ordered_results(c(1:2, 15:16), 'plasmid', 'Rep', 'CL403', 'CL048', 'CL403-CL048_dox.csv')
make_volcano_plot(c(1:2, 15:16), 'plasmid', 'Rep', 'CL403', 'CL048')
get_GO_up(c(1:2, 15:16), 'plasmid', 'Rep', 'CL403', 'CL048', 'CL403-CL048_dox_GO_up.csv')
get_GO_down(c(1:2, 15:16), 'plasmid', 'Rep', 'CL403', 'CL048', 'CL403-CL048_dox_GO_down.csv')

# EBNA2_A-A-A vs Citrine (dox)
get_DESeq(c(1:2, 17:18), 'plasmid', 'Rep', 'CL404', 'CL048')
get_ordered_results(c(1:2, 17:18), 'plasmid', 'Rep', 'CL404', 'CL048', 'CL404-CL048_dox.csv')
make_volcano_plot(c(1:2, 17:18), 'plasmid', 'Rep', 'CL404', 'CL048')
get_GO_up(c(1:2, 17:18), 'plasmid', 'Rep', 'CL404', 'CL048', 'CL404-CL048_dox_GO_up.csv')
get_GO_down(c(1:2, 17:18), 'plasmid', 'Rep', 'CL404', 'CL048', 'CL404-CL048_dox_GO_down.csv')

# EBNA2_B-A-B vs Citrine (dox)
get_DESeq(c(1:2, 19:20), 'plasmid', 'Rep', 'CL405', 'CL048')
get_ordered_results(c(1:2, 19:20), 'plasmid', 'Rep', 'CL405', 'CL048', 'CL405-CL048_dox.csv')
make_volcano_plot(c(1:2, 19:20), 'plasmid', 'Rep', 'CL405', 'CL048')
get_GO_up(c(1:2, 19:20), 'plasmid', 'Rep', 'CL405', 'CL048', 'CL405-CL048_dox_GO_up.csv')
get_GO_down(c(1:2, 19:20), 'plasmid', 'Rep', 'CL405', 'CL048', 'CL405-CL048_dox_GO_down.csv')

# EBNA2_A-B-A vs Citrine (dox)
get_DESeq(c(1:2, 21:22), 'plasmid', 'Rep', 'CL406', 'CL048')
get_ordered_results(c(1:2, 21:22), 'plasmid', 'Rep', 'CL406', 'CL048', 'CL406-CL048_dox.csv')
make_volcano_plot(c(1:2, 21:22), 'plasmid', 'Rep', 'CL406', 'CL048')
get_GO_up(c(1:2, 21:22), 'plasmid', 'Rep', 'CL406', 'CL048', 'CL406-CL048_dox_GO_up.csv')
get_GO_down(c(1:2, 21:22), 'plasmid', 'Rep', 'CL406', 'CL048', 'CL406-CL048_dox_GO_down.csv')

# RL5A vs Citrine (dox)
get_DESeq(c(1:2, 23:24), 'plasmid', 'Rep', 'CL418', 'CL048')
get_ordered_results(c(1:2, 23:24), 'plasmid', 'Rep', 'CL418', 'CL048', 'CL418-CL048_dox.csv')
make_volcano_plot(c(1:2, 23:24), 'plasmid', 'Rep', 'CL418', 'CL048')

# EBNA2_B-B-B vs EBNA2_A-A-A (dox)
get_DESeq(15:18, 'plasmid', 'Rep', 'CL403', 'CL404')
get_ordered_results(15:18, 'plasmid', 'Rep', 'CL403', 'CL404', 'CL403-CL404_dox.csv')
make_volcano_plot(15:18, 'plasmid', 'Rep', 'CL403', 'CL404')
get_GO_up(15:18, 'plasmid', 'Rep', 'CL403', 'CL404', 'CL403-CL404_dox_GO_up.csv')
get_GO_down(15:18, 'plasmid', 'Rep', 'CL403', 'CL404', 'CL403-CL404_dox_GO_down.csv')

# EBNA2_B-B-B vs EBNA2_B-A-B (dox)
get_DESeq(c(15:16, 19:20), 'plasmid', 'Rep', 'CL403', 'CL405')
get_ordered_results(c(15:16, 19:20), 'plasmid', 'Rep', 'CL403', 'CL405', 'CL403-CL405_dox.csv')
make_volcano_plot(c(15:16, 19:20), 'plasmid', 'Rep', 'CL403', 'CL405')

# EBNA2_B-B-B vs EBNA2_A-B-A (dox)
get_DESeq(c(15:16, 21:22), 'plasmid', 'Rep', 'CL403', 'CL406')
get_ordered_results(c(15:16, 21:22), 'plasmid', 'Rep', 'CL403', 'CL406', 'CL403-CL406_dox.csv')
make_volcano_plot(c(15:16, 21:22), 'plasmid', 'Rep', 'CL403', 'CL406')

# EBNA2_A-A-A vs EBNA2_B-A-B (dox)
get_DESeq(17:20, 'plasmid', 'Rep', 'CL404', 'CL405')
get_ordered_results(17:20, 'plasmid', 'Rep', 'CL404', 'CL405', 'CL404-CL405_dox.csv')
make_volcano_plot(17:20, 'plasmid', 'Rep', 'CL404', 'CL405')

# EBNA2_A-A-A vs EBNA2_A-B-A (dox)
get_DESeq(c(17:18, 21:22), 'plasmid', 'Rep', 'CL404', 'CL406')
get_ordered_results(c(17:18, 21:22), 'plasmid', 'Rep', 'CL404', 'CL406', 'CL404-CL406_dox.csv')
make_volcano_plot(c(17:18, 21:22), 'plasmid', 'Rep', 'CL404', 'CL406')

# EBNA2_B-A-B vs EBNA2_A-B-A (dox)
get_DESeq(19:22, 'plasmid', 'Rep', 'CL405', 'CL406')
get_ordered_results(19:22, 'plasmid', 'Rep', 'CL405', 'CL406', 'CL405-CL406_dox.csv')
make_volcano_plot(19:22, 'plasmid', 'Rep', 'CL405', 'CL406')




