### Set working directory
setwd("/Users/claudiamagagnini/Desktop/Magagnini_1884218")

### Produce a count matrix with the gene lenght
library(GenomicFeatures)
library(tibble)


## Creation of the lenght table of all the genome genes

# import the GTF file
txdb <- makeTxDbFromGFF("/Users/claudiamagagnini/Desktop/Magagnini_1884218/gencode.vM36.basic.annotation.gtf", format="gtf")

# extract all the exonic ranges fo each gene
exons.list.per.gene <- exonsBy(txdb,by="gene")

# compute the total exonic lenght for each gene (x) of the exons.list.per.gene excluding overlapping regions. 
exonic.gene.sizes <- lapply(exons.list.per.gene,function(x){sum(width(reduce(x)))})

# convert the exonic.gene.sizes list into a dataframe
length.table <- as.data.frame(do.call(rbind, exonic.gene.sizes))

# moves the gene IDs into a proper column called "Gene ID"
length.table <- rownames_to_column(length.table, var = "Gene ID")

# rename the second column to "Length"
colnames(length.table)[2] <- "Length"



### Creation of the count table
library(rtracklayer)
library(GenomicRanges)
library(dplyr)

## Load the count matrix
raw.counts <-read.table("/Users/claudiamagagnini/Desktop/Magagnini_1884218/GSE255401_TroyKOvsHet_Raw_count_matrix.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
raw.counts <- raw.counts[, !colnames(raw.counts) %in% "X"]

# sum duplicated genes
raw.counts.sum <- raw.counts %>%
  group_by(Geneid, external_gene_name) %>%
  summarise(across(where(is.numeric), sum), .groups = "drop")

### Merge the length table and the count table

## remove suffixes from Ensembl gene IDs
length.table$Geneid <- sub("\\.\\d+$", "", length.table$`Gene ID`)

## merge the count table with the gene length table
counts.length <-merge(raw.counts.sum, length.table[, c("Geneid", "Length")], by = "Geneid")

## set the row names of the final matrix to the gene IDs
rownames(counts.length) <- counts.length$Geneid

## remove the Geneid excess column 
counts.length <- counts.length[, -which(names(counts.length) == "Geneid")]

### Load the sample annotation
sample_info <- read.table("//Users/claudiamagagnini/Desktop/Magagnini_1884218/GSE255401_TroyKOvsHet_sample_annotation.txt",
                          header = TRUE, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
## rename columns
colnames(sample_info) <- c("Sample", "Genotype", "Celltype", "Reporter", "Mouse_ID", "Sample_ID")

## replace "L001" with "L002" in sample names to have consistency with the matrix
sample_info$Sample <- gsub("L001", "L002", sample_info$Sample)

## match the sample order to the count matrix
sample_cols <- intersect(colnames(counts.length), sample_info$Sample)
counts <- counts.length[, sample_cols]

## prepare a gene lengths vector for FPKM calculation
gene.lengths <- counts.length$Length
names(gene.lengths) <- rownames(counts.length)


### DE analysis using edgeR
library(edgeR)
library(corrplot)
library(limma)
library(pheatmap)
library(biomaRt)

## separate raw counts and gene lengths from counts.length
sizes.all <- counts.length[, "Length", drop = FALSE]
counts.all <- counts.length %>%
  select(-external_gene_name, -Length)
  
## filter out low expressed genes
sample_info$group <- factor(paste(sample_info$Genotype, sample_info$Reporter, sep = "_"))
keep <- filterByExpr(counts.all, group = sample_info$group)
counts <- counts.all[keep, ]

## keep only the length of expressed genes
sizes <- sizes.all[keep, , drop=FALSE]$Length
names(sizes) <- rownames(sizes.all[keep, , drop=FALSE])

## create DGE object
dge <- DGEList(counts = counts)

## calculate normalization factors
dge <- calcNormFactors(dge)

## calculate fpkms and cpms
fpkms <- rpkm(dge, gene.length = sizes) 
cpms <- cpm(dge)

## write fpkms and cpms in a file
write.table(fpkms,"/Users/claudiamagagnini/Desktop/Magagnini_1884218/EdgeR/fpkms.txt")
write.table(cpms,"/Users/claudiamagagnini/Desktop/Magagnini_1884218/EdgeR/cpms.txt")

## calculate logCPM for subsequent clustering and heatmap
logcpm <- log2( cpms + 1)

## create a table of average fpkm values for each condition
samples_pos <- which(sample_info$Reporter == "pos")
samples_neg <- which(sample_info$Reporter == "neg")
fpkm.table <- data.frame(
  TroyPos = rowMeans(fpkms[, samples_pos]),
  TroyNeg = rowMeans(fpkms[, samples_neg])
)
rownames(fpkm.table) <- rownames(fpkms)

## create the design matrix
design <- model.matrix(~ 0 + group, data = sample_info)
colnames(design) <- levels(sample_info$group)

## estimate dispersion
dge <- estimateDisp(dge, design)

## fit a negative binomial generalized linear model
fit <- glmFit(dge, design)

### Open a pdf file for plotting
pdf("/Users/claudiamagagnini/Desktop/Magagnini_1884218/EdgeR/plots.pdf")

## create shorter labels for samples
short_labels <- gsub(".*S(\\d+)_.*", "S\\1", colnames(logcpm))
colnames(dge) <- short_labels
colnames(logcpm) <- short_labels
colnames(fpkms) <- short_labels
sample_info$Sample <- short_labels
rownames(sample_info) <- short_labels
colnames(cpms) <- short_labels


## MDS plot with distances computed based on logFC
plotMDS(dge)

## MDS plot with distances computed based on biological coefficient of variation
plotMDS(dge, method="bcv")

## correlation plot
M <- cor(logcpm)
corrplot(M, method = "number", number.cex = 0.7)

## hierarchical clustering of samples
distCor <- as.dist(1-M)
hc <- hclust(distCor)
plot(hc)

## create contrasts
my.contrasts <- makeContrasts(
  KOvsHet_Pos = Ko_pos - Het_pos,
  KOvsHet_Neg = Ko_neg - Het_neg,
  levels = design
)

## perform likelihood ratio test
lrt_KOvsHet_Pos <- glmLRT(fit, contrast = my.contrasts[,"KOvsHet_Pos"])
lrt_KOvsHet_Neg <- glmLRT(fit, contrast = my.contrasts[,"KOvsHet_Neg"])

## extract the DEGs tables for each contrast
tp_KOvsHet_Pos <- topTags(lrt_KOvsHet_Pos, n = Inf)$table
tp_KOvsHet_Neg <- topTags(lrt_KOvsHet_Neg, n = Inf)$table

## annotate DGE results

# create a tx2gene object: 
mart <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")
tx2gene <- getBM(attributes = c("ensembl_transcript_id", "ensembl_gene_id", "external_gene_name", "gene_biotype", "transcript_biotype"),
                 mart = mart)  
# rename columns
colnames(tx2gene) <- c("ensembl_transcript_id","ensembl_gene_id","gene_symbol","gene_type","transcript_type")

# remove duplicates
tx2gene.dedup <- tx2gene[!duplicated(tx2gene$ensembl_gene_id), ] 

# define function
annotateDE <- function(rTable, tx2gene.dedup) {
  rTable$gene_id <- rownames(rTable)
  
# merge DE results with the annotation table
  rTable <- merge(rTable, tx2gene.dedup, by.x = "gene_id", by.y = "ensembl_gene_id", all.x = TRUE)
  rownames(rTable) <- rTable$gene_id
  
# keep only relevant columns
  stat_cols <- c("gene_symbol", "ensembl_transcript_id", "transcript_type", "ensembl_gene_id", "gene_type", "logFC", "logCPM", "LR", "PValue", "FDR")
  stat_cols <- stat_cols[stat_cols %in% colnames(rTable)]
  rTable <- rTable[, stat_cols]
  
# sort by p-value
  rTable$PValue <- as.numeric(rTable$PValue)
  rTable <- rTable[order(rTable$PValue), ]
  
# filter out genes with NA FDR 
  rTable <- rTable[complete.cases(rTable$FDR), ]
  
  return(rTable)
}

# apply annotation
res_KOvsHet_Pos <- annotateDE(tp_KOvsHet_Pos, tx2gene.dedup)
res_KOvsHet_Neg <- annotateDE(tp_KOvsHet_Neg, tx2gene.dedup)

## save results on file
write.table(cbind(ensembl_gene_id = rownames(res_KOvsHet_Pos), res_KOvsHet_Pos), file = "/Users/claudiamagagnini/Desktop/Magagnini_1884218/EdgeR/results/KOvsHet_Pos.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(cbind(ensembl_gene_id = rownames(res_KOvsHet_Neg), res_KOvsHet_Neg), file = "/Users/claudiamagagnini/Desktop/Magagnini_1884218/EdgeR/results/KOvsHet_Neg.tsv", sep = "\t", quote = FALSE, row.names = FALSE)


### Plots

## Volcano Plots

# define function
plotVolcano<- function(vTable, comparison, lVal, pCut) {
with(vTable, plot(logFC, -log10(PValue), pch=20, xlim=c(-lVal, lVal),xlab="", ylab="")) 
with(subset(vTable, FDR< pCut), points(logFC, -log10(PValue), pch=20, col="red"))
title(comparison, xlab="logFC",ylab="-log10 p-value")
}

# apply function
plotVolcano(res_KOvsHet_Pos, "KO vs Het (Reporter Pos)", lVal = 5, pCut = 0.05)
plotVolcano(res_KOvsHet_Neg, "KO vs Het (Reporter Neg)", lVal = 5, pCut = 0.05)

## merge significantly expressed genes
DEGs_Pos <- rownames(subset(tp_KOvsHet_Pos, FDR < 0.05))
DEGs_Neg <- rownames(subset(tp_KOvsHet_Neg, FDR < 0.05))
 
all_DEGs <- union(DEGs_Pos, DEGs_Neg)

## Heatmap

#define function
expressionHeatmap <- function(DE1,DE2,exprMatrix,sample_info,pCut,...) {
  all_DEGs <- union( rownames(subset(DE1, FDR<pCut)), rownames(subset(DE2, FDR<pCut ))   )
  all_DEGs <- intersect(all_DEGs, rownames(exprMatrix))
  expr_subset <- exprMatrix[all_DEGs, , drop=FALSE]
  annotation_col <- sample_info[, c("Genotype", "Reporter")]
  rownames(annotation_col) <- sample_info$Sample
  pHM <- pheatmap(expr_subset, annotation_col= annotation_col,
                  show_rownames = FALSE,
                  scale = "row",
                  ...)
  print(pHM)
}

#apply function
expressionHeatmap(res_KOvsHet_Pos, res_KOvsHet_Neg,
                        exprMatrix = logcpm,
                        sample_info = sample_info,
                        pCut = 0.05)

### close pdf file
dev.off()



### Functional Enrichment Analysis


## ORA
library("WebGestaltR")

## KOvsHet_Pos

projectName_KOvsHet_Pos <- "KOvsHet_Pos"

# create a geneList with up-regulated and down-regulated genes
up_regulated_KOvsHet_Pos <- rownames(subset(res_KOvsHet_Pos, FDR < 0.05 & logFC >= 0))
down_regulated_KOvsHet_Pos <-rownames(subset(res_KOvsHet_Pos, FDR < 0.05 & logFC <= 0))

geneList_KOvsHet_Pos <- list("up_regulated" = up_regulated_KOvsHet_Pos,
                        "down_regulated" = down_regulated_KOvsHet_Pos)
geneListNames_KOvsHet_Pos <- names(geneList_KOvsHet_Pos)

# define the refList
refList <- rownames(counts.length)

# print the list of organisms
listOrganism(hostName = "http://www.webgestalt.org/")

# print the list of supported id types
listIdType(organism = "mmusculus",hostName = "http://www.webgestalt.org/")

# print the list of gene sets
genesets <- listGeneSet(organism = "mmusculus",
                        hostName = "http://www.webgestalt.org/")


# define categories of interest
categories <- c("geneontology_Biological_Process_noRedundant", "geneontology_Cellular_Component_noRedundant", 
                "geneontology_Molecular_Function_noRedundant", "pathway_Reactome")


# create output directory for each category
outputDirs <- file.path("/Users/claudiamagagnini/Desktop/Magagnini_1884218/webgestalt/ORA",categories)

# create runORA function
runORA_category_Pos <- function(category, geneSetName, geneSet) {
  outputDir <- file.path("/Users/claudiamagagnini/Desktop/Magagnini_1884218/webgestalt/ORA",
    category,projectName_KOvsHet_Pos,geneSetName)
  dir.create(outputDir, recursive = TRUE, showWarnings = FALSE)
  
WebGestaltR(
    enrichMethod = "ORA",
    organism = "mmusculus",
    enrichDatabase = category,
    interestGene = geneSet,
    referenceGene = refList,
    interestGeneType = "ensembl_gene_id",
    referenceGeneType = "ensembl_gene_id",
    outputDirectory = outputDir,
    isOutput = TRUE,
    projectName = projectName_KOvsHet_Pos, 
    nThreads = 2,
    hostName = "https://www.webgestalt.org/"
  )
}

# run webgestalt ORA
ORA_results_Pos <- lapply(geneListNames_KOvsHet_Pos, function(geneSetName) {
  geneSet <- geneList_KOvsHet_Pos[[geneSetName]]
  mapply(function(category, outputDir) {
    runORA_category_Pos(category, geneSetName, geneSet)
  }, categories, outputDirs, SIMPLIFY = FALSE)
})


###GSEA 

## KOvsHet_Pos

# remove rows with NA values
res_KOvsHet_Pos.GSEA <- res_KOvsHet_Pos[complete.cases(res_KOvsHet_Pos),]

# extract logFC sign
res_KOvsHet_Pos.GSEA$FCsign <- sign(res_KOvsHet_Pos.GSEA$logFC)

# extract PValue and calculate -log10(Pvalue),
res_KOvsHet_Pos.GSEA$logP <- -log10(res_KOvsHet_Pos.GSEA$PValue)

# apply logFC sign to transformed PValue and add it as a metric column to the table
res_KOvsHet_Pos.GSEA$metric <- res_KOvsHet_Pos.GSEA$logP/res_KOvsHet_Pos.GSEA$FCsign

# create the ranked list for GSEA by extracting gene ID and metric
res_KOvsHet_Pos.GSEA$gene_id <- rownames(res_KOvsHet_Pos.GSEA)
ranklist_KOvsHet_Pos <- res_KOvsHet_Pos.GSEA[, c("gene_id", "metric")]

# create a GSEA function
runGSEA_KOvsHet_Pos <- function(category) {
  outputDirectory <- file.path("/Users/claudiamagagnini/Desktop/Magagnini_1884218/webgestalt/GSEA", category)
  dir.create(outputDirectory, showWarnings = FALSE)

  WebGestaltR(enrichMethod="GSEA", 
              organism = "mmusculus", 
              enrichDatabase = category, 
              interestGene = ranklist_KOvsHet_Pos,
              interestGeneType = "ensembl_gene_id",
              outputDirectory = outputDirectory, 
              isOutput = TRUE,
              minNum = 10, maxNum = 1000, 
              nThreads = 2, 
              sigMethod = "top",
              topThr = 20,
              projectName = projectName_KOvsHet_Pos,
              hostName = "https://www.webgestalt.org/")}

# run webgestalt GSEA 
GSEAResults_KOvsHet_Pos <- lapply(categories, runGSEA_KOvsHet_Pos)
