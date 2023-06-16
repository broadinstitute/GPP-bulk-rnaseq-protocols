library(DESeq2)
library(edgeR)
library(ggpubr)
library(ggplot2)
library(data.table)
library(biomaRt) # for importing ensembl data

mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
annot <- read.csv(paste(analysis_dir, "TP53-BE-tiling-annotation.csv", sep = ""))

####################
# HELPER FUNCTIONS #
####################
DESEQ.makeResults = function(reads, condition, reference) {
  design <- relevel(condition, reference)
  # Rows of colData correspond to columns of countData
  info <- data.frame(row.names = colnames(reads), design)
  data <- DESeqDataSetFromMatrix(countData = reads, 
                                 colData = info, 
                                 design= ~design)
  return(DESeq(data))
}

DESEQ.resultWithGeneNames = function(res) {
  res <- as.data.frame(res)
  res$geneID <- rownames(res)
  # Remove the version number
  res$geneID <- lapply(res$geneID, function(x) unlist(strsplit(x, "\\."))[1])
  res$geneID<-as.character(res$geneID)
  genelist <- getBM(filters=c("ensembl_gene_id"), attributes= c("ensembl_gene_id","hgnc_symbol"), 
                    values=res$geneID, mart=mart)
  return(merge(res, genelist, by.x="geneID", by.y="ensembl_gene_id"))
}

DESEQ.meanDifferencePlot = function(deseq.result, title) {
  return(
    ggplot(deseq.result, aes(x=refAvg, y=log2FoldChange, color=log2FoldChange < 0 )) + 
      geom_point() + 
      theme_classic() + 
      scale_y_continuous(limits = c(-8, 8)) + 
      labs(title = title, x = "Average Expression", color = "Expression Change") + 
      scale_color_manual(labels = c("Upregulated", "Downregulated"), values = c("pink", "lightblue")) + 
      geom_text(aes(label=ifelse(hgnc_symbol=="TP53", hgnc_symbol,'')), color="black") + 
      geom_text(aes(label=ifelse(hgnc_symbol=="MDM2", hgnc_symbol,'')), color="black") 
  )
}



# Remove 569-Nutlin-RepA from analysis 
annot.minus.outlier <- annot[annot$Outlier=="FALSE",]
reads_df.minus.outlier <- subset(reads_df, select=-c(`RDA569-Nutlin-RepA`))

# Make DGE object
reads.dge.minus.outlier <- DGEList(counts=reads_df.minus.outlier, samples = annot.minus.outlier)
reads.dge.minus.outlier$samples$group <- as.factor(c(annot.minus.outlier$Sample))

# Remove lowly expressed genes
keep.exprs.reads <- filterByExpr(reads.dge.minus.outlier, group = reads.dge.minus.outlier$samples$group)
reads.dge.minus.outlier <- reads.dge.minus.outlier[keep.exprs.reads,, keep.lib.sizes=FALSE]
reads.matrix.filtered <- round(as.data.frame(reads.dge.minus.outlier$counts)) # must be non-negative integers


# Cannot use contrast because it is not compatible with apeglm, 
# we must do this the Means Reference way and re-level as we go 
condition.nutlin <- as.factor(reads.dge.minus.outlier$samples$group)
deseq.refList <- c("nutlin_ctrl", "untreated_ctrl", "untreated_sg1", "untreated_sg9", "nutlin_sg1")
deseq.res.list <- list()
for(ref in deseq.refList){
  print(ref)
  deseq.res.list[[ref]] <- DESEQ.makeResults(reads.matrix.filtered, condition.nutlin, ref)
}

resultsNames(deseq.res.list[["untreated_sg9"]])

# Get all pairwise comparisons, 6 choose 2 is 15 
# Do shrink LFC using the default apeglm algorithm
deseq.res.list <- list("Untreated sg1 vs Untreated Control" = 
                       lfcShrink(deseq.res.list[["untreated_ctrl"]], coef="condition_untreated_sg1_vs_untreated_ctrl"), 
                     "Untreated sg9 vs Untreated Control" = 
                       lfcShrink(deseq.res.list[["untreated_ctrl"]], coef="condition_untreated_sg9_vs_untreated_ctrl"), 
                     "Nutlin Control vs Untreated Control" = 
                       lfcShrink(deseq.res.list[["untreated_ctrl"]], coef="condition_nutlin_ctrl_vs_untreated_ctrl"), 
                     "Nutlin sg1 vs Untreated Control" = 
                       lfcShrink(deseq.res.list[["untreated_ctrl"]], coef="condition_nutlin_sg1_vs_untreated_ctrl"), 
                     "Nutlin sg9 vs Untreated Control" = 
                       lfcShrink(deseq.res.list[["untreated_ctrl"]], coef="condition_nutlin_sg9_vs_untreated_ctrl"), 
                     
                     "Untreated sg9 vs Untreated sg1" = 
                       lfcShrink(deseq.res.list[["untreated_sg1"]], coef="condition_untreated_sg9_vs_untreated_sg1"), 
                     "Nutlin Control vs Untreated sg1" = 
                       lfcShrink(deseq.res.list[["untreated_sg1"]], coef="condition_nutlin_ctrl_vs_untreated_sg1"), 
                     "Nutlin sg1 vs Untreated sg1" = 
                       lfcShrink(deseq.res.list[["untreated_sg1"]], coef="condition_nutlin_sg1_vs_untreated_sg1"), 
                     "Nutlin sg9 vs Untreated sg1" = 
                       lfcShrink(deseq.res.list[["untreated_sg1"]], coef="condition_nutlin_sg9_vs_untreated_sg1"), 
                     
                     "Nutlin Control vs Untreated sg9" = 
                       lfcShrink(deseq.res.list[["untreated_sg9"]], coef="condition_nutlin_ctrl_vs_untreated_sg9"), 
                     "Nutlin sg1 vs Untreated sg9" = 
                       lfcShrink(deseq.res.list[["untreated_sg9"]], coef="condition_nutlin_sg1_vs_untreated_sg9"), 
                     "Nutlin sg9 vs Untreated sg9" = 
                       lfcShrink(deseq.res.list[["untreated_sg9"]], coef="condition_nutlin_sg9_vs_untreated_sg9"), 
                     
                     "Nutlin sg1 vs Nutlin Control" = 
                       lfcShrink(deseq.res.list[["nutlin_ctrl"]], coef="condition_nutlin_sg1_vs_nutlin_ctrl"), 
                     "Nutlin sg9 vs Nutlin Control" = 
                       lfcShrink(deseq.res.list[["nutlin_ctrl"]], coef="condition_nutlin_sg9_vs_nutlin_ctrl"), 
                     
                     "Nutlin sg9 vs Nutlin sg1" = 
                       lfcShrink(deseq.res.list[["nutlin_sg1"]], coef="condition_nutlin_sg9_vs_nutlin_sg1")
                     )


##############################
# LFC vs Reference Mean Plot #
##############################

# the baseMean column of the result dataframe is the average expression of ALL samples
# We need to get the mean expression of reference-only
refAvg.df <- data.frame(matrix(ncol=4, nrow=nrow(reads.dge.lcpm), 
                               dimnames=list(rownames(reads.dge.lcpm), 
                                             c("untreated_ctrl_mean", "untreated_sg1_mean", "untreated_sg9_mean",
                                               "nutlin_ctrl_mean", "nutlin_sg1_mean"))))
rownames(refAvg.df) <- lapply(rownames(refAvg.df), function(x) unlist(strsplit(x, "\\."))[1])
refAvg.df$untreated_ctrl_mean <- rowMeans(reads.dge.lcpm[, annot$Sample == "untreated_ctrl"])
refAvg.df$nutlin_ctrl_mean <- rowMeans(reads.dge.lcpm[, annot$Sample == "nutlin_ctrl"])
refAvg.df$untreated_sg1_mean <- rowMeans(reads.dge.lcpm[, annot$Sample == "untreated_sg1"])
refAvg.df$untreated_sg9_mean <- rowMeans(reads.dge.lcpm[, annot$Sample == "untreated_sg9"])
refAvg.df$nutlin_sg1_mean <- rowMeans(reads.dge.lcpm[, annot$Sample == "nutlin_sg1"])

refMap <- c(
  "Untreated sg1 vs Untreated Control" = "untreated_ctrl_mean", 
  "Untreated sg9 vs Untreated Control" = "untreated_ctrl_mean", 
  "Nutlin Control vs Untreated Control" = "untreated_ctrl_mean",
  "Nutlin sg1 vs Untreated Control" = "untreated_ctrl_mean",
  "Nutlin sg9 vs Untreated Control" = "untreated_ctrl_mean",
  
  "Untreated sg9 vs Untreated sg1" = "untreated_sg1_mean", 
  "Nutlin Control vs Untreated sg1" = "untreated_sg1_mean",
  "Nutlin sg1 vs Untreated sg1" = "untreated_sg1_mean",
  "Nutlin sg9 vs Untreated sg1" = "untreated_sg1_mean",
  
  "Nutlin Control vs Untreated sg9" = "untreated_sg9_mean",
  "Nutlin sg1 vs Untreated sg9" = "untreated_sg9_mean",
  "Nutlin sg9 vs Untreated sg9" = "untreated_sg9_mean",
  
  "Nutlin sg1 vs Nutlin Control" = "nutlin_ctrl_mean", 
  "Nutlin sg9 vs Nutlin Control" = "nutlin_ctrl_mean", 
  
  "Nutlin sg9 vs Nutlin sg1" = "nutlin_sg1_mean", 
)

deseq.plotList <- list()
for (name in names(deseq.res.list)) {
  print(name)
  #.Machine$double.xmin gives the value of the smallest positive number whose 
  # representation meets the requirements of IEEE 754 technical standard for floating point computation.
  deseq.res.list[[name]]$padj[deseq.res.list[[name]]$padj == 0] <- .Machine$double.xmin
  
  # Add the reference mean as a new column
  deseq.res.list[[name]]$refAvg <- refAvg.df[refMap[name]][deseq.res.list[[name]]$geneID,]
  
  # Add gene names
  deseq.res.list[[name]] <- DESEQ.resultWithGeneNames(deseq.res.list[[name]])
  
  # make plot
  deseq.plotList[[name]] <- DESEQ.meanDifferencePlot(deseq.res.list[[name]], name)
}

pdf(paste(analysis_dir, "Figures/DESeq/DESeq-MeanModel-mean-reference.pdf", sep=""))
ggarrange(plotlist=deseq.plotList, ncol=2, nrow=2, common.legend = TRUE, legend="bottom")
dev.off()
