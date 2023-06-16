library(limma)
library(edgeR)
library(ggplot2)
library(ggpmisc)
library(ggpubr) # ggarrange
library(biomaRt)

mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
annot <- read.csv(paste(analysis_dir, "TP53-BE-tiling-annotation.csv", sep = ""))

####################
# HELPER FUNCTIONS #
####################
makeTopTableWithGeneNames = function(fit, colName){
  toptable <- topTable(fit, coef=colName, number=Inf, adjust.method="BH")
  toptable$geneID <- rownames(toptable)
  toptable$geneID <- lapply(toptable$geneID, function(x) unlist(strsplit(x, "\\."))[1])
  toptable$geneID<-as.character(toptable$geneID)
  genelist <- getBM(filters=c("ensembl_gene_id"), attributes= c("ensembl_gene_id","hgnc_symbol"), values=toptable$geneID, mart=mart)
  return(merge(toptable, genelist, by.x="geneID", by.y="ensembl_gene_id"))
}

meanDifferencePlot = function(topTable, title) {
  return(
    ggplot(topTable, aes(x=RefAvg, y=logFC, color=logFC < 0 )) + 
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
reads.dge.minus.outlier <- calcNormFactors(reads.dge.minus.outlier, method = "TMM")


# No intercept, Means Model
# https://bioconductor.org/packages/release/workflows/vignettes/RNAseq123/inst/doc/designmatrices.html#means-model-for-factors
group <- reads.dge.minus.outlier$samples$group
design.no_intercept <- model.matrix(~0+group)
colnames(design.no_intercept) <- gsub("group", "", colnames(design.no_intercept))

# Make contrasts for the comparisons we care about
contrast.no_intercept <- makeContrasts(
  "Untreated sg1 vs Untreated Control" = untreated_sg1 - untreated_ctrl, 
  "Untreated sg9 vs Untreated Control" = untreated_sg9 - untreated_ctrl, 
  "Nutlin Control vs Untreated Control" = nutlin_ctrl - untreated_ctrl,
  "Nutlin sg1 vs Untreated Control" = nutlin_sg1 - untreated_ctrl,
  "Nutlin sg9 vs Untreated Control" = nutlin_sg9 - untreated_ctrl,
  
  "Untreated sg9 vs Untreated sg1" = untreated_sg9 - untreated_sg1, 
  "Nutlin Control vs Untreated sg1" = nutlin_ctrl - untreated_sg1,
  "Nutlin sg1 vs Untreated sg1" = nutlin_sg1 - untreated_sg1,
  "Nutlin sg9 vs Untreated sg1" = nutlin_sg9 - untreated_sg1,
  
  "Nutlin Control vs Untreated sg9" = nutlin_ctrl - untreated_sg9,
  "Nutlin sg1 vs Untreated sg9" = nutlin_sg1 - untreated_sg9,
  "Nutlin sg9 vs Untreated sg9" = nutlin_sg9 - untreated_sg9,
  
  "Nutlin sg1 vs Nutlin Control" = nutlin_sg1 - nutlin_ctrl, 
  "Nutlin sg9 vs Nutlin Control" = nutlin_sg9 - nutlin_ctrl, 
  
  "Nutlin sg9 vs Nutlin sg1" = nutlin_sg9 - nutlin_sg1, 
  levels = colnames(design.no_intercept))
contrast.no_intercept

################
# REGULAR VOOM #
################
v.no_intercept <- voom(reads.dge.minus.outlier, design.no_intercept, plot=TRUE)
fit.no_intercept <- lmFit(v.no_intercept, design.no_intercept)
fit.no_intercept <- contrasts.fit(fit.no_intercept, contrasts=contrast.no_intercept)
fit.no_intercept <- eBayes(fit.no_intercept)
summary(decideTests(fit.no_intercept))

toptables <- list()
for(name in colnames(fit.no_intercept)){
  toptables[[name]] <- makeTopTableWithGeneNames(fit.no_intercept, colName=name)
}

##############################
# LFC vs Reference Mean Plot #
##############################
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

lfc.plotList <- list()
for(name in names(toptables.qw)){
  toptables.qw[[name]]$RefAvg <- refAvg.df[refMap[name]][toptables.qw[[name]]$geneID,]
  lfc.plotList[[name]] <- meanDifferencePlot(toptables.qw[[name]], name) 
}

# Condense the legend into one at the bottom of the page instead of one for every graph
pdf(paste(analysis_dir, "Figures/Limma-MeanModel-mean-reference.pdf", sep=""))
ggpubr::ggarrange(plotlist=lfc.plotList, ncol=2, nrow=2, common.legend = TRUE, legend="bottom")
dev.off()