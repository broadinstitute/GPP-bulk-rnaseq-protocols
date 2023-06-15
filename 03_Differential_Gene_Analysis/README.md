Differential Gene Expression Analysis
================
Dany Gould
2023-06-14

There are numerous tools to test for differentially expressed genes
(DEGs), all different in their approach and results. Even within the
same library there are multiple ways to carry out the analysis. This
folder contains example code for two of the most popular libraries
currently -
[Limma](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4402510/) and
[DESeq2](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-014-0550-8).
We will default to using DESeq2 but it may be necessary to adopt Limma
(or other alternatives) depending on experimental specifics.

### DESeq2

------------------------------------------------------------------------

#### 1. Import packages

``` r
library(DESeq2)
library(edgeR)
library(ggpubr)
library(ggplot2)
library(data.table)
library(biomaRt) # for importing ensembl data
```

#### 2. Import Read files

Read files are the `.rsem.genes.results` or `.rsem.genes.results.gz`
files downloaded from Terra

``` r
reads_file_dir <- "/Users/fu/Library/CloudStorage/GoogleDrive-fu@broadinstitute.org/Shared drives/GPP RNA Seq/2023/TP53/Terra/RSEM_Results/All-reads/"
read_files <- list.files(path=reads_file_dir, full.names=TRUE)
# Create a dataframe of reads from each file
reads_list <- lapply(read_files, function(x) {
  r <- fread(x, select = c("gene_id", "expected_count"))
  file_name <- unlist(strsplit(x, "\\."))[2]
  # Change the expected_count column to the sample name
  setnames(r, c("gene_id", sub(".*TP53-A549-", "", file_name)))
})

# Merge all the reads into one dataframe
reads_df <- Reduce(function(x, y) merge(x, y, by = "gene_id", all=TRUE), reads_list)
reads_df<- as.data.frame(reads_df)
# Gene IDs are in the first column, make them rownames and then drop the column
rownames(reads_df) <- reads_df[,1]
reads_df[,1] <- NULL
```

#### 3. Pre-process the reads

``` r
# Make sure the samples in the annotation file are in the same order as the columns names of the dataframe!
annot <- read.csv("/Users/fu/Library/CloudStorage/GoogleDrive-fu@broadinstitute.org/Shared drives/GPP Cloud /R&D/People/Dany/RNAseq Analysis/TP53 Base Editing/TP53-BE-tiling-annotation.csv")

# Remove 569-Nutlin-RepA from analysis (see 02_Quality_Control for reasoning behind this)
annot.minus.outlier <- annot[annot$Outlier=="FALSE",]
reads_df.minus.outlier <- subset(reads_df, select=-c(`RDA569-Nutlin-RepA`))

# Make DGE object
reads.dge.minus.outlier <- edgeR::DGEList(counts=reads_df.minus.outlier, samples = annot.minus.outlier)
reads.dge.minus.outlier$samples$group <- as.factor(c(annot.minus.outlier$Sample))

# Remove lowly expressed genes
# There are many valid ways to do so - for example, remove any row where the sum of the row < 10
# Here we use the filterByExpr from the edgeR package described by Chen et al (2016)
keep.exprs.reads <- edgeR::filterByExpr(reads.dge.minus.outlier, group = reads.dge.minus.outlier$samples$group)
reads.dge.minus.outlier <- reads.dge.minus.outlier[keep.exprs.reads,, keep.lib.sizes=FALSE]

# DESeq2 requires all data to be non-negative integers
reads.matrix.filtered <- round(as.data.frame(reads.dge.minus.outlier$counts)) 
```

Note: we do not need to normalize the reads prior to using the DESeq2
library.

#### 4. Differentially Expressed Genes (DEG) Analysis

##### 4.1 Create Design Matrix

Refer to this
[guide](https://bioconductor.org/packages/release/workflows/vignettes/RNAseq123/inst/doc/designmatrices.html)
to understand design matrix terminology and how to create a suitable one
for your experiment.

Note: In an orthogonal data set such as this, the Means model and
Means-Reference model will give the same results (see Limma code
example). Personally, I find the Means model (and the use of contrasts)
to be more intuitive and easier to code. However, the following is an
example of the Means-Reference model because the default shrinkage
algorithm (`apeglm`) used in a later step is not compatible with
`contrasts` (see step 4.5).

``` r
# In this experiment we combine guide and treatment into one factor so that the data set is orthogonal 
design <- as.factor(reads.dge.minus.outlier$samples$group)
design
```

    ##  [1] nutlin_ctrl    nutlin_ctrl    nutlin_ctrl    untreated_ctrl untreated_ctrl
    ##  [6] untreated_ctrl nutlin_sg1     nutlin_sg1     untreated_sg1  untreated_sg1 
    ## [11] untreated_sg1  nutlin_sg9     nutlin_sg9     nutlin_sg9     untreated_sg9 
    ## [16] untreated_sg9  untreated_sg9 
    ## 6 Levels: nutlin_ctrl nutlin_sg1 nutlin_sg9 untreated_ctrl ... untreated_sg9

##### 4.2 Set Reference

By default, R will choose the first element in the list (by alphabetical
order) as the reference. We can reset the reference with the `relevel`
function.

``` r
design <- relevel(design, "untreated_ctrl")
# The first "level" is the new reference but the order of the list will not change.
design
```

    ##  [1] nutlin_ctrl    nutlin_ctrl    nutlin_ctrl    untreated_ctrl untreated_ctrl
    ##  [6] untreated_ctrl nutlin_sg1     nutlin_sg1     untreated_sg1  untreated_sg1 
    ## [11] untreated_sg1  nutlin_sg9     nutlin_sg9     nutlin_sg9     untreated_sg9 
    ## [16] untreated_sg9  untreated_sg9 
    ## 6 Levels: untreated_ctrl nutlin_ctrl nutlin_sg1 nutlin_sg9 ... untreated_sg9

##### 4.3 Create table of sample information

The first column of the dataframe contains the name of samples *which
must match* the order of samples in the matrix that contains the reads.
The number of subsequent columns must match the number of factors in
your design.

``` r
# In this experiment we only have one factor
info <- data.frame(row.names = colnames(reads.matrix.filtered), design)
info
```

    ##                               design
    ## RDA429-Nutlin-RepA       nutlin_ctrl
    ## RDA429-Nutlin-RepB       nutlin_ctrl
    ## RDA429-Nutlin-RepC       nutlin_ctrl
    ## RDA429-untreated-RepA untreated_ctrl
    ## RDA429-untreated-RepB untreated_ctrl
    ## RDA429-untreated-RepC untreated_ctrl
    ## RDA569-Nutlin-RepB        nutlin_sg1
    ## RDA569-Nutlin-RepC        nutlin_sg1
    ## RDA569-untreated-RepA  untreated_sg1
    ## RDA569-untreated-RepB  untreated_sg1
    ## RDA569-untreated-RepC  untreated_sg1
    ## RDA692-Nutlin-RepA        nutlin_sg9
    ## RDA692-Nutlin-RepB        nutlin_sg9
    ## RDA692-Nutlin-RepC        nutlin_sg9
    ## RDA692-untreated-RepA  untreated_sg9
    ## RDA692-untreated-RepB  untreated_sg9
    ## RDA692-untreated-RepC  untreated_sg9

##### 4.3 Create DEseq2 Data Object

``` r
# The ~ in front of design indicates that this is a Mean-Reference model
data <- DESeqDataSetFromMatrix(countData = reads.matrix.filtered, 
                               colData = info, 
                               design= ~design)
```

##### 4.4 Run DESeq2

``` r
dds <- DESeq(data)
# Default estimates for DESeq2
res_dds <- results(dds)
head(res_dds)
```

    ## log2 fold change (MLE): design untreated sg9 vs untreated ctrl 
    ## Wald test p-value: design untreated sg9 vs untreated ctrl 
    ## DataFrame with 6 rows and 6 columns
    ##                     baseMean log2FoldChange     lfcSE      stat      pvalue
    ##                    <numeric>      <numeric> <numeric> <numeric>   <numeric>
    ## ENSG00000000003.14 1592.3889     -0.1453116 0.0844991 -1.719682 0.085490172
    ## ENSG00000000419.12  759.0205     -0.0932702 0.0757468 -1.231341 0.218195507
    ## ENSG00000000457.13  128.0968     -0.0551149 0.1743322 -0.316149 0.751889480
    ## ENSG00000000460.16  378.3852      0.0620898 0.0967179  0.641968 0.520894118
    ## ENSG00000000938.12   11.9082      2.3057122 0.6305278  3.656797 0.000255386
    ## ENSG00000000971.15 2998.2009     -0.0928059 0.1433031 -0.647620 0.517230987
    ##                          padj
    ##                     <numeric>
    ## ENSG00000000003.14 0.31436393
    ## ENSG00000000419.12 0.52597333
    ## ENSG00000000457.13 0.90801719
    ## ENSG00000000460.16 0.79112038
    ## ENSG00000000938.12 0.00448122
    ## ENSG00000000971.15 0.78842176

##### 4.5 Shrink Log Fold Change (LFC)

``` r
# Use the resultNames function to get all the coefficients with the reference that we set in step 4.2
resultsNames(dds)
```

    ## [1] "Intercept"                             
    ## [2] "design_nutlin_ctrl_vs_untreated_ctrl"  
    ## [3] "design_nutlin_sg1_vs_untreated_ctrl"   
    ## [4] "design_nutlin_sg9_vs_untreated_ctrl"   
    ## [5] "design_untreated_sg1_vs_untreated_ctrl"
    ## [6] "design_untreated_sg9_vs_untreated_ctrl"

As mentioned previously, the `apeglm` algorithm is not compatible with
the use of `contrasts`, but `ashr` is. See [extended section on
shrinkage
estimators](http://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#extended-section-on-shrinkage-estimators).
We did not explore whether the result is the same between `apeglm` vs
`ashr`.

``` r
resLFC <- lfcShrink(dds, coef="design_untreated_sg1_vs_untreated_ctrl", type="apeglm")
resLFC
```

    ## log2 fold change (MAP): design untreated sg1 vs untreated ctrl 
    ## Wald test p-value: design untreated sg1 vs untreated ctrl 
    ## DataFrame with 16301 rows and 5 columns
    ##                     baseMean log2FoldChange     lfcSE      pvalue        padj
    ##                    <numeric>      <numeric> <numeric>   <numeric>   <numeric>
    ## ENSG00000000003.14 1592.3889     -0.4355068 0.0861720 9.39897e-08 1.16550e-06
    ## ENSG00000000419.12  759.0205     -0.2195115 0.0762236 2.21082e-03 1.29680e-02
    ## ENSG00000000457.13  128.0968     -0.0325750 0.1423997 7.74533e-01 9.01925e-01
    ## ENSG00000000460.16  378.3852     -0.3318970 0.1004589 2.99287e-04 2.20961e-03
    ## ENSG00000000938.12   11.9082      0.0920927 0.2465692 2.45210e-01 4.90727e-01
    ## ...                      ...            ...       ...         ...         ...
    ## ERCC-00096           95.7556     0.16733766  0.194161   0.2034789    0.438682
    ## ERCC-00113           20.4528     0.14314654  0.273172   0.1360344    0.338108
    ## ERCC-00130          214.4978     0.24107872  0.148350   0.0462672    0.158133
    ## ERCC-00136           12.4840     0.00101547  0.221664   0.9916441    0.997197
    ## ERCC-00171           14.0237     0.04334880  0.218856   0.6408350    0.825552

### Visualization

------------------------------------------------------------------------

### Notes

------------------------------------------------------------------------

It is difficult to put into a single guide how RNASeq analysis can be
accomplished, as there are many decisions to be made throughout the
process. For example, we analyzed the data set as a whole in this
experiment, allowing the linear model to “borrow” information across all
the samples. Subsetting the data *will* change the results, and while
the creators of Limma recommend analyzing the whole experiment together,
this is not a hard and fast rule. See the [discussion on
Bioconductor](https://support.bioconductor.org/p/p132527/) and the
answer from the creator of DESeq2 on this topic. We have also included a
[code sample](DESeq2/DESeq2-data-subset.md) using DESeq2 where the data
was subsetted.
