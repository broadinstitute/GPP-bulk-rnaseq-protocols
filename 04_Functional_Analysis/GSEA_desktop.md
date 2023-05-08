## GSEA - Desktop
Unlike the `fgsea` package, using the Broad desktop version does not require that you run differential gene analysis prior. Instead, you can upload a normalized read counts file and use their supplied ranking metrics. The option to run GSEAPreranked in the software is also available. 

1. Install the [application](https://www.gsea-msigdb.org/gsea/downloads.jsp)

2. Create a normalized read counts file   
The `reads.dge.lcpm` object referred to here is the log2 normalized reads. Steps to creating this file can be found in Step #3 of the Read Count Metrics section in the Quality Control folder. 
```
lcpm.gsea <- as.data.frame(reads.dge.lcpm)
lcpm.gsea$geneID <- rownames(lcpm.gsea)
lcpm.gsea$geneID <- lapply(lcpm.gsea$geneID, function(x) unlist(strsplit(x, "\\."))[1])
lcpm.gsea$geneID<-as.character(lcpm.gsea$geneID)

# Check for duplicates and decide what to do with them before GSEA 
lcpm.gsea[duplicated(lcpm.gsea$geneID),] 

# Properly format for GSEA
lcpm.gsea$Description <- ""
lcpm.gsea <- lcpm.gsea[,c(19, 20, 1:18)] # Rearrange the order
colnames(lcpm.gsea)[1] <- "Name"

write.table(lcpm.gsea, file=paste(analysis_dir, "lcpm.txt"), quote=FALSE, sep='\t', row.names=FALSE)
```

3. Create a Phenotype Label file     
You can manually create this file based on your experiment. Details of how to do so is outlined in the Broad's [Data Format wiki](https://software.broadinstitute.org/cancer/software/gsea/wiki/index.php/Data_formats#Phenotype_Data_Formats). The order of the samples MUST be the same as the columns in the expression dataset! An example `cls` file is included here for the TP53 Base Editor experiment.

4. Load data to GSEA Desktop     
![GSEA Load File](https://user-images.githubusercontent.com/7750862/236883733-7bd32196-3c65-457a-879f-3c3862fe7a2e.png)

5. Run GSEA 
* **Expression dataset** - the normalized read counts loaded in step 4
* **Gene sets database** - preloaded MSigDB gene sets
* **Number of permutations** - for assessing statistical significance, default is 1000
* **Phenotype labels** - Comparison that you want to make (eg test vs control) using the labels from the cls file uploaded in step 4. 
* **Collapse / Remap to gene symbols** - use Collapse if you supplied gene IDs (recommended)
* **Permutation type**  - for assessing statistical significance
	* Phenotype permutation - randomly shuffles sample labels and ranks genes and calculates ES score, 7+ replicates per phenotype recommended
	* Gene set permutation - how likely it is that a random gene set of a given size was to be enriched within the given dataset, recommended for <7 replicates (this is the `fGsea` default)
* **Chip platform** - for mapping gene IDs to symbols, use `Human_ENSEMBL_Gene_ID_MSigDB.vX.chip`. See details on [wiki](https://software.broadinstitute.org/cancer/software/gsea/wiki/index.php/RNA-Seq_Data_and_Ensembl_CHIP_files)

The other parameter worth noting is **Metric for ranking genes** in Basic fields. The Broad default is to use Signal2Noise ratio. 

![GSEA Parameters](https://user-images.githubusercontent.com/7750862/236887212-1ae9ca00-63ce-484a-9028-bafda603f9cd.png)
