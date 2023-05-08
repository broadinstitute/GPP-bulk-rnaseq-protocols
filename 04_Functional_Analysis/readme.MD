# Functional Analysis

## Gene Set Enrichment Analysis
GSEA refers to a type of analysis that groups together genes to find statistically significant pathways effected in a given phenotype. While there are various libraries that will accomplish this, the GSEA protocols in this folder focuses on the `fgsea` R package and a Desktop tool, both of which implement the original Broad algorithm by [Subramanian et al](https://www.pnas.org/doi/10.1073/pnas.0506580102).

While these guides are specific to the Broad version, they will help you understand the parameters and algorithms of GSEA in general. 
* GSEA [User Guide](https://www.gsea-msigdb.org/gsea/doc/GSEAUserGuideFrame.html)
* GSEA [Wiki](https://software.broadinstitute.org/cancer/software/gsea/wiki/index.php/Main_Page)

### GSEA Desktop
Although the [R version](https://github.com/GSEA-MSigDB/GSEA_R) of the Broad software exists, it is not actively maintained by the developers and there are untested portions of the code. See GSEA_desktop.md for instructions on how to use the desktop version maintained by the Broad.

### fgsea (pre-ranked)
The [fast-GSEA]((https://www.biorxiv.org/content/10.1101/060012v3)) R package is not the exact same implementation as the Broad algorithm, but the enrichment scores and normalized enrichment scores are calculated in the same way. See GSEA.md for the R protocol.