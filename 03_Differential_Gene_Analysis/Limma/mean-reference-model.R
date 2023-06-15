library("limma")

##########################
# GUIDE EFFECT UNTREATED #
##########################
group.untreated <- relevel(group, "untreated_ctrl")
design.intercept.untreated <- model.matrix(~group.untreated)
colnames(design.intercept.untreated) <- gsub("group.untreated", "", colnames(design.intercept.untreated))

v.intercept.untreated <- voom(reads.dge, design.intercept.untreated, plot=TRUE)
fit.intercept.untreated <- lmFit(v.intercept.untreated,design.intercept.untreated)
fit.intercept.untreated <- eBayes(fit.intercept.untreated)
summary(decideTests(fit.intercept.untreated))

#(Intercept) treated_ctrl treated_ex1 treated_ex9 untreated_ex1 untreated_ex9
#Down          1143         3647        1586         748          2090           963
#NotSig        1865         8298       13266       14951         12454         14629
#Up           12921         3984        1077         230          1385           337

toptable.intercept.untreated_ex1 <- makeTopTableWithGeneNames(fit.intercept.untreated, colNum=5)
toptable.intercept.untreated_ex9 <- makeTopTableWithGeneNames(fit.intercept.untreated, colNum=6)



########################
# GUIDE EFFECT TREATED #
########################
design.intercept.treated <- model.matrix(~group)
colnames(design.intercept.treated) <- gsub("group", "", colnames(design.intercept.treated))

v.intercept.treated <- voom(reads.dge, design.intercept.treated, plot=TRUE)
fit.intercept.treated <- lmFit(v.intercept.treated,design.intercept.treated)
fit.intercept.treated <- eBayes(fit.intercept.treated)
summary(decideTests(fit.intercept.treated))
#(Intercept) treated_ex1 treated_ex9 untreated_ctrl untreated_ex1 untreated_ex9
#Down           669        3791        4287           3984          4157          4773
#NotSig        1940        8449        7648           8298          7315          6924
#Up           13320        3689        3994           3647          4457          4232
toptable.intercept.treated_ex1 <- makeTopTableWithGeneNames(fit.intercept.treated, colNum=2)
toptable.intercept.treated_ex9 <- makeTopTableWithGeneNames(fit.intercept.treated, colNum=3)


pdf(file = paste(analysis_dir, "DEG-MeanReference.pdf"))
p1<-meanDifferencePlot(toptable.intercept.untreated_ex1, "Untreated Ex 1 vs Untreated Control") 
p2<-meanDifferencePlot(toptable.intercept.untreated_ex9, "Untreated Ex 9 vs Untreated Control")  
p3<-meanDifferencePlot(toptable.intercept.treated_ex1, "Nutlin Ex 9 vs Nutlin Control") 
p4<-meanDifferencePlot(toptable.intercept.treated_ex9, "Nutlin Ex 9 vs Nutlin Control")
ggarrange(p1, p2, p3, p4, ncol=2, nrow=2, common.legend = TRUE, legend="bottom")
dev.off() 




