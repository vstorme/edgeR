#' ---
#' title: "case study 4: one group with 4 levels, intercept versus no-intercept model"
#' author: "Veronique Storme"
#' date: "14/11/2018"
#' output: html_document
#' ---
#' 
#' 
#' bioconductor must have been installed as well as the edgeR package
#' source("http://bioconductor.org/biocLite.R")
#' biocLite("BiocUpgrade")
#' biocLite("edgeR")
#' biocLite("limma")
#' 
#' The locfit package is also needed
#' biocLite("locfit")
#' 
#' other packages are needed for plotting and clustering, eg ggplo2
#' install.packages("ggplot2")
#' 
#' load the necessary libraries
#' 
## ----libraries, results = "hide", message = FALSE------------------------
library(knitr)
library(limma)
library(edgeR)
library(locfit)
library(statmod)
library(Rcpp)
library(pheatmap)
library(gplots)

#' 
#' getting the version and correct citations
#' 
## ----versions, results = "hide"------------------------------------------
sessionInfo()
citation("base")
citation("edgeR")

#' 
#' ### reading the data
#' 
#' start by setting the working directory
#' 
## ------------------------------------------------------------------------
setwd("//client/C$/MYCOURSES/RNAseq_v2/DEMO/data4")

#' 
#' all the counts are available in one file: data4.txt
#' create a DGEList object from it
#' 
## ------------------------------------------------------------------------
x <- read.delim("data4.txt", row.names="geneID")
head(x)
group <- factor(c(1,1,1,2,2,2,3,3,3,4,4,4))
y <- DGEList(counts=x,group=group)
y$samples
head(y$counts)

#' 
#' ### explore the data
#' 
## ------------------------------------------------------------------------
hist(asinh(y$counts))

#' 
## ------------------------------------------------------------------------
barplot(y$samples$lib.size*1e-6, names=1:12, ylab="Library size (millions)",cex.names=0.5)

#' 
#' 
#' ### filter low expression tags (little power to detect DE for filtered tags)
#' 
#' decide on a ad hoc cutoff of 5 read counts
#' the expression must exceed 5 read counts in at least 3 samples  (3 replicate samples in each group, hence 2)
#' libsize divided by 1 million gives an idea about the relationship between raw count and cpm
#' thus 5 read counts correspond to 5 divided by the libsize and multiplied by a million
#' 
## ------------------------------------------------------------------------
prop=y$samples$lib.size/1000000
prop
cutoff = 5/prop
cutoff

#' 
#' filter on 0.4 cpm
#' 
## ------------------------------------------------------------------------
keep = rowSums(cpm(y)> 0.4) >=3
y = y[keep,]
dim(y)

#' 
#' reset library sizes
#' 
## ------------------------------------------------------------------------
y$samples$lib.size = colSums(y$counts)
y$samples

#' 
#' ### apply TMM normalisation (TMM is the default method)
#' 
## ------------------------------------------------------------------------
y = calcNormFactors(y)
y$samples

#' 
#' multidimensional scaling plot
#' 
## ------------------------------------------------------------------------
plotMDS(y)

#' 
#' ### construct the design matrix (assuming no batch effect)
#' 
#' Before we fit GLMs, we need to define the design matrix based on the experimental design.
#' There are always many equivalent ways to define this matrix. The traditional way is to include an intercept term that represents the first level of the factor.
#' 
## ------------------------------------------------------------------------
design <- model.matrix(~ group, data=y$samples)
rownames(design) <- colnames(y)
design

#' 
#' the regression equation becomes:
#' log(E(Y)) ~ b0 + b1*A + b2*B + b3*C
#' 
#' WT is the ref level
#' 
#' predicted mean for WT: log(E(Y|WT)) = b0
#' predicted mean for mutant A: log(E(Y|A)) = b0 + b1*A, hence b1 = log(E(Y|A)) - log(E(Y|WT)) = logFC(A/WT)
#' predicted mean for mutant B: log(E(Y|B)) = b0 + b2*B, hence b2 = log(E(Y|B)) - log(E(Y|WT)) = logFC(B/WT)
#' predicted mean for mutant C: log(E(Y|C)) = b0 + b3*C, hence b3 = log(E(Y|C)) - log(E(Y|WT)) = logFC(C/WT)
#' 
#' 
#' estimating the dispersions
#' 
## ------------------------------------------------------------------------
y = estimateDisp(y, design, robust=TRUE)

#' 
#' Plot the tagwise BCV against its abundance (in log2 counts per million logCPM)
#' 
## ------------------------------------------------------------------------
plotBCV(y,cex=0.4)

#' 
#' ### fit the GLM
#' 
## ------------------------------------------------------------------------
fit.data4 = glmQLFit(y,design, robust=TRUE)


### Save the y object and fitted object for later use
#-----------------------------------------------------
save(y, fit.data4, file="data4_objects.RData")

#' 
#' ### Differential expression: 
#' 
#' research question is here to find DE between each mutant line and the ref line
#' 
#' For mutant A:
#' 
#' H0: log(E(Y/efflib)|A) - log(E(Y/efflib)|WT) = 0
#' 
#' log(E(Y/efflib)|A) = b0 + b1
#' log(E(Y/efflib)|WT) = b0
#' 
#' thus test whether b1=0
#'   
## ------------------------------------------------------------------------
AvsWT=glmQLFTest(fit.data4,coef=2)
results.AvsWT =topTags (AvsWT,n=Inf, sort.by="none")
save(results.AvsWT, file="data4_AvsWT.RData")
write.csv(results.AvsWT$table, file="results_AvsWT.csv")

#' 
#' table of the top differentially expressed tags (default n=10)
#' default adjustment for multiple testing with the BH method
#' sorted by (unadjusted) pvalue (default) or by absolute log-fold change
#' 
#' some stat info for the interested reader:
#' F is the value of the test-statistic, this test-statistic follows an F distribution
#' empirical Bayes quasi-likelihood F-tests
#'   
#' The p-values from glmQLFTest are always greater than or equal to those that would be obtained from glmLRT using the same negative binomial dispersions.
#'   
#' The lrt test compares the full model with the reduced model (b1=0)
#' since the difference between the full model and the reduced model is only 1df (the model differs in 1 parameter), 
#' the corresponding p-values are calculated as "pchisq(q=LR,df=1, lower.tail=FALSE)"
#' 
#' logFC is log(E(A)) - log(E(WT))
#' 
#' comparing B vs WT proceeds the same way
#' 
## ------------------------------------------------------------------------
BvsWT=glmQLFTest(fit.data4,coef=3)
results.BvsWT =topTags (BvsWT,n=Inf, sort.by="none")
write.csv(results.BvsWT$table, file="results_BvsWT.csv")

#' 
#' comparing C vs WT proceeds the same way
#' 
## ------------------------------------------------------------------------
CvsWT=glmQLFTest(fit.data4,coef=4)
results.CvsWT =topTags (CvsWT,n=Inf, sort.by="none")
write.csv(results.CvsWT$table, file="results_CvsWT.csv")

#' 
#' ### Venndiagrams{limma} based on decideTestsDGE (default adj BH and selection on FDR 0.05)
#' 
## ------------------------------------------------------------------------
de.AvsWT = decideTestsDGE(AvsWT)
de.BvsWT = decideTestsDGE(BvsWT)
de.CvsWT = decideTestsDGE(CvsWT)
sum(de.AvsWT == 1)
sum(de.BvsWT == 1)
sum(de.CvsWT == 1)
sum(de.AvsWT == -1)
sum(de.BvsWT == -1)
sum(de.CvsWT == -1)
de = as.matrix(cbind(de.AvsWT,de.BvsWT,de.CvsWT))
rownames(de) = rownames(y$counts)
head(de)

#' 
## ------------------------------------------------------------------------
vennDiagram(de, names=c("A vs WT", "B vs WT", "C vs WT"))

#' 
#' 
#' ### selecting those that have a significant two-fold change between C and ref
#' 
## ------------------------------------------------------------------------
twofold.CvsWT = glmTreat(fit.data4, coef=4, lfc=1)

DE.treat = topTags(twofold.CvsWT, n=Inf,sort.by="none")

select.treat = DE.treat$table[DE.treat$table$FDR<0.05,]
dim(select.treat)

#' 
#' ! Note that multiple testing is done on each list separataley.
#' 
#' classical ANOVA approach:
#' evaluate first F-test
#' only apply pairwise comparisons where F-test was sign (post-hoc tests)
#' 
#' ### EXTRA: fitted values
#' 
#' The fitted values differ between the 3 biological replicates because of the offset
#' To obtain fitted values normalised to the effective library size and resized to 1 million cpm:
#' 
## ------------------------------------------------------------------------
names(fit.data4)
fit.data4$coefficients[1:3,]
fit.data4$fitted.values[1:3,]

#' 
#' 
## ------------------------------------------------------------------------
N <- dim(y$counts)[[1]]
gene.fitted <- matrix(rep(NA,N*12),nrow=N)
for (i in 1:N)
  { 
  beta <- as.matrix(fit.data4$coefficients[i,])
  gene.fitted[i,] = exp(t(design %*% beta))*1000000
  }
  
row.names(gene.fitted) <- row.names(y$counts)
colnames(gene.fitted)  <- colnames(y$counts)
head(gene.fitted)

