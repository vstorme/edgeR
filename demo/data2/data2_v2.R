#' ---
#' title: "case study 2: two group comparison and testing for batch effect"
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
#' this RNAseq experiment involves a control-treatment case
#' there are 3 independent biological samples for each treatment group
#' 
#' ### reading the data
#' 
#' start by setting the working directory
#' 
## ------------------------------------------------------------------------
setwd("//client/C$/MYCOURSES/RNAseq_v2/DEMO/data2")

#' 
#' counts  are stored in 6 separate plain text files (with no headers)
#' 
## ------------------------------------------------------------------------
dir()

#' 
#' the file target_data2.txt gives the filename, the group and a brief description for each sample
#' 
## ------------------------------------------------------------------------
targets = readTargets("target_data2.txt")
targets

#' 
#' read the tables of counts, calculate the sizes of the libraries and produce a DGElist object
#' 
## ------------------------------------------------------------------------
d = readDGE(targets, header=FALSE)
d$samples
dim(d)

#' 
#' ### explore the data
#' 
## ------------------------------------------------------------------------
hist(asinh(d$counts))

#' 
## ------------------------------------------------------------------------
barplot(d$samples$lib.size*1e-6, names=1:6, ylab="Library size (millions)")

#' 
#' 
#' ### filter low expression tags
#' 
#' decide on a ad hoc cutoff of 5 read counts
#' 
## ------------------------------------------------------------------------
prop=d$samples$lib.size/1000000
prop
cutoff = 5/prop
cutoff

#' 
#' filter on 0.2 cpm
#' 
## ------------------------------------------------------------------------
keep = rowSums(cpm(d)> 0.2) >=3
d = d[keep,]
dim(d)

#' 
#' reset library sizes
#' 
## ------------------------------------------------------------------------
d$samples$lib.size = colSums(d$counts)
d$samples

#' 
#' ### apply TMM normalisation (TMM is the default method)
#' 
## ------------------------------------------------------------------------
d = calcNormFactors(d)
d$samples

#' 
#' MD plots
#' 
## ------------------------------------------------------------------------
for (i in 1:6)
  {
  jpeg(file=paste("MDplot_sample", i, ".jpg", sep=""),width = 10, height = 10, 
       units = "cm", res = 200) 
  plotMD(d, column = i)
  abline(h=0, col="red", lty=2, lwd=2)
  dev.off()
  }

#' 
#' ### multidimensional scaling plot
#' 
## ------------------------------------------------------------------------
plotMDS(d,labels=d$samples$description, cex=0.5)

#' 
#' ### construct the design matrix
#' 
#' Before we fit GLMs, we need to define the design matrix based on the experimental design.
#' There are always many equivalent ways to define this matrix. The traditional way is to include an intercept term that represents the first level of the factor.
#' 
#' include also a batch effect
#' 
## ------------------------------------------------------------------------
d$samples
d$samples$time = c("day1","day2","day3","day1","day2","day3")

design <- model.matrix(~ time + group, data=d$samples)
rownames(design) <- d$samples$group
design

#' 
#' the regression equation becomes:
#' log(E(Y)) ~ b0 + b1*timeday2 + b2*timeday3 +b3*grouptreated
#' 
#' The intercept is the expression level of the control treatment at time1
#' 
#' ### estimating the dispersions
#' 
#' > robust
#' 
## ------------------------------------------------------------------------
dr = estimateDisp(d, design, robust=TRUE)
quantile(dr$prior.df,c(0,0.001, 0.005, 0.01,0.05,0.5,0.95,0.99,0.995, 1))

#' 
#' outliers with a prior.df value < 1
#' 
## ------------------------------------------------------------------------
outliers = dr$counts[dr$prior.df < 1,]
dim(outliers)

#' 
#' few outliers detected
#' 
#' Plot the tagwise BCV against its abundance (in log2 counts per million logCPM)
#' 
## ------------------------------------------------------------------------
plotBCV(dr,cex=0.4)

#' 
#' ### fitting a quasi-likelihood model
#' 
## ------------------------------------------------------------------------
fit.data2 = glmQLFit(dr, design, robust=TRUE)
head(fit.data2$coefficients)

#' 
#' note, all fold changes are stored internally in natural log, meaning that fit$coefficient values are natural log, but all fold changes that are presented in output tables or stored in fit$table components are log2.
#' 
#' the coefficient b3 (grouptreated) represents the change in ln counts normalised to the effective library size between the control group and the treated group when time is held constant,
#' exp(b3) is the change in counts.
#' (log2(exp(b3))) is thus the log2FC
#' 
#' ### testing significance batch effect
#' 
#' test whether b1=b2=0 (ie is there need to adjust for timeday3 )
#' LR: compares the full model log(E(Y)) ~ b0 + b1*timeday2 + b2*timeday3 +b3*grouptreated
#' with the model log(E(Y)) ~ b0 +b3*grouptreated
#' 
## ------------------------------------------------------------------------
DE.batch=glmQLFTest(fit.data2,coef=2:3)

#' 
## ------------------------------------------------------------------------
DE.batch.fdr = topTags(DE.batch,n=Inf,sort.by="none")
head(DE.batch.fdr$table)

#' 
#' number of tags where FDR < 0.05
#' 
## ------------------------------------------------------------------------
sigbatch <- DE.batch.fdr$table[DE.batch.fdr$table$FDR < 0.05,]
sum(sigbatch$FDR < 0.05)

#' 
#' Conclusion: there is a need to correct for the batch (time) effects
#' 
#' ### Differential expression due to treatment effect: 
#' 
#' research question is here to find DE between each mutant line and the ref line
#' 
#' test whether b1=0 ie H0: log(E(Y|control)) = log(E(Y|treated))
#' LR: compares the full model log(E(Y)) ~ b0 + b1*timeday2 + b2*timeday3 +b3*grouptreated
#' with the reduced model log(E(Y)) ~ b0 + b1*timeday2 + b2*timeday3
#'   
## ------------------------------------------------------------------------
DE.TvsC=glmQLFTest(fit.data2,coef=4)

#' 
#' 
## ------------------------------------------------------------------------
DE.fdr = topTags(DE.TvsC,n=Inf,sort.by="none")

#' 
#' export as csv file
#' 
## ------------------------------------------------------------------------
write.csv(DE.fdr$table,file="de_fdr_logFC.csv")

#' 
#' save the object
#' 
## ------------------------------------------------------------------------
save(DE.fdr, file="data2_treated_vs_control.RData")

#' 
#' > use of decideTestsDGE function
#' 
#' total nr of DE genes at FDR < 0.05 is and FC of at least 2
#' 
## ------------------------------------------------------------------------
de <- decideTestsDGE(DE.TvsC, p=0.05, adjust="BH", lfc=1)
summary(de)

#' 
#' > test which genes show a DE significantly above a FC of 2
#' 
## ------------------------------------------------------------------------
tr = glmTreat(fit.data2, coef=4, lfc=1)

DE.treat = topTags(tr, n=Inf,sort.by="none")

select.treat = DE.treat$table[DE.treat$table$FDR<0.05,]
dim(select.treat)

#' 
#' compare with an ad hoc selection
#' 
## ------------------------------------------------------------------------
select.adhoc = DE.fdr$table[DE.fdr$table$FDR < 0.05 & abs(DE.fdr$table$logFC) > 1,]
dim(select.adhoc)

#' 
#' export as csv file
#' 
## ------------------------------------------------------------------------
write.csv(DE.fdr$table,file="de_fdr_logFC.csv")
write.csv(DE.treat$table,file="de_treat_logFC.csv")

#' 
