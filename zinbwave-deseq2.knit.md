---
title: "ZINB-WaVE + DESeq2 integration for single-cell RNA-seq"
author: "Michael Love"
output: html_document
---

Here we use the *splatter* package to simulate single-cell RNA-seq
data.

* Zappia, Phipson, and Oshlack "Splatter: simulation of single-cell RNA
sequencing data" *Genome Biology* (2017)
[doi: 10.1186/s13059-017-1305-0](https://doi.org/10.1186/s13059-017-1305-0)

We then use the methods defined in the following paper to combine
*zinbwave* observation weights with *DESeq2* modeling of negative
binomial counts.

* Van den Berge & Perraudeau *et al* "Observation weights unlock bulk
RNA-seq tools for zero inflation and single-cell applications" *Genome Biology* (2018)
[doi: 10.1186/s13059-018-1406-4](https://doi.org/10.1186/s13059-018-1406-4)

> It is important to note that while methods such as ZINB-WaVE and
> ZINGER can successfully identify excess zeros, they cannot, however,
> readily discriminate between their underlying causes, i.e., between
> technical (e.g., dropout) and biological (e.g., bursting) zeros. 

The above note implies that the zero-inflation weighting approach
outlined below can be used when the interesting signal is not in the
zero component. That is, if you wanted to find biological differences
in transcriptional bursting across groups of cells, the below approach
would not help you find these differences. It instead helps to uncover
differences in counts besides the zero component (whether those zeros
be biological or technical). 

### Simulate single-cell count data with *splatter*


```r
suppressPackageStartupMessages(library(splatter))
params <- newSplatParams()
#params
#slotNames(params)
# note: these DE params are natural log scale
params <- setParam(params, "de.facLoc", 1) 
params <- setParam(params, "de.facScale", .25)
# add a lot more dropout - see if ZI weighting works
params <- setParam(params, "dropout.type", "experiment")
params <- setParam(params, "dropout.mid", 3)
```


```r
set.seed(1)
sim <- splatSimulate(params, group.prob=c(.5,.5), method="groups")
```

```
## Getting parameters...
```

```
## Creating simulation object...
```

```
## Simulating library sizes...
```

```
## Simulating gene means...
```

```
## Simulating group DE...
```

```
## Simulating cell means...
```

```
## Simulating BCV...
```

```
## Simulating counts...
```

```
## Simulating dropout (if needed)...
```

```
## Sparsifying assays...
```

```
## Automatically converting to sparse matrices, threshold = 0.95
```

```
## Skipping 'BatchCellMeans': estimated sparse size 1.5 * dense matrix
```

```
## Skipping 'BaseCellMeans': estimated sparse size 1.5 * dense matrix
```

```
## Skipping 'BCV': estimated sparse size 1.5 * dense matrix
```

```
## Skipping 'CellMeans': estimated sparse size 1.49 * dense matrix
```

```
## Skipping 'TrueCounts': estimated sparse size 1.62 * dense matrix
```

```
## Skipping 'DropProb': estimated sparse size 1.5 * dense matrix
```

```
## Warning in sparsifyMatrices(assays(sim), auto = TRUE, verbose = verbose): matrix 'Dropout' is class
## 'matrixarray', unable to estimate size reduction factor
```

```
## Converting 'Dropout' to sparse matrix: estimated sparse size NA * dense matrix
```

```
## Converting 'counts' to sparse matrix: estimated sparse size 0.34 * dense matrix
```

```
## Done!
```

```r
# hist(rowSums(counts(sim) >= 1))
```


```r
plot(log10(rowMeans(assays(sim)[["TrueCounts"]])), rowMeans(assays(sim)[["Dropout"]]))
```

<img src="zinbwave-deseq2_files/figure-html/dropout-1.png" width="672" />


```r
# note: each group gets it's own DE genes -- meaning some will be "doubly DE"
#z <- rowData(sim)$DEFacGroup1
#hist(log(z[z > 1]), breaks=30, col="grey", freq=FALSE, ylim=c(0,5))
rowData(sim)$log2FC <- with(rowData(sim), log2(DEFacGroup2/DEFacGroup1))
```


```r
rowData(sim)$trueDisp <- rowMeans(assays(sim)[["BCV"]])^2
gridlines <- c(1e-2,1e-1,1); cols <- c("blue","red","darkgreen")
with(rowData(sim)[rowData(sim)$GeneMean> 1,],
     plot(GeneMean, trueDisp, log="xy", xlim=c(1,300), ylim=c(.01,5)))
abline(h=gridlines, col=cols)
text(300, gridlines, labels=gridlines, col=cols, pos=3)
```

<img src="zinbwave-deseq2_files/figure-html/trueDisp-1.png" width="672" />

### Model zero component using *zinbwave*


```r
library(zinbwave)
library(BiocParallel)
# low count filter - at least 10 with count of 5 or more
keep <- rowSums(counts(sim) >= 5) >= 10
table(keep)
```

```
## keep
## FALSE  TRUE 
##  7891  2109
```

```r
zinb <- sim[keep,]
zinb$condition <- factor(zinb$Group)
# we need to reorganize the assays in the SumExp from splatter
nms <- c("counts", setdiff(assayNames(zinb), "counts"))
assays(zinb) <- assays(zinb)[nms]
assay(zinb) <- as.matrix(assay(zinb))
# epsilon setting as recommended by the ZINB-WaVE integration paper
system.time({
  zinb <- zinbwave(zinb, K=0, observationalWeights=TRUE,
                   BPPARAM=SerialParam(), epsilon=1e12)
})
```

```
##    user  system elapsed 
## 177.506   3.066 195.533
```

### Estimate size factors


```r
suppressPackageStartupMessages(library(DESeq2))
dds <- DESeqDataSet(zinb, design=~condition)
```

```
## converting counts to integer mode
```

```r
dds <- estimateSizeFactors(dds, type="poscounts")
library(scran)
scr <- computeSumFactors(dds)
dat <- data.frame(true=dds$ExpLibSize,
                  pos=sizeFactors(dds),
                  sum=sizeFactors(scr))
dat$true <- dat$true / exp(mean(log(dat$true)))
panel.scatter <- function(x,y,...) {
  points(x,y,...)
  abline(0,1,col="red",lwd=2)
  legend("topleft", legend=round(cor(x,y),3))
}
pairs(dat, panel=panel.scatter)
```

<img src="zinbwave-deseq2_files/figure-html/unnamed-chunk-5-1.png" width="672" />

### Estimate dispersion and DE using *DESeq2*

Van den Berge and Perraudeau and others have shown the LRT may perform
better for null hypothesis testing, so we use the LRT. In order to use
the Wald test, it is recommended to set `useT=TRUE`.


```r
# use scran's sum factors:
sizeFactors(dds) <- sizeFactors(scr)
# run DESeq:
system.time({
  dds <- DESeq(dds, test="LRT", reduced=~1,
               minmu=1e-6, minRep=Inf)
})
```

```
## using pre-existing size factors
```

```
## estimating dispersions
```

```
## gene-wise dispersion estimates
```

```
## Warning in getAndCheckWeights(object, modelMatrix, weightThreshold = weightThreshold): for 1 row(s), the weights as supplied won't allow parameter estimation, producing a
##   degenerate design matrix. These rows have been flagged in mcols(dds)$weightsFail
##   and treated as if the row contained all zeros (mcols(dds)$allZero set to TRUE).
##   If you are blocking for donors/organisms, consider design = ~0+donor+condition.
```

```
## mean-dispersion relationship
```

```
## final dispersion estimates
```

```
## fitting model and testing
```

```
##    user  system elapsed 
##  16.422   0.592  18.723
```

### Plot dispersion estimates

It is recommended to plot the dispersion estimates for *DESeq2* on
single-cell data. As discussed in the *DESeq2* paper, it becomes
difficult to accurately estimate the dispersion when the counts are
very small, because the Poisson component of the variance is
dominant. Therefore we see some very low dispersion estimates here,
although the trend is still accurately capturing the upper proportion.
So here everything looks good.


```r
plotDispEsts(dds)
```

<img src="zinbwave-deseq2_files/figure-html/plotDispEsts-1.png" width="672" />

If the parametric trend fails to fit (there would be a warning in this
case), one should check the dispersion plot as above. If it looks like
the dispersion fit is being thrown off by the low count genes with low
dispersion estimates at the bottom of the plot, there is a relatively
easy solution: one can filter out more of the low count genes only for
the dispersion estimation step, so that the trend still captures the upper
portion. This is pretty easy to do in *DESeq2*, to filter genes solely
for the dispersion trend estimation, but to use a larger set for the
rest of the analysis. An example of how this can be done:


```r
keepForDispTrend <- rowSums(counts(dds) >= 10) >= 25
dds2 <- estimateDispersionsFit(dds[keepForDispTrend,])
plotDispEsts(dds2, ylim=c(1e-3,1))
```

<img src="zinbwave-deseq2_files/figure-html/plotDispEsts2-1.png" width="672" />

One would then assign the dispersion function to the original dataset,
re-estimate final dispersions, check `plotDispEsts`, and then either
re-run the Wald or LRT function (this chunk not evaluated):


```r
dispersionFunction(dds) <- dispersionFunction(dds2)
dds <- estimateDispersionsMAP(dds)
dds <- nbinomLRT(dds, reduced=~1, minmu=1e-6)
```

### Evaluate how well we did on simulated data

Compare dispersion on the non-zero-component counts to the true value
used for simulation. 


```r
with(mcols(dds), plot(trueDisp, dispMAP, log="xy"))
abline(0,1,col="red")
```

<img src="zinbwave-deseq2_files/figure-html/trueDispVsMAP-1.png" width="672" />

Extract results table.


```r
# we already performed low count filtering
res <- results(dds, independentFiltering=FALSE)
plot(mcols(dds)$log2FC, res$log2FoldChange, ylim=c(-4,4)); abline(0,1,col="red")
```

<img src="zinbwave-deseq2_files/figure-html/trueLFCVsMLE-1.png" width="672" />

Below we show that the "simple" LFC does not work - it over-estimates
the true DE LFC because of the dropout zeros in the group with the
smaller mean. It also has a lot of noise for the null genes.


```r
ncts <- counts(dds, normalized=TRUE)
simple.lfc <- log2(rowMeans(ncts[,dds$condition == "Group2"])/
                   rowMeans(ncts[,dds$condition == "Group1"]))
plot(mcols(dds)$log2FC, simple.lfc, ylim=c(-4,4)); abline(0,1,col="red")
```

<img src="zinbwave-deseq2_files/figure-html/trueLFCVsSimple-1.png" width="672" />

How well do we do in null hypothesis testing:


```r
tab <- table(DE.status=mcols(dds)$log2FC != 0, sig=res$padj < .05)
tab
```

```
##          sig
## DE.status FALSE TRUE
##     FALSE  1607   26
##     TRUE    172  303
```

```r
round(prop.table(tab, 2), 3)
```

```
##          sig
## DE.status FALSE  TRUE
##     FALSE 0.903 0.079
##     TRUE  0.097 0.921
```


```r
session_info()
```

```
## ─ Session info ───────────────────────────────────────────────────────────────────────────────────
##  setting  value                       
##  version  R version 4.0.3 (2020-10-10)
##  os       macOS Catalina 10.15.6      
##  system   x86_64, darwin17.0          
##  ui       X11                         
##  language (EN)                        
##  collate  en_US.UTF-8                 
##  ctype    en_US.UTF-8                 
##  tz       America/New_York            
##  date     2021-01-06                  
## 
## ─ Packages ───────────────────────────────────────────────────────────────────────────────────────
##  package              * version  date       lib source        
##  annotate               1.68.0   2020-10-27 [1] Bioconductor  
##  AnnotationDbi          1.52.0   2020-10-27 [1] Bioconductor  
##  assertthat             0.2.1    2019-03-21 [1] CRAN (R 4.0.0)
##  backports              1.2.1    2020-12-09 [1] CRAN (R 4.0.2)
##  beachmat               2.6.4    2020-12-20 [1] Bioconductor  
##  Biobase              * 2.50.0   2020-10-27 [1] Bioconductor  
##  BiocGenerics         * 0.36.0   2020-10-27 [1] Bioconductor  
##  BiocNeighbors          1.8.2    2020-12-07 [1] Bioconductor  
##  BiocParallel         * 1.24.1   2020-11-06 [1] Bioconductor  
##  BiocSingular           1.6.0    2020-10-27 [1] Bioconductor  
##  bit                    4.0.4    2020-08-04 [1] CRAN (R 4.0.2)
##  bit64                  4.0.5    2020-08-30 [1] CRAN (R 4.0.2)
##  bitops                 1.0-6    2013-08-17 [1] CRAN (R 4.0.0)
##  blob                   1.2.1    2020-01-20 [1] CRAN (R 4.0.0)
##  bluster                1.0.0    2020-10-27 [1] Bioconductor  
##  callr                  3.5.1    2020-10-13 [1] CRAN (R 4.0.2)
##  checkmate              2.0.0    2020-02-06 [1] CRAN (R 4.0.0)
##  cli                    2.2.0    2020-11-20 [1] CRAN (R 4.0.2)
##  colorspace             2.0-0    2020-11-11 [1] CRAN (R 4.0.2)
##  crayon                 1.3.4    2017-09-16 [1] CRAN (R 4.0.0)
##  DBI                    1.1.0    2019-12-15 [1] CRAN (R 4.0.0)
##  DelayedArray           0.16.0   2020-10-27 [1] Bioconductor  
##  DelayedMatrixStats     1.12.1   2020-11-24 [1] Bioconductor  
##  desc                   1.2.0    2018-05-01 [1] CRAN (R 4.0.0)
##  DESeq2               * 1.30.0   2020-10-27 [1] Bioconductor  
##  devtools             * 2.3.2    2020-09-18 [1] CRAN (R 4.0.2)
##  digest                 0.6.27   2020-10-24 [1] CRAN (R 4.0.2)
##  dplyr                  1.0.2    2020-08-18 [1] CRAN (R 4.0.2)
##  dqrng                  0.2.1    2019-05-17 [1] CRAN (R 4.0.0)
##  edgeR                  3.32.0   2020-10-27 [1] Bioconductor  
##  ellipsis               0.3.1    2020-05-15 [1] CRAN (R 4.0.2)
##  evaluate               0.14     2019-05-28 [1] CRAN (R 4.0.0)
##  fansi                  0.4.1    2020-01-08 [1] CRAN (R 4.0.0)
##  fs                     1.5.0    2020-07-31 [1] CRAN (R 4.0.2)
##  genefilter             1.72.0   2020-10-27 [1] Bioconductor  
##  geneplotter            1.68.0   2020-10-27 [1] Bioconductor  
##  generics               0.1.0    2020-10-31 [1] CRAN (R 4.0.2)
##  GenomeInfoDb         * 1.26.2   2020-12-08 [1] Bioconductor  
##  GenomeInfoDbData       1.2.4    2020-11-09 [1] Bioconductor  
##  GenomicRanges        * 1.42.0   2020-10-27 [1] Bioconductor  
##  ggplot2                3.3.3    2020-12-30 [1] CRAN (R 4.0.3)
##  glue                   1.4.2    2020-08-27 [1] CRAN (R 4.0.2)
##  gtable                 0.3.0    2019-03-25 [1] CRAN (R 4.0.0)
##  htmltools              0.5.0    2020-06-16 [1] CRAN (R 4.0.2)
##  httr                   1.4.2    2020-07-20 [1] CRAN (R 4.0.2)
##  igraph                 1.2.6    2020-10-06 [1] CRAN (R 4.0.2)
##  IRanges              * 2.24.1   2020-12-12 [1] Bioconductor  
##  irlba                  2.3.3    2019-02-05 [1] CRAN (R 4.0.0)
##  knitr                  1.30     2020-09-22 [1] CRAN (R 4.0.2)
##  lattice                0.20-41  2020-04-02 [1] CRAN (R 4.0.3)
##  lifecycle              0.2.0    2020-03-06 [1] CRAN (R 4.0.0)
##  limma                  3.46.0   2020-10-27 [1] Bioconductor  
##  locfit                 1.5-9.4  2020-03-25 [1] CRAN (R 4.0.0)
##  magrittr               2.0.1    2020-11-17 [1] CRAN (R 4.0.2)
##  Matrix                 1.3-0    2020-12-22 [1] CRAN (R 4.0.2)
##  MatrixGenerics       * 1.2.0    2020-10-27 [1] Bioconductor  
##  matrixStats          * 0.57.0   2020-09-25 [1] CRAN (R 4.0.2)
##  memoise                1.1.0    2017-04-21 [1] CRAN (R 4.0.0)
##  munsell                0.5.0    2018-06-12 [1] CRAN (R 4.0.0)
##  pillar                 1.4.7    2020-11-20 [1] CRAN (R 4.0.2)
##  pkgbuild               1.2.0    2020-12-15 [1] CRAN (R 4.0.2)
##  pkgconfig              2.0.3    2019-09-22 [1] CRAN (R 4.0.0)
##  pkgload                1.1.0    2020-05-29 [1] CRAN (R 4.0.2)
##  prettyunits            1.1.1    2020-01-24 [1] CRAN (R 4.0.0)
##  processx               3.4.5    2020-11-30 [1] CRAN (R 4.0.2)
##  ps                     1.5.0    2020-12-05 [1] CRAN (R 4.0.2)
##  purrr                  0.3.4    2020-04-17 [1] CRAN (R 4.0.0)
##  R6                     2.5.0    2020-10-28 [1] CRAN (R 4.0.2)
##  RColorBrewer           1.1-2    2014-12-07 [1] CRAN (R 4.0.0)
##  Rcpp                   1.0.5    2020-07-06 [1] CRAN (R 4.0.2)
##  RCurl                  1.98-1.2 2020-04-18 [1] CRAN (R 4.0.0)
##  remotes                2.2.0    2020-07-21 [1] CRAN (R 4.0.2)
##  rlang                  0.4.9    2020-11-26 [1] CRAN (R 4.0.2)
##  rmarkdown            * 2.6      2020-12-14 [1] CRAN (R 4.0.2)
##  rprojroot              2.0.2    2020-11-15 [1] CRAN (R 4.0.2)
##  RSQLite                2.2.1    2020-09-30 [1] CRAN (R 4.0.2)
##  rsvd                   1.0.3    2020-02-17 [1] CRAN (R 4.0.0)
##  S4Vectors            * 0.28.1   2020-12-09 [1] Bioconductor  
##  scales                 1.1.1    2020-05-11 [1] CRAN (R 4.0.0)
##  scran                * 1.18.3   2020-12-21 [1] Bioconductor  
##  scuttle                1.0.4    2020-12-17 [1] Bioconductor  
##  sessioninfo            1.1.1    2018-11-05 [1] CRAN (R 4.0.0)
##  SingleCellExperiment * 1.12.0   2020-10-27 [1] Bioconductor  
##  softImpute             1.4      2015-04-08 [1] CRAN (R 4.0.0)
##  sparseMatrixStats      1.2.0    2020-10-27 [1] Bioconductor  
##  splatter             * 1.14.1   2020-12-01 [1] Bioconductor  
##  statmod                1.4.35   2020-10-19 [1] CRAN (R 4.0.2)
##  stringi                1.5.3    2020-09-09 [1] CRAN (R 4.0.2)
##  stringr                1.4.0    2019-02-10 [1] CRAN (R 4.0.0)
##  SummarizedExperiment * 1.20.0   2020-10-27 [1] Bioconductor  
##  survival               3.2-7    2020-09-28 [1] CRAN (R 4.0.3)
##  testthat             * 3.0.1    2020-12-17 [1] CRAN (R 4.0.2)
##  tibble                 3.0.4    2020-10-12 [1] CRAN (R 4.0.2)
##  tidyselect             1.1.0    2020-05-11 [1] CRAN (R 4.0.2)
##  usethis              * 2.0.0    2020-12-10 [1] CRAN (R 4.0.2)
##  vctrs                  0.3.6    2020-12-17 [1] CRAN (R 4.0.2)
##  withr                  2.3.0    2020-09-22 [1] CRAN (R 4.0.2)
##  xfun                   0.19     2020-10-30 [1] CRAN (R 4.0.2)
##  XML                    3.99-0.5 2020-07-23 [1] CRAN (R 4.0.2)
##  xtable                 1.8-4    2019-04-21 [1] CRAN (R 4.0.0)
##  XVector                0.30.0   2020-10-28 [1] Bioconductor  
##  yaml                   2.2.1    2020-02-01 [1] CRAN (R 4.0.0)
##  zinbwave             * 1.12.0   2020-10-28 [1] Bioconductor  
##  zlibbioc               1.36.0   2020-10-28 [1] Bioconductor  
## 
## [1] /Library/Frameworks/R.framework/Versions/4.0/Resources/library
```
