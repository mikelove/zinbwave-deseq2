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
## Simulating counts..
```

```
## Simulating dropout (if needed)...
```

```
## Done!
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
# low count filter - at least 25 samples with count of 5 or more
keep <- rowSums(counts(sim) >= 5) >= 25
table(keep)
```

```
## keep
## FALSE  TRUE 
##  9026   974
```

```r
zinb <- sim[keep,]
zinb$condition <- factor(zinb$Group)
# we need to reorganize the assays in the SumExp from splatter
nms <- c("counts", setdiff(assayNames(zinb), "counts"))
assays(zinb) <- assays(zinb)[nms]
# epsilon setting as recommended by the ZINB-WaVE integration paper
system.time({
  zinb <- zinbwave(zinb, K=0, BPPARAM=SerialParam(), epsilon=1e12)
})
```

```
##    user  system elapsed 
##  14.189   0.072  14.261
```

### Estimate dispersion and DE using *DESeq2*

Van De Berge and Perraudeau and others have shown the LRT may perform
better for null hypothesis testing, so we use the LRT. In order to use
the Wald test, it is recommended to set `useT=TRUE`.


```r
suppressPackageStartupMessages(library(DESeq2))
dds <- DESeqDataSet(zinb, design=~condition)
system.time({
  dds <- DESeq(dds, test="LRT", reduced=~1,
               sfType="poscounts", minmu=1e-6, minRep=Inf)
})
```

```
## estimating size factors
```

```
## estimating dispersions
```

```
## gene-wise dispersion estimates
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
##  19.722   0.003  19.725
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
plotDispEsts(dds2)
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
##     FALSE   715   18
##     TRUE     14  227
```

```r
round(prop.table(tab, 2), 3)
```

```
##          sig
## DE.status FALSE  TRUE
##     FALSE 0.981 0.073
##     TRUE  0.019 0.927
```


```r
session_info()
```

```
## Session info -------------------------------------------------------------
```

```
##  setting  value                                             
##  version  R Under development (unstable) (2017-12-25 r73962)
##  system   x86_64, linux-gnu                                 
##  ui       X11                                               
##  language (EN)                                              
##  collate  en_US.UTF-8                                       
##  tz       America/New_York                                  
##  date     2018-04-27
```

```
## Packages -----------------------------------------------------------------
```

```
##  package              * version   date       source        
##  acepack                1.4.1     2016-10-29 CRAN (R 3.5.0)
##  ADGofTest              0.3       2011-12-28 CRAN (R 3.5.0)
##  annotate               1.57.2    2017-12-27 Bioconductor  
##  AnnotationDbi          1.41.4    2017-12-27 Bioconductor  
##  backports              1.1.2     2017-12-13 CRAN (R 3.5.0)
##  base                 * 3.5.0     2017-12-27 local         
##  base64enc              0.1-3     2015-07-28 CRAN (R 3.5.0)
##  Biobase              * 2.39.2    2018-02-03 Bioconductor  
##  BiocGenerics         * 0.25.3    2018-03-20 Bioconductor  
##  BiocInstaller        * 1.29.6    2018-04-12 Bioconductor  
##  BiocParallel         * 1.13.1    2018-01-17 Bioconductor  
##  bit                    1.1-12    2014-04-09 CRAN (R 3.5.0)
##  bit64                  0.9-7     2017-05-08 CRAN (R 3.5.0)
##  bitops                 1.0-6     2013-08-17 CRAN (R 3.5.0)
##  blob                   1.1.0     2017-06-17 CRAN (R 3.5.0)
##  checkmate              1.8.5     2017-10-24 CRAN (R 3.5.0)
##  cluster                2.0.6     2017-03-10 CRAN (R 3.5.0)
##  codetools              0.2-15    2016-10-05 CRAN (R 3.5.0)
##  colorspace             1.3-2     2016-12-14 CRAN (R 3.5.0)
##  compiler               3.5.0     2017-12-27 local         
##  copula                 0.999-18  2017-09-01 CRAN (R 3.5.0)
##  data.table             1.10.4-3  2017-10-27 CRAN (R 3.5.0)
##  datasets             * 3.5.0     2017-12-27 local         
##  DBI                    0.8       2018-03-02 CRAN (R 3.5.0)
##  DelayedArray         * 0.5.31    2018-04-20 Bioconductor  
##  DESeq2               * 1.19.51   2018-04-27 Bioconductor  
##  devtools             * 1.13.5    2018-02-18 CRAN (R 3.5.0)
##  digest                 0.6.15    2018-01-28 CRAN (R 3.5.0)
##  edgeR                  3.21.9    2018-03-20 Bioconductor  
##  evaluate               0.10.1    2017-06-24 CRAN (R 3.5.0)
##  foreach                1.4.4     2017-12-12 CRAN (R 3.5.0)
##  foreign                0.8-70    2017-11-28 CRAN (R 3.5.0)
##  Formula                1.2-2     2017-07-10 CRAN (R 3.5.0)
##  genefilter             1.61.1    2018-02-03 Bioconductor  
##  geneplotter            1.57.0    2017-12-27 Bioconductor  
##  GenomeInfoDb         * 1.15.5    2018-03-20 Bioconductor  
##  GenomeInfoDbData       1.1.0     2017-12-27 Bioconductor  
##  GenomicRanges        * 1.31.22   2018-03-20 Bioconductor  
##  ggplot2                2.2.1     2016-12-30 CRAN (R 3.5.0)
##  glmnet                 2.0-13    2017-09-22 CRAN (R 3.5.0)
##  graphics             * 3.5.0     2017-12-27 local         
##  grDevices            * 3.5.0     2017-12-27 local         
##  grid                   3.5.0     2017-12-27 local         
##  gridExtra              2.3       2017-09-09 CRAN (R 3.5.0)
##  gsl                    1.9-10.3  2017-01-05 CRAN (R 3.5.0)
##  gtable                 0.2.0     2016-02-26 CRAN (R 3.5.0)
##  Hmisc                  4.1-1     2018-01-03 CRAN (R 3.5.0)
##  htmlTable              1.11.2    2018-01-20 CRAN (R 3.5.0)
##  htmltools              0.3.6     2017-04-28 CRAN (R 3.5.0)
##  htmlwidgets            1.0       2018-01-20 CRAN (R 3.5.0)
##  IRanges              * 2.13.28   2018-03-20 Bioconductor  
##  iterators              1.0.9     2017-12-12 CRAN (R 3.5.0)
##  knitr                  1.20      2018-02-20 CRAN (R 3.5.0)
##  lattice                0.20-35   2017-03-25 CRAN (R 3.5.0)
##  latticeExtra           0.6-28    2016-02-09 CRAN (R 3.5.0)
##  lazyeval               0.2.1     2017-10-29 CRAN (R 3.5.0)
##  limma                  3.35.13   2018-03-20 Bioconductor  
##  locfit                 1.5-9.1   2013-04-20 CRAN (R 3.5.0)
##  magrittr             * 1.5       2014-11-22 CRAN (R 3.5.0)
##  Matrix                 1.2-12    2017-11-20 CRAN (R 3.5.0)
##  matrixStats          * 0.53.1    2018-02-11 CRAN (R 3.5.0)
##  memoise                1.1.0     2017-04-21 CRAN (R 3.5.0)
##  methods              * 3.5.0     2017-12-27 local         
##  munsell                0.4.3     2016-02-13 CRAN (R 3.5.0)
##  mvtnorm                1.0-7     2018-01-26 CRAN (R 3.5.0)
##  nnet                   7.3-12    2016-02-02 CRAN (R 3.5.0)
##  numDeriv               2016.8-1  2016-08-27 CRAN (R 3.5.0)
##  parallel             * 3.5.0     2017-12-27 local         
##  pcaPP                  1.9-73    2018-01-14 CRAN (R 3.5.0)
##  pillar                 1.2.1     2018-02-27 CRAN (R 3.5.0)
##  plyr                   1.8.4     2016-06-08 CRAN (R 3.5.0)
##  pspline                1.0-18    2017-06-12 CRAN (R 3.5.0)
##  R6                     2.2.2     2017-06-17 CRAN (R 3.5.0)
##  RColorBrewer           1.1-2     2014-12-07 CRAN (R 3.5.0)
##  Rcpp                   0.12.16   2018-03-13 CRAN (R 3.5.0)
##  RCurl                  1.95-4.10 2018-01-04 CRAN (R 3.5.0)
##  rlang                  0.2.0     2018-02-20 CRAN (R 3.5.0)
##  rmarkdown            * 1.9       2018-03-01 CRAN (R 3.5.0)
##  rpart                  4.1-13    2018-02-23 CRAN (R 3.5.0)
##  rprojroot              1.3-2     2018-01-03 CRAN (R 3.5.0)
##  RSQLite                2.0       2017-06-19 CRAN (R 3.5.0)
##  rstudioapi             0.7       2017-09-07 CRAN (R 3.5.0)
##  S4Vectors            * 0.17.37   2018-03-20 Bioconductor  
##  scales                 0.5.0     2017-08-24 CRAN (R 3.5.0)
##  SingleCellExperiment * 1.1.3     2018-04-20 Bioconductor  
##  softImpute             1.4       2015-04-08 CRAN (R 3.5.0)
##  splatter             * 1.3.5     2018-04-27 Bioconductor  
##  splines                3.5.0     2017-12-27 local         
##  stabledist             0.7-1     2016-09-12 CRAN (R 3.5.0)
##  stats                * 3.5.0     2017-12-27 local         
##  stats4               * 3.5.0     2017-12-27 local         
##  stringi                1.1.7     2018-03-12 CRAN (R 3.5.0)
##  stringr                1.3.0     2018-02-19 CRAN (R 3.5.0)
##  SummarizedExperiment * 1.9.15    2018-03-20 Bioconductor  
##  survival               2.41-3    2017-04-04 CRAN (R 3.5.0)
##  testthat             * 2.0.0     2017-12-13 CRAN (R 3.5.0)
##  tibble                 1.4.2     2018-01-22 CRAN (R 3.5.0)
##  tools                  3.5.0     2017-12-27 local         
##  utils                * 3.5.0     2017-12-27 local         
##  withr                  2.1.2     2018-03-15 CRAN (R 3.5.0)
##  XML                    3.98-1.10 2018-02-19 CRAN (R 3.5.0)
##  xtable                 1.8-2     2016-02-05 CRAN (R 3.5.0)
##  XVector                0.19.9    2018-03-20 Bioconductor  
##  yaml                   2.1.18    2018-03-08 CRAN (R 3.5.0)
##  zinbwave             * 1.1.7     2018-04-27 Bioconductor  
##  zlibbioc               1.25.0    2017-12-27 Bioconductor
```
