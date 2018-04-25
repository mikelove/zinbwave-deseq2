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
##  9070   930
```

```r
zinb <- sim[keep,]
zinb$condition <- factor(zinb$Group)
# we need to reorganize the assays in the SumExp from splatter
nms <- c("counts", setdiff(assayNames(zinb), "counts"))
assays(zinb) <- assays(zinb)[nms]
# ~15 sec for ~900 genes and 100 cells
# epsilon setting as recommended by Koen and Fanny
system.time({
  zinb <- zinbwave(zinb, K=0, BPPARAM=SerialParam(), epsilon=1e12)
})
```

```
##    user  system elapsed 
##  16.783   1.873  18.865
```

### Estimate dispersion and DE using *DESeq2*

Below we use the default Wald test, although Van De Berge and
Perraudeau and others have shown the LRT may perform better for null
hypothesis testing.


```r
suppressPackageStartupMessages(library(DESeq2))
dds <- DESeqDataSet(zinb, design=~condition)
# ~15 sec
system.time({
  dds <- DESeq(dds, sfType="poscounts", useT=TRUE, minmu=1e-6, minRep=Inf)
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
##  14.191   0.750  15.289
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

<img src="zinbwave-deseq2_files/figure-html/unnamed-chunk-6-1.png" width="672" />

One would then assign the dispersion function to the original dataset,
re-estimate final dispersions, check `plotDispEsts`, and then either
re-run the Wald or LRT function (this chunk not evaluated):


```r
dispersionFunction(dds) <- dispersionFunction(dds2)
dds <- estimateDispersionsMAP(dds)
# for the Wald:
dds <- nbinomWaldTest(dds, useT=TRUE, minmu=1e-6)
# or for the LRT:
dds <- nbinomLRT(dds, minmu=1e-6)
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
##     FALSE   676   29
##     TRUE      6  219
```

```r
round(prop.table(tab, 2), 3)
```

```
##          sig
## DE.status FALSE  TRUE
##     FALSE 0.991 0.117
##     TRUE  0.009 0.883
```


```r
session_info()
```

```
## Session info -------------------------------------------------------------
```

```
##  setting  value                                             
##  version  R Under development (unstable) (2017-12-07 r73859)
##  system   x86_64, darwin15.6.0                              
##  ui       X11                                               
##  language (EN)                                              
##  collate  en_US.UTF-8                                       
##  tz       America/New_York                                  
##  date     2018-04-25
```

```
## Packages -----------------------------------------------------------------
```

```
##  package              * version  date      
##  acepack                1.4.1    2016-10-29
##  ADGofTest              0.3      2011-12-28
##  annotate               1.57.2   2017-12-07
##  AnnotationDbi          1.41.4   2018-01-10
##  backports              1.1.2    2017-12-13
##  base                 * 3.5.0    2017-12-07
##  base64enc              0.1-3    2015-07-28
##  Biobase              * 2.39.2   2018-02-07
##  BiocGenerics         * 0.25.3   2018-04-25
##  BiocInstaller        * 1.29.6   2018-04-10
##  BiocParallel         * 1.13.3   2018-04-25
##  bit                    1.1-12   2014-04-09
##  bit64                  0.9-7    2017-05-08
##  bitops                 1.0-6    2013-08-17
##  blob                   1.1.0    2017-06-17
##  checkmate              1.8.5    2017-10-24
##  cluster                2.0.6    2017-03-10
##  codetools              0.2-15   2016-10-05
##  colorspace             1.3-2    2016-12-14
##  compiler               3.5.0    2017-12-07
##  copula                 0.999-18 2017-09-01
##  data.table             1.10.4-3 2017-10-27
##  datasets             * 3.5.0    2017-12-07
##  DBI                    0.7      2017-06-18
##  DelayedArray         * 0.5.32   2018-04-25
##  DESeq2               * 1.19.51  2018-04-23
##  devtools             * 1.13.4   2017-11-09
##  digest                 0.6.15   2018-01-28
##  edgeR                  3.21.10  2018-04-25
##  evaluate               0.10.1   2017-06-24
##  foreach                1.4.4    2017-12-12
##  foreign                0.8-70   2017-11-28
##  Formula                1.2-2    2017-07-10
##  genefilter             1.61.1   2018-02-07
##  geneplotter            1.57.0   2017-12-07
##  GenomeInfoDb         * 1.15.5   2018-02-07
##  GenomeInfoDbData       1.1.0    2018-01-10
##  GenomicRanges        * 1.31.23  2018-04-25
##  ggplot2                2.2.1    2016-12-30
##  glmnet                 2.0-16   2018-04-02
##  graphics             * 3.5.0    2017-12-07
##  grDevices            * 3.5.0    2017-12-07
##  grid                   3.5.0    2017-12-07
##  gridExtra              2.3      2017-09-09
##  gsl                    1.9-10.3 2017-01-05
##  gtable                 0.2.0    2016-02-26
##  Hmisc                  4.0-3    2017-05-02
##  htmlTable              1.11.2   2018-01-20
##  htmltools              0.3.6    2017-04-28
##  htmlwidgets            1.0      2018-01-20
##  IRanges              * 2.13.28  2018-04-25
##  iterators              1.0.9    2017-12-12
##  knitr                  1.17     2017-08-10
##  lattice                0.20-35  2017-03-25
##  latticeExtra           0.6-28   2016-02-09
##  lazyeval               0.2.1    2017-10-29
##  limma                  3.35.15  2018-04-25
##  locfit                 1.5-9.1  2013-04-20
##  magrittr               1.5      2014-11-22
##  Matrix                 1.2-12   2017-11-16
##  matrixStats          * 0.53.1   2018-02-11
##  memoise                1.1.0    2017-04-21
##  methods              * 3.5.0    2017-12-07
##  munsell                0.4.3    2016-02-13
##  mvtnorm                1.0-7    2018-01-26
##  nnet                   7.3-12   2016-02-02
##  numDeriv               2016.8-1 2016-08-27
##  parallel             * 3.5.0    2017-12-07
##  pcaPP                  1.9-73   2018-01-14
##  pillar                 1.1.0    2018-01-14
##  plyr                   1.8.4    2016-06-08
##  pspline                1.0-18   2017-06-12
##  R6                     2.2.2    2017-06-17
##  RColorBrewer           1.1-2    2014-12-07
##  Rcpp                   0.12.15  2018-01-20
##  RCurl                  1.95-4.8 2016-03-01
##  rlang                  0.1.6    2017-12-21
##  rmarkdown            * 1.8      2017-11-17
##  rpart                  4.1-12   2018-01-10
##  rprojroot              1.3-2    2018-01-03
##  RSQLite                2.0      2017-06-19
##  rstudioapi             0.7      2017-09-07
##  S4Vectors            * 0.17.42  2018-04-25
##  scales                 0.5.0    2017-08-24
##  SingleCellExperiment * 1.1.3    2018-04-20
##  softImpute             1.4      2015-04-08
##  splatter             * 1.3.5    2018-04-25
##  splines                3.5.0    2017-12-07
##  stabledist             0.7-1    2016-09-12
##  stats                * 3.5.0    2017-12-07
##  stats4               * 3.5.0    2017-12-07
##  stringi                1.1.6    2017-11-17
##  stringr                1.2.0    2017-02-18
##  SummarizedExperiment * 1.9.17   2018-04-25
##  survival               2.41-3   2017-04-04
##  testthat             * 2.0.0    2017-12-13
##  tibble                 1.4.2    2018-01-22
##  tools                  3.5.0    2017-12-07
##  utils                * 3.5.0    2017-12-07
##  withr                  2.1.1    2017-12-19
##  XML                    3.98-1.9 2017-06-19
##  xtable                 1.8-2    2016-02-05
##  XVector                0.19.9   2018-04-25
##  yaml                   2.1.16   2017-12-12
##  zinbwave             * 1.1.7    2018-04-20
##  zlibbioc               1.25.0   2017-12-07
##  source                           
##  CRAN (R 3.5.0)                   
##  CRAN (R 3.5.0)                   
##  Bioconductor                     
##  Bioconductor                     
##  cran (@1.1.2)                    
##  local                            
##  CRAN (R 3.5.0)                   
##  Bioconductor                     
##  cran (@0.25.3)                   
##  Bioconductor                     
##  Bioconductor                     
##  CRAN (R 3.5.0)                   
##  CRAN (R 3.5.0)                   
##  CRAN (R 3.5.0)                   
##  CRAN (R 3.5.0)                   
##  CRAN (R 3.5.0)                   
##  CRAN (R 3.5.0)                   
##  CRAN (R 3.5.0)                   
##  CRAN (R 3.5.0)                   
##  local                            
##  CRAN (R 3.5.0)                   
##  CRAN (R 3.5.0)                   
##  local                            
##  CRAN (R 3.5.0)                   
##  cran (@0.5.32)                   
##  Bioconductor                     
##  CRAN (R 3.5.0)                   
##  cran (@0.6.15)                   
##  cran (@3.21.10)                  
##  CRAN (R 3.5.0)                   
##  CRAN (R 3.5.0)                   
##  CRAN (R 3.5.0)                   
##  CRAN (R 3.5.0)                   
##  Bioconductor                     
##  Bioconductor                     
##  Bioconductor                     
##  Bioconductor                     
##  cran (@1.31.23)                  
##  CRAN (R 3.5.0)                   
##  CRAN (R 3.5.0)                   
##  local                            
##  local                            
##  local                            
##  CRAN (R 3.5.0)                   
##  CRAN (R 3.5.0)                   
##  CRAN (R 3.5.0)                   
##  CRAN (R 3.5.0)                   
##  CRAN (R 3.5.0)                   
##  CRAN (R 3.5.0)                   
##  CRAN (R 3.5.0)                   
##  cran (@2.13.28)                  
##  CRAN (R 3.5.0)                   
##  CRAN (R 3.5.0)                   
##  CRAN (R 3.5.0)                   
##  CRAN (R 3.5.0)                   
##  CRAN (R 3.5.0)                   
##  cran (@3.35.15)                  
##  CRAN (R 3.5.0)                   
##  CRAN (R 3.5.0)                   
##  CRAN (R 3.5.0)                   
##  CRAN (R 3.5.0)                   
##  CRAN (R 3.5.0)                   
##  local                            
##  CRAN (R 3.5.0)                   
##  CRAN (R 3.5.0)                   
##  CRAN (R 3.5.0)                   
##  CRAN (R 3.5.0)                   
##  local                            
##  CRAN (R 3.5.0)                   
##  CRAN (R 3.5.0)                   
##  CRAN (R 3.5.0)                   
##  CRAN (R 3.5.0)                   
##  CRAN (R 3.5.0)                   
##  CRAN (R 3.5.0)                   
##  cran (@0.12.15)                  
##  CRAN (R 3.5.0)                   
##  CRAN (R 3.5.0)                   
##  CRAN (R 3.5.0)                   
##  CRAN (R 3.5.0)                   
##  cran (@1.3-2)                    
##  CRAN (R 3.5.0)                   
##  CRAN (R 3.5.0)                   
##  cran (@0.17.42)                  
##  CRAN (R 3.5.0)                   
##  Bioconductor                     
##  CRAN (R 3.5.0)                   
##  Github (Oshlack/splatter@d5ba92c)
##  local                            
##  CRAN (R 3.5.0)                   
##  local                            
##  local                            
##  CRAN (R 3.5.0)                   
##  CRAN (R 3.5.0)                   
##  cran (@1.9.17)                   
##  CRAN (R 3.5.0)                   
##  CRAN (R 3.5.0)                   
##  CRAN (R 3.5.0)                   
##  local                            
##  local                            
##  CRAN (R 3.5.0)                   
##  CRAN (R 3.5.0)                   
##  CRAN (R 3.5.0)                   
##  cran (@0.19.9)                   
##  CRAN (R 3.5.0)                   
##  Bioconductor                     
##  Bioconductor
```
