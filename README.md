# rVictis
R interface to https://github.com/stuchly/vaevictis with naive importance sampling

## Installation
November 12, 2020 - reticulate version 1.16 or devel version - install_github("rstudio/reticulate") - is needed

``` r
devtools::install_github("stuchly/rVictis")
```
If you did not install vaevictis, restart session and load library - you will be prompted to install (type Yes); Vaevictis version 0.3.1 and higher is required.
## Example
``` r
library(rVictis)

res <- vaevictis(iris[,1:4],upsample=list(N=100,cluster=10)) ## 10 clusters by clara, 100 examples from each cluster

plot(res[[1]],col=as.factor(iris[,5])) ## plot reduced data

rediris <- res[[3]](as.matrix(iris[,1:4])) ##apply train reduction on "new" data, now we pass the data to python function - must be matrix 

plot(rediris,col=as.factor(iris[,5]))
```


