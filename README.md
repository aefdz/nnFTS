[![License](https://img.shields.io/badge/license-GPL%20v3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

Depth-based prediction bands for functional data
================================================

Enveloping without prediction or forecasting
--------------------------------------------

``` r
focal <- '1'
dist <- 'l2' # dist<-'supremum'
plotting <- TRUE 

resultsBand <- banddpeeling(data, focal, dist, plotting)
```

![](README_files/figure-markdown_github/unnamed-chunk-1-1.png)

    ## Press <Enter> to continue...

![](README_files/figure-markdown_github/unnamed-chunk-1-2.png)

    ## Press <Enter> to continue...

![](README_files/figure-markdown_github/unnamed-chunk-1-3.png)

    ## Press <Enter> to continue...

![](README_files/figure-markdown_github/unnamed-chunk-1-4.png)

    ## Press <Enter> to continue...

``` r
resultsBand #Envelope 
```

    ## $subsample
    ##  [1] "1"  "43" "29" "84" "22" "88" "28" "52" "63" "15" "41" "86" "97" "4"  "46"
    ## [16] "79" "14" "10"

Curve Extension
---------------

``` r
cut <- 25 # number of points observed of the partially observed function
kcurves <- 10 # number of curves of the envelope involved in the band

results <- extension(data, focal, cut, dist)

pl <- plotBand(data, cut, results$Jordered, kcurves, focal)
```

![](README_files/figure-markdown_github/unnamed-chunk-2-1.png)

To explore different values of cut and kcurves (only running in Rstudio)
------------------------------------------------------------------------

``` r
manipulate(
  {
  plotBand(data, cut, results$Jordered, kcurves, focal)
  },
  kcurves = slider(min = 1, max = length(results$Jordered), step = 1, ticks = TRUE),
  cut = slider(1, 99, initial = 50, step = 1)
)
```
