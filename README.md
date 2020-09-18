# Turbofuns

The goal of Turbofuns is to estimate polychoric correlations and their asymptotic covariance matrix (ACM).



This package can be installed directly from GitHub:

install.packages("devtools")

devtools::install_github("pdand/Turbofuns")



## Example

Examples using the data sets included in the packages:

data("BFI228")                   # Big-five inventory (N = 228)
#For ordinal data, estimating the polychoric correlation and its ACM
#with 5 cores and 1/(nc*nr) added to all cells

polyACM = PolychoricRM(BFI228,NCore=5, IAdjust=1, estimate.acm=TRUE)
```R
...
```
