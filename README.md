# venneR
There are several amazing packages to create a Venn diagram in R. 
However, they all have a different structure and there is a steep learning curve for majority of them. 
**venneR** has the motto of "Easy Venn Diagrams with R" and it creates Venn diagrams only using the a list object.

## Install
Please install `devtools` if you haven't yet.
```r
install.packages("devtools")
```
Then, you basically install the `venneR` package by using:
```r
library(devtools)
install_github("altintasali/venneR")
```

## Usage
```r
library(venneR)

# Generate artificial sets
x <- list(A=c(letters[1:5]),
          B=c(letters[3:10]),
          C=c(letters[10:15]))
universe <- length(letters)

# Plot 'all' Venn diagrams
res <- venneR(x, universe, stats = TRUE, pairwise = FALSE)
res$stat
res$plot

# Plot 'pairwise' Venn diagrams
res <- venneR(x, universe, stats = TRUE, pairwise = TRUE)
res
```
