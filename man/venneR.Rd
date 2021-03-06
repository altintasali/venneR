% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/venneR.R, R/venneR_package.R
\docType{package}
\name{venneR}
\alias{venneR}
\title{Plot Venn diagrams and check if the pairwise intersection (overlap) is statistically significant}
\usage{
venneR(
  x,
  universe,
  stats = TRUE,
  pairwise = TRUE,
  cols = NULL,
  alpha = 0.5,
  cex = 1.5,
  p.adjust = "BH",
  ...
)
}
\arguments{
\item{x}{List object with character vectors of Venn (set) elements. Each element of list is a set. List names will be used as set names.}

\item{universe}{Numeric defining the background for Fisher's Exact Test. See \code{\link{newGeneOverlap}} for details.}

\item{stats}{Logical for plotting the statistics (default TRUE). It uses Fisher's exact test with alternative hypothesis "greater". See \code{\link{fisher.test}} for details.}

\item{pairwise}{Logical for plotting Venn diagrams in pairwise fashion (default FALSE). If TRUE, pairwise Venn diagrams will be generated.}

\item{cols}{Vector of colors. The length should be equal to the number of sets. If NULL, automatics colors will be generated.}

\item{alpha}{Color Transparency parameter as integer between 0 (fully transparent) and 1 (fully opaque). Defaults to 0.5.}

\item{cex}{Scaling parameter for font sizes. Defaults to 1.5.}

\item{p.adjust}{Method for adjusting p-values for multiple testing. Default is "BH". See \code{\link{p.adjust}} for details.}

\item{...}{Other parameters for see \code{\link{venn.diagram}} function.}
}
\value{
A list object containing
\itemize{
    \item \code{plot} A grid object ready to be plotted. Can be directly plotted by calling the object with preferred graphics output device (e.g. \code{pdf})
    \item \code{stat} A data.table with Fisher's Exact Test statistics
    \itemize{
         \item \code{Group1} The name of the first Venn diagram
         \item \code{Group2} The name of the second Venn diagram
         \item \code{p} p-value
         \item \code{p.adj} Adjusted p-value
         \item \code{OR} Odds ratio
    }
}
}
\description{
Generates Venn diagrams and checks whether pairwise set intersections are statistically significant
}
\examples{
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
}
\author{
Ali Altintas \email{altintas.ali@gmail.com}
}
