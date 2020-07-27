#'Plot Venn diagrams and check if the pairwise intersection (overlap) is statistically significant
#'
#'@param x List object with character vectors of Venn (set) elements. Each element of list is a set. List names will be used as set names.
#'@param universe Numeric defining the background for Fisher's Exact Test. See \code{\link{newGeneOverlap}} for details.
#'@param stats Logical for plotting the statistics (default TRUE). It uses Fisher's exact test with alternative hypothesis "greater". See \code{\link{fisher.test}} for details.
#'@param pairwise Logical for plotting Venn diagrams in pairwise fashion (default FALSE). If TRUE, pairwise Venn diagrams will be generated.
#'@param cols Vector of colors. The length should be equal to the number of sets. If NULL, automatics colors will be generated.
#'@param alpha Color Transparency parameter as integer between 0 (fully transparent) and 1 (fully opaque). Defaults to 0.5.
#'@param cex Scaling parameter for font sizes. Defaults to 1.5.
#'@param p.adjust Method for adjusting p-values for multiple testing. Default is "BH". See \code{\link{p.adjust}} for details.
#'@param ... Other parameters for see \code{\link{venn.diagram}} function.
#'@return A list object containing
#' \itemize{
#'     \item \code{plot} A grid object ready to be plotted. Can be directly plotted by calling the object with preferred graphics output device (e.g. \code{pdf})
#'     \item \code{stat} A data.table with Fisher's Exact Test statistics
#'     \itemize{
#'          \item \code{Group1} The name of the first Venn diagram
#'          \item \code{Group2} The name of the second Venn diagram
#'          \item \code{p} p-value
#'          \item \code{p.adj} Adjusted p-value
#'          \item \code{OR} Odds ratio
#'     }
#'}
#'@examples
#'# Generate artificial sets
#'x <- list(A=c(letters[1:5]),
#'          B=c(letters[3:10]),
#'          C=c(letters[10:15]))
#'universe <- length(letters)
#'
#'# Plot 'all' Venn diagrams
#'res <- venneR(x, universe, stats = TRUE, pairwise = FALSE)
#'res$stat
#'res$plot
#'
#'# Plot 'pairwise' Venn diagrams
#'res <- venneR(x, universe, stats = TRUE, pairwise = TRUE)
#'res
#'@export
venneR <- function(x,
                  universe,
                  stats = TRUE,
                  pairwise = TRUE,
                  cols = NULL,
                  alpha = 0.5,
                  cex = 1.5,
                  p.adjust = "BH",
                  ...
)
{
  ##------------------------------------------------------------------------
  ## Prepare input and global variables
  ##------------------------------------------------------------------------
  ##-----------------------------------------------
  ## Check input parameters
  ##-----------------------------------------------
  if(!is.list(x)){stop("Input parameter is not a list object.")}
  if(is.null(names(x))){
    warning("'x' is not a named list. Using 'Set' as the default name for Venn Sets")
    names(x) <- paste("Set", 1:length(x))
  }

  # turn-off logger
  flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger") #futile.logger

  ##-----------------------------------------------
  ## Plot parameters
  ##-----------------------------------------------
  cex <- cex
  fonts <- "sans"
  if(is.null(cols)){cols <- brewer.pal(length(x), "Set1")[1:length(x)] %>% rev}

  ##-----------------------------------------------
  ## split into pairwise sets
  ##-----------------------------------------------
  sets <- combn(length(x), 2)
  setcombs <- split(sets, rep(1:ncol(sets), each = nrow(sets)))

  ##------------------------------------------------------------------------
  ## Calculate Stats
  ##------------------------------------------------------------------------
  stat.out <- lapply(setcombs, function(a){
    go <- newGeneOverlap(listA = x[[a[1]]] %>% as.character %>% unique,
                         listB = x[[a[2]]] %>% as.character %>% unique,
                         genome.size = universe)
    go <- testGeneOverlap(go)
    go.out <- data.table(Group1 = names(x)[a[1]],
                         Group2 = names(x)[a[2]],
                         p = go@pval,
                         OR = go@odds.ratio)
    return(go.out)
  })

  statsMatrix <- rbindlist(stat.out)
  statsMatrix[, p.adj := p.adjust(p, method = p.adjust)]
  statsMatrix <- statsMatrix[, c(1,2,3,5,4)]

  statsGrob <- statsMatrix
  statsGrob[, p := scientific(p, 3)]
  statsGrob[, p.adj := scientific(p.adj, 3)]
  statsGrob[, OR := round(OR, 2)]
  statsGrob <- tableGrob(statsGrob, rows = NULL)

  ##------------------------------------------------------------------------
  ## Plot Venn Diagrams
  ##------------------------------------------------------------------------
  ##-----------------------------------------------
  ## Subtitles
  ##-----------------------------------------------
  if(stats){

    id <- (sets[1,] %in% setcombs[[1]][1]) & (sets[2,] %in% setcombs[[1]][2])
    id <- which(id)

    subTitle <- paste0("p = ", stat.out[[id]]$p %>% scientific(., 3),
                       ", OR = ", stat.out[[id]]$OR %>% round(., 2))
  }else{
    subTitle <- NULL
  }

  ##-----------------------------------------------
  ## Plot pair-wise combinations
  ##-----------------------------------------------
  if(pairwise){

    plots <- lapply(setcombs, function(a){

      # mainTitle <- paste(names(x)[a], collapse = " vs ")

      if(stats){

        id <- (sets[1,] %in% a[1]) & (sets[2,] %in% a[2])
        id <- which(id)

        subTitle <- paste0("p = ", stat.out[[id]]$p %>% scientific(., 3),
                           ", OR = ", stat.out[[id]]$OR %>% round(., 2))
      }else{
        subTitle <- NULL
      }

      p <- venn.diagram(x=list(A=x[[a[1]]] %>% as.character %>% unique,
                               B=x[[a[2]]] %>% as.character %>% unique),
                        category.names = names(x)[a],
                        alpha = 0.5,
                        fill = cols[a],
                        filename = NULL,
                        main = NULL,
                        sub = subTitle,
                        sub.pos = c(0.5, 0.1),
                        cex = cex,
                        main.cex = cex,
                        cat.cex = cex,
                        fontfamily = fonts,
                        cat.fontfamily = fonts,
                        main.fontfamily = fonts,
                        sub.fontfamily = fonts,
                        cat.pos = 0 + c(-35,35), #TODO: For some reason, 0 degrees either start at 12 or 6 o'clock. Needs a fix!
                        ...
      )

      return(p)
    })

    pOut <- do.call(plot_grid, plots)

    ##-----------------------------------------------
    ## Plot NON pair-wise combinations
    ##-----------------------------------------------
  }else{
    p <- venn.diagram(x=x,
                      category.names = names(x),
                      alpha = 0.5,
                      fill = cols,
                      filename = NULL,
                      main = NULL,
                      sub = NULL,
                      cex = cex,
                      main.cex = cex,
                      cat.cex = cex,
                      fontfamily = fonts,
                      cat.fontfamily = fonts,
                      main.fontfamily = fonts,
                      sub.fontfamily = fonts,
                      #cat.pos = 0 + c(-35,35),
                      ...
    )

    if(stats){
      pOut <- plot_grid(p, statsGrob, ncol = 1, rel_heights = c(3,1))
    }else{
      pOut <- plot_grid(p)
    }
  }

  ##------------------------------------------------------------------------
  ## Generate output
  ##------------------------------------------------------------------------
  if(stats){
    # TODO: Weird bug: the numbers do not round if I do not recalculate the 'statsMatrix'
    statsMatrix <- rbindlist(stat.out)
    statsMatrix[, p.adj := p.adjust(p, method = "BH")]
    statsMatrix <- statsMatrix[, c(1,2,3,5,4)]

    out <- list(plot = pOut, stat = statsMatrix)
  }else{
    out <- list(plot = pOut)
  }

  return(out)
}
