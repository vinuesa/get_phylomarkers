### Copyright (C) 2014 -- Grady Weyenberg ###
##' Summarize a kdetrees object in human-readable form.
##'
##' Pretty-prints the results of a kdetrees ananlysis to console.
##'
##' @param x object to be printed
##' @param ... unused, required for generic compatability
##' @return invisible(x)
##' @author Grady Weyenberg
##' @method print kdetrees
##' @export
print.kdetrees <- function(x,...){
  outliers <- if(is.null(names(x$outliers))) x$i else names(x$outliers)
  cat("Call: ")
  print(attr(x,"call"))
  cat("Density estimates:\n")
  print(summary(x$density))
  cat("Cutoff: ",attr(x,"c"),"\n")
  cat("\nOutliers detected:\n")
  print(outliers,quote=FALSE)
  invisible(x)
}

##' Convert kdetrees object to data.frame
##'
##' Converts a kdetrees object to a data.frame suitable for saving as
##' output. It contains the density estimates for each tree, a Boolean
##' value indicating if the tree was selected as an outlier, and
##' optionally the newick string corresponding to the tree.
##' 
##' @param x kdetrees object to be converted
##' @param row.names ignored
##' @param optional ignored
##' @param trees If given the original list of trees, will convert to
##' newick and add a column to the output
##' @param ... unused
##' @return a data.frame
##' @author Grady Weyenberg
##' @method as.data.frame kdetrees
##' @export
##' @examples
##' result <- kdetrees(apicomplexa)
##' as.data.frame(result)
as.data.frame.kdetrees <- function(x, row.names, optional, trees=NULL, ...){
    out <- data.frame(density = x$density,
                      outlier=names(x$density) %in% names(x$outliers),
                      stringsAsFactors=FALSE)
    if (inherits(trees,"multiPhylo")) out$newick <- write.tree(trees,digits=4L)
    out
}
  
##' Plot the unnormalized density estimates for each tree.
##'
##' @param x kdetrees object to be plotted
##' @param ... additional arguments passed to ggplot
##' @return a ggplot object
##' @author Grady Weyenberg
##' @export
##' @method plot kdetrees
##' @examples
##' result <- kdetrees(apicomplexa)
##' plot(result)
plot.kdetrees <- function(x,...){
  df <- with(x,data.frame(density=unname(density),
                          index=seq_along(density),
                          outlier=seq_along(density) %in% i))
  ylab <- "Non-normalized Density"
  xlab <- "Tree Index"
  main <- paste(length(x$outliers),"Outliers Removed")
  
  ggplot(df,aes(index,density,color=outlier),...) + geom_point() +
    labs(title=main, x=xlab,y=ylab) + theme(legend.position="top")
}


##' Create a histogram of tree density estimates
##' 
##' @param x kdetrees object to plot
##' @param ... additional arguments passed to ggplot
##' @return a ggplot object
##' @author Grady Weyenberg
##' @export
##' @method hist kdetrees
##' @examples
##' result <- kdetrees(apicomplexa)
##' hist(result)
hist.kdetrees <- function(x,...){
  df <- with(x,data.frame(density=unname(density),
                          index=seq_along(density),
                          outlier=seq_along(density) %in% i))
  bw <- with(x,diff(range(density))/nclass.FD(density))
  main <- paste("Histogram of Estimates:",length(x$outliers),"Outliers Removed")
  xlab <- "Non-normalized Density"
  ylab <- "Count"

  ggplot(df,aes(density,fill=outlier),...) + geom_histogram(binwidth=bw) +
    labs(title=main,x=xlab,y=ylab) + theme(legend.position="top")
}


## Suppress the NOTES from R CMD check about undefined variables (ggplot calls)
globalVariables(c("outlier","index"))

