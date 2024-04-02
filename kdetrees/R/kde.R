### Copyright (C) 2014 -- Grady Weyenberg ###

##' Analyze a set of phylogenetic trees and attempt to identify trees
##' which are significantly discordant with other trees in the sample
##' (outlier trees).
##' 
##' If bw is a single number, it will be used as a single constant
##' bandwidth. It can also be a vector, in which case it will be used
##' as variable bandwidths for each tree, repectively. Finally, if it
##' is a list (default), the list will be passed as arguments to the bw.nn
##' adaptive bandwith function.
##'
##' ... Is passed to either \code{distory::dist.multiPhylo} or
##' \code{dist.diss}, as appropriate. See the help for these functions
##' for more details.
##'
##' @title Identify discordant trees in a sample
##' @param trees multiPhylo object
##' @param k IQR multiplier for outlier detection
##' @param distance Select "geodesic" or "dissimilarity" distance
##' calculation method
##' @param outgroup if a character, reroot all trees with this species
##' as outgroup. The geodesic distance method requires rooted trees.
##' @param topo.only set all branch lengths to 1 before analyzing?
##' @param bw see Details
##' @param greedy greedy outlier detection?
##' @param ... additional arguments for distance calculation function, see details
##' @return a kdetrees object; list(density,outliers)
##' @author Grady Weyenberg
##' @export
##' @examples
##' kdeobj <- kdetrees(apicomplexa)
##' print(kdeobj)
##' kdeobj$outliers
##'
##' kdetrees(apicomplexa, k=2.0, distance="dissimilarity",topo.only=FALSE)
kdetrees <- function(trees,k=1.5,distance=c("geodesic","dissimilarity"), outgroup=NULL,
                     topo.only=FALSE,bw=list(),greedy=FALSE,...) {
  distance <- match.arg(distance)
  if (!inherits(trees,"multiPhylo") && all(sapply(trees,inherits,"phylo"))) class(trees) <- "multiPhylo"

  if (topo.only) {
    trees <- lapply(trees,compute.brlen, method = 1)
    class(trees) <- "multiPhylo"
  }
  
  if (is.character(outgroup)) {
    trees <- lapply(trees,root,outgroup,resolve.root=TRUE)
    trees <- lapply(trees,"[<-","node.label", NULL)
    class(trees) <- "multiPhylo"
  }

  dm <- switch(distance,
               geodesic = as.matrix(dist.multiPhylo(trees,...)),
               dissimilarity = as.matrix(dist.diss(trees,...)))
  dimnames(dm) <- list(names(trees),names(trees))

  cutoff <- function(x, c = 1.5){
    qs <- quantile(x, c(0.25,0.75))
    unname(diff(qs) * -c + qs[1])
  }
  
  if(is.list(bw)) bw <- do.call(bw.nn,c(list(dm),bw))
  km <- normkern(dm,bw)
  x <- estimate(km)
  c <- cutoff(x, k)
  i <- which( x < c )
  if (greedy) {
    while(TRUE){
      i <- which( x < c )
      if (length(i) < 1) break
      x <- estimate(km,i)
      c2 <- cutoff(x[-i], k)
      ## if(is.na(c2)) browser()
      if(c2 > c) c <- c2 else break
    }
  }
  structure(list(density=x, i=i, outliers=trees[i]), class="kdetrees",
            call=match.call(), c=c)
}


##' Performs a complete kdetrees analysis, starting with reading trees
##' from a newick file on disk, and writing result files to the
##' working directory. Names and location of output files may be
##' controlled by optional arguments.
##'
##' @title Complete kdetrees analysis convenience function
##' @param infile newick file with trees
##' @param ... additional parameters for kdetrees
##' @param treeoutfile write outlier trees in newick format to this file
##' @param csvfile write density results to this file
##' @param plotfile print scatterplot of results to this file
##' @param histfile print histogram of density estimates to this file
##' @return result of kdetrees call
##' @author Grady Weyenberg
##' @export
kdetrees.complete <- function(infile,...,treeoutfile="outliers.tre",
                              csvfile="results.csv",plotfile="plot.png",
                              histfile="hist.png"){
  trees <- read.tree(infile)
  if (is.null(names(trees))) names(trees) <- paste("tree",seq_along(trees),sep="")
  if (!inherits(trees,"multiPhylo")) stop("Could not read tree file")
  
  res <- kdetrees(trees,...)
  if (is.character(plotfile)) ggsave(plotfile,plot(res)) #scatterplot
  if (is.character(histfile)) ggsave(histfile,hist(res)) #histogram
  if (is.character(csvfile)) write.csv(as.data.frame(res),csvfile) #csv
  if (is.character(treeoutfile) && length(res$outliers) > 0)
    write.tree(res$outliers, treeoutfile, tree.names=TRUE,digits=5) #out-trees
  
  res
}

##' estimate densities from kernel matrix
##'
##' @param x matrix of kernel contributions
##' @param i vector of columns to exclude from calculation
##' @return vector of density estimates for each tree
##' @author Grady Weyenberg
estimate <- function(x,i=integer()){
  if(length(i) > 0)
    rowSums(x[,-i]) - (!(1:nrow(x) %in% i)) * diag(x)
  else
    rowSums(x) - diag(x)
}



