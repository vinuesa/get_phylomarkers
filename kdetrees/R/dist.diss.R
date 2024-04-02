### Copyright (C) 2014 -- Grady Weyenberg ###

##' dissimilarity map tree vectorization
##'
##' Dissimilarity maps convert trees to vectors using tip-to-tip path
##' lengths. Branch length information may be optionally discarded
##' (the default), resulting in vectors based solely on tree
##' topology.
##' 
##' @param x an ape::multiPhylo object.
##' @param ... additional options for \code{ape::cophenetic.phylo}
##' @return a row matrix of tree vectors
##' @author Grady Weyenberg
##' @method as.matrix multiPhylo
##' @export
##' @examples
##' as.matrix(apicomplexa[1:5])
as.matrix.multiPhylo <- function(x,...){
  tip.labels <- unique(unlist(lapply(x,"[[","tip.label")))

  ##Unroot the trees (is there any reason not to do this?)
  ##if(unroot){ x <- lapply(x,unroot) }
  ##Set branch lengths to 1
  #if(!use.blen){ x <- lapply(x, compute.brlen, method = 1) }
  ##find tip-tip distances for a tree, if tips are missing, add NA
  tip2tip <- function(y){
    missing.tips <- setdiff(tip.labels,y$tip.label)
    o <- cophenetic(y,...)
    o <- rbind(o,matrix(nrow=length(missing.tips),ncol=ncol(o),dimnames=list(missing.tips)))
    o <- cbind(o,matrix(nrow=nrow(o),ncol=length(missing.tips),dimnames=list(rownames(o),missing.tips)))
    o[tip.labels,tip.labels][upper.tri(o)]
  }
  cnames <- outer(tip.labels, tip.labels, paste, sep="-")
  ##convert multiPhylo to a row matrix of tip distances
  out <- t(sapply(x,tip2tip))
  colnames(out) <- cnames[upper.tri(cnames)]
  if(is.null(rownames(out))) rownames(out) <- paste("tree", 1:nrow(out), sep="")
  out
}

##' Compute pairwise tree distances
##'
##' 
##'
##' @param x either a row matrix of tree vectors, or a multiPhylo object
##' @param ... additional arguments passed to as.matrix.multiPhylo
##' @param method option passed to dist
##' @param p option passed to dist
##' @seealso dist
##' @return a dist object with tree-to-tree distances
##' @author Grady Weyenberg
##' @export
##' @examples
##' dist.diss(apicomplexa[1:5])
dist.diss <- function(x,...,method="euclidean",p=2){
  ##pairwise tree distances
  if(inherits(x,"multiPhylo")) d <- as.matrix.multiPhylo(x,...)
  ##missing tip imputation
  if (any(is.na(d))) {
    cm <- apply(d,2,median,na.rm=TRUE)
    if (any(is.na(cm)))
      stop("There are some species which never appear in the same tree: ",
           paste(names(cm)[is.na(cm)], collapse=", "))
    for (j in 1:ncol(d)) { d[is.na(d[,j]),j] <- cm[j] }
    warning("Tip labels were not the same for all trees, missing values have been imputed.")
  }
  dist(d,method=method,p=p)
} 
