### Copyright (C) 2014 -- Grady Weyenberg ###

##' nearest-neighbor adaptive bandwidth selection
##' 
##' For each row in pairwise distance matrix find the distance to the
##' closest prop fraction of trees. 
##' @param x pairwise distance matrix
##' @param prop fraction of data to define the local neighborhood
##' @param tol tolerance for zero-bandwidth check
##' @return a vector of bandwidths for each tree (row) in x
##' @author Grady Weyenberg
##' @export
##' @examples
##' dm <- as.matrix(dist.diss(apicomplexa[1:20]))
##' bw.nn(dm)
bw.nn <- function(x,prop=0.2,tol=1e-6){
  out <- apply(x,1,function(y) quantile(y,prop))
  is.zero <- out < tol
  if(sum(is.zero)>0) out[is.zero] <- apply(x[is.zero,],1,function(y) min(y[y>tol]))
  out
}

