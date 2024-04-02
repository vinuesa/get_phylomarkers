####
#### This is a modified version of the Phylo-MCOA program for use in
#### the simulation routine. The main modification is the removal of
#### numerous progress bars, diagnostic messages, and other extraneous
#### code.
####

######################################################### 
##                   Damien M. de Vienne               ##
##                   Phylo-MCOA functions              ##
##                    13th of July 2012                ##
##                      Version 1.1                    ##
##                                                     ##
##         Go to http://phylomcoa.cgenomics.org/       ##
##        for tutorial, help for installing, and       ##
##           description of all the functions          ##
##                                                     ##
## If you use this tool for a publication, please cite ##
##   de Vienne, DM., Ollier, S and Aguileta, G. 2012.  ##
##           Molecular Biology and Evolution           ## 
#########################################################

require(ape)
require(ade4)
require(utils)
require(lattice)

trees2matrices <- function(trees,distance=c("nodal","patristic"),bvalue=0){
  distance <- match.arg(distance)
  collapse.nodes <- function(y){
    if(is.null(y$node.label))
      return(di2multi(y,bvalue))
    ##else the node labels have bootstrap(?) values
    idx <- which(as.numeric(y$node.label) < bvalue) + Ntip(y)
    y$edge.length[y$edge[,1] == idx] <- 1e-10
    return(di2multi(y,1e-9))
  }
  if(bvalue != 0)
    trees <- lapply(trees,collapse.nodes)
  if(distance == "nodal")
    trees <- lapply(trees,compute.brlen,1.0)
  lapply(trees,cophenetic)
}

##' Prepares list of tree matrices for analysis. Inserts missing
##' rows/cols into each matrix. Reorders all rows and columns to match
##' each other. Fills any missing numbers with elemetwise average
##' across matrices.
##' @param matrices a list of tree matrices
##' @return a list of processed matrices
gestion.mat<-function(matrices) {
  listsp <- Reduce(union,lapply(matrices,colnames))
### if all trees have same # of tips then only  re-order the columns
  if ( all(sapply(matrices,ncol) == length(listsp)) )
    return(lapply(matrices,"[", rownames(matrices[[1]]), colnames(matrices[[1]])))
### else we insert missing columns
  pad.matrix <- function(y){
    new.names <- setdiff(listsp, colnames(y))
    new.dim <- nrow(y) + length(new.names)
    res <- matrix(NA, new.dim, new.dim, dimnames = rep(list(c(rownames(y), new.names)), 2))
    res[1:nrow(y), 1:ncol(y)] <- y
  }
  matrices <- lapply(matrices, pad.matrix)
### rearrange columns to match
  matrices <- lapply(matrices,"[", rownames(matrices[[1]]), colnames(matrices[[1]]))
### find elementwise averages
  meanmat <- apply(simplify2array(matrices), 1:2, mean, na.rm=TRUE)
  meanmat[is.nan(meanmat)] <- mean(meanmat,na.rm=TRUE)
### and fill in NAs with mean values
  lapply(matrices,function(y) y[is.na(y)] <- meanmat[is.na(y)])
}

mat2mcoa<-function(matrices, wtts=NULL, scannf=TRUE, nf="auto") {
  dist.all <- lapply(matrices, function(x) suppressWarnings(cailliez(sqrt(as.dist(x)))))
  dist.all.ktab <- kdist2ktab(kdist(dist.all))

  if (nf=="auto") {
    scannf <- FALSE
    nf <- min(20, ncol(matrices[[1]]) - 1)
  }
  
  if (is.null(wtts))
    return(mcoa(dist.all.ktab, scannf=scannf, nf=nf))
  else
    dist.all.ktab$tabw <- abs(as.numeric(wtts))
  return(mcoa(dist.all.ktab,option="internal", scannf=scannf,nf=nf))
}

mcoa2WRmat <- function(mcoa) {
  result <- with(mcoa, sqrt(rowSums((Tl1 - SynVar[TL[ , 2], ])^2)))
  dim(result) <- with(mcoa, c(nrow(SynVar), nrow(lambda)))
  rownames(result) <- rownames(mcoa$SynVar)
  colnames(result) <- rownames(mcoa$lambda)
  result
}

pMCOA <- function(trees, distance="nodal", bvalue=0, wtts=NULL, scannf=TRUE, nf="auto"){
  if (is.list(trees) && class(trees[[1]]) != "phylo")
    stop("The trees should be in the \"phylo\" format!")
  if (class(trees)=="character")
    trees <- read.tree(trees)
  if (is.null(names(trees)))
    names(trees) <- paste("gene", (1:length(trees)))
  if (length(unique(names(trees))) < length(trees))
    warning("Gene tree names not unique")
  b<-trees2matrices(trees, distance=distance, bvalue=bvalue)
  c<-gestion.mat(b)
  d<-mat2mcoa(c, wtts, scannf, nf)
  mcoa2WRmat(d)
}






pMCOA.complete <- function(trees, distance="nodal", bvalue=0, wtts=NULL, scannf=TRUE, nf="auto", k=1.5, thres=0.2) {
  res1 <- pMCOA(trees, distance, bvalue, wtts, scannf, nf)
  total.outl <- detect.complete.outliers(res1, k, thres)
  if (with(total.outl, sum(TFgn, TFsp)) > 0) {
    trees <- rm.gene.and.species(total.outl,trees)
    res2 <- pMCOA(trees, distance, bvalue, wtts, scannf, nf)
  }
  else
    res2 <- res1
  res3 <- detect.cell.outliers(res2, k)
  res <- list(trees = trees, step1 = res1, outcompl = total.outl, step2 = res2, outcell = res3)
  invisible(res)
}

rm.gene.and.species <- function(obj, trees){
  if (length(obj$outsp) > 0)
    trees <- lapply(trees, drop.tip, obj$outsp)
  structure(trees[!obj$TFgn], class="multiPhylo")
}

normalize<-function(mat, what=c("none","species","genes")) {
  what <- match.arg(what)
  fun <- function(x) x / mean(x)
  switch(what, species = apply(mat, 2, fun), genes = t(apply(mat, 1, fun)), mat)
}

outl.sub <- function(x,k) return(x > quantile(x,0.75) + k * IQR(x))


detect.complete.outliers <- function(mat2WR, k=1.5, thres=0.5) {

  tabgn <- normalize(mat2WR, "genes")
  tabgn.TF <- t(apply(tabgn,1,outl.sub, k))

  tabsp <- normalize(mat2WR, "species")
  tabsp.TF <- apply(tabsp,2, outl.sub, k)

  score.genes <- colMeans(tabgn.TF)
  tf.gn <- score.genes > thres

  score.species <- rowMeans(tabsp.TF)
  tf.sp <- score.species > thres
  
  invisible(list(scoregn = score.genes, scoresp = score.species,
                 TFgn = tf.gn, TFsp = tf.sp,
                 outgn = names(tf.gn[tf.gn]), outsp = names(tf.sp[tf.sp])))
}



detect.island <- function(arr) {
  spi <- 1:length(arr)
  names(spi) <- names(arr)
  
  true.names <- names(arr)[arr == TRUE]
  
  if (length(true.names) == 1)
    return(list(true.names))

  if (length(true.names) > 1) {
    true.i <- spi[true.names]
    res <- dist(true.i)
    table.i <- cbind(t(combn(attributes(res)$Labels,2)), array(res))
    in.island <- NULL
    if (length(table.i[table.i[,3] == "1", 3]) == 0) {
      in.island <- "nopair"
      list.i <- NULL
    }
    if (length(table.i[table.i[,3]=="1",3])==1) {
      in.island <- table.i[table.i[,3]=="1",c(1,2)]
      list.i <- list(in.island)
    }
    if (is.null(in.island)) {
      table.small <- table.i[table.i[,3]=="1",c(1,2)]        
      list.i <- list()
      for (i in 1:nrow(table.small))
        list.i[[i]] <- table.small[i,]
      for (i in 1:(length(list.i)-1)) {
        for (j in (i+1):length(list.i)) {
          if (length(intersect(list.i[[i]], list.i[[j]]))>0) {
            list.i[[i]] <- c(list.i[[i]], list.i[[j]])
            list.i[[j]] <- "out"
            list.i[[i]] <- unique(list.i[[i]])
          }
        }
      }
      list.i2 <- list()
      w <- 0
      for (i in 1:length(list.i)) {
        if ((length(list.i[[i]])>1)&&(list.i[[i]][1]!="out")) {
          w <- w+1
          list.i2[[w]] <- list.i[[i]]
          in.island <- c(in.island, list.i[[i]])
        }
      }
      list.i <- list.i2
    }
    out.island <- as.list(setdiff(true.names,in.island))
    return(c(list.i,out.island))
  }
  return(NULL)
}

detect.cell.outliers <- function(mat2WR, k=3) {
  MATspgn <- normalize(mat2WR, "genes") * normalize(mat2WR, "species")
  testspgn1 <- apply(MATspgn,2,outl.sub, k)
  testspgn2 <- t(apply(MATspgn,1,outl.sub, k))
  testspgn <- testspgn1 * testspgn2

  if (sum(testspgn) > 0) {
    out.list <- apply(testspgn,2, detect.island)
    genes <- colnames(testspgn)
    res <- c(NA,NA)
    for (i in 1:length(out.list)) {
      if (!is.null(out.list[[i]])) {
        for (j in 1:length(out.list[[i]])) {
          if (length(out.list[[i]][[j]])==1) res <- rbind(res, c(out.list[[i]][[j]],genes[i]))
          if (length(out.list[[i]][[j]])>1) {
            vals <- MATspgn[out.list[[i]][[j]], genes[i]]
            res <- rbind(res, c(names(vals)[vals==max(vals)], genes[i]))
          }
        }
      }
    }
    colnames(res) <- c("Species", "Genes")
    ##we construct the MATfinal
    MATfinal <- testspgn
    MATfinal[,] <- 0
    for (w in 2:nrow(res)) MATfinal[res[w,1], res[w,2]] <- 1

    return(list(mat2WR = mat2WR, matspgn = MATspgn, matfinal = MATfinal, testFALSE = testspgn, outcell = res[2:nrow(res), ]))
  }
  else return(NULL)
}
