### R code from vignette source 'kdetrees.Rnw'

###################################################
### code chunk number 1: preamble
###################################################
ops <- options(width=60)
library(kdetrees)


###################################################
### code chunk number 2: kdetrees.complete (eval = FALSE)
###################################################
## kdetrees.complete("trees.tre")


###################################################
### code chunk number 3: read.tree (eval = FALSE)
###################################################
## library(ape)
## apicomplexa <- read.tree("apicomplexa.tre")


###################################################
### code chunk number 4: kdetrees
###################################################
result <- kdetrees(apicomplexa)
result


###################################################
### code chunk number 5: kdetrees.diss (eval = FALSE)
###################################################
## kdetrees(apicomplexa, k=1.25, distance="diss", topo.only=TRUE)


###################################################
### code chunk number 6: plot
###################################################
plot(result)


###################################################
### code chunk number 7: hist
###################################################
hist(result)


###################################################
### code chunk number 8: outlierplot
###################################################
plot(result$outliers[[1]])


###################################################
### code chunk number 9: write.tree (eval = FALSE)
###################################################
## write.tree(result$outliers,file="outliers.tre")
## result.df <- as.data.frame(result)
## write.csv(result.df, file="scores.csv")


###################################################
### code chunk number 10: adv1 (eval = FALSE)
###################################################
## kdetrees(apicomplexa, bw=list(prop=0.5))


###################################################
### code chunk number 11: adv2 (eval = FALSE)
###################################################
## kdetrees(apicomplexa, bw=6)


###################################################
### code chunk number 12: postamble
###################################################
options(ops)


