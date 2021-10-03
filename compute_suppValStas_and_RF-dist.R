#!/usr/bin/env Rscript

# compute_suppValStas_and_RF-dist_of_geneTrees2concat_tree_automatically.R

VERSION <- 'Version: 0.8 - May 16, 2017' # improved graphical summary of SH-supp vals and RF-distances 2 concat tree
AUTHOR <- "Authors: Pablo Vinuesa [CCG-UNAM], Bruno Contreras Moreira [EEAD-CSIC]; "
REPOS <- "https://cloud.r-project.org" 

# find script path 
cmd.args <- commandArgs()
m <- regexpr("(?<=^--file=).+", cmd.args, perl=TRUE)
script.dir <- dirname(regmatches(cmd.args, m))

LOCAL_LIB = paste(script.dir,"/lib/R",sep = "")
.libPaths( c( .libPaths(), LOCAL_LIB) )

# load/download required packages
# see require() vs. library() discussions in:
# http://yihui.name/en/2014/07/library-vs-require/
# http://stackoverflow.com/questions/5595512/what-is-the-difference-between-require-and-library

for (package in c('ape', 'gplots',  'seqinr')) {
  if (!require(package, character.only=T, quietly=T)) {
    sprintf("unable to locate %s, will install it for you from %s", package, REPOS)
    install.packages(package, repos=REPOS, lib=LOCAL_LIB)
    # library(package, character.only=T)
  }
  #else{ library(package, character.only=T, quietly=T) ) } # don't load all packages to avoid namespace pollution/collision
}

# library(package, character.only=T) # don't load all packages to avoid namespace pollution/collision
# load the required seqinr apt and gplots::functions individually
for ( pkg in  c("plyr", "stringr", "ggplot2") ){
  library( pkg, character.only=T, quietly=T) 
}

#-----------------------------
#>>> FUNCTION DEFINITIONS <<<#
#-----------------------------
print_help <- function(){
  cat("", "AIM: Computes and plots tree support value stats and gene2concatTree RF-distances",
        VERSION, AUTHOR,
        "",
        ">>> Usage: ~/R_code/scripts/compute_phyloStats_and_geneTrees2concatTree_RF-dist_autom.R <full_path2_ph_files> <Runmode> <aln_ext> <tree_ext> <genoFlag[0|1]> [basegraphs|ggplots]",
        "",
        ">>> Notes:",
        "# 1. ARGUMENTS:",
        "runmodes: 1 ==> get only tree support stats; 2 ==> compute also RF genetree2concat-tree distances",
        "genoFlag: 0 ==> a standard MLSA/MLST set, with < 20 loci/trees; 1 ==> a genomic dataset with >> 20 loci/trees",
        "",
        "# ASSUMPTIONS ON DATA:",
        " The sequences should be collapsed to haplotypes",
        " but beware that the R code below complains when the tree labels",
        " contain '#' symbols, as introduced by collapse2haplotypes.pl",
        " so use a line like the following to change # by - in the collapsed fastas (or on the trees)",
        " for file in *FAS; do collapse2haplotypes.pl $file | awk '{print $1, $2}' | perl -pe 'if(/^>/){ s/\\h/_/ }' | sed 's/#/-/' > ${file}UNIQ; done",
        " or: for file in *ph; do sed 's/#/-/ $file > ${file}ed; done",
        " IMPORTANT: remove trees collapsed with very few leaves: ls -l *ph | awk '{if ($5 < 250) print $9}' | tee -a filelog | xargs rm",
        "",
        "TODO:",
        "1. add progress messages to STDOUT ...",
        "2. Need to check trees for no. of leaves; skip those with lt 4 leaves from within the script",
        "3. Need to check why ggplots are not working when invoked as a script, but work when executing lines from rstudio",
        "4. Integrate code from run_kdetrees",
        "",
      sep ="\n")
}
#-----------------------------------------------------------------------

# see http://www.inside-r.org/r-doc/base/file.copy
# for details on file manipulation from R
checkFileCreated <- function(F){
  if( file.exists(F) ){
      sprintf("File %s was created", F)
  }else{
    sprintf("WARNING: File %s could not be written to disk!", F)
  }
} 
#-----------------------------------------------------------------------


#-------------
#>>> MAIN <<<#
#-------------

# Notes on passing argument to the script
# add something like getops to pass the working dir to setwd()
# see: http://tuxette.nathalievilla.org/?p=1696
# https://cran.r-project.org/doc/manuals/R-intro.html#Scripting-with-R

# see also
# https://cwcode.wordpress.com/2013/04/16/the-joys-of-rscript/
# http://stackoverflow.com/questions/14167178/passing-command-line-arguments-to-r-cmd-batch
# https://shihho.wordpress.com/2012/11/30/r-how-to-run-r-scripts-in-batch-mode-with-arguments/

# get the usr arguments
args = commandArgs(trailingOnly=TRUE)

# test if there is at least one argument: if not, return an error
if (length(args) < 5) {
  print_help()
  stop("requires at least 5 arguments: <full_path2_ph_files> <Runmode> <aln_ext> <tree_ext> <genoFlag[0|1]> [basegraphs|ggplots]\n", call.=FALSE)
}  else if (length(args) >= 5 & length(args) < 7 ) {
  sprintf("Will run with args: %s, %s, %s, %s, %s, %s", args[1], args[2], args[3], args[4], args[5], args[6])
}

# set the working dir containing the tree files
cat("will setwd as: ", args[1], "\n")
setwd(args[1])
runmode <- args[2]
aln_ext <- args[3]
tree_ext <- args[4]
genoFlag <- args[5]
if (length(args) == 5) {
   GraphSyst <- "basegraphs"
}else if (length(args) == 6){
  GraphSyst <- as.character(args[6])
}  


#-------------------------------------------------#
#>>> PARALLEL BOXPLOTS FOR TREE SUPPORT VALUES <<<#
#-------------------------------------------------#

# 1. get all tree files with extension .ph
#  IMPORTANT NOTE: the sequences should be collapsed to haplotypes
#  but beware that the R code below complains when the tree labels 
#  contain '#' symbols, as introduced by collapse2haplotypes.pl
#  so use a line like the following to change # by - in the collapsed fastas (or on the trees)
# for file in *FAS; do collapse2haplotypes.pl $file | awk '{print $1, $2}' | perl -pe 'if(/^>/){ s/\h/_/ }' | sed 's/#/-/ > ${file}UNIQ; done
# for file in *ph; do sed 's/#/-/ $file > ${file}ed; done

# 1. get the list of tree files
# >>> pass the tree_ext arg to list.files() funct as a regex
tree_rgx <- paste("\\.", tree_ext, "$", sep = "")
files <- list.files(pattern = tree_rgx)
# check there are tree files with tree_ext extension in the working directory
if(length(files) == 0) stop("There are no tree files with", tree_ext, " extension in the working directory! Will stop now ...")

# 2. initialize vars
supp.dfr <- c()
aln.len.dfr <- c()
supp.vals <- c()
br.lens <- c()
br.dfr <- c()
br.lens.dfr <- c()
fac <- c()
aln_files <- c()
aln_rgx <- c()

# 2.1. get the gene names out of the tree file names. strsplit produces a list; 
#      get the 1st elements of each vector out of the list
gene.names <- strsplit(files, '_')
genes <- as.character(lapply(gene.names, "[[", 1))

# 2.2 GET ALIGNMENT LENGTHs for each input file and save them into aln.len.dfr
aln_rgx <- paste("\\.", aln_ext, "$", sep = "")
aln_files <- list.files(pattern = aln_rgx)
for (i in 1:length(aln_files) ){ 
    aln <- seqinr::read.alignment(file = aln_files[i], format = "fasta")
    aln.lens <- as.numeric(nchar(aln$seq[[1]])); aln.dfr <- data.frame(genes[i], aln.lens)
    aln.len.dfr <- rbind(aln.len.dfr, aln.dfr) 
}

# 2.3 add colnames to the cleaned dfr
colnames(aln.len.dfr) <- c("loci", "aln_length")

# 3. loop over files, read.trees and get support values; add those and corresponding factor level to intermediate dfr; finally rbind(supp.dfr, dfr)
#for (i in 1:length(files) ){ tr <- read.tree(files[i]); supp.vals <- as.numeric(tr$node.label); fac <- rep(i, length(supp.vals)); dfr <- data.frame(fac, supp.vals); supp.dfr <- rbind(supp.dfr, dfr) }
for (i in 1:length(files) ){ tr <- ape::read.tree(files[i]); supp.vals <- as.numeric(tr$node.label); fac <- rep(genes[i], length(supp.vals)); dfr <- data.frame(fac, supp.vals); supp.dfr <- rbind(supp.dfr, dfr) }

# 3.1 remove any rows containing NAs in supp.dfr
supp.dfr <- supp.dfr[complete.cases(supp.dfr), ]

# 3.2 add colnames to the cleaned dfr
colnames(supp.dfr) <- c("loci", "support.values")

# 3.3 use the aggregate function to plot the median support values
# notice the use of plyr::each()
aggr.supp.val.stats.dfr <- aggregate(support.values ~ loci, supp.dfr, each(median, mean, var, sd, length) )
aggr.supp.val.median.dfr <- aggregate(support.values ~ loci, supp.dfr, median)
aggr.supp.val.mean.dfr <- aggregate(support.values ~ loci, supp.dfr, mean)
aggr.supp.val.mean.dfr <- merge(aggr.supp.val.mean.dfr, aln.len.dfr, by = "loci")
aggr.supp.val.mean.dfr$signal <-  c(100 * (aggr.supp.val.mean.dfr$support.values / aggr.supp.val.mean.dfr$aln_length) )

med_supp_val <- median(aggr.supp.val.mean.dfr$support.values)
mean_supp_val <- mean(aggr.supp.val.mean.dfr$support.values)
quant070_supp_val <- quantile(aggr.supp.val.mean.dfr$support.values, probs = 0.70)
quant030_supp_val <- quantile(aggr.supp.val.mean.dfr$support.values, probs = 0.30)
stdevs_supp_val <- mean_supp_val + c(-2, -1, +1, +2)*sd(aggr.supp.val.mean.dfr$support.values)

write.table(aggr.supp.val.mean.dfr[order(aggr.supp.val.mean.dfr$support.values, decreasing = T), ], file="sorted_aggregated_support_values4loci.tab", row.names = FALSE, sep ="\t")
checkFileCreated("sorted_aggregated_support_values4loci.tab")

# save orignal par() params to opar
# no.readonly = TRUE; only parameters are returned which can be set by a subsequent par() call _on the same device_.
opar <- par(no.readonly = TRUE) 


# 3.4 order the dfr by support.values
if(genoFlag == 1){
  aggr.supp.val.mean.decr.ord.dfr <- aggr.supp.val.mean.dfr[order(aggr.supp.val.mean.dfr$support.values, decreasing = TRUE), ]
  no.trees <- dim(aggr.supp.val.stats.dfr)[1]

  # Make scatterplots; here need to make a faceted plot, to include boxplot and distribution wiht density
  if(GraphSyst == "basegraphs"){
      svg("scatterplot_for_gene_tree_support_values.svg")
      par(mfrow=c(2,1))
      cols <- ifelse(aggr.supp.val.mean.dfr$support.values >= quant070_supp_val, "black", "grey")
      plot(aggr.supp.val.mean.dfr$support.values, xaxt='n', main=sprintf("Mean SH-support values for %s gene trees (+-1 & 2 SD)", no.trees), xlab="loci", ylab="SH-support values", col=cols )
      abline(h=mean_supp_val, lwd = 2, col = "black") # "lty = 1# solid"
     # abline(h=med_supp_val, lty = "dotted", lwd = 2, col = "black") # "lty = "dotted" # solid"
      abline(h=stdevs_supp_val, lty = 2, lwd = 2, col = "grey") # "dashed"
      hist(aggr.supp.val.mean.dfr$support.values, 30, main ="Histogram of mean SH-support values", xlab="SH-support values", prob = TRUE)
      lines(density(aggr.supp.val.mean.dfr$support.values))
      dev.off()
      par(opar)
  }else if (GraphSyst == "ggplots"){
    # see http://www.cookbook-r.com/Graphs/Axes_%28ggplot2%29/#setting-tick-mark-labels
    # http://www.cookbook-r.com/Graphs/Titles_%28ggplot2%29/
    # ggplot(data=aggr.supp.val.mean.dfr, aes(y=support.values, x = loci) ) + geom_point() + theme(axis.title.x = element_blank(),  axis.text.x = element_blank() ) + labs(title="Median SH-support values for no.trees gene trees")
    svg("scatterplot_for_gene_tree_support_values.svg")
    print("GrahSyst=ggplot ==> scatterplot_for_gene_tree_support_values.svg")
    g <- ggplot(data = aggr.supp.val.mean.dfr, aes(y = support.values, x = loci) )
    g + geom_point() + theme(axis.title.x = element_blank(),  axis.text.x = element_blank() ) + ggtitle("Median SH-support values for gene trees") + theme(plot.title = element_text(lineheight=.8, face="bold") )
    dev.off()
  }
  checkFileCreated("scatterplot_for_gene_tree_support_values.svg")
  
    aggr.supp.val.median.top100.dfr <- head(aggr.supp.val.mean.decr.ord.dfr, 100)
    sink(file = "top100_median_support_values4loci.tab", type=c("output") )
    print(aggr.supp.val.median.top100.dfr)
    sink()
    checkFileCreated("top100_median_support_values4loci.tab")
}

# 4. AOV and plot group means with confidence inetaval (p = 0.95)
if(genoFlag < 1){
  fit.aov <- aov(support.values ~ loci, supp.dfr)
  sink(file = "summary_fit_AOV_support_values.out", type=c("output") )
  summary(fit.aov)
  sink()
  checkFileCreated("summary_fit_AOV_support_values.out")
}

# 5. prepare diverse plots
if(genoFlag < 1){
  svg(file="meanplot_for_bipartiton_support_values_with_095CI.svg")
  gplots::plotmeans(support.values ~ loci, supp.dfr, main=c("Mean tree support values/locus and 95% CIs") )
  dev.off()
  checkFileCreated("meanplot_for_bipartiton_support_values_with_095CI.svg")

  # 4. Plot the parallel boxplots; note the use of notch and names
  svg(file="parallel_boxplots_for_bipartiton_support_values.svg")
  if(GraphSyst == "basegraphs"){
    boxplot(support.values ~ loci, data = supp.dfr, notch = TRUE, names=genes, las = 2, main="SH-like support values", ylab="SH-support values")
    dev.off()
  }else if (GraphSyst == "ggplots"){
    ggplot(data=supp.dfr, aes(y = support.values, x = loci)) + geom_boxplot(aes(group = loci), notch = TRUE)
    dev.off()
  }
  checkFileCreated("parallel_boxplots_for_bipartiton_support_values.svg")

  # 5. loop over files, read.trees and get support values; add those and corresponding factor level to intermediate dfr; finally rbind(supp.dfr, dfr)
  for (i in 1:length(files) ){ tr <- read.tree(files[i]); br.lens <- tr$edge.length; fac <- rep(i, length(br.lens)); br.dfr <- data.frame(fac, br.lens); br.lens.dfr <- rbind(br.lens.dfr, br.dfr) }
  
  svg(file="parallel_boxplots_for_branch_lengths.svg")
  boxplot(br.lens ~ fac, data = br.lens.dfr, notch = FALSE, names=genes, las = 2, main="Tree branch lengths boxplot analysis", ylab="tree branch length")
  dev.off()
  checkFileCreated("parallel_boxplots_for_branch_lengths.svg")
}

# exit here if runmode == 1
if(runmode == 1) q(save = "no")

#---------------------------------------------------------------------------------------------------------------#
#>>> Compute Robinson-Foulds distances between the concatenated phylo and the individual (source) gene trees <<<#
#---------------------------------------------------------------------------------------------------------------#

# 0. get the list of "\\.ph$" files
tree_rgx <- paste("\\.", tree_ext, "$", sep = "")
files <- list.files(pattern = tree_rgx)
# check there are tree files with .ph extension in the working directory
if(length(files) == 0) stop("There are no tree files with .ph extension in the working directory! Will stop now ...")

# 1. filter out the concatenated and the gene trees and make sure it exists or stop
concat_tree_rgx <- paste("^concat.*", "\\.", tree_ext, "$", sep = "")
concat.tree.file <- str_extract(files, concat_tree_rgx)
concat.tree.file <- concat.tree.file[!is.na(concat.tree.file)]
if( length(concat.tree.file) < 1 ) stop("There is no concat file in the working directory! Will exit now ...") 

concat_tree <- read.tree(concat.tree.file)
gene.tree.files <- files[grep("concat.*", files, invert = TRUE)]

# 2. extract the gene names from the file names; see p. 103-104 R in action
gene.names <- strsplit(gene.tree.files, '_')
#genes <- as.character(lapply(gene.names, "[[", 1))
genes <- sapply(gene.names, "[", 1)

# 3. Compute the RF-scores 
# print the out once to screen
if(genoFlag == 0){
for ( i in 1:length(gene.tree.files)){ gene.tree <- read.tree(gene.tree.files[i]); print(dist.topo(concat_tree, gene.tree, "score")) }
}

if(genoFlag == 0){
  # 3.1 capture the output to file
  concat_vs_gene_scoreMat <- capture.output(for ( i in 1:length(gene.tree.files)){ gene.tree <- read.tree(gene.tree.files[i]); print(gene.tree.files[i]); print(dist.topo(concat_tree, gene.tree, "score")) })
  sink("gene_vs_concat_scoreMat.out", type=c("output") )
   # print(genes)
  print(as.data.frame(concat_vs_gene_scoreMat))
  sink()
  checkFileCreated("gene_vs_concat_scoreMat.out")
}

# 3.2 capture the RF-scores in a vector to plot a barplot; take the gene names from genes 
topo.vec <-c()
for ( i in 1:length(gene.tree.files)){ gene.tree <- read.tree(gene.tree.files[i]); top.dist<-(dist.topo(concat_tree, gene.tree, "score")); topo.vec <-c(topo.vec, top.dist) }

RF.dfr <- data.frame(genes,topo.vec)
colnames(RF.dfr) <- c("loci", "RF_dist")

# 3.3 print summary_stats_RFdist summary_stats_aggr_suppVal_medians to file
sink( file = "summary_stats_RFdist.out", type=c("output") )
summary(RF.dfr$RF_dist)
checkFileCreated("summary_stats_RFdist.out")

sink(file = "summary_stats_aggr_suppVal_medians.out", type=c("output") )
summary(aggr.supp.val.mean.dfr$support.values)
checkFileCreated("summary_stats_aggr_suppVal_medians.out")

# 3.4 compute the quantiles RF.30quant aggr.supp.val.mean.070quant to filter out best markers
#RF.30quant <- quantile(RF.dfr$RF_dist, probs = 0.10)
RF.30quant <- quantile(RF.dfr$RF_dist, probs = 0.3)
RFdist_mean <- mean(RF.dfr$RF_dist)
RFdist_stdevs <- RFdist_mean + c(-2, -1, +1, +2)*sd(RF.dfr$RF_dist)

#aggr.supp.val.mean.070quant <- quantile(aggr.supp.val.mean.dfr$support.values, probs = 0.70)
#aggr.supp.val.median.signal.075quant <- quantile(aggr.supp.val.mean.dfr$signal, probs = 0.75)

# 4. merge support and RFdist dfrs.
support_plus_RFdist.dfr <- merge(aggr.supp.val.mean.dfr, RF.dfr, by="loci")
support_plus_RFdist_ord.dfr <- support_plus_RFdist.dfr[order(support_plus_RFdist.dfr$RF_dist, decreasing = FALSE), ]
write.table(support_plus_RFdist_ord.dfr, file = "gene_trees2_concat_tree_RF_distances.tab", row.names = FALSE, sep = "\t" )
checkFileCreated("gene_trees2_concat_tree_RF_distances.tab")

# 5. filter dfr for high qual loci ( supp.val >= 0.75 quantile &  RF dist <= 0.10 quantile)
if(genoFlag == 1){
  # & support_plus_RFdist_ord.dfr$signal >= aggr.supp.val.median.signal.075quant
  top_markers.dfr <- support_plus_RFdist_ord.dfr[support_plus_RFdist_ord.dfr$support.values >= mean_supp_val & support_plus_RFdist_ord.dfr$RF_dist <= RF.30quant & support_plus_RFdist_ord.dfr$aln_length > 180, ]
  no.mark <- nrow(top_markers.dfr)
  #file_name <- sprintf("top %s markers suppVal gt075quant RFdist lt 010quant.out", no.mark)
  #file_name <- str_replace_all(string = file_name, pattern = "\\s", replacement = "_") 
  #sink(file = file_name)
  #print(top_markers.dfr)
  #sink()
  #checkFileCreated(file_name)
  
  tab_file_name <- sprintf("top %s markers suppVal gt70quant RFdist lt30quant.tab", no.mark)
  tab_file_name <- str_replace_all(string = tab_file_name, pattern = "\\s", replacement = "_") 
  write.table(top_markers.dfr, file=tab_file_name, row.names = FALSE, sep = "\t")
  checkFileCreated( tab_file_name ) # does not check it! neither with or withour print; complains with cat!?
}

# 6. barplots/scatterplots of RF-dist_of_gene_trees2concat_phylo
if(genoFlag == 0){
  svg(file="RF-dist_of_gene_trees2concat_phylo.svg")
  if(GraphSyst == "basegraphs"){
    barplot(topo.vec, names=genes, las = 2, main="RF-dist of gene-trees to concat-phylogeny", ylab="RF-dist")
    dev.off()
  }else if(GraphSyst == "ggplots"){
    ggplot(data=RF.dfr, aes(x=loci, y = RF_dist)) + geom_bar(stat = "identity")
    dev.off()
  }
  dev.off()
  checkFileCreated("RF-dist_of_gene_trees2concat_phylo.svg")
}else if(genoFlag == 1){
  svg("scatterplot_RF-dist_of_gene_trees2concat_phylo.svg")
  par(mfrow=c(2,1))
  cols <- ifelse(support_plus_RFdist.dfr$RF_dist <= RF.30quant, "black", "grey")
  plot(support_plus_RFdist.dfr$RF_dist, xaxt='n', main = sprintf("RF-dist. of %s gene trees to concat. tree (mean-+1&2SD)", no.trees), xlab="loci", ylab="RF-distance", col=cols )
  abline(h=RFdist_mean, lwd = 2, col = "black") # "lty = 1# solid"
  abline(h=RFdist_stdevs, lty = 2, lwd = 2, col = "grey") # "dashed"
  hist(RF.dfr$RF_dist, 30, main ="Histogram and density estimate of RF distances", xlab="RF-distance", prob = TRUE)
  lines(density(RF.dfr$RF_dist))
  dev.off()
  par(opar)
  checkFileCreated("scatterplot_RF-dist_of_gene_trees2concat_phylo.svg")
}

# exit without saving workspace
q(save = "no")

