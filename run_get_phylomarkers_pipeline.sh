#!/usr/bin/env bash

#: PROGRAM: run_get_phylomarkers_pipeline.sh
#: AUTHORS: Pablo Vinuesa, Center for Genomic Sciences, UNAM, Mexico
#:          http://www.ccg.unam.mx/~vinuesa/
#           Bruno Contreras Moreira, EEAD-CSIC, Zaragoza, Spain
#           https://digital.csic.es/cris/rp/rp02661/
#
#: PROJECT START: April 2017; This is a wrapper script to automate the whole process of marker selection and downstream analyses.
#
#: AIM: select optimal molecular markers for phylogenomics and population genomics from orthologous gene clusters computed by get_homologues
#
#: OUTPUT: multiple sequence alignments (of protein and DNA sequences) of selected markers, gene and supermatrix phylogenies, 
#          along with graphics and tables summarizing the filtering procedure. 
#          

progname=${0##*/} # run_get_phylomarkers_pipeline.pl
VERSION='1.9.9.2_31Dic17' # grep '^Gamma20LogLk' instead of ML_Lengths2 to print the FastTree lnL score to STDOUT
    # '1.9.9.1_23Dic17'  # 1.9.9.1_23Dic17: renamed the iqtree binary to iqtree-omp to be explicit about the multicore version
    # 1.9.9.0_22Dic17: added IQ-tree searching option for the concatenated alignment, controlled with new options -A, -N and -S
    # 1.9.8.4_17Nov17: improved/expanded -h help message; thorough and consistent tidying of directories; cleanup of code comments
    # 1.9.8.3_17Nov17: added get_homologues manual url to ERROR message to better assist users
    # 1.9.8.2_17Nov17: another sanity check: make sure there are equal number of fna and faa files to start working on
    # 1.9.8.1_17Nov17: fixed name of the add_labels2tree.pl in one of the calls with -R 1 -t PROT and code cleanup
    # 1.9.8_17Nov17' to avoid problems with old versions of locally insalled binaries or scripts, we set $bindir and $distrodir in front of PATH
    # 1.9.7.5_15Nov17 added check for minimal versions of clustalo FastTree parallel and paup
    # 1.9.7.5_15Nov17 fixed check for $HOME/bin in $PATH
    # 1.9.7.4_15Nov17 matched paths of results dir and log
    # 1.9.7.3_14Nov17 fixed bug/typo in get_script_PID()
    # 1.9.7.2_14Nov17 popGen dir cleanup + top_X_markers dir for -R 1 -t PROT
    #1.9.7.1_14Nov17: added strain composition check on f?aed files to make sure each one contains a single instance for the same number of strains
    #		      This is a critical check to avoid problems with inparalogues in some fastas if get_homologues.pl was run with -t X and without -e.
    #         Trees could be mislabeled in that case, and some alignments will most likely contain a different strain composition, 
    #         generating a chimaeric concatenated file; A useful ERROR message is printed (run compare_clusters.pl with -t NUM_GENOMES.
    #		      Also calls run_parallel_cmmds.pl with parallel --gnu; some code cleanup
    # 1.9.7_14Nov17: now using parallel instead of pexec binary, which failed in CentOS  
    # 1.9.6.4_12Nov17: prints total number of trees with < 1 internal branches; 
    # corrected regex so that all trees and alns with < 1 int br are removed (do not pass to downstream analyses)
    # 1.9.6.3_11Nov17:added -R${runmode} to dir_suffix
    # 1.9.6.2_11Nov17: all DEBUG|VERBOSE messages written to logfile
    # v1.9.6_11nov17: fixed code that writes Phi-test results to Phi_results_11nov17.tsv; 
    #		      prints the Warning: will remove file* because it has < 5 branches! to log file (no only STDOUT)
    # v1.9.5_8Nov17: the number of sequences in input FASTA files is checked upfront
    # v1.9.5_7Nov17: prints a logfile for FastTree runs and captures the lnL value of the corresponding tree; 
    #		     more detailed info reported on alignments with too few informative sites for the Phi test to be work on
    # v1.9.4_6Nov17: added error checking for Phi run when sequences lack enough polymorphism for the test to work
    # v1.9.3_10Oct17: improved error message when kdetrees fail
    # v1.9.2_1Jun17: code cleanup > removed comments and print_development_notes(); append the $distrodir/lib/perl to PERL5LIB and export
    # v1.9.1_1Jun17: more sensible directory cleanup ...
    # v1.9_31May17: improved progress messages and directory cleanup ...
    # 1.8.4_30May17: prepended $ditrodir/ to perl scripts that use FindBin; so that it can find the required libs in $ditrodir/lib/perl
    # Added -n $n_cores flag, which is passed to run_parallel_cmmds.pl '' $n_cores, so that it runs on MacOSX!!! <<< Thanks Alfredo!
    #	 automatically set n_cores=no_proc if [ -z $n_cores ]
    # v1.8.1_24May17 fixed problmes with @INC searching of rename.pl by prepending $distrodir/rename.pl
    #v1.8_23May17. Added R code in count_tree_branches() to use local_lib if ape is not installed systemwide; 
    #	 searches and prints the number of available cores on HOSTNAME 
    #	 exports R_LIBS="$R_LIBS:$distrodir/lib/R" to fix issues with library paths in R scripts
    # 1.7_17May17. Changed exit for Warning when pal2nal does not return a codon alignment! so that the pipeline can proceed
    # v1.6_15May17 added extensive debugging messages throughout the code for easier debugging; activated the -V flag
    # v1.5 fixed set_bindirs and check_homebinpath(), to export to PATH; 
    #	   fully tested on a new linux account without $HOME/bin dir using freshly cloned distro; Note that R and Perl libs were already in ENV
    # v1.3 further refinement in set_bindirs() and check_homebinpath(), validated on yaxche; minor code cleanup
    # v1.2_13May17 refined the logic of set_bindirs(); added get_start_time(); improved error checking code, including get_script_PID()
    # fixed a bug in -t PROT 
    
    #v1.1_13May17; Major update to facilitate installation by users: added set_pipeline_environment(), check_homebinpath(), set_bindirs(), 
    #		   which ckeck and setup the ENVIRONMENT for the pipeline. Added new bin/ and pop_gen_tables/ dirs for consistency. 
    #		   The bin/linux dir contains the 2cnd party binaries compiled for 64bit linux machines. Need to compile for darwin
    
    #v1.0_11May17; git commited; first version on GitHub: https://github.com/vinuesa/get_phylomarkers
    # v0.9_10May17 important speed improvement due to running pal2nal.pl and run_parallel_molecClock_test_with_paup.sh
    #		    in parallel with run_parallel_cmmds.pl 
    # v0.8_9May17, added -R 1|2, for Phylo|popGen, based on popGen_summStats.pl and pre-computed Tajima's D and Fu-Li crit value tables
   #	    with missing values predicted by lm() in R, based on the original critical values in Tajima's 89 and  Fu and Li 1993 papers
   #	    The critical CI values are gathered by get_critical_TajD_values() and get_critical_FuLi_values()
   # v0.7_2May17 added -e min_no_ext_branches -c $codontable and -C print_codontable() for pal2nal.pl, which is now run from within a loop
   #		 runs compute_suppValStas_and_RF-dist.R in the top_X_markers dir
   # v0.6_1May17 added get_script_PID() and count_tree_branches(); added -t BOTH
   # v0.5_29April17 Added -t PROT
   # v0.4_27April17 Added print_usage_notes()
# Set GLOBALS
DEBUG=0
wkdir=$(pwd) #echo "# working in $wkdir"

DATEFORMAT_SHORT="%d%b%y"
TIMESTAMP_SHORT=$(date +${DATEFORMAT_SHORT})

DATEFORMAT_HMS="%H:%M:%S"
TIMESTAMP_HMS=$(date +${DATEFORMAT_HMS})

TIMESTAMP_SHORT_HMS=$(date +${DATEFORMAT_SHORT}-${DATEFORMAT_HMS})

#>>> set colors in bash 
# ANSI escape codes
# Black        0;30     Dark Gray     1;30
# Red          0;31     Light Red     1;31
# Green        0;32     Light Green   1;32
# Brown/Orange 0;33     Yellow        1;33
# Blue         0;34     Light Blue    1;34
# Purple       0;35     Light Purple  1;35
# Cyan         0;36     Light Cyan    1;36
# Light Gray   0;37     White         1;37

RED='\033[0;31m'
LRED='\033[1;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
LBLUE='\033[1;34m'
CYAN='\033[0;36m'
NC='\033[0m' # No Color => end color
#printf "I ${RED}like${NC} ${GREEN}Stack Overflow${NC}\n"


#---------------------------------------------------------------------------------#
#>>>>>>>>>>>>>>>>>>>>>>>>>>>> FUNCTION DEFINITIONS <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<#
#---------------------------------------------------------------------------------#

function print_start_time()
{
   echo -n "[$(date +%T)] "
}
#----------------------------------------------------------------------------------------- 

function set_pipeline_environment()
{
  if [[ "$OSTYPE" == "linux-gnu" ]]
  then
    scriptdir=$(readlink -f ${BASH_SOURCE[0]})
    distrodir=$(dirname $scriptdir) #echo "scriptdir: $scriptdir|basedir:$distrodir|OSTYPE:$OSTYPE"
    bindir=$distrodir/bin/linux
    OS='linux'
    no_cores=$(awk '/^processor/{n+=1}END{print n}' /proc/cpuinfo)
  elif [[ "$OSTYPE" == "darwin"* ]]
  then
    distrodir=$(cd "$(dirname "$0")"; pwd)
    bindir=$distrodir/bin/macosx-intel
    OS='darwin'
    no_cores=$(sysctl -n hw.ncpu)
  fi
  echo "$distrodir $bindir $OS $no_cores"
}
#----------------------------------------------------------------------------------------- 

function check_dependencies()
{
    
    # check if scripts are in path; if not, set flag
    for prog in bash R perl awk cut grep sed sort uniq Rscript
    do
       bin=$(type -P $prog)
       if [ -z $bin ]; then
          echo
          printf "${RED}# ERROR: system binary $prog is not in \$PATH!${NC}\n"
	  printf "${LRED}  >>>  you will need to install $prog for $progname to run on $HOSTNAME ${NC}\n"
	  printf "${LRED}  >>>  will exit now ... ${NC}\n"
	  exit 1
       fi
    done	  
}    
#----------------------------------------------------------------------------------------- 

function check_scripts_in_path()
{
    distrodir=$1
    not_in_path=0
    homebinflag=0
    homebinpathflag=0

    [ $DEGUB ] && echo "check_scripts_in_path() distrodir:$distrodir"
    
    bash_scripts=(run_parallel_molecClock_test_with_paup.sh )
    perl_scripts=( run_parallel_cmmds.pl add_nos2fasta_header.pl pal2nal.pl rename.pl concat_alignments.pl add_labels2tree.pl convert_aln_format_batch_bp.pl popGen_summStats.pl convert_aln_format_batch_bp.pl )
    R_scripts=( run_kdetrees.R compute_suppValStas_and_RF-dist.R )
    
    # check if scripts are in path; if not, set flag
    for prog in "${bash_scripts[@]}" "${perl_scripts[@]}" "${R_scripts[@]}"
    do
       bin=$(type -P $prog)
       if [ -z $bin ]; then
          echo
          printf "${LRED}# WARNING: script $prog is not in \$PATH!${NC}\n"
	        printf "${CYAN}  >>>  Will generate a symlink from $HOME/bin or add it to \$PATH ${NC}\n"
	        not_in_path=1
       fi
    done	  
    
    # if flag $not_in_path -eq 1, then either generate symlinks into $HOME/bin (if in $PATH) or export $distrodir to PATH
    if [ $not_in_path -eq 1 ]
    then
       if [ ! -d $HOME/bin ]
       then
            printf "${CYAN} found no $HOME/bin directory for $USER ...${NC}\n"
	    printf "${CYAN} ... will update PATH=$distrodir:$PATH ${NC}\n"
	    export PATH="${distrodir}:${PATH}" # prepend $ditrodir to $PATH
       else
           homebinflag=1
       fi
    
       # check if $HOME/bin is in $PATH
       echo $PATH | sed 's/:/\n/g'| grep "$HOME/bin$" &> /dev/null
       if [ $? -eq 0 ]
       then
             homebinpathflag=1

             printf "${CYAN} Found the $HOME/bin for $USER in \$PATH ...${NC}\n"
             printf "${CYAN} ... will generate symlinks in $HOME/bin to all scripts in $distrodir ...${NC}"
             ln -s $distrodir/*.sh $HOME/bin &> /dev/null
             ln -s $distrodir/*.R $HOME/bin &> /dev/null
             ln -s $distrodir/*.pl $HOME/bin &> /dev/null
             #ln -s $distrodir/rename.pl $HOME/bin &> /dev/null
       else
           printf "${CYAN} Found the $HOME/bin for $USER, but it is NOT in \$PATH ...${NC}\n"
           printf "${CYAN} updating PATH=$PATH:$distrodir ${NC}\n"
	   export PATH="${distrodir}:${PATH}" # prepend $distrodir to $PATH
       fi
    fi
    #echo "$homebinflag $homebinpathflag"
}
#----------------------------------------------------------------------------------------- 

function set_bindirs()
{  
    # receives: $bindir $homebinpathflag
    bindir=$1
#   not_in_path=1

#   # prepend $bindir to $PATH in any case, to ensure that the script runs with the distributed binaries in $bindir
    export PATH="${bindir}:${PATH}" 

#    bins=( clustalo FastTree parallel Phi paup consense )
#
#    for prog in "${bins[@]}" 
#    do
#       bin=$(type -P $prog)
#       if [ -z $bin ]
#       then
#          echo
#          printf "${LRED}# $prog not found in \$PATH ... ${NC}\n"
#	        not_in_path=1
#       fi	  
#   done	  
# 
#    # check if scripts are in path; if not, set flag 
#   if [ $not_in_path -eq 1 ]
#   then
#   	   printf "${CYAN} updating PATH=$PATH:$bindir ${NC}\n"
#   	   #export PATH=$PATH:$bindir # append $HOME/bin to $PATH, (at the end, to not interfere with the system PATH)  
#	   # we do not export, so that this PATH update lasts only for the run of the script, 
#	   # avoiding a longer-lasting alteration of $ENV; 
#	   export PATH="${bindir}:${PATH}" # prepend $bindir to $PATH to ensure that the script runs with the distributed binaries
#   fi
   #echo $setbindir_flag
}
#----------------------------------------------------------------------------------------- 

function print_usage_notes()
{
   cat <<USAGE
   
   $progname v$VERSION important usage notes.
   
   1. Start the run from within the directory holding core gene clusters generated by get_homologues.pl -e 
      or compare_clusters.pl
   
      NOTE: Both .faa and .fna files are required to generate codon alignments from DNA fasta files. This
       means that two runs of compare_clusters.pl (from the get_homologues package) are required, one of them 
       using the -n flag.
	    
   2. $progname is intended to run on a collection of single-copy sequence clusters from different species or strains. 
   
      NOTES: An absolute minimum of 4 distinct genomes are required. 
       However, the power of the pipeline for selecting optimal genome loci 
	     for phylogenomics improves when a larger number of genomes are available 
	     for analysis. Reasonable numbers lie in the range of 10 to 100 clearly
	     distinct genomes from multiple species of a genus, family, order or phylum.
	     The pipeline may not perform satisfactorily with too distant genome sequences,
	     particularly when sequences with significantly distinct nucleotide or aminoacid
	     compositions are used. This type of sequence heterogeneity is well known to 
	     cause systematic bias in phylogenetic inference. Otherwise, too distantly related
	     organisms, such as those from different phyla or even domains, are also not
	     properly handled by ${progname}.

   3. On the locus filtering criteria. $progname uses a hierarchical filtering scheme, as follows:
   
      i) Detection of recombinant loci. Codon or protein alignments (depending on runmode) 
       are first screened with Phi for the presence of potential recombinant sequences. 
	    It is a well established fact that recombinant sequences negatively impact 
	    phylogenetic inference when using algorithms that do not account for the effects 
	    of this evolutionary force. The permutation test with 1000 permutations is used
	    to compute the p-values. These are considerd significant if < 0.05.
    
      ii) Detection of trees deviating from expectations of the (multispecies) coalescent.
       The next filtering step is provided by the kdetrees test, which checks the distribution of
       topologies, tree lengths and branch lenghts. kdetrees is a non-parametric method for 
	    estimating distributions of phylogenetic trees, with the goal of identifying trees that 
	    are significantly different from the rest of the trees in the sample. Such "outlier" 
	    trees may arise for example from horizontal gene transfers or gene duplication 
	    (and subsequent neofunctionalization) followed by differential loss of paralogues among
	    lineages. Such processes will cause the affected genes to exhibit a history distinct 
	    from those of the majority of genes, which are expected to be generated by the 
	    (multispecies) coalescent as species or populations diverge. Alignments producing 
	    significantly deviating trees in the kdetrees test are discarded.
	 
	    * Parameter for controlling kdetrees stingency:
	    -k <real> kde stringency (0.7-1.6 are reasonable values; less is more stringent)
	              [default: $kde_stringency]
 
      iii) Phylogenetic signal content and congruence. The alignments passing the two previous
       filters are subjected to maximum likelihood (ML) tree searches with FastTree to 
	    infer the corresponding ML gene trees. The phylogenetic signal of these trees is 
	    computed from the Shimodaria-Hasegawa-like likelihood ratio test of branch support
	    values, which vary between 0-1. The support values of each internal branch or 
	    bipartition are parsed to compute the mean support value for each tree. Trees 
	    with a mean support value below a cutoff threshold are discarded. 
	    In addition, a consensus tree is computed from the collection of trees that 
	    passed filters i and ii, and the Robinson-Fould distance is computed between 
	    each gene tree and the consensus tree.  

         * Parameters controlling filtering based on mean support values.
         -m <real> min. average support value (0.7-0.8 are reasonable values) 
	           for trees to be selected [default: $min_supp_val]


      iv) On tree searching: $progname performs tree searches using the FastTree ML algorithm.
       This program meets an excellent compromise between speed and accuracy, runnig both
	    with DNA and protein sequence alignments. It computes the above-mentioned 
	    Shimodaria-Hasegawa-like likelihood ratio test of branch support values.
	    A limitation though, is that it implements only very few substitution models. 
	    However, for divergent sequences of different species within a bacterial taxonomic
	    genus or family, our experience has shown that almost invariably the GTR+G model
	    is selected by jmodeltest2, particularly when there is base frequency heterogeneity.
	    The GTR+G+CAT is the substitution model used by $progname calls of FastTree on codon
	    alignments. The gene trees are computed with high accuracy by performing a thorough
	    tree search, as hardcoded in the following FastTree call:
	 
	   -nt -gtr -gamma -bionj -slownni -mlacc 3 -spr 8 -sprlength 8 
	   
	   For concatenated codon alignments, which may take a considerable time (up to several hours)
	   for large datasets (~ 100 taxa and > 300 concatenated genes) the user can choose to run 
	   FastTree with at different levels of tree-search thoroughness: high|medium|low|lowest 
	 
	   high:   -nt -gtr -bionj -slownni -gamma -mlacc 3 -spr $spr -sprlength $spr_length
	   medium: -nt -gtr -bionj -slownni -gamma -mlacc 2 -spr $spr -sprlength $spr_length 
	   low:    -nt -gtr -bionj -slownni -gamma -spr $spr -sprlength $spr_length 
      lowest: -nt -gtr -gamma -mlnni 4
	 
	   where -s \$spr and -l \$spr_length can be set by the user. 
	   The lines above show their default values.
	 
	   For protein alignments, the search parameters are the same, only the model changes to lg
	 
	   high: -lg -bionj -slownni -gamma -mlacc 3 -spr $spr -sprlength $spr_length
	 
	   Please refer to the FastTree manual for the details. 

USAGE

exit 0

}
#----------------------------------------------------------------------------------------- 

function print_codontables()
{
  cat <<CODONTBL
    1  Universal code (default)
    2  Vertebrate mitochondrial code
    3  Yeast mitochondrial code
    4  Mold, Protozoan, and Coelenterate Mitochondrial code
       and Mycoplasma/Spiroplasma code
    5  Invertebrate mitochondrial
    6  Ciliate, Dasycladacean and Hexamita nuclear code
    9  Echinoderm and Flatworm mitochondrial code
   10  Euplotid nuclear code
   11  Bacterial, archaeal and plant plastid code
   12  Alternative yeast nuclear code
   13  Ascidian mitochondrial code
   14  Alternative flatworm mitochondrial code
   15  Blepharisma nuclear code
   16  Chlorophycean mitochondrial code
   21  Trematode mitochondrial code
   22  Scenedesmus obliquus mitochondrial code
   23  Thraustochytrium mitochondrial code

CODONTBL

exit 0

}
#----------------------------------------------------------------------------------------- 

function get_script_PID()
{
    # returns the PID of the script, run by USER
    prog=${0%.*} # remove the script's .extension_name
    #proc_ID=$(ps -eaf|grep "$prog"|grep -v grep|grep '-'|grep $USER|awk '{print $2}')
    proc_ID=$(ps aux|grep "$prog"|grep -v grep|grep '-'|grep $USER|awk '{print $2}')
    echo $proc_ID
    [ $DEBUG -eq 1 ] && echo "$progname PID is: $proc_ID"
}
#----------------------------------------------------------------------------------------- 

function install_Rlibs_msg()
{
   outfile=$1
   Rpackage=$2

   printf "${RED} ERROR: the expected outfile $outfile was not produced\n"  
   printf "       This may be because R package(s) $Rpackage are not properly installed.\n"
   printf "       Please run the script './install_R_deps.R' from within $distrodir to install them.\n"
   printf "       Further installation tips can be found in file Rscript install_R_deps.R\n ${NC}"  
   
   R --no-save --quiet <<RCMD 2> /dev/null
   
       print("Your R installation currently searches for packages in :")
       print(.libPaths())
RCMD

exit 1
  
}

#----------------------------------------------------------------------------------------- 

function count_tree_labels()
{
   # call R pacake 'ape' to count the number of labels in a tree
   ext=$1
   outfile=$2

   R --no-save --quiet <<RCMD

   pkg <- c("ape")
   
   local_lib <- c("$distrodir/lib/R")
   distrodir <-c("$distrodir")

   .libPaths( c( .libPaths(), "$distrodir/lib/R") )
   
   library("ape")
   
   if (!require(pkg, character.only=T, quietly=T)) {
       sprintf("# cannot load %s. Install it in %s using the command 'Rscript install_R_deps.R' from within %s",package,local_lib,distrodir)
       print("Your R installation currently searches for packages in :\n")
       print(.libPaths())
   }

   trees <- list.files(pattern = "$ext\$")
   sink("$outfile")

   for( i in 1:length(trees)){
      tr <- read.tree(trees[i])
      labels <- tr\$tip.label
      no_labels <- length(labels)
      cat(trees[i], "\t", no_labels, "\n")
   }
   sink()
RCMD

}
#----------------------------------------------------------------------------------------- 

function count_tree_branches()
{
   ext=$1
   outfile=$2

   [ $DEBUG -eq 1 ] && echo "# running count_tree_branches $ext $outfile"
    
   R --no-save --quiet <<RCMD &> /dev/null
  
   local_lib <- c("$distrodir/lib/R")
   distrodir <- c("$distrodir")

   .libPaths( c(.libPaths(), local_lib) )
   required_packages <- c("ape")
   
   for (pkg in required_packages) { 
     if (!require("ape", character.only=T, quietly=T)) {
       sprintf("# cannot load %s. Install it in %s using the command 'Rscript install_R_deps.R' from within %s",package,local_lib,distrodir)
       print("Your R installation currently searches for packages in :\n")
       print(.libPaths())
     }
   }
   library("ape")

   trees <- list.files(pattern = "$ext\$")
   sink("$outfile")
   cat("#Tree","\t", "n_leaf_lab", "\t", "n_zero_len_br" , "\t", "n_nodes", "\t", "n_br", "\t", "n_int_br", "\t", "n_ext_br", "\t", "is_binary_tree", "\n")

   no_int_br_counter <- 0

   for( i in 1:length(trees)){
      tr <- read.tree(trees[i])
      is_binary_tree <- is.binary.tree(tr)
      #no_branches <- length(tr\$edge.length) # equivalent to the line below
      no_branches <- dim(tr\$edge)[1]
      
      # we need to count the number of branches with lenght == 0 
      # to remove those from the total count of no_branches computed above
      n_zero_length_branches <- length(which(tr\$edge.length == 0)) 
      
      # this is the real num. of branches on our tree
      no_branches <- no_branches - n_zero_length_branches
      
      Nnodes <- tr\$Nnode  # these are the internal nodes
      Ntips <- length(tr\$tip.label) # this is the number of tree labels or taxa    
      
      # this is the number of nodes in a binary tree for n taxa/tip.label
      branches_in_binary_tree <- (2*Ntips - 3) 
      #int_branches_in_binary_tree <- branches_in_binary_tree - Ntips
      
      # Using the eq. above, we can now estimate the real num. of ext branches
      # taking into account the real num of branches, after removing those with zelo length 
      no_ext_branches <- floor((no_branches + 3) / 2) # floor() in case there is no real internal branch, we get -0.5
      no_int_branches <- ceiling(no_branches - no_ext_branches) # ceiling() in case there is no real internal branch, we get -0.5
      if ( no_int_branches < 1 ) no_int_br_counter <- no_int_br_counter + 1
      cat(trees[i], "\t", Ntips, "\t", n_zero_length_branches, "\t", Nnodes, "\t", no_branches, "\t", no_int_branches, "\t", no_ext_branches, "\t", is_binary_tree, "\n")
    }
   sink()
   
   # return this number from function to print in main script
   cat(no_int_br_counter)
   
RCMD

}
#----------------------------------------------------------------------------------------- 

function concat_alns()
{
    aln_ext=$1
    pPID=$2
    
    ls -1 *${aln_ext} > list2concat
    check_output list2concat $parent_PID
    
    if [ "$mol_type" == "DNA" ]
    then
        concat_file="concat_cdnAlns.fna"
    fi

    if [ "$mol_type" == "PROT" ]
    then
        concat_file="concat_protAlns.faa"
    fi
    
    $distrodir/concat_alignments.pl list2concat > $concat_file
    #concat_alignments.pl list2concat > $concat_file
    check_output $concat_file $pPID
    perl -ne 'if (/^#|^$/){ next }else{print}' $concat_file > ed && mv ed $concat_file
    check_output $concat_file $pPID 
}
#----------------------------------------------------------------------------------------- 

function fix_fastaheader()
{
   # extract the relevant fields from the long '|'-delimited fasta headers generated by get_homologues
   # to construct a shorte one, suitable for display as tree labels
   
   file=$1
   awk 'BEGIN {FS = "|"}{print $1, $2, $3}' $file|perl -pe 'if(/^>/){s/>\S+/>/; s/>\h+/>/; s/\h+/_/g; s/,//g; s/;//g; s/://g; s/\(//g; s/\)//g}' > ${file}ed;
   
   check_output ${file}ed $parent_PID

}
#----------------------------------------------------------------------------------------- 

function get_critical_TajD_values()
{
   no_seq=$1
   TajD_low=$(awk -v n=$no_seq '$1 == n' ${distrodir}/pop_gen_tables/TajD_predicted_CIs.tsv|cut -d' ' -f2)
   TajD_up=$(awk -v n=$no_seq '$1 == n' ${distrodir}/pop_gen_tables/TajD_predicted_CIs.tsv|cut -d' ' -f3)
   echo "$TajD_low $TajD_up"
}
#----------------------------------------------------------------------------------------- 

function get_critical_FuLi_values()
{
   no_seq=$1
   FuLi_low=$(awk -v n=$no_seq 'BEGIN{FS="\t"}$1 == n' ${distrodir}/pop_gen_tables/Fu_Li_predicted_CIs.tsv|cut -f2)
   FuLi_up=$(awk -v n=$no_seq 'BEGIN{FS="\t"}$1 == n'  ${distrodir}/pop_gen_tables/Fu_Li_predicted_CIs.tsv|cut -f3)
   echo "$FuLi_low $FuLi_up"
}
#----------------------------------------------------------------------------------------- 

function check_output()
{
    outfile=$1
    pPID=$2

    #>>> set color in bash 
    #  SEE: echo http://stackoverflow.com/questions/5947742/how-to-change-the-output-color-of-echo-in-linux
    #  And very detailed: http://misc.flogisoft.com/bash/tip_colors_and_formatting
    # ANSI escape codes
    # Black        0;30     Dark Gray     1;30
    # Red          0;31     Light Red     1;31
    # Green        0;32     Light Green   1;32
    # Brown/Orange 0;33     Yellow        1;33
    # Blue         0;34     Light Blue    1;34
    # Purple       0;35     Light Purple  1;35
    # Cyan         0;36     Light Cyan    1;36
    # Light Gray   0;37     White         1;37
    RED='\033[0;31m'
    GREEN='\033[0;32m'
    YELLOW='\033[1;33m'
    BLUE='\033[0;34m'
    CYAN='\033[0;36m'
    NC='\033[0m' # No Color => end color
    #printf "I ${RED}love${NC} ${GREEN}Stack Overflow${NC}\n"
    
    if [ -s $outfile ]
    then
         printf "${GREEN} >>> wrote file $outfile ...${NC}\n"
	 return 0
    else
        echo
	printf "${RED} >>> ERROR! The expected output file $outfile was not produced, will exit now!${NC}\n"
        echo
	exit 9
	[ $DEBUG -eq 1 ] && echo "check_output running: kill -9 $pPID"
	kill -9 $pPID
    fi
}
#----------------------------------------------------------------------------------------- 

function print_help()
{
   cat <<EOF
   $progname v.$VERSION usage:
   
   REQUIRED:
    -R <integer> RUNMODE
       1 select optimal markers for phylogenetics/phylogenomics (genomes from different species)
       2 select optimal markers for population genetics (genomes from the same species)
    -t <string> type of input sequences: DNA|PROT
    
   OPTIONAL:
     -h flag to print this short help notes
     -A <string> Algorithm for tree searching: <F|I> [FastTree|IQ-TREE]                            [default:$search_algorithm]                            
     -H flag to print additional usage Notes
     -c <integer> NCBI codontable number (1-23) for pal2nal.pl to generate codon alignment         [default:$codontable] 
     -C flag to print codontables
     -D flag to print debugging messages; please use if you encounter problems executing the code  [default: $DEBUG]
     -e <integer> select gene trees with at least (min. = 4) external branches                     [default: $min_no_ext_branches]
     -k <real> kde stringency (0.7-1.6 are reasonable values; less is more stringent)              [default: $kde_stringency]
     -K <integer> run molecular clock test on codon alignments                                     [default: $eval_clock]
     -l <integer> max. spr length (7-12 are reasonable values)                                     [default: $spr_length]
     -m <real> min. average support value (0.7-0.8 are reasonable values) for trees to be selected [default: $min_supp_val]
     -M <string> base Model for clock test (use one of: GTR|TrN|HKY|K2P|F81); uses +G in all cases [default: $base_mod]
     -n <integer> number of cores/threads to use                                                   [default: all cores]
     -N <integer> number of IQ-TREE searches to run [only active with -T high]                     [default: $nrep]
     -q <real> quantile (0.95|0.99) of Chi-square distribution for computing molec. clock p-value  [default: $q]
     -r <string> root method (midpoint|outgroup)                                                   [default: $root_method]
     -s <integer> number of spr rounds (4-20 are reasonable values) for FastTree tree searching    [default: $spr]
     -S <string> quoted 'comma-separated list' of base models to be evaluated 
              by ModelFinder for IQ-TREE, when -A I. IQ-TREE runs with the concatenated alignments
              <'JC,F81,K2P,HKY,TrN,TNe,K3P,K81u,TPM2,TPM2u,TPM3,TPM3u,
	      TIM,TIMe,TIM2,TIM2e,TIM3,TIM3e,TVM,TVMe,SYM,GTR'>              for DNA alignments    [default: $IQT_DNA_models] 
              <'BLOSUM62,cpREV,Dayhoff,DCMut,FLU,HIVb,HIVw,JTT,JTTDCMut,LG,
                mtART,mtMAM,mtREV,mtZOA,Poisson,PMB,rtREV,VT,WAG'>           for PROT alignments   [default: $IQT_PROT_models]       
     -T <string> tree search Thoroughness: high|medium|low|lowest                                  [default: $search_thoroughness]
     -V flag to activate verbose command execution lines                                           [default: $VERBOSITY]
     
   Invocation examples:
     1. default on DNA sequences: 
          $progname -R 1 -t DNA
     2. thorough FastTree searching and molecular clock analysis on DNA sequences:
          $progname -R 1 -t DNA -k 1.2 -m 0.7 -s 8 -l 10 -T high -K 1 -M HKY -q 0.95
     3. FastTree searching on a huge protein dataset
          $progname -R 1 -t PROT -m 0.6 -k 1.0 -T lowest
     4. To run the pipeline on a remote server, we recommend using the nohup command upfront, as shown below:
        nohup $progname -R 1 -t DNA -A I -S 'TNe,TVM,TVMe,GTR' -k 1.0 -m 0.7 -T high -N 5 &> /dev/null &	  
     
   NOTES
     1: run from within the directory holding core gene clusters generated by get_homologues.pl -e or
          compare_clusters.pl with -t no_genome_gbk files (core genome: all clusters with a single gene copy per genome)
     2: If you encounter any problems, please run the script with the -D -V flags added at the end of the command line, 
        redirect STOUT to a file and send us the output, so that we can better diagnose the problem.
	e.g. $progname -R 1 -t DNA -k 1.0 -m 0.7 -s 8 -l 10 -T high -K 1 -D -V &> ERROR.log
	       
EOF

exit

}
#----------------------------------------------------------------------------------------- 

#------------------------------------#
#----------- GET OPTIONS ------------#
#------------------------------------#

# Required
runmode=
mol_type=

# Optional, with defaults
search_thoroughness=medium
kde_stringency=1.5
min_supp_val=0.75
min_no_ext_branches=4
n_cores=
VERBOSITY=0
spr_length=8
spr=4
codontable=11 # bacterial by default, sorry for the bias ;)
base_mod=GTR
eval_clock=0
root_method=midpoint
tree_prefix=concat 
q=0.99

search_algorithm=F
IQT_DNA_models=GTR
IQT_PROT_models=LG
IQT_models=
nrep=10

while getopts 'c:e:k:K:l:m:M:n:N:p:q:r:s:t:A:T:R:S:hHCDV' OPTIONS
do
   case $OPTIONS in
   h)   print_help
        ;;
   A)   search_algorithm=$OPTARG
        ;;
   k)   kde_stringency=$OPTARG
        ;;
   c)   codontable=$OPTARG
        ;;
   C)	print_codontables
	;;
   D)   DEBUG=1
        ;;
   e)   min_no_ext_branches=$OPTARG
	;;
   K)   eval_clock=$OPTARG
        ;;
   l)   spr_length=$OPTARG
        ;;
   m)   min_supp_val=$OPTARG
        ;;
   M)   base_mod=$OPTARG
        ;;
   n)   n_cores=$OPTARG
        ;;
   N)   nrep=$OPTARG
        ;;
   t)   mol_type=$OPTARG
        ;;
   T)   search_thoroughness=$OPTARG
        ;;
   p)   tree_prefix=$OPTARG
        ;;
   q)   q=$OPTARG
        ;;
   r)   root_method=$OPTARG
        ;;
   s)   spr=$OPTARG
        ;;
   S)   IQT_models=$OPTARG
        ;;
   R)   runmode=$OPTARG
        ;;
   H)   print_usage_notes
        ;;
   V)   VERBOSITY=1
        ;;
   \:)   printf "argument missing from -%s option\n" $OPTARG
   	 print_help
     	 exit 2 
     	 ;;
   \?)   echo "need the following args: "
   	 print_help
         exit 3
	 ;;
    *)   echo "An  unexpected parsing error occurred"
         echo
         print_help
	 exit 4
	 ;;	 
   esac >&2   # print the ERROR MESSAGES to STDERR
done

shift $(($OPTIND - 1))


#-------------------------------------------------------#
# >>>BLOCK 0.1 SET THE ENVIRONMENT FOR THE PIPELINE <<< #
#-------------------------------------------------------#

# 0. Set the distribution base directory and OS-specific (linux|darwin) bindirs
env_vars=$(set_pipeline_environment) # returns: $distrodir $bindir $OS $no_proc
[ $DEGUG ] && echo "env_vars:$env_vars"
distrodir=$(echo $env_vars|awk '{print $1}')
bindir=$(echo $env_vars|awk '{print $2}')
OS=$(echo $env_vars|awk '{print $3}')
no_proc=$(echo $env_vars|awk '{print $4}')

[ $DEBUG -eq 1 ] && echo "distrodir:$distrodir|bindir:$bindir|OS:$OS|no_proc:$no_proc"

# 0.1 Determine if pipeline scripts are in $PATH; 
# if not, add symlinks from ~/bin, if available
check_scripts_in_path $distrodir

# 0.2  Determine the bindir and add prepend it to PATH and export
set_bindirs $bindir

# 0.3 append the $distrodir/lib/R to R_LIBS and export
export R_LIBS="$R_LIBS:$distrodir/lib/R"

# 0.4 append the $distrodir/lib/perl to PERL5LIB and export
export PERL5LIB="${PERL5LIB}:${distrodir}/lib/perl:${distrodir}/lib/perl/bioperl-1.5.2_102"

#--------------------------------------#
# >>> BLOCK 0.2 CHECK USER OPTIONS <<< #
#--------------------------------------#

# check for bare minimum dependencies: bash R perl awk cut grep sed sort uniq Rscript

logdir=$(pwd)

check_dependencies

if [ -z $runmode ]
then
       echo "# ERROR: no runmode defined!"
       print_help
       exit 1    
fi

if [ "$search_algorithm" != "I" -a "$search_algorithm" != "F" ]
then
       echo "# ERROR: search_algorithm "$search_algorithm" is not recognized!"
       print_help
       exit 1    
fi

if [ $min_no_ext_branches -lt 4 ]
then
    printf "${RED}>>> ERROR: -e has to be >= 4\n\n${NC}"
    print_help
    exit 1
fi

if [ -z $n_cores ]
then
     n_cores=$no_proc   
fi

if [ "$mol_type" != "DNA" -a "$mol_type" != "PROT" ] # "$mol_type" == "BOTH" not implemented yet
then
     printf "\n${RED}ERROR: -t must be DNA or PROT${NC}\n"
     print_help
     exit 1
fi

if [ -z $IQT_models ]
then
   [ "$mol_type" == "DNA" ] &&  IQT_models=$IQT_DNA_models
   [ "$mol_type" == "PROT" ] && IQT_models=$IQT_PROT_models
fi

if [ "$search_thoroughness" != "high" -a "$search_thoroughness" != "medium" -a "$search_thoroughness" != "low" -a "$search_thoroughness" != "lowest" ]
then
     printf "\n${RED}ERROR: -T must be lowest|low|medium|high${NC}\n"
     print_help
     exit 1
fi

if [ "$base_mod" != "GTR" -a "$base_mod" != "TrN" -a "$base_mod" != "HKY" -a "$base_mod" != "K2P" -a "$base_mod" != "F81" ]
then
     printf "\n${RED}ERROR: -M must be one of: GTR|TrN|HKY|K2P|F81${NC}\n"
     print_help
     exit 1
fi

if [ $eval_clock -gt 0 -a "$mol_type" != "DNA" ] # MolClock currently only with DNA
then
     printf "\n${RED}ERROR: -K 1 (evaluate clock) must be run on codon alignments with -t DNA${NC}\n"
     print_help
     exit 1
fi

if [ $runmode -gt 1 -a "$mol_type" != "DNA" ] # PopGen analysis currently only with DNA
then
     printf "\n${RED}ERROR: runmode $runmode must be run on codon alignments with -t DNA${NC}\n"
     print_help
     exit 1
fi

#---------------------#
# >>>> MAIN CODE <<<< #
#---------------------#

[ $DEBUG -eq 1 ] && echo "running on $OSTYPE" && echo "path contains: " && echo $PATH|sed 's/:/\n/g' 

start_time=$(date +%s)

parent_PID=$(get_script_PID $progname)
[ $DEBUG -eq 1 ] && echo "parent_PID:$parent_PID"

dir_suffix=R${runmode}t${mol_type}_k${kde_stringency}_m${min_supp_val}_s${spr}_l${spr_length}_T${search_thoroughness}_K${eval_clock}

printf "
 ${CYAN}>>> $(basename $0) vers. $VERSION run with the following parameters:${NC}
 ${YELLOW}Run started on $TIMESTAMP_SHORT_HMS under $OSTYPE on $HOSTNAME with $n_cores cores
 wkdir=$wkdir
 distrodir=$distrodir
 bindir=$bindir

 > General run settings:
      runmode=$runmode|mol_type=$mol_type|search_algorithm=$search_algorithm 
 > Filtering parameters: 
     kde_stringency=$kde_stringency|min_supp_val=$min_supp_val
 > FastTree parameters: 
     spr=$spr|spr_length=$spr_length|search_thoroughness=$search_thoroughness
 > IQ-TREE parameters: 
     IQT_models=$IQT_models|search_thoroughness=$search_thoroughness|nrep=$nrep
 > Molecular Clock parmeters: 
     eval_clock=$eval_clock|root_method=$root_method|base_model=$base_mod|ChiSq_quantile=$q

 DEBUG=$DEBUG|VERBOSITY=$VERBOSITY${NC}

" | tee -a ${logdir}/get_phylomarkers_run_${dir_suffix}_${TIMESTAMP_SHORT}.log 

#----------------------------------------------------------------------------------------------------------------
#>>>BLOCK 1. make a new subdirectory within the one holding core genome clusters generated by compare_clusters.pl
#    and generate symlinks to the faa and fna files (NOTE: BOTH REQUIRED). Fix fasta file names and headers.
#    Check that we have fna faa input FASTA file pairs and that they contain the same number of sequences and
#    instances of each taxon. Mark dir as top_dir
#----------------------------------------------------------------------------------------------------------------

if [ -d get_phylomarkers_run_${dir_suffix}_${TIMESTAMP_SHORT} ]
then
    printf "${RED} >>> Found and older get_phylomarkers_run_${dir_suffix}_${TIMESTAMP_SHORT}/ directory. Please remove or rename and re-run!${NC}\n\n"| \
    tee -a ${logdir}/get_phylomarkers_run_${dir_suffix}_${TIMESTAMP_SHORT}.log
    exit 2
fi 

# make sure we have *.faa and *.fna file pairs to work on
nfna=$(ls *.fna | wc -l) 
if [ "$?" -gt "0" ]
then
   printf "\n${RED} >>> ERROR: there are no input fna files to work on!\n\tPlease check input FASTA files: [you may need to run compare_clusters.pl with -t NUM_OF_INPUT_GENOMES -n]\n\tPlease check the GET_HOMOLOGUES manual${NC}\n${LBLUE}http://eead-csic-compbio.github.io/get_homologues/manual/${NC}\n" | tee -a ${logdir}/get_phylomarkers_run_${dir_suffix}_${TIMESTAMP_SHORT}.log 
   exit 2
fi

nfaa=$(ls *.faa | wc -l) 
if [ "$?" -gt "0" ]
then
   printf "\n${RED} >>> ERROR: there are no input faa files to work on!\n\tPlease check input FASTA files: [you may need to run compare_clusters.pl with -t NUM_OF_INPUT_GENOMES]\n\tPlease check the GET_HOMOLOGUES manual${NC}\n${LBLUE}http://eead-csic-compbio.github.io/get_homologues/manual/${NC}\n" | tee -a ${logdir}/get_phylomarkers_run_${dir_suffix}_${TIMESTAMP_SHORT}.log   
   exit 2 
fi

if [ "$nfna" -ne "$nfaa" ]
then
 printf "\n${RED} >>> ERROR: there are no equal numbers of fna and faa input files to work on!\n\tPlease check input FASTA files: [you may need to run compare_clusters.pl with -t NUM_OF_INPUT_GENOMES -n; and a second time: run compare_clusters.pl with -t NUM_OF_INPUT_GENOMES]\n\tPlease check the GET_HOMOLOGUES manual${NC}\n${LBLUE}http://eead-csic-compbio.github.io/get_homologues/manual/${NC}\n" | tee -a ${logdir}/get_phylomarkers_run_${dir_suffix}_${TIMESTAMP_SHORT}.log 
 exit 3
fi 

mkdir get_phylomarkers_run_${dir_suffix}_${TIMESTAMP_SHORT} && cd get_phylomarkers_run_${dir_suffix}_${TIMESTAMP_SHORT}
top_dir=$(pwd)

print_start_time && printf "${LBLUE}# processing source fastas in directory get_phylomarkers_run_${dir_suffix}_${TIMESTAMP_SHORT} ...${NC}\n"| \
tee -a ${logdir}/get_phylomarkers_run_${dir_suffix}_${TIMESTAMP_SHORT}.log

ln -s ../*.faa .
ln -s ../*.fna .

# fix fasta file names with two and three dots
$distrodir/rename.pl 's/\.\.\./\./g' *.faa
$distrodir/rename.pl 's/\.\.\./\./g' *.fna

# 1.0 check that all fasta files contain the same number of sequences
FASTASIZES=$(grep -c "^>" *.f[na]a | cut -d":" -f 2 | sort | uniq)
NSEQSFASTA=$(grep -c "^>" *.f[na]a | cut -d":" -f 2 | sort | uniq | wc -l)
[ $NSEQSFASTA -gt 1 ] && printf "\n${RED} >>> ERROR: Input FASTA files do not contain the same number of sequences...${NC}\n$FASTASIZES\n" && exit 4

# 1.1 fix fastaheaders of the source protein and DNA fasta files
for file in *faa; do awk 'BEGIN {FS = "|"}{print $1, $2, $3}' $file|perl -pe 'if(/^>/){s/>\S+/>/; s/>\h+/>/; s/\h+/_/g; s/,//g; s/;//g; s/://g; s/\(//g; s/\)//g}' > ${file}ed; done
for file in *fna; do awk 'BEGIN {FS = "|"}{print $1, $2, $3}' $file|perl -pe 'if(/^>/){s/>\S+/>/; s/>\h+/>/; s/\h+/_/g; s/,//g; s/;//g; s/://g; s/\(//g; s/\)//g}' > ${file}ed; done

print_start_time && printf "${BLUE}# Performing strain composition check on f?aed files ...${NC}\n"
faaed_strain_intersection_check=$(grep '>' *faaed | cut -d: -f2 | sort | uniq -c | awk '{print $1}' | sort | uniq -c | wc -l)
fnaed_strain_intersection_check=$(grep '>' *fnaed | cut -d: -f2 | sort | uniq -c | awk '{print $1}' | sort | uniq -c | wc -l)

# check that each file has the same number of strains and a single instance for each strain
if [ $faaed_strain_intersection_check -eq 1 -a $fnaed_strain_intersection_check -eq 1 ]
then 
   printf "${GREEN} >>> Strain check OK: each f?aed file has the same number of strains and a single instance for each strain${NC}\n" | \
   tee -a ${logdir}/get_phylomarkers_run_${dir_suffix}_${TIMESTAMP_SHORT}.log
else 
     printf "\n${RED} >>> ERROR: Input f?aed files do not contain the same number of strains and a single instance for each strain...\n\tPlease check input FASTA files: [you may need to run compare_clusters.pl with -t NUM_OF_INPUT_GENOMES]\n\tPlease check the GET_HOMOLOGUES manual${NC}\n${LBLUE}http://eead-csic-compbio.github.io/get_homologues/manual/${NC}\n" | \
     tee -a ${logdir}/get_phylomarkers_run_${dir_suffix}_${TIMESTAMP_SHORT}.log 
     exit 5
fi


# 1.2 add_nos2fasta_header.pl to avoid problems with duplicate labels
[ $DEBUG -eq 1 -o $VERBOSITY -eq 1 ] && echo " > ${distrodir}/run_parallel_cmmds.pl faaed 'add_nos2fasta_header.pl $file > ${file}no' $n_cores &> /dev/null" | \
tee -a ${logdir}/get_phylomarkers_run_${dir_suffix}_${TIMESTAMP_SHORT}.log
${distrodir}/run_parallel_cmmds.pl faaed 'add_nos2fasta_header.pl $file > ${file}no' $n_cores &> /dev/null

[ $DEBUG -eq 1 -o $VERBOSITY -eq 1 ] && echo " > ${distrodir}/run_parallel_cmmds.pl fnaed 'add_nos2fasta_header.pl $file > ${file}no' $n_cores &> /dev/null" | \
tee -a ${logdir}/get_phylomarkers_run_${dir_suffix}_${TIMESTAMP_SHORT}.log
${distrodir}/run_parallel_cmmds.pl fnaed 'add_nos2fasta_header.pl $file > ${file}no' $n_cores &> /dev/null

no_alns=$(ls *.fnaedno | wc -l)

[ $no_alns -eq 0 ] && printf "\n${RED} >>> ERROR: There are no codon alignments to work on! Something went wrong. Please check input and settings ...${NC}\n\n" | \
 tee -a ${logdir}/get_phylomarkers_run_${dir_suffix}_${TIMESTAMP_SHORT}.log && exit 4
print_start_time && printf "${BLUE}# Total number of alignments to be computed $no_alns ${NC}\n" | \
 tee -a ${logdir}/get_phylomarkers_run_${dir_suffix}_${TIMESTAMP_SHORT}.log

# 1.3 generate a tree_labels.list file for later tree labeling
print_start_time && printf "${BLUE}# generating the labels file for tree-labeling ...${NC}\n" | \
 tee -a ${logdir}/get_phylomarkers_run_${dir_suffix}_${TIMESTAMP_SHORT}.log
tree_labels_dir=$(pwd)
grep '>' $(ls *fnaedno|head -1) > tree_labels.list

[ $DEBUG -eq 1 -o $VERBOSITY -eq 1 ] && echo " > perl -pe '$c++; s/>/$c\t/; s/\h\[/_[/' tree_labels.list > ed && mv ed tree_labels.list" | \
tee -a ${logdir}/get_phylomarkers_run_${dir_suffix}_${TIMESTAMP_SHORT}.log
perl -pe '$c++; s/>/$c\t/; s/\h\[/_[/' tree_labels.list > ed && mv ed tree_labels.list

#------------------------------------------------------------------------------------------------------#
# >>>BLOCK 2. Generate cdnAlns with with pal2nal, maintaining the input-order of the source fastas <<< #
#------------------------------------------------------------------------------------------------------#

# 2.1 generate the protein alignments using clustalo
print_start_time &&  printf "${BLUE}# generating $no_alns codon alignments ...${NC}\n"|tee -a ${logdir}/get_phylomarkers_run_${dir_suffix}_${TIMESTAMP_SHORT}.log
[ $DEBUG -eq 1 -o $VERBOSITY -eq 1 ] && echo " > '${distrodir}/run_parallel_cmmds.pl faaedno clustalo -i $file -o ${file%.*}_cluo.faaln --output-order input-order' $n_cores &> /dev/null" | \
tee -a ${logdir}/get_phylomarkers_run_${dir_suffix}_${TIMESTAMP_SHORT}.log
${distrodir}/run_parallel_cmmds.pl faaedno 'clustalo -i $file -o ${file%.*}_cluo.faaln --output-order input-order' $n_cores &> clustalo.log

if grep -q "Thread creation failed" clustalo.log; then
	printf "\n${RED} >>> ERROR: This system cannot launch too many threads, please use option -n and re-run ...${NC}\n" | \
		tee -a ${logdir}/get_phylomarkers_run_${dir_suffix}_${TIMESTAMP_SHORT}.log && exit 4
fi

# 2.2 generate the codon alignments (files with *_cdnAln.fasta extension) using pal2nal.pl, 
#     excluding gaps, and mismatched codons, assuming a bacterial genetic code
# NOTE: to execute run_parallel_cmmds.pl with a customized command, resulting from the interpolation of multiple varialbles,
#       we have to first construct the command line in a variable and pipe its content into bash for execution
faaln_ext=faaln
command="${distrodir}/run_parallel_cmmds.pl $faaln_ext '${distrodir}/pal2nal.pl \$file \${file%_cluo.faaln}.fnaedno -output fasta -nogap -nomismatch -codontable $codontable > \${file%_cluo.faaln}_cdnAln.fasta' $n_cores"

# now we can execute run_parallel_cmmds.pl with a customized command, resulting from the interpolation of multiple varialbles
[ $DEBUG -eq 1 -o $VERBOSITY -eq 1 ] && echo " > $command | bash &> /dev/null" | \
tee -a ${logdir}/get_phylomarkers_run_${dir_suffix}_${TIMESTAMP_SHORT}.log
echo "$command" | bash &> /dev/null

# check we got non-empty *cdnAln.fasta files
for f in *cdnAln.fasta
do
     if [ ! -s $f ] 
     then
           printf "\n${LRED} >>> Warning: produced empty codon alignment $f!\n     ... Will skip this locus and move it to problematic_alignments/ ...\n\n${NC}" | \
	   tee -a ${logdir}/get_phylomarkers_run_${dir_suffix}_${TIMESTAMP_SHORT}.log
	   [ ! -d problematic_alignments ] && mkdir problematic_alignments 
	   locus_base=${f%_cdnAln.fasta}
	   mv $f problematic_alignments
	   mv ${locus_base}* problematic_alignments
     fi	   
done

# 2.3 cleanup: remove the source faa, fna, fnaed and faaed files; make numbered_fna_files.tgz and numbered_faa_files.tgz; rm *aedno
rm *fnaed *faaed *faa *fna
[ "$DEBUG" -eq "1" -o "$VERBOSITY" -eq "1" ] && echo " > tar -czf numbered_fna_files.tgz *fnaedno" | tee -a ${logdir}/get_phylomarkers_run_${dir_suffix}_${TIMESTAMP_SHORT}.log
tar -czf numbered_fna_files.tgz *fnaedno
[ "$DEBUG" -eq "1" -o "$VERBOSITY" -eq "1" ] && echo " > tar -czf numbered_fna_files.tgz *faaedno" | tee -a ${logdir}/get_phylomarkers_run_${dir_suffix}_${TIMESTAMP_SHORT}.log
tar -czf numbered_faa_files.tgz *faaedno
rm ./*aedno

#----------------------------------------------------------------------------------------------------------#
#>>> BLOCK 3. run Phi-test to identify recombinant codon alignments on all *_cdnAln.fasta source files <<< #
#----------------------------------------------------------------------------------------------------------#
# 3.1 make a new PhiPack subdirectory to work in. generate symlinks to ../*fasta files
#     Mark dir as phipack_dir
mkdir PhiPack && cd PhiPack
phipack_dir=$(pwd)
ln -s ../*fasta .

# 3.1.2 check that we have codon alignments before proceeding
no_fasta_files=$(ls *.fasta | wc -l)

[ $no_fasta_files -lt 1 ] && print_start_time && printf "\n${RED} >>> ERROR: there are no codon alignments to run Phi on. Will exit now!${NC}\n\n" | \
tee -a ${logdir}/get_phylomarkers_run_${dir_suffix}_${TIMESTAMP_SHORT}.log && exit 3

# 3.2 run Phi from the PhiPack in parallel
print_start_time && printf "${LBLUE}# running Phi recombination test in PhiPack dir ...${NC}\n" | tee -a ${logdir}/get_phylomarkers_run_${dir_suffix}_${TIMESTAMP_SHORT}.log
[ $DEBUG -eq 1 -o $VERBOSITY -eq 1 ] && echo " > ${distrodir}/run_parallel_cmmds.pl fasta 'Phi -f $file -p 1000 > ${file%.*}_Phi.log' $n_cores &> /dev/null"
${distrodir}/run_parallel_cmmds.pl fasta 'Phi -f $file -p 1000 > ${file%.*}_Phi.log' $n_cores &> /dev/null

# 3.3 process the *_Phi.log files generated by Phi to write a summary table and print short overview to STDOUT
declare -a nonInfoAln
COUNTERNOINFO=0

[ -s Phi.log ] && rm Phi.log

for f in *_Phi.log
do 
   grep '^Too few' $f &> /dev/null # if exit status is not 0
   if [ $? -eq 0 ]
   then 
       # if there are "Too few informative sites"; assign dummy, non significative values
       # and report this in the logfile!
       let COUNTERNOINFO++
       norm=1
       perm=1
       echo -e "$f\t$norm\t$perm\tTOO_FEW_INFORMATIVE_SITES"
       nonInfoAln[$COUNTERNOINFO]="$f"
   else
      perm=$(grep Permut $f | awk '{print $3}')
      norm=$(grep Normal $f | awk '{print $3}')
      
      [ "$perm" == "--" ] && perm=1
      [ "$norm" == "--" ] && norm=1
      #nonInfoAln[$COUNTERNOINFO]=""
    
      echo -e "$f\t$norm\t$perm"
   fi
done > Phi_results_${TIMESTAMP_SHORT}.tsv

check_output Phi_results_${TIMESTAMP_SHORT}.tsv $parent_PID | tee -a ${logdir}/get_phylomarkers_run_${dir_suffix}_${TIMESTAMP_SHORT}.log

no_non_recomb_alns_perm_test=$(awk '$2 > 5e-02 && $3 > 5e-02' Phi_results_${TIMESTAMP_SHORT}.tsv | wc -l)
total_no_cdn_alns=$(ls *_cdnAln.fasta | wc -l)

# 3.4 Check the number of non-recombinant alignments remaining
if [ ${#nonInfoAln[@]} == 0 ]
then 
  print_start_time && printf "${GREEN} >>> Phi test result: there are $no_non_recomb_alns_perm_test non-recombinant alignments out of $total_no_cdn_alns input alignments${NC}\n" | \
  tee -a ${logdir}/get_phylomarkers_run_${dir_suffix}_${TIMESTAMP_SHORT}.log
fi

if [ ${#nonInfoAln[@]} -gt 0 ]
then
  print_start_time && printf "${LRED} >>> Phi test WARNING: there ${#nonInfoAln[@]} alignments with too few informative sites for the Phi test to work on ... ${NC} \n"| \
  tee -a ${logdir}/get_phylomarkers_run_${dir_suffix}_${TIMESTAMP_SHORT}.log

  # print the names of the alignments with too few informative sites for the Phi test to be work on
  if [ "$VERBOSITY" -eq "1" -o  "$DEBUG" -eq "1" ]
  then
      printf "${LRED} >>> The alignments with too few informative sites for the Phi test to be work on are:${NC} \n"| \
      tee -a ${logdir}/get_phylomarkers_run_${dir_suffix}_${TIMESTAMP_SHORT}.log
      for f in "${nonInfoAln[@]}"
      do
           printf "${LRED} >>> ${f}${NC}\n"| tee -a ${logdir}/get_phylomarkers_run_${dir_suffix}_${TIMESTAMP_SHORT}.log
      done
  fi
fi

[ $no_non_recomb_alns_perm_test -lt 1 ] && print_start_time && printf "\n${LRED} >>> Warning: All alignments seem to have recombinant sequences. will exit now!${NC}\n\n" | \
tee -a ${logdir}/get_phylomarkers_run_${dir_suffix}_${TIMESTAMP_SHORT}.log && exit 3

#3.5 cleanup dir
tar -czf Phi_test_log_files.tgz *Phi.log
[ -s Phi_test_log_files.tgz ] && rm *Phi.log Phi.inf*

# 3.5.1 mv non-recombinant codon alignments and protein alignments to their own directories:
#     non_recomb_cdn_alns/ and  non_recomb_cdn_alns/
#     Mark dir as non_recomb_cdn_alns
mkdir non_recomb_cdn_alns

for base in $(awk '$2 > 5e-02 && $3 > 5e-02' Phi_results_${TIMESTAMP_SHORT}.tsv | awk '{print $1}'|sed 's/_Phi\.log//')
do 
   cp ${base}.fasta non_recomb_cdn_alns
done 

mkdir non_recomb_FAA_alns
for base in $(awk '$2 > 5e-02 && $3 > 5e-02' Phi_results_${TIMESTAMP_SHORT}.tsv | awk '{print $1}'|sed 's/_cdnAln_Phi\.log//')
do
   cp ../${base}*.faaln non_recomb_FAA_alns
done

# 3.5.2 cleanup phipack_dir; rm *fasta, which are the same as in topdir
rm *cdnAln.fasta

#------------------------------------------------------------------------------------------------
#>>>BLOCK 4. Compute individual ML gene trees from codon or protein alignments 
#            and run the kdetrees test to identify deviant trees.
#            Keep the ones passing the test (non-significant tests) for supermatrix analysis.
#------------------------------------------------------------------------------------------------
#
# 4.1 cd into non_recomb_cdn_alns and run FastTree in parallel on all codon alignments, 
#     using exhaustive tree-rearrangements: turn off NNI heuristics, and to always optimize 
#     all 5 branches at each NNI, optimizing all 5 branches in as many as 3 rounds. 
#     This is hardcoded because we want optimal gene trees at this stage,
#     which are rapidly computed by FastTree due to the low number of columns in CDS/product
#     alignments. 
#     NOTES: 
#         1. at this point the code is divided into two large if blocks to run -t DNA|PROT
#                    [ "$mol_type" == "DNA" ]
#                    [ "$mol_type" == "PROT" ]
#         2. blocks 4 and 5 are highly repetitive; refactor into a subroutine, passing proper suffixes
if [ "$mol_type" == "DNA" ]
then
    cd non_recomb_cdn_alns
    non_recomb_cdn_alns_dir=$(pwd)

    print_start_time && printf "${LBLUE}# working in dir non_recomb_cdn_alns ...${NC}\n" | \
    tee -a ${logdir}/get_phylomarkers_run_${dir_suffix}_${TIMESTAMP_SHORT}.log

    print_start_time && printf "${BLUE}# estimating $no_non_recomb_alns_perm_test gene trees from non-recombinant sequences ...${NC}\n"| \
    tee -a ${logdir}/get_phylomarkers_run_${dir_suffix}_${TIMESTAMP_SHORT}.log
    [ $DEBUG -eq 1 -o $VERBOSITY -eq 1 ] && echo " > ${distrodir}/run_parallel_cmmds.pl fasta 'FastTree -quiet -nt -gtr -gamma -bionj -slownni -mlacc 3 -spr 8 -sprlength 8 < $file > ${file%.*}_allFTGTRG.ph' $n_cores &> /dev/null"  | \
    tee -a ${logdir}/get_phylomarkers_run_${dir_suffix}_${TIMESTAMP_SHORT}.log   
    ${distrodir}/run_parallel_cmmds.pl fasta 'FastTree -quiet -nt -gtr -gamma -bionj -slownni -mlacc 3 -spr 8 -sprlength 8 < $file > ${file%.*}_allFTGTRG.ph' $n_cores &> /dev/null
    
    no_gene_trees=$(ls *allFTGTRG.ph|wc -l)
    [ $no_gene_trees -lt 1 ] && print_start_time && printf "\n${LRED} >>> Warning: There are no gene tree to work on in non_recomb_cdn_alns/. will exit now!${NC}\n\n" | \
    tee -a ${logdir}/get_phylomarkers_run_${dir_suffix}_${TIMESTAMP_SHORT}.log && exit 3
    
    #remove trees with < 5 branches
    print_start_time && printf "${BLUE}# counting branches on $no_non_recomb_alns_perm_test gene trees ...${NC}\n" | \
    tee -a ${logdir}/get_phylomarkers_run_${dir_suffix}_${TIMESTAMP_SHORT}.log
    [ $DEBUG -eq 1 -o $VERBOSITY -eq 1 ] && echo " > count_tree_branches ph no_tree_branches.list &> /dev/null" | \
    tee -a ${logdir}/get_phylomarkers_run_${dir_suffix}_${TIMESTAMP_SHORT}.log

    count_tree_branches ph no_tree_branches.list # &> /dev/null
    [ ! -s no_tree_branches.list ] && install_Rlibs_msg no_tree_branches.list ape

    check_output no_tree_branches.list $parent_PID | tee -a ${logdir}/get_phylomarkers_run_${dir_suffix}_${TIMESTAMP_SHORT}.log
    
    # remove trees with < 5 external branches (leaves)
     [ $DEBUG -eq 1 -o $VERBOSITY -eq 1 ] && echo " >  removing trees with < 5 external branches (leaves)" | \
     tee -a ${logdir}/get_phylomarkers_run_${dir_suffix}_${TIMESTAMP_SHORT}.log
    
    no_tree_counter=0
    for phy in $(grep -v '^#Tree' no_tree_branches.list | awk -v min_no_ext_branches=$min_no_ext_branches 'BEGIN{FS="\t"; OFS="\t"}$7 < min_no_ext_branches' |cut -f1)
    do
         base=$(echo $phy | sed 's/_allFT.*ph//')
	 print_start_time && printf "${LRED} >>> will remove ${base}* because it has < 5 branches!${NC}\n" | \
	 tee -a ${logdir}/get_phylomarkers_run_${dir_suffix}_${TIMESTAMP_SHORT}.log 
	 rm ${base}*
	 let no_tree_counter++
    done

    printf "${LRED} >>> WARNING: there are $no_tree_counter trees with < 1 internal branches (no real trees) that will be discarded ...${NC}\n" | tee -a ${logdir}/get_phylomarkers_run_${dir_suffix}_${TIMESTAMP_SHORT}.log
    
    # 4.1 generate the all_GTRG_trees.tre holding all source trees, which is required by kdetrees
    #     Make a check for the existence of the file to interrupt the pipeline if something has gone wrong
    [ $DEBUG -eq 1 -o $VERBOSITY -eq 1 ] && echo " > cat *allFTGTRG.ph > all_GTRG_trees.tre" | \
    tee -a ${logdir}/get_phylomarkers_run_${dir_suffix}_${TIMESTAMP_SHORT}.log
    cat *allFTGTRG.ph > all_GTRG_trees.tre
    check_output all_GTRG_trees.tre $parent_PID | tee -a ${logdir}/get_phylomarkers_run_${dir_suffix}_${TIMESTAMP_SHORT}.log
    [ ! -s all_GTRG_trees.tre ] && exit 3

    # 4.2 run_kdetrees.R at desired stringency 
    print_start_time && printf "${BLUE}# running kde test ...${NC}\n" | tee -a ${logdir}/get_phylomarkers_run_${dir_suffix}_${TIMESTAMP_SHORT}.log
    [ $DEBUG -eq 1 -o $VERBOSITY -eq 1 ] && echo " > run_kdetrees.R ph all_GTRG_trees.tre $kde_stringency &> /dev/null" | \
    tee -a ${logdir}/get_phylomarkers_run_${dir_suffix}_${TIMESTAMP_SHORT}.log
    run_kdetrees.R ph all_GTRG_trees.tre $kde_stringency &> /dev/null
    [ ! -s kde_dfr_file_all_GTRG_trees.tre.tab ] && install_Rlibs_msg kde_dfr_file_all_GTRG_trees.tre.tab kdetrees,ape
    check_output kde_dfr_file_all_GTRG_trees.tre.tab $parent_PID | tee -a ${logdir}/get_phylomarkers_run_${dir_suffix}_${TIMESTAMP_SHORT}.log

    # 4.3 mv outliers to kde_outliers
    no_kde_outliers=$(grep -c outlier kde_dfr_file_all_GTRG_trees.tre.tab)
    no_kde_ok=$(grep -v outlier kde_dfr_file_all_GTRG_trees.tre.tab|grep -vc '^file')

    # 4.4 Check how many cdnAlns passed the test and separate into two subirectories those passing and failing the test
    if [ $no_kde_outliers -gt 0 ]
    then
        print_start_time && printf "${BLUE}# making dir kde_outliers/ and moving $no_kde_outliers outlier files into it ...${NC}\n" | \
	tee -a ${logdir}/get_phylomarkers_run_${dir_suffix}_${TIMESTAMP_SHORT}.log
        mkdir kde_outliers
        for f in $(grep outlier kde_dfr_file_all_GTRG_trees.tre.tab|cut -f1|sed 's/_allFTGTRG.ph//'); do mv ${f}* kde_outliers; done
    else
        print_start_time && printf "${GREEN} >>> there are no kde-test outliers ...${NC}\n" | \
	tee -a ${logdir}/get_phylomarkers_run_${dir_suffix}_${TIMESTAMP_SHORT}.log
    fi

    if [ $no_kde_ok -gt 0 ]
    then
        print_start_time && printf "${BLUE}# making dir kde_ok/ and linking $no_kde_ok selected files into it ...${NC}\n" | \
	tee -a ${logdir}/get_phylomarkers_run_${dir_suffix}_${TIMESTAMP_SHORT}.log
        mkdir kde_ok
        cd kde_ok
        ln -s ../*ph .
   
        print_start_time && printf "${BLUE}# labeling $no_kde_ok gene trees in dir kde_ok/ ...${NC}\n" | \
	tee -a ${logdir}/get_phylomarkers_run_${dir_suffix}_${TIMESTAMP_SHORT}.log
        [ $DEBUG -eq 1 -o $VERBOSITY -eq 1 ] && echo " > ${distrodir}/run_parallel_cmmds.pl ph 'add_labels2tree.pl ../../../tree_labels.list $file' $n_cores &> /dev/null" | \
	tee -a ${logdir}/get_phylomarkers_run_${dir_suffix}_${TIMESTAMP_SHORT}.log
        ${distrodir}/run_parallel_cmmds.pl ph 'add_labels2tree.pl ../../../tree_labels.list $file' $n_cores &> /dev/null
        
	# remove symbolic links to cleanup kde_ok/
        for f in $(ls *ph|grep -v '_ed\.ph'); do rm $f; done
    
        cd ..
    else
        print_start_time &&  printf "${RED}# There are $no_kde_ok gene trees producing non-significant kde-test results! Increase the actual -k $kde_stringency value. Will stop here. ...${NC}\n" | \
	tee -a ${logdir}/get_phylomarkers_run_${dir_suffix}_${TIMESTAMP_SHORT}.log
	 exit 5
    fi

#===============================================#
# >>> ---------- 5. RUNMODES -------------- <<< #
#===============================================#
# >>> RUNMODE 1: PHYLOGENETICS <<< #
#----------------------------------#    
    if [ $runmode -eq 1 ] 
    then
        # 5.1 compute average bipartition support values for each gene tree
        wkdir=$(pwd) 
        
        print_start_time && printf "${BLUE}# computing tree support values ...${NC}\n" | tee -a ${logdir}/get_phylomarkers_run_${dir_suffix}_${TIMESTAMP_SHORT}.log
        [ $DEBUG -eq 1 -o $VERBOSITY -eq 1 ] && echo " > compute_suppValStas_and_RF-dist.R $wkdir 1 fasta ph 1 &> /dev/null" | \
	tee -a ${logdir}/get_phylomarkers_run_${dir_suffix}_${TIMESTAMP_SHORT}.log
        compute_suppValStas_and_RF-dist.R $wkdir 1 fasta ph 1 &> /dev/null
	
        print_start_time && printf "${BLUE}# writing summary tables ...${NC}\n" | tee -a ${logdir}/get_phylomarkers_run_${dir_suffix}_${TIMESTAMP_SHORT}.log
        min_supp_val_perc=${min_supp_val#0.}
        no_digits=${#min_supp_val_perc}
        [ $no_digits -eq 1 ] && min_supp_val_perc=${min_supp_val_perc}0
    
        awk -v min_supp_val=$min_supp_val '$2 >= min_supp_val' sorted_aggregated_support_values4loci.tab > sorted_aggregated_support_values4loci_ge${min_supp_val_perc}perc.tab
        check_output sorted_aggregated_support_values4loci_ge${min_supp_val_perc}perc.tab $parent_PID | \
	tee -a ${logdir}/get_phylomarkers_run_${dir_suffix}_${TIMESTAMP_SHORT}.log
    
        no_top_markers=$(perl -lne 'END{print $.}' sorted_aggregated_support_values4loci_ge${min_supp_val_perc}perc.tab)
        top_markers_dir="top_${no_top_markers}_markers_ge${min_supp_val_perc}perc"
        top_markers_tab=$(ls sorted_aggregated_support_values4loci_ge${min_supp_val_perc}perc.tab)
    
        print_start_time && printf "${LBLUE}# making dir $top_markers_dir and moving $no_top_markers top markers into it ...${NC}\n" | \
	tee -a ${logdir}/get_phylomarkers_run_${dir_suffix}_${TIMESTAMP_SHORT}.log
        mkdir $top_markers_dir && cd $top_markers_dir
        top_markers_dir=$(pwd)
        ln -s ../$top_markers_tab .
        for base in $(awk '{print $1}' $top_markers_tab|grep -v loci|sed 's/"//g'); do ln -s ../${base}* .; done
   
        [ $no_top_markers -lt 2 ] && print_start_time && printf "\n${LRED} >>> Warning: There are less than 2 top markers. Relax your filtering thresholds. will exit now!${NC}\n\n" | tee -a ${logdir}/get_phylomarkers_run_${dir_suffix}_${TIMESTAMP_SHORT}.log && exit 3


        # 5.2 generate supermatrix (concatenated alignment) 
        print_start_time && printf "${BLUE}# concatenating $no_top_markers top markers into supermatrix ...${NC}\n" | \
	tee -a ${logdir}/get_phylomarkers_run_${dir_suffix}_${TIMESTAMP_SHORT}.log
        [ $DEBUG -eq 1 -o $VERBOSITY -eq 1 ] && echo " > concat_alns fasta $parent_PID &> /dev/null" | \
	tee -a ${logdir}/get_phylomarkers_run_${dir_suffix}_${TIMESTAMP_SHORT}.log
        concat_alns fasta $parent_PID &> /dev/null
	
        # 5.3 remove uninformative sites from the concatenated alignment to speed up computation
        print_start_time && printf "${BLUE}# removing uninformative sites from concatenated alignment ...${NC}\n"| \
	tee -a ${logdir}/get_phylomarkers_run_${dir_suffix}_${TIMESTAMP_SHORT}.log
        [ $DEBUG -eq 1 -o $VERBOSITY -eq 1 ] && echo " > ${distrodir}/remove_uninformative_sites_from_aln.pl < concat_cdnAlns.fna > concat_cdnAlns.fnainf" | \
	tee -a ${logdir}/get_phylomarkers_run_${dir_suffix}_${TIMESTAMP_SHORT}.log
        #$distrodir/remove_uninformative_sites_from_aln.pl < concat_cdnAlns.fna > concat_cdnAlns.fnainf
	${distrodir}/remove_uninformative_sites_from_aln.pl < concat_cdnAlns.fna > concat_cdnAlns.fnainf
        check_output concat_cdnAlns.fnainf $parent_PID | tee -a ${logdir}/get_phylomarkers_run_${dir_suffix}_${TIMESTAMP_SHORT}.log
	
	[ ! -s concat_cdnAlns.fnainf ] && print_start_time && printf "\n${RED} >>> ERROR: The expected file concat_cdnAlns.fnainf was not produced! will exit now!${NC}\n\n" | \
	tee -a ${logdir}/get_phylomarkers_run_${dir_suffix}_${TIMESTAMP_SHORT}.log && exit 3

        
	#if [ "$search_algorithm" == "F" ] # always run FastTree, ok
	#then
	  # 5.4 run FasTree under the GTR+G model 
          print_start_time && printf "${BLUE}# running FastTree on the concatenated alignment with $search_thoroughness thoroughness. This may take a while ...${NC}\n" | \
	  tee -a ${logdir}/get_phylomarkers_run_${dir_suffix}_${TIMESTAMP_SHORT}.log

          if [ "$search_thoroughness" == "high" ]
          then
            FastTree -quiet -nt -gtr -bionj -slownni -gamma -mlacc 3 -spr $spr -sprlength $spr_length -log ${tree_prefix}_nonRecomb_KdeFilt_cdnAlns_FTGTRG.log < concat_cdnAlns.fnainf > ${tree_prefix}_nonRecomb_KdeFilt_cdnAlns_FTGTRG.ph
          fi

          if [ "$search_thoroughness" == "medium" ]
          then
            FastTree -quiet -nt -gtr -bionj -slownni -gamma -mlacc 2 -spr $spr -sprlength $spr_length -log ${tree_prefix}_nonRecomb_KdeFilt_cdnAlns_FTGTRG.log < concat_cdnAlns.fnainf > ${tree_prefix}_nonRecomb_KdeFilt_cdnAlns_FTGTRG.ph
          fi

          if [ "$search_thoroughness" == "low" ]
          then
            FastTree -quiet -nt -gtr -bionj -gamma -spr $spr -sprlength $spr_length -log ${tree_prefix}_nonRecomb_KdeFilt_cdnAlns_FTGTRG.log < concat_cdnAlns.fnainf > ${tree_prefix}_nonRecomb_KdeFilt_cdnAlns_FTGTRG.ph
          fi

          if [ "$search_thoroughness" == "lowest" ]
          then
            FastTree -quiet -nt -gtr -gamma -mlnni 4 -log ${tree_prefix}_nonRecomb_KdeFilt_cdnAlns_FTGTRG.log < concat_cdnAlns.fnainf > ${tree_prefix}_nonRecomb_KdeFilt_cdnAlns_FTGTRG.ph
          fi

          check_output ${tree_prefix}_nonRecomb_KdeFilt_cdnAlns_FTGTRG.ph $parent_PID | tee -a ${logdir}/get_phylomarkers_run_${dir_suffix}_${TIMESTAMP_SHORT}.log
	
	  if [ -s "${tree_prefix}_nonRecomb_KdeFilt_cdnAlns_FTGTRG.log" -a -s "${tree_prefix}_nonRecomb_KdeFilt_cdnAlns_FTGTRG.ph" ]
	  then
	    #lnL=$(grep ML_Lengths2 "${tree_prefix}_nonRecomb_KdeFilt_cdnAlns_FTGTRG.log" | grep TreeLogLk | sed 's/TreeLogLk[[:space:]]ML_Lengths2[[:space:]]//')
	    lnL=$(grep '^Gamma20LogLk' "${tree_prefix}_nonRecomb_KdeFilt_cdnAlns_FTGTRG.log" |awk '{print $2}'
	    printf "${GREEN} >>> lnL for ${tree_prefix}_nonRecomb_KdeFilt_cdnAlns_FTGTRG.ph = $lnL ${NC}\n" | \
	    tee -a ${logdir}/get_phylomarkers_run_${dir_suffix}_${TIMESTAMP_SHORT}.log
	  else
	    printf "${LRED} >>> WARNING: ${tree_prefix}_nonRecomb_KdeFilt_cdnAlns_FTGTRG.log could not be produced!${NC}\n" | \
	    tee -a ${logdir}/get_phylomarkers_run_${dir_suffix}_${TIMESTAMP_SHORT}.log
	  fi
    
          print_start_time && printf "${BLUE}# Adding labels back to tree ...${NC}\n" | tee -a ${logdir}/get_phylomarkers_run_${dir_suffix}_${TIMESTAMP_SHORT}.log
          [ $DEBUG -eq 1 -o $VERBOSITY -eq 1 ] && echo " > ${distrodir}/add_labels2tree.pl ${tree_labels_dir}/tree_labels.list ${tree_prefix}_nonRecomb_KdeFilt_cdnAlns_FTGTRG.ph &> /dev/null"
          ${distrodir}/add_labels2tree.pl ${tree_labels_dir}/tree_labels.list ${tree_prefix}_nonRecomb_KdeFilt_cdnAlns_FTGTRG.ph &> /dev/null

          if [ -s ${tree_prefix}_nonRecomb_KdeFilt_cdnAlns_FTGTRG_ed.ph ]
          then
            mv ${tree_prefix}_nonRecomb_KdeFilt_cdnAlns_FTGTRG_ed.ph ${tree_prefix}_nonRecomb_KdeFilt_${no_top_markers}cdnAlns_FTGTRG_ed.sptree # for compute_suppValStats_and_RF-dist.R
            check_output ${tree_prefix}_nonRecomb_KdeFilt_${no_top_markers}cdnAlns_FTGTRG_ed.sptree $parent_PID | \
	    tee -a ${logdir}/get_phylomarkers_run_${dir_suffix}_${TIMESTAMP_SHORT}.log
            printf "${GREEN} >>> found in dir $top_markers_dir ...${NC}\n\n" | tee -a ${logdir}/get_phylomarkers_run_${dir_suffix}_${TIMESTAMP_SHORT}.log
          else
              printf "${LRED} >>> WARNING: ${tree_prefix}_nonRecomb_KdeFilt_${no_top_markers}cdnAlns_FTGTRG_ed.sptree could not be produced!${NC}"
          fi 
	
          print_start_time && printf "${BLUE}# computing the mean support values and RF-distances of each gene tree to the concatenated tree   ...${NC}\n" | \
	  tee -a ${logdir}/get_phylomarkers_run_${dir_suffix}_${TIMESTAMP_SHORT}.log
	  [ $DEBUG -eq 1 -o $VERBOSITY -eq 1 ] && echo " > compute_suppValStas_and_RF-dist.R $top_markers_dir 2 fasta ph 1 &> /dev/null" | \
	  tee -a ${logdir}/get_phylomarkers_run_${dir_suffix}_${TIMESTAMP_SHORT}.log
          compute_suppValStas_and_RF-dist.R $top_markers_dir 2 fasta ph 1 &> /dev/null
	
	  # top100_median_support_values4loci.tab should probably not be written in the first instance
	  [ -s top100_median_support_values4loci.tab -a "${no_top_markers}" -lt 101 ] && rm top100_median_support_values4loci.tab
    
       #fi # [ "$search_algorithm" == "F" ]
       
       if [ "$search_algorithm" == "I" ]
       then
	  # 5.5 run IQ-tree in addition to FastTree, if requested
          echo
	  printf "${YELLOW} >>>>>>>>>>>>>>> ModelFinder + IQ-TREE run <<<<<<<<<<<<<<< ${NC}\n" | \
	  tee -a ${logdir}/get_phylomarkers_run_${dir_suffix}_${TIMESTAMP_SHORT}.log
	  echo 

	  print_start_time && printf "${BLUE}# running ModelFinder on the concatenated alignment with $IQT_models. This will take a while ...${NC}\n" | \
	  tee -a ${logdir}/get_phylomarkers_run_${dir_suffix}_${TIMESTAMP_SHORT}.log
           
          iqtree-omp -s concat_cdnAlns.fnainf -st DNA -mset "$IQT_models" -m MF -nt AUTO &> /dev/null 
	  
	  check_output concat_cdnAlns.fnainf.log $parent_PID | tee -a ${logdir}/get_phylomarkers_run_${dir_suffix}_${TIMESTAMP_SHORT}.log
	  
	  best_model=$(grep '^Best-fit model' concat_cdnAlns.fnainf.log | cut -d' ' -f 3)
	  printf "${GREEN} >>> Best-fit model: ${best_model} ...${NC}\n" | \
	  tee -a ${logdir}/get_phylomarkers_run_${dir_suffix}_${TIMESTAMP_SHORT}.log
	  
	  mkdir iqtree_abayes && cd iqtree_abayes
	  ln -s ../concat_cdnAlns.fnainf .
	  
	  if [ "$search_thoroughness" == "high" ]
	  then
	     print_start_time && printf "${BLUE}# running IQ-TREE on the concatenated alignment with best model ${best_model} -abayes -bb 1000, starting from $nrep random trees!. This will take a while ...${NC}\n" | \
	     tee -a ${logdir}/get_phylomarkers_run_${dir_suffix}_${TIMESTAMP_SHORT}.log
	     
	     # run nrep IQ-TREE searches under the best-fit model found
	     for ((rep=1;rep<=$nrep;rep++))
	     do 
	         print_start_time && printf "${LBLUE} > iqtree-omp -s concat_cdnAlns.fnainf -st DNA -m "$best_model" -abayes -bb 1000 -nt AUTO -pre abayes_run${rep} &> /dev/null${NC}\n" | \
	         tee -a ${logdir}/get_phylomarkers_run_${dir_suffix}_${TIMESTAMP_SHORT}.log
	        
		 iqtree-omp -s concat_cdnAlns.fnainf -st DNA -m "$best_model" -abayes -bb 1000 -nt AUTO -pre abayes_run${rep} &> /dev/null 
	     done
             
	     grep '^BEST SCORE' *log | sort -nrk5 > sorted_IQ-TREE_searches.out
	     
	     check_output sorted_IQ-TREE_searches.out $parent_PID | tee -a ${logdir}/get_phylomarkers_run_${dir_suffix}_${TIMESTAMP_SHORT}.log
	     best_search=$(head -1 sorted_IQ-TREE_searches.out)
	     best_search_base_name=$(head -1 sorted_IQ-TREE_searches.out | cut -d\. -f 1)
	     
	     printf "${GREEN}# >>> Best IQ-TREE run was: $best_search ...${NC}\n" | \
	     tee -a ${logdir}/get_phylomarkers_run_${dir_suffix}_${TIMESTAMP_SHORT}.log
	     
	     best_tree_file=${tree_prefix}_${best_search_base_name}_nonRecomb_KdeFilt_iqtree_${best_model}.ph
	     cp ${best_search_base_name}.treefile ${best_tree_file}
	     
	     print_start_time && printf "${BLUE}# Adding labels back to ${best_tree_file} ...${NC}\n" | tee -a ${logdir}/get_phylomarkers_run_${dir_suffix}_${TIMESTAMP_SHORT}.log
             [ $DEBUG -eq 1 -o $VERBOSITY -eq 1 ] && echo " > ${distrodir}/add_labels2tree.pl ${tree_labels_dir}/tree_labels.list ${best_tree_file} &> /dev/null"
             ${distrodir}/add_labels2tree.pl ${tree_labels_dir}/tree_labels.list ${best_tree_file} &> /dev/null
             
	     check_output ${best_tree_file} $parent_PID | tee -a ${logdir}/get_phylomarkers_run_${dir_suffix}_${TIMESTAMP_SHORT}.log
	     
	     cp *_ed.ph sorted_IQ-TREE_searches.out $top_markers_dir
	     cd $top_markers_dir
	     rm -rf iqtree_abayes concat_cdnAlns.fnainf.treefile concat_cdnAlns.fnainf.uniqueseq.phy concat_cdnAlns.fnainf.ckp.gz
	  else
	     print_start_time && printf "${BLUE}# running IQ-tree on the concatenated alignment with best model ${best_model} -abayes -bb 1000. This will take a while ...${NC}\n" | \
	     tee -a ${logdir}/get_phylomarkers_run_${dir_suffix}_${TIMESTAMP_SHORT}.log

	     iqtree-omp -s concat_cdnAlns.fnainf -st DNA -m "$best_model" -abayes -bb 1000 -nt AUTO -pre iqtree_abayes &> /dev/null 
	    
	     best_tree_file=${tree_prefix}_nonRecomb_KdeFilt_cdnAlns_iqtree_${best_model}.ph
	     cp iqtree_abayes.treefile ${best_tree_file}
	     
	     print_start_time && printf "${BLUE}# Adding labels back to ${best_tree_file} ...${NC}\n" | tee -a ${logdir}/get_phylomarkers_run_${dir_suffix}_${TIMESTAMP_SHORT}.log
             [ $DEBUG -eq 1 -o $VERBOSITY -eq 1 ] && echo " > ${distrodir}/add_labels2tree.pl ${tree_labels_dir}/tree_labels.list ${best_tree_file} &> /dev/null"
             ${distrodir}/add_labels2tree.pl ${tree_labels_dir}/tree_labels.list ${best_tree_file} &> /dev/null
             
	     check_output ${best_tree_file} $parent_PID | tee -a ${logdir}/get_phylomarkers_run_${dir_suffix}_${TIMESTAMP_SHORT}.log
	     cp *_ed.ph $top_markers_dir
	     cd $top_markers_dir
	     rm -rf iqtree_abayes concat_cdnAlns.fnainf.treefile concat_cdnAlns.fnainf.uniqueseq.phy concat_cdnAlns.fnainf.ckp.gz
	  fi

       fi # if [ "$search_algorithm" == "I" ]
       
        # NOTE: after v0.9 this process is prallelized with run_parallel_cmmds.pl
	if [ $eval_clock -gt 0 ]
        then 
	     echo
	     printf "${YELLOW} >>>>>>>>>>>>>>> TESTING THE MOLECULAR CLOCK HYPOTHESIS <<<<<<<<<<<<<<< ${NC}\n" | \
	     tee -a ${logdir}/get_phylomarkers_run_${dir_suffix}_${TIMESTAMP_SHORT}.log
	     echo
	     
 	     # 1. convert fasta2nexus
             print_start_time && printf "${BLUE}# converting fasta files to nexus files${NC}\n"| \
	     tee -a ${logdir}/get_phylomarkers_run_${dir_suffix}_${TIMESTAMP_SHORT}.log
	     [ $DEBUG -eq 1 -o $VERBOSITY -eq 1 ] && echo " > convert_aln_format_batch_bp.pl fasta fasta nexus nex &> /dev/null" | \
	     tee -a ${logdir}/get_phylomarkers_run_${dir_suffix}_${TIMESTAMP_SHORT}.log
             $distrodir/convert_aln_format_batch_bp.pl fasta fasta nexus nex &> /dev/null
	     
	     # FIX the nexus file format produced by bioperl: (recent paup version error message provided below)
	     # User-defined symbol 'A' conflicts with predefined DNA state symbol.
             # If you are using a predefined format ('DNA', 'RNA', 'nucleotide', or 'protein'), 
             # you may not specify predefined states for this format as symbols in  the Format command.
	     for nexusf in *.nex
	     do
	       perl -pe 'if(/^format /){ s/symbols.*$/;/}' $nexusf > ed && mv ed $nexusf
	     done
	     
	     print_start_time && printf "${BLUE}# Will test the molecular clock hypothesis for $no_top_markers top markers. This will take some time ...${NC}\n" | \
	     tee -a ${logdir}/get_phylomarkers_run_${dir_suffix}_${TIMESTAMP_SHORT}.log
             #run_molecClock_test_jmodeltest2_paup.sh -R 1 -M $base_mod -t ph -e fasta -b molec_clock -q $q &> /dev/null
	     
	     
	    # 2. >>> print table header and append results to it
            no_dot_q=$(echo $q | sed 's/\.//' )
            results_table=mol_clock_M${base_mod}G_r${rooting_method}_o${outgroup_OTU_nomol}_q${no_dot_q}_ClockTest.tab
             
            echo -e "#nexfile\tlnL_unconstr\tlnL_clock\tLRT\tX2_crit_val\tdf\tp-val\tmol_clock" > $results_table
	    
	     cmd="${distrodir}/run_parallel_cmmds.pl nex 'run_parallel_molecClock_test_with_paup.sh -R 1 -f \$file -M $base_mod -t ph -b global_mol_clock -q $q' $n_cores"
	     [ $DEBUG -eq 1 -o $VERBOSITY -eq 1 ] && echo "run_parallel_molecClock.cmd: $cmd" | \
	     tee -a ${logdir}/get_phylomarkers_run_${dir_suffix}_${TIMESTAMP_SHORT}.log
	     echo $cmd | bash &> /dev/null
	     
	     mol_clock_tab=$(ls *_ClockTest.tab)

       	     if [ -s $mol_clock_tab ] 
	     then
	        printf "${GREEN} >>> generated the molecular clock results file $mol_clock_tab ...${NC}\n" | \
		tee -a ${logdir}/get_phylomarkers_run_${dir_suffix}_${TIMESTAMP_SHORT}.log
	    
               # Paste filoinfo and clocklikeness tables
               cut -f1 gene_trees2_concat_tree_RF_distances.tab | grep -v loci | sed 's/"//; s/"$/_/' > list2grep.tmp
               head -1 "$mol_clock_tab" > header.tmp
           
	       # sort lines in molecular clock output file according to order in gene_trees2_concat_tree_RF_distances.tab
	       while read line
	       do
	           grep "$line" $mol_clock_tab 
	       done < list2grep.tmp >> ${mol_clock_tab}sorted
	   
	       cat header.tmp ${mol_clock_tab}sorted > ed && mv ed ${mol_clock_tab}sorted
           
               paste gene_trees2_concat_tree_RF_distances.tab ${mol_clock_tab}sorted > phylogenetic_attributes_of_top${no_top_markers}_gene_trees.tab
	   
	       check_output phylogenetic_attributes_of_top${no_top_markers}_gene_trees.tab $parent_PID | \
	       tee -a ${logdir}/get_phylomarkers_run_${dir_suffix}_${TIMESTAMP_SHORT}.log
           
	       printf "${GREEN} >>> Top markers and associated stats are found in:\n$top_markers_dir ...${NC}\n\n" | \
	       tee -a ${logdir}/get_phylomarkers_run_${dir_suffix}_${TIMESTAMP_SHORT}.log
            else
	          printf "${RED} >>> ${mol_clock_tab} not found" | tee -a ${logdir}/get_phylomarkers_run_${dir_suffix}_${TIMESTAMP_SHORT}.log
	    fi

	 # 5.4 cleanup
         tar -czf molClock_PAUP_files.tgz *_paup.block *.nex  *_clockTest.log *tre *clock.scores *critical_X2_val.R $mol_clock_tab ${mol_clock_tab}sorted mol_clock_MGTRG_r_o_q099_ClockTest.ta*
         [ -s molClock_PAUP_files.tgz ] && rm *_paup.block *.nex  *_clockTest.log *tre *clock.scores *critical_X2_val.R
         [ "$DEBUG" -eq "0" ] && rm list2concat Rplots.pdf header.tmp list2grep.tmp concat_nonRecomb_KdeFilt_cdnAlns_FTGTRG.ph 

	 tar -czf concatenated_alignment_files.tgz concat_cdnAlns.fna concat_cdnAlns.fnainf 
         [ -s concatenated_alignment_files.tgz ] && rm concat_cdnAlns.fna concat_cdnAlns.fnainf 
	 [ "$DEBUG" -eq "0" ] && rm mol_clock_MGTRG_r_o_q099_ClockTest.ta* gene_trees2_concat_tree_RF_distances.tab *cdnAln_allFTGTRG.ph *cdnAln.fasta
	 [ "$DEBUG" -eq "0" ] && rm ../Rplots.pdf ../sorted*perc.tab sorted_aggregated_*tab ../all_*trees.tre

        cd $non_recomb_cdn_alns_dir	
	tar -czf non_recombinant_kdeOK_codon_alignments.tgz ./*_cdnAln.fasta
	[ -s non_recombinant_kdeOK_codon_alignments.tgz ] && rm ./*_cdnAln.fasta
	tar -czf gene_trees_from_non_recombinant_kdeOK_codon_alignments.tgz *_cdnAln*.ph
	[ -s gene_trees_from_non_recombinant_kdeOK_codon_alignments.tgz ] && rm ./*_cdnAln*.ph

	
	cd $top_dir
	tar -czf codon_alignments.tgz ./*_cdnAln.fasta
        [ -s codon_alignments.tgz ] && rm ./*_cdnAln.fasta clustalo.log
        tar -czf protein_alignments.tgz ./*.faaln
        [ -s protein_alignments.tgz ] && rm ./*.faaln
       
      fi
    fi # if [ $runmode -eq 1 ]; then run phylo pipeline on DNA seqs
    

#---------------------------------#    
# >>>>> RUNMODE 2: PopGen <<<<<<< #
#---------------------------------#    

    if [ $runmode -eq 2 ]
    then
        mkdir popGen && cd popGen
	popGen_dir=$(pwd)

        print_start_time && printf "${LBLUE}# Moved into dir popGen ...${NC}\n" | \
	tee -a ${logdir}/get_phylomarkers_run_${dir_suffix}_${TIMESTAMP_SHORT}.log
	
	ln -s ../*fasta .
	no_top_markers=$(ls *fasta | wc -l)
	tmpf=$(ls -1 *fasta | head -1)
	no_seqs=$(grep -c '>' $tmpf)
	[ $DEBUG -eq 1 ] && echo "no_seqs:$no_seqs"  | tee -a ${logdir}/get_phylomarkers_run_${dir_suffix}_${TIMESTAMP_SHORT}.log

        print_start_time && printf "${BLUE}# Will run descriptive DNA polymorphism statistics for $no_top_markers top markers. This will take some time ...${NC}\n"| \
	tee -a ${logdir}/get_phylomarkers_run_${dir_suffix}_${TIMESTAMP_SHORT}.log
	
	TajD_crit_vals=$(get_critical_TajD_values $no_seqs)
	TajD_l=$(echo $TajD_crit_vals|awk '{print $1}')
	TajD_u=$(echo $TajD_crit_vals|awk '{print $2}')

	FuLi_crit_vals=$(get_critical_FuLi_values $no_seqs)
	FuLi_l=$(echo "$FuLi_crit_vals"|awk '{print $1}')
	FuLi_u=$(echo "$FuLi_crit_vals"|awk '{print $2}')
	
	[ $DEBUG -eq 1 ] && echo "TajD_crit_vals:$TajD_crit_vals|TajD_l:$TajD_l|TajD_u:$TajD_u|FuLi_crit_vals:$FuLi_crit_vals|FuLi_l:$FuLi_l|FuLi_u:$FuLi_u" | \
	tee -a ${logdir}/get_phylomarkers_run_${dir_suffix}_${TIMESTAMP_SHORT}.log
	
        print_start_time && printf "${BLUE}# converting $no_top_markers fasta files to nexus format ...${NC}\n" | \
	tee -a ${logdir}/get_phylomarkers_run_${dir_suffix}_${TIMESTAMP_SHORT}.log
        [ $DEBUG -eq 1 -o $VERBOSITY -eq 1 ] && echo " > convert_aln_format_batch_bp.pl fasta fasta nexus nex &> /dev/null" | \
	tee -a ${logdir}/get_phylomarkers_run_${dir_suffix}_${TIMESTAMP_SHORT}.log
	#$distrodir/
	${distrodir}/convert_aln_format_batch_bp.pl fasta fasta nexus nex &> /dev/null 
	  
        print_start_time && printf "${BLUE}# Running popGen_summStats.pl ...${NC}\n" | \
       	tee -a ${logdir}/get_phylomarkers_run_${dir_suffix}_${TIMESTAMP_SHORT}.log
	[ $DEBUG -eq 1 -o $VERBOSITY -eq 1 ] && echo " > popGen_summStats.pl -R 2 -n nex -f fasta -F fasta -H -r 100 -t $TajD_l -T $TajD_u -s $FuLi_l -S $FuLi_u &> popGen_summStats_hs100.log"
	#$distrodir/popGen_summStats.pl -R 2 -n nex -f fasta -F fasta -H -r 100 -t $TajD_l -T $TajD_u -s $FuLi_l -S $FuLi_u &> popGen_summStats_hs100.log
	${distrodir}/popGen_summStats.pl -R 2 -n nex -f fasta -F fasta -H -r 100 -t $TajD_l -T $TajD_u -s $FuLi_l -S $FuLi_u &> popGen_summStats_hs100.log
	
	check_output polymorphism_descript_stats.tab $parent_PID | tee -a ${logdir}/get_phylomarkers_run_${dir_suffix}_${TIMESTAMP_SHORT}.log
	
	printf "${GREEN} >>> descriptive DNA polymorphism stats are found in:\n$popGen_dir ...${NC}\n\n" | \
	tee -a ${logdir}/get_phylomarkers_run_${dir_suffix}_${TIMESTAMP_SHORT}.log
	
	
	#>>> CLEANUP <<<#
	tar -czf clean_cdnAlns.tgz *_cdnAln_clean.fasta *.nex
	[ -s clean_cdnAlns.tgz ] && rm -rf ./*_cdnAln_clean.fasta ./*.nex ./paup.cmd ./popGen_summStats_*.log ./*_cdnAln.fasta
	
	cd $non_recomb_cdn_alns_dir	
	tar -czf non_recombinant_kdeOK_codon_alignments.tgz ./*_cdnAln.fasta
	[ -s non_recombinant_kdeOK_codon_alignments.tgz ] && rm ./*_cdnAln.fasta ./all_*trees.tre ./Rplots.pdf 
	tar -czf gene_trees_from_non_recombinant_kdeOK_codon_alignments.tgz *_cdnAln*.ph
	[ -s gene_trees_from_non_recombinant_kdeOK_codon_alignments.tgz ] && rm ./*_cdnAln*.ph

	cd $top_dir
	tar -czf codon_alignments.tgz ./*_cdnAln.fasta
        [ -s codon_alignments.tgz ] && rm ./*_cdnAln.fasta clustalo.log
        tar -czf protein_alignments.tgz ./*.faaln
        [ -s protein_alignments.tgz ] && rm ./*.faaln
    fi
fi # if [ "$mol_type" == "DNA"    

#----------------------------------------------------------------------------------------------------------------
#>>>BLOCK 5. Compute individual gene trees from codon or protein alignments
# 5.1 run FastTree in parallel on all codon alignments, using exhausitve tree-rearrangements, 
#            turn off NNI heuristics, and to always optimize all 5 branches at each NNI, 
#            optimizing all 5 branches in as many as 3 rounds
#----------------------------------------------------------------------------------------------------------------

if [ "$mol_type" == "PROT" ]
then
    cd non_recomb_FAA_alns
    non_recomb_FAA_alns_dir=$(pwd)

    print_start_time && printf "${BLUE}# estimating $no_non_recomb_alns_perm_test gene trees from non-recombinant sequences ...${NC}\n" | \
    tee -a ${logdir}/get_phylomarkers_run_${dir_suffix}_${TIMESTAMP_SHORT}.log
    [ "$DEBUG" -eq 1 -o "$VERBOSITY" -eq 1 ] && echo " > ${distrodir}/run_parallel_cmmds.pl faaln 'FastTree -quiet -lg -gamma -bionj -slownni -mlacc 3 -spr 8 -sprlength 8 < $file > ${file%.*}_allFTlgG.ph'  $n_cores &> /dev/null"
    ${distrodir}/run_parallel_cmmds.pl faaln 'FastTree -quiet -lg -gamma -bionj -slownni -mlacc 3 -spr 8 -sprlength 8 < $file > ${file%.*}_allFTlgG.ph' $n_cores &> /dev/null
    
    #remove trees with < 5 branches
    print_start_time && printf "${BLUE}# counting branches on $no_non_recomb_alns_perm_test gene trees ...${NC}\n" | \
    tee -a ${logdir}/get_phylomarkers_run_${dir_suffix}_${TIMESTAMP_SHORT}.log
    [ "$DEBUG" -eq 1 -o "$VERBOSITY" -eq 1 ] && echo " > count_tree_branches ph no_tree_branches.list &> /dev/null"  | \
    tee -a ${logdir}/get_phylomarkers_run_${dir_suffix}_${TIMESTAMP_SHORT}.log
    count_tree_branches ph no_tree_branches.list #&> /dev/null
   
    check_output no_tree_branches.list "$parent_PID" | tee -a ${logdir}/get_phylomarkers_run_${dir_suffix}_${TIMESTAMP_SHORT}.log
    
    # remove trees with < 5 external branches (leaves)
    [ $DEBUG -eq 1 -o $VERBOSITY -eq 1 ] && echo " > removing trees with < 5 external branches (leaves)" | \
    tee -a ${logdir}/get_phylomarkers_run_${dir_suffix}_${TIMESTAMP_SHORT}.log
    for phy in $(grep -v '^#Tree' no_tree_branches.list|awk -v min_no_ext_branches=$min_no_ext_branches 'BEGIN{FS="\t"; OFS="\t"}$7 < min_no_ext_branches'|cut -f1)
    do
         base=$(echo $phy|sed 's/_allFTlgG\.ph//')
	 printf "${LRED} >>> will remove ${base}* because it has < 5 branches!${NC}\n" | tee -a ${logdir}/get_phylomarkers_run_${dir_suffix}_${TIMESTAMP_SHORT}.log
	 rm ${base}*
    done

    # 6.1 generate the all_GTRG_trees.tre holding all source trees. This is required for run_kdetrees.R ph
    #     Make a check for the existence of the file to interrupt the pipeline if something has gone wrong
    [ "$DEBUG" -eq 1 -o "$VERBOSITY" -eq 1 ] && echo " > cat *_allFTlgG.ph > all_FTlgG_trees.tre" | \
    tee -a ${logdir}/get_phylomarkers_run_${dir_suffix}_${TIMESTAMP_SHORT}.log
    cat *_allFTlgG.ph > all_FTlgG_trees.tre
    check_output all_FTlgG_trees.tre "$parent_PID" | tee -a ${logdir}/get_phylomarkers_run_${dir_suffix}_${TIMESTAMP_SHORT}.log
    [ ! -s all_FTlgG_trees.tre ] && exit 5
    
    # 6.2 run_kdetrees.R at desired stringency 
    print_start_time && printf "${BLUE}# running kde test ...${NC}\n" | tee -a ${logdir}/get_phylomarkers_run_${dir_suffix}_${TIMESTAMP_SHORT}.log
    [ $DEBUG -eq 1 -o $VERBOSITY -eq 1 ] && echo " > run_kdetrees.R ph all_FTlgG_trees.tre $kde_stringency &> /dev/null" | \
    tee -a ${logdir}/get_phylomarkers_run_${dir_suffix}_${TIMESTAMP_SHORT}.log
    run_kdetrees.R ph all_FTlgG_trees.tre $kde_stringency &> /dev/null
    check_output kde_dfr_file_all_FTlgG_trees.tre.tab $parent_PID | tee -a ${logdir}/get_phylomarkers_run_${dir_suffix}_${TIMESTAMP_SHORT}.log

    # 6.3 mv outliers to kde_outliers
    no_kde_outliers=$(grep -c outlier kde_dfr_file_all_FTlgG_trees.tre.tab)
    no_kde_ok=$(grep -vc outlier kde_dfr_file_all_FTlgG_trees.tre.tab)

    if [ "$no_kde_outliers" -gt "0" ]
    then
        print_start_time && printf "${BLUE}# making dir kde_outliers/ and moving $no_kde_outliers outlier files into it ...${NC}\n" | \
	tee -a ${logdir}/get_phylomarkers_run_${dir_suffix}_${TIMESTAMP_SHORT}.log
        mkdir kde_outliers
        for f in $(grep outlier kde_dfr_file_all_FTlgG_trees.tre.tab|cut -f1|sed 's/_allFTlgG.ph//'); do mv ${f}* kde_outliers; done
    else
        printf "${GREEN} >>> there are no kde-test outliers ...${NC}\n" | tee -a ${logdir}/get_phylomarkers_run_${dir_suffix}_${TIMESTAMP_SHORT}.log
    fi

    if [ $no_kde_ok -gt 0 ]
    then
        print_start_time && printf "${BLUE}# making dir kde_ok/ and linking $no_kde_ok selected files into it ...${NC}\n" | \
	tee -a ${logdir}/get_phylomarkers_run_${dir_suffix}_${TIMESTAMP_SHORT}.log
        mkdir kde_ok
        cd kde_ok
        ln -s ../*ph .
   
        print_start_time && printf "${BLUE}# labeling $no_kde_ok gene trees in dir kde_ok/ ...${NC}\n" | \
	tee -a ${logdir}/get_phylomarkers_run_${dir_suffix}_${TIMESTAMP_SHORT}.log
        [ $DEBUG -eq 1 -o $VERBOSITY -eq 1 ] && echo " > ${distrodir}/run_parallel_cmmds.pl ph 'add_labels2tree.pl ../../../tree_labels.list $file' $n_cores &> /dev/null"  | \
	tee -a ${logdir}/get_phylomarkers_run_${dir_suffix}_${TIMESTAMP_SHORT}.log
        ${distrodir}/run_parallel_cmmds.pl ph 'add_labels2tree.pl ../../../tree_labels.list $file' $n_cores &> /dev/null
    
        # remove symbolic links to cleanup kde_ok/
        for f in $(ls *ph|grep -v '_ed\.ph'); do rm $f; done
    
        cd ..
    else
         printf "${RED}# There are $no_kde_ok gene trees producing non-significant kde-test results! Increase the actual -k $kde_stringency value. Will stop here. ...${NC}\n" | \
	 tee -a ${logdir}/get_phylomarkers_run_${dir_suffix}_${TIMESTAMP_SHORT}.log
	 exit 5
    fi
    
    # 6.4 compute MJRE consensus tree with consense; does not work with trees with different number of leaves
    # cat *ph > intree
    # echo "y" > consense.cmd

    # consense < consense.cmd
    # check_output outtree $parent_PID
    
    #  6.5 concatenated prefix required by run compute_suppValStas_and_RF-dist.R
    #  ln -s outtree concatenated.ph 

    # 6.6 compute average bipartition support values for each gene tree
    #     and their RF-fistance to the consensus tree computed with consense
    wkdir=$(pwd) 
        
    print_start_time && printf "${BLUE}# computing tree support values ...${NC}\n" | tee -a ${logdir}/get_phylomarkers_run_${dir_suffix}_${TIMESTAMP_SHORT}.log
    [ $DEBUG -eq 1 -o $VERBOSITY -eq 1 ] && echo " > compute_suppValStas_and_RF-dist.R $wkdir 1 faaln ph 1 &> /dev/null"  | \
    tee -a ${logdir}/get_phylomarkers_run_${dir_suffix}_${TIMESTAMP_SHORT}.log
    compute_suppValStas_and_RF-dist.R $wkdir 1 faaln ph 1 &> /dev/null
    
    print_start_time && printf "${BLUE}# writing summary tables ...${NC}\n" | tee -a ${logdir}/get_phylomarkers_run_${dir_suffix}_${TIMESTAMP_SHORT}.log
    min_supp_val_perc=${min_supp_val#0.}
    no_digits=${#min_supp_val_perc}
    [ $no_digits -eq 1 ] && min_supp_val_perc=${min_supp_val_perc}0
    
    awk -v min_supp_val=$min_supp_val '$2 >= min_supp_val' sorted_aggregated_support_values4loci.tab > sorted_aggregated_support_values4loci_ge${min_supp_val_perc}perc.tab
    check_output sorted_aggregated_support_values4loci_ge${min_supp_val_perc}perc.tab $parent_PID | tee -a ${logdir}/get_phylomarkers_run_${dir_suffix}_${TIMESTAMP_SHORT}.log
    
    no_top_markers=$(wc -l sorted_aggregated_support_values4loci_ge${min_supp_val_perc}perc.tab|cut -d' ' -f1)
    top_markers_dir="top_${no_top_markers}_markers_ge${min_supp_val_perc}perc"
    top_markers_tab=$(ls sorted_aggregated_support_values4loci_ge${min_supp_val_perc}perc.tab)
    
    print_start_time && printf "${LBLUE}# making dir $top_markers_dir and moving $no_top_markers top markers into it ...${NC}\n" | \
    tee -a ${logdir}/get_phylomarkers_run_${dir_suffix}_${TIMESTAMP_SHORT}.log
    mkdir $top_markers_dir && cd $top_markers_dir
    ln -s ../$top_markers_tab .
    for base in $(awk '{print $1}' $top_markers_tab|grep -v loci|sed 's/"//g'); do ln -s ../${base}* .; done

    # 6.7 generate supermatrix (concatenated alignment) 
    print_start_time && printf "${BLUE}# concatenating $no_top_markers top markers into supermatrix ...${NC}\n" | tee -a ${logdir}/get_phylomarkers_run_${dir_suffix}_${TIMESTAMP_SHORT}.log
    [ $DEBUG -eq 1 -o $VERBOSITY -eq 1 ] && echo " > concat_alns faaln $parent_PID &> /dev/null" | \
    tee -a ${logdir}/get_phylomarkers_run_${dir_suffix}_${TIMESTAMP_SHORT}.log
    concat_alns faaln $parent_PID &> /dev/null

    # 6.8 remove uninformative sites from the concatenated alignment to speed up computation
        print_start_time && printf "${BLUE}# removing uninformative sites from concatenated alignment ...${NC}\n" | tee -a ${logdir}/get_phylomarkers_run_${dir_suffix}_${TIMESTAMP_SHORT}.log
    $distrodir/remove_uninformative_sites_from_aln.pl < concat_protAlns.faa > concat_protAlns.faainf

    check_output concat_protAlns.faainf $parent_PID|tee -a ${logdir}/get_phylomarkers_run_${dir_suffix}_${TIMESTAMP_SHORT}.log

    # 6.9 run FasTree under the GTR+G model 
    print_start_time && printf "${BLUE}# running FastTree on the concatenated alignment with $search_thoroughness thoroughness. This may take a while ...${NC}\n" | \
    tee -a ${logdir}/get_phylomarkers_run_${dir_suffix}_${TIMESTAMP_SHORT}.log
    
    if [ "$search_thoroughness" == "high" ]
    then
        FastTree -quiet -lg -bionj -slownni -gamma -mlacc 3 -spr $spr -sprlength $spr_length -log ${tree_prefix}_nonRecomb_KdeFilt_protAlns_FTlgG.log < concat_protAlns.faainf > ${tree_prefix}_nonRecomb_KdeFilt_protAlns_FTlgG.ph
    fi

    if [ "$search_thoroughness" == "medium" ]
    then
        FastTree -quiet -lg -bionj -slownni -gamma -mlacc 2 -spr $spr -sprlength $spr_length -log ${tree_prefix}_nonRecomb_KdeFilt_protAlns_FTlgG.log < concat_protAlns.faainf > ${tree_prefix}_nonRecomb_KdeFilt_protAlns_FTlgG.ph
    fi

    if [ "$search_thoroughness" == "low" ]
    then
        FastTree -quiet -lg -bionj -gamma -spr $spr -sprlength $spr_length -log ${tree_prefix}_nonRecomb_KdeFilt_protAlns_FTlgG.log < concat_protAlns.faainf > ${tree_prefix}_nonRecomb_KdeFilt_protAlns_FTlgG.ph
    fi

    if [ "$search_thoroughness" == "lowest" ]
    then
        FastTree -quiet -lg -gamma -mlnni 4 -log ${tree_prefix}_nonRecomb_KdeFilt_protAlns_FTlgG.log < concat_protAlns.faainf > ${tree_prefix}_nonRecomb_KdeFilt_protAlns_FTlgG.ph
    fi

    check_output ${tree_prefix}_nonRecomb_KdeFilt_protAlns_FTlgG.ph $parent_PID | tee -a ${logdir}/get_phylomarkers_run_${dir_suffix}_${TIMESTAMP_SHORT}.log
    
    
   if [ -s "${tree_prefix}_nonRecomb_KdeFilt_protAlns_FTlgG.log" -a -s "${tree_prefix}_nonRecomb_KdeFilt_protAlns_FTlgG.ph" ]
   then
       lnL=$(grep ML_Lengths2 "${tree_prefix}_nonRecomb_KdeFilt_protAlns_FTlgG.log" | grep TreeLogLk | sed 's/TreeLogLk[[:space:]]ML_Lengths2[[:space:]]//')
       printf "${GREEN} >>> lnL for ${tree_prefix}_nonRecomb_KdeFilt_protAlns_FTlgG.ph = $lnL ${NC}\n" | tee -a ${logdir}/get_phylomarkers_run_${dir_suffix}_${TIMESTAMP_SHORT}.log
   else
       printf "${LRED} >>> WARNING: ${tree_prefix}_nonRecomb_KdeFilt_protAlns_FTlgG.log could not be produced!${NC}\n" | tee -a ${logdir}/get_phylomarkers_run_${dir_suffix}_${TIMESTAMP_SHORT}.log
   fi

    
    print_start_time && printf "${BLUE}# Adding labels back to tree ...${NC}\n"|tee -a ${logdir}/get_phylomarkers_run_${dir_suffix}_${TIMESTAMP_SHORT}.log
    [ $DEBUG -eq 1 -o $VERBOSITY -eq 1 ] && echo " > add_labels2tree.pl ${tree_labels_dir}/tree_labels.list ${tree_prefix}_nonRecomb_KdeFilt_protAlns_FTlgG.ph &> /dev/null" | \
    tee -a ${logdir}/get_phylomarkers_run_${dir_suffix}_${TIMESTAMP_SHORT}.log
    add_labels2tree.pl ${tree_labels_dir}/tree_labels.list ${tree_prefix}_nonRecomb_KdeFilt_protAlns_FTlgG.ph &> /dev/null
    
    [ -s ${tree_prefix}_nonRecomb_KdeFilt_protAlns_FTlgG_ed.ph ] && \
    mv ${tree_prefix}_nonRecomb_KdeFilt_protAlns_FTlgG_ed.ph ${tree_prefix}_${no_top_markers}nonRecomb_KdeFilt_protAlns_FTlgG.spTree
    [ -s ${tree_prefix}_nonRecomb_KdeFilt_protAlns_FTlgG_ed.ph ] && rm ${tree_prefix}_nonRecomb_KdeFilt_protAlns_FTlgG.ph
    
    check_output ${tree_prefix}_${no_top_markers}nonRecomb_KdeFilt_protAlns_FTlgG.spTree $parent_PID | tee -a ${logdir}/get_phylomarkers_run_${dir_suffix}_${TIMESTAMP_SHORT}.log
     
    printf "${GREEN} >>> found in dir $wkdir ...${NC}\n\n" | tee -a ${logdir}/get_phylomarkers_run_${dir_suffix}_${TIMESTAMP_SHORT}.log
       
    wkdir=$(pwd)
    [ $DEBUG -eq 1 -o $VERBOSITY -eq 1 ] && echo " > compute_suppValStas_and_RF-dist.R $wkdir 2 faaln ph 1 &> /dev/null" | \
    tee -a ${logdir}/get_phylomarkers_run_${dir_suffix}_${TIMESTAMP_SHORT}.log
    ${distrodir}/compute_suppValStas_and_RF-dist.R $wkdir 2 faaln ph 1 &> /dev/null
    
    
       if [ "$search_algorithm" == "I" ]
       then
	  # 5.5 run IQ-tree in addition to FastTree, if requested
          echo
	  printf "${YELLOW} >>>>>>>>>>>>>>> IQ-TREE + ModelFinder run <<<<<<<<<<<<<<< ${NC}\n" | \
	  tee -a ${logdir}/get_phylomarkers_run_${dir_suffix}_${TIMESTAMP_SHORT}.log
	  echo 

	  print_start_time && printf "${BLUE}# running ModelFinder on the concatenated alignment with $IQT_models. This will take a while ...${NC}\n" | \
	  tee -a ${logdir}/get_phylomarkers_run_${dir_suffix}_${TIMESTAMP_SHORT}.log
           
          iqtree-omp -s concat_protAlns.faainf -st PROT -mset "$IQT_models" -m MF -nt AUTO &> /dev/null 
	  
	  check_output concat_protAlns.faainf.log $parent_PID | tee -a ${logdir}/get_phylomarkers_run_${dir_suffix}_${TIMESTAMP_SHORT}.log
	  
	  best_model=$(grep '^Best-fit model' concat_protAlns.faainf.log | cut -d' ' -f 3)
	  printf "${GREEN} >>> Best-fit model: ${best_model} ...${NC}\n" | \
	  tee -a ${logdir}/get_phylomarkers_run_${dir_suffix}_${TIMESTAMP_SHORT}.log
	  
	  mkdir iqtree_abayes && cd iqtree_abayes
	  ln -s ../concat_protAlns.faainf .
	  
	  if [ "$search_thoroughness" == "high" ]
	  then
	     print_start_time && printf "${BLUE}# running IQ-TREE on the concatenated alignment with best model ${best_model} -abayes -bb 1000, starting from $nrep random trees!. This will take a while ...${NC}\n" | \
	     tee -a ${logdir}/get_phylomarkers_run_${dir_suffix}_${TIMESTAMP_SHORT}.log
	     
	     # run nrep IQ-TREE searches under the best-fit model found
	     for ((rep=1;rep<=$nrep;rep++))
	     do 
	         print_start_time && printf "${LBLUE} > iqtree-omp -s concat_protAlns.faainf -st DNA -m "$best_model" -abayes -bb 1000 -nt AUTO -pre abayes_run${rep} &> /dev/null${NC}\n" | \
	         tee -a ${logdir}/get_phylomarkers_run_${dir_suffix}_${TIMESTAMP_SHORT}.log
	        
		  iqtree-omp -s concat_protAlns.faainf -st PROT -m "$best_model" -abayes -bb 1000 -nt AUTO -pre abayes_run${rep} &> /dev/null 
	     done
             
	     grep '^BEST SCORE' *log | sort -nrk5 > sorted_IQ-TREE_searches.out
	     
	     check_output sorted_IQ-TREE_searches.out $parent_PID | tee -a ${logdir}/get_phylomarkers_run_${dir_suffix}_${TIMESTAMP_SHORT}.log
	     best_search=$(head -1 sorted_IQ-TREE_searches.out)
	     best_search_base_name=$(head -1 sorted_IQ-TREE_searches.out | cut -d\. -f 1)
	     
	     printf "${GREEN}# >>> Best IQ-TREE run was: $best_search ...${NC}\n" | \
	     tee -a ${logdir}/get_phylomarkers_run_${dir_suffix}_${TIMESTAMP_SHORT}.log
	     
	     best_tree_file=${tree_prefix}_${best_search_base_name}_nonRecomb_KdeFilt_iqtree_${best_model}.ph
	     cp ${best_search_base_name}.treefile ${best_tree_file}
	     
	     print_start_time && printf "${BLUE}# Adding labels back to ${best_tree_file} ...${NC}\n" | tee -a ${logdir}/get_phylomarkers_run_${dir_suffix}_${TIMESTAMP_SHORT}.log
             [ $DEBUG -eq 1 -o $VERBOSITY -eq 1 ] && echo " > ${distrodir}/add_labels2tree.pl ${tree_labels_dir}/tree_labels.list ${best_tree_file} &> /dev/null"
             ${distrodir}/add_labels2tree.pl ${tree_labels_dir}/tree_labels.list ${best_tree_file} &> /dev/null
             
	     check_output ${best_tree_file} $parent_PID | tee -a ${logdir}/get_phylomarkers_run_${dir_suffix}_${TIMESTAMP_SHORT}.log
	     cp *_ed.ph sorted_IQ-TREE_searches.out $wkdir
	     cd $wkdir
	     rm -rf iqtree_abayes concat_protAlns.faainf.treefile concat_protAlns.faainf.uniqueseq.phy *ckp.gz
	  else
	     print_start_time && printf "${BLUE}# running IQ-tree on the concatenated alignment with best model ${best_model} -abayes -bb 1000. This will take a while ...${NC}\n" | \
	     tee -a ${logdir}/get_phylomarkers_run_${dir_suffix}_${TIMESTAMP_SHORT}.log

	     print_start_time && printf "${BLUE}# running: iqtree-omp -s concat_protAlns.faainf -st PROT -m $best_model -abayes -bb 1000 -nt AUTO -pre iqtree_abayes &> /dev/null  ...${NC}\n" | \
	     tee -a ${logdir}/get_phylomarkers_run_${dir_suffix}_${TIMESTAMP_SHORT}.log
	     iqtree-omp -s concat_protAlns.faainf -st PROT -m "$best_model" -abayes -bb 1000 -nt AUTO -pre iqtree_abayes &> /dev/null 
	    
	     best_tree_file=${tree_prefix}_nonRecomb_KdeFilt_protAlns_iqtree_${best_model}.ph
	     cp iqtree_abayes.treefile ${best_tree_file}
	     
	     print_start_time && printf "${BLUE}# Adding labels back to ${best_tree_file} ...${NC}\n" | tee -a ${logdir}/get_phylomarkers_run_${dir_suffix}_${TIMESTAMP_SHORT}.log
             [ $DEBUG -eq 1 -o $VERBOSITY -eq 1 ] && echo " > ${distrodir}/add_labels2tree.pl ${tree_labels_dir}/tree_labels.list ${best_tree_file} &> /dev/null"
             ${distrodir}/add_labels2tree.pl ${tree_labels_dir}/tree_labels.list ${best_tree_file} &> /dev/null
	                  
	     check_output ${best_tree_file} $parent_PID | tee -a ${logdir}/get_phylomarkers_run_${dir_suffix}_${TIMESTAMP_SHORT}.log
	     cp *_ed.ph $wkdir
	     cd $wkdir
	     rm -rf iqtree_abayes concat_protAlns.faainf.treefile concat_protAlns.faainf.uniqueseq.phy *ckp.gz
	  fi

       fi # if [ "$search_algorithm" == "I" ]
       

   # >>> 6.9 CLEANUP <<< #
   [ "$DEBUG" -eq "0" ] && rm list2concat Rplots.pdf sorted*perc.tab concat_nonRecomb_KdeFilt_protAlns_FT*.ph | \
    tee -a ${logdir}/get_phylomarkers_run_${dir_suffix}_${TIMESTAMP_SHORT}.log
   
   rm ./*allFT*.ph ./*faaln

   tar -czf concatenated_alignment_files.tgz concat_protAlns.faa concat_protAlns.faainf
   [ -s concatenated_alignment_files.tgz ] && rm concat_protAlns.faa concat_protAlns.faainf
   [ "$DEBUG" -eq "0" ] && rm  ../Rplots.pdf ../all_FT*trees.tre ../sorted*perc.tab sorted_aggregated_*tab | \
    tee -a ${logdir}/get_phylomarkers_run_${dir_suffix}_${TIMESTAMP_SHORT}.log
    
   cd $non_recomb_FAA_alns_dir
   tar -czf non_recomb_kdeOK_FAA_alignments.tgz ./*_cluo.faaln
   [ -s non_recomb_kdeOK_FAA_alignments.tgz ] && rm ./*_cluo.faaln
   tar -czf non_recomb_kdeOK_prot_trees.tgz ./*_cluo_*.ph
   [ -s non_recomb_kdeOK_prot_trees.tgz ] && rm ./*_cluo_*.ph
   
   cd $top_dir
   tar -czf codon_alignments.tgz *_cdnAln.fasta
   [ -s codon_alignments.tgz ] && rm *_cdnAln.fasta clustalo.log
   tar -czf protein_alignments.tgz *.faaln
   [ -s protein_alignments.tgz ] && rm *.faaln
fi

# compute the elapsed time since the script was fired
end_time=`date +%s`
secs=$(($end_time-$start_time))

#printf '%dh:%dm:%ds\n' $(($secs/3600)) $(($secs%3600/60)) $(($secs%60))
printf "\n${LBLUE} >>> Total runtime of $progname: ${NC}" 
printf '%dh:%dm:%ds\n' $(($secs/3600)) $(($secs%3600/60)) $(($secs%60))
echo

cat <<REF

* PROVISIONAL CITATION: 

If you find the code useful for your academic work, please use the following citation: 
Pablo Vinuesa and Bruno Contreras-Moreira 2017. Get_PhyloMarkers, a pipeline to select 
optimal markers for microbial phylogenomics, population genetics and genomic taxomy. 
Available at https://github.com/vinuesa/get_phylomarkers 

A publication is in preparation. The abstract was accepted in Frontiers in Microbiology, 
for the Research topic on "microbial taxonomy, phylogeny and biodiversity" 
http://journal.frontiersin.org/researchtopic/5493/microbial-taxonomy-phylogeny-and-biodiversity

* NOTES: 
  1. The links to the the corresponding manuscript will be provided here 
      as soon as it is available at bioRxiv, and latter, to the paper.

  2. If you encounter problems or bugs while trying to run the pipeline
      please report them through the issues page at
      https://github.com/vinuesa/get_phylomarkers/issues
      
      Alternatively, see our contact details at:
      http://www.ccg.unam.mx/~vinuesa/
      https://digital.csic.es/cris/rp/rp02661/
      
      Please run the script with the -D flag added at the end of the 
      command line and send us the output, so that we can better diagnose 
      the problem.
      
      Thanks!

REF



