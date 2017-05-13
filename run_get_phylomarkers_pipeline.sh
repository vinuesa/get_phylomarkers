#!/usr/bin/env bash

#: PROGRAM: run_get_phylomarkers_pipeline.sh
#: AUTHOR: Pablo Vinuesa, Center for Genomic Sciences, UNAM, Mexico
#:         http://www.ccg.unam.mx/~vinuesa/
#
#: PROJECT START: April 2017; This is this is a wrapper script to automate the whole process of marker selection and downstream analyses.
#
#: AIM: select optimal molecular markers for phylogenomics and population genomics from orthologous gene clusters computed by get_homologues
#
#: OUTPUT: multiple sequence alignments (of protein and DNA sequences) of selected markers, gene and supermatrix phylogenies, 
#          along with summary graphics and tables summarizing the filtering procedure. If requested. 
#          

progname=${0##*/} # run_get_phylomarkers_pipeline.pl
VERSION='1.2_13May17' #v1.2_13May17 refined the logic of set_bindirs(); added get_start_time(); improved error checking code, including get_script_PID()
                      # fixed a bug in -t PROT 
                      
		      #v1.1_13May17; Major update to facilitate installation by users: added set_pipeline_environment(), check_homebinpath(), set_bindirs(), 
                      #              which ckeck and setup the ENVIRONMENT for the pipeline. Added new bin/ and pop_gen_tables/ dirs for consistency. 
		      #              The bin/linux dir contains the 2cnd party binaries compiled for 64bit linux machines. Need to compile for darwin
		      
                      #v1.0_11May17; git commited; first version on GitHub: https://github.com/vinuesa/get_phylomarkers
                      # v0.9_10May17 important speed improvement due to running pal2nal.pl and run_parallel_molecClock_test_with_paup.sh
                      #               in parallel with run_pexec_cmmds.sh 
                     # v0.8_9May17, added -R 1|2, for Phylo|popGen, based on popGen_summStats.pl and pre-computed Tajima's D and Fu-Li crit value tables
                     #        with missing values predicted by lm() in R, based on the original critical values in Tajima's 89 and  Fu and Li 1993 papers
		     #        The critical CI values are gathered by get_critical_TajD_values() and get_critical_FuLi_values()
                     # v0.7_2May17 added -e min_no_ext_branches -c $codontable and -C print_codontable() for pal2nal.pl, which is now run from within a loop
                     #             runs compute_suppValStas_and_RF-dist.R in the top_X_markers dir
                     # v0.6_1May17 added get_script_PID() and count_tree_branches(); added -t BOTH
                     # v0.5_29April17 Added -t PROT
                     # v0.4_27April17 Added print_usage_notes() and make_labels_4_trees()
# Set GLOBALS

wkdir=$(pwd) #echo "# working in $wkdir"

DATEFORMAT_SHORT="%d%b%y"
TIMESTAMP_SHORT=$(date +${DATEFORMAT_SHORT})

DATEFORMAT_HMS="%H:%M:%S"
TIMESTAMP_HMS=$(date +${DATEFORMAT_HMS})

TIMESTAMP_SHORT_HMS=$(date +${DATEFORMAT_SHORT}-${DATEFORMAT_HMS})

#>>> set color in bash 
#  SEE: echo http://stackoverflow.com/questions/5947742/how-to-change-the-output-color-of-echo-in-linux
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
    elif [[ "$OSTYPE" == "darwin"* ]]
    then
        scriptdir=$(greadlink -f ${BASH_SOURCE[0]})
	distrodir=$(dirname $scriptdir) #echo "scriptdir: $scriptdir|basedir:$distrodir|OSTYPE:$OSTYPE"
        bindir=$distrodir/bin/MacOSX
	OS='darwin'
    fi
    echo "$distrodir $bindir $OS"
}
#----------------------------------------------------------------------------------------- 

function set_bindirs()
{  
    # $bindir $homebinflag $homebinpathflag
    bindir=$1
    homebinpathflag=$2
    
    setbindir_flag=0
    bins=( clustalo FastTree pexec Phi paup consense )

    for prog in "${bins[@]}" 
    do
       bin=$(type -P $prog)
       if [ -z $bin ]; then
          echo
          printf "${RED}# $prog not found in \$PATH ... ${NC}\n"
	  if [ $homebinpathflag -eq 1 ] 
	  then
          printf " >>> ${CYAB}# will generate a softlink in $HOME/bin to $prog ${NC}\n"
	        ln -s $bindir/$prog $HOME/bin 
	  else
                printf " >>> ${CYAN} will append $bindir to the \$PATH variable${NC}\n"
                PATH=$PATH:$bindir/$prog # append $HOME/bin to $PATH 
	  fi    
       fi
    done
}
#----------------------------------------------------------------------------------------- 

function check_homebinpath()
{
   distrodir=$1
   
   homebinflag=
   homebinpathflag=
   
   if [ -d $HOME/bin ]
   then
       homebinflag=1
   else
      printf "${RED} No $HOME/bin directory found, will append $distrodir to the \$PATH variable${NC}"
      PATH=$PATH:$HOME/$distrodir # append $HOME/bin to $PATH 
   fi
   
   if [ -d $(echo $PATH | sed 's/:/\n/g' | grep "$HOME/bin") ] 
   then
        homebinpathflag=1
	
	printf "${CYAN} Will generate symlinks in $HOME/bin to Bash, Perl and R scripts in $distrodir ...${NC}"
        ln -s $distrodir/*.sh $HOME/bin &> /dev/null
	ln -s $distrodir/*.R $HOME/bin &> /dev/null
	ln -s $distrodir/*.pl $HOME/bin &> /dev/null
   fi
   echo "$homebinflag $homebinpathflag"
}
#----------------------------------------------------------------------------------------- 

function check_scripts_in_path()
{
    distrodir=$1
    
    #bins=( clustalo FastTree pexec Phi paup consense )
    bash_scripts=( run_pexec_cmmds.sh run_parallel_molecClock_test_with_paup.sh )
    perl_scripts=( add_nos2fasta_header.pl pal2nal.pl rename concat_alignments.pl add_labels2tree.pl convert_aln_format_batch_bp.pl popGen_summStats.pl convert_aln_format_batch_bp.pl )
    R_scripts=( run_kdetrees.R compute_suppValStas_and_RF-dist.R )
    
    for prog in "${bash_scripts[@]}" "${perl_scripts[@]}" "${R_scripts[@]}"
    do
       bin=$(type -P $prog)
       if [ -z $bin ]; then
          echo
          printf "${RED}# WARNING: script $prog is not in \$PATH!${NC}\n"
	  printf "${CYAN}  >>>  Will add it to \$PATH ${NC}\n"
	  homebinpath_flags=$(check_homebinpath $distrodir)
       fi
    done
    echo "$homebinpath_flags"
}
#----------------------------------------------------------------------------------------- 

#function check_binaries_in_path()
#{
#    bindir=$1
#    
#    bins=( clustalo FastTree pexec Phi paup consense )
#    
#    for prog in "${bins[@]}"
#    do
#       bin=$(type -P $prog)
#       if [ -z $bin ]; then
#          echo
#          printf "${RED}# WARNING: binary $prog is not in \$PATH!\n${NC}"
#	  printf "${CYAN} >>> Will append $bindir to \$PATH ${NC}\n"
#	  #check_homebinpath $distrodir
#          #echo "# ... you will need to install \"$prog\" first or include it in \$PATH"
#          #echo "# ... exiting"
##          exit 1
#       fi
#    done
#}
#----------------------------------------------------------------------------------------- 

#run_pexec_cmmds()
#{
#  # run pexec to parallelize processes on all available cores, or user-defined no. of cores
#  #   requires 3 or max 4 args
#  #  <file extension name> <'command'> [no_of_cores]
#
#   ext=$1
#   command=$2
#   suffix=$3
#   no_of_cores=$4
#
#   total_files=$(ls *$ext | wc -l)
#
#   if [ -z $no_of_cores ]
#   then
#       # use all available cores
#       pexec -r *.$ext -e file -c -o - -- "$command; [ $VERBOSITY -eq 1 ] && echo \$file \=\=\> processed!; done"  
#   else
#       # use defined cores
#       pexec -n $no_of_cores -r *.$ext -e file -c -o - -- "$command; [ $VERBOSITY -eq 1 ] && echo \$file \=\=\> processed!; done"
#   fi
#
#   echo 
#   printf  "${GREEN}run_pexec_cmmds is done executing $command on $total_files input files ...${NC}"
#   echo 
#}
#----------------------------------------------------------------------------------------- 

function print_usage_notes()
{
   cat <<USAGE
   
   $progname v$VERSION important usage notes.
   
   1. Start the run from within the directory holding core genome clusters generated by compare_clusters.pl
   
      NOTE: both faa and fna files are required to generate codon alignments from DNA fasta files. This
            means that two runs of compare_clusters.pl (from the get_homologues package) are required,
	    one of them using the -n flag.
	    
   2. $progname is intended to run on a collection of genomes from different species. 
   
      NOTES: an absolute minimum of 4 distinct genomes are required. 
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
	 ...
	 
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
    prog=${1%.*} # remove the script's .extension_name
    #proc_ID=$(ps -eaf | grep "$prog" | grep -v grep | grep '-' | grep $USER | awk '{print $2}')
    proc_ID=$(ps aux | grep "$prog" | grep -v grep | grep '-' | grep $USER | awk '{print $2}')
    echo $proc_ID
    [ $DEBUG ] && echo "$progname PID is: $proc_ID"
}
#----------------------------------------------------------------------------------------- 

function count_tree_labels()
{
   # call R pacake 'ape' to count the number of labels in a tree
   ext=$1
   outfile=$2

   R --no-save <<RCMD

   library("ape")

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

   R --no-save <<RCMD &> /dev/null

   library("ape")

   trees <- list.files(pattern = "$ext\$")
   sink("$outfile")
   cat("#Tree","\t", "n_leaf_lab", "\t", "n_zero_len_br" , "\t", "n_nodes", "\t", "n_br", "\t", "n_int_br", "\t", "n_ext_br", "\t", "is_binary_tree", "\n")

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
      cat(trees[i], "\t", Ntips, "\t", n_zero_length_branches, "\t", Nnodes, "\t", no_branches, "\t", no_int_branches, "\t", no_ext_branches, "\t", is_binary_tree, "\n")
   }
   sink()
RCMD

check_output $outfile 
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
    
    concat_alignments.pl list2concat > $concat_file
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
   awk 'BEGIN {FS = "|"}{print $1, $2, $3}' $file | perl -pe 'if(/^>/){s/>\S+/>/; s/>\h+/>/; s/\h+/_/g; s/,//g; s/;//g; s/://g; s/\(//g; s/\)//g}' > ${file}ed;
   
   check_output ${file}ed $parent_PID

}
#----------------------------------------------------------------------------------------- 

function make_labels_4_trees()
{
    #label_dir=$1
    printf "${BLUE}# Adding labels back to tree ...${NC}\n"
    grep '>' $(ls ${label_dir}/*fnaedno | head -1) > tree_labels.list
    perl -pe '$c++; s/>/$c\t/; s/\h\[/_[/' tree_labels.list > ed && mv ed tree_labels.list
    add_labels2tree.pl tree_labels.list ${tree_prefix}_nonRecomb_KdeFilt_cdnAlns_FTGTRG.ph &> /dev/null
    check_output ${tree_prefix}_nonRecomb_KdeFilt_cdnAlns_FTGTRG_ed.ph $parent_PID
}
#----------------------------------------------------------------------------------------- 

function get_critical_TajD_values()
{
   no_seq=$1
   TajD_low=$(awk -v n=$no_seq '$1 == n' ${distrodir}/pop_gen_tables/TajD_predicted_CIs.tsv | cut -d' ' -f2)
   TajD_up=$(awk -v n=$no_seq '$1 == n' ${distrodir}/pop_gen_tables/TajD_predicted_CIs.tsv | cut -d' ' -f3)
   echo "$TajD_low $TajD_up"
}
#----------------------------------------------------------------------------------------- 

function get_critical_FuLi_values()
{
   no_seq=$1
   FuLi_low=$(awk -v n=$no_seq 'BEGIN{FS="\t"}$1 == n' ${distrodir}/pop_gen_tables/Fu_Li_predicted_CIs.tsv | cut -f2)
   FuLi_up=$(awk -v n=$no_seq 'BEGIN{FS="\t"}$1 == n'  ${distrodir}/pop_gen_tables/Fu_Li_predicted_CIs.tsv | cut -f3)
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
	[ $DEBUG ] && echo "check_output running: kill -9 $pPID"
	kill -9 $pPID
	
    fi
}
#----------------------------------------------------------------------------------------- 

function print_development_notes()
{
   cat <<DEV
   $progname v.$VERSION usage:
 
   TODO: (last review: May 14th, 2017)
    I. CRITICAL/Important:
    
    0. Need to make sure that the concatenated tree gets concat prefix in file name for compute_suuValStats_and_RF-dist.R to work!
       Rename the final *ed.ph to *ed.tre to avoid problems with compute_suppValStats_and_RF-dist.R to work
       
    0.1 Compile external binaries on Mac OS X (64 bits) to incude in $distrodir/bin/darwin 
    
    1. run_kdetrees.R and compute_suppValStas_and_RF-dist.R run_molecClock_test_with_paup.sh may be refactored into functions run from within $progname
    1.1 verify if kdetrees or RF-dist require haplotyes and/or rooted trees! 
    1.2 rationalize names of output files (*out) and improve graphics; think of using ggplot2 graphics-
    1.3 verify the output generated by runmodes 1 and 2 of compute_suppValStas_and_RF-dist.R; we are calling it with runmode 2 in both occasions
    1.4 add array dependencies lists for perl modules/packages and write install_get_phylomarkers.sh script
        to make install as easy as possible for users; incude required R and Perl libs.
    1.5 Write documentation/Manual    
    1.6. refactor -t PROT and -t DNA codes into subs using appropriate suffixes; make_labels_4_trees() calibrated for DNA!!!
         Make also sure that we use a concat_prefix in the concatenated tree file name, because this is expected by 
	 compute_suuValStats_and_RF-dist.R
    1.7 Need to check the check_output function	and related checking code; add more checking code
    1.8 Need to think of most suitable install procedure, including directory structure
  
    II. DESIRABLE:
    2.1 implement functions to run AU tests and to compute RF distances against the concatenated tree of top markers
    2.2 Think if implementing here the pop_genetics or in a separte script; remove -R if not used! A reasonable place would be
       to run the evaluation pipeline in the non_recomb_cdn_aln/ dir, and run it after the kdetrees test. 
    2.3 finish/test the run_pexec_cmmds() function to run from within the script (minimize the usage of external scripts)
    2.4 Add the functionality of get_TajD_critical_values() and get_FuLi_critical_values() to popGen_summStats.pl	 

    III: CODE CLEANUP
    3.1 Cleanup all companion scripts
    3.2 Rename some scripts to something less pathological, particularly add_labels2tree.pl
        and run_molecClock_test_jmodeltest2_paup.sh

    NOTES:
      1. read the descriptions of code blocks with: grep -A 300 BLOCK $progname | egrep '^#|^[[:space:]]+#'
      
      
    VERSION HISTORY:
      
DEV

exit 0
}
#----------------------------------------------------------------------------------------- 

function print_help()
{
   cat <<EOF
   $progname v.$VERSION usage:
   
   REQUIRED:
    -R <integer> RUNMODE
          1 select optimal markers for phylogenetics/phylogenomics (genomes form different species).
	  2 select optimal markers for population genetics (genomes form the same species).
    -t <string> type of input sequences: DNA|PROT
    
   OPTIONAL:
     -h flag to print this short help notes
     -H flag to print additional usage Notes
     -c <integer> NCBI codontable number (1-23) for pal2nal.pl to generate codon alignment;        [default:$codontable] 
     -C flag to print codontables
     -d flag to print debugging messages                                                          [default: $DEBUG]
     -D flag to print development notes and TODOs                                            
     -e <integer> select gene trees with at least (min. = 4) external branches                     [default: $min_no_ext_branches]
     -k <real> kde stringency (0.7-1.6 are reasonable values; less is more stringent)              [default: $kde_stringency]
     -K <integer> run molecular clock test on codon alignments                                     [default: $eval_clock]
     -l <integer> max. spr length (7-12 are reasonable values)                                     [default: $spr_length]
     -m <real> min. average support value (0.7-0.8 are reasonable values) for trees to be selected [default: $min_supp_val]
     -M <string> base Model for clock test (use one of: GTR|TrN|HKY|K2P|F81); uses +G in all cases [default: $base_mod]
     -q <real> quantile (0.95|0.99) of Chi-square distribution for computing molec. clock p-value  [default: $q]
     -r <string> root method (midpoint|outgroup)                                                   [default: $root_method]
     -s <integer> number of spr rounds (4-20 are reasonable values) for FastTree tree searching    [default: $spr]
     -T <string> tree search Thoroughness: high|medium|low|lowest                                  [default: $search_thoroughness]
     -V <integer> Verbosity level (how much rubish sent to STDOUT)                                 [default: $VERBOSITY]
     
   Invocation examples:
     1. Default: $progname -R 1 -t DNA
     2. thorough searching and molecular clock analysis on DNA sequences:
          $progname -R 1 -t DNA -k 1.3 -m 0.8 -spr 8 -l 10 -p my_Xgenomes -T high -K 1 -M HKY -q 0.95
     3. fastest searching on a huge protein dataset
          $progname -R 1 -t PROT -m 0.7 -k 1.0 -T lowest
     
   NOTES
     1: run from within the directory holding core genome clusters generated by 
          compare_clusters.pl with -t no_genome_gbk files (core genome: all clusters with a single copy per genome)
     
EOF

   exit 2  
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
VERBOSITY=0
DEBUG=
spr_length=8
spr=4
codontable=11 # bacterial by default
base_mod=GTR
eval_clock=0
root_method=midpoint
tree_prefix=concat 
q=0.99


# See bash cookbook 13.1 and 13.2
while getopts ':c:e:k:K:l:m:M:p:q:r:s:t:T:R:V:hHdCD?:' OPTIONS
do
   case $OPTIONS in
   h)   print_help
        ;;
   k)   kde_stringency=$OPTARG
        ;;
   c)   codontable=$OPTARG
        ;;
   C)	print_codontables
	;;
   d)   DEBUG=1
        ;;
   D)	print_development_notes
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
   R)   runmode=$OPTARG
        ;;
   H)   print_usage_notes
        ;;
   V)   VERBOSITY=$OPTARG
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
env_vars=$(set_pipeline_environment) # returns $distrodir $bindir $OS
distrodir=$(echo $env_vars | awk '{print $1}')
bindir=$(echo $env_vars | awk '{print $2}')
OS=$(echo $env_vars | awk '{print $3}')
bindir="$distrodir/bin/$OS"

[ $DEBUG ] && echo "distrodir:$distrodir"
[ $DEBUG ] && echo "bindir:$bindir"

# 0.1 Determine if pipeline scripts are in $PATH; 
# if not, add them

# returns $homebinflag $homebinpathflag; if $homebinpathflag add symlinks to scripts
homebinpath_flags=$(check_scripts_in_path $distrodir) 
homebinflag=$(echo $homebinpath_flags | awk '{print $1}')
homebinpathflag=$(echo $homebinpath_flags | awk '{print $2}')
[ $DEBUG ] && echo "homebinpath_flags:$homebinpath_flags"

# 0.2  Determine if second-party binaries are in $PATH; 
#  if they are not in $PATH then:
# i) will generate a symlink from $HOME/bin (if in path)
# to the binaries provided in $bindir.
# ii) If no $HOME/bin exists, or it is not in $PATH,
# then $bindir will be added to $PATH
  set_bindirs $bindir $homebinpathflag

#-------------------------------------#
# >>>BLOCK 0.2 CHECK USER OPTIONS <<< #
#-------------------------------------#

if [ -z $runmode ]
then
       echo "# ERROR: no runmode defined!"
       print_help
       exit 1    
fi

if [ $min_no_ext_branches -lt 4 ]
then
    printf "${RED}>>> ERROR: -e has to be >= 4\n\n${NC}"
    print_help
    exit 1
fi

if [ "$mol_type" != "DNA" -a "$mol_type" != "PROT" ] # "$mol_type" == "BOTH" not implemented yet
then
     printf "\n${RED}ERROR: -t must be DNA or PROT${NC}\n"
     print_help
     exit 1
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

start_time=$(date +%s)

parent_PID=$(get_script_PID $progname)
[ $DEBUG ] && echo "parent_PID:$parent_PID"

logdir=$(pwd)

dir_suffix=t${mol_type}_k${kde_stringency}_m${min_supp_val}_s${spr}_l${spr_length}_T${search_thoroughness}

printf "
 ${CYAN}>>> $(basename $0) vers. $VERSION run with the following parameters:${NC}
 ${YELLOW}Run time stamp:$TIMESTAMP_SHORT_HMS
 wkdir=$wkdir
 setbindir_flag=$setbindir_flag|bindir=$bindir
 runmode=$runmode|mol_type=$mol_type|eval_clock=$eval_clock|root_method=$root_method|base_model=$base_mod|ChiSq_quantile=$q
 kde_stringency=$kde_stringency|min_supp_val=$min_supp_val|spr=$spr|spr_length=$spr_length|search_thoroughness=$search_thoroughness${NC}

" | tee ${logdir}/get_phylomarkers_run_${dir_suffix}_${TIMESTAMP_SHORT}.log 

#----------------------------------------------------------------------------------------------------------------
#>>>BLOCK 1. make a new subdirectory within the one holding core genome clusters generated by compare_clusters.pl
#    and generate symlinks to the faa and fna files (NOTE: BOTH REQUIRED). Fix fasta file names and headers
#    Mark dir as top_dir
#----------------------------------------------------------------------------------------------------------------

if [ -d get_phylomarkers_run_${dir_suffix}_${TIMESTAMP_SHORT} ]
then
    printf "${RED} >>> Found and older get_phylomarkers_run_${dir_suffix}_${TIMESTAMP_SHORT}/ directory. Please remove or rename and re-run!${NC}\n\n" | \
    tee -a ${logdir}/get_phylomarkers_run_${dir_suffix}_${TIMESTAMP_SHORT}.log
    exit 2
fi 

mkdir get_phylomarkers_run_${dir_suffix}_${TIMESTAMP_SHORT} && cd get_phylomarkers_run_${dir_suffix}_${TIMESTAMP_SHORT}
top_dir=$(pwd)

print_start_time && printf "${BLUE}# processing source fastas in directory get_phylomarkers_run_${dir_suffix}_${TIMESTAMP_SHORT} ...${NC}\n" | \
tee -a ${logdir}/get_phylomarkers_run_${dir_suffix}_${TIMESTAMP_SHORT}.log

ln -s ../*faa .
ln -s ../*fna .

# fix fasta file names with two and three dots
rename 's/\.\.\./\./g' *.faa
rename 's/\.\.\./\./g' *.fna

# 1.1 fix fastaheaders of the source protein and DNA fasta files
for file in *faa; do awk 'BEGIN {FS = "|"}{print $1, $2, $3}' $file | perl -pe 'if(/^>/){s/>\S+/>/; s/>\h+/>/; s/\h+/_/g; s/,//g; s/;//g; s/://g; s/\(//g; s/\)//g}' > ${file}ed; done
for file in *fna; do awk 'BEGIN {FS = "|"}{print $1, $2, $3}' $file | perl -pe 'if(/^>/){s/>\S+/>/; s/>\h+/>/; s/\h+/_/g; s/,//g; s/;//g; s/://g; s/\(//g; s/\)//g}' > ${file}ed; done

# 1.2 add_nos2fasta_header.pl to avoid problems with duplicate labels
run_pexec_cmmds.sh faaed 'add_nos2fasta_header.pl $file > ${file}no' &> /dev/null
run_pexec_cmmds.sh fnaed 'add_nos2fasta_header.pl $file > ${file}no' &> /dev/null

no_alns=$(ls *.fnaedno | wc -l)

[ $no_alns -eq 0 ] && printf "\n${RED} >>> ERROR: There are no codon alignments to work on! Something went wrong. Please check input and settings ...${NC}\n" && exit 4

# 1.3 generate a tree_labels.list file for later tree labeling
print_start_time && printf "${BLUE}# generating the labels file for tree-labeling ...${NC}\n" | tee -a ${logdir}/get_phylomarkers_run_${dir_suffix}_${TIMESTAMP_SHORT}.log
tree_labels_dir=$(pwd)
grep '>' $(ls *fnaedno | head -1) > tree_labels.list
perl -pe '$c++; s/>/$c\t/; s/\h\[/_[/' tree_labels.list > ed && mv ed tree_labels.list

#------------------------------------------------------------------------------------------------
#>>>BLOCK 2. Generate cdnAlns with with pal2nal, maintaining the input-order of the source fastas
#------------------------------------------------------------------------------------------------

# 2.1 generate the protein alignments using clustalo
print_start_time &&  printf "${BLUE}# generating $no_alns codon alignments ...${NC}\n" | tee -a ${logdir}/get_phylomarkers_run_${dir_suffix}_${TIMESTAMP_SHORT}.log
run_pexec_cmmds.sh faaedno 'clustalo -i $file -o ${file%.*}_cluo.faaln --output-order input-order' &> /dev/null

# 2.2 generate the codon alignments (files with *_cdnAln.fasta extension) using pal2nal.pl, 
#     excluding gaps, and mismatched codons, assuming a bacterial genetic code
#for file in *faaln
#do
#   pal2nal.pl $file ${file%_cluo.faaln}.fnaedno -output fasta -nogap -nomismatch -codontable $codontable > ${file%_cluo.faaln}_cdnAln.fasta
#done

# NOTE: to execute run_pexec_cmmds.sh with a customized command, resulting from the interpolation of multiple varialbles,
#       we have to first construct the command line in a variable and pipe its content into bash for execution
faaln_ext=faaln
command="run_pexec_cmmds.sh $faaln_ext 'pal2nal.pl \$file \${file%_cluo.faaln}.fnaedno -output fasta -nogap -nomismatch -codontable $codontable > \${file%_cluo.faaln}_cdnAln.fasta'"

# now we can execute run_pexec_cmmds.sh with a customized command, resulting from the interpolation of multiple varialbles
echo "$command" | bash &> /dev/null

# check we got the expected *cdnAln.fasta files or die!
for f in *cdnAln.fasta
do
     [ ! -s $f ] && printf "\n${RED} >>> ERROR: produced empty alignment! Check your input. Will stop now ...\n\n${NC}" && exit 5
done

# 2.3 cleanup: remove the source fnaed and faaed files; make numbered_fna_files.tgz and numbered_faa_files.tgz; rm *aedno
rm *fnaed *faaed
tar -czf numbered_fna_files.tgz *fnaedno
tar -czf numbered_faa_files.tgz *faaedno
rm *aedno

#----------------------------------------------------------------------------------------------------------#
#>>> BLOCK 3. run Phi-test to identify recombinant codon alignments on all *_cdnAln.fasta source files <<< #
#----------------------------------------------------------------------------------------------------------#
# 3.1 make a new PhiPack subdirectory to work in. generate symlinks to ../*fasta files
#     Mark dir as phipack_dir
mkdir PhiPack && cd PhiPack
phipack_dir=$(pwd)

ln -s ../*fasta .

no_fasta_files=$(ls *.fasta | wc -l)
[ $no_fasta_files -gt 1 ] && print_start_time && printf "${BLUE}# running Phi in PhiPack dir on $no_fasta_files codon alignments ...${NC}\n" | \
tee -a ${logdir}/get_phylomarkers_run_${dir_suffix}_${TIMESTAMP_SHORT}.log
[ $no_fasta_files -lt 1 ] && print_start_time && printf "\n${RED} >>> ERROR: there are no codon alignments to run Phi on. Will exit now!${NC}\n\n" | \
tee -a ${logdir}/get_phylomarkers_run_${dir_suffix}_${TIMESTAMP_SHORT}.log && exit 3

# 3.2 run Phi from the PhiPack in parallel
print_start_time && printf "${BLUE}# running Phi test in PhiPack dir ...${NC}\n"| tee -a ${logdir}/get_phylomarkers_run_${dir_suffix}_${TIMESTAMP_SHORT}.log
run_pexec_cmmds.sh fasta 'Phi -f $file -p 1000 > ${file%.*}_Phi.log' &> /dev/null

# 3.3 process the *_Phi.log files generated by Phi to write a summary table and print short overview to STDOUT
for f in *_Phi.log
do 
   perm=$(grep Permut $f | awk '{print $3}')
   norm=$(grep Normal $f | awk '{print $3}')
   echo -e "$f\t$norm\t$perm"
done > Phi_results_${TIMESTAMP_SHORT}.tsv

check_output Phi_results_${TIMESTAMP_SHORT}.tsv $parent_PID | tee -a ${logdir}/get_phylomarkers_run_${dir_suffix}_${TIMESTAMP_SHORT}.log

no_non_recomb_alns_perm_test=$(awk '$2 > 5e-02 && $3 > 5e-02' Phi_results_${TIMESTAMP_SHORT}.tsv | wc -l)
total_no_cdn_alns=$(ls *_cdnAln.fasta | wc -l)

printf "${GREEN} >>> Phi test result: there are $no_non_recomb_alns_perm_test non-recomb alignments out of $total_no_cdn_alns total alignments${NC}\n" | \
tee -a ${logdir}/get_phylomarkers_run_${dir_suffix}_${TIMESTAMP_SHORT}.log


# 3.4 mv non-recombinant codon alignments and protein alignments to their own directories:
#     non_recomb_cdn_alns/ and  non_recomb_cdn_alns/
#     Mark dir as non_recomb_cdn_alns
mkdir non_recomb_cdn_alns
non_recomb_cdn_alns_dir=$(pwd)

print_start_time && printf "${BLUE}# working in dir non_recomb_cdn_alns ...${NC}\n" | \
tee -a ${logdir}/get_phylomarkers_run_${dir_suffix}_${TIMESTAMP_SHORT}.log

for base in $(awk '$2 > 5e-02 && $3 > 5e-02' Phi_results_${TIMESTAMP_SHORT}.tsv | awk '{print $1}' | sed 's/_Phi\.log//')
do 
   cp ${base}.fasta non_recomb_cdn_alns
done 

mkdir non_recomb_FAA_alns
for base in $(awk '$2 > 5e-02 && $3 > 5e-02' Phi_results_${TIMESTAMP_SHORT}.tsv | awk '{print $1}' | sed 's/_cdnAln_Phi\.log//')
do
   cp ../${base}*.faaln non_recomb_FAA_alns
done

#------------------------------------------------------------------------------------------------
#>>>BLOCK 4. Compute individual ML gene trees from codon or protein alignments 
#            and run the kdetrees test to identify deviant trees
#            and keep the ones passing the test (non-significant tests) for supermatrix analysis.
#------------------------------------------------------------------------------------------------
#
# 4.1 cd into non_recomb_cdn_alns and run FastTree in parallel on all codon alignments, 
#     using exhausitve tree-rearrangements: turn off NNI heuristics, and to always optimize 
#     all 5 branches at each NNI, optimizing all 5 branches in as many as 3 rounds. 
#     This is hardcoded in run_pexec_cmmd because we want optimal gene trees at this stage,
#     which are rapidly computed by FastTree due to the low number of columns in CDS/product
#     alignments. 
#     NOTES: 
#         1. at this point the code is divided into two large if blocks to run -t DNA|PROT|BOTH
#                    [ "$mol_type" == "DNA" ]
#                    [ "$mol_type" == "PROT" ]
#         2. blocks 4 and 5 are highly repetitive; refactor into a subroutine, passing proper suffixes
if [ "$mol_type" == "DNA" ]
then
    cd non_recomb_cdn_alns

    #  IMPORTANT NOTE: the sequences should be collapsed to haplotypes
    #  but beware that the R code below complains when the tree labels 
    #  contain '#' symbols, as introduced by collapse2haplotypes.pl
    #  so use a line like the following to change # by - in the collapsed fastas (or on the trees)

    #printf "${BLUE}# collapsing fastas to haplotypes ...${NC}\n"
    #for file in *fasta; do collapse2haplotypes.pl $file | awk '{print $1, $2}' | perl -pe 'if(/^>/){ s/\h\#\d+// }' > ${file}UNIQ; done &> /dev/null

    #printf "${BLUE}# running FastTree on fastas collapsed to haplotypes ...${NC}\n"
    #run_pexec_cmmds.sh fastaUNIQ 'FastTree -quiet -nt -gtr -gamma -bionj -slownni -mlacc 3 -spr 8 -sprlength 8 < $file > ${file%.*}_haploFTGTRG.ph' &> /dev/null

    print_start_time && printf "${BLUE}# estimating $no_non_recomb_alns_perm_test gene trees from non-recombinant sequences ...${NC}\n" | \
    tee -a ${logdir}/get_phylomarkers_run_${dir_suffix}_${TIMESTAMP_SHORT}.log
    run_pexec_cmmds.sh fasta 'FastTree -quiet -nt -gtr -gamma -bionj -slownni -mlacc 3 -spr 8 -sprlength 8 < $file > ${file%.*}_allFTGTRG.ph' &> /dev/null
    
    #remove trees with < 5 branches
    print_start_time && printf "${BLUE}# counting branches on $no_non_recomb_alns_perm_test gene trees ...${NC}\n" | \
    tee -a ${logdir}/get_phylomarkers_run_${dir_suffix}_${TIMESTAMP_SHORT}.log
    count_tree_branches ph no_tree_branches.list &> /dev/null
   
    check_output no_tree_branches.list $parent_PID | tee -a ${logdir}/get_phylomarkers_run_${dir_suffix}_${TIMESTAMP_SHORT}.log
    
    # remove trees with < 5 external branches (leaves)
    for phy in $(grep -v '^#Tree' no_tree_branches.list | awk -v min_no_ext_branches=$min_no_ext_branches 'BEGIN{FS="\t"; OFS="\t"}$7 < min_no_ext_branches' | cut -f1)
    do
         base=$(echo $phy | sed 's/_allFTlgG\.ph//')
	 print_start_time && printf "${RED} >>> will remove ${base}* because it has < 5 branches!${NC}\n"
	 rm ${base}*
    done

    # 4.1 generate the all_GTRG_trees.tre holding all source trees, which is required by kdetrees
    #     Make a check for the existence of the file to interrupt the pipeline if something has gone wrong
    cat *allFTGTRG.ph > all_GTRG_trees.tre
    check_output all_GTRG_trees.tre $parent_PID | tee -a ${logdir}/get_phylomarkers_run_${dir_suffix}_${TIMESTAMP_SHORT}.log

    # 4.2 run_kdetrees.R at desired stringency 
    print_start_time && printf "${BLUE}# running kde test ...${NC}\n" | tee -a ${logdir}/get_phylomarkers_run_${dir_suffix}_${TIMESTAMP_SHORT}.log
    run_kdetrees.R ph all_GTRG_trees.tre $kde_stringency &> /dev/null
    check_output kde_dfr_file_all_GTRG_trees.tre.tab $parent_PID | tee -a ${logdir}/get_phylomarkers_run_${dir_suffix}_${TIMESTAMP_SHORT}.log

    # 4.3 mv outliers to kde_outliers
    no_kde_outliers=$(grep -c outlier kde_dfr_file_all_GTRG_trees.tre.tab)
    no_kde_ok=$(grep -vc outlier kde_dfr_file_all_GTRG_trees.tre.tab)

    if [ $no_kde_outliers -gt 0 ]
    then
        print_start_time && printf "${BLUE}# making dir kde_outliers/ and moving $no_kde_outliers outlier files into it ...${NC}\n" | \
	tee -a ${logdir}/get_phylomarkers_run_${dir_suffix}_${TIMESTAMP_SHORT}.log
        mkdir kde_outliers
        for f in $(grep outlier kde_dfr_file_all_GTRG_trees.tre.tab | cut -f1 | sed 's/_allFTGTRG.ph//'); do mv ${f}* kde_outliers; done
    else
        print_start_time && printf "${GREEN} >>> there are no kde-test outliers ...${NC}\n" | tee -a ${logdir}/get_phylomarkers_run_${dir_suffix}_${TIMESTAMP_SHORT}.log
    fi

    if [ $no_kde_ok -gt 0 ]
    then

        print_start_time && printf "${BLUE}# making dir kde_ok/ and linking $no_kde_ok selected files into it ...${NC}\n" | tee -a ${logdir}/get_phylomarkers_run_${dir_suffix}_${TIMESTAMP_SHORT}.log
        mkdir kde_ok
        cd kde_ok
        ln -s ../*ph .
   
        print_start_time && printf "${BLUE}# labeling $no_kde_ok gene trees in dir kde_ok/ ...${NC}\n" | tee -a ${logdir}/get_phylomarkers_run_${dir_suffix}_${TIMESTAMP_SHORT}.log
        run_pexec_cmmds.sh ph 'add_labels2tree.pl ../../../tree_labels.list $file' &> /dev/null
    
        # remove symbolic links to cleanup kde_ok/
        for f in $(ls *ph | grep -v '_ed\.ph'); do rm $f; done
    
        cd ..
    else
        print_start_time &&  printf "${RED}# There are $no_kde_ok gene trees producing non-significant kde-test results! Increase the actual -k $kde_stringency value. Will stop here. ...${NC}\n" | \
	 tee -a ${logdir}/get_phylomarkers_run_${dir_suffix}_${TIMESTAMP_SHORT}.log
	 exit 5
    fi
       

    if [ $runmode -eq 1 ]
    then
        # 4.6 compute average bipartition support values for each gene tree
        wkdir=$(pwd) 
        
        print_start_time && printf "${BLUE}# computing tree support values ...${NC}\n" | tee -a ${logdir}/get_phylomarkers_run_${dir_suffix}_${TIMESTAMP_SHORT}.log
        ~/R_code/scripts/compute_suppValStas_and_RF-dist.R $wkdir 1 fasta ph 1 &> /dev/null
    
        print_start_time && printf "${BLUE}# writing summary tables ...${NC}\n" | tee -a ${logdir}/get_phylomarkers_run_${dir_suffix}_${TIMESTAMP_SHORT}.log
        min_supp_val_perc=${min_supp_val#0.}
        no_digits=${#min_supp_val_perc}
        [ $no_digits -eq 1 ] && min_supp_val_perc=${min_supp_val_perc}0
    
        awk -v min_supp_val=$min_supp_val '$2 >= min_supp_val' sorted_aggregated_support_values4loci.tab > sorted_aggregated_support_values4loci_ge${min_supp_val_perc}perc.tab
        check_output sorted_aggregated_support_values4loci_ge${min_supp_val_perc}perc.tab $parent_PID | tee -a ${logdir}/get_phylomarkers_run_${dir_suffix}_${TIMESTAMP_SHORT}.log
    
        no_top_markers=$(wc -l sorted_aggregated_support_values4loci_ge${min_supp_val_perc}perc.tab | cut -d' ' -f1)
        top_markers_dir="top_${no_top_markers}_markers_ge${min_supp_val_perc}perc"
        top_markers_tab=$(ls sorted_aggregated_support_values4loci_ge${min_supp_val_perc}perc.tab)
    
        print_start_time && printf "${BLUE}# making dir $top_markers_dir and moving $no_top_markers top markers into it ...${NC}\n" | \
        tee -a ${logdir}/get_phylomarkers_run_${dir_suffix}_${TIMESTAMP_SHORT}.log
        mkdir $top_markers_dir && cd $top_markers_dir
        top_markers_dir=$(pwd)
        ln -s ../$top_markers_tab .
        for base in $(awk '{print $1}' $top_markers_tab | grep -v loci | sed 's/"//g'); do ln -s ../${base}* .; done


        # 4.7 generate supermatrix (concatenated alignment) 
        print_start_time && printf "${BLUE}# concatenating $no_top_markers top markers into supermatrix ...${NC}\n" | \
	tee -a ${logdir}/get_phylomarkers_run_${dir_suffix}_${TIMESTAMP_SHORT}.log
        concat_alns fasta $parent_PID &> /dev/null

        # 4.8 remove uninformative sites from the concatenated alignment to speed up computation
        print_start_time && printf "${BLUE}# removing uninformative sites from concatenated alignment ...${NC}\n" | \
	tee -a ${logdir}/get_phylomarkers_run_${dir_suffix}_${TIMESTAMP_SHORT}.log
        remove_uninformative_sites_from_aln.pl < concat_cdnAlns.fna > concat_cdnAlns.fnainf

        check_output concat_cdnAlns.fnainf $parent_PID | tee -a ${logdir}/get_phylomarkers_run_${dir_suffix}_${TIMESTAMP_SHORT}.log

        # 4.9 run FasTree under the GTR+G model 
        print_start_time && printf "${BLUE}# running FastTree on the concatenated alignment with $search_thoroughness thoroughness. This may take a while ...${NC}\n" | \
        tee -a ${logdir}/get_phylomarkers_run_${dir_suffix}_${TIMESTAMP_SHORT}.log
    
        if [ "$search_thoroughness" == "high" ]
        then
            FastTree -quiet -nt -gtr -bionj -slownni -gamma -mlacc 3 -spr $spr -sprlength $spr_length < concat_cdnAlns.fnainf > ${tree_prefix}_nonRecomb_KdeFilt_cdnAlns_FTGTRG.ph
        fi

        if [ "$search_thoroughness" == "medium" ]
        then
            FastTree -quiet -nt -gtr -bionj -slownni -gamma -mlacc 2 -spr $spr -sprlength $spr_length < concat_cdnAlns.fnainf > ${tree_prefix}_nonRecomb_KdeFilt_cdnAlns_FTGTRG.ph
        fi

        if [ "$search_thoroughness" == "low" ]
        then
            FastTree -quiet -nt -gtr -bionj -gamma -spr $spr -sprlength $spr_length < concat_cdnAlns.fnainf > ${tree_prefix}_nonRecomb_KdeFilt_cdnAlns_FTGTRG.ph
        fi

        if [ "$search_thoroughness" == "lowest" ]
        then
            FastTree -quiet -nt -gtr -gamma -mlnni 4 < concat_cdnAlns.fnainf > ${tree_prefix}_nonRecomb_KdeFilt_cdnAlns_FTGTRG.ph
        fi

        check_output ${tree_prefix}_nonRecomb_KdeFilt_cdnAlns_FTGTRG.ph $parent_PID | tee -a ${logdir}/get_phylomarkers_run_${dir_suffix}_${TIMESTAMP_SHORT}.log
    
        print_start_time && printf "${BLUE}# Adding labels back to tree ...${NC}\n" | tee -a ${logdir}/get_phylomarkers_run_${dir_suffix}_${TIMESTAMP_SHORT}.log
        add_labels2tree.pl ${tree_labels_dir}/tree_labels.list ${tree_prefix}_nonRecomb_KdeFilt_cdnAlns_FTGTRG.ph &> /dev/null
    
        if [ -s ${tree_prefix}_nonRecomb_KdeFilt_cdnAlns_FTGTRG_ed.ph ]
        then
            mv ${tree_prefix}_nonRecomb_KdeFilt_cdnAlns_FTGTRG_ed.ph ${tree_prefix}_nonRecomb_KdeFilt_cdnAlns_FTGTRG_ed.sptree # for compute_suppValStats_and_RF-dist.R
            check_output ${tree_prefix}_nonRecomb_KdeFilt_cdnAlns_FTGTRG_ed.sptree $parent_PID | tee -a ${logdir}/get_phylomarkers_run_${dir_suffix}_${TIMESTAMP_SHORT}.log
            printf "${GREEN} >>> found in dir $top_markers_dir ...${NC}\n\n" | tee -a ${logdir}/get_phylomarkers_run_${dir_suffix}_${TIMESTAMP_SHORT}.log
        else
              printf "${RED} >>> WARNING: ${tree_prefix}_nonRecomb_KdeFilt_cdnAlns_FTGTRG_ed.sptree could not be produced!${NC}"
        fi 
    
    
        print_start_time && printf "${BLUE}# computing the mean support values and RF-distances of each gene tree to the concatenated tree   ...${NC}\n" | \
        tee -a ${logdir}/get_phylomarkers_run_${dir_suffix}_${TIMESTAMP_SHORT}.log
        compute_suppValStas_and_RF-dist.R $top_markers_dir 2 fasta ph 1 &> /dev/null
    
        # NOTE: after v0.9 this process is prallelized with run_pexec_cmmds.sh
	if [ $eval_clock -gt 0 ]
        then 
 	     # 1. convert fasta2nexus
             print_start_time && printf "${BLUE}# converting fasta files to nexus files${NC}\n" | \
       	     tee -a ${logdir}/get_phylomarkers_run_${dir_suffix}_${TIMESTAMP_SHORT}.log
             convert_aln_format_batch_bp.pl fasta fasta nexus nex &> /dev/null
	     
	     print_start_time && printf "${BLUE}# Will test the molecular clock hypothesis for $no_top_markers top markers. This will take some time ...${NC}\n" | \
       	     tee -a ${logdir}/get_phylomarkers_run_${dir_suffix}_${TIMESTAMP_SHORT}.log
             #run_molecClock_test_jmodeltest2_paup.sh -R 1 -M $base_mod -t ph -e fasta -b molec_clock -q $q &> /dev/null
	     
	     
	    # 2. >>> print table header and append results to it
            no_dot_q=$(echo $q | sed 's/\.//' )
            results_table=mol_clock_M${base_mod}G_r${rooting_method}_o${outgroup_OTU_nomol}_q${no_dot_q}_ClockTest.tab
             
            echo -e "#nexfile\tlnL_unconstr\tlnL_clock\tLRT\tX2_crit_val\tdf\tp-val\tmol_clock" > $results_table
	    
	     cmd="run_pexec_cmmds.sh nex 'run_parallel_molecClock_test_with_paup.sh -R 1 -f \$file -M $base_mod -t ph -b global_mol_clock -q $q'"
	     [ $DEBUG ] && echo "run_parallel_molecClock.cmd: $cmd"
	     echo $cmd | bash &> /dev/null
	     
	     mol_clock_tab=$(ls *_ClockTest.tab)

       	     if [ -s $mol_clock_tab ] 
	     then
	        printf "${GREEN} >>> the molecular clock test file can be found in:\n$top_markers_dir/$mol_clock_tab ...${NC}\n\n" | \
		tee -a ${logdir}/get_phylomarkers_run_${dir_suffix}_${TIMESTAMP_SHORT}.log
	    
               # Paste filoinfo and clocklikeness tables
               cut -f1 gene_trees2_concat_tree_RF_distances.tab | grep -v loci | sed 's/"//; s/"$/_/'  > list2grep.tmp
               head -1 "$mol_clock_tab" > header.tmp
           
	       # sort lines in molecular clock output file according to order in gene_trees2_concat_tree_RF_distances.tab
	       while read line
	       do
	           grep "$line" $mol_clock_tab 
	       done < list2grep.tmp >> ${mol_clock_tab}sorted
	   
	       cat header.tmp ${mol_clock_tab}sorted > ed && mv ed ${mol_clock_tab}sorted
           
               paste gene_trees2_concat_tree_RF_distances.tab ${mol_clock_tab}sorted > phylogenetic_attributes_of_top${no_top_markers}_gene_trees.tab
	   
	       check_output phylogenetic_attributes_of_top${no_top_markers}_gene_trees.tab $parent_PID | tee -a ${logdir}/get_phylomarkers_run_${dir_suffix}_${TIMESTAMP_SHORT}.log
           
	       printf "${GREEN} >>> Top markers and associated stats are found in:\n$top_markers_dir ...${NC}\n\n" | \
	       tee -a ${logdir}/get_phylomarkers_run_${dir_suffix}_${TIMESTAMP_SHORT}.log
            else
	          printf "${RED} >>> ${mol_clock_tab} not found"
	    fi

	 # 4. cleanup
         tar -czf molClock_PAUP_files.tgz *_paup.block *.nex  *_clockTest.log *tre *clock.scores
         [ -s molClock_PAUP_files.tgz ] && rm *_paup.block *.nex  *_clockTest.log *tre *clock.scores
         rm list2concat Rplots.pdf
        fi
    fi # if [ $runmode -eq 1 ]; then run phylo pipeline on DNA seqs
    
    if [ $runmode -eq 2 ]
    then
        mkdir popGen && cd popGen
        popGen_dir=$(pwd)

	ln -s ../*fasta .
	no_top_markers=$(ls *fasta | wc -l)
	tmpf=$(ls -1 *fasta | head -1)
	no_seqs=$(grep -c '>' $tmpf)
	[ $DEBUG -eq 1 ] && echo "no_seqs:$no_seqs"

        print_start_time && printf "${BLUE}# Will run descriptive DNA polymorphism statistics for $no_top_markers top markers. This will take some time ...${NC}\n" | \
       	tee -a ${logdir}/get_phylomarkers_run_${dir_suffix}_${TIMESTAMP_SHORT}.log
	
	TajD_crit_vals=$(get_critical_TajD_values $no_seqs)
	TajD_l=$(echo $TajD_crit_vals | awk '{print $1}')
	TajD_u=$(echo $TajD_crit_vals | awk '{print $2}')

	FuLi_crit_vals=$(get_critical_FuLi_values $no_seqs)
	FuLi_l=$(echo "$FuLi_crit_vals" | awk '{print $1}')
	FuLi_u=$(echo "$FuLi_crit_vals" | awk '{print $2}')
	
	[ $DEBUG -eq 1 ] && echo "TajD_crit_vals:$TajD_crit_vals | TajD_l:$TajD_l | TajD_u:$TajD_u | FuLi_crit_vals:$FuLi_crit_vals | FuLi_l:$FuLi_l | FuLi_u:$FuLi_u"
	
        print_start_time && printf "${BLUE}# converting $no_top_markers fasta files to nexus format ...${NC}\n" | \
       	tee -a ${logdir}/get_phylomarkers_run_${dir_suffix}_${TIMESTAMP_SHORT}.log
	convert_aln_format_batch_bp.pl fasta fasta nexus nex &> /dev/null 
	  
        print_start_time && printf "${BLUE}# Running popGen_summStats.pl -R 2 -n nex -f fasta -F fasta -H -r 100 -t $TajD_l -T $TajD_u -s $FuLi_l -S $FuLi_u &> popGen_summStats_hs100.log ...${NC}\n" | \
       	tee -a ${logdir}/get_phylomarkers_run_${dir_suffix}_${TIMESTAMP_SHORT}.log
	popGen_summStats.pl -R 2 -n nex -f fasta -F fasta -H -r 100 -t $TajD_l -T $TajD_u -s $FuLi_l -S $FuLi_u &> popGen_summStats_hs100.log
	
	check_output polymorphism_descript_stats.tab $parent_PID
	
	printf "${GREEN} >>> descriptive DNA polymorphism stats are found in:\n$popGen_dir ...${NC}\n\n" | \
	tee -a ${logdir}/get_phylomarkers_run_${dir_suffix}_${TIMESTAMP_SHORT}.log
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
    
    print_start_time && printf "${BLUE}# estimating $no_non_recomb_alns_perm_test gene trees from non-recombinant sequences ...${NC}\n" | \
    tee -a ${logdir}/get_phylomarkers_run_${dir_suffix}_${TIMESTAMP_SHORT}.log
    run_pexec_cmmds.sh faaln 'FastTree -quiet -lg -gamma -bionj -slownni -mlacc 3 -spr 8 -sprlength 8 < $file > ${file%.*}_allFTlgG.ph' &> /dev/null
    
    #remove trees with < 5 branches
    print_start_time && printf "${BLUE}# counting branches on $no_non_recomb_alns_perm_test gene trees ...${NC}\n" | \
    tee -a ${logdir}/get_phylomarkers_run_${dir_suffix}_${TIMESTAMP_SHORT}.log
    count_tree_branches ph no_tree_branches.list &> /dev/null
   
    check_output no_tree_branches.list $parent_PID | tee -a ${logdir}/get_phylomarkers_run_${dir_suffix}_${TIMESTAMP_SHORT}.log
    
    # remove trees with < 5 external branches (leaves)
    for phy in $(grep -v '^#Tree' no_tree_branches.list | awk -v min_no_ext_branches=$min_no_ext_branches 'BEGIN{FS="\t"; OFS="\t"}$7 < min_no_ext_branches' | cut -f1)
    do
         base=$(echo $phy | sed 's/_allFTlgG\.ph//')
	 printf "${RED} >>> will remove ${base}* because it has < 5 branches!${NC}\n"
	 rm ${base}*
    done

    # 5.1 generate the all_GTRG_trees.tre holding all source trees. This is required for run_kdetrees.R ph
    #     Make a check for the existence of the file to interrupt the pipeline if something has gone wrong
    cat *_allFTlgG.ph > all_FTlgG_trees.tre
    check_output all_FTlgG_trees.tre $parent_PID | tee -a ${logdir}/get_phylomarkers_run_${dir_suffix}_${TIMESTAMP_SHORT}.log

    # 5.2 run_kdetrees.R at desired stringency 
    print_start_time && printf "${BLUE}# running kde test ...${NC}\n" | tee -a ${logdir}/get_phylomarkers_run_${dir_suffix}_${TIMESTAMP_SHORT}.log
    run_kdetrees.R ph all_FTlgG_trees.tre $kde_stringency &> /dev/null
    check_output kde_dfr_file_all_FTlgG_trees.tre.tab $parent_PID | tee -a ${logdir}/get_phylomarkers_run_${dir_suffix}_${TIMESTAMP_SHORT}.log

    # 5.3 mv outliers to kde_outliers
    no_kde_outliers=$(grep -c outlier kde_dfr_file_all_FTlgG_trees.tre.tab)
    no_kde_ok=$(grep -vc outlier kde_dfr_file_all_FTlgG_trees.tre.tab)

    if [ $no_kde_outliers -gt 0 ]
    then
        print_start_time && printf "${BLUE}# making dir kde_outliers/ and moving $no_kde_outliers outlier files into it ...${NC}\n" | \
	tee -a ${logdir}/get_phylomarkers_run_${dir_suffix}_${TIMESTAMP_SHORT}.log
        mkdir kde_outliers
        for f in $(grep outlier kde_dfr_file_all_FTlgG_trees.tre.tab | cut -f1 | sed 's/_allFTlgG.ph//'); do mv ${f}* kde_outliers; done
    else
        printf "${GREEN} >>> there are no kde-test outliers ...${NC}\n" | tee -a ${logdir}/get_phylomarkers_run_${dir_suffix}_${TIMESTAMP_SHORT}.log
    fi

    if [ $no_kde_ok -gt 0 ]
    then
        print_start_time && printf "${BLUE}# making dir kde_ok/ and linking $no_kde_ok selected files into it ...${NC}\n" | tee -a ${logdir}/get_phylomarkers_run_${dir_suffix}_${TIMESTAMP_SHORT}.log
        mkdir kde_ok
        cd kde_ok
        ln -s ../*ph .
   
        print_start_time && printf "${BLUE}# labeling $no_kde_ok gene trees in dir kde_ok/ ...${NC}\n" | tee -a ${logdir}/get_phylomarkers_run_${dir_suffix}_${TIMESTAMP_SHORT}.log
        run_pexec_cmmds.sh ph 'faalnheader2treetag_PhyML_Consense_Topol_V03.pl ../../../tree_labels.list $file' &> /dev/null
    
        # remove symbolic links to cleanup kde_ok/
        for f in $(ls *ph | grep -v '_ed\.ph'); do rm $f; done
    
        cd ..
    else
         printf "${RED}# There are $no_kde_ok gene trees producing non-significant kde-test results! Increase the actual -k $kde_stringency value. Will stop here. ...${NC}\n" | \
	 tee -a ${logdir}/get_phylomarkers_run_${dir_suffix}_${TIMESTAMP_SHORT}.log
	 exit 5
    fi
    
    # 5.4 compute MJRE consensus tree with consense; does not work with trees with different number of leaves
    # cat *ph > intree
    # echo "y" > consense.cmd

    # consense < consense.cmd
    # check_output outtree $parent_PID
    
    #  4.5 concatenated prefix required by run compute_suppValStas_and_RF-dist.R
    #  ln -s outtree concatenated.ph 

    # 5.6 compute average bipartition support values for each gene tree
    #     and their RF-fistance to the consensus tree computed with consense
    wkdir=$(pwd) 
        
    print_start_time && printf "${BLUE}# computing tree support values ...${NC}\n" | tee -a ${logdir}/get_phylomarkers_run_${dir_suffix}_${TIMESTAMP_SHORT}.log
    ~/R_code/scripts/compute_suppValStas_and_RF-dist.R $wkdir 1 faaln ph 1 &> /dev/null
    
    print_start_time && printf "${BLUE}# writing summary tables ...${NC}\n" | tee -a ${logdir}/get_phylomarkers_run_${dir_suffix}_${TIMESTAMP_SHORT}.log
    min_supp_val_perc=${min_supp_val#0.}
    no_digits=${#min_supp_val_perc}
    [ $no_digits -eq 1 ] && min_supp_val_perc=${min_supp_val_perc}0
    
    awk -v min_supp_val=$min_supp_val '$2 >= min_supp_val' sorted_aggregated_support_values4loci.tab > sorted_aggregated_support_values4loci_ge${min_supp_val_perc}perc.tab
    check_output sorted_aggregated_support_values4loci_ge${min_supp_val_perc}perc.tab $parent_PID | tee -a ${logdir}/get_phylomarkers_run_${dir_suffix}_${TIMESTAMP_SHORT}.log
    
    no_top_markers=$(wc -l sorted_aggregated_support_values4loci_ge${min_supp_val_perc}perc.tab | cut -d' ' -f1)
    top_markers_dir="top_${no_top_markers}_markers_ge${min_supp_val_perc}perc"
    top_markers_tab=$(ls sorted_aggregated_support_values4loci_ge${min_supp_val_perc}perc.tab)
    
    print_start_time && printf "${BLUE}# making dir $top_markers_dir and moving $no_top_markers top markers into it ...${NC}\n" | \
    tee -a ${logdir}/get_phylomarkers_run_${dir_suffix}_${TIMESTAMP_SHORT}.log
    mkdir $top_markers_dir && cd $top_markers_dir
    ln -s ../$top_markers_tab .
    for base in $(awk '{print $1}' $top_markers_tab | grep -v loci | sed 's/"//g'); do ln -s ../${base}* .; done

    # 5.7 generate supermatrix (concatenated alignment) 
    print_start_time && printf "${BLUE}# concatenating $no_top_markers top markers into supermatrix ...${NC}\n" | tee -a ${logdir}/get_phylomarkers_run_${dir_suffix}_${TIMESTAMP_SHORT}.log
    concat_alns faaln $parent_PID &> /dev/null

    # 5.8 remove uninformative sites from the concatenated alignment to speed up computation
        print_start_time && printf "${BLUE}# removing uninformative sites from concatenated alignment ...${NC}\n" | tee -a ${logdir}/get_phylomarkers_run_${dir_suffix}_${TIMESTAMP_SHORT}.log
    remove_uninformative_sites_from_aln.pl < concat_protAlns.faa > concat_protAlns.faainf

    check_output concat_protAlns.faainf $parent_PID | tee -a ${logdir}/get_phylomarkers_run_${dir_suffix}_${TIMESTAMP_SHORT}.log

    # 5.9 run FasTree under the GTR+G model 
    print_start_time && printf "${BLUE}# running FastTree on the concatenated alignment with $search_thoroughness thoroughness. This may take a while ...${NC}\n" | \
    tee -a ${logdir}/get_phylomarkers_run_${dir_suffix}_${TIMESTAMP_SHORT}.log
    
    if [ "$search_thoroughness" == "high" ]
    then
        FastTree -quiet -lg -bionj -slownni -gamma -mlacc 3 -spr $spr -sprlength $spr_length < concat_protAlns.faainf > ${tree_prefix}_nonRecomb_KdeFilt_protAlns_FTlgG.ph
    fi

    if [ "$search_thoroughness" == "medium" ]
    then
        FastTree -quiet -lg -bionj -slownni -gamma -mlacc 2 -spr $spr -sprlength $spr_length < concat_protAlns.faainf > ${tree_prefix}_nonRecomb_KdeFilt_protAlns_FTlgG.ph
    fi

    if [ "$search_thoroughness" == "low" ]
    then
        FastTree -quiet -lg -bionj -gamma -spr $spr -sprlength $spr_length < concat_protAlns.faainf > ${tree_prefix}_nonRecomb_KdeFilt_protAlns_FTlgG.ph
    fi

    if [ "$search_thoroughness" == "lowest" ]
    then
        FastTree -quiet -lg -gamma -mlnni 4 < concat_protAlns.faainf > ${tree_prefix}_nonRecomb_KdeFilt_protAlns_FTlgG.ph
    fi

    check_output ${tree_prefix}_nonRecomb_KdeFilt_protAlns_FTlgG.ph $parent_PID | tee -a ${logdir}/get_phylomarkers_run_${dir_suffix}_${TIMESTAMP_SHORT}.log
    
    print_start_time && printf "${BLUE}# Adding labels back to tree ...${NC}\n" | tee -a ${logdir}/get_phylomarkers_run_${dir_suffix}_${TIMESTAMP_SHORT}.log
    add_labels2tree.pl ${tree_labels_dir}/tree_labels.list ${tree_prefix}_nonRecomb_KdeFilt_protAlns_FTlgG.ph &> /dev/null
    
    check_output ${tree_prefix}_nonRecomb_KdeFilt_protAlns_FTlgG_ed.ph $parent_PID | tee -a ${logdir}/get_phylomarkers_run_${dir_suffix}_${TIMESTAMP_SHORT}.log
     
    printf "${GREEN} >>> found in dir $wkdir ...${NC}\n\n" | tee -a ${logdir}/get_phylomarkers_run_${dir_suffix}_${TIMESTAMP_SHORT}.log
    
    check_output ${logdir}/get_phylomarkers_run_${dir_suffix}_${TIMESTAMP_SHORT}.log $parent_PID | tee -a ${logdir}/get_phylomarkers_run_${dir_suffix}_${TIMESTAMP_SHORT}.log
    
    wkdir=$(pwd)
    compute_suppValStas_and_RF-dist.R $wkdir 2 faaln ph 1 &> /dev/null
fi


# compute the elapsed time since the script was fired
end_time=`date +%s`
secs=$(($end_time-$start_time))

#printf '%dh:%dm:%ds\n' $(($secs/3600)) $(($secs%3600/60)) $(($secs%60))
printf "\n${LBLUE} >>> Total runtime of $progname: ${NC}" 
printf '%dh:%dm:%ds\n' $(($secs/3600)) $(($secs%3600/60)) $(($secs%60))
echo
