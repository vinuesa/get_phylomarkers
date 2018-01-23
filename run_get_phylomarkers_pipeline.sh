#!/usr/bin/env bash

#: PROGRAM: run_get_phylomarkers_pipeline.sh
#: AUTHORS: Pablo Vinuesa, Center for Genomic Sciences, UNAM, Mexico
#:          http://www.ccg.unam.mx/~vinuesa/
#           Bruno Contreras Moreira, EEAD-CSIC, Zaragoza, Spain
#           https://digital.csic.es/cris/rp/rp02661/
#
#: DISCLAIMER: programs of the GET_PHYLOMARKERS package are distributed in the hope that it will be useful, 
#              but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
#              See the GNU General Public License for more details. 
#
#: LICENSE: This software is freely available under the GNU GENERAL PUBLIC LICENSE v.3.0
#           see https://github.com/vinuesa/get_phylomarkers/blob/master/LICENSE
#
#: AVAILABILITY: freely available from GitHub @ https://github.com/vinuesa/get_phylomarkers
#
#: PROJECT START: April 2017; This is a wrapper script to automate the whole process of marker selection and downstream analyses.
#
#: AIM: select optimal molecular markers for phylogenomics and population genomics from orthologous gene clusters computed by GET_HOMOLOGUES
#           which is freely available from GitHub @ https://github.com/eead-csic-compbio/get_homologues
#
#: OUTPUT: multiple sequence alignments (of protein and DNA sequences) of selected markers, gene trees and species tree 
#              inferred from concatenated supermatrix of top-ranking markers, 
#              along with graphics and tables summarizing the results of the pipeline obtained at different levels.
#
progname=${0##*/} # run_get_phylomarkers_pipeline.sh
VERSION='2.0.1_22Jan18'

# Set GLOBALS
DEBUG=0
wkdir=$(pwd) #echo "# working in $wkdir"

dir_suffix=
gene_tree_ext=
sp_tree_ext=

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
#printf "${RED}%s${NC} ${GREEN}%s${NC}\n" "I like" "GitHub"


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
  [ $DEBUG -eq 1 ] && msg " => working in $FUNCNAME ..." DEBUG NC
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
  [ $DEBUG -eq 1 ] && msg " <= exiting $FUNCNAME ..." DEBUG NC
}
#-----------------------------------------------------------------------------------------

function check_dependencies()
{

    # check if scripts are in path; if not, set flag
    [ $DEBUG -eq 1 ] && msg " => working in $FUNCNAME ..." DEBUG NC
    for prog in bash R perl awk bc cut grep sed sort uniq Rscript
    do
       bin=$(type -P $prog)
       if [ -z $bin ]; then
          echo
          printf "${RED}%s${NC}\n"  "# ERROR: system binary $prog is not in \$PATH!"
	  printf "${LRED}%s${NC}\n" " >>> you will need to install $prog for $progname to run on $HOSTNAME"
	  printf "${LRED}%s${NC}\n" " >>> will exit now ..."
	  exit 1
       fi
    done
}
#-----------------------------------------------------------------------------------------

function check_scripts_in_path()
{
    [ $DEBUG -eq 1 ] && msg " => working in $FUNCNAME ..." DEBUG NC
    distrodir=$1
    not_in_path=0
    homebinflag=0
    homebinpathflag=0

    [ $DEGUB ] && msg "check_scripts_in_path() distrodir:$distrodir" DEBUG NC

    bash_scripts=(run_parallel_molecClock_test_with_paup.sh )
    perl_scripts=( run_parallel_cmmds.pl add_nos2fasta_header.pl pal2nal.pl rename.pl concat_alignments.pl add_labels2tree.pl convert_aln_format_batch_bp.pl popGen_summStats.pl convert_aln_format_batch_bp.pl )
    R_scripts=( run_kdetrees.R compute_suppValStas_and_RF-dist.R )

    # check if scripts are in path; if not, set flag
    for prog in "${bash_scripts[@]}" "${perl_scripts[@]}" "${R_scripts[@]}"
    do
       bin=$(type -P $prog)
       if [ -z $bin ]; then
          echo
          msg "# WARNING: script $prog is not in \$PATH!" WARNING LRED
	        msg " >>>  Will generate a symlink from $HOME/bin or add it to \$PATH" WARNING CYAN
	        not_in_path=1
       fi
    done

    # if flag $not_in_path -eq 1, then either generate symlinks into $HOME/bin (if in $PATH) or export $distrodir to PATH
    if [ $not_in_path -eq 1 ]
    then
       if [ ! -d $HOME/bin ]
       then
          msg "Could not find a $HOME/bin directory for $USER ..."  WARNING CYAN
	  msg " ... will update PATH=$distrodir:$PATH"  WARNING CYAN
	  export PATH="${distrodir}:${PATH}" # prepend $ditrodir to $PATH
       else
           homebinflag=1
       fi

       # check if $HOME/bin is in $PATH
       echo $PATH | sed 's/:/\n/g'| grep "$HOME/bin$" &> /dev/null
       if [ $? -eq 0 ]
       then
          homebinpathflag=1

          msg "Found dir $HOME/bin for $USER in \$PATH ..." WARNING CYAN
          msg " ... will generate symlinks in $HOME/bin to all scripts in $distrodir ..." WARNING CYAN
          ln -s $distrodir/*.sh $HOME/bin &> /dev/null
          ln -s $distrodir/*.R $HOME/bin &> /dev/null
          ln -s $distrodir/*.pl $HOME/bin &> /dev/null
          #ln -s $distrodir/rename.pl $HOME/bin &> /dev/null
       else
          msg " Found dir $HOME/bin for $USER, but it is NOT in \$PATH ..." WARNING CYAN
          msg " ... updating PATH=$PATH:$distrodir" WARNING CYAN
	  export PATH="${distrodir}:${PATH}" # prepend $distrodir to $PATH
       fi
    fi
    #echo "$homebinflag $homebinpathflag"
    [ $DEBUG -eq 1 ] && msg " <= exiting $FUNCNAME ..." DEBUG NC
}
#-----------------------------------------------------------------------------------------

function set_bindirs()
{
    [ $DEBUG -eq 1 ] && msg " => working in $FUNCNAME ..." DEBUG NC
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
   [ $DEBUG -eq 1 ] && msg " <= exiting $FUNCNAME ..." DEBUG NC
}
#-----------------------------------------------------------------------------------------

function print_usage_notes()
{
   cat <<USAGE

   $progname v$VERSION extensive Help and details on the search modes and models.

   1. Start the run from within the directory holding core gene clusters generated by 
      get_homologues.pl -e -t number_of_genomes or compare_clusters.pl -t number_of_genomes

      NOTE: Both .faa and .fna files are required to generate codon alignments from DNA fasta files. This
        means that two runs of compare_clusters.pl (from the get_homologues package) are required, one of them
        using the -n flag. See GET_HOMOLOGUES online help http://eead-csic-compbio.github.io/get_homologues/manual/ 

   2. $progname is intended to run on a collection of single-copy sequence clusters from different species or strains.

      NOTES: An absolute minimum of 4 distinct genomes are required.
       However, the power of the pipeline for selecting optimal genome loci
	  for phylogenomics improves when a larger number of genomes are available
	  for analysis. Reasonable numbers lie in the range of 10 to 200 distinct genomes
	  from multiple species of a genus, family, order or phylum.
	  The pipeline may not perform satisfactorily with very distant genome sequences,
	  particularly when sequences with significantly distinct nucleotide or aminoacid
	  compositions are used. This type of sequence heterogeneity is well known to
	  cause systematic bias in phylogenetic inference.

   3. On the locus filtering criteria. $progname uses a hierarchical filtering scheme, as follows:

      i) Detection of recombinant loci. Codon or protein alignments (depending on runmode)
          are first screened with phi(w) for the presence of potential recombinant sequences.
	  It is a well established fact that recombinant sequences negatively impact
	  phylogenetic inference when using algorithms that do not account for the effects
	  of this evolutionary force. The permutation test with 1000 permutations is used
	  to compute the p-values. These are considerd significant if < 0.05.

      ii) Detection of trees deviating from the expectation of the (multispecies) coalescent.
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

      iii) Phylogenetic signal content. 
            The alignments passing the two previous filters are subjected to maximum likelihood 
	    (ML) tree searches with FastTree or IQ-TREE to infer the corresponding ML gene trees. 
	    The phylogenetic signal of these trees is computed from the Shimodair-Hasegawa-like 
	    likelihood ratio test (SH-alrt) of branch support values, which vary between 0-1. 
	    The support values of each internal branch or bipartition are parsed to compute the 
	    mean support value for each tree. Trees with a mean support value below a cutoff 
	    threshold are discarded.

         * Parameters controlling filtering based on mean support values.
         -m <real> min. average support value (0.7-0.8 are reasonable values)
	           for trees to be selected [default: $min_supp_val]

      iv) On tree searching: From version 2.0 onwards, $progname performs tree searches using
          either the FastTree (FT) or IQ-TREE (IQT) fast ML tree search algorithms,
	  controlled with the -A <F|I> option [default: $search_algorithm]

       a) FT searches:
	  FT meets a good compromise between speed and accuracy, runnig both
	  with DNA and protein sequence alignments. It computes the above-mentioned
	  Shimodaria-Hasegawa-like likelihood ratio test of branch support values.
	  A limitation though, is that it implements only tow substitution models.
	  However, for divergent sequences of different species within a bacterial taxonomic
	  genus or family, our experience has shown that almost invariably the GTR+G model
	  is selected by jmodeltest2, particularly when there is base frequency heterogeneity.
	  The GTR+G+CAT is the substitution model used by $progname calls of FastTree on codon
	  alignments. From version 2.0 onwards, $progname can compute gene trees wit FT using
	  differetn levels of tree search intensity as defined with the -T parameter, both
	  during the gene-tree and species-tree estimation phases.

	   high:   -nt -gtr -bionj -slow -slownni -gamma -mlacc 3 -spr $spr -sprlength $spr_length
	   medium: -nt -gtr -bionj -slownni -gamma -mlacc 2 -spr $spr -sprlength $spr_length
	   low:    -nt -gtr -bionj -gamma -spr $spr -sprlength $spr_length
           lowest: -nt -gtr -gamma -mlnni 4

	   where -s \$spr and -l \$spr_length can be set by the user.
	   The lines above show their default values.

	   The same applies for tree searches on concatenated codon alignments.
	   Note that these may take a considerable time (up to several hours)
	   for large datasets (~ 100 taxa and > 300 concatenated genes).

	   For protein alignments, the search parameters are the same, only the model changes to lg

	   Please refer to the FastTree manual for the details.

       b) IQT searches:
          Our benchmark analyses have shown that IQ-TREE (v1.6.1) runs quickly enough when the '-fast'
	  flag is passed to make it feasible to include a modelselection step withouth incurring in 
	  prohibitively long computation times. Combined with its superior tree-searching algorithm, 
	  makes IQT the clear winner in our benchmarks. Therefore, from version 2.0 (22Jan18) onwards,
	  $progname uses IQT as its default tree searching algorithm, both for the estimation of
	  gene-trees and the species-tree (from the concatenated, top-scoring alignments), now using
	  model selection in both cases (v 1.9 only used model selection and IQT for supermatrix search).
	  
	  However, the number of models evaluated by ModelFinder (integrated in IQ-TREE) differ for the
	  gene-tree and species-tree search phases, as detailed below: 

	-IQT gene-tree searches (hard-coded): -T <high|medium|low|lowest> [default: $search_thoroughness]
	  high:   -m MFP -nt 1 -alrt 1000 -fast [ as jModelTest]
	  medium: -mset K2P,HKY,TN,TNe,TIM,TIMe,TIM2,TIM2e,TIM3,TIM3e,TVM,TVMe,GTR -m MFP -nt 1 -alrt 1000 -fast
	  low:	  -mset K2P,HKY,TN,TNe,TVM,TVMe,TIM,TIMe,GTR -m MFP -nt 1 -alrt 1000 -fast
	  lowest: -mset K2P,HKY,TN,TNe,TVM,TIM,GTR -m MFP -nt 1 -alrt 1000 -fast

	  all gene trees are run in parallel under the modelset with the following parameters: 
	    -mset X,Y,Z -m MFP -nt 1 -alrt 1000 -fast

	-IQT species-tree searches on the supermatrix: 
	      -S <string> comma-separated list of base models to be evaluated by IQ-TREE
	         when estimating the species tree from the concatenated supermatrix.
		 If no -S is passed, then sinlge default models are used, as shown below
              <'JC,F81,K2P,HKY,TrN,TNe,K3P,K81u,TPM2,TPM2u,TPM3,TPM3u,
	      TIM,TIMe,TIM2,TIM2e,TIM3,TIM3e,TVM,TVMe,SYM,GTR'>              for DNA alignments    [default: $IQT_DNA_models]
              <'BLOSUM62,cpREV,Dayhoff,DCMut,FLU,HIVb,HIVw,JTT,JTTDCMut,LG,
                mtART,mtMAM,mtREV,mtZOA,Poisson,PMB,rtREV,VT,WAG'>           for PROT alignments   [default: $IQT_PROT_models]
		
          In addition, if -T high, $progname will launch -N <integer> [default: $nrep_IQT_searches] independent IQT searches
	     on the supermatrix of concatenated top-scoring markers.

    v) Running the pipeline in population-genetics mode (-R 2 -t DNA): 
       When invoked in popGen mode (-R 2), the pipeline will perform the same 4 initial steps as in phylo mode (-R 1): 
          1. generate codon alginments
          2. chech them for the presence of recombinant sequences
          3. estimate phylogenetic trees from non-recombinant alignments
          4. filter gene trees for outliers with kdetrees test
	  
	  The filtered alignments (non-recombinant and non-outlier) will then enter into the DNA polymorphims analysis,
	    which involves computing basic descriptive statistics from the alignments as well as performing the popular
	    Tajima\'s D and Fu and Li\'s D* neutrality tests. The results are summarized in a table with the following 
	    fifteen columns:
      
              Alignment_name
              no_seqs
              aln_len
              avg_perc_identity
              pars_info_sites
              consistency_idx
              homoplasy_idx
              segregating_sites
              singletons
              pi_per_gene
              pi_per_site
              theta_per_gene
              theta_per_site
              tajimas_D
              fu_and_li_D_star

        This allows the user to identify neutral markers with desirable properties for standard population genetic analyses.

   INVOCATION EXAMPLES:
     1. default on DNA sequences (uses IQ-TREE evaluating a subset of models specified in the detailed help)
          $progname -R 1 -t DNA
     2. thorough FastTree searching and molecular clock analysis on DNA sequences using 10 cores:
          $progname -R 1 -t DNA -A F -k 1.2 -m 0.7 -s 20 -l 12 -T high -K -M HKY -q 0.95 -n 10
     3. FastTree searching on a huge protein dataset for fast inspection
          $progname -R 1 -A F -t PROT -m 0.6 -k 1.0 -T lowest
     4. To run the pipeline on a remote server, we recommend using the nohup command upfront, as shown below:
        nohup $progname -R 1 -t DNA -S 'TNe,TVM,TVMe,GTR' -k 1.0 -m 0.75 -T high -N 5 &> /dev/null &
     5. Run in population-genetics mode (generates a table with descritive statistics for DNA-polymorphisms 
          and the results of diverse neutrality tests)
	  $progname -R 2 -t DNA

   NOTES
     1: run from within the directory holding core gene clusters generated by get_homologues.pl -e or
          compare_clusters.pl with -t no_genome_gbk files (core genome: all clusters with a single gene copy per genome)
     2: If you encounter any problems, please run the script with the -D -V flags added at the end of the command line,
          redirect STOUT to a file and send us the output, so that we can better diagnose the problem.
	  e.g. $progname -R 1 -t DNA -k 1.0 -m 0.7 -s 8 -l 10 -T high -K -D &> ERROR.log
      
USAGE

exit 0

}


#-----------------------------------------------------------------------------------------

function print_help()
{
#     -V flag to activate verbose command execution lines                                           [default: $VERBOSITY]

   cat <<EOF
   $progname v.$VERSION OPTIONS:

   REQUIRED:
    -R <integer> RUNMODE
       1 select optimal markers for phylogenetics/phylogenomics (genomes from different species)
       2 select optimal markers for population genetics (genomes from the same species)
    -t <string> type of input sequences: DNA|PROT

   OPTIONAL:
     -h flag, print this short help notes
     -H flag, print extensive Help and details about search modes and models
     -A <string> Algorithm for tree searching: <F|I> [FastTree|IQ-TREE]                            [default:$search_algorithm]
     -c <integer> NCBI codontable number (1-23) for pal2nal.pl to generate codon alignment         [default:$codontable]
     -C flag to print codontables
     -D flag to print debugging messages; please use if you encounter problems executing the code  [default: $DEBUG]
#    -e <integer> select gene trees with at least (min. = 4) external branches                     [default: $min_no_ext_branches]
     -k <real> kde stringency (0.7-1.6 are reasonable values; less is more stringent)              [default: $kde_stringency]
     -K flag to run molecular clock test on codon alignments                                       [default: $eval_clock]
     -l <integer> max. spr length (7-12 are reasonable values)                                     [default: $spr_length]
     -m <real> min. average support value (0.7-0.8 are reasonable values) for trees to be selected [default: $min_supp_val]
     -M <string> base Model for clock test (use one of: GTR|TrN|HKY|K2P|F81); uses +G in all cases [default: $base_mod]
     -n <integer> number of cores/threads to use                                                   [default: all cores]
     -N <integer> number of IQ-TREE searches to run [only active with -T high]                     [default: $nrep_IQT_searches]
     -q <real> quantile (0.95|0.99) of Chi-square distribution for computing molec. clock p-value  [default: $q]
     -r <string> root method (midpoint|outgroup)                                                   [default: $root_method]
     -s <integer> number of spr rounds (4-20 are reasonable values) for FastTree tree searching    [default: $spr]
     -S <string> quoted 'comma-separated list' of base models to be evaluated by IQ-TREE 
                 when estimating the species tree from the concatenated supermatrix  (see -H for details). 
		 If no -S is passed, then sinlge default models are used, as shown below
              <'JC,F81,K2P,HKY,TrN,TNe,K3P,K81u,TPM2,TPM2u,TPM3,TPM3u,
	      TIM,TIMe,TIM2,TIM2e,TIM3,TIM3e,TVM,TVMe,SYM,GTR'>              for DNA alignments    [default: $IQT_DNA_models]
              <'BLOSUM62,cpREV,Dayhoff,DCMut,FLU,HIVb,HIVw,JTT,JTTDCMut,LG,
                mtART,mtMAM,mtREV,mtZOA,Poisson,PMB,rtREV,VT,WAG'>           for PROT alignments   [default: $IQT_PROT_models]
     -T <string> tree search Thoroughness: high|medium|low|lowest (see -H for details)             [default: $search_thoroughness]

   INVOCATION EXAMPLES:
     1. default on DNA sequences (uses IQ-TREE evaluating a subset of models specified in the detailed help)
          $progname -R 1 -t DNA
     2. thorough FastTree searching and molecular clock analysis on DNA sequences using 10 cores:
          $progname -R 1 -t DNA -A F -k 1.2 -m 0.7 -s 20 -l 12 -T high -K -M HKY -q 0.95 -n 10
     3. FastTree searching on a huge protein dataset for fast inspection
          $progname -R 1 -A F -t PROT -m 0.6 -k 1.0 -T lowest
     4. To run the pipeline on a remote server, we recommend using the nohup command upfront, as shown below:
        nohup $progname -R 1 -t DNA -S 'TNe,TVM,TVMe,GTR' -k 1.0 -m 0.75 -T high -N 5 &> /dev/null &
     5. Run in population-genetics mode (generates a table with descritive statistics for DNA-polymorphisms 
          and the results of diverse neutrality test)
	  $progname -R 2 -t DNA

   NOTES
     1: run from within the directory holding core gene clusters generated by get_homologues.pl -e or
          compare_clusters.pl with -t no_genome_gbk files (core genome: all clusters with a single gene copy per genome)
     2: If you encounter any problems, please run the script with the -D -V flags added at the end of the command line,
          redirect STOUT to a file and send us the output, so that we can better diagnose the problem.
	  e.g. $progname -R 1 -t DNA -k 1.0 -m 0.7 -s 8 -l 10 -T high -K -D &> ERROR.log

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
search_thoroughness='medium'
kde_stringency=1.5
min_supp_val=0.7
min_no_ext_branches=4
n_cores=
#VERBOSITY=0
spr_length=10
spr=4
codontable=11 # bacterial by default, sorry for the bias ;)
base_mod=GTR
eval_clock=0
root_method=midpoint
tree_prefix=concat
q=0.99

search_algorithm=I
IQT_DNA_models=GTR
IQT_PROT_models=LG
IQT_models=
nrep_IQT_searches=5

while getopts ':c:e:k:l:m:M:n:N:p:q:r:s:t:A:T:R:S:hCDHK' OPTIONS
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
   K)   eval_clock=1
        ;;
   l)   spr_length=$OPTARG
        ;;
   m)   min_supp_val=$OPTARG
        ;;
   M)   base_mod=$OPTARG
        ;;
   n)   n_cores=$OPTARG
        ;;
   N)   nrep_IQT_searches=$OPTARG
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
    :)   sprintf "argument missing from -%s option\n" "-$OPTARG" 
   	 print_help
     	 exit 2
     	 ;;
   \?)   echo "invalid option: -$OPTARG"
   	 print_help
         exit 3
	 ;;
    *)   msg "An unexpected parsing error occurred" ERROR RED
         echo
         print_help
	 exit 4
	 ;;
   esac >&2   # print the ERROR MESSAGES to STDERR
done

shift $((OPTIND - 1))


#-------------------------------------------------------#
# >>>BLOCK 0.1 SET THE ENVIRONMENT FOR THE PIPELINE <<< #
#-------------------------------------------------------#

# 0. Set the distribution base directory and OS-specific (linux|darwin) bindirs
env_vars=$(set_pipeline_environment) # returns: $distrodir $bindir $OS $no_proc
[ $DEBUG ] && echo "env_vars:$env_vars"
distrodir=$(echo $env_vars|awk '{print $1}')
bindir=$(echo $env_vars|awk '{print $2}')
OS=$(echo $env_vars|awk '{print $3}')
no_proc=$(echo $env_vars|awk '{print $4}')

# source get_phylomarkers_fun_lib into the script to get access to reminder of the functions
. ${distrodir}/lib/get_phylomarkers_fun_lib

#-----------------------------------------------------------------------------------------

[ $DEBUG -eq 1 ] && msg "distrodir:$distrodir|bindir:$bindir|OS:$OS|no_proc:$no_proc" DEBUG LBLUE

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
       msg "# ERROR: no runmode defined!" HELP RED
       print_help
       exit 1
fi

if [ "$search_algorithm" != "I" -a "$search_algorithm" != "F" ]
then
       msg "# ERROR: search_algorithm $search_algorithm is not recognized!" ERROR RED
       print_help
       exit 1
fi

if [ $min_no_ext_branches -lt 4 ]
then
    msg ">>> ERROR: -e has to be >= 4" HELP RED
    print_help
    exit 1
fi

if [ -z $n_cores ]
then
     n_cores=$no_proc
fi

if [ "$mol_type" != "DNA" -a "$mol_type" != "PROT" ] # "$mol_type" == "BOTH" not implemented yet
then
     msg "ERROR: -t must be DNA or PROT" ERROR RED
     print_help
     exit 1
fi

if [ -z $IQT_models ]
then
   [ "$mol_type" == "DNA" ] &&  IQT_models=$IQT_DNA_models
   [ "$mol_type" == "PROT" ] && IQT_models=$IQT_PROT_models
fi

if [ "$search_algorithm" == "I" -a "$mol_type" == "DNA" ]
then
     check_IQT_DNA_models "$IQT_models"
fi

if [ "$search_algorithm" == "I" -a "$mol_type" == "PROT" ]
then
     check_IQT_PROT_models "$IQT_models"
fi

if [ "$search_thoroughness" != "high" -a "$search_thoroughness" != "medium" -a "$search_thoroughness" != "low" -a "$search_thoroughness" != "lowest" ]
then
     msg "ERROR: -T must be lowest|low|medium|high" HELP RED
     print_help
     exit 1
fi

if [ "$base_mod" != "GTR" -a "$base_mod" != "TrN" -a "$base_mod" != "HKY" -a "$base_mod" != "K2P" -a "$base_mod" != "F81" ]
then
     msg "ERROR: -M must be one of: GTR|TrN|HKY|K2P|F81" HELP RED
     print_help
     exit 1
fi

if [ $eval_clock -gt 0 -a "$mol_type" != "DNA" ] # MolClock currently only with DNA
then
     msg "-K 1 (evaluate clock) must be run on codon alignments with -t DNA" HELP RED
     print_help
     exit 1
fi

if [ $runmode -gt 1 -a "$mol_type" != "DNA" ] # PopGen analysis currently only with DNA
then
     msg "ERROR: runmode $runmode must be run on codon alignments with -t DNA" HELP RED
     print_help
     exit 1
fi

#---------------------#
# >>>> MAIN CODE <<<< #
#---------------------#

[ $DEBUG -eq 1 ] && msg "running on $OSTYPE" && echo "path contains: " DEBUG LBLUE && echo $PATH|sed 's/:/\n/g'

start_time=$(date +%s)

parent_PID=$(get_script_PID $progname)
[ $DEBUG -eq 1 ] && msg "parent_PID:$parent_PID" DEBUG LBLUE

if [ $eval_clock -eq 1 -a $search_algorithm == "F" ]
then
    dir_suffix=A${search_algorithm}R${runmode}t${mol_type}_k${kde_stringency}_m${min_supp_val}_s${spr}_l${spr_length}_T${search_thoroughness}_K
elif [ $eval_clock -ne 1 -a $search_algorithm == "F" ]
then
    dir_suffix=A${search_algorithm}R${runmode}t${mol_type}_k${kde_stringency}_m${min_supp_val}_s${spr}_l${spr_length}_T${search_thoroughness}
elif [ $eval_clock -eq 1 -a $search_algorithm == "I" ]
then
    dir_suffix=A${search_algorithm}R${runmode}t${mol_type}_k${kde_stringency}_m${min_supp_val}_T${search_thoroughness}_K
elif [ $eval_clock -ne 1 -a $search_algorithm == "I" ]
then
    dir_suffix=A${search_algorithm}R${runmode}t${mol_type}_k${kde_stringency}_m${min_supp_val}_T${search_thoroughness}
fi

lmsg=" >>> $(basename $0) vers. $VERSION run with the following parameters:"
msg "$lmsg" PROGR CYAN

lmsg="Run started on $TIMESTAMP_SHORT_HMS under $OSTYPE on $HOSTNAME with $n_cores cores
 wkdir:$wkdir
 distrodir:$distrodir
 bindir:$bindir

 > General run settings:
      runmode:$runmode|mol_type:$mol_type|search_algorithm:$search_algorithm
 > Filtering parameters:
     kde_stringency:$kde_stringency|min_supp_val:$min_supp_val
 > FastTree parameters:
     spr:$spr|spr_length:$spr_length|search_thoroughness:$search_thoroughness
 > IQ-TREE parameters:
     IQT_models:$IQT_models|search_thoroughness:$search_thoroughness|nrep_IQT_searches:$nrep_IQT_searches
 > Molecular Clock parmeters:
     eval_clock:$eval_clock|root_method:$root_method|base_model:$base_mod|ChiSq_quantile:$q

 DEBUG=$DEBUG"
 
msg "$lmsg" PROGR YELLOW

#printf "
# ${CYAN}>>> $(basename $0) vers. $VERSION run with the following parameters:${NC}
# ${YELLOW}Run started on $TIMESTAMP_SHORT_HMS under $OSTYPE on $HOSTNAME with $n_cores cores
# wkdir:$wkdir
# distrodir:$distrodir
# bindir:$bindir
#
# > General run settings:
#      runmode:$runmode|mol_type:$mol_type|search_algorithm:$search_algorithm
# > Filtering parameters:
#     kde_stringency:$kde_stringency|min_supp_val:$min_supp_val
# > FastTree parameters:
#     spr:$spr|spr_length:$spr_length|search_thoroughness:$search_thoroughness
# > IQ-TREE parameters:
#     IQT_models:$IQT_models|search_thoroughness:$search_thoroughness|nrep_IQT_searches:$nrep_IQT_searches
# > Molecular Clock parmeters:
#     eval_clock:$eval_clock|root_method:$root_method|base_model:$base_mod|ChiSq_quantile:$q
#
# DEBUG=$DEBUG${NC}
#
#" | tee -a ${logdir}/get_phylomarkers_run_${dir_suffix}_${TIMESTAMP_SHORT}.log

#----------------------------------------------------------------------------------------------------------------
#>>>BLOCK 1. make a new subdirectory within the one holding core genome clusters generated by compare_clusters.pl
#    and generate symlinks to the faa and fna files (NOTE: BOTH REQUIRED). Fix fasta file names and headers.
#    Check that we have fna faa input FASTA file pairs and that they contain the same number of sequences and
#    instances of each taxon. Mark dir as top_dir
#----------------------------------------------------------------------------------------------------------------

if [ -d get_phylomarkers_run_${dir_suffix}_${TIMESTAMP_SHORT} ]
then
    msg "  >>> ERROR Found and older get_phylomarkers_run_${dir_suffix}_${TIMESTAMP_SHORT}/ directory. Please remove or rename and re-run!" ERROR RED
    exit 2
fi

msg "" PROGR NC
msg " >>>>>>>>>>>>>>> running input data sanity checks <<<<<<<<<<<<<<< " PROGR YELLOW
msg "" PROGR NC


# make sure we have *.faa and *.fna file pairs to work on
nfna=$(find . -name "*.fna" | wc -l)
if [ "$?" -ne "0" ]
then
   msg " >>> ERROR: there are no input fna files to work on!\n\tPlease check input FASTA files: [you may need to run compare_clusters.pl with -t NUM_OF_INPUT_GENOMES -n]\n\tPlease check the GET_HOMOLOGUES manual" ERROR RED
msg "http://eead-csic-compbio.github.io/get_homologues/manual/" ERROR BLUE
   exit 2
fi

nfaa=$(find . -name "*.faa" | wc -l)
if [ "$?" -ne "0" ]
then
   msg " >>> ERROR: there are no input faa files to work on!\n\tPlease check input FASTA files: [you may need to run compare_clusters.pl with -t NUM_OF_INPUT_GENOMES]\n\tPlease check the GET_HOMOLOGUES manual" ERROR RED
   msg "http://eead-csic-compbio.github.io/get_homologues/manual/" ERROR BLUE
   exit 2
fi

if [ "$nfna" -ne "$nfaa" ]
then
  lmsg=" >>> ERROR: there are no equal numbers of fna and faa input files to work on!\n\tPlease check input FASTA files: [you may need to run compare_clusters.pl with -t NUM_OF_INPUT_GENOMES -n; and a second time: run compare_clusters.pl with -t NUM_OF_INPUT_GENOMES]"
  msg "$lmsg" ERROR RED
  msg "  Please check the GET_HOMOLOGUES manual at: " ERROR RED
  msg "  http://eead-csic-compbio.github.io/get_homologues/manual/" ERROR LBLUE
fi

mkdir get_phylomarkers_run_${dir_suffix}_${TIMESTAMP_SHORT} && cd get_phylomarkers_run_${dir_suffix}_${TIMESTAMP_SHORT}
top_dir=$(pwd)

print_start_time && msg "# processing source fastas in directory get_phylomarkers_run_${dir_suffix}_${TIMESTAMP_SHORT} ..." PROGR BLUE

ln -s ../*.faa .
ln -s ../*.fna .

# fix fasta file names with two and three dots
$distrodir/rename.pl 's/\.\.\./\./g' ./*.faa
$distrodir/rename.pl 's/\.\.\./\./g' ./*.fna

# 1.0 check that all fasta files contain the same number of sequences
FASTASIZES=$(grep -c "^>" ./*.f[na]a | cut -d":" -f 2 | sort | uniq)
NSEQSFASTA=$(grep -c "^>" ./*.f[na]a | cut -d":" -f 2 | sort | uniq | wc -l)
[ $NSEQSFASTA -gt 1 ] && msg " >>> ERROR: Input FASTA files do not contain the same number of sequences...${NC}\n$FASTASIZES\n" ERROR RED && exit 4

# 1.1 fix fastaheaders of the source protein and DNA fasta files
for file in ./*faa; do awk 'BEGIN {FS = "|"}{print $1, $2, $3}' $file|perl -pe 'if(/^>/){s/>\S+/>/; s/>\h+/>/; s/\h+/_/g; s/,//g; s/;//g; s/://g; s/\(//g; s/\)//g}' > ${file}ed; done
for file in ./*fna; do awk 'BEGIN {FS = "|"}{print $1, $2, $3}' $file|perl -pe 'if(/^>/){s/>\S+/>/; s/>\h+/>/; s/\h+/_/g; s/,//g; s/;//g; s/://g; s/\(//g; s/\)//g}' > ${file}ed; done

print_start_time && printf "${BLUE}# Performing strain composition check on f?aed files ...${NC}\n"
faaed_strain_intersection_check=$(grep '>' ./*faaed | cut -d: -f2 | sort | uniq -c | awk '{print $1}' | sort | uniq -c | wc -l)
fnaed_strain_intersection_check=$(grep '>' ./*fnaed | cut -d: -f2 | sort | uniq -c | awk '{print $1}' | sort | uniq -c | wc -l)

# check that each file has the same number of strains and a single instance for each strain
if [ $faaed_strain_intersection_check -eq 1 -a $fnaed_strain_intersection_check -eq 1 ]
then
   msg " >>> Strain check OK: each f?aed file has the same number of strains and a single instance for each strain" PROGR GREEN
else
     msg " >>> ERROR: Input f?aed files do not contain the same number of strains and a single instance for each strain...\n\tPlease check input FASTA files: [you may need to run compare_clusters.pl with -t NUM_OF_INPUT_GENOMES]\n\tPlease check the GET_HOMOLOGUES manual" ERROR RED
     msg "http://eead-csic-compbio.github.io/get_homologues/manual" ERROR BLUE
     exit 5
fi


# 1.2 add_nos2fasta_header.pl to avoid problems with duplicate labels
[ $DEBUG -eq 1 ] && msg " > ${distrodir}/run_parallel_cmmds.pl faaed 'add_nos2fasta_header.pl $file > ${file}no' $n_cores &> /dev/null" DEBUG NC
${distrodir}/run_parallel_cmmds.pl faaed 'add_nos2fasta_header.pl $file > ${file}no' $n_cores &> /dev/null

[ $DEBUG -eq 1 ] && msg " > ${distrodir}/run_parallel_cmmds.pl fnaed 'add_nos2fasta_header.pl $file > ${file}no' $n_cores &> /dev/null" DEBUG NC
${distrodir}/run_parallel_cmmds.pl fnaed 'add_nos2fasta_header.pl $file > ${file}no' $n_cores &> /dev/null

no_alns=$(find . -name "*.fnaedno" | wc -l)

[ $no_alns -eq 0 ] && msg " >>> ERROR: There are no codon alignments to work on! Something went wrong. Please check input and settings ... " ERROR RED && exit 4
print_start_time && msg "# Total number of alignments to be computed $no_alns" PROGR BLUE

# 1.3 generate a tree_labels.list file for later tree labeling
print_start_time && msg "# generating the labels file for tree-labeling ..." PROGR BLUE

tree_labels_dir=$(pwd)
grep '>' "$(find . -name "*fnaedno" | head -1)" > tree_labels.list

[ $DEBUG -eq 1 ] && msg " > perl -pe '$c++; s/>/$c\t/; s/\h\[/_[/' tree_labels.list > ed && mv ed tree_labels.list" DEBUG NC
perl -pe '$c++; s/>/$c\t/; s/\h\[/_[/' tree_labels.list > ed && mv ed tree_labels.list

#------------------------------------------------------------------------------------------------------#
# >>>BLOCK 2. Generate cdnAlns with with pal2nal, maintaining the input-order of the source fastas <<< #
#------------------------------------------------------------------------------------------------------#

# 2.1 generate the protein alignments using clustalo
msg "" PROGR NC
msg " >>>>>>>>>>>>>>> parallel clustalo and pal2nal runs to generate protein and codon alignments <<<<<<<<<<<<<<< " PROGR YELLOW
msg "" PROGR NC


print_start_time &&  msg "# generating $no_alns protein alignments ..." PROGR BLUE
[ $DEBUG -eq 1 ] && msg " > '${distrodir}/run_parallel_cmmds.pl faaedno clustalo -i $file -o ${file%.*}_cluo.faaln --output-order input-order' $n_cores &> /dev/null" DEBUG NC
${distrodir}/run_parallel_cmmds.pl faaedno 'clustalo -i $file -o ${file%.*}_cluo.faaln --output-order input-order' $n_cores &> clustalo.log

if grep -q "Thread creation failed" clustalo.log; then
   msg " >>> ERROR: This system cannot launch too many threads, please use option -n and re-run ..." ERROR RED && exit 4
fi

# 2.2 generate the codon alignments (files with *_cdnAln.fasta extension) using pal2nal.pl,
#     excluding gaps, and mismatched codons, assuming a bacterial genetic code
# NOTE: to execute run_parallel_cmmds.pl with a customized command, resulting from the interpolation of multiple varialbles,
#       we have to first construct the command line in a variable and pipe its content into bash for executio

print_start_time && msg "# running pal2nal to generate codon alignments ..." PROGR LBLUE

faaln_ext=faaln
command="${distrodir}/run_parallel_cmmds.pl $faaln_ext '${distrodir}/pal2nal.pl \$file \${file%_cluo.faaln}.fnaedno -output fasta -nogap -nomismatch -codontable $codontable > \${file%_cluo.faaln}_cdnAln.fasta' $n_cores"

# now we can execute run_parallel_cmmds.pl with a customized command, resulting from the interpolation of multiple varialbles
[ $DEBUG -eq 1 ] && msg " > $command | bash &> /dev/null" DEBUG NC
echo "$command" | bash &> /dev/null

# check we got non-empty *cdnAln.fasta files
for f in ./*cdnAln.fasta
do
     if [ ! -s $f ]
     then
           msg " >>> Warning: produced empty codon alignment $f!" WARNING LRED
	   msg "    ... Will skip this locus and move it to problematic_alignments ..." WARNING LRED
	   [ ! -d problematic_alignments ] && mkdir problematic_alignments
	   locus_base=${f%_cdnAln.fasta}
	   mv $f problematic_alignments
	   mv ${locus_base}* problematic_alignments
     fi
done

# 2.3 cleanup: remove the source faa, fna, fnaed and faaed files; make numbered_fna_files.tgz and numbered_faa_files.tgz; rm *aedno
rm ./*fnaed ./*faaed ./*faa ./*fna
[ "$DEBUG" -eq "1" ] && msg " > tar -czf numbered_fna_files.tgz ./*fnaedno" DEBUG NC
tar -czf numbered_fna_files.tgz ./*fnaedno
[ "$DEBUG" -eq "1" ] && msg " > tar -czf numbered_fna_files.tgz ./*faaedno" DEBUG NC
tar -czf numbered_faa_files.tgz ./*faaedno
rm ./*aedno

#---------------------------------------------------------------------------------------------------------#
#>>>BLOCK 3. run Phi-test to identify recombinant codon alignments on all *_cdnAln.fasta source files <<< #
#---------------------------------------------------------------------------------------------------------#
# 3.1 make a new PhiPack subdirectory to work in. generate symlinks to ../*fasta files
#     Mark dir as phipack_dir
mkdir PhiPack && cd PhiPack
phipack_dir=$(pwd)
ln -s ../*fasta .

# 3.1.2 check that we have codon alignments before proceeding
no_fasta_files=$(find . -name "*.fasta" | wc -l)

[ $no_fasta_files -lt 1 ] && print_start_time && msg " >>> ERROR: there are no codon alignments to run Phi on. Will exit now!" ERROR RED && exit 3

msg "" PROGR NC
msg " >>>>>>>>>>>>>>> parallel phi(w) runs to identify alignments with recombinant sequences <<<<<<<<<<<<<<< " PROGR YELLOW
msg "" PROGR NC

# 3.2 run Phi from the PhiPack in parallel
print_start_time && msg "# running Phi recombination test in PhiPack dir ..." PROGR LBLUE
[ $DEBUG -eq 1 ] && msg " > ${distrodir}/run_parallel_cmmds.pl fasta 'Phi -f $file -p 1000 > ${file%.*}_Phi.log' $n_cores &> /dev/null" DEBUG NC
${distrodir}/run_parallel_cmmds.pl fasta 'Phi -f $file -p 1000 > ${file%.*}_Phi.log' $n_cores &> /dev/null

# 3.3 process the *_Phi.log files generated by Phi to write a summary table and print short overview to STDOUT
declare -a nonInfoAln
COUNTERNOINFO=0

[ -s Phi.log ] && rm Phi.log

for f in ./*_Phi.log
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
total_no_cdn_alns=$(ls ./*_cdnAln.fasta | wc -l)

# 3.4 Check the number of non-recombinant alignments remaining
if [ ${#nonInfoAln[@]} == 0 ]
then
  lmsg=" >>> Phi test result: there are $no_non_recomb_alns_perm_test non-recombinant alignments out of $total_no_cdn_alns input alignments"
  print_start_time && msg "$lmsg" PROGR GREEN
fi

if [ ${#nonInfoAln[@]} -gt 0 ]
then
  lmsg=" >>> Phi test WARNING: there ${#nonInfoAln[@]} alignments with too few informative sites for the Phi test to work on ... "
  print_start_time && msg "$lmsg" WARNING LRED

  # print the names of the alignments with too few informative sites for the Phi test to be work on
  if [ "$DEBUG" -eq "1" ]
  then
      msg " >>> The alignments with too few informative sites for the Phi test to be work on are:" WARNING LRED
      for f in "${nonInfoAln[@]}"
      do
           msg " >>> ${f}" WARNING LRED
      done
  fi
fi

[ $no_non_recomb_alns_perm_test -lt 1 ] && print_start_time && msg " >>> ERROR: All alignments seem to have recombinant sequences. will exit now!" ERROR RED && exit 3

#3.5 cleanup dir
tar -czf Phi_test_log_files.tgz ./*Phi.log
[ -s Phi_test_log_files.tgz ] && rm ./*Phi.log Phi.inf*

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
rm ./*cdnAln.fasta


#=====================================#
# >>>  Block 4. -t DNA RUNMODES   <<< #
#=====================================#
#
#  NOTES:
#   * cd into non_recomb_cdn_alns and run FastTree/IQ-TREE in parallel on all codon alignments.
#   * For FastTree -T controls desired search thoroughness      (see details with -H)
#   * For IQT -T controls the number of models to be evaluated  (see details with -H)
# 
# 	 1. The if statement below divides the flow into two large if blocks to run -t DNA|PROT
# 		    [ "$mol_type" == "DNA" ]
# 		    [ "$mol_type" == "PROT" ]
#   BLOCK 4 will work on DNA, either in phylogenetics (-R 1) or population genetics (-R 2 runmodes)
#   BLOCK 5 works on protein sequences in -R 1 exclusively


#:::::::::::::::::::::::::::#
# >>> Filter gene trees <<< #
#:::::::::::::::::::::::::::#
if [ "$mol_type" == "DNA" ]
then
    cd non_recomb_cdn_alns
    non_recomb_cdn_alns_dir=$(pwd)

    print_start_time && msg "# working in dir non_recomb_cdn_alns ..." PROGR LBLUE
    print_start_time && msg "# estimating $no_non_recomb_alns_perm_test gene trees from non-recombinant sequences ..." PROGR LBLUE

    # 4.1 >>> estimate_IQT_gene_trees | estimate_FT_gene_trees
    if [ "$search_algorithm" == "F" ]
    then
	msg "" PROGR NC
	msg " >>>>>>>>>>>>>>> parallel FastTree runs to estimate gene trees <<<<<<<<<<<<<<< " PROGR YELLOW
	msg "" PROGR NC

        gene_tree_ext="ph"
	lmsg=" > running estimate_FT_gene_trees $mol_type $search_thoroughness ..."
        [ $DEBUG -eq 1 ] && msg "$lmsg" DEBUG NC
	estimate_FT_gene_trees $mol_type $search_thoroughness

	# 4.1.1 check that FT computed the expected gene trees
        no_gene_trees=$(find . -name "*.ph" | wc -l)
        [ $no_gene_trees -lt 1 ] && print_start_time && msg " >>> ERROR: There are no gene tree to work on in non_recomb_cdn_alns/. will exit now!" ERROR RED && exit 3

	# 4.1.2 generate computation-time and lnL stats
	[ $DEBUG -eq 1 ] && msg "compute_FT_gene_tree_stats $mol_type $search_thoroughness" DEBUG NC
	compute_FT_gene_tree_stats $mol_type $search_thoroughness
    else
        gene_tree_ext="treefile"
	msg "" PROGR NC
	msg " >>>>>>>>>>>>>>> parallel IQ-TREE runs to estimate gene trees <<<<<<<<<<<<<<< " PROGR YELLOW
	msg "" PROGR NC

	estimate_IQT_gene_trees $mol_type $search_thoroughness $IQT_models

	# 4.1.1 check that IQT computed the expected gene trees
	no_gene_trees=$(find . -name "*.treefile" | wc -l)
        [ $no_gene_trees -lt 1 ] && print_start_time && msg " >>> ERROR: There are no gene tree to work on in non_recomb_cdn_alns/. will exit now!" ERROR RED && exit 3

	# 4.1.2 generate computation-time, lnL and best-model stats
	compute_IQT_gene_tree_stats $mol_type $search_thoroughness
    fi

    #remove trees with < 5 branches
    print_start_time && msg "# counting branches on $no_non_recomb_alns_perm_test gene trees ..." PROGR BLUE
    [ $DEBUG -eq 1 ] && msg " > count_tree_branches $gene_tree_ext no_tree_branches.list &> /dev/null" DEBUG NC
    [ $DEBUG -eq 1 ] && msg " search_thoroughness: ${search_thoroughness}" DEBUG NC

    count_tree_branches $gene_tree_ext no_tree_branches.list # &> /dev/null
    [ ! -s no_tree_branches.list ] && install_Rlibs_msg no_tree_branches.list ape

    check_output no_tree_branches.list $parent_PID | tee -a ${logdir}/get_phylomarkers_run_${dir_suffix}_${TIMESTAMP_SHORT}.log

    # remove trees with < 5 external branches (leaves)
     [ $DEBUG -eq 1 ] && msg " >  removing trees with < 5 external branches (leaves)" DEBUG NC
     
    no_tree_counter=0
    for phy in $(grep -v '^#Tree' no_tree_branches.list | awk -v min_no_ext_branches=$min_no_ext_branches 'BEGIN{FS="\t"; OFS="\t"}$7 < min_no_ext_branches' |cut -f1)
    do
         [ "$search_algorithm" == "F" ] && base=$(echo ${phy//_FTGTR\.ph/})
	 [ "$search_algorithm" == "I" ] && base=$(echo ${phy//\.treefile/})
	 print_start_time && msg " >>> will remove ${base}* because it has < 5 branches" WARNING LRED
	 rm ${base}*
	 let no_tree_counter++
    done

    msg " >>> WARNING: there are $no_tree_counter trees with < 1 internal branches (no real trees) that will be discarded ..." WARNING LRED

    msg "" PROGR NC
    msg " >>>>>>>>>>>>>>> filter gene trees for outliers with kdetrees test <<<<<<<<<<<<<<< " PROGR YELLOW
    msg "" PROGR NC


    if [ "$search_algorithm" == "F" ]
    then
        # 4.1 generate the all_GTRG_trees.tre holding all source trees, which is required by kdetrees
        #     Make a check for the existence of the file to interrupt the pipeline if something has gone wrong
        [ $DEBUG -eq 1 ] && msg " > cat ./*.ph > all_gene_trees.tre" DEBUG NC
	[ $DEBUG -eq 1 ] && msg " search_thoroughness: ${search_thoroughness}" DEBUG NC
        cat ./*.ph > all_gene_trees.tre
        check_output all_gene_trees.tre $parent_PID | tee -a ${logdir}/get_phylomarkers_run_${dir_suffix}_${TIMESTAMP_SHORT}.log
        [ ! -s all_gene_trees.tre ] && exit 3
    else
        # 4.1 generate the all_IQT_trees.tre holding all source trees, which is required by kdetrees
        #     Make a check for the existence of the file to interrupt the pipeline if something has gone wrong
        [ $DEBUG -eq 1 ] && msg " > cat ./*.treefile > all_gene_trees.tre" DEBUG NC
	[ $DEBUG -eq 1 ] && msg " search_thoroughness:$search_thoroughness" DEBUG NC
	cat ./*.treefile > all_gene_trees.tre
        check_output all_gene_trees.tre $parent_PID | tee -a ${logdir}/get_phylomarkers_run_${dir_suffix}_${TIMESTAMP_SHORT}.log
        [ ! -s all_gene_trees.tre ] && exit 3
    fi

    # 4.2 run_kdetrees.R at desired stringency
    print_start_time && msg "# running kde test ..." PROGR BLUE
    [ $DEBUG -eq 1 ] && msg " > ${distrodir}/run_kdetrees.R ${gene_tree_ext} all_gene_trees.tre $kde_stringency &> /dev/null" DEBUG NC
    ${distrodir}/run_kdetrees.R ${gene_tree_ext} all_gene_trees.tre $kde_stringency &> /dev/null
    [ ! -s kde_dfr_file_all_gene_trees.tre.tab ] && install_Rlibs_msg kde_dfr_file_all_gene_trees.tre.tab kdetrees,ape
    check_output kde_dfr_file_all_gene_trees.tre.tab $parent_PID | tee -a ${logdir}/get_phylomarkers_run_${dir_suffix}_${TIMESTAMP_SHORT}.log

    # 4.3 mv outliers to kde_outliers
    no_kde_outliers=$(grep -c outlier kde_dfr_file_all_gene_trees.tre.tab)
    no_kde_ok=$(grep -v outlier kde_dfr_file_all_gene_trees.tre.tab|grep -vc '^file')

    # 4.4 Check how many cdnAlns passed the test and separate into two subirectories those passing and failing the test
    if [ $no_kde_outliers -gt 0 ]
    then
        print_start_time && msg "# making dir kde_outliers/ and moving $no_kde_outliers outlier files into it ..." PROGR BLUE
        mkdir kde_outliers
        [ $search_algorithm == "F" ] && for f in $(grep outlier kde_dfr_file_all_gene_trees.tre.tab|cut -f1|sed 's/\.ph//'); do mv ${f}* kde_outliers; done
	[ $search_algorithm == "I" ] && for f in $(grep outlier kde_dfr_file_all_gene_trees.tre.tab|cut -f1|sed 's/\.treefile//'); do mv ${f}* kde_outliers; done
    else
        print_start_time && msg " >>> there are no kde-test outliers ..." PROGR GREEN
    fi

    if [ $no_kde_ok -gt 0 ]
    then
        print_start_time && msg "# making dir kde_ok/ and linking $no_kde_ok selected files into it ..." PROGR BLUE
        mkdir kde_ok
        cd kde_ok
        ln -s ../*.${gene_tree_ext} .

        print_start_time && msg "# labeling $no_kde_ok gene trees in dir kde_ok/ ..." PROGR BLUE
        [ $DEBUG -eq 1 ] && msg " > ${distrodir}/run_parallel_cmmds.pl ${gene_tree_ext} 'add_labels2tree.pl ../../../tree_labels.list $file' $n_cores &> /dev/null" DEBUG NC
        ${distrodir}/run_parallel_cmmds.pl ${gene_tree_ext} 'add_labels2tree.pl ../../../tree_labels.list $file' $n_cores &> /dev/null

	# remove symbolic links to cleanup kde_ok/
        for f in $(ls ./*.${gene_tree_ext}|grep -v "_ed\.${gene_tree_ext}"); do rm $f; done

        cd ..
    else
        print_start_time && msg "# ERROR There are $no_kde_ok gene trees producing non-significant kde-test results! Increase the actual -k $kde_stringency value. Will stop here!" ERROR RED
	exit 5
    fi

#----------------------------------#
# >>> BLOCK 4.2: PHYLOGENETICS <<< #
#----------------------------------#
    if [ $runmode -eq 1 ]
    then
        # >>> 5.1 compute average bipartition support values for each gene tree
	#         and write them to sorted_aggregated_support_values4loci_ge${min_supp_val_perc}perc.tab
        wkdir=$(pwd)

        msg "" PROGR NC
        msg " >>>>>>>>>>>>>>> filter gene trees by phylogenetic signal content <<<<<<<<<<<<<<< " PROGR YELLOW
        msg "" PROGR NC


	print_start_time && msg "# computing tree support values ..." PROGR BLUE
        [ $DEBUG -eq 1 ] && msg " > compute_suppValStas_and_RF-dist.R $wkdir 1 fasta ${gene_tree_ext} 1 &> /dev/null" DEBUG NC
        compute_suppValStas_and_RF-dist.R $wkdir 1 fasta ${gene_tree_ext} 1 &> /dev/null

        print_start_time && msg "# writing summary tables ..." PROGR BLUE
        min_supp_val_perc=${min_supp_val#0.}
        no_digits=${#min_supp_val_perc}
        [ $no_digits -eq 1 ] && min_supp_val_perc=${min_supp_val_perc}0


	# NOTE: IQ-TREE -alrt 1000 provides support values in 1-100 scale, not as 0-1 as FT!
	#       therefore we convert IQT-bases SH-alrt values to 0-1 scale for consistency
	[ $search_algorithm == "I" ] && min_supp_val=$(echo "$min_supp_val * 100" | bc)

        awk -v min_supp_val=$min_supp_val '$2 >= min_supp_val' sorted_aggregated_support_values4loci.tab > sorted_aggregated_support_values4loci_ge${min_supp_val_perc}perc.tab
        check_output sorted_aggregated_support_values4loci_ge${min_supp_val_perc}perc.tab $parent_PID | \
	tee -a ${logdir}/get_phylomarkers_run_${dir_suffix}_${TIMESTAMP_SHORT}.log

        no_top_markers=$(perl -lne 'END{print $.}' sorted_aggregated_support_values4loci_ge${min_supp_val_perc}perc.tab)
        top_markers_dir="top_${no_top_markers}_markers_ge${min_supp_val_perc}perc"
        top_markers_tab=$(ls sorted_aggregated_support_values4loci_ge${min_supp_val_perc}perc.tab)

        # >>> 5.2 move top-ranking markers to $top_markers_dir
	print_start_time && msg "# making dir $top_markers_dir and moving $no_top_markers top markers into it ..." PROGR LBLUE
        mkdir $top_markers_dir && cd $top_markers_dir
        top_markers_dir=$(pwd)
        ln -s ../$top_markers_tab .
        for base in $(awk '{print $1}' $top_markers_tab|grep -v loci|sed 's/"//g'); do ln -s ../${base}* .; done

        [ $no_top_markers -lt 2 ] && print_start_time && msg " >>> Warning: There are less than 2 top markers. Relax your filtering thresholds. will exit now!" ERROR LRED && exit 3

        msg "" PROGR NC
        msg " >>>>>>>>>>>>>>> generate supermatrix from concatenated, top-ranking alignments <<<<<<<<<<<<<<< " PROGR YELLOW
        msg "" PROGR NC

        # >>> 5.3 generate supermatrix (concatenated alignment)
        print_start_time && msg "# concatenating $no_top_markers top markers into supermatrix ..." PROGR BLUE
        [ $DEBUG -eq 1 ] && msg " > concat_alns fasta $parent_PID &> /dev/null" DEBUG NC
        concat_alns fasta $parent_PID &> /dev/null

        # >>> 5.4 remove uninformative sites from the concatenated alignment to speed up computation
        print_start_time && msg "# removing uninformative sites from concatenated alignment ..." PROGR BLUE
        [ $DEBUG -eq 1 ] && msg " > ${distrodir}/remove_uninformative_sites_from_aln.pl < concat_cdnAlns.fna > concat_cdnAlns.fnainf" DEBUG NC
        #$distrodir/remove_uninformative_sites_from_aln.pl < concat_cdnAlns.fna > concat_cdnAlns.fnainf
	${distrodir}/remove_uninformative_sites_from_aln.pl < concat_cdnAlns.fna > concat_cdnAlns.fnainf
        check_output concat_cdnAlns.fnainf $parent_PID | tee -a ${logdir}/get_phylomarkers_run_${dir_suffix}_${TIMESTAMP_SHORT}.log

	[ ! -s concat_cdnAlns.fnainf ] && print_start_time && msg " >>> ERROR: The expected file concat_cdnAlns.fnainf was not produced! will exit now!" ERROR RED && exit 3


        #:::::::::::::::::::::::::::::::#
        # >>> FastTree species tree <<< #
        #:::::::::::::::::::::::::::::::#

	if [ "$search_algorithm" == "F" ]	
        then
	   msg "" PROGR NC
	   msg " >>>>>>>>>>>>>>> FastTree run on supermatrix to estimate the species tree <<<<<<<<<<<<<<< " PROGR YELLOW
	   msg "" PROGR NC

	  # 5.4 run FasTree under the GTR+G model
          print_start_time && msg "# running FastTree on the supermatrix with $search_thoroughness thoroughness. This may take a while ..." PROGR BLUE
	  
          [ $DEBUG -eq 1 ] && msg " search_thoroughness: $search_thoroughness" DEBUG NC
	  if [ "$search_thoroughness" == "high" ]
          then
            $bindir/FastTree -quiet -nt -gtr -gamma -bionj -slow -slownni -mlacc 3 -spr $spr -sprlength $spr_length \
	    -log ${tree_prefix}_nonRecomb_KdeFilt_cdnAlns_FTGTRG.log < concat_cdnAlns.fnainf > \
	     ${tree_prefix}_nonRecomb_KdeFilt_cdnAlns_FTGTRG.ph
          fi

          if [ "$search_thoroughness" == "medium" ]
          then
            $bindir/FastTree -quiet -nt -gtr -gamma -bionj -slownni -mlacc 2 -spr $spr -sprlength $spr_length \
	    -log ${tree_prefix}_nonRecomb_KdeFilt_cdnAlns_FTGTRG.log < concat_cdnAlns.fnainf > \
	    ${tree_prefix}_nonRecomb_KdeFilt_cdnAlns_FTGTRG.ph
          fi

          if [ "$search_thoroughness" == "low" ]
          then
            $bindir/FastTree -quiet -nt -gtr -gamma -bionj -spr $spr -sprlength $spr_length \
	    -log ${tree_prefix}_nonRecomb_KdeFilt_cdnAlns_FTGTRG.log < concat_cdnAlns.fnainf > \
	    ${tree_prefix}_nonRecomb_KdeFilt_cdnAlns_FTGTRG.ph
          fi

          if [ "$search_thoroughness" == "lowest" ]
          then
            $bindir/FastTree -quiet -nt -gtr -gamma -mlnni 4 -log ${tree_prefix}_nonRecomb_KdeFilt_cdnAlns_FTGTRG.log \
	    < concat_cdnAlns.fnainf > ${tree_prefix}_nonRecomb_KdeFilt_cdnAlns_FTGTRG.ph
          fi

          check_output ${tree_prefix}_nonRecomb_KdeFilt_cdnAlns_FTGTRG.ph $parent_PID | tee -a ${logdir}/get_phylomarkers_run_${dir_suffix}_${TIMESTAMP_SHORT}.log

	  if [ -s "${tree_prefix}_nonRecomb_KdeFilt_cdnAlns_FTGTRG.log" -a -s "${tree_prefix}_nonRecomb_KdeFilt_cdnAlns_FTGTRG.ph" ]
	  then
	    #lnL=$(grep ML_Lengths2 "${tree_prefix}_nonRecomb_KdeFilt_cdnAlns_FTGTRG.log" | grep TreeLogLk | sed 's/TreeLogLk[[:space:]]ML_Lengths2[[:space:]]//')
	    lnL=$(grep '^Gamma20LogLk' "${tree_prefix}_nonRecomb_KdeFilt_cdnAlns_FTGTRG.log" |awk '{print $2}')
	    msg " >>> lnL for ${tree_prefix}_nonRecomb_KdeFilt_cdnAlns_FTGTRG.ph = $lnL" PROGR GREEN
	  else
	    msg " >>> ERROR: ${tree_prefix}_nonRecomb_KdeFilt_cdnAlns_FTGTRG.log could not be produced, will stop here" ERROR LRED
	    exit 5
	  fi

          print_start_time && msg "# Adding labels back to tree ..." PROGR BLUE

	  longmsg=" > ${distrodir}/add_labels2tree.pl ${tree_labels_dir}/tree_labels.list ${tree_prefix}_nonRecomb_KdeFilt_cdnAlns_FTGTRG.ph &> /dev/null"
          [ $DEBUG -eq 1 ] && msg "$longmsg" DEBUG NC
          ${distrodir}/add_labels2tree.pl ${tree_labels_dir}/tree_labels.list ${tree_prefix}_nonRecomb_KdeFilt_cdnAlns_FTGTRG.ph &> /dev/null

          if [ -s ${tree_prefix}_nonRecomb_KdeFilt_cdnAlns_FTGTRG_ed.ph ]
          then
	     # for compute_suppValStats_and_RF-dist.R
             mv ${tree_prefix}_nonRecomb_KdeFilt_cdnAlns_FTGTRG_ed.ph ${tree_prefix}_nonRecomb_KdeFilt_${no_top_markers}cdnAlns_FTGTRG_ed.sptree
             check_output ${tree_prefix}_nonRecomb_KdeFilt_${no_top_markers}cdnAlns_FTGTRG_ed.sptree $parent_PID | \
	     tee -a ${logdir}/get_phylomarkers_run_${dir_suffix}_${TIMESTAMP_SHORT}.log
             msg " >>> found in dir $top_markers_dir ..." PROGR GREEN
          else
             msg " >>> WARNING: ${tree_prefix}_nonRecomb_KdeFilt_${no_top_markers}cdnAlns_FTGTRG_ed.sptree could not be produced" WARNING LRED
          fi

          print_start_time && msg "# computing the mean support values and RF-distances of each gene tree to the concatenated tree ..." PROGR BLUE
	  [ $DEBUG -eq 1 ] && msg " > compute_suppValStas_and_RF-dist.R $top_markers_dir 2 fasta ph 1 &> /dev/null" DEBUG NC
          $distrodir/compute_suppValStas_and_RF-dist.R $top_markers_dir 2 fasta ph 1 &> /dev/null

	  # top100_median_support_values4loci.tab should probably not be written in the first instance
	  [ -s top100_median_support_values4loci.tab -a "${no_top_markers}" -lt 101 ] && rm top100_median_support_values4loci.tab

       fi # [ "$search_algorithm" == "F" ]


       #::::::::::::::::::::::::::::::#
       # >>> IQ-TREE species tree <<< #
       #::::::::::::::::::::::::::::::#
       if [ "$search_algorithm" == "I" ]
       then
	  # 5.5 run IQ-tree in addition to FastTree, if requested
          msg "" PROGR NC
	  msg " >>>>>>>>>>>>>>> ModelFinder + IQ-TREE run on supermatrix to estimate the species tree <<<<<<<<<<<<<<< " PROGR YELLOW
	  msg "" PROGR NC

	  print_start_time && msg "# running ModelFinder on the concatenated alignment with $IQT_models. This will take a while ..." PROGR BLUE

          iqtree -s concat_cdnAlns.fnainf -st DNA -mset "$IQT_models" -m MF -nt AUTO -fast &> /dev/null

	  check_output concat_cdnAlns.fnainf.log $parent_PID | tee -a ${logdir}/get_phylomarkers_run_${dir_suffix}_${TIMESTAMP_SHORT}.log

	  best_model=$(grep '^Best-fit model' concat_cdnAlns.fnainf.log | cut -d' ' -f 3)
	  msg " >>> Best-fit model: ${best_model} ..." PROGR GREEN

	  mkdir iqtree_abayes && cd iqtree_abayes
	  ln -s ../concat_cdnAlns.fnainf .

	  if [ "$search_thoroughness" == "high" ]
	  then
	     lmsg="# Will launch $nrep_IQT_searches IQ-TREE searches on the supermatrix with best model ${best_model} -abayes -bb 1000!.
	                This will take a while ..."
	     print_start_time && msg "$lmsg" PROGR BLUE

	     # run nrep_IQT_searches IQ-TREE searches under the best-fit model found
	     for ((rep=1;rep<=nrep_IQT_searches;rep++))
	     do
	         print_start_time && msg " > iqtree -s concat_cdnAlns.fnainf -st DNA -m $best_model -abayes -bb 1000 -nt AUTO -pre abayes_run${rep} &> /dev/null" PROGR LBLUE

		 iqtree -s concat_cdnAlns.fnainf -st DNA -m $best_model -abayes -bb 1000 -nt AUTO -pre abayes_run${rep} &> /dev/null
	     done

	     grep '^BEST SCORE' ./*log | sed 's#./##' | sort -nrk5 > sorted_IQ-TREE_searches.out

	     check_output sorted_IQ-TREE_searches.out $parent_PID | tee -a ${logdir}/get_phylomarkers_run_${dir_suffix}_${TIMESTAMP_SHORT}.log
	     best_search=$(head -1 sorted_IQ-TREE_searches.out)
	     best_search_base_name=$(head -1 sorted_IQ-TREE_searches.out | cut -d\. -f 1)

	     msg "# >>> Best IQ-TREE run was: $best_search ..." PROGR GREEN
	     best_tree_file=${tree_prefix}_${best_search_base_name}_nonRecomb_KdeFilt_iqtree_${best_model}.treefile
	     
	     # Note: this function works within iqtree_abayes/ and takes care of:
	     # 1. labeling the species tree and moving it to top_markers_dir 
	     # 2. making a cleanup in iqtree_abayes/
	     # 3. moves to top_markers_dir to compute RF-dist of gene-trees to species-tree
	     # 4. removes the double extension name *.fasta.treefile and changes treefile for ph to make it paup-compatible for clock-test
             process_IQT_species_trees_for_molClock $best_search_base_name $best_tree_file
	  else
	     print_start_time && msg "# running IQ-tree on the concatenated alignment with best model ${best_model} -abayes -bb 1000. This will take a while ..." PROGR BLUE

	     iqtree -s concat_cdnAlns.fnainf -st DNA -m "$best_model" -abayes -bb 1000 -nt AUTO -pre iqtree_abayes &> /dev/null

	     grep '^BEST SCORE' ./*log | sed 's#./##' | sort -nrk5 > sorted_IQ-TREE_searches.out

	     check_output sorted_IQ-TREE_searches.out $parent_PID | tee -a ${logdir}/get_phylomarkers_run_${dir_suffix}_${TIMESTAMP_SHORT}.log
	     best_search=$(head -1 sorted_IQ-TREE_searches.out)
	     best_search_base_name=$(head -1 sorted_IQ-TREE_searches.out | cut -d\. -f 1)

	     msg "# >>> Best IQ-TREE run was: $best_search ..." PROGR GREEN
	     best_tree_file=${tree_prefix}_nonRecomb_KdeFilt_iqtree_${best_model}.treefile
             process_IQT_species_trees_for_molClock iqtree_abayes $best_tree_file
	  fi
       fi # if [ "$search_algorithm" == "I" ]


        # NOTE: after v0.9 this process is prallelized with run_parallel_cmmds.pl
	if [ $eval_clock -gt 0 ]
        then
	     msg "" PROGR NC
	     msg " >>>>>>>>>>>>>>> TESTING THE MOLECULAR CLOCK HYPOTHESIS <<<<<<<<<<<<<<< " PROGR YELLOW
	     msg "" PROGR NC

 	     # 1. convert fasta2nexus
             print_start_time && msg "# converting fasta files to nexus files" PROGR BLUE
	     [ $DEBUG -eq 1 ] && msg " > $distrodir/convert_aln_format_batch_bp.pl fasta fasta nexus nex &> /dev/null" DEBUG NC
             $distrodir/convert_aln_format_batch_bp.pl fasta fasta nexus nex &> /dev/null

	     # FIX the nexus file format produced by bioperl: (recent paup version error message provided below)
	     # User-defined symbol 'A' conflicts with predefined DNA state symbol.
             # If you are using a predefined format ('DNA', 'RNA', 'nucleotide', or 'protein'),
             # you may not specify predefined states for this format as symbols in  the Format command.
	     for nexusf in ./*.nex
	     do
	       perl -pe 'if(/^format /){ s/symbols.*$/;/}' $nexusf > ed && mv ed $nexusf
	     done

	     print_start_time && msg "# Will test the molecular clock hypothesis for $no_top_markers top markers. This will take some time ..." PROGR BLUE
             #run_molecClock_test_jmodeltest2_paup.sh -R 1 -M $base_mod -t ph -e fasta -b molec_clock -q $q &> /dev/null

	     # 2. >>> print table header and append results to it
             no_dot_q=$(echo ${q//\./})
             results_table=mol_clock_M${base_mod}G_r${rooting_method}_o${outgroup_OTU_nomol}_q${no_dot_q}_ClockTest.tab
             echo -e "#nexfile\tlnL_unconstr\tlnL_clock\tLRT\tX2_crit_val\tdf\tp-val\tmol_clock" > $results_table

	     cmd="${distrodir}/run_parallel_cmmds.pl nex '${distrodir}/run_parallel_molecClock_test_with_paup.sh -R 1 -f \$file -M $base_mod -t ph -b global_mol_clock -q $q' $n_cores"
	     [ $DEBUG -eq 1 ] && msg "run_parallel_molecClock.cmd: $cmd" DEBUG NC
	     echo $cmd | bash &> /dev/null

	     mol_clock_tab=$(ls ./*_ClockTest.tab)

       	     if [ -s $mol_clock_tab ]
	     then
	        msg " >>> generated the molecular clock results file $mol_clock_tab ..." PROGR GREEN

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

	       msg " >>> Top markers and associated stats are found in: $top_markers_dir ..." PROGR GREEN
            else
	       msg " >>> ${mol_clock_tab} not found" ERROR RED
	    fi
        fi # if [ $eval_clock -gt 0 ]

 #---------------------#
 # >>> 5.4 cleanup <<< #
 #---------------------#

	if [ $search_algorithm == "F" ]
	then
	    if [ $eval_clock -eq 1 ]
	    then
	        tar -czf molClock_PAUP_files.tgz ./*_paup.block ./*.nex  ./*_clockTest.log ./*tre ./*clock.scores ./*critical_X2_val.R \
	        $mol_clock_tab ${mol_clock_tab}sorted mol_clock_MGTRG_r_o_q099_ClockTest.ta*
                [ $DEBUG -eq 0 ] && [ -s molClock_PAUP_files.tgz ] && rm ./*_paup.block ./*.nex  ./*_clockTest.log ./*tre ./*clock.scores ./*critical_X2_val.R
                [ $DEBUG -eq 0 ] && rm list2concat Rplots.pdf header.tmp list2grep.tmp concat_nonRecomb_KdeFilt_cdnAlns_FTGTRG.ph
	    fi
	    
	    tar -czf concatenated_alignment_files.tgz concat_cdnAlns.fna concat_cdnAlns.fnainf
            [ -s concatenated_alignment_files.tgz ] && rm concat_cdnAlns.fna concat_cdnAlns.fnainf
	    [ $DEBUG -eq 0 ] && rm gene_trees2_concat_tree_RF_distances.tab ./*cdnAln_FTGTR.ph ./*cdnAln.fasta
	    [ $DEBUG -eq 0 ] && rm ../Rplots.pdf ../sorted*perc.tab sorted_aggregated_*tab ../all_*trees.tre

            cd $non_recomb_cdn_alns_dir
	    [ $DEBUG -eq 0 ] && tar -czf non_recombinant_kdeOK_codon_alignments.tgz ./*_cdnAln.fasta
	    [ -s non_recombinant_kdeOK_codon_alignments.tgz ] && rm ./*_cdnAln.fasta
	    tar -czf gene_trees_from_non_recombinant_kdeOK_codon_alignments.tgz ./*_cdnAln*.ph
	    [ $DEBUG -eq 0 ] && [ -s gene_trees_from_non_recombinant_kdeOK_codon_alignments.tgz ] && rm ./*_cdnAln*.ph

	    cd $top_dir
	    tar -czf codon_alignments.tgz ./*_cdnAln.fasta
            [ $DEBUG -eq 0 ] && [ -s codon_alignments.tgz ] && rm ./*_cdnAln.fasta clustalo.log
            tar -czf protein_alignments.tgz ./*.faaln
            [ $DEBUG -eq 0 ] && [ -s protein_alignments.tgz ] && rm ./*.faaln
	else
	    if [ $eval_clock -eq 1 ]
	    then
	        tar -czf IQT_molClock_PAUP_files.tgz ./*_paup.block ./*.nex  ./*_clockTest.log ./*tre ./*clock.scores ./*critical_X2_val.R \
	        $mol_clock_tab ${mol_clock_tab}sorted mol_clock_MGTRG_r_o_q099_ClockTest.ta*
                [ -s IQT_molClock_PAUP_files.tgz ] && rm ./*_paup.block ./*.nex ./*_clockTest.log ./*tre ./*clock.scores ./*critical_X2_val.R
                [ "$DEBUG" -eq "0" ] && rm list2concat Rplots.pdf header.tmp list2grep.tmp ./*cdnAln.ph 
            fi
	    
	    tar -czf concatenated_alignment_files.tgz concat_cdnAlns.fna concat_cdnAlns.fnainf
            [ -s concatenated_alignment_files.tgz ] && rm concat_cdnAlns.fna concat_cdnAlns.fnainf
	    [ $DEBUG -eq 0 ] && rm gene_trees2_concat_tree_RF_distances.tab ./*.ph ./*cdnAln.fasta ./*.log sorted_aggregated_support_values4loci_*.tab

            cd $non_recomb_cdn_alns_dir
	    [ $DEBUG -eq 0 ] && rm ./Rplots.pdf sorted_aggregated_*tab ./all_*trees.tre kde_*out
	    tar -czf non_recombinant_kdeOK_codon_alignments.tgz ./*_cdnAln.fasta
	    [ $DEBUG -eq 0 ] && [ -s non_recombinant_kdeOK_codon_alignments.tgz ] && rm ./*_cdnAln.fasta
	    tar -czf IQT_gene_trees_from_non_recombinant_kdeOK_codon_alignments.tgz ./*_cdnAln*.treefile
	    [ $DEBUG -eq 0 ] && [ -s IQT_gene_trees_from_non_recombinant_kdeOK_codon_alignments.tgz ] && rm ./*_cdnAln*.treefile
	    tar -czf IQT_gene_tree_logfiles.tgz ./*fasta.log
	    [ $DEBUG -eq 0 ] && [ -s IQT_gene_tree_logfiles.tgz ] && rm ./*fasta.log
	    
	    cd $top_dir
	    tar -czf codon_alignments.tgz ./*_cdnAln.fasta
            [ $DEBUG -eq 0 ] && [ -s codon_alignments.tgz ] && rm ./*_cdnAln.fasta clustalo.log
            tar -czf protein_alignments.tgz ./*.faaln
            [ $DEBUG -eq 0 ] && [ -s protein_alignments.tgz ] && rm ./*.faaln
	fi
    fi # if [ $runmode -eq 1 ]; then run phylo pipeline on DNA seqs

#----------------------------------------#
# >>> BLOCK 4.3: POPULATION GENETICS <<< #
#----------------------------------------#

   if [ $runmode -eq 2 ]
    then

        msg "" PROGR NC
        msg " >>>>>>>>>>>>>>> run descriptive DNA polymorphism statistics and neutrality tests <<<<<<<<<<<<<<< " PROGR YELLOW
        msg "" PROGR NC

        mkdir popGen && cd popGen
	popGen_dir=$(pwd)

        print_start_time && msg "# Moved into dir popGen ..." PROGR LBLUE

	ln -s ../*fasta .
	no_top_markers=$(find . -name "*.fasta" | wc -l)
	tmpf=$(find . -name "*.fasta" | head -1)
	no_seqs=$(grep -c '>' $tmpf)
	[ $DEBUG -eq 1 ] && msg "no_seqs:$no_seqs" DEBUG NC

        print_start_time && msg "# Will run descriptive DNA polymorphism statistics for $no_top_markers top markers. This will take some time ..." PROGR BLUE

	TajD_crit_vals=$(get_critical_TajD_values $no_seqs)
	TajD_l=$(echo $TajD_crit_vals|awk '{print $1}')
	TajD_u=$(echo $TajD_crit_vals|awk '{print $2}')

	FuLi_crit_vals=$(get_critical_FuLi_values $no_seqs)
	FuLi_l=$(echo "$FuLi_crit_vals"|awk '{print $1}')
	FuLi_u=$(echo "$FuLi_crit_vals"|awk '{print $2}')

	[ $DEBUG -eq 1 ] && msg "TajD_crit_vals:$TajD_crit_vals|TajD_l:$TajD_l|TajD_u:$TajD_u|FuLi_crit_vals:$FuLi_crit_vals|FuLi_l:$FuLi_l|FuLi_u:$FuLi_u" DEBUG NC

        print_start_time && msg "# converting $no_top_markers fasta files to nexus format ..." PROGR BLUE
        [ $DEBUG -eq 1 ] && msg " > convert_aln_format_batch_bp.pl fasta fasta nexus nex &> /dev/null" DEBUG NC
	#$distrodir/
	${distrodir}/convert_aln_format_batch_bp.pl fasta fasta nexus nex &> /dev/null

        print_start_time && msg "# Running popGen_summStats.pl ..." PROGR BLUE
	lmsg=" > popGen_summStats.pl -R 2 -n nex -f fasta -F fasta -H -r 100 -t $TajD_l -T $TajD_u -s $FuLi_l -S $FuLi_u &> popGen_summStats_hs100.log"
	[ $DEBUG -eq 1 ] && msg "$lmsg" DEBUG NC
	#$distrodir/popGen_summStats.pl -R 2 -n nex -f fasta -F fasta -H -r 100 -t $TajD_l -T $TajD_u -s $FuLi_l -S $FuLi_u &> popGen_summStats_hs100.log
	${distrodir}/popGen_summStats.pl -R 2 -n nex -f fasta -F fasta -H -r 100 -t $TajD_l -T $TajD_u -s $FuLi_l -S $FuLi_u &> popGen_summStats_hs100.log

	check_output polymorphism_descript_stats.tab $parent_PID | tee -a ${logdir}/get_phylomarkers_run_${dir_suffix}_${TIMESTAMP_SHORT}.log

	msg " >>> descriptive DNA polymorphism stats are found in: $popGen_dir ..." PROGR GREEN


	#>>> CLEANUP <<<#
	tar -czf clean_cdnAlns.tgz ./*_cdnAln_clean.fasta ./*.nex
	[ -s clean_cdnAlns.tgz ] && rm -rf ./*_cdnAln_clean.fasta ./*.nex ./paup.cmd ./popGen_summStats_*.log ./*_cdnAln.fasta

	cd $non_recomb_cdn_alns_dir
	tar -czf non_recombinant_kdeOK_codon_alignments.tgz ./*_cdnAln.fasta
	[ -s non_recombinant_kdeOK_codon_alignments.tgz ] && rm ./*_cdnAln.fasta ./all_*trees.tre ./Rplots.pdf
	tar -czf gene_trees_from_non_recombinant_kdeOK_codon_alignments.tgz ./*_cdnAln*.ph
	[ -s gene_trees_from_non_recombinant_kdeOK_codon_alignments.tgz ] && rm ./*_cdnAln*.ph

	cd $top_dir
	tar -czf codon_alignments.tgz ./*_cdnAln.fasta
        [ -s codon_alignments.tgz ] && rm ./*_cdnAln.fasta clustalo.log
        tar -czf protein_alignments.tgz ./*.faaln
        [ -s protein_alignments.tgz ] && rm ./*.faaln
    fi
fi # if [ "$mol_type" == "DNA"


#::::::::::::::::: END -t DNA RUNMODES :::::::::::::::::#


#======================================================#
# >>>  Block 5. -t PROTEIN GENE AND SPECIES TREES  <<< #
#======================================================#

if [ "$mol_type" == "PROT" ]
then
    cd non_recomb_FAA_alns
    non_recomb_FAA_alns_dir=$(pwd)

    print_start_time && msg "# working in dir $non_recomb_FAA_alns_dir ..." PROGR LBLUE
    print_start_time && msg "# estimating $no_non_recomb_alns_perm_test gene trees from non-recombinant sequences ..." PROGR LBLUE


    # 5.1 >>> estimate_IQT_gene_trees | estimate_FT_gene_trees
    if [ "$search_algorithm" == "F" ]
    then
	msg "" PROGR NC
	msg " >>>>>>>>>>>>>>> parallel FastTree runs to estimate gene trees <<<<<<<<<<<<<<< " PROGR YELLOW
	msg "" PROGR NC

        gene_tree_ext="ph"
	lmsg=" > running estimate_FT_gene_trees $mol_type $search_thoroughness ..."
        [ $DEBUG -eq 1 ] && msg "$lmsg" DEBUG NC
	estimate_FT_gene_trees $mol_type $search_thoroughness

	# 4.1.1 check that FT computed the expected gene trees
        no_gene_trees=$(find . -name "*.ph" | wc -l)
        [ $no_gene_trees -lt 1 ] && print_start_time && msg " >>> ERROR: There are no gene tree to work on in non_recomb_cdn_alns/. will exit now!" ERROR RED && exit 3

	# 4.1.2 generate computation-time and lnL stats
	compute_FT_gene_tree_stats $mol_type $search_thoroughness
    else
        gene_tree_ext="treefile"
	msg "" PROGR NC
	msg " >>>>>>>>>>>>>>> parallel IQ-TREE runs to estimate gene trees <<<<<<<<<<<<<<< " PROGR YELLOW
	msg "" PROGR NC

	estimate_IQT_gene_trees $mol_type $search_thoroughness $IQT_models

	# 4.1.1 check that IQT computed the expected gene trees
	no_gene_trees=$(ls ./*.treefile | wc -l)
        [ $no_gene_trees -lt 1 ] && print_start_time && msg " >>> ERROR: There are no gene tree to work on in non_recomb_cdn_alns/. will exit now!" ERROR RED && exit 3

	# 4.1.2 generate computation-time, lnL and best-model stats
	compute_IQT_gene_tree_stats $mol_type $search_thoroughness
    fi

    #remove trees with < 5 branches
    print_start_time && msg "# counting branches on $no_non_recomb_alns_perm_test gene trees ..." PROGR BLUE
    [ $DEBUG -eq 1 ] && msg " > count_tree_branches $gene_tree_ext no_tree_branches.list &> /dev/null" DEBUG NC

    count_tree_branches $gene_tree_ext no_tree_branches.list # &> /dev/null
    [ ! -s no_tree_branches.list ] && install_Rlibs_msg no_tree_branches.list ape

    check_output no_tree_branches.list $parent_PID | tee -a ${logdir}/get_phylomarkers_run_${dir_suffix}_${TIMESTAMP_SHORT}.log

    # remove trees with < 5 external branches (leaves)
     [ $DEBUG -eq 1 ] && msg " >  removing trees with < 5 external branches (leaves)" DEBUG NC

    no_tree_counter=0
    for phy in $(grep -v '^#Tree' no_tree_branches.list | awk -v min_no_ext_branches=$min_no_ext_branches 'BEGIN{FS="\t"; OFS="\t"}$7 < min_no_ext_branches' |cut -f1)
    do
         [ "$search_algorithm" == "F" ] && base=$(echo ${phy//_allFT\.ph/}) 
	 [ "$search_algorithm" == "I" ] && base=$(echo ${phy//\.treefile/})
	 print_start_time && msg " >>> will remove ${base}* because it has < 5 branches" WARNING LRED
	 rm ${base}*
	 let no_tree_counter++
    done

    msg " >>> WARNING: there are $no_tree_counter trees with < 1 internal branches (no real trees) that will be discarded ..." WARNING LRED

    msg "" PROGR NC
    msg " >>>>>>>>>>>>>>> filter gene trees for outliers with kdetrees test <<<<<<<<<<<<<<< " PROGR YELLOW
    msg "" PROGR NC


    if [ "$search_algorithm" == "F" ]
    then
        # 4.1 generate the all_GTRG_trees.tre holding all source trees, which is required by kdetrees
        #     Make a check for the existence of the file to interrupt the pipeline if something has gone wrong
        [ $DEBUG -eq 1 ] && msg " > cat ./*.ph > all_gene_trees.tre" DEBUG NC
        cat ./*.ph > all_gene_trees.tre
        check_output all_gene_trees.tre $parent_PID | tee -a ${logdir}/get_phylomarkers_run_${dir_suffix}_${TIMESTAMP_SHORT}.log
        [ ! -s all_gene_trees.tre ] && exit 3
    else
        # 4.1 generate the all_IQT_trees.tre holding all source trees, which is required by kdetrees
        #     Make a check for the existence of the file to interrupt the pipeline if something has gone wrong
        [ $DEBUG -eq 1 ] && msg " > cat ./*treefile > all_gene_trees.tre" DEBUG NC
	cat ./*.treefile > all_gene_trees.tre
        check_output all_gene_trees.tre $parent_PID | tee -a ${logdir}/get_phylomarkers_run_${dir_suffix}_${TIMESTAMP_SHORT}.log
        [ ! -s all_gene_trees.tre ] && exit 3
    fi

    # 5.2 run_kdetrees.R at desired stringency
    print_start_time && msg "# running kde test ..." PROGR BLUE
    [ $DEBUG -eq 1 ] && msg " > ${distrodir}/run_kdetrees.R ${gene_tree_ext} all_gene_trees.tre $kde_stringency &> /dev/null" DEBUG NC
    ${distrodir}/run_kdetrees.R ${gene_tree_ext} all_gene_trees.tre $kde_stringency &> /dev/null
    [ ! -s kde_dfr_file_all_gene_trees.tre.tab ] && install_Rlibs_msg kde_dfr_file_all_gene_trees.tre.tab kdetrees,ape
    check_output kde_dfr_file_all_gene_trees.tre.tab $parent_PID | tee -a ${logdir}/get_phylomarkers_run_${dir_suffix}_${TIMESTAMP_SHORT}.log

    # 5.3 mv outliers to kde_outliers
    no_kde_outliers=$(grep -c outlier kde_dfr_file_all_gene_trees.tre.tab)
    no_kde_ok=$(grep -v outlier kde_dfr_file_all_gene_trees.tre.tab|grep -vc '^file')

    # 5.4 Check how many cdnAlns passed the test and separate into two subirectories those passing and failing the test
    if [ $no_kde_outliers -gt 0 ]
    then
        print_start_time && msg "# making dir kde_outliers/ and moving $no_kde_outliers outlier files into it ..." PROGR BLUE
        mkdir kde_outliers
        [ $search_algorithm == "F" ] && for f in $(grep outlier kde_dfr_file_all_gene_trees.tre.tab|cut -f1|sed 's/_allFTlgG\.ph//'); do mv ${f}* kde_outliers; done
	[ $search_algorithm == "I" ] && for f in $(grep outlier kde_dfr_file_all_gene_trees.tre.tab|cut -f1|sed 's/\.treefile//'); do mv ${f}* kde_outliers; done
    else
        print_start_time && msg " >>> there are no kde-test outliers ..." PROGR GREEN
    fi

    if [ $no_kde_ok -gt 0 ]
    then
        print_start_time && msg "# making dir kde_ok/ and linking $no_kde_ok selected files into it ..." PROGR BLUE
        mkdir kde_ok
        cd kde_ok
        ln -s ../*.${gene_tree_ext} .

        print_start_time && msg "# labeling $no_kde_ok gene trees in dir kde_ok/ ..." PROGR BLUE
        [ $DEBUG -eq 1 ] && msg " > ${distrodir}/run_parallel_cmmds.pl ${gene_tree_ext} 'add_labels2tree.pl ../../../tree_labels.list $file' $n_cores &> /dev/null" DEBUG NC
        ${distrodir}/run_parallel_cmmds.pl ${gene_tree_ext} 'add_labels2tree.pl ../../../tree_labels.list $file' $n_cores &> /dev/null

	# remove symbolic links to cleanup kde_ok/
        for f in $(ls ./*.${gene_tree_ext}|grep -v "_ed\.${gene_tree_ext}"); do rm $f; done

        cd ..
    else
        print_start_time && msg "# ERROR There are $no_kde_ok gene trees producing non-significant kde-test results! Increase the actual -k $kde_stringency value. Will stop here!" ERROR RED
	exit 5
    fi

#----------------------------------#
# >>> BLOCK 5.2: PHYLOGENETICS <<< #
#----------------------------------#

    wkdir=$(pwd)

    msg "" PROGR NC
    msg " >>>>>>>>>>>>>>> filter gene trees by phylogenetic signal content <<<<<<<<<<<<<<< " PROGR YELLOW
    msg "" PROGR NC


    print_start_time && msg "# computing tree support values ..." PROGR BLUE
    [ $DEBUG -eq 1 ] && msg " > compute_suppValStas_and_RF-dist.R $wkdir 1 fasta ${gene_tree_ext} 1 &> /dev/null" DEBUG NC
    compute_suppValStas_and_RF-dist.R $wkdir 1 faaln ${gene_tree_ext} 1 &> /dev/null

    print_start_time && msg "# writing summary tables ..." PROGR BLUE
    min_supp_val_perc=${min_supp_val#0.}
    no_digits=${#min_supp_val_perc}
    [ $no_digits -eq 1 ] && min_supp_val_perc=${min_supp_val_perc}0


    # NOTE: IQ-TREE -alrt 1000 provides support values in 1-100 scale, not as 0-1 as FT!
    #	    therefore we convert IQT-bases SH-alrt values to 0-1 scale for consistency
    [ $search_algorithm == "I" ] && min_supp_val=$(echo "$min_supp_val * 100" | bc)

    awk -v min_supp_val=$min_supp_val '$2 >= min_supp_val' sorted_aggregated_support_values4loci.tab > sorted_aggregated_support_values4loci_ge${min_supp_val_perc}perc.tab
    check_output sorted_aggregated_support_values4loci_ge${min_supp_val_perc}perc.tab $parent_PID | \
    tee -a ${logdir}/get_phylomarkers_run_${dir_suffix}_${TIMESTAMP_SHORT}.log

    no_top_markers=$(perl -lne 'END{print $.}' sorted_aggregated_support_values4loci_ge${min_supp_val_perc}perc.tab)
    top_markers_dir="top_${no_top_markers}_markers_ge${min_supp_val_perc}perc"
    top_markers_tab=$(ls sorted_aggregated_support_values4loci_ge${min_supp_val_perc}perc.tab)

    # >>> 5.2 move top-ranking markers to $top_markers_dir
    print_start_time && msg "# making dir $top_markers_dir and moving $no_top_markers top markers into it ..." PROGR LBLUE
    mkdir $top_markers_dir && cd $top_markers_dir
    top_markers_dir=$(pwd)
    ln -s ../$top_markers_tab .
    for base in $(awk '{print $1}' $top_markers_tab|grep -v loci|sed 's/"//g'); do ln -s ../${base}* .; done

    [ $no_top_markers -lt 2 ] && print_start_time && msg " >>> Warning: There are less than 2 top markers. Relax your filtering thresholds. will exit now!" ERROR LRED && exit 3

    msg "" PROGR NC
    msg " >>>>>>>>>>>>>>> generate supermatrix from $no_top_markers concatenated, top-ranking alignments <<<<<<<<<<<<<<< " PROGR YELLOW
    msg "" PROGR NC

    # >>> 5.3 generate supermatrix (concatenated alignment)
    print_start_time && msg "# concatenating $no_top_markers top markers into supermatrix ..." PROGR BLUE
    [ $DEBUG -eq 1 ] && msg " > concat_alns faaln $parent_PID &> /dev/null" DEBUG NC
    concat_alns faaln $parent_PID &> /dev/null

    # >>> 5.4 remove uninformative sites from the concatenated alignment to speed up computation
    print_start_time && msg "# removing uninformative sites from concatenated alignment ..." PROGR BLUE
    [ $DEBUG -eq 1 ] && msg " > ${distrodir}/remove_uninformative_sites_from_aln.pl < concat_protAlns.faa > concat_protAlns.faainf" DEBUG NC
    ${distrodir}/remove_uninformative_sites_from_aln.pl < concat_protAlns.faa > concat_protAlns.faainf
    check_output concat_protAlns.faainf $parent_PID | tee -a ${logdir}/get_phylomarkers_run_${dir_suffix}_${TIMESTAMP_SHORT}.log

    [ ! -s concat_protAlns.faainf ] && print_start_time && msg " >>> ERROR: The expected file concat_protAlns.faainf was not produced! will exit now!" ERROR RED && exit 3


    #:::::::::::::::::::::::::::::::#
    # >>> FastTree species tree <<< #
    #:::::::::::::::::::::::::::::::#
    if [ "$search_algorithm" == "F" ]
    then
        msg "" PROGR NC
        msg " >>>>>>>>>>>>>>> FastTree run on $no_top_markers concatenated $mol_type alignments to estimate the species tree <<<<<<<<<<<<<<< " PROGR YELLOW
        msg "" PROGR NC

        # 5.4 run FasTree under the LG+G model
        print_start_time && msg "# running FastTree on $no_top_markers concatenated $mol_type alignments with $search_thoroughness thoroughness. This may take a while ..." PROGR BLUE

        if [ "$search_thoroughness" == "high" ]
        then
            $bindir/FastTree -quiet -lg -bionj -slow -slownni -gamma -mlacc 3 -spr $spr -sprlength $spr_length -log ${tree_prefix}_nonRecomb_KdeFilt_protAlns_FTlgG.log \
	    < concat_protAlns.faainf > ${tree_prefix}_nonRecomb_KdeFilt_protAlns_FTlgG.ph
        fi

        if [ "$search_thoroughness" == "medium" ]
        then
            $bindir/FastTree -quiet -lg -bionj -slownni -gamma -mlacc 2 -spr $spr -sprlength $spr_length -log ${tree_prefix}_nonRecomb_KdeFilt_protAlns_FTlgG.log \
	    < concat_protAlns.faainf > ${tree_prefix}_nonRecomb_KdeFilt_protAlns_FTlgG.ph
        fi

        if [ "$search_thoroughness" == "low" ]
        then
            $bindir/FastTree -quiet -lg -bionj -gamma -spr $spr -sprlength $spr_length -log ${tree_prefix}_nonRecomb_KdeFilt_protAlns_FTlgG.log \
	    < concat_protAlns.faainf > ${tree_prefix}_nonRecomb_KdeFilt_protAlns_FTlgG.ph
        fi

        if [ "$search_thoroughness" == "lowest" ]
        then
            $bindir/FastTree -quiet -lg -gamma -mlnni 4 -log ${tree_prefix}_nonRecomb_KdeFilt_protAlns_FTlgG.log < concat_protAlns.faainf > \
	    ${tree_prefix}_nonRecomb_KdeFilt_protAlns_FTlgG.ph
        fi

        check_output ${tree_prefix}_nonRecomb_KdeFilt_protAlns_FTlgG.ph $parent_PID | tee -a ${logdir}/get_phylomarkers_run_${dir_suffix}_${TIMESTAMP_SHORT}.log


       if [ -s "${tree_prefix}_nonRecomb_KdeFilt_protAlns_FTlgG.log" -a -s "${tree_prefix}_nonRecomb_KdeFilt_protAlns_FTlgG.ph" ]
       then
         # lnL=$(grep ML_Lengths2 "${tree_prefix}_nonRecomb_KdeFilt_protAlns_FTlgG.log" | grep TreeLogLk | sed 's/TreeLogLk[[:space:]]ML_Lengths2[[:space:]]//')
           lnL=$(grep '^Gamma20LogLk' "${tree_prefix}_nonRecomb_KdeFilt_protAlns_FTlgG.log" |awk '{print $2}')
           msg " >>> lnL for ${tree_prefix}_nonRecomb_KdeFilt_protAlns_FTlgG.ph = $lnL" PROGR GREEN
       else
           msg " >>> WARNING: ${tree_prefix}_nonRecomb_KdeFilt_protAlns_FTlgG.log could not be produced!" WARNING LRED
       fi

        print_start_time && msg "# Adding labels back to tree ..." PROGR LBLUE
        lmsg=" > add_labels2tree.pl ${tree_labels_dir}/tree_labels.list ${tree_prefix}_nonRecomb_KdeFilt_protAlns_FTlgG.ph &> /dev/null"
        [ $DEBUG -eq 1 ] && msg "$lmsg" DEBUG NC
        add_labels2tree.pl ${tree_labels_dir}/tree_labels.list ${tree_prefix}_nonRecomb_KdeFilt_protAlns_FTlgG.ph &> /dev/null

        [ -s ${tree_prefix}_nonRecomb_KdeFilt_protAlns_FTlgG_ed.ph ] && \
        mv ${tree_prefix}_nonRecomb_KdeFilt_protAlns_FTlgG_ed.ph ${tree_prefix}_${no_top_markers}nonRecomb_KdeFilt_protAlns_FTlgG.spTree
        mv ${tree_prefix}_nonRecomb_KdeFilt_protAlns_FTlgG.ph ${tree_prefix}_${no_top_markers}nonRecomb_KdeFilt_protAlns_FTlgG_numbered.tre

        check_output ${tree_prefix}_${no_top_markers}nonRecomb_KdeFilt_protAlns_FTlgG.spTree $parent_PID | tee -a ${logdir}/get_phylomarkers_run_${dir_suffix}_${TIMESTAMP_SHORT}.log

       msg " >>> found in dir $top_markers_dir ..." PROGR GREEN

        [ $DEBUG -eq 1 ] && msg " > compute_suppValStas_and_RF-dist.R $top_markers_dir 2 faaln ph 1 &> /dev/null" DEBUG NC
        ${distrodir}/compute_suppValStas_and_RF-dist.R $top_markers_dir 2 faaln ph 1 &> /dev/null
   fi # [ "$search_algorithm" == "F" ]


    #::::::::::::::::::::::::::::::#
    # >>> IQ-TREE species tree <<< #
    #::::::::::::::::::::::::::::::#

   if [ "$search_algorithm" == "I" ]
   then
       # 5.5 run IQ-tree in addition to FastTree, if requested
       msg "" PROGR NC
       msg " >>>>>>>>>>>>>>> IQ-TREE + ModelFinder run on $no_top_markers concatenated $mol_type alignments to estimate the species tree <<<<<<<<<<<<<<< " PROGR YELLOW
       msg "" PROGR NC

       #wkdir=$(pwd)
       
       print_start_time && msg "# running ModelFinder on the $no_top_markers concatenated $mol_type alignments with $IQT_models. This will take a while ..." PROGR BLUE

       iqtree -s concat_protAlns.faainf -st PROT -mset "$IQT_models" -m MF -nt AUTO -fast &> /dev/null

       check_output concat_protAlns.faainf.log $parent_PID | tee -a ${logdir}/get_phylomarkers_run_${dir_suffix}_${TIMESTAMP_SHORT}.log

       best_model=$(grep '^Best-fit model' concat_protAlns.faainf.log | cut -d' ' -f 3)
       msg " >>> Best-fit model: ${best_model} ..." PROGR GREEN

       mkdir iqtree_abayes && cd iqtree_abayes
       ln -s ../concat_protAlns.faainf .

       if [ "$search_thoroughness" == "high" ]
       then
    	  lmsg="# will launch $nrep_IQT_searches independent IQ-TREE searches on the supermatrix with best model ${best_model} -abayes -bb 1000! 
	            This will take a while ..."
    	  print_start_time && msg "$lmsg" PROGR BLUE

    	  # run nrep_IQT_searches IQ-TREE searches under the best-fit model found
    	  for ((rep=1;rep<=nrep_IQT_searches;rep++))
    	  do
    	     lmsg=" > iqtree -s concat_protAlns.faainf -st PROT -m $best_model -abayes -bb 1000 -nt AUTO -pre abayes_run${rep} &> /dev/null"
    	     print_start_time && msg "$lmsg" PROGR LBLUE

    	      iqtree -s concat_protAlns.faainf -st PROT -m $best_model -abayes -bb 1000 -nt AUTO -pre abayes_run${rep} &> /dev/null
    	  done

    	  grep '^BEST SCORE' ./*log | sed 's#./##' | sort -nrk5 > sorted_IQ-TREE_searches.out

    	  check_output sorted_IQ-TREE_searches.out $parent_PID | tee -a ${logdir}/get_phylomarkers_run_${dir_suffix}_${TIMESTAMP_SHORT}.log
    	  best_search=$(head -1 sorted_IQ-TREE_searches.out)
    	  best_search_base_name=$(head -1 sorted_IQ-TREE_searches.out | cut -d\. -f 1)

    	  msg "# >>> Best IQ-TREE run was: $best_search ..." PROGR GREEN

    	  best_tree_file=${tree_prefix}_${best_search_base_name}_nonRecomb_KdeFilt_${no_top_markers}concat_protAlns_iqtree_${best_model}.spTree
	  numbered_nwk=${tree_prefix}_${best_search_base_name}_nonRecomb_KdeFilt_${no_top_markers}concat_protAlns_iqtree_${best_model}_numbered.nwk
    	  cp ${best_search_base_name}.treefile $best_tree_file
	  cp ${best_search_base_name}.treefile $numbered_nwk

    	  print_start_time && msg "# Adding labels back to ${best_tree_file} ..." PROGR BLUE
   	  [ $DEBUG -eq 1 ] && msg " > ${distrodir}/add_labels2tree.pl ${tree_labels_dir}/tree_labels.list ${best_tree_file} &> /dev/null" DEBUG NC
   	  ${distrodir}/add_labels2tree.pl ${tree_labels_dir}/tree_labels.list $best_tree_file &> /dev/null
	  
	  sp_tree=$(ls ./*ed.spTree)

    	  check_output "$sp_tree" $parent_PID | tee -a ${logdir}/get_phylomarkers_run_${dir_suffix}_${TIMESTAMP_SHORT}.log
    	  cp $sp_tree $numbered_nwk $top_markers_dir
    	  cd $top_markers_dir
    	  rm -rf iqtree_abayes concat_protAlns.faainf.treefile concat_protAlns.faainf.uniqueseq.phy ./*ckp.gz
       else
    	  print_start_time && msg "# running IQ-tree on the concatenated alignment with best model ${best_model} -abayes -bb 1000. This will take a while ..." PROGR BLUE

    	  print_start_time && msg "# running: iqtree -s concat_protAlns.faainf -st PROT -m $best_model -abayes -bb 1000 -nt AUTO -pre iqtree_abayes &> /dev/null  ..." PROGR BLUE
    	  iqtree -s concat_protAlns.faainf -st PROT -m $best_model -abayes -bb 1000 -nt AUTO -pre iqtree_abayes &> /dev/null

    	  best_tree_file=${tree_prefix}_nonRecomb_KdeFilt_${no_top_markers}concat_protAlns_iqtree_${best_model}.spTree
	  numbered_nwk=${tree_prefix}_nonRecomb_KdeFilt_${no_top_markers}concat_protAlns_iqtree_${best_model}_numbered.nwk
    	  cp iqtree_abayes.treefile ${best_tree_file}
	  cp iqtree_abayes.treefile $numbered_nwk

    	  print_start_time && msg "# Adding labels back to ${best_tree_file} ..." PROGR BLUE
   	  [ $DEBUG -eq 1 ] && msg " > ${distrodir}/add_labels2tree.pl ${tree_labels_dir}/tree_labels.list ${best_tree_file} &> /dev/null" DEBUG NC
   	  ${distrodir}/add_labels2tree.pl ${tree_labels_dir}/tree_labels.list ${best_tree_file} &> /dev/null
	  
	  sp_tree=$(ls ./*ed.spTree)

    	  check_output $sp_tree $parent_PID | tee -a ${logdir}/get_phylomarkers_run_${dir_suffix}_${TIMESTAMP_SHORT}.log
    	  cp $sp_tree $numbered_nwk $top_markers_dir
    	  cd $top_markers_dir
    	  rm -rf iqtree_abayes concat_protAlns.faainf.treefile concat_protAlns.faainf.uniqueseq.phy ./*ckp.gz
       fi
    fi # if [ "$search_algorithm" == "I" ]

    # >>> 5.9 CLEANUP <<< #

    if [ $search_algorithm == "F" ]
    then
        [ $DEBUG -eq 0 ] && rm list2concat Rplots.pdf sorted*perc.tab ./*.ph ./*faaln ./*log top*tab
        tar -czf concatenated_alignment_files.tgz concat_protAlns.faa concat_protAlns.faainf
        [ -s concatenated_alignment_files.tgz ] && rm concat_protAlns.faa concat_protAlns.faainf 
        [ $DEBUG -eq 0 ] && rm ../sorted*perc.tab sorted_aggregated_*tab | \
        tee -a ${logdir}/get_phylomarkers_run_${dir_suffix}_${TIMESTAMP_SHORT}.log

        cd $non_recomb_FAA_alns_dir
        tar -czf non_recomb_kdeOK_FAA_alignments.tgz ./*_cluo.faaln
        [ -s non_recomb_kdeOK_FAA_alignments.tgz ] && rm ./*_cluo.faaln
        tar -czf non_recomb_kdeOK_prot_trees.tgz ./*_cluo_*.ph ./*.log
        [ -s non_recomb_kdeOK_prot_trees.tgz ] && rm ./*_cluo_*.ph ./*.log top100* kde*out
	rm Rplots.pdf no_tree_branches.list

        cd $top_dir
        tar -czf codon_alignments.tgz ./*_cdnAln.fasta
        [ -s codon_alignments.tgz ] && rm ./*_cdnAln.fasta clustalo.log
        tar -czf protein_alignments.tgz ./*.faaln
        [ -s protein_alignments.tgz ] && rm ./*.faaln
    else
        [ $DEBUG -eq 0 ] && rm list2concat sorted*perc.tab
        rm ./*faaln ./*faaln.treefile ./*faaln.log 

        tar -czf concatenated_alignment_files.tgz concat_protAlns.faa concat_protAlns.faainf
        [ -s concatenated_alignment_files.tgz ] && rm concat_protAlns.faa concat_protAlns.faainf
        [ $DEBUG -eq 0 ] && rm  ../Rplots.pdf ../sorted*perc.tab 

        cd $non_recomb_FAA_alns_dir
        tar -czf non_recomb_kdeOK_FAA_alignments.tgz ./*_cluo.faaln
        [ -s non_recomb_kdeOK_FAA_alignments.tgz ] && rm ./*_cluo.faaln
        tar -czf non_recomb_kdeOK_prot_trees.tgz ./*faaln.treefile ./*.faaln.log all_gene_trees.tre
        [ -s non_recomb_kdeOK_prot_trees.tgz ] && rm ./*faaln.treefile all_gene_trees.tre ./*.faaln.log kde*.out

        cd $top_dir
        tar -czf codon_alignments.tgz ./*_cdnAln.fasta
        [ -s codon_alignments.tgz ] && rm ./*_cdnAln.fasta clustalo.log
        tar -czf protein_alignments.tgz ./*.faaln
        [ -s protein_alignments.tgz ] && rm ./*.faaln
    fi

fi # [ "$mol_type" == "PROT" ]

# compute the elapsed time since the script was fired
end_time=$(date +%s)
secs=$((end_time-start_time))

#printf '%dh:%dm:%ds\n' $(($secs/3600)) $(($secs%3600/60)) $(($secs%60))
msg "" PROGR NC
msg " >>> Total runtime of $progname:" PROGR LBLUE
printf '%dh:%dm:%ds\n' "$((secs/3600))" "$((secs%3600/60))" "$((secs%60))" | tee -a ${logdir}/get_phylomarkers_run_${dir_suffix}_${TIMESTAMP_SHORT}.log
echo

cat <<REF

* PROVISIONAL CITATION:

If you find the code useful for your academic work, please use the following citation:

Pablo Vinuesa, Luz-Edith Ochoa-Sanchez and Bruno Contreras-Moreira 2018.
GET_PHYLOMARKERS, a pipeline to select optimal phylogenetic markers for phylogenomics
and inference of pan-genome phylogenies: identification of cryptic species in
the Stenotrophomonas maltophilia complex.
Available at https://github.com/vinuesa/get_phylomarkers

The manuscript was submitted to Frontiers in Microbiology on January 15th 2018,
to be considered for its publication in the Research Topic on "Microbial Taxonomy, Phylogeny and Biodiversity"
http://journal.frontiersin.org/researchtopic/5493/microbial-taxonomy-phylogeny-and-biodiversity

* NOTES:
  1. The links to the the corresponding manuscript will be provided here
      as soon as it is available at bioRxiv, and latter, to the paper.

  2. If you encounter problems or bugs while running the pipeline
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


