#!/usr/bin/env bash

#: AUTHOR: Pablo Vinuesa; CCG-UNAM, Mexico. http://www.ccg.unam.mx/~vinuesa/
#: $0; run_molecClock_test_paup.sh (up to v. 0.3) -> run_molecClock_test_jmodeltest2_paup.sh

#: See print_help function for documentation

# PERFORMING MOLECULAR CLOCK TEST IN PAUP
# http://andrelevy.net/bioinfo/paup/molecular.clock.nex
# my PAUP* notes for LCG

# Common PAUP analysis commands by Peter Unmack
# see: http://peter.unmack.net/molecular/programs/paup.command.blocks.html

#--------------------------------------------------------------------------

progname=$(basename "$0")
VERSION='0.8_29May17' # v0.8_29May1: improved the regex to remove the symbol=""; part, which worked on Linux but not on MacOSX
                      # v0.7_25May17 1. added code to remove the symbols statement from PAUP's data block, 
                      #     as it seems to conflict with predefined DNA state symbol (on buluc!)
		      #     2. fixed regexes in run_paup_clock_test_with_user_tree that capture the -lnL values
                      # v0.6_3May17' Major upgrade version: added -R 1, which calles the function described in 1 below.
                     #    1. run_paup_clock_test_with_user_tree() to run without jmodeltest2, which is now -R 2
                     #    by using pre-computed userteres and passing either -M GTR|TrN|HKY|K2P|F81 as base models to adjust
		     #    model params on the user tree, with gamma correction of rates. 
		     #    2. Also uses -r to pass rootmethod midpoint|outgroup, levaing midpoint as default! much better
		     #   TODO: The logic of run_paup_clock_test() is not good, as it runs the for loop within the function;
		     #        Recode that function as in run_paup_clock_test_with_user_tree() and add the rootmethod option to itW
            # 0.5_28April17set default no_subst_schemes for jmodeltest2 to the minimum=3, 
                        #      used as default by jmodeltest2 to speed up genomic analyses.
	                #      In jmodeltest2 changed -S BEST to -t BIONJ to speed up!
            # v0.4, July 21st, 2014; Major re-write: changed modeltest3.7 for jmodeltest-2.1.5-20140405!
            #        This has many advantages: 1) parallel execution; 2) run using different schemes of substitution types
	    #        3) Get the PhyML tree under best-fit model and use that one for the mol. clock test, instead of
	    #           the original NJ tree used in previous versoins.
            # v0.3, July 19th, 2014; Fully integrated with run_amplicon_evalutation_pipeline.sh
            # v0.2, July 19th, 2014; Major upgrade: added get_opts, various subs and R code to compute X2 stats
            # 0.1,  July 2013; first prototype. Depended on usr provied df and crit X2 val
#-------------------------------------------------------------------------

function check_dependencies()
{
    check_dep_flag=$1
  
    # 1) check other binaries are in PATH
    echo
    echo "# >>> checking binaries and scripts in \$PATH..."
    for programname in R perl add_nos2fasta_header.pl convert_aln_format_batch_bp.pl
    do
       #if which $programname >/dev/null; then <== avoid which
       # see: http://stackoverflow.com/questions/592620/check-if-a-program-exists-from-a-bash-script
       bin=$(type -P $programname)
       if [ -z "$bin" ]; then
          echo
          echo " ERROR: $programname not in place!"
          echo " ... you will need to install \"$programname\" first or include it in \$PATH"
          echo " ... exiting"
          exit 1
       else
           echo " -> OK: $programname is in place ..."  
       fi
     done
   
   # 1.1 check jmodeltest2 is in place
   modeltest_home=$(echo "$JMODELTEST_HOME")
   if(( $? ))
   then
         echo "# jmodeltest2 seems not be in place. Check  echo \$JMODELTEST_HOME"
   else
         echo -n " -> OK: jmodeltest2_home=$modeltest_home"
   fi

   # 2) check all required perl modules are in @INC  
   echo
   echo "# >>> checking Perl modules ..."
   #for mod in Bio::SeqIO
   #do
      perldoc "Bio::SeqIO" > module_check.tmp
      grep '^No docum' module_check.tmp
      if (( $? ))
      then
          echo " -> OK: Perl module $mod is in place ..."  
      else
          echo '----------------------------------------------------------------------------------'
	  echo " >>> ERROR: Perl module $mod is was not found in @INC"
	  echo "      you will need to install it first by issuing the command: sudo cpan -i $mod"
	  echo "      Will exit now!"
          echo '----------------------------------------------------------------------------------'
	  exit 1
      fi
   #done

   echo
   echo '>>>>>>>>>>>>>>>>>>>>> DEPENDENCIES CHECK RESULT: OK <<<<<<<<<<<<<<<<<<<<<<<<<<'
   echo
  
  [ -s module_check.tmp ] && rm module_check.tmp

  [ "$check_dep_flag" -eq 1 ] && exit 0
}
#-------------------------------------------------------------------------

function set_script_environment()
{
    if [[ "$OSTYPE" == "linux-gnu" ]]
    then
         scriptdir=$(readlink -f "${BASH_SOURCE[0]}")
    	 distrodir=$(dirname "$scriptdir") #echo "scriptdir: $scriptdir|basedir:$distrodir|OSTYPE:$OSTYPE"
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

function set_bindirs()
{  
    # receives: $bindir $homebinpathflag
    bindir=$1

    not_in_path=0

    bins=( paup )

    for prog in "${bins[@]}" 
    do
       bin=$(type -P "$prog")
       if [ -z "$bin" ]
       then
          echo
          printf "${RED}# $prog not found in \$PATH ... ${NC}\n"
	        not_in_path=1
       fi	  
   done	  
 
   if [ $not_in_path -eq 1 ]
   then
   	   printf "${CYAN} updating PATH=$PATH:$bindir ${NC}\n"
   	   #export PATH=$PATH:$bindir # append $HOME/bin to $PATH, (at the end, to not interfere with the system PATH)  
	     # we do not export, so that this PATH update lasts only for the run of the script, 
	     # avoiding a longer alteration of $ENV; by appending to the end of PATH, no user-preferences should be altered 
	     PATH=$PATH:$bindir # append $HOME/bin to $PATH, (at the end, to not interfere with the system PATH)
   fi
   #echo $setbindir_flag
}
#----------------------------------------------------------------------------------------- 

function check_no_seqs_in_fasta()
{
    # find and quantify the classes of alignments by sequence numbers they conatin
    # sort aln_classes by decreasing no. of containing sequences and 
    # get the most frequent class to filter out those amps that do not match this no.
   
    fasta_file_ext=$1    

    no_seqNo_classes=$(grep -c '>' *."${fasta_file_ext}" | cut -d: -f2 | sort | uniq -c | wc -l)
    modal_no_seqs_in_amps=$(grep -c '>' *."${fasta_file_ext}" | cut -d: -f2 | sort | uniq -c | sort -nrk1 | awk '{print $2}' | head -1)
    
    if [ "$no_seqNo_classes" -ne 1 ]
    then
        echo 
	echo '-------------------------------------------------------------------------------------------------------------------'
	echo "# >>> function check_no_seqs_in_fasta WARNING: found $no_seqNo_classes sequence no. classes among fasta files"
	echo "   the fasta files with extension ${fasta_file_ext} have different nos. of sequences!"
	echo "   The modal number of sequences in alignments is = $modal_no_seqs_in_amps"
	echo "   Will set CHECK_NO_SEQS_IN_NEXUS = 1 ..."
	CHECK_NO_SEQS_IN_NEXUS=1
    else
        echo 
	echo "# >>> function check_no_seqs_in_fasta found that all *.${fasta_file_ext} files have $modal_no_seqs_in_amps sequences ..."

        no_seq=$modal_no_seqs_in_amps
	df=$(($modal_no_seqs_in_amps - 2))
	CHECK_NO_SEQS_IN_NEXUS=0

        echo "#  function check_no_seqs_in_fasta will set no_seq=$no_seq and df=$df"
    fi
}
#-------------------------------------------------------------------------

function check_is_numbered_fasta_file()
{
    fasta_file=$1
    
    fasta_id=$(head -1 "$fasta_file" | cut -d' ' -f1 | sed 's/>//')
    
    # http://stackoverflow.com/questions/806906/how-do-i-test-if-a-variable-is-a-number-in-bash
    re='^[0-9]+$'
    if ! [[ $fasta_id =~ $re ]]
    then
        echo "# Warning of function check_is_numbered_fasta_file: $fasta_file has not numbered header!" >&2
	echo "#  Will run add_nos2fasta_header.pl! and process the file" >&2
	add_nos2fasta_header.pl "$fasta_file"
	if [ -s my_"${fasta_file}" ]
	then
	    awk '{print $1}' my_"${fasta_file}" > ed
	    mv ed "my_${fasta_file}"
	    [ -s "my_${fasta_file}" ] && mv "${fasta_file}" src_fasta_files
	fi
    fi
}
#-------------------------------------------------------------------------

function run_X2_stats()
{
    # En R caluclamos el valor critico usando nivel de certidumbre (0.95 o 0.99) y df) asi:
    # qchisq(0.95, df); 
    # rechazamos H0 si X2 observada > valor critico
    # La p la calculamos asi: 1 - pchisq(X2,df)
       # $nex_basename $q $df $lnL_unconstr $lnL_clock
       nex_basename=$1
       q=$2
       df=$3
       lnL_unconstr=$4
       lnL_clock=$5
       
       crit_X2_val_file=${nex_basename}_critical_X2_val.txt
       p_ChiSq_test_val_file=${nex_basename}_p_ChiSq_test_val.txt
       LRT_file=${nex_basename}_LRT.txt
 
       R --no-save -q <<RCMD &> ${nex_basename}_compute_critical_X2_val.R

       LRT <- 2*($lnL_unconstr - $lnL_clock)
       
       sink(file="$LRT_file")
       print(LRT)
       sink()

       sink(file="$crit_X2_val_file")
       print(round(qchisq($q, $df)))
       sink()

       sink(file="$p_ChiSq_test_val_file")
       print( 1 - pchisq(LRT, $df))
       sink()
RCMD

    if [ -s "$crit_X2_val_file" ]
    then
	# need to get rid of R's vector notation in file [1] 54
	cut -d' ' -f2 "$crit_X2_val_file" > "${crit_X2_val_file}ed"
	mv "${crit_X2_val_file}ed" "$crit_X2_val_file"
    else
        echo
	echo "# ERROR: function compute_critical_X2_val did not reutrn file $crit_X2_val_file"
	echo "# ... will exit now!"
	exit 1
    fi

    if [ -s "$p_ChiSq_test_val_file" ]
    then
	# need to get rid of R's vector notation in file [1] 1.571452e-07
	cut -d' ' -f2 "$p_ChiSq_test_val_file" > "${p_ChiSq_test_val_file}ed"
	mv "${p_ChiSq_test_val_file}ed" "$p_ChiSq_test_val_file"
    else
        echo
	echo "# ERROR: function compute_critical_X2_val did not reutrn file $p_ChiSq_test_val_file"
	echo "# ... will exit now!"
	exit 2
    fi

    if [ -s "$LRT_file" ]
    then
	# need to get rid of R's vector notation in file [1] 1.571452e-07
	cut -d' ' -f2 "$LRT_file" > "${LRT_file}ed"
	mv "${LRT_file}ed" "$LRT_file"
    else
        echo
	echo "# ERROR: function compute_critical_X2_val did not reutrn file $LRT_file"
	echo "# ... will exit now!"
	exit 3
    fi
}
#-------------------------------------------------------------------------

function newick2nexus()
{
    # this function converts a newick-formatted tree
    # to a simple nexus-formatted one
    # NO translation block included
    intree=$1
    intree_string=$(cat $intree)
    outtree=${intree%.*}.tre
    
    [ -s $outtree ] && rm $outtree
    

    # convert newick tree to simple Nexus
    echo -e "#NEXUS\nBEGIN TREES;" >> $outtree
    echo "tree PAUP_1 = [&U] $intree_string" >> $outtree
    echo "END;" >> $outtree
}
#-------------------------------------------------------------------------

function run_paup_clock_test_with_user_tree()
{
    nexus=$1
    utree=$2
    basemod=$3
    res_table=$4
    root_method=$5
    outgroup=$6

   [ "$basemod" == "GTR" ] && lset_line1='nst=6 rmat=est rates=gamma shape=est' && lset_line2='nst=6 rmat=prev rates=gamma shape=prev'
   [ "$basemod" == "TrN" ] && lset_line1='nst=6 rmat=est rclass=(abaaea) rates=gamma shape=est' && lset_line2='nst=6 rmat=prev rclass=prev rates=gamma shape=prev'
   [ "$basemod" == "HKY" ] && lset_line1='nst=2 trat=est rates=gamma shape=est' && lset_line2='nst=2 trat=prev rates=gamma shape=prev'
   [ "$basemod" == "K2P" ] && lset_line1='nst=2 trat=est rates=gamma shape=est basefreq=equal' && lset_line2='nst=2 trat=prev rates=gamma shape=prev basefreq=equal'
   [ "$basemod" == "F81" ] && lset_line1='nst=1 rates=gamma shape=est basefreq=eq' && lset_line2='nst=1 rates=gamma shape=prev basefreq=eq'

   nex_basename=${nexus%.*}
   paupcmdfile=${nexus%.*}_paup.block
   logfile=${nexus%.*}_clockTest.log;
  
   scorefile_unconstr=${nexus%.*}_M${basemod}_r${root_method}_o${outgroup}_unconstrained_clock.scores
   scorefile_constr=${nexus%.*}_M${basemod}_r${root_method}_o${outgroup}_constrained_clock.scores
   
# 1. write the paup commands file for each nexus input file

if [ "$root_method" == "outgroup" ]
then
cat << PAUPBLOCK > $paupcmdfile
#nexus;
begin paup;
  set autoclose=y warntree=no warnreset=no;
  log file=$logfile start replace;
  outgroup $outgroup;
  GetTrees file=$utree StoreBrLens=Y StoreTreeWts=Y;    
  
  set criterion=likelihood;
  lset clock=no $lset_line1;
  lscores /scorefile=$scorefile_unconstr replace=yes;

  roottrees rootmethod=outgroup;
  lset clock=yes $lset_line2;
  lscores /scorefile=$scorefile_constr replace=yes;
  tstatus;
  log stop;
  end;
  quit;

PAUPBLOCK

fi

if [ "$root_method" == "midpoint" ]
then
cat << PAUPBLOCK > $paupcmdfile
#nexus;
begin paup;
  set autoclose=y warntree=no warnreset=no;
  log file=$logfile start replace;
  GetTrees file=$utree StoreBrLens=Y StoreTreeWts=Y;    
  
  set criterion=likelihood;
  lset clock=no $lset_line1;
  lscores /scorefile=$scorefile_unconstr replace=yes;

  roottrees rootmethod=midpoint;
  lset clock=yes $lset_line2;
  lsc /scorefile=$scorefile_constr replace=yes; 
  tstatus;
  log stop;
  end;
  quit;
  
PAUPBLOCK

fi
   
   # remove the symbols statement from PAUP's data block, as it seems to conflict with predefined DNA state symbol
   perl -pe 'if(/^format /){ s/\h+symbols=.*?;/;/}' $nexus > ${nexus}ed && mv ${nexus}ed $nexus
   [ $DEBUG -eq 1 ] && head $nexus

   # 2. run paup* with the file-specific cmd file
   [ $DEBUG -eq 1 ] && echo "# running: cat $paupcmdfile | paup -n $nexus &> /dev/null"
   cat $paupcmdfile | paup -n $nexus &> /dev/null
   
   # get the lnL values for the unconstrained and constrained trees: watch out the specific regexes
   lnL_unconstr=$(egrep -A 3 '^Tree' $logfile | egrep '^-ln' | head -1 | perl -pe 's/-ln L\s+//' )
   lnL_clock=$(egrep -A 3 '^Tree' $logfile | egrep '^-ln' | tail -1 | perl -pe 's/-ln L\s+//')
   
   # see: https://www.shell-tips.com/2010/06/14/performing-math-calculation-in-bash/
   #ChiSq=$(echo "2*($lnL_unconstr - $lnL_clock)" | bc)
   #ChiSqRound=$(echo $ChiSq | cut -d\. -f1 )
   
   if [ $CHECK_NO_SEQS_IN_NEXUS -eq 1 ]
   then
        # dimensions ntax=47 nchar=871;
	no_seq=$(grep 'dimensions ntax=' $nexus | cut -d\= -f2 | cut -d' ' -f1)
	df=$(($no_seq - 2))
   fi
   
   # 3. run R for the ChiSQ stats
   [ $DEBUG -eq 1 ] && echo "# running: run_X2_stats $nex_basename $q $df $lnL_unconstr $lnL_clock"
   run_X2_stats $nex_basename $q $df $lnL_unconstr $lnL_clock

   critical_X2_val=$(cat ${nex_basename}_critical_X2_val.txt)
   p_val=$(cat ${nex_basename}_p_ChiSq_test_val.txt)
   LRT=$(cat ${nex_basename}_LRT.txt)
   
   LRTround=$(echo $LRT | cut -d\. -f1)
   
   # 4. calc the mol_clock flag
   if [ $LRTround -ge $critical_X2_val ]
   then
      mol_clock='NO'
   else
      mol_clock='YES' 
   fi     
   
   # 5. concatenate to the results table. The header was written by the main program, befor calling this funcition
   echo -e "$nexus\t$lnL_unconstr\t$lnL_clock\t$LRT\t$critical_X2_val\t$df\t$p_val\t$mol_clock" >> $res_table
}
#-------------------------------------------------------------------------

function print_help()
{
   cat << HELP

 $progname v.$VERSION usage:

 Ia)  REQUIRED arguments:  
      -f input nexus file alignment

 Ib)  Conditionally REQUIRED arguments:
      -o outgroup OTU number (use numbered fasta headers, as in my_ files) to be used as outgroup for tree rooting 
            required if -r outgroup

 IIa) OPTIONAL arguments for data parsing and program execution modes:  
      -b basename for output table holding the analysis results                       [default: $basename_output_table]
      -m <string> modeltest output file suffix                                        [default: $modeltest_outfile_suffix]
      -M <string> base Model for -R 1 (use one of: GTR|TrN|HKY|K2P|F81)               [default: $base_mod]
      -q <real> quantile for computing critical and p-value of X2 square distribution [default $q]
      -r <string> rooting method: [outgroup|midpoint]                                 [default: $rooting_method]
      -t <string> extension name of pre-computed user-trees
      -D print DEBUGGING messages                                                     [flag, def $DEBUG]
      -V print verbose output to screen                                               [flag, def $VERBOSE]
      -K run check_dependencies                                                       [flag, def $CHECK_DEPENDENCIES]
      -k check number of sequences in fasta                                           [flag, def $CHECK_NO_SEQS_IN_NEXUS]

 IIb) OPTIONAL arguments for running jmodeltest2 
      -g no. of categories to discretizize the gamma distribution                    [def $no_cat]
      -I run also models with prop. inv. sites                                       [flag, def $pInv]
         (jmodeltest is called with -f -g 4 to use also uneq_freq and +G models
	  Only -I is ignored by default; can be added with -I)
      -s no. of substitution schemes for jmodeltest2 <3,5,7,11,203>              [default $no_subst_schemes]
      -S search algorithm <NNI, SPR, BEST>                                       [default $search_alg]
        
 EXAMPLE invocation lines
    1) standard run:
    $progname -R 1 -M HKY -b test_clock -f nexus_file
    $progname -R 1 -M HKY -r outgroup -o 23 -b test_clock -f nexus_file -q 0.95 
   
    2) default: [ -q 0.99 -s 5 -g 4]
           $progname -R 2 -b test_with_fastas -e dna_amp -s 3 -g 5 -I   
    
 AIM:
    Used to run global molec. clock test in PAUP*  (-R 1) using as input files the outfile generated by jmodeltest2. 
    The script takes the newick string from the jmodelstes2.out files, converts them to nexus-formatted tree and
    writes a paup block to run lscores on this tree without clock enforcement, under the best-fit model selected under AIC.
    Modeltest runs are launched for -R 1 -a -R 2. df, no-seqs and LRT stats are automatically computed by the script, 
    which uses R code for the ChiSQ stats.

 INPUT:
      aligned fasta files; providing their extension name, runs the molecular clock analysis on all of them
      
 OUTPUT:
    The script prints results to file test_with_fastas_q099_molClockTest.tab
       #nexfile lnL_unconstr    lnL_clock       LRT     X2_crit_val     df      p-val   mol_clock
       my_219956_thrA_cluO_lowH_341_1255.nex    -2160.56291     -2232.72785     144.3299        80      45      2.43372e-12     NO
       my_219967_dnaK_cluO_lowH_533_1394.nex    -1550.45096     -1582.35056     63.7992 80      45      0.0339092       YES     
 
 NOTES:
      1) qchisq(0.95, df);  H0 is rejected if X2 observed > critical value. df = no_tax -2
         The p-value is computed liake so: 1 - pchisq(X2,df)

 TODO:
     *  The logic of run_paup_clock_test() is not good, as it runs the for loop within the function;
          Recode that function as in run_paup_clock_test_with_user_tree() and add the rootmethod option -r to it
     * Explore the MPI parallelization possibilities of jmodeltes2 (read documentation)
     * implement clock analysis using dnaml|proml from the PAML package, particularly for proteins

HELP

exit

}

#>>>>>>>>>>>>>>>>>>>> END FUNCTIONS <<<<<<<<<<<<<<<<<<<<<


#-------------------------------------------------------#
#-------------------- GET OPTIONS ----------------------#
#-------------------------------------------------------#

# initialize globals
nex_file=
basename_output_table=global_molecular_ClockTest.tab
modeltest_outfile_suffix=
tree_extension=
rooting_method=midpoint
outgroup_OTU_no=
base_mod=GTR
runmode=
no_seq=
df=
q=0.99

# flags
VERBOSE=0
DEBUG=0
CHECK_DEPENDENCIES=0
RUN_PEXEC_CMDS=0
CONVERT_FASTA2NEXUS=0
CHECK_NO_SEQS_IN_NEXUS=1

# See bash cookbook 13.1 and 13.2
while getopts ':b:f:g:k:m:M:o:q:r:R:s:t:S:hHIVDK?:' OPTIONS
do
   case $OPTIONS in
   h)   print_help
        ;;
   H)   print_documentation
        ;;
   b)   basename_output_table=$OPTARG
        ;;
   f)   nex_file=$OPTARG
        ;;
   g)   no_cat=$OPTARG
        ;;
   k)	CHECK_NO_SEQS_IN_NEXUS=$OPTARG
        ;;
   m)   modeltest_outfile_suffix=$OPTARG
        ;;
   M)   base_mod=$OPTARG
        ;;
   o)   outgroup_OTU_no=$OPTARG
        ;;
   q)   q=$OPTARG
        ;;
   r)   rooting_method=$OPTARG
        ;;
   s)   no_subst_schemes=$OPTARG
        ;;
   t)   tree_extension=$OPTARG
        ;;
   S)   search_alg=$OPTARG
        ;;
   I)   pInv=1
        ;;
   K)   CHECK_DEPENDENCIES=1
        ;;
   n)	no_top_markers2keep=$OPTARG
        ;;
   R)   runmode=$OPTARG
        ;;
   D)   DEBUG=1
        ;;
   V)   VERBOSE=1
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

[ $CHECK_DEPENDENCIES -eq 1 ] && check_dependencies $CHECK_DEPENDENCIES

if [ -z $nex_file ]
then
      echo "# ERROR: need to provide the name of an alignment in nexus format"
      print_help
      exit 5
fi

[ -z $nex_file ] && CHECK_NO_SEQS_IN_NEXUS=1


if [ -z $no_cores ]
then
    no_cores=1
fi

if [ "$rooting_method" == "outgroup" ] && [ -z $outgroup_OTU_no ]
then
     echo "# ERROR: rooting method $rooting_method requires an outgroup sequence number to be provied!"
     print_help
     exit 7
fi

if [ -z $runmode ]
then
       echo "# ERROR: no runmode defined!"
       print_help
       exit 8   
fi

if [ $runmode -eq 1 ] && [ -z $tree_extension ]
then
       echo "# ERROR: runmode $runmode requires user tree extension file names!"
       print_help
       exit 8   
fi

###>>> Set the script's environment
#env_vars=$(set_script_environment) # returns: $distrodir $bindir $OS $no_proc
#[ $DEGUG ] && echo "env_vars:$env_vars"
#distrodir=$(echo $env_vars | awk '{print $1}')
#bindir=$(echo $env_vars | awk '{print $2}')
#OS=$(echo $env_vars | awk '{print $3}')
#no_proc=$(echo $env_vars | awk '{print $4}')
#
#[ $DEBUG -eq 1 ] && echo "distrodir:$distrodir|bindir:$bindir|OS:$OS|no_proc:$no_proc"

# 0.1 Determine if pipeline scripts are in $PATH; 
# if not, add them
#check_scripts_in_path $distrodir

# 0.2  Determine if second-party binaries are in $PATH; 
#  if they are not in $PATH then:
# i) will generate a symlink from $HOME/bin (if in path)
# to the binaries provided in $bindir.
# ii) If no $HOME/bin exists, or it is not in $PATH,
# then $bindir will be added to $PATH
set_bindirs $bindir

# 0.3 append the $distrodir/lib/R to R_LIBS and export
export R_LIBS="$R_LIBS:$distrodir/lib/R"


#-----------------------#
# >>>>> MAIN CODE <<<<< #
#-----------------------#

wkdir=$(pwd)

date_F=$(date +%F |sed 's/-/_/g')
date_T=$(date +%T |sed 's/:/./g')
start_time="$date_F$date_T"

echo
echo '----------------------------------------------------------------------------'
echo '>>>>>>>>>>>>>>>>>>>>>>>>> PROGRAM RUN PARAMETERS <<<<<<<<<<<<<<<<<<<<<<<<<<<'
echo '----------------------------------------------------------------------------'
echo "### $progname vers. $VERSION run on $start_time with the following parameters:"
echo "# work_directory=$wkdir"
echo "# nex_file=$nex_file | tree_extension=$tree_extension | basename_output_table=$basename_output_table"
echo "# runmode=$runmode | CONVERT_FASTA2NEXUS=$CONVERT_FASTA2NEXUS | root_method=$root_method | outgroup_OTU_no=$outgroup_OTU_no | q=$q"
echo "# jmodeltest2 params: no_subst_schemes=$no_subst_schemes | base_mod=$base_mod | no_cat=$no_cat | pInv=$pInv | search_alg=$search_alg"
echo "# DEBUG=$DEBUG | VERBOSE=$VERBOSE | COMPRESS_SRC_FASTAS=$COMPRESS_SRC_FASTAS" 
echo '----------------------------------------------------------------------------'
echo

#-------------------#
# >>> runmode 1 <<< #
#-------------------#
# NOTE: this is a slowet process in a genomic pipeline; re-write run_molecClock_test_jmodeltest2_paup.sh to parallelize with run_pexec_cmmds.sh
#cmd="run_pexec_cmmds.sh nex 'run_molecClock_test_jmodeltest2_paup.sh -R 1 -M $base_mod -t ph -e fasta -b molec_clock -q $q'"
#echo $cmd | bash

   # 1. convert newick trees to nexus-trees
   nex_basename=${nex_file%.*}
   [ $DEBUG -eq 1 ] && echo "# nex_basename:$nex_basename"
   tree=$(ls ${nex_basename}* | egrep "${tree_extension}$")
   [ $DEBUG -eq 1 ] && echo "# Converting tree $tree"
   newick2nexus $tree
   nextree=$(ls ${nex_basename}* | egrep "tre$")
   [ $DEBUG -eq 1 ] && echo "# nextree:$nextree"
   
   #results_table=$(ls mol_clock_M*_ClockTest.tab)
   [ -s mol_clock_M*_ClockTest.tab ] && results_table=$(ls mol_clock_M*_ClockTest.tab)
   [ -z $results_table ] && results_table=

   [ $DEBUG -eq 1 ] && echo "# Will test molecular clock with base model: $base_mod, root method: $rooting_method and outgroup: $outgroup_OTU_no"
   [ $DEBUG -eq 1 ] && echo "# running: run_paup_clock_test_with_user_tree() $nex_file $nextree $base_mod $results_table $rooting_method $outgroup_OTU_no"
   run_paup_clock_test_with_user_tree $nex_file $nextree $base_mod $results_table $rooting_method $outgroup_OTU_no
 
   # 5) cleanup and compress intermediary files
   [ -s ${nex_basename}_critical_X2_val.txt ] &&  rm ${nex_basename}_critical_X2_val.txt
   [ -s ${nex_basename}_LRT.txt ] && rm ${nex_basename}_LRT.txt
   [ -s ${nex_basename}_p_ChiSq_test_val.txt ] && rm ${nex_basename}_p_ChiSq_test_val.txt  	 



# 6) report exit status
date_F=$(date +%F |sed 's/-/_/g')
date_T=$(date +%T |sed 's/:/./g')
end_time="$date_F$date_T"


if [ -s $results_table ]
then
     echo
     echo '**********************************************************************************'
     echo " >>> $progname v.$VERSION exit status: SUCCESS!"
     echo "	  * Run finished at $end_time"         
     echo "	  * The results table is in $wkdir/$results_table"
     echo '----------------------------------------------------------------------------------' 
     echo
else
     echo
     echo '**********************************************************************************'
     echo " >>> $progname v.$VERSION exit status: ERROR!"
     echo "	  * The expected results table $results_table is not found in:"
     echo "	    $wkdir"
     echo '----------------------------------------------------------------------------------' 
     echo
fi


