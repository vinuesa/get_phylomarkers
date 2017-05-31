#!/usr/bin/env perl

# Author: Pablo Vinuesa, CCG-UNAM, Cuernavaca, Mexico. http://www.ccg.unam.mx/~vinuesa/
# popGen_summStats.pl, project started on April 16th, 2010.

### AIM:
# Used to compute basic population genetics descriptive statistics (pi, theta, tajima's_D) from aligned DNA sequences; 
#   codon alignments can be computed from unaligned DNA sequences if the requiered dependencies are fulfilled
#   PAUP* is used in -R 2|3 to compute the no. of parsimony informative sites, CI and HI

### DEPENDENCIES
# * Binaries: needs muscle in path for -R [1|3], and PAUP* for -R [2|3]
# * Requires BioPerl modules: Bio::AlignIO, Bio::PopGen::Utilities,Bio::PopGen::Statistics; 
# * Auxilliary scripts: assumes that run_muscle.pl and munge_sequence_files_bp.pl are found in $PATH


### TO DO:
# Compare results with those from DNAsp!!! ===> Oh no, the results differ notably!
# 0) VERY IMPORTANT: compute no. of haplotypes to derive Hd (haplotype diversity); Also compute the pi per gene and S per gene (not only per site)
# 0.1) tadapt code from http://www.pasteur.fr/recherche/unites/sis/formation/bioperl/ch03s02.html to firlter out X,N,R,Y,K,W,S,M,B,D,H,V, from MSAs!
# 0.2) include a sub that reads the significance values of Tajimas_D and Fu_Lis_D* and computes the significance for the results, printing * or **
#      or pass it as args!
# 1) include the subs from run_muscle.pl and munge_sequence_files_bp.pl in this script to reduce dependencies! NOTE: clustalw profile alns are more robust!
# 2) write a run_popGen_summStats.pl script to run it on the cluster or locally
# 3) write subs for additional analyses, like Fst, reading different populations and multiple loci ...
# 4) integrate the code with the Fasta2primers pipeline and the corresponding parsing scripts
#    to get the info on the same final output table
# 5) Finish/Re-think -R 7 
# 6) Write POD

# IMPORTANT NOTES: 
# 1) The bp methods seem to work using the segregating sites! (as deduced by comparing with the variscan output using total mutations or segr. sites!)
# 2) It is very important to scan alginments to eliminate columns containing unresolved or gapped sites!!! 
#     (see variscan manual pg. 10 CompleteDeletion = 0) as these affect the computation of Tajima's D and Fu and Li's statistics 
#     (see the discrepancies for the aroC aln, which contains an R and a Y, produced by the bp-based and variscan (with CompleteDeletion = 1) results!!!

#pablo@Tenerife:~/Projects/marfil/Escherichia_Abril10/variscan_test$ cat polymorphism_descript_stats.tab
#Alignment_name	no_seqs	aln_len	avg_perc_identity	pars_info_sites	consistency_idx	homoplasy_idx	segregating_sites	singletons	pi	theta	tajimas_D	fu_and_li_D_star
#1224_purB_cdnaln.fna	35	1371	98.49	61	0.6859	0.3141	103.00	45	0.01510	0.01824	-0.64	-1.48
#2381_aroC_cdnaln.fna	35	1089	96.06	126	0.4889	0.5111	161.00	47	0.03932	0.03590	0.36	-0.34
#
#pablo@Tenerife:~/Projects/marfil/Escherichia_Abril10/variscan_test$ cat variscan_sorted_summary.tab |grep -v 'phy.vs'
#        S        Eta      Eta_E         Pi      Theta   Tajima_D FuLi_Dstar FuLi_Fstar
#my_1224_purB_cdnaln.vs	103	107	45	0.0151040	0.0182428	-0.6425409	-1.4750553	-1.4098177	
#my_2381_aroC_cdnaln.vs	159	175	44	0.0393004	0.0356172	0.3897912	-0.2151974	-0.0056957

#------------------------------------------------------------------------------------------------#

# Code development history: April 16th 2010 -> 
my $VERSION ="v1.3.2"; # v1.3.1 some code cleanup; remove the symbols=.*;/;/ line form nexu data block
                    # v1.3 May 14th, 2017. Added portable shebang line for get_phylomarkers
                    # v1.2 Aug 9st, 2011: included an important test 'unless ($paup_data[2]){}' # i.e. next if the alignment has no variable sites! 
                    #    the alignment won't be processed by get_pop_summary_stats() if it has no variable sites; otherwise it will die (see comments in the unless block)
                    # v1.1 Jul 21st, 2011: simply changed method no_sequences by num_sequences, since the former is deprecated in BioPerl > 1.5
                    # v1.0 Nov. 10th, 2010; just corrected the order of table headers pi_gene pi_site ... printed to the *tab output file
                    # v0.99 Nov 7th, 2010; incudes -tT and -sS to provide the script with the lower and upper CI for Tajima'D and Fu & Li's D*
                    #                      Prints also the per-gene pi and theta
                    # v0.92; Oct 7th, 2010; fixed a subtle problem in  sub fas2nex(); added a system "sed 's/symbols.*/;/' statement
                    # to remove symbols="AcTagCGt"; from the nexus file, originated due to upper and lowercase letters of CDS and IGS regions
		    # from IGS amplicons !!! This resulted in PAUP* failing to run the analyses, leaving the following fields in polymorphism_descript_stats.tab empty:
		    # consistency_idx[5], homoplasy_idx[6], segregating_sites[7], singletons[8],

                    # v0.9 Jul 31st, uses an extended read_FASTA_sequence() sub, that removes cols with gaps and abiguous sites
                    # which cause problems in the estimation of Tajimas_D and Fu_and_Li's statistics.

                    # v0.8 Jul 27th, 2010, corrected a bug in computation of Pi/site and Theta/site; 
                    #  should use $stats->pi($pop) / @{$stats_aref}[1] (aln_length) and not @{$stats_aref}[0], no of seqs
		    #  See also the IMPORTANT NOTES section above. 
		    
                    ## v0.7 May 13th 2010. 
                    ## v0.6 May 11th 2010. Further debugging and simplification/documentation of getpots,  
                     #          and generalization of code for -R 2 for getting the alignment extension  -f cdnaln.fna or fna 
                    ## v0.5, May 1st 2010 ==> major debugging and code cleanup; apparently no bugs left; tested in all runmodes

use strict;
use warnings;
use File::Basename;
use Getopt::Std;

use FindBin '$Bin'; # added by Bruno May2017
use lib "$Bin/lib/perl/bioperl-1.5.2_102/";

use Bio::AlignIO;
use Bio::PopGen::IO;
use Bio::PopGen::Utilities;
use Bio::PopGen::Statistics;

my $progname = basename($0); # popGen_summStats.pl

my $ALN_THRESHOLD_PERCENTAGE = 95; # To be used in -R 7
#my $ENV{'TRAILINGZEROS'} = '';
my $alignment_ext ='';
my $unaln_dna_file_ext = '';
my $alignment_format = '';
my $run_mode = '';
my $threshold_percent = '';
my $nexus_aln_ext = '';
my $fasta2nexus = '';
my $hs = '';
my $nrep  = '';
my %opts =();
my ($tajimas_lowerCI, $tajimas_upperCI, $fu_lis_Dstar_lowerCI, $fu_lis_Dstar_upperCI);

getopts('e:f:n:r:F:R:t:T:s:S:hHND', \%opts);  # opt p could be a=alphabet [dna|prot]

if(($opts{'h'})||(scalar(keys(%opts))==0)) 
{ 
	print   "\nusage: $progname version $VERSION [options]\n";
	print   "-h \t this message\n";
	print   "-e \t dna file (unaligned) extension                [fna|fas|fasta]\n";
	print   "-f \t fasta alignment extension name                [fna|fas|fasta] (default fna for -R 3)\n";
	print   "-n \t nexus alignment extension name                [nex|nxs] (only for PAUP* searches)\n";
	print   "-F \t alignment file format                         [fasta|nexus|phylip|clustalw ... for Bio::PopGen::*]\n";
	print   "-H \t run heuristic parsimony search in PAUP*       (optional, requires nexus-formatted alns; default: NJ; describe;)\n";
	print   "   \t  (hs nrep=\$nrep start=step add=rand)\n";
	print   "-N \t generate nexus files from fasta alignments    (optional, requires -f)\n";
	print   "-r \t nrep for heuristic search in PAUP*            (optional, default 10)\n";
	print   "-D \t print code Development history and notes      (optional, default NO)\n";
	print   "-t \t provide the lower CI for tajima's D (-R 4)    (optional, default '')\n";
	print   "-T \t provide the upper CI for tajima's D (-R 4)    (optional, default '')\n";
	print   "-s \t provide the lower CI for Fu-Li's D* (-R 5)    (optional, default '')\n";
	print   "-S \t provide the upper CI for Fu-Li's D* (-R 5)    (optional, default '')\n";

	#print   "-T \t Threshold percentage for consensus alignments   (optional, default 95)\n";

	print STDOUT <<EOF;
-R  <int>   invokes pre-defined runmodes, as follows:
   1 -> generate codon alignments, starting from unaligned dna files [requires -e]
   2 -> run the basic popGen stats: Pi and Tajima's D on the cdnAlns [requires -f -F -n | -f -F -N] 
   3 -> run both, steps 1 and 2 [requires -e]
   4 -> print TajimasD CI table for given sample sizes
   5 -> print Fu and Li's D* CI table for given sample sizes
   6 -> print notes on processing run_Fasta2primers.pl output files for this script
   7 -> alignment stats to tab file (NOT IMPLEMENTED YET)

Run examples:
     popGen_summStats.pl -R 1 -e fas [assumes fasta-formatted input file for muscle]
     popGen_summStats.pl -R 2 -N -n nex -F fasta -f fna -H -r 100 || -R 2 -n nex -f fna -F fasta -H -r 100 (assumes nexus-formatted alns are present!!!) 
     # to evaluate a run_Fasta2primers.pl output need only the *dna_amp files; use -N to make new clean.nex files:
          nohup popGen_summStats.pl -f dna_amp -F fasta -N -R 2 -H -r 100 &> popGen_summStats_hs100.log & 
     popGen_summStats.pl -R 3 -e fna
     popGen_summStats.pl -R 4 
     popGen_summStats.pl -R 5 
     popGen_summStats.pl -R 6 

TODO: (Oct. 2012)
      1. include also the 99% confidence interval values for Tajima's D (only Fu-Li's D* with both 95% & 99% CIs).
      2. add Simpson's discriminatory power from compute_discriminatory_power.pl
      3. Need to read the significance talbles automatically, for example, from a __DATA__ section!
      4. Automatically remove ? from the nexus data symbols='A?GCT' line generated by convert_aln_format_batch_bp.pl
      5. Write a CGI version of the script to run on a web server
      6. Revise new Bio::Perl popGen HOWTO, to include further tests, including coalescent
      7. Revise the other bioperl code included herein, for further updates/optimizations
      

EOF

exit;
}

if(defined($opts{'D'})){ print_code_devel_history(); }
if(defined($opts{'R'})){ $run_mode = $opts{'R'}; }
else{ die "# $progname : need a run mode -R <1|2|3|4|5|6>; see help.\n";}

if(defined($opts{'N'}) and defined($opts{'f'}) ){ $fasta2nexus = 1; $alignment_ext = $opts{'f'} }
elsif(defined($opts{'N'}) and !defined($opts{'f'}) ){ die "# $progname : need to define -f to make fasta2nexus convesion; see help.\n";}
elsif($run_mode == 2 and !defined($opts{'N'}) and !defined($opts{'n'})){ die "# $progname : need to define -n nexus_aln_ext; see help.\n";}

if(defined($opts{'e'}) and $run_mode == 1 ){ $unaln_dna_file_ext = $opts{'e'}}
elsif(!defined($opts{'e'}) and $run_mode == 1 ){ die "# $progname : 0RUNMODE $run_mode needs [-e]\n"; }

if($run_mode == 2 and defined($opts{'F'}) and defined($opts{'n'}) and defined($opts{'f'}) and !defined($opts{'N'}))
{ 
    $alignment_format = $opts{'F'}; $nexus_aln_ext = $opts{'n'}; $alignment_ext = $opts{'f'};
}
elsif($run_mode == 2 and defined($opts{'f'}) and defined($opts{'F'}) and defined($opts{'N'}) )
{
   $alignment_ext = $opts{'f'}; $alignment_format = $opts{'F'}; $nexus_aln_ext = $opts{'n'};
}
elsif($run_mode == 2 and !defined($opts{'f'}) and defined($opts{'F'}) and defined($opts{'N'}) )
{ 
    die "# $progname: 1RUNMODE $run_mode needs -f and -F and -N; see help \n";
}
elsif( $run_mode == 2 and (!defined($opts{'F'}) or !defined($opts{'f'})) )
{ 
   die "# $progname: 1RUNMODE $run_mode needs -F and -f; see help \n";
}

if(defined($opts{'e'}) and $run_mode == 3){ $unaln_dna_file_ext = $opts{'e'}; $alignment_ext = 'cdnaln_clean.fna';  $alignment_format = 'fasta'; }
elsif($run_mode == 3 and !defined($opts{'e'})){ die "# $progname : 2RUNMODE $run_mode needs [-e]; see help\n"; }

if(defined($opts{'T'})){ $threshold_percent = $opts{'T'}; }
else{ $threshold_percent = $ALN_THRESHOLD_PERCENTAGE; }

if(defined($opts{'H'})){ $hs = 1; }
else { $hs = ''; }
if(defined($opts{'r'})){ $nrep = $opts{'r'} }
else { $nrep = 10; }

if(defined($opts{'t'})){ $tajimas_lowerCI = $opts{'t'}; }
else { $tajimas_lowerCI = ''; }
if(defined($opts{'T'})){ $tajimas_upperCI = $opts{'T'}; }
else { $tajimas_upperCI = ''; }
if(defined($opts{'t'})){ $fu_lis_Dstar_lowerCI = $opts{'s'}; }
else { $fu_lis_Dstar_lowerCI = ''; }
if(defined($opts{'T'})){ $fu_lis_Dstar_upperCI = $opts{'S'}; }
else { $fu_lis_Dstar_upperCI = ''; }

print "# $progname vers. $VERSION running in RUNMODE: $run_mode with the following parameters:\n"
     ."# -e $unaln_dna_file_ext -F $alignment_format -H $hs -f $alignment_ext -N $fasta2nexus -r $nrep\n"
     ."# -t $tajimas_lowerCI -T $tajimas_upperCI -s $fu_lis_Dstar_lowerCI -S $fu_lis_Dstar_upperCI\n\n";
   

if($run_mode == 1)
{
    make_cdn_aln($unaln_dna_file_ext);
    
    while( my $cdn_aln = <*cdnaln.fna>)
    {
          my %cleanFasta = read_FASTA_sequence($cdn_aln,0,0,0,1);
    }
}
elsif($run_mode == 2)
{
    ######### THIS SHOULD BE BETTER WRITTEN AS IN -R 3, using fas2nex() and providing a cdnAln extension name to pass to fas2nex(); 
    ### don't call system to run convert_aln_format_batch_bp.pl
    # convert file to nexus format for paup* to run
    
    while( my $cdn_aln = <*$alignment_ext*> )
    {
    	 #print "# working on cdn_aln: $cdn_aln\n"; exit;
	 my %cleanFasta = read_FASTA_sequence($cdn_aln,0,0,0,1);
    }	  
    
    if ($fasta2nexus)
    {  
	print "# Running fas2nex() with: $alignment_ext\n";
	my $aln2process = "clean.".$alignment_ext;
        fas2nex($aln2process); #$alignment_ext  'clean.fna'
	$nexus_aln_ext = 'clean.nex';
    }
    
    # remove the "symbols=.*;" ending from the format interleave datatype=dna   gap=- symbols="GCTA"; nexus line
    #my $nexf="";
    #system("for nexf in *nex; do sed 's/ symbols=.*;/;/' $nexf > ${nexf}ed; mv ${nexf}.ed $nexf; done");
    
    open PAUP, "> paup.cmd", or die "can't write file paup.cmd: $!\n";
    if( $hs ){ print PAUP "hs nrep=$nrep start=step add=rand; describe; "; }
    else { print PAUP "nj; describe; "; }
    close (PAUP);
       
    open OUT, ">polymorphism_descript_stats.tab" or die "can't write to file polymorphism_descript_stats.tab: $!\n";
    print OUT "#Alignment_name\tno_seqs\taln_len\tavg_perc_identity\tpars_info_sites\tconsistency_idx\thomoplasy_idx\tsegregating_sites\tsingletons\tpi_per_gene\tpi_per_site\ttheta_per_gene\ttheta_per_site\ttajimas_D\tfu_and_li_D_star\n";
    close OUT;
    print     "#Alignment_name\tno_seqs\taln_len\tavg_perc_identity\tpars_info_sites\tconsistency_idx\thomoplasy_idx\tsegregating_sites\tsingletons\tpi_per_gene\tpi_per_site\ttheta_per_gene\ttheta_per_site\ttajimas_D\tfu_and_li_D_star\n";
    
    foreach my $nexus_file (<*$nexus_aln_ext>)   
    { 
         #print "# working on nexus file $nexus_file\n";
	 my $alignment_file = '';
	 if($fasta2nexus and $alignment_ext =~ /\w+\.(\w+)$/)  # for example to eliminate _cdnaln from _cdnaln.fna 
	 {
	     $alignment_ext = $1;
	     $alignment_file = (split /\./, $nexus_file)[0] . ".$alignment_ext"; #print "# working on nexus_file $nexus_file alignment_file $alignment_file\n"; 
	 }
	 elsif($alignment_ext =~ /\w+\.(\w+)$/)  # for example to eliminate _cdnaln from _cdnaln.fna 
	 {
	     $alignment_ext = $1;
	     $alignment_file = (split /\./, $nexus_file)[0] . ".$alignment_ext"; #print "# working on nexus_file $nexus_file alignment_file $alignment_file\n"; 
	 }
	 else
	 {
	     $alignment_file = (split /\./, $nexus_file)[0] . ".$alignment_ext"; #print "# working on nexus_file $nexus_file alignment_file $alignment_file\n";
	 }
	 my @paup_data = paup_parsimony($nexus_file, $hs);  # print "# paup_data are: [@paup_data]: $nexus_file, $hs\n"; 
	 
	 
	 # Check the alignments are not invariant; if they are, we would get the following ERROR MSGs from Bio::Perl modules
	 #   Use of uninitialized value $nms[0] in pattern match (m//) at /usr/local/share/perl/5.10.1/Bio/PopGen/Population.pm line 400.
         #   Illegal division by zero at /usr/local/share/perl/5.10.1/Bio/PopGen/Statistics.pm line 388.
	 unless ($paup_data[2]) # i.e. next if the alignment has no variable sites! 
	 {
           # THIS IS AN IMPORTANT TEST, as of version 
	   # paup_data are: [0]: M_avium_hsp65_65_aislados_BER_classified_clean.nex, 1  <== no segregating sites => no data for consitency idx nor for homoplasy idx
           # paup_data are: [0 1.0000 0.0000]: M_avium_recA_65_aislados_BER_classified_clean.nex, 1
           # paup_data are: [0 1.0000 0.0000]: M_avium_rpoB1_65_aislados_BER_classified_clean.nex, 1
           # paup_data are: [1 1.0000 0.0000]: M_avium_rpoB2_65_aislados_BER_classified_clean.nex, 1

           open INVARIANT ,'>>', 'invariant_alignments.txt' or warn "can't write to file monomorphic_alignments.txt $!\n\n";
	   print INVARIANT "# found invariant alignment: $nexus_file\n";
	   close INVARIANT;
	   next; # <== next if the alignment has no variable sites!
	 } 
	    my $pop = get_pop_object_from_aln($alignment_file, $alignment_format);
            get_pop_summary_stats($pop, $alignment_file,$alignment_format,$tajimas_lowerCI,$tajimas_upperCI,$fu_lis_Dstar_lowerCI,$fu_lis_Dstar_upperCI,@paup_data);

    }
}
elsif($run_mode == 3)
{
    $alignment_format = 'fasta';
    make_cdn_aln($unaln_dna_file_ext);
    while( my $cdn_aln = <*cdnaln.fna> )
    {
         read_FASTA_sequences($cdn_aln,0,0,0,1); # $infile, $remove_gap_cols, $skipbadCDSs, $skipidentical, $skip_ambiguous_and_gap_cols
    }
    
    fas2nex('cdnaln_clean.fna');

    open PAUP, "> paup.cmd", or die "can't write file paup.cmd: $!\n";
    if( $hs ){ print PAUP "hs nrep=$nrep start=step add=rand; describe; "; }
    else { print PAUP "nj; describe; "; }
    close (PAUP);

    open OUT, ">polymorphism_descript_stats.tab" or die "can't write to file polymorphism_descript_stats.tab: $!\n";
    print OUT "#Alignment_name\tno_seqs\taln_len\tavg_perc_identity\tpars_info_sites\tconsistency_idx\thomoplasy_idx\tsegregating_sites\tsingletons\tpi_per_gene\tpi_per_site\ttheta_per_gene\ttheta_per_site\ttajimas_D\tfu_and_li_D_star\n";
    
    print "#Alignment_name\tno_seqs\taln_len\tavg_perc_identity\tpars_info_sites\tconsistency_idx\thomoplasy_idx\tsegregating_sites\tsingletons\tpi_per_gene\tpi_per_site\ttheta_per_gene\ttheta_per_site\ttajimas_D\tfu_and_li_D_star\n";
    foreach my $nexus_file (<*clean.nex>)
    { 
         my $alignment_file = (split /\./, $nexus_file)[0] . ".$unaln_dna_file_ext"; #print "# working on nexus_file $nexus_file alignment_file $alignment_file\n"; 
	 my @paup_data = paup_parsimony($nexus_file, $hs);                             #print "# paup_data are: @paup_data\n"; 
	 my $pop = get_pop_object_from_aln($alignment_file, $alignment_format);
         get_pop_summary_stats($pop, $alignment_file, $alignment_format, $tajimas_lowerCI, $tajimas_upperCI, $fu_lis_Dstar_lowerCI, $fu_lis_Dstar_upperCI, @paup_data, );
    }
}
elsif($run_mode == 4)
{
    print_TajimasD_CI_table();
}
elsif($run_mode == 5)
{
    print_Fu_Li_Dstar_CI_table();
}
elsif($run_mode == 6)
{
	print STDOUT <<EOF;
       
Notes: 
 
 1.1) # After running run_Fasta2primers.pl, get the best amplicons from the parse_Fasta2primers.pl output file by running:
   mkdir selected; for file in \$(grep '>' parsed_F2P_* | perl -ne '\$alrt = (split /aLRT=/)[1]; print if \$alrt >= 0.85' \
     |cut -d' ' -f13 |cut -d= -f2); do cp \$file selected/; done
 1.2) # untar and extract the *dna_amp files from the selected tgzs; then remove all source tgz files
    cd selected/; for file in \$(ls *tgz); do tar -xzf \$file --wildcards --no-anchored '*dna_amp'; rm -f \$file; done

 2)   # Run popGen_summStats.pl
 2.1) convert_aln_format_batch_bp.pl fasta dna_amp nexus nex
 2.2) nohup popGen_summStats.pl -e dna_amp -F fasta -n nex -R 2 -H -r 100 &> popGen_summStats_hs100.log & # to evaluate a run_Fasta2primers.pl output
 
 3) # Parse the polymorphism_descript_stats.tab output file generated by \$progname using the following 1liner 
    (and the confidence values indicated by Tajima 1989, which are displayed by \$progname -R 4):
     
      grep 'my_' polymorphism_descript_stats.tab |sort -k9nr | \
      perl -ne '\$TajimasD=(split)[11]; print if(\$TajimasD >= -1.086 && \$TajimasD <= 2.027);' |sort -k1
      > sorted_233_neutral_markers4Ecoli_TajimasD.out  

4) # grep out the neutral markers from parsed_F2P_Tm35P80.out_stats.tab
cut -f1 sorted_233_neutral_markers4Ecoli_TajimasD.out | sed 's/my_//; s/\.dna_amp//' > neutral_markers2grep.txt 
grep -wf neutral_markers2grep.txt parsed_F2P_Tm35P80.out_stats.tab | sort -k1 > sorted_selected_parsed_F2P_Tm35P80.out_stats.tab
paste sorted_F2P_alns_sorted_by_pi.out sorted_selected_parsed_F2P_Tm35P80.out_stats.tab |less
paste sorted_233_neutral_markers4Ecoli_TajimasD.out sorted_selected_parsed_F2P_Tm35P80.out_stats.tab | sort -k9nr > SORTED_PASTED_FILES.out
wc -l SORTED_PASTED_FILES.out 

EOF
exit;
}
## Runmode 7 has to be finished / re-designed / thought-over
elsif($run_mode == 7)
{
    open OUT, ">>alignment_stats.tab" or die "can't write to file alignment_stats.tab: $!\n";
    print OUT "#Alignment_name\tpi\tsegregating_sites\tsingleton\ttheta\ttajimas_D\tfu_and_li_D_star\n";
    
    print "#Alignment_name\tpi\tsegregating_sites\tsingletons\ttheta\ttajimas_D\tfu_and_li_D_star\n";
    foreach my $alignment_file (<*$alignment_ext>)
    {
        print_alignment_stats($alignment_file,$alignment_format,$threshold_percent,);
    }	
}



## Make cleaup

print "\n";
for my $file (<*.*>)
{
    if (! -s $file)
    {
       print "# will remove empty file: $file\n";
       unlink $file;
    }
}

#------------------------------------------------------------------------------------------------#
##############################>>>>>>>>>> SUBROUTINES <<<<<<<<<<###################################
#------------------------------------------------------------------------------------------------#

sub print_code_devel_history
{
    print STDOUT <<EOF;

    Code development history for $progname $VERSION 
     v0.7 May 13th 2010. 
     v0.6 May 11th 2010. Further debugging and simplification/documentation of getpots,  
    	      and generalization of code for -R 2 for getting the alignment extension  -f cdnaln.fna or fna 
     v0.5, May 1st 2010 ==> major debugging and code cleanup; apparently no bugs left; tested in all runmodes


# Author: Pablo Vinuesa, CCG-UNAM, Cuernavaca, Mexico. http://www.ccg.unam.mx/~vinuesa/
# popGen_summStats.pl, project started on April 16th, 2010.

### AIM:
# Used to compute basic population genetics descriptive statistics (pi, theta, tajima's_D) from aligned DNA sequences; 
#   codon alignments can be computed from unaligned DNA sequences if the requiered dependencies are fulfilled
#   PAUP* is used in -R 2|3 to compute the no. of parsimony informative sites, CI and HI

### DEPENDENCIES
# * Binaries: needs muscle in path for -R [1|3], and PAUP* for -R [2|3]
# * Requires BioPerl modules: Bio::AlignIO, Bio::PopGen::Utilities,Bio::PopGen::Statistics; 
# * Auxilliary scripts: assumes that run_muscle.pl and munge_sequence_files_bp.pl are found in \$PATH


### TO DO:
# Compare results with those from DNAsp!!! ===> Oh no, the results differ notably!
# 0) include a sub that reads the significance values of Tajimas_D and Fu_Lis_D* and computes the significance for the results, printing * or **
# 1) include the subs from run_muscle.pl and munge_sequence_files_bp.pl in this script to reduce dependencies!
# 2) write a run_popGen_summStats.pl script to run it on the cluster or locally
# 3) write subs for additional analyses, like Fst, reading different populations and multiple loci ...
# 4) integrate the code with the Fasta2primers pipeline and the corresponding parsing scripts
#    to get the info on the same final output table
# 5) Finish/Re-think -R 7 
# 6) Write POD

EOF
exit;
}

#------------------------------------------------------------------------------------------------#

sub get_pop_object_from_aln
{
    my $stats = Bio::PopGen::Statistics->new();
    my ($alignment_file, $alignment_format) = @_; #print "# get_pop_object_from_aln: alignment_file $alignment_file, alignment_format  $alignment_format\n";
    my $in  = Bio::AlignIO->new(-format => $alignment_format, -file => $alignment_file);
    my $aln = $in->next_aln;
    my $pop = Bio::PopGen::Utilities->aln_to_population(-alignment=>$aln);
    return $pop;							
}
#------------------------------------------------------------------------------------------------#
sub get_pop_summary_stats
{       
     my ($pop, $alignment_file, $alignment_format, $tajimas_lowerCI, $tajimas_upperCI, $fu_lis_Dstar_lowerCI,$fu_lis_Dstar_upperCI, @paup_data ) = @_;
     #print "# get_pop_summary_stats() received:\n $pop, $alignment_file, $alignment_format, $tajimas_lowerCI, $tajimas_upperCI, $fu_lis_Dstar_lowerCI,$fu_lis_Dstar_upperCI, @paup_data\n";
     my $stats = Bio::PopGen::Statistics->new();
     my $stats_aref = get_alignment_stats($alignment_file,$alignment_format); # \@stats = $aln_length, $avg_percent_identity
     
     my $pi_gene = ($stats->pi($pop));
     my $pi_site = ($stats->pi($pop) / @{$stats_aref}[1]); # <- print $pi_site per site #print "#", $stats->pi($pop), " / ", @{$stats_aref}[1], "\n"; 
     my $tajima_D = $stats->tajima_D($pop);
     my $theta_gene = ($stats->theta($pop));
     my $theta_site = ($stats->theta($pop)/ @{$stats_aref}[1]); # <- print $theta_site per site
     my $singleton_count = $stats->singleton_count($pop);
     my $segregating_sites_count = $stats->segregating_sites_count($pop);
     my $fu_and_li_D_star = $stats->fu_and_li_D_star($pop);
     
     $tajima_D = sprintf("%.3f", $tajima_D);
     if ( ($tajimas_lowerCI && $tajimas_upperCI) && ( $tajima_D < $tajimas_lowerCI || $tajima_D > $tajimas_upperCI ))
     {
         $tajima_D = $tajima_D . '*';
     }    
     
     $fu_and_li_D_star = sprintf("%.3f", $fu_and_li_D_star);
     if (($fu_lis_Dstar_lowerCI && $fu_lis_Dstar_upperCI) && ($fu_and_li_D_star < $fu_lis_Dstar_lowerCI || $fu_and_li_D_star > $fu_lis_Dstar_upperCI)) 
     {
         $fu_and_li_D_star = $fu_and_li_D_star . '*';
     }

     open OUT, ">>polymorphism_descript_stats.tab" or die "can't write to file polymorphism_descript_stats.tab: $!\n";
     print OUT $alignment_file, "\t", @{$stats_aref}[0],"\t", @{$stats_aref}[1],"\t", @{$stats_aref}[2],"\t", $paup_data[0], "\t", $paup_data[1], "\t", $paup_data[2], "\t", sprintf("%.2f",$segregating_sites_count), "\t", $singleton_count, "\t", sprintf("%.2f", $pi_gene), "\t", sprintf("%.5f",$pi_site), "\t", sprintf("%.2f", $theta_gene), "\t", sprintf("%.5f",$theta_site), "\t", $tajima_D, "\t", $fu_and_li_D_star, "\n";
     print $alignment_file, "\t", @{$stats_aref}[0],"\t", @{$stats_aref}[1],"\t", @{$stats_aref}[2],"\t", $paup_data[0], "\t", $paup_data[1], "\t", $paup_data[2], "\t", $segregating_sites_count, "\t", $singleton_count, "\t", sprintf("%.2f", $pi_gene), "\t", sprintf("%.5f",$pi_site), "\t", sprintf("%.2f", $theta_gene), "\t", sprintf("%.5f",$theta_site), "\t", $tajima_D, "\t", $fu_and_li_D_star, "\n";
     close OUT;
}
#------------------------------------------------------------------------------------------------#

sub make_cdn_aln
{
    my ($unaln_dna_file_ext) = @_;
    my @ualn_dna_files = <*$unaln_dna_file_ext>;
    #print "#will work on the following ualn_dna_files: @ualn_dna_files\n"; exit;
    #print "# running munge_sequence_files_bp.pl -R 5 -e $unaln_dna_file_ext\n";

    system "munge_sequence_files_bp.pl -R 5 -e $unaln_dna_file_ext";
    
    print "# running run_muscle.pl faa ...\n";
    system ("run_muscle.pl faa > /dev/null 2>&1");
       
    my $counter = 0;
    foreach my $unaln_dna_file ( @ualn_dna_files)
    { 
        $counter++;
	my $basename = (split /\./, $unaln_dna_file)[0];
	my $prot_aln = $basename . '_translated_ref2.faa';
	print "# running AAaln2cndAln.pl [$unaln_dna_file] [$prot_aln]\n"; #exit;
        system "AAaln2cndAln.pl $unaln_dna_file $prot_aln";
    }
    #"for file in $(ls *.$unaln_dna_file_ext); do cdnAln=${file%.$unaln_dna_file_ext}_translated_ref2.faa; echo $cdnAln; AAaln2cndAln.pl $file $cdnAln; done"
    print "# make_cdn_aln() processed $counter alignments ...\n";
}

#------------------------------------------------------------------------------------------------#

sub read_FASTA_sequence
{
	# if $remove_gap_col assumes FASTA file is actually a multiple alignment
	# in FASTA format
	# if $skipbadCDSs ignores nt sequences with length%3 > 0

    	 my ( $infile, $remove_gap_cols, $skipbadCDSs, $skipidentical, $skip_ambiguous_and_gap_cols ) = @_;
	 
	 my (%FASTA,$name,$seq,$n_of_sequences,$length,$maxlength,$pos,$isgap,$seqid,$isambig);
	 
	 my ($fasta_base, $ext) = (split /\./, $infile)[0,1];
	 my $clean_fasta_file = $fasta_base . '_clean' . ".$ext";
	 my $infile_stats = $fasta_base . '_gap_and_amb_sites.stats';
	 $n_of_sequences = $maxlength = 0;
	 open(FASTA,$infile) || die "# read_FASTA_sequence: cannot read $infile $!:\n";
	 while(<FASTA>)
	 {
	 	  if(/^\>/)
		  {
		   	$name = $_; 
			$n_of_sequences++;
			#$seqid = sprintf("%0$ENV{'TRAILINGZEROS'}d",$n_of_sequences); #print "#$seqid#\n";
			$seqid = $n_of_sequences;
			#print "# seqid: $seqid\n";
			$FASTA{$seqid}{'NAME'} = $name;
		  }    	 	   	
		  else
		  {
		   	$FASTA{$seqid}{'SEQ'} .= $_;
			$FASTA{$seqid}{'SEQ'} =~ s/[\s|\n]//g;
				
			$length = length($FASTA{$seqid}{'SEQ'});
			if($length > $maxlength){ $maxlength = $length; }
		  }
	 }
	 close(FASTA);

	if($skipbadCDSs)
	{
		# check length is divisible by 3 and that no inframe STOP codons exist
		my %goodCDSs_FASTA;
		foreach $seq (keys(%FASTA))
                {
			my $ntseq = Bio::Seq->new( -display_id => 'tmp', -seq => $FASTA{$seq}{'SEQ'} );
                        my $protseq = $ntseq->translate()->seq(); chop $protseq;
			if($protseq =~ /\*/)
                        {
                                print "# read_FASTA_sequence : skipped CDS sequence (inframe STOP codon) $FASTA{$seq}{'NAME'}\n";
                                next;
                        }
		
                	if($FASTA{$seq}{'SEQ'} !~ /QWERYIPSDFHKLMNV/ && length($FASTA{$seq}{'SEQ'})%3)
			{
				my $choppednts = 0;
				while(length($FASTA{$seq}{'SEQ'})%3){ chop $FASTA{$seq}{'SEQ'}; $choppednts++; }
				print "# read_FASTA_sequence : chopped CDS sequence (3' $choppednts nts) $FASTA{$seq}{'NAME'}\n";
			}
			
			#if(length($protseq) < 100){ print "mira: $protseq\n"; }
			$goodCDSs_FASTA{$seq} = $FASTA{$seq};
		}
		%FASTA = %goodCDSs_FASTA;
	}

	if($skipidentical)
	{
		my (%nrFASTA,%identical,$seq2);
                foreach $seq (keys(%FASTA))
                {
			$FASTA{$seq}{'IDENTICALS'} = '';
			next if($identical{$seq});
			foreach $seq2 (keys(%FASTA))
			{
				next if($seq == $seq2);
				if($FASTA{$seq}{'SEQ'} eq $FASTA{$seq2}{'SEQ'} ||        # same length
					$FASTA{$seq}{'SEQ'} =~ /$FASTA{$seq2}{'SEQ'}/ || # different length 
					$FASTA{$seq2}{'SEQ'} =~ /$FASTA{$seq}{'SEQ'}/)
				{
					$identical{$seq2} = $seq;
					$FASTA{$seq}{'IDENTICALS'} .= "$seq2,"; #print "mira $seq $seq2\n";
					$FASTA{$seq}{'NUMBER_IDENTICALS'}++;
				}
			}
		}

		foreach $seq (keys(%FASTA))
                {
			if($identical{$seq})
			{
				print "# read_FASTA_sequence : skipped sequence identical to $identical{$seq}: $FASTA{$seq}{'NAME'}\n";
				next;
			}
			$nrFASTA{$seq} = $FASTA{$seq};

			# keep track of identical sequences
			if($FASTA{$seq}{'NUMBER_IDENTICALS'})
			{ 
				chomp($nrFASTA{$seq}{'NAME'});
				$nrFASTA{$seq}{'NAME'} .= " | identical sequences=$FASTA{$seq}{'NUMBER_IDENTICALS'}\n";
			}
		}
		%FASTA = %nrFASTA;
	}

	#if($remove_gap_cols && $remove_gap_cols == 1) 
	if($remove_gap_cols)
	{
		for($pos=0;$pos < $length; $pos++)
		{	
			$isgap=0;
			foreach $seq (keys(%FASTA))
			{
				if(substr($FASTA{$seq}{'SEQ'},$pos,1))
				{
					if(substr($FASTA{$seq}{'SEQ'},$pos,1) eq '-'){ $isgap=1; last; }
				}
				else{ $isgap=1; last }
			}
			
			if($isgap == 0)
			{
				foreach $seq (keys(%FASTA))
				{
					$FASTA{$seq}{'NOGAPSEQ'} .= substr($FASTA{$seq}{'SEQ'},$pos,1);
				}
			}
		}
		
		
		foreach $seq (keys(%FASTA))
		{
			$FASTA{$seq}{'SEQ'} = $FASTA{$seq}{'NOGAPSEQ'}
		}
	}

        ### Pablo Jul 29th, 2010
	my $skip_col_counter = 0;
	if($skip_ambiguous_and_gap_cols)
	{
		for($pos=0;$pos < $length; $pos++)
		{	
			$isambig=0;
			foreach $seq (keys(%FASTA))
			{
				if(substr($FASTA{$seq}{'SEQ'},$pos,1))
				{
					if(substr($FASTA{$seq}{'SEQ'},$pos,1) !~ /[ACGT]/i){ $isambig=1; $skip_col_counter++; last; }
				}
				else{ $isambig=1; last }
			}
			
			if($isambig == 0)
			{
				foreach $seq (keys(%FASTA))
				{
					$FASTA{$seq}{'NOGAPAMBSEQ'} .= substr($FASTA{$seq}{'SEQ'},$pos,1);
				}
			}
		}
		
		
		foreach $seq (keys(%FASTA))
		{
			$FASTA{$seq}{'SEQ'} = $FASTA{$seq}{'NOGAPAMBSEQ'}
		}
	}
        
	open OUT, ">$clean_fasta_file" or die "can't write to output file $clean_fasta_file: $!\n";
	foreach my $seq (sort keys %FASTA)
        {
            print OUT ">$seq\n", $FASTA{$seq}{SEQ}, "\n";
        } 
        close(OUT);
	return %FASTA;
}

#------------------------------------------------------------------------------------------------#

sub get_alignment_stats
{
        # See documentation in Bio::SimpleAlign
        # See perldoc Bio::Align::AlignI
        my($input_file,$inputformat,$threshold_percent) = @_;   #print "# sub get_alignment_stats input_file: $input_file\n";
        my @stats   = ();
        my $stream  = Bio::AlignIO->new(-file => $input_file,   -format => $inputformat);
        my    $aln  = $stream->next_aln();

        my $aln_length = $aln->length;  #print "# sub get_alignment_stats  alignment_length : $aln_length\n"; exit;
        my $avg_percent_identity = sprintf ("%.2f",$aln->average_percentage_identity);
        my $no_of_seqs = $aln->num_sequences;
        push @stats, $no_of_seqs, $aln_length, $avg_percent_identity;
        return \@stats;
}

#------------------------------------------------------------------------------------------------#
sub paup_parsimony
{
     my ($nexus_file, $hs) = @_ ;
     my @trash;
     my $paup_data_string = (split /\./, $nexus_file)[0] . '_paup_data.tmp';
     my @paup_data = ();
     if(!$hs)
     {
         system "cat paup.cmd | paup -n $nexus_file | egrep 'parsimony-informative characters|Consistency index|Homoplasy index' |cut -d= -f2 | perl -pe 'chomp;' | sed 's/ /\t/g'> $paup_data_string";
     }
     elsif ($hs)
     {
         system "cat paup.cmd | paup -n $nexus_file | egrep 'parsimony-informative characters|Consistency index|Homoplasy index' |cut -d= -f2 | perl -pe 'chomp;' | sed 's/ /\t/g' | cut -f 3-5 > $paup_data_string";
     } 
      
     # print "# paup_parsimony() -> paup string for $nexus_file is: $paup_data_string\n";
     open DATA, $paup_data_string or die "can't open file $paup_data_string: $!\n";
     while (<DATA>)
     {
        @paup_data = split;
     }
     push @trash, $paup_data_string;
     unlink @trash;
     return @paup_data;
}
#------------------------------------------------------------------------------------------------#

sub fas2nex
{
    my ($ext) = @_;
    my $counter = 0;
    foreach my $infile ( <*$ext> )
    {  
	my $basename = (split(/\./, $infile))[0]; 
		
	my $in  = Bio::AlignIO->new(-file => $infile,   -format => 'fasta');
	my $out = Bio::AlignIO->new(-file => ">$basename.nex", -format => 'nexus');

	while ( my $aln = $in->next_aln() ) 
	{
		    $out->write_aln($aln);
		    
		    # this is to remove symbols="AcTagCGt"; from the nexus file,
		    # originated due to upper and lowercase letters of CDS and IGS regions from IGS amplicons
		    system "sed 's/ symbols=.*;/;/' $basename.nex > $basename.nexed; mv $basename.nexed $basename.nex"
	}
	$counter++
    }
    print "# fas2nex(): $counter input files with were converted to nexus format\n\n";
}
#------------------------------------------------------------------------------------------------#

sub print_TajimasD_CI_table
{
  # Taken from Table2 of Tajima 1989. Genetics 123:585-595
  # @sample_size = qw(4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48 49 50 55 60 65 70 75 80 85 90 95 100 110 120 130 140 150 200 250 300 350 400 450 500 600 800 1000)
  # @upper_lower_95percCI = qw(-0.876 - 2.232 -1.269 - 1.834 -1.478 - 1.999 - 1.608 - 1.932 -1.663 - 1.975 -1.713 - 1.954 -1.733 - 1.975 -1.757 - 1.966 - 1.765 - 1.979 -1.779 - 1.976 -1.783 - 1.985 -1.791 - 1.984 -1.793 - 1.990 -1.798 - 1.990 -1.799 - 1.996 -1.802 - 1.996 -1.803 - 2.001 -1.805 - 2.001 -1.804 - 2.005 -1.806 - 2.006 -1.806 - 2.009 -1.807 - 2.010 -1.807 - 2.013 -1.807 - 2.014 -1.807 - 2.017 -1.807 - 2.018 -1.807 - 2.020 -1.807 - 2.021 -1.806 - 2.023 -1.806 - 2.024 -1.806 - 2.026 -1.806 - 2.027 -1.805 - 2.029 -1.805 - 2.030 -1.804 - 2.031 -1.804 - 2.032 -1.804 - 2.033 -1.803 - 2.034 -1.803 - 2.036 -1.803 - 2.037 -1.802 - 2.038 -1.802 -, 2.039 -1,801 - 2.040 -1.801 - 2.041 -1 .SO0 - 2.042 -1.800 - 2.042 -1.800 - 2.044 -1.797 - 2.048 -1.795 - 2.052 -1.793 - 2.055 -1.791 - 2.058 -1.790 - 2.061 - 1.788 5 2.064 -1.786 - 2.066 - 1.784 - 2.069 -1.783 - 2.071 -1.781 - 2.073 -1.779 - 2.077 -1.776 - 2.080 -1.774 - 2.084 -1.771 - 2.086 -1.769 - 2.089 -1.765 - 2.095 -1.760 - 2.100 -1.754 - 2.107 -1.748 - 2.114 -1.744-2.119 -1.740 - 2.123 -1.737 - 2.127 -1.734 - 2.130 -1.728 - 2.135 -1.721 - 2.143 -1.715 - 2.150 -1.760 - 2.100 -1.754 - 2.107 -1.748 - 2.114 -1.744-2.119 -1.740 - 2.123 -1.737 - 2.127 -1.734 - 2.130 -1.728 - 2.135 -1.721 - 2.143 -1.715 - 2.150)

print STDOUT <<EOF;
# Confidence intervals for TajimasD given a sample size, assuming a beta distribution
#n	lower	upper
4	-0.876	2.232 
5	-1.269	1.834 
6	-1.478	1.999
7	-1.608	1.932 
8	-1.663	1.975 
9	-1.713	1.954 
10	-1.733	1.975 
11	-1.757	1.966
12	-1.765	1.979 
13	-1.779	1.976 
14	-1.783	1.985 
15	-1.791	1.984 
16	-1.793	1.990 
17	-1.798	1.990 
18	-1.799	1.996 
19	-1.802	1.996 
20	-1.803	2.001 
21	-1.805	2.001 
22	-1.804	2.005 
23	-1.806	2.006 
24	-1.806	2.009 
25	-1.807	2.010 
26	-1.807	2.013 
27	-1.807	2.014 
28	-1.807	2.017 
29	-1.807	2.018 
30	-1.807	2.020 
31	-1.807	2.021 
32	-1.806	2.023 
33	-1.806	2.024 
34	-1.806	2.026 
35	-1.806	2.027 
36	-1.805	2.029 
37	-1.805	2.030 
38	-1.804	2.031 
39	-1.804	2.032 
40	-1.804	2.033 
41	-1.803	2.034 
42	-1.803	2.036 
43	-1.803	2.037 
44	-1.802	2.038 
45	-1.802	2.039 
46	-1.801	2.040 
47	-1.801	2.041 
48	-1.800	2.042 
49	-1.800	2.042 
50	-1.800	2.044 
55	-1.797	2.048 
60	-1.795	2.052 
65	-1.793	2.055 
70	-1.791	2.058 
75	-1.790	2.061
80	-1.788 	2.064
85	-1.786	2.066  
90	-1.784	2.069 
95	-1.783	2.071 
100	-1.781	2.073 
110	-1.779	2.077 
120	-1.776	2.080 
130	-1.774	2.084 
140	-1.771	2.086 
150	-1.769	2.089 
175	-1.765	2.095 
200	-1.760	2.100 
250	-1.754	2.107 
300	-1.748	2.114
350	-1.744	2.119 
400	-1.740	2.123 
450	-1.737	2.127 
500	-1.734	2.130 
600	-1.728	2.135 
800	-1.721	2.143 
1000	-1.715	2.150 

EOF
}
#------------------------------------------------------------------------------------------------#

sub print_Fu_Li_Dstar_CI_table
{
# Taken from Table 2 of Fu and Li 1993. Genetics 133:693-709
# Statistical tests of neutrality of mutations

print STDOUT <<EOF;
# Table 2 of Fu and Li 1993. Genetics 133:693-709
# Percentage points of statistic D* as functions of sample size
# If n=10 and D*obs = -2.11 then the result is significant at the  0.01 sign. level
#       <<<Left tail vals>>>  <<<right tail vals>>>
#n	0.01	0.025	0.05	0.95 0.975 0.99
4	-0.87	-0.87	-0.87	1.89 2.08 2.19
5	-1.26	-1.26	-1.20	1.57 1.68 1.77
6	-1.54	-1.54	-1.43	1.46 1.55 1.62
7	-1.75	-1.75	-1.57	1.37 1.46 1.56
8	-1.93	-1.93	-1.67	1.34 1.43 1.51
9	-2.07	-2.07	-1.74	1.32 1.40 1.49
10	-2.19	-2.19	-1.79	1.30 1.38 1.48
11	-2.30	-2.30	-1.86	1.27 1.37 1.47
12	-2.39	-2.39	-1.87	1.26 1.36 1.47
13	-2.49	-2.49	-1.91	1.29 1.37 1.47
14	-2.54	-2.54	-1.92	1.28 1.36 1.47
15	-2.61	-2.61	-1.93	1.27 1.36 1.47
16	-2.68	-2.68	-1.96	1.27 1.35 1.48
17	-2.75	-2.75	-1.98	1.26 1.35 1.47
18	-2.79	-2.79	-1.97	1.25 1.36 1.49
19	-2.84	-2.84	-1.97	1.25 1.35 1.49
20	-2.87	-2.87	-2.02	1.29 1.37 1.50
21	-2.93	-2.93	-1.99	1.29 1.37 1.50
22	-2.99	-2.99	-1.96	1.29 1.37 1.50
23	-3.02	-3.02	-1.95	1.29 1.37 1.50
24	-3.04	-3.04	-1.96	1.28 1.37 1.50
25	-3.08	-3.08	-1.95	1.28 1.38 1.51
26	-3.09	-3.09	-1.94	1.28 1.38 1.52
27	-3.11	-3.11	-1.92	1.27 1.38 1.52
28	-3.17	-3.17	-1.95	1.27 1.38 1.52
29	-3.17	-3.17	-1.96	1.27 1.38 1.54
30	-3.18	-3.18	-1.91	1.27 1.39 1.54
32	-3.25	-3.25	-1.94	1.32 1.40 1.54
34	-3.23	-3.23	-1.96	1.31 1.40 1.55
36	-3.28	-3.28	-2.00	1.31 1.41 1.55
38	-3.29	-3.29	-2.05	1.31 1.40 1.57
40	-3.33	-3.33	-1.86	1.31 1.42 1.58
42	-3.34	-3.34	-1.88	1.30 1.42 1.57
44	-3.35	-3.35	-1.86	1.30 1.42 1.59
46	-3.40	-3.40	-1.84	1.30 1.44 1.59
48	-3.37	-3.37	-1.87	1.30 1.44 1.60
50	-3.38	-3.38	-1.88	1.30 1.44 1.61
55	-3.34	-3.34	-1.87	1.31 1.46 1.62
60	-3.41	-3.41	-1.90	1.34 1.46 1.63
65	-3.39	-3.39	-1.87	1.34 1.47 1.64
70	-3.27	-3.27	-1.19	1.33 1.48 1.66
75	-3.32	-3.32	-1.89	1.33 1.49 1.67
80	-3.22	-3.22	-1.91	1.33 1.50 1.68
85	-3.40	-3.40	-1.88	1.33 1.50 1.68
90	-3.27	-3.27	-1.91	1.33 1.51 1.70
95	-3.19	-3.19	-1.94	1.37 1.52 1.70
100	-3.27	-3.27	-1.90	1.37 1.53 1.71

EOF

}
#------------------------------------------------------------------------------------------------#

sub print_alignment_stats
{ 
	# See documentation in Bio::SimpleAlign
	# better see perldoc Bio::Align::AlignI
	my($input_file,$inputformat,$threshold_percent) = @_;
	
	my $stream  = Bio::AlignIO->new(-file => $input_file,   -format => $inputformat);
	my    $aln  = $stream->next_aln();

	# Describe
	 print "## Stats for alignment $input_file with format $inputformat and threshold percent = $threshold_percent\n";
	 print "\tthe alignment length is: ", $aln->length, "\n";
         print "\tis the alignment is flush? ", $aln->is_flush, "\n";
         print "\tthe alignment has ", $aln->no_residues, " residues\n";
         print "\tthe no. of seqs in the alingment ",  $aln->num_sequences, "\n";
         my $percent_identity = sprintf ("%.2f",$aln->percentage_identity);
	 print "\tthe pecentage identity is: $percent_identity\n";
	 my $avg_percent_identity = sprintf ("%.2f",$aln->average_percentage_identity);
	 print "\tthe average pecentage identity is: $avg_percent_identity\n";
	 
         # Analyze
         print $aln->consensus_string($threshold_percent), "\n";
         print $aln->match_line(), "\n";
	 print $aln->gap_line(), "\n";
         # print $aln->cigar_line(), "\n";
	 print $aln->consensus_iupac, "\n"; # if (defined -p) then $prot_seqs == 1;
	 #return $aln->length;
}
#------------------------------------------------------------------------------------------------#
