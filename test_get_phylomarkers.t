# version 2024-04-16
use strict;
use warnings;
use Test::More tests => 24;

use lib "lib";
use lib "lib/perl/bioperl-1.5.2_102/";
use lib "lib/perl/File::Rename";

### Test Perl scripts
# test 1
ok( eval{ `perl ./add_labels2tree.pl 2>&1` } =~ /two arguments/ , 'add_labels2tree.pl' );

# test 2
ok( eval{ `perl ./add_nos2fasta_header.pl 2>&1` } =~ /fasta file/ , 'add_nos2fasta_header.pl' );

# test 3
ok( eval{ `perl ./concat_alignments.pl 2>&1` } =~ /list of alignment/ , 'concat_alignments.pl' );

# test 4
ok( eval{ `perl ./convert_aln_format_batch_bp.pl` } =~ /Usage:/ , 'convert_aln_format_batch_bp.pl [Bio::AlignIO]' );

# test 5
ok( eval{ `perl ./pal2nal.pl 2>&1` } =~ /Usage:/ , 'pal2nal.pl' );

# test 6
ok( eval{ `perl ./popGen_summStats.pl` } =~ /usage/ , 'popGen_summStats.pl' );

# test 7
ok( eval{ `perl ./remove_uninformative_sites_from_aln.pl -h 2>&1` } =~ /Removes/ , 'remove_uninformative_sites_from_aln.pl [Bio::SeqIO]' );

# test 8
ok( eval{ `perl ./rename.pl 2>&1` } =~ /Usage/ , 'rename.pl' );

# test 9
ok( eval{ `perl ./run_parallel_cmmds.pl` } =~ /usage/ , 'run_parallel_cmmds.pl' );

### Test R scripts
# test 10
ok( eval{ `Rscript ./compute_suppValStas_and_RF-dist.R 2>&1` } =~ /Usage/ , 'compute_suppValStas_and_RF-dist.R [ape,gplots,seqinr,plyr,stringr,ggplot2]' ); 

# test 11
ok( eval{ `Rscript ./run_kdetrees.R 2>&1` } =~ /Usage/ , 'run_kdetrees.R [stringr,ape,kdetrees,vioplot]' );


### Test Bash scripts
# test 12
ok( eval{ `bash ./estimate_pangenome_phylogenies.sh 2>&1` } =~ /PARS PARALLELIZATION/ , 'estimate_pangenome_phylogenies.sh [seqboot,pars,consense,nw_reroot,nw_support]' ); 

# test 13
ok( eval{ `bash ./run_parallel_molecClock_test_with_paup.sh 2>&1` } =~ /OUTPUT/, 'run_parallel_molecClock_test_with_paup.sh [Bio::SeqIO,paup]' );

# test 14
ok( eval{ `bash ./run_get_phylomarkers_pipeline.sh 2>&1` } =~ /INVOCATION/, 'run_get_phylomarkers_pipeline.sh [clustalo,FastTree,parallel,Phi,paup,consense,iqtree]' );

# test 15
ok( eval{ `bash ./run_pexec_cmmds.sh 2>&1` } =~ /example4/, 'run_pexec_cmmds.sh [pexec]' );

# test 16
ok( eval{ `bash ./run_test_suite.sh 2>&1` } =~ /LAUNCHING THE CONTAINER/, 'run_test_suite.sh' );

# test 17; which calls check_bash_version, to ensure that we are running with bash >= 4.3
# this test fails while building the docker image, but not when running the container; 
# This test does not fail during Travis CI
# ok( eval{ `bash ./run_get_phylomarkers_pipeline.sh -V 2>&1` } =~ /clustalo/, 'run_get_phylomarkers_pipeline.sh' );


### Test runs of main scripts
## NOTE: use -I 2 (but not higher, or a core may be dumped in some occasions by iqtree -T IQT_threads) instead of AUTO to speed-up tests 
# test 17 default run on DNA sequence ... Builds on docker and locally, fails on travis 
ok( eval{ `cd test_sequences/core_genome && ../../run_get_phylomarkers_pipeline.sh -R 1 -t DNA -I 2 -N2 | grep "markers into supermatrix" && rm -rf get_phylomarkers_run_*` }, 'run_get_phylomarkers_pipeline.sh -R 1 -t DNA -I 2 -N 2' ); 

# test 18 IQT run on proteins sequence ...
ok( eval{ `cd test_sequences/core_genome && ../../run_get_phylomarkers_pipeline.sh -R 1 -t PROT -k 1.5 -m 0.3 -I 2 -N 2 | grep "wrote file" && rm -rf get_phylomarkers_run_*` }, 'run_get_phylomarkers_pipeline.sh -R 1 -t PROT -k 1.5 -m 0.3 -I 2 -N 2' ); 

# test 19 Thorough FastTree searching and molecular clock analysis on DNA sequences using 10 cores and increasing k stringency
ok( eval{  `cd test_sequences/core_genome && ../../run_get_phylomarkers_pipeline.sh -R 1 -t DNA -A F -k 1.2 -m 0.7 -s 20 -l 12 -T high -K -M HKY -q 0.95 | grep "running kde test"` }, 'run_get_phylomarkers_pipeline.sh -R 1 -t DNA -A F -k 1.2 ...' );

# test 20 FastTree thorough searching on a protein dataset with moderate average bipartition support
ok( eval{  `cd test_sequences/core_genome && ../../run_get_phylomarkers_pipeline.sh -R 1 -t PROT -A F -T high -m 0.2 | grep "writing summary tables"` }, 'run_get_phylomarkers_pipeline.sh -R 1 -t PROT -A F ...' );

# test 21  Run in population-genetics mode under K2P model, using IQ-Tree
ok( eval{ `cd test_sequences/core_genome && ../../run_get_phylomarkers_pipeline.sh -R 2 -t DNA -S K2P -k 1.3 -I 2 | grep "wrote file concat_cdnAlns_SNPs.fasta"` }, 'run_get_phylomarkers_pipeline.sh -R 2 -t DNA -S K2P -k 1.3 -I 2' );

# test 22  Run in population-genetics mode, estimating population tree using FastTree
ok( eval{ `cd test_sequences/core_genome && ../../run_get_phylomarkers_pipeline.sh -R 2 -t DNA -A F -k 1.2 | grep "wrote file concat_cdnAlns_SNPs.fasta"` }, 'run_get_phylomarkers_pipeline.sh -R 2 -t DNA -A F -k 1.2' );


# test 23 estimate a ML pan-genome tree from the pan-genome matrix, using 2 independent IQT runs and UFBoot
my $testOK = ok( eval{ `cd test_sequences/pan_genome && ../../estimate_pangenome_phylogenies.sh -f pangenome_matrix_t0.fasta -r 2 -S UFBoot -I 2 | grep "done!"` }, 'estimate_pangenome_phylogenies.sh -r 2 -S UFBoot ... [libnw.so]' );

if(!$testOK) {
  print "\n# Note: this test requires setting up libnw.so . Type the following:\n\n"; 
  print "sudo cp /PATH/TO/get_phylomarkers/lib/libnw.so /usr/local/lib && sudo echo \"export LD_LIBRARY_PATH=/usr/local/lib\" && sudo ldconfig && make clean'\n\n";
  print "# Read more at https://github.com/eead-csic-compbio/get_phylomarkers/blob/master/INSTALL.md\n\n";
}

# test 24 estimate a PARS  pan-genome tree with bootstrapping; 50 bootstrap replicates divided on 10 core (5 reps / core)
ok( eval{ `cd test_sequences/pan_genome && ../../estimate_pangenome_phylogenies.sh -c PARS -R 3 -i pangenome_matrix_t0.phylip -b 5 -j 1 -t 1 | grep "seq_key"` }, 'estimate_pangenome_phylogenies.sh -c PARS -R 3 ...' );
