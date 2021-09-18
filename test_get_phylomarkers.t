# version 2021-09-16; vivaMX
use strict;
use warnings;
use Test::More tests => 22;

use lib "lib";
use lib "lib/perl/bioperl-1.5.2_102/";
use lib "lib/perl/File::Rename";

### Test Perl scripts
# test 1
ok( eval{ `perl ./add_labels2tree.pl 2>&1` } =~ /two args/ , 'add_labels2tree.pl' );

# test 2
ok( eval{ `perl ./add_nos2fasta_header.pl 2>&1` } =~ /fasta file/ , 'add_nos2fasta_header.pl' );

# test 3
ok( eval{ `perl ./concat_alignments.pl 2>&1` } =~ /list of alignment/ , 'concat_alignments.pl' );

# test 4
ok( eval{ `perl ./convert_aln_format_batch_bp.pl` } =~ /Usage:/ , 'convert_aln_format_batch_bp.pl' );

# test 5
ok( eval{ `perl ./popGen_summStats.pl` } =~ /usage/ , 'popGen_summStats.pl' );

# test 6
ok( eval{ `perl ./pal2nal.pl 2>&1` } =~ /Usage:/ , 'pal2nal.pl' );

# test 7
ok( eval{ `perl ./remove_uninformative_sites_from_aln.pl -h 2>&1` } =~ /Removes/ , 'remove_uninformative_sites_from_aln.pl' );

# test 8
ok( eval{ `perl ./rename.pl 2>&1` } =~ /Usage/ , 'rename.pl' );

# test 9
ok( eval{ `perl ./run_parallel_cmmds.pl` } =~ /usage/ , 'run_parallel_cmmds.pl' );

### Test R scripts
# test 10
ok( eval{ `Rscript ./compute_suppValStas_and_RF-dist.R 2>&1` } =~ /Usage/ , 'compute_suppValStas_and_RF-dist.R' ); 


### Test Bash scripts
# test 11
ok( eval{ `bash ./estimate_pangenome_phylogenies.sh 2>&1` } =~ /NOTES ON/ , 'estimate_pangenome_phylogenies.sh' ); 

# test 12
ok( eval{ `bash ./hcluster_pangenome_matrix.sh 2>&1` } =~ /IMPORTANT NOTES/, 'hcluster_pangenome_matrix.sh' );

# test 13
ok( eval{ `bash ./run_get_phylomarkers_pipeline.sh 2>&1` } =~ /INVOCATION EXAMPLES/, 'run_get_phylomarkers_pipeline.sh' );

# test 14
ok( eval{ `bash ./run_parallel_molecClock_test_with_paup.sh 2>&1` } =~ /OUTPUT/, 'run_parallel_molecClock_test_with_paup.sh' );

# test 15
ok( eval{ `bash ./run_pexec_cmmds.sh 2>&1` } =~ /example4/, 'run_pexec_cmmds.sh' );

# test 16
ok( eval{ `bash ./run_test_suite.sh 2>&1` } =~ /LAUNCHING THE CONTAINER/, 'run_test_suite.sh' );

### Test runs of main scripts
# test 17 default run on DNA sequence ...
ok( eval{ `cd test_sequences/core_genome && ../../run_get_phylomarkers_pipeline.sh -R 1 -t DNA | grep "wrote file gene_trees" && rm -rf get_phylomarkers_run_*` }, 'run_get_phylomarkers_pipeline.sh -R 1 -t DNA' ); 

# test 18 Thorough FastTree searching and molecular clock analysis on DNA sequences using 10 cores and increasing k stringency
ok( eval{  `cd test_sequences/core_genome && ../../run_get_phylomarkers_pipeline.sh -R 1 -t DNA -A F -k 1.2 -m 0.7 -s 20 -l 12 -T high -K -M HKY -q 0.95 -n 10 | grep "wrote file phylogenetic_attributes_of_top"` }, 'run_get_phylomarkers_pipeline.sh -R 1 -t DNA -A F -k 1.2 ...' );

# test 19 FastTree thorough searching on a protein dataset with moderate average bipartition support
ok( eval{  `cd test_sequences/core_genome && ../../run_get_phylomarkers_pipeline.sh  -R 1 -t PROT -A F -T high -m 0.6 | grep "wrote file concat_nonRecomb_KdeFilt_protAlns_FTlgG"` }, 'run_get_phylomarkers_pipeline.sh -R 1 -t PROT -A F ...' );

# test 20  Run in population-genetics mode (generates a table with descritive statistics for DNA-polymorphisms) with K2P model
ok( eval{ `cd test_sequences/core_genome && ../../run_get_phylomarkers_pipeline.sh -R 2 -t DNA -S K2P | grep "wrote file polymorphism_descript_stats.tab"` }, 'run_get_phylomarkers_pipeline.sh -R 2 -t DNA -S K2P' );

# test 21 estimate a ML pan-genome tree from the pan-genome matrix, using 2 independent IQT runs and UFBoot
ok( eval{ `cd test_sequences/pan_genome && ../../estimate_pangenome_phylogenies.sh -f pangenome_matrix_t0.fasta -r 2 -S UFBoot | grep "done!"` }, 'estimate_pangenome_phylogenies.sh -r 2 -S UFBoot ...' );

# test 22 estimate a PARS  pan-genome tree with bootstrapping; 50 bootstrap replicates divided on 10 core (5 reps / core)
ok( eval{ `cd test_sequences/pan_genome && ../../estimate_pangenome_phylogenies.sh -c PARS -R 3 -i pangenome_matrix_t0.phylip -n 10 -b 5 -j 1 -t 1 | grep "wrote file full_pars_tree_rooted_withBoot.ph"` }, 'estimate_pangenome_phylogenies.sh -c PARS -R 3 ...' );
