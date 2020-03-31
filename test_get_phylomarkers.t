use strict;
use warnings;
use Test::More tests => 10;
use lib "lib";
use lib "lib/perl/bioperl-1.5.2_102/";
use lib "lib/perl/File::Rename";

ok( eval{ `perl ./add_labels2tree.pl 2>&1` } =~ /two args/ , 'add_labels2tree.pl' );

ok( eval{ `perl ./add_nos2fasta_header.pl 2>&1` } =~ /fasta file/ , 'add_nos2fasta_header.pl' );

ok( eval{ `perl ./concat_alignments.pl 2>&1` } =~ /list of alignment/ , 'concat_alignments.pl' );

ok( eval{ `perl ./convert_aln_format_batch_bp.pl` } =~ /Usage:/ , 'convert_aln_format_batch_bp.pl' );

ok( eval{ `perl ./popGen_summStats.pl` } =~ /usage/ , 'popGen_summStats.pl' );

ok( eval{ `perl ./pal2nal.pl 2>&1` } =~ /Usage:/ , 'pal2nal.pl' );

ok( eval{ `perl ./remove_uninformative_sites_from_aln.pl -h 2>&1` } =~ /Removes/ , 'remove_uninformative_sites_from_aln.pl' );

ok( eval{ `perl ./rename.pl 2>&1` } =~ /Usage/ , 'rename.pl' );

ok( eval{ `perl ./run_parallel_cmmds.pl` } =~ /usage/ , 'run_parallel_cmmds.pl' );

ok( eval{ `cd test_sequences/core_genome && ../../run_get_phylomarkers_pipeline.sh -R 1 -t DNA | grep "wrote file gene_trees"` }, 'run_get_phylomarkers_pipeline.sh' ); 