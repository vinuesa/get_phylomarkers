#!/usr/bin/env perl

# Use BioPerl's Bio::AlignIO to convert between multiple sequence alignment formats
# Pablo Vinuesa July 2007. vinuesa@ccg.unam.mx

use strict;
use warnings;
use File::Basename;

my $progname = basename($0);
my $VERSION = 0.2;  # May 14th, 2017. Added portable shebang; some cleanup

use FindBin '$Bin'; # added by Bruno May2017
use lib "$Bin/lib/perl/bioperl-1.5.2_102/";
use Bio::AlignIO;           

# 1) Declare variable, get input arguments from the command line and print help if needed

if ($#ARGV < 3)
{
	print "Usage: $progname inputformat infile_ext outputformat outputfile_ext\n
	          Supported formats include:

              bl2seq      Bl2seq Blast output
              clustalw    clustalw (.aln) format
              emboss      EMBOSS water and needle format
              fasta       FASTA format
              maf         Multiple Alignment Format
              mase        mase (seaview) format
              mega        MEGA format
              meme        MEME format
              msf         msf (GCG) format
              nexus       Swofford et al NEXUS format
              pfam        Pfam sequence alignment format
              phylip      Felsenstein PHYLIP format
              prodom      prodom (protein domain) format
              psi         PSI?BLAST format
              selex       selex (hmmer) format
              stockholm   stockholm format\n\n\n";
	exit;
}    

my($inputformat, $infile_ext, $outputformat, $outputfile_ext)=@ARGV;
my @infiles = <*$infile_ext>;
my($basename,$counter);


# 2) Process all files in cwd having infile_ext, converting them from inputformat to outputformat 
#    with the outputfile_ext provided at the command line

foreach my $infile ( @infiles )
{  
	$basename = (split(/\./, $infile))[0]; 
		
	my $in  = Bio::AlignIO->new(-file => $basename . ".$infile_ext",   -format => $inputformat);
	my $out = Bio::AlignIO->new(-file => ">$basename.$outputfile_ext", -format => $outputformat);

	while ( my $aln = $in->next_aln() ) 
	{
		    $out->write_aln($aln);
	}
	$counter++;
}
print "\n\t# I'm done: $counter $inputformat input files with $infile_ext extension were converted to $outputformat format with $outputfile_ext extension\n\n";
