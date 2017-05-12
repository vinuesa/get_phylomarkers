#!/usr/bin/perl -w

# Use BioPerl's Bio::AlignIO to convert between multiple sequence alignment formats
# Pablo Vinuesa July 2007. vinuesa@ccg.unam.mx

use strict;

chomp(my $HOSTNAME = `hostname`);
if ($HOSTNAME =~ /kayab/){use lib '/state/partition1/space12/PIG/vinuesa/marfil/CVS_code/BioPerl-1.6.1'};

use Bio::AlignIO;           

# 1) Declare variable, get input arguments from the command line and print help if needed

if ($#ARGV < 3)
{
	print "Usage: $0 inputformat infile_ext outputformat outputfile_ext\n
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
