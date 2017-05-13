#!/usr/bin/env perl
# Pablo Vinuesa, Bruno Contreras Moreira
# 2007 UNAM, Mexico

use strict;
use warnings;
use File::Basename;

my $progname = basename($0);
my $REMOVETMPFILES = 1;
my $VERSION = 1.0; # May 13th, 2017. Stripped down version, that does not depend on phyTools

######################################################

my ( $listfile, @aln_files, $file, $seq, %concat_aln );
my ( $length);

######################################################

if(!$ARGV[0]){ die "\n\n# usage: $progname <list of alignment filenames>\n\n"; }
if($ARGV[1])
{ 
	$listfile = $ARGV[0]; 
}
else{ ($listfile) = @ARGV;  }

if(! -s $listfile){ die "#$progname : need a valid list file\n"; }

#print "# $progname params: input file $listfile\n";

######################################################

# 1) read list file
open(LIST,$listfile) || die "# $progname : cannot read $listfile\n";
while(<LIST>)
{
	next if(/^#/);
	$file = (split)[0];
	push(@aln_files, $file);
}
close(LIST);  

# 2) loop reading aln files and concatenate sequences
foreach $file (@aln_files)
{
    my %FASTAaln;
    %FASTAaln = read_FASTA_sequence( $file );

    foreach $seq (keys(%FASTAaln))
    {
    	    if(!$concat_aln{$seq}{'LENGTH'}){ $concat_aln{$seq}{'LENGTH'} = 0; }

    	    # concat headers ommited in concat_aln_local.pl 
    	    #$concat_aln{$seq}{'NAME'} .= $FASTAaln{$seq}{'NAME'};
    	    $concat_aln{$seq}{'NAME'} = $FASTAaln{$seq}{'NAME'}; #concatenation ommited
    	    chomp $concat_aln{$seq}{'NAME'};

    	    # concat sequence
    	    $concat_aln{$seq}{'SEQ'} .= $FASTAaln{$seq}{'SEQ'};
    	    chomp $concat_aln{$seq}{'SEQ'};
    	    $concat_aln{$seq}{'SEQ'} =~ s/\*//g;
    	    
    	    # store coordinates 
    	    $length = length($FASTAaln{$seq}{'SEQ'});
    	    $concat_aln{$seq}{'COORDS'} .= basename($file) . " = " . ($concat_aln{$seq}{'LENGTH'}+1) . "-" . ($length + $concat_aln{$seq}{'LENGTH'}) . "\n";
    	    
    	    # update lengths
    	    $concat_aln{$seq}{'LENGTH'} += $length;
    }
}	

# 3) print concatenated alignment with full names to STDOUT
my $coords_file="concatenation_coordinates.txt";
open COORDS, ">", "$coords_file" or die "can't write to $coords_file: $!\n\n";
print COORDS "# concatenation coordinates:\n$concat_aln{1}{'COORDS'}\n"; 
close COORDS;

foreach $seq (sort {$a<=>$b} (keys(%concat_aln)))  # numeric sort of seq IDs
{
	print "$concat_aln{$seq}{'NAME'}\n$concat_aln{$seq}{'SEQ'}\n"; # $concat_aln{$seq}{'NAME'} carries ">" along
}

#---------------------#
# >>> SUBROUTINES <<< #
#---------------------#

sub read_FASTA_sequence
{
    my ($infile) = @_;
	 
    my (%FASTA,$name,$seq,$n_of_sequences);

    $n_of_sequences = 0;
    open(FASTA,$infile) || die "# read_FASTA_sequence: cannot read $infile\n";
    while(<FASTA>)
    {
    	 if(/^\>/)
	 {
	       $name = $_; 
	       $n_of_sequences++;
	       $FASTA{$n_of_sequences}{'NAME'} = $name;
	 }		       
	 else
	 {
	       $FASTA{$n_of_sequences}{'SEQ'} .= $_;
	       $FASTA{$n_of_sequences}{'SEQ'} =~ s/[\s|\n]//g;
	 }
    }
    close(FASTA);

    return %FASTA;
}
