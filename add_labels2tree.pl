#!/usr/bin/env perl

# add_labels2tree.pl
# by Pablo Vinuesa, CCG-UNAM, Mexico. 

# Used to edit the taxon labels on newick phylogenies generated by PhyML or CONSENSE
# converting numbers to full taxon labels. The latter are read from the 
# fasta headers of the original *.faa files, or from the my_*.faa files
# in the formats of GABO, Neurogadgets and IH-PVs

# usage: run from within the directory containing the tree files,
# providing two args: 1.- an *.faa file and 2.- the extension of the topologies to edit phyml_tree.txt

# perl add_labels2tree.pl my_rpoA_NCTC11168_aln.faa phyml_tree.txt

use warnings;
use strict;
use File::Basename;
my $progname = basename($0); # add_labels2tree.pl
my $VERSION = 1.7; # v1.7 2022-01-08; fixed regex for ASTRAL
          # v1.6 June 12th, 2022; fixed regexes for ASTRAL (topology-only) tree; updated POD and HELP
          # May 14th, 2017. added portable shebang line

          # Sept 15th, 2013; added general case (/^>?(\d+)[\t\s+](.*$); again allows >ID [bichococcus ..] type labels to work!
          # Jul 10th, 2013; added the general case: (/(\S+)[\t\s+](\S+)/)
          # Aug. 4th, 2012; 
          # 9 Feb 2010 ahora lee tambi�n elsif(/^>?(\d+)\s+(.*)/)
          # v 1.0 ahora lee if(/^>(\d+)\s+\[(.*?)\]/) y  elsif(/^(\d+)\s+\[(.*?)\]/) ... es decir >ID [blah] y ID [blah]

die "# $progname v. $VERSION requires two argyments: 
     1) the name of the 2-column, tab-separated ID\tfull_name conversion table or an alignment file with full names
     2) the tree file name to edit, or an extension name of tree files on which to edit the taxon tags
     for example:
     $0 rpoC-Rhizobiales_extended_aln.faa rpoC-Rhizobiales_extended_aln_nj.ph
     or
     $0 tree_labels.list astral.tre\n" unless @ARGV ==2;
     
my($aln_file, $tree_file_ext) = (@ARGV);
my (@tree_files, $tree_file, $species, $strain,$tax_label, $seq_key, %seq_ID, $basename, $outfile, $ext, $gi);

open ALNFILE, $aln_file or die "$0 can't open aln_file $aln_file: $!\n";
print "the seq_key => seq_ID correspondences for aln $aln_file are:\n";
while(<ALNFILE>){
    chomp; 
    if(/^>(\d+)\h+\[(.*?)\]/){ 
        #$GI_number = $1;
	$seq_key = $1;
	$strain = $2;
	$seq_ID{$seq_key}=$strain;
	print "$seq_key => $seq_ID{$seq_key}\n";
    }
    elsif(/^(\d+)\h+(\S+)/)
    {
        #$GI_number = $1;
	$seq_key = $1;
	$strain = $2;
	$seq_ID{$seq_key}=$strain;
	print "$seq_key => $seq_ID{$seq_key}\n";	
     }
    elsif(/^(\d+)\h+\[(.*?)\]/)
    {
        #$GI_number = $1;
	$seq_key = $1;
	$strain = $2;
	$seq_ID{$seq_key}=$strain;
	print "$seq_key => $seq_ID{$seq_key}\n";	
     }
    elsif(/^>?(\S+)[\t\s+](\S+)/)
    {
        #$GI_number = $1;
	$seq_key = $1;
	$strain = $2;
	$seq_ID{$seq_key}=$strain;
	print "$seq_key => $seq_ID{$seq_key}\n";	
     }
    elsif(/^>?(\d+)[\t\s+](.*$)/)
    {
        #$GI_number = $1;
	$seq_key = $1;
	$strain = $2;
	$strain =~ s/[\s\|;-]+/_/g;
	$seq_ID{$seq_key}=$strain;
	print "$seq_key => $seq_ID{$seq_key}\n";	

     }
    elsif(/^>?(\d+)\h+\|\d+\|\[(.*)\]/)
    {
        #$GI_number = $1;
	$seq_key = $1;
	$strain = $2;
	$strain =~ s/[\s\|;-]+/_/g; #print "# strain label is $strain\n"; exit;
	$seq_ID{$seq_key}=$strain;
	print "$seq_key => $seq_ID{$seq_key}\n";	

     }
     elsif(/^>(\d+) \d+\|\[(.*)\]\|(.*)\|/)
     {
        #$GI_number = $1;
	$seq_key = $1;
	$species = $2;
	$strain = $3;
	$tax_label = $species . '_' . $strain;
	$tax_label =~ s/[\s+|;|-]/_/g;
	$seq_ID{$seq_key}=$tax_label;
	print "$seq_key => $seq_ID{$seq_key}\n"; # format of our extract_CDSfromGBK_batch.pl
     }elsif(/^>(\d+)\s+\S+\s+\|\[(.*?)\]/){ # Fasta header for get_homologues.pl
        $seq_key=$1;
	$strain = $2;
	$seq_ID{$seq_key}=$strain;
	print "$seq_key => $seq_ID{$seq_key}\n";	
    }elsif(/^>(\d+) \d+_\[(\S+)\]/){
        #$GI_number = $1;
	$seq_key = $1;
	$strain = $2;
	$seq_ID{$seq_key}=$strain;
	print "$seq_key => $seq_ID{$seq_key}\n";
    }elsif(/^>\d+_(\S+)/){
        #$GI_number = $1;
	$seq_key++;
	$strain = $1;
	$seq_ID{$seq_key}=$strain;
	print "$seq_key => $seq_ID{$seq_key}\n";
    }elsif(/^>(\w+)\|/){
        #$GI_number = $1;
	$seq_key++;
	$strain = $1;
	$seq_ID{$seq_key}=$strain;
	print "$seq_key => $seq_ID{$seq_key}\n";
    }elsif(/^>\d+\s+\S+\s+\[(\S+)\]/){ #GABO's fasta header
        $seq_key++;
	$strain = $1;
	$seq_ID{$seq_key}=$strain;
	print "$seq_key => $seq_ID{$seq_key}\n";
     }elsif(/^>(\d+) \d+ \S* (\d+)\|\[(\w+\s*\w*\s*\w*)\]/){ #Marfil's fasta header after filtering redundant sequences
        $seq_key = $1; #print "saw seq_key $seq_key\n";
	$gi      = $2;
	$strain  = $3; #print "saw strain $strain\n";
	$strain  =~ s/\s+/_/g;
	$strain  = $strain . '_' . $gi;
	$seq_ID{$seq_key}=$strain;
	print "$seq_key => $seq_ID{$seq_key}\n";
    }elsif (/^>(\w+) (.*)/){  # RDPII fasta header; most general fasta header ID w+ w+ w+ ...
        #$GI_number = $1;
	$seq_key = $1;
	$strain = $2;
	$strain =~ s/\s+/_/g;
	$strain =~ s/,/_/g;
	$strain =~ s/:/_/g;
	$strain =~ s/".*"/_/g;
	$strain =~ s/\(//g;
	$strain =~ s/\)//g;
	$seq_ID{$seq_key}=$strain;
	print "$seq_key => $seq_ID{$seq_key}\n";	
    }else{
        next;
    }
}
print "\n";


@tree_files = <*$tree_file_ext>;
foreach $tree_file(@tree_files){
    ($basename, $ext) = (split(/\./, $tree_file))[0,1];
    $outfile = $basename . "_ed." . "$ext";
    open TREEFILE, $tree_file or die "$0 can't open aln_file $aln_file: $!\n";
    while(<TREEFILE>){
	foreach $seq_key (sort keys %seq_ID){
	    open OUT, ">$outfile" or die "$0 can't write to outfile $outfile: $!\n";
	    if ($_ =~ /\($seq_key:/){
	         print STDERR "the seq_key is $seq_key\n";
	        s/\($seq_key:/\($seq_ID{$seq_key}:/g;
		 print STDERR "the seq_key is $seq_key\n";
	    }elsif ($_ =~ /,$seq_key:/){
	       s/,$seq_key:/,$seq_ID{$seq_key}:/g;
	       print STDERR "the seq_key is $seq_key\n";
	    }elsif  ($_ =~ /\(($seq_key)(\__\d+):/){    #if we used Fasta2cluster, with tree labels such as 000000001__3
	        s/\($1$2:/\($seq_ID{$seq_key}$2:/g;     #if we used Fasta2cluster, with tree labels such as 000000001__3
	    }elsif ($_ =~ /,($seq_key)(\__\d+):/){      #if we used Fasta2cluster, with tree labels such as 000000001__3
	       s/,$1$2:/,$seq_ID{$seq_key}$2:/g;        #if we used Fasta2cluster, with tree labels such as 000000001__3
	    }
	    
	    # for topology-only trees, such as those produced by ASTRAL
	    elsif ($_ =~ /\(($seq_key,)/){                 
	       s/\($seq_key,/\($seq_ID{$seq_key},/g; 
	        print STDERR "the seq_key is $seq_key\n"; 
	    }
	    elsif ($_ =~ /(,$seq_key\))/){ 
	        s/,$seq_key\)/,$seq_ID{$seq_key}\)/g;
		 print STDERR "the seq_key is $seq_key\n";
	    }
#	   elsif ($_ =~ /\)(,$seq_key)/) {
#	      s/$1/,$seq_ID{$seq_key}/g;
#	   }
#	   elsif ($_ =~ /(,$seq_key)\)/) {
#	      s/$1/,$seq_ID{$seq_key}/g;
#	   }  
 
	}
        print OUT;
	print "the edited topology for treefile $tree_file is:\n $_\n";		
    }
    
    close TREEFILE;
    close OUT;
}

=head1 NAME
 
<add_labels2tree.pl> - <Used for editing of the taxon labels on a newick phylogeny, converting seqID numbers to full taxon labels>
 
 
=head1 VERSION
 

This documentation refers to <add_labels2tree.pl> version 1.7
 
 
=head1 USAGE
 
    usage: run from within the directory containing the tree files, providing two args: 
    1.- tsv conversion/correspondence table | my_alignment | Fasta_File.faa  
    2.- the file name | extension file name of the topologies to edit, like phyml_tree.txt or aln.ph 
        The latter should complement the basename of the aln files
    
 Usage example:  
  path/to/add_labels2tree.pl my_list.tsv tree.ph
 
 
=head1 REQUIRED ARGUMENTS
 
	1) name of tab-separated, 2 column correspondences file | alignment file
	2) tree name or tree extension name, like phyml_tree.txt or aln.ph
 
 
=head1 DESCRIPTION
 
The program will process all trees in a directory, changing the numeric sequence identifiers to full taxon labels using the alignment provided

=head1 DIAGNOSTICS
 
Just standard diagnostics
 
 
=head1 DEPENDENCIES
 
No external dependiencies.

=head1 BUGS AND LIMITATIONS
 
There are no known bugs in this module. 
Please report problems to <Maintainer name(s)>  (<vinuesa@ccg.unam.mx>)
Patches are welcome.
 
=head1 AUTHOR
 
<Pablo Vinuesa>  (<vinuesa@ccg.unam.mx>)
 

=head1 LICENCE AND COPYRIGHT
 
Copyright (c) <2023> <Pablo Vinuesa> (<vinuesa@ccg.unam.mx>). All rights reserved.
 
This program is free software; you can redistribute it and/or
modify it under the same terms as Perl itself. See L<perlartistic>.
 
This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
