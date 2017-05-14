#!/usr/bin/env perl 

# Pablo Vinuesa
# Centro de Ciencias Genomicas, UNAM, Mexico
# Jul 26th, 2010

# add_nos2fasta_header.pl, adds consecutive numbers after '>' in fasta files
# translated from the onliner: 
# perl -pe 'if(/^>/){$counter++ ; s/>/>$counter /}' infile | grep '>'

use File::Basename;

$progname = basename($0); #add_nos2fasta_header.pl
$version = 1.0;

die "\n# $progname v.$version needs a fasta file name provided as single argument\n"
  unless @ARGV == 1;
while(<>)
{
   if(/^>/)
   {
       $counter++; 
       s/>/>$counter /
   }
}
continue
{
   print
}

