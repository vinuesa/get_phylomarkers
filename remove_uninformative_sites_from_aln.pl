#!/usr/bin/env perl

# Removes uninformative sites from a fasta multiple sequence alignment: Ns, gaps, and invariant sites
# Author: Lee Katz
# https://github.com/lskatz/lyve-SET/blob/master/scripts/removeUninformativeSites.pl

use strict;
use warnings;
use Data::Dumper;
use Getopt::Long;

use FindBin '$Bin'; # added by Bruno May2017
use lib "$Bin/lib/perl/bioperl-1.5.2_102/";

use Bio::SeqIO;

sub logmsg{$|++;print STDERR "@_\n"; $|--;}
exit main();
sub main{
  my $settings={};
  GetOptions($settings,qw(help verbose ambiguities-allowed gaps-allowed sort=s)) or die $!;
  die usage() if($$settings{help});
  $$settings{"ambiguities-allowed"} ||=0;
  $$settings{"gaps-allowed"} ||=0;
  $$settings{sort}||="";
  $$settings{sort}=lc($$settings{sort});

  ## read in the fasta file into @seq and %seq, and keep the deflines
  my $in;
  if(defined($ARGV[0]) && -f $ARGV[0]){
    $in=Bio::SeqIO->new(-file=>$ARGV[0]);
  } else {
    $in=Bio::SeqIO->new(-fh=>\*STDIN,-format=>"fasta");
  }
  my($length,$defline,@defline,@seq);
  while(my $seqObj=$in->next_seq){
    push(@seq,$seqObj->seq);
    push(@defline,$seqObj->id);
    $length=$seqObj->length;
  }
  $in->close;

  ## read informative positions into @pos
  # compare each sequence to the reference sequence (the first sequence)
  my $removeAmbiguities=!$$settings{"ambiguities-allowed"};
  my $removeGaps=!$$settings{"gaps-allowed"};
  my (%aln,@pos);
  #my $numOtherSeq=@seq;
  my $numSeq=@seq;
  my $informativeCount=0;
  POSITION:for(my $j=0;$j<$length;$j++){ 
    # Find a good reference base to compare against
    my $refNt='N'; # a dummy reference base to start
    for(my $i=0;$i<$numSeq;$i++){
      $refNt=uc(substr($seq[$i],$j,1)); # uppercase for comparison
      last if($refNt ne 'N');
    }

    # Start checking if the site is informative or not
    next if($removeAmbiguities && $refNt =~ /[NRYKMSWBHDV]/); # if it's informative, but the ref base isn't good, then skip
    next if($removeGaps && $refNt eq '-');        # if it's informative, but the ref base isn't good, then skip

    # See if the rest of this column shows that it is informative.
    # Cycle through the rest of the members of this column to check.
    my $informative=0; # guilty until proven innocent
    for(my $i=0;$i<$numSeq;$i++){
      my $nt=uc(substr($seq[$i],$j,1)); # compare everything on the uppercase level
      die "ERROR: Sequence $i does not have a nucleotide at position $j! Is the MSA flush?" if(!$nt);
      # It's informative, but you have to continue reviewing
      # the other sequences if you are ignoring N or gap columns.
      # Also, do not consider ambiguities when deciding if this is an invariant site or not
      if($nt ne 'N' && $nt ne $refNt){
        $informative=1;
      }

      # Check to make sure it doesn't have an N anywhere
      next POSITION if($removeAmbiguities && $nt =~ /[NRYKMSWBHDV]/);
      # Check to make sure there is no gap at all
      next POSITION if($removeGaps && $nt eq '-');
    } 
    next if(!$informative);
    push(@pos,$j); # retain this informative position
    
    $informativeCount++;
    logmsg $j if($$settings{verbose});
  }

  ## Print all nucleotides found at informative positions.
  my %seq;
  @seq{@defline}=@seq;
  # sort the sequences by identifier
  my @sortedId;
  if($$settings{sort} eq 'alpha'){
    @sortedId=sort {$a cmp $b} @defline;
  } elsif($$settings{sort} eq 'num'){
    @sortedId=sort {$a <=> $b} @defline;
  } else {
    @sortedId=@defline;
  }
  for (my $i=0;$i<@sortedId;$i++){
    my $id=$sortedId[$i];
    my $sequence=$seq{$id};
    print ">$id\n";
    for my $pos(@pos){
      print substr($sequence,$pos,1);
    }
    print "\n";
  }
  return 0;


  for (my $i=0;$i<@seq;$i++){
    my $sequence=$seq[$i];
    my $id=$defline[$i];
    print ">$id\n";
    for my $pos(@pos){
      print substr($sequence,$pos,1);
    }
    print "\n";
  }

  return 0;
}

sub usage{
  "Removes all the uninformative sites in a multiple sequence alignment fasta file: Ns, gaps, and invariant sites
  An invariant site will also be defined as a site with only ambiguities and no other variations
  Usage: $0 < aln.fasta > informative.fasta
  -v for verbose (technically makes the script slower)
  Using the following two options will allow you to keep a master MSA list at-hand if you are converting all VCFs in a project
  --gaps-allowed Allow gaps in the alignment
  --ambiguities-allowed Allow ambiguous bases in the alignment
  --sort alpha,num Sort the sequences by their deflines.  Values for sort are either ALPHA or NUM (case insensitive)
  "
}
