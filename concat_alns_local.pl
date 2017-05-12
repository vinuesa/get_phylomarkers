#!/usr/bin/perl -w
# Pablo Vinuesa, Bruno Contreras Moreira
# 2007 UNAM, Mexico


use strict;
use File::Basename;
use lib '/home/vinuesa/bin';
use phyTools;
set_phyTools_env();

my $progname = "concat_alns_local.pl";
my $REMOVETMPFILES = 1;
my $EVOLRATE_CUTOFF = $ENV{'ABS_MAX_SITE_RATE'};

######################################################

my ( $listfile, @aln_files, $file, $seq, %concat_aln );
my ( $length , $strip_gaps, %models, $removeFES, $r4s_info, $info );

######################################################

$removeFES=$strip_gaps=0;

if(!$ARGV[0]){ die "\n\n# usage: $progname <list of alignment filenames [tab model, optional, only for remove_fast_sites()], one per line> <options: [g=remove gaps OR r=remove fast evolving sites]>\n\n\n"; }
if($ARGV[1])
{ 
	$listfile = $ARGV[0]; 
	if($ARGV[1] eq 'g'){ $strip_gaps = 1; }
	elsif($ARGV[1] eq 'r'){ $removeFES = 1; }
}
else{ ($listfile) = @ARGV;  }

if(! -s $listfile){ die "#$progname : need a valid list file\n"; }

if( $removeFES == 1){ $strip_gaps = 1 }

print "# $progname params: input file $listfile g=$strip_gaps r=$removeFES rate=$EVOLRATE_CUTOFF\n";

######################################################

# 1) read list file
open(LIST,$listfile) || die "# $progname : cannot read $listfile\n";
while(<LIST>)
{
	next if(/^#/);
	$file = (split)[0];
	push(@aln_files, $file);
	if((split)[1]){ $models{$file} = (split)[1]; }
}
close(LIST);  

# 2) loop reading aln files and concatenate sequences
foreach $file (@aln_files)
{
	my %FASTAaln;

	if($removeFES)
	{
		my %gFASTAaln = read_FASTA_sequence( $file , 1 );
		($info,%FASTAaln) = filter_aln_rate4site( $file, $models{$file}, \%gFASTAaln , $EVOLRATE_CUTOFF );
		$r4s_info .= $info;
	}
	elsif($strip_gaps)
        {
                %FASTAaln = read_FASTA_sequence( $file , 1 );
	}
	else
	{
	 	%FASTAaln = read_FASTA_sequence( $file );
	} 
	
	foreach $seq (keys(%FASTAaln))
	{
		if(!$concat_aln{$seq}{'LENGTH'})	{ $concat_aln{$seq}{'LENGTH'} = 0; }

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
		if($models{$file}){ $concat_aln{$seq}{'COORDS'} .= $models{$file} .", " }
		$concat_aln{$seq}{'COORDS'} .= basename($file) . " = " . ($concat_aln{$seq}{'LENGTH'}+1) . "-" . ($length + $concat_aln{$seq}{'LENGTH'}) . "\n";
		
		# update lengths
		$concat_aln{$seq}{'LENGTH'} += $length;
	}
}	

# 3) print concatenation information for removed FFES (fast-evolving sites)
if( $removeFES )
{
	print "$r4s_info\n";
}

print "# concatenation coordinates:\n$concat_aln{1}{'COORDS'}\n"; 
# print "# concat_aln_seq_lengts are: \n$concat_aln{1}{'LENGTH'} \n"; ???? doesn't print

## 4) print concatenated alignment with names collapsed to IDs
#foreach $seq (sort {$a<=>$b} (keys(%concat_aln)))
#{
#	print ">$seq\n$concat_aln{$seq}{'SEQ'}\n";
#}

# 4) print concatenated alignment with full names
# my $outfile = $concat_aln{1}{'LENGTH'}_concat
#open(OUT,">$outfile") || die "$0 can't write to outputfile $outfile :$!\n";
foreach $seq (sort {$a<=>$b} (keys(%concat_aln)))  # numeric sort of seq IDs
{
	print "$concat_aln{$seq}{'NAME'}\n$concat_aln{$seq}{'SEQ'}\n"; # $concat_aln{$seq}{'NAME'} carries ">" along
	# print OUT $concat_aln{$seq}{'NAME'},"\n"$concat_aln{$seq}{'SEQ'},"\n";
}
close(OUT);
