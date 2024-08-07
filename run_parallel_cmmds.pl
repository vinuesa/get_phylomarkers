#!/usr/bin/env perl

# AUTHORS: Bruno Contreras Moreira, Pablo Vinuesa 
# AIM: run commands in parallel to process multiple files using GNU parallel

use FindBin '$Bin';

my $VERBOSE = 0; 

# check whether an attached parallel binary is available to avoid errors caused by older, 
# non-compatible parallel binaries in the system
my $parallelEXE = "$Bin/bin/";
if($ENV{'OSTYPE'} =~ /darwin/){ $parallelEXE .= "macosx-intel/parallel" }
else{ $parallelEXE .= "linux/parallel" }

if(! -s $parallelEXE){
  $parallelEXE = `which parallel`;
  chomp $parallelEXE;
  if($parallelEXE eq ''){ 
    print "# ERROR: parallel not in place!\n";
    print "# ... you will need to install \"parallel\" first or include it in \$PATH\n";
    print "# ... exiting\n";
    exit(1)
  }
}

if($VERBOSE == 1){ print "# parallelEXE=$parallelEXE\n" }

if(!$ARGV[1]){
	print << "HELP";
  
	$0 usage: requires 2 or max 3 args
	$0 <file extension name> <'command with \$file interpolations'> [no_of_cores/threads]

	example1: $0 faa   'muscle < \$file > \${file%.faa}_musAln.FAA'
	example2: $0 faa   'mafft-linsi \$file > \${file%.faa}_mlinsi.FAA'
	example3: $0 FAA   'FastTree -slow -sprlength 12 -log \${file%.FAA}_FT.log < \$file > \${file%.FAA}_FT.ph'

HELP
	exit(2);
}

my ( $ext, $command ) = ($ARGV[0], $ARGV[1]);

my $total_files = `ls -1 *$ext | wc -l`;
chomp($total_files);

# prepare command for parallel syntax, which is different to pexec's
$command =~ s/\$file/{}/g;
$command =~ s/\$\{/{/g;
$command =~ s/file//g;
# eliminate indicated extensions
$command =~ s/\%\.\*/./g;
$command =~ s/\%([\.\-\_\w]+)/=s\/$1\/\/=/g;

# added --gnu flag for compatibility
if($ARGV[2] && $ARGV[2] > 0){ 
	my $n_of_cores = $ARGV[2];
	$command = "ls -1 *$ext | $parallelEXE --gnu --willcite -j $n_of_cores \"$command\"";
}
else{
	$command = "ls -1 *$ext | $parallelEXE --gnu --willcite \"$command\"";
}

warn "# $command" if($VERBOSE);
open(RUN,"$command |") || die "# cannot run $command\n";
while(<RUN>){
	if(/Thread creation failed/){
		print "\n >>> ERROR: Thread creation failed, stop ...\n";
		exit(3);
	}
}
close(RUN);

print "\ndone processing $total_files files ...\n";
exit(4);
