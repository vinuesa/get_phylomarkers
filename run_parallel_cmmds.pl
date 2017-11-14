#!/usr/bin/env perl

# AUTHORS: Bruno Contreras Moreira, Pablo Vinuesa 
# AIM: run commands in parallel to process multiple files using GNU parallel

my $VERBOSE = 0; 

my $parallelEXE = `which parallel`;
if($parallelEXE eq ''){ 
	print "# ERROR: parallel not in place!\n";
	print "# ... you will need to install \"parallel\" first or include it in \$PATH";
   print "# ... exiting";
   exit(1)
}
elsif(!$ARGV[1]){
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

if($ARGV[2] && $ARGV[2] > 0){ 
	my $n_of_cores = $ARGV[2];
	$command = "ls -1 *$ext | parallel -j $n_of_cores \"$command\"";
}
else{
	$command = "ls -1 *$ext | parallel \"$command\"";
}

warn "# $command" if($VERBOSE);
system("$command");

print "\ndone processing $total_files files ...\n";
exit(3);
