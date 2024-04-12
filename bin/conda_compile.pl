use strict;
use FindBin qw($Bin);
use Config;

my $move_binaries = 0;
if($ARGV[0]) { 
	$move_binaries = 1;
}

my ($targetdir,$cmd) = ('');

if($Config{osname} =~ /linux/) {
	$targetdir = $Bin.'/linux/';
} #elsif($Config{osname} =~ /darwin/) { $targetdir = $Bin.'/macosx-intel/' }

## phylip binaries (3)
chdir("$Bin/phylip-3.695/src");
if($Config{osname} =~ /linux/) {
	$cmd = "make -f Makefile.unx consense pars seqboot";
} #elsif($Config{osname} =~ /darwin/) { $cmd = "make -f Makefile.osx consense pars seqboot" }

system("$cmd");
if ($? == 0 && $move_binaries == 1) {
	system("mv consense pars seqboot $targetdir");
}
chdir($Bin);

## FastTree (1)
chdir("$Bin/FastTree_v2.1.11");
if($Config{osname} =~ /linux/) {
	$cmd ="gcc -DUSE_DOUBLE -O3 -finline-functions -funroll-loops -Wall -o FastTree FastTree-2.1.11.c -lm"
} 

system("$cmd");
if ($? == 0 && $move_binaries == 1) {
	system("mv FastTree $targetdir")
}
chdir($Bin);

# PhiPack (1)
chdir("$Bin/PhiPack/src");
if($Config{osname} =~ /linux/) {
	$cmd ="make Phi"
}

system("$cmd");
if ($? == 0 && $move_binaries == 1) {
        system("mv Phi $targetdir")
}
chdir($Bin);

# ASTER (2)
chdir("$Bin/ASTERcommit1840895");
if($Config{osname} =~ /linux/) {
	$cmd ="make astral wastral"
}

system("$cmd");
if ($? == 0 && $move_binaries == 1) {
	system("mv wastral $targetdir");
	system("mv astral4 $targetdir/astral");
}
chdir($Bin);


