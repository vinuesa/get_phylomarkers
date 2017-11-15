# Installation and execution notes for the get_phylomarkers pipeline

Version: Nov, 14 2017
 
This file lists the software components of the *get_phylomarkers* pipeline and briefly describes how to install them.

The pipeline runs on Linux (Ubuntu and RedHat distros) and Mac OS X environments.

It assumes that recent versions of Bash, Perl and R are installed on a multicore (64bit) machine, ideally a server running Linux.
The pipeline is designed to take advantage of modern multicore machines to parallelize all repetitive tasks that have to be performed on each entry sequence, like generating multiple sequence alignments, deriving codon alignments, inferring maximum likelihood gene phylogenies and testing their clock-likeness. Therefore, if your intention is to select optimal genome markers to infer genome phylogenies, you should run the pipeline on a multiprocessor/multicore server to speed up computations. 


## Quick install and test notes

1. clone the repository into a suitable directory (e.g. $HOME/src/gitHub/) using the command 'git clone https://github.com/vinuesa/get_phylomarkers.git'

2. Make sure R is configured in your system. If not, please install R package (r-base in linux). See CRAN packages for OSX [here](https://cran.r-project.org/bin/macosx/).

3. cd into get_phylomarkers/ and run './install_R_deps.R', which will install R packages into get_phylomarkers/lib/R

4. Copy the test_sequences directory into a suitable place (e.g. 'cp -r test_sequences $HOME)

5. cd into the test_sequences dir (e.g. '$HOME/test_sequences')

6. Issue the following command to test if the distro is working on your system: '/path/to/get_phylomarkers/run_get_phylomarkers_pipeline.sh -R 1 -t DNA -K 1', which will run in phylogenomics mode (-R 1), on DNA sequences (-t DNA), and will test the molecular clock hypothesis (-K 1). 
 
7. Check it now on the protein level: 'run_get_phylomarkers_pipeline.sh -R 1 -t PROT'. Note that for this second invocation, you will probably not need to prepend the full path to the script anymore (see NOTES below).

8. Explore the help menu of the main script to see the options available for customization of the run. It is printed to STDOUT when issuing run_get_phylomarkers_pipeline.sh -h or simply run_get_phylomarkers_pipeline.sh

That's it, enjoy. 


## Scripts distributed through GitHub

NOTES: 

1. The main script run_get_phylomarkers_pipeline.sh will automatically identify the directory
on the local machine where the get_phylomarkers package was downloaded. It will also check if the host machine has a 
\$HOME/bin dir included in \$PATH. If so, the main script will automatically generate 
symlinks in \$HOME/bin to the package scripts, so that they become visible system-wide.
Otherwise, it will append the distribution directory holding the scripts to the \$PATH variable, which is set from within the script, but not written down to .bash_profile to avoid any interference with user settings.

2. The auxiliary scripts are called sequentially by the main script run_get_phylomarkers_pipeline.sh 
according to predefined runmodes and on DNA or protein sequences. However, the auxiliary scripts all have
their own user interface and may be useful to perform specific computations without having to run the pipeline.
All scripts display usage instructions and describe their aims.

### Bash scripts

* run_get_phylomarkers_pipeline.sh (the main script to run the pipeline)
* run_parallel_molecClock_test_with_paup.sh

### Perl scripts
* run_parallel_cmmds.pl
* add_labels2tree.pl 
* add_nos2fasta_header.pl 
* concat_alignments.pl 
* convert_aln_format_batch_bp.pl
* pal2nal.pl 
* popGen_summStats.pl
* rename

#### Perl modules
  File::Rename
 From the BioPerl suite: 
  Bio::AlignIO;
  Bio::PopGen::IO;
  Bio::PopGen::Utilities;
  Bio::PopGen::Statistics;
  Bio::SeqIO;

These are bundled with the software and should work out-of-the-box. Should this fail, 
they can be easily installed in Ubuntu with: 'sudo apt-get install libbio-perl-perl'
For alternative installation options see [bioperl.org INSTALL](http://bioperl.org/INSTALL.html)

### R scripts
* compute_suppValStas_and_RF-dist.R
* run_kdetrees.R consense 

#### R packages
The dependencies can be easily installed in local folder lib/R by calling script *./install_R_deps.R* .
Instead, if you wish these packages to be installed on a system-wide basis, then you should call R with 
superuser privileges and within it execute the following command: 

install.packages( c("ape", "kdetrees", "stingr", "vioplot", "ggplot2", "gplots", "plyr", "seqinr"), dep=T)

Please see examples in the source code of *install_R_deps.R* to solve problems that might arise when
old versions of the, particularly *Rcpp*, are already in the system. Tips are also provided to install
"ape" in Mac systems.

## External dependencies: second party binaries to be installed by the user. 

NOTES: 

1. the corresponding binaries, after installation, shoud be found in on one of the directories listed in the $PATH variable.  On Linux machines this is typically /usr/local/bin (if you have superuser privileges) or \$HOME/bin if not. 

2. If the required binaries are not found in \$PATH, the main script will automatically try to use the statically compiled versions packaged in the distribution under bin/linux or bin/macosx-intel directories, as required for the local environment. 

3. It is strongly recommended, howerver, that the user downloads the latest versions of the binaries from the links provided below. This is particularly important for [paup*](https://people.sc.fsu.edu/~dswofford/paup_test/), since the current version is a test version that automatically expires every 6 months:

* [clustal omega](http://www.clustal.org/omega/). Multiple sequence alignment software. [Sievers et al. 2011](http://msb.embopress.org/content/7/1/539.long). On Ubuntu try: 'sudo apt-get install clustalo'
* [parallel](https://www.gnu.org/software/parallel/). Executes processes in parallel on multicore machines. On Ubuntu try: 'sudo apt-get install parallel'
<!--* [pexec](https://www.gnu.org/software/pexec/). Execute processes in parallel on multicore machines. On Ubuntu try: 'sudo apt-get install pexec'-->
* [Phi test](https://www.maths.otago.ac.nz/~dbryant/software/PhiPack.tar.gz). Recombination test software. [Bruen et al. 2006](http://www.genetics.org/content/172/4/2665.long)
* [FastTree](http://microbesonline.org/fasttree/). Fast maximum-likelihood tree searching program. [Price et al. 2010](http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0009490). On Ubuntu try: 'sudo apt-get install fasttree'
* [paup*](https://people.sc.fsu.edu/~dswofford/paup_test/). Multipurpose phylogenetics software package developed by David Swofford and colleagues. NOTE: This is a test version that expires every 6 months! So please update regularly.

4. Source code and manual compilation instructions are also provided in case bundled binaries fail.
