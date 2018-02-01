# Installation and execution notes for the get_phylomarkers pipeline

Version: Jan. 31, 2018
 
This file lists the software components of the *get_phylomarkers* pipeline and briefly describes how to install them.

The pipeline runs on Linux (Ubuntu and RedHat distros) and Mac OS X environments.

It assumes that recent versions of Bash, Perl and R are installed on a multicore (64bit) machine, ideally a server running Linux.
The pipeline is designed to take advantage of modern multicore machines to parallelize all repetitive tasks that have to be performed on each entry sequence, like tagging sequences, generating multiple sequence alignments, deriving codon alignments, runnig the Phi test on them and inferring maximum likelihood phylogenies from each alignment. We recommend running the pipeline on a multiprocessor/multicore server to speed up these  computations.


## Quick install and test notes

- We highly recommend that you download the [**Docker image**](https://hub.docker.com/r/csicunam/get_homologues), which bundles GET_PHYLOMARKERS with [**GET_HOMOLOGUES**](https://github.com/eead-csic-compbio/get_homologues), ready to use. This is the easiest way to get the full pipeline up and running, avoiding potential architecture-specific configuration and installation problems of the diverse dependencies (Perl modules, R packages, binaries ...). 

- To get a quick impression of the capabilities of the [**GET_HOMOLOGUES**](https://github.com/eead-csic-compbio/get_homologues) + **GET_PHYLOMARKERS** combo, We recommend following the [**tutorial**](docs/GET_PHYLOMARKERS_manual.md#get_phylomarkers-tutorial) with the test sequences provided in the test_sequences/ directory (and subdirectories contained therein). Read the [**manual**](docs/GET_PHYLOMARKERS_manual.md) for the implementation and advanced usage details.

Alternatively, you can try to perform a manual install, as follows:

1. Download the [latest release](https://github.com/vinuesa/get_phylomarkers/releases) or clone the repository into a suitable directory (e.g. $HOME/src/gitHub/). To clone the repo, issue the command 'git clone https://github.com/vinuesa/get_phylomarkers.git' from within $HOME/src/gitHub/

2. Make sure a recent version of R is configured in your system. If not, please install R package (r-base in linux). See CRAN packages for OSX [here](https://cran.r-project.org/bin/macosx/). Your system will also need recent versions of Bash and Perl installed.

3. cd into get_phylomarkers/ and run './install_R_deps.R', which will install R packages into get_phylomarkers/lib/R

4. Copy the test_sequences directory into a suitable place (e.g. 'cp -r test_sequences $HOME)

5. cd into the test_sequences dir (e.g. 'cd $HOME/test_sequences/core_genome')

6. Issue the following command to test if the distro is working on your system: '/path/to/get_phylomarkers/run_get_phylomarkers_pipeline.sh -R 1 -t DNA', which will run in phylogenomics mode (-R 1), on DNA sequences (-t DNA). 
 
7. Check it now on the protein level: 'run_get_phylomarkers_pipeline.sh -R 1 -t PROT'. Note that for this second invocation, you will probably not need to prepend the full path to the script anymore, as symlinks are created to the scripts from your $HOME/bin dir (see NOTES below).

8. Explore the help menu of the main script to see the options available for customization of the run. It is printed to STDOUT when issuing run_get_phylomarkers_pipeline.sh -h or simply run_get_phylomarkers_pipeline.sh

9.  Read the [**manual**](docs/GET_PHYLOMARKERS_manual.md) for the implementation details.

That's it, enjoy. 


## Scripts distributed through GitHub

NOTES: 

1. The main script run_get_phylomarkers_pipeline.sh will automatically identify the directory
on the local machine where the get_phylomarkers package was downloaded on its second invocation. 
It will also check if the host machine has a \$HOME/bin dir included in \$PATH. If so, the main script will automatically generate 
symlinks in \$HOME/bin to the package scripts, so that they become visible system-wide. **We highly encourage users to generate the \$HOME/bin directory, if not available**.
Otherwise, it will prepend the distribution directory holding the scripts to the \$PATH variable, which is set from within the script, but not written down to .bash_profile to avoid any interference with user settings.

2. The auxiliary scripts are called sequentially by the main script run_get_phylomarkers_pipeline.sh 
according to predefined runmodes and on DNA or protein sequences. However, the auxiliary scripts all have
their own user interface and may be useful to perform specific computations without having to (re-)run the full pipeline.
All scripts display usage instructions and describe their aims.

### Bash scripts

* run_get_phylomarkers_pipeline.sh (the main script to run the pipeline)
* run_parallel_molecClock_test_with_paup.sh
* estimate_pangenome_phylogenies.sh

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

install.packages( c("ape", "kdetrees", "stingr", "vioplot", "ggplot2", "gplots", "dplyr", "seqinr"), dep=T)

Please see examples in the source code of *install_R_deps.R* to solve problems that might arise when
old versions of the, particularly *Rcpp*, are already in the system. Tips are also provided to install
"ape" in Mac systems.

## External dependencies: second party binaries called by GET_PHYLOMARKERS scripts. 

NOTES: 

1. The required binaries are provided as part of the distribution. The main script will determine the location of the statically compiled versions packaged in the distribution under bin/linux or bin/macosx-intel directories, as required for the local environment. The corresponding \$bindir is prepended to \$PATH and exported only for the duration of the run of the script to avoid polluting the user's ENVIRONMENT. This also ensures that the pipeline always gets access to tested versions of the binaries. Older, pre-installed versions available on the system may not work properly.

2. The required second-party binaries required by GET_PHYLOMARKERS are all freely available for download from the links provided below. However, the *run_get_phylomarkers_pipeline.sh* script will call the binaries provided with the distribution in the bin/$OSTYPE/ directory, which avoids problems with old versions that may be installed on the user's system. If you

* [clustal omega](http://www.clustal.org/omega/). Multiple sequence alignment software. [Sievers et al. 2011](http://msb.embopress.org/content/7/1/539.long). On Ubuntu try: 'sudo apt-get install clustalo'
* [parallel](https://www.gnu.org/software/parallel/). Executes processes in parallel on multicore machines. On Ubuntu try: 'sudo apt-get install parallel'
* [Phi test](https://www.maths.otago.ac.nz/~dbryant/software/PhiPack.tar.gz). Recombination test software. [Bruen et al. 2006](http://www.genetics.org/content/172/4/2665.long)
* [FastTree](http://microbesonline.org/fasttree/): Fast maximum-likelihood tree searching program. [Price et al. 2010](http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0009490). On Ubuntu try: 'sudo apt-get install fasttree'. Note that generally more recent versions are available at the [FastTree](http://microbesonline.org/fasttree/) distribution page. In order to achieve the highest likelihood scores possible, the binary should be compiled with the double precission flag as shown below:

```
   gcc -DUSE_DOUBLE -O3 -finline-functions -funroll-loops -Wall -o FastTree FastTree.c -lm

```

* [ModelFinder](http://www.iqtree.org/ModelFinder/): Fast model selection for accurate phylogenetic estimates. [(Kalyaanamoorthy et al. 2017)](https://www.nature.com/articles/nmeth.4285)
* [IQ-TREE](http://www.iqtree.org/). Highly accurate maximum-likelihood tree searching program. [Nguyen et. al (2015)](https://academic.oup.com/mbe/article/32/1/268/2925592). You will need to install the latest version 1.6.\*, not the old 1.5.\* version installed on ubuntu with the 'sudo apt install iqtree' option.
* [paup*](https://people.sc.fsu.edu/~dswofford/paup_test/). Multipurpose phylogenetics software package developed by David Swofford and colleagues. NOTE: This is a test version that expires every 6 months! So please update regularly. Shipped version will expire on April 1st, 2018.
* pars, seqboot and consense from Joe Felsenstein's [PHYLIP](http://evolution.genetics.washington.edu/phylip.html) package.
* nw_reroot and nw_support from the [Newick utilities](http://bioinformatics.oxfordjournals.org/cgi/content/abstract/btq243v1) package.
* you may also need to install **bc**,an arbitrary-precision language for performing math calculations with Bash and other shells

3. Source code and manual compilation instructions are also provided in the corresponding \$bindir, in case bundled binaries fail.
