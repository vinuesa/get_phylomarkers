# Install notes for get_phylomarkers
This file lists the software components of the get_phylomarkers pipeline and briefly describes how to install them.

## Scripts distributed through GitHub
### Bash scripts
* run_pexec_cmmds.sh 
* run_parallel_molecClock_test_with_paup.sh

### Perl scripts
* add_labels2tree.pl 
* add_nos2fasta_header.pl 
* concat_alns_local.pl 
* convert_aln_format_batch_bp.pl
* pal2nal.pl 
* popGen_summStats.pl
* rename

### R scripts
* compute_suppValStas_and_RF-dist.R
* run_kdetrees.R consense 

## External dependencies: second party binaries to be installed by the user. 

NOTE: the corresponding binaries, after installation, shoud be found in the list of directories hold in the $PATH variable.


* [clustal omega](http://www.clustal.org/omega/). Multiple sequence alignmet software. [Sievers et al. 2011](http://msb.embopress.org/content/7/1/539.long). On Ubuntu try: 'sudo apt-get install clustalo'
* [pexec](https://www.gnu.org/software/pexec/). Execute processes in parallel on multicore machines. On Ubuntu try: 'sudo apt-get install pexec'
* [Phi test](https://www.maths.otago.ac.nz/~dbryant/software/PhiPack.tar.gz). Recombination test software. [Bruen et al. 2006](http://www.genetics.org/content/172/4/2665.long)
* [FastTree](http://microbesonline.org/fasttree/). Fast maximum-likelihood tree searching program. [Price et al. 2010](http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0009490). On Ubuntu try: 'sudo apt-get install fasttree'
* [paup*](https://people.sc.fsu.edu/~dswofford/paup_test/). Multipurpose phylogenetics software developed by David Swofford and colleagues. 
