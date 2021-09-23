# Installation and execution notes for GET_PHYLOMARKERS

This file lists the software components of the **GET_PHYLOMARKERS** software package and its dependencies, briefly describing how to install them. It also provides instructions on **how to setup your Docker environment** to run local instances (containers) of the [**GET_PHYLOMARKERS Docker image**](https://hub.docker.com/r/vinuesa/get_phylomarkers) or [**GET_HOMOLOGUES + GET_PHYLOMARKERS Docker image**](https://hub.docker.com/r/csicunam/get_homologues/), both available from the [Docker registry](https://hub.docker.com/). We highly recommend running the containerized, ready-to-use version, which will free you from struggling with system-specific installation and configuration issues.  

The pipeline has been developed and extensively tested on Linux (Ubuntu and RedHat distros). It should also run on Mac OS X machines, but it has been less tested in that environment.

A standard local install pulling from GitHub assumes that recent versions of Bash, Perl and R are installed on a multicore (64bit) machine, ideally a server running Linux.
**The pipeline is designed to take advantage of modern multicore machines to parallelize all repetitive tasks** that have to be performed on each cluster of homologous sequences, like tagging sequences, generating multiple sequence alignments, deriving codon alignments, runnig the Phi test on them and inferring maximum likelihood phylogenies from each alignment. For large genomic datasets, the pipeline should be run on a multiprocessor/multicore server to speed up these computations.


## Quick install and test notes
We provide two ways of installing and running GET_PHYLOMARKERS and GET_HOMOLOGUES on your machine - Through Docker containers or using a standard local install.

### Pulling and running the Docker image
We highly recommend that you download the [**Docker image**](https://hub.docker.com/r/csicunam/get_homologues), which bundles **GET_PHYLOMARKERS** with [**GET_HOMOLOGUES**](https://github.com/eead-csic-compbio/get_homologues), both ready to use. This is probably the easiest way to get the full pipeline up and running with minimal fuzz, avoiding potential architecture-specific configuration and installation problems of the diverse dependencies (Perl modules, R packages, binaries ...). For your convenience, we've also packaged a [**Docker image for GET_PHLYLOMARKERS**](https://hub.docker.com/r/vinuesa/get_phylomarkers) alone. Detailed instructions are [provided below](#installing-and-running-docker-on-your-machine).

### Cloning the repository and manual installation
Alternatively, you can try to perform a manual install, as follows:

1. **GET_PHYLOMARKERS** requires recent versions of Bash and Perl installed on your system. In addition, make sure a recent version of R is installed in your system. If this should not be the case, please visit [The Comprehensive R Archive Network - CRAN](https://cran.r-project.org/) and follow the detailed instructions provided there to install R on different host architectures and operating systems.  

2. Download the [latest release](https://github.com/vinuesa/get_phylomarkers/releases) or clone the repository into a suitable directory (e.g. $HOME/GitHub/). To clone the repo, issue the command 'git clone https://github.com/vinuesa/get_phylomarkers.git' from within $HOME/GitHub/

3. cd into get_phylomarkers/ and run './install_R_deps.R', which will install R packages into get_phylomarkers/lib/R

4. from within the get_phylomarkers distribution directory, as regular user type:
```
rlibs=`for p in $(R -q -e 'print(.libPaths())'); do if [[ "$p" =~ '/' ]]; then echo -n "$p:"; fi; done; echo -n "$wkd"/"$distrodir/lib/R"` && echo "export R_LIBS_SITE=$rlibs >> $HOME/.Rprofile"
```
  - or as sudo append ':/PATH/TO/get_phylomarkers/lib/R' to <code>R_LIBS_SITE=${R_LIBS_SITE-'/usr/local/lib/R/site-library:/usr/lib/R/site-library:/usr/lib/R/library'}</code> in /etc/R/Profile for system-wide permanent changes for all users and session types.
   
5. To setup the *libnw* library required by *estimate_pangenome_phylogenies.sh* in parsimony mode (-R 3), type the following <code>'sudo cp /PATH/TO/get_phylomarkers/lib/libnw.so /usr/local/lib && sudo echo "export LD_LIBRARY_PATH=/usr/local/lib" && sudo ldconfig'</code>

6. cd into the test_sequences/core_genome directory, or copy that directory into a suitable place (e.g. 'cp -r test_sequences $HOME && cd $HOME/test_sequences/core_genome')

7. If you want to perform a system-wide install, you will have to become the superuser (e.g. 'sudo su').

8. Issue the following command from within /path/to/test_sequences/core_genome to test if the distro is working on your system: '/path/to/get_phylomarkers/run_get_phylomarkers_pipeline.sh -R 1 -t DNA', which will run in phylogenomics mode (-R 1), on DNA sequences (-t DNA). 
 
9. Check it now on the protein level: 'run_get_phylomarkers_pipeline.sh -R 1 -t PROT'. Note that for this second invocation, you will probably not need to prepend the full path to the script anymore, as symlinks were created to the scripts from your $HOME/bin dir, or if you run the lines above with root privileges, from /usr/local/bin.

10. Explore the help menu of the master script to see the options available for customizing the runs. It is printed to STDOUT when issuing run_get_phylomarkers_pipeline.sh -h or simply run_get_phylomarkers_pipeline.sh

That's it, enjoy!

## Test the pipeline by running the tutorial execises

To get a quick impression of the capabilities of the **GET_HOMOLOGUES + GET_PHYLOMARKERS** combo, we recommend that you run the [**tutorial**](docs/GET_PHYLOMARKERS_manual.md#get_phylomarkers-tutorial) with the test sequences provided in the test_sequences/ directory (and subdirectories contained therein). This will provide you with a quick grasp of its possibilities.

Read the [**manual**](docs/GET_PHYLOMARKERS_manual.md) for the implementation details and additional functionality not described in the [tutorial](https://github.com/vinuesa/get_phylomarkers/blob/master/docs/GET_PHYLOMARKERS_manual.md#get_phylomarkers-tutorial).


## Installing and running Docker on your machine
Before you can run [docker containers](https://www.docker.com/resources/what-container) on your machine, you will need to install the [Docker engine](https://docs.docker.com/engine/) on your host. You can freely download and install the [Docker Community Edition (CE)](https://docs.docker.com/get-docker/), available for [Linux](https://docs.docker.com/engine/install/), [Mac OS](https://docs.docker.com/desktop/mac/install/) and [Windows](https://docs.docker.com/desktop/windows/install/).

### What is this "Docker" thing you are recommending to install?
- Docker images are a standard way to build/package/ship software applications in a portable format, making it easy for developers to deploy complex software applications as **Docker images** and for users to run on their machines. Images can then be easily pulled (downloaded) by users from an image registry like the official [Docker Hub](https://www.docker.com/products/docker-hub). Once pulled, users run them on their host as an instance of the image known in the jargon as a [**Docker container**](https://docs.docker.com/get-started/part2/).

- After a successful installation of the Docker engine, you will have a Docker daemon running on your machine as well as a Docker client. You will use the latter to issue commands like <code> docker pull; docker run</code> to pull and interact with fecthed **Docker images** like those for[GET_PHYLOMARKERS](https://hub.docker.com/r/vinuesa/get_phylomarkers) or [GET_HOMOLOGUES + GET_PHYLOMARKER image](https://hub.docker.com/r/csicunam/get_homologues/). Have a look at the official [get-started tutorial](https://docs.docker.com/get-started/) if you'd like more details.

#### Recap of basic Docker concepts
- **A container** is launched by running an image with Docker client command like the following: <code>docker run --rm -it <image_name> <command_to_execute_by_started_container></code>. 
- **An image** is an executable package that includes everything needed to run an application–the code, a runtime, libraries, environment variables, and configuration files.

A container is a runtime instance of an image –what the image becomes in memory when executed (that is, an image with state, or a user process). You can see a list of your running containers with the command, <code>docker ps</code>, just as you would in Linux (see below).


### Setting up proper permissions on your host machine to manage Docker as a non-root user
Once Docker is installed on your machine (the host), it is highly recommended to generate a docker group and add your user to it, as explained here: [docker linux-postinstall](https://docs.docker.com/engine/install/linux-postinstall/) and summarized below -

```
# From your Docker host terminal type the following commands:

sudo groupadd docker
sudo usermod -aG docker $USER
```

### Pulling and running the latest version of the [GET_PHYLOMARKERS](https://hub.docker.com/r/vinuesa/get_phylomarkers) and/or [GET_HOMOLOGUES + GET_PHYLOMARKERS](https://hub.docker.com/r/csicunam/get_homologues/) Docker images.
Once Docker is set up on your machine (the host), you can *pull* (fetch) the latest image version available at Dockerhub by typing the following command <code>docker pull <image:version_tag></code> in your terminal: 

```
docker pull vinuesa/get_phylomarkers:latest
```

You can both *pull and run* (execute) an *image* as shown below for the [GET_HOMOLOGUES + GET_PHYLOMARKERS](https://hub.docker.com/r/csicunam/get_homologues/) Docker image

```
# If you are pulling the distro for the first time, the following command will take a while to download everything from Dockerhub ...
docker run --rm -it csicunam/get_homologues:latest /bin/bash

# you will get a bash command prompt within the container that will look something like:
you@ce7444a1ffbd:~$

# make sure GET_HOMOLOGUES is properly setup
get_homologues.pl -v

# do the same for GET_PHYLOMARKERS
run_get_phylomarkers_pipeline.sh -v

# exit this docker container
exit

```

Note: The <code>docker run --rm -it</code> launches an interactive container (*-i* keeps STDIN open even if not attached; *-t* allocates a pseudo-TTY) and removes it after exiting the session (*--rm* Automatically remove the container when it exits). Once the session starts, you will be logged in as **user 'you'** in directory '/home/you/'. The software packages are installed in the root directory, at '/get_homologues' and '/get_phylomarkers'.

#### Update your local images regularly with docker pull
Every once in a time, in order to have the **latest image version** running on your machine, you need to [pull](https://docs.docker.com/engine/reference/commandline/pull/) it again from Dockerhub as shown below

```
# pull the latest version of the image available on Dockerhub
docker pull csicunam/get_homologues:latest

# this will be much faster as the first install, as only modified files are updated.
# if your image corresponds to the last version available, you'll see output like:
latest: Pulling from csicunam/get_homologues
Digest: sha256:e33fce348e28d4544c726961024b8450f48ac7ed1a30c0199247e1d29e0874ae
Status: Image is up to date for csicunam/get_homologues:latest

```

### Final setup - *bind mount* a host directory as a volume into the container for presistent data sharing between the host file system and Docker containers
**IMPORTANT**: By default, Docker containers cannot access data on the host system. This means the container won't be able to access the hosts file system and when stopped, any files written to the thin container's rw storage layer will be lost. Hence, input files and results should be stored in a local, persistent folder, which has to be made visible to the container. The ideal option to achieve this is to [bind mount a host file or directory into the container](https://docs.docker.com/storage/bind-mounts/), using the following general syntax: <code>docker run ... -v /path/in/host:/path/in/container ...</code>. Writes to one will affect the other. Note that both paths have to be absolute paths. Lets perform such a binding mount of a directory containing clusters of homologous sequences (protein and DNA FASTA files) generated by GET_HOMOLOGUES, as detailed in the tutorial. For testing purposes, make a '~/data/genomes/' directory on your host machine.

```
# make a ~/data/genomes/ directory on your host machine
mkdir -p $HOME/data/genomes
```

To *bind-mount* that directory in a get_phylomarker container under '/home/you/data' we would use a command like:

```
# bind mount the host directory ~/data/genomes on the container's /home/you/data dir
docker run --rm -dit -P --name get_phylo -v ~/data/genomes:/home/you/data vinuesa/get_phylomarkers:latest /bin/bash
01417f649c728193bee4e6631b6974037ba41ba46c29d69d9df89103b1fa9266
```

The command launches a detached (*-d*) and interactive (*-it*) *get_phylomarkers:latest* container, bind-mounting the host's ~/data/genomes directory under the container's /home/you/data directory. This allows the container access to the ~/data/genomes directory on the Docker host and sub-directories below it containing results of a previous GET_HOMOLOGUES run to be analyze with GET_PHYLOMARKERS. If you do not have such data ready for GET_PHYLOMARKERS yet, don't worry, the container has a small data set for you to use, as we'll demonstrate next.

To get access to the get_phylo container running in the background, you need to attach your terminal to the container with this command:

```
docker attach <ID> # use the crypto ID printed to STDOUT by the previous command, i.e. 1417f649c728193bee4e6631b6974037ba41ba46c29d69d9df89103b1fa9266

# now your terminal is attached to the Docker container named get_phylo. Listing contents will show the bind-mounted data directory
you@01417f649c72:~$ ls
data
you@01417f649c72:~$ cd data/

# listing directory contents reveals that it is empty, if you just created it.
you@01417f649c72:~/data$ ls
# copy the small test data packaged in the image and available in the running container instance
you@01417f649c72:~/data$ cp -r  cp -r /get_phylomarkers/test_sequences/ .
you@01417f649c72:~/data/$ ls
test_sequences

# cd into the test_sequences directory and list contents
you@01417f649c72:~/data/$ cd test_sequences && ls
core_genome  pan_genome  pIncAC  pIncAC_homologues

# cd into the core_genome directory, which contains pre-computed core-genome sequences for 12 pIncA/C plasmids
you@01417f649c72:~/data/test_sequences$ cd core_genome/

# listing contents reveals pairs of faa and fna files with single-copy orthologous genes and their translations
you@01417f649c72:~/data/test_sequences/core_genome$ ls
1961_hypothetical_protein.faa    1971_ParA.fna                    1999_hypothetical_protein.faa    2014_KfrA_protein.fna
1961_hypothetical_protein.fna    1984_hypothetical_protein.faa    1999_hypothetical_protein.fna    2015_DNA_replication_term...faa
1962_DNA_topoisomerase_II...faa  1984_hypothetical_protein.fna    2003_RepA.faa                    2015_DNA_replication_term...fna
1962_DNA_topoisomerase_II...fna  1989_hypothetical_protein.faa    2003_RepA.fna                    2016_hypothetical_protein.faa
1964_hypothetical_protein.faa    1989_hypothetical_protein.fna    2005_StbA_family_protein.faa     2016_hypothetical_protein.fna
```

Once you access the data in the mounted '~/data/core_genome' volume, you can launch a simple GET_PHYLOMARKERS run as follows:

```
you@01417f649c72:~/data/test_sequences/core_genome$ run_get_phylomarkers_pipeline.sh -R 1 -t DNA
```

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

* run_get_phylomarkers_pipeline.sh (the master script to run the pipeline)
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
* [paup*](https://people.sc.fsu.edu/~dswofford/paup_test/). Multipurpose phylogenetics software package developed by David Swofford and colleagues. NOTE: This is a test version that changes quickly and expires every 6 months! So please update regularly.
* pars, seqboot and consense from Joe Felsenstein's [PHYLIP](http://evolution.genetics.washington.edu/phylip.html) package.
* nw_reroot and nw_support from the [Newick utilities](http://bioinformatics.oxfordjournals.org/cgi/content/abstract/btq243v1) package.
* you may also need to install **bc**,an arbitrary-precision language for performing math calculations with Bash and other shells

3. Source code and manual compilation instructions are also provided in the corresponding \$bindir, in case bundled binaries fail.
