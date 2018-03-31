---
output:
  html_document: default
  pdf_document: default
---
# Installation and execution notes for GET_PHYLOMARKERS

This file lists the software components of the **GET_PHYLOMARKERS** software package and its dependencies, briefly describing how to install them. It also provides instructions on how to prepare your Docker environment to run local instances (containers) of the [**GET_HOMOLOGUES + GET_PHYLOMARKERS Docker image**](https://hub.docker.com/r/csicunam/get_homologues/) available on Docker hub. We highly recommend installing the platform-independent, ready-to-use Docker version, which will free you from struggling with system-specific configuration issues.  

The pipeline has been developed and extensively tested on Linux (Ubuntu and RedHat distros). It should also run on Mac OS X machines, but it has been less tested in that environment.

It assumes that recent versions of Bash, Perl and R are installed on a multicore (64bit) machine, ideally a server running Linux.
**The pipeline is designed to take advantage of modern multicore machines to parallelize all repetitive tasks** that have to be performed on each cluster of homologous sequences, like tagging sequences, generating multiple sequence alignments, deriving codon alignments, runnig the Phi test on them and inferring maximum likelihood phylogenies from each alignment. For large genomic datasets, the pipeline should be run on a multiprocessor/multicore server to speed up these computations.


## Quick install and test notes

### Installing the Docker image
We highly recommend that you download the [**Docker image**](https://hub.docker.com/r/csicunam/get_homologues), which bundles **GET_PHYLOMARKERS** with [**GET_HOMOLOGUES**](https://github.com/eead-csic-compbio/get_homologues), ready to use. This is probably the easiest way to get the full pipeline up and running with minimal fuzz, avoiding potential architecture-specific configuration and installation problems of the diverse dependencies (Perl modules, R packages, binaries ...). Detailed instructions are provided below.

### Cloning the repository
Alternatively, you can try to perform a manual install, as follows:

1. **GET_PHYLOMARKERS** requires recent versions of Bash and Perl installed on your system. In addition, make sure a recent version of R is configured in your system. If not, please install base R (r-base in linux) and recommended packages. See CRAN packages for MACOSX [here](https://cran.r-project.org/bin/macosx/).

2. Download the [latest release](https://github.com/vinuesa/get_phylomarkers/releases) or clone the repository into a suitable directory (e.g. $HOME/GitHub/). To clone the repo, issue the command 'git clone https://github.com/vinuesa/get_phylomarkers.git' from within $HOME/GitHub/

3. cd into get_phylomarkers/ and run './install_R_deps.R', which will install R packages into get_phylomarkers/lib/R

4. cd into the test_sequences/core_genome directory, or copy that directory into a suitable place (e.g. 'cp -r test_sequences $HOME && cd $HOME/test_sequences/core_genome')

5. If you want to perform a system-wide install, you will have to become the superuser (e.g. 'sudo su').

6. Issue the following command from within /path/to/test_sequences/core_genome to test if the distro is working on your system: '/path/to/get_phylomarkers/run_get_phylomarkers_pipeline.sh -R 1 -t DNA', which will run in phylogenomics mode (-R 1), on DNA sequences (-t DNA). 
 
7. Check it now on the protein level: 'run_get_phylomarkers_pipeline.sh -R 1 -t PROT'. Note that for this second invocation, you will probably not need to prepend the full path to the script anymore, as symlinks were created to the scripts from your $HOME/bin dir, or if you run the lines above with root privileges, from /usr/local/bin.

8. Explore the help menu of the master script to see the options available for customizing the runs. It is printed to STDOUT when issuing run_get_phylomarkers_pipeline.sh -h or simply run_get_phylomarkers_pipeline.sh

That's it, enjoy!

## Test the pipeline by running the tutorial execises

To get a quick impression of the capabilities of the **GET_HOMOLOGUES + GET_PHYLOMARKERS** combo, we recommend that you run the [**tutorial**](docs/GET_PHYLOMARKERS_manual.md#get_phylomarkers-tutorial) with the test sequences provided in the test_sequences/ directory (and subdirectories contained therein). This will provide you with a quick grasp of its possibilities.

Read the [**manual**](docs/GET_PHYLOMARKERS_manual.md) for the implementation details and additional functionality not described in the [tutorial](https://github.com/vinuesa/get_phylomarkers/blob/master/docs/GET_PHYLOMARKERS_manual.md#get_phylomarkers-tutorial).


## Installing the Docker environment on your machine
Here we provide brief notes and provide the relevant links to the official and up-to-date instructions to install the freely available [Docker Community Edition (CE)](https://docs.docker.com/install/).

### What is this "Docker" thing you are recommending to install?
- Docker is a standard way to build/package/ship software applications in a portable format, making it easy for developers to deploy complex software applications developed in a Linux environment as **Docker images**. These can then be very easily downloaded by users, who run them as an image instance known in the jargon as a [**Docker container**](https://docs.docker.com/get-started/part2/) on essentially any modern computer running recent OS versions of Linux, MacOS or Windows. 
- You will need to download and install the matching version for your platform [supported platforms](https://docs.docker.com/install/#supported-platforms). 
- After a successful install, you will have a Docker daemon running on your machine and you will be able to interact with **Docker images** you download from [**Dockerhub**](https://hub.docker.com/) like the [GET_HOMOLOGUES + GET_PHYLOMARKER image](https://hub.docker.com/r/csicunam/get_homologues/) (or built yourself on your system) by issuing commands from the **docker client**.
- Read the official Docker [Get Started, Part 1: Orientation and setup](https://docs.docker.com/get-started/) documentation if you'd like further background information.
- There are thousands of images available for you to download. You can search them at [Docker store](https://store.docker.com/)

#### Recap of basic Docker concepts
**A container** is launched by running an image. **An image** is an executable package that includes everything needed to run an application–the code, a runtime, libraries, environment variables, and configuration files.

A container is a runtime instance of an image –what the image becomes in memory when executed (that is, an image with state, or a user process). You can see a list of your running containers with the command, docker ps, just as you would in Linux (see below).

#### Installing Docker on a Linux box
To install Docker CE, you need the 64-bit version of one of these Ubuntu versions:

    Artful 17.10 (Docker CE 17.11 Edge and higher only)
    Zesty 17.04
    Xenial 16.04 (LTS)
    Trusty 14.04 (LTS)

You can install either from:

- the [Docker's repository](https://docs.docker.com/install/linux/docker-ce/ubuntu/#install-docker-ce) 
- or from a [.deb package](https://docs.docker.com/install/linux/docker-ce/ubuntu/#install-from-a-package)
    
#### Install Docker for Mac
The [Docker for Mac install package](https://docs.docker.com/docker-for-mac/install/) includes everything you need to run Docker CE on a Mac.

#### Install Docker for Windows
The [Docker for Windows install package](https://docs.docker.com/docker-for-windows/install/) includes everything you need to run Docker on a Windows system. 

### Fetching and running the latest version of the [GET_HOMOLOGUES + GET_PHYLOMARKERS Docker image](https://hub.docker.com/r/csicunam/get_homologues/) 
Once Docker is set up on your machine (the host), you can run get (fetch) the latest image version available at Dockerhub by typing the following commands in your terminal: 

#### Fetch and run a container instance of the downloaded image

```
# If you are fetching the distro for the first time, the following command will take a while to download everything from Dockerhub ...
docker run --rm -it csicunam/get_homologues:latest /bin/bash

# you will get a container command prompt that will look something like:
you@ce7444a1ffbd:~$

# lets make sure GET_HOMOLOGUES is properly setup
get_homologues.pl -v

# lets do the same for GET_PHYLOMARKERS
run_get_phylomarkers_pipeline.sh -v

# lets exit this docker session
exit

```

#### Update your local image of the GET_HOMOLOGUES + GET_PHYLOMARKERS Docker distro with docker pull
In order to have the **latest version of the image** running on your machine, you need to [pull](https://docs.docker.com/engine/reference/commandline/pull/) it from [Dockerhub](https://hub.docker.com/r/csicunam/get_homologues/)

```
# pull the latest version of the image available on Dockerhub to your machine
docker pull csicunam/get_homologues:latest

# this will be much faster as the first install, as only modified files are updated.
# if your image corresponds to the last version available, you'll see output like:
latest: Pulling from csicunam/get_homologues
Digest: sha256:e33fce348e28d4544c726961024b8450f48ac7ed1a30c0199247e1d29e0874ae
Status: Image is up to date for csicunam/get_homologues:latest

```

#### A few basic docker commands
You can search the full official [reference documentation](https://docs.docker.com/reference/). There are many good basic Docker tutorials on the Internet. You may like [Rominirami's tutorial for beginners](https://rominirani.com/docker-tutorial-series-part-2-basic-commands-baaf70807fd3) to quickly learn more of the basics.

What follows is just a tiny bit to get you quickly up and running.

- Get the installed docker version, search for images on Docker hub and general help + command-specific help

```
docker --version
# Docker version 17.12.0-ce, build c97c6d6

# >>> Display system-wide information
docker info

docker --help # output not shown

docker pull --help
#Usage:	docker pull [OPTIONS] NAME[:TAG|@DIGEST]
#
#Pull an image or a repository from a registry
#
#Options:
#  -a, --all-tags                Download all tagged images in the repository
#      --disable-content-trust   Skip image verification (default true)

# >>> search for an image on Docker hub 
docker search get_phylomarkers

#NAME                             DESCRIPTION                                     STARS               OFFICIAL            AUTOMATED
#eeadcsiccompbio/get_homologues   Ubuntu-based image with GET_HOMOLOGUES and G…   0                                       
#csicunam/get_homologues          Ubuntu-based image with GET_HOMOLOGUES and G…   0

```

- Listing installed images

```
docker images
# or docker image ls

#REPOSITORY                TAG                     IMAGE ID            CREATED             SIZE
#csicunam/get_homologues   01022018-2.1.1_2Feb18   9626fda936c7        14 hours ago        1.37GB
#csicunam/get_homologues   latest                  9626fda936c7        14 hours ago        1.37GB
#ubuntu                    16.04                   0458a4468cbc        9 days ago          112MB
#r-base                    3.4.3                   524f705b5ed1        7 weeks ago         647MB
#hello-world               latest                  f2a91732366c        2 months ago        1.85kB

```
- Listing running containers

```
## List Docker containers (running, all, all in quiet mode)
docker container ls
docker container ls -all
docker container ls -a -q

# example output
docker ps -a
#CONTAINER ID        IMAGE               COMMAND             CREATED             STATUS                    PORTS               NAMES
#c0fe986eb5d8        d61f3b9a7428        "/bin/bash"         41 hours ago        Exited (0) 37 hours ago                       test-container

```

- Removing images and containers

```
# remove and image (think twice befor performing)
rmi IMAGE

# remove one or more containers (good household habit)

```

- How to access the sequence data on the host (your machine) from a container instance
This is an essential operation that is documented in the [**GET_PHYLOMARKERS tutorial**](https://github.com/vinuesa/get_phylomarkers/blob/master/docs/GET_PHYLOMARKERS_manual.md#get_phylomarkers-tutorial).

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
