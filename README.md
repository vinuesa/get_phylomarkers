# GET_PHYLOMARKERS

[![Build Status](https://app.travis-ci.com/vinuesa/get_phylomarkers.svg?branch=master)](https://travis-ci.com/vinuesa/get_phylomarkers)
[![Publication](https://img.shields.io/badge/DOI-10.3389/fmicb.2018.00771-blue)](https://doi.org/10.3389/fmicb.2018.00771)
[![GPLv3-like license](https://img.shields.io/badge/License-GPLv3-blue.svg)](./LICENSE.txt)
[![DockerHub](https://badgen.net/badge/icon/docker?icon=docker&label)](https://hub.docker.com/r/vinuesa/get_phylomarkers)


**GET_PHYLOMARKERS** ([Vinuesa et al., 2018](https://www.frontiersin.org/articles/10.3389/fmicb.2018.00771/full)) is primarily an open-source software package to select optimal markers for microbial phylogenomics and species tree estimation using the multispecies coalescent and concatenation approaches. It implements a [**bioinformatics pipeline**](https://vinuesa.github.io/get_phylomarkers/#brief-presentation-and-graphical-overview-of-the-pipeline) to filter orthologous gene clusters computed by the companion package [**GET_HOMOLOGUES**](https://github.com/eead-csic-compbio/get_homologues) to select those with optimal attributes for phylogenetic inference using maximum likelihood (ML). Multiple sequence alignments of loci passing the filters are concatenated into a supermatrix to estimate a species tree under the ML criterion using the state-of-the-art fast ML tree searching algorithms [FastTree](http://www.microbesonline.org/fasttree) or [IQ-TREE](http://www.iqtree.org). We have also tested it with plant coding sequences (details [here](https://github.com/vinuesa/get_phylomarkers?tab=readme-ov-file#manual-and-tutorials)).

Starting with **release 2.0.0** (2022-11-20), a **concatenation-independent species tree** is computed from the ML gene trees estimated from top-scoring alignments using [ASTRAL-III](https://github.com/smirarab/ASTRAL). 

 **Release 2.1.0** (2024-03-31) implemented the **maximal matched-pairs tests** to asses violations of the data to the Stationarity, Reversibility, and Homogeneity (SRH) assumptions made by maximum-likelihood phylogenetic models, as implemented in [IQ-TREE](http://www.iqtree.org/). In addition, [**ASTRAL-IV**](https://github.com/chaoszhang/ASTER/blob/master/tutorial/astral4.md) is used to estimate the species tree directly from the filtered ML source gene trees, which computes terminal and internal branch lengths in substitution-per-site units.

 **GET_PHYLOMARKERS** can also estimate **ML and parsimony trees from the pan-genome matrix**, including unsupervised learning methods to determine the optimal number of clusters from the pan-genome and average genomic distance matrices. A detailed [**manual**](https://vinuesa.github.io/get_phylomarkers/#get_phylomarkers-manual) and step-by-step [**tutorials**](https://vinuesa.github.io/get_phylomarkers/#get_phylomarkers-tutorial) document the software and help the user to get up and running quickly. For convenience, [**html**](https://vinuesa.github.io/get_phylomarkers/) and [**markdown**](https://github.com/vinuesa/get_phylomarkers/blob/master/docs/GET_PHYLOMARKERS_manual.md) versions of the documentation material are available.


## Installation, dependencies and Docker image

For detailed instructions and dependencies please check [**INSTALL.md**](INSTALL.md).

A [**GET_PHYLOMARKERS Docker image**](https://hub.docker.com/r/vinuesa/get_phylomarkers) is available, as well as an [**image bundling GET_PHYLOMARKERS + 
 GET_HOMOLOGUES**](https://github.com/eead-csic-compbio/get_homologues), ready to use. Detailed instructions for setting up the Docker environment are provided in [**INSTALL.md**](INSTALL.md). How to run container instances with the test sequences distributed with **GET_PHYLOMARKERS** is described in the [**tutorial**](https://vinuesa.github.io/get_phylomarkers/#get_phylomarkers-tutorial).

## Aim
**GET_PHYLOMARKERS** ([Vinuesa et al. 2018](https://www.frontiersin.org/articles/10.3389/fmicb.2018.00771/full)) implements a series of sequential filters (detailed below) to selects markers from the homologous gene clusters produced by [GET_HOMOLOGUES](https://github.com/eead-csic-compbio/get_homologues) with optimal attributes for phylogenomic inference. It estimates **gene-trees** and **species-trees** under the **maximum likelihood (ML) optimality criterion** using state-of-the-art fast ML tree searching algorithms. The species tree is estimated from the supermatrix of concatenated, top-scoring alignments that passed the quality filters outlined in the figures below and explained in detail in the [**manual**](https://vinuesa.github.io/get_phylomarkers/#get_phylomarkers-manual) and [publication](https://www.frontiersin.org/articles/10.3389/fmicb.2018.00771/full).

![**flowchart**](./pics/fmicb-09-00771-g001.jpg) ![**Filtering actions.** GET_PHYLOMARKERS](./pics/fmicb-09-00771-g003.jpg) 

<b>Figure 1A.</b> Simplified flow-chart of the GET_PHYLOMARKERS pipeline showing only those parts used and described in this work. The left branch, starting at the top of the diagram, is fully under control of the master script run_get_phylomarkes_pipeline.sh. The names of the worker scripts called by the master program are indicated on the relevant points along the flow, as detailed in the [**manual**](https://vinuesa.github.io/get_phylomarkers/#get_phylomarkers-manual). The image corresponds to [**Fig. 1 of Vinuesa et al. 2018**](https://www.frontiersin.org/files/Articles/351767/fmicb-09-00771-HTML/image_m/fmicb-09-00771-g001.jpg).

<b>Figure 1B.</b> Combined filtering actions performed by GET_HOMOLOGUES and GET_PHYLOMARKERS to select top-ranking phylogenetic markers to be concatenated for phylogenomic analyses, and benchmark results of the performance of the FastTree (FT) and IQ-TREE (IQT) maximum-likelihood (ML) phylogeny inference programs. The image corresponds to [**Fig. 3 of Vinuesa et al. 2018**](https://www.frontiersin.org/files/Articles/351767/fmicb-09-00771-HTML/image_m/fmicb-09-00771-g003.jpg).


**GET_HOMOLOGUES** is a genome-analysis software package for microbial pan-genomics and comparative genomics originally described in the following publications: 

- [Contreras-Moreira and Vinuesa, AEM 2013](https://www.ncbi.nlm.nih.gov/pubmed/24096415)
- [Vinuesa and Contreras-Moreira, Meth. Mol. Biol. 2015](https://www.ncbi.nlm.nih.gov/pubmed/25343868) 

More recently we developed [GET_HOMOLOGUES-EST](https://github.com/eead-csic-compbio/get_homologues), 
which can be used to cluster eukaryotic genes and transcripts, as described in [Contreras-Moreira et al, Front. Plant Sci. 2017](http://journal.frontiersin.org/article/10.3389/fpls.2017.00184/full). 

If GET_HOMOLOGUES_EST is fed both .fna and .faa files of CDS sequences it will produce **identical output to that of GET_HOMOLOGUES and thus can be analyzed with GET_PHYLOMARKERS all the same**.

* * *

**GET_PHYLOMARKERS** is primarily tailored towards selecting CDSs (gene markers) to infer DNA-level phylogenies of different species of the same genus or family. It can also select optimal markers for population genetics, when the source genomes belong to the same species ([Vinuesa et al. 2018](https://www.frontiersin.org/articles/10.3389/fmicb.2018.00771/full)).
For more divergent genome sequences, classified in different genera, families, orders or higher taxa,
the pipeline should be run using protein instead of DNA sequences.

![ML core-genome phylogeny of Stenotrophomonas](./pics/fmicb-09-00771-g005.jpg) ![ML core-genome phylogeny of Stenotrophomonas](./pics/fmicb-09-00771-g006.jpg)

<b>Figure 2A</b>. Best maximum-likelihood **core-genome phylogeny** for the genus <i>Stenotrophomonas</i> found in the IQ-TREE search, based on the supermatrix obtained by concatenation of 55 top-ranking alignments. The image corresponds to [**Fig. 5 of Vinuesa et al. 2018**](https://www.frontiersin.org/files/Articles/351767/fmicb-09-00771-HTML/image_m/fmicb-09-00771-g005.jpg).

<b>Figure 2B</b>. Maximum-likelihood **pan-genome phylogeny** estimated with IQ-TREE from the consensus pan-genome clusters displayed in the Venn diagram. Clades of lineages belonging to the *S. maltophilia* complex are collapsed and are labeled as in Figure 2A. Numbers on the internal nodes represent the approximate Bayesian posterior probability/UFBoot2 bipartition support values (see methods). The tabular inset shows the results of fitting either the binary (GTR2) or morphological (MK) models implemented in IQ-TREE, indicating that the former has an overwhelmingly better fit. The scale bar represents the number of expected substitutions per site under the binary GTR2+F0+R4 substitution model.  The image corresponds to [**Fig. 6 of Vinuesa et al. 2018**](https://www.frontiersin.org/files/Articles/351767/fmicb-09-00771-HTML/image_m/fmicb-09-00771-g006.jpg).

* * * 

## Manual and tutorials

Please, follow the links for a detailed [**manual**](https://vinuesa.github.io/get_phylomarkers/#get_phylomarkers-manual) and [**tutorials**](https://vinuesa.github.io/get_phylomarkers/#get_phylomarkers-tutorial), including a [**graphical flowchart**](https://vinuesa.github.io/get_phylomarkers/#brief-presentation-and-graphical-overview-of-the-pipeline) of the pipeline and explanations of the implementation details. 
See also our [plant tutorial](http://eead-csic-compbio.github.io/get_homologues/plant_pangenome/protocol.html#downstream-phylogenomic-analyses).

## Citation.

Pablo Vinuesa, Luz-Edith Ochoa-Sanchez and Bruno Contreras-Moreira (2018).
GET_PHYLOMARKERS, a software package to select optimal orthologous clusters for phylogenomics 
and inferring pan-genome phylogenies, used for a critical geno-taxonomic revision of the 
genus *Stenotrophomonas*. [Front. Microbiol. | doi: 10.3389/fmicb.2018.00771](https://doi.org/10.3389/fmicb.2018.00771) 

Published in the Research Topic on "Microbial Taxonomy, Phylogeny and Biodiversity"
http://journal.frontiersin.org/researchtopic/5493/microbial-taxonomy-phylogeny-and-biodiversity

## Code

- Source sode is freely available from [GitHub](https://github.com/vinuesa/get_phylomarkers) and released under a [GPLv3-like license](./LICENSE.txt).
- Docker images ready to pull
    - [GET_PHYLOMARKERS Docker image](https://hub.docker.com/repository/docker/vinuesa/get_phylomarkers)
    - [GET_HOMOLOGUES+GET_PHYLOMARKERS Docker image](https://hub.docker.com/r/csicunam/get_homologues)

## Developers

The code is developed and maintained by [Pablo Vinuesa](https://www.ccg.unam.mx/~vinuesa/) 
at [CCG-UNAM, Mexico](https://www.ccg.unam.mx) and 
Bruno Contreras-Moreira  at [EEAD-CSIC, Spain](https://www.eead.csic.es/compbio). 

## Acknowledgements

We thank Alfredo J. Hernández and Víctor del Moral at CCG-UNAM for technical support with server administration.

### Funding
We gratefully acknowledge the funding provided over the years by [DGAPA-PAPIIT/UNAM](https://dgapa.unam.mx/index.php/impulso-a-la-investigacion/papiit) (grants IN201806-2, IN211814, IN206318, and IN216424) and [CONAHCyT-Mexico](https://conahcyt.mx/) (grants P1-60071, 179133 and A1-S-11242) to [Pablo Vinuesa](https://www.ccg.unam.mx/~vinuesa/), as well as the Fundación ARAID,Consejo  Superior  de Investigaciones Científicas (grant 200720I038 and Spanish MINECO (AGL2013-48756-R) to Bruno Contreras-Moreira.
