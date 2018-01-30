# GET_PHYLOMARKERS MANUAL
<!--knitr::opts_chunk$set(fig.align = "center", out.height = '800px', out.width = '0.8\\textheigh'); # to avoid excessively large figures-->
<!--Version Dec. 29cnd, 2017.-->

## Brief presentation and graphical overview of the pipeline

This manual provides the usage details for **GET_PHYLOMARKERS**, a software package designed to select "well-behaved" phylogenetic markers to estimate a **maximum likelihoood (ML) species tree** from the supermatrix of concatenated, top-scoring alignments. These are identified through a series of sequential filters that operate on orthologous gene/protein clusters computed by **GET_HOMOLOGUES** to exclude:

1. alignments with evidence for **recombinant sequences**
2. sequences that yield "**outlier gene trees**" in the context of the distributions of topologies and tree-lengths expected under the multispecies coalescent
3. **poorly resolved gene trees** 

**Figure 1** provides a graphical overview of the **GET_PHYLOMARKERS** pipeline. The Manual will describe in detail each of these steps along with the options available to the user to control the behaviour fo the pipeline, the stringency of the filters and the tree-search thoroughness.

![Figure 1](pics/getphylo_flowchart_FINAL.png)

<!--
<img src="pics/getphylo_flowchart_FINAL.png" width="100%" height="1200px" style="display: block; margin: auto;" />
-->

<!--<img src="pics/getphylo_flowchart_FINAL.png" alt="GET_PHYLOMARKERS pipeline workflow"
     width="600px" height="1500px"></img>-->

In addition, the script *estimate_pangenome_phylogenies.sh* can search for ML and parsimony **pan-genome phylogenies** using the pan-genome matrix computed by *compare_clusters.pl* from the *GET_HOMOLOGUES* suite.

## Installation, dependencies and Docker image

The GET_HOMOLOGUES pagacke can be downloaded as a [GitHub distribution](https://github.com/vinuesa/get_phylomarkers/releases). For detailed instructions and dependencies please check [INSTALL.md](https://github.com/vinuesa/get_phylomarkers/blob/master/INSTALL.md).

A [Docker image](https://hub.docker.com/r/csicunam/get_homologues) is available with GET_PHYLOMARKERS
bundled with [GET_HOMOLOGUES](https://github.com/eead-csic-compbio/get_homologues), ready to use. We highly recommend installing the docker image to avoid potential problems with the installation of the many second-party dependencies.

## Aim
**GET_PHYLOMARKERS** implements a series of sequential filters (**Fig. 1** and explained below) to selects markers from the homologous gene clusters produced by [GET_HOMOLOGUES](https://github.com/eead-csic-compbio/get_homologues) with optimal attributes for phylogenomic inference. It estimates **gene-trees** and **species-trees** under the **maximum likelihood (ML) optimality criterion** using *state-of-the-art* fast ML tree searching algorithms. The species tree is estimated from the supermatrix of concatenated, top-scoring alignments that passed the quality filters. 

**GET_HOMOLOGUES** is a genome-analysis software package for microbial pan-genomics and comparative genomics originally described in the following publications: 

- [Contreras-Moreira and Vinuesa, AEM 2013](https://www.ncbi.nlm.nih.gov/pubmed/24096415)
- [Vinuesa and Contreras-Moreira, 2015](https://www.ncbi.nlm.nih.gov/pubmed/25343868) 

More recently we have developed [GET_HOMOLOGUES-EST](https://github.com/eead-csic-compbio/get_homologues), 
which can be used to cluster eukaryotic genes and transcripts, as described in [Contreras-Moreira et al., 2017](http://journal.frontiersin.org/article/10.3389/fpls.2017.00184/full).

**GET_PHYLOMARKERS** is primarily tailored towards selecting CDSs (gene markers) to infer DNA-level phylogenies of different species of the same genus or family. It can also select optimal markers for population genetics, when the source genomes belong to the same species.
For more divergent genome sequences, belonging for example to different genera, families, orders or higher taxa,
the pipeline should be run using protein instead of DNA sequences.

## Usage synopsis

1. The pipeline is run by executing the main script *run_get_phylomarkers_pipeline.sh* inside a folder containing twin \*.fna and \*.faa FASTA files for orthologous **single-copy** CDSs and translation products. <!-- computed by the *get_homologues.pl -e* or *compare_clusters.pl* -t number_of_genomes scripts of the **GET_HOMOLOGUES** suite.-->
2. There are two **runmodes**: **-R 1** (for phylogenetics) and **-R 2** (for population genetics).
3. The pipeline can be run on two **molecular types**: **DNA** or **protein** sequences (**-t DNA|PROT**). The latter is intended for the analysis of more divergent genome sequences, typically above the genus level.
4. From version 2.0 onwards, **GET_PHYLOMARKERS** uses either the **FastTree (FT)** or **IQ-TREE (IQT)** fast ML tree search algorithms, controlled with the **-A [F|I]** option [default: I], respectively. Note that previous versions used FastTree as the default. This change is based on our benchmark analyses (Vinuesa et. al, 2018. Submitted), which showed that IQ-TREE systematically finds better-scoring trees which allowed for finer-grained filtering of markers based on their phylogentic attributes. This is in line with a recent publication that concluded that IQT is the best ML tree-searching algorithm available to date for datasets with no more than 100-200 sequences, as demonstrated in a large benchmark study with empirical phylogenomic datasets (Zhou et al. 2017. Mol Biol Evol. Nov 21. doi: 10.1093/molbev/msx302.) [PMID:29177474](https://www.ncbi.nlm.nih.gov/pubmed/29177474).
5. As of version 2.0 (January, 22cnd, 2018), GET_PHYLOMARKERS uses IQ-TREE version 1.6.1 (released Dec. 23rd, 2017) which implements a fast search option, which almost matches the speed of FastTree, but rataining the accuracy of IQ-TREE 1.5.* Model-selection is performed with the **-fast flag** to for maximal speed both for gene- and species tree searches.
6. The **global molecular-clock hypothesis** can be evaluated for DNA (codon) alignments (-R 1 -t DNA -K 1). It is not yet implemented for protein sequences.
7. **GET_PHYLOMARKERS** can compute **basic population genetics descritive statistics and neutrality tests** when run in popGen mode (**-R 2**). 
8. A **toy sequence dataset is provided** with the distribution in the test_sequences/ directory for easy and fast testing (~16-60 seconds on a commodity GNU/Linux desktop machine with 4 cores; see [INSTALL.md](INSTALL.md), depending on the search thoroughness and/or nuber of models to be evaluated). 

### Basic usage examples

```
 # default run: launches the IQT-based filtering pipeline, selecting best models for gene trees and fitting GTR+RATE_AND_F_PARAMS models for species tree search based on the supermatrix.
 run_get_phylomarkers_pipeline.sh -R 1 -t DNA

 # Same as above, but adding molecular-clock analysis assuming a HKY85+G substitution model                
 run_get_phylomarkers_pipeline.sh -R 1 -t DNA -K 1 -M HKY

 # population-genetics mode   
 run_get_phylomarkers_pipeline.sh -R 2 -t DNA                
 
 # protein alignments, user-defined kdetrees & mean branch support cutoff values
 run_get_phylomarkers_pipeline.sh -R 1 -t PROT -k 1.2 -m 0.7 

 # To run the pipeline on a remote server, we recommend using the nohup command upfront, as shown below:
 #   in this case, calling also IQ-TREE, which will select among the (TNe,TVM,TVMe,GTR)+RATE models and do
 #   5 independent tree searches under the best-fit model, computing ultrafast bootstrapp and aproximate Bayes
 #   branch support values 
 nohup run_get_phylomarkers_pipeline.sh -R 1 -t DNA -S 'HKY,TN,TVM,TIM,SYM,GTR' -k 1.0 -m 0.7 -T high -N 5 &> /dev/null &
```

## Usage and design details

1. Start the run from within the directory holding core gene clusters generated by either *get_homologues.pl -e* or 
subsequent intersection (OMCL,COGS,BDBH) clusters produced with *compare_clusters.pl -t number_of_genomes* from the 
[GET_HOMOLOGUES](https://github.com/eead-csic-compbio/get_homologues) package.
   
  NOTE: **both .faa and .fna files are required** to generate codon alignments from DNA fasta files. This
            means that two runs of *compare_clusters.pl* (from the **GET_HOMOLOGUES** package) are required,
	          one of them using the -n flag. 
	    
2. *run_get_phylomarkers_pipeline.sh* is intended to run on a collection of **single-copy** sequence clusters from 
different species or strains.

   NOTES: an absolute minimum of 4 distinct haplotypes per cluster are required
        for the cluster to be evaluated. Clustes with < 4 haplotypes are automatically
	      discarded. This means that at least 4 distinct genomes should be used as input.
	      However, the power of the pipeline for selecting optimal genome loci 
     	  for phylogenomics improves when a larger number of genomes are available 
     	  for analysis. Reasonable numbers lie in the range of 10 to 200 clearly
     	  distinct genomes from multiple species of a genus, family, order or phylum.
     	  The pipeline may not perform satisfactorily with too distant genome sequences,
     	  particularly when sequences with significantly distinct nucleotide or aminoacid
     	  compositions are used. This type of sequence heterogeneity is well known to 
     	  cause systematic bias in phylogenetic inference.

## On the filtering criteria. 

*run_get_phylomarkers_pipeline.sh* uses a **hierarchical filtering scheme**, as follows:

### Detection of recombinant loci. 

Codon or protein alignments (depending on runmode) are first screened with **Phi-test** 
([Bruen et al. 2006](http://www.genetics.org/content/172/4/2665.long)) for the presence of potential recombinant sequences. It is a well established fact that recombinant sequences negatively impact phylogenetic inference when using algorithms that do not account for the effects of this evolutionary force. The permutation test with 1000 permutations is used to compute the *p*-values. These are considered significant if *p* < 0.05. Some loci may not contain sufficient polymorphisms for the test to
work. In that case, the main script assumes that the locus does not contain recombinant sites.
 
### ii) Detection of trees deviating from expectations of the (multispecies) coalescent.

The next filtering step is provided by the **kdetrees-test**, which checks the distribution of topologies, tree lengths and branch lengths. *kdetrees* is a non-parametric method for estimating distributions of phylogenetic trees 
([Weyenberg et al. 2014](https://academic.oup.com/bioinformatics/article-lookup/doi/10.1093/bioinformatics/btu258)), 
with the goal of identifying trees that are significantly different from the rest of the trees in the sample, based on the analysis of topology and branch length distributions. Such "outlier" trees may arise for example from horizontal gene transfers or gene duplication (and subsequent neofunctionalization) followed by differential loss of paralogues among lineages. Such processes will cause the affected genes to exhibit a history distinct from those of the majority of genes, which are expected to be generated by the  (multispecies) coalescent as species or populations diverge. Alignments producing significantly deviating trees in the kdetrees test are identified and saved in the kde_outliers/ directory. The corresponding alignments are not used in downstream analyses.

```      
      * Parameter for controlling kdetrees stingency:
      -k <real> kde stringency (0.7-1.6 are reasonable values; less is more stringent)
     			       [default: 1.5]
```

### Phylogenetic signal content. 

The alignments passing the two previous filters are subjected to **maximum likelihood (ML) tree searches** with *FastTree* or *IQ-TREE* to 
infer the corresponding ML gene trees. Their **phylogenetic signal is computed from the Shimodaria-Hasegawa-like likelihood ratio test branch support values**, which vary between 0-1, as we have reported previously ([Vinuesa et al. 2008](http://aem.asm.org/content/74/22/6987.long)). The support values of each internal branch or bipartition are parsed to compute the mean support value for each tree.
**Alignments/Trees with a mean support value below a cutoff threshold are discarded**.

```
      * Parameters controlling filtering based on mean support values.
      -m <real> min. average support value (0.7-0.8 are reasonable values) 
     		for trees/loci to be selected as informative [default: 0.75]
```

### Evaluating the global molecular clock hypothesis.

*run_get_phylomarkers_pipeline.sh* calls the auxiliary script *run_parallel_molecClock_test_with_paup.sh*
to evaluate the **global molecular clock hypothesis** on the top markers, selected according to the criteria explained in the three previous
points. The script calls [paup*](https://people.sc.fsu.edu/~dswofford/paup_test/) 
to evaluate the free-rates and clock hypothesis using likelihood ratio tests using R. Currently this is only performed on codon alignments. Future versions will implement the global clock hypothesis test also for protein alignments.

### On tree searching: 
#### ModelFinder + IQ-TREE searches:
As of version 1.9.9.0_22Dec17, **GET_PHYLOMERKERS** implements the **-A I** option, which calls [IQ-TREE](http://www.iqtree.org/) for ML tree searching [Nguyen et. al (2015)](https://academic.oup.com/mbe/article/32/1/268/2925592). This is the most recent fast ML software on the scene, which was developed with the aim of escaping from early local maxima encountered during "hill-climbing" by generating multiple seed trees to intiate tree searches, maintaining a pool of candidate trees during the entire run. Overall, it was the best-performing ML tree search algorithm among those evaluated by [Zhou et al. (2017)](https://www.ncbi.nlm.nih.gov/pubmed/29177474) in their above-mentioned benchmarking paper, at least for datasets < 200 taxa. For larger datasets (several hudreds to thousands of sequencs), current implementations of algorithms like RAxML and ExaML, which make hevay use of SPR-moves, will most likey outperform IQ-TREE, which makes more intensive use of local NNI-moves. 

Our benchmark analyses have shown that [IQ-TREE](http://www.iqtree.org/) (v1.6.1 and above) [Nguyen et. al (2015)](https://academic.oup.com/mbe/article/32/1/268/2925592) runs quickly enough when the **-fast** flag is passed to make it feasible to include a model selection step using  [ModelFinder](http://www.iqtree.org/ModelFinder/) [(Kalyaanamoorthy et al. 2017)](https://www.nature.com/articles/nmeth.4285) withouth incurring in prohibitively long computation times. 

Combined with its superior tree-searching algorithm, makes IQT the clear winner in our benchmarks. Therefore, **from version 2.0 onwards, *GET_PHYLOMARKERS* uses IQT as its default tree searching algorithm**, both for the estimation of gene-trees and the species-tree (from the concatenated, top-scoring alignments).
	  
However, the number of models evaluated by ModelFinder (integrated in IQ-TREE) differ for the gene-tree and species-tree search phases, as shown below: 
```
IQT gene tree searches (hard-coded): -T <high|medium|low|lowest> [default: medium]

- high:   -m MFP -nt 1 -alrt 1000 -fast

- medium: -mset K2P,HKY,TN,TNe,TIM,TIMe,TIM2,TIM2e,TIM3,TIM3e,TVM,TVMe,GTR

- low:	   -mset K2P,HKY,TN,TNe,TVM,TVMe,TIM,TIMe,GTR

- lowest: -mset K2P,HKY,TN,TNe,TVM,TIM,GTR

```

All gene trees are run in parallel under the modelset with the following parameters: -mset XXX -m MFP -nt 1 -alrt 1000 -fast

- IQT species-tree searches on the supermatrix: 

	      -S <string> quoted 'comma-separated list' of base models to be evaluated by IQ-TREE
	         when estimating the species tree from the concatenated supermatrix.
	         
If no -S is passed, then sinlge default models are used, as shown below:

```
 - for DNA alignments [default: GTR]

  <'JC,F81,K2P,HKY,TrN,TNe,K3P,K81u,TPM2,TPM2u,TPM3,TPM3u,TIM,TIMe,TIM2,TIM2e,TIM3,TIM3e,TVM,TVMe,SYM,GTR'>  
     
 - for PROT alignments [default: LG]

   <'BLOSUM62,cpREV,Dayhoff,DCMut,FLU,HIVb,HIVw,JTT,JTTDCMut,LG,mtART,mtMAM,mtREV,mtZOA,Poisson,PMB,rtREV,VT,WAG'>          
```
                
In addition, if **-T high**, *run_get_phylomarkers_pipeline.sh* will launch -N <integer> [default: 10] independent IQT searches on the supermatrix of concatenated top-scoring markers.

After selecting the best substitution model, which includes taking care of among-site rate variation, IQ-TREE will search for the best tree, including bootstrapping with the **UFBoot2 ultrafast bootstrapping algorithm** [(Hoang et al. 2017)](https://www.ncbi.nlm.nih.gov/pubmed/29077904) and estimation of approximate Bayes support values. 

The following are example *run_get_phylomarkers_pipeline.sh* invocations to perform IQ-TREE searches. Note that as of version 1.9.10 (January 1st, 2018), the script calls IQ-TREE 1.6.1 with the -fast flag enabled for maximal speed.


```    
# 1. Default IQ-TREEE run (-T medium), evaluating the corresponding set models during gene-tree searches and evaluating the base GTR model  on the concatenated DNA supermatrix, using a single search
run_get_phylomarkers_pipeline.sh -R 1 -t DNA

# 2. Run 10 independent IQ-TREEE runs on a concatenated DNA supermatrix, evaluating the TN,TIM,TVM,GTR base models
run_get_phylomarkers_pipeline.sh -R 1 -t DNA  -S 'TN,TIM,TVM,GTR' -k 0.9 -m 0.8 -T high -N 10 &> /dev/null &

# 3. Run 5 independent IQ-TREEE runs on a concatenated PROT supermatrix, evaluating the LG,WAG,JTT matrices 
run_get_phylomarkers_pipeline.sh -R 1 -t PROT -A I -S 'LG,WAG,JTT,VT' -k 1.0 -m 0.7 -T high -N 5 &> /dev/null &
 
# 4. To run the pipeline on a remote server, we recommend using the nohup command upfront, as shown below:
nohup run_get_phylomarkers_pipeline.sh -R 1 -t DNA -S 'HKY,TN,TVM,TIM,SYM,GTR' -k 1.0 -m 0.7 -T high -N 5 &> /dev/null &	  

```

#### FastTree searches:
- *run_get_phylomarkers_pipeline.sh* performs tree searches using the [FastTree](http://microbesonline.org/fasttree/) program ([Price et al. 2010](http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0009490)) when invoked with the -A F parameter. 
- This program implements many heuristics tailored towards improving the cpu time and memory usage, making it the fastest ML tree searching program that is currently available. However, this comes at a price: in a recent comprehensive benchmark analysis of fast ML phylogenetic programs using large and diverse phylogeneomic datasets, FastTree was shown to be the worst-scoring one in regard to lnL values and topological accuracy [(Zhou et al. 2017)](https://www.ncbi.nlm.nih.gov/pubmed/29177474). 
- Given its outstanding speed, low RAM requirements, acceptance of DNA and protein sequence alignments and on-the-fly computation ofthe above-mentioned Shimodaria-Hasegawa-like likelihood ratio test of branch support, we recommend it for fast exploration of large datasets, but not for final analyses.
- A clear limitation is that it implements only very few substitution models. However, for divergent sequences of different species within a bacterial taxonomic genus or family, our experience has shown that almost invariably the GTR+G model is selected by jmodeltest2, particularly when there is base frequency heterogeneity. The GTR+G+CAT is the substitution model used by *run_get_phylomarkers_pipeline.sh* 
calls of FastTree on codon alignments. 
- To maximize accuracy, it is important to compile FastTree with double precission enabled in order to obtain the highest *lnL* scores possible. This is particularly critical when highly similar sequences are present in the dataset.
- To further enhance accuracy, the gene trees are computed by performing a thorough tree search, as hardcoded in the following FastTree call, which performs a significantly more intense tree search than the default setting used by [Zhou et al. (2017)](https://www.ncbi.nlm.nih.gov/pubmed/29177474). 

```      
     	-nt -gtr -gamma -bionj -slownni -mlacc 3 -spr 8 -sprlength 10 
```    

For concatenated codon alignments, which may take a considerable time (up to several hours) or
for large datasets (~ 100 taxa and > 300 concatenated genes) the user can choose to run FastTree with at different **levels of tree-search thoroughness**: high|medium|low|lowest 
```      
      high:   -nt -gtr -gamma -bionj -slow -slownni -mlacc 3 -spr 8 -sprlength 10
      medium: -nt -gtr -gamma -bionj -slownni -mlacc 2 -spr 4 -sprlength 10 
      low:    -nt -gtr -gamma -bionj -slownni -spr 4 -sprlength 10 
      lowest: -nt -gtr -gamma -mlnni 4
```      
where '-s spr' and '-l spr_length' can be set by the user. The lines above show their default values.
      
For protein alignments, the search parameters are the same, only the model changes to lg

```      
      high: -lg -gamma -bionj -slow -slownni -mlacc 3 -spr 4 -sprlength 10
```

Example invokations for FastTree searches

```
# FastTree searching on a huge protein dataset for fast inspection
run_get_phylomarkers_pipeline.sh -R 1 -A F -t PROT -m 0.6 -k 1.0 -T lowest


# FastTree searching on DNA dataset using the most through possible search and extensive spr rounds
run_get_phylomarkers_pipeline.sh -R 1 -A F -t DNA -m 0.7 -k 1.0 -T high -s 20 -l 12

```

# GET_PHYLOMARKERS TUTORIAL

## Test datasets
The **GET_PHYLOMARKERS** distribution provides *test_sequences* directory which holds the subdirectories core_genome and pan_genome. The first one contains \*.fna and \*.faa FASTA files with the consensus (BDBH, COGtriangles and OMCL) core-genome clusters computed with **GET_HOMOLOGUES** from a set of 12 GenBank-formatted pIncA/C plasmids. The latter holds the pan-genome matrix computed by *compare_clusters.pl* from the **GET_HOMOLOGUES** suite in tabular (\*.tab), FASTA (\*.fasta) and phylip (\*.phy) formats. The pIncAC_gbk directory holds the source \*.gbk GenBank files. This directory has a README.txt file that briefly describes the GenBank files.

The following sections will provide code examples on how to run the full GET_HOMOLOGUES + GET_PHYLOMARKERS pipelines using the test sequences.

## Computing a consensus core-genome with GET_HOMOLOGUES

Go (cd) into the distribution directory and cd into the subidrectory test_sequences and issue the commands shown below.

```
# 1. cd into the directory holding the test_sequences
cd test_sequences/
ls

# 2. run get_homologues to compute the set of homologous clusters using the BDBH, COGtriangles and OMCL clustering algorithms 
get_homologues.pl -d pIncAC -t 12 -e -n 4  # BDBH clusters containing sequences for the 12 genomes, excluding those with inparalogues (-e)
get_homologues.pl -d pIncAC -G -t 0        # COGtriangles, computing clusters of all sizes (-t 0)
get_homologues.pl -d pIncAC/ -M -t 0       # OMCL

```

## Computing a consensus core-genome with GET_HOMOLOGUES

```
# Compute consensus core-genome clusters using compare_clusters.pl of the GET_HOMOLOGUES suite.
# Note that we run the script twice, once with the -n flag (to compute the consesus clusters at the DNA level, *.fna files)
# and a second instance without the flag, to get the protein clusters (*.faa files)

cd pIncAC_homologues/
find . -type d

compare_clusters.pl -d ./KlebsiellapneumoniaeplasmidpNDM-KNNC019153_f0_12taxa_algBDBH_e1_,./KlebsiellapneumoniaeplasmidpNDM-KNNC019153_f0_0taxa_algCOG_e0_,./KlebsiellapneumoniaeplasmidpNDM-KNNC019153_f0_0taxa_algCOG_e0_ -o core_BCM -n

compare_clusters.pl -d ./KlebsiellapneumoniaeplasmidpNDM-KNNC019153_f0_12taxa_algBDBH_e1_,./KlebsiellapneumoniaeplasmidpNDM-KNNC019153_f0_0taxa_algCOG_e0_,./KlebsiellapneumoniaeplasmidpNDM-KNNC019153_f0_0taxa_algCOG_e0_ -o core_BCM

```

## Computing a consensus pan-genome with GET_HOMOLOGUES

```
# Compute consensus pan-genome clusters and matrix using compare_clusters.pl of the GET_HOMOLOGUES suite.
# Note that we run the script twice, once with the -n and -m flags (to compute the consesus clusters at the DNA level, *.fna files
#  and the pan-genome matrix) and a second instance without these flags, to get the protein clusters (*.faa files). 
# Note also that we exclude the directory holding the BDBH clusters, 
#  as these are not suitable to compute a proper pan-genome matrix, since the BDBH clusters always contain the reference sequence.

compare_clusters.pl -d ./KlebsiellapneumoniaeplasmidpNDM-KNNC019153_f0_0taxa_algCOG_e0_,./KlebsiellapneumoniaeplasmidpNDM-KNNC019153_f0_0taxa_algCOG_e0_ -o pan_CM -n -m
compare_clusters.pl -d ./KlebsiellapneumoniaeplasmidpNDM-KNNC019153_f0_0taxa_algCOG_e0_,./KlebsiellapneumoniaeplasmidpNDM-KNNC019153_f0_0taxa_algCOG_e0_ -o pan_CM

```


## Searching for the best core-genome phylogeny of the pIncA/C plasmids with GET_PHYLOMARKERS

```
# assumes that you are within pIncAC_homologues/
cd core_BCM

# 1. make sure we have the same nuber of faa and fna cluster files
find . -name '*.faa' | wc -l
find . -name '*.fna' | wc -l

# 2. run the pipeline under default values
run_get_phylomarkers_pipeline.sh -R 1 -t DNA 

```

## Estimating the ML pan-genome phylogeny of the pIncA/C plasmids with GET_PHYLOMARKERS

```
# assumes that you are within pIncAC_homologues/
cd pan_CM

# 1. find the pangenome_matrix
ls pangenome_matrix*


# 2. run estimate_pangenome_phylogenies.sh launching 10 independent iqtree searches, fitting binary (BIN) models
estimate_pangenome_phylogenies.sh -f pangenome_matrix_t0.fasta -r 10

```

## Estimating the pan-genome phylogeny of the pIncA/C plasmids under the parsimony criterion with GET_PHYLOMARKERS

```
# assumes that you are within pIncAC_homologues/
cd pan_CM

# 1. find the pangenome_matrix
ls pangenome_matrix*


# 2. run estimate_pangenome_phylogenies.sh launching 10 independent iqtree searches, fitting binary (BIN) models

estimate_pangenome_phylogenies.sh -c PARS -R 3 -i pangenome_matrix_t0.phylip -n 4 -b 25 -j 10 -t 1

```

# Developers
The code is developed and maintained by [Pablo Vinuesa](http://www.ccg.unam.mx/~vinuesa/) 
at [CCG-UNAM, Mexico](http://www.ccg.unam.mx/) and 
[Bruno Contreras-Moreira](https://digital.csic.es/cris/rp/rp02661/) 
 at [EEAD-CSIC, Spain](http://www.eead.csic.es/). It is released to the public domain under the GNU GPLv3 [license](./LICENSE).
 


# Citation.

On Jaunary 15th, 2018 we submitted a manuscript describing the implementation get_phylomarkers pipeline and its use in genomic taxonomy and population genomics of *Stenotrophomonas* bacteria to the [Research Topic of Frontiers in Microbiology: Microbial Taxonomy, Phylogeny and Biodiversity](http://journal.frontiersin.org/researchtopic/5493/microbial-taxonomy-phylogeny-and-biodiversity).

Meanwhile, if you find the code useful for your academic work, please use the following citation:
Pablo Vinuesa and Bruno Contreras-Moreira 2018. Get_PhyloMarkers, a pipeline to select optimal markers for microbial phylogenomics, systematics and genomic taxomy. Available at https://github.com/vinuesa/get_phylomarkers and released under the GNU GPLv3 license.

# Acknowledgements

## Personal
We thank Alfredo J. Hernández and Víctor del Moral at CCG-UNAM for technical support.

## Funding
We gratefully acknowledge the funding provided by [DGAPA-PAPIIT/UNAM](http://dgapa.unam.mx/index.php/impulso-a-la-investigacion/papiit) (grants IN201806-2, IN211814 and IN206318) and [CONACyT-Mexico](http://www.conacyt.mx/) (grants P1-60071, 179133 and FC-2015-2-879) to [Pablo Vinuesa](http://www.ccg.unam.mx/~vinuesa/), as well as the Fundación ARAID,Consejo  Superior  de Investigaciones Científicas (grant 200720I038 and Spanish MINECO (AGL2013-48756-R) to [Bruno Contreras-Moreira](https://digital.csic.es/cris/rp/rp02661).
