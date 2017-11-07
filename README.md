# get_phylomarkers

<!--Version June 1st, 2017.-->

This repository hosts the code for the *get_phylomarkers* pipeline. This file describes its
aim and basic usage notes. <!--See [INSTALL.md](INSTALL.md) for installation instructions.-->
The code is developed and maintained by [Pablo Vinuesa](http://www.ccg.unam.mx/~vinuesa/) 
at [CCG-UNAM, Mexico](http://www.ccg.unam.mx/) and 
[Bruno Contreras-Moreira](https://digital.csic.es/cris/rp/rp02661/) 
 at [EEAD-CSIC, Spain](http://www.eead.csic.es/). It is released to the public domain under the GNU GPLv3 [license](./LICENSE).
 
## Installation

For installation instructions please check [INSTALL.md](INSTALL.md).

## Aim
The pipeline selects markers with optimal phylogenetic attributes from the homologous gene 
clusters produced by [GET_HOMOLOGUES](https://github.com/eead-csic-compbio/get_homologues). 
This is a genome-analysis pipeline for microbial pan-genomics and comparative genomics originally 
described in the following publications: 
[Contreras-Moreira and Vinuesa, AEM 2013](https://www.ncbi.nlm.nih.gov/pubmed/24096415) and 
[Vinuesa and Contreras-Moreira, 2015](https://www.ncbi.nlm.nih.gov/pubmed/25343868). More recently
we have developed [GET_HOMOLOGUES-EST](https://github.com/eead-csic-compbio/get_homologues), 
which can be used to cluster eukaryotic genes and transcripts, as described in [Contreras-Moreira et al., 2017](http://journal.frontiersin.org/article/10.3389/fpls.2017.00184/full).

The homologous gene/protein clusters selected by the *get_phylomarkers* pipeline are optimally 
suited for genome phylogenies. The pipeline is primarily tailored towards selecting 
CDSs (gene markers) to infer DNA-level phylogenies of different species of the same genus or family. 
It can also select optimal markers for population genetics, when the source genomes belong to the same species.
For more divergent genome sequences, belonging to different taxonomic genera, families or even orders,
the pipeline can be also run using protein instead of DNA sequences.

## Usage synopsis

1. The pipeline is run by executing the main script *run_get_phylomarkers_pipeline.sh* inside a folder containing twin \*.fna and \*.faa FASTA files for orthologous **single-copy** CDSs and translation products. <!-- computed by the *get_homologues.pl -e* or *compare_clusters.pl* scripts of the **GET_HOMOLOGUES** suite.-->
2. There are two **runmodes**: -R 1 (for phylogenetics) and -R 2 (for population genetics).
3. The pipeline can be run on two **molecular types**: DNA or protein sequences (-t DNA|PROT). The latter is intended for the analysis of more divergent genome sequences, above the genus level.
4. The **global molecular-clock hypothesis** can be evaluated for DNA (codon) alignments. It is not yet implemented for protein sequences.
5. A **toy sequence dataset is provided** with the distribution in the test_sequences/ directory for easy and fast testing (~14 seconds on a commodity GNU/Linux desktop machine with 4 cores; see [INSTALL.md](INSTALL.md)). 

### Basic usage examples

```
 run_get_phylomarkers_pipeline.sh -R 1 -t DNA                # default run
 run_get_phylomarkers_pipeline.sh -R 1 -t DNA -K 1 -M HKY    # add molecular-clock analysis assuming a HKY85+G substitution model
 run_get_phylomarkers_pipeline.sh -R 2 -t DNA                # population-genetics mode
 run_get_phylomarkers_pipeline.sh -R 1 -t PROT -k 1.2 -m 0.7 # protein alignments, user-defined kdetrees & mean branch support cutoff values
```

## Usage and design details

1. Start the run from **within the directory** holding core gene clusters generated by either *get_homologues.pl -e* or 
subsequent intersection (OMCL,COGS,BDBH) clusters produced with *compare_clusters.pl* from the 
[GET_HOMOLOGUES](https://github.com/eead-csic-compbio/get_homologues) package.
   
  NOTE: **both .faa and .fna files are required** to generate codon alignments from DNA fasta files. This
            means that two runs of *compare_clusters.pl* (from the **GET_HOMOLOGUES** package) are required,
	          one of them using the -n flag. 
	    
2. *run_get_phylomarkers_pipeline.sh* is intended to run on a collection of **single-copy** sequence clusters from 
different species or strains.

   NOTES: an absolute minimum of 4 distinct genomes are required. 
	  However, the power of the pipeline for selecting optimal genome loci 
     	  for phylogenomics improves when a larger number of genomes are available 
     	  for analysis. Reasonable numbers lie in the range of 10 to 100 clearly
     	  distinct genomes from multiple species of a genus, family, order or phylum.
     	  The pipeline may not perform satisfactorily with too distant genome sequences,
     	  particularly when sequences with significantly distinct nucleotide or aminoacid
     	  compositions are used. This type of sequence heterogeneity is well known to 
     	  cause systematic bias in phylogenetic inference. In general, very distantly related
     	  organisms, such as those from different phyla or even domains, may not be
     	  properly handled by run_get_phylomarkers_pipeline.sh.

## On the filtering criteria. 

*run_get_phylomarkers_pipeline.sh* uses a **hierarchical filtering scheme**, as follows:

###   i) Detection of recombinant loci. 

Codon or protein alignments (depending on runmode) are first screened with **Phi-test** 
([Bruen et al. 2006](http://www.genetics.org/content/172/4/2665.long)) for the presence of potential recombinant sequences. It is a well established fact that recombinant sequences negatively impact phylogenetic inference when using algorithms that do not account for the effects of this evolutionary force. The permutation test with 1000 permutations is used to compute the *p*-values. These are considerd significant if *p* < 0.05.
 
### ii) Detection of trees deviating from expectations of the (multispecies) coalescent.

The next filtering step is provided by the **kdetrees-test**, which checks the distribution of topologies, tree lengths and branch lenghts. *kdetrees* is a non-parametric method for estimating distributions of phylogenetic trees 
([Weyenberg et al. 2014](https://academic.oup.com/bioinformatics/article-lookup/doi/10.1093/bioinformatics/btu258)), 
with the goal of identifying trees that are significantly different from the rest of the trees in the sample, based on the analysis of topology and branch lenght distributions. Such "outlier" trees may arise for example from horizontal gene transfers or gene duplication (and subsequent neofunctionalization) followed by differential loss of paralogues among lineages. Such processes will cause the affected genes to exhibit a history distinct from those of the majority of genes, which are expected to be generated by the  (multispecies) coalescent as species or populations diverge. Alignments producing significantly deviating trees in the kdetrees test are identified and saved in the kde_outliers/ directory. The corresponding alignments are not used in downstream analyses.
      
      * Parameter for controlling kdetrees stingency:
      -k <real> kde stringency (0.7-1.6 are reasonable values; less is more stringent)
     			       [default: 1.5]

### iii) Phylogenetic signal content and congruence. 

The alignments passing the two previous filters are subjected to **maximum likelihood (ML) tree searches with FastTree** to 
infer the corresponding ML gene trees. Their **phylogenetic signal is computed from the Shimodaria-Hasegawa-like likelihood ratio test branch support values**, which vary between 0-1, as we have reported previously ([Vinuesa et al. 2008](http://aem.asm.org/content/74/22/6987.long)). The support values of each internal branch or bipartition are parsed to compute the mean support value for each tree. 
**Alignments/Trees with a mean support value below a cutoff threshold are discarded**. In addition, a consensus tree is computed from the collection of trees that passed filters *i* and *ii*, and the Robinson-Fould distance (**RF**) is computed between each gene tree and the consensus tree.

      * Parameters controlling filtering based on mean support values.
      -m <real> min. average support value (0.7-0.8 are reasonable values) 
     		for trees to be selected [default: 0.75]

### iv) Evaluating the global molecular clock hypothesis.

*run_get_phylomarkers_pipeline.sh* calls the auxiliary script *run_parallel_molecClock_test_with_paup.sh*
to evaluate the **global molecular clock hypothesis** on the topo markers, selected according to the criteria explained in the three previous
points. The script calls [paup*](https://people.sc.fsu.edu/~dswofford/paup_test/) 
to evaluate the free-rates and clock hypothesis using likelihood ratio tests using R. Currently this is only performed on codon alignments. Future versions will implement the global clock hypothesis test also for protein alignments.

### v) On tree searching: 

*run_get_phylomarkers_pipeline.sh* performs **tree searches under the maximum-likelihood criterion** 
using the [FastTree](http://microbesonline.org/fasttree/) program ([Price et al. 2010](http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0009490)). This program meets an excellent compromise between speed and accuracy, using modest ammounts of RAM. It accepts DNA and protein sequence alignments. It computes the above-mentioned Shimodaria-Hasegawa-like likelihood ratio test of branch support. 
A limitation though, is that it implements only very few substitution models. However, for divergent sequences of different species within a bacterial taxonomic genus or family, our experience has shown that almost invariably the GTR+G model is selected by jmodeltest2, particularly when there is base frequency heterogeneity. The GTR+G+CAT is the substitution model used by *run_get_phylomarkers_pipeline.sh* 
calls of FastTree on codon alignments. The gene trees are computed with high accuracy by performing a thorough tree search, as hardcoded in the following FastTree call:
      
     	-nt -gtr -gamma -bionj -slownni -mlacc 3 -spr 8 -sprlength 8 
     	
For concatenated codon alignments, which may take a considerable time (up to several hours) or
for large datasets (~ 100 taxa and > 300 concatenated genes) the user can choose to run FastTree with at different **levels of tree-search thoroughness**: high|medium|low|lowest 
      
      high:   -nt -gtr -bionj -slownni -gamma -mlacc 3 -spr 4 -sprlength 8
      medium: -nt -gtr -bionj -slownni -gamma -mlacc 2 -spr 4 -sprlength 8 
      low:    -nt -gtr -bionj -slownni -gamma -spr 4 -sprlength 8 
      lowest: -nt -gtr -gamma -mlnni 4
      
where -s $spr and -l $spr_length can be set by the user. The lines above show their default values.
      
For protein alignments, the search parameters are the same, only the model changes to lg
      
      high: -lg -bionj -slownni -gamma -mlacc 3 -spr 4 -sprlength 8
      
Please refer to the FastTree manual for the details.
      
## Citation.

We are preparing a publication describing the implementation get_phylomarkers pipeline and its use in genomic taxonomy and population genomics of *Stenotrophomonas* bacteria, which will be submitted to the [Research Topic of Frontiers in Microbiology: Microbial Taxonomy, Phylogeny and Biodiversity](http://journal.frontiersin.org/researchtopic/5493/microbial-taxonomy-phylogeny-and-biodiversity).

Meanwhile, if you find the code useful for your academic work, please use the following citation:
Pablo Vinuesa and Bruno Contreras-Moreira 2017. Get_PhyloMarkers, a pipeline to select optimal markers for microbial phylogenomics, systematics and genomic taxomy. Available at https://github.com/vinuesa/get_phylomarkers and released under the GNU GPLv3 license.

## Acknowledgements

### Personal
We  thank  Romualdo  Zayas,  Víctor  del  Moral,  and  Alfredo  J. Hernández at CCG-UNAM for technical support.

### Funding
We gratefully acknowledge the funding provided by [DGAPA-PAPIIT/UNAM](http://dgapa.unam.mx/index.php/impulso-a-la-investigacion/papiit) (grants IN201806-2 and IN211814) and [CONACyT-Mexico](http://www.conacyt.mx/) (grants P1-60071 and 179133) to [Pablo Vinuesa](http://www.ccg.unam.mx/~vinuesa/), as well as the Fundación ARAID, and Consejo  Superior  de Investigaciones Científicas (grant 200720I038) to [Bruno Conteras-Moreira](http://161.111.227.80/compbio/staff/bruno_contreras_moreira.html).
