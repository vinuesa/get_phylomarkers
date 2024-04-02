#!/usr/bin/env python
import dendropy, optparse, sys
from dendropy import treesim

def mylabel(x,i):
    "returns labels from taxon x"
    return x.label

##comand line parsing
parser = optparse.OptionParser()
parser.add_option("-n", "--pop-size", type="float", help="effective population size, N_e", metavar="NUMBER",default=1e4)
parser.add_option("-N", "--number-trees", type="int", help="number of trees to generate for each species tree", metavar="NUMBER",default=100)
parser.add_option("-o", "--output-file", help="output coalescent trees to FILE", metavar="FILE")
parser.add_option("-s", "--species-file", help="read species trees from FILE", metavar="FILE", default="species.nex")
(options,args) = parser.parse_args()

if options.output_file==None:
    genehandle = sys.stdout
else:
    genehandle = open(options.output_file,"w")

species = dendropy.TreeList.get_from_path(options.species_file,"nexus",as_rooted=True)
taxon_map = dendropy.TaxonSetMapping.create_contained_taxon_mapping(species.taxon_set, 1,contained_taxon_label_func=mylabel)

tlist = dendropy.TreeList()

for tree in species:
    for i in range(0,options.number_trees):
        gene_tree = treesim.contained_coalescent(tree,taxon_map,default_pop_size=options.pop_size)
        tlist.append(gene_tree)
        ##gene_tree.write_to_stream(genehandle,"newick",suppress_rooting=True)

tlist.write_to_stream(genehandle,"newick",suppress_rooting=True)
genehandle.close()

