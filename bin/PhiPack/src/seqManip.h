/* This file is part of PhiPack.
 Copyright (c)2005, Trevor Bruen
 Foobar is free software: you can redistribute it and/or modify
 it under the terms of the GNU Lesser General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.
 
 Foobar is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU Lesser General Public License for more details.
 
 You should have received a copy of the GNU General Public License
 along with PhiPack.  If not, see <http://www.gnu.org/licenses/>.
 */
#include "global.h"
#ifndef SEQMANIP
#define SEQMANIP 


/* Check validity of alignment */
cbool validate_alignment(align_type **alignment, alignmentClass alignKind,int num_taxa, int num_sites);

void get_polymorphic( align_type **alignment,  site* site_desc,  int *site_states, int num_taxa, int num_sites, int *new_num_sites, align_type ***new_alignment, int **new_site_states, site **new_site_desc);

void get_seq_div(align_type **alignment, alignmentClass alignKind, int num_taxa, int num_sites, double *div);

/*void get_seq_div(align_type **alignment, int num_taxa, int num_sites, double *div);*/


void extract( align_type **alignment, int num_taxa, int num_sites, int start, int end, align_type ***new_alignment);

void get_unambig( align_type **alignment,  site* site_desc,  int *site_states, int num_taxa, int num_sites, int *new_num_sites, align_type ***new_alignment, int **new_site_states, site **new_site_desc);

/* In:   alignment       - alignment ordered by taxa
         site_desc       - a description of each site
         site_states     - the number of states at each site
         num_taxa.. num_inf  - number of taxa char, etc..
   Out:  new_alignment   - alignment ordered by taxa with only informative
         new_site_states - number of states at each informative site
	 new_site_desc   - description of states at each informative site
*/
void get_informative( align_type **alignment,  site* site_desc,  int *site_states, int num_taxa, int num_sites, int *new_num_sites, align_type ***new_alignment, int **new_site_states, site **new_site_desc);

/* In:  alignment       - ordered by taxa
   Out: site_states     - number of states at every site
        site_desc       - type of every site
	num_inf         - number of informative sites
*/
void find_states( align_type **alignment, alignmentClass alignKind, cbool ignoreMissing, int num_taxa, int num_sites,  site** site_desc, int** site_states);

     /*void find_states( align_type **alignment, int num_taxa, int num_sites,  site** site_desc, int** site_states);*/



/* In: alignment     - ordered by taxa
       site_states   - number of states at each site
       num_sites     - number of sites
       num_taxa      - number of taxa...
Out:
       new_alignment - ordered by character 
                       where each character goes from 
		       0..k, for k states, with NEGATIVE
		       states for missing/gaps...
		      
*/
/*void reorder_chars( align_type **alignment,  int num_sites, int num_taxa, align_type ***new_alignment);*/
void reorder_chars( align_type **alignment,  alignmentClass alignKind, cbool ignoreMissing,int num_sites, int num_taxa, align_type ***new_alignment);

#endif
