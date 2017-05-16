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
#ifndef MAXCHI
#define MAXCHI

void maxchi(double windowScale,cbool verbose,FILE *logfile,fileType inFile,align_type **alignment, site *site_desc, int *site_states, char **taxa_names,int num_taxa, int num_sites, int num_trials, int *emp_MAXCHI);

void all_pairs_maxchi( align_type **alignment, int num_taxa, int num_sites, int windowsize, int *permutation, site *site_desc,double *maxvalue,char **taxa_names,cbool report_breakpoint,FILE *logfile);

#endif
