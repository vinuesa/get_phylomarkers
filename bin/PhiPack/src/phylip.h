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

#ifndef PHYLIP
#define PHYLIP

/* char skip_all_space(FILE *in_file); */


/* char skip_newlines(FILE *in_file); */

/* char skip_non_newline_space(FILE *in_file); */



/* void read_strict_name(FILE *in_file,char *name, int capacity); */




/* void read_relaxed_name(FILE *in_file,char *name, int capacity); */

/* int read_seq_bit(FILE *in_file, int base_limit, align_type* states,int index, int *num_bases); */



void read_phylip(const char* file_name, cbool strict, char ***taxa_names, align_type ***alignment, int *num_taxa, int *num_sites);

void write_phylip(FILE* cur_stream, cbool strict, char **taxa_names,align_type **alignment, int num_taxa,int num_sites);

#endif
