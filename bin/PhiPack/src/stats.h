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
#ifndef STATS
#define STATS

void identity_permutation(int *permutation, int n);


void sample_permutation(int *permutation, int n);


/* Print out histogram of incompatibilities */
void inc_stats(int max_size,int num_chars,cbool verbose,inc_type **inc_matrix,int **counts, FILE *logfile);

double PHI(inc_type **inc_matrix,int *num_states,int *permutation, int num_chars, int k);


double NSS(inc_type **inc_matrix, int *permutation, int num_chars);


void get_f_and_g(inc_type **inc_matrix,int *num_states,int num_chars, double *f_values,double *g_values);

void write_diag(const char *file_name,int k,double small_variance, double mean, inc_type** inc_matrix,int *num_states, int num_chars);

void analytic_mean_variance(inc_type** inc_matrix, int* num_states, int num_chars,int k_val ,double *mean, double *varnce) ;
#endif
