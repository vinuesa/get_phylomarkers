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
/*
   Copyright (c)2005, Trevor Bruen
   All rights reserved.                          

   Any feedback is very welcome.
   email: david.bryant@otago.ac.nz
*/

#include "global.h"
#include "stats.h"
#include "math.h"
#include "mem.h"

void identity_permutation(int *permutation, int n)
{
  int i;
 
  i=0;
  for(i=0;i<n;i++)
    permutation[i]=i;
}


void sample_permutation(int *permutation, int n)
{
  int i,temp,val;

  i=0;
  identity_permutation(permutation,n);
  
  for(i=0;i<n-1;i++)
    {
      val=rand() % (n-i) +i;
      //fprintf(stdout, "Val is between %d & %d - got %d\n",i,(n-1),val);
      temp=permutation[i];
      permutation[i]=permutation[val];
      permutation[val]=temp;
    }
  
}

/* Print out histogram of incompatibilities */
void inc_stats(int max_size,int num_chars,cbool verbose,inc_type **inc_matrix,int **counts, FILE *logfile)
{
  int i,j,total=0,act_max=0;
  double prop;
  
  *counts=(int *)mmalloc(max_size * sizeof(int));
  
  for(i=0;i<max_size;i++)
    (*counts)[i]=0;
  
  for(i=0;i<num_chars;i++)
    {
      for(j=i+1;j<num_chars;j++)
	{
	  (*counts)[inc_matrix[i][j]]=(*counts)[inc_matrix[i][j]] +1;
	  if(inc_matrix[i][j] > act_max)
	    act_max=inc_matrix[i][j];
	}
    }
  
  total=(num_chars*(num_chars-1))/2;
  
  
  /* Print out values  */
  if(verbose)
    {
      fprintf(stdout,"Distribution of scaled incompatibility scores:\n");
      fprintf(stdout,"Score (%c):\n",'%');
      for(i=0;i<=act_max;i++)
	{
	  prop=((double) ((*counts)[i]) )/((double)total) * 50.0;
	  
	  fprintf(stdout,"%2d   (%4.1f): ",i,prop*2);
	  fprintf(logfile,"%2d   (%4.1f): ",i,prop*2);
	  for(j=0;j<prop;j++)
	    {
	      fprintf(stdout, "o");
	      fprintf(logfile, "o");
	    }
	  fprintf(stdout,"\n");
	  fprintf(logfile,"\n");
	}
      fprintf(stdout,"\n");
      fprintf(logfile,"\n");
    }

  
  return;
}
  

double PHI(inc_type **inc_matrix,int *num_states,int *permutation, int num_chars, int k)
{
  int i,j;
  int p_i,p_j,states_i,states_j;
  double val;
  double score=0.0;
  int terms;
  
  for(i=0;i<(num_chars-1);i++)
    {
      p_i=permutation[i];
      states_i=num_states[p_i];
      
      for(j=1;j<=k;j++)
	{	  
	  if(i<(num_chars-j))
	    {
	      p_j=permutation[i+j];
	      states_j=num_states[p_j];
	      val=((double)((int)inc_matrix[p_i][p_j]));

	      score=score+val;
	    }
	}
    }
  terms=(k*(2*num_chars-k-1))/2;
  return score/((double)(terms));
}


double NSS(inc_type **inc_matrix, int *permutation, int num_chars)
{
  int i,j,count=0,val,next_val;
  int p_i,p_j,p_j_1;
  double score=0.0;
    
  for(i=0;i<num_chars;i++)
    {
      p_i=permutation[i];
      for(j=0;j<(num_chars-1);j++)
	{
	  p_j=permutation[j];
	  p_j_1=permutation[j+1];
	  val=(int)inc_matrix[p_i][p_j];
	  next_val=(int)inc_matrix[p_i][p_j_1];
	  if(((val>0) && (next_val >0)) || ((val ==0) && (next_val ==0)))
	    score++;
	  count++;
	}
    }
  return score/((double)count);
} 


void get_f_and_g(inc_type **inc_matrix,int *num_states,int num_chars, double *f_values,double *g_values)
 {
   int i,j;
   double fscore=0.0,gscore=0.0;
   double val;
   for(i=0;i<num_chars;i++)
     {
       fscore=0.0;
       gscore=0.0;
       f_values[i]=0.0;
       g_values[i]=0.0;

       for(j=0;j<num_chars;j++)
	 {
	   val=((double)((int)inc_matrix[i][j]));
	   fscore=fscore+val;
	   gscore=gscore+val*val;

	 }
       f_values[i]=fscore;
       g_values[i]=gscore;
     }
 }

void analytic_mean_variance(inc_type** inc_matrix, int* num_states, int num_chars,int k_val, double *mean, double *varnce)
{
  double coeff_1,coeff_2,coeff_3,sum1,sum2,sum3;
  int i;
  
  double n,top,bottom;
  double f_vals[num_chars],g_vals[num_chars];

  
  double k=(double)k_val;
  n=(double)num_chars;
 
  
  get_f_and_g(inc_matrix,num_states,num_chars,f_vals,g_vals);

  
  /* Figure out mean */
  sum1=0;
  for(i=0;i<num_chars;i++)
    {
      sum1=sum1+f_vals[i];
    }
  
  (*mean)=sum1/((double)(num_chars * (num_chars -1)));
 
  /* Figure out variance - one term at a time */
  
  /* First coeff */
  top=(double)(27*k*n-18*pow(k,2)+28*pow(k,2)*n-21*k*pow(n,2)-9*k+5*n-9*pow(k,3)-11*pow(n,2)+6*pow(n,3)+6*pow(k,3)*n-4*pow(k,2)*pow(n,2));
  bottom=(double)(k*pow((k+1-2*n),2)*pow((n-1),2)*(n-2)*(n-3)*pow(n,2));
  
  coeff_1=((2.0)/(3.0))*(top/bottom);
  
  /* Second coeff */
  top=(double)(8*pow(k,2)*n-14*pow(k,2)+39*k*n+19*n-21*k+3*pow(k,3)-15*k*pow(n,2)+6*pow(n,3)-21*pow(n,2)-4);  
  bottom=(double)(k*pow((2*n-k-1),2)*n*(n-1)*(n-2)*(n-3));
  
  coeff_2=((2.0)/(3.0))*(top/bottom);
  
  /* Third coeff */
  top=(double)(-18*k*n-2*pow(k,2)*n+16*pow(k,2)+6*pow(n,2)-10*n+2+15*k+3*pow(k,3));
  bottom=(double)(k*pow((k+1-2*n),2)*n*(n-1)*(n-2)*(n-3));
  coeff_3=(-4.0/3.0)*(top/bottom);
  
  /* Now figure out (sum(f))^2, sumg(g) and sum(f^2) */
  sum1=sum1*sum1;
  
  sum2=0;
  for(i=0;i<num_chars;i++)
    {
      sum2=sum2+(g_vals[i]);
    }
  
  sum3=0;
  for(i=0;i<num_chars;i++)
    {
      sum3=sum3+(f_vals[i])*(f_vals[i]);
    }

  *varnce=coeff_1*sum1+coeff_2*sum2+coeff_3*sum3;
}
