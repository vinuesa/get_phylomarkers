/*  
   Copyright (c)2005, Trevor Bruen
   All rights reserved.                          

   Any feedback is very welcome.
   email: david.bryant@otago.ac.nz
*/
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
#include "mem.h"

#include "seqManip.h"
//#include "misc.h"
#include "stats.h"
#include "phylip.h"
#include "fasta.h"
#include "queue.h"
#include "maxChi.h"


 
/* Calculate best breakpoint considering all pairs of sequences*/
void all_pairs_maxchi( align_type **alignment, int num_taxa, int num_sites, int windowsize, int *permutation, site *site_desc,double *maxvalue,char **taxa_names,cbool report_breakpoint,FILE *logfile)
{
  int i,j,l;
  queue all_poly;

  int n=windowsize;
  int k=(int)(n/2);
  
  int r,s;
  queue_elem diff;
  
  double best_val=0,cur_val=0;

  /* Remove ? */
  /* Sequences and breakpoint */
  int best_i=0;
  int best_j=0;
  int best_l=0;
  
  int best_s=0;
 int  best_r=0;

  /* For efficient calculation use queue */
  all_poly.max_size=n;
   
  init_queue(&all_poly);
 
  for(i=0;i<(num_taxa-1);i++)
    {
      for(j=(i+1);j<num_taxa;j++)
	{
	  /* Set up scan */
	  s=0;
	  r=0;
	  clear_queue(&all_poly);
	
	  /* First initialize values queue */
	  for(l=0;l<n;l++)
	    {
	      if(alignment[i][permutation[l]] == alignment[j][permutation[l]])
		diff=0;
	      else
		{
		  diff=1;
		}
	      s=s+(int)diff;
	      enqueue(&all_poly,diff);
	      
	      if(l<k)
		{
		  r=r+(int)diff;
		}
	    }

	  cur_val=((double)n*(k*s-n*r)*(k*s-n*r))/((double)(k*s)*(n-k)*(n-s));
	  if(cur_val > best_val)
	    {
	      best_s=s;
	      best_r=r;
	
	      best_i=i;
	      best_j=j;
	      best_l=l;
	      best_val=cur_val;
	    }
	  // fprintf(stdout,"i and j are %d and %d. r and s are %d and %d. val is %lf\n",i,j,r,s,cur_val);
	  for(;l<num_sites;l++)
	    {
	      dequeue_front(&all_poly,&diff); 
	      s=s-(int)diff;
	      r=r-(int)diff;
	      
	      if(alignment[i][permutation[l]] == alignment[j][permutation[l]])
		diff=0;
	      else
		diff=1;

	      s=s+(int)diff;
 	      enqueue(&all_poly,diff);
	      
	      /* Update value of r as well */
	      element_at(&all_poly,&diff,(k-1));
	      r=r+(int)diff;
	      cur_val=((double)n*(k*s-n*r)*(k*s-n*r))/((double)(k*s)*(n-k)*(n-s));
	      if(cur_val > best_val)
		{
		  best_s=s;
		  best_r=r;
		  best_i=i;
		  best_j=j;
		  best_l=l;
		  best_val=cur_val;
		}
	    }
	}
    }

  
  if(report_breakpoint)
    {

      
      fprintf(stdout,"\nBest breakpoint for Max Chi found with sequences %s and %s. r and s are %d and %d\n",taxa_names[best_i],taxa_names[best_j],best_r,best_s);  
      fprintf(logfile,"\nBest breakpoint for Max Chi found with sequences %s and %s.r and s are %d and %d\n",taxa_names[best_i],taxa_names[best_j],best_r,best_s);  
      
      fprintf(stdout,"Value of maximum breakpoint is:    %3.1lf\n\n",best_val);
      fprintf(logfile,"Value of maximum breakpoint is:    %3.1lf\n\n",best_val);
      
      
      fprintf(stdout,"Coordinates of breakpoint with only polymorphic sites (start,breakpoint,end) = (%d, %d, %d)\n",best_l-n,best_l-k,best_l);
      fprintf(stdout,"Coordinates of breakpoint with all sites (start,breakpoint,end)=(%d, %d, %d)\n",site_desc[best_l-n].orig_index,site_desc[best_l-k].orig_index,site_desc[best_l].orig_index);
      
      fprintf(logfile,"Coordinates of breakpoint with only polymorphic sites (start,breakpoint,end) = (%d, %d, %d)\n",best_l-n,best_l-k,best_l);
      fprintf(logfile,"Coordinates of breakpoint with all sites (start,breakpoint,end)=(%d, %d, %d)\n",site_desc[best_l-n].orig_index,site_desc[best_l-k].orig_index,site_desc[best_l].orig_index);

    }
  
  *maxvalue=best_val;
}

void maxchi(double windowScale,cbool verbose,FILE *logfile,fileType inFile,align_type **alignment, site *site_desc, int *site_states, char **taxa_names,int num_taxa, int num_sites, int num_trials, int *emp_MAXCHI)
{
  int i;
  double orig_MAXCHI,cur_MAXCHI;
  int maxChiWindow;
  int *permutation;
  int num_poly,num_poly_unambig;
  align_type **poly_alignment,**poly_alignment_unambig;
  int *poly_sites,*poly_sites_unambig;
  site*poly_site_desc,*poly_site_desc_unambig;
  FILE *alignfile;
  fprintf(stdout,"\nDoing Permutation test for MAXCHI\n\n");
  fprintf(logfile,"\nDoing Permutation test for MAXCHI\n\n");
  
  get_polymorphic( alignment,  site_desc, site_states, num_taxa, num_sites, &num_poly, &poly_alignment, &poly_sites, &poly_site_desc);
  
  get_unambig( poly_alignment,  poly_site_desc,poly_sites, num_taxa, num_poly, &num_poly_unambig, &poly_alignment_unambig, &poly_sites_unambig, &poly_site_desc_unambig);
  
  if(verbose)
    {
  fprintf(stdout, "Number of umabiguous polymorphic sites is %d\n",num_poly_unambig);
  fprintf(logfile,"Number of umabiguous polymorphic sites is %d\n",num_poly_unambig);
    }
  
  alignfile=ffopen("Phi.poly.unambig.sites","w");
  
  fprintf(stdout  ,"Writing  alignment of polymorphic unambig sites to: Phi.poly.sites\n");
  fprintf(logfile," Writing  alignment of polymorphic unambig sites to: Phi.poly.sites\n");
  
  switch(inFile)
    {
    case strict:
      write_phylip(alignfile,TRUE, taxa_names,poly_alignment_unambig,num_taxa,num_poly_unambig);
      break;
    case relaxed:
      write_phylip(alignfile,FALSE, taxa_names,poly_alignment_unambig,num_taxa,num_poly_unambig);
      break;
    case fasta:
      write_fasta(alignfile, taxa_names,poly_alignment_unambig,num_taxa,num_poly_unambig );
      break;
    }
  fclose(alignfile);
  
  permutation=(int *)mmalloc(num_poly_unambig * sizeof(int));
  identity_permutation(permutation,num_poly_unambig);
  
  maxChiWindow=(int)(((double)num_poly_unambig) * windowScale);
  fprintf(stdout,"Window size is %d polymorphic sites\n",maxChiWindow);
  fprintf(logfile,"Window size is %d polymorphic sites\n",maxChiWindow);

  if(verbose)
    {
  all_pairs_maxchi(poly_alignment_unambig, num_taxa, num_poly_unambig,maxChiWindow,permutation,poly_site_desc_unambig,&orig_MAXCHI,taxa_names,TRUE,logfile);
    }
  else
    {
       all_pairs_maxchi(poly_alignment_unambig, num_taxa, num_poly_unambig,maxChiWindow,permutation,poly_site_desc_unambig,&orig_MAXCHI,taxa_names,
FALSE,logfile);
    }
  *emp_MAXCHI=0;
  
  for(i=0;i<num_trials;i++)
    {
      sample_permutation(permutation,num_poly_unambig);
      all_pairs_maxchi(poly_alignment_unambig, num_taxa, num_poly_unambig,maxChiWindow,permutation,poly_site_desc_unambig,&cur_MAXCHI,taxa_names,FALSE,logfile);
      if(cur_MAXCHI >= orig_MAXCHI)
	(*emp_MAXCHI)++;
    }
  return;
}
