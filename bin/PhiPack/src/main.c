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

#include <time.h>

#include "global.h"
#include "options.h"
#include "misc.h"
#include  "pairScore.h"


#include "normal.h"
#include "stats.h"
#include "mem.h"
#include "phylip.h"
#include "fasta.h"
#include "seqManip.h"
#include "maxChi.h"
#include "graphCode.h"


void print_usage()
{
  
  fprintf(stderr, "Input file not specified!!! (use -f|-r|-s)\n");
  
  fprintf(stderr, "Usage: Phi   [-f|-s|-r] Filename  [-t] AlignmentType\n"); 
  fprintf(stderr, "             [-p [#]] [-o] [-v] [-g [#]] [-b [#]]\n\n");
  
  
  fprintf(stderr, "Options:\n"); 
  fprintf(stderr, "  -f: Filename = FASTA format\n");
  fprintf(stderr, "  -s: Filename = Strict phylip file\n");
  fprintf(stderr, "  -r: Filename = Relaxed phylip file\n");
  fprintf(stderr, "  -t: AlignmentType = D|A|O where D=DNA\n"); 
  fprintf(stderr, "                      A=AA and O=OTHER [default DNA]\n");
  
  fprintf(stderr, "  -p: [#] = PHI permutation test [default = FALSE, #=1000]\n");
  fprintf(stderr, "  -w: # = Change default window size [default w = 100]\n");
  
 
  fprintf(stderr, "  -o: Report other statistics (NSS and Max Chi^2) [default = FALSE]\n");

  fprintf(stderr, "  -v: Verbose [default = FALSE]\n");
 
  fprintf(stderr, "  -g: [i] = Print color (scaled) incompatibility matrix (graph.ppm)\n");
  fprintf(stderr, "          i - Image only (no ticks...) [default = FALSE]\n");
 
}


void get_params(int argc, char**argv, options *opt) 
{ 
  char *cur,ch,nextch;
  char temp[MAX_SIZE+1];
  int i;
  cbool inFileFound=FALSE;
  
  opt->otherStats=FALSE;
  opt->doPerm=FALSE; 
  opt->winSize=100;
  opt->k=0;
  opt->ntrials=1000;
  opt->alignKind=DNA;
  opt->printMatrix=FALSE;
  opt->graphType=4;

  opt->printBreakpoints=FALSE;
  opt->breakWindow=100;
  opt->embellish=TRUE;

  opt->verbose=FALSE;
  
  if(argc < 3)
    {
      print_usage();
      exit(1);
    }
  
  
  for(i=1;i<argc;i++)
    {
      cur=argv[i];
      ch=cur[0];
      
      //cout<<"Have "<<cur<<" and "<<ch<<endl;
      if(ch == '-')
	{
	  ch=toupper(cur[1]);
	  switch(ch)
	    {
	      /* Print other stats */
	    case 'O':
	      opt->otherStats=TRUE;
	      break;
	      
	      /* Type of sequence */
	    case 'T':
	      i++;
	      cur=argv[i];
	      nextch=toupper(cur[0]);
	      if(nextch == 'D')
		opt->alignKind=DNA;
	      else if(nextch == 'A')
		opt->alignKind=AA;
	      else if(nextch == 'O')
		opt->alignKind=OTHER;
	      else
		opt->alignKind=OTHER;
	      //Error
	      break;
	      
	      /* Fasta file */ 
	    case 'F':
	      if(!inFileFound)
		{
		  inFileFound=TRUE;
		  opt->inFile=fasta;
		  i++;
		  strsub(opt->seqFile,argv[i],0);
		  break;

		}
	      else
		{
		  //Error
		}
	      break;
	      
	      /* Strict phylip */
	    case 'S':
	      if(!inFileFound)
		{
		  inFileFound=TRUE;
		  opt->inFile=strict;
		  i++;
		  strsub(opt->seqFile,argv[i],0);
		  break;

		}
	      else
		{
		  //Error
		}
	      break;
	      
	      /* Relaxed phylip */
	    case 'R':
	       if(!inFileFound)
		{
		  inFileFound=TRUE;
		  opt->inFile=relaxed;
		  i++;
		  strsub(opt->seqFile,argv[i],0);
		  break;

		}
	      else
		{
		  //Error
		}
	       break;
	    
	    case 'V':
	      opt->verbose=TRUE;
	      break;
	      

	      
	      /* Or window size */
	    case 'W':
	      i++;
	      strsub(temp,argv[i],0);
	      opt->winSize=atof(temp);
	      break;
	   
	      /* Do permutations ? */
	    case 'P':
	      opt->doPerm=TRUE;
	      i++;
	      if(i<argc && argv[i][0] != '-')    
		{
		  strsub(temp,argv[i],0);
		  opt->ntrials=atoi(temp);
		}
	      else
		i--;
	      
	      break;

	    case 'G':
	      opt->printMatrix=TRUE;
	      i++;
	      if(i<argc && argv[i][0] == 'i')
		opt->embellish=FALSE;
	      else
		i--;
	      break;
	      
	    case 'B':
	      opt->printBreakpoints=TRUE;
	      i++;
	      if(i<argc && argv[i][0] != '-')    
		{
		  strsub(temp,argv[i],0);
		  opt->breakWindow=atoi(temp);
		  
		}
	      else 
		i--;
	      break;
	      
	    }
	}
      else
	{
	  //Error
	}
    }
  // Error if no file
  if(!(inFileFound))
    {
      print_usage();
      exit(1);
    }
  
} 
void print_vals(FILE *logfile,cbool print_val_a,cbool print_val_b,double val_a, double val_b)
{
  if(print_val_a)
     {
       fprintf(stdout,"%4.2le          ",val_a);
       fprintf(logfile,"%4.2le          ",val_a);
     }
  else
    {
      fprintf(stdout,"     --          ");
      fprintf(logfile,"     --          ");
    }
  if(print_val_b)
    {
      fprintf(stdout,"%4.2le\n",val_b);
      fprintf(logfile,"%4.2le\n",val_b);
    }
  else
    {
      fprintf(stdout,"--\n");
         fprintf(logfile,"--\n");
     }
 
}


int main(int argc, char* argv[])
{
  align_type **alignment;
  char **taxa_names;
  int num_taxa,num_sites;
  int num_inf;
  int pair_inc;
  double divg,val;

  /*Parameters */
  options opt;
  
  /* Number of states at each character and description of each site */
  int  *site_states;
  site *site_desc;

  /* Counts of incompatibilities */
  int *counts,max_state,max_inc;

  /* Alignment for informative states */
  align_type **inf_alignment;
  int *inf_states;
  site *inf_site_desc;

  /* Alignment by characters */
  align_type **char_alignment;
  
  /* Matrix of incompatibilities */
  int i,j;
  inc_type **inc_matrix;
  
  /* Log file */
  FILE *alignfile,*logfile;
 
  
  /* Values of statistics */
  cbool valid_normal_approx;
  double orig_PHI,orig_NSS, cur_PHI, cur_NSS;
  double sum_PHI=0,sum_sq_PHI=0,obs_mean=0,obs_varnce=0;
  double mean=0.0,variance=0.0;
  double difference,normal_p_val=0.0;
  int emp_PHI=0,emp_NSS=0, emp_MAXCHI=0;
 
  int *permutation;

  /* For maxChi */
  double windowScale= 2/3;
  
  /* Seed random number generator */
  srand(time(NULL));
     
  get_params(argc,argv,&opt);
  
  /* Open log file */

  logfile=ffopen("Phi.log","w");
  
  fprintf(stdout,"Reading sequence file %s\n",opt.seqFile);
  fprintf(logfile,"Reading sequence file %s\n",opt.seqFile);
  
  switch(opt.inFile)
    {
    case strict:
      read_phylip(opt.seqFile,TRUE, &taxa_names,&alignment,&num_taxa,&num_sites);
      break;
    case relaxed:
      read_phylip(opt.seqFile,FALSE, &taxa_names,&alignment,&num_taxa,&num_sites);
      break;
    case fasta:
      read_fasta(opt.seqFile, &taxa_names,&alignment,&num_taxa,&num_sites );
      break;
    }
  
  fprintf(logfile,"Allocating space for %d taxa and %d sites\n",num_taxa,num_sites);
  
  if(num_taxa<=0 || num_sites<=0)
    error("No sequences or sitesin alignment");
  
  
  /* Validate alignment type */
  validate_alignment(alignment,opt.alignKind,num_taxa,num_sites);
  fprintf(stdout,"Alignment looks like a valid ");
  fprintf(logfile,"Alignment looks like a valid ");
  switch(opt.alignKind)
    {
    case DNA:
      fprintf(stdout,"DNA alignment.\n");
      fprintf(logfile,"DNA alignment.\n");
      break;
    case AA:
      fprintf(stdout,"AA alignment.\n");
      fprintf(logfile,"AA alignment.\n");
      break;
    case OTHER:
      fprintf(stdout,"OTHER alignment.\n");
      fprintf(logfile,"OTHER alignment.\n");
      break;
    }

  /* FIX - Add complete deletion ?  */
  
  get_seq_div(alignment, opt.alignKind, num_taxa, num_sites, &divg);
  fprintf(stdout,"Estimated diversity is (pairwise deletion - ignoring missing/ambig): %4.1lf%%\n",(divg*100.0));
  fprintf(logfile,"Estimated diversity is (pairwise deletion - ignoring missing/ambig): %4.1lf%%\n",(divg*100.0));

  /* Get informative sites */
  find_states(alignment, opt.alignKind,FALSE,num_taxa, num_sites, &site_desc,  &site_states);
  
  get_informative(alignment, site_desc, site_states, num_taxa, num_sites, &num_inf, &inf_alignment, &inf_states, &inf_site_desc);
  
  fprintf(stdout,"Found %d informative sites.\n",num_inf);
  fprintf(logfile,"Found %d informative sites.\n",num_inf);
  
  /* Open seq file */

  alignfile=ffopen("Phi.inf.sites","w");
  
  fprintf(stdout  ,"Writing alignment of informative sites to: Phi.inf.sites\n");
  fprintf(logfile,"Writing alignment of informative sites to: Phi.inf.sites\n");
  
    
   switch(opt.inFile)
    {
    case strict:
      write_phylip(alignfile,TRUE, taxa_names,inf_alignment,num_taxa,num_inf);
      break;
    case relaxed:
      write_phylip(alignfile,FALSE, taxa_names,inf_alignment,num_taxa,num_inf);
      break;
    case fasta:
      write_fasta(alignfile, taxa_names,inf_alignment,num_taxa,num_inf );
      break;
    }
   fclose(alignfile);
   
   /* Write where those sites occured */
   alignfile=ffopen("Phi.inf.list","w");
  
   fprintf(stdout  ,"Writing list of informative sites to:      Phi.inf.list\n");
   fprintf(logfile, "Writing list of informative sites to:      Phi.inf.list\n");
   for(i=0;i<num_inf;i++)
     fprintf(alignfile,"%4d %4d\n",i+1,(inf_site_desc[i]).orig_index);
   
   fclose(alignfile);

   /* Allocate possibly sparsely - big matrix sequential stride */
  inc_matrix=(inc_type **)mcalloc(num_inf , sizeof(inc_type *) );
  for(i=0;i<num_inf;i++)
      inc_matrix[i]=(inc_type *)mcalloc( num_inf , sizeof(inc_type));
  
  /* Reorder by character */
  reorder_chars(inf_alignment, opt.alignKind,FALSE,num_inf, num_taxa, &char_alignment);
  
  
  switch(opt.alignKind)
    {
    case DNA:
      max_state=4;
      break;
    case AA:
      max_state=20;
      break;
    case OTHER:
      max_state=MAX_STATE;
      break;
    default:
      max_state=0;
    }
  max_inc=max_state*max_state-2*max_state+1;
  
  fprintf(stdout,"Calculating all pairwise incompatibilities...\n");
  /* Now get incompatibilities between all pairs... */
  for(i=0;i<num_inf;i++)
    {
      /* Fix / remove this... */
      if(i % 100 == 0)
	{
	  if(i != 0)
	    fprintf(stdout,"\b\b\b\b\b\b");
	  else
	    fprintf(stdout,"Done: ");
	  
	  val=((double)i)*((double)num_inf)-((double)i)*((double)(i+1))/2.0;
	  val=val/(((double)num_inf)*((double)num_inf-1)/2);
	  val=val*100.0;
	  fprintf(stdout,"%5.1f%%",val);
	  
	  fflush(stdout);
	}
      
      for(j=i;j<num_inf;j++)
 	{
	  if(i == j)
	    inc_matrix[i][j]=(inc_type)0;
	  else
	    {
	      pair_inc=pair_score(char_alignment, inf_states,i, j, num_inf, num_taxa);
	      
	      
	      inc_matrix[i][j]=(inc_type)pair_inc;
	      inc_matrix[j][i]=(inc_type)pair_inc;
	    }
 	}
    }
  
  fprintf(stdout,"\b\b\b\b\b\b");
  fprintf(stdout,"%5.1f%%",100.0);
  fflush(stdout);
  fprintf(stdout,"\n\n");

  inc_stats(max_inc,num_inf,opt.verbose,inc_matrix,&counts,logfile);

  if(opt.k == 0)
    {
      val=((double)num_inf)/((double)num_sites);
      val=val*opt.winSize;
      
      opt.k=(int)(val+0.5);
      
      if(opt.k == 0)
	opt.k++;
    }
  else
    {
      val=((double)num_sites)/((double)num_inf);
      val=val*opt.k;
      opt.winSize=val;
    }
  
  fprintf(stdout, "Using a window size of %3.0lf with k as %d\n",opt.winSize,opt.k);
  fprintf(logfile, "Using a window size of %3.0lf with k as %d\n",opt.winSize,opt.k);

   
   permutation=(int *)mmalloc(num_inf * sizeof(int));
   identity_permutation(permutation,num_inf);
   orig_PHI=PHI(inc_matrix,inf_states,permutation,num_inf,opt.k);
   
   if(num_inf <= 2*opt.k)
     {
       valid_normal_approx=FALSE;
       fprintf(stdout, "Too few informative sites to use normal approximation.\nTry doing a permutation test or increasing alignment length\nCan also try decreasing windowsize.\n\n");
       fprintf(logfile, "Too few informative sites to use normal approximation.\nTry doing a permutation test or increasing alignment length\nCan also try decreasing windowsize.\n\n");
     }
   else
     {
       valid_normal_approx=TRUE;
       fprintf(stdout,"\nCalculating analytical mean and variance\n");
       fprintf(logfile,"\nCalculating analytical mean and variance\n");
       
       /* Fix argument with 1 */
      analytic_mean_variance(inc_matrix,inf_states,num_inf,opt.k,&mean,&variance);
      
      difference=mean-orig_PHI;
      normal_p_val=normal_01_cdf((-difference)/sqrt(variance));
     }
   
   if(opt.doPerm)
     {
       if(opt.k >= num_inf)
	 {
	   error("Too few informative sites to test significance.  Try decreasing windowsize or increasing alignment length\n\n");
	 }
       fprintf(stdout,"\nDoing permutation test for PHI\n");
       fprintf(logfile,"\nDoing permutation test for PHI\n");
       for(i=0;i<opt.ntrials;i++)
	{
	  sample_permutation(permutation,num_inf);
	  cur_PHI=PHI(inc_matrix,inf_states,permutation,num_inf,opt.k);
	  if(cur_PHI <= orig_PHI)
	    emp_PHI++;
	  sum_PHI=sum_PHI+cur_PHI;
	  sum_sq_PHI=sum_sq_PHI+cur_PHI*cur_PHI;
	  
	}
       obs_mean=sum_PHI/((double)opt.ntrials);
       obs_varnce=(sum_sq_PHI-opt.ntrials*obs_mean*obs_mean)/(opt.ntrials-1);
       
     }
   
   if(opt.otherStats)
     {
       /* For NSS */
       fprintf(stdout,"\nDoing permutation test for NSS\n");
       fprintf(logfile,"\nDoing permutation test for NSS\n");
       
         identity_permutation(permutation,num_inf);
 
       orig_NSS=NSS(inc_matrix,permutation,num_inf);
       if(opt.verbose)
	 {
	   fprintf(stdout,"\nThe Neighbour Similarity score is %6.4le\n",orig_NSS);
	   fprintf(logfile,"\nThe Neighbour Similarity score is %6.4le\n",orig_NSS);
	 }

       for(i=0;i<opt.ntrials;i++)
	 {
	   sample_permutation(permutation,num_inf);
	   cur_NSS=NSS(inc_matrix,permutation,num_inf);
	   if(cur_NSS >= orig_NSS)
	    emp_NSS++;
	 }
       
       /* For Max Chi^2 */
       windowScale=2.0/3.0;
       maxchi(windowScale,opt.verbose,logfile,opt.inFile,alignment,site_desc,site_states,taxa_names,num_taxa,num_sites,opt.ntrials,&emp_MAXCHI);
    
     }
   
   if(opt.verbose)
     {
       fprintf(stdout,"\n                      PHI Values\n");
       fprintf(stdout,"                      ----------\n");
       
       fprintf(logfile,"\n                      PHI Values\n");
       fprintf(logfile,"                      ----------\n");
       
       if(!opt.doPerm)
	 i=0;
       else
	 i=opt.ntrials;
       
       fprintf(stdout,"              Analytical    (%d) Permutations\n\n",i);
       fprintf(logfile,"              Analytical    (%d) Permutations\n\n",i);
       
       fprintf(stdout,"Mean:          ");
       fprintf(logfile,"Mean:          ");
       
       print_vals(logfile,valid_normal_approx,opt.doPerm,mean,obs_mean);
       
       fprintf(stdout,"Variance:      ");
       fprintf(logfile,"Variance:      ");
       
       print_vals(logfile,valid_normal_approx,opt.doPerm,variance,obs_varnce);
       
       fprintf(stdout,"Observed:      %4.2le          %4.2le   \n\n",orig_PHI,orig_PHI);
        fprintf(logfile,"Observed:      %4.2le          %4.2le   \n\n",orig_PHI,orig_PHI);
     }
  
   if(opt.printMatrix)
     {
       fprintf(stdout, "\nOutputting incompatibility matrix to file matrix.ppm\n");
       fprintf(logfile, "\nOutputting incompatibility matrix to file matrix.ppm\n");
       
       create_incompat_pic(max_state,max_inc,inc_matrix,opt.embellish,num_inf);
       

     }

   fprintf(stdout,"\n     **p-Value(s)**     \n");
   fprintf(stdout,"       ----------\n\n");

  
   fprintf(logfile,"\n       p-Value(s)\n");
   fprintf(logfile,"       ----------\n\n");

   
   if(opt.otherStats)
     {
       val=((double)emp_NSS)/((double)(opt.ntrials));
       
       fprintf(stdout,"NSS:                 %4.2le  (%d permutations)\n",val,opt.ntrials);
       fprintf(logfile,"NSS:                 %4.2le  (%d permutations)\n",val,opt.ntrials);
       
       val=((double)emp_MAXCHI)/((double)(opt.ntrials));
       
       fprintf(stdout,"Max Chi^2:           %4.2le  (%d permutations)\n",val,opt.ntrials);
       fprintf(logfile,"Max Chi^2:           %4.2le  (%d permutations)\n",val,opt.ntrials);
       
     }
   if(opt.doPerm)
     {
       val=((double)emp_PHI)/((double)(opt.ntrials));
       fprintf(stdout,"PHI (Permutation):   %4.2le  (%d permutations)\n",val,opt.ntrials);
       fprintf(logfile,"PHI (Permutation):   %4.2le  (%d permutations)\n",val,opt.ntrials);
     }
   
   if(valid_normal_approx)
     {
       fprintf(stdout,"PHI (Normal):        %4.2le\n",normal_p_val);
       fprintf(logfile,"PHI (Normal):        %4.2le\n",normal_p_val);
     }
   else
     {
       fprintf(stdout,"PHI (Normal):        --\n");
       fprintf(logfile,"PHI (Normal):        --\n");
     }
     

   fprintf(stdout,"\n");
   fprintf(logfile,"\n");
   fclose(logfile);
   return 0;

}
