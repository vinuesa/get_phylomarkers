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
#include "mem.h"
#include "seqManip.h"

void get_polymorphic( align_type **alignment,  site* site_desc,  int *site_states, int num_taxa, int num_sites, int *new_num_sites, align_type ***new_alignment, int **new_site_states, site **new_site_desc)
{
  int i,j;
  int site_count=0,polymorphic_count=0;
  for(i=0;i<num_sites;i++)
    {
      if((site_desc[i].poly)==polymorphic)
	polymorphic_count++;
    }
  *new_num_sites=polymorphic_count;
  
  *new_site_states=(int *)mcalloc(polymorphic_count , sizeof(int));
  *new_site_desc=(site *)mcalloc(polymorphic_count , sizeof(site));
  *new_alignment=(align_type **)mcalloc(num_taxa , sizeof(align_type *));
  

  for(i=0;i<num_taxa;i++)
    ((*new_alignment)[i])=(align_type *)mcalloc(polymorphic_count , sizeof(align_type) );
  
  for(j=0;j<num_sites;j++)
    {

      if((site_desc[j].poly)==polymorphic)
	{
	  
	  for(i=0;i<num_taxa;i++)
	    {
	      (*new_alignment)[i][site_count]=alignment[i][j];
	    }
	  (*new_site_states)[site_count]=site_states[j];
	  (*new_site_desc)[site_count].inf=site_desc[j].inf;
	  (*new_site_desc)[site_count].poly=site_desc[j].poly;
	  (*new_site_desc)[site_count].gap=site_desc[j].gap;
	  (*new_site_desc)[site_count].orig_index=site_desc[j].orig_index;
	  (*new_site_desc)[site_count].num_missing=site_desc[j].num_missing;
	  site_count++;
	}
    }
 
}

void extract( align_type **alignment, int num_taxa, int num_sites, int start, int end, align_type ***new_alignment)
{
  int i,j;
  int site_count=0,new_count=end-start;

  *new_alignment=(align_type **)mcalloc(num_taxa , sizeof(align_type *));
  

  for(i=0;i<num_taxa;i++)
    ((*new_alignment)[i])=(align_type *)mcalloc(new_count , sizeof(align_type) );
  
  for(j=start;j<end;j++)
    {
      
      for(i=0;i<num_taxa;i++)
	{
	  (*new_alignment)[i][site_count]=alignment[i][j];
	}
      site_count++;
    }
  
}


void get_unambig( align_type **alignment,  site* site_desc,  int *site_states, int num_taxa, int num_sites, int *new_num_sites, align_type ***new_alignment, int **new_site_states, site **new_site_desc)
{
  int i,j;
  int site_count=0,nongapped_count=0;
  for(i=0;i<num_sites;i++)
    {
      if((site_desc[i].gap)==nongapped  && (site_desc[i].num_missing == 0))
	nongapped_count++;
    }
  *new_num_sites=nongapped_count;
  
  *new_site_states=(int *)mcalloc(nongapped_count , sizeof(int));
  *new_site_desc=(site *)mcalloc(nongapped_count , sizeof(site));
  *new_alignment=(align_type **)mcalloc(num_taxa , sizeof(align_type *));
  

  for(i=0;i<num_taxa;i++)
    ((*new_alignment)[i])=(align_type *)mcalloc(nongapped_count , sizeof(align_type) );
  
  for(j=0;j<num_sites;j++)
    {

      if(((site_desc[j].gap)==nongapped) && (site_desc[j].num_missing == 0))
	{
	  
	  for(i=0;i<num_taxa;i++)
	    {
	      (*new_alignment)[i][site_count]=alignment[i][j];
	    }
	  (*new_site_states)[site_count]=site_states[j];
	  (*new_site_desc)[site_count].inf=site_desc[j].inf;
	  (*new_site_desc)[site_count].poly=site_desc[j].poly;
	  (*new_site_desc)[site_count].gap=site_desc[j].gap;
	  (*new_site_desc)[site_count].orig_index=site_desc[j].orig_index;
	  (*new_site_desc)[site_count].num_missing=site_desc[j].num_missing;
	  site_count++;
	}
    }
  
}




/* In:   alignment           - alignment ordered by taxa
         site_desc           - a description of each site
         site_states         - the number of states at each site
         num_taxa.. num_inf  - number of taxa char, etc..
   Out:  new_num_sites       - the new number of sites
         new_alignment       - alignment ordered by taxa with only informative
         new_site_states     - number of states at each informative site
	 new_site_desc       - description of states at each informative site
*/

void get_informative( align_type **alignment,  site* site_desc,  int *site_states, int num_taxa, int num_sites, int* new_num_sites, align_type ***new_alignment, int **new_site_states, site **new_site_desc)
{
  int i,j;
  int site_count=0,informative_count=0;
   for(i=0;i<num_sites;i++)
    {
      if((site_desc[i].inf)==informative)
	informative_count++;
    }
  *new_num_sites=informative_count;

  *new_site_states=(int *)mcalloc(informative_count , sizeof(int));
  *new_site_desc=(site *)mcalloc(informative_count , sizeof(site));
  *new_alignment=(align_type **)mcalloc(num_taxa , sizeof(align_type *));
  
 
  for(i=0;i<num_taxa;i++)
    ((*new_alignment)[i])=(align_type *)mcalloc(informative_count , sizeof(align_type) );
  
  for(j=0;j<num_sites;j++)
    {
  
      if((site_desc[j].inf)==informative)
	{
		
	  for(i=0;i<num_taxa;i++)
	    {
	      (*new_alignment)[i][site_count]=alignment[i][j];
	    }
	  (*new_site_states)[site_count]=site_states[j];
	  (*new_site_desc)[site_count].inf=site_desc[j].inf;
	  (*new_site_desc)[site_count].poly=site_desc[j].poly;
	  (*new_site_desc)[site_count].gap=site_desc[j].gap;
	  (*new_site_desc)[site_count].orig_index=site_desc[j].orig_index;
	  (*new_site_desc)[site_count].num_missing=site_desc[j].num_missing;
	   
	  site_count++;
	}
    }
  
}

/* Check validity of alignment */
cbool validate_alignment(align_type **alignment, alignmentClass alignKind,int num_taxa, int num_sites)
{
  int i,j;
  char ch;
  char s[MAX_SIZE+1];

  for(i=0;i<num_taxa;i++)
    {
      for(j=0;j<num_sites;j++)
	{
	  ch=alignment[i][j];
	  if(!(validState(alignKind,ch) || missing_ambig_State(alignKind,ch) || (ch == (GLOBAL_GAP-CHAR_START))))
	    {
	      sprintf(s,"Error %c (%d) is not a valid state for this type of sequence\n",ch+CHAR_START,(int)ch+CHAR_START);
	      error(s);
	    }
	  
	}
      
    }
  return TRUE;

}
/* Estimate diversity using pairwise deletion */
void get_seq_div(align_type **alignment, alignmentClass alignKind, int num_taxa, int num_sites, double *div)
{
  int i,j,l,count=0,num=0;
  double val=0;
  
  /* fprintf(stdout,"Looking at alignment %d \n",alignKind);*/
  for(i=0;i<(num_taxa-1);i++)
    {
      for(j=(i+1);j<num_taxa;j++)
	{
	  count=0;
	  num=0;
	  for(l=0;l<num_sites;l++)
	    {
	      if(validState(alignKind,alignment[i][l]) && validState(alignKind,alignment[j][l]))
		/*   if((alignment[i][l] < MAX_STATE) && (alignment[j][l] <MAX_STATE)) */
		{
		  if(alignment[i][l]!=alignment[j][l])
		    count++;
		  num++;
		}
	    }
	  val=val+((double)count)/((double)num);
	  
	}
      
    }
  
  *div=(val*2)/(num_taxa*(num_taxa-1));
  //fprintf(stdout,"Found %lf\n",*div);
}


/* In:  alignment       - ordered by taxa
        num_taxa        - number of taxa
	num_sites       - number of sites
   Out: site_states     - number of states at every site
        site_desc       - type of every site

*/
void find_states( align_type **alignment, alignmentClass alignKind, cbool stateMissing, int num_taxa, int num_sites,  site** site_desc, int** site_states)
{
  int state_count,big_states, char_A,i,j;
  int state_map[MAX_STATE];
  int missing_ambig_count=0;
  char s[250];

  cbool gap_char;
  
  *site_desc=(site *)mcalloc(num_sites , sizeof(site)); 
  *site_states=(int *)mcalloc(num_sites , sizeof(int));

  for(j=0;j<num_sites;j++)
    {
      for(i=0;i<MAX_STATE;i++)
	state_map[i]=0;
      
      state_count=0;
      gap_char=FALSE;
      missing_ambig_count=0;


      for(i=0;i<num_taxa;i++)
	{
	  char_A=(int)alignment[i][j];

	  if(validState(alignKind,char_A))
	    /*   if(char_A <=MAX_STATE) */
	    {
	      if(state_map[char_A]<=0)
		{
		  state_count++;
		  state_map[char_A]++;
		  
		}
	      else
		state_map[char_A]++;
	    }
	  else if(missing_ambig_State(alignKind,char_A))
	    {
	      missing_ambig_count++;
	      
	    }
	  else if(valid_gap(alignKind,char_A))
	    {
	      gap_char=TRUE;
	    }
	  else
	    {
	      sprintf(s,"Invalid state %c\n",(char_A+CHAR_START));
	      error(s);
	    }
	  
	}
      /* Add current index to description */
      ((*site_desc)[j]).orig_index=j;
       
      if(stateMissing && (missing_ambig_count > 0))
	state_count++;
      
      /*Add state count */
      (*site_states)[j]=state_count;
      
      ((*site_desc)[j]).num_states=state_count;
      
      /* Add missing count */
      ((*site_desc)[j]).num_missing=missing_ambig_count;
      
      
      /* Check if gapped */
      if(gap_char)
	{
	  ((*site_desc)[j]).gap=gapped;
	}
      else
	{
	  ((*site_desc)[j]).gap=nongapped;
	}

      /* Check if informative */
      big_states=0;
      
      if(state_count >= 2)
	{
	  for(i=0;i<MAX_STATE;i++)
	    if(state_map[i]>= 2)
	      big_states++;
	  if(stateMissing && (missing_ambig_count >=2))
	    big_states++;
  
	}


      if(state_count>=2)
	{
	  ((*site_desc)[j]).poly=polymorphic;
	  if(big_states>=2)
	    {
	
	      ((*site_desc)[j]).inf=informative;
	    
	
	    }
	  else
	    {
	      
	      ((*site_desc)[j]).inf=uninformative;
	      
	      
	    }
	    
	}
      else
	{
	  ((*site_desc)[j]).inf=uninformative;

	  (*site_desc)[j].poly=constant;


	}
    }
 

 
  
}


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
void reorder_chars( align_type **alignment,  alignmentClass alignKind, cbool stateMissing,int num_sites, int num_taxa, align_type ***new_alignment)
{
  int state_map[MAX_STATE];
  int state_count,cur_state;
  int i,j,k,char_A;
  *new_alignment=(align_type **)mcalloc(num_sites , sizeof(align_type *));
  
  for(i=0;i<num_sites;i++)
    ((*new_alignment)[i])=(align_type *)mcalloc(num_taxa , sizeof(align_type));
 
  

  for(j=0;j<num_sites;j++)
    {
      for(k=0;k<MAX_STATE;k++)
	state_map[k]=-1;
      
      state_count=0;
      for(i=0;i<num_taxa;i++)
	{
	  char_A=(int)alignment[i][j];
	  
	  /* Figure out which positive state 0..k */
	  /*	  if(char_A <=MAX_STATE)*/
	  if(validState(alignKind,char_A))
	    {
	      if(state_map[char_A]<0)
		{
		  state_map[char_A]=state_count;
		  cur_state=state_count;
		  state_count++;
		  
		}
	      else
		{
		  cur_state=state_map[char_A];
		}
	    }
	  else if(missing_ambig_State(alignKind,char_A) && (stateMissing == TRUE))
	    {
	      if(state_map[MISSING_STATE]<0)
		{
		  state_map[MISSING_STATE]=state_count;
		  cur_state=state_count;
		  state_count++;
		  
		}
	      else
		{
		  cur_state=state_map[MISSING_STATE];
		}
	    }
	  /* Otherwise - state that should be ignored - use ignore state */
	  else
	    {
	      cur_state=IGNORE_STATE;
	    }
	  /* Add taxa */
	  (*new_alignment)[j][i]=(align_type)cur_state;
	}
      
    } 
      
}
