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
struct node_struct{
  
  int neighbourindex;
  struct node_struct *next;
};

typedef struct node_struct node;



/* Input:  alignment       - ordered by site (each site [0,k] + negative for gap)
           site_states     - number of states at each state
           char_a          - a site index
	   char_b          - another site index
	   num_sites       - total number of sites 
	   num_taxa        - number of taxa
   Output: 
           incompatibility - (pars_score - #states)+2
 */
int pair_score( align_type **alignment,  int *site_states,int char_a, int char_b, int num_sites, int num_taxa)
{
  /* Keep both list & matrix - one for quick adding, other for quick DFS */
  node **adjacency_list;
  cbool **adjacency_matrix;
  node *cur_node,*new_node,*next_node;
  int char_a_states = site_states[char_a];
  int char_b_states = site_states[char_b];
  int total_states = char_a_states+char_b_states;
  int char_a_val,char_b_val,i,j,edge_count;

  /* For DFS */
  node **DFS_adjacency;
  int* array_stack;
  cbool* marked;
  int potential_neighbour;
  cbool has_valid_neighbour;
  int comp_count;
  int top=0,cur_vertex=0;
  
  /* For score */
  
  int inc_score=-1;
  
  
  /* Initialize list and matrix */
  adjacency_list=(node **)mmalloc(total_states * sizeof(node *));
  for(i=0;i<total_states;i++)
    adjacency_list[i]=NULL;
  
  adjacency_matrix=(cbool **)mmalloc(total_states * sizeof(cbool *) );
  for(i=0;i<total_states;i++)
    {
      adjacency_matrix[i]=(cbool *)mmalloc(total_states * sizeof(cbool));
    }
  
  for(i=0;i<total_states;i++)
    {
      for(j=0;j<total_states;j++)
	{
	  adjacency_matrix[i][j]=FALSE;
	}
    }

  
  /* Initialize stuff for DFS... */
  DFS_adjacency=(node **)mmalloc(total_states * sizeof(node *) );
  array_stack=(int *)mmalloc(  (total_states) * sizeof(int));
  marked=(cbool *)mmalloc(total_states * sizeof(cbool) );


  /* Build adjacency list */
  edge_count=0;


  for(i=0;i<num_taxa;i++)
    {
      char_a_val=(int)(alignment[char_a][i]);
      /* Number vertices [0...char_a_states-1] then [char_a_states..total_states] */
      char_b_val=(int)(alignment[char_b][i]);
      
      /* Add the edge - if necessary */
      
      if((char_a_val <= MAX_STATE) && (char_b_val <= MAX_STATE))
	{
	  /* Increase index to "global index" */
	  char_b_val=char_b_val+char_a_states;
	

	  if(adjacency_matrix[char_a_val][char_b_val]==FALSE)
	    {
	
	      /* Update symmetric adjacency matrix (undirected graph)*/
	      adjacency_matrix[char_a_val][char_b_val]=TRUE;
	      adjacency_matrix[char_b_val][char_a_val]=TRUE;
	      edge_count++;
	      
	      /* Add to adjacency lists */
	      cur_node=adjacency_list[char_a_val];
	      new_node=(node *)mmalloc(sizeof(node) );
	      new_node->neighbourindex=char_b_val;
	      new_node->next=cur_node;
	      adjacency_list[char_a_val]=new_node;

	      /* And other list */
	      cur_node=adjacency_list[char_b_val];
	      new_node=(node *)mmalloc(sizeof(node) );
	      new_node->neighbourindex=char_a_val;
	      new_node->next=cur_node;
	      adjacency_list[char_b_val]=new_node;
	      
	      
	    }
	}
    }
  /* Now do DFS to count components */
  
  for(i=0;i<total_states;i++)
    {
      marked[i]=FALSE;
      DFS_adjacency[i]=adjacency_list[i];
    }

  top=-1;

  comp_count=0;

  for(i=0;i<total_states;i++)
    {
      if(!marked[i])
	{

	  comp_count++;
	  /* "push" index onto stack */
	  array_stack[++top]=i;

	  while(top >= 0)
	    {
	      cur_vertex=array_stack[top];
	      marked[cur_vertex]=TRUE;
	      has_valid_neighbour=FALSE;
	      while((DFS_adjacency[cur_vertex] != NULL) && (has_valid_neighbour == FALSE))
		{
		  potential_neighbour=(DFS_adjacency[cur_vertex])->neighbourindex;
		  if(marked[potential_neighbour]==FALSE)
		    {

		      array_stack[++top]=potential_neighbour;
		      has_valid_neighbour=TRUE;
		    }
		  else
		    {
		      DFS_adjacency[cur_vertex]=DFS_adjacency[cur_vertex]->next;
		    }
		}
	      
	      if(has_valid_neighbour == FALSE)
		{
		  top--;
		}

	    }
	}
    }


  /* De allocate adjacency list */
  for(i=0;i<total_states;i++)
    {

      cur_node=adjacency_list[i];
      
      while(cur_node!=NULL)
	  {
	    next_node=cur_node->next;
	    free(cur_node);
	    cur_node=next_node;
	  }
      /*   free(cur_node); ? */
      free(adjacency_matrix[i]); 
    }
  
  free(adjacency_list);
  free(DFS_adjacency);
  free(adjacency_matrix);
  
  free(marked);
  free(array_stack);
  
  /* For the pairwise incompatibility */
  
  inc_score=edge_count-total_states+comp_count;

  return inc_score;
}


 
