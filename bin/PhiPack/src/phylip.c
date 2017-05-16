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
#include "phylip.h"
#include "misc.h"

/* Read strict phylip name.  in_file should point to beginning of line */
void read_strict_name(FILE *in_file,char *name, int capacity)
{
  int i;
  char ch='\0';
  for(i=0;(i<PHYLIP_SIZE) && (ch!=EOF);i++)
    {
      ch=fgetc(in_file);
      name[i]=ch;
    }

  if(ch==EOF)
    error("Taxa name appears to be cut off\n");

  /* Insert many null chars to be safe */ 
  for(;i<=capacity;i++)
    name[i]='\0';
  
}


/* Read relaxed phylip name.  in_file should point to beginning of line */
void read_relaxed_name(FILE *in_file,char *name, int capacity)
{
  int i=0;
  char ch='\0',next_ch;
  cbool double_space=FALSE;
  
  while((i<capacity) && (ch!=EOF) && (double_space == FALSE))
    {
      ch=fgetc(in_file);
      
      /* Check for end of string.. messy */
      if(isspace(ch))
	{
	  
	  next_ch=fgetc(in_file);
	  if(next_ch != EOF)
	    {
	      if(isspace(next_ch))
		double_space=TRUE;
	      else if(next_ch != EOF)
		ungetc(next_ch,in_file);
	    }
	}
      if(double_space==FALSE)
	name[i++]=ch;
    }
  
  
  
  if(ch==EOF)
    error("Taxa name appears to be cut off\n");

  /* Insert many null chars to be safe */ 
  for(;i<=capacity;i++)
    name[i]='\0';
  
}


void write_phylip(FILE* cur_stream, cbool strict, char **taxa_names,align_type **alignment, int num_taxa,int num_sites)
{
  int i;
  int j;
  
  fprintf(cur_stream,"%d %d\n",num_taxa,num_sites);
  for(i=0;i<num_taxa;i++)
    {
      if(strict == TRUE)
	{
	  for(j=0;j<10;j++)
	    {
	      if(taxa_names[i][j] != '\0')
		fprintf(cur_stream,"%c",taxa_names[i][j]);
	      else
		fprintf(cur_stream," ");
	    }
	  fprintf(cur_stream," ");
	}
      else
	fprintf(cur_stream,"%s\n",taxa_names[i]);
      
      for(j=0;j<num_sites;j++)
	{
/* 	  if((int)alignment[i][j] <=MAX_STATE ) */
	  fprintf(cur_stream,"%c",((char)(alignment[i][j] + CHAR_START)));
/* 	  else if((int)alignment[i][j] == MISSING_STATE) */
/* 	    fprintf(cur_stream,"%c",MISSING); */
/* 	  else if((int)alignment[i][j]== GAP_STATE) */
/* 	    fprintf(cur_stream,"%c",GAP); */
	  
	}
      fprintf(cur_stream,"\n");
    }
  
  
}


int read_seq_bit(FILE *in_file, int base_limit, align_type* states, int index, int *num_bases)
{
  int bases_read=num_bases[index];
  char ch;
  char s[250];
  
  ch=fgetc(in_file);

 
  while(ch != EOF && ch != '\n' && bases_read<base_limit)
    {
      if(!isspace(ch))
	{
	 /*  if(ch == MISSING) */
/* 	    { */
/* 	      states[bases_read]=(align_type)MISSING_STATE; */
/* 	    } */
/* 	  else if(ch == GAP) */
/* 	    { */
/* 	      states[bases_read]=(align_type)GAP_STATE; */
/* 	    } */
/* 	  else */
/* 	    { */
	  if((ch < CHAR_START) || (ch > CHAR_END))
	    {
	      sprintf(s,"Illegal state encountered: %c\n",ch);
	      error(s);
	    }
	  else
	    {
	      states[bases_read]=(align_type)(toupper(ch)-CHAR_START);
	    }
	  /*    }*/
	  bases_read++;
	}
      else
	{
	  ch=skip_non_newline_space(in_file);
	}
      ch=fgetc(in_file);
    }
  
  if(bases_read >= base_limit)
    {
     
      if((ch != EOF) && !isspace(ch))
	{

	  sprintf(s,"Seems like too many bases in seq %d.\nMake sure file is in Phylip file format.\nTry also using the strict and relaxed  options (-s and -r)",(index+1));
	  error(s);
	}
     
    }
  
  num_bases[index]=bases_read;
  return ch;
  
}

/* In: file_name   - name of file with sequences
       strict      - taxa names are either 10 chars or 1--MAX_SIZE chars
                     in first case first 10 chars name
		     in second case first x chars without 2 spaces are name 
		     (cannot be longer than MAX_SIZE)
  Out: taxa_name   - name of all the taxa
       alignment   - ordered by taxa
       num_taxa   
       num_sites

  Note: Line with
  10 200
  seq1234567 ADSFS

  is treated same as:
  10 2000
  [spaces]seq1234567 ADSFS
*/   


void read_phylip(const char* file_name, cbool strict, char ***taxa_names, align_type ***alignment, int *num_taxa, int *num_sites)
{
  char s[250];
  FILE *in_file;
  /* Keep track of how far we are in each sequence... */
  int *seq_counter;
  int i;
  int cur_seq=0;
  char ch;

  in_file=fopen(file_name,"r");
  if(in_file == NULL)
    {
      sprintf(s,"Could not open file %s",file_name);
      error(s);
    }
    

    
  if(fscanf(in_file,"%d%d",num_taxa,num_sites)==EOF)
    error("Could not read number of taxa and/or number of sites.  Make sure file is in phylip format\n");
  
  fprintf(stdout,"Allocating space for %d taxa and %d sites\n",*num_taxa,*num_sites);
  
  /* Allocate space for sequence and names*/
  (*taxa_names)=(char **)mcalloc((*num_taxa) , sizeof(char *) );

  (*alignment)=(align_type **)mcalloc((*num_taxa) , sizeof(align_type *) );
  seq_counter=(int *)mmalloc((*num_taxa) * sizeof(int)  );

  
  for(i=0;i<(*num_taxa);i++)
    {
      (*taxa_names)[i]=(char *)mcalloc((MAX_SIZE+1) , sizeof(char)  );
      (*alignment)[i]=(align_type *)mcalloc((*num_sites) , sizeof(align_type)  );
      seq_counter[i]=0;
    }
  
 
  /* Read in name + first part of sequence */
  for(i=0;i<(*num_taxa);i++)
   {
     /* Advance to first non--space */
     if(skip_non_newline_space(in_file)== EOF)
       error("Cannot read sequences.  Make sure file is in phylip format\n");
     if(skip_newlines(in_file)== EOF)
       error("Cannot read sequences.  Make sure file is in phylip format\n");
     
     if(strict == TRUE)
       {
	 read_strict_name(in_file, (*taxa_names)[i],MAX_SIZE);
       }
     else
       {
	 read_relaxed_name(in_file,(*taxa_names)[i],MAX_SIZE);
       }
     
     /*fprintf(stderr,"Read taxa:%s.\n",(*taxa_names)[i]);*/
     ch=skip_all_space(in_file);
     
     read_seq_bit(in_file, (*num_sites), (*alignment)[i],i,seq_counter);
     
   }

  /*fprintf(stderr,"Checking for more...\n");*/
  /* For interleaved format ... */
  ch=skip_all_space(in_file);
  while(ch!= EOF)
    {
     
      read_seq_bit(in_file, (*num_sites),(*alignment)[cur_seq],cur_seq,seq_counter);
     
      cur_seq++;
      cur_seq=cur_seq % (*num_taxa);
      ch=skip_all_space(in_file);
    
    }

  /* Check size of sequences */
  for(i=0;i<(*num_taxa);i++)
    if(seq_counter[i] != (*num_sites))
      {
	sprintf(s,"Number of sites read does not match number of sites expected for seq %s",(*taxa_names)[cur_seq]);
      }
  
  free(seq_counter);
}

