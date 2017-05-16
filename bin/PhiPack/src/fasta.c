
/*  
   Copyright (c)2005, Trevor Bruen

   Any feedback is very welcome.
   email: david.bryant@otago.ac.nz
*/
/* This file is part of PhiPack.

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
#include "fasta.h"
#include "misc.h"


/* Read fasta sequence name.  in_file should point to beginning of line */
/* Format: 
   >TAXA_NAME COMMENT 
*/
void read_sequence_name(FILE *in_file,char *name, int capacity)
{
  int i=0;
  char ch='\0';

  ch=fgetc(in_file);
  if(ch != '>')
    error("Sequence name should begin with >\n");
 
  ch=fgetc(in_file);
  
 
  while((i<capacity) && (ch!=EOF) && (ch!= '\n'))
    {
      name[i++]=ch;
      ch=fgetc(in_file);
    }
  
  if(ch == '\n')
    ungetc(ch,in_file);
  if(i==capacity)
    error("Taxa names too large for this routine\n");
  
  if(ch==EOF)
    error("Taxa name appears to be cut off\n");

  /* Insert many null chars to be safe */ 
  for(;i<=capacity;i++)
    name[i]='\0';
  // fprintf(stdout,"Read %s\n",name);
  //fflush(stdout);
}


void write_fasta(FILE* cur_stream,  char **taxa_names,align_type **alignment, int num_taxa,int num_sites)
{
  int i;
  int j;
  
  for(i=0;i<num_taxa;i++)
    {
      fprintf(cur_stream,">%s \n",taxa_names[i]);
      
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


char read_seq_part(FILE *in_file,  int base_limit, align_type*** alignment, int index, int *num_bases)
{
  int bases_read=0;
  int new_limit;
  char ch;
  char s[250];
  
  ch=fgetc(in_file);

 
  while(ch != EOF && ch != '>')
    {
      
      if(!isspace(ch))
	{
	  /* Re allocate memory */
	  if(bases_read == base_limit)
	    {
	      // fprintf(stderr,"Read %d bases - reallocating\n",bases_read);
	      new_limit=base_limit+base_limit;
	      (*alignment)[index]=(align_type *)mrealloc((*alignment)[index],(base_limit) * sizeof(align_type),new_limit*sizeof(align_type)  );
	      base_limit=new_limit;
	    }
/* 	  if(ch == MISSING) */
/* 	    { */
/* 	      (*alignment)[index][bases_read]=(align_type)MISSING_STATE; */
/* 	    } */
/* 	  else if(ch == GAP) */
/* 	    { */
/* 	      (*alignment)[index][bases_read]=(align_type)GAP_STATE; */
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
		  (*alignment)[index][bases_read]=(align_type)(toupper(ch)-CHAR_START);
		}
	/*     } */
	  bases_read++;
	}
      else
	{
	  ch=skip_all_space(in_file);
	}
      ch=fgetc(in_file);
    }
  if(ch != EOF)
    ungetc(ch,in_file);


  num_bases[index]=bases_read;
  //  fprintf(stderr,"Final tally is %d\n",bases_read);
  (*alignment)[index]=(align_type *)mrealloc((*alignment)[index],(base_limit) * sizeof(align_type),bases_read*sizeof(align_type)  );
  base_limit=bases_read;
  return ch;
  
}



/* In: file_name   - name of fasta file with sequences
   Out: taxa_name   - name of all the taxa
       alignment   - ordered by taxa
       num_taxa   
       num_sites
*/   


void read_fasta(const char* file_name,  char ***taxa_names, align_type ***alignment, int *num_taxa, int *num_sites)
{
  char s[250];
  FILE *in_file;

  /* Keep track of how many bases in each sequence... */
  int *seq_counter;
  int i;
 
  char ch;

  int cur_seqs;
  int num_seqs=1,new_num;
  int length=99;

  in_file=fopen(file_name,"r");
  if(in_file == NULL)
    {
      sprintf(s,"Could not open file %s",file_name);
      error(s);
    }
  
  
  
  /* Allocate space for sequence and names - 1 at a time*/
  (*taxa_names)=(char **)mcalloc(num_seqs , sizeof(char *) );
  
  (*alignment)=(align_type **)mcalloc(num_seqs , sizeof(align_type *) );
  seq_counter=(int *)mmalloc((num_seqs * sizeof(int)  ));
			     
  
  for(i=0;i<(num_seqs);i++)
    {
      (*taxa_names)[i]=(char *)mcalloc((MAX_SIZE+1) , sizeof(char)  );
      (*alignment)[i]=(align_type *)mcalloc(length , sizeof(align_type)  );
      seq_counter[i]=0;
    }
  
 
			     
  ch=skip_non_newline_space(in_file);
  cur_seqs=0;
  while(ch != EOF)
    {
      if(cur_seqs==num_seqs)
	{
	  new_num=num_seqs+1;
	  (*taxa_names)=(char **)mrealloc(*taxa_names,num_seqs *sizeof(char *),new_num*sizeof(char *) );
	  (*alignment)=(align_type **)mrealloc(*alignment,num_seqs * sizeof(align_type *),new_num*sizeof(align_type *));
	  seq_counter=(int *)mrealloc(seq_counter,(num_seqs * sizeof(int)  ),new_num*sizeof(int));
	  
	  for(i=num_seqs;i<(new_num);i++)
	    {
	      (*taxa_names)[i]=(char *)mcalloc((MAX_SIZE+1) , sizeof(char)  );
	      (*alignment)[i]=(align_type *)mcalloc(length , sizeof(align_type)  );
	      seq_counter[i]=0;
	    }
	  
	  num_seqs=new_num;

	}
      ch=skip_all_space(in_file);
      read_sequence_name(in_file,(*taxa_names)[cur_seqs],MAX_SIZE);
      /* Might have comment - remove */
      ch=skip_non_newline(in_file);
      /* Advance to sequence */
      ch=skip_newlines(in_file);
      /* Read sequence */
      ch=read_seq_part(in_file,  length,alignment, cur_seqs, seq_counter);
      
      cur_seqs++;
  
    }

  
  *num_taxa=cur_seqs;
  new_num=cur_seqs;
  


/*   /\* Fix allocated size *\/ */
/*   (*taxa_names)=(char **)mrealloc(*taxa_names,num_seqs *sizeof(char *),new_num*sizeof(char *) ); */
/*   (*alignment)=(align_type **)mrealloc(*alignment,num_seqs * sizeof(align_type *),new_num*sizeof(align_type *)); */
/*   seq_counter=(int *)mrealloc(seq_counter,(num_seqs * sizeof(int)  ),new_num*sizeof(int)); */
  
  /* Assume all have equal length - use first taxa */
  *num_sites=seq_counter[0];
  /* Check size of sequences */
  for(i=0;i<(*num_taxa);i++)
    if(seq_counter[i] != *num_sites)
      {
	sprintf(s,"Number of sites read does not match number of sites expected for seq %s",(*taxa_names)[i]);
	error(s);
      }
  fprintf(stdout,"Found %d sequences of length %d\n",*num_taxa,*num_sites);
  free(seq_counter);
}

