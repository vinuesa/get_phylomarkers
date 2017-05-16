/*****************************************************************************/
/* iof.c 								     */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/* Tools for high-level file manipulation (I/O with f*() functions).	     */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/* See function prototypes and the usage of the functions in iof.h	     */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/* Copyright (C) 2004, Pal, A. (apal@szofi.elte.hu). 			     */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*  This library is free software: you can redistribute it and/or modify     */
/*  it under the terms of the GNU General Public License as published by     */
/*  the Free Software Foundation, either version 3 of the License, or	     */
/*  (at your option) any later version.					     */
/*									     */
/*  This program is distributed in the hope that it will be useful,	     */
/*  but WITHOUT ANY WARRANTY; without even the implied warranty of	     */
/*  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the	     */
/*  GNU General Public License for more details.			     */
/*									     */
/*  You should have received a copy of the GNU General Public License	     */
/*  along with the program.  If not, see <http://www.gnu.org/licenses/>.     */
/*****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifdef HAVE_LIBINTL_H
#include <libintl.h>
#define _(string)	gettext(string)
#else
#define _(string)	string
#endif

#include "iof.h"

#define	IOF_ABORT_ON_MALLOC
#ifdef	IOF_ABORT_ON_MALLOC
#define	realloc_check(ptr,size) \
	do \
	 {	if ( ptr==NULL && size>0 ) \
		 {	fprintf(stderr,"iof.c: %s.\n",_("memory exhausted"));\
			abort(); \
		 } \
	 } while(0)
#else
#define	realloc_check(ptr,size) 
#endif
#define	malloc_check(ptr)	realloc_check(ptr,1)

/*****************************************************************************/

char *freadline(FILE *fr)
{
 char	*ret,buff[256];
 ret=NULL;
 if ( feof(fr) )	return(NULL);
 while ( 1 )
  {	if ( fgets(buff,255,fr)==NULL )	break;
	if ( ret==NULL )	ret=(char *)malloc(strlen(buff)+1),ret[0]=0;
	else			ret=(char *)realloc(ret,strlen(ret)+strlen(buff)+1);
	malloc_check(ret);
	strcat(ret,buff);
	if ( ret[strlen(ret)-1]==10 )	break;
  }
 return(ret);
}

char *freadline_bs(FILE *fr)
{
 char	*ret,*wr;
 int	n;
 ret=NULL;
 while ( 1 ) 
  {	wr=freadline(fr);
	if ( wr==NULL )	return(ret);
	else
	 {	if ( ret==NULL )	ret=wr;
		else
		 {	ret=realloc(ret,strlen(ret)+strlen(wr)+1);
			malloc_check(ret);
			strcat(ret,wr);
			free(wr);
		 }
		n=strlen(ret);
		if ( n>=2 )
		 {	if ( ret[n-1]==0x0A && ret[n-2]=='\\' )
				ret[n-2]=0;	
			else	
				return(ret);
		 }
		else	return(ret);
	 }
  };
 return(NULL);			
}

/*****************************************************************************/

FILE *fopenread(char *name)
{
 FILE	*fr;
 if ( name==NULL )			fr=stdin;
 else if ( strcmp(name,"-")==0 )	fr=stdin;
 else					fr=fopen(name,"rb");
 return(fr);
}
FILE *fopenwrite(char *name)
{
 FILE	*fw;
 if ( name==NULL )			fw=stdout;
 else if ( strcmp(name,"-")==0 )	fw=stdout;
 else					fw=fopen(name,"wb");
 return(fw);
}

int fcloseread(FILE *fr)
{
 if ( fileno(fr) != fileno(stdin) )	fclose(fr);
 return(0);
}        

int fclosewrite(FILE *fw)
{
 if ( fileno(fw) != fileno(stdout) )	fclose(fw);
 return(0);
}        

/*****************************************************************************/
      
