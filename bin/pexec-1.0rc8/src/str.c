/*****************************************************************************/
/* str.c								     */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/* Some functions related to (dynamic) string handling.			     */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/* Copyright (C) 2006; Pal, A. (apal@szofi.elte.hu)			     */
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

#include <string.h>
#include <stdlib.h>
#include <stdarg.h>
#include <stdio.h>

#include "str.h"

#ifdef HAVE_LIBINTL_H
#include <libintl.h>
#define _(string)	gettext(string)
#else
#define _(string)	string
#endif

#define	STR_ABORT_ON_MALLOC
#ifdef	STR_ABORT_ON_MALLOC
#define	realloc_check(ptr,size) \
	do \
	 {	if ( ptr==NULL && size>0 ) \
		 {	fprintf(stderr,"str.c: %s.\n",_("memory exhausted"));\
			abort(); \
		 } \
	 } while(0)
#else
#define	realloc_check(ptr,size) 
#endif
#define	malloc_check(ptr)	realloc_check(ptr,1)

/*****************************************************************************/

char * strkcpy(char *out,char *in,int size)
{
 strncpy(out,in,size);
 out[size-1]=0;
 return(out);
}

/*****************************************************************************/

int strappend(char **str,char *cat)
{
 int	l1,l2;

 if ( str==NULL )
	return(-1);
 else if ( *str==NULL )
  {	*str=strdup(cat);
	malloc_check(*str);
	return(0);
  }
 else
  {	l1=strlen(*str);
	l2=strlen(cat);
	*str=realloc(*str,l1+l2+1);
	malloc_check(*str);
	strcpy((*str)+l1,cat);
	return(0);
  }
}

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

int vstrappendf(char **str,char *fmt,va_list ap)
{
 int		n,l,size;

 if ( str==NULL )	return(0);

 if ( *str==NULL )	l=0;
 else			l=strlen(*str);
 size=128;

 *str=realloc(*str,l+size);
 realloc_check(*str,l+size);
 if ( *str==NULL )	return(-1);
 while ( 1 )
  {	n=vsnprintf((*str)+l,size,fmt,ap);
	if ( n>-1 && n<size )
		return(0);
	else if ( n>-1 )
		size=n+1;	
	else
		size=size*2;	
	if ( (*str=realloc(*str,l+size))==NULL )
		return(-1);
  };
 return(0);	
}

int strappendf(char **str,char *fmt,...)
{
 int		ret;
 va_list	ap;
 
 va_start(ap,fmt);
 ret=vstrappendf(str,fmt,ap);
 va_end(ap);

 return(ret);
}

/*****************************************************************************/
      
