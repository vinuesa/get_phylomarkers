/*****************************************************************************/
/* format.c								     */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/* Simple *printf()-like implementation for arbitrary string formatting	     */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/* Copyright (C) 2007; Pal, A. (apal@szofi.elte.hu)			     */
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
#include <string.h>
#include <stdlib.h>
#include <ctype.h>
#include <stdarg.h>

#ifdef HAVE_LIBINTL_H
#include <libintl.h>
#define	_(string)	gettext(string)
#else
#define	_(string)	string
#endif

#include "format.h"

#define	FORMAT_ABORT_ON_MALLOC
#ifdef	FORMAT_ABORT_ON_MALLOC
#define	realloc_check(ptr,size) \
	do \
	 {	if ( ptr==NULL && size>0 ) \
		 {	fprintf(stderr,"format.c: %s.\n",_("memory exhausted"));\
			abort(); \
		 } \
	 } while(0)
#else
#define	realloc_check(ptr,size) 
#endif
#define	malloc_check(ptr)	realloc_check(ptr,1)

/*****************************************************************************/

#define	BLOCK_FORMAT		256

#define	FORMAT_ABORT_ON_MALLOC

/*****************************************************************************/

static int format_int_realloc(char **out,int *len,int *alen,int req)
{
 if ( (*len)+req > (*alen) )
  {	req+=(*len);
	*alen = BLOCK_FORMAT*((req+BLOCK_FORMAT-1)/BLOCK_FORMAT);
	*out  = realloc(*out,*alen);
	realloc_check(*out,*alen);
  }
 return(0);
}

static int format_get_parameters(char *format,int *d1,int *d2,int *fmt,char **next)
{
 int	set;
 char	*fstart;

 fstart=format;

 set=0;
 while ( isspace((int)(*format)) )	format++;
 if ( *format && sscanf(format,"%d",d1)==1 )
  {	set|=1;
	while ( *format=='-' || isdigit((int)(*format)) )	format++;
	while ( isspace((int)(*format)) )	format++;
  }
 if ( *format=='.' )			format++;
 while ( isspace((int)(*format)) )	format++;
 if ( *format && sscanf(format,"%d",d2)==1 )
  {	set|=2;
	while ( *format=='-' || isdigit((int)(*format)) )	format++;
	while ( isspace((int)(*format)) )	format++;
  }
 
 if ( islower((int)(*format)) )
  {	*fmt=(int)(*format);
	format++;
	if ( next != NULL )	*next=format;
	return(set);
  } 
 else
  {	if ( next != NULL )	*next=format;
	return(-1);
  }

}

/* call: format_check_if_formatted(fmt,"<c1><c2>..."); */
int format_check_if_formatted(char *format,char *fchars)
{
 int		fmt,d1,d2,set;
 char		*next;

 while ( *format )
  {	if ( *format=='%' && *(format+1) )
	 {	format++;
		set=format_get_parameters(format,&d1,&d2,&fmt,&next);
		if ( set>=0 )
			format=next;
		else if ( *format=='%' )
			format++;

		if ( strchr(fchars,fmt) != NULL )
			return(1);
	 }
	else
		format++;
  }

 return(0);
}	

static int format_int_add_string(char **out,int *len,int *alen,int set,int d1,char *vs)
{
 int	slen,req;

 slen=strlen(vs);
 if ( (set & 1) && abs(d1)>slen )
	req=abs(d1);
 else
	req=slen;

 format_int_realloc(out,len,alen,req);

 if ( (set & 1) && d1>slen )
  {	memset((*out)+(*len),32,d1-slen);
	strcpy((*out)+(*len)+(d1-slen),vs);
	*len+=req;
  }
 else if ( (set & 1) && -d1>slen )
  {	strcpy((*out)+(*len),vs);
	memset((*out)+(*len)+slen,32,(-d1)-slen);
	(*out)[(*len)+req]=0;
	*len+=req;
  }
 else
  { 	strcpy((*out)+(*len),vs);
	*len+=slen;
  }

 return(0);
}

/* call: format_replace(fmt,esc,'c1',TYPE1,value1,'c2',TYPE2,value2,...,0); */
char *format_replace(char *format,int is_escape,...)
{ 
 char		*out,*fstart,*next;
 int		len,alen,set,d1,d2,fmt,chr,type;
 va_list	ap;
 
 alen=BLOCK_FORMAT;
 out=(char *)malloc(alen);
 malloc_check(out);

 out[0]=0;
 len=0;

 while ( *format )
  {	if ( is_escape && *format=='\\' )
	 {	format++;
		format_int_realloc(&out,&len,&alen,1);
		switch ( *format )
		 {   case 'n':
			out[len]='\n',len++;
			break;
		     case 't':
			out[len]='\t',len++;
			break;
		     case 'r':
			out[len]='\r',len++;
			break;
		     case 0:
			out[len]='\\',len++;
			break;
		     default:
			out[len]=*format,len++;
		 }
		if ( *format )	format++;
	 }
	else if ( *format=='%' && *(format+1) )
	 {	fstart=format;
		format++;
		set=format_get_parameters(format,&d1,&d2,&fmt,&next);
		/* valid printf-like format */
		if ( set>=0 )
		 {	va_start(ap,is_escape);
			while ( (chr=va_arg(ap,int))>0 )
			 {	type=va_arg(ap,int);
				if ( chr==fmt )
					break;
				else
				 {	switch ( type )
					 {   case FORMAT_INTEGER:
						(void)va_arg(ap,int);
						break;	
					     case FORMAT_STRING:
						(void)va_arg(ap,char *);
						break;
					 }
				 }
			 }
			/* valid printf-like format and found in the list */
			if ( chr>0 )
			 {	int	vi;
				char	*vs;
				char	buff[256],sbuff[16];
				switch ( type )
				 {   case FORMAT_INTEGER:
					vi=va_arg(ap,int);
					if ( (set & 2) && d2>0 )
					 {	if ( d2>128 )	d2=128;
						sprintf(sbuff,"%%.%dd",d2);
						sprintf(buff,sbuff,vi);
					 }
					else
						sprintf(buff,"%d",vi);
					format_int_add_string(&out,&len,&alen,set,d1,buff);
					break;	
				     case FORMAT_STRING:
					vs=va_arg(ap,char *);
					if ( vs == NULL )
						vs="(null)";
					format_int_add_string(&out,&len,&alen,set,d1,vs);
					break;
				 }
				format=next;
			 }
			/* valid but not in the list: */
			else
			 {	format--; /* go back to the percent mark */
				format_int_realloc(&out,&len,&alen,next-format);
				memcpy(out+len,format,next-format);
				len+=next-format;
				format=next;
			 }
			va_end(ap);
		 }
		else if ( (*format)=='%' )
		 {	format++;
			format_int_realloc(&out,&len,&alen,1);
			out[len]='%';
			len++;
		 }
		else
		 {	format_int_realloc(&out,&len,&alen,2);
			out[len]='%';
			out[len+1]=(*format);
			format++;
			len+=2;
		 }
	 }
	else
	 {	format_int_realloc(&out,&len,&alen,1);
		out[len]=*format,len++;
		format++;
	 }
  };

 format_int_realloc(&out,&len,&alen,1);
 out[len]=0;

 return(out);
}

/*****************************************************************************/
                        
