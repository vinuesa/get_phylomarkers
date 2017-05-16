/*****************************************************************************/
/* tokenize.c 								     */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/* Standalone library for basic text processing.			     */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/* Copyright (C) 2001, 2004, Pal, A. (apal@szofi.elte.hu). 		     */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/* See function prototypes and the usage of the functions in tokenize.h	     */
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

#include "tokenize.h"

#ifdef HAVE_LIBINTL_H
#include <libintl.h>
#define _(string)	gettext(string)
#else
#define _(string)	string
#endif

#define	TOKENIZE_ABORT_ON_MALLOC
#ifdef	TOKENIZE_ABORT_ON_MALLOC
#define	realloc_check(ptr,size) \
	do \
	 {	if ( ptr==NULL && size>0 ) \
		 {	fprintf(stderr,"tokenize.c: %s.\n",_("memory exhausted"));\
			abort(); \
		 } \
	 } while(0)
#else
#define	realloc_check(ptr,size) 
#endif
#define	malloc_check(ptr)	realloc_check(ptr,1)

/*****************************************************************************/

void remove_newlines_and_comments(char *buff)
{
 int k;
 while ( *buff )
  {	if ( *buff=='#' )	*buff=0;
	else 
	 {	for ( k=0 ; buff[k]==10 || buff[k]==13 ; )	k++;
		if ( k )	memmove(buff,buff+k,strlen(buff)+1-k);
		else		buff++;
	 }
  }
	
}

void remove_spaces_and_comments(char *buff)
{
 int k;
 while ( *buff )
  {	if ( *buff=='#' )	*buff=0;
	else 
	 {	for ( k=0 ; buff[k]==9 || buff[k]==32 || buff[k]==10 || buff[k]==13 ; )	k++;
		if ( k )	memmove(buff,buff+k,strlen(buff)+1-k);
		else		buff++;
	 }
  }
	
}

void remove_spaces(char *buff)
{
 int k;
 while ( *buff )
  {	for ( k=0 ; buff[k]==9 || buff[k]==32 || buff[k]==10 || buff[k]==13 ; )	k++;
	if ( k )	memmove(buff,buff+k,strlen(buff)+1-k);
	else		buff++;
  }
}

/*****************************************************************************/

void remove_quotes(char *buff)
{
 int k;
 while ( *buff )
  {	for ( k=0 ; buff[k]=='"' ; )	k++;
	if ( k )	memmove(buff,buff+k,strlen(buff)+1-k);
	else		buff++;
  }
	
}

/*****************************************************************************/

int char_is_space(int c)
{
 if ( c==32 || c==13 || c==10 || c==9 )	return(1);
 else					return(0);
}

/*****************************************************************************/

int tokenize_spaces(char *buff,char **tokens,int max)
{
 int	intoken,inquota,n;
 char 	**tsave;

 tsave=tokens;

 intoken=0,inquota=0;n=0;
 while ( *buff && n<max )
  {	if ( ( ! char_is_space(*buff) ) && ! intoken )
	 {	*tokens=buff;
		intoken=!0,inquota=0;n++;
		if ( *buff=='"' )	inquota=!0;
		tokens++,buff++;
	 }
	else if ( intoken && ( (char_is_space(*buff) && inquota) || (!char_is_space(*buff)) ) )
	 {	if ( *buff=='"' )	inquota=!inquota;
		buff++;
	 }
	else if ( intoken && ! inquota && char_is_space(*buff) )
	 {	*buff=0,buff++;
		intoken=0;
	 }
	else	buff++;
  };
 *tokens=NULL;

 while ( *tsave != NULL )
  {	remove_quotes(*tsave);
	tsave++;
  };

 return(n);
}

char **tokenize_spaces_dyn(char *buff)
{
 int	intoken,inquota,i,n,nm;
 char 	**rtokens;

 nm=16;
 rtokens=(char **)malloc(sizeof(char *)*nm);
 malloc_check(rtokens);
 if ( rtokens==NULL )
	return(NULL);

 intoken=0,inquota=0;n=0;
 while ( *buff )
  {	if ( ( ! char_is_space(*buff) ) && ! intoken )
	 {	rtokens[n]=buff;
		intoken=!0,inquota=0;n++;
		if ( *buff=='"' )	inquota=!0;
		buff++;
		if ( n>=nm-1 )
		 {	nm+=16;
			rtokens=(char **)realloc(rtokens,sizeof(char *)*nm);
			realloc_check(rtokens,sizeof(char *)*nm);
		 }
	 }
	else if ( intoken && ( (char_is_space(*buff) && inquota) || (!char_is_space(*buff)) ) )
	 {	if ( *buff=='"' )	inquota=!inquota;
		buff++;
	 }
	else if ( intoken && ! inquota && char_is_space(*buff) )
	 {	*buff=0,buff++;
		intoken=0;
	 }
	else	buff++;
  };

 rtokens[n]=NULL;

 for ( i=0 ; i<n ; i++ )
  {	remove_quotes(rtokens[i]);	}

 return(rtokens);
}

int tokenize_char(char *buff,char **tokens,int tchar,int max)
{
 int	n;

 if ( *buff==0 )
  {	*tokens=NULL;return(0);		}
	
 n=1;*tokens=buff,tokens++;
 while ( *buff && n<max )
  {	if ( *buff != tchar )	buff++;
	else			*buff=0,buff++,*tokens=buff,tokens++,n++;
  };
 *tokens=NULL;
 return(n);
}

char **tokenize_char_dyn_wwt(char *buff,int tchar,int is_terminate)
{
 int	n;
 char	**tokens;

 if ( buff==NULL )
  {	return(NULL);				}

 tokens=(char **)malloc(sizeof(char *));
 if ( tokens==NULL )
	return(NULL);

 if ( *buff==0 )
  {	*tokens=NULL;return(tokens);		}
	
 n=0;tokens[n]=buff,n++;
 while ( *buff )
  {	if ( *buff != tchar )	buff++;
	else
	 {	if ( is_terminate )	*buff=0;
		buff++;
		tokens=(char **)realloc(tokens,sizeof(char *)*(n+1));
		malloc_check(tokens);
		tokens[n]=buff,n++;
	 }
  };
 tokens=(char **)realloc(tokens,sizeof(char *)*(n+1));
 malloc_check(tokens);
 tokens[n]=NULL;
 return(tokens);
}
char **tokenize_char_dyn(char *buff,int tchar)
{
 char	**ret;
 ret=tokenize_char_dyn_wwt(buff,tchar,1);
 return(ret);
}

/*****************************************************************************/
                            
                                                        

