/*****************************************************************************/
/* linebuffer.c								     */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/* Copyright (C) 2007; Pal, A. (apal@szofi.elte.hu)			     */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/* Functions related to buffered line reading, from arbitrary streams	     */
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
#include <unistd.h>
#include <sys/select.h>
#include <errno.h>

#include "linebuffer.h"

#ifdef HAVE_LIBINTL_H
#include <libintl.h>
#define _(string)	gettext(string)
#else
#define _(string)	string
#endif

#define	LINEBUFFER_ABORT_ON_MALLOC
#ifdef	LINEBUFFER_ABORT_ON_MALLOC
#define	realloc_check(ptr,size) \
	do \
	 {	if ( ptr==NULL && size>0 ) \
		 {	fprintf(stderr,"linebuffer.c: %s.\n",_("memory exhausted"));\
			abort(); \
		 } \
	 } while(0)
#else
#define	realloc_check(ptr,size) 
#endif
#define	malloc_check(ptr)	realloc_check(ptr,1)

/*****************************************************************************/

int linebuffer_reset(linebuffer *lb)
{
 lb->buffer=NULL;
 lb->length=0;
 return(0);
}
int linebuffer_free(linebuffer *lb)
{
 if ( lb->buffer != NULL )
	free(lb->buffer);
 linebuffer_reset(lb);
 return(0);
}
int linebuffer_concatenate(linebuffer *lb,char *buff,size_t size)
{
 if ( buff==NULL || size<=0 )
	return(0);
 lb->buffer=(char *)realloc(lb->buffer,lb->length+size+1);
 malloc_check(lb->buffer);
 memcpy(lb->buffer+lb->length,buff,size);
 lb->length+=size;
 lb->buffer[lb->length]=0;
 return(0);
}

char * linebuffer_read_line(int fd,linebuffer *lb,int timeout)
{
 char	*line,*eoc,buff[256];
 int	length,n;
 struct	timeval	tv,*ptv,tn;
 fd_set	set;

 if ( timeout>0 )
  {	tv.tv_sec=timeout;
	tv.tv_usec=0;
	ptv=&tv;
  }
 else
	ptv=NULL;

 while ( 1 )
  {	if ( lb->buffer != NULL )
		eoc=memchr(lb->buffer,10,lb->length);
	else
		eoc=NULL;

	if ( eoc != NULL )
	 {	length=(int)(eoc-lb->buffer);
		line=(char *)malloc(1+length);
		malloc_check(line);
		memcpy(line,lb->buffer,length);
		line[length]=0;
		length++;
		if ( length>=lb->length )
		 {	free(lb->buffer);
			lb->buffer=NULL;
			lb->length=0;
		 }
		else
		 {	memmove(lb->buffer,lb->buffer+length,lb->length-length);
			lb->length -= length;
		 }
		return(line);
	 }

	FD_ZERO(&set);
	FD_SET(fd,&set);
	tn.tv_sec=0;
	tn.tv_usec=0;
	select(fd+1,&set,NULL,NULL,&tn);
	if ( FD_ISSET(fd,&set) )
	 {	n=read(fd,buff,256);
		if ( n>0 )
			linebuffer_concatenate(lb,buff,n);
		else if ( ! ( n<0 && errno==EINTR ) )
			return(NULL);
	 }
	else
	 {	FD_ZERO(&set);
		FD_SET(fd,&set);
		select(fd+1,&set,NULL,NULL,ptv);
		if ( FD_ISSET(fd,&set) )
		 {	n=read(fd,buff,256);
			if ( n>0 )
				linebuffer_concatenate(lb,buff,n);
			else if ( ! ( n<0 && errno==EINTR ) )
				return(NULL);
		 }
		else
			return(NULL);
	 }
  }
 return(NULL);	/* unreachable */ 
}

char *linebuffer_fetch(linebuffer *lb)
{
 char	*eoc,*line;
 size_t	length;

 if ( lb->buffer != NULL && lb->length>0 )
	eoc=memchr(lb->buffer,10,lb->length);
 else
 	eoc=NULL;

 if ( eoc != NULL )
  {	length=(size_t)(eoc-lb->buffer);
	line=(char *)malloc(1+length);
	malloc_check(line);
	memcpy(line,lb->buffer,length);
	line[length]=0;
	length++;
	if ( length>=lb->length )
	 {	free(lb->buffer);
		lb->buffer=NULL;
		lb->length=0;
	 }
	else
	 {	memmove(lb->buffer,lb->buffer+length,lb->length-length);
		lb->length -= length;
	 }
	return(line);
  }
 else
	return(NULL);
 
}

char * linebuffer_flush(linebuffer *lb)
{
 char	*line;
 size_t	length;

 if ( lb->buffer != NULL && lb->length>0 )
  {	line=(char *)malloc(lb->length+1);
	malloc_check(line);
	length=lb->length;
	memcpy(line,lb->buffer,length);
	line[length]=0;
	free(lb->buffer);
	linebuffer_reset(lb);
	return(line);
  }
 else
  {	if ( lb->buffer != NULL )
		free(lb->buffer);
	linebuffer_reset(lb);
	return(NULL);
  }
	
}

/*****************************************************************************/

