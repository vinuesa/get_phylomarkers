/*****************************************************************************/
/* fifo.c								     */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/* Simple FIFO (first-in-first-out) implementation.			     */
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
#include <stdlib.h>
#include <string.h>

#ifdef HAVE_LIBINTL_H
#include <libintl.h>
#define _(string)       gettext(string)
#else
#define _(string)       string
#endif

#include "fifo.h"

#define	FIFO_ABORT_ON_MALLOC
#ifdef	FIFO_ABORT_ON_MALLOC
#define	realloc_check(ptr,size) \
	do \
	 {	if ( ptr==NULL && size>0 ) \
		 {	fprintf(stderr,"fifo.c: %s.\n",_("memory exhausted"));\
			abort(); \
		 } \
	 } while(0)
#else
#define	realloc_check(ptr,size) 
#endif
#define	malloc_check(ptr)	realloc_check(ptr,1)

/*****************************************************************************/

int fifo_init(fifo *f)
{
 f->buffer=NULL;
 f->size=0;
 f->rpnt=0;
 f->wrts=0;
 return(0);
}

int fifo_free(fifo *f)
{
 if ( f->buffer != NULL )
	free(f->buffer);
 fifo_init(f);
 return(0);
}

/*****************************************************************************/

#define		FIFO_FRAGSIZE		256

int fifo_write(fifo *f,void *vbuffer,size_t size)
{
 unsigned char	*buffer=(unsigned char *)vbuffer;
 size_t		msize,wpnt,osize;

 if ( size > f->size-f->wrts )
  {	osize=f->size;
	f->size=f->wrts+size+FIFO_FRAGSIZE-1;
	f->size=FIFO_FRAGSIZE*(f->size/FIFO_FRAGSIZE);
	f->buffer=realloc(f->buffer,f->size);
	realloc_check(f->buffer,f->size);

	if ( f->rpnt+f->wrts > osize )
	 {	memmove(f->buffer+f->rpnt+(f->size-osize),f->buffer+f->rpnt,osize-f->rpnt);
		f->rpnt+=(f->size-osize);
	 }
  }

 while ( size>0 )
  {	wpnt=f->rpnt+f->wrts;
	if ( wpnt>=f->size )	wpnt-=f->size;
	msize=f->size-wpnt;
	if ( size<msize )	msize=size;
	if ( buffer != NULL )
	 {	memcpy(f->buffer+wpnt,buffer,msize);
		buffer+=msize;
	 }
	f->wrts+=msize;
	size-=msize;
  };

 /*
 if ( f->rpnt <= f->wpnt )
  {	if ( f->size-f->wpnt >= size )
	 {	if ( buffer != NULL )
			memcpy(f->buffer+f->wpnt,buffer,size);
		f->wpnt+=size;
	 }
	else if ( f->rpnt > (size-(f->size-f->wpnt)) )
	 {	if ( buffer != NULL )
		 {	if ( f->wpnt < f->size )
				memcpy(f->buffer+f->wpnt,buffer,f->size-f->wpnt);
			memcpy(f->buffer,buffer+(f->size-f->wpnt),size-(f->size-f->wpnt));
		 }
		f->wpnt=size-(f->size-f->wpnt);
	 }
	else
	 {	f->size=(f->wpnt-f->rpnt)+size+FIFO_FRAGSIZE;
		f->size=FIFO_FRAGSIZE*(f->size/FIFO_FRAGSIZE);
		f->buffer=realloc(f->buffer,f->size);
		if ( f->size-f->wpnt >= size )
		 {	if ( buffer != NULL )
				memcpy(f->buffer+f->wpnt,buffer,size);
			f->wpnt+=size;
		 }
		else 
		 {	if ( buffer != NULL )
			 {	if ( f->wpnt < f->size )
					memcpy(f->buffer+f->wpnt,buffer,f->size-f->wpnt);
				memcpy(f->buffer,buffer+(f->size-f->wpnt),size-(f->size-f->wpnt));
			 }
			f->wpnt=size-(f->size-f->wpnt);
		 }
	 }
  }
 else if ( f->rpnt-f->wpnt > size )
  {	if ( buffer != NULL )
		memcpy(f->buffer+f->wpnt,buffer,size);
	f->wpnt+=size;
  }
 else
  {	size_t	oldsize,diff;
	oldsize=f->size;
	f->size=(f->size-(f->rpnt-f->wpnt))+size+FIFO_FRAGSIZE;
	f->size=FIFO_FRAGSIZE*(f->size/FIFO_FRAGSIZE);
	f->buffer=realloc(f->buffer,f->size);
	diff=f->size-oldsize;
	memmove(f->buffer+f->rpnt+diff,f->buffer+f->rpnt,oldsize-f->rpnt);
	f->rpnt+=diff;
	if ( buffer != NULL )
		memcpy(f->buffer+f->wpnt,buffer,size);
	f->wpnt+=size;
  }
 */

 return(0);
}

size_t fifo_available(fifo *f)
{
 return(f->wrts);
}

size_t fifo_read(fifo *f,void *vbuffer,size_t size)
{
 unsigned char	*buffer=(unsigned char *)vbuffer;
 size_t		rsize,rsize0,msize;

 if ( size<f->wrts )	rsize=size;
 else			rsize=f->wrts;

 rsize0=rsize;
 while ( rsize>0 )
  {	msize=f->size-f->rpnt;
	if ( rsize<msize )	msize=rsize;
	if ( buffer != NULL )
	 {	memcpy(buffer,f->buffer+f->rpnt,msize);
		buffer+=msize;
	 }
	f->rpnt+=msize;
	f->wrts-=msize;
	if ( f->rpnt >= f->size )
		f->rpnt=0;
	rsize-=msize;
  }
 return(rsize0);
}

size_t fifo_skip(fifo *f,size_t size)
{
 return(fifo_read(f,NULL,size));
}

size_t fifo_peek(fifo *f,void *vbuffer,size_t size)
{
 size_t	rpnt,wrts,ret;
 rpnt=f->rpnt;
 wrts=f->wrts;
 ret=fifo_read(f,vbuffer,size);
 f->rpnt=rpnt;
 f->wrts=wrts;
 return(ret);
}

size_t fifo_flush(fifo *f,void *vbuffer)
{
 size_t	size;

 size=fifo_available(f);
 if ( size>0 )	fifo_read(f,vbuffer,size);
 return(size);
} 

/*****************************************************************************/
                              
