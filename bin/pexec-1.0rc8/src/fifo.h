/*****************************************************************************/
/* fifo.h								     */
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

#ifndef	__FIFO_H_INCLUDED
#define	__FIFO_H_INCLUDED	1

/*****************************************************************************/

typedef struct
 {	unsigned char	*buffer;
	size_t		size;
	size_t		rpnt,wrts;
 } fifo;

/*****************************************************************************/

/* fifo_init():
   This function initializes the FIFO 'f'.				     */
int	fifo_init(fifo *f);

/* fifo_write():
   This function pushes 'size' bytes from the buffer 'buffer' to the FIFO 'f'.
   The function returns 0 if the call was successful, otherwise it returns
   a nonzero value.							     */
int	fifo_write(fifo *f,void *buffer,size_t size);

/* fifo_available():
   This function returns the number of available bytes to read in the FIFO.  */
size_t	fifo_available(fifo *f);

/* fifo_read():
   This function tries to read 'size' bytes from the FIFO 'f' into the buffer 
   'buffer'. The number of bytes actually read is returned, which can be
   less than 'size' if the desired amount of data are not available 
   (i.e. fifo_available() is less than 'size').				     */
size_t	fifo_read(fifo *f,void *buffer,size_t size);

/* fifo_skip():
   Skips the next 'size' bytes from the FIFO. Equivalent to 
   fifo_read(f,NULL,size);						     */
size_t	fifo_skip(fifo *f,size_t size);

/* fifo_peek():
   This function tries to read 'size' bytes from the FIFO 'f' into the buffer 
   'buffer' but does not remove the data from the FIFO.		             */
size_t	fifo_peek(fifo *f,void *buffer,size_t size);

/* fifo_flush():
   This function copies all available data from FIFO to 'buffer'. The number
   of bytes copied is returned (which can be zero if the FIFO is empty). 
   Note that there should be enough space in 'buffer' to store the appropriate
   amount of data (i.e. which is fifo_available()).			     */
size_t	fifo_flush(fifo *f,void *buffer);

/* fifo_reduce():
   This function reduces the memory usage of the FIFO to the minimal.        */
size_t	fifo_reduce(fifo *f);

/* fifo_free():
   This function releases and resets the FIFO 'f'. All data in the FIFO
   are lost. The FIFO object 'f' can be re-used as an empty one without
   initialization: there is no need to call fifo_init() before pushing
   new data into the FIFO.						     */
int	fifo_free(fifo *f);

/*****************************************************************************/

#endif

/*****************************************************************************/
                    
