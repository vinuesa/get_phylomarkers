/*****************************************************************************/
/* linebuffer.h								     */
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

#ifndef	__LINEBUFFER_H_INCLUDED
#define	__LINEBUFFER_H_INCLUDED	1

/*****************************************************************************/

typedef struct
 {	char	*buffer;
	size_t	length;
 } linebuffer;

/*****************************************************************************/

int	linebuffer_reset(linebuffer *lb);
int	linebuffer_free(linebuffer *lb);

int	linebuffer_concatenate(linebuffer *lb,char *buff,size_t size);

char *	linebuffer_read_line(int fd,linebuffer *lb,int timeout);

char *	linebuffer_fetch(linebuffer *lb);
char *	linebuffer_flush(linebuffer *lb);

/*****************************************************************************/

#endif

/*****************************************************************************/
   
