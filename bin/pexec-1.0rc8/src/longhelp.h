/*****************************************************************************/
/* longhelp.h								     */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/* Copyright (C) 2008; Pal, A. (apal@szofi.elte.hu)			     */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/* Functions for creating nice ``long helps''				     */
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

#ifndef	__LONGHELP_H_INCLUDED
#define	__LONGHELP_H_INCLUDED	1

/*****************************************************************************/

/* these are specific terms which are supported by the ``help2man'' utility. */
#define		LONGHELP_OPTIONS	{ "Options:", NULL }
#define		LONGHELP_EXAMPLES	{ "Examples:", NULL }

/*****************************************************************************/

typedef struct
 {	char	*options;
	char	*description;
 } longhelp_entry;

/*****************************************************************************/

/* 
   longhelp_fprint():
   Prints the (NULL,NULL)-terminated longhelp_entry array `entries` to the
   file stream referred by `fw`. The output is formatted for `width` width,
   if it is positive, ignored at all if `width` is zero and calculated
   automatically if it is negative (in this case, if `fw` refers to a terminal,
   the width is derived automatically from the terminal window size argument,
   otherwise the formatting is ignored).
*/
int	longhelp_fprint(FILE *fw,longhelp_entry *entries,int flags,int width);

/*****************************************************************************/

#endif

/*****************************************************************************/
                             
