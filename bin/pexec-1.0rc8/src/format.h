/*****************************************************************************/
/* format.h								     */
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

#ifndef	__FORMAT_H_INCLUDED
#define	__FORMAT_H_INCLUDED	1

/*****************************************************************************/

#define	FORMAT_INTEGER		1
#define	FORMAT_STRING		2

/*****************************************************************************/

/* call: format_check_if_formatted(fmt,"<c1><c2>..."); */
int	format_check_if_formatted(char *format,char *fchars);

/* call: format_replace(fmt,esc,'c1',TYPE1,value1,'c2',TYPE2,value2,...,0); */
char *	format_replace(char *format,int is_escape,...);

/*****************************************************************************/

#endif	/* #ifndef __FORMAT_H_INCLUDED */

/*****************************************************************************/
                               
