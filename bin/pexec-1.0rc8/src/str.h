/*****************************************************************************/
/* str.h								     */
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

#ifndef	__STR_H_INCLUDED
#define	__STR_H_INCLUDED	1

/*****************************************************************************/

/* strkcpy():
   It works like strncpy() but ensures that 'out' is zero-terminated 
   (therefore, the length of 'out' is always _less_ than 'max').	     */
char *	strkcpy(char *out,char *in,int max);

/* strappend():
   Appends the string 'cat' to the dynamically allocated string 'string'.  
   'string' can be NULL, in this case, the function works like strdup().     */
int	strappend(char **string,char *cat);

/* strappendf():
   It is a combination of realloc(), strcat() and sprintf(). Appends the 
   printf-formatted argument 'format' to the dynamically allocated string
   'string'. It is a bit primitive, therefore use only if long-long strings
   could not occur (otherwise it also works, but can be slow).		     */
int	strappendf(char **string,char *format,...);

/* vstrappendf(): */
int	vstrappendf(char **string,char *format,va_list ap);

/*****************************************************************************/

#endif

/*****************************************************************************/
                  
