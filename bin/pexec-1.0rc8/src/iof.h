/*****************************************************************************/
/* iof.h 								     */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/* Tools for high-level file manipulation (I/O with f*() functions).	     */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/* Copyright (C) 2004, Pal, A. (apal@szofi.elte.hu). 			     */
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

#ifndef	__IOF_H_INCLUDED
#define	__IOF_H_INCLUDED	1

/* freadline():
   Reads one line (terminated by EOL, 012) from the stream 'fr' and stores it
   in a newly allocated string which is returned by the function. If the
   file stream is over, the function returns NULL. The returned string is 
   zero-terminated and if the line read was not the last one, it includes the 
   trailing EOL character. The string returned can be released with free().  */
char *	freadline(FILE *fr);

/* freadline_bs(): 
   Reads one or more lines from the stream 'fr', and concates them if the 
   ending is a backslash-EOL ("\\\n", 0114 0012) character pair. The function 
   returns a pointer pointing to the newly allocated string with the line read.
   This string is zero-terminated, includes the last trailing EOL character and
   all internal backslash-EOLs are removed. The string returned can be released
   with free().								     */
char *	freadline_bs(FILE *fr);

/* fopenread(), fopenwrite():
   Opens the file 'name' for reading or writing in binary mode. If the name 
   was "-" or NULL, the functions return stdin or stdout respectively.	     */
FILE *	fopenread(char *name);
FILE *	fopenwrite(char *name);

/* fcloseread(), fclosewrite():
   Clses the stream 'fr' if it is not stdin or stdout. Can be used with 
   the functions fopenread()/fopenwrite().				     */
int	fcloseread(FILE *fr);
int	fclosewrite(FILE *fw);

/*****************************************************************************/

#endif

/*****************************************************************************/
                                                                        
