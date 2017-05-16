/*****************************************************************************/
/* tokenize.h 								     */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/* Standalone library for basic text processing.			     */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/* Copyright (C) 2001, 2004, Pal, A. (apal@szofi.elte.hu). 		     */
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

#ifndef	__TOKENIZE_H_INCLUDED
#define	__TOKENIZE_H_INCLUDED	1

/*****************************************************************************/

/* remove_newlines_and_comments():
   Removes newlines (012 and 015) and comments (all data after a #) from the
   string 'buff'.							     */
void	remove_newlines_and_comments(char *buff);
/* remove_spaces_and_comments():
   Removes whitespaces (011, 012, 015 and 040) and comments from the string 
   'buff'.								     */
void	remove_spaces_and_comments(char *buff);
/* remove_spaces():
   Removes whitespaces (011, 012, 015 and 040) from the string 'buff'.	     */
void	remove_spaces(char *buff);

/* remove_quote():
   Removes all double quotation marks (042) from the string 'buff'.	     */
void	remove_quotes(char *buff);

/* char_is_space():
   Returns a nonzero value if 'c' is a whitespace character, otherwise 0.    */
int	char_is_space(int c);

/* tokenize_spaces():
   Tokenizes the string 'buff', where tokens are delimited by one or more 
   whitespaces. The string 'buff' are modified and the pointers in the 
   array 'tokens' will point to the tokens stored in the newly modified string.
   After each token, a trailing zero (000) are put. The number of tokens are 
   maximalized by the argument 'max'. The function returns the number of 
   tokens and the last element of the array 'tokens' is set to NULL: 
   tokens[max]==NULL (so, note that the array 'tokens' should be declared 
   somehow like char *tokens[max+1]). All parts of buff are treated as one 
   token if it is put between double quotation marks (042) and after the 
   whole process, the quotation marks are removed (from the begining 
   and the end of the tokens).						     */
int	tokenize_spaces(char *buff,char **tokens,int max);
char ** tokenize_spaces_dyn(char *buff);

/* tokenize_char():
   Tokenizes the string 'buff', where tokens are delimited by a single 
   character given in the argument 'tchar'. There's no effect of double-
   quotation marks, spaces or any other special characters. The string 'buff' 
   is modified and the pointers in the array 'tokens' point to the tokens 
   in the newly modified string. After each token, a trailing zero (000) are
   put. The number of tokens are maximalized by the argument 'max'.	     */
int	tokenize_char(char *buff,char **tokens,int tchar,int max);
char ** tokenize_char_dyn_wwt(char *buff,int tchar,int is_terminate);
char ** tokenize_char_dyn(char *buff,int tchar);

/*****************************************************************************/

#endif

/*****************************************************************************/
                                    
