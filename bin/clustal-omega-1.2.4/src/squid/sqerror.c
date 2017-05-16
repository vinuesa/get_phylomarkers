/*****************************************************************
 * SQUID - a library of functions for biological sequence analysis
 * Copyright (C) 1992-2002 Washington University School of Medicine
 * 
 *     This source code is freely distributed under the terms of the
 *     GNU General Public License. See the files COPYRIGHT and LICENSE
 *     for details.
 *****************************************************************/

/* sqerror.c
 * 
 * error handling for the squid library
 * RCS $Id: sqerror.c 217 2011-03-19 10:27:10Z andreas $ (Original squid RCS Id: sqerror.c,v 1.4 1999/07/15 22:28:31 eddy Exp)
 */

				/* a global errno equivalent */
int squid_errno;

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>

#ifdef MEMDEBUG
#include "dbmalloc.h"
#endif

/* Function: Die()
 * 
 * Purpose:  Print an error message and die. The arguments
 *           are formatted exactly like arguments to printf().
 *           
 * Return:   None. Exits the program.
 */          
/* VARARGS0 */
void
Die(char *format, ...)
{
  va_list  argp;
				/* format the error mesg */
  fprintf(stderr, "\nFATAL: ");
  va_start(argp, format);
  vfprintf(stderr, format, argp);
  va_end(argp);
  fprintf(stderr, "\n");
  fflush(stderr);
				/* exit  */
  exit(1);
}



/* Function: Warn()
 * 
 * Purpose:  Print an error message and return. The arguments
 *           are formatted exactly like arguments to printf().
 *           
 * Return:   (void)
 */          
/* VARARGS0 */
void
#ifdef CLUSTALO
Warning(char *format, ...)
#else
Warn(char *format, ...)
#endif
{
  va_list  argp;
				/* format the error mesg */
  fprintf(stderr, "WARNING: ");
  va_start(argp, format);
  vfprintf(stderr, format, argp);
  va_end(argp);
  fprintf(stderr, "\n");
  fflush(stderr);
}

/* Function: Panic()
 * 
 * Purpose:  Die from a lethal error that's not my problem,
 *           but instead a failure of a StdC/POSIX call that
 *           shouldn't fail. Call perror() to get the
 *           errno flag, then die.
 *           
 *           Usually called by the PANIC macro which adds
 *           the __FILE__ and __LINE__ information; see
 *           structs.h.
 *           
 *           Inspired by code in Donald Lewine's book, _POSIX 
 *           Programmer's Guide_.
 */
void
Panic(char *file, int line)
{
  (void) fprintf(stderr, "\nPANIC [%s line %d] ", file, line);
  (void) perror("Unusual error");
  exit(EXIT_FAILURE);
}

