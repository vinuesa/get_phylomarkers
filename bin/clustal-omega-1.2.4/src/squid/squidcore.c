/************************************************************
 * SQUID - a library of functions for biological sequence analysis
 * Copyright (C) 1992-2002 Washington University School of Medicine
 * 
 *     This source code is freely distributed under the terms of the
 *     GNU General Public License. See the files COPYRIGHT and LICENSE
 *     for details.
 ************************************************************/

/* squidcore.c
 * SRE, Sun Jun 20 17:19:04 1999 [Graeme's kitchen]
 * 
 * Core functions for SQUID library.
 * RCS $Id: squidcore.c 217 2011-03-19 10:27:10Z andreas $ (Original squid RCS Id: squidcore.c,v 1.1 1999/07/15 22:24:41 eddy Exp)
 */

#include <stdio.h>
#include "version.h"

/* Function: Banner()
 * Date:     SRE, Sun Jun 20 17:19:41 1999 [Graeme's kitchen]
 *
 * Purpose:  Print a package version and copyright banner.
 *           Used by all the main()'s.
 *           
 *           Expects to be able to pick up defined macros:
 *           macro         example
 *           ------        --------------  
 *           PACKAGE       "HMMER"
 *           RELEASE       "2.0.42"
 *           RELEASEDATE   "April 1 1999"
 *           COPYRIGHT     "Copyright (C) 1992-1999 Washington University School of Medicine"
 *           LICENSE       "HMMER is freely distributed under the GNU General Public License (GPL)."
 *           
 *           This gives us a general mechanism to update release information
 *           without changing multiple points in the code; we can also override
 *           SQUID release data with another package's release data (e.g.
 *           HMMER) just by redefining macros.
 * 
 * Args:     fp     - where to print it
 *           banner - one-line program description, e.g.:
 *                    "foobar - make bars from foo with elan" 
 * Returns:  (void)
 */
void
Banner(FILE *fp, char *banner)
{
  fprintf(fp, "%s\n%s %s (%s)\n%s\n%s\n", banner, PACKAGE, RELEASE, RELEASEDATE, COPYRIGHT, LICENSE);
  fprintf(fp, "- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -\n");
}


