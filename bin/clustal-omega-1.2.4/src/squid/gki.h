/*****************************************************************
 * SQUID - a library of functions for biological sequence analysis
 * Copyright (C) 1992-2002 Washington University School of Medicine
 * 
 *     This source code is freely distributed under the terms of the
 *     GNU General Public License. See the files COPYRIGHT and LICENSE
 *     for details.
 *****************************************************************/

#ifndef SQUID_GKI_INCLUDED
#define SQUID_GKI_INCLUDED

/* gki.h
 * SRE, Sat May  1 15:07:22 1999
 * 
 * Declarations of structures, functions for generic key index
 * module: emulation of Perl hashes. See gki.c.
 * 
 * RCS $Id: gki.h 217 2011-03-19 10:27:10Z andreas $ (Original squid RCS Id: gki.h,v 1.2 1999/07/15 22:30:45 eddy Exp)
 */

/* gki_elem:
 *    key, array index pairs are kept in linked list structures.
 */
struct gki_elem {
  char            *key;
  int              idx;
  struct gki_elem *nxt;
};

/* gki:
 *    a dynamically resized hash structure; 
 *    contains a hash table and associated data
 */
typedef struct {
  struct gki_elem **table;
  
  int primelevel;
  int nhash;
  int nkeys;
} GKI;

GKI *GKIInit(void);
void GKIFree(GKI *hash);
int  GKIHashValue(GKI *hash, char *key);
int  GKIStoreKey(GKI *hash, char *key);
int  GKIKeyIndex(GKI *hash, char *key);
void GKIStatus(GKI *hash);

#endif /* SQUID_GKI_INCLUDED */
