/*********************************************************************
 * Clustal Omega - Multiple sequence alignment
 *
 * Copyright (C) 2010 University College Dublin
 *
 * Clustal-Omega is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License as
 * published by the Free Software Foundation; either version 2 of the
 * License, or (at your option) any later version.
 *
 * This file is part of Clustal-Omega.
 *
 ********************************************************************/

/*
 *  RCS $Id: util.h 230 2011-04-09 15:37:50Z andreas $
 */

#include <limits.h>
#include <string.h>
#include <strings.h>
#include <stdarg.h>
#include <stdbool.h>

#ifndef CLUSTALO_UTIL_H
#define CLUSTALO_UTIL_H


#define CKMALLOC(b) CkMalloc((b), __FUNCTION__, __LINE__)
#define CKCALLOC(c, s) CkCalloc((c), (s), __FUNCTION__, __LINE__)
#define CKREALLOC(p, b) CkRealloc((p), (b), __FUNCTION__, __LINE__)
#define CKFREE(b) ((b)=CkFree((b), __FUNCTION__, __LINE__))

#ifndef MAX
#define MAX(a,b) ((a)>(b)?(a):(b))
#endif
#ifndef MIN
#define MIN(a,b) ((a)<(b)?(a):(b))
#endif

/* STR_EQ: strings are equal, case sensitive */
#define STR_EQ(a,b)         (strcmp((a),(b)) == 0)
/*  STR_NC_EQ: strings are equal, ignoring case */
#define STR_NC_EQ(a,b)      (strcasecmp((a),(b)) == 0)


/* type boolean and false and true defined in stdbool.h */
#ifndef TRUE
#define TRUE true
#endif
#ifndef FALSE
#define FALSE false
#endif

/* clashes with hhalign
#define FAIL -1
#define OK 0
*/



/* don't use the following directly; use macros provided above instead
 */
void *CkMalloc(size_t size, const char *function, const int line);
void *CkCalloc(size_t count, size_t size, const char *function, const int line);
void *CkRealloc(void *ptr, size_t bytes, const char *function, const int line);
void *CkFree(void *ptr, const char *function, const int line);
char *CkStrdup(const char *s);
void PermutationArray(int **array, const int len);
void RandomUniqueIntArray(int *array, const int array_len, const int max_value);
int IntCmp(const void *a, const void *b);
bool FileIsWritable(char *pcFileName);
void QSortAndTrackIndex(int *piSortedIndices, int *piArrayToSort,
                        const int uArrayLen, const char cOrder, const bool bOverwriteArrayToSort);

#endif
