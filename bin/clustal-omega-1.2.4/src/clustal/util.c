/* -*- mode: c; tab-width: 4; c-basic-offset: 4;  indent-tabs-mode: nil -*- */

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
 *  RCS $Id: util.c 235 2011-04-13 14:13:19Z andreas $
 */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <time.h>

#include "log.h"
#include "util.h"




/* struct for QSortAndTrackIndex and SortAndTrackIndexCmp[Asc|Desc] */
typedef struct {
    int piIndex;
    int piValue;
} sortwithindex_t;




/**
 * @brief Copy of squid's FileExists(). Copied here to make squid independent.
 */
int
CheckIfFileExists(char *pcFilename)
{
  FILE *prFile;
  if ((prFile = fopen(pcFilename, "r"))) { 
      fclose(prFile);
      return TRUE; 
  }
  return FALSE;
}
/* end of CheckIfFileExists */


/**
 * @brief Allocates memory (malloc)
 *
 * @param[in] bytes
 * bytes to allocated
 * @param[in] function
 * calling function name
 * @param[in] line
 * calling function line
 *
 * @return void pointer to the newly allocated memory
 *
 * @note use provided macro CKMALLOC() which automatically adds
 * function name and line number
 *
 */
void *
CkMalloc(size_t bytes, const char *function, const int line)
{
    void *ret;
    
    if(NULL == (ret = malloc(bytes * sizeof(char)))) {
        Log(&rLog, LOG_FATAL, "Out of memory (requested from %s:%d)\n", function, line);
    }

    return ret; 
}
/***   end: ckmalloc   ***/



/**
 * @brief Allocates memory (calloc). Memory will be
 * set to zero.
 *
 * @param[in] count
 * Allocate space for count objects
 * @param[in] size
 * Objects are of this size
 * @param[in] function
 * calling function name
 * @param[in] line
 * calling function line
 *
 * @return void pointer to the newly allocated and zeroed memory (calloc).
 *
 * @note use provided macro CKCALLOC() which automatically adds
 * function name and line number
 *
 *
 */
void *
CkCalloc(size_t count, size_t size, const char *function, const int line)
{
    void *ret;
    
    if(NULL == (ret = calloc(count, size))) {
        Log(&rLog, LOG_FATAL, "Out of memory (requested from %s:%d)\n",
                function, line);
        exit(EXIT_FAILURE);
    }
    
    return ret; 
}
/***   end: CkCalloc   ***/


/**
* @brief Reallocates memory
 *
 * @param[in] ptr
 * Pointer to memory to be reallocated
 * @param[in] bytes
 * bytes to allocated
 * @param[in] function
 * calling function name
 * @param[in] line
 * calling function line
 *
 * @return void pointer to the newly allocated memory
 *
 * @note use provided macro CKREALLOC() which automatically adds
 * function name and line number   
 *
 */
void *
CkRealloc(void *ptr, size_t bytes, const char *function, const int line)
{
    void *ret=NULL;

    if(NULL == (ret = realloc(ptr, bytes))) {
        Log(&rLog, LOG_FATAL, "FATAL: Out of memory (requested from %s:%d)\n",
                function, line);
    }

    return ret; 
}
/***   end: ckrealloc   ***/



/**
 *
 * @brief Frees memory
 *
 * @param[in] ptr
 * Pointer to memory to be freed
 * @param[in] function
 * calling function name
 * @param[in] line
 * calling function line
 *
 * @return void pointer to the now zeroed memory
 *
 * @note use provided macro CKFREE()
 *
 */
void *
CkFree(void *ptr, const char *function, const int line)
{
    if (ptr == NULL) {
        Log(&rLog, LOG_WARN, "Bad call to CkFree from %s:%d (pointer was NULL)\n", function, line);
    } else {
        free(ptr);
        ptr = NULL;
    }
    return ptr;
}
/***   end: CkFree   ***/




/**
 * @brief safe version of strdup
 *
 * @param[in] src
 * String to copy from
 *
 * @return copy of string
 *
 * @note src is not allowed to be NULL. 
 *
 */
char *
CkStrdup(const char *src)
{
    char *cp;
    assert(NULL!=src);
    
    /*cp = strdup(src); always makes trouble... */
    cp = (char *) CKMALLOC(strlen(src) +1);
    strcpy(cp, src);

    return cp;
}
/***   end: CkStrdup   ***/




/**
 * @brief Creates an int array of size len with len-1 random but unique
 * integers with values 0--len-1. This is "a modified version of
 * Fisher-Yates known as Durstenfeld-Fisher-Yates or
 * Knuth-Fisher-Yates". See
 * http://stackoverflow.com/questions/196017/unique-random-numbers-in-o1
 * 
 * @param[in] perm
 * The permutation array. Has to be preallocated
 * @param[out] len
 * Length of the permutation array
 *
 */    
void
PermutationArray(int **perm, const int len)
{
    int max, i;

    assert(len>0);
    
    srand((unsigned int) time(0));
    (*perm) = (int *) CKMALLOC(len * sizeof(int));
    max = len-1;
    for (i=0; i<len; i++) {
        (*perm)[i] = i;
    }

    while (max>=0) {
        int tmp, randno;
        /* randno = addrand((unsigned long) len); / returns 0--n-1 */
        randno = (rand()%len);
        
        /* swap */
        tmp = (*perm)[randno];
        (*perm)[randno] = (*perm)[max];
        (*perm)[max] = tmp;

        max--;
    }
    return;
}
/***   end: PermutationArray()   ***/



/**
 * @brief Creates an array filled with random, but unique ints of
 * given range. Implementation of the Floyd algorithm. See
 * http://stackoverflow.com/questions/1608181/unique-random-numbers-in-an-integer-array-in-the-c-programming-language
 * Assuming M is the length of the desired array and N is the numeric
 * range: "This approach has the complexity of about O(M) (depends on
 * the search structure used), so it is better suitable when M << N.
 * This approach keeps track of already generated random numbers, so
 * it requires extra memory. However, the beauty of it is that it does
 * not make any of those horrendous "trial and error" iterations,
 * trying to find an unused random number. This algorithm is
 * guaranteed to generate one unique random number after each single
 * call to the random number generator."
 *
 * @warning This won't work if max_value<=array_len. Only use for
 * cases where array_len<<max_value
 *
 * @param[in] array
 * Preallocated int array whose values will be set
 * @param[in] array_len
 * Size of array
 * @param[in] max_value
 * 
 */    
void
RandomUniqueIntArray(int *array, const int array_len, const int max_value)
{
    bool *is_used;
    int in, im;

    assert(array_len<max_value);
    srand((unsigned int) time(0));
    is_used = (bool *) CKCALLOC(max_value, sizeof(bool));
    
    im = 0;

    for (in = max_value - array_len; in < max_value && im < array_len; ++in) {
        int r = rand() % (in + 1);
        if (is_used[r]==TRUE) {
            r = in; /* use 'in' instead of generated number */
        }
        assert(! is_used[r]);
        array[im++] = r;
        is_used[r] = TRUE;
    }
    
    assert(im == array_len);
    
    free(is_used);
    
    return; 
}
/***   end: RandomUniqueIntArray()   ***/



/**
 * @brief int comparison function for qsort
 *
 */
int
IntCmp(const void *a, const void *b)
{
    const int *ia = (const int *)a;
    const int *ib = (const int *)b;
    return *ia  - *ib; 
}
/***   end: IntCmp()   ***/





/**
 * @brief Compare two sortwithindex_t pointers and return the
 * difference between the int value of the 1st sortwithindex_t and the
 * 2nd. Used for ascending sort order in QSortWithIndexes()/
 *
 * @see SortAndTrackIndexCmpDesc() and QSortAndTrackIndex()
 */
int
SortAndTrackIndexCmpAsc(const void *a, const void *b)
{
    const sortwithindex_t *a_t = (const sortwithindex_t *)a;
    const sortwithindex_t *b_t = (const sortwithindex_t *)b;
    
    const int ia = (const int) a_t->piValue;
    const int ib = (const int) b_t->piValue;
    return ia  - ib;     
}
/***   end: SortAndTrackIndexCmpAsc   ***/




/**
 * @brief Compare two sortwithindex_t pointers and return the
 * difference between the int value of the 2nd sortwithindex_t and the
 * 1st. Used for descending sort order in QSortWithIndexes()
 *
 * @see SortAndTrackIndexCmpDesc() and QSortAndTrackIndex()
 */
int
SortAndTrackIndexCmpDesc(const void *a, const void *b)
{
    const sortwithindex_t *a_t = (const sortwithindex_t *)a;
    const sortwithindex_t *b_t = (const sortwithindex_t *)b;
    
    const int ia = (const int) a_t->piValue;
    const int ib = (const int) b_t->piValue;
    return ib - ia;
}
/***   end: SortAndTrackIndexCmpDesc   ***/




/**
 * @brief Sort a given int array in ascending or descending order,
 * while keeping track of the element order.
 *
 * @param[out] piSortedIndices
 * Will contain the indices of the sorted elements. Has to be preallocated.
 * @param[out] piArrayToSort
 * Array with values to sort. Will only be overwritten if
 * bOverwriteArrayToSort it true.
 * @param[in] iArrayLen
 * Number of elements in piArrayToSort.
 * @param[in] cOrder
 * Sort order. 'a' for ascending, 'd' for descending.
 * @param[in] bOverwriteArrayToSort
 * If false do not overwrite the array to sort.
 *
 */    
void
QSortAndTrackIndex(int *piSortedIndices, int *piArrayToSort,
                   const int iArrayLen, const char cOrder,
                   const bool bOverwriteArrayToSort)
{
    int iCtr; /**< aux */
 
    sortwithindex_t *prSort;

    
    assert(NULL!=piSortedIndices);
    assert(iArrayLen>0);
    assert(NULL!=piArrayToSort);
    assert('a'==cOrder || 'd'==cOrder);

    prSort = (sortwithindex_t *) CKMALLOC(iArrayLen * sizeof(sortwithindex_t));
    
    for (iCtr=0; iCtr<iArrayLen; iCtr++) {
        prSort[iCtr].piIndex = iCtr;
        prSort[iCtr].piValue = piArrayToSort[iCtr];
#if 0
        LOG_DEBUG("b4 sort: prSort idx %d = val %d",
                  prSort[iCtr].piIndex, prSort[iCtr].piValue);
#endif
    }

    if ('a'==cOrder) {
        qsort(prSort, iArrayLen, sizeof(sortwithindex_t), SortAndTrackIndexCmpAsc);
    } else if ('d'==cOrder) {
        qsort(prSort, iArrayLen, sizeof(sortwithindex_t), SortAndTrackIndexCmpDesc);
    } else {
        Log(&rLog, LOG_FATAL, "Internal error: unknown order %c", cOrder);
    }

    for (iCtr=0; iCtr<iArrayLen; iCtr++) {
        piSortedIndices[iCtr] = prSort[iCtr].piIndex;
        
        if (bOverwriteArrayToSort) {
            piArrayToSort[iCtr] = prSort[iCtr].piValue;
        }
#if 0
        LOG_DEBUG("after sort: prSort idx %d = val %d",
                  prSort[iCtr].piIndex, prSort[iCtr].piValue);
#endif
    }
    free(prSort);
    
    return; 
}
/***   end: QSortWithIndexes()   ***/




/**
 * @brief Test if file is writable. File may or may not exist.
 *
 * @param[in] pcFileName
 * Filename to check
 *
 * @return True if file is writable at the time of calling
 *
 */    
bool
FileIsWritable(char *pcFileName)
{
    bool bFileAlreadyExisted;
    FILE *prFilePointer;
    bool bIsWritable;
    
    if (0 != CheckIfFileExists(pcFileName)) {
        bFileAlreadyExisted = TRUE;
    } else {
        bFileAlreadyExisted = FALSE;
    }

    if (NULL == (prFilePointer=fopen(pcFileName, "a"))) {
        bIsWritable = FALSE;
    } else {
        bIsWritable = TRUE;
        if (0 != fclose(prFilePointer)) {
            Log(&rLog, LOG_ERROR, "Couldn't close temporily created file %s. Expect trouble...");
        }
    }
    
#if 0
    LOG_DEBUG("%s existed?=%d writable=%d", pcFileName,
              (TRUE==bFileAlreadyExisted? 1 : 0),
              (TRUE==bIsWritable? 1:0));
#endif
    
    /* delete if file didn't exist before and was created here
     * temporarily
     */
    if (FALSE==bFileAlreadyExisted && TRUE==bIsWritable) {
        if (0 != remove(pcFileName)) {
            Log(&rLog, LOG_ERROR, "Removing of temporarily created file %s failed. Expect trouble...");
        }
    }
    return bIsWritable; 
}
/***   end: FileIsWritable()   ***/


