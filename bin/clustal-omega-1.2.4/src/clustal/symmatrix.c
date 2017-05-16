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
 *  RCS $Id: symmatrix.c 303 2016-06-13 13:37:47Z fabian $
 *
 *
 * Functions for symmetric (square) matrices including diagonal.
 * supports the notion of non-square sub-matrices of a symmetric
 * matrix, i.e. where |rows|<|cols|.
 *
 * FIXME Allocating one big chunk of memory is probably
 * much faster and also easier to maintain.
 */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <stdlib.h>
#include <stdio.h>
#include <ctype.h>
#include <string.h>
#include <assert.h>
#include <stdbool.h>
#include "symmatrix.h"

#include "util.h"

#if 0
#define DEBUG
#endif

#define MAX_BUF_SIZE 65536

/**
 * @brief Allocates symmat and its members and initialises them. Data
 * will be calloced, i.e. initialised with zeros.
 *
 * @param[out] symmat
 * newly allocated and initialised symmatrix instance
 * @param[in] nrows
 * number of rows
 * @param[in]
 * ncols number of columns
 *
 * @return: non-zero on error
 *
 * @see FreeSymMatrix()
 *
 * @note: symmat data will be of fake shape nrows x ncols
 *
 */
int
NewSymMatrix(symmatrix_t **symmat, int nrows, int ncols)
{
    int i; /* aux */
    
    assert(nrows>0 && ncols>0 && ncols>=nrows);
    assert(ncols>0 && ncols>=nrows);

    (*symmat) = (symmatrix_t *) malloc(1*sizeof(symmatrix_t));
    if (NULL == (*symmat)) {
        fprintf(stderr, "Couldn't allocate memory (%s|%s)\n",
                __FILE__, __FUNCTION__);
        return -1;
    } 
    else {
        (*symmat)->nrows = nrows;
        (*symmat)->ncols = ncols;
    }
    
    (*symmat)->data = (double **) malloc(nrows * sizeof(double *));
    if (NULL == (*symmat)->data) {
        fprintf(stderr, "Couldn't allocate memory (%s|%s)\n",
                __FILE__, __FUNCTION__);
        free(*symmat);
        *symmat = NULL;
        return -1;
    }

#ifdef HAVE_OPENMP
    /* one malloc for the full matrix data structure
       assumes ncols >= nrows (asserted above) */
    double *mat_data;
    int elements_to_date=0;
    if ((mat_data = (double *)malloc(((nrows+1)*nrows/2 + (ncols-nrows)*nrows) * sizeof(double))) == NULL) {
        fprintf(stderr, "Couldn't allocate MPI memory (%s|%s)\n", __FILE__, __FUNCTION__); 
        free((*symmat)->data);
        free(*symmat);
        *symmat = NULL;
        return -1;
    }
    for (i=0; i<nrows; i++) {
        (*symmat)->data[i] = &mat_data[elements_to_date];
        elements_to_date += ncols-i;
#ifdef TRACE
        fprintf(stderr, "DEBUG(%s|%s():%d) initialising symmat with the number of the beast\n",
                __FILE__, __FUNCTION__, __LINE__);
#endif
        {
            int j;
            for (j=0; j<ncols-i; j++) {
#ifdef TRACE
                (*symmat)->data[i][j] = -666.0;
#else
                (*symmat)->data[i][j] = 0.0;
#endif
            }
        }
    }
#else  /* not HAVE_OPENMP */
    for (i=0; i<nrows; i++) {
        (*symmat)->data[i] = (double *) calloc((ncols-i), sizeof(double));
        if (NULL == (*symmat)->data[i]) {
            fprintf(stderr, "Couldn't allocate memory (%s|%s)\n",
                    __FILE__, __FUNCTION__); 
            while (0!=--i) {
                free((*symmat)->data[i]);
            }           
            free((*symmat)->data);
            free(*symmat);
            *symmat = NULL;
            return -1;
        }
#ifdef TRACE
        fprintf(stderr, "DEBUG(%s|%s():%d) initialising symmat with the number of the beast\n",
                __FILE__, __FUNCTION__, __LINE__);
        {
            int j;
            for (j=0; j<ncols-i; j++) {
                (*symmat)->data[i][j] = -666.0;
            }

        }
#endif
    }
#endif /* HAVE_OPENMP */
    
    return 0;
}
/***   end: new_symmatrix   ***/


/**
 * @brief Sets symmat data of given index to given value
 *
 * @param[in] symmat
 * symmatrix_t whose data is to be set
 * @param[in] i
 * first index 
 * @param[in] j
 * second index 
 * @param[in] value
 * value used to set data point
 *
 * @see SymMatrixGetValue()
 *
 * @note This is a convenience function that checks index order.
 *
 */
void
SymMatrixSetValue(symmatrix_t *symmat, const int i, const int j, const double value)
{
    assert(NULL != symmat);

    if (i<=j) {
        assert(i < symmat->nrows && j < symmat->ncols);
        symmat->data[i][j-i] = value;
    } else {
        assert(j < symmat->nrows && i < symmat->ncols);
        symmat->data[j][i-j] = value;
    }
}
/***   end: symmatrix_setvalue   ***/



/**
 * @brief Returns element of symmat corresponding to given indices
 *
 * @param[in] symmat
 * symmatrix_t of question
 * @param[in] i
 * index i
 * @param[in] j
 * index j
 *
 * @return requested value
 *
 * @see SymMatrixSetValue()
 *
 * @note This is a convenience function that checks index order.
 */
double
SymMatrixGetValue(symmatrix_t *symmat, const int i, const int j)
{
    assert(NULL != symmat);

    if (i<=j) {
        assert(i < symmat->nrows && j < symmat->ncols);
        return symmat->data[i][j-i];
    } else {
        assert(j < symmat->nrows && i < symmat->ncols);
        return symmat->data[j][i-j];
    }
}
/***   end: symmatrix_getvalue   ***/





/**
 * @brief Returns a pointer to an element of symmat corresponding to
 * given indices
 *
 * @param[out] val
 * Value to be set
 * @param[in] symmat
 * symmatrix_t of question
 * @param[in] i
 * index i
 * @param[in] j
 * index j
 *
 * @return pointer to value
 *
 * @see SymMatrixGetValue()
 *
 * @note This is a convenience function that checks index order.
 *
 */
void
SymMatrixGetValueP(double **val, symmatrix_t *symmat,
                   const int i, const int j)
{
    assert(NULL != symmat);
    
    if (i<=j) {
        assert(i < symmat->nrows && j < symmat->ncols);
        *val = &(symmat->data[i][j-i]);
    } else {
        assert(j < symmat->nrows && i < symmat->ncols);
        *val = &(symmat->data[j][i-j]);
    }
}
/***   end: symmatrix_getvaluep   ***/


/**
 * @brief Frees memory allocated by data members of symmat and symmat
 * itself. 
 *
 * @param[in] symmat
 * symmatrix_t to be freed
 *
 * @note Use in conjunction with NewSymMatrix()
 *
 * @see NewSymMatrix()
 *
 */
void
FreeSymMatrix(symmatrix_t **symmat)
{
#ifndef HAVE_OPENMP
    int i;
#endif
    if (NULL != (*symmat)) {
        if (NULL != (*symmat)->data) {
#ifdef HAVE_OPENMP
            free((*symmat)->data[0]);
#else
            for (i=0; i<(*symmat)->nrows; i++) {
                free((*symmat)->data[i]);
            }
#endif
            free((*symmat)->data);
        }
    }
    free(*symmat);
    *symmat = NULL;
}
/***   end: free_symmatrix   ***/



/**
 * @brief Print out a symmat in phylip style. Distances are printed on
 * one line per sequence/object. Since we also support matrices with
 * rows\<cols, the first line can also be nrows by ncolumns instead
 * of just nrows
 *
 * @param[in] symmat
 * the symmatrix_t to print
 * @param[in] labels
 * sequence/objects labels/names. must be at least of length symmat nrows
 * @param[in] path
 * filename or NULL. If NULL stdout will be used.
 *
 */
void
SymMatrixPrint(symmatrix_t *symmat, char **labels,  const char *path, bool bPercID)
{
    FILE *fp = NULL;
    int i, j; /* aux */
    int max_label_len = 0;
    
    if (NULL==symmat || NULL==labels) {
        fprintf(stderr,
                "One of the provided arguments is empty or NULL (print_matrix)\n");
        return;
    }
    if (NULL==path) {
        fp = stdout;
    } else {
        if (NULL==(fp=fopen(path, "w"))) {
            fprintf(stderr, "Couldn't open %s for writing.", path);
            return;
        }
    }

    /* find maximum label length for pretty printing
     */
    for (i=0; i<symmat->nrows; i++) {
        int this_len;
        assert(NULL != labels[i]);
        this_len = strlen(labels[i]);
        if (this_len>max_label_len)
            max_label_len = this_len;
    }

    if (symmat->ncols==symmat->nrows) {
        fprintf(fp, "%u\n", symmat->ncols);
    } else {
        /* this is not standard phylip, but necessary to indicate
           seed matrices */
        fprintf(fp, "%u x %u\n", symmat->nrows, symmat->ncols);
    }
    for (i=0; i<symmat->nrows; i++) {
        /* phylip restriction is 10 characters. we don't care here
         */
        fprintf(fp, "%-*s", max_label_len, labels[i]);
        /* we are lazy and do it like fastdist: write all distances
         *  for one seq to one line
         */
        if (TRUE == bPercID){
            for (j=0; j<symmat->ncols; j++) {
                fprintf(fp, " %f", 100*(1.00-SymMatrixGetValue(symmat, i, j)));
            }
        }
        else{
            for (j=0; j<symmat->ncols; j++) {
                fprintf(fp, " %f", SymMatrixGetValue(symmat, i, j));
            }
        }
        fprintf(fp, "%s", "\n");
    }
    if (NULL != path) {
        (void) fclose(fp);
    } else {
        (void) fflush(fp);
    }
}
/***   end: SymMatrixPrint   ***/



/**
 *
 * @brief Read a distance matrix in phylip format
 *
 * @param[in] pcFileIn
 * distance matrix filename 
 * @param[in] prMSeq
 * optional mseq pointer. only used for consistency checking of names/labels 
 * @param[out] prSymMat_p
 * the symmatrix_t. will be allocated here.
 * @return:
 * non-zero on error
 *
 * @note:
 * FIXME untested code
 */
int
SymMatrixRead(char *pcFileIn, symmatrix_t **prSymMat_p, mseq_t *prMSeq)
{    
    FILE *prFilePointer;
    char *buf;
    /* number of parsed sequences */
    int iNParsedSeq = 0; 
    /* number of parsed distances per sequence */
    int iNParsedDists = 0;
    /* total number of sequences/objects */
    int iTotalNSeq = 0; 
    int iRetCode = 0;

    assert(NULL!=pcFileIn);

    /* could pass down prMSeq for label checking if non NULL */
    if (NULL == (buf = (char *) malloc(MAX_BUF_SIZE * sizeof(char)))) {
        fprintf(stderr, "ERROR: couldn't allocate memory at %s:%s:%d\n",
                __FILE__, __FUNCTION__, __LINE__);
        return -1;
    }
    
    if (NULL == (prFilePointer = fopen(pcFileIn, "r"))) {
        fprintf(stderr, "ERROR: Couldn't open %s for reading\n", pcFileIn);
        free(buf);
        return -1;
    }
    
    /* get number of sequences from first line and allocate memory for
     * distance matrix
     *
     */
    if (NULL == fgets(buf, MAX_BUF_SIZE, prFilePointer) ) {
        fprintf(stderr, "Couldn't read first line from %s\n", pcFileIn);
        iRetCode = -1;
        goto closefile_and_freebuf;
    }
    if (strlen(buf)==MAX_BUF_SIZE-1) {
        fprintf(stderr, "%s\n", "Looks like I couldn't read complete line. Wrong format (or too small MAX_BUF_SIZE)");
        iRetCode = -1;
        goto closefile_and_freebuf;
    }
    if (sscanf(buf, "%d", &iTotalNSeq)!=1) {
        fprintf(stderr, "ERROR: couldn't parse number of sequences from first line of %s\n", pcFileIn); 
        iRetCode = -1;
        goto closefile_and_freebuf;
    }

#if TRACE
    Log(&rLog, LOG_FORCED_DEBUG, "iTotalNSeq parsed from %s is %d\n", pcFileIn, iTotalNSeq);
#endif
    
    if (NewSymMatrix(prSymMat_p, iTotalNSeq, iTotalNSeq)) {
        fprintf(stderr, "FATAL %s", "Memory allocation for distance matrix failed");
        iRetCode = -1;
        goto closefile_and_freebuf;
    }


    /* parse file line by line
     *
     */
    while (NULL != fgets(buf, MAX_BUF_SIZE, prFilePointer)) {
        char *szToken;
        int is_newseq;

        if (MAX_BUF_SIZE-1 == strlen(buf)) {
            fprintf(stderr, "%s\n", "Looks like I couldn't read complete line. Wrong format (or too small MAX_BUF_SIZE)");
            iRetCode = -1;
            goto closefile_and_freebuf;
        }

#ifdef DEBUG
        Log(&rLog, LOG_FORCED_DEBUG, "Got line: %s\n", buf);
#endif

        /* new sequence label at beginning of line?
         */
        if (isblank(buf[0])) {
            is_newseq = 0;
        } else {
            is_newseq = 1;
        }

        /* tokenise line and treat new sequence specially
         */
        szToken = strtok(buf, " \t");
        if (is_newseq==1) {
            iNParsedSeq++;
            iNParsedDists=0;

            /* if format is lower dimensional matrix then first
             * sequence has no distances but might have newline
             * character at it's end.
             */
            while (isspace(szToken[strlen(szToken)-1])) {
                szToken[strlen(szToken)-1]='\0';
            }
            if (strcmp(szToken, prMSeq->sqinfo[iNParsedSeq-1].name) != 0) {
                fprintf(stderr, "Sequence ordering in mseq and distmat differ (expected %s and got %s from distmat %s)n", 
                        prMSeq->sqinfo[iNParsedSeq-1].name, szToken, pcFileIn);
                iRetCode = -1;
                goto closefile_and_freebuf;                
            }
            szToken = strtok(NULL, " \t");
        }
        /* from here on it's only parsing of distances */
        while (szToken != NULL) {
            double dist;
            iNParsedDists++;

            /* only parse what's needed */
            if (iNParsedDists!=iNParsedSeq) {
                /* parse and store distance
                 */
                if (sscanf(szToken, "%lf", &dist)!=1) {
                    fprintf(stderr, "Couldn't parse float from entry '%s'\n", szToken);
                    iRetCode = -1;
                    goto closefile_and_freebuf;
                }
#if TRACE
                Log(&rLog, LOG_FORCED_DEBUG, "Parsed distance %d for seq %d = %f\n", iNParsedDists-1, iNParsedSeq-1, dist);
#endif
                SymMatrixSetValue(*prSymMat_p, iNParsedSeq-1, iNParsedDists-1, dist);
                SymMatrixSetValue(*prSymMat_p, iNParsedDists-1, iNParsedSeq-1, dist);
            }
            szToken = strtok(NULL, " \t");
        }        
    }

    if (iTotalNSeq!=iNParsedSeq) {
        fprintf(stderr, "expected %d seqs, but only parsed %d\n", iTotalNSeq, iNParsedSeq);
        iRetCode = -1;
        goto closefile_and_freebuf;
    }
#if TRACE
    for (i=0; i<iNParsedSeq; i++) {
        int j;
        for (j=0; j<iNParsedSeq; j++) {
            Log(&rLog, LOG_FORCED_DEBUG, "prSymMat_p[%d][%d]=%f\n", i, j, (*prSymMat_p)[i][j]);
        }
    }
#endif

closefile_and_freebuf:
    fclose(prFilePointer);
    free(buf);
    
    return iRetCode;
}
/***   end: SymMatrixRead   ***/

