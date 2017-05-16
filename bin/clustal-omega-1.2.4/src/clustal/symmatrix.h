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
 *  RCS $Id: symmatrix.h 288 2013-07-29 13:15:50Z andreas $
 */

/**
 * Functions for symmetric (square) matrices including diagonal.
 *
 * Supports the notion of non-square sub-matrices of a symmetric
 * matrix, i.e. where |rows|<|cols| and the corresponding full matrix
 * would be |cols|x|cols|
 *
 * Instead of making this one big chunk of memory we keep pointers to
 * pointers, so that we can easily realloc (the project where this file
 * originated from needed this for growing a "seed" matrix).
 *
 * FIXME Allocating one big chunk of memory is probably
 * much faster and also easier to maintain.
 * 
 * 
 */

#ifndef CLUSTALO_SYMMATRIX_H
#define CLUSTALO_SYMMATRIX_H

/* ugly, but necessary for cross-checking of names in distmat and
 * mseq */
#include "seq.h"



/**
 * @brief symmetric matrix structure
 */
typedef struct {
    int nrows; /**< number of rows */
    int ncols; /**< number of columns */
    /**
     * stored data
     *
     * @note indices range: [i][j-i] i<=j. use getvalue() and
     * setvalue() instead of accessing directly
     *
     * @see SymMatrixGetValue(), SymMatrixSetValue()
     */
    double **data; 
} symmatrix_t;



extern int
NewSymMatrix(symmatrix_t **symmat, const int nrows, const int ncols);

extern void
SymMatrixSetValue(symmatrix_t *symmat, const int i, const int j, const double value);

extern double
SymMatrixGetValue(symmatrix_t *symmat, const int i, const int j);

extern void
SymMatrixGetValueP(double **value, symmatrix_t *symmat, const int i, const int j);

extern void
FreeSymMatrix(symmatrix_t **symmat);

extern void
SymMatrixPrint(symmatrix_t *symmat, char **labels,  const char *path, bool bPercID);

extern int
SymMatrixRead(char *pcFileIn, symmatrix_t **prSymMat_p, mseq_t *prMSeq);

#endif
