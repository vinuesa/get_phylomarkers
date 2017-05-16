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
 *  RCS $Id: ktuple_pair.h 193 2011-02-07 15:45:21Z andreas $
 */

/* K-Tuple code for pairwise alignment (Wilbur and Lipman (1983)
 * Most code taken from showpair (Clustal 1.83)
 */


#ifndef CLUSTALO_KTUPLE_PAIR_H
#define CLUSTALO_KTUPLE_PAIR_H

#include "seq.h"
#include "symmatrix.h"
#include "progress.h"

typedef struct {
    int ktup;
    int window;
    int wind_gap;
    int signif;
} ktuple_param_t;


extern void
KTuplePairDist(symmatrix_t *tmat, mseq_t *mseq,
               int istart, int iend,
               int jstart, int jend,
               ktuple_param_t *aln_param,
               progress_t *prProgress, 
			   unsigned long int *ulStepNo, unsigned long int ulTotalStepNo);

#endif
