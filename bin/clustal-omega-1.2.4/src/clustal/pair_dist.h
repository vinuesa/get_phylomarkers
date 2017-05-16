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
 *  RCS $Id: pair_dist.h 302 2016-06-13 13:35:50Z fabian $
 */


#ifndef CLUSTALO_PAIR_DIST_H
#define CLUSTALO_PAIR_DIST_H

#define PAIRDIST_UNKNOWN 0
/* k-tuple distances: Wilbur and Lipman (1983) */
#define PAIRDIST_KTUPLE 1
/* fractional identity between aligned sequences. denominator is
 * minimum seq len (see squid:aligneval.c) */
#define PAIRDIST_SQUIDID 2
/* SQUIDID + Kimura correction */
#define PAIRDIST_SQUIDID_KIMURA 3

#include "seq.h"
#include "symmatrix.h"

extern int
PairDistances(symmatrix_t **distmat, mseq_t *mseq, const int pairdist_type, bool bPercID, 
              const int istart, const unsigned long int iend,
              const int jstart, const unsigned long int jend,
              char *fdist_in, char *fdist_out);

#endif

