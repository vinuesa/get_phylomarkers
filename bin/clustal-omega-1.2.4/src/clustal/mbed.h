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
 *  RCS $Id: mbed.h 278 2013-05-16 15:53:45Z fabian $
 */

#ifndef CLUSTALO_MBED_H
#define CLUSTALO_MBED_H

#include "muscle_tree.h"
#include "seq.h"


extern int
Mbed(tree_t **tree, mseq_t *prMSeq,
     const int iPairDistType, const char *pcGuidetreeOutfile,
     int iClustersizes, const char *pcClusterFile);

#endif
