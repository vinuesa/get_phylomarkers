/* This the fast UPGMA algorithm (O(N^2)) as implemented in Bob Edgar's
 * Muscle (UPGMA2.cpp; version 3.7) ported to pure C.
 *
 * Muscle's code is public domain and so is this code here.
 *
 * From http://www.drive5.com/muscle/license.htm:
 * """
 * MUSCLE is public domain software
 *
 * The MUSCLE software, including object and source code and
 * documentation, is hereby donated to the public domain.
 *
 * Disclaimer of warranty
 * THIS SOFTWARE IS PROVIDED "AS IS" WITHOUT WARRANTY OF ANY KIND,
 * EITHER EXPRESS OR IMPLIED, INCLUDING WITHOUT LIMITATION IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * """
 *
 */

/*
 *  RCS $Id: muscle_upgma.h 193 2011-02-07 15:45:21Z andreas $
 */

#ifndef CLUSTALO_UPGMA_H
#define CLUSTALO_UPGMA_H

#include "symmatrix.h"
#include "muscle_tree.h"

enum linkage_e {
    LINKAGE_MIN,
    LINKAGE_AVG,
    LINKAGE_MAX,
    LINKAGE_NEIGHBORJOINING,
    LINKAGE_BIASED
};
typedef enum linkage_e linkage_t;

void MuscleUpgma2(tree_t *tree, symmatrix_t *distmat,
                   linkage_t linkage, char **names);

#endif
