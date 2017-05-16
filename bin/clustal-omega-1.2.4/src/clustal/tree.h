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
 *  RCS $Id: tree.h 193 2011-02-07 15:45:21Z andreas $
 */

#ifndef CLUSTALO_TREE_H
#define CLUSTALO_TREE_H

#include "symmatrix.h"
#include "muscle_tree.h"
#include "seq.h"

enum {LEFT_NODE = 0, RGHT_NODE, PRNT_NODE, DIFF_NODE};

extern void
GuideTreeUpgma(tree_t **tree,
               char **labels, symmatrix_t *tmat, char *ftree);

extern int
GuideTreeFromFile(tree_t **tree,
                  mseq_t *mseq, char *ftree);
    
extern void
TraverseTree(int **piOrderLR_p, 
	      tree_t *tree, mseq_t *mseq);

#endif
