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
 *  RCS $Id: tree.c 278 2013-05-16 15:53:45Z fabian $
 */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <stdlib.h>
#include <string.h>
#include <assert.h>

#include "util.h"
#include "log.h"
#include "muscle_upgma.h"
#include "tree.h"

/**
 *
 * @brief Creates a UPGMA guide tree. This is a frontend function to
 * the ported Muscle UPGMA code ().
 *
 * @param[out] tree
 * created upgma tree. will be allocated here. use FreeMuscleTree()
 * to free
 * @param[in] labels
 * pointer to nseq sequence names
 * @param[in] distmat
 * distance matrix
 * @param[in] ftree
 * optional: if non-NULL, tree will be written to this files
 *
 * @see FreeMuscleTree()
 * @see MuscleUpgma2()
 *
 */
void
GuideTreeUpgma(tree_t **tree, char **labels,
                symmatrix_t *distmat, char *ftree)
{
    linkage_t linkage = LINKAGE_AVG;
    FILE *fp = NULL;

    if (NULL != ftree) {
        if (NULL == (fp=fopen(ftree, "w"))) {
            Log(&rLog, LOG_ERROR, "Couldn't open tree-file '%s' for writing. Skipping", ftree);
        }
        /* fp NULL is handled later */
    }

    (*tree) = (tree_t *) CKMALLOC(1 * sizeof(tree_t));
    MuscleUpgma2((*tree), distmat, linkage, labels);

    if (rLog.iLogLevelEnabled <= LOG_DEBUG) {
        Log(&rLog, LOG_DEBUG, "tree logging...");
        LogTree((*tree), LogGetFP(&rLog, LOG_DEBUG));
    }
    
    if (NULL != fp) {
        MuscleTreeToFile(fp, (*tree));
        Log(&rLog, LOG_INFO, "Guide tree written to %s", ftree);
        fclose(fp);
    }
}
/***   end: guidetree_upgma   ***/



/**
 *
 * @brief
 *
 * @param[out] tree
 * created upgma tree. will be allocated here. use FreeMuscleTree()
 * to free
 * @param[in] mseq
 * @param[in] ftree
 *
 * @return non-zero on error
 *
 */    
int
GuideTreeFromFile(tree_t **tree, mseq_t *mseq, char *ftree)
{
    int iNodeCount;
    int iNodeIndex;
    
    (*tree) = (tree_t *) CKMALLOC(1 * sizeof(tree_t));
    if (MuscleTreeFromFile((*tree), ftree)!=0) {
        Log(&rLog, LOG_ERROR, "%s", "MuscleTreeFromFile failed");
        return -1;
    }

    /* Make sure tree is rooted */
    if (!IsRooted((*tree))) {
        Log(&rLog, LOG_ERROR, "User tree must be rooted");
        return -1;
    }
    
    if ((int)GetLeafCount((*tree)) != mseq->nseqs) {
        Log(&rLog, LOG_ERROR, "User tree does not match input sequences");
        return -1;
    }

    /* compare tree labels and sequence names and set leaf-ids */
    iNodeCount = GetNodeCount((*tree));
    for (iNodeIndex = 0; iNodeIndex < iNodeCount; ++iNodeIndex) {
        char *LeafName;
        int iSeqIndex;
        
        if (!IsLeaf(iNodeIndex, (*tree)))
            continue;
        LeafName = GetLeafName(iNodeIndex, (*tree));

        if ((iSeqIndex=FindSeqName(LeafName, mseq))==-1) {
            Log(&rLog, LOG_ERROR, "Label '%s' in tree could not be found in sequence names", LeafName);
            return -1;
        }
        
        SetLeafId((*tree), iNodeIndex, iSeqIndex);
    }

    if (rLog.iLogLevelEnabled <= LOG_DEBUG) {
        Log(&rLog, LOG_DEBUG, "tree logging...");
        LogTree((*tree),  LogGetFP(&rLog, LOG_DEBUG));
    }
    
    return 0;
}
/***   end: GuideTreeFromFile()   ***/



/**
 *
 * @brief Depth first traversal of tree, i.e. leaf nodes (sequences)
 * will be visited first. Order can be used to guide progressive
 * alignment order.
 * 
 * @param[out] piOrderLR_p
 * order in which left/right nodes (profiles) are to be aligned.
 * allocated here; caller must free.
 * @param[in] tree
 * The tree to traverse; has to be rooted
 * @param[in] mseq
 * corresponding multiple sequence structure
 *
 */    
void
TraverseTree(int **piOrderLR_p, 
              tree_t *tree, mseq_t *mseq)
{
    int tree_nodeindex = 0;
    int order_index = 0;
    int iLeafCount = 0;

    assert(NULL!=tree);
    assert(NULL!=mseq);    
    assert(IsRooted(tree));
    
    /* allocate memory for node/profile alignment order;
     * for every node allocate DIFF_NODE (3) int (1 left, 1 right, 1 parent)
     */
    *piOrderLR_p = (int *)CKCALLOC(DIFF_NODE * GetNodeCount(tree), sizeof(int));
  
    /* Log(&rLog, LOG_FORCED_DEBUG, "print tree->m_iNodeCount=%d", tree->m_iNodeCount); */
  
  
    tree_nodeindex = FirstDepthFirstNode(tree);
    /*LOG_DEBUG("Starting with treenodeindex = %d", tree_nodeindex);*/
  
    order_index = 0;
  
    do {
        if (IsLeaf(tree_nodeindex, tree)) {
            int leafid = GetLeafId(tree_nodeindex, tree);
            if (leafid >= mseq->nseqs){
                Log(&rLog, LOG_FATAL, "Sequence index out of range during tree traversal (leafid=%d nseqs=%d)",
                    leafid, mseq->nseqs);
            }
            if (NULL != mseq->tree_order){
                mseq->tree_order[iLeafCount] = leafid;
                iLeafCount++;
            }

            /* this is a leaf node, 
             * indicate this by registering same leafid for left/right
             */
      
            (*piOrderLR_p)[DIFF_NODE*order_index+LEFT_NODE] = leafid;
            (*piOrderLR_p)[DIFF_NODE*order_index+RGHT_NODE] = leafid;
            (*piOrderLR_p)[DIFF_NODE*order_index+PRNT_NODE] = tree_nodeindex;
      
            Log(&rLog, LOG_DEBUG, "Tree traversal: Visited leaf-node %d (leaf-id %d = Seq '%s')",
                 tree_nodeindex, leafid, mseq->sqinfo[leafid].name);
      
        } else {
            int merge_nodeindex;
            int left;
            int right;
      
            merge_nodeindex = tree_nodeindex;
            left  = GetLeft(tree_nodeindex, tree);
            right = GetRight(tree_nodeindex, tree);
      
            /* this is not a leaf node but a merge node, 
             * register left node (even) and right node (odd)
             */
            (*piOrderLR_p)[DIFF_NODE*order_index+LEFT_NODE] = left;
            (*piOrderLR_p)[DIFF_NODE*order_index+RGHT_NODE] = right;
            (*piOrderLR_p)[DIFF_NODE*order_index+PRNT_NODE] = merge_nodeindex;
      
            Log(&rLog, LOG_DEBUG, "Tree traversal: Visited non-leaf node %d with siblings %d (L) and %d (R)",
                 merge_nodeindex, left, right);
        }
        tree_nodeindex = NextDepthFirstNode(tree_nodeindex, tree);
    
        order_index++;
    
    } while (NULL_NEIGHBOR != tree_nodeindex);
  
    return;
}
/***   end: TraverseTree   ***/
