/* This a mix of tree functions and data-structures from
 * Bob Edgar's Muscle (version 3.7) ported to pure C.
 *
 * Used files: phy.cpp, tree.h, phytofile.cpp and phyfromclust.cpp
 *
 * Muscle's code is public domain and so is this code here.

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
 *  RCS $Id: muscle_tree.h 230 2011-04-09 15:37:50Z andreas $
 */


#ifndef CLUSTALO_MUSCLE_CLUSTALO_TREE_H
#define CLUSTALO_MUSCLE_CLUSTALO_TREE_H

#include <stdio.h>
#include "util.h"

#ifndef uint
/* limit use of uint (see coding_style_guideline.txt) */
typedef unsigned int uint;
#endif

static const uint NULL_NEIGHBOR = UINT_MAX;

/**
 * @brief guide-tree structure
 *     
 * @note We kept the original variable names here, to make it easy to
 * search through Muscle's source code.
 * From phy.cpp:
 *    Node has 0 to 3 neighbors:
 *    0 neighbors: singleton root
 *    1 neighbor:  leaf, neighbor is parent
 *    2 neigbors:  non-singleton root
 *    3 neighbors: internal node (other than root)
 *    
 *    Minimal rooted tree is single node.
 *    Minimal unrooted tree is single edge.
 *    Leaf node always has nulls in neighbors 2 and 3, neighbor 1 is parent.
 *    When tree is rooted, neighbor 1=parent, 2=left, 3=right.
 *
 */
typedef struct {
    uint m_uNodeCount;/**< number of nodes */
    uint m_uCacheCount;/**< reserved memory */
    
    uint *m_uNeighbor1;/**< parent node */
    uint *m_uNeighbor2;/**< left node */
    uint *m_uNeighbor3;/**< right node */
    
    /* do we have edge lengths info stored (m_dEdgeLength[123]) */
    bool *m_bHasEdgeLength1;
    bool *m_bHasEdgeLength2;
    bool *m_bHasEdgeLength3;

    double *m_dEdgeLength1;
    double *m_dEdgeLength2;
    double *m_dEdgeLength3;
    
#if USE_HEIGHT
    /* unused in our version of the code. we might need it at some
     * stage so keep it in here, but disable via USE_HEIGHT throughout
     * the code */
    double *m_dHeight;
    bool *m_bHasHeight;
#endif

    /**
     * leaf labels.
     * index range: 0 -- (m_uNodeCount+1)/2
     */
    char **m_ptrName;
    
    /**
     * node id.
     * index range: 0 -- m_uNodeCount
     */
    uint *m_Ids;
    
    bool m_bRooted; /**< tree is rooted */
    uint m_uRootNodeIndex;
} tree_t;


extern void
MuscleTreeCreate(tree_t *tree, uint uLeafCount, uint uRoot, const uint *Left,
           const uint  *Right, const float *LeftLength, const float* RightLength,
           const uint *LeafIds, char **LeafNames);

extern void
MuscleTreeToFile(FILE *fp, tree_t *tree);

extern int
MuscleTreeFromFile(tree_t *tree, char *ftree);

extern void
FreeMuscleTree(tree_t *tree);

extern void
LogTree(tree_t *tree, FILE *fp);

extern bool
IsRooted(tree_t *tree);

extern uint
GetNodeCount(tree_t *tree);

extern uint
GetLeafCount(tree_t *tree);
        
extern uint
FirstDepthFirstNode(tree_t *tree);

extern uint
NextDepthFirstNode(uint nodeindex, tree_t *tree);

extern bool
IsLeaf(uint nodeindex, tree_t *tree);

extern void
SetLeafId(tree_t *tree, uint uNodeIndex, uint uId);
    
extern uint
GetLeafId(uint nodeindex, tree_t *tree);

extern char *
GetLeafName(unsigned uNodeIndex, tree_t *tree);

extern uint
GetLeft(uint nodeindex, tree_t *tree);

extern uint
GetRight(uint nodeindex, tree_t *tree);

extern uint
GetRootNodeIndex(tree_t *tree);

extern bool
IsRoot(uint uNodeIndex, tree_t *tree);

extern uint
GetParent(unsigned uNodeIndex, tree_t *tree);

extern double
GetEdgeLength(uint uNodeIndex1, uint uNodeIndex2, tree_t *tree);

extern uint
LeafIndexToNodeIndex(uint uLeafIndex, tree_t *prTree);

extern void
AppendTree(tree_t *prDstTree,
          uint uDstTreeNodeIndex, tree_t *prSrcTree);

extern void
TreeValidate(tree_t *tree);

#endif
