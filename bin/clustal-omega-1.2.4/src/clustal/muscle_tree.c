/* -*- mode: c; tab-width: 4; c-basic-offset: 4;  indent-tabs-mode: nil -*- */

/* This a mix of tree functions and data-structures from
 * Bob Edgar's Muscle (version 3.7) ported to pure C, topped up with
 * some of our own stuff.
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
 *  RCS $Id: muscle_tree.c 230 2011-04-09 15:37:50Z andreas $
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <limits.h>
#include <assert.h>
#include <ctype.h>

#include "util.h"
#include "log.h"
#include "muscle_tree.h"


static const double VERY_NEGATIVE_DOUBLE = -9e29;
/*const double dInsane = VERY_NEGATIVE_DOUBLE;*/
static const double dInsane = -9e29;
static const unsigned uInsane = 8888888;

typedef enum 
{
    NTT_Unknown,

    /* Returned from Tree::GetToken: */
    NTT_Lparen,
    NTT_Rparen,
    NTT_Colon,
    NTT_Comma,
    NTT_Semicolon,
    NTT_String,
    
    /* Following are never returned from Tree::GetToken: */
    NTT_SingleQuotedString,
    NTT_DoubleQuotedString,
    NTT_Comment
} NEWICK_TOKEN_TYPE;


static void
InitCache(uint uCacheCount, tree_t *tree);
static void
TreeZero(tree_t *tree);
static uint
GetNeighborCount(unsigned uNodeIndex, tree_t *tree);
static bool
IsEdge(unsigned uNodeIndex1, unsigned uNodeIndex2, tree_t *tree);
static bool
HasEdgeLength(uint uNodeIndex1, uint uNodeIndex2, tree_t *tree);
static void
TreeToFileNodeRooted(tree_t *tree, uint m_uRootNodeIndex, FILE *fp);
static void
ValidateNode(uint uNodeIndex, tree_t *tree);
static void
AssertAreNeighbors(unsigned uNodeIndex1, unsigned uNodeIndex2, tree_t *tree);
static void
ExpandCache(tree_t *tree);
static void
TreeCreateRooted(tree_t *tree);
static bool
GetGroupFromFile(FILE *fp, uint uNodeIndex, double *ptrdEdgeLength, tree_t *tree);
static NEWICK_TOKEN_TYPE
GetToken(FILE *fp, char szToken[], uint uBytes);
/* stuff from textfile.cpp */
static void
FileSkipWhite(FILE *fp);
static bool
FileSkipWhiteX(FILE *fp);
static void
SetLeafName(uint uNodeIndex, const char *ptrName, tree_t *tree);
uint
AppendBranch(tree_t *tree, uint uExistingLeafIndex);
static void
SetEdgeLength(uint uNodeIndex1, uint uNodeIndex2,
              double dLength, tree_t *tree);
static uint
UnrootFromFile(tree_t *tree);
uint
GetNeighbor(uint uNodeIndex, uint uNeighborSubscript, tree_t *tree);
static void
InitNode(tree_t *prTree, uint uNodeIndex);


/**
 * @param[in] uNodeIndex
 * node to examine
 * @param[in] tree
 * tree to examine
 * @return id of left node
 * @note called GetRight in Muscle3.7
 */
uint
GetLeft(uint uNodeIndex, tree_t *tree) 
{
    assert(NULL != tree);
    assert(tree->m_bRooted && uNodeIndex < tree->m_uNodeCount);
    return tree->m_uNeighbor2[uNodeIndex];
}
/***   end: GetLeft   ***/



/**
 * @param[in] uNodeIndex
 * node to examine
 * @param[in] tree
 * tree to examine
 * @return id of right node
 * @note called GetRight in Muscle3.7
 */
uint
GetRight(uint uNodeIndex, tree_t *tree)
{
    assert(NULL != tree);
    assert(tree->m_bRooted && uNodeIndex < tree->m_uNodeCount);
    return tree->m_uNeighbor3[uNodeIndex];
}
/***   end: GetRight   ***/



/**
 * @param[in] uNodeIndex
 * node to examine
 * @param[in] tree
 * tree to examine
 * @return leaf id of current node
 */
uint
GetLeafId(uint uNodeIndex, tree_t *tree)
{
    assert(NULL != tree);
    assert(uNodeIndex < tree->m_uNodeCount);
    assert(IsLeaf(uNodeIndex, tree));
    
    return tree->m_Ids[uNodeIndex];
}
/***   end: GetLeafId   ***/



/**
 * @note originally called GetLeafName
 *
 */
char *
GetLeafName(unsigned uNodeIndex, tree_t *tree)
{
    assert(NULL != tree);
    assert(uNodeIndex < tree->m_uNodeCount);
    assert(IsLeaf(uNodeIndex, tree));
    return tree->m_ptrName[uNodeIndex];
}
/***   end: GetLeafName   ***/



/**
 * @brief returns first leaf node for a depth-first traversal of tree
 *
 * @param[in] tree
 * tree to traverse
 *
 * @return node index of first leaf node for depth-first traversal
 *
 * @note called FirstDepthFirstNode in Muscle3.7
 *
 */
uint
FirstDepthFirstNode(tree_t *tree)
{
    uint uNodeIndex;
        
    assert(NULL != tree);
    assert(IsRooted(tree));    
    
    /* Descend via left branches until we hit a leaf */
    uNodeIndex =  tree->m_uRootNodeIndex;
    while (!IsLeaf(uNodeIndex, tree))
        uNodeIndex = GetLeft(uNodeIndex, tree);
    return uNodeIndex;
}
/***   end: FirstDepthFirstNode   ***/


/**
 * @brief returns next leaf node index for depth-first traversal of
 * tree
 *
 * @param[in] tree
 * tree to traverse
 * @param[in] uNodeIndex
 * Current node index
 *
 * @return node index of next leaf node during depth-first traversal
 *
 * @note called NextDepthFirstNode in Muscle3.7
 */
uint
NextDepthFirstNode(uint uNodeIndex, tree_t *tree)
{
    uint uParent;
    
    assert(NULL != tree);
    assert(IsRooted(tree));
    assert(uNodeIndex < tree->m_uNodeCount);
    
    if (IsRoot(uNodeIndex, tree))
    {
        /* Node uNodeIndex is root, end of traversal */
        return NULL_NEIGHBOR;
    }
    
    uParent = GetParent(uNodeIndex, tree);
    if (GetRight(uParent, tree) == uNodeIndex) {
        /* Is right branch, return parent=uParent */
        return uParent;
    }
    
    uNodeIndex = GetRight(uParent, tree);
    /* Descend left from right sibling uNodeIndex */
    while (!IsLeaf(uNodeIndex, tree)) {
        uNodeIndex = GetLeft(uNodeIndex, tree);
    }
    
    /*  bottom out at leaf uNodeIndex */
    return uNodeIndex;
}
/***   end: NextDepthFirstNode   ***/



/**
 * @brief check if tree is a rooted tree
 * @param[in] tree
 * tree to check
 * @return TRUE if given tree is rooted, FALSE otherwise
 *
 */
bool
IsRooted(tree_t *tree) 
{
    assert(NULL != tree);
    return tree->m_bRooted;
}
/***   end: IsRooted   ***/


/**
 *
 */
void
FreeMuscleTree(tree_t *tree)
{
    uint i;

    assert(tree!=NULL);
    

    /* FIXME use GetLeafNodeIndex? or
       for (unsigned uNodeIndex = 0; uNodeIndex < m_uNodeCount; ++uNodeIndex)
       {
       if (tree.IsLeaf(uNodeIndex))
       {
       const char *ptrName =
       tree.GetLeafName(uNodeIndex);
    */
    /* IsLeaf needs m_uNodeCount and all m_uNeighbor's
     * so free first
     */
    for (i=0; i<tree->m_uNodeCount; i++) {
        /* IsLeaf needs neighbour count, so free those guys later */
        if (IsLeaf(i, tree)) {
            CKFREE(tree->m_ptrName[i]);
        }
    }
    CKFREE(tree->m_ptrName);

    CKFREE(tree->m_uNeighbor1);
    CKFREE(tree->m_uNeighbor2);
    CKFREE(tree->m_uNeighbor3);
    
    CKFREE(tree->m_Ids);
    
    CKFREE(tree->m_dEdgeLength1);
    CKFREE(tree->m_dEdgeLength2);
    CKFREE(tree->m_dEdgeLength3);
#if USE_HEIGHT
    CKFREE(tree->m_dHeight);
    CKFREE(tree->m_bHasHeight);
#endif 
    CKFREE(tree->m_bHasEdgeLength1);
    CKFREE(tree->m_bHasEdgeLength2);
    CKFREE(tree->m_bHasEdgeLength3);
    
    TreeZero(tree);

    free(tree);
}
/***   end: FreeMuscleTree   ***/



/**
 * @brief create a muscle tree
 *
 * @note Original comment in Muscle code: "Create rooted tree from a
 * vector description. Node indexes are 0..N-1 for leaves, N..2N-2 for
 * internal nodes. Vector subscripts are i-N and have values for
 * internal nodes only, but those values are node indexes 0..2N-2. So
 * e.g. if N=6 and Left[2]=1, this means that the third internal node
 * (node index 8) has the second leaf (node index 1) as its left
 * child. uRoot gives the vector subscript of the root, so add N to
 * get the node index."
 *
 * @param[out] tree
 * newly created tree
 * @param[in] uLeafCount
 * number of leaf nodes
 * @param[in] uRoot
 * Internal node index of root node
 * @param[in] Left
 * Node index of left sibling of an internal node.
 * Index range: 0--uLeafCount-1
 * @param[in] Right
 * Node index of right sibling of an internal node.
 * Index range: 0--uLeafCount-1
 * @param[in] LeftLength
 * Branch length of left branch of an internal node.
 * Index range: 0--uLeafCount-1
 * @param[in] RightLength
 * Branch length of right branch of an internal node.
 * Index range: 0--uLeafCount-1
 * @param[in] LeafIds
 * Leaf id. Index range: 0--uLeafCount
 * @param[in] LeafNames
 * Leaf label. Index range: 0--uLeafCount
 *
 */
void
MuscleTreeCreate(tree_t *tree,
                  uint uLeafCount, uint uRoot,
                  const uint *Left, const uint *Right,
                  const float *LeftLength, const float* RightLength,
                  const uint *LeafIds, char **LeafNames)
{
    uint uNodeIndex;

	TreeZero(tree);
    tree->m_uNodeCount = 2*uLeafCount - 1;
	InitCache(tree->m_uNodeCount, tree);
    
	for (uNodeIndex = 0; uNodeIndex < uLeafCount; ++uNodeIndex) {
		tree->m_Ids[uNodeIndex] = LeafIds[uNodeIndex];
		tree->m_ptrName[uNodeIndex] = CkStrdup(LeafNames[uNodeIndex]);
    }

	for (uNodeIndex = uLeafCount; uNodeIndex < tree->m_uNodeCount; ++uNodeIndex) {
		uint v = uNodeIndex - uLeafCount;
		uint uLeft = Left[v];
		uint uRight = Right[v];
		float fLeft = LeftLength[v];
		float fRight = RightLength[v];
        
		tree->m_uNeighbor2[uNodeIndex] = uLeft;
		tree->m_uNeighbor3[uNodeIndex] = uRight;
        
		tree->m_bHasEdgeLength2[uNodeIndex] = TRUE;
		tree->m_bHasEdgeLength3[uNodeIndex] = TRUE;

		tree->m_dEdgeLength2[uNodeIndex] = fLeft;
		tree->m_dEdgeLength3[uNodeIndex] = fRight;

		tree->m_uNeighbor1[uLeft] = uNodeIndex;
		tree->m_uNeighbor1[uRight] = uNodeIndex;

		tree->m_dEdgeLength1[uLeft] = fLeft;
		tree->m_dEdgeLength1[uRight] = fRight;

		tree->m_bHasEdgeLength1[uLeft] = TRUE;
		tree->m_bHasEdgeLength1[uRight] = TRUE;
    }

	tree->m_bRooted = TRUE;
	tree->m_uRootNodeIndex = uRoot + uLeafCount;
#ifndef NDEBUG
	TreeValidate(tree);
#endif
}
/***   end: MuscleTreeCreate   ***/



/**
 * @param[in] tree
 * tree to write
 * @param[out] fp
 * filepointer to write to2
 * 
 * @brief write a muscle tree to a file in newick format (distances
 * and all names)
 *
 * 
 */
void
MuscleTreeToFile(FILE *fp, tree_t *tree)
{
    assert(NULL != tree);
    if (IsRooted(tree)) {
        TreeToFileNodeRooted(tree, tree->m_uRootNodeIndex, fp);
        fprintf(fp, ";\n");
        return;
    } else {
        Log(&rLog, LOG_FATAL, "FIXME: output of unrooted muscle trees not implemented");
    }
}
/***   end: MuscleTreeToFile   ***/



/**
 * @brief check if given node is a leaf node
 *
 * @param[in] uNodeIndex
 * the node index
 * @param tree
 * the tree
 *
 * @return TRUE if given node is a leaf, FALSE otherwise
 *
 * @note called IsLeaf in Muscle3.7. See tree.h in original code
 */
bool
IsLeaf(uint uNodeIndex, tree_t *tree)
{
    assert(NULL != tree);
    assert(uNodeIndex < tree->m_uNodeCount);
    if (1 == tree->m_uNodeCount)
        return TRUE;
    return 1 == GetNeighborCount(uNodeIndex, tree);
}
/***   end: IsLeaf   ***/



/**
 */
bool
IsRoot(uint uNodeIndex, tree_t *tree)
{
    assert(NULL != tree);
    return IsRooted(tree) && tree->m_uRootNodeIndex == uNodeIndex;
}
/***   end: IsRoot   ***/



/**
 */
uint
GetNeighborCount(uint uNodeIndex, tree_t *tree)
{
    uint n1, n2, n3;
    assert(uNodeIndex < tree->m_uNodeCount);
    assert(NULL != tree);
    assert(NULL != tree->m_uNeighbor1);
    assert(NULL != tree->m_uNeighbor2);
    assert(NULL != tree->m_uNeighbor3);
    n1 = tree->m_uNeighbor1[uNodeIndex];
    n2 = tree->m_uNeighbor2[uNodeIndex];
    n3 = tree->m_uNeighbor3[uNodeIndex];
    return (NULL_NEIGHBOR != n1) + (NULL_NEIGHBOR != n2) + (NULL_NEIGHBOR != n3);
}
/***   end: GetNeighborCount   ***/


/**
 */
uint
GetParent(unsigned uNodeIndex, tree_t *tree)
{
    assert(NULL != tree);
    assert(tree->m_bRooted && uNodeIndex < tree->m_uNodeCount);
    return tree->m_uNeighbor1[uNodeIndex];
}
/***   end: GetParent   ***/



/**
 */
bool
IsEdge(unsigned uNodeIndex1, unsigned uNodeIndex2, tree_t *tree)
{
    assert(uNodeIndex1 < tree->m_uNodeCount && uNodeIndex2 < tree->m_uNodeCount);
    assert(NULL != tree);
    
    return tree->m_uNeighbor1[uNodeIndex1] == uNodeIndex2 ||
        tree->m_uNeighbor2[uNodeIndex1] == uNodeIndex2 ||
        tree->m_uNeighbor3[uNodeIndex1] == uNodeIndex2;
}
/***   end: IsEdge   ***/



/**
 */
bool
HasEdgeLength(uint uNodeIndex1, uint uNodeIndex2, tree_t *tree)
{
    assert(NULL != tree);
    assert(uNodeIndex1 < tree->m_uNodeCount);
    assert(uNodeIndex2 < tree->m_uNodeCount);
    assert(IsEdge(uNodeIndex1, uNodeIndex2, tree));
    
    if (tree->m_uNeighbor1[uNodeIndex1] == uNodeIndex2)
        return tree->m_bHasEdgeLength1[uNodeIndex1];
    else if (tree->m_uNeighbor2[uNodeIndex1] == uNodeIndex2)
        return tree->m_bHasEdgeLength2[uNodeIndex1];
    assert(tree->m_uNeighbor3[uNodeIndex1] == uNodeIndex2);
    return tree->m_bHasEdgeLength3[uNodeIndex1];
}
/***   end:   ***/


/**
 */
double
GetEdgeLength(uint uNodeIndex1, uint uNodeIndex2, tree_t *tree)
{
    assert(NULL != tree);
    assert(uNodeIndex1 < tree->m_uNodeCount && uNodeIndex2 < tree->m_uNodeCount);
    if (!HasEdgeLength(uNodeIndex1, uNodeIndex2, tree))
    {
        Log(&rLog, LOG_FATAL, "Missing edge length in tree %u-%u", uNodeIndex1, uNodeIndex2);
    }
    
    if (tree->m_uNeighbor1[uNodeIndex1] == uNodeIndex2)
        return tree->m_dEdgeLength1[uNodeIndex1];
    else if (tree->m_uNeighbor2[uNodeIndex1] == uNodeIndex2)
        return tree->m_dEdgeLength2[uNodeIndex1];
    assert(tree->m_uNeighbor3[uNodeIndex1] == uNodeIndex2);
    return tree->m_dEdgeLength3[uNodeIndex1];
}
/***   end: GetEdgeLength   ***/


/**
 *
 */
void
InitNode(tree_t *prTree, uint uNodeIndex)
{
    prTree->m_uNeighbor1[uNodeIndex] = NULL_NEIGHBOR;
    prTree->m_uNeighbor2[uNodeIndex] = NULL_NEIGHBOR;
    prTree->m_uNeighbor3[uNodeIndex] = NULL_NEIGHBOR;
    prTree->m_bHasEdgeLength1[uNodeIndex] = FALSE;
    prTree->m_bHasEdgeLength2[uNodeIndex] = FALSE;
    prTree->m_bHasEdgeLength3[uNodeIndex] = FALSE;
#if USE_HEIGHT
    prTree->m_bHasHeight[uNodeIndex] = FALSE;
    prTree->m_dHeight[uNodeIndex] = dInsane;
#endif
    prTree->m_dEdgeLength1[uNodeIndex] = dInsane;
    prTree->m_dEdgeLength2[uNodeIndex] = dInsane;
    prTree->m_dEdgeLength3[uNodeIndex] = dInsane;
    
    prTree->m_ptrName[uNodeIndex] = NULL;
    prTree->m_Ids[uNodeIndex] = uInsane;
}
/***   end: InitNode   ***/


/**
 *
 */
void
InitCache(uint uCacheCount, tree_t *tree)
{
    uint uNodeIndex;
    
    tree->m_uCacheCount = uCacheCount;
    
    tree->m_uNeighbor1 = (uint *) CKMALLOC(sizeof(uint) * tree->m_uCacheCount);
    tree->m_uNeighbor2 = (uint *) CKMALLOC(sizeof(uint) * tree->m_uCacheCount);
    tree->m_uNeighbor3 = (uint *) CKMALLOC(sizeof(uint) * tree->m_uCacheCount);
    
    tree->m_Ids = (uint *) CKMALLOC(sizeof(uint) * tree->m_uCacheCount);
    
    tree->m_dEdgeLength1 = (double *) CKMALLOC(sizeof(double) * tree->m_uCacheCount);
    tree->m_dEdgeLength2 = (double *) CKMALLOC(sizeof(double) * tree->m_uCacheCount);
    tree->m_dEdgeLength3 = (double *) CKMALLOC(sizeof(double) * tree->m_uCacheCount);
#if USE_HEIGHT
    tree->m_dHeight = (double *) CKMALLOC(sizeof(double) * tree->m_uCacheCount);
    tree->m_bHasHeight = (bool *) CKMALLOC(sizeof(bool) * tree->m_uCacheCount);
#endif    
    tree->m_bHasEdgeLength1 = (bool *) CKMALLOC(sizeof(bool) * tree->m_uCacheCount);
    tree->m_bHasEdgeLength2 = (bool *) CKMALLOC(sizeof(bool) * tree->m_uCacheCount);
    tree->m_bHasEdgeLength3 = (bool *) CKMALLOC(sizeof(bool) * tree->m_uCacheCount);

    tree->m_ptrName = (char **) CKMALLOC(sizeof(char *) * tree->m_uCacheCount);
    
    for (uNodeIndex = 0; uNodeIndex < tree->m_uNodeCount; ++uNodeIndex) {
        InitNode(tree, uNodeIndex);
    }
}
/***   end: InitCache   ***/




/**
 *
 * @note Replacing Tree::Clear but no freeing of memory! Just setting
 * everything to 0/NULL
 */
void
TreeZero(tree_t *tree)
{
    assert(NULL != tree);
    tree->m_uNodeCount = 0;
    tree->m_uCacheCount = 0;

    tree->m_uNeighbor1 = NULL;
    tree->m_uNeighbor2 = NULL;
    tree->m_uNeighbor3 = NULL;
    
    tree->m_dEdgeLength1 = NULL;
    tree->m_dEdgeLength2 = NULL;
    tree->m_dEdgeLength3 = NULL;
    
#if USE_HEIGHT
    tree->m_dHeight = NULL;
    tree->m_bHasHeight = NULL;
#endif
    tree->m_bHasEdgeLength1 = NULL;
    tree->m_bHasEdgeLength2 = NULL;
    tree->m_bHasEdgeLength3 = NULL;

    tree->m_ptrName = NULL;
    tree->m_Ids = NULL;
    
    tree->m_bRooted = FALSE;
    tree->m_uRootNodeIndex = 0;
}
/* end: TreeZero */



/**
 *
 */
void
TreeToFileNodeRooted(tree_t *tree, uint uNodeIndex, FILE *fp)
{
    bool bGroup;
    
    assert(IsRooted(tree));
    assert(NULL != tree);
    bGroup = !IsLeaf(uNodeIndex, tree) || IsRoot(uNodeIndex, tree);
    
    if (bGroup) 
        fprintf(fp, "(\n");
    

    if (IsLeaf(uNodeIndex, tree)) {
        /* File.PutString(GetName(uNodeIndex)); */
        fprintf(fp, "%s", tree->m_ptrName[uNodeIndex]);
    } else {
        TreeToFileNodeRooted(tree, GetLeft(uNodeIndex, tree), fp);
        fprintf(fp, ",\n");
        TreeToFileNodeRooted(tree, GetRight(uNodeIndex, tree), fp);
    }

    if (bGroup)
        fprintf(fp, ")");

    if (!IsRoot(uNodeIndex, tree)) {
        uint uParent = GetParent(uNodeIndex, tree);
        if (HasEdgeLength(uNodeIndex, uParent, tree))
            fprintf(fp, ":%g", GetEdgeLength(uNodeIndex, uParent, tree));
    }
    fprintf(fp, "\n");
}
/***   end: TreeToFileNodeRooted   ***/


/**
 *
 *
 */
void
TreeValidate(tree_t *tree)
{
    uint uNodeIndex;
    assert(NULL != tree);
    for (uNodeIndex = 0; uNodeIndex < tree->m_uNodeCount; ++uNodeIndex) {
        ValidateNode(uNodeIndex, tree);
    }
    /* FIXME: maybe set negative length to zero? What impact would
     * that have? */
}
/***   end TreeValidate   ***/



/**
 *
 *
 */
void
ValidateNode(uint uNodeIndex, tree_t *tree)
{
    uint uNeighborCount;
    uint n1, n2, n3;
    assert(NULL != tree);
    
    if (uNodeIndex >= tree->m_uNodeCount)
        Log(&rLog, LOG_FATAL, "ValidateNode(%u), %u nodes", uNodeIndex, tree->m_uNodeCount);

    uNeighborCount = GetNeighborCount(uNodeIndex, tree);
    

    if (2 == uNeighborCount) {
        if (!tree->m_bRooted) {
            Log(&rLog, LOG_FATAL, "Tree::ValidateNode: Node %u has two neighbors, tree is not rooted",
                  uNodeIndex);
        }
        if (uNodeIndex != tree->m_uRootNodeIndex) {
            Log(&rLog, LOG_FATAL, "Tree::ValidateNode: Node %u has two neighbors, but not root node=%u",
                  uNodeIndex, tree->m_uRootNodeIndex);
        }
    }

    n1 = tree->m_uNeighbor1[uNodeIndex];
    n2 = tree->m_uNeighbor2[uNodeIndex];
    n3 = tree->m_uNeighbor3[uNodeIndex];
    
    if (NULL_NEIGHBOR == n2 && NULL_NEIGHBOR != n3) {
        Log(&rLog, LOG_FATAL, "Tree::ValidateNode, n2=null, n3!=null", uNodeIndex);
    }
    if (NULL_NEIGHBOR == n3 && NULL_NEIGHBOR != n2) {
        Log(&rLog, LOG_FATAL, "Tree::ValidateNode, n3=null, n2!=null", uNodeIndex);
    }
    
    if (n1 != NULL_NEIGHBOR)
        AssertAreNeighbors(uNodeIndex, n1, tree);
    if (n2 != NULL_NEIGHBOR)
        AssertAreNeighbors(uNodeIndex, n2, tree);
    if (n3 != NULL_NEIGHBOR)
        AssertAreNeighbors(uNodeIndex, n3, tree);


    
    if (n1 != NULL_NEIGHBOR && (n1 == n2 || n1 == n3)) {
        Log(&rLog, LOG_FATAL, "Tree::ValidateNode, duplicate neighbors in node %u", uNodeIndex);
    }
    if (n2 != NULL_NEIGHBOR && (n2 == n1 || n2 == n3)) {
        Log(&rLog, LOG_FATAL, "Tree::ValidateNode, duplicate neighbors in node %u", uNodeIndex);
    }
    if (n3 != NULL_NEIGHBOR && (n3 == n1 || n3 == n2)) {
        Log(&rLog, LOG_FATAL, "Tree::ValidateNode, duplicate neighbors in node %u", uNodeIndex);
    }


    if (IsRooted(tree)) {
        if (NULL_NEIGHBOR == GetParent(uNodeIndex, tree)) {
            if (uNodeIndex != tree->m_uRootNodeIndex) {
                Log(&rLog, LOG_FATAL, "Tree::ValiateNode(%u), no parent", uNodeIndex);
            }
        } else if (GetLeft(GetParent(uNodeIndex, tree), tree) != uNodeIndex &&
                   GetRight(GetParent(uNodeIndex, tree), tree) != uNodeIndex) {
            Log(&rLog, LOG_FATAL, "Tree::ValidateNode(%u), parent / child mismatch", uNodeIndex);
        }
    }
}
/***   end: ValidateNode   ***/


/**
 *
 *
 */
void
AssertAreNeighbors(unsigned uNodeIndex1, unsigned uNodeIndex2, tree_t *tree)
{
    bool Has12, Has21;
    assert(NULL != tree);

    if (uNodeIndex1 >= tree->m_uNodeCount || uNodeIndex2 >= tree->m_uNodeCount)
        Log(&rLog, LOG_FATAL, "AssertAreNeighbors(%u,%u), are %u nodes",
              uNodeIndex1, uNodeIndex2, tree->m_uNodeCount);

    if (tree->m_uNeighbor1[uNodeIndex1] != uNodeIndex2 &&
        tree->m_uNeighbor2[uNodeIndex1] != uNodeIndex2 &&
        tree->m_uNeighbor3[uNodeIndex1] != uNodeIndex2) {
        Log(&rLog, LOG_FATAL, "AssertAreNeighbors(%u,%u) failed", uNodeIndex1, uNodeIndex2);
    }

    
    if (tree->m_uNeighbor1[uNodeIndex2] != uNodeIndex1 &&
        tree->m_uNeighbor2[uNodeIndex2] != uNodeIndex1 &&
        tree->m_uNeighbor3[uNodeIndex2] != uNodeIndex1) {
        Log(&rLog, LOG_FATAL, "AssertAreNeighbors(%u,%u) failed", uNodeIndex1, uNodeIndex2);
    }


    Has12 = HasEdgeLength(uNodeIndex1, uNodeIndex2, tree);
    Has21 = HasEdgeLength(uNodeIndex2, uNodeIndex1, tree);
    if (Has12 != Has21) {
        HasEdgeLength(uNodeIndex1, uNodeIndex2, tree);
        HasEdgeLength(uNodeIndex2, uNodeIndex1, tree);
        Log(&rLog, LOG_ERROR, "HasEdgeLength(%u, %u)=%c HasEdgeLength(%u, %u)=%c\n",
              uNodeIndex1,
              uNodeIndex2,
              Has12 ? 'T' : 'F',
              uNodeIndex2,
              uNodeIndex1,
              Has21 ? 'T' : 'F');
        Log(&rLog, LOG_FATAL, "Tree::AssertAreNeighbors, HasEdgeLength not symmetric");
    }

    
    if (Has12) {
        double d12 = GetEdgeLength(uNodeIndex1, uNodeIndex2, tree);
        double d21 = GetEdgeLength(uNodeIndex2, uNodeIndex1, tree);
        if (d12 != d21) {
            Log(&rLog, LOG_FATAL, "Tree::AssertAreNeighbors, Edge length disagrees %u-%u=%.3g, %u-%u=%.3g",
                  uNodeIndex1, uNodeIndex2, d12,
                  uNodeIndex2, uNodeIndex1, d21);
        }
    }
}
/***   end: AssertAreNeighbors   ***/



/**
 *
 * @note see phyfromfile.cpp in Muscle3.7. Tree has to be a pointer to
 * an already allocated tree structure.
 *
 * return non-Zero on failure
 *
 * leafids will not be set here (FIXME:CHECK if true)
 */
int
MuscleTreeFromFile(tree_t *tree, char *ftree)
{
    double dEdgeLength;
    bool bEdgeLength;
    char szToken[16];
    NEWICK_TOKEN_TYPE NTT;
    unsigned uThirdNode;
    FILE *fp = NULL;

    assert(tree!=NULL);
    assert(ftree!=NULL);

    if (NULL == (fp=fopen(ftree, "r"))) {
        Log(&rLog, LOG_ERROR, "Couldn't open tree-file '%s' for reading. Skipping", ftree);
        return -1;
    }
    
    /* Assume rooted.
     * If we discover that it is unrooted, will convert on the fly.
     */
    TreeCreateRooted(tree);

    bEdgeLength = GetGroupFromFile(fp, 0, &dEdgeLength, tree);


    /* Next token should be either ';' for rooted tree or ',' for
     * unrooted.
     */
    NTT = GetToken(fp, szToken, sizeof(szToken));

    /* If rooted, all done. */
    if (NTT_Semicolon == NTT) {
        if (bEdgeLength)
            Log(&rLog, LOG_WARN, " *** Warning *** edge length on root group in Newick file %s\n", ftree);
        TreeValidate(tree);
        fclose(fp);
        return 0;
    }

    if (NTT_Comma != NTT)
        Log(&rLog, LOG_FATAL, "Tree::FromFile, expected ';' or ',', got '%s'", szToken);

    uThirdNode = UnrootFromFile(tree);
    bEdgeLength = GetGroupFromFile(fp, uThirdNode, &dEdgeLength, tree);
    if (bEdgeLength)
        SetEdgeLength(0, uThirdNode, dEdgeLength, tree);
    TreeValidate(tree);

    fclose(fp);
    return 0;
}
/***   end MuscleTreeFromFile   ***/



/**
 *
 */
void
ExpandCache(tree_t *tree)
{
    const uint uNodeCount = 100;
    uint uNewCacheCount;
    uint *uNewNeighbor1, *uNewNeighbor2, *uNewNeighbor3;
    uint *uNewIds;
    double *dNewEdgeLength1, *dNewEdgeLength2, *dNewEdgeLength3;
#if USE_HEIGHT
    double *dNewHeight;
    bool *bNewHasHeight;
#endif
    bool *bNewHasEdgeLength1, *bNewHasEdgeLength2, *bNewHasEdgeLength3;
    char **ptrNewName;

    assert(NULL != tree);
    uNewCacheCount = tree->m_uCacheCount + uNodeCount;
    uNewNeighbor1 = (uint *) CKMALLOC(
        uNewCacheCount * sizeof(uint));
    uNewNeighbor2 = (uint *) CKMALLOC(
        uNewCacheCount * sizeof(uint));
    uNewNeighbor3 = (uint *) CKMALLOC(
        uNewCacheCount * sizeof(uint));

    uNewIds = (uint *) CKCALLOC(
        uNewCacheCount, sizeof(uint));
    
    dNewEdgeLength1 = (double *) CKMALLOC(
        uNewCacheCount * sizeof(double));
    dNewEdgeLength2 =  (double *) CKMALLOC(
        uNewCacheCount * sizeof(double));
    dNewEdgeLength3 = (double *) CKMALLOC(
        uNewCacheCount * sizeof(double));
#if USE_HEIGHT
    dNewHeight = (double *) CKMALLOC(
        uNewCacheCount * sizeof(double));
    bNewHasHeight = (bool *) CKMALLOC(
        uNewCacheCount * sizeof(bool));
#endif
     
    bNewHasEdgeLength1 = (bool *) CKMALLOC(
        uNewCacheCount * sizeof(bool));
    bNewHasEdgeLength2 = (bool *) CKMALLOC(
        uNewCacheCount * sizeof(bool));
    bNewHasEdgeLength3 = (bool *) CKMALLOC(
        uNewCacheCount * sizeof(bool));
    ptrNewName = (char **) CKCALLOC(uNewCacheCount, sizeof(char*));

    if (tree->m_uCacheCount > 0) {
        uint uUnsignedBytes, uEdgeBytes;
        uint uBoolBytes, uNameBytes;

        uUnsignedBytes = tree->m_uCacheCount*sizeof(uint);
        
        memcpy(uNewNeighbor1, tree->m_uNeighbor1, uUnsignedBytes);
        memcpy(uNewNeighbor2, tree->m_uNeighbor2, uUnsignedBytes);
        memcpy(uNewNeighbor3, tree->m_uNeighbor3, uUnsignedBytes);

        memcpy(uNewIds, tree->m_Ids, uUnsignedBytes);

        uEdgeBytes = tree->m_uCacheCount*sizeof(double);
        memcpy(dNewEdgeLength1, tree->m_dEdgeLength1, uEdgeBytes);
        memcpy(dNewEdgeLength2, tree->m_dEdgeLength2, uEdgeBytes);
        memcpy(dNewEdgeLength3, tree->m_dEdgeLength3, uEdgeBytes);
#if USE_HEIGHT
        memcpy(dNewHeight, tree->m_dHeight, uEdgeBytes);
#endif        
          
        uBoolBytes = tree->m_uCacheCount*sizeof(bool);
        memcpy(bNewHasEdgeLength1, tree->m_bHasEdgeLength1, uBoolBytes);
        memcpy(bNewHasEdgeLength2, tree->m_bHasEdgeLength2, uBoolBytes);
        memcpy(bNewHasEdgeLength3, tree->m_bHasEdgeLength3, uBoolBytes);
#if USE_HEIGHT
        memcpy(bNewHasHeight, tree->m_bHasHeight, uBoolBytes);
#endif
        uNameBytes = tree->m_uCacheCount*sizeof(char *);
        memcpy(ptrNewName, tree->m_ptrName, uNameBytes);

        /* similiar to FreeMuscleTree
         */
        
        /* IsLeaf needs m_uNodeCount and all m_uNeighbor's
         * so free first
         */
#if 0
        for (i=0; i<tree->m_uNodeCount; i++) {
            if (IsLeaf(i, tree)) {
#ifndef NDEBUG
                if (NULL==tree->m_ptrName[i]) {
                    Log(&rLog, LOG_WARN, "FIXME tree->m_ptrName[%d] is already NULL", i);
                }
#endif
                CKFREE(tree->m_ptrName[i]);
            }
        }
#endif
        CKFREE(tree->m_ptrName);

        CKFREE(tree->m_uNeighbor1);
        CKFREE(tree->m_uNeighbor2);
        CKFREE(tree->m_uNeighbor3);

        CKFREE(tree->m_Ids);

        CKFREE(tree->m_dEdgeLength1);
        CKFREE(tree->m_dEdgeLength2);
        CKFREE(tree->m_dEdgeLength3);

        CKFREE(tree->m_bHasEdgeLength1);
        CKFREE(tree->m_bHasEdgeLength2);
        CKFREE(tree->m_bHasEdgeLength3);
#if USE_HEIGHT
        CKFREE(tree->m_bHasHeight);
        CKFREE(tree->m_dHeight);
#endif
    }
    
    tree->m_uCacheCount = uNewCacheCount;
    tree->m_uNeighbor1 = uNewNeighbor1;
    tree->m_uNeighbor2 = uNewNeighbor2;
    tree->m_uNeighbor3 = uNewNeighbor3;
    tree->m_Ids = uNewIds;
    tree->m_dEdgeLength1 = dNewEdgeLength1;
    tree->m_dEdgeLength2 = dNewEdgeLength2;
    tree->m_dEdgeLength3 = dNewEdgeLength3;
        
#ifdef USE_HEIGHT
    tree->m_dHeight = dNewHeight;
    tree->m_bHasHeight = bNewHasHeight;
#endif
    tree->m_bHasEdgeLength1 = bNewHasEdgeLength1;
    tree->m_bHasEdgeLength2 = bNewHasEdgeLength2;
    tree->m_bHasEdgeLength3 = bNewHasEdgeLength3;
        
    tree->m_ptrName = ptrNewName;

}
/***   end: ExpandCache   ***/



/**
 *
 * Tree must be pointer to an already allocated tree structure
 *
 */
void
TreeCreateRooted(tree_t *tree)
{
    TreeZero(tree);
    ExpandCache(tree);
    tree->m_uNodeCount = 1;

    tree->m_uNeighbor1[0] = NULL_NEIGHBOR;
    tree->m_uNeighbor2[0] = NULL_NEIGHBOR;
    tree->m_uNeighbor3[0] = NULL_NEIGHBOR;
    
    tree->m_bHasEdgeLength1[0] = FALSE;
    tree->m_bHasEdgeLength2[0] = FALSE;
    tree->m_bHasEdgeLength3[0] = FALSE;
#if USE_HEIGHT
    tree->m_bHasHeight[0] = FALSE;
#endif
    
    tree->m_uRootNodeIndex = 0;
    tree->m_bRooted = TRUE;
    
#ifndef NDEBUG
    TreeValidate(tree);
#endif
}
/***   end: TreeCreateRooted   ***/


/**
 *
 */
uint
UnrootFromFile(tree_t *tree)
{
    uint uThirdNode;
#if TRACE
    Log("Before unroot:\n");
    LogMe();
#endif
    
    if (!tree->m_bRooted)
        Log(&rLog, LOG_FATAL, "Tree::Unroot, not rooted");
    
    /* Convention: root node is always node zero */
    assert(IsRoot(0, tree));
    assert(NULL_NEIGHBOR == tree->m_uNeighbor1[0]);

    uThirdNode = tree->m_uNodeCount++;
    
    tree->m_uNeighbor1[0] = uThirdNode;
    tree->m_uNeighbor1[uThirdNode] = 0;
    
    tree->m_uNeighbor2[uThirdNode] = NULL_NEIGHBOR;
    tree->m_uNeighbor3[uThirdNode] = NULL_NEIGHBOR;
    
    tree->m_dEdgeLength1[0] = 0;
    tree->m_dEdgeLength1[uThirdNode] = 0;
    tree->m_bHasEdgeLength1[uThirdNode] = TRUE;
    
    tree->m_bRooted = FALSE;
    
#if TRACE
    Log("After unroot:\n");
    LogMe();
#endif

    return uThirdNode;
}
/***   end: UnrootFromFile   ***/



/**
 *
 *
 */
bool
GetGroupFromFile(FILE *fp, uint uNodeIndex, double *ptrdEdgeLength, tree_t *tree)
{
    char szToken[1024];
    NEWICK_TOKEN_TYPE NTT = GetToken(fp, szToken, sizeof(szToken));
    bool bRightLength;
    bool bEof;
    char c;

    /* Group is either leaf name or (left, right). */
    if (NTT_String == NTT) {
        SetLeafName(uNodeIndex, szToken, tree);
#if TRACE
        Log(&rLog, LOG_FORCED_DEBUG, "Group is leaf '%s'", szToken);
#endif
    } else if (NTT_Lparen == NTT) {
        const unsigned uLeft = AppendBranch(tree, uNodeIndex);
        const unsigned uRight = uLeft + 1;
        double dEdgeLength;
        bool bLeftLength = GetGroupFromFile(fp, uLeft, &dEdgeLength, tree);
        
        /* Left sub-group...
         */
#if TRACE
        Log(&rLog, LOG_FORCED_DEBUG, "%s", "Got '(', group is compound, expect left sub-group");

        if (bLeftLength) {
            Log(&rLog, LOG_FORCED_DEBUG, "Edge length for left sub-group: %.3g", dEdgeLength);
        } else {
            Log(&rLog, LOG_FORCED_DEBUG, "%s", "No edge length for left sub-group");
        }
#endif
        if (bLeftLength)
            SetEdgeLength(uNodeIndex, uLeft, dEdgeLength, tree);

        /* ... then comma ...
         */
#if TRACE
        Log(&rLog, LOG_FORCED_DEBUG, "%s", "Expect comma");
#endif
        NTT = GetToken(fp, szToken, sizeof(szToken));
        if (NTT_Comma != NTT)
            Log(&rLog, LOG_FATAL, "Tree::GetGroupFromFile, expected ',', got '%s'", szToken);
        
        /* ...then right sub-group...
         */
#if TRACE
        Log(&rLog, LOG_FORCED_DEBUG, "%s", "Expect right sub-group");
#endif
        bRightLength = GetGroupFromFile(fp, uRight, &dEdgeLength, tree);
        if (bRightLength)
            SetEdgeLength(uNodeIndex, uRight, dEdgeLength, tree);
        
#if TRACE
        if (bRightLength)
            Log(&rLog, LOG_FORCED_DEBUG, "Edge length for right sub-group: %.3g", dEdgeLength);
        else
            Log(&rLog, LOG_FORCED_DEBUG, "%s", "No edge length for right sub-group");
#endif

        /* ... then closing parenthesis.
         */
#if TRACE
        Log(&rLog, LOG_FORCED_DEBUG, "%s", "Expect closing parenthesis (or comma if > 2-ary)");
#endif
        NTT = GetToken(fp, szToken, sizeof(szToken));
        if (NTT_Rparen == NTT)
            ;
        else if (NTT_Comma == NTT) {
            if (ungetc(',', fp)==EOF)
                Log(&rLog, LOG_FATAL, "%s" "ungetc failed");
            return FALSE;
        } else
            Log(&rLog, LOG_FATAL, "Tree::GetGroupFromFile, expected ')' or ',', got '%s'", szToken);
    } else {
        Log(&rLog, LOG_FATAL, "Tree::GetGroupFromFile, expected '(' or leaf name, got '%s'",
              szToken);
    }
    
    /* Group may optionally be followed by edge length.
     */
    bEof = FileSkipWhiteX(fp);
    if (bEof)
        return FALSE;
    if ((c = fgetc(fp))==EOF) /* GetCharX */
        Log(&rLog, LOG_FATAL, "%s", "fgetc reached end of file");
#if TRACE
    Log(&rLog, LOG_FORCED_DEBUG, "Character following group, could be colon, is '%c'", c);
#endif
    if (':' == c) {
        NTT = GetToken(fp, szToken, sizeof(szToken));
        if (NTT_String != NTT)
            Log(&rLog, LOG_FATAL, "Tree::GetGroupFromFile, expected edge length, got '%s'", szToken);
        *ptrdEdgeLength = atof(szToken);
        return TRUE;
    }
    if (ungetc(c, fp)==EOF)
        Log(&rLog, LOG_FATAL, "%s" "ungetc failed");
    
    return FALSE;
}
/***   end: GetGroupFromFile   ***/




/**
 *
 */
void
FileSkipWhite(FILE *fp)
{
    bool bEof = FileSkipWhiteX(fp);
    if (bEof)
        Log(&rLog, LOG_FATAL, "%s", "End-of-file skipping white space");
}
/***   end: FileSkipWhite   ***/




/**
 *
 */
bool
FileSkipWhiteX(FILE *fp)
{
    for (;;) {
        int c;
        bool bEof;
        
        /* GetChar */
        if ((c = fgetc(fp))==EOF) {
            bEof = TRUE;
        } else {
            bEof = FALSE;
        }

        if (bEof)
            return TRUE;
        if (!isspace(c)) {
            if (ungetc(c, fp)==EOF)
                Log(&rLog, LOG_FATAL, "%s" "ungetc failed");
            break;
        }
    }
    return FALSE;
}
/***   end: FileSkipWhiteX   ***/




/**
 *
 */
NEWICK_TOKEN_TYPE
GetToken(FILE *fp, char szToken[], uint uBytes)
{
    char c;
    unsigned uBytesCopied = 0;
    NEWICK_TOKEN_TYPE TT;

    /* Skip leading white space */
    FileSkipWhite(fp);

    if ((c = fgetc(fp))==EOF) /* GetCharX */
        Log(&rLog, LOG_FATAL, "%s", "fgetc reached end of file");
    
    /* In case a single-character token */
    szToken[0] = c;
    szToken[1] = 0;

    switch (c) {

    case '(':
        return NTT_Lparen;
        
    case ')':
        return NTT_Rparen;
        
    case ':':
        return NTT_Colon;
        
    case ';':
        return NTT_Semicolon;
        
    case ',':
        return NTT_Comma;
        
    case '\'':
        TT = NTT_SingleQuotedString;
        if ((c = fgetc(fp))==EOF) /* GetCharX */
            Log(&rLog, LOG_FATAL, "%s", "fgetc reached end of file");
        break;
        
    case '"':
        TT = NTT_DoubleQuotedString;
        if ((c = fgetc(fp))==EOF) /* GetCharX */
            Log(&rLog, LOG_FATAL, "%s", "fgetc reached end of file");
        break;
        
    case '[':
        TT = NTT_Comment;
        break;
        
    default:
        TT = NTT_String;
        break;
    }
    
    for (;;)
    {
        bool bEof;
        if (TT != NTT_Comment) {
            if (uBytesCopied < uBytes - 2)  {
                szToken[uBytesCopied++] = c;
                szToken[uBytesCopied] = 0;
            } else {
                Log(&rLog, LOG_FATAL, "Tree::GetToken: input buffer too small, token so far='%s'", szToken);
            }
        }
        c = fgetc(fp); /* GetChar */
        bEof = (c==EOF ? TRUE : FALSE);
        if (bEof)
            return TT;
        
        switch (TT) {

        case NTT_String:
            if (0 != strchr("():;,", c)) {
                if (ungetc(c, fp)==EOF)
                    Log(&rLog, LOG_FATAL, "%s" "ungetc failed");
                return NTT_String;
            }
            if (isspace(c))
                return NTT_String;
            break;
            
        case NTT_SingleQuotedString:
            if ('\'' == c)
                return NTT_String;
            break;
            
        case NTT_DoubleQuotedString:
            if ('"' == c)
                return NTT_String;
            break;
            
        case NTT_Comment:
            if (']' == c)
                return GetToken(fp, szToken, uBytes);
            break;
            
        default:
            Log(&rLog, LOG_FATAL, "Tree::GetToken, invalid TT=%u", TT);
        }
    }
}
/***   end: GetToken   ***/



/***   SetLeafName
 *
 */
void
SetLeafName(unsigned uNodeIndex, const char *ptrName, tree_t *tree)
{
    assert(uNodeIndex < tree->m_uNodeCount);
    assert(IsLeaf(uNodeIndex, tree));
    free(tree->m_ptrName[uNodeIndex]);
    /*LOG_DEBUG("Setting tree->m_ptrName[uNodeIndex=%d] to %s", uNodeIndex, ptrName);*/
    tree->m_ptrName[uNodeIndex] = CkStrdup(ptrName);
}
/***   end: SetLeafName   ***/




/**
 *
 * Append a new branch. This adds two new nodes and joins them to an
 * existing leaf node. Return value is k, new nodes have indexes k and
 * k+1 respectively.
 *
 */
uint
AppendBranch(tree_t *tree, uint uExistingLeafIndex)
{
    uint uNewLeaf1;
    uint uNewLeaf2;
    
    assert(tree!=NULL);
    if (0 == tree->m_uNodeCount) {
        Log(&rLog, LOG_FATAL, "%s(): %s", __FUNCTION__, "tree has not been created");
    }
    assert(NULL_NEIGHBOR == tree->m_uNeighbor2[uExistingLeafIndex]);
    assert(NULL_NEIGHBOR == tree->m_uNeighbor3[uExistingLeafIndex]);
    assert(uExistingLeafIndex < tree->m_uNodeCount);
#ifndef NDEBUG
    if (!IsLeaf(uExistingLeafIndex, tree)) {
        Log(&rLog, LOG_FATAL, "AppendBranch(%u): not leaf", uExistingLeafIndex);
    }
#endif

    if (tree->m_uNodeCount >= tree->m_uCacheCount - 2) {
        ExpandCache(tree);
    }
    uNewLeaf1 = tree->m_uNodeCount;
    uNewLeaf2 = tree->m_uNodeCount + 1;

    tree->m_uNodeCount += 2;
    
    tree->m_uNeighbor2[uExistingLeafIndex] = uNewLeaf1;
    tree->m_uNeighbor3[uExistingLeafIndex] = uNewLeaf2;
    
    tree->m_uNeighbor1[uNewLeaf1] = uExistingLeafIndex;
    tree->m_uNeighbor1[uNewLeaf2] = uExistingLeafIndex;
    
    tree->m_uNeighbor2[uNewLeaf1] = NULL_NEIGHBOR;
    tree->m_uNeighbor2[uNewLeaf2] = NULL_NEIGHBOR;
    
    tree->m_uNeighbor3[uNewLeaf1] = NULL_NEIGHBOR;
    tree->m_uNeighbor3[uNewLeaf2] = NULL_NEIGHBOR;
    
    tree->m_dEdgeLength2[uExistingLeafIndex] = 0;
    tree->m_dEdgeLength3[uExistingLeafIndex] = 0;
    
    tree->m_dEdgeLength1[uNewLeaf1] = 0;
    tree->m_dEdgeLength2[uNewLeaf1] = 0;
    tree->m_dEdgeLength3[uNewLeaf1] = 0;
    
    tree->m_dEdgeLength1[uNewLeaf2] = 0;
    tree->m_dEdgeLength2[uNewLeaf2] = 0;
    tree->m_dEdgeLength3[uNewLeaf2] = 0;
    
    tree->m_bHasEdgeLength1[uNewLeaf1] = FALSE;
    tree->m_bHasEdgeLength2[uNewLeaf1] = FALSE;
    tree->m_bHasEdgeLength3[uNewLeaf1] = FALSE;
    
    tree->m_bHasEdgeLength1[uNewLeaf2] = FALSE;
    tree->m_bHasEdgeLength2[uNewLeaf2] = FALSE;
    tree->m_bHasEdgeLength3[uNewLeaf2] = FALSE;
    
#if USE_HEIGHT
    tree->m_bHasHeight[uNewLeaf1] = FALSE;
    tree->m_bHasHeight[uNewLeaf2] = FALSE;
#endif 
    tree->m_Ids[uNewLeaf1] = uInsane;
    tree->m_Ids[uNewLeaf2] = uInsane;
    
    return uNewLeaf1;
}
/***   end: AppendBranch   ***/


/**
 *
 *
 */
void
SetEdgeLength(uint uNodeIndex1, uint uNodeIndex2,
              double dLength, tree_t *tree)
{
    assert(uNodeIndex1 < tree->m_uNodeCount && uNodeIndex2 < tree->m_uNodeCount);
    assert(IsEdge(uNodeIndex1, uNodeIndex2, tree));
    
    if (tree->m_uNeighbor1[uNodeIndex1] == uNodeIndex2) {
        tree->m_dEdgeLength1[uNodeIndex1] = dLength;
        tree->m_bHasEdgeLength1[uNodeIndex1] = TRUE;
    } else if (tree->m_uNeighbor2[uNodeIndex1] == uNodeIndex2) {
        tree->m_dEdgeLength2[uNodeIndex1] = dLength;
        tree->m_bHasEdgeLength2[uNodeIndex1] = TRUE;
        
    } else {
        assert(tree->m_uNeighbor3[uNodeIndex1] == uNodeIndex2);
        tree->m_dEdgeLength3[uNodeIndex1] = dLength;
        tree->m_bHasEdgeLength3[uNodeIndex1] = TRUE;
    }
    
    if (tree->m_uNeighbor1[uNodeIndex2] == uNodeIndex1) {
        tree->m_dEdgeLength1[uNodeIndex2] = dLength;
        tree->m_bHasEdgeLength1[uNodeIndex2] = TRUE;
    } else if (tree->m_uNeighbor2[uNodeIndex2] == uNodeIndex1) {
        tree->m_dEdgeLength2[uNodeIndex2] = dLength;
        tree->m_bHasEdgeLength2[uNodeIndex2] = TRUE;
    } else {
        assert(tree->m_uNeighbor3[uNodeIndex2] == uNodeIndex1);
        tree->m_dEdgeLength3[uNodeIndex2] = dLength;
        tree->m_bHasEdgeLength3[uNodeIndex2] = TRUE;
    }
}
/***   end: SetEdgeLength   ***/



/**
 *
 * Debug output
 *
 * LogMe in phy.cpp
 *
 */
void
LogTree(tree_t *tree, FILE *fp)
{
    uint uNodeIndex;
    uint n1, n2, n3;
    char *ptrName;
    
    fprintf(fp, "This is a tree with %u nodes, which is ", tree->m_uNodeCount);
    
    if (IsRooted(tree)) {
        fprintf(fp, "rooted:\n");
        fprintf(fp, "Index  Parnt  LengthP  Left   LengthL  Right  LengthR     Id  Name\n");
        fprintf(fp, "-----  -----  -------  ----   -------  -----  -------  -----  ----\n");

    } else {
        fprintf(fp, "unrooted;\n");
        fprintf(fp, "Index  Nbr_1  Length1  Nbr_2  Length2  Nbr_3  Length3     Id  Name\n");
        fprintf(fp, "-----  -----  -------  -----  -------  -----  -------  -----  ----\n");
    }
    
    for (uNodeIndex = 0; uNodeIndex < tree->m_uNodeCount; ++uNodeIndex) {
        fprintf(fp, "%5u  ", uNodeIndex);
        n1 = tree->m_uNeighbor1[uNodeIndex];
        n2 = tree->m_uNeighbor2[uNodeIndex];
        n3 = tree->m_uNeighbor3[uNodeIndex];
        
        if (NULL_NEIGHBOR != n1) {
            fprintf(fp, "%5u  ", n1);
            if (tree->m_bHasEdgeLength1[uNodeIndex])
                fprintf(fp, "%7.3g  ", tree->m_dEdgeLength1[uNodeIndex]);
            else
                fprintf(fp, "      *  ");
        } else {
            fprintf(fp, "                ");
        }
        
        if (NULL_NEIGHBOR != n2) {
            fprintf(fp, "%5u  ", n2);
            if (tree->m_bHasEdgeLength2[uNodeIndex])
                fprintf(fp, "%7.3g  ", tree->m_dEdgeLength2[uNodeIndex]);
            else
                fprintf(fp, "      *  ");
        } else {
            fprintf(fp, "                ");
        }
        
        if (NULL_NEIGHBOR != n3) {
            fprintf(fp, "%5u  ", n3);
            if (tree->m_bHasEdgeLength3[uNodeIndex])
                fprintf(fp, "%7.3g  ", tree->m_dEdgeLength3[uNodeIndex]);
            else
                fprintf(fp, "      *  ");
        } else {
            fprintf(fp, "                ");
        }
        
        if (tree->m_Ids != 0 && IsLeaf(uNodeIndex, tree)) {
            unsigned uId = tree->m_Ids[uNodeIndex];
            if (uId == uInsane)
                fprintf(fp, "    *");
            else
                fprintf(fp, "%5u", uId);
        } else {
            fprintf(fp, "     ");
        }
        
        if (tree->m_bRooted && uNodeIndex == tree->m_uRootNodeIndex)
            fprintf(fp, "  [ROOT] ");
        ptrName = tree->m_ptrName[uNodeIndex];
        if (ptrName != 0)
            fprintf(fp, "  %s", ptrName);
        
        fprintf(fp, "\n");
    }     
}
/***   end: LogTree   ***/



/**
 *
 * replaces m_uLeafCount
 */
uint
GetLeafCount(tree_t *tree)
{
    assert(tree!=NULL);
    
    return (tree->m_uNodeCount+1)/2;
}
/***   GetLeafCount   ***/



/**
 *
 */
uint
GetNodeCount(tree_t *tree)
{
    assert(tree!=NULL);
    
    return 2*(GetLeafCount(tree)) - 1;
}
/***   end: GetNodeCount   ***/


/**
 *
 */
uint
GetNeighbor(uint uNodeIndex, uint uNeighborSubscript, tree_t *prTree)
{
    assert(uNodeIndex < prTree->m_uNodeCount);
    switch (uNeighborSubscript)
    {
    case 0:
        return prTree->m_uNeighbor1[uNodeIndex];
    case 1:
        return prTree->m_uNeighbor2[uNodeIndex];
    case 2:
        return prTree->m_uNeighbor3[uNodeIndex];
    }
    Log(&rLog, LOG_FATAL, "Internal error in %s: sub=%u", __FUNCTION__, uNeighborSubscript);
    return NULL_NEIGHBOR;
}
/***   end: GetNeighbor   ***/





/**
 *
 */
void
SetLeafId(tree_t *tree, uint uNodeIndex, uint uId)
{
    assert(uNodeIndex < tree->m_uNodeCount);
    assert(IsLeaf(uNodeIndex, tree));
    tree->m_Ids[uNodeIndex] = uId;
}
/***   end: SetLeafId    ***/


/**
 *
 */
uint
GetRootNodeIndex(tree_t *tree)
{
    assert(NULL!=tree);
    return tree->m_uRootNodeIndex;
}
/***   end: GetRootNodeIndex   ***/



/**
 * @note avoid usage if you want to iterate over all indices, because
 * this will be slow
 *
 */
uint
LeafIndexToNodeIndex(uint uLeafIndex, tree_t *prTree) {
    uint uLeafCount = 0;
    unsigned uNodeCount = GetNodeCount(prTree);
    uint uNodeIndex;
    
    for (uNodeIndex = 0; uNodeIndex < uNodeCount; uNodeIndex++) {
        if (IsLeaf(uNodeIndex, prTree)) {
            if (uLeafCount == uLeafIndex) {
                return uNodeIndex;
            } else {
                uLeafCount++;
            }
        }
    }
    Log(&rLog, LOG_FATAL, "Internal error: node index out of range");
    return 0;
}
/***   end: LeafIndexToNodeIndex   ***/




/**
 * @brief Append a (source) tree to a (dest) tree to a given node
 * which will be replaced. All other nodes in that tree will stay the
 * same.
 *
 * @param[out] prDstTree
 * The tree to append to
 * @param[in] uDstTreeReplaceNodeIndex
 * Dest tree node to which source tree will be appended
 * @param[in] prSrcTree
 * The tree to append
 * 
 * @note No nodes inside prDstTree will change except
 * uDstTreeReplaceNodeIndex
 *
 *
 * @warning: Function won't check or touch the m_Ids/leaf-indices!
 * That means if you want to join two trees with leaf indices 1-10 and
 * 1-10 your m_Ids/leaf-indices won't be unique anymore and the
 * association between your sequences and the tree are broken. Make
 * sure m_Ids are unique before calling me.
 *
 * The easiest would have been to do this by recursively calling
 * AppendBranch() (after adding uSrcTreeNodeIndex as extra argument to
 * this function). But recursion is evil. Yet another version would be
 * to setup all the data and call MuscleTreeCreate() to create a third
 * tree, which seems complicated and wasteful. The approach taken here
 * is the following: increase dest tree memory, copy over each src
 * tree node data and update the indices and counters. This is tricky
 * and has a lot of potential for bugs if tree interface is changed
 * (and even if it isn't).
 *
 */
void
AppendTree(tree_t *prDstTree, uint uDstTreeReplaceNodeIndex, tree_t *prSrcTree)
{
    uint uSrcTreeNodeIndex;
    uint uOrgDstTreeNodeCount;

    assert(NULL!=prDstTree);
    assert(NULL!=prSrcTree);
    assert(uDstTreeReplaceNodeIndex < prDstTree->m_uNodeCount);
    assert(IsLeaf(uDstTreeReplaceNodeIndex, prDstTree));
    assert(IsRooted(prDstTree) && IsRooted(prSrcTree));

    
    uOrgDstTreeNodeCount = prDstTree->m_uNodeCount;

    
    /* increase dest tree memory
     */
    while (prDstTree->m_uCacheCount
           <
           (GetNodeCount(prDstTree) + GetNodeCount(prSrcTree))) {
        ExpandCache(prDstTree);
    }


    /* copy all src tree nodes
     *
     */
    for (uSrcTreeNodeIndex=0;
         uSrcTreeNodeIndex<GetNodeCount(prSrcTree); uSrcTreeNodeIndex++) {
        uint uNewDstNodeIndex = prDstTree->m_uNodeCount;
        
        /* distinguish 4 cases for src nodes to copy:
         *
         * 1. src node is the only node, i.e. root & leaf
         *
         * 2. src node is root: set only left & right, but not parent
         * and just replace the given dest index
         *
         * 3. src node is leaf: set only parent
         *
         * 4. src node is internal node: update all three neighbours
         *
         * FIXME: this is messy. Is there a cleaner way to do this by
         * merging all cases?
         *
         */
        if (IsRoot(uSrcTreeNodeIndex, prSrcTree) && IsLeaf(uSrcTreeNodeIndex, prSrcTree)) {
            /* special case: if this is the only member in
             * tree, i.e. it's root and leaf. Copy leaf name and leaf
             * id. No neighbours to update
             */

            /* free dst node name if set */
            if (NULL != prDstTree->m_ptrName[uDstTreeReplaceNodeIndex]) {
                CKFREE(prDstTree->m_ptrName[uDstTreeReplaceNodeIndex]);
            }

            prDstTree->m_ptrName[uDstTreeReplaceNodeIndex] =
                CkStrdup(GetLeafName(uSrcTreeNodeIndex, prSrcTree));

            prDstTree->m_Ids[uDstTreeReplaceNodeIndex] =
                prSrcTree->m_Ids[uSrcTreeNodeIndex];

            /* no updated of uNodeCount, because we used the replace node */

#if TRACE
            Log(&rLog, LOG_FORCED_DEBUG, "Updated dst rpl node %d with the only src leaf node: parent=%d (%f)",
                      uDstTreeReplaceNodeIndex,
                      prDstTree->m_uNeighbor1[uDstTreeReplaceNodeIndex], prDstTree->m_dEdgeLength1[uDstTreeReplaceNodeIndex]);
#endif
            
        } else if (IsRoot(uSrcTreeNodeIndex, prSrcTree)) {
            /* src node is root: replace uDstTreeReplaceNodeIndex
             * (not uNewDstNodeIndex) with src root, i.e. the
             * uDstTreeReplaceNodeIndex becomes an internal node now.
             *
             * We only have two neighbours 2 & 3 (no parent). Keep old
             * parent info (neighbor 1).
             */

            /* free dst node name if set */
            if (NULL != prDstTree->m_ptrName[uDstTreeReplaceNodeIndex]) {
                CKFREE(prDstTree->m_ptrName[uDstTreeReplaceNodeIndex]);
            }
            
            prDstTree->m_uNeighbor2[uDstTreeReplaceNodeIndex] =
                prSrcTree->m_uNeighbor2[uSrcTreeNodeIndex] + uOrgDstTreeNodeCount;
            prDstTree->m_uNeighbor3[uDstTreeReplaceNodeIndex] =
                prSrcTree->m_uNeighbor3[uSrcTreeNodeIndex] + uOrgDstTreeNodeCount;          
            
            prDstTree->m_bHasEdgeLength2[uDstTreeReplaceNodeIndex] =
                prSrcTree->m_bHasEdgeLength2[uSrcTreeNodeIndex];
            prDstTree->m_bHasEdgeLength3[uDstTreeReplaceNodeIndex] =
                prSrcTree->m_bHasEdgeLength3[uSrcTreeNodeIndex];            

            prDstTree->m_dEdgeLength2[uDstTreeReplaceNodeIndex] =
                prSrcTree->m_dEdgeLength2[uSrcTreeNodeIndex];
            prDstTree->m_dEdgeLength3[uDstTreeReplaceNodeIndex] =
                prSrcTree->m_dEdgeLength3[uSrcTreeNodeIndex];

            /* make Id invalid */
            prDstTree->m_Ids[uDstTreeReplaceNodeIndex] = uInsane;

            /* no updated of uNodeCount, because we used the replace node */

#if TRACE
            Log(&rLog, LOG_FORCED_DEBUG, "Updated dst rpl node %d with the src root node: (untouched) parent=%d (%f) left=%d (%f) right=%d (%f)",
                      uDstTreeReplaceNodeIndex,
                      prDstTree->m_uNeighbor1[uDstTreeReplaceNodeIndex], prDstTree->m_dEdgeLength1[uDstTreeReplaceNodeIndex],
                      prDstTree->m_uNeighbor2[uDstTreeReplaceNodeIndex], prDstTree->m_dEdgeLength2[uDstTreeReplaceNodeIndex],
                      prDstTree->m_uNeighbor3[uDstTreeReplaceNodeIndex], prDstTree->m_dEdgeLength3[uDstTreeReplaceNodeIndex]);
#endif
            
        } else if (IsLeaf(uSrcTreeNodeIndex, prSrcTree)) {
            /* src node is a leaf, which means we only have one
             * neighbour, and that is its parent, i.e. n1
             *
             */

            /* initialise/zero new node to default values
             */
            InitNode(prDstTree, uNewDstNodeIndex);

        
            /*  update m_ptrName/leaf name
             */
            prDstTree->m_ptrName[uNewDstNodeIndex] =
                CkStrdup(GetLeafName(uSrcTreeNodeIndex, prSrcTree));

            /* update parent node (beware of special case: parent was
               src tree root */
            if (IsRoot(prSrcTree->m_uNeighbor1[uSrcTreeNodeIndex], prSrcTree)) {
                    prDstTree->m_uNeighbor1[uNewDstNodeIndex] =
                        uDstTreeReplaceNodeIndex;
            } else {
                prDstTree->m_uNeighbor1[uNewDstNodeIndex] =
                    prSrcTree->m_uNeighbor1[uSrcTreeNodeIndex] + uOrgDstTreeNodeCount;
            }

            /* update edge length info to parent
             */
            prDstTree->m_bHasEdgeLength1[uNewDstNodeIndex] =
                prSrcTree->m_bHasEdgeLength1[uSrcTreeNodeIndex];
            prDstTree->m_dEdgeLength1[uNewDstNodeIndex] =
                prSrcTree->m_dEdgeLength1[uSrcTreeNodeIndex];

            /* update sequence/object id
             */
            prDstTree->m_Ids[uNewDstNodeIndex] =
                prSrcTree->m_Ids[uSrcTreeNodeIndex];            

            /* we used a new node so increase their count */
            prDstTree->m_uNodeCount += 1;
            
#if TRACE
            Log(&rLog, LOG_FORCED_DEBUG, "Updated dst node %d with a src leaf node: parent=%d (%f)",
                      uNewDstNodeIndex,
                      prDstTree->m_uNeighbor1[uNewDstNodeIndex], prDstTree->m_dEdgeLength1[uNewDstNodeIndex]);
#endif
            
        } else  {
            /* src node is not root neither leaf, means we have an
             * internal node. Update all neighbour info
             * 
             */

            /* initialise/zero node values to default values
             */
            InitNode(prDstTree, uNewDstNodeIndex);
          
            /* update neigbours
             */
            /* parent: special case if parent was src tree root */
            if (IsRoot(prSrcTree->m_uNeighbor1[uSrcTreeNodeIndex], prSrcTree)) {
                    prDstTree->m_uNeighbor1[uNewDstNodeIndex] =
                        uDstTreeReplaceNodeIndex;
            } else {
                prDstTree->m_uNeighbor1[uNewDstNodeIndex] =
                    prSrcTree->m_uNeighbor1[uSrcTreeNodeIndex] + uOrgDstTreeNodeCount;
            }
            /* left */
            prDstTree->m_uNeighbor2[uNewDstNodeIndex] =
                prSrcTree->m_uNeighbor2[uSrcTreeNodeIndex] + uOrgDstTreeNodeCount;
            /* right */
            prDstTree->m_uNeighbor3[uNewDstNodeIndex] =
                prSrcTree->m_uNeighbor3[uSrcTreeNodeIndex] + uOrgDstTreeNodeCount;

            /* update edge length info
             */
            /* parent */
            prDstTree->m_bHasEdgeLength1[uNewDstNodeIndex] =
                prSrcTree->m_bHasEdgeLength1[uSrcTreeNodeIndex];
            prDstTree->m_dEdgeLength1[uNewDstNodeIndex] =
                prSrcTree->m_dEdgeLength1[uSrcTreeNodeIndex];
            /* left */
            prDstTree->m_bHasEdgeLength2[uNewDstNodeIndex] =
                prSrcTree->m_bHasEdgeLength2[uSrcTreeNodeIndex];
            prDstTree->m_dEdgeLength2[uNewDstNodeIndex] =
                prSrcTree->m_dEdgeLength2[uSrcTreeNodeIndex];
            /* right */
            prDstTree->m_bHasEdgeLength3[uNewDstNodeIndex] =
                prSrcTree->m_bHasEdgeLength3[uSrcTreeNodeIndex];
            prDstTree->m_dEdgeLength3[uNewDstNodeIndex] =
                prSrcTree->m_dEdgeLength3[uSrcTreeNodeIndex];

            /* we used a new node so increase their count */
            prDstTree->m_uNodeCount += 1;
            
#if TRACE
            Log(&rLog, LOG_FORCED_DEBUG, "Updated dst node %d with an internal src node: parent=%d (%f) left=%d (%f) right=%d (%f)",
                      uNewDstNodeIndex,
                      prDstTree->m_uNeighbor1[uNewDstNodeIndex], prDstTree->m_dEdgeLength1[uNewDstNodeIndex],
                      prDstTree->m_uNeighbor2[uNewDstNodeIndex], prDstTree->m_dEdgeLength2[uNewDstNodeIndex],
                      prDstTree->m_uNeighbor3[uNewDstNodeIndex], prDstTree->m_dEdgeLength3[uNewDstNodeIndex]);
#endif
        }

    }
    /* end for each src tree node */

    
    /*
     * m_uRootNodeIndex stays the same.
     *
     * No need to touch m_uCacheCount.
     *
     */    
#if USE_HEIGHT
    Log(&rLog, LOG_FATAL, "Internal error: Height usage not implemented in %s", __FUNCTION__);
#endif 

    
#ifndef NDEBUG
    TreeValidate(prDstTree);
#endif
        
    return;
}
/***   end: AppendTree()   ***/

