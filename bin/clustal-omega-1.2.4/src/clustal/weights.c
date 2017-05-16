/* -*- mode: c; tab-width: 4; c-basic-offset: 4;  indent-tabs-mode: nil -*- */

/* Module for deriving sequence weights from a tree. Largely based on
 * Bob Edgar's Muscle (mainly clwwt.cpp; version 3.7). Ported to pure
 * C. Most functions where apparently based on Clustal 1.8 anyway.
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
 * RCS $Id: weights.c 231 2011-04-09 17:13:06Z andreas $
 */

/*
 * Documentation from Muscle
 *
 * """
 *  Compute weights by the CLUSTALW method.
 *  Thompson, Higgins and Gibson (1994), CABIOS (10) 19-29;
 *  see also CLUSTALW paper.
 *   
 *  Weights are computed from the edge lengths of a rooted tree.
 *   
 *  Define the strength of an edge to be its length divided by the number
 *  of leaves under that edge. The weight of a sequence is then the sum
 *  of edge strengths on the path from the root to the leaf.
 *   
 *  Example.
 *   
 *          0.2
 *         -----A     0.1
 *           -x         ------- B     0.7
 *             --------y           ----------- C
 *              0.3     ----------z
 *                      0.4    -------------- D
 *                                   0.8
 *   
 *  Edge    Length  Leaves  Strength
 *  ----    -----   ------  --------
 *  xy              0.3             3               0.1
 *  xA              0.2             1               0.2
 *  yz              0.4             2               0.2
 *  yB              0.1             1               0.1
 *  zC              0.7             1               0.7
 *  zD              0.8             1               0.8
 *   
 *  Leaf    Path            Strengths                       Weight
 *  ----    ----            ---------                       ------
 *  A               xA                      0.2                                     0.2
 *  B               xy-yB           0.1 + 0.1                       0.2
 *  C               xy-yz-zC        0.1 + 0.2 + 0.7         1.0
 *  D               xy-yz-zD        0.1 + 0.2 + 0.8         1.1
 * """
 *
 */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <assert.h>
#include "log.h"
#include "muscle_tree.h"
#include "weights.h"


/* 
   #undef DEBUG 
*/


/**
 * @brief FIXME
 *
 * @param[out] puLeavesUnderNode
 * FIXME
 * @param[in] prTree
 * FIXME
 * @param[in] uNodeIndex
 * FIXME
 *
 * @return The return value
 *
 * @note see Muscle3.7:clwwt.cpp
 * 
 */    
uint
CountLeaves(uint *puLeavesUnderNode, tree_t *prTree, uint uNodeIndex)
{
    uint uLeft;
    uint uRight;
    uint uRightCount;
    uint uLeftCount;
    uint uCount;


    if (IsLeaf(uNodeIndex, prTree)) {
		puLeavesUnderNode[uNodeIndex] = 1;
		return 1;
    }

    uLeft = GetLeft(uNodeIndex, prTree);
    uRight = GetRight(uNodeIndex, prTree);
    uRightCount = CountLeaves(puLeavesUnderNode, prTree, uRight);
    uLeftCount = CountLeaves(puLeavesUnderNode, prTree, uLeft);
    uCount = uRightCount + uLeftCount;

    puLeavesUnderNode[uNodeIndex] = uCount;
        
    return uCount;    
}
/***   end: CountLeaves()   ***/




/**
 * @brief Normalise values in a double array to values between 0 and 1.
 *
 * @param[out] p
 * double array with n elements
 * @param[in] n
 * number of elements in p
 *
 * @note From Muscle3.7: intmath.cpp:Normalize()
 * 
 */    
void
Normalise(double *p, uint n) {
    unsigned i;
    double dSum = 0.0;
    for (i = 0; i < n; ++i) {
        dSum += p[i];
    }
    if (0.0 == dSum) {
        Log(&rLog, LOG_FATAL, "Normalise, sum=0");
    }
    for (i = 0; i < n; ++i) {
        p[i] /= dSum;
    }
}
/***   end: Normalise()   ***/



/**
 * @brief Calculate "Clustal" weights from a tree.
 *
 * FIXME see doc in muscle:clwwt.cpp
 *
 * @param[out] pdWeights_p
 * Will contain a weight for each leaf/sequence. Allocated here. User
 * has to free
 * @param[in] prTree
 * Tree to derive weights from
 *
 * @return 0 on success, non-zero otherwise
 *
 * @note Largely copied from Muscle3.7: clwwt.cpp:CalcClustalWWeights()
 * 
 * @warning FIXME Not sure if Muscle routines are most efficient here.
 * Couldn't we do all this while traversing the tree and thereby safe
 * time?
 *
 */    
int
CalcClustalWeights(double **pdWeights_p, tree_t *prTree)
{
    int i; /* aux */
    uint uLeafCount;
    uint uNodeCount;
	uint *puLeavesUnderNode;
	uint uLeavesUnderRoot;
    uint uRootNodeIndex;
    double *pdStrengths;
    uint uNodeIndex;
    bool bLogWeights = FALSE; /* verbose output of weights */

    
    assert(NULL != pdWeights_p);
    assert(NULL != prTree);

    if (rLog.iLogLevelEnabled <= LOG_DEBUG) {
        bLogWeights = TRUE;
    }
    
    uLeafCount = GetLeafCount(prTree);
    uNodeCount = GetNodeCount(prTree);


    (*pdWeights_p) = (double *) CKMALLOC(uNodeCount * sizeof(double));

    if (0 == uLeafCount) {
		return 0;
    } else if (1 == uLeafCount) {
		(*pdWeights_p)[0] = 1.0;
		return 0;
    } else if (2 == uLeafCount) {
		(*pdWeights_p)[0] = 0.5;
		(*pdWeights_p)[1] = 0.5;
		return 0;
    }
    
    if (!IsRooted(prTree)) {
        Log(&rLog, LOG_ERROR, "Tree must be rooted to get weights");
        CKFREE(pdWeights_p);
        return -1;
    }


#ifdef TRACE
        Log(&rLog, LOG_FORCED_DEBUG, "%s", "Weights follow");
        fprintf(stderr, "Node  Leaves    Length  Strength\n");
        fprintf(stderr, "----  ------  --------  --------\n");
        /*    1234  123456  12345678  12345678 */
#endif
    
    uRootNodeIndex = GetRootNodeIndex(prTree);
    puLeavesUnderNode = (uint *) CKCALLOC(uNodeCount, sizeof(uint));

    uLeavesUnderRoot = CountLeaves(puLeavesUnderNode, prTree, uRootNodeIndex);
	if (uLeavesUnderRoot != uLeafCount) {
		Log(&rLog, LOG_FATAL, "Internal error, root count %u %u",
              uLeavesUnderRoot, uLeafCount);
    }
#if 0
    for (uNodeIndex=0; uNodeIndex<uNodeCount; uNodeIndex++) {
        Log(&rLog, LOG_FORCED_DEBUG, "LeavesUnderNode[%d]=%d", uNodeIndex, puLeavesUnderNode[uNodeIndex]);
    }
#endif

    pdStrengths = (double *) CKMALLOC(uNodeCount * sizeof(double));
                                      
    for (uNodeIndex=0; uNodeIndex < uNodeCount; uNodeIndex++) {
        uint uParent;
        double dLength;
        uint uLeaves;
        double dStrength;
        
        if (IsRoot(uNodeIndex, prTree)) {
            pdStrengths[uNodeIndex] = 0.0;
            continue;
        }
        
        uParent = GetParent(uNodeIndex, prTree);
        dLength = GetEdgeLength(uNodeIndex, uParent, prTree);
        uLeaves = puLeavesUnderNode[uNodeIndex];
        dStrength = dLength / (double) uLeaves;
        pdStrengths[uNodeIndex] = dStrength;
        
#ifdef TRACE
        fprintf(stderr, "%4u  %6u  %8g  %8g\n", uNodeIndex, uLeaves, dLength, dStrength);
#endif        
    }




    
    if (bLogWeights){
        fprintf(stderr, "\n");
        fprintf(stderr, "                 Seq  Path..Weight\n");
        fprintf(stderr, "--------------------  ------------\n");
    }
	for (i=0; i<uLeafCount; i++) {
		double dWeight = 0.0;
		unsigned uLeafNodeIndex;
		unsigned uNode;

        uLeafNodeIndex = LeafIndexToNodeIndex(i, prTree);
        uNode = uLeafNodeIndex;

        if (bLogWeights){
            fprintf(stderr, "%20.20s  %4u ", GetLeafName(uLeafNodeIndex, prTree), uLeafNodeIndex);
        }
		if (! IsLeaf(uLeafNodeIndex, prTree)) {
			Log(&rLog, LOG_FATAL, 
                "Internal error: non-leaf-node %d", uLeafNodeIndex);
        }
            
        /*LOG_DEBUG("dWeight = %f", dWeight);*/
		while (! IsRoot(uNode, prTree)) {
			dWeight += pdStrengths[uNode];
            /*LOG_DEBUG("dWeight +== %f", pdStrengths[uNode]);*/
			uNode = GetParent(uNode, prTree);
            if (bLogWeights){
                fprintf(stderr, "->%u(%g)", uNode, pdStrengths[uNode]);
            }
        }
        /* AW: no idea what this is, but it's done like this in Muscle */
		if (dWeight < 0.0001) {
#ifdef TRACE
			fprintf(stderr, "zero->one");
#endif
			dWeight = 1.0;
        }

        /* @note: the only difference to the muscle code is here: we
         * use the input index for storing weights, instead of the
         * tree leaf index
         */
        (*pdWeights_p)[GetLeafId(uLeafNodeIndex, prTree)] = dWeight;
        if (bLogWeights){
            fprintf(stderr, " = %g\n", dWeight);
        }
    }

#if 0
    for (i=0; i<uLeafCount; i++) {
        Log(&rLog, LOG_FORCED_DEBUG, "Weights before normalisation: pdWeights_p[%d]=%f", i, (*pdWeights_p)[i]);
        /*LOG_DEBUG("Should be %d", GetLeafId(LeafIndexToNodeIndex(i, prTree), prTree));*/
    }
#endif

	Normalise((*pdWeights_p), uLeafCount);
    

    CKFREE(puLeavesUnderNode);
    CKFREE(pdStrengths);
    
    return 0;
}
/***   end: CalcWeights()   ***/
