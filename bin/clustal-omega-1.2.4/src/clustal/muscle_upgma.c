/* -*- mode: c; tab-width: 4; c-basic-offset: 4;  indent-tabs-mode: nil -*- */

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
 *  RCS $Id: muscle_upgma.c 230 2011-04-09 15:37:50Z andreas $
 *
 *
 * Notes:
 * ------
 * LINKAGE become linkage_t here
 *
 * Replaced the the following member functions for DistCalc DC:
 * DC.GetId = sequence id as int
 * DC.GetName = sequence name
 * DC.GetCount = matrix dim
 * DC.DistRange = vector / matrix row for object i with index j<i
 *
 * Log() has been replaced with Clustal's Info(), Quiet() with Log(&rLog, LOG_FATAL)
 *
 * Made TriangleSubscript() and g_ulTriangleSize ulong to prevent overflow for many sequences
 */

#ifndef ulint
/* limit use of unsigned vars (see coding_style_guideline.txt) */
typedef unsigned long int ulong;
#endif



#include <stdlib.h>
#include <stdio.h>
#include <assert.h>

#include "util.h"
#include "log.h"
#include "symmatrix.h"

#include "muscle_tree.h"
#include "muscle_upgma.h"

/* from distcalc.h */
typedef float dist_t;
static const dist_t BIG_DIST = (dist_t) 1e29;
/* from muscle.h */
static const unsigned uInsane = 8888888;




/*static inline*/
ulong TriangleSubscript(uint uIndex1, uint uIndex2);




#define TRACE   0

#ifndef MIN
#define MIN(x, y)   ((x) < (y) ? (x) : (y))
#endif
#ifndef MIN
#define MAX(x, y)   ((x) > (y) ? (x) : (y))
#endif
#define AVG(x, y)   (((x) + (y))/2)

static uint g_uLeafCount;
static ulong g_ulTriangleSize;
static uint g_uInternalNodeCount;
static uint g_uInternalNodeIndex;

/* Triangular distance matrix is g_Dist, which is allocated
 * as a one-dimensional vector of length g_ulTriangleSize.
 * TriangleSubscript(i,j) maps row,column=i,j to the subscript
 * into this vector.
 * Row / column coordinates are a bit messy.
 * Initially they are leaf indexes 0..N-1.
 * But each time we create a new node (=new cluster, new subtree),
 * we re-use one of the two rows that become available (the children
 * of the new node). This saves memory.
 * We keep track of this through the g_uNodeIndex vector.
 */
static dist_t *g_Dist;

/* Distance to nearest neighbor in row i of distance matrix.
 * Subscript is distance matrix row.
 */
static dist_t *g_MinDist;

/* Nearest neighbor to row i of distance matrix.
 * Subscript is distance matrix row.
 */
static uint *g_uNearestNeighbor;

/* Node index of row i in distance matrix.
 * Node indexes are 0..N-1 for leaves, N..2N-2 for internal nodes.
 * Subscript is distance matrix row.
 */
static uint *g_uNodeIndex;

/* The following vectors are defined on internal nodes,
 * subscripts are internal node index 0..N-2.
 * For g_uLeft/Right, value is the node index 0 .. 2N-2
 * because a child can be internal or leaf.
 */
static uint *g_uLeft;
static uint *g_uRight;
static dist_t *g_Height;
static dist_t *g_LeftLength;
static dist_t *g_RightLength;


/***   CalcDistRange
 *
 * Imitation of DistCalc.DistRange
 *
 * Sets values of row (vector / matrix row) to distances for object i with index j<i
 *
 * row must be preallocated
 */
void CalcDistRange(symmatrix_t *distmat, uint i, dist_t *row)
{
    uint j;
    for (j = 0; j < i; ++j) {
        row[j] = SymMatrixGetValue(distmat, i, j);
    }
}
/* end of CalcDistRange */



/*static inline*/
ulong
TriangleSubscript(uint uIndex1, uint uIndex2)
{
    ulong v;
#ifndef NDEBUG
    if (uIndex1 >= g_uLeafCount || uIndex2 >= g_uLeafCount)
        Log(&rLog, LOG_FATAL, "TriangleSubscript(%u,%u) %u", uIndex1, uIndex2, g_uLeafCount);
#endif
    if (uIndex1 >= uIndex2)
        v = uIndex2 + (uIndex1*(uIndex1 - 1))/2;
    else
        v = uIndex1 + (uIndex2*(uIndex2 - 1))/2;
    assert(v < (g_uLeafCount*(g_uLeafCount - 1))/2);
    return v;
}

#ifdef UNUSED
static void ListState()
{
    uint i, j;
    Info("Dist matrix\n");
    Info("     ");
    for (i = 0; i < g_uLeafCount; ++i)
    {
        if (uInsane == g_uNodeIndex[i])
            continue;
        Info("  %5u", g_uNodeIndex[i]);
    }
    Info("\n");

    for (i = 0; i < g_uLeafCount; ++i)
    {
        if (uInsane == g_uNodeIndex[i])
            continue;
        Info("%5u  ", g_uNodeIndex[i]);
        for (j = 0; j < g_uLeafCount; ++j)
        {
            if (uInsane == g_uNodeIndex[j])
                continue;
            if (i == j)
                Info("       ");
            else
            {
                ulong v = TriangleSubscript(i, j);
                Info("%5.2g  ", g_Dist[v]);
            }
        }
        Info("\n");
    }

    Info("\n");
    Info("    i   Node   NrNb      Dist\n");
    Info("-----  -----  -----  --------\n");
    for (i = 0; i < g_uLeafCount; ++i)
    {
        if (uInsane == g_uNodeIndex[i])
            continue;
        Info("%5u  %5u  %5u  %8.3f\n",
             i,
             g_uNodeIndex[i],
             g_uNearestNeighbor[i],
             g_MinDist[i]);
    }

    Info("\n");
    Info(" Node      L      R  Height  LLength  RLength\n");
    Info("-----  -----  -----  ------  -------  -------\n");
    for (i = 0; i <= g_uInternalNodeIndex; ++i)
        Info("%5u  %5u  %5u  %6.2g  %6.2g  %6.2g\n",
             i,
             g_uLeft[i],
             g_uRight[i],
             g_Height[i],
             g_LeftLength[i],
             g_RightLength[i]);
}
#endif
/* ifdef UNUSED */

/**
 * @brief Creates a UPGMA in O(N^2) tree from given distmat
 *
 * @param[out] tree
 * newly created rooted UPGMA tree
 * @param[in] distmat
 * distance matrix to be clustered
 * @param[in] linkage
 * linkage type
 * @param[in] names
 * leaf names, will be copied
 *
 * @note called UPGMA2() in Muscle3.7.
 * caller has to free with FreeMuscleTree()
 *
 * @see FreeMuscleTree()
 */
void
MuscleUpgma2(tree_t *tree, symmatrix_t *distmat, linkage_t linkage, char **names)
{
    int i, j;
    uint *Ids;

    /* only works on full symmetric matrices */
    assert (distmat->nrows==distmat->ncols);
   
    g_uLeafCount = distmat->ncols;    
    g_ulTriangleSize = (g_uLeafCount*(g_uLeafCount - 1))/2;
    g_uInternalNodeCount = g_uLeafCount - 1;

    g_Dist = (dist_t *) CKMALLOC(g_ulTriangleSize * sizeof(dist_t));

    g_uNodeIndex = (uint*) CKMALLOC(sizeof(uint) * g_uLeafCount);
    g_uNearestNeighbor = (uint*) CKMALLOC(sizeof(uint) * g_uLeafCount);
    g_MinDist = (dist_t *) CKMALLOC(sizeof(dist_t) * g_uLeafCount);
    Ids = (uint*) CKMALLOC(sizeof(uint) * g_uLeafCount);
    /* NOTE: we replaced Names with argument names */

    /**
     * left and right node indices, as well as left and right
     * branch-length and height for for internal nodes
     */
    g_uLeft =  (uint*) CKMALLOC(sizeof(uint) * g_uInternalNodeCount);
    g_uRight =  (uint*) CKMALLOC(sizeof(uint) * g_uInternalNodeCount);
    g_Height =  (dist_t*) CKMALLOC(sizeof(dist_t) * g_uInternalNodeCount);
    g_LeftLength =  (dist_t*) CKMALLOC(sizeof(dist_t) * g_uInternalNodeCount);
    g_RightLength =  (dist_t*) CKMALLOC(sizeof(dist_t) * g_uInternalNodeCount);
    
    for (i = 0; i < g_uLeafCount; ++i) {
        g_MinDist[i] = BIG_DIST;
        g_uNodeIndex[i] = i;
        g_uNearestNeighbor[i] = uInsane;
        Ids[i] = i;
    }
    
    for (i = 0; i < g_uInternalNodeCount; ++i) {
        g_uLeft[i] = uInsane;
        g_uRight[i] = uInsane;
        g_LeftLength[i] = BIG_DIST;
        g_RightLength[i] = BIG_DIST;
        g_Height[i] = BIG_DIST;
    }
    
/* Compute initial NxN triangular distance matrix.
 * Store minimum distance for each full (not triangular) row.
 * Loop from 1, not 0, because "row" is 0, 1 ... i-1,
 * so nothing to do when i=0.
 */
    for (i = 1; i < g_uLeafCount; ++i) {
        dist_t *Row = g_Dist + TriangleSubscript(i, 0);
        CalcDistRange(distmat, i, Row);
        for (j = 0; j < i; ++j) {
            const dist_t d = Row[j];
            if (d < g_MinDist[i]) {
                g_MinDist[i] = d;
                g_uNearestNeighbor[i] = j;
            }
            if (d < g_MinDist[j]) {
                g_MinDist[j] = d;
                g_uNearestNeighbor[j] = i;
            }
        }
    }

#if TRACE
    Info("Initial state:\n");
    ListState();
#endif

    for (g_uInternalNodeIndex = 0;
         g_uInternalNodeIndex < g_uLeafCount - 1;
         ++g_uInternalNodeIndex) {

        dist_t dtNewMinDist = BIG_DIST;
        uint uNewNearestNeighbor = uInsane;

#if TRACE
        Info("\n");
        Info("Internal node index %5u\n", g_uInternalNodeIndex);
        Info("-------------------------\n");
#endif

        /* Find nearest neighbors */
        uint Lmin = uInsane;
        uint Rmin = uInsane;
        dist_t dtMinDist = BIG_DIST;
        for (j = 0; j < g_uLeafCount; ++j) {
            dist_t d;
            if (uInsane == g_uNodeIndex[j])
                continue;

            d = g_MinDist[j];
            if (d < dtMinDist) {
                dtMinDist = d;
                Lmin = j;
                Rmin = g_uNearestNeighbor[j];
                assert(uInsane != Rmin);
                assert(uInsane != g_uNodeIndex[Rmin]);
            }
        }

        assert(Lmin != uInsane);
        assert(Rmin != uInsane);
        assert(dtMinDist != BIG_DIST);

#if TRACE
        Info("Nearest neighbors Lmin %u[=%u] Rmin %u[=%u] dist %.3g\n",
             Lmin,
             g_uNodeIndex[Lmin],
             Rmin,
             g_uNodeIndex[Rmin],
             dtMinDist);
#endif

        /* Compute distances to new node
         * New node overwrites row currently assigned to Lmin
         */
        for ( j = 0; j < g_uLeafCount; ++j) {
            ulong vL, vR;
            dist_t dL, dR;
            dist_t dtNewDist;
            
            if (j == Lmin || j == Rmin)
                continue;
            if (uInsane == g_uNodeIndex[j])
                continue;

            vL = TriangleSubscript(Lmin, j);
            vR = TriangleSubscript(Rmin, j);
            dL = g_Dist[vL];
            dR = g_Dist[vR];
            dtNewDist = 0.0;

            switch (linkage) {
            case LINKAGE_AVG:
                dtNewDist = AVG(dL, dR);
                break;

            case LINKAGE_MIN:
                dtNewDist = MIN(dL, dR);
                break;

            case LINKAGE_MAX:
                dtNewDist = MAX(dL, dR);
                break;
/* couldn't be arsed to figure out proper usage of g_dSUEFF */
#if 0
            case LINKAGE_BIASED:
                dtNewDist = g_dSUEFF*AVG(dL, dR) + (1 - g_dSUEFF)*MIN(dL, dR);
                break;
#endif
            default:
                Log(&rLog, LOG_FATAL, "UPGMA2: Invalid LINKAGE_%u", linkage);
            }

            /* Nasty special case.
             * If nearest neighbor of j is Lmin or Rmin, then make the new
             * node (which overwrites the row currently occupied by Lmin)
             * the nearest neighbor. This situation can occur when there are
             * equal distances in the matrix. If we don't make this fix,
             * the nearest neighbor pointer for j would become invalid.
             * (We don't need to test for == Lmin, because in that case
             * the net change needed is zero due to the change in row
             * numbering).
             */
            if (g_uNearestNeighbor[j] == Rmin)
                g_uNearestNeighbor[j] = Lmin;

#if TRACE
            Info("New dist to %u = (%u/%.3g + %u/%.3g)/2 = %.3g\n",
                 j, Lmin, dL, Rmin, dR, dtNewDist);
#endif
            g_Dist[vL] = dtNewDist;
            if (dtNewDist < dtNewMinDist) {
                dtNewMinDist = dtNewDist;
                uNewNearestNeighbor = j;
            }
        }

        assert(g_uInternalNodeIndex < g_uLeafCount - 1 || BIG_DIST != dtNewMinDist);
        assert(g_uInternalNodeIndex < g_uLeafCount - 1 || uInsane != uNewNearestNeighbor);

        const ulong v = TriangleSubscript(Lmin, Rmin);
        const dist_t dLR = g_Dist[v];
        const dist_t dHeightNew = dLR/2;
        const uint uLeft = g_uNodeIndex[Lmin];
        const uint uRight = g_uNodeIndex[Rmin];
        const dist_t HeightLeft =
            uLeft < g_uLeafCount ? 0 : g_Height[uLeft - g_uLeafCount];
        const dist_t HeightRight =
            uRight < g_uLeafCount ? 0 : g_Height[uRight - g_uLeafCount];

        g_uLeft[g_uInternalNodeIndex] = uLeft;
        g_uRight[g_uInternalNodeIndex] = uRight;
        g_LeftLength[g_uInternalNodeIndex] = dHeightNew - HeightLeft;
        g_RightLength[g_uInternalNodeIndex] = dHeightNew - HeightRight;
        g_Height[g_uInternalNodeIndex] = dHeightNew;

        /* Row for left child overwritten by row for new node */
        g_uNodeIndex[Lmin] = g_uLeafCount + g_uInternalNodeIndex;
        g_uNearestNeighbor[Lmin] = uNewNearestNeighbor;
        g_MinDist[Lmin] = dtNewMinDist;

        /* Delete row for right child */
        g_uNodeIndex[Rmin] = uInsane;

#if TRACE
        Info("\nInternalNodeIndex=%u Lmin=%u Rmin=%u\n",
             g_uInternalNodeIndex, Lmin, Rmin);
        ListState();
#endif
    }

    uint uRoot = g_uLeafCount - 2;

#if TRACE
    Log(&rLog, LOG_FORCED_DEBUG, "uRoot=%d g_uLeafCount=%d g_uInternalNodeCount=%d", uRoot, g_uLeafCount, g_uInternalNodeCount);
    for (i=0; i<g_uInternalNodeCount; i++) {
        Log(&rLog, LOG_FORCED_DEBUG, "internal node=%d:  g_uLeft=%d g_uRight=%d g_LeftLength=%f g_RightLength=%f g_Height=%f",
                  i, g_uLeft[i], g_uRight[i],
                  g_LeftLength[i], g_RightLength[i],
                  g_Height[i]);
    }
    for (i=0; i<g_uLeafCount; i++) {
        Log(&rLog, LOG_FORCED_DEBUG, "leaf node=%d:  Ids=%d names=%s",
                  i, Ids[i], names[i]);
    }
#endif
    
    MuscleTreeCreate(tree, g_uLeafCount, uRoot,
                      g_uLeft, g_uRight,
                      g_LeftLength, g_RightLength,
                      Ids, names);
#if TRACE
    tree.LogMe();
#endif

    free(g_Dist);

    free(g_uNodeIndex);
    free(g_uNearestNeighbor);
    free(g_MinDist);
    free(g_Height);

    free(g_uLeft);
    free(g_uRight);
    free(g_LeftLength);
    free(g_RightLength);

    /* NOTE: Muscle's "Names" variable is here the argument "names" */
    free(Ids);
}
/***   end of UPGMA2   ***/
