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
 *  RCS $Id: mbed.c 316 2016-12-16 16:14:39Z fabian $
 *
 *
 * Reimplementation from scratch of mBed (Blackshields et al., 2010;
 * PMID 20470396) and addition of bisecting K-Means which allows
 * control over subcluster size
 *
 */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <assert.h>
#include <float.h>

#include "mbed.h"
#include "pair_dist.h"
#include "symmatrix.h"
#include "ktuple_pair.h"
#include "tree.h"
#include "util.h"
#include "progress.h"
#include "queue.h"

#include "log.h"
#include "kmpp/KMeans.h"
#include "mbed.h"

#define TIMING 0
#if TIMING
#include "squid/stopwatch.h"
#endif


/* If FULL_WITHIN_CLUSTER_DISTANCES is not 0, distances within each
 * bisecting kmeans subcluster are not estimated using the vectors,
 * but calculated normally (using ktuple or kimura). Surprisingly this
 * results in 3% loss on a Homfam p24-h2010-08-09 subset (100-5000
 * seqs in test, at least 5 ref seqs; MAX_SEQ 100 vs 10000; NUM_SEEDS
 * log2 instead of log2^2). And of course it slows things down.
 */
#define FULL_WITHIN_CLUSTER_DISTANCES 1

#define COMPUTE_WITHIN_SUBCLUSTER_AVERAGE 0

/* Cluster size limits. Maximum is a soft limit, which might be
 * exceeded if a K-Means split was unsuccesful, where unsuccesful
 * might also mean that the minimum required number seqs. was not
 * reached */
#if FULL_WITHIN_CLUSTER_DISTANCES
static const int MAX_ALLOWED_SEQ_PER_PRECLUSTER = 100;
static const int MIN_REQUIRED_SEQ_PER_PRECLUSTER = 1;
#else
static const int MAX_ALLOWED_SEQ_PER_PRECLUSTER = 10000;
static const int MIN_REQUIRED_SEQ_PER_PRECLUSTER = 100;
#endif

/* How many restarts per bisecting kmeans split. 10 give 0.5% increase
 * in quality on original HOMFAM over just 1. It also increases kmeans
 * time by a factor of ~3, but overall time is insignificant
 * compared to pairdist/progressive alignment part.
 */
static const int RESTARTS_PER_SPLIT = 10;


/* Use standard kmeans (lloyds) or kmeans++. Both are almost
 * indistinguishable here, although kmeans++ is slightly ahead the
 * fewer seeds you pick (and it should be better in theory)
 */
#define USE_KMEANS_LLOYDS 0


//#ifndef HAVE_LOG2
#ifndef CLUSTAL_OMEGA_HAVE_LOG2
//#define log2(x)  (log(x) / 0.69314718055994530942)
#define log2(x)  ( x<=0 ? (float)(-100000) : 1.442695041*log(x) )
#endif
#define NUMBER_OF_SEEDS(n) pow(log2(((double)n)), 2)


/* Seed selection method: SAMPLE_SEEDS_BY_LENGTH is the original mBed
 * approach: Sample iSeeds sequence with constant stride from length-sorted X.
 * It might be better to pick the seeds randomly, because length sorting might
 * be more prone to including outliers (e.g. very long and very short seqs).
 * However, we could never observer such a thing. So just stick to the
 * original version as this also removes the random element.
 */
enum SEED_SELECTION_TYPE {
    SELECT_SEEDS_RANDOMLY,
    SELECT_SEEDS_BY_LENGTH
};
#define SEED_SELECTION SELECT_SEEDS_BY_LENGTH


/* Tests on BAliBase (RV11,12,20,30,40,50; 10 runs each) show there is
 * no difference between mbed-trees created from cosine or euclidean
 * distances (simple version, just using disparities).
 */
#define USE_EUCLIDEAN_DISTANCE 1


/* print some mbed pre-cluster usage to screen */
#define PRINT_CLUSTER_DISTRIBUTION 0


#define TRACE 0


typedef struct {
    /* Number of final clusters
     */
    int iNClusters;

    /* Coordinates (columns) for each cluster (rows)
     * valid indices: [0...iNClusters][0...dim-1]
     */
    double **ppdClusterCenters;

    /* Dimensionality of input data and cluster center coordinates
     */
    int iDim;

    /* Number of objects per cluster
     */
    int *piNObjsPerCluster;

    /* Objects (indices) for each cluster. [i][j] will point to (index
     * of) object j in cluster i
     */
    int **ppiObjIndicesPerCluster;
} bisecting_kmeans_result_t;



/**
 * @brief Free KMeans result structure.
 *
 * @param[out]  prResult_p
 * K-Means result to free
 *
 * @see NewKMeansResult()
 */
void
FreeKMeansResult(bisecting_kmeans_result_t **prResult_p)
{
    int iAux;
    CKFREE((*prResult_p)->piNObjsPerCluster);
    for (iAux=0; iAux<(*prResult_p)->iNClusters; iAux++) {
        CKFREE((*prResult_p)->ppiObjIndicesPerCluster[iAux]);
        CKFREE((*prResult_p)->ppdClusterCenters[iAux]);
    }
    CKFREE((*prResult_p)->ppiObjIndicesPerCluster);
    CKFREE((*prResult_p)->ppdClusterCenters);
    (*prResult_p)->iNClusters = 0;
    (*prResult_p)->iDim = 0;
    CKFREE(*prResult_p);
}
/***   end: FreeKMeansResult()   ***/



/**
 * @brief Allocate new KMeans result
 *
 * @param[out] prKMeansResult_p
 * K-Means result to free
 *
 * @see NewKMeansResult()
 */
void
NewKMeansResult(bisecting_kmeans_result_t **prKMeansResult_p)
{
    (*prKMeansResult_p) = (bisecting_kmeans_result_t *)
        CKCALLOC(1, sizeof(bisecting_kmeans_result_t));
    (*prKMeansResult_p)->iNClusters = 0;
    (*prKMeansResult_p)->iDim = 0;
    (*prKMeansResult_p)->ppdClusterCenters = NULL;
    (*prKMeansResult_p)->piNObjsPerCluster = NULL;
    (*prKMeansResult_p)->ppiObjIndicesPerCluster = NULL;

}
/***   end: NewKMeansResult()   ***/



/**
 * @brief Calculate the euclidean distance between two vectors
 *
 * @param[in] v1
 * First vector with dim dimensions
 * @param[in] v2
 * Second vector with dim dimensions
 * @param[in] dim
 * Dimensionality of v1 and v2
 *
 * @return euclidean distance as double
 *
 * @note Could probably be inlined. Also perfect for SSE
 */
double
EuclDist(const double *v1, const double *v2, const int dim)
{
    int i; /* aux */
    double dist; /* returned distance */

    assert(v1!=NULL);
    assert(v2!=NULL);

    dist = 0.0;
    for (i=0; i<dim; i++) {
        dist += pow(v1[i]-v2[i], 2.0);
    }

    return sqrt(dist);
}
/***   end: EuclDist()   ***/



/**
 * @brief Calculate the cosine distance between two vectors
 *
 * @param[in] v1
 * First vector with dim dimensions
 * @param[in] v2
 * Second vector with dim dimensions
 * @param[in] dim
 * Dimensionality of v1 and v2
 *
 * @return cosine distance as double
 *
 * @note Could probably be inlined. Also perfect for SSE
 */
double
CosDist(const double *v1, const double *v2, const int dim)
{
    int i; /* aux */
    double dist; /* returned distance */
    double s, sq1, sq2;

    assert(v1!=NULL);
    assert(v2!=NULL);

    s = 0.0;
    sq1 = 0.0;
    sq2 = 0.0;
    for (i=0; i<dim; i++) {
        s += (v1[i]*v2[i]);
        sq1 += pow(v1[i], 2.0);
        sq2 += pow(v2[i], 2.0);
    }
    sq1 = sqrt(sq1);
    sq2 = sqrt(sq2);

    if ((sq1 * sq2) < DBL_EPSILON) {
        dist = 1 - s / (sq1 * sq2);
    } else {
        /* if the denominator gets small, the fraction gets big, hence dist
           gets small: */
        dist = 0.0;
    }
    return dist;
}
/***   end: CosDist()   ***/



/**
 * @brief convert sequences into mbed-like (distance) vector
 * representation. Seeds (prMSeq sequence indices) have to be picked before
 *
 * @param[out] ppdSeqVec
 * Pointer to preallocated matrix of size nseqs x iSeeds
 * @param[in] prMSeq
 * Sequences which are to be converted
 * @param[in] piSeeds
 * Array of sequences in prMSeq which are to be used as seeds.
 * @param[in] iNumSeeds
 * Number of seeds/elements in piSeeds
 * @param[in] iPairDistType
 * Type of pairwise distance comparison
 *
 */
int
SeqToVec(double **ppdSeqVec, mseq_t *prMSeq, 
         int *piSeeds, const int iNumSeeds,
         const int iPairDistType)
{
    symmatrix_t *prDistmat = NULL;
    int iSeqIndex;
    int iSeedIndex;
   /* indices for restoring order */
    int *restore = NULL; 
    /* sorted copy of piSeeds */
    int *piSortedSeeds = NULL; 

#if TIMING
    Stopwatch_t *stopwatch = StopwatchCreate();
    StopwatchZero(stopwatch);
    StopwatchStart(stopwatch);
#endif
    
    assert(NULL != ppdSeqVec);
    assert(iNumSeeds>0 && iNumSeeds<=prMSeq->nseqs);

    piSortedSeeds = CKMALLOC(iNumSeeds * sizeof(int));
    memcpy(piSortedSeeds, piSeeds, iNumSeeds*sizeof(int));

    /* need them sorted, otherwise the swapping approach below might
     * break 
     */
    qsort(piSortedSeeds, iNumSeeds, sizeof(int), IntCmp);


    /* rearrange mseq order so that we can use ktuple_pairdist code as
     * is. That is, swap seeds and non-seeds such that the seeds end
     * up in front of mseq. This way we can use the KTuplePairDist()
     * code, without making a copy of mseq. Also, keep an array of
     * indices which serves to restore the original sequence order
     * after putting the seeds in front
     *
     */
    restore = (int *) CKMALLOC(prMSeq->nseqs * sizeof(int));
    for (iSeqIndex=0; iSeqIndex<prMSeq->nseqs; iSeqIndex++) {
        restore[iSeqIndex] = iSeqIndex;
    }
    for (iSeedIndex=0; iSeedIndex<iNumSeeds; iSeedIndex++) {
#if TRACE
        Log(&rLog, LOG_FORCED_DEBUG, "Swapping seqs %d and %u", piSortedSeeds[iSeedIndex], iSeedIndex);
#endif
        if (piSortedSeeds[iSeedIndex]!=iSeedIndex) {
            int swap;
            SeqSwap(prMSeq, piSortedSeeds[iSeedIndex], iSeedIndex);

            swap = restore[iSeedIndex];
            restore[iSeedIndex] = restore[piSortedSeeds[iSeedIndex]];
            restore[piSortedSeeds[iSeedIndex]] = swap;
        }
    }
#if TRACE
    printf("DEBUG(%s|%s():%d): restore indices =",
           __FILE__, __FUNCTION__, __LINE__);
    for (iSeqIndex=0; iSeqIndex<prMSeq->nseqs; iSeqIndex++) {
        printf(" %u:%u", iSeqIndex, restore[iSeqIndex]);
    }
    printf("\n");

    for (iSeqIndex=0; iSeqIndex<prMSeq->nseqs; iSeqIndex++) {
        Log(&rLog, LOG_FORCED_DEBUG, "swapped seq no %d now:  seq = %s",
                  iSeqIndex, prMSeq->sqinfo[iSeqIndex].name);
    }
#endif


    /* convert sequences into vectors of distances by calling pairwise
     * distance function for each seq/seed pair
     *
     */
    /* 4th argument (FALSE) is bPercID, in mBed mode never use percent-identities */
    if (0 != PairDistances(&prDistmat, prMSeq, iPairDistType, FALSE, 
                           0, iNumSeeds, 0, prMSeq->nseqs,
                           NULL, NULL)) {
        Log(&rLog, LOG_ERROR, "Could not compute pairwise distances for mbed.");
        FreeSymMatrix(&prDistmat);
        CKFREE(piSortedSeeds);
        CKFREE(restore);
        return -1;
    }
#if TRACE
    {
        char **labels;
        labels = (char **) CKMALLOC(prMSeq->nseqs * sizeof(char *));
        for (iSeqIndex=0; iSeqIndex<prMSeq->nseqs; iSeqIndex++) {
            labels[iSeqIndex] = prMSeq->sqinfo[iSeqIndex].name;
        }
        SymMatrixPrint(prDistmat, labels, NULL, FALSE);
        CKFREE(labels);
    }
#endif


    /* update ppdSeqVec according to prDistmat. keep in mind that we
     * changed mseq order
     *
     */
    for (iSeqIndex=0; iSeqIndex<prMSeq->nseqs; iSeqIndex++) {
        for (iSeedIndex=0; iSeedIndex<iNumSeeds; iSeedIndex++) {
            ppdSeqVec[restore[iSeqIndex]][iSeedIndex] =
                SymMatrixGetValue(prDistmat, iSeqIndex, iSeedIndex);
        }
    }


    /* restore mseq order by reverse swapping
     *
     * need reverse order now, so start from top index
     */
    iSeedIndex=iNumSeeds-1;
    do {
#if TRACE
        Log(&rLog, LOG_FORCED_DEBUG, "Swapping seqs %d and %d", piSortedSeeds[iSeedIndex], iSeedIndex);
#endif
        if (piSortedSeeds[iSeedIndex]!=iSeedIndex) {
            int swap;

            SeqSwap(prMSeq, piSortedSeeds[iSeedIndex], iSeedIndex);

            swap = restore[iSeedIndex];
            restore[iSeedIndex] = restore[piSortedSeeds[iSeedIndex]];
            restore[piSortedSeeds[iSeedIndex]] = swap;
        }
    } while (0 != iSeedIndex--);
#if TRACE
     for (iSeqIndex=0; iSeqIndex<prMSeq->nseqs; iSeqIndex++) {
         Log(&rLog, LOG_FORCED_DEBUG, "restored seq no %d:  seq = %s %s",
                   iSeqIndex, prMSeq->sqinfo[iSeqIndex].name, prMSeq->seq[iSeqIndex]);
     }
#endif

#if TRACE
    for (iSeqIndex=0; iSeqIndex<prMSeq->nseqs; iSeqIndex++) {
        printf("DEBUG: seq %-4u as distance vector =", iSeqIndex);
        for (iSeedIndex=0; iSeedIndex<iNumSeeds; iSeedIndex++) {
            printf(" %f", ppdSeqVec[iSeqIndex][iSeedIndex]);
        }
        printf("\n");
    }
#endif

    
    FreeSymMatrix(&prDistmat);
    CKFREE(restore);
    CKFREE(piSortedSeeds);
#if TIMING
    StopwatchStop(stopwatch);
    StopwatchDisplay(stdout, "Total time for SeqToVec(): ", stopwatch);
    StopwatchFree(stopwatch);
#endif


    return 0;
}
/***   end: SeqToVec()   ***/



/**
 * @brief Select seeds to be used from an prMSeq
 *
 * @param[out] piSeeds
 * Will store the indices of prMSeqs seqs used to be as seeds here. Must be preallocated.
 * @param[in] iNumSeeds
 * Number of seeds to be picked
 * @param[in] iSelectionMethod
 * Seed selection method to be used
 * @param[in] prMSeq
 * The prMSeq structure to pick sequences from
 *
 * @return: Non-zero on failure
 *
 */    
int
SeedSelection(int *piSeeds, int iNumSeeds, int iSelectionMethod,  mseq_t *prMSeq)
{
    /* seed iterator */
    int iSeedIndex;
    /* sequence iterator */
    int iSeqIndex;

    if (SELECT_SEEDS_RANDOMLY == iSelectionMethod) {
        int *piPermArray;

        Log(&rLog, LOG_INFO,
             "Using %d seeds (randomly chosen) for mBed (from a total of %d sequences)",
             iNumSeeds, prMSeq->nseqs);

        PermutationArray(&piPermArray, iNumSeeds);
        /* copy to piSeeds */
        for (iSeedIndex=0; iSeedIndex<iNumSeeds; iSeedIndex++) {
            piSeeds[iSeedIndex] = piPermArray[iSeedIndex];
        }
        CKFREE(piPermArray);
        
    }  else if (SELECT_SEEDS_BY_LENGTH == iSelectionMethod) {

        /* step size for picking with constant stride */
        int iStep;
        int *piSeqLen =  CKMALLOC(prMSeq->nseqs * sizeof(int));
        int *piOrder = CKMALLOC(prMSeq->nseqs * sizeof(int));

        Log(&rLog, LOG_INFO,
             "Using %d seeds (chosen with constant stride from length sorted seqs) for mBed (from a total of %d sequences)",
             iNumSeeds, prMSeq->nseqs);

        iStep = prMSeq->nseqs / iNumSeeds; /* iStep will never get too big
                                            * due to rounding */
        /* first get an array of seq indices order according to
           corresponding sequence length: piOrder */
        for (iSeqIndex=0; iSeqIndex<prMSeq->nseqs; iSeqIndex++) {
            piSeqLen[iSeqIndex] = prMSeq->sqinfo[iSeqIndex].len;
        }
        QSortAndTrackIndex(piOrder, piSeqLen, prMSeq->nseqs, 'd', FALSE);
#if 0
        for (iSeqIndex=0; iSeqIndex<prMSeq->nseqs; iSeqIndex++) {
            Log(&rLog, LOG_FORCED_DEBUG, "Descending order (no %d): seq %d has len %d)",
                      iSeqIndex, piOrder[iSeqIndex], piSeqLen[piOrder[iSeqIndex]]);
        }
#endif
        CKFREE(piSeqLen);
        iSeedIndex = 0;
        for (iSeqIndex=0; iSeedIndex<iNumSeeds; iSeqIndex+=iStep) {
            piSeeds[iSeedIndex++] = piOrder[iSeqIndex];
        }
        CKFREE(piOrder);

    } else {
        
        Log(&rLog, LOG_ERROR, "Internal error: unknown seed selection type");
        return -1;
    }


#ifdef PICK_SEEDS_FROM_KMEANS_CLUSTERS
    /* do initial mbedding and kmeans. then pick seeds from each cluster and
       do full mbed with these seeds. idea is that those seeds represent the
       sequence space better than random seeds. however, this reduces the
       quality slightly. random is almost always better.
    */

    if (1) {
        /* seqs represented as distance vectors */
        double **ppdSeqVec;
        double dCost = -1.0;
        /* assignments for each seq to the K newly created clusters */
        int *piKMeansClusterAssignments = CKMALLOC(prMSeq->nseqs * sizeof(int));
        double *pdKMeansClusterCenters = CKCALLOC(iNumSeeds * iNumSeeds, sizeof(double));
        /* create a copy of ppdSeqVec suitable for KMeans */
        double *pdKMeansVectors = CKMALLOC(prMSeq->nseqs * iNumSeeds * sizeof(double));
        bool *pbClusterUsed = CKCALLOC(iNumSeeds, sizeof(bool));

        Log(&rLog, LOG_FORCED_DEBUG, "%s", "FIXME Experimental feature: K-means on first round embedding to get better seeds");
        Log(&rLog, LOG_FORCED_DEBUG, "%s", "FIXME Reuse seeds from clusters");
        Log(&rLog, LOG_FORCED_DEBUG, "%s", "FIXME Could try log(n) instead of log(n)^2 in first round");
        Log(&rLog, LOG_FORCED_DEBUG, "%s", "FIXME hardcoded iPairDistType PAIRDIST_KTUPLE");
        Log(&rLog, LOG_FORCED_DEBUG, "%s", "FIXME Show that this is fast *and* good");

        /*
          Log(&rLog, LOG_FORCED_DEBUG, "%s", "Overriding iNumSeeds");
          iNumSeeds = (int) pow(log2((double)prMSeq->nseqs), 2);
        */

        ppdSeqVec = (double **) CKMALLOC(prMSeq->nseqs * sizeof(double *));
        for (iSeqIndex=0; iSeqIndex<prMSeq->nseqs; iSeqIndex++) {
            ppdSeqVec[iSeqIndex] = (double *) CKMALLOC(iNumSeeds * sizeof(double));
        }
        if (0 != SeqToVec(ppdSeqVec, prMSeq, piSeeds, iNumSeeds, PAIRDIST_KTUPLE)) {
            Log(&rLog, LOG_ERROR, "Could not convert sequences into vectors for mbed");
            return -1;
        }

        for (iSeqIndex=0; iSeqIndex<prMSeq->nseqs; iSeqIndex++) {
            int iDim;
            for (iDim=0; iDim<iNumSeeds; ++iDim) {
                pdKMeansVectors[iDim*iSeqIndex + iDim] = ppdSeqVec[iSeqIndex][iDim];
            }
        }


        Log(&rLog, LOG_FORCED_DEBUG, "%s\n", "FIXME hardcoded RESTARTS_PER_SPLIT");
        dCost = KMeans(prMSeq->nseqs, iNumSeeds, iNumSeeds,
                       pdKMeansVectors,
                       RESTARTS_PER_SPLIT, USE_KMEANS_LLOYDS,
                       pdKMeansClusterCenters, piKMeansClusterAssignments);
        Log(&rLog, LOG_FORCED_DEBUG, "Best split cost = %f", dCost);
        
        Log(&rLog, LOG_FORCED_DEBUG, "%s", "FIXME Check for Nan in cluster centers");

#if TRACE
        Log(&rLog, LOG_FORCED_DEBUG, "%s", "K-Means output:");
        for (iSeqIndex=0; iSeqIndex<prMSeq->nseqs; iSeqIndex++) {
            Log(&rLog, LOG_FORCED_DEBUG, " Raw assignment: Seq %u (%s) to cluster %d",
                      iSeqIndex,  prMSeq->sqinfo[iSeqIndex].name,
                      piKMeansClusterAssignments[iSeqIndex]);
        }
#endif

        Log(&rLog, LOG_FORCED_DEBUG, "FIXME %s", "proof of concept implementation: Pick first sequences from each clusters instead of reusing seeds");
        iSeedIndex = 0;
        for (iSeqIndex=0; iSeqIndex<prMSeq->nseqs; iSeqIndex++) {
            int iAssignedCluster = piKMeansClusterAssignments[iSeqIndex];
            if (pbClusterUsed[iAssignedCluster]) {
                continue;
            } else {
                /*LOG_DEBUG("Picked seed %d from cluster %d", iSeqIndex,iAssignedCluster);*/
                piSeeds[iSeedIndex++] = iSeqIndex;
                pbClusterUsed[iAssignedCluster] = TRUE;
            }
        }
        CKFREE(pbClusterUsed);

        CKFREE(pdKMeansVectors);
        CKFREE(pdKMeansClusterCenters);
        CKFREE(piKMeansClusterAssignments);
        
    }
    /* if (1) */
#endif


    if (LOG_DEBUG <= rLog.iLogLevelEnabled) {
        for (iSeedIndex=0; iSeedIndex<iNumSeeds; iSeedIndex++) {
            Log(&rLog, LOG_DEBUG, "Picked sequence %d (%s) as seed no %d",
                 piSeeds[iSeedIndex], prMSeq->sqinfo[piSeeds[iSeedIndex]].name, iSeedIndex);
        }
    }
    
    return 0; 
}
/* end of SeedSelection() */




/**
 * @brief Bisecting K-Means clustering. Repeatedly calls K-Means with
 * a K of 2 until no cluster has more than iMaxAllowedObjsPerCluster.
 *
 * @param[out] prKMeansResult_p
 * Result of Bisecting KMeans. Will be allocated here.
 * Caller has to free. See @see FreeKMeansResult()
 * @param[in] iNObjs
 * Number of objects/sequences to cluster
 * @param[in] iDim
 * Dimensionality of input data
 * @param[in] ppdVectors
 * each row holds iDim points for this object's coordinates
 * @param[in] iMinRequiredObjsPerCluster
 * Minimum number of objects per Cluster (inclusive)/
 * @param[in] iMaxAllowedObjsPerCluster
 * Maximum number of objects per Cluster (inclusive). Function returns once no
 * cluster contains more then this number of objects. Soft limit!
 *
 * @note Convoluted code. Could use some restructuring. My apologies.
 * AW
 *
 */
void
BisectingKmeans(bisecting_kmeans_result_t **prKMeansResult_p,
                const int iNObjs, const int iDim, double **ppdVectors,
                const int iMinRequiredObjsPerCluster,
                const int iMaxAllowedObjsPerCluster, char ***ppcClusterSplits_p)
{
    int iN, iD;
    /* cluster centers for each cluster created at each split */
    double *pdKClusterCenters;
    /* keep track of updated object indices per newly created
     * cluster */
    int *piCurObjToUpdate;
    /* number of times we called 2-means */
    int iNRounds;
    /* flag for unsuccessful cluster split */
    bool bNaNDetected = FALSE;
    /* flag for detected small cluster after split  */
    bool bSmallClusterDetected;
    /* queue of clusters which are to be split */
    int_queue_t rClusterSplitQueue;
    /* */
    int iLog2N2 = 0;
    
#if TIMING
    Stopwatch_t *stopwatch = StopwatchCreate();
#endif


    piCurObjToUpdate = (int *) CKMALLOC(2 * sizeof(int));

    NewKMeansResult(prKMeansResult_p);

    /* new cluster centers created at each split/KMeans run
     */
    pdKClusterCenters = (double *) CKCALLOC(2 * iDim, sizeof(double));

    /* init results by setting a first cluster that contains all objects
     */
    (*prKMeansResult_p)->iNClusters = 1;
    (*prKMeansResult_p)->iDim = iDim;
    /* fake center coordinates of first cluster */
    (*prKMeansResult_p)->ppdClusterCenters =
        (double **) CKMALLOC(1 * sizeof(double *));
    (*prKMeansResult_p)->ppdClusterCenters[0] =
        (double *) CKMALLOC(iDim * sizeof(double));
    /* objects per cluster */
    (*prKMeansResult_p)->piNObjsPerCluster =
        (int *) CKMALLOC(1 * sizeof(int));
    (*prKMeansResult_p)->piNObjsPerCluster[0] = iNObjs;
    /* object indices per cluster */
    (*prKMeansResult_p)->ppiObjIndicesPerCluster =
        (int **) CKMALLOC(1 * sizeof(int *));
    (*prKMeansResult_p)->ppiObjIndicesPerCluster[0] = (int *)
        CKMALLOC(iNObjs * sizeof(int));
    for (iN=0; iN<iNObjs; iN++) {
        (*prKMeansResult_p)->ppiObjIndicesPerCluster[0][iN] = iN;
    }


    /*
     * this is an array that encodes where a sequences fell at a bisecting k-means
     * ppcClusterSplits_p can consume quite a bit of memory, 
     but it is only needed if clustering info is written to file.
     Use ppcClusterSplits_p as a flag, initialise in the calling function (Mbed) 
     to NULL, if no info written out then leave ppcClusterSplits_p=NULL.
     If info is to be written out then malloc one char*, so that 
     ppcClusterSplits_p!=NULL. In that case realloc proper amount.
     Only expect a certain number of splits, iLog2N2. 
     Worst case would be iNObjs but then ppcClusterSplits_p 
     would be quadratic, which is too big for many sequences.
     */
    if (NULL != *ppcClusterSplits_p){
        /* (square) of logarithm (base 2) of number of objects */
        double dLog2N = log10(iNObjs) / log10(2);
        iLog2N2 = (int)(dLog2N*dLog2N);

        *ppcClusterSplits_p      = (char **)CKREALLOC(*ppcClusterSplits_p, iNObjs*sizeof(char *));
        (*ppcClusterSplits_p)[0] =  (char *)CKMALLOC(iNObjs*iLog2N2);
        (*ppcClusterSplits_p)[0][0] = '\0';
        for (iN = 1; iN < iNObjs; iN++){
            (*ppcClusterSplits_p)[iN] = (*ppcClusterSplits_p)[iN-1] + iLog2N2;
            (*ppcClusterSplits_p)[iN][0] = '\0';
        }
    }


    /* Starting with the first cluster that now contains all the
     * sequences keep splitting until no cluster contains more than
     * iMaxAllowedObjsPerCluster
     *
     * Keep a queue of clusters (rClusterSplitQueue) to split
     *
     * At each split values/memory of the just split cluster will be
     * reused and exactly one new only allocated.
     *
     * Example:
     * Given the following cluster assignments
     * 0 0 0 0 1 1 2 2 2 3 3 3 3 3 3 3 4 4
     * and a K-Means split in cluster 3 at |:
     * 0 0 0 0 1 1 2 2 2 3 3 3 | 3 3 3 3 4 4
     * The result should be this:
     * 0 0 0 0 1 1 2 2 2 3 3 3 | 5 5 5 5 4 4
     *
     *
     */
    INT_QUEUE_INIT(&rClusterSplitQueue);
    if (iNObjs>iMaxAllowedObjsPerCluster) {
        /* pish fake first cluster index */
        INT_QUEUE_PUSH(&rClusterSplitQueue, 0); 
    }
    iNRounds = 0;
    while (! INT_QUEUE_EMPTY(&rClusterSplitQueue)) {
        /* assignments for each seq to the K newly created clusters */
        int *piKClusterAssignments;
        /* number of objects in cluster that is to be split */
        int iNObjsInClusterToSplit;
        /* coordinates of points in cluster that is to be split
         * array of size n*d where [d*i + j] gives coordinate j of point i
         */
        double *pdVectorsInClusterToSplit;
        /* copy of object indices in split cluster */
        int *piObjIndicesOfSplitCluster;
        /* best cost of kmeans rounds */
        double dCost = -1.0;

        /*  indices for the two created clusters
         */
        /* index of cluster to split */
        int iClusterToSplot;
        /* index of newly created cluster at each round. data for the
           other created cluster goes to the one just split,
           i.e. (*piClusterToSplot) */
        int iNewClusterIdx;

#if TIMING
        StopwatchZero(stopwatch);
        StopwatchStart(stopwatch);
#endif

        INT_QUEUE_POP(&rClusterSplitQueue, &iClusterToSplot);
        
        iNObjsInClusterToSplit = (*prKMeansResult_p)->piNObjsPerCluster[iClusterToSplot];
        piKClusterAssignments = (int *)
            CKMALLOC(iNObjsInClusterToSplit * sizeof(int));

        pdVectorsInClusterToSplit = (double *)
            CKMALLOC(iNObjsInClusterToSplit * iDim * sizeof(double));
        for (iN=0; iN<iNObjsInClusterToSplit; iN++) {
            for (iD=0; iD<iDim; ++iD) {
                int iThisObjIdx =
                    (*prKMeansResult_p)->ppiObjIndicesPerCluster[iClusterToSplot][iN];
                pdVectorsInClusterToSplit[iDim*iN + iD] = ppdVectors[iThisObjIdx][iD];
            }
        }

#if TRACE
        Log(&rLog, LOG_FORCED_DEBUG, "Round %d: Will split cluster %d which has %d objects",
            iNRounds, iClusterToSplot, iNObjsInClusterToSplit);
        fprintf(stderr, "DEBUG(%s|%s():%d): Object indices in cluster to split are:",
                __FILE__, __FUNCTION__, __LINE__);
        for (iN=0; iN<iNObjsInClusterToSplit; iN++) {
            fprintf(stderr, " %u",
                    (*prKMeansResult_p)->ppiObjIndicesPerCluster[iClusterToSplot][iN]);
        }
        fprintf(stderr, "\n");
        (void) fflush(stderr);
#endif

#if TRACE
        for (iN=0; iN<iNObjsInClusterToSplit; iN++) {
            fprintf(stderr,
                    "DEBUG(%s|%s():%d): Coordinate of object %u (real index %u) in cluster to split:",
                    __FILE__, __FUNCTION__, __LINE__, iN,
                    (*prKMeansResult_p)->ppiObjIndicesPerCluster[iClusterToSplot][iN]);
            for (iD=0; iD<iDim; iD++) {
                fprintf(stderr, " %f", pdVectorsInClusterToSplit[iDim*iN + iD]);
            }
            fprintf(stderr, "\n");
            (void) fflush(stderr);
        }
#endif
        
       /* KMeans(1 "The number of points in the data set",
        *        2 "The number of clusters to look for",
        *        3 "The number of dimensions that the data set lives in",
        *        4 "points: An array of size n*d where points[d*i + j] gives coordinate j of point i",
        *        5 "attempts: The number of times to independently run k-means",
        *        6 "use_lloyds_method: uses kmpp if false, otherwise lloyds method",
        *        7 "centers: This can either be null or an array of size k*d.
        *           In the latter case, centers[d*i + j] will give coordinate j of center i.
        *           If the cluster is unused, it will contain NaN instead.",
        *        8 "assignments: This can either be null or an array of size n.
        *           In the latter case, it will be filled with the cluster that each point is assigned to
        *           (an integer between 0 and k-1 inclusive).");
        */
        dCost = KMeans(iNObjsInClusterToSplit, 2, iDim,
                       pdVectorsInClusterToSplit,
                       RESTARTS_PER_SPLIT, USE_KMEANS_LLOYDS,
                       pdKClusterCenters, piKClusterAssignments);

#if TRACE
        Log(&rLog, LOG_FORCED_DEBUG, "%s", "Raw K-Means output:");
        for (iN=0; iN<2; iN++) {
            fprintf(stderr, "DEBUG(%s|%s():%d): Cluster Center %u =",
                    __FILE__, __FUNCTION__, __LINE__, iN);
            for (iD=0; iD<iDim; iD++) {
                fprintf(stderr, " %f", pdKClusterCenters[iN*iDim+iD]);
            }
            fprintf(stderr, "\n");
            (void) fflush(stderr);
        }
        for (iN=0; iN<iNObjsInClusterToSplit; iN++) {
            Log(&rLog, LOG_FORCED_DEBUG, " Raw assignment: Seq %u to cluster %d (of #%u)",
                      iN,  piKClusterAssignments[iN], 2);
        }
#endif
        
        /* real index of one of the newly created clusters. the other
         * one is iClusterToSplot */
        iNewClusterIdx = (*prKMeansResult_p)->iNClusters;

        
        /* We don't want Clusters which are too small. Check here if a
         * split created a small cluster and if yes, discard the
         * solution. Because the cluster has already been removed from
         * the queue, we can just continue.
         */
        bSmallClusterDetected = FALSE;
        if (iMinRequiredObjsPerCluster>1) {
            int iNObjsCluster[2];
            iNObjsCluster[0] = 0; /* first cluster */
            iNObjsCluster[1] = 0; /* second cluster */
            for (iN=0; iN<iNObjsInClusterToSplit; iN++) {
                iNObjsCluster[piKClusterAssignments[iN]]+=1;
            }
            
            if (iNObjsCluster[0]<iMinRequiredObjsPerCluster
                ||
                iNObjsCluster[1]<iMinRequiredObjsPerCluster) {
                bSmallClusterDetected = TRUE;
                Log(&rLog, LOG_FORCED_DEBUG, "Skipping this split because objs in 1st/2nd cluster = %d/%d < %d",
                          iNObjsCluster[0], iNObjsCluster[1], iMinRequiredObjsPerCluster);
            }
        }
        
        /* Same logic as for small clusters applies if KMeans couldn't
         * split the cluster. In this case one of its center
         * coordinates will be NaN, which we check for in the
         * following.
         */
        if (! bSmallClusterDetected) {
            for (iN=0; iN<2; iN++) {
                bNaNDetected = FALSE;
                for (iD=0; iD<iDim; iD++) {
                    if (0 != isnan(pdKClusterCenters[iN*iDim+iN])) {
                        /* Got NaN as coordinate after splitting */
                        Log(&rLog, LOG_WARN, "%s(): Can't split cluster no. %d which has %d objects any further. %s",
                             __FUNCTION__,
                             iClusterToSplot, iNObjsInClusterToSplit,
                            "Hope it's not too big and doesn't slow things down."
                            );
                        bNaNDetected = TRUE;
                        break;
                    }
                }
                if (bNaNDetected) {
                    break;
                }
            }
        }
        
        /* Discarding split of this cluster. It has been removed from the
         * queue already. so just continue
         */
        if (bNaNDetected || bSmallClusterDetected) {
            CKFREE(piKClusterAssignments);
            CKFREE(pdVectorsInClusterToSplit);
            continue;
        }


        /* update cluster centers: pdClusterCenters
         *
         */
        /* reuse memory of existing/old/split cluster
         */
        for (iN=0; iN<iDim; iN++) {
            double dCoord = pdKClusterCenters[0*iDim+iN];
            (*prKMeansResult_p)->ppdClusterCenters[iClusterToSplot][iN] = dCoord;
        }
        /* realloc and set new one
         */
        (*prKMeansResult_p)->ppdClusterCenters = (double **)
            CKREALLOC((*prKMeansResult_p)->ppdClusterCenters,
                      ((*prKMeansResult_p)->iNClusters+1)  * sizeof(double *));
        (*prKMeansResult_p)->ppdClusterCenters[iNewClusterIdx] = (double *)
            CKMALLOC(iDim * sizeof(double));
        for (iD=0; iD<iDim; iD++) {
            double dCoord = pdKClusterCenters[1*iDim+iD];
            (*prKMeansResult_p)->ppdClusterCenters[iNewClusterIdx][iD] = dCoord;
        }
#if TRACE
        Log(&rLog, LOG_FORCED_DEBUG, "%s", "* Update of cluster centers done. Cluster centers so far:");
        for (iN=0; iN<(*prKMeansResult_p)->iNClusters+1; iN++) {
            fprintf(stderr, "DEBUG(%s|%s():%d): Center %u =",
                    __FILE__, __FUNCTION__, __LINE__, iN);
            for (iD=0; iD<iDim; iD++) {
                fprintf(stderr, " %f", (*prKMeansResult_p)->ppdClusterCenters[iN][iD]);
            }
            fprintf(stderr, "\n");
            (void) fflush(stderr);
        }
#endif


        /* update #seqs per cluster: piNObjsPerCluster
         *
         */
        (*prKMeansResult_p)->piNObjsPerCluster = (int *)
            CKREALLOC((*prKMeansResult_p)->piNObjsPerCluster,
                      ((*prKMeansResult_p)->iNClusters+1) * sizeof(int));
        /* init new and old one to zero */
        (*prKMeansResult_p)->piNObjsPerCluster[iClusterToSplot] = 0;
        (*prKMeansResult_p)->piNObjsPerCluster[iNewClusterIdx] = 0;
        /* now update values */
        for (iN=0; iN<iNObjsInClusterToSplit; iN++) {
            if (0 == piKClusterAssignments[iN]) {
                (*prKMeansResult_p)->piNObjsPerCluster[iClusterToSplot] += 1;
            } else if (1 == piKClusterAssignments[iN]) {
                (*prKMeansResult_p)->piNObjsPerCluster[iNewClusterIdx] += 1;
            } else {
                /* there used to be code for iK>=2 in r101 */
                Log(&rLog, LOG_FATAL, "Internal error: split into more than two clusters (got assignment %d)",
                      piKClusterAssignments[iN]);
            }
        }
#if TRACE
        Log(&rLog, LOG_FORCED_DEBUG, "%s", "* Update of NObjs per cluster done:");
        for (iN=0; iN<(*prKMeansResult_p)->iNClusters+1; iN++) {
            Log(&rLog, LOG_FORCED_DEBUG, "Cluster %d contains %d sequences",
                      iN, (*prKMeansResult_p)->piNObjsPerCluster[iN]);
        }
#endif
        /* queue clusters if still they are still too big
         */
        if ((*prKMeansResult_p)->piNObjsPerCluster[iClusterToSplot] > iMaxAllowedObjsPerCluster) {
            INT_QUEUE_PUSH(&rClusterSplitQueue, iClusterToSplot);
        }
        if ((*prKMeansResult_p)->piNObjsPerCluster[iNewClusterIdx] > iMaxAllowedObjsPerCluster) {
            INT_QUEUE_PUSH(&rClusterSplitQueue, iNewClusterIdx);
        }
        
        

        /* update which objects belong to which cluster:
         *
         * note:  piNObjsPerCluster needs to be already updated
         *
         */
        /* create a copy of the object indices in the split cluster
         * (original will be overwritten)
         */
        piObjIndicesOfSplitCluster = (int *) CKMALLOC(iNObjsInClusterToSplit * sizeof(int));
        memcpy(piObjIndicesOfSplitCluster,
               (*prKMeansResult_p)->ppiObjIndicesPerCluster[iClusterToSplot],
               iNObjsInClusterToSplit * sizeof(int));

        (*prKMeansResult_p)->ppiObjIndicesPerCluster = (int **)
            CKREALLOC((*prKMeansResult_p)->ppiObjIndicesPerCluster,
                      ((*prKMeansResult_p)->iNClusters+1) * sizeof(int *));


        (*prKMeansResult_p)->ppiObjIndicesPerCluster[iClusterToSplot] =
            (int *) CKREALLOC((*prKMeansResult_p)->ppiObjIndicesPerCluster[iClusterToSplot],
                              (*prKMeansResult_p)->piNObjsPerCluster[iClusterToSplot] * sizeof(int));

        (*prKMeansResult_p)->ppiObjIndicesPerCluster[iNewClusterIdx] =
            (int *) CKMALLOC((*prKMeansResult_p)->piNObjsPerCluster[iNewClusterIdx] * sizeof(int));


        /* now reassign the object indices to the assigned cluster
         */
        piCurObjToUpdate[0] = 0;
        piCurObjToUpdate[1] = 0;
        for (iN=0; iN<iNObjsInClusterToSplit; iN++) {
            int iThisObjIdx = piObjIndicesOfSplitCluster[iN];
            int iThisClusterAssignment =  piKClusterAssignments[iN];
            int iThisOffset = piCurObjToUpdate[iThisClusterAssignment];
            int iThisClusterIndex = 0;

            if (0 == iThisClusterAssignment) {
                iThisClusterIndex = iClusterToSplot;
                if (*ppcClusterSplits_p && strlen((*ppcClusterSplits_p)[iThisObjIdx])<iLog2N2){
                    strcat((*ppcClusterSplits_p)[iThisObjIdx], "0");
                }
            } else if (1 == iThisClusterAssignment) {
                iThisClusterIndex = iNewClusterIdx;
                if (*ppcClusterSplits_p && strlen((*ppcClusterSplits_p)[iThisObjIdx])<iLog2N2){
                    strcat((*ppcClusterSplits_p)[iThisObjIdx], "1");
                }
            } else {
                /* there used to be code for iK>=2 in r101 */
                Log(&rLog, LOG_FATAL, "Internal error: split into more than two clusters (got assignment %d)",
                      piKClusterAssignments[iN]);
            }
#if 0
            Log(&rLog, LOG_FORCED_DEBUG, "Setting (*prKMeansResult_p)->ppiObjIndicesPerCluster[%d][%d] = %d",
                      iThisClusterIndex, iThisOffset, iThisObjIdx);
#endif
            (*prKMeansResult_p)->ppiObjIndicesPerCluster[iThisClusterIndex][iThisOffset] = iThisObjIdx;
            piCurObjToUpdate[iThisClusterAssignment]+=1;
        }
        CKFREE(piObjIndicesOfSplitCluster);
#if TRACE
        for (iN=0; iN<(*prKMeansResult_p)->iNClusters+1; iN++) {
            int iObj;
            fprintf(stderr, "DEBUG(%s|%s():%d): Objects in cluster %u: ",
                    __FILE__, __FUNCTION__, __LINE__, iN);
            for (iObj=0; iObj<(*prKMeansResult_p)->piNObjsPerCluster[iN]; iObj++) {
                fprintf(stderr, " %u", (*prKMeansResult_p)->ppiObjIndicesPerCluster[iN][iObj]);
            }
            fprintf(stderr, "\n");
            (void) fflush(stderr);
        }
#endif


#if TIMING
        StopwatchStop(stopwatch);
        StopwatchDisplay(stdout, "Total time after next round in Bisecting-KMeans: ", stopwatch);
#endif


        /* finally: increase number of clusters
         */
        (*prKMeansResult_p)->iNClusters += 1;
        iNRounds += 1;
        CKFREE(piKClusterAssignments);
        CKFREE(pdVectorsInClusterToSplit);

    } /* while */
    INT_QUEUE_DESTROY(&rClusterSplitQueue);


    Log(&rLog, LOG_DEBUG,
         "Bisecting K-means finished after %d rounds (no more clusters to split)",
         iNRounds);
         
#if TIMING
    StopwatchFree(stopwatch);
#endif

    /* @note could use progress/timer */

    CKFREE(pdKClusterCenters);
    CKFREE(piCurObjToUpdate);

    return;
}
/***   end: BisectingKmeans()   ***/



/**
 *
 * @brief From scratch reimplementation of mBed: Blackshields et al.
 * (2010); PMID 20470396.
 *
 * Idea is a follows:
 * - convert sequences into vectors of distances
 * - cluster the vectors using k-means
 * - cluster each of the k clusters using upgma (used cached distances
 *    from above?)
 * - join the sub-clusters to create on tree (use UPGMA on k-means
 *   medoids)
 *
 *
 * @param[out] prMbedTree_p
 * Created upgma tree. will be allocated here. use FreeMuscleTree()
 * to free
 * @param[in] prMSeq
 * Multiple sequences
 * @param[in] iPairDistType
 * Distance measure for pairwise alignments
 * @param[in] pcGuidetreeOut
 * Passed down to GuideTreeUpgma()
 *
 * @note: if the number of sequences is smaller than
 * MAX_ALLOWED_SEQ_PER_PRECLUSTER then there's no need to do the subclustering
 * first. In fact it costs some extra time. However, it's insignificant and
 * for simplicities sake we don't do any special checks

 * 
 * @return Zero on success, non-zero on error
 *
 */
int
Mbed(tree_t **prMbedTree_p, mseq_t *prMSeq, const int iPairDistType,
     const char *pcGuidetreeOut, int iClustersizes, const char *pcClusterFile)
{
    /* number of seeds */
    int iNumSeeds = 0;
    /* seed indices matching prMSeq */
    int *piSeeds = NULL; 
    /* seqs represented as distance vectors */
    double **ppdSeqVec = NULL;
    /* kmeans result */
    bisecting_kmeans_result_t *prKMeansResult = NULL;
    /* distance matrix of kmeans (pre-)cluster centers */
    symmatrix_t *prPreClusterDistmat = NULL;
    /* auxiliary for symmetric matrix output; tree routines etc */
    char **ppcLabels = NULL;
    int iNodeIndex = 0;
    /* mapping of cluster-center tree node indices to corresponding
       cluster */
    int *piClusterToTreeNode = NULL;
    int iClusterIndex = 0;
    int iI, iJ;
    FILE *pfOut = NULL;
    progress_t *prSubClusterDistanceProgress = NULL;
    bool bPrintCR = (rLog.iLogLevelEnabled<=LOG_VERBOSE) ? FALSE : TRUE;
    char **ppcClusterSplits = NULL;

#if FULL_WITHIN_CLUSTER_DISTANCES
    Log(&rLog, LOG_DEBUG, "Computing real distances within subclusters for mBed.");
#else
    Log(&rLog, LOG_DEBUG, "Computing vector distances within subclusters for mBed.");
#endif
    
#if MBED_TIMING
    Stopwatch_t *stopwatch = StopwatchCreate();
#endif

    assert(NULL != prMbedTree_p);
    assert(NULL != prMSeq);


    iNumSeeds = (int) NUMBER_OF_SEEDS(prMSeq->nseqs);
    if (iNumSeeds >= prMSeq->nseqs) {
        /* -1 is condition for RandomUniqueIntArray */
        iNumSeeds = prMSeq->nseqs-1; 
        Log(&rLog, LOG_DEBUG,
             "Automatically determined number of seeds is bigger (or equal) the number of sequences. Will set it to %d",
             iNumSeeds);
    }

       
    /* Turn sequences into vectors of distances to the seeds
     *
     */
    piSeeds = (int *) CKMALLOC(iNumSeeds * sizeof(int));
    if (0 != SeedSelection(piSeeds, iNumSeeds, SEED_SELECTION, prMSeq)) {
        Log(&rLog, LOG_ERROR, "Something went wrong during seed selection for mbed");
        return -1;
    }
    ppdSeqVec = (double **) CKMALLOC(prMSeq->nseqs * sizeof(double *));
    for (iI=0; iI<prMSeq->nseqs; iI++) {
        ppdSeqVec[iI] = (double *) CKMALLOC(iNumSeeds * sizeof(double));
    }
    if (0 != SeqToVec(ppdSeqVec, prMSeq, piSeeds, iNumSeeds, iPairDistType)) {
        Log(&rLog, LOG_ERROR, "Could not convert sequences into vectors for mbed");
        return -1;
    }
    CKFREE(piSeeds);

    
    /* Calculate (pre-)clusters of sequence vectors by applying
     * bisecting kmeans
     *
     */
#if MBED_TIMING
    /* start clock only here, to make sure we don't include pairwise
     * distance computation */
    StopwatchZero(stopwatch);
    StopwatchStart(stopwatch);
#endif

    /* ppcClusterSplits can consume quite a bit of memory,
       however, it is only needed if cluster information is 
       to be written to file. If no pcClusterFile is specified, 
       then ppcClusterSplits does not have to be allocated.
       Use ppcClusterSplits as a flag, initialise to NULL, 
       if NULL is passed into BisectingKmeans then do NOT malloc, 
       if !NULL is passed in then realloc the appropriate memory */
    if (NULL != pcClusterFile){
        ppcClusterSplits = (char **)malloc(sizeof(char*));
    }
    else {
        ppcClusterSplits = NULL;
    }

    BisectingKmeans(&prKMeansResult, prMSeq->nseqs, iNumSeeds, ppdSeqVec,
                    MIN_REQUIRED_SEQ_PER_PRECLUSTER,
                    iClustersizes/*MAX_ALLOWED_SEQ_PER_PRECLUSTER*/, &ppcClusterSplits);
    Log(&rLog, LOG_INFO,
         "mBed created %u cluster/s (with a minimum of %d and a soft maximum of %d sequences each)",
         prKMeansResult->iNClusters,
         MIN_REQUIRED_SEQ_PER_PRECLUSTER,
        iClustersizes/*MAX_ALLOWED_SEQ_PER_PRECLUSTER*/);

    if (NULL != pcClusterFile){ /* print clustering: FS, r275 -> */
        FILE *pfClust = NULL;
        if (NULL == (pfClust = fopen(pcClusterFile, "w"))){
            Log(&rLog, LOG_FATAL, "Could not open file %s for writing", pcClusterFile);
        }
        for (iI = 0; iI < prKMeansResult->iNClusters; iI++) {
            for (iJ=0; iJ<prKMeansResult->piNObjsPerCluster[iI]; iJ++) {
                int iRealIndex = prKMeansResult->ppiObjIndicesPerCluster[iI][iJ];
                fprintf(pfClust, "Cluster %u: object %u has index %u (=seq %s %d~len)\t %s\n",
                        iI, iJ,  iRealIndex, prMSeq->sqinfo[iRealIndex].name, prMSeq->sqinfo[iRealIndex].len, ppcClusterSplits[iRealIndex]);
            }
        }
        fclose(pfClust); pfClust = NULL;
        /* if there is no pcClusterFile then no memory will have been allocated */
        CKFREE(ppcClusterSplits[0]);
        CKFREE(ppcClusterSplits);
    } /* there was a request to write out clustering */

#if PRINT_CLUSTER_DISTRIBUTION
    Log(&rLog, LOG_FORCED_DEBUG, "Bisecting Kmeans returned %d clusters", prKMeansResult->iNClusters);
    for (iI=0; iI<prKMeansResult->iNClusters; iI++) {
#if TRACE
        int iD;
        Log(&rLog, LOG_FORCED_DEBUG, "Diagnostic output for cluster %d follows:", iI);
        fprintf(stderr, "DEBUG(%s|%s():%d): center coordinates =",
                __FILE__, __FUNCTION__, __LINE__);
        for (iD=0; iD<iNumSeeds; iD++) {
            fprintf(stderr, " %f", prKMeansResult->ppdClusterCenters[iI][iD]);
        }
        fprintf(stderr, "\n");
        fflush(stderr);
#endif
        Log(&rLog, LOG_FORCED_DEBUG, "Cluster %d has %d objects assigned",
                  iI, prKMeansResult->piNObjsPerCluster[iI]);
#if TRACE
        for (iJ=0; iJ<prKMeansResult->piNObjsPerCluster[iI]; iJ++) {
            int iRealIndex = prKMeansResult->ppiObjIndicesPerCluster[iI][iJ];
            
            Log(&rLog, LOG_FORCED_DEBUG, "Cluster %u: object %u has index %u (= seq %s)",
                      iI, iJ,  iRealIndex, prMSeq->sqinfo[iRealIndex].name);
        }
#endif
    }
#endif

    
    /* Cluster pre-clusters produced by k-means.
     *
     * Do this by calculating the vector distances of the cluster
     * centers and applying UPGMA.
     *
     * @note could try to force-balance the tree here
     *
     */
    if (0 != NewSymMatrix(&prPreClusterDistmat,
                     prKMeansResult->iNClusters, prKMeansResult->iNClusters)) {
        Log(&rLog, LOG_FATAL, "%s", "Memory allocation for pre-cluster distance-matrix failed");
    }
    for (iI=0; iI<prKMeansResult->iNClusters; iI++) {
        for (iJ=iI+1; iJ<prKMeansResult->iNClusters; iJ++) {
            double dDist;
            dDist = EuclDist(prKMeansResult->ppdClusterCenters[iI],
                             prKMeansResult->ppdClusterCenters[iJ],
                             iNumSeeds);
            SymMatrixSetValue(prPreClusterDistmat, iI, iJ, dDist);
            /* Log(&rLog, LOG_FORCED_DEBUG, "Euclidean distance between clusters %d and %d = %f",
               iI, iJ, dDist); */
        }
    }
    

    /* labels needed for the guide tree building routine only */
    ppcLabels = (char **) CKMALLOC(prKMeansResult->iNClusters * sizeof(char*));
    for (iI=0; iI<prKMeansResult->iNClusters; iI++) {
        ppcLabels[iI] = (char *) CKMALLOC(32 * sizeof(char));
        (void) snprintf(ppcLabels[iI], 32, "Subcluster-%u", iI);
    }
#if TRACE
    Log(&rLog, LOG_FORCED_DEBUG, "%s", "Distance matrix for pre-cluster centers:");
    SymMatrixPrint(prPreClusterDistmat, ppcLabels, NULL, FALSE);
#endif

    GuideTreeUpgma(prMbedTree_p,
                   ppcLabels, prPreClusterDistmat, NULL);

    for (iI=0; iI<prKMeansResult->iNClusters; iI++) {
        CKFREE(ppcLabels[iI]);
    }
    CKFREE(ppcLabels);
#if TRACE
    Log(&rLog, LOG_FORCED_DEBUG, "%s", "Cluster-center guide-tree:");
    LogTree(*prMbedTree_p, stderr);
#endif



    /* Now replace each leaf in the pre-cluster-center tree
     * appropriately, i.e. with the corresponding sub-cluster.
     *
     * For each leaf/sub-cluster, create a distance matrix for the
     * corresponding sequences. Use distances between vectors as
     * approximated distances between sequences. Then create a
     * guide-tree by applying UPGMA.
     *
     */

    /* Get a mapping of (pre)cluster number and leaf-node indices in
     * the cluster-center tree. We can add trees to prMbedTree_p
     * because AppendTrees() guarantees that no other than the node to
     * append to changes.
     */
    piClusterToTreeNode = (int*)
        CKMALLOC(prKMeansResult->iNClusters * sizeof(int));
    iNodeIndex = FirstDepthFirstNode(*prMbedTree_p);
    do {
        if (IsLeaf(iNodeIndex, *prMbedTree_p)) {
            int iLeafId = GetLeafId(iNodeIndex, *prMbedTree_p);
            piClusterToTreeNode[iLeafId] = iNodeIndex;
        }
        iNodeIndex =  NextDepthFirstNode(iNodeIndex, *prMbedTree_p);
    } while (NULL_NEIGHBOR != iNodeIndex);




    /* Now step through all the leafs and replace them with the
     * corresponding sub-trees
     */
    NewProgress(&prSubClusterDistanceProgress, LogGetFP(&rLog, LOG_INFO),
                "Distance calculation within sub-clusters", bPrintCR);
    /* for each cluster */
    for (iClusterIndex=0;
         iClusterIndex < prKMeansResult->iNClusters; iClusterIndex++) {
        /* distance matrix for the sub-cluster */
        symmatrix_t *prWithinClusterDistances = NULL;
        int iNSeqInCluster;
        tree_t *prSubClusterTree = NULL;
        
        ProgressLog(prSubClusterDistanceProgress,
                    iClusterIndex, prKMeansResult->iNClusters, FALSE);

#if FULL_WITHIN_CLUSTER_DISTANCES
        mseq_t *prSubClusterMSeq;
        int iPairDistType;
        int iSeqIndex;
        int iOldLogLevel;

        Log(&rLog, LOG_DEBUG, 
            "%s\n", "Calling new Mbed use makes only sense if nseq>MAX_ALLOWED_SEQ_PER_PRECLUSTER");

        if (TRUE == prMSeq->aligned)  {
            if (SEQTYPE_PROTEIN == prMSeq->seqtype){
                iPairDistType = PAIRDIST_SQUIDID_KIMURA;
            }
            else {
                iPairDistType = PAIRDIST_SQUIDID;
            }
        } else {
            iPairDistType = PAIRDIST_KTUPLE;
        }
#endif
            
        iNSeqInCluster = prKMeansResult->piNObjsPerCluster[iClusterIndex];
#if TRACE
        Log(&rLog, LOG_FORCED_DEBUG, "#seqs in subcluster no %d = %d",
                  iClusterIndex, iNSeqInCluster);
#endif

#if FULL_WITHIN_CLUSTER_DISTANCES
        /* create an mseq structure for sequences in this cluster
         * don't need most of the members (e.g. orig_seq) but copy
         * them anyway for the sake of completeness
         */
        NewMSeq(&prSubClusterMSeq);
    
        prSubClusterMSeq->nseqs = iNSeqInCluster;
        prSubClusterMSeq->seqtype = prMSeq->seqtype;
        if (NULL!=prMSeq->filename) {
            prSubClusterMSeq->filename = CkStrdup(prMSeq->filename);
        }
        prSubClusterMSeq->aligned = prMSeq->aligned; /* FS, r252 */
        prSubClusterMSeq->seq =  (char **)
            CKMALLOC(prSubClusterMSeq->nseqs * sizeof(char *));
        prSubClusterMSeq->orig_seq =  (char **)
            CKMALLOC(prSubClusterMSeq->nseqs * sizeof(char *));
        prSubClusterMSeq->sqinfo =  (SQINFO *)
            CKMALLOC(prSubClusterMSeq->nseqs * sizeof(SQINFO));
        
        for (iSeqIndex=0; iSeqIndex<iNSeqInCluster; iSeqIndex++) {
            int iRealSeqIndex = prKMeansResult->ppiObjIndicesPerCluster[iClusterIndex][iSeqIndex];
            prSubClusterMSeq->seq[iSeqIndex] = CkStrdup(prMSeq->seq[iRealSeqIndex]);
            prSubClusterMSeq->orig_seq[iSeqIndex] = CkStrdup(prMSeq->orig_seq[iRealSeqIndex]);
            SeqinfoCopy(&prSubClusterMSeq->sqinfo[iSeqIndex], &prMSeq->sqinfo[iRealSeqIndex]);
#if TRACE
            Log(&rLog, LOG_DEBUG, "seq no %d in cluster %d is %s (real index = %d)",
                iSeqIndex, iClusterIndex, prSubClusterMSeq->sqinfo[iSeqIndex].name,
                iRealSeqIndex);
#endif
        }
#endif

        
        /* Create a distance matrix for this sub-cluster
         * (prWithinClusterDistances) by using the vector distances or
         * ktuple distances.
         *
         * Then apply UPGMA to get a subcluster tree
         * (prSubClusterTree) and append created tree to the
         * pre-cluster-tree (prMbedTree_p)
         *
         */
#if FULL_WITHIN_CLUSTER_DISTANCES
        
        iOldLogLevel = rLog.iLogLevelEnabled;
        rLog.iLogLevelEnabled = LOG_WARN;
        /* compute distances, but be quiet */
        /* 4th argument (FALSE) is bPercID, in mBed mode never use percent-identities */
        if (PairDistances(&prWithinClusterDistances, prSubClusterMSeq, iPairDistType, FALSE, 
                          0, prSubClusterMSeq->nseqs, 0, prSubClusterMSeq->nseqs,
                          NULL, NULL)) {
            Log(&rLog, LOG_ERROR, "Couldn't compute pair distances");
            return -1;
        }
        rLog.iLogLevelEnabled = iOldLogLevel;

#if COMPUTE_WITHIN_SUBCLUSTER_AVERAGE
        {
            double dSum = 0.0;
            for (iI=0; iI<prSubClusterMSeq->nseqs; iI++) {
                for (iJ=iI+1; iJ<prSubClusterMSeq->nseqs; iJ++) {
                    dSum += SymMatrixGetValue(prWithinClusterDistances, iI, iJ);
                }
            }
            Log(&rLog, LOG_FORCED_DEBUG,
                "mean pair-wise distance within subcluster %d of %d = %f", 
                iClusterIndex, prKMeansResult->iNClusters, dSum/prSubClusterMSeq->nseqs);
        }
#endif

#else
        if (NewSymMatrix(&prWithinClusterDistances, iNSeqInCluster, iNSeqInCluster)!=0) {
            Log(&rLog, LOG_FATAL, "%s", "Memory allocation for disparity matrix failed");
        }
        for (iI=0; iI<iNSeqInCluster; iI++) {
            int iRealIndexI = prKMeansResult->ppiObjIndicesPerCluster[iClusterIndex][iI];

            for (iJ=iI+1; iJ<iNSeqInCluster; iJ++) {
                int iRealIndexJ = prKMeansResult->ppiObjIndicesPerCluster[iClusterIndex][iJ];
                double dist;

                /* Log(&rLog, LOG_FORCED_DEBUG, "Cluster %d: compute distance between %d:%s and %d:%s",
                   iClusterIndex, i, prMSeq->sqinfo[iRealIndexI].name,
                   iJ, prMSeq->sqinfo[iRealIndexJ].name); */

                if (1 == USE_EUCLIDEAN_DISTANCE) {
                    dist = EuclDist(ppdSeqVec[iRealIndexI],
                                    ppdSeqVec[iRealIndexJ], iNumSeeds);
                } else {
                    dist = CosDist(ppdSeqVec[iRealIndexI],
                                   ppdSeqVec[iRealIndexJ], iNumSeeds);
                }
                SymMatrixSetValue(prWithinClusterDistances, iI, iJ, dist);
            }
        }
#endif

        /* labels needed for the guide tree building routine only */
        ppcLabels = (char**) CKMALLOC(iNSeqInCluster * sizeof(char*));
        for (iI=0; iI<iNSeqInCluster; iI++) {
#if FULL_WITHIN_CLUSTER_DISTANCES
            ppcLabels[iI] = prSubClusterMSeq->sqinfo[iI].name;
#else       
            int iRealIndex = prKMeansResult->ppiObjIndicesPerCluster[iClusterIndex][iI];
            ppcLabels[iI] = prMSeq->sqinfo[iRealIndex].name;
#endif
        }        
#if TRACE
        Log(&rLog, LOG_FORCED_DEBUG, "Distance matrix for seqs within sub cluster %d/%d",
                  iClusterIndex, prKMeansResult->iNClusters);
        SymMatrixPrint(prWithinClusterDistances, ppcLabels, NULL, FALSE);
#endif

        
        GuideTreeUpgma(&prSubClusterTree, ppcLabels,
                       prWithinClusterDistances, NULL);

        CKFREE(ppcLabels); /* don't free members, they just point */
#if 0
        Log(&rLog, LOG_FORCED_DEBUG, "Cluster %d guide-tree:", iClusterIndex);
        LogTree(prSubClusterTree);
#endif

        
        /* The guide tree id's (that point to the sequences) now start
         * from 0, i.e. the association with the prMSeq numbering is
         * broken and fixed in the following
         */
        for (iNodeIndex = 0; iNodeIndex < (int)GetNodeCount(prSubClusterTree); iNodeIndex++) {
            if (IsLeaf(iNodeIndex, prSubClusterTree)) {
                int iLeafId = GetLeafId(iNodeIndex, prSubClusterTree);
                int iRealId = prKMeansResult->ppiObjIndicesPerCluster[iClusterIndex][iLeafId];
#if 0
                Log(&rLog, LOG_FORCED_DEBUG, "Correcting leaf node %d which has (wrong) id %d and name %s to id %d (prMSeq name %s)",
                          iNodeIndex, iLeafId,
                          GetLeafName(iNodeIndex, prSubClusterTree),
                          iRealId, prMSeq->sqinfo[iRealId].name);
#endif
                SetLeafId(prSubClusterTree, iNodeIndex, iRealId);
            }
        }


        /* Append the newly created tree (prSubClusterTree) to the
         * corresponding node index of prMbedTree_p.
         */
#if TRACE
        Log(&rLog, LOG_FORCED_DEBUG, "Will join trees at leaf node %d = %s",
                  piClusterToTreeNode[iClusterIndex],
                  GetLeafName(piClusterToTreeNode[iClusterIndex], *prMbedTree_p));
#endif

        AppendTree(*prMbedTree_p,
                   piClusterToTreeNode[iClusterIndex],
                   prSubClusterTree);
        /* Note: piClusterToTreeNode is still valid, because
         * AppendTrees() guarantees that no other than the node to
         * append to changes. */

#if 0
        Log(&rLog, LOG_FORCED_DEBUG, "%s", "prMbedTree_p after cluster %d has appended:", iClusterIndex);
        LogTree(*prMbedTree_p);

        if (0) {
            char fname[] = "mbed-joined-tree.dnd";
            FILE *pfOut;
            if (NULL == (pfOut = fopen(fname, "w"))) {
                Log(&rLog, LOG_FATAL, "Couldn't open %s for writing", fname);
            }
            MuscleTreeToFile(pfOut, *prMbedTree_p);
            Log(&rLog, LOG_FORCED_DEBUG, "Joined tree written to %s", fname);
            fclose(pfOut);
            Log(&rLog, LOG_FATAL, "DEBUG EXIT");
        }
#endif

        /* cleanup
         */
        FreeMuscleTree(prSubClusterTree);
        FreeSymMatrix(&prWithinClusterDistances);
#if FULL_WITHIN_CLUSTER_DISTANCES
        FreeMSeq(&prSubClusterMSeq);
#endif
    } /* end for each cluster */
    ProgressDone(prSubClusterDistanceProgress);
    FreeProgress(&prSubClusterDistanceProgress);

    
    if (NULL != pcGuidetreeOut) {
        if (NULL == (pfOut = fopen(pcGuidetreeOut, "w"))) {
            Log(&rLog, LOG_ERROR, "Couldn't open %s for writing", pcGuidetreeOut);
        } else {
            MuscleTreeToFile(pfOut, *prMbedTree_p);
            Log(&rLog, LOG_INFO, "Guide tree written to %s", pcGuidetreeOut);
            (void) fclose(pfOut);
        }
    }

    
    /* cleanup
     *
     */
#if MBED_TIMING
    StopwatchStop(stopwatch);
    StopwatchDisplay(stdout, "mBed time (without pairwise distance computation): ", stopwatch);
    StopwatchFree(stopwatch);
#endif

    FreeKMeansResult(&prKMeansResult);
    FreeSymMatrix(&prPreClusterDistmat);
    for (iI=0; iI<prMSeq->nseqs; iI++) {
        CKFREE(ppdSeqVec[iI]);
    }
    CKFREE(ppdSeqVec);
    CKFREE(piClusterToTreeNode);

#ifndef NDEBUG
    TreeValidate(*prMbedTree_p);
#endif
    
    return 0;
}
/***   end: Mbed()   ***/
