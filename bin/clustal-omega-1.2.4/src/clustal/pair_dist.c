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
 *  RCS $Id: pair_dist.c 301 2016-06-13 13:32:55Z fabian $
 */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <stdlib.h>
#include <ctype.h>
#include <assert.h>
#include <time.h>

/* only neededfor iNumberOfThreads */
#include "clustal-omega.h"

#include "ktuple_pair.h"
#include "pair_dist.h"
#include "progress.h"
#include "util.h"

/* Made iend/jend const unsigned long int (originally just int), FS, 2016-04-04
 */


/* Up to rev 173 we had a USE_SYM_KTUPLE switch implemented here. When active
 * ktuple distances were computed twice for each pair and averaged. Idea was
 * to avoid assymmetries in the pairwise scores (score(a, b) is often not the
 * same as score(b, a)). Results on BAliBASE indicate that this is overkill:
 *
 * r92_default            core columns: avg-sp=0.800656 avg-tc=0.47711  (of total 218)
 * r93-mod--norm-ktuple/  core columns: avg-sp=0.800656 avg-tc=0.47711  (of total 218)
 * r93-mod--sym-ktuple/   core columns: avg-sp=0.801083 avg-tc=0.476544 (of total 217)
 * r93-mod--rand-ktuple-1 core columns: avg-sp=0.799289 avg-tc=0.468028 (of total 218)
 * r93-mod--rand-ktuple-2 core columns: avg-sp=0.801654 avg-tc=0.47659  (of total 217)
 * r93-mod--rand-ktuple-3 core columns: avg-sp=0.800234 avg-tc=0.474908 (of total 218)
 * r93-mod--rand-ktuple-4 core columns: avg-sp=0.800573 avg-tc=0.476514 (of total 218)
 * r93-mod--rand-ktuple-5 core columns: avg-sp=0.799679 avg-tc=0.468716 (of total 218)
 *
 */

static double
KimuraCorrection(double frac_id);

static int
SquidIdPairDist(symmatrix_t *tmat, mseq_t *mseq,
                int istart, const unsigned long int iend,
                int jstart, const unsigned long int jend,
                bool use_KimuraCorrection, progress_t *prProgress,
                unsigned long int *ulStepNo, unsigned long int ulTotalStepNo);

/* Taken from Muscle's msadistkimura.cpp */
static int DAYHOFF_PAMS[]={
  195,   /* 75.0% observed d; 195 PAMs estimated = 195% estimated d */
  196,   /* 75.1% observed d; 196 PAMs estimated */
                  197,    198,    199,    200,    200,    201,    202,  203,
  204,    205,    206,    207,    208,    209,    209,    210,    211,  212,
  213,    214,    215,    216,    217,    218,    219,    220,    221,  222,
  223,    224,    226,    227,    228,    229,    230,    231,    232,  233,
  234,    236,    237,    238,    239,    240,    241,    243,    244,  245,
  246,    248,    249,    250,    /* 250 PAMs = 80.3% observed d */
                                  252,    253,    254,    255,    257,  258,
  260,    261,    262,    264,    265,    267,    268,    270,    271,  273,
  274,    276,    277,    279,    281,    282,    284,    285,    287,  289,
  291,    292,    294,    296,    298,    299,    301,    303,    305,  307,
  309,    311,    313,    315,    317,    319,    321,    323,    325,  328,
  330,    332,    335,    337,    339,    342,    344,    347,    349,  352,
  354,    357,    360,    362,    365,    368,    371,    374,    377,  380,
  383,    386,    389,    393,    396,    399,    403,    407,    410,  414,
  418,    422,    426,    430,    434,    438,    442,    447,    451,  456,
  461,    466,    471,    476,    482,    487,    493,    498,    504,  511,
  517,    524,    531,    538,    545,    553,    560,    569,    577,  586,
  595,    605,    615,    626,    637,    649,    661,    675,    688,  703,
  719,    736,    754,    775,    796,    819,    845,    874,    907,  945,
         /* 92.9% observed; 945 PAMs */
  988    /* 93.0% observed; 988 PAMs */
};
static int DAYHOFF_TABLE_ENTRIES = sizeof(DAYHOFF_PAMS)/sizeof(DAYHOFF_PAMS[0]);



/**
 *
 * @brief Compute Kimura corrected distance.
 *
 * Original Muscle documentation following:
 * """
 * This is defined to be:
 *     log_e(1 - p - p*p/5)
 * where p is the fraction of residues that differ, i.e.:
 *     p = (1 - fractional_conservation)
 * This measure is infinite for p = 0.8541 and is considered
 * unreliable for p >= 0.75 (according to the ClustalW docs).
 * ClustalW uses a table lookup for values > 0.75. The following table
 * was copied from the ClustalW file dayhoff.h.
 * """
 *
 * @note copied from Muscle's msadistkimura.cpp:KimuraDist()
 *
 * @warning For protein only (uses Dayhoff substitution parameters)
 *
 * @param[in] p
 * distance, e.g. 1.0 - fractional/relative identity
 *
 * @return The Kimura corrected distance
 *
 */
double
KimuraCorrection(double p)
{
    int table_index;

    /* Typical case: use Kimura's empirical formula */
    if (p < 0.75)
        return -log(1 - p - (p*p)/5);

    /* Per ClustalW, return 10.0 for anything over 93% */
    if (p > 0.93)
        return 10.0;

    /* If 0.75 >= p <= 0.93, use table lookup */
    table_index = (int) ((p - 0.75)*1000 + 0.5);
    if (table_index < 0 || table_index >= DAYHOFF_TABLE_ENTRIES)
        Log(&rLog, LOG_FATAL, "Internal error in %s:%s", __FILE__, __FUNCTION__);

    return DAYHOFF_PAMS[table_index] / 100.0;
}
/***   end: KimuraCorrection()   ***/




/**
 * @brief Compute distances between all aligned sequence pairs using
 * squid's PairwiseIdentity, which is: idents / MIN(len1, len2)
 *
 * @param[out] tmat
 * Where to store the computed distances
 * @param[in] mseq
 * The aligned sequences
 * @param[in] istart
 * For distances [i][j] i>=istart, i<j
 * @param[in] iend
 * For distances [i][j] i<iend, i<j
 * @param[in] jstart
 * For distances [i][j] j>=jstart, i<j
 * @param[in] jend
 * For distances [i][j] i<j<jend, i<j
 * @param[in] use_kimura
 * Use Kimura corrected values (Proteins only)
 *
 * @return Non-zero on error
 *
 */
int
SquidIdPairDist(symmatrix_t *tmat, mseq_t *mseq,
                int istart, const unsigned long int iend,
                int jstart, const unsigned long int jend,
                bool use_kimura, progress_t *prProgress,
                unsigned long int *ulStepNo, unsigned long int ulTotalStepNo)
{
    int i, j; /* aux */
    /* progress_t *prProgress; */
    bool bPrintCR = (rLog.iLogLevelEnabled<=LOG_VERBOSE) ? FALSE : TRUE;
    /* unsigned long int ulStepNo;
    unsigned long ulTotalStepNo; */

    assert(NULL != tmat);
    assert(NULL != mseq);

    if (TRUE != mseq->aligned) {
        Log(&rLog, LOG_ERROR, "Sequences need to be aligned (%s)", __FUNCTION__);
        return -1;
    }
    if (SEQTYPE_PROTEIN != mseq->seqtype && TRUE == use_kimura) {
        Log(&rLog, LOG_WARN, "Using Kimura distance corretion which includes Dayhoff substitution table lookup for non-protein sequences");
    }

    NewProgress(&prProgress, LogGetFP(&rLog, LOG_INFO),
                "Pairwise distance calculation progress", bPrintCR);
    /* estimation of total number of steps (if istart and jstart are
     * both 0)
     */
    /* ulTotalStepNo = iend*jend - iend*iend/2 + iend/2;
    ulStepNo = 0; */
    /*LOG_DEBUG("istart=%d iend=%d jstart=%d jend=%d", istart, iend, jstart, jend);*/
    for (i=istart; i<iend; ++i) {
        /* by definition a sequence compared to itself should give a
           score of 0 */
        SymMatrixSetValue(tmat, i, i, 0.0);
#ifdef HAVE_OPENMP
        #pragma omp critical(squidid)
#endif
        {
            ProgressLog(prProgress, *ulStepNo, ulTotalStepNo, FALSE);
        }

        for (j=MAX(i+1, jstart); j<jend; ++j) {
            float dist;
            dist = 1.0 - PairwiseIdentity(mseq->seq[i], mseq->seq[j]);


#ifdef HAVE_OPENMP
        #pragma omp atomic
#endif
            (*ulStepNo)++;
            /*LOG_DEBUG("%d:%d raw dist = %f", i, j, dist);*/
            if (use_kimura) {
                dist = KimuraCorrection(dist);
                /*LOG_DEBUG("cor dist = %f", dist);*/
            }
            SymMatrixSetValue(tmat, i, j, dist);
#ifdef HAVE_OPENMP
            #pragma omp critical(squidid)
#endif
            {
                Log(&rLog, LOG_DEBUG, "Aligned distance for sequence pair %d:%d= %lg",
                     i+1, j+1, dist);
            }
        }
    }

    return 0;
}
/***   end: SquidIdPairDist()   ***/




/**
 * @brief compute or read precomputed distances for given sequences
 *
 * @param[out] distmat
 * Distances will be written to this matrix. will be allocated here as
 * well. Caller must free with FreeSymMatrix()
 * @param[in] mseq
 * Distances will be computed for these sequences
 * @param[in] pairdist_type
 * Type of pairwise distance comparison
 * @param[in] fdist_in
 * If not NULL, sequences will be written from this file instead of
 * computing them
 * @param[in] istart
 * Compute distances for sequences i:j, i>=istart, i<j.
 * Usually 0.
 * @param[in] iend
 * Compute distances for sequences i:j, i<iend, i<j
 * Usually mseq->nseqs.
 * @param[in] jstart
 * Compute distances for sequences i:j, j>=jstart, i<j
 * Usually 0.
 * @param[in] jend
 * Compute distances for sequences i:j, j<iend, i<j
 * Usually mseq->nseqs.
 * @param[in] fdist_out
 * If not NULL, distances will be written to this files
 *
 *
 */
int
PairDistances(symmatrix_t **distmat, mseq_t *mseq, int pairdist_type, bool bPercID, 
              int istart, const unsigned long int iend,
              int jstart, const unsigned long int jend,
              char *fdist_in, char *fdist_out)
{
    int uSeqIndex;
    unsigned long int ulStepNo = 0, ulTotalStepNo; /* DD: moved from SquidIdPairDist so progress bar works multithreaded */
    int iChunk, iChunkStart, iChunkEnd;
    int iChunkStarts[iNumberOfThreads];
    int iChunkEnds[iNumberOfThreads];
    progress_t *prProgress = NULL;
    int iSquidSuccess = 0;
    bool bPrintCR = (rLog.iLogLevelEnabled<=LOG_VERBOSE) ? FALSE : TRUE;

    assert(NULL!=distmat);
    assert(NULL!=mseq);
    assert(istart<iend);
    assert(jstart<jend);


    /* compute pairwise distances or read from file
     *
     */
#if 0
#include "random-dist.h"
#else
    if (NULL != fdist_in) {
        Log(&rLog, LOG_WARN,
            "Please use distance matrix input only, if you know exactly what you're doing!");
        if (SymMatrixRead(fdist_in, distmat, mseq)) {
            Log(&rLog, LOG_FATAL, "%s", "Reading distance matrix failed");
        }

    } else {

        if (NewSymMatrix(distmat, iend, jend)!=0) {
            Log(&rLog, LOG_FATAL, "%s", "Memory allocation for distance matrix failed");
        }

        /* break into chunks, one for each thread
           matrix is a triangle, not a square
           hence making even chunk sizes is slightly fiddlier
           */
        ulTotalStepNo = iend*jend - iend*iend/2 + iend/2;
        
        /* FIXME: can get rid of iChunkStart, iChunkEnd now that we're using the arrays */
        iChunkStart = iend;
        for(iChunk = 0; iChunk <= iNumberOfThreads; iChunk++)
        {
            iChunkEnd = iChunkStart;
            if (iChunk == iNumberOfThreads - 1){
                iChunkStart = 0;
            }
            else if (iend == jend){
                iChunkStart = iend - ((double)(iend - istart) * sqrt(((double)iChunk + 1.0)/(double)iNumberOfThreads));
            }
            else {
                iChunkStart = iend - (iend - istart) * (iChunk + 1) / (double)(iNumberOfThreads);
            }
            iChunkStarts[iChunk] = iChunkStart;
            iChunkEnds[iChunk] = iChunkEnd;
            /*printf("%s:%d: C=%d, ie=%d, is=%d, je=%d, js=%d, Cstart=%d, Cend=%d, diff=%d\n", 
                   __FILE__, __LINE__, iChunk, iend, istart, jend, jstart, iChunkStart, iChunkEnd, iChunkEnd-iChunkStart);*/
        }

        if (PAIRDIST_KTUPLE == pairdist_type) {

            Log(&rLog, LOG_INFO, "Calculating pairwise ktuple-distances...");
            
            NewProgress(&prProgress, LogGetFP(&rLog, LOG_INFO),
                        "Ktuple-distance calculation progress", bPrintCR); 
#ifdef HAVE_OPENMP
            #pragma omp parallel for private(iChunk) schedule(dynamic)
#endif
            for(iChunk = 0; iChunk < iNumberOfThreads; iChunk++)
            {
                KTuplePairDist((*distmat), mseq, iChunkStarts[iChunk], 
                    iChunkEnds[iChunk], jstart, jend, NULL, prProgress, 
                    &ulStepNo, ulTotalStepNo);
            }

#if 0
            printf("total ops %d\n", ulStepNo);
#endif
            /* old format:
            KTuplePairDist((*distmat), mseq,
                istart, iend,
                jstart, jend, NULL); */

        } else if (PAIRDIST_SQUIDID == pairdist_type) {
            Log(&rLog, LOG_INFO, "Calculating pairwise aligned identity distances...");

            NewProgress(&prProgress, LogGetFP(&rLog, LOG_INFO),
                        "Pairwise identity calculation progress", bPrintCR);
#ifdef HAVE_OPENMP
            #pragma omp parallel for private(iChunk) schedule(dynamic)
#endif
            for(iChunk = 0; iChunk < iNumberOfThreads; iChunk++)
            {
                 iSquidSuccess = SquidIdPairDist((*distmat), mseq,
                                    iChunkStarts[iChunk], iChunkEnds[iChunk],
                                    jstart, jend, FALSE, prProgress,
                                    &ulStepNo, ulTotalStepNo);
            }
            if(iSquidSuccess != 0)
                return -1;

        } else if (PAIRDIST_SQUIDID_KIMURA == pairdist_type) {
            Log(&rLog, LOG_INFO, "Calculating Kimura-corrected pairwise aligned identity distances...");
            NewProgress(&prProgress, LogGetFP(&rLog, LOG_INFO),
                        "Pairwise identity calculation progress", bPrintCR);
#ifdef HAVE_OPENMP
            #pragma omp parallel for private(iChunk) schedule(dynamic)
#endif
            for(iChunk = 0; iChunk < iNumberOfThreads; iChunk++)
            {
                iSquidSuccess = SquidIdPairDist((*distmat), mseq,
                                    iChunkStarts[iChunk], iChunkEnds[iChunk],
                                    jstart, jend, TRUE, prProgress,
                                    &ulStepNo, ulTotalStepNo);
            }
            if(iSquidSuccess != 0)
                return -1;
        } else {
            Log(&rLog, LOG_FATAL, "INTERNAL ERROR: don't know about pairdist_type %d",
                  pairdist_type);
        }
    }
#endif /* random/proper distance calculation */

    
    /* optional printing of matrix to file
     */
    if (NULL != fdist_out) {
        /* need a copy of sequence names for printing */
        char **names;
        names = (char **)CKMALLOC(mseq->nseqs * sizeof(char*));
        for (uSeqIndex=0; uSeqIndex<mseq->nseqs; uSeqIndex++) {
            names[uSeqIndex] = mseq->sqinfo[uSeqIndex].name;
        }

        SymMatrixPrint((*distmat), names, fdist_out, bPercID);

        Log(&rLog, LOG_INFO, "Pairwise distance matrix written to %s",
             fdist_out);
        CKFREE(names);
    }

#if 0
#include "distance-distrib.h" 
#endif

    if (NULL != prProgress) {
        ProgressDone(prProgress);
        FreeProgress(&prProgress);
    }
    
    return 0;
}
/***   end: PairDistances()   ***/



