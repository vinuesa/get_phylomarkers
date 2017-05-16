/* -*- mode: c; tab-width: 4; c-basic-offset: 4; indent-tabs-mode: nil -*- */

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
 *  RCS $Id: hhalign_wrapper.c 306 2016-06-13 13:49:04Z fabian $
 */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <ctype.h>
#include <stdbool.h>

#include "seq.h"
#include "tree.h"
#include "progress.h"
#include "hhalign/general.h"
#include "hhalign/hhfunc.h"
#include "hhalign/hhalign.h"

/* up to this level (from leaf) will background HMM info be applied */
#define APPLY_BG_HMM_UP_TO_TREE_DEPTH 10

#define TIMING 0

#define TRACE 0

/**
 * @brief FIXME
 * 
 * @note prHalignPara has to point to an already allocated instance
 *
 */
void SetDefaultHhalignPara(hhalign_para *prHhalignPara)
{
    prHhalignPara->iMacRamMB = 8000;  /* default|give just under 8GB to MAC algorithm, FS, 2016-04-18, went from 2GB to 8GB */
    prHhalignPara->bIsDna = false; /* protein mode unless we say otherwise */	
    prHhalignPara->bIsRna = false;	
    prHhalignPara->pca = -UNITY;
    prHhalignPara->pcb = -UNITY;
    prHhalignPara->pcc = -UNITY;
    prHhalignPara->pcw = -UNITY;
    prHhalignPara->gapb = -UNITY;
    prHhalignPara->gapd = -UNITY;
    prHhalignPara->gape = -UNITY;
    prHhalignPara->gapf = -UNITY;
    prHhalignPara->gapg = -UNITY;
    prHhalignPara->gaph = -UNITY;
    prHhalignPara->gapi = -UNITY;
    prHhalignPara->pcaV = -UNITY;
    prHhalignPara->pcbV = -UNITY;
    prHhalignPara->pccV = -UNITY;
    prHhalignPara->pcwV = -UNITY;
    prHhalignPara->gapbV = -UNITY;
    prHhalignPara->gapdV = -UNITY;
    prHhalignPara->gapeV = -UNITY;
    prHhalignPara->gapfV = -UNITY;
    prHhalignPara->gapgV = -UNITY;
    prHhalignPara->gaphV = -UNITY;
    prHhalignPara->gapiV = -UNITY;

}
/*** end: SetDefaultHhalignPara()  ***/



/**
 * @brief get rid of unknown residues
 *
 * @note HHalignWrapper can be entered in 2 different ways: (i) all
 * sequences are un-aligned (ii) there are 2 (aligned) profiles.  in
 * the un-aligned case (i) the sequences come straight from Squid,
 * that is, they have been sanitised, all non-alphabetic residues 
 * have been rendered as X's. In profile mode (ii) one profile may 
 * have been produced internally. In that case residues may have 
 * been translated back into their 'native' form, that is, they may
 * contain un-sanitised residues. These will cause trouble  during
 * alignment
 * FS, r213->214
 */
void
SanitiseUnknown(mseq_t *mseq)
{

    int iS; /* iterator for sequence */
    /*int iR;*/ /* iterator for residue */
    /*int iLen;*/ /* length of sequence */
    char *pcRes = NULL;


    for (iS = 0; iS < mseq->nseqs; iS++){

        for (pcRes = mseq->seq[iS]; '\0' != *pcRes; pcRes++){

            if (isgap(*pcRes)){
                /* don't like MSF gap characters ('~'), 
                   sanitise them (and '.' and ' '); FS, r258 -> r259 */
                *pcRes = '-';
                continue;
            }

            if (mseq->seqtype==SEQTYPE_PROTEIN) {
                if (NULL == strchr(AMINO_ALPHABET, toupper(*pcRes))) {
                    *pcRes = AMINOACID_ANY;
                }
            } else if (mseq->seqtype==SEQTYPE_DNA) {
                if (NULL == strchr(DNA_ALPHABET, toupper(*pcRes))) {
                    *pcRes = NUCLEOTIDE_ANY;
                }
            } else if (mseq->seqtype==SEQTYPE_RNA) {
                if (NULL == strchr(RNA_ALPHABET, toupper(*pcRes))) {
                    *pcRes = NUCLEOTIDE_ANY;
                }
            }

        } /* !EO String */

    } /* 0 <= iS < mseq->nseqs */

    return;

} /*** end:  SanitiseUnknown()  ***/

/**
 * @brief translate unknown residues back to ambiguity codes;
 * hhalign translates ambiguity codes (B,Z) into unknown residue (X).
 * we still have the original (un-aligned) residue information,
 * by iterating along the original and aligned sequences we can
 * reconstruct where codes have been changed and restore them
 * to their original value
 *
 * @param[in,out] mseq
 * sequence/profile data, mseq->seq [in,out] is changed to conform
 * with mseq->orig_seq [in]
 *
 */
void
TranslateUnknown2Ambiguity(mseq_t *mseq)
{

    int iS; /* iterator for sequence */
    int iR, iRo; /* iterator for residue (original) */
    int iChange, iCase, iAmbi; /* counts how many replacements */
    static int siOffset = 'a' - 'A';

    for (iS = 0; iS < mseq->nseqs; iS++){

        iR = iRo = 0;
        iChange = iCase = iAmbi = 0;

        while(('\0' != mseq->seq[iS][iR]) &&
              ('\0' != mseq->orig_seq[iS][iRo])) {

            /* skip gaps in aligned sequences */
            while(isgap(mseq->seq[iS][iR])) {
                iR++;
            } /* was gap in unaligned seq
               * this should probably not happen */
            while(isgap(mseq->orig_seq[iS][iRo])) {
                iRo++;
            } /* was gap in aligned seq */

            /* check if we reached the end of the sequence after
             * skipping the gaps */
            if ( ('\0' == mseq->seq[iS][iR]) ||
                 ('\0' == mseq->orig_seq[iS][iRo]) ){
                break;
            }


            if (mseq->seq[iS][iR] != mseq->orig_seq[iS][iRo]){
                /* FIXME: count replacements, discard case changes */
                iChange++;
                if ( (mseq->seq[iS][iR] == mseq->orig_seq[iS][iRo]+siOffset) ||
                     (mseq->seq[iS][iR] == mseq->orig_seq[iS][iRo]-siOffset) ){
                    iCase++;
                }
                else {
                    iAmbi++;
                }
#if TRACE
                Log(&rLog, LOG_FORCED_DEBUG, "seq=%d, pos=(%d:%d), (%c:%c)",
                     iS, iR, iRo,
                     mseq->seq[iS][iR], mseq->orig_seq[iS][iRo]);
#endif
                mseq->seq[iS][iR] = mseq->orig_seq[iS][iRo];
            }

            iR++;
            iRo++;

        } /* !EO seq */

        Log(&rLog, LOG_DEBUG,
             "in seq %d re-translated %d residue codes (%d true, %d case)",
             iS, iChange, iAmbi, iCase);

        /* assert that both sequences (un/aligned) have terminated */
        /* skip gaps in aligned sequences */
        while(isgap(mseq->seq[iS][iR])) {
            iR++;
        } /* was gap in unaligned seq
           * this should probably not happen */
        while(isgap(mseq->orig_seq[iS][iRo])) {
            iRo++;
        } /* was gap in aligned seq */


        if ( ('\0' != mseq->seq[iS][iR]) ||
             ('\0' != mseq->orig_seq[iS][iRo]) ){

            Log(&rLog, LOG_FATAL, "inconsistency in un/aligned sequence %d\n>%s\n>%s\n",
                  iS, mseq->seq[iS], mseq->orig_seq[iS]);
        }

    } /* 0 <= iS < mseq->nseqs */

} /*** end: TranslateUnknown2Ambiguity() ***/


/**
 * @brief re-attach leading and trailing gaps to alignment
 *
 * @param[in,out] prMSeq
 * alignment structure (at this stage there should be no un-aligned sequences)
 * @param[in] iProfProfSeparator
 * gives sizes of input profiles, -1 if no input-profiles but un-aligned sequences
 *
 * @note leading and tailing profile columns 
 * that only contain gaps have no effect on the alignment 
 * and are removed during the alignment. if they are 
 * encountered a warning message is printed to screen.
 * some users like to preserve these gap columns
 * FS, r213->214
 */
void
ReAttachLeadingGaps(mseq_t *prMSeq, int  iProfProfSeparator)
{

    int i, j;
    int iSize1 = 0; /* #seqs in 1st profile */
    int iSize2 = 0; /* #seqs in 2nd profile */
    int iPPS   = iProfProfSeparator;
    int iLeadO1  = 0; /* leading  gaps in 1st seq of 1st profile */
    int iLeadO2  = 0; /* leading  gaps in 1st seq of 2nd profile */
    int iLeadA1  = 0; /* leading  gaps in 1st seq of final alignment */
    int iLeadA2  = 0; /* leading  gaps in PPS seq of final alignment */
    int iTrailO1 = 0; /* trailing gaps in 1st seq of 1st profile */
    int iTrailO2 = 0; /* trailing gaps in 1st seq of 2nd profile */
    int iTrailA1 = 0; /* trailing gaps in 1st seq of final alignment */
    int iTrailA2 = 0; /* trailing gaps in PPS seq of final alignment */
    int iLen  = 0; /* length of final alignment */
    int iLen1 = 0; /* length of 1st profile */
    int iLen2 = 0; /* length of 2nd profile */
    int iCutHead = 0; /* make up truncation at head */
    int iCutTail = 0; /* make up truncation at tail */
    char *pcIter = NULL;

    if (-1 == iProfProfSeparator){
        return;
    }
    else {
        assert(prMSeq->aligned);
        assert(prMSeq->nseqs > iProfProfSeparator);
    }
    iSize1 = iProfProfSeparator;
    iSize2 = prMSeq->nseqs - iProfProfSeparator;
    iLen  = strlen(prMSeq->seq[0]);
    iLen1 = strlen(prMSeq->orig_seq[0]);
    iLen2 = strlen(prMSeq->orig_seq[iPPS]);

    /* count leading/trailing gaps in 1st sequence of 1st/2nd profile and final alignmant */
    for (iLeadO1  = 0, pcIter =  prMSeq->orig_seq[0];             isgap(*pcIter); pcIter++, iLeadO1++);
    for (iLeadO2  = 0, pcIter =  prMSeq->orig_seq[iPPS];          isgap(*pcIter); pcIter++, iLeadO2++);
    for (iLeadA1  = 0, pcIter =  prMSeq->seq[0];                  isgap(*pcIter); pcIter++, iLeadA1++);
    for (iLeadA2  = 0, pcIter =  prMSeq->seq[iPPS];               isgap(*pcIter); pcIter++, iLeadA2++);
    for (iTrailO1 = 0, pcIter = &prMSeq->orig_seq[0][iLen1-1];    isgap(*pcIter); pcIter--, iTrailO1++);
    for (iTrailO2 = 0, pcIter = &prMSeq->orig_seq[iPPS][iLen2-1]; isgap(*pcIter); pcIter--, iTrailO2++);
    for (iTrailA1 = 0, pcIter = &prMSeq->seq[0][iLen-1];          isgap(*pcIter); pcIter--, iTrailA1++);
    for (iTrailA2 = 0, pcIter = &prMSeq->seq[iPPS][iLen-1];       isgap(*pcIter); pcIter--, iTrailA2++);


    /* turn leading/trailing gaps into truncation */
    iLeadO1  = iLeadO1  > iLeadA1  ? iLeadO1-iLeadA1   : 0;
    iLeadO2  = iLeadO2  > iLeadA2  ? iLeadO2-iLeadA2   : 0;
    iTrailO1 = iTrailO1 > iTrailA1 ? iTrailO1-iTrailA1 : 0;
    iTrailO2 = iTrailO2 > iTrailA2 ? iTrailO2-iTrailA2 : 0;
    iCutHead = iLeadO1  > iLeadO2  ? iLeadO1  : iLeadO2;
    iCutTail = iTrailO1 > iTrailO2 ? iTrailO1 : iTrailO2;

    /* re-allocate and shift memory */
    if ( (iCutHead > 0) || (iCutTail > 0) ){ /* skip if no re-attachment, FS, r244 -> r245 */
        for (i = 0; i < prMSeq->nseqs; i++){
            
            CKREALLOC(prMSeq->seq[i], iLen+iCutHead+iCutTail+2);
            if (iCutHead > 0){ /* skip if no re-attachment, FS, r244 -> r245 */
                memmove(prMSeq->seq[i]+iCutHead, prMSeq->seq[i], iLen);
            }
            for (j = 0; j < iCutHead; j++){
                prMSeq->seq[i][j] = '-';
            }
            for (j = iLen+iCutHead; j < iLen+iCutHead+iCutTail; j++){
                prMSeq->seq[i][j] = '-';
            }
            prMSeq->seq[i][j] = '\0';
        }
    } /* (iCutHead > 0) || (iCutTail > 0) */

} /***  end: ReAttachLeadingGaps()  ***/

/**
 * @brief reallocate enough memory for alignment and
 * attach sequence pointers to profiles
 *
 * @param[in,out] mseq
 * sequence/profile data, increase memory for sequences in profiles
 * @param[out] ppcProfile1
 * pointers to sequencese in 1st profile
 * @param[out] ppcProfile2
 * pointers to sequencese in 2nd profile
 * @param[out] pdWeightsL
 * weights (normalised to 1.0) for sequences in left  profile
 * @param[out] pdWeightsR
 * weights (normalised to 1.0) for sequences in right profile
 * @param[in] pdSeqWeights
 * weights for _all_ sequences in alignment
 * @param[in] iLeafCountL
 * number of sequences in 1st profile
 * @param[in] piLeafListL
 * array of integer IDs of sequences in 1st profile
 * @param[in] iLeafCountR
 * number of sequences in 2nd profile
 * @param[in] piLeafListR
 * array of integer IDs of sequences in 2nd profile
 *
 */
void
PrepareAlignment(mseq_t *mseq, char **ppcProfile1, char **ppcProfile2,
                 double *pdWeightsL, double *pdWeightsR, double *pdSeqWeights,
                 int iLeafCountL, int *piLeafListL,
                 int iLeafCountR, int *piLeafListR)
{

    int iLenL = 0; /* length of 1st profile */
    int iLenR = 0; /* length of 2nd profile */
    int iMaxLen = 0; /* maximum possible length of alignment */
    int i; /* aux */
    double dWeight = 0.00;
    double dWeightInv = 0.00;

    assert(NULL!=mseq);
    assert(NULL!=ppcProfile1);
    assert(NULL!=ppcProfile2);
    assert(NULL!=piLeafListL);
    assert(NULL!=piLeafListR);

    /* get length of profiles,
     * in a profile all sequences should have same length so only look at 1st
     */
    iLenL = strlen(mseq->seq[piLeafListL[0]]);
    iLenR = strlen(mseq->seq[piLeafListR[0]]);
    iMaxLen = iLenL + iLenR + 1;

    /* reallocate enough memory for sequences in alignment (for adding
     * gaps)
     */
    for (i = 0; i < iLeafCountL; i++){
        mseq->seq[piLeafListL[i]] =
            CKREALLOC(mseq->seq[piLeafListL[i]], iMaxLen);
    }
    for (i = 0; i < iLeafCountR; i++){
        mseq->seq[piLeafListR[i]] =
            CKREALLOC(mseq->seq[piLeafListR[i]], iMaxLen);
    }

    /* attach sequences to profiles
     */
    for (i = 0; i < iLeafCountL; i++){
        ppcProfile1[i] = mseq->seq[piLeafListL[i]];
    }
    ppcProfile1[i] = NULL;

    for (i = 0; i < iLeafCountR; i++){
        ppcProfile2[i] = mseq->seq[piLeafListR[i]];
    }
    ppcProfile2[i] = NULL;


    /* remove terminal 'X' for single sequences:
     * it is a quirk of hhalign() to delete 2 individual sequences
     * if 2 terminal X's meet during alignment,
     * just replace (one of) them.
     * this can be undone at the end.
     * profiles -consisting of more than 1 sequence-
     * appear to be all-right.
     * there seems to be no problem with B's and Z's
     */
    if ( (1 == iLeafCountL) && (1 == iLeafCountR) ){

        if ( ('X' == ppcProfile1[0][0]) && ('X' == ppcProfile2[0][0]) ){
#define NOTX 'N'
            ppcProfile1[0][0] = NOTX; /* FIXME: arbitrary assignment */
            ppcProfile2[0][0] = NOTX; /* FIXME: arbitrary assignment */
        }
        if ( ('X' == ppcProfile1[0][iLenL-1]) && ('X' == ppcProfile2[0][iLenR-1]) ){
            ppcProfile1[0][iLenL-1] = NOTX; /* FIXME: arbitrary assignment */
            ppcProfile2[0][iLenR-1] = NOTX; /* FIXME: arbitrary assignment */
        }
    }

    /* obtain sequence weights
     */
    if (NULL != pdSeqWeights){

        dWeight = 0.00;
        for (i = 0; i < iLeafCountL; i++){
            register double dAux = pdSeqWeights[piLeafListL[i]];
#ifndef NDEBUG
            if (dAux <= 0.00){
                Log(&rLog, LOG_DEBUG, "seq-weight %d = %f", piLeafListL[i], dAux);
            }
#endif
            pdWeightsL[i] = dAux;
            dWeight      += dAux;
        } /* 0 <= i < iLeafCountL */
        dWeightInv = 1.00 / dWeight;
        for (i = 0; i < iLeafCountL; i++){
            pdWeightsL[i] *= dWeightInv;
        }

        dWeight = 0.00;
        for (i = 0; i < iLeafCountR; i++){
            register double dAux = pdSeqWeights[piLeafListR[i]];
#ifndef NDEBUG
            if (dAux <= 0.00){
                Log(&rLog, LOG_DEBUG, "seq-weight %d = %f", piLeafListR[i], dAux);
            }
#endif
            pdWeightsR[i] = dAux;
            dWeight      += dAux;
        } /* 0 <= i < iLeafCountL */
        dWeightInv = 1.00 / dWeight;
        for (i = 0; i < iLeafCountR; i++){
            pdWeightsR[i] *= dWeightInv;
        }
    } /* (NULL != pdSeqWeights) */
    else {
        pdWeightsL[0] = pdWeightsR[0] = -1.00;
    }

#if TRACE
    for (i = 0; i < iLeafCountL; i++){
        Log(&rLog, LOG_FORCED_DEBUG, "ppcProfile1[%d/%d] pointing to mseq %d = %s",
                  i, iLeafCountR, piLeafListL[i], ppcProfile1[i]);
    }
    for (i = 0; i < iLeafCountR; i++){
        Log(&rLog, LOG_FORCED_DEBUG, "ppcProfile2[%d/%d] pointing to mseq %d = %s",
                  i, iLeafCountR, piLeafListR[i], ppcProfile2[i]);
    }
#endif
    
    return;

} /*** end: PrepareAlignment() ***/


/**
 * @brief PosteriorProbabilities() calculates posterior probabilities 
 *  of aligning a single sequences on-to an alignment containing this sequence
 *
 * @param[in] prMSeq
 * holds the aligned sequences [in] 
 * @param[in] rHMMalignment
 * HMM of the alignment in prMSeq
 * @param[in] rHhalignPara
 * various parameters read from commandline, needed by hhalign()
 * @param[in] pcPosteriorfile
 * name of file into which posterior probability information is written
 *
 * @return score of the alignment FIXME what is this?
 *
 * @note the PP-loop can be parallelised easily FIXME
 */

int
PosteriorProbabilities(mseq_t *prMSeq, hmm_light rHMMalignment, hhalign_para rHhalignPara, char *pcPosteriorfile)
{

    double dScore = 0.0;
    int i;
    int iS; /* iterator for sequences */
    int iNseq = prMSeq->nseqs; /* number of sequences */
    int iLenHMM = rHMMalignment.L; /* length of alignment */
    int iViterbiCount = 0; /* counts how often Viterbi is triggered */
    char **ppcCopy = NULL;
    char **ppcRepresent = NULL;
    char zcAux[10000] = {0};
    char zcError[10000] = {0};
    hhalign_scores *prHHscores = NULL;
    FILE *pfPP = NULL;

    pfPP = fopen(pcPosteriorfile, "w");
    fprintf(pfPP, "#1.i\t2.name\t\t3.L1\t4.L2\t5.sum\t\t6.sum/L1\t7.HH\n");


    prHHscores = CKCALLOC(iNseq, sizeof(hhalign_scores));
    for (iS = 0; iS < iNseq; iS++){
        memset(&(prHHscores[iS]), 0, sizeof(hhalign_scores));
    }

    ppcCopy         = CKCALLOC(1, sizeof(char *));
    ppcCopy[0]      = CKCALLOC(2*iLenHMM, sizeof(char));
    ppcRepresent    = CKCALLOC(1, sizeof(char *));
    ppcRepresent[0] = CKCALLOC(2*iLenHMM, sizeof(char));

    for (i = 0; i < iLenHMM; i++){
        ppcRepresent[0][i] = rHMMalignment.seq[rHMMalignment.ncons][i+1];
    }

    /* FIXME: this loop can be parallelised, FS r288 */
    iViterbiCount = 0;
    for (iS = 0; iS < iNseq; iS++){

        /* single sequences may very well contain leading/trailing gaps,
         * this will trigger warning in Transfer:hhalignment-C.h */
        char *pcIter = NULL;
        for (pcIter = prMSeq->orig_seq[iS]; ('-' == *pcIter) || ('.' == *pcIter); pcIter++);
        strcpy(ppcCopy[0], pcIter/*prMSeq->orig_seq[iS]*/);
        for (pcIter = &ppcCopy[0][strlen(ppcCopy[0])-1]; ('-' == *pcIter) || ('.' == *pcIter); pcIter--);
        pcIter++;
        *pcIter = '\0';
        
        zcError[0] = '\0';
        hhalign(ppcCopy, 1, NULL, 
                ppcRepresent, 0, NULL, 
                &dScore, &rHMMalignment, &rHMMalignment, 
                NULL, NULL, NULL, NULL,
                rHhalignPara, &prHHscores[iS], 0/* debug */, FALSE /* Log-Level*/, 
                zcAux, zcError);
        if (NULL != strstr(zcError, "Viterbi")){
            iViterbiCount++;
        }

    } /* 0 <= iS < iNseq */
    Log(&rLog, LOG_INFO, "Viterbi algorithm triggered %d times (out of %d)", iViterbiCount, iNseq-1);

    for (iS = 0; iS < iNseq; iS++){

        fprintf(pfPP, "%d\t%10s\t%3d"
               , iS, prMSeq->sqinfo[iS].name, (int)(strlen(prMSeq->orig_seq[iS])));
        fprintf(pfPP, "\t%3d\t%f\t%f\t%f", 
                prHHscores[iS].L, prHHscores[iS].sumPP, prHHscores[iS].sumPP/strlen(prMSeq->orig_seq[iS]), prHHscores[iS].hhScore);
        fprintf(pfPP, "\n");
    } /* 0 <= iS < iNseq */



    fclose(pfPP); pfPP = NULL;
    free(ppcRepresent[0]); ppcRepresent[0] = NULL;
    free(ppcCopy[0]); ppcCopy[0] = NULL;
    free(ppcRepresent); ppcRepresent = NULL;
    free(ppcCopy); ppcCopy = NULL;
    free(prHHscores); prHHscores = NULL;

    return 1;

} /* this is the end of PosteriorProbabilities() */


/**
 * @brief sequentially align (chain) sequences
 *
 * @param[in,out] prMSeq
 * holds the un-aligned sequences (in) and the final alignment (out)
 * @param[in] rHhalignPara
 * various parameters needed by hhalign()
 * @param[in] iClustersize
 * parameter that controls how often HMM is updated
 *
 * @note chained alignment takes much longer than balanced alignment
 * because at every step ClustalO has to scan all previously aligned residues.
 * for a balanced tree this takes N*log(N) time but for a chained tree 
 * it takes N^2 time.
 * This function has a short-cut, that the HMM need not be updated 
 * for every single alignment step, but the HMM from the previous 
 * step(s) can be re-cycled. The HMM is updated (i) at the very first step, 
 * (ii) if a gap has been inserted into the HMM during alignment or 
 * (iii) if the HMM has been used for too many steps without having 
 * been updated. This update-frequency is controlled by the input 
 * parameter iClustersize. iClustersize is the number of sequences used 
 * to build a HMM to allow for one non-updating step. For example, 
 * if iClustersize=100 and a HMM has been build from 100 sequences, 
 * then this HMM can be used once without updating. If the HMM has 
 * been built from 700 sequences (and iClustersize=100) then this 
 * HMM can be used 7-times without having to be updated, etc. 
 * For this reason the initial iClustersize sequences are always 
 * aligned with fully updated HMMs. 
 */
double 
PileUp(mseq_t *prMSeq, hhalign_para rHhalignPara, int iClustersize)
{
    /* @<variables local to PileUp> */
    int iI; /* iterator for sequences */
    int iCountPro = 0; /* count how many sequences are already in profile */
    int *piLeafListSeq = NULL; /* lists required by PrepareAlignment() for attaching sequences pointers to ppcProfileSeq */
    int *piLeafListPro = NULL; /* -"-  ppcProfilePro */
    double *pdWeightL = NULL; /* argument used by hhalign, not used here */
    double *pdWeightR = NULL; /* -"- */
    double *pdWeightP = NULL; /* -"- */
    char **ppcProfileSeq = NULL; /* pointer to sequences in profile */
    char **ppcProfilePro = NULL; /* pointer to sequences in profile */
    double dScore = -1.0; /* alignment score, calculated by hhalign() */
    hhalign_scores rHHscores = {0}; /* more explicit scores than dScore */
    hmm_light *prHMM = NULL; /* HMM against which to align single sequence */
    hmm_light rHMMprofile = {0}; /* HMM calculated externally from profile */
    int iHHret = -1; /* return value of hhalign() */
    char zcAux[10000] = {0}; /* auxilliary print string used by hhalign() */
    char zcError[10000] = {0}; /* error message printed by hhalign() */
    char *pcConsensPro = NULL, /* auxilliary sequence strings required by hhalign() */
        *pcReprsntPro = NULL, /* not used here */
        *pcConsensSeq = NULL, *pcReprsntSeq = NULL;
    bool bHMM_stale = YES; /* flag if HMM has to be updated */
    int iSinceLastUpdate = -1; /* counts how often HMM has been used since last update */
    progress_t *prProgress; /* structure required for progress report */
    bool bPrintCR = (rLog.iLogLevelEnabled<=LOG_VERBOSE) ? FALSE : TRUE; /* progress report */
    char **ppcReprsnt = NULL; /* string used to represent HMM */
    ppcReprsnt    = CKCALLOC(1, sizeof(char *));

    /* weights are not really used, but hhalign() expects them */
    pdWeightL = (double *)CKMALLOC(prMSeq->nseqs * sizeof(double));
    pdWeightR = (double *)CKMALLOC(prMSeq->nseqs * sizeof(double));
    pdWeightP = (double *)CKMALLOC(prMSeq->nseqs * sizeof(double));

    /* lists used for attaching ppcProfileSeq/ppcProfilePro to prMSeq */
    piLeafListSeq = (int *)CKMALLOC(1 * sizeof(int));
    piLeafListPro = (int *)CKMALLOC(prMSeq->nseqs * sizeof(int));

    ppcProfileSeq    = (char **)CKMALLOC(1 * sizeof(char *));
    ppcProfileSeq[0] = NULL;
    ppcProfilePro    = (char **)CKMALLOC(prMSeq->nseqs * sizeof(char *));
    for (iI = 0; iI < prMSeq->nseqs; iI++){
        ppcProfilePro[iI] = NULL;
    }
    
    piLeafListPro[0] = 0; /* first sequences in profile is simply un-aligned sequence */
    iCountPro = 1; /* 'profile' now contains 1 sequence */

    NewProgress(&prProgress, LogGetFP(&rLog, LOG_INFO),
                "Progressive alignment progress", bPrintCR);

    /* build up the initial alignment:
     * sequences are chained but hhalign() is used normally, 
     * that is, one sequence is aligned to a profile */
    for (iI = 1; iI < MIN(prMSeq->nseqs, iClustersize); iI++){

        piLeafListSeq[0] = iI;

        /* PrepareAlignment() connects ppcProfileSeq/ppcProfilePro and prMSeq */
        PrepareAlignment(prMSeq, ppcProfilePro, ppcProfileSeq,
                         pdWeightL, pdWeightR, pdWeightP, 
                         iI, piLeafListPro,
                         1, piLeafListSeq);

        /* align 1 (5th argument) sequence to the growing profile 
         * (2nd argument number of sequences is iI), 
         * external HMM not used at this stage (prHMM=NULL) */
        iHHret = hhalign(ppcProfilePro, iI, NULL/*pdWeightL*/,
                         ppcProfileSeq,  1, NULL/*pdWeightR*/,
                         &dScore, prHMM, prHMM,
                         pcConsensPro, pcReprsntPro,
                         pcConsensSeq, pcReprsntSeq,
                         rHhalignPara, &rHHscores,
                         CALL_FROM_REGULAR/* DEBUG ARGUMENT */, rLog.iLogLevelEnabled,
                         zcAux, zcError);

        piLeafListPro[iI] = iI;
        iCountPro = iI+1;
        ProgressLog(prProgress, iI, prMSeq->nseqs-1, FALSE); /* FIXME: test! FS, r288 */

    } /* 1 <= iI < MIN(prMSeq->nseqs, iClustersize) */

    /* now align single sequences not to a profile but to a HMM of a potentially old profile */
    while (iI < prMSeq->nseqs){

        /* if HMM not up-to-date then create new HMM for profile:
         * HMM is stale (i) at the beginning, 
         * (ii) if a gap has been inserted into the HMM during previous step,
         * (iii) too many steps after the last up-date */
        if ( (YES == bHMM_stale) ){
            FreeHMMstruct(&rHMMprofile);
            if (OK !=
                AlnToHMM2(&rHMMprofile, rHhalignPara, prMSeq->seq, iI)
                ) {
                Log(&rLog, LOG_ERROR, "Couldn't convert alignment to HMM. Will not proceed with chained alignment, step %d", iI);
            }
            bHMM_stale = NO;
            iSinceLastUpdate = 0;
        }

        /* align single sequence to HMM */
        piLeafListSeq[0] = iI;

        PrepareAlignment(prMSeq, ppcProfilePro, ppcProfileSeq,
                         pdWeightL, pdWeightR, pdWeightP, 
                         iI, piLeafListPro,
                         1, piLeafListSeq);

        /* use representative sequence to track where gaps are inserted into HMM */
        ppcReprsnt[0] = CKREALLOC(ppcReprsnt[0], strlen(ppcProfilePro[0])+strlen(ppcProfileSeq[0])+1);
        memcpy(ppcReprsnt[0], rHMMprofile.seq[rHMMprofile.ncons]+1, rHMMprofile.L);
        ppcReprsnt[0][rHMMprofile.L] = '\0';

        /* align 1 (5th argument) sequence to a HMM (2nd argument number of sequences is '0', 
         * not iI as for a profile or '1' for one representative sequence), 
         * external HMM (rHMMprofile calculated by AlnToHMM2()) is used at this stage */
        prHMM = &rHMMprofile;
        hhalign(ppcReprsnt, 0, NULL,
                ppcProfileSeq, 1, NULL, 
                &dScore, prHMM, prHMM,
                pcConsensPro, pcReprsntPro,
                pcConsensSeq, pcReprsntSeq,
                rHhalignPara, &rHHscores,
                CALL_FROM_ALN_HMM/* DEBUG ARGUMENT */, rLog.iLogLevelEnabled,
                zcAux, zcError);
        iSinceLastUpdate++; /* HMM has been used one more time */
        ProgressLog(prProgress, iI, prMSeq->nseqs-1, FALSE); /* FIXME: test! FS, r288 */

        /* check if gaps have been inserted into HMM.
         * if this is the case then HMM is stale 
         * and the profile used to create the HMM must get gaps at the same positions.
         * gaps are shown in the representative sequence */
        if (rHMMprofile.L != strlen(ppcReprsnt[0])){

            int iSeq; /* iterate through sequences of profile */

            /* manually insert gaps into existing profile, 
             * there should already be enough memory */
#ifdef HAVE_OPENMP
            /* all sequences in profile obtain gaps at the same position, 
             * can do this in parallel */
	        #pragma omp parallel for 
#endif
            for (iSeq = 0; iSeq < iI; iSeq++){

                int iPtrRep, /* points to 'residue' in representative sequence */
                    iPtrPrf; /* points to the corresponding position in the profile */

                for (iPtrRep = strlen(ppcReprsnt[0])-1, iPtrPrf = strlen(ppcProfilePro[iSeq])-1; 
                     iPtrRep >= 0; iPtrRep--){

                    if ('-' == ppcReprsnt[0][iPtrRep]){ /* gaps are newly introduced into profile */
                        ppcProfilePro[iSeq][iPtrRep] = '-';
                    }
                    else { /* non-gap characters (residues) are shifted towards the back */
                        ppcProfilePro[iSeq][iPtrRep] = ppcProfilePro[iSeq][iPtrPrf];
                        iPtrPrf--;
                    }
                    if (iPtrRep == iPtrPrf){ 
                        /* if profile and representative seq point to same position then no more gaps */
                        break; 
                    }
                } /* strlen(ppcReprsnt[0])-1 >= iPtrRep >= 0 */
                ppcProfilePro[iSeq][strlen(ppcReprsnt[0])] = '\0';

            } /* 0 <= iSeq < iI */
            
            /* make HMM stale */
            bHMM_stale = YES;

        } /* strlen(represent) != HMM.L (gaps inserted into HMM) */

        /* check if maximum number of sequences exceeded */
        if ( (iClustersize <= 1) || (iSinceLastUpdate > (double)(rHMMprofile.N_in)/(double)(iClustersize)) ){
            bHMM_stale = YES;
        }


        piLeafListPro[iI] = iI;
        iCountPro = iI+1;

        iI++;

    } /* iI < prMSeq->nseqs */
    ProgressDone(prProgress);
    FreeProgress(&prProgress);

    FreeHMMstruct(&rHMMprofile);
    CKFREE(pdWeightL); CKFREE(pdWeightR); CKFREE(pdWeightP); 
    if (NULL != ppcReprsnt){ 
        if (NULL != ppcReprsnt[0]){ 
            CKFREE(ppcReprsnt[0]); 
        }
        CKFREE(ppcReprsnt); 
    }
    CKFREE(piLeafListSeq); CKFREE(piLeafListPro);

    return 0.0;

} /* this is the end of PileUp() */





/**
 * @brief wrapper for hhalign. This is a frontend function to
 * the ported hhalign code.
 *
 * @param[in,out] prMSeq
 * holds the unaligned sequences [in] and the final alignment [out]
 * @param[in] piOrderLR
 * holds order in which sequences/profiles are to be aligned,
 * even elements specify left nodes, odd elements right nodes,
 * if even and odd are same then it is a leaf
 * @param[in] pdSeqWeights
 * Weight per sequence. No weights used if NULL
 * @param[in] iNodeCount
 * number of nodes in tree, piOrderLR has 2*iNodeCount elements
 * @param[in] prHMMList
 * List of background HMMs (transition/emission probabilities)
 * @param[in] iHMMCount
 * Number of input background HMMs
 * @param[in] iProfProfSeparator
 * Gives the number of sequences in the first profile, if in
 * profile/profile alignment mode (iNodeCount==3). That assumes mseqs
 * holds the sequences of profile 1 and profile 2.
 * @param[in] rHhalignPara
 * various parameters read from commandline
 *
 * @return score of the alignment FIXME what is this?
 *
 * @note complex function. could use some simplification, more and
 * documentation and a struct'uring of piOrderLR
 * 
 * @note HHalignWrapper can be entered in 2 different ways:
 * (i) all sequences are un-aligned (ii) there are 2 (aligned) profiles.
 * in the un-aligned case (i) the sequences come straight from Squid, 
 * that is, they have been sanitised, all non-alphabetic residues 
 * have been rendered as X's. In profile mode (ii) one profile may 
 * have been produced internally. In that case residues may have 
 * been translated back into their 'native' form, that is, they 
 * may contain un-sanitised residues. These will cause trouble 
 * during alignment
 *
 * @note: introduced argument hhalign_para rHhalignPara; FS, r240 -> r241
 * @note: if hhalign() fails then try with Viterbi by setting MAC-RAM=0; FS, r241 -> r243
 */


double
HHalignWrapper(mseq_t *prMSeq, int *piOrderLR,
               double *pdSeqWeights, int iNodeCount,
               hmm_light *prHMMList, int iHMMCount,
               int iProfProfSeparator, hhalign_para rHhalignPara)
{
    int iN; /* node iterator */
    int *piLeafCount = NULL; /* number of leaves beneath a certain node */
    int **ppiLeafList = NULL; /* list of leaves beneath a certain node */
    char **ppcProfile1 = NULL; /* pointer to sequences in profile */
    char **ppcProfile2 = NULL; /* pointer to sequences in profile */
    char *pcReprsnt1 = NULL; /* representative of HMM aligned to left */
    char *pcReprsnt2 = NULL; /* representative of HMM aligned to right */
    char **ppcReprsnt1 = &pcReprsnt1; /* representative of HMM aligned to L */
    char **ppcReprsnt2 = &pcReprsnt2; /* representative of HMM aligned to R */
    char *pcConsens1 = NULL; /* copy of  left sequence */
    char *pcConsens2 = NULL; /* copy of right sequence */
    char **ppcCopy1 = /*&pcCopy1*/NULL; /* copy of  left sequences */
    char **ppcCopy2 = /*&pcCopy2*/NULL; /* copy of right sequences */
    double *pdScores = NULL; /* alignment scores (seq/HMM) */
    double dScore = 0.0; /* alignment score (seq/HMM) */
    int iAux_FS = 0;
    char zcAux[10000] = {0};
    char zcError[10000] = {0};
    int i; /* aux */
    progress_t *prProgress;
    int iAlnLen; /* alignment length */
    double *pdWeightsL = NULL; /* sequence weights of left  profile */
    double *pdWeightsR = NULL; /* sequence weights of right profile */
    int iMergeNodeCounter = 0;
    hmm_light *prHMM = NULL; 
    hmm_light *prHMMleft = NULL;
    hmm_light *prHMMrght = NULL;
    hmm_light *prHMMnull = NULL;
    bool bPrintCR = (rLog.iLogLevelEnabled<=LOG_VERBOSE) ? FALSE : TRUE;
#if TIMING
    char zcStopwatchMsg[1024];
    Stopwatch_t *stopwatch = StopwatchCreate();
    StopwatchZero(stopwatch);
    StopwatchStart(stopwatch);
#endif
    hhalign_scores rHHscores = {0};

#if 0 /* DEVEL 291 */
    if (NULL != prHMMList) { /* FIXME DEVEL 291: replace this outer test with iHMMCount>0*/
        if (iHMMCount>1) {
            Log(&rLog, LOG_WARN, "FIXME: Using only first of %u HMMs (needs implementation)", iHMMCount);
        }
        prHMM = &(prHMMList[0]); /* FIXME DEVEL 291: check for NULL */
    } else {
        /* FIXME: prHMM not allowed to be NULL and therefore pseudo allocated here */
        prHMM = (hmm_light *) CKCALLOC(1, sizeof(hmm_light));
    }
#else
    prHMMnull = (hmm_light *) CKCALLOC(1, sizeof(hmm_light));
    if ( (iHMMCount > 0) && (NULL != prHMMList) ){
        prHMM = &(prHMMList[0]);
        if (iHMMCount > 1) {
            Log(&rLog, LOG_WARN, "FIXME: Using only first of %u HMMs (needs implementation)", iHMMCount);
        }
    }
    else {
        prHMM = prHMMnull; /* prHMM not allowed to be NULL and therefore assigned to pseudo allocated  */
    }
#endif
    
    assert(NULL!=prMSeq);
    
    if (NULL==piOrderLR) {
        assert(3==iNodeCount);
    }
    SanitiseUnknown(prMSeq);

    /* hhalign was not made for DNA/RNA. So warn if sequences are not
     * protein
     */
    if (SEQTYPE_PROTEIN != prMSeq->seqtype) {
        /*Log(&rLog, LOG_WARN, "%s alignment is still experimental.",
          SeqTypeToStr(prMSeq->seqtype));*/
		if(prMSeq->seqtype == SEQTYPE_DNA)
			rHhalignPara.bIsDna = true;
		if(prMSeq->seqtype == SEQTYPE_RNA)
			rHhalignPara.bIsRna = true;
    }

    /* hhalign produces funny results if sequences contain gaps, so
     * dealign. Only way to use alignment info is to use it as a
     * background HMM
     */
    if (TRUE == prMSeq->aligned) {
        Log(&rLog, LOG_DEBUG, "Dealigning aligned sequences (inside %s)", __FUNCTION__);
        DealignMSeq(prMSeq);
    }

#if TRACE
    Log(&rLog, LOG_FORCED_DEBUG, "iNodeCount = %d", iNodeCount);
#endif

    
    /* allocate top-level memory for leaf tracking arrays and profiles,
     * and sequence weights*/
    piLeafCount = CKCALLOC(iNodeCount, sizeof(int));
    ppiLeafList = CKCALLOC(iNodeCount, sizeof(int *));
    ppcProfile1 = CKCALLOC(prMSeq->nseqs*2-1, sizeof(char *));
    ppcProfile2 = CKCALLOC(prMSeq->nseqs*2-1, sizeof(char *));
    pdScores    = CKCALLOC(iNodeCount, sizeof(double));
    pdWeightsL  = CKCALLOC(iNodeCount, sizeof(double));
    pdWeightsR  = CKCALLOC(iNodeCount, sizeof(double));

    NewProgress(&prProgress, LogGetFP(&rLog, LOG_INFO),
                "Progressive alignment progress", bPrintCR);


    
    /* Profile-profile alignment? Then setup piLeafCount and
     * piLeafList here. FIXME this is just an awful haaaack
     */
    if (iNodeCount==3 && NULL==piOrderLR) {
        int iSizeProf1 = iProfProfSeparator;
        int iSizeProf2 = prMSeq->nseqs - iProfProfSeparator;

        piLeafCount[0] = iSizeProf1;
        ppiLeafList[0] = (int *)CKMALLOC(iSizeProf1 * sizeof(int));
        for (i=0;i<iSizeProf1;i++) {
            ppiLeafList[0][i] = i;
        }
        piLeafCount[1] = iSizeProf2;
        ppiLeafList[1] = (int *)CKMALLOC(iSizeProf2 * sizeof(int));
        for (i=0;i<iSizeProf2;i++) {
            ppiLeafList[1][i] = i+iSizeProf1;
        }
        
        /* awful hack inside an awful hack: we had to setup piLeafCount
         * and piLeafList outside the node iteration. this which is
         * normally done at leaf level inside the node iteration. to
         * avoid overwriting the already setup vars set...
         */
        iNodeCount=1;

        piOrderLR =  (int *)CKCALLOC(DIFF_NODE, sizeof(int));
        piOrderLR[LEFT_NODE] = 0;
        piOrderLR[RGHT_NODE] = 1;
        piOrderLR[PRNT_NODE] = 2;
    }
    


    iMergeNodeCounter = 0;
    for (iN = 0; iN < iNodeCount; iN++){

        register int riAux = DIFF_NODE * iN;

        /*LOG_DEBUG("node %d ", iN);*/

        if (piOrderLR[riAux+LEFT_NODE] == piOrderLR[riAux+RGHT_NODE]){

            register int riLeaf = piOrderLR[riAux+LEFT_NODE];
#if TRACE
            if (NULL == pdSeqWeights) {
                Log(&rLog, LOG_FORCED_DEBUG, "node %d is a leaf with entry %d (seq %s)",
                          iN, riLeaf, prMSeq->sqinfo[riLeaf].name);
            } else {
                Log(&rLog, LOG_FORCED_DEBUG, "node %d is a leaf with entry %d  (seq %s) and weight %f",
                          iN, riLeaf, prMSeq->sqinfo[riLeaf].name, pdSeqWeights[riLeaf]);
            }
#endif
            /* left/right entry same, this is a leaf */
            piLeafCount[piOrderLR[riAux+PRNT_NODE]] = 1; /* number of leaves is '1' */
            ppiLeafList[piOrderLR[riAux+PRNT_NODE]] = (int *)CKMALLOC(1 * sizeof(int));
            ppiLeafList[piOrderLR[riAux+PRNT_NODE]][0] = riLeaf;

        } /* was a leaf */
        else {
            int iL, iR, iP; /* ID of left/right nodes, parent */
            int i, j; /* aux */

            Log(&rLog, LOG_DEBUG,
                "merge profiles at node %d", iN, piOrderLR[riAux]);

            /* iNodeCount - prMSeq->nseqs = total # of merge-nodes 
             * unless in profile/profile alignment mode
             */
            if (1 == iNodeCount) {
                ProgressLog(prProgress, ++iMergeNodeCounter,
                            1, FALSE);
            } else {
                ProgressLog(prProgress, ++iMergeNodeCounter,
                            iNodeCount - prMSeq->nseqs, FALSE);                
            }

            /* left/right entry are not same, this is a merge node */
            iL = piOrderLR[riAux+LEFT_NODE];
            iR = piOrderLR[riAux+RGHT_NODE];
            iP = piOrderLR[riAux+PRNT_NODE];
            piLeafCount[iP] = piLeafCount[iL] + piLeafCount[iR];
            ppiLeafList[iP] = (int *)CKMALLOC(piLeafCount[iP] * sizeof(int));

            for (i = j = 0; i < piLeafCount[iL]; i++, j++){
                ppiLeafList[iP][j] = ppiLeafList[iL][i];
            }
            for (i = 0; i < piLeafCount[iR]; i++, j++){
                ppiLeafList[iP][j] = ppiLeafList[iR][i];
            }

            /* prepare simulation arena:
             * - make sure enough memory in sequences
             * - attach sequence pointers to profiles
             */
            /* idea: switch template and query according to nseq? */

            PrepareAlignment(prMSeq, ppcProfile1, ppcProfile2,
                             pdWeightsL, pdWeightsR, pdSeqWeights,
                             piLeafCount[iL], ppiLeafList[iL],
                             piLeafCount[iR], ppiLeafList[iR]);
            if (rLog.iLogLevelEnabled <= LOG_DEBUG){
                int i;
                FILE *fp = LogGetFP(&rLog, LOG_DEBUG);
                Log(&rLog, LOG_DEBUG, "merging profiles %d & %d", iL, iR);
                for (i = 0; i < piLeafCount[iL]; i++){
                    fprintf(fp, "L/#=%3d (ID=%3d, w=%f): %s\n",
                           i, ppiLeafList[iL][i], pdWeightsL[i], ppcProfile1[i]);
                }
                for (i = 0; i < piLeafCount[iR]; i++){
                    fprintf(fp, "R/#=%3d (ID=%3d, w=%f): %s\n",
                           i, ppiLeafList[iR][i], pdWeightsR[i], ppcProfile2[i]);
                }
            }

#if 1 /* DEVEL 291 */
            /*
             * Note: if there is a HMM-batch file, then prMSeq->ppiHMMBindex is not NULL;
             *       ppiLeafList[iL/iR][0] is the 'lead' sequence in a profile;
             *       prMSeq->ppiHMMBindex[ppiLeafList[iL][0]] are the HMMs associated with the lead sequence;
             *         this could be NULL if there are no HMMs associated with this particular sequence
             *       at the moment only 1st HMM can be used, prMSeq->ppiHMMBindex[ppiLeafList[iL][0]][0];
             *         the index of this HMM can be '-1' if the specified HMM file does not exist
             * Note: we only use prHMMleft/prHMMrght, even if global HMM (--hmm-in) is specified
             **/
            if ( (NULL != prMSeq->ppiHMMBindex) && (NULL != prMSeq->ppiHMMBindex[ppiLeafList[iL][0]]) && 
                 (prMSeq->ppiHMMBindex[ppiLeafList[iL][0]][0] > -1) ){
                prHMMleft = &(prHMMList[prMSeq->ppiHMMBindex[ppiLeafList[iL][0]][0]]);
            }
            else if (iHMMCount > 0){
                prHMMleft = prHMM;
            }
            else {
                prHMMleft = prHMMnull;
            }

            if ( (NULL != prMSeq->ppiHMMBindex) && (NULL != prMSeq->ppiHMMBindex[ppiLeafList[iR][0]]) &&
                 (prMSeq->ppiHMMBindex[ppiLeafList[iR][0]][0] > -1) ){
                prHMMrght = &(prHMMList[prMSeq->ppiHMMBindex[ppiLeafList[iR][0]][0]]);
            }
            else if (iHMMCount > 0){
                prHMMrght = prHMM;
            }
            else {
                prHMMrght = prHMMnull;
            }
#endif
            
            /* align individual sequences to HMM;
             * - use representative sequence to get gapping
             * - create copies of both, individual/representative sequences
             *   as we don't want to introduce gaps into original
             *
             * FIXME: representative sequence is crutch, should use
             * full HMM but that does not seem to work at all
             * -- try harder! Fail better!
             */
            /*if ( (piLeafCount[iL] <= APPLY_BG_HMM_UP_TO_TREE_DEPTH) && (0 != prHMM->L) ){*/
            if (0 != prHMMleft->L){
                int i, j;
                pcReprsnt1 = CKCALLOC(prHMMleft->L+strlen(ppcProfile1[0])+1, sizeof(char));
                for (i = 0; i < prHMMleft->L; i++){
                    pcReprsnt1[i] = prHMMleft->seq[prHMMleft->ncons][i+1];
                }
                ppcCopy1 = CKCALLOC(piLeafCount[iL], sizeof(char *));
                for (j = 0; j < piLeafCount[iL]; j++){
                    ppcCopy1[j] = CKCALLOC(prHMMleft->L+strlen(ppcProfile1[0])+1, sizeof(char));
                    for (i = 0; i < (int) strlen(ppcProfile1[0]); i++){
                        ppcCopy1[j][i] = ppcProfile1[j][i];
                    }
                }

                {
                    /* the size of the elements in the forward/backward matrices 
                       depends very much on the lengths of the profiles _and_ 
                       in which position (1st/2nd) the longer/shorter profile/HMM is.
                       the matrix elements can easily exceed the size of a (long?) double
                       if the longer profile/HMM is associated with the query (q) and the 
                       shorter with the target (t). 
                       FIXME: however, pseudo-count adding may also depend on position, 
                       this is only really tested for the HMM being in the 1st position (q)
                       MUST TEST THIS MORE THOROUGHLY

                       this switch appears to be most easily (although unelegantly) 
                       effected here. Don't want to do it (upstairs) in PrepareAlignment() 
                       because it might jumble up the nodes. Don't want to do it in hhalign() 
                       either because ppcProfile1/2 and q/t may be used independently.
                       FS, r236 -> r237
                    */
                    int iLenA = strlen(ppcCopy1[0]);
                    int iLenH = prHMMleft->L;
                    int iHHret = 0;
                    
                    if (iLenH < iLenA){
                        iHHret = hhalign(ppcReprsnt1, 0/* only one representative seq*/, NULL,
                                         ppcCopy1, piLeafCount[iL], pdWeightsL,
                                         &dScore, prHMMleft, prHMMleft,
                                         NULL, NULL, NULL, NULL,
                                         rHhalignPara, &rHHscores, 
                                         iAux_FS++, /* DEBUG ARGUMENT */ rLog.iLogLevelEnabled,
                                         zcAux, zcError);
                    }
                    else {
                        iHHret = hhalign(ppcCopy1, piLeafCount[iL], pdWeightsL,
                                         ppcReprsnt1, 0/* only one representative seq*/, NULL,
                                         &dScore, prHMMleft, prHMMleft,
                                         NULL, NULL, NULL, NULL,
                                         rHhalignPara, &rHHscores, 
                                         iAux_FS++, /* DEBUG ARGUMENT */ rLog.iLogLevelEnabled,
                                         zcAux, zcError);
                    }
                    if ( (0 != iHHret) && (rLog.iLogLevelEnabled <= LOG_VERBOSE) ){ /* FS, r255 -> */
                        fprintf(stderr, "%s:%s:%d: (not essential) HMM pre-alignment failed, error %d, \n"
                                "\t#=%d (len=%d), lead-seq=%s, len(HMM)=%d\n%s\nCARRY ON REGARDLESS\n", 
                                __FUNCTION__, __FILE__, __LINE__, iHHret, 
                                piLeafCount[iL], (int)strlen(ppcCopy1[0]), prMSeq->sqinfo[ppiLeafList[iL][0]].name, 
                                (int)strlen(ppcReprsnt1[0]), zcError);
                    }
                }
                pdScores[ppiLeafList[iL][0]] = dScore;
#if 0
                printf("score: %f\nL: %s\nH: %s\n",
                       dScore, ppcCopy1[0], ppcReprsnt1[0]);
#endif
                /* assemble 'consensus';
                 * this is not a real consensus, it is more a gap indicator,
                 * for each position it consists of residues/gaps in the 1st sequences,
                 * or a residue (if any) of the other sequences.
                 * it only contains a gap if all sequences of the profile
                 * have a gap at this position
                 */
                pcConsens1 = CKCALLOC(prHMMleft->L+strlen(ppcProfile1[0])+1, sizeof(char));
                for (i = 0; i < prHMMleft->L; i++){
                    for (j = 0, pcConsens1[i] = '-'; (j < piLeafCount[iL]) && ('-' == pcConsens1[i]); j++){
                        pcConsens1[i] = ppcCopy1[j][i];
                    }
                }
#if 0
                for (j = 0; (j < piLeafCount[iL]); j++){
                    printf("L%d:%s\n", j, ppcCopy1[j]);
                }
                printf("LC:%s\n", pcConsens1);
#endif
            } /* ( (1 == piLeafCount[iL]) && (0 != prHMM->L) ) */

            /*if ( (piLeafCount[iR] <= APPLY_BG_HMM_UP_TO_TREE_DEPTH) && (0 != prHMM->L) ){*/
            if (0 != prHMMrght->L){
                int i, j;

                pcReprsnt2 = CKCALLOC(prHMMrght->L+strlen(ppcProfile2[0])+1, sizeof(char));
                for (i = 0; i < prHMMrght->L; i++){
                    pcReprsnt2[i] = prHMMrght->seq[prHMMrght->ncons][i+1];
                }
                ppcCopy2 = CKCALLOC(piLeafCount[iR], sizeof(char *));
                for (j = 0; j < piLeafCount[iR]; j++){
                    ppcCopy2[j] = CKCALLOC(prHMMrght->L+strlen(ppcProfile2[0])+1, sizeof(char));
                    for (i = 0; i < (int) strlen(ppcProfile2[0]); i++){
                        ppcCopy2[j][i] = ppcProfile2[j][i];
                    }
                }

                {
                    /* the size of the elements in the forward/backward matrices 
                       depends very much on the lengths of the profiles _and_ 
                       in which position (1st/2nd) the longer/shorter profile/HMM is.
                       the matrix elements can easily exceed the size of a (long?) double
                       if the longer profile/HMM is associated with the query (q) and the 
                       shorter with the target (t). 
                       FIXME: however, pseudo-count adding may also depend on position, 
                       this is only really tested for the HMM being in the 1st position (q)
                       MUST TEST THIS MORE THOROUGHLY
                       
                       this switch appears to be most easily (although unelegantly) 
                       effected here. Don't want to do it (upstairs) in PrepareAlignment() 
                       because it might jumble up the nodes. Don't want to do it in hhalign() 
                       either because ppcProfile1/2 and q/t may be used independently.
                       FS, r236 -> r237
                    */
                    int iLenA = strlen(ppcCopy2[0]);
                    int iLenH = prHMMrght->L;
                    int iHHret = 0;

                    if (iLenH < iLenA){
                        iHHret = hhalign(ppcReprsnt2, 0/* only one representative seq */, NULL,
                                         ppcCopy2,    piLeafCount[iR], pdWeightsR,
                                         &dScore, prHMMrght, prHMMrght,
                                         NULL, NULL, NULL, NULL,
                                         rHhalignPara, &rHHscores, 
                                         iAux_FS++, /* DEBUG ARGUMENT */ rLog.iLogLevelEnabled,
                                         zcAux, zcError);
                    }
                    else {
                        iHHret = hhalign(ppcCopy2,    piLeafCount[iR], pdWeightsR,
                                         ppcReprsnt2, 0/* only one representative seq */, NULL,
                                         &dScore, prHMMrght, prHMMrght,
                                         NULL, NULL, NULL, NULL,
                                         rHhalignPara, &rHHscores, 
                                         iAux_FS++, /* DEBUG ARGUMENT */ rLog.iLogLevelEnabled,
                                         zcAux, zcError);
                    }
                    if ( (0 != iHHret) && (rLog.iLogLevelEnabled <= LOG_VERBOSE) ){ /* FS, r255 -> */
                        fprintf(stderr, "%s:%s:%d: (not essential) HMM pre-alignment failed, error %d, \n"
                                "\t#=%d (len=%d), lead-seq=%s, len(HMM)=%d\n%s\nCARRY ON REGARDLESS\n", 
                                __FUNCTION__, __FILE__, __LINE__, iHHret, 
                                piLeafCount[iR], (int)strlen(ppcCopy2[0]), prMSeq->sqinfo[ppiLeafList[iR][0]].name, 
                                (int)strlen(ppcReprsnt2[0]), zcError);
                    }
                }
                pdScores[ppiLeafList[iR][0]] = dScore;
#if 0
                printf("H: %s\nR: %s\nscore: %f\n",
                       ppcReprsnt2[0], ppcCopy2[0], dScore);
#endif
                /* assemble 'consensus';
                 * this is not a real consensus, it is more a gap indicator,
                 * for each position it consists of residues/gaps in the 1st sequences,
                 * or a residue (if any) of the other sequences.
                 * it only contains a gap if all sequences of the profile
                 * have a gap at this position
                 */
                pcConsens2 = CKCALLOC(prHMMrght->L+strlen(ppcProfile2[0])+1, sizeof(char));
                for (i = 0; i < prHMMrght->L; i++){
                    for (j = 0, pcConsens2[i] = '-'; (j < piLeafCount[iR]) && ('-' == pcConsens2[i]); j++){
                        pcConsens2[i] = ppcCopy2[j][i];
                    }
                }
#if 0
                for (j = 0; (j < piLeafCount[iR]); j++){
                    printf("R%d:%s\n", j, ppcCopy2[j]);
                }
                printf("RC:%s\n", pcConsens2);
#endif
            } /*  ( (1 == piLeafCount[iR]) && (0 != prHMM->L) ) */

            

            /* do alignment here (before free)
             */
            {
                /* the size of the elements in the forward/backward matrices 
                   depends very much on the lengths of the profiles _and_ 
                   in which position (1st/2nd) the longer/shorter profile is.
                   the matrix elements can easily exceed the size of a (long?) double
                   if the longer profile is associated with the query (q) and the 
                   shorter with the target (t). 
                   this switch appears to be most easily (although unelegantly) 
                   effected here. Don't want to do it (upstairs) in PrepareAlignment() 
                   because it might jumble up the nodes. Don't want to do it in hhalign() 
                   either because ppcProfile1/2 and q/t may be used independently.
                   FS, r228 -> 229
                 */
                int iLen1 = strlen(ppcProfile1[0]);
                int iLen2 = strlen(ppcProfile2[0]);
                /* potential problem with empty profiles, FS, r249 -> r250 */
                if ( (0 == iLen1) || (0 == iLen2) ){
                    Log(&rLog, LOG_FATAL, "strlen(prof1)=%d, strlen(prof2)=%d -- nothing to align\n", 
                        iLen1, iLen2);
                }

                if (iLen1 < iLen2){
                    int iHHret = 0;
                    int iOldMacRam = rHhalignPara.iMacRamMB;
                    iHHret = hhalign(ppcProfile1, piLeafCount[iL], pdWeightsL,
                                     ppcProfile2, piLeafCount[iR], pdWeightsR,
                                     &dScore, prHMMleft, prHMMrght,
                                     pcConsens1, pcReprsnt1,
                                     pcConsens2, pcReprsnt2,
                                     rHhalignPara, &rHHscores, 
                                     iAux_FS++, /* DEBUG ARGUMENT */ rLog.iLogLevelEnabled,
                                     zcAux, zcError);

                    if (RETURN_OK != iHHret){ /* FS, r241 -> */
                        /*fprintf(stderr, "%s:%d: emergency EXIT\n", __FILE__, __LINE__); exit(-1);*/
                        fprintf(stderr, 
                                "%s:%s:%d: problem in alignment (profile sizes: %d + %d) (%s + %s), forcing Viterbi\n"
                                "\thh-error-code=%d (mac-ram=%d)\n%s",
                                __FUNCTION__, __FILE__, __LINE__, piLeafCount[iL], piLeafCount[iR],
                                prMSeq->sqinfo[ppiLeafList[iL][0]].name, prMSeq->sqinfo[ppiLeafList[iR][0]].name,
                                iHHret, rHhalignPara.iMacRamMB, zcError);
                        /* at this stage hhalign() has failed, 
                           the only thing we can do (easily) is to re-run it in Viterbi mode, 
                           for this set MAC-RAM=0, set it back to its original value after 2nd try. 
                           FS, r241 -> r243 */

                        if (RETURN_FROM_MAC == iHHret){
                            /* Note: the default way to run hhalign() is to initially select MAC 
                               by giving it all the memory it needs. MAC may fail due to overflow (repeats?).
                               alternatively, the problem may be (genuinely) too big for MAC.
                               in thses cases it is legitimate to switch to Viterbi. 
                               However, selecting Viterbi from the outset is an abuse (abomination!), 
                               should this 1st invocation of Viterbi fail, then we (FS) will overrule 
                               the user and hammer the system with a massive memory request. 
                               (Jos 2:19) If anyone goes outside your house into the street, 
                               his blood will be on his own head; we will not be responsible. FS, r246 -> r247 */
                            rHhalignPara.iMacRamMB = 0;
                        }
                        else {
                            rHhalignPara.iMacRamMB = REALLY_BIG_MEMORY_MB;
                        }

                        iHHret = hhalign(ppcProfile1, piLeafCount[iL], pdWeightsL,
                                         ppcProfile2, piLeafCount[iR], pdWeightsR,
                                         &dScore, prHMMleft, prHMMrght,
                                         pcConsens1, pcReprsnt1,
                                         pcConsens2, pcReprsnt2,
                                         rHhalignPara, &rHHscores, 
                                         iAux_FS++, /* DEBUG ARGUMENT */ rLog.iLogLevelEnabled,
                                         zcAux, zcError);
                        if (RETURN_OK != iHHret){ /* at this stage hhalign() has failed twice, 
                                             1st time MAC, 2nd time Viterbi, don't know what to do else. FS, r241 -> r243 */
                            fprintf(stderr, "%s:%s:%d: problem in alignment, Viterbi did not work\n"
                                    "\thh-error-code=%d (mac-ram=%d)\n%s",
                                    __FUNCTION__, __FILE__, __LINE__, iHHret, rHhalignPara.iMacRamMB, zcError);
                            Log(&rLog, LOG_FATAL, "could not perform alignment -- bailing out\n");
                        }
                        else {
                            fprintf(stderr, "%s:%s:%d: 2nd attempt worked", __FUNCTION__, __FILE__, __LINE__);
                        }
                        rHhalignPara.iMacRamMB = iOldMacRam; 
                    } /* 1st invocation failed */

                } /* 1st profile was shorter than 2nd */
                else {
                    int iHHret = 0;
                    int iOldMacRam = rHhalignPara.iMacRamMB;
                    iHHret = hhalign(ppcProfile2, piLeafCount[iR], pdWeightsR,
                                     ppcProfile1, piLeafCount[iL], pdWeightsL,
                                     &dScore, prHMMrght, prHMMleft,
                                     pcConsens2, pcReprsnt2,
                                     pcConsens1, pcReprsnt1,
                                     rHhalignPara, &rHHscores, 
                                     iAux_FS++, /* DEBUG ARGUMENT */ rLog.iLogLevelEnabled,
                                     zcAux, zcError);

                    if (RETURN_OK != iHHret){ /* FS, r241 -> r243 */
                        /*fprintf(stderr, "%s:%d: emergency EXIT\n", __FILE__, __LINE__); exit(-1);*/
                        fprintf(stderr, 
                                "%s:%s:%d: problem in alignment (profile sizes: %d + %d) (%s + %s), forcing Viterbi\n"
                                "\thh-error-code=%d (mac-ram=%d)\n%s",
                                __FUNCTION__, __FILE__, __LINE__, piLeafCount[iL], piLeafCount[iR],
                                prMSeq->sqinfo[ppiLeafList[iL][0]].name, prMSeq->sqinfo[ppiLeafList[iR][0]].name,
                                iHHret, rHhalignPara.iMacRamMB, zcError);
                        /* at this stage hhalign() has failed, 
                           the only thing we can do (easily) is to re-run it in Viterbi mode, 
                           for this set MAC-RAM=0, set it back to its original value after 2nd try. 
                           FS, r241 -> r243 */

                        if (RETURN_FROM_MAC == iHHret){
                            /* see above */
                            rHhalignPara.iMacRamMB = 0;
                        }
                        else {
                            rHhalignPara.iMacRamMB = REALLY_BIG_MEMORY_MB;
                        }

                        iHHret = hhalign(ppcProfile2, piLeafCount[iR], pdWeightsR,
                                         ppcProfile1, piLeafCount[iL], pdWeightsL,
                                         &dScore, prHMMrght, prHMMleft,
                                         pcConsens2, pcReprsnt2,
                                         pcConsens1, pcReprsnt1,
                                         rHhalignPara, &rHHscores, 
                                         iAux_FS++, /* DEBUG ARGUMENT */ rLog.iLogLevelEnabled,
                                         zcAux, zcError);
                        if (RETURN_OK != iHHret){ /* at this stage hhalign() has failed twice, 
                                             1st time MAC, 2nd time Viterbi, don't know what to do else. FS, r241 -> r243 */
                            fprintf(stderr, "%s:%s:%d: problem in alignment, Viterbi did not work\n"
                                    "\thh-error-code=%d (mac-ram=%d)\n%s",
                                    __FUNCTION__, __FILE__, __LINE__, iHHret, rHhalignPara.iMacRamMB, zcError);
                            Log(&rLog, LOG_FATAL, "could not perform alignment -- bailing out\n");
                        }
                        else {
                            fprintf(stderr, "%s:%s:%d: 2nd attempt worked", __FUNCTION__, __FILE__, __LINE__);
                        }
                        rHhalignPara.iMacRamMB = iOldMacRam; 
                    } /* 1st invocation failed */

                } /* 2nd profile was shorter than 1st */

                /* 
                 * at this stage have performed alignment of 2 profiles/sequences.
                 * if HMM batch information had been used then have choices:
                 * (i) if HMM info only intended for initial alignment (of sequences) then make both HMMs stale;
                 * (iia) if alignment of 2 profiles/sequences where same HMM used, then retain;
                 * (iib) if alignment of 2 profiles/sequences where different HMMs used, then make both stale;
                 * (iii) some majority voting
                 */
#if 0  /* always make HMM batch stale (after 1st invocation) */
                if ( (NULL != prMSeq->ppiHMMBindex) && (NULL != prMSeq->ppiHMMBindex[ppiLeafList[iL][0]]) ){
                    prMSeq->ppiHMMBindex[ppiLeafList[iL][0]][0] = -1;
                }
                if ( (NULL != prMSeq->ppiHMMBindex) && (NULL != prMSeq->ppiHMMBindex[ppiLeafList[iR][0]]) ){
                    prMSeq->ppiHMMBindex[ppiLeafList[iR][0]][0] = -1;
                }
#else /* retain HMMs if they were the same for both profiles */
                if (NULL != prMSeq->ppiHMMBindex) {
                    int i;

                    if ( (NULL != prMSeq->ppiHMMBindex[ppiLeafList[iL][0]]) && (NULL != prMSeq->ppiHMMBindex[ppiLeafList[iR][0]]) ){
                        if ( prMSeq->ppiHMMBindex[ppiLeafList[iL][0]][0] == -1){
                            prMSeq->ppiHMMBindex[ppiLeafList[iR][0]][0] = -1; /* this is conservative, could do H[iL] = H[iR] */
                            for (i = 0; i < piLeafCount[iR]; i++){
                                prMSeq->ppiHMMBindex[ppiLeafList[iR][i]][0] = -1;
                            }
                        }
                        else if ( prMSeq->ppiHMMBindex[ppiLeafList[iR][0]][0] == -1){
                            prMSeq->ppiHMMBindex[ppiLeafList[iL][0]][0] = -1; /* this is conservative, could do H[iR] = H[iL] */
                            for (i = 0; i < piLeafCount[iL]; i++){
                                prMSeq->ppiHMMBindex[ppiLeafList[iL][i]][0] = -1;
                            }
                        }
                        else if (prMSeq->ppiHMMBindex[ppiLeafList[iL][0]][0] != prMSeq->ppiHMMBindex[ppiLeafList[iR][0]][0]){
                            prMSeq->ppiHMMBindex[ppiLeafList[iL][0]][0] = -1; /* this is NOT conservative, mandatory */
                            for (i = 0; i < piLeafCount[iL]; i++){
                                prMSeq->ppiHMMBindex[ppiLeafList[iL][i]][0] = -1;
                            }
                            prMSeq->ppiHMMBindex[ppiLeafList[iR][0]][0] = -1; /* this is NOT conservative, mandatory */
                            for (i = 0; i < piLeafCount[iR]; i++){
                                prMSeq->ppiHMMBindex[ppiLeafList[iR][i]][0] = -1;
                            }
                        }
                        else {
                            /* void, HMMs should be same */
                        }
                    }
                } /* there was a HMM batch */
#endif
                
                if (rLog.iLogLevelEnabled <= LOG_DEBUG){
                    int i;

                    printf("@@iL=%d, #(iL)=%d, iR=%d, #(iR)=%d\n", iL, piLeafCount[iL], iR, piLeafCount[iR]);
                    for (i = 0; i < piLeafCount[iL]; i++){
                        char *pc = ppcProfile1[i];
                        printf("@@>%s\n", prMSeq->sqinfo[ppiLeafList[iL][i]].name);
                        printf("@@");
                        while('\0' != *pc){
                            printf("%c", toupper(*pc));
                            pc++;
                        }
                        printf("\n");
                    }
                    for (i = 0; i < piLeafCount[iR]; i++){
                        char *pc = ppcProfile2[i];
                        printf("@@>%s\n", prMSeq->sqinfo[ppiLeafList[iR][i]].name);
                        printf("@@");
                        while('\0' != *pc){
                            printf("%c", toupper(*pc));
                            pc++;
                        }
                        printf("\n");
                    }
                    printf("\n");
                } /* LOG_DEBUG */    

            }
            /* free left/right node lists,
             * after alignment left/right profiles no longer needed
             */
            if (NULL != ppcCopy1){
                int i;
                for (i = 0; i <  piLeafCount[iL]; i++){
                    CKFREE(ppcCopy1[i]);
                }
                CKFREE(ppcCopy1);
                CKFREE(pcReprsnt1);
                CKFREE(pcConsens1);
            }
            if (NULL != ppcCopy2){
                int i;
                for (i = 0; i <  piLeafCount[iR]; i++){
                    CKFREE(ppcCopy2[i]);
                }
                CKFREE(ppcCopy2);
                CKFREE(pcReprsnt2);
                CKFREE(pcConsens2);
            }
            ppiLeafList[iL] = CKFREE(ppiLeafList[iL]);
            ppiLeafList[iR] = CKFREE(ppiLeafList[iR]);
            piLeafCount[iL] = piLeafCount[iR] = 0;
            
        } /* was a merge node */

        if (rLog.iLogLevelEnabled <= LOG_DEBUG){
            int i, j;
            FILE *fp = LogGetFP(&rLog, LOG_DEBUG);
            for (i = 0; i < iNodeCount; i++){
                if (0 == piLeafCount[i]){
                    continue;
                }
                fprintf(fp, "node %3d, #leaves=%d:\t", i, piLeafCount[i]);
                for (j = 0; ppiLeafList && (j < piLeafCount[i]); j++){
                    fprintf(fp, "%d,", ppiLeafList[i][j]);
                }
                fprintf(fp, "\n");
            }
        }


    } /* 0 <= iN < iNodeCount */
    ProgressDone(prProgress);


    /* check length and set length info
     */
    iAlnLen = strlen(prMSeq->seq[0]);
    for (i=0; i<prMSeq->nseqs; i++) {
#if 0
        Log(&rLog, LOG_FORCED_DEBUG, "seq no %d: name %s; len %d; %s",
                  i, prMSeq->sqinfo[i].name, strlen(prMSeq->seq[i]), prMSeq->seq[i]);
#endif
        
#ifndef NDEBUG
        assert(iAlnLen == strlen(prMSeq->seq[i]));
#endif
        prMSeq->sqinfo[i].len = iAlnLen;
    }
    prMSeq->aligned = TRUE;


    if (rLog.iLogLevelEnabled <= LOG_DEBUG){
        if (0 != prHMM->L){
            int i;
            Log(&rLog, LOG_DEBUG, "Alignment scores with HMM:");
            for (i = 0; /*pdScores[i] > 0.0*/i < prMSeq->nseqs; i++){
                Log(&rLog, LOG_DEBUG, "%2d:\t%f\n", i, pdScores[i]);
            }
        }
    }


    /** translate back ambiguity residues
     * hhalign translates ambiguity codes (B,Z) into unknown residues (X).
     * as we still have the original input we can substitute them back
     */
    TranslateUnknown2Ambiguity(prMSeq);
    ReAttachLeadingGaps(prMSeq, iProfProfSeparator);

#if 0 /* DEVEL 291 */
    if (NULL == prHMMList){
        CKFREE(prHMM);
    }
#else
    CKFREE(prHMMnull);
#endif
    CKFREE(ppcProfile2);
    CKFREE(ppcProfile1);
    CKFREE(ppiLeafList[piOrderLR[DIFF_NODE*(iNodeCount-1)+PRNT_NODE]]);
    CKFREE(ppiLeafList);
    CKFREE(piLeafCount);
    CKFREE(pdScores);
    FreeProgress(&prProgress);
    CKFREE(pdWeightsL);
    CKFREE(pdWeightsR);

#if TIMING
    StopwatchStop(stopwatch);
    StopwatchDisplay(stdout, "Total time for HHalignWrapper():" , stopwatch);
    StopwatchFree(stopwatch);
#endif

    return dScore; /* FIXME alternative: return averaged pdScores */

}
/***   end: HHalignWrapper() ***/
