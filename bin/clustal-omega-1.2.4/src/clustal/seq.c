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
 *  RCS $Id: seq.c 298 2014-11-07 12:18:36Z fabian $
 *
 *
 * Module for sequence/alignment IO and misc.
 *
 * This depends heavily on Sean Eddy's squid library, which is obsoleted by
 * HMMER3's Easel. However, easel doesn't support that many non-aligned input
 * formats.
 *
 */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <assert.h>
#include "squid/squid.h"
#include <ctype.h>

#include "util.h"
#include "log.h"
#include "seq.h"
/*#include "../../mymemmonitor.h"*/

#define ALLOW_ONLY_PROTEIN 0 // DD



/**
 * @brief Stripped down version of squid's alistat
 *
 *
 * @param[in] prMSeq
 * The alignment to analyse
 * @param[in] bSampling
 * For many sequences: samples from pool
 * @param[in] bReportAll
 * Report identities for all sequence pairs
 *
 * Don't have to worry about sequence case because our version of PairwiseIdentity is case insensitive
 */
void
AliStat(mseq_t *prMSeq, bool bSampling, bool bReportAll) {

    /*
     * bSampling = squid's do_fast
     * bReportAll = squid's allreport
     */
    float  **ppdIdentMx;  /* identity matrix (squid: imx) */
    const int iNumSample = 1000; /* sample size (squid: nsample) */


    MSA *msa; /* squid's alignment structure */
    int small, large;	
    int bestj, worstj;
    float sum;
    float worst_worst, worst_best, best_best;
    float avgid;
    int i, j;
    int nres; /* number of residues */

    if (bSampling && bReportAll) {
        Log(&rLog, LOG_WARN,
            "Cannot report all and sample at the same time. Skipping %s()", __FUNCTION__);
        return;
    }
    if (FALSE == prMSeq->aligned) {
        Log(&rLog, LOG_WARN,
            "Sequences are not aligned. Skipping %s()", __FUNCTION__);
        return;
    }

    /* silence gcc warnings about uninitialized variables
     */
    worst_worst = worst_best = best_best = 0.0;
    bestj = worstj = -1;


    /** mseq to squid msa
     *
     * FIXME code overlap with WriteAlignment. Make it a function and take
     * code there (contains more comments) as template
     *
     */
    msa  = MSAAlloc(prMSeq->nseqs, 
                    /* derive alignment length from first seq */
                    strlen(prMSeq->seq[0]));
    for (i=0; i<prMSeq->nseqs; i++) {
        int key; /* MSA struct internal index for sequence */
        char *this_name = prMSeq->sqinfo[i].name; /* prMSeq sequence name */
        char *this_seq = prMSeq->seq[i]; /* prMSeq sequence */
        SQINFO *this_sqinfo = &prMSeq->sqinfo[i]; /* prMSeq sequence name */

        key = GKIStoreKey(msa->index, this_name);
        msa->sqname[key] = sre_strdup(this_name, strlen(this_name));
        /* setting msa->sqlen[idx] and msa->aseq[idx] */
        msa->sqlen[key] = sre_strcat(&(msa->aseq[key]), msa->sqlen[key],
                                     this_seq, strlen(this_seq));
        if (this_sqinfo->flags & SQINFO_DESC) {
            MSASetSeqDescription(msa, key, this_sqinfo->desc);
        }           
        msa->nseq++;
    }    
    


    nres = 0;
    small = large = -1;
    for (i = 0; i < msa->nseq; i++) {
        int rlen;		/* raw sequence length           */
        rlen  = DealignedLength(msa->aseq[i]);
        nres +=  rlen;
        if (small == -1 || rlen < small) small = rlen;
        if (large == -1 || rlen > large) large = rlen;
    }


    if (bSampling) {
        avgid = AlignmentIdentityBySampling(msa->aseq, msa->alen, 
                                            msa->nseq, iNumSample);

	} else {
        float best, worst;

        /* this might be slow...could use openmp inside squid */
        MakeIdentityMx(msa->aseq, msa->nseq, &ppdIdentMx);
        if (bReportAll) {
            printf("  %-15s %5s %7s %-15s %7s %-15s\n",
                   "NAME", "LEN", "HIGH ID", "(TO)", "LOW ID", "(TO)");
            printf("  --------------- ----- ------- --------------- ------- ---------------\n");
        }

        sum = 0.0;
        worst_best  = 1.0;
        best_best   = 0.0;
        worst_worst = 1.0;
        for (i = 0; i < msa->nseq; i++) {
            worst = 1.0;
            best  = 0.0;
            for (j = 0; j < msa->nseq; j++) {
                /* closest seq to this one = best */
                if (i != j && ppdIdentMx[i][j] > best)  {
                    best  = ppdIdentMx[i][j]; bestj = j; 
                }
                if (ppdIdentMx[i][j] < worst) {
                    worst = ppdIdentMx[i][j]; worstj = j; 
                }
            }

            if (bReportAll)  {
                printf("* %-15s %5d %7.1f %-15s %7.1f %-15s\n",
                       msa->sqname[i], DealignedLength(msa->aseq[i]),
                       best * 100.,  msa->sqname[bestj],
                       worst * 100., msa->sqname[worstj]);
            }
            if (best > best_best)    best_best = best;
            if (best < worst_best)   worst_best = best;
            if (worst < worst_worst) worst_worst = worst;
            for (j = 0; j < i; j++)
                sum += ppdIdentMx[i][j];
	    }
        avgid = sum / (float) (msa->nseq * (msa->nseq-1)/2.0);
        if (bReportAll)
            puts("");
        FMX2Free(ppdIdentMx);
    } /* else bSampling */


    
    /* Print output
     */
    if (msa->name != NULL)
        printf("Alignment name:      %s\n", msa->name); 
    /*printf("Format:              %s\n",     SeqfileFormat2String(afp->format));*/
    printf("Number of sequences: %d\n", msa->nseq);
    printf("Total # residues:    %d\n", nres);
    printf("Smallest:            %d\n", small);
    printf("Largest:             %d\n", large);
    printf("Average length:      %.1f\n", (float) nres / (float) msa->nseq);
    printf("Alignment length:    %d\n", msa->alen);
    printf("Average identity:    %.2f%%\n", 100.*avgid);

    if (! bSampling) {
        printf("Most related pair:   %.2f%%\n", 100.*best_best);
        printf("Most unrelated pair: %.2f%%\n", 100.*worst_worst);
        printf("Most distant seq:    %.2f%%\n", 100.*worst_best);
    }

    /*
	char *cs;
	cs = MajorityRuleConsensus(msa->aseq, msa->nseq, msa->alen);
    printf cs;
    */

    MSAFree(msa);
}
/* end of AliStat() */





/**
 * @brief Shuffle mseq order
 *
 * @param[out] mseq
 * mseq structure to shuffle
 *
 */
void
ShuffleMSeq(mseq_t *mseq)
{
    int iSeqIndex;
    int *piPermArray;


    /* quick and dirty: create an array with permuted indices and
     * swap accordingly (which amounts to shuffling twice)
     */

    PermutationArray(&piPermArray, mseq->nseqs);

    for (iSeqIndex=0; iSeqIndex<mseq->nseqs; iSeqIndex++) {    
        SeqSwap(mseq, piPermArray[iSeqIndex], iSeqIndex);
    }

    CKFREE(piPermArray);
}
/***   end: ShuffleMSeq()   ***/


/**
 * @brief Swap two sequences in an mseq_t structure.
 *
 * @param[out] prMSeq
 * Multiple sequence struct
 * @param[in] i
 * Index of first sequence
 * @param[in] j
 * Index of seconds sequence
 *
 * 
 */    
void
SeqSwap(mseq_t *prMSeq, int i, int j)
{
    char *pcTmp;
    SQINFO rSqinfoTmp;

    assert(NULL!=prMSeq);
    assert(i<prMSeq->nseqs && j<prMSeq->nseqs);

    if (i==j) {
        return;
    }
    
    pcTmp = prMSeq->seq[i];
    prMSeq->seq[i] = prMSeq->seq[j];
    prMSeq->seq[j] = pcTmp;

    pcTmp = prMSeq->orig_seq[i];
    prMSeq->orig_seq[i] = prMSeq->orig_seq[j];
    prMSeq->orig_seq[j] = pcTmp;

    rSqinfoTmp = prMSeq->sqinfo[i];
    prMSeq->sqinfo[i] = prMSeq->sqinfo[j];
    prMSeq->sqinfo[j] = rSqinfoTmp;
    
    return; 
}
/***   end: SeqSwap()   ***/




/**
 * @brief Dealigns all sequences in mseq structure, updates the
 * sequence length info and sets aligned to FALSE
 *
 * @param[out] mseq
 * The mseq structure to dealign
 * 
 */    
void
DealignMSeq(mseq_t *mseq)
{
    int i; /* aux */
    for (i=0; i<mseq->nseqs; i++) {
        DealignSeq(mseq->seq[i]);
        mseq->sqinfo[i].len = strlen(mseq->seq[i]);
    }
    mseq->aligned = FALSE;
    
    return; 
}
/***   end: DealinMSeq()   ***/



/**
 * @brief debug output of sqinfo struct
 *
 * @param[in] sqinfo
 * Squid's SQINFO struct for a certain seqeuence
 *
 * @note useful for debugging only
 * 
 */    
void
LogSqInfo(SQINFO *sqinfo)
{
    
    /*LOG_DEBUG("sqinfo->flags = %d", sqinfo->flags);*/
    if (sqinfo->flags & SQINFO_NAME)
        Log(&rLog, LOG_FORCED_DEBUG, "sqinfo->name = %s", sqinfo->name);

    if (sqinfo->flags & SQINFO_ID)
        Log(&rLog, LOG_FORCED_DEBUG, "sqinfo->id = %s", sqinfo->id);
        
    if (sqinfo->flags & SQINFO_ACC)
        Log(&rLog, LOG_FORCED_DEBUG, "sqinfo->acc = %s", sqinfo->acc);
        
    if (sqinfo->flags & SQINFO_DESC)
        Log(&rLog, LOG_FORCED_DEBUG, "sqinfo->desc = %s", sqinfo->desc);

    if (sqinfo->flags & SQINFO_LEN)
        Log(&rLog, LOG_FORCED_DEBUG, "sqinfo->len = %d", sqinfo->len);

    if (sqinfo->flags & SQINFO_START)
        Log(&rLog, LOG_FORCED_DEBUG, "sqinfo->start = %d", sqinfo->start);
        
    if (sqinfo->flags & SQINFO_STOP)
        Log(&rLog, LOG_FORCED_DEBUG, "sqinfo->stop = %d", sqinfo->stop);

    if (sqinfo->flags & SQINFO_OLEN)
        Log(&rLog, LOG_FORCED_DEBUG, "sqinfo->olen = %d", sqinfo->olen);

    if (sqinfo->flags & SQINFO_TYPE)
        Log(&rLog, LOG_FORCED_DEBUG, "sqinfo->type = %d", sqinfo->type);
        
    if (sqinfo->flags & SQINFO_SS) 
        Log(&rLog, LOG_FORCED_DEBUG, "sqinfo->ss = %s", sqinfo->ss);
    
    if (sqinfo->flags & SQINFO_SA)
        Log(&rLog, LOG_FORCED_DEBUG, "sqinfo->sa = %s", sqinfo->sa);
}
/***   end: log_sqinfo   ***/


/**
 * @brief convert int-encoded iSeqType to string
 *
 * @param[in] iSeqType int-encoded sequence type
 *
 * @return character pointer describing the sequence type
 *
 */
const char*
SeqTypeToStr(int iSeqType)
{
    switch (iSeqType)  {
    case SEQTYPE_DNA:
        return "DNA";
    case SEQTYPE_RNA:
        return "RNA";
    case SEQTYPE_PROTEIN:
        return "Protein";
    case SEQTYPE_UNKNOWN:
        return "UNKNOWN";
    default:
        Log(&rLog, LOG_FATAL, "Internal error in %s", __FUNCTION__);
    }
    return "Will never get here";
}
/***   end: SeqTypeToStr   ***/



/**
 * @brief reads sequences from file
 *
 * @param[out] prMSeq
 * Multiple sequence struct. Must be preallocated.
 * FIXME: would make more sense to allocate it here.
 * @param[in] seqfile
 * Sequence file name. If '-' sequence will be read from stdin. 
 * @param[in] iSeqType
 * int-encoded sequence type. Set to
 * SEQTYPE_UNKNOWN for autodetect (guessed from first sequence)
 * @param[in] iMaxNumSeq
 * Return an error, if more than iMaxNumSeq have been read
 * @param[in] iMaxSeqLen 
 * Return an error, if a seq longer than iMaxSeqLen has been read
 *
 * @return 0 on success, -1 on error
 *
 * @note
 *  - Depends heavily on squid
 *  - Sequence file format will be guessed
 *  - If supported by squid, gzipped files can be read as well.
 */
int
ReadSequences(mseq_t *prMSeq, char *seqfile, 
              int iSeqType, int iSeqFmt, bool bIsProfile, bool bDealignInputSeqs, 
              int iMaxNumSeq, int iMaxSeqLen, char *pcHMMBatch)
{
    SQFILE *dbfp; /* sequence file descriptor */
    char *cur_seq = NULL;
    SQINFO cur_sqinfo = {0};
    int iSeqIdx; /* sequence counter */
    int iSeqPos; /* sequence string position counter */
    
    assert(NULL!=seqfile);

    
    /* Try to work around inability to autodetect from a pipe or .gz:
     * assume FASTA format
     */
    if (SQFILE_UNKNOWN == iSeqFmt  &&
        (Strparse("^.*\\.gz$", seqfile, 0) || strcmp(seqfile, "-") == 0)) {
        iSeqFmt = SQFILE_FASTA;
    }
  
    /* Using squid routines to read input. taken from seqstat_main.c. we don't
     * know if input is aligned, so we use SeqfileOpen instead of MSAFileOpen
     * etc. NOTE this also means we discard some information, e.g. when
     * reading from and writing to a stockholm file, all extra MSA
     * info/annotation will be lost.
     *
     */

    if (NULL == (dbfp = SeqfileOpen(seqfile, iSeqFmt, NULL))) {
        Log(&rLog, LOG_ERROR, "Failed to open sequence file %s for reading", seqfile);
        return -1;
    }

    
    /* FIXME squid's ReadSeq() will exit with fatal error if format is
     * unknown. This will be a problem for a GUI. Same is true for many squid
     * other functions.
     *
     * The original squid:ReadSeq() dealigns sequences on input. We
     * use a patched version.
     *
     */
    while (ReadSeq(dbfp, dbfp->format,
                   &cur_seq,
                   &cur_sqinfo)) { 
        
        if (prMSeq->nseqs+1>iMaxNumSeq) {
            Log(&rLog, LOG_ERROR, "Maximum number of sequences (=%d) exceeded after reading sequence '%s' from '%s'",
                  iMaxNumSeq, cur_sqinfo.name, seqfile);
            return -1;
        }
        if ((int)strlen(cur_seq)>iMaxSeqLen) {
            Log(&rLog, LOG_ERROR, "Sequence '%s' has %d residues and is therefore longer than allowed (max. sequence length is %d)",
                  cur_sqinfo.name, strlen(cur_seq), iMaxSeqLen);
            return -1;
        }
        if ((int)strlen(cur_seq)==0) {
            Log(&rLog, LOG_ERROR, "Sequence '%s' has 0 residues",
                cur_sqinfo.name);
            return -1;
        }
            
        /* FIXME: use modified version of AddSeq() that allows handing down SqInfo
         */

        prMSeq->seq =  (char **)
            CKREALLOC(prMSeq->seq, (prMSeq->nseqs+1) * sizeof(char *));
        prMSeq->seq[prMSeq->nseqs] = CkStrdup(cur_seq);
        

        prMSeq->sqinfo =  (SQINFO *)
            CKREALLOC(prMSeq->sqinfo, (prMSeq->nseqs+1) * sizeof(SQINFO)); 
        memset(&prMSeq->sqinfo[prMSeq->nseqs], 0, sizeof(SQINFO));
        SeqinfoCopy(&prMSeq->sqinfo[prMSeq->nseqs], &cur_sqinfo);

#ifdef TRACE
        Log(&rLog, LOG_FORCED_DEBUG, "seq no %d: seq = %s", prMSeq->nseqs, prMSeq->seq[prMSeq->nseqs]);
        LogSqInfo(&prMSeq->sqinfo[prMSeq->nseqs]);
#endif
        /* always guess type from first seq. use squid function and
         * convert value
         */
        if (0 == prMSeq->nseqs) {
            int type = Seqtype(prMSeq->seq[prMSeq->nseqs]);
            switch (type)  {
            case kDNA:
                prMSeq->seqtype = SEQTYPE_DNA;
                break;
            case kRNA:
                prMSeq->seqtype = SEQTYPE_RNA;
                break;
            case kAmino:
                prMSeq->seqtype = SEQTYPE_PROTEIN;
                break;
            case kOtherSeq:
                prMSeq->seqtype = SEQTYPE_UNKNOWN;
                break;
            default:
                Log(&rLog, LOG_FATAL, "Internal error in %s", __FUNCTION__);
            }

            /* override with given sequence type but check with
             * automatically detected type and warn if necessary
             */
            if (SEQTYPE_UNKNOWN != iSeqType) {
                if (prMSeq->seqtype != iSeqType) { 
                    Log(&rLog, LOG_WARN, "Overriding automatically determined seq-type %s to %s as requested",
                         SeqTypeToStr(prMSeq->seqtype), SeqTypeToStr(iSeqType));
                    prMSeq->seqtype = iSeqType;
                }
            }
            /* if type could not be determined and was not set return error */
            if (SEQTYPE_UNKNOWN == iSeqType && SEQTYPE_UNKNOWN == prMSeq->seqtype) {
                Log(&rLog, LOG_ERROR, "Couldn't guess sequence type from first sequence");
                FreeSequence(cur_seq, &cur_sqinfo);        
                SeqfileClose(dbfp);
                return -1;
            }
        }

        Log(&rLog, LOG_DEBUG, "seq-no %d: type=%s name=%s len=%d seq=%s",
             prMSeq->nseqs, SeqTypeToStr(prMSeq->seqtype),
             prMSeq->sqinfo[prMSeq->nseqs].name, prMSeq->sqinfo[prMSeq->nseqs].len,
             prMSeq->seq[prMSeq->nseqs]);
        
        /* FIXME IPUAC and/or case conversion? If yes see
         * corresponding squid functions. Special treatment of
         * Stockholm tilde-gaps for ktuple code?
         */

        prMSeq->nseqs++;

        FreeSequence(cur_seq, &cur_sqinfo); 
    }
    SeqfileClose(dbfp);

/*#if ALLOW_ONLY_PROTEIN
    if (SEQTYPE_PROTEIN != prMSeq->seqtype) {
        Log(&rLog, LOG_FATAL, "Sequence type is %s. %s only works on protein.",
              SeqTypeToStr(prMSeq->seqtype), PACKAGE_NAME);
    }
#endif*/
    
    /* Check if sequences are aligned */
    prMSeq->aligned = SeqsAreAligned(prMSeq, bIsProfile, bDealignInputSeqs);

    
    /* keep original sequence as copy and convert "working" sequence
     *
     */
    prMSeq->orig_seq = (char**) CKMALLOC(prMSeq->nseqs * sizeof(char *));
    for (iSeqIdx=0; iSeqIdx<prMSeq->nseqs; iSeqIdx++) {

        prMSeq->orig_seq[iSeqIdx] = CkStrdup(prMSeq->seq[iSeqIdx]);

        
        /* convert unknown characters according to set seqtype
         * be conservative, i.e. don't allow any fancy ambiguity
         * characters to make sure that ktuple code etc. works.
         */

        /* first on the fly conversion between DNA and RNA
         */
        if (prMSeq->seqtype==SEQTYPE_DNA)
            ToDNA(prMSeq->seq[iSeqIdx]);
        if (prMSeq->seqtype==SEQTYPE_RNA)
            ToRNA(prMSeq->seq[iSeqIdx]);
        
        /* then check of each character
         */
        for (iSeqPos=0; iSeqPos<(int)strlen(prMSeq->seq[iSeqIdx]); iSeqPos++) {
            char *res = &(prMSeq->seq[iSeqIdx][iSeqPos]);
            if (isgap(*res))
                continue;
            
            if (prMSeq->seqtype==SEQTYPE_PROTEIN) {
                if (NULL == strchr(AMINO_ALPHABET, toupper(*res))) {
                    *res = AMINOACID_ANY;
                }
            } else if (prMSeq->seqtype==SEQTYPE_DNA) {                    
                if (NULL == strchr(DNA_ALPHABET, toupper(*res))) {
                    *res = NUCLEOTIDE_ANY;
                }
            } else if (prMSeq->seqtype==SEQTYPE_RNA) {
                if (NULL == strchr(RNA_ALPHABET, toupper(*res))) {
                    *res = NUCLEOTIDE_ANY;
                }
            }
        }
    }

    /* order in which sequences appear in guide-tree 
     * only allocate if different output-order desired */
    prMSeq->tree_order = NULL;

    prMSeq->filename = CkStrdup(seqfile);
    Log(&rLog, LOG_INFO, "Read %d sequences (type: %s) from %s",
         prMSeq->nseqs, SeqTypeToStr(prMSeq->seqtype), prMSeq->filename);

    prMSeq->pppcHMMBNames = NULL;
    prMSeq->ppiHMMBindex = NULL;

    /* read HMM-batch file if existent */
    if (NULL != pcHMMBatch) {

        enum {MAXLINE=10000};
        FILE *pfHMMBatch = NULL;
        char zcScanline[MAXLINE] = {0};
        char *pcToken = NULL;
        char *pcSeqName = NULL;
        int iSeq = 0;

        /* check that file exists */
        if (NULL == (pfHMMBatch = fopen(pcHMMBatch, "r"))){
            Log(&rLog, LOG_ERROR, "Failed to open HMM-batch file %s for reading", pcHMMBatch);
            return -1;
        }

        /* initialise names and indices */
        prMSeq->pppcHMMBNames = (char ***)CKMALLOC(prMSeq->nseqs * sizeof(char **));
        for (iSeq = 0; iSeq < prMSeq->nseqs; iSeq++){
            prMSeq->pppcHMMBNames[iSeq] = NULL;
        }
        prMSeq->ppiHMMBindex = (int **)CKMALLOC(prMSeq->nseqs * sizeof(int *));
        for (iSeq = 0; iSeq < prMSeq->nseqs; iSeq++){
            prMSeq->ppiHMMBindex[iSeq] = (int *)CKMALLOC(1 * sizeof(int));
            prMSeq->ppiHMMBindex[iSeq][0] = -1;
        }

        /* read batch file line-by-line */
        while (NULL != fgets(zcScanline, MAXLINE, pfHMMBatch)){

            pcToken = strtok(zcScanline, " \040\t\n");
            if (NULL == pcToken){
                continue;
            }
            else {
                pcSeqName = pcToken;
            }
            /* identify sequence label from batch file in labels read from sequence file */
            for (iSeq = 0; iSeq < prMSeq->nseqs; iSeq++){
                int iHMM = 0;

                if (0 == strcmp(pcSeqName, prMSeq->sqinfo[iSeq].name)){

                    while (NULL != (pcToken = strtok(NULL, " \040\t\n"))){
                        prMSeq->pppcHMMBNames[iSeq] = (char **)CKREALLOC(prMSeq->pppcHMMBNames[iSeq], 
                                                                         (iHMM+2) * sizeof(char *));
                        prMSeq->pppcHMMBNames[iSeq][iHMM] = CkStrdup(pcToken);
                        prMSeq->ppiHMMBindex[iSeq] = (int *)CKREALLOC(prMSeq->ppiHMMBindex[iSeq], 
                                                                      (iHMM+2) * sizeof(int));
                        prMSeq->ppiHMMBindex[iSeq][iHMM] = 0;
                        iHMM++;
                        prMSeq->pppcHMMBNames[iSeq][iHMM] = NULL;
                        prMSeq->ppiHMMBindex[iSeq][iHMM] = 0;
                    }
                    break;
                }

            } /* 0 <= iSeq < prMSeq->nseqs */
            if (iSeq >= prMSeq->nseqs) {
                Log(&rLog, LOG_WARN,
                    "sequence %s not found in input sequences (%s), will be ignored",
                    pcSeqName, seqfile);
            }

        } /* !EOF */

        fclose(pfHMMBatch); pfHMMBatch = NULL;

    } /* there was a HMM batch file */
    else {
        prMSeq->pppcHMMBNames = NULL;
        prMSeq->ppiHMMBindex = NULL;
    } /* there was no HMM batch file */

    return 0;
}
/***   end: ReadSequences   ***/


/**
 * @brief allocate and initialise new mseq_t
 *
 * @param[out] prMSeq
 * newly allocated and initialised mseq_t
 *
 * @note caller has to free by calling FreeMSeq()
 *
 * @see FreeMSeq
 *
 */
void
NewMSeq(mseq_t **prMSeq)
{
    *prMSeq = (mseq_t *) CKMALLOC(1 * sizeof(mseq_t));

    (*prMSeq)->nseqs = 0;
	(*prMSeq)->seq = NULL;
	(*prMSeq)->orig_seq = NULL;
	(*prMSeq)->seqtype = SEQTYPE_UNKNOWN;
	(*prMSeq)->sqinfo = NULL;
	(*prMSeq)->filename = NULL;
	(*prMSeq)->tree_order = NULL;
    (*prMSeq)->pppcHMMBNames = NULL;
    (*prMSeq)->ppiHMMBindex = NULL;
}
/***   end: NewMSeq   ***/



/**
 * @brief copies an mseq structure
 *
 * @param[out] prMSeqDest_p
 * Copy of mseq structure
 * @param[in]  prMSeqSrc
 * Source mseq structure to copy
 *
 * @note caller has to free copy by calling FreeMSeq()
 *
 */
void
CopyMSeq(mseq_t **prMSeqDest_p, mseq_t *prMSeqSrc)
{
    int i;
    assert(prMSeqSrc != NULL && prMSeqDest_p != NULL);
    
    NewMSeq(prMSeqDest_p);
    
    (*prMSeqDest_p)->nseqs = prMSeqSrc->nseqs;
    (*prMSeqDest_p)->seqtype = prMSeqSrc->seqtype;
    if (prMSeqSrc->filename!=NULL) {
        (*prMSeqDest_p)->filename = CkStrdup(prMSeqSrc->filename);
    }

    (*prMSeqDest_p)->seq =  (char **)
        CKMALLOC((*prMSeqDest_p)->nseqs * sizeof(char *));
    (*prMSeqDest_p)->orig_seq =  (char **)
        CKMALLOC((*prMSeqDest_p)->nseqs * sizeof(char *));
    (*prMSeqDest_p)->sqinfo =  (SQINFO *)
        CKMALLOC((*prMSeqDest_p)->nseqs * sizeof(SQINFO));
    

        
    for (i=0; i<(*prMSeqDest_p)->nseqs; i++) {
        (*prMSeqDest_p)->seq[i] = CkStrdup(prMSeqSrc->seq[i]);
        (*prMSeqDest_p)->orig_seq[i] = CkStrdup(prMSeqSrc->orig_seq[i]);
        SeqinfoCopy(&(*prMSeqDest_p)->sqinfo[i], &prMSeqSrc->sqinfo[i]);
    }
}
/***   end: CopyMSeq   ***/



/**
 * @brief
 *
 * @param[in] seqname
 * The sequence name to search for
 * @param[in] mseq
 * The multiple sequence structure to search in
 *
 * @return -1 on failure, sequence index of matching name otherwise
 *
 * @warning If sequence name happens to be used twice, only the first
 * one will be reported back
 * 
 */    
int
FindSeqName(char *seqname, mseq_t *mseq)
{
    int i; /* aux */

    assert(NULL!=mseq);

    for (i=0; i<mseq->nseqs; i++) {
        if (STR_EQ(mseq->sqinfo[i].name, seqname)) {
            return i;
        }
    }
    
    return -1;
}
/***   end: FindSeqName()   ***/

    
/**
 * @brief Frees an mseq_t and it's members and zeros all members
 *
 * @param[in] mseq mseq_to to free
 *
 * @note use in conjunction with NewMSeq()
 * @see new_mseq
 */
void
FreeMSeq(mseq_t **mseq)
{
    int i;
    
    if (NULL==(*mseq)) {
        return;
    }
        
	if ((*mseq)->filename) {
        (*mseq)->filename = CKFREE((*mseq)->filename);
    }

    for (i=0; i<(*mseq)->nseqs; i++) {
        FreeSequence((*mseq)->seq[i], &(*mseq)->sqinfo[i]);
        CKFREE((*mseq)->orig_seq[i]);
    }
    if ((*mseq)->seq) {
        CKFREE((*mseq)->seq);
    }
    if ((*mseq)->orig_seq) {  /* FIXME (FS): only ptr to ptr freed, actual sequences NOT freed*/
        CKFREE((*mseq)->orig_seq);
    }
    if ((*mseq)->sqinfo) {
        CKFREE((*mseq)->sqinfo);
    }

    /* FIXME (FS): problem with freeing tree_order */
    if ((*mseq)->tree_order){
        CKFREE((*mseq)->tree_order);
    }

    if (NULL != (*mseq)->pppcHMMBNames){ /* FS, r291 -> */
        for (i = 0; (*mseq)->pppcHMMBNames[i] && (i < (*mseq)->nseqs); i++){
            int iIter = 0;
            for (iIter = 0; NULL != (*mseq)->pppcHMMBNames[i][iIter]; iIter++){
                CKFREE((*mseq)->pppcHMMBNames[i][iIter]);
            }
        }
    }
	(*mseq)->seqtype = SEQTYPE_UNKNOWN;
	(*mseq)->nseqs = 0;

    CKFREE((*mseq));
}
/***   end: FreeMSeq   ***/


/**
 * @brief Write alignment to file.
 *
 * @param[in] mseq
 * The mseq_t struct containing the aligned sequences
 * @param[in] pcAlnOutfile
 * The name of the output file
 * @param[in] outfmt
 * The alignment output format (defined in squid.h)
 * @param[in] iWrap
 * length of line for Clustal/Fasta format
 *
 * @return Non-zero on error
 *
 * @note We create a temporary squid MSA struct in here because we never
 * use it within clustal. We might be better of using the old clustal
 * output routines instead.
 * 
 */    
int
WriteAlignment(mseq_t *mseq, const char *pcAlnOutfile, int outfmt, int iWrap, bool bResno)
{
    int i; /* aux */
    MSA *msa; /* squid's alignment structure */
    FILE *pfOut = NULL;
    int key; /* MSA struct internal index for sequence */
    int alen; /* alignment length */
    bool use_stdout;
    
    assert(mseq!=NULL);

    if (MSAFILE_UNKNOWN == outfmt) {
        Log(&rLog, LOG_ERROR, "Unknown output format chosen");
        return -1;
    }

    if (NULL == pcAlnOutfile) {
        pfOut = stdout;
        use_stdout = TRUE;
    } else {
        use_stdout = FALSE;
        if (NULL == (pfOut = fopen(pcAlnOutfile, "w"))) {
            Log(&rLog, LOG_ERROR, "Could not open file %s for writing", pcAlnOutfile);
            return -1;
        }
    }
    

    /* derive alignment length from first seq */
    alen = strlen(mseq->seq[0]);
    
    msa  = MSAAlloc(mseq->nseqs, alen);

    /* basic structure borrowed code from squid-1.9g/a2m.c:ReadA2M()
     * we actually create a copy of mseq. keeping the pointers becomes
     * messy when calling MSAFree()
     */
    for (i=0; i<mseq->nseqs; i++) {
        char *this_name = NULL; /* mseq sequence name */
        char *this_seq = NULL; /* mseq sequence */
        SQINFO *this_sqinfo = NULL; /* mseq sequence name */
        int iI;

        /* mseq->tree_order encodes to order in which sequences are listed in the guide-tree,
           if the user wants the sequence output in the input-order then mseq->tree_order==NULL, 
           otherwise mseq->tree_order!=NULL, containing the indices of the sequences, FS, r274 ->  */
        iI = (NULL == mseq->tree_order) ? i : mseq->tree_order[i];

        this_name = mseq->sqinfo[iI].name; /* mseq sequence name */
        this_seq = mseq->seq[iI]; /* mseq sequence */
        this_sqinfo = &mseq->sqinfo[iI]; /* mseq sequence name */

        key = GKIStoreKey(msa->index, this_name);
        msa->sqname[key] = sre_strdup(this_name, strlen(this_name));

        /* setting msa->sqlen[idx] and msa->aseq[idx] */
        msa->sqlen[key] = sre_strcat(&(msa->aseq[key]), msa->sqlen[key],
                                     this_seq, strlen(this_seq));
        
        if (this_sqinfo->flags & SQINFO_DESC) {
            /* FIXME never get here ... */
            MSASetSeqDescription(msa, key, this_sqinfo->desc);
        }           
        /* FIXME extend this by copying more stuff according to flags.
         * See MSAFileRead() in msa.c and used functions there
         *
         * Problem is that we never parse MSA information as we use squid'sSeqFile
         */

        msa->nseq++;

    } /* 0 <= i < mseq->nseqs */


    /* FIXME Would like to, but can't use MSAVerifyParse(msa) here, as it
     * will die on error. Need to implement our own version
     */
#if 0
    MSAVerifyParse(msa);
#endif

    /* The below is copy of MSAFileWrite() which originally only writes to stdout.
     */

    /* Be sloppy and make a2m and fasta the same. same for vienna (which is
       the same). same same. can can. boleh boleh */
    if (outfmt==SQFILE_FASTA)
        outfmt = MSAFILE_A2M;
    if (outfmt==SQFILE_VIENNA)
        outfmt = MSAFILE_VIENNA;
    
    switch (outfmt) {
    case MSAFILE_A2M:
        /*WriteA2M(pfOut, msa, 0);*/
        WriteA2M(pfOut, msa, iWrap);
        break;
    case MSAFILE_VIENNA:
        /*WriteA2M(pfOut, msa, 1);*/
        WriteA2M(pfOut, msa, INT_MAX);
        break;
    case MSAFILE_CLUSTAL:
        WriteClustal(pfOut, msa, iWrap, TRUE==bResno ? 1 : 0, mseq->seqtype);
        break;
    case MSAFILE_MSF:
        WriteMSF(pfOut, msa);
        break;
    case MSAFILE_PHYLIP:
        WritePhylip(pfOut, msa);
        break;
    case MSAFILE_SELEX:
        WriteSELEX(pfOut, msa);
        break;
    case MSAFILE_STOCKHOLM:
        WriteStockholm(pfOut, msa);
        break;
    default:
        Log(&rLog, LOG_FATAL, "internal error: %s",
              "invalid output format should have been detected before");
    }

    if (use_stdout == FALSE) {
        (void) fclose(pfOut);
        Log(&rLog, LOG_INFO,
             "Alignment written to %s", pcAlnOutfile);
    }
    MSAFree(msa);
    
    return 0; 
}
/***   end of WriteAlignment()   ***/


/**
 * @brief Removes all gap-characters from a sequence.
 *
 * @param[out] seq
 * Sequence to dealign
 *
 * @note seq will not be reallocated
 */    
void
DealignSeq(char *seq)
{
    int aln_pos;
    int dealn_pos;

    assert(seq!=NULL);

    dealn_pos=0;
    for (aln_pos=0; aln_pos<(int)strlen(seq); aln_pos++) {
        if (! isgap(seq[aln_pos])) {
            seq[dealn_pos++] = seq[aln_pos];
        }
    }
    seq[dealn_pos] = '\0';

    return; 
}
/***   end: DealignSeq()   ***/



/**
 * @brief Sort sequences by length
 *
 * @param[out] prMSeq
 * mseq to sort by length
 * @param[out] cOrder
 * Sorting order. 'd' for descending, 'a' for ascending.
 *
 * 
 */    
void
SortMSeqByLength(mseq_t *prMSeq, const char cOrder)
{
    int *piSeqLen;
    int *piOrder;
    int iSeqIndex;
    mseq_t *prMSeqCopy = NULL;

    assert('a'==cOrder || 'd'==cOrder);

    Log(&rLog, LOG_WARN, 
        "FIXME: This modifies sequence ordering. Might not be what user wants. Will change output order as well");

    piSeqLen = (int *) CKMALLOC(prMSeq->nseqs * sizeof(int));
    piOrder = (int *) CKMALLOC(prMSeq->nseqs * sizeof(int));
    for (iSeqIndex=0; iSeqIndex<prMSeq->nseqs; iSeqIndex++) {
        piSeqLen[iSeqIndex] = prMSeq->sqinfo[iSeqIndex].len;
    }
    QSortAndTrackIndex(piOrder, piSeqLen, prMSeq->nseqs, cOrder, FALSE);
    
    CopyMSeq(&prMSeqCopy, prMSeq);
    for (iSeqIndex=0; iSeqIndex<prMSeq->nseqs; iSeqIndex++) {    
        /* copy mseq entry
         */
        CKFREE(prMSeq->seq[iSeqIndex]);
        prMSeq->seq[iSeqIndex] = CkStrdup(prMSeqCopy->seq[piOrder[iSeqIndex]]);
        
        CKFREE(prMSeq->orig_seq[iSeqIndex]);
        prMSeq->orig_seq[iSeqIndex] = CkStrdup(prMSeqCopy->orig_seq[piOrder[iSeqIndex]]);
        
        SeqinfoCopy(&prMSeq->sqinfo[iSeqIndex], &prMSeqCopy->sqinfo[piOrder[iSeqIndex]]);
    }

    CKFREE(piSeqLen);
    CKFREE(piOrder);
    FreeMSeq(&prMSeqCopy);
    
    return; 
}
/***   end: SortMSeqByLength()   ***/



/**
 * @brief Checks if sequences in given mseq structure are aligned. By
 * definition this is only true, if sequences are of the same length
 * and at least one gap was found
 *
 * @param[in] prMSeq
 * Sequences to check
 *
 * @return TRUE if sequences are aligned, FALSE if not
 *
 * 
 */    
bool
SeqsAreAligned(mseq_t *prMSeq, bool bIsProfile, bool bDealignInputSeqs)
{
    bool bGapFound, bSameLength;
    int iSeqIdx; /* sequence counter */
    int iSeqPos; /* sequence string position counter */

    /* Special case of just one sequence:
     * it is arguable that a single sequence qualifies as a profile, 
     * however, this is what we do at the first stage of MSA anyway. 
     * So, if there is only 1 sequence it is a 1-profile 
     * and it is (defined to be) aligned (with itself). FS, r240 -> 241
     */
    if (1 == prMSeq->nseqs) {
        return TRUE;
    }


    /* Check if sequences are aligned. For being aligned, the
     * sequences have to be of same length (bSameLength) and at least
     * one of them has to contain at least one gap (bGapFound)
     */
    bGapFound = FALSE;
    bSameLength = TRUE;
    for (iSeqIdx=0; (iSeqIdx < prMSeq->nseqs); iSeqIdx++) {
        if ( (FALSE == bGapFound) ){
            for (iSeqPos=0;
                 iSeqPos<prMSeq->sqinfo[iSeqIdx].len && false==bGapFound;
                 iSeqPos++) {
                if  (isgap(prMSeq->seq[iSeqIdx][iSeqPos])) {
                    bGapFound = TRUE;
                    /* skip rest of sequence */
                    break;
                }
            }
        } /* gap not (yet) found */

        if (iSeqIdx>0) {
            if (prMSeq->sqinfo[iSeqIdx].len != prMSeq->sqinfo[iSeqIdx-1].len) {
                bSameLength = FALSE;
                /* no need to continue search, bSameLength==FALSE is
                 * sufficient condition */
                break;
            }
        }
    } /* 0 <= iSeqIdx < prMSeq->nseqs */
#if 0
    Log(&rLog, LOG_FORCED_DEBUG, "bSameLength=%d bGapFound=%d", bSameLength, bGapFound);
#endif

#if 0
    if ( (TRUE == bSameLength) && ((TRUE == bGapFound) || (TRUE == bIsProfile)) ) {
        return TRUE;
    } else {
        if ((FALSE == bSameLength) && (TRUE == bGapFound) && (FALSE == bDealignInputSeqs)){
            Log(&rLog, LOG_FORCED_DEBUG, "Potential Problem: Gaps encountered but not all sequences have same length, consider using --dealign");
        }
        return FALSE;
    }   
#else
    if (FALSE == bSameLength){
        /* if sequences don't have same lengths they can never be profile */
        if (TRUE == bGapFound){
            Log(&rLog, LOG_FORCED_DEBUG, "Potential Problem: sequences (N=%d) don't have same lengths but contain gaps, consider using --dealign", prMSeq->nseqs);
        }
        return FALSE;
    }
    else { /* here all sequences have same lengths */
        if (TRUE == bGapFound){
            /* if at least one sequence contains gaps (and all have the same lengths) 
               then we can be sure it is a profile */
            return TRUE;
        }
        /* here all sequences have same lengths but no sequences contain any gaps */
        else if (TRUE == bIsProfile){
            /* if the user says it is a profile then it is */
            return TRUE;
        }
        else {
            return FALSE;
        }
    }
#endif

} 
/***   end: SeqsAreAligned()   ***/



/**
 * @brief Creates a new sequence entry and appends it to an existing mseq
 * structure.
 *
 * @param[out] prMSeqDest_p
 * Already existing and initialised mseq structure
 * @param[in] pcSeqName
 * sequence name of the sequence to add
 * @param[in] pcSeqRes
 * the actual sequence (residues) to add
 * 
 * @note Don't forget to update the align and type flag if necessary!
 *
 * FIXME allow adding of more features
 *
 */    
void
AddSeq(mseq_t **prMSeqDest_p, char *pcSeqName, char *pcSeqRes)
{
    int iSeqIdx = 0;
    SQINFO sqinfo;

    assert(NULL != prMSeqDest_p);
    assert(NULL != pcSeqName);
    assert(NULL != pcSeqRes);

    iSeqIdx = (*prMSeqDest_p)->nseqs;

    (*prMSeqDest_p)->seq =  (char **)
        CKREALLOC((*prMSeqDest_p)->seq, (iSeqIdx+1) * sizeof(char *));
    (*prMSeqDest_p)->orig_seq =  (char **)
        CKREALLOC((*prMSeqDest_p)->orig_seq, (iSeqIdx+1) * sizeof(char *));
    (*prMSeqDest_p)->sqinfo =  (SQINFO *)
        CKREALLOC((*prMSeqDest_p)->sqinfo, (iSeqIdx+1) * sizeof(SQINFO));


    (*prMSeqDest_p)->seq[iSeqIdx] = CkStrdup(pcSeqRes);
    (*prMSeqDest_p)->orig_seq[iSeqIdx] = CkStrdup(pcSeqRes);

    /* should probably get ri of SqInfo altogether in the long run and just
       transfer the intersting members into our own struct
     */
    sqinfo.flags = 0; /* init */

    sqinfo.len = strlen(pcSeqRes);
    sqinfo.flags |= SQINFO_LEN;

    /* name is an array of SQINFO_NAMELEN length */
    strncpy(sqinfo.name, pcSeqName, SQINFO_NAMELEN-1);
    sqinfo.name[SQINFO_NAMELEN-1] = '\0';
    sqinfo.flags |= SQINFO_NAME;
    
    SeqinfoCopy(&(*prMSeqDest_p)->sqinfo[iSeqIdx],
                & sqinfo);

    (*prMSeqDest_p)->nseqs++;

    return; 
}
/* end of  AddSeq() */




/**
 * @brief Appends an mseq structure to an already existing one.
 * filename will be left untouched.
 *
 * @param[in] prMSeqDest_p
 * MSeq structure to which to append to
 * @param[out] prMSeqToAdd
 * MSeq structure which is to append
 * 
 */    
void
JoinMSeqs(mseq_t **prMSeqDest_p, mseq_t *prMSeqToAdd)
{
    int iSrcSeqIndex;
    int iNewNSeq;
    
    assert(NULL != prMSeqDest_p && NULL != (*prMSeqDest_p));
    assert(NULL != prMSeqToAdd);
    
    if (0 == prMSeqToAdd->nseqs) {
        Log(&rLog, LOG_WARN, "Was asked to add 0 sequences");
        return;
    }
    
    /* warn on seqtype mismatch and keep original seqtype */
    if ((*prMSeqDest_p)->seqtype != prMSeqToAdd->seqtype) {
        Log(&rLog, LOG_WARN, "Joining sequences of different type");
    }
    
    /* leave filename as it is */

    /*
     * copy new seq/s, orig_seq/s, sqinfo/s
     */
    iNewNSeq = (*prMSeqDest_p)->nseqs + prMSeqToAdd->nseqs;
    
    (*prMSeqDest_p)->seq =  (char **)
        CKREALLOC((*prMSeqDest_p)->seq, iNewNSeq * sizeof(char *));
    
    (*prMSeqDest_p)->orig_seq =  (char **)
        CKREALLOC((*prMSeqDest_p)->orig_seq, iNewNSeq * sizeof(char *));
    
    (*prMSeqDest_p)->sqinfo =  (SQINFO *)
        CKREALLOC((*prMSeqDest_p)->sqinfo, iNewNSeq * sizeof(SQINFO));
    
    
    for (iSrcSeqIndex=0; iSrcSeqIndex < prMSeqToAdd->nseqs; iSrcSeqIndex++) {
        int iDstSeqIndex = (*prMSeqDest_p)->nseqs++;
        
        (*prMSeqDest_p)->seq[iDstSeqIndex] =
            CkStrdup(prMSeqToAdd->seq[iSrcSeqIndex]);
        
        (*prMSeqDest_p)->orig_seq[iDstSeqIndex] =
            CkStrdup(prMSeqToAdd->orig_seq[iSrcSeqIndex]);
        
        SeqinfoCopy(&(*prMSeqDest_p)->sqinfo[iDstSeqIndex],
                    & prMSeqToAdd->sqinfo[iSrcSeqIndex]);
    }

    (*prMSeqDest_p)->nseqs = iNewNSeq;

#if 0
    /* 2nd arg is bIsProfile, which when set TRUE skips 
     * the check for gaps. here always check for gaps, 
     * so set FALSE (main reason is that it is easier), FS, r282 -> */    
    /* had a problem at this stage, therefore dispense with gap check, FS, r290 -> */
    /* 3rd argument is dealignment flag, do not dealign profiles */
    (*prMSeqDest_p)->aligned = SeqsAreAligned(*prMSeqDest_p, TRUE/*FALSE*/, FALSE);
#else 
	(*prMSeqDest_p)->aligned = TRUE;
#endif
    
    return; 
}
/***   end: JoinMSeqs()   ***/
