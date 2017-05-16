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
 *  RCS $Id: ktuple_pair.c 305 2016-06-13 13:46:02Z fabian $
 *
 *
 * K-Tuple code for pairwise alignment (Wilbur and Lipman, 1983; PMID
 * 6572363). Most code taken from showpair.c (Clustal 1.83)
 * DD: some functions now have lots of parameters as static variables
 * were removed to make code OpenMP-friendly
 *
 */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif


#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>

#ifdef HAVE_OPENMP
#include <omp.h>
#endif

#include "squid/squid.h"
#include "util.h"
#include "symmatrix.h"
#include "ktuple_pair.h"
#include "log.h"
#include "progress.h"

#define END_MARK -3 /* see interface.c in 1.83 */
#define NUMRES 32 /* max size of comparison matrix */

/* see notes below */
#undef SORT_LAST_ELEMENT_AS_WELL

/* gap_pos1 = NUMRES-2; /@ code for gaps inserted by clustalw @/ */
static const int GAP_POS2 = NUMRES-1; /* code for gaps already in alignment */
static bool DNAFLAG = FALSE;

static const char *AMINO_ACID_CODES = "ABCDEFGHIKLMNPQRSTUVWXYZ-";
static const char *NUCLEIC_ACID_CODES = "ACGTUN-";
/* As far as I understand the gap symbol should not be necessary here,
 * because we use isgap for testing later anyway. But changing this,
 * will affect max_res_code and max_nuc as well. So I leave it for now
 * as it is. AW
 */

static bool percent = TRUE;

static void make_ptrs(int *tptr, int *pl, const int naseq, const int l, const int ktup, const int max_res_code, char **seq_array);
static void put_frag(const int fs, const int v1, const int v2, const int flen, const int curr_frag, int *next, int *maxsf, int **accum);
static bool frag_rel_pos(int a1, int b1, int a2, int b2, int ktup);
static void des_quick_sort(int *array1, int *array2, const int array_size);
static void pair_align(int seq_no, int l1, int l2, int max_res_code, ktuple_param_t *aln_param,
    char **seq_array, int *maxsf, int **accum, int max_aln_length,
    int *zza, int *zzb, int *zzc, int *zzd);
static void encode(char *seq, char *naseq, int l, const char *res_codes);
static int res_index(const char *lookup, char c);


typedef struct {
    int i1;
    int i2;
} two_ints_t;



/* default ktuple pairwise alignment parameters
 *
 */
/* protein
 */
/* designated initializer */
const ktuple_param_t default_protein_param = {
     .ktup = 1,
     .wind_gap = 3,
     .signif = 5,
     .window = 5,
};
/* dna
 */
/* designated initializer */
const ktuple_param_t default_dna_param = {
    .ktup = 2,
    .wind_gap = 5,
    .signif = 4,
    .window = 4,
};


/**
 * note: naseq should be unit-offset
 */
static void
encode(char *seq, char *naseq, int l, const char *res_codes)
{
    /* code seq as ints .. use GAP_POS2 for gap */
    register int i;
    bool seq_contains_unknown_char = FALSE;
    /*LOG_DEBUG("seq=%s naseq=%p l=%d", &(seq[1]), naseq, l); */


    for (i=1; i<=l; i++) {
        char res = toupper(seq[i]);
        if (isgap(res)) {
            naseq[i] = GAP_POS2; /* gap in input */
        } else {
            naseq[i] = res_index(res_codes, res);
        }

        /*LOG_DEBUG("Character '%c' at pos %d", res, i);*/
        if (-1 == naseq[i]) {
            seq_contains_unknown_char = TRUE;
            /*LOG_DEBUG("Unknown character '%c' at pos %d", res, i);*/
        }
        /*LOG_DEBUG("na_seq[%d]=%d", i, naseq[i]);*/
    }

    if (TRUE == seq_contains_unknown_char)
        Log(&rLog, LOG_WARN, "Unknown character in seq '%s'", &(seq[1]));

    naseq[i] = END_MARK;

    return;
}
/* end of encode */


/**
 *
 */
static int
res_index(const char *t, char c)
{
    register int i;
    for (i=0; t[i] && t[i] != c; i++)
        ;
    if (t[i]) {
        return (i);
    } else {
        return -1;
    }
}
/* end of res_index */


/**
 *
 */
static void
make_ptrs(int *tptr, int *pl, const int naseq, const int l, const int ktup, const int max_res_code, char **seq_array)
{
    /* FIXME make 10 a constant and give it a nice name */
    static int a[10];
    int i, j, code, flag;
    char residue;
    const int limit = (int) pow((double)(max_res_code+1),(double)ktup);

    for (i=1;i<=ktup;i++)
        a[i] = (int) pow((double)(max_res_code+1),(double)(i-1));

    for (i=1; i<=limit; ++i)
        pl[i]=0;
    for (i=1; i<=l; ++i)
        tptr[i]=0;

    for (i=1; i<=(l-ktup+1); ++i) {
        code=0;
        flag=FALSE;
        for (j=1; j<=ktup; ++j) {
            /* Log(&rLog, LOG_FORCED_DEBUG, "naseq=%d i=%d j=%d seq_array[naseq]=%p",
             * naseq, i, j, seq_array[naseq]);
             */
            residue = seq_array[naseq][i+j-1];
            /* Log(&rLog, LOG_FORCED_DEBUG, "residue = %d", residue); */
            if ((residue<0) || (residue > max_res_code)){
                flag=TRUE;
                break;
            }
            code += ((residue) * a[j]);
        }
        if (flag)
            continue;
        ++code;
        if (0 != pl[code])
            tptr[i] =pl[code];
        pl[code] = i;
    }

    return;
}
/* end of make_ptrs */


/**
 *
 * FIXME Why hardcoding of 5?
 */
static void
put_frag(const int fs, const int v1, const int v2, const int flen, const int curr_frag, int *next, int *maxsf, int **accum)
{
    int end;
    accum[0][curr_frag]=fs;
    accum[1][curr_frag]=v1;
    accum[2][curr_frag]=v2;
    accum[3][curr_frag]=flen;

    if (!*maxsf) {
        *maxsf=1;
        accum[4][curr_frag]=0;
        return;
    }

    if (fs >= accum[0][*maxsf]) {
        accum[4][curr_frag]=*maxsf;
        *maxsf=curr_frag;
        return;
    } else {
        *next=*maxsf;
        while (TRUE) {
            end=*next;
            *next=accum[4][*next];
            if (fs>=accum[0][*next])
                break;
        }
        accum[4][curr_frag]=*next;
        accum[4][end]=curr_frag;
    }

    return;
}
/* end of put_frag */


/**
 *
 */
static bool
frag_rel_pos(int a1, int b1, int a2, int b2, int ktup)
{
    if (a1-b1 == a2-b2) {
        if (a2<a1) {
            return TRUE;
        }
    } else {
        if (a2+ktup-1<a1 && b2+ktup-1<b1) {
            return TRUE;
        }
    }
    return FALSE;
}
/* end of frag_rel_pos */




/**
 *
 * @note: This is together with des_quick_sort most time consuming
 * routine according to gprof on r110. Tried to replace it with qsort
 * and/or QSortAndTrackIndex(), which is always slower! So we keep the
 * original.
 *
 * Original doc: Quicksort routine, adapted from chapter 4, page 115
 * of software tools by Kernighan and Plauger, (1986). Sort the
 * elements of array1 and sort the elements of array2 accordingly
 *
 * There might be a bug here. The original function apparently never
 * touches the last element and keeps it as is. Tried to fix this (see
 * SORT_LAST_ELEMENT_AS_WELL) which gives slightly worse performance
 * (-0.5% on BB). My fix might not be working or it's not a bug at
 * all...
 *
 *
 *
 */
static void
des_quick_sort(int *array1, int *array2, const int array_size)
{
    int temp1, temp2;
    int p, pivlin;
    int i, j;
    int lst[50], ust[50];       /* the maximum no. of elements must be*/
                                /* < log(base2) of 50 */

#if 0
    for (i=1; i<=array_size; i++) {
        Log(&rLog, LOG_FORCED_DEBUG, "b4 sort array1[%d]=%d array2[%d]=%d", i, array1[i], i, array2[i]);
    }
#endif
    lst[1] = 1;

#ifdef SORT_LAST_ELEMENT_AS_WELL
    ust[1] = array_size;
#else
    /* original */
    ust[1] = array_size-1;
#endif
    p = 1;


    while (p > 0) {
        if (lst[p] >= ust[p]) {
            p--;
        } else {
            i = lst[p] - 1;
            j = ust[p];
            pivlin = array1[j];
            while (i < j) {
                for (i=i+1; array1[i] < pivlin; i++)
                    ;
                for (j=j-1; j > i; j--)
                    if (array1[j] <= pivlin) break;
                if (i < j) {
                    temp1     = array1[i];
                    array1[i] = array1[j];
                    array1[j] = temp1;

                    temp2     = array2[i];
                    array2[i] = array2[j];
                    array2[j] = temp2;
                }
            }

            j = ust[p];

            temp1     = array1[i];
            array1[i] = array1[j];
            array1[j] = temp1;

            temp2     = array2[i];
            array2[i] = array2[j];
            array2[j] = temp2;

            if (i-lst[p] < ust[p] - i) {
                lst[p+1] = lst[p];
                ust[p+1] = i - 1;
                lst[p]   = i + 1;

            } else {
                lst[p+1] = i + 1;
                ust[p+1] = ust[p];
                ust[p]   = i - 1;
            }
            p = p + 1;
        }
    }

#if 0
    for (i=1; i<=array_size; i++) {
        Log(&rLog, LOG_FORCED_DEBUG, "after sort array1[%d]=%d array2[%d]=%d", i, array1[i], i, array2[i]);
    }
#endif

    return;
}
/* end of des_quick_sort */



/**
 *
 * FIXME together with des_quick_sort most time consuming routine
 * according to gprof on r110
 *
 */
static void
pair_align(int seq_no, int l1, int l2, int max_res_code, ktuple_param_t *aln_param,
    char **seq_array, int *maxsf, int **accum, int max_aln_length,
    int *zza, int *zzb, int *zzc, int *zzd)
{
    int next; /* forrmerly static */
    int pot[8],i, j, l, m, flag, limit, pos, vn1, vn2, flen, osptr, fs;
    int tv1, tv2, encrypt, subt1, subt2, rmndr;
    char residue;
    int *diag_index;
    int *displ;
    char *slopes;
    int curr_frag;
    const int tl1 = (l1+l2)-1;

    assert(NULL!=aln_param);

    /*
      Log(&rLog, LOG_FORCED_DEBUG, "DNAFLAG=%d seq_no=%d l1=%d l2=%d window=%d ktup=%d signif=%d wind_gap=%d",
      DNAFLAG, seq_no, l1, l2, window, ktup, signif,
      wind_gap);
    */

    slopes = (char *) CKCALLOC(tl1+1, sizeof(char));
    displ = (int *) CKCALLOC(tl1+1, sizeof(int));
    diag_index = (int *) CKMALLOC((tl1+1) * sizeof(int));

    for (i=1; i<=tl1; ++i) {
        /* unnecessary, because we calloced: slopes[i] = displ[i] = 0; */
        diag_index[i] = i;
    }

    for (i=1;i<=aln_param->ktup;i++)
        pot[i] = (int) pow((double)(max_res_code+1),(double)(i-1));
    limit = (int) pow((double)(max_res_code+1),(double)aln_param->ktup);



    /* increment diagonal score for each k_tuple match */

    for (i=1; i<=limit; ++i) {
        vn1=zzc[i];
        while (TRUE) {
            if (!vn1) break;
            vn2 = zzd[i];
            while (0 != vn2) {
                osptr = vn1-vn2+l2;
                ++displ[osptr];
                vn2 = zzb[vn2];
            }
            vn1=zza[vn1];
        }
    }


    /* choose the top SIGNIF diagonals
     */

#ifdef QSORT_REPLACEMENT
    /* This was an attempt to replace des_quick_sort with qsort(),
     * which turns out to be much slower, so don't use this
     */

    /* FIXME: if we use this branch, we don't need to init diag_index
     * before, because that is done in QSortAndTrackIndex()
     * automatically.
     */
#if 0
    for (i=1; i<=tl1; i++) {
        Log(&rLog, LOG_FORCED_DEBUG, "b4 sort disp[%d]=%d diag_index[%d]=%d", i, diag_index[i], i, displ[i]);
    }
#endif

    QSortAndTrackIndex(&(diag_index[1]), &(displ[1]), tl1, 'a', TRUE);

#if 0
    for (i=1; i<=tl1; i++) {
        Log(&rLog, LOG_FORCED_DEBUG, "after sort disp[%d]=%d diag_index[%d]=%d", i, diag_index[i], i, displ[i]);
    }
#endif

#else

    des_quick_sort(displ, diag_index, tl1);

#endif

    j = tl1 - aln_param->signif + 1;

    if (j < 1) {
        j = 1;
    }

    /* flag all diagonals within WINDOW of a top diagonal */

    for (i=tl1; i>=j; i--)  {
        if (displ[i] > 0) {
            pos = diag_index[i];
            l = (1   > pos - aln_param->window) ?
                1 :  pos - aln_param->window;
            m = (tl1 < pos + aln_param->window) ?
                tl1 : pos + aln_param->window;
            for (; l <= m; l++)
                slopes[l] = 1;
        }
    }

    for (i=1; i<=tl1; i++) {
        displ[i] = 0;
    }

    curr_frag=*maxsf=0;

    for (i=1; i<=(l1-aln_param->ktup+1); ++i) {
        encrypt=flag=0;
        for (j=1; j<=aln_param->ktup; ++j) {
            residue = seq_array[seq_no][i+j-1];
            if ((residue<0) || (residue>max_res_code)) {
                flag=TRUE;
                break;
            }
            encrypt += ((residue)*pot[j]);
        }
        if (flag) {
            continue;
        }
        ++encrypt;

        vn2=zzd[encrypt];

        flag=FALSE;
        while (TRUE) {
            if (!vn2) {
                flag=TRUE;
                break;
            }
            osptr=i-vn2+l2;
            if (1 != slopes[osptr]) {
                vn2=zzb[vn2];
                continue;
            }
            flen=0;
            fs=aln_param->ktup;
            next=*maxsf;

            /*
             * A-loop
             */

            while (TRUE) {
                if (!next) {
                    ++curr_frag;
                    if (curr_frag >= 2*max_aln_length) {
                        Log(&rLog, LOG_VERBOSE, "(Partial alignment)");
                        goto free_and_exit; /* Yesss! Always wanted to
                                             * use a goto (AW) */
                    }
                    displ[osptr]=curr_frag;
                    put_frag(fs, i, vn2, flen, curr_frag, &next, maxsf, accum);

                } else {
                    tv1=accum[1][next];
                    tv2=accum[2][next];

                    if (frag_rel_pos(i, vn2, tv1, tv2, aln_param->ktup)) {
                        if (i-vn2 == accum[1][next]-accum[2][next]) {
                            if (i > accum[1][next]+(aln_param->ktup-1)) {
                                fs = accum[0][next]+aln_param->ktup;
                            } else {
                                rmndr = i-accum[1][next];
                                fs = accum[0][next]+rmndr;
                            }
                            flen=next;
                            next=0;
                            continue;

                        } else {
                            if (0 == displ[osptr]) {
                                subt1=aln_param->ktup;
                            } else {
                                if (i > accum[1][displ[osptr]]+(aln_param->ktup-1)) {
                                    subt1=accum[0][displ[osptr]]+aln_param->ktup;
                                } else {
                                    rmndr=i-accum[1][displ[osptr]];
                                    subt1=accum[0][displ[osptr]]+rmndr;
                                }
                            }
                            subt2=accum[0][next] - aln_param->wind_gap + aln_param->ktup;
                            if (subt2>subt1) {
                                flen=next;
                                fs=subt2;
                            } else {
                                flen=displ[osptr];
                                fs=subt1;
                            }
                            next=0;
                            continue;
                        }
                    } else {
                        next=accum[4][next];
                        continue;
                    }
                }
                break;
            }
            /*
             * End of Aloop
             */

            vn2=zzb[vn2];
        }
    }

free_and_exit:
    CKFREE(displ);
    CKFREE(slopes);
    CKFREE(diag_index);

    return;
}
/* end of pair_align */



/**
 *
 * Will compute ktuple scores and store in tmat
 * Following values will be set: tmat[i][j], where
 * istart <= i <iend
 * and
 * jstart <= j < jend
 * i.e. zero-offset
 * tmat data members have to be preallocated
 *
 * if ktuple_param_t *aln_param == NULL defaults will be used
 */
void
KTuplePairDist(symmatrix_t *tmat, mseq_t *mseq,
               int istart, int iend,
               int jstart, int jend,
               ktuple_param_t *param_override,
               progress_t *prProgress, 
               unsigned long int *ulStepNo, unsigned long int ulTotalStepNo)
{
    /* this first group of variables were previously static
       and hence un-parallelisable */
    char **seq_array;
    int maxsf;
    int **accum;
    int max_aln_length;
    /* divide score with length of smallest sequence */
    int *zza, *zzb, *zzc, *zzd;
    int private_step_no = 0;

    int i, j, dsr;
    double calc_score;
    int max_res_code = -1;

    int max_seq_len;
    int *seqlen_array;
    /* progress_t *prProgress; */
    /* int uStepNo, uTotalStepNo; */
    ktuple_param_t aln_param = default_protein_param;
    bool bPrintCR = (rLog.iLogLevelEnabled<=LOG_VERBOSE) ? FALSE : TRUE;
    

    if(prProgress == NULL) {
        NewProgress(&prProgress, LogGetFP(&rLog, LOG_INFO), 
                    "Ktuple-distance calculation progress", bPrintCR);
    }

    /* conversion to old style data types follows
     *
     */

    seqlen_array = (int*) CKMALLOC((mseq->nseqs+1) * sizeof(int));
    for (i=0; i<mseq->nseqs; i++) {
        seqlen_array[i+1] = mseq->sqinfo[i].len;
    }

    /* setup alignment parameters
     */
    if (SEQTYPE_PROTEIN == mseq->seqtype) {
        DNAFLAG = FALSE;
        max_res_code = strlen(AMINO_ACID_CODES)-2;
        aln_param = default_protein_param;

    } else if (SEQTYPE_RNA == mseq->seqtype || SEQTYPE_DNA == mseq->seqtype) {
        DNAFLAG = TRUE;
        max_res_code = strlen(NUCLEIC_ACID_CODES)-2;
        aln_param = default_dna_param;

    } else {
        Log(&rLog, LOG_FATAL, "Internal error in %s: Unknown sequence type.", __FUNCTION__);
    }

    if (NULL!=param_override) {
        aln_param.ktup = param_override->ktup;
        aln_param.wind_gap = param_override->wind_gap;
        aln_param.signif = param_override->signif;
        aln_param.window = param_override->window;
    }

    /*LOG_DEBUG("DNAFLAG = %d max_res_code = %d", DNAFLAG, max_res_code);*/

    /* convert mseq to clustal's old-style int encoded sequences (unit-offset)
     */
    max_aln_length = 0;
    max_seq_len = 0;
    seq_array =  (char **) CKMALLOC((mseq->nseqs+1) * sizeof(char *));
    seq_array[0] = NULL;
    /* FIXME check that non of the seqs is smaller than ktup (?).
     * Otherwise segfault occurs
     */
    for (i=0; i<mseq->nseqs; i++) {
        seq_array[i+1] = (char *) CKMALLOC((seqlen_array[i+1]+2) * sizeof (char));;
    }
    for (i=0; i<mseq->nseqs; i++) {
        /*LOG_DEBUG("calling encode with seq_array[%d+1] len=%d and seq=%s",
          i, seqlen_array[i+1], mseq->seq[i]);*/
        if (TRUE == DNAFLAG) {
            encode(&(mseq->seq[i][-1]), seq_array[i+1],
                   seqlen_array[i+1], NUCLEIC_ACID_CODES);
        } else  {
            encode(&(mseq->seq[i][-1]), seq_array[i+1],
                   seqlen_array[i+1], AMINO_ACID_CODES);
        }

        if (seqlen_array[i+1]>max_seq_len) {
            max_seq_len = seqlen_array[i+1];
        }
    }
    max_aln_length = max_seq_len * 2;
    /* see sequence.c in old source */

    /* FIXME: short sequences can cause seg-fault 
     * because max_aln_length can get shorter 
     * than (max_res_code+1)^k 
     * FS, r222->r223 */
    max_aln_length = max_aln_length > pow((max_res_code+1), aln_param.ktup)+1 ? 
        max_aln_length : pow((max_res_code+1), aln_param.ktup)+1;

    /*
     *
     * conversion to old style clustal done (in no time) */


    accum = (int **) CKCALLOC(5, sizeof (int *));
    for (i=0;i<5;i++) {
        accum[i] = (int *) CKCALLOC((2*max_aln_length+1), sizeof(int));
    }
    zza = (int *) CKCALLOC( (max_aln_length+1), sizeof(int));
    zzb = (int *) CKCALLOC( (max_aln_length+1), sizeof(int));
    zzc = (int *) CKCALLOC( (max_aln_length+1), sizeof(int));
    zzd = (int *) CKCALLOC( (max_aln_length+1), sizeof(int));

    /* estimation of total number of steps (if istart and jstart are
     * both 0) (now handled in the calling routine)
     */
    /* uTotalStepNo = iend*jend - iend*iend/2 + iend/2;
    uStepNo = 0; */
    /*LOG_DEBUG("istart=%d iend=%d jstart=%d jend=%d", istart, iend, jstart, jend);*/

    for (i=istart+1; i<=iend; ++i) {
        /* by definition a sequence compared to itself should give
           a score of 0. AW */
        SymMatrixSetValue(tmat, i-1, i-1, 0.0);
        make_ptrs(zza, zzc, i, seqlen_array[i], aln_param.ktup, max_res_code, seq_array);

#ifdef HAVE_OPENMP
        #pragma omp critical(ktuple)
#endif
        {
            ProgressLog(prProgress, *ulStepNo, ulTotalStepNo, FALSE);
        }

        for (j=MAX(i+1, jstart+1); j<=jend; ++j) {
            (*ulStepNo)++;
            private_step_no++;
            /*LOG_DEBUG("comparing pair %d:%d", i, j);*/

            make_ptrs(zzb, zzd, j, seqlen_array[j], aln_param.ktup, max_res_code, seq_array);
            pair_align(i, seqlen_array[i], seqlen_array[j], max_res_code, &aln_param,
                seq_array, &maxsf, accum, max_aln_length, zza, zzb, zzc, zzd);

            if (!maxsf) {
                calc_score=0.0;
            } else {
                calc_score=(double)accum[0][maxsf];
                if (percent) {
                    dsr=(seqlen_array[i]<seqlen_array[j]) ?
                        seqlen_array[i] : seqlen_array[j];
                    calc_score = (calc_score/(double)dsr) * 100.0;
                }
            }

            /* printf("%d %d %d\n", i-1, j-1, (100.0 - calc_score)/100.0); */
            SymMatrixSetValue(tmat, i-1, j-1, (100.0 - calc_score)/100.0);

            /* the function allows you not to compute the full matrix.
             * here we explicitely make the resulting matrix a
             * rectangle, i.e. we always set full rows. in other
             * words, if we don't complete the full matrix then we
             * don't have a full symmetry. so only use the defined
             * symmetric part. AW
             */
            /*LOG_DEBUG("setting %d : %d = %f", j, i, tmat[i][j]);*/
            /* not needed anymore since we use symmatrix_t
               if (j<=iend) {
               tmat[j][i] = tmat[i][j];
               }
            */
#ifdef HAVE_OPENMP
            #pragma omp critical(ktuple)
#endif
            {
                Log(&rLog, LOG_DEBUG, "K-tuple distance for sequence pair %d:%d = %lg",
                    i, j, SymMatrixGetValue(tmat, i-1, j-1));
            }
        }
    }
    /*
      Log(&rLog, LOG_FORCED_DEBUG, "uTotalStepNo=%d for istart=%d iend=%d jstart=%d jend=%d", uStepNo, istart, iend, jstart, jend);
      Log(&rLog, LOG_FORCED_DEBUG, "Fabian = %d", iend*jend - iend*iend/2 + iend/2);
    */

/*    printf("\n\n%d\t%d\t%d\t%d\n\n", omp_get_thread_num(), uStepNo, istart, iend); */

    for (i=0;i<5;i++) {
        CKFREE(accum[i]);
    }
    CKFREE(accum);

#ifdef HAVE_OPENMP
    #pragma omp critical(ktuple)
#if 0
    {
        int tid;
        tid = omp_get_thread_num();
        printf("%s:%d: tid %d: steps %d\n", __FILE__, __LINE__, tid, private_step_no);
    }
#endif
#endif

    CKFREE(zza);
    CKFREE(zzb);
    CKFREE(zzc);
    CKFREE(zzd);

    free(seqlen_array);

    for (i=1; i<=mseq->nseqs; i++) {
        CKFREE(seq_array[i]);
    }
    CKFREE(seq_array);
}
/* end of KTuplePairDist */
