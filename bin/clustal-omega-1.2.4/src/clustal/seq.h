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
 *  RCS $Id: seq.h 296 2014-10-07 12:15:41Z fabian $
 */

#ifndef CLUSTALO_SEQ_H
#define CLUSTALO_SEQ_H

#include "squid/squid.h"

#include "util.h"


/**
 * int-encoded sequence types.
 * these are in sync with squid's seqtypes and only used for
 * convenience here
 */
#define SEQTYPE_UNKNOWN kOtherSeq
#define SEQTYPE_DNA kDNA
#define SEQTYPE_RNA kRNA
#define SEQTYPE_PROTEIN kAmino

/* Alphabets are defined in squid.h: AMINO_ALPHABET, DNA_ALPHABET,
 * RNA_ALPHABET (all uppercase)
 */
#define AMINOACID_ANY 'X'
#define NUCLEOTIDE_ANY 'N'

/**
 * @brief structure for storing multiple sequences
 *
 */
typedef struct {
    int nseqs; /**< number of sequences */
    int seqtype; /**< sequence type */
    char *filename; /**< input file / source of sequences */
    bool aligned; /**< true if all seqs are same length **/
    
    /** (working) sequence residues as char pointer.
     * range for first index: 0--nseq-1.
     * changes during alignment.
     */
    char **seq;
    
    /** original sequence residues as char pointer.
     * range for first index: 0--nseq-1.
     * only set during input
     */
    char **orig_seq;
    
  /** order in which sequences appear in guide-tree
   */
  int *tree_order;

    /**
     * @brief Squid's sequence info structure.
     * Index range: 0--nseq-1.
     *
     * extra data are available:
     * int flags;
     *
     * name:
     * char name[SQINFO_NAMELEN];
     *
     * database identifier:
     * char id[SQINFO_NAMELEN];
     *
     * database accession no:
     * char acc[SQINFO_NAMELEN];
     *
     * description:      
     * char desc[SQINFO_DESCLEN];
     *
     * length of this seq, incl gaps in our case!:
     * int len;
     *
     * start position on source seq (valid range: 1..len):
     * int start;
     *
     * end position on source seq (valid range: 1..len):
     * int stop;
     *
     * original length of source seq:
     * int olen;
     *
     * kRNA, kDNA, kAmino, or kOther:
     * int type;
     *
     * secondary structure string (index range: 0..len-1):
     * char *ss;
     *
     * percent side chain surface access (index range: 0..len-1):
     * char *sa;                  
     * 
     * @see squid.h
     * @see LogSqInfo()
     *
     */
    SQINFO *sqinfo; 

  /* HMM batch information */
  char ***pppcHMMBNames;
  int **ppiHMMBindex;
} mseq_t;

extern void
AliStat(mseq_t *prMSeq, bool bSampling, bool bReportAll);

extern void
AddSeq(mseq_t **prMSeqDest_p, char *pcSeqName, char *pcSeqRes);

extern void
SeqSwap(mseq_t *mseq, int i, int j);

extern void
DealignMSeq(mseq_t *mseq);

extern const char *
SeqTypeToStr(int seqtype);

extern int
ReadSequences(mseq_t *prMSeq_p, char *pcSeqFile, 
              int iSeqType,  int iSeqFmt, bool bIsProfile, bool bDealignInputSeqs, 
              int iMaxNumSeq, int iMaxSeqLen, char *pcHMMBatch);

extern void
NewMSeq(mseq_t **mseq);

extern void
FreeMSeq(mseq_t **mseq);

extern void
CopyMSeq(mseq_t **prMSeqDest_p, mseq_t *prMSeqSrc);

extern void
LogSqInfo(SQINFO *sqinfo);

extern int
FindSeqName(char *seqname, mseq_t *mseq);

extern int
WriteAlignment(mseq_t *mseq, const char *aln_outfile, int msafile_format, int iWrap, bool bResno);

extern void
DealignSeq(char *seq);

extern void
ShuffleMSeq(mseq_t *prMSeq);

extern void
SortMSeqByLength(mseq_t *prMSeq, const char cOrder);

void
JoinMSeqs(mseq_t **prMSeqDest_p, mseq_t *prMSeqToAdd);

bool
SeqsAreAligned(mseq_t *prMSeq, bool bIsProfile, bool bDealignInputSeqs);

#endif
