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
 *  RCS $Id: clustal-omega.h 212 2011-03-10 15:09:46Z andreas $
 */

#ifndef CLUSTALO_H
#define CLUSTALO_H



#ifdef HAVE_OPENMP
#include <omp.h>
#endif
#include <stdbool.h>

#include "clustal-omega-config.h"

/* the following needs to be kept in sync with library_include_HEADERS of all
 * subdir Makefile.am's 
 */

/* hhalign */
#include "hhalign/general.h"
#include "hhalign/hhfunc.h"


/* clustal */
#include "clustal/log.h"
#include "clustal/util.h"
#include "clustal/symmatrix.h"
#include "clustal/tree.h"
#include "clustal/seq.h"
#include "clustal/mbed.h"
#include "clustal/weights.h"
#include "clustal/pair_dist.h"
#include "clustal/hhalign_wrapper.h"



#define CLUSTERING_UNKNOWN 0
#define CLUSTERING_UPGMA 1

/* weights will be computed if 1. but are not really used for now and they
 * might slow things down. also, mbed's screws up branch lengths which will
 * have a negative effect on weights 
*/
#define USE_WEIGHTS 0

extern int iNumberOfThreads;

/* output order */
enum {INPUT_ORDER = 0, TREE_ORDER};


/** user/commandline options
 *
 * changes here will have to be reflected in ParseCommandLine()
 * and during setup of the default opts
 *
 */
typedef struct {
    /* auto: Clustal (know what) is good for you
     */
    bool bAutoOptions;

    /* Distance matrix
     */
    /** distance matrix input file */
    char *pcDistmatInfile;
    /** distance matrix output file */
    char *pcDistmatOutfile;
    
    /* Clustering / guide-tree
     */
    /** clustering type (from cmdline arg) */
    int iClusteringType;
    /** number of sequences in cluster */
    int iClustersizes; 
    /** file with clustering information */
    char *pcClustfile; 
    /** use transitivity */
    int iTransitivity; 
    /** file with posterior probability information */
    char *pcPosteriorfile; 
    /** pairwise distance method */
    int iPairDistType;
    /** use mbed-like clustering */
    bool bUseMbed;
    /** use mbed-like clustering also during iteration */
    bool bUseMbedForIteration;
    /** pile-up flag */
    bool bPileup;
    /** guidetree output file */
    char *pcGuidetreeOutfile;
    /** guidetree input file */
    char *pcGuidetreeInfile;
    /** use Kimura corrected distance */
    bool bUseKimura;
    /** print percentage identity */
    bool bPercID;

    /* HMMs
     */
    /** HMM input files. index range: 0..iHMMInputFiles */
    char **ppcHMMInput;
    /** number of provided HMM input files. not really a user
       option but need for ppcHMMInput */
    int iHMMInputFiles;
    /** HMM batch-file, specify HMMs for individual sequences. FS, r291 -> */
    char *pcHMMBatch;

    /* Iteration
     */
    /** number of iterations */
    int iNumIterations;
    /** determine number of iterations automatically */
	bool bIterationsAuto;
    /** maximum number of hmm iterations */
    int iMaxHMMIterations;
    /** max number of guidetree iterations */
	int iMaxGuidetreeIterations;
    
    hhalign_para rHhalignPara;

    /* changes here will have to be reflected in FreeAlnOpts(),
	 * SetDefaultAlnOpts(), AlnOptsLogicCheck() etc 
	 */
} opts_t;



extern void 
PrintLongVersion(char *pcStr, int iSize);

extern void
SetDefaultAlnOpts(opts_t *opts);

extern void
FreeAlnOpts(opts_t *aln_opts);

extern void
AlnOptsLogicCheck(opts_t *opts);

extern void
PrintAlnOpts(FILE *prFile, opts_t *opts);

extern void
InitClustalOmega(int iNumThreadsToUse);

extern void
SequentialAlignmentOrder(int **piOrderLR_p, int iNumSeq);

extern int
AlignmentOrder(int **piOrderLR_p, double **pdSeqWeights_p, mseq_t *prMSeq,
               int iPairDistType, char *pcDistmatInfile, char *pcDistmatOutfile,
               int iClusteringType, int iClustersizes,
               char *pcGuidetreeInfile, char *pcGuidetreeOutfile, char *pcClusterFile,
               bool bUseMBed, bool bPercID);

extern int
Align(mseq_t *prMSeq, 
      mseq_t *prMSeqProfile,
	  opts_t *prOpts);

extern int
AlignProfiles(mseq_t *prMSeqProfile1, 
			  mseq_t *prMSeqProfile2, hhalign_para rHhalignPara);

#endif
extern int
ReadPseudoCountParams(hhalign_para *rHhalignPara_p, char *pcPseudoFile);
