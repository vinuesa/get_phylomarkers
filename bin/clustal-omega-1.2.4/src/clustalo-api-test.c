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
 *  RCS $Id: clustalo-api-test.c 280 2013-05-16 16:04:12Z fabian $
 */


#include <stdio.h>

/* see clustal-omega.c for documentation */

/* Include clustal-omega's header. That's all you need 
 *
 * If you developing in C++, use the following instead:
 * extern "C" {
 * #include "clustal-omega.h"
 * }
 */
#include "clustal-omega.h"


int
main(int argc, char **argv)
{
    /* the multiple sequence structure */
    mseq_t *prMSeq = NULL;
    /* for openmp: number of threads to use */
    int iThreads = 1;
    /* alignment options to use */
    opts_t rAlnOpts;
    /* an input file */
    char *pcSeqInfile;
    int iAux;


    /* Must happen first: setup logger */
    LogDefaultSetup(&rLog);

    SetDefaultAlnOpts(&rAlnOpts);

    InitClustalOmega(iThreads);

    /* Get sequence input file name from command line
     */
    if (argc!=2) {
        Log(&rLog, LOG_FATAL, "Need sequence file as argument");
    }
    pcSeqInfile = argv[1];

    /* Read sequence file
     */
    NewMSeq(&prMSeq);
    if (ReadSequences(prMSeq, pcSeqInfile,
                      SEQTYPE_UNKNOWN,
                      INT_MAX, INT_MAX)) {
        Log(&rLog, LOG_FATAL, "Reading sequence file '%s' failed", pcSeqInfile);
    }

    /* Dump some info about the sequences
     */
    for (iAux=0; iAux<prMSeq->nseqs; iAux++) {
        Log(&rLog, LOG_INFO, 
             "Sequence no %d has the following name: %s",
             iAux, prMSeq->sqinfo[iAux].name);
        Log(&rLog, LOG_INFO, 
             "Sequence no %d has the following residues: %s",
             iAux, prMSeq->seq[iAux]);
        /* more info can be found in prMSeq->sqinfo[iAux] */
    }


    /* Align the sequences without a profile (NULL)
     */
    if (Align(prMSeq, NULL, &rAlnOpts)) {
        Log(&rLog, LOG_FATAL, "A fatal error happended during the alignment process");
    }


    /* Output of final alignment to stdout (NULL) as aligned fasta/a2m
     */
#define LINE_WRAP 60
    if (WriteAlignment(prMSeq, NULL, MSAFILE_A2M, LINE_WRAP)) {
        Log(&rLog, LOG_FATAL, "Could not save alignment");
    } 

    FreeMSeq(&prMSeq);

    Log(&rLog, LOG_INFO, "Successfull program exit");

    return EXIT_SUCCESS;
}
/***   end of main()   ***/
