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
 *  RCS $Id: progress.h 193 2011-02-07 15:45:21Z andreas $
 */


#ifndef CLUSTALO_PROGRESS_H
#define CLUSTALO_PROGRESS_H

#include "squid/stopwatch.h"

typedef struct {
    /* where to write to */
    FILE *prFile;
    /* prefix printed before each step */
    char *pcPrefix;
    bool bPrintCR;
    char pcLastLogMsg[1024];
    Stopwatch_t *prStopwatch;
} progress_t;


extern void
NewProgress(progress_t **pprProgress, FILE *prFile, char *pcPrefix, bool bPrintCR);

extern void
FreeProgress(progress_t **pprProgress);

extern void
ProgressLog(progress_t *prProgress, 
			unsigned long int iStep, unsigned long int iTotalSteps,
			bool bForceOutput);

extern void
ProgressDone(progress_t *pprProgress);

#endif
