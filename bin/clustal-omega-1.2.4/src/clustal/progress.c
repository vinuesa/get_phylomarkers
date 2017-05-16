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
 *  RCS $Id: progress.c 230 2011-04-09 15:37:50Z andreas $
 */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <stdlib.h>
#include <stdio.h>
#include <assert.h>

#include "util.h"
#include "log.h"
#include "progress.h"

#define LOGLEVEL_THRESHOLD LOG_INFO

/**
 * @brief Allocates a new progress structure and initialises its members.
 * Free with FreeProgress()
 *
 * @note Starts the internal stopwatch immediatly!
 *
 * @see FreeProgress()
 *
 * @param[out] pprProgress
 * Pointer pointer to progress structure. Progress structure will be
 * allocated here. 
 * @param[in] prFile
 * Where to log messages to
 * @param[in] pcPrefix
 * What prefix to use for messages
 * @param[in] bPrintCR
 * If TRUE carriage return instead of newline will be printed between log messages
 */    
void
NewProgress(progress_t **pprProgress, FILE *prFile,
            char *pcPrefix, bool bPrintCR)
{
    assert(NULL!=pprProgress);
    assert(NULL!=prFile);
    assert(NULL!=pcPrefix);
    
    (*pprProgress) = (progress_t *) CKMALLOC(1*sizeof(progress_t));
    (*pprProgress)->prFile = prFile;
    (*pprProgress)->bPrintCR = bPrintCR;    
    (*pprProgress)->pcPrefix = CkStrdup(pcPrefix);
    strcpy((*pprProgress)->pcLastLogMsg, "\0");
    (*pprProgress)->prStopwatch = StopwatchCreate();
    StopwatchZero((*pprProgress)->prStopwatch);
    StopwatchStart((*pprProgress)->prStopwatch);
    
    return; 
}
/***   end: NewProgress()   ***/



/**
 * @brief Frees progress structure and its members
 *
 * @param[out] pprProgress
 * Pointer pointer to progress structure
 *
 * @see NewProgress()
 * 
 */    
void
FreeProgress(progress_t **pprProgress)
{
    (*pprProgress)->prFile = NULL;
    CKFREE((*pprProgress)->pcPrefix);
    StopwatchFree((*pprProgress)->prStopwatch);

    CKFREE(*pprProgress);
    return; 
}
/***   end: FreeProgress()   ***/




/**
 * @brief Prints a progress update (and a carriage return)
 *
 * @param[in] prProgress
 * Pointer to the progress structure
 * @param[in] iStep
 * Current step number
 * @param[in] iTotalSteps
 * Total step number
 * @param[in] bForceOutput
 * If percentage hasn't changed output is normally supressed
 * normally. Output can be forced with this flag.
 *
 */    
void
ProgressLog(progress_t *prProgress, 
            unsigned long int iStep, unsigned long int iTotalSteps, 
            bool bForceOutput)
{
    char pcLogMsg[1024];
    assert(0!=iTotalSteps);
    
    if (rLog.iLogLevelEnabled>LOGLEVEL_THRESHOLD) {
        return;
    }
    
    (void) snprintf(pcLogMsg, sizeof(pcLogMsg), "%s: %lu %%",
                    prProgress->pcPrefix, (unsigned long int)(iStep/(float)iTotalSteps*100.0));

    if (! bForceOutput) {
        /* Skip logging, if we've just logged the same message */
        if (STR_EQ(pcLogMsg, prProgress->pcLastLogMsg)) {
            return;
        }
    }

    strncpy(prProgress->pcLastLogMsg, pcLogMsg,
            sizeof(prProgress->pcLastLogMsg));
    
    fprintf(prProgress->prFile, "%s (%lu out of %lu)", pcLogMsg, iStep, iTotalSteps);
    if (prProgress->bPrintCR) {
        fprintf(prProgress->prFile, "\r");
    } else {
        fprintf(prProgress->prFile, "\n");

    }
    (void) fflush(prProgress->prFile);

    return; 
}
/***   end: ProgressLog()   ***/


/**
 * @brief Finishes progress output by printing the elapsed time
 *
 * @param[in] prProgress
 * Pointer to the progress structure
 *
 */    
void
ProgressDone(progress_t *prProgress)
{
    char pcBuf[1024];

    if (rLog.iLogLevelEnabled>LOGLEVEL_THRESHOLD) {
        return;
    }

    (void) snprintf(pcBuf, sizeof(pcBuf), "%s done. CPU time: ", 
                    prProgress->pcPrefix);
    StopwatchStop(prProgress->prStopwatch);
    StopwatchDisplay(prProgress->prFile, pcBuf, prProgress->prStopwatch);
    (void) fflush(prProgress->prFile);

    return; 
}
/***   end: ProgressDone()   ***/
