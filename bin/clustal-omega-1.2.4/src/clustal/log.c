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
 *  RCS $Id$
 */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif


#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <assert.h>

#include "log.h"

/* a default standard logger */
log_t rLog;



void
LogVfprintf(FILE *prFP, char *pcFmt, va_list rVArgList);
void
LogWarn(FILE *prFP, char *pcFmt, va_list rVArgList);
void
LogError(FILE *prFP, char *pcFmt, va_list rVArgList);
void
LogCritical(FILE *prFP, char *pcFmt, va_list rVArgList);
void
LogFatal(FILE *prFP, char *pcFmt, va_list rVArgList);
void
LogForcedDebug(FILE *prFP, char *pcFmt, va_list rVArgList);




/**
 * @brief Plain, default print function
 *
 * Newline character is automatically appended to message.
 *
 */
void
LogVfprintf(FILE *prFP, char *pcFmt, va_list rVArgList)
{
	/* print prefix */
    vfprintf(prFP, pcFmt, rVArgList);
    fprintf(prFP,"\n");

#ifndef NDEBUG
    fflush(prFP);
#endif
}
/* end of LogVfprintf() */



/**
 * @brief 
 *
 */
void
LogWarn(FILE *prFP, char *pcFmt, va_list rVArgList)
{
	fprintf(prFP, "WARNING: ");
    LogVfprintf(prFP, pcFmt, rVArgList);
}
/* end of LogWarn() */



/**
 * @brief 
 *
 */
void
LogError(FILE *prFP, char *pcFmt, va_list rVArgList)
{
	fprintf(prFP, "ERROR: ");
    LogVfprintf(prFP, pcFmt, rVArgList);
}
/* end of LogError() */



/**
 * @brief 
 *
 */
void
LogCritical(FILE *prFP, char *pcFmt, va_list rVArgList)
{
	fprintf(prFP, "CRITICAL ERROR: ");
    LogVfprintf(prFP, pcFmt, rVArgList);
}
/* end of LogCritical() */



/**
 * @brief Will also exit!
 */
void
LogFatal(FILE *prFP, char *pcFmt, va_list rVArgList)
{
	fprintf(prFP, "FATAL: ");
    LogVfprintf(prFP, pcFmt, rVArgList);
    
    exit(EXIT_FAILURE);
}
/* end of LogFatal() */



/**
 * @brief 
 *
 */
void
LogForcedDebug(FILE *prFP, char *pcFmt, va_list rVArgList)
{
	fprintf(prFP, "FORCED DEBUG: ");
    LogVfprintf(prFP, pcFmt, rVArgList);
}
/* end of LogForcedDebug() */



/**
 *
 * @brief Sets up default function pointers
 *
 */
void
LogDefaultSetup(log_t *log)
{
    log->iLogLevelEnabled = LOG_WARN;

	log->prFP[LOG_DEBUG] = stdout;
	log->prFP[LOG_INFO] = stdout;
	log->prFP[LOG_WARN] = stderr;
	log->prFP[LOG_ERROR] = stderr;
	log->prFP[LOG_CRITICAL] = stderr;
	log->prFP[LOG_FATAL] = stderr;
	log->prFP[LOG_FORCED_DEBUG] = stderr;

	log->prFunc[LOG_DEBUG] = &LogVfprintf;
	log->prFunc[LOG_INFO] = &LogVfprintf;
	log->prFunc[LOG_WARN] = &LogWarn;
	log->prFunc[LOG_ERROR] = &LogError;
	log->prFunc[LOG_CRITICAL] = &LogCritical;
	log->prFunc[LOG_FATAL] = &LogFatal;
	log->prFunc[LOG_FORCED_DEBUG] = &LogForcedDebug;
}
/* end of LogDefaultSetup() */



/**
 * @brief Log to certain level
 *
 * See also comp.lang.c FAQ list · Question 15.12
 * http://c-faq.com/varargs/handoff.html How can I write a function which
 * takes a variable number of arguments and passes them to some other function
 * (which takes a variable number of arguments)?
 *
 */
void
Log(log_t *prLog, int iLevel, char *pcFmt, ...) 
{
	va_list rVArgList;
    void (*prFunc) (FILE *prFP, char *pcFormat, va_list rVArgList);
	FILE *prFP;

	assert(iLevel>=0 && iLevel<=LOG_NUM_LEVELS);

    /* fprintf(stderr, "DEBUG: iLevel=%d and iLogLevelEnabled=%d\n", iLevel, prLog->iLogLevelEnabled); */

	/* return if below current loglevel */
	if (iLevel < prLog->iLogLevelEnabled) {
		return;
	}

	prFunc = prLog->prFunc[iLevel];
	prFP = prLog->prFP[iLevel];

	/* return if muted */
	if (NULL == prFunc) {
		return;
	}

	va_start(rVArgList, pcFmt);
	prFunc(prFP, pcFmt, rVArgList);
	va_end(rVArgList);
}
/* end of Log() */



/**
 *
 * @brief Change file pointer for certain level
 *
 */
void
LogSetFP(log_t *prLog, int iLevel, FILE *prFP)
{
	assert(iLevel>=0 && iLevel<=LOG_NUM_LEVELS);

	prLog->prFP[iLevel] = prFP;
}
/* end of LogSetFP() */



/**
 *
 * @brief Return file pointer for certain level
 *
 */
FILE *
LogGetFP(log_t *prLog, int iLevel)
{
	assert(iLevel>=0 && iLevel<=LOG_NUM_LEVELS);

	return prLog->prFP[iLevel];
}
/* end of LogGetFP() */





/**
 *
 * @brief Change file pointer for all levels
 *
 */
void
LogSetFPForAll(log_t *prLog, FILE *prFP)
{
	int iAux;

	for (iAux=0; iAux<LOG_NUM_LEVELS; iAux++) {
		prLog->prFP[iAux] = prFP;
	}
}
/* end of LogSetFP() */



/**
 *
 * @brief Mute certain level (i.e set the corresponding function to NULL)
 *
 */
void
LogMute(log_t *prLog, int iLevel)
{
	assert(iLevel>=0 && iLevel<=LOG_NUM_LEVELS);

	prLog->prFunc[iLevel] = NULL;
}
/* end of LogMute() */


/**
 *
 * @brief Mute all channels
 *
 */
void
LogMuteAll(log_t *prLog)
{
	int iAux;

	for (iAux=0; iAux<LOG_NUM_LEVELS; iAux++) {
		LogMute(prLog, iAux);
	}
}
/* end of LogMuteAll() */



/**
 * @brief
 *
 */
void
LogFuncOverwrite(log_t *prLog, int iLevel, 
                 void (*Func) (FILE *prFP, char *pcFormat, va_list rVArgList))
{
	assert(iLevel>=0 && iLevel<=LOG_NUM_LEVELS);

    prLog->prFunc[iLevel] = Func;

}
/* end of LogFuncOverwrite() */



#ifdef LOG_TEST


#define TEXT "Lorem ipsum dolor sit amet"

void
PrintSomeTextToAll(log_t *prLog) {
	int iAux;
	for (iAux=0; iAux<LOG_NUM_LEVELS; iAux++) {
		Log(prLog, iAux, TEXT);
	}
}


int
main(int argc, char**argv) {
	log_t prLog;
	log_t pr2ndLog;
	FILE *prFP;
	char *pcTmpFileName = "schmock.txt";

	prFP = fopen(pcTmpFileName, "w");

	 
	LogDefaultSetup(&prLog);
	LogDefaultSetup(&pr2ndLog);

	printf("Printing to log:\n");
	PrintSomeTextToAll(&prLog);
	printf("---\n");

	printf("Printing to pr2ndLog:\n");
	PrintSomeTextToAll(&prLog);
	printf("---\n");


	
	printf("Changing second log's FP to write to '%s'\n", pcTmpFileName);
	LogSetFPForAll(&pr2ndLog, prFP);

	printf("Printing to pr2ndLog:\n");
	PrintSomeTextToAll(&pr2ndLog);
	printf("---\n");


	printf("Changing Info() to new function (Fatal()) in log:\n");
    LogFuncOverwrite(&prLog, LOG_INFO, &LogFatal);

	printf("Printing to log:\n");
	PrintSomeTextToAll(&prLog);
	printf("---\n");


   
	fclose(prFP);

	return 0;
}


#endif
