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


#include <stdio.h>
#include <stdarg.h>

#ifndef LOG_H
#define LOG_H


#define LOG_DEBUG 0
#define LOG_VERBOSE 1
#define LOG_INFO 2
#define LOG_WARN 3
#define LOG_FORCED_DEBUG 4
#define LOG_ERROR 5
#define LOG_CRITICAL 6
#define LOG_FATAL 7

#define LOG_NUM_LEVELS 8


typedef struct {
    /* the higher the level, the more priority it has. numbers must be
     *  sequential 
     */
    
    /* array of function pointers */
    void (*prFunc[LOG_NUM_LEVELS]) (FILE *prFP, char *pcFormat, va_list rVArgList);
    FILE *prFP[LOG_NUM_LEVELS];
    char *prPrefix[LOG_NUM_LEVELS];

    /* everything above this level will be printed */
    int iLogLevelEnabled;
} log_t;



/* a standard logger */
extern log_t rLog;



void
LogDefaultSetup(log_t *log);
void
Log(log_t *prLog, int iLevel, char *pcFmt, ...);
void
LogSetFP(log_t *log, int level, FILE *fp);
void
LogSetFPForAll(log_t *log, FILE *fp);
FILE *
LogGetFP(log_t *prLog, int iLevel);
void
LogMute(log_t *log, int level);
void
LogMuteAll(log_t *log);
void
LogFuncOverwrite(log_t *prLog, int iLevel, 
                 void (*Func) (FILE *prFP, char *pcFormat, va_list rVArgList));


#endif
