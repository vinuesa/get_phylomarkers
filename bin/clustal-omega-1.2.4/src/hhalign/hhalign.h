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
 * RCS $Id: hhalign.h 296 2014-10-07 12:15:41Z fabian $
 */


extern int 
hhalign(char **ppcFirstProf, int iFirstCnt, double *pdWeightsL, 
        char **ppcSecndProf, int iSecndCnt, double *pdWeightsR, 
        double *dScore_p, hmm_light *prHMM, hmm_light *prHMM2, 
        char *pcPrealigned1, char *pcRepresent1,
        char *pcPrealigned2, char *pcRepresent2,
        hhalign_para rHhalignPara, hhalign_scores *rHHscores, 
        int iFlag, int iVerbosity,
        char zcAux[], char zcError[]);
