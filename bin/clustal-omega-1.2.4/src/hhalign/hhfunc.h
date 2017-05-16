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
 *  RCS $Id: hhfunc.h 290 2013-09-20 15:18:12Z fabian $
 */


int readHMMWrapper(hmm_light *rHMM_p,
		     char *pcHMM_input);
void FreeHMMstruct(hmm_light *prHMM);
int AlnToHMM2(hmm_light *rHMM_p, hhalign_para rHhalignPara, 
              char **ppcSeq, int iN);
int HHMake_Wrapper(char *tmp_aln, char *hmm_out);
