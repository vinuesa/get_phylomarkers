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
 * RCS $Id: general.h 290 2013-09-20 15:18:12Z fabian $
 */


#ifndef GENERAL_H
#define GENERAL_H


#include "../clustal/log.h"
#include <stdbool.h>

/*
***** Omega definitions************************************
 FS, 2010-02-19
*/
enum {NO = 0, YES};
enum {BASE10 = 10};
enum {AMINOACIDS = 20, STATE_TRANSITIONS = 7};
enum {MAXWORD = 100, MAXLEN = 10000};
enum {OVER_ALLOCATE = 2};
enum {FAILURE = -1, OK};
enum {RETURN_OK = 0, RETURN_FROM_MAC, RETURN_FROM_VITERBI, RETURN_FROM_PRINT_ALI, RETURN_FROM_RNP, RETURN_UNKNOWN};
enum {REALLY_BIG_MEMORY_MB = 64000};
enum {F_OFFSET = 1};
enum {INTERN_ALN_2_HMM = 0, READ_ALN_2_HMM, READ_HMM_2_HMM, INTERN_HMM_2_HMM};
#define UNITY 1.00
enum {CALL_FROM_DONT_KNOW = -1, CALL_FROM_ALN_HMM, CALL_FROM_REGULAR};

/*#define MIN(a,b) ((a)<(b)?(a):(b))*/


/* parameters passed from Clustal-Omega to hhalign; FS, r240 -> */
typedef struct {

    int iMacRamMB; /* dedicated amount of RAM for Maximum Accuracy (in MB) */
	bool bIsDna; /* indicates we're in nucleotide mode */
	bool bIsRna; /* indicates we're in nucleotide mode */
    double pca;
    double pcb;
    double pcc;
    double pcw;
    double gapb;
    double gapd;
    double gape;
    double gapf;
    double gapg;
    double gaph;
    double gapi;
    double pcaV;
    double pcbV;
    double pccV;
    double pcwV;
    double gapbV;
    double gapdV;
    double gapeV;
    double gapfV;
    double gapgV;
    double gaphV;
    double gapiV;

} hhalign_para;

/* 'scores' passed out from hhalign */
typedef struct{

	double forwardProb;
    double backwardProb;
    double hhScore;
    double sumPP;
    double *PP;
    int L;

} hhalign_scores;


typedef struct {
  /***public***/
  int n_display;
  char **sname;
  char **seq;
  int ncons;
  int nfirst;
  int nss_dssp;
  int nsa_dssp;
  int nss_pred;
  int nss_conf;
  int L;
  int N_in;
  int N_filtered;
  float *Neff_M;
  float *Neff_I;
  float *Neff_D;
  float Neff_HMM;
  char *longname;
  char name[511];
  char file[511];
  char fam[511];
  char sfam[511];
  char fold[511];
  char cl[511];
  float lamda;
  float mu;
  /***private***/
  float **f;
  float **g;
  float **p;
  float **tr;
  float **linTr;
  char trans_lin;
  char *ss_dssp;
  char *sa_dssp;
  char *ss_pred;
  char *ss_conf;
  char *Xcons;
  float pav[20];
  float pnul[20];
  int *l;

} hmm_light;


#endif

