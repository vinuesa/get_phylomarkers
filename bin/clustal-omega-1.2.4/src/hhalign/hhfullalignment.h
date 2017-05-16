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
 * RCS $Id: hhfullalignment.h 284 2013-06-12 10:10:11Z fabian $
 */

//////////////////////////////////////////////////////////////////////////////
// Class for output alignment of query against template sequences
//////////////////////////////////////////////////////////////////////////////

class FullAlignment
{
public:
  FullAlignment(int maxseqdis=MAXSEQDIS);
  ~FullAlignment();
  void FreeMemory();
  int Build(HMM& q, Hit& hit, char zcError[]);
  void PrintHeader(FILE* outf, HMM& q, Hit& hit);
  void PrintHHR(FILE* outf, Hit& hit);
  void PrintA2M(FILE* outf, Hit& hit);
  void PrintFASTA(FILE* outf, Hit& hit);
  void PrintA3M(FILE* outf, Hit& hit);
  void OverWriteSeqs(char **ppcFirstProf, char **ppcSecndProf);
  int identities;      // number of identical residues in query and template sequence
  float score_sim;     // substitution matrix similarity score between query and template

private:
  HalfAlignment* qa; //query and template parts of the alignment
  HalfAlignment* ta; //query and template parts of the alignment
  char symbol[LINELEN];         //symbol[h] = symbol (= - . + |) indicating match score for col h of alignment    
  void ClearSymbols()      {for (int h=0; h<LINELEN-1; h++) symbol[h]=' ';}
  void AddColumns(int i, int j, char prev_state, char state, float S);
  void AddGaps();
  int ScoreChr(float S) {return (S<-1.5?'=':(S<-0.5?'-':(S<0.5?'.':(S<1.5?'+':'|'))));}
};
