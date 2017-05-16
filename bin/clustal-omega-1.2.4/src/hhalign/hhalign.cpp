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
 *  RCS $Id: hhalign.cpp 309 2016-06-13 14:10:11Z fabian $
 */

/* hhalign.C: 
 * Align a multiple alignment to an alignment or HMM 
 * Print out aligned input sequences in a3m format
 * Compile:              g++ hhalign.C -o hhalign -I/usr/include/ -L/usr/lib -lpng -lz -O3 -fno-strict-aliasing 
 * Compile with efence:  g++ hhalign.C -o hhalign -I/usr/include/ -lefence -L/usr/lib -lpng -lz -O -g  
 * 
 * Error codes: 0: ok  1: file format error  2: file access error  
 *              3: memory error  4: internal numeric error  
 *              5: command line error
 */
#undef PNG           /* include options for making png files? 
			(will need the png library) */
#define MAIN
#include <iostream>   // cin, cout, cerr
#include <fstream>    // ofstream, ifstream 
#include <string.h>     // strcmp, strstr
#include <stdio.h>    // printf
#include <stdlib.h>   // exit
#include <math.h>     // sqrt, pow
#include <limits.h>   // INT_MIN
#include <float.h>    // FLT_MIN
#include <time.h>     // clock
#include <ctype.h>    // islower, isdigit etc

using std::cout;
using std::cerr;
using std::endl;
using std::ios;
using std::ifstream;
using std::ofstream;

int iAux_GLOBAL;

#include "general.h"
#include "util-C.h"     /* imax, fmax, iround, iceil, ifloor, strint, strscn, 
			 strcut, substr, uprstr, uprchr, Basename etc. */
#include "list-C.h"     // list data structure
#include "hash-C.h"     // hash data structure

#include "hhdecl-C.h"      // Constants, global variables, struct Parameters
#include "hhutil-C.h"      /* MatchChr, InsertChr, aa2i, i2aa, log2, 
			    fast_log2, WriteToScreen, */
#include "hhmatrices-C.h"  // BLOSUM50, GONNET, HSDM
#include "hhhmm.h"       // class HMM
#include "hhhit.h"       // class Hit
#include "hhalignment.h" // class Alignment
#include "hhhalfalignment.h" // class HalfAlignment
#include "hhfullalignment.h" // class FullAlignment
#include "hhhitlist.h"   // class Hit
#include "hhhmm-C.h"       // class HMM
#include "hhalignment-C.h" // class Alignment
#include "hhhit-C.h"       // class Hit
#include "hhhalfalignment-C.h" // class HalfAlignment
#include "hhfullalignment-C.h" // class FullAlignment
#include "hhhitlist-C.h"   // class HitList 
#include "hhfunc-C.h"      // some functions common to hh programs

#ifdef PNG
#include "pngwriter.h"   //PNGWriter (http://pngwriter.sourceforge.net/)
#include "pngwriter.cc"  //PNGWriter (http://pngwriter.sourceforge.net/)
#endif	    

////////////////////////////////////////////////////////////////////////
// Global variables 
////////////////////////////////////////////////////////////////////////
HMM *q;                 // Create query HMM with maximum of MAXRES match states
HMM *t;                 /* Create template HMM with maximum of MAXRES 
			  match states  */
Alignment *qali;        /* (query alignment might be needed outside of hhfunc.C 
			  for -a option) */
Hit hit;               // Ceate new hit object pointed at by hit
HitList hitlist;       /* list of hits with one Hit object for each 
			  pairwise comparison done */
char aliindices[256];  /* hash containing indices of all alignments 
			  which to show in dot plot */
char* dmapfile=NULL;   /* where to write the coordinates for the HTML map file 
			  (to click the alignments) */
char* alitabfile=NULL; // where to write pairs of aligned residues
char* strucfile=NULL;  // where to read structure scores
char* pngfile=NULL;    // pointer to pngfile
char* tcfile=NULL;     // TCoffee output file name
float probmin_tc=0.05; /* 5% minimum posterior probability for printing 
			  pairs of residues for TCoffee */

int dotW=10;           // average score of dot plot over window [i-W..i+W]
float dotthr=0.5;      // score threshold for dot plot
int dotscale=600;      // size scale of dotplot
char dotali=0;         // show no alignments in dotplot
float dotsat=0.3;      // saturation of grid and alignments in dot plot
float pself=0.001;     // maximum p-value of 2nd and following self-alignments
int Nstochali=0;       // # of stochastically traced alignments in dot plot
float** Pstruc=NULL;   /* structure matrix which can be multiplied to prob 
			  ratios from aa column comparisons in Forward() */
float** Sstruc=NULL;   /* structure matrix which can be added to log odds 
			  from aa column comparisons in Forward() */
bool nucleomode = false;

////////////////////////////////////////////////////////////////////////////
// Help functions
////////////////////////////////////////////////////////////////////////////
/**
 * @brief  general help function, not accessible in Clustal-Omega
 */
#if 0
void 
help()
{
  printf("\n");
  printf("HHalign %s\n",VERSION_AND_DATE);
  printf("Align a query alignment/HMM to a template alignment/HMM by HMM-HMM alignment\n");
  printf("If only one alignment/HMM is given it is compared to itself and the best\n");
  printf("off-diagonal alignment plus all further non-overlapping alignments above \n");
  printf("significance threshold are shown.\n");
  printf("%s",REFERENCE);
  printf("%s",COPYRIGHT);
  printf("\n");
  printf("Usage: %s -i query [-t template] [options]  \n",program_name);
  printf(" -i <file>     input query alignment  (fasta/a2m/a3m) or HMM file (.hhm)\n");
  printf(" -t <file>     input template alignment (fasta/a2m/a3m) or HMM file (.hhm)\n");
#ifdef PNG
  printf(" -png <file>   write dotplot into PNG-file (default=none)           \n");
#endif
  printf("\n");         
  printf("Output options:                                                           \n");
  printf(" -o <file>     write output alignment to file\n"); 
  printf(" -ofas <file>  write alignments in FASTA, A2M (-oa2m) or A3M (-oa3m) format   \n"); 
  printf(" -a <file>     write query alignment in a3m format to file (default=none)\n");
  printf(" -aa <file>    append query alignment in a3m format to file (default=none)\n");
  printf(" -atab <file>  write alignment as a table (with posteriors) to file (default=none)\n");
  printf(" -v <int>      verbose mode: 0:no screen output  1:only warings  2: verbose\n");
  printf(" -seq  [1,inf[ max. number of query/template sequences displayed  (def=%i)  \n",par.nseqdis);
  printf(" -nocons       don't show consensus sequence in alignments (default=show) \n");
  printf(" -nopred       don't show predicted 2ndary structure in alignments (default=show) \n");
  printf(" -nodssp       don't show DSSP 2ndary structure in alignments (default=show) \n");
  printf(" -aliw int     number of columns per line in alignment list (def=%i)\n",par.aliwidth);
  printf(" -P <float>    for self-comparison: max p-value of alignments (def=%.2g\n",pself);
  printf(" -p <float>    minimum probability in summary and alignment list (def=%G) \n",par.p);
  printf(" -E <float>    maximum E-value in summary and alignment list (def=%G)     \n",par.E);
  printf(" -Z <int>      maximum number of lines in summary hit list (def=%i)       \n",par.Z);
  printf(" -z <int>      minimum number of lines in summary hit list (def=%i)       \n",par.z);
  printf(" -B <int>      maximum number of alignments in alignment list (def=%i)    \n",par.B);
  printf(" -b <int>      minimum number of alignments in alignment list (def=%i)    \n",par.b);
  printf(" -rank int     specify rank of alignment to write with -a or -aa option (default=1)\n");
  printf("\n");         
#ifdef PNG
  printf("Dotplot options:\n");
  printf(" -dwin <int>   average score in dotplot over window [i-W..i+W] (def=%i)      \n",dotW);
  printf(" -dthr <float> score threshold for dotplot (default=%.2f)                    \n",dotthr);
  printf(" -dsca <int>   if value <= 20: size of dot plot unit box in pixels           \n");
  printf("               if value > 20: maximum dot plot size in pixels (default=%i)   \n",dotscale);
  printf(" -dali <list>  show alignments with indices in <list> in dot plot            \n");
  printf("               <list> = <index1> ... <indexN>  or  <list> = all              \n");
  printf("\n");         
#endif
  printf("Filter input alignment (options can be combined):                         \n");
  printf(" -id   [0,100] maximum pairwise sequence identity (%%) (def=%i)   \n",par.max_seqid);
  printf(" -diff [0,inf[ filter most diverse set of sequences, keeping at least this    \n");
  printf("               many sequences in each block of >50 columns (def=%i)\n",par.Ndiff);
  printf(" -cov  [0,100] minimum coverage with query (%%) (def=%i) \n",par.coverage);
  printf(" -qid  [0,100] minimum sequence identity with query (%%) (def=%i) \n",par.qid);
  printf(" -qsc  [0,100] minimum score per column with query  (def=%.1f)\n",par.qsc);
//   printf(" -csc  [0,100] minimum score per column with core alignment (def=%-.2f)\n",par.coresc);
//   printf(" -qscc [0,100] minimum score per column of core sequence with query (def=%-.2f)\n",par.qsc_core);
  printf("\n");         
  printf("Input alignment format:                                                     \n");
  printf(" -M a2m        use A2M/A3M (default): upper case = Match; lower case = Insert;\n");         
  printf("               '-' = Delete; '.' = gaps aligned to inserts (may be omitted)   \n");
  printf(" -M first      use FASTA: columns with residue in 1st sequence are match states\n");
  printf(" -M [0,100]    use FASTA: columns with fewer than X%% gaps are match states   \n");
  printf("\n");    
  printf("HMM-HMM alignment options:                                                  \n");
  printf(" -glob/-loc    global or local alignment mode (def=local)         \n");
  printf(" -alt <int>    show up to this number of alternative alignments (def=%i)    \n",par.altali);
  printf(" -vit          use Viterbi algorithm for alignment instead of MAC algorithm \n");
  printf(" -mac          use Maximum Accuracy (MAC) alignment (default)  \n");
  printf(" -mact [0,1[   posterior probability threshold for MAC alignment (def=%.3f) \n",par.mact);
  printf("               A threshold value of 0.0 yields global alignments.\n");
  printf(" -sto <int>    use global stochastic sampling algorithm to sample this many alignments\n");
  printf(" -excl <range> exclude query positions from the alignment, e.g. '1-33,97-168'\n");
  printf(" -shift [-1,1] score offset (def=%-.3f)                                      \n",par.shift);
  printf(" -corr [0,1]   weight of term for pair correlations (def=%.2f)               \n",par.corr);
  printf(" -ssm  0-4     0:no ss scoring [default=%i]               \n",par.ssm);
  printf("               1:ss scoring after alignment                                  \n");
  printf("               2:ss scoring during alignment                                 \n");
  printf(" -ssw  [0,1]   weight of ss score  (def=%-.2f)                               \n",par.ssw);
  printf("\n");
  printf(" -def          read default options from ./.hhdefaults or <home>/.hhdefault. \n");
  printf("\n");
  printf("Example: %s -i T0187.a3m -t d1hz4a_.hhm -png T0187pdb.png \n",program_name);
  cout<<endl;
//   printf("More help:                                                         \n");
//   printf(" -h out        output options                                      \n");
//   printf(" -h hmm        options for building HMM from multiple alignment    \n");
//   printf(" -h gap        options for setting gap penalties                   \n");
//   printf(" -h ali        options for HMM-HMM alignment                       \n");
//   printf(" -h all        all options \n");
}
#endif

/**
 * @brief  helpt for output, not accessible in Clustal-Omega
 */
#if 0
void 
help_out()
{
  printf("\n");
  printf("Output options: \n");
  printf(" -v            verbose mode (default: show only warnings)                 \n");
  printf(" -v 0          suppress all screen output                                 \n");
  printf(" -nocons       don't show consensus sequence in alignments (default=show) \n");
  printf(" -nopred       don't show predicted 2ndary structure in alignments (default=show) \n");
  printf(" -nodssp       don't show DSSP SS 2ndary structure in alignments (default=show) \n");
  printf(" -seq  [1,inf[ max. number of query/template sequences displayed  (def=%i)  \n",par.nseqdis);
  printf(" -aliw [40,..[ number of columns per line in alignment list (def=%i)\n",par.aliwidth);
  printf(" -P <float>    for self-comparison: max p-value of alignments (def=%.2g\n",pself);
  printf(" -p <float>    minimum probability in summary and alignment list (def=%G) \n",par.p);
  printf(" -E <float>    maximum E-value in summary and alignment list (def=%G)     \n",par.E);
  printf(" -Z <int>      maximum number of lines in summary hit list (def=%i)       \n",par.Z);
  printf(" -z <int>      minimum number of lines in summary hit list (def=%i)       \n",par.z);
  printf(" -B <int>      maximum number of alignments in alignment list (def=%i)    \n",par.B);
  printf(" -b <int>      minimum number of alignments in alignment list (def=%i)    \n",par.b);
  printf(" -rank int     specify rank of alignment to write with -a or -aa option (def=1)\n");
  printf(" -tc <file>    write a TCoffee library file for the pairwise comparison   \n");         
  printf(" -tct [0,100]  min. probobability of residue pairs for TCoffee (def=%i%%)\n",iround(100*probmin_tc));         
  printf("\n");         
#ifdef PNG
  printf("Dotplot options:\n");
  printf(" -dwin int     average score in dotplot over window [i-W..i+W] (def=%i)   \n",dotW);
  printf(" -dthr float   score threshold for dotplot (default=%.2f)                 \n",dotthr);
  printf(" -dsca int     size of dot plot box in pixels  (default=%i)               \n",dotscale);
  printf(" -dali <list>  show alignments with indices in <list> in dot plot\n");
  printf("               <list> = <index1> ... <indexN>  or  <list> = all              \n");
  printf(" -dmap <file>  print list of coordinates in png plot  \n");
#endif
}
#endif

/**
 * @brief  help hit HMM options, not accessible in Clustal-Omega
 */
#if 0
void 
help_hmm()
{
  printf("\n");
  printf("Options to filter input alignment (options can be combined):              \n");
  printf(" -id   [0,100] maximum pairwise sequence identity (%%) (def=%i)   \n",par.max_seqid);
  printf(" -diff [0,inf[ filter most diverse set of sequences, keeping at least this    \n");
  printf("               many sequences in each block of >50 columns (def=%i)\n",par.Ndiff);
  printf(" -cov  [0,100] minimum coverage with query (%%) (def=%i) \n",par.coverage);
  printf(" -qid  [0,100] minimum sequence identity with query (%%) (def=%i) \n",par.qid);
  printf(" -qsc  [0,100] minimum score per column with query  (def=%.1f)\n",par.qsc);
//   printf(" -csc  [0,100] minimum score per column with core alignment (def=%-.2f)\n",par.coresc);
//   printf(" -qscc [0,100] minimum score per column of core sequence with query (def=%-.2f)\n",par.qsc_core);
  printf("                                                                          \n");
  printf("HMM-building options:                                                     \n");
  printf(" -M a2m        use A2M/A3M (default): upper case = Match; lower case = Insert;\n");         
  printf("               '-' = Delete; '.' = gaps aligned to inserts (may be omitted)   \n");
  printf(" -M first      use FASTA: columns with residue in 1st sequence are match states\n");
  printf(" -M [0,100]    use FASTA: columns with fewer than X%% gaps are match states   \n");
  printf(" -tags         do NOT neutralize His-, C-myc-, FLAG-tags, and \n");
  printf("               trypsin recognition sequence to background distribution    \n");
  printf("                                                                          \n");
  printf("Pseudocount options:                                                      \n");
  printf(" -Gonnet       use the Gonnet substitution matrix (default)               \n");
  printf(" -Blosum50     use the Blosum50 substitution matrix                       \n");
  printf(" -Blosum62     use the Blosum62 substitution matrix                       \n");
  printf(" -HSDM         use the structure-derived HSDM substitution matrix         \n");
  printf(" -pcm  0-3     Pseudocount mode (default=%-i)                             \n",par.pcm);
  printf("               tau = substitution matrix pseudocount admixture            \n");
  printf("               0: no pseudo counts:     tau = 0                           \n");
  printf("               1: constant              tau = a                           \n");
  printf("               2: divergence-dependent: tau = a/(1 + ((Neff-1)/b)^c)       \n");
  printf("                  Neff=( (Neff_q^d+Neff_t^d)/2 )^(1/d)                       \n");
  printf("                  Neff_q = av number of different AAs per column in query  \n");
  printf("               3: column-specific:      tau = \'2\' * (Neff(i)/Neff)^(w*Neff/20)\n");
  printf(" -pca  [0,1]   set a (overall admixture) (def=%-.1f)                      \n",par.pca);
  printf(" -pcb  [1,inf[ set b (threshold for Neff) (def=%-.1f)                      \n",par.pcb);
  printf(" -pcc  [0,3]   set c (extinction exponent for tau(Neff))  (def=%-.1f)      \n",par.pcc);
  printf(" -pcw  [0,3]   set w (weight of pos-specificity for pcs) (def=%-.1f)      \n",par.pcw);
}
#endif

/**
 * @brief  help with gap handling, not accessible in Clustal-Omega
 */
#if 0
void 
help_gap()
{
  printf("\n");
  printf("Gap cost options:                                                         \n");
  printf(" -gapb [0,inf[ transition pseudocount admixture (def=%-.2f)               \n",par.gapb);
  printf(" -gapd [0,inf[ Transition pseudocount admixture for opening gap (default=%-.2f)\n",par.gapd);
  printf(" -gape [0,1.5] Transition pseudocount admixture for extending gap (def=%-.1f)\n",par.gape);
  printf(" -gapf ]0,inf] factor for increasing/reducing the gap open penalty for deletes (def=%-.2f)\n",par.gapf);
  printf(" -gapg ]0,inf] factor for increasing/reducing the gap open penalty for deletes (def=%-.2f)\n",par.gapg);
  printf(" -gaph ]0,inf] factor for increasing/reducing the gap extension penalty for deletes(def=%-.2f)\n",par.gaph);
  printf(" -gapi ]0,inf] factor for increasing/reducing the gap extension penalty for inserts(def=%-.2f)\n",par.gapi);
  printf(" -egq  [0,inf[ penalty (bits) for end gaps aligned to query residues (def=%-.2f)\n",par.egq);
  printf(" -egt  [0,inf[ penalty (bits) for end gaps aligned to template residues (def=%-.2f)\n",par.egt);
}
#endif

/**
 * @brief  help with alignment options, not accessible in Clustal-Omega
 */
#if 0
void 
help_ali()
{
  printf("\n");
  printf("Alignment options:  \n");
  printf(" -glob/-loc    global or local alignment mode (def=global)                \n");
  printf(" -mac          use Maximum Accuracy (MAC) alignment instead of Viterbi\n");
  printf(" -mact [0,1]   posterior prob threshold for MAC alignment (def=%.3f)      \n",par.mact);
  printf(" -sto <int>    use global stochastic sampling algorithm to sample this many alignments\n");
  printf(" -sc   <int>   amino acid score         (tja: template HMM at column j) (def=%i)\n",par.columnscore);
  printf("        0      = log2 Sum(tja*qia/pa)   (pa: aa background frequencies)    \n");
  printf("        1      = log2 Sum(tja*qia/pqa)  (pqa = 1/2*(pa+ta) )               \n");
  printf("        2      = log2 Sum(tja*qia/ta)   (ta: av. aa freqs in template)     \n");
  printf("        3      = log2 Sum(tja*qia/qa)   (qa: av. aa freqs in query)        \n");
   printf(" -corr [0,1]   weight of term for pair correlations (def=%.2f)            \n",par.corr);
  printf(" -shift [-1,1] score offset (def=%-.3f)                                   \n",par.shift);
  printf(" -r            repeat identification: multiple hits not treated as independent\n");
  printf(" -ssm  0-2     0:no ss scoring [default=%i]               \n",par.ssm);
  printf("               1:ss scoring after alignment                               \n");
  printf("               2:ss scoring during alignment                              \n");
  printf(" -ssw  [0,1]   weight of ss score compared to column score (def=%-.2f)    \n",par.ssw);
  printf(" -ssa  [0,1]   ss confusion matrix = (1-ssa)*I + ssa*psipred-confusion-matrix [def=%-.2f)\n",par.ssa);
}
#endif

/**
 * @brief  general help menu, not accessible in Clustal-Omega
 */
#if 0
void 
help_all()
{
  help();
  help_out();
  help_hmm();
  help_gap();
  help_ali();
  printf("\n");
  printf("Default options can be specified in './.hhdefaults' or '~/.hhdefaults'\n");
}
#endif

/////////////////////////////////////////////////////////////////////////////
//// Processing input options from command line and .hhdefaults file
/////////////////////////////////////////////////////////////////////////////
/**
 * @brief  process arguments from commandline, not accessible from Clustal-Omega
 */
#if 0
void 
ProcessArguments(int argc, char** argv)
{
  //Processing command line input
  for (int i=1; i<argc; i++)
    { 
      if (v>=4) cout<<i<<"  "<<argv[i]<<endl; //PRINT
      if (!strcmp(argv[i],"-i"))
	{
	  if (++i>=argc || argv[i][0]=='-') 
	    {help(); cerr<<endl<<"Error in "<<program_name<<": no query file following -i\n"; exit(4);}
	  else strcpy(par.infile,argv[i]);
	}
      else if (!strcmp(argv[i],"-t"))
	{
	  if (++i>=argc || argv[i][0]=='-') 
	    {help(); cerr<<endl<<"Error in "<<program_name<<": no template file following -d\n"; exit(4);}
	  else strcpy(par.tfile,argv[i]); /* FS, ProcessArguments */
	}
      else if (!strcmp(argv[i],"-o"))
	{
	  if (++i>=argc) 
	    {help(); cerr<<endl<<"Error in "<<program_name<<": no filename following -o\n"; exit(4);}
	  else strcpy(par.outfile,argv[i]);
	}
      else if (!strcmp(argv[i],"-ofas"))
	{
	  par.outformat=1;
	  if (++i>=argc || argv[i][0]=='-') 
	    {help() ; cerr<<endl<<"Error in "<<program_name<<": no output file following -o\n"; exit(4);}
	  else strcpy(par.pairwisealisfile,argv[i]);
	}
      else if (!strcmp(argv[i],"-oa2m"))
	{
	  par.outformat=2;
	  if (++i>=argc || argv[i][0]=='-') 
	    {help() ; cerr<<endl<<"Error in "<<program_name<<": no output file following -o\n"; exit(4);}
	  else strcpy(par.pairwisealisfile,argv[i]);
	}
      else if (!strcmp(argv[i],"-oa3m"))
	{
	  par.outformat=3;
	  if (++i>=argc || argv[i][0]=='-') 
	    {help() ; cerr<<endl<<"Error in "<<program_name<<": no output file following -o\n"; exit(4);}
	  else strcpy(par.pairwisealisfile,argv[i]);
	}
      else if (!strcmp(argv[i],"-rank") && (i<argc-1)) par.hitrank=atoi(argv[++i]); 
      else if (!strcmp(argv[i],"-Oa3m"))
	{
	  par.append=0;
	  if (++i>=argc || argv[i][0]=='-') 
	    {help() ; cerr<<endl<<"Error in "<<program_name<<": no output file following -Oa3m\n"; exit(4);}
	  else strcpy(par.alnfile,argv[i]);
	}
      else if (!strcmp(argv[i],"-Aa3m"))
	{
	  par.append=1;
	  if (++i>=argc || argv[i][0]=='-') 
	    {help() ; cerr<<endl<<"Error in "<<program_name<<": no output file following -Aa3m\n"; exit(4);}
	  else strcpy(par.alnfile,argv[i]);
	}
      else if (!strcmp(argv[i],"-Ohhm"))
	{
	  par.append=0;
	  if (++i>=argc || argv[i][0]=='-') 
	    {help() ; cerr<<endl<<"Error in "<<program_name<<": no output file following -Ohhm\n"; exit(4);}
	  else strcpy(par.hhmfile,argv[i]);
	}
      else if (!strcmp(argv[i],"-Ahhm"))
	{
	  par.append=1;
	  if (++i>=argc || argv[i][0]=='-') 
	    {help() ; cerr<<endl<<"Error in "<<program_name<<": no output file following -Ahhm\n"; exit(4);}
	  else strcpy(par.hhmfile,argv[i]);
	}
      else if (!strcmp(argv[i],"-Opsi"))
	{
	  par.append=0;
	  if (++i>=argc || argv[i][0]=='-') 
	    {help() ; cerr<<endl<<"Error in "<<program_name<<": no output file following -Opsi\n"; exit(4);}
	  else strcpy(par.psifile,argv[i]);
	}
      else if (!strcmp(argv[i],"-Apsi"))
	{
	  par.append=1;
	  if (++i>=argc || argv[i][0]=='-') 
	    {help() ; cerr<<endl<<"Error in "<<program_name<<": no output file following -Apsi\n"; exit(4);}
	  else strcpy(par.psifile,argv[i]);
	}
      else if (!strcmp(argv[i],"-png"))
	{
	  if (++i>=argc) 
	    {help(); cerr<<endl<<"Error in "<<program_name<<": no filename following -png\n"; exit(4);}
	  else 
	    {
	      pngfile = new(char[strlen(argv[i])+1]);
	      strcpy(pngfile,argv[i]);
	    }
	}
      else if (!strcmp(argv[i],"-Struc"))
	{
	  if (++i>=argc || argv[i][0]=='-') 
	    {help(); cerr<<endl<<"Error in "<<program_name<<": no query file following -Struc\n"; exit(4);}
	  else 
	    {
	      strucfile = new(char[strlen(argv[i])+1]);
	      strcpy(strucfile,argv[i]);
	    }
	}
      else if (!strcmp(argv[i],"-atab") || !strcmp(argv[i],"-Aliout"))
	{
	  if (++i>=argc || argv[i][0]=='-') 
	    {help(); cerr<<endl<<"Error in "<<program_name<<": no query file following -Struc\n"; exit(4);}
	  else 
	    {
	      alitabfile = new(char[strlen(argv[i])+1]);
	      strcpy(alitabfile,argv[i]);
	    }
	}
      else if (!strcmp(argv[i],"-tc"))
	{
	  if (++i>=argc || argv[i][0]=='-') 
	    {help() ; cerr<<endl<<"Error in "<<program_name<<": no output file following -Opsi\n"; exit(4);}
	  else 
	    {
	      tcfile = new(char[strlen(argv[i])+1]);
	      strcpy(tcfile,argv[i]);
	    }
	}
      else if (!strcmp(argv[i],"-h")|| !strcmp(argv[i],"--help"))
	{
	  if (++i>=argc) {help(); exit(0);} 
	  if (!strcmp(argv[i],"out")) {help_out(); exit(0);} 
	  if (!strcmp(argv[i],"hmm")) {help_hmm(); exit(0);} 
	  if (!strcmp(argv[i],"gap")) {help_gap(); exit(0);} 
	  if (!strcmp(argv[i],"ali")) {help_ali(); exit(0);} 
	  if (!strcmp(argv[i],"all")) {help_all(); exit(0);} 
	  else {help(); exit(0);}
	}
      else if (!strcmp(argv[i],"-excl"))
	{
	  if (++i>=argc) {help(); exit(4);} 
	  par.exclstr = new(char[strlen(argv[i])+1]);
	  strcpy(par.exclstr,argv[i]);
	}
      else if (!strcmp(argv[i],"-v") && (i<argc-1) && argv[i+1][0]!='-' ) v=atoi(argv[++i]);
      else if (!strcmp(argv[i],"-v"))  v=2;
      else if (!strcmp(argv[i],"-v0")) v=0;
      else if (!strcmp(argv[i],"-v1")) v=1;
      else if (!strcmp(argv[i],"-v2")) v=2;
      else if (!strcmp(argv[i],"-v3")) v=3;
      else if (!strcmp(argv[i],"-v4")) v=4;
      else if (!strcmp(argv[i],"-v5")) v=5;
      else if (!strcmp(argv[i],"-P") && (i<argc-1)) pself=atof(argv[++i]);
      else if (!strcmp(argv[i],"-p") && (i<argc-1)) par.p = atof(argv[++i]);
      else if (!strcmp(argv[i],"-e") && (i<argc-1)) par.E = atof(argv[++i]);
      else if (!strcmp(argv[i],"-E") && (i<argc-1)) par.E = atof(argv[++i]);
      else if (!strcmp(argv[i],"-b") && (i<argc-1)) par.b = atoi(argv[++i]);
      else if (!strcmp(argv[i],"-B") && (i<argc-1)) par.B = atoi(argv[++i]);
      else if (!strcmp(argv[i],"-z") && (i<argc-1)) par.z = atoi(argv[++i]);
      else if (!strcmp(argv[i],"-Z") && (i<argc-1)) par.Z = atoi(argv[++i]);
      else if (!strncmp(argv[i],"-nocons",7)) par.showcons=0;
      else if (!strncmp(argv[i],"-nopred",7)) par.showpred=0;
      else if (!strncmp(argv[i],"-nodssp",7)) par.showdssp=0;
      else if (!strncmp(argv[i],"-mark",7)) par.mark=1;
      else if (!strcmp(argv[i],"-seq") && (i<argc-1))  par.nseqdis=atoi(argv[++i]); 
      else if (!strcmp(argv[i],"-aliw") && (i<argc-1)) par.aliwidth=atoi(argv[++i]); 
      else if (!strcmp(argv[i],"-id") && (i<argc-1))   par.max_seqid=atoi(argv[++i]); 
      else if (!strcmp(argv[i],"-tct") && (i<argc-1))  probmin_tc=atoi(argv[++i]); 
      else if (!strcmp(argv[i],"-dwin") && (i<argc-1)) dotW=atoi(argv[++i]); 
      else if (!strcmp(argv[i],"-dsca") && (i<argc-1)) dotscale=atoi(argv[++i]); 
      else if (!strcmp(argv[i],"-dthr") && (i<argc-1)) dotthr=atof(argv[++i]); 
      else if (!strcmp(argv[i],"-dali") && (i<argc-1))  
	{
	  dotali=1; 
	  for (int index=0; index<256; index++) aliindices[index]=0;
	  while (i+1<argc && argv[i+1][0]!='-') // adds index to hash aliindices
	    {
	      i++;
	      if (strcmp(argv[i],"all")) aliindices[atoi(argv[i])]=1;	      
	      else dotali=2;
	    }
	}
      else if (!strcmp(argv[i],"-dmap"))  
	{
	  if (++i>=argc) 
	    {help(); cerr<<endl<<"Error in "<<program_name<<": no filename following -o\n"; exit(4);}
	  else 
	    {
	      dmapfile = new(char[strlen(argv[i])+1]);
	      strcpy(dmapfile,argv[i]);
	    }
	}
      else if (!strcmp(argv[i],"-dsat") && (i<argc-1)) dotsat=atof(argv[++i]); 
      else if (!strcmp(argv[i],"-qid") && (i<argc-1))  par.qid=atoi(argv[++i]); 
      else if (!strcmp(argv[i],"-qsc") && (i<argc-1))  par.qsc=atof(argv[++i]); 
      else if (!strcmp(argv[i],"-cov") && (i<argc-1))  par.coverage=atoi(argv[++i]); 
      else if (!strcmp(argv[i],"-diff") && (i<argc-1)) par.Ndiff=atoi(argv[++i]); 
      else if (!strcmp(argv[i],"-qscc") && (i<argc-1))    par.qsc_core=atof(argv[++i]); 
      else if (!strcmp(argv[i],"-csc") && (i<argc-1))     par.coresc=atof(argv[++i]); 
      else if (!strcmp(argv[i],"-Gonnet")) par.matrix=0; 
      else if (!strcmp(argv[i],"-HSDM")) par.matrix=1; 
      else if (!strcmp(argv[i],"-BLOSUM50")) par.matrix=2; 
      else if (!strcmp(argv[i],"-Blosum50")) par.matrix=2; 
      else if (!strcmp(argv[i],"-B50")) par.matrix=2; 
      else if (!strcmp(argv[i],"-BLOSUM62")) par.matrix=3; 
      else if (!strcmp(argv[i],"-Blosum62")) par.matrix=3; 
      else if (!strcmp(argv[i],"-B62")) par.matrix=3; 
      else if (!strcmp(argv[i],"-pcm") && (i<argc-1)) par.pcm=atoi(argv[++i]); 
      else if (!strcmp(argv[i],"-pca") && (i<argc-1)) par.pca=atof(argv[++i]); 
      else if (!strcmp(argv[i],"-pcb") && (i<argc-1)) par.pcb=atof(argv[++i]); 
      else if (!strcmp(argv[i],"-pcc") && (i<argc-1)) par.pcc=atof(argv[++i]);  
      else if (!strcmp(argv[i],"-pcw") && (i<argc-1)) par.pcw=atof(argv[++i]); 
      else if (!strcmp(argv[i],"-gapb") && (i<argc-1)) { par.gapb=atof(argv[++i]); if (par.gapb<=0.01) par.gapb=0.01;} 
      else if (!strcmp(argv[i],"-gapd") && (i<argc-1)) par.gapd=atof(argv[++i]); 
      else if (!strcmp(argv[i],"-gape") && (i<argc-1)) par.gape=atof(argv[++i]); 
      else if (!strcmp(argv[i],"-gapf") && (i<argc-1)) par.gapf=atof(argv[++i]); 
      else if (!strcmp(argv[i],"-gapg") && (i<argc-1)) par.gapg=atof(argv[++i]); 
      else if (!strcmp(argv[i],"-gaph") && (i<argc-1)) par.gaph=atof(argv[++i]); 
      else if (!strcmp(argv[i],"-gapi") && (i<argc-1)) par.gapi=atof(argv[++i]); 
      else if (!strcmp(argv[i],"-egq") && (i<argc-1)) par.egq=atof(argv[++i]); 
      else if (!strcmp(argv[i],"-egt") && (i<argc-1)) par.egt=atof(argv[++i]); 
      else if (!strcmp(argv[i],"-ssgap")) par.ssgap=1;
      else if (!strcmp(argv[i],"-ssgapd") && (i<argc-1)) par.ssgapd=atof(argv[++i]); 
      else if (!strcmp(argv[i],"-ssgape") && (i<argc-1)) par.ssgape=atof(argv[++i]); 
      else if (!strcmp(argv[i],"-ssgapi") && (i<argc-1)) par.ssgapi=atoi(argv[++i]); 
      else if (!strcmp(argv[i],"-ssm") && (i<argc-1)) par.ssm=atoi(argv[++i]); 
      else if (!strcmp(argv[i],"-ssw") && (i<argc-1)) par.ssw=atof(argv[++i]); 
      else if (!strcmp(argv[i],"-ssa") && (i<argc-1)) par.ssa=atof(argv[++i]); 
      else if (!strncmp(argv[i],"-glo",3)) {par.loc=0; if (par.mact>0.3 && par.mact<0.301) {par.mact=0;} }
      else if (!strncmp(argv[i],"-loc",3)) par.loc=1;
      else if (!strncmp(argv[i],"-alt",4) && (i<argc-1)) par.altali=atoi(argv[++i]); 
      else if (!strcmp(argv[i],"-map") || !strcmp(argv[i],"-MAP")) par.forward=2; 
      else if (!strcmp(argv[i],"-mac") || !strcmp(argv[i],"-MAC")) par.forward=2; 
      else if (!strcmp(argv[i],"-vit")) par.forward=0; 
      else if (!strcmp(argv[i],"-sto") && (i<argc-1))  {Nstochali=atoi(argv[++i]); par.forward=1;}
      else if (!strcmp(argv[i],"-r")) par.repmode=1; 
      else if (!strcmp(argv[i],"-M") && (i<argc-1)) 
	if (!strcmp(argv[++i],"a2m") || !strcmp(argv[i],"a3m"))  par.M=1; 
	else if(!strcmp(argv[i],"first"))  par.M=3; 
	else if (argv[i][0]>='0' && argv[i][0]<='9') {par.Mgaps=atoi(argv[i]); par.M=2;}
	else cerr<<endl<<"WARNING: Ignoring unknown argument: -M "<<argv[i]<<"\n";
      else if (!strcmp(argv[i],"-shift") && (i<argc-1)) par.shift=atof(argv[++i]); 
      else if (!strcmp(argv[i],"-mact") && (i<argc-1)) {par.mact=atof(argv[++i]); par.forward=2;}
      else if (!strcmp(argv[i],"-wstruc") && (i<argc-1)) par.wstruc=atof(argv[++i]); 
      else if (!strcmp(argv[i],"-opt") && (i<argc-1)) par.opt=atoi(argv[++i]); 
      else if (!strcmp(argv[i],"-sc") && (i<argc-1)) par.columnscore=atoi(argv[++i]); 
      else if (!strcmp(argv[i],"-corr") && (i<argc-1)) par.corr=atof(argv[++i]); 
      else if (!strcmp(argv[i],"-def")) par.readdefaultsfile=1; 
      else if (!strcmp(argv[i],"-ovlp") && (i<argc-1)) par.min_overlap=atoi(argv[++i]);
      else if (!strcmp(argv[i],"-tags")) par.notags=0;
      else if (!strcmp(argv[i],"-notags")) par.notags=1;
      else cerr<<endl<<"WARNING: Ignoring unknown option "<<argv[i]<<" ...\n";
      if (v>=4) cout<<i<<"  "<<argv[i]<<endl; //PRINT
    } // end of for-loop for command line input
}
#endif




/////////////////////////////////////////////////////////////////////////////
//// MAIN PROGRAM
/////////////////////////////////////////////////////////////////////////////
/**
 *
 * hhalign()
 *
 * @param[in,out] ppcFirstProf
 * first profile to be aligned
 * @param[in] iFirstCnt
 * number of sequences in 1st profile
 * @param[in,out] ppcSecndProf
 * second profile to be aligned
 * @param[in] iSecndCnt
 * number of sequences in 2nd profile
 * @param[out] dScore_p
 * score of alignment
 * @param[in] prHMM
 * HMM info of external HMM (background)
 * @param[in] pcPrealigned1
 * (gapped) 1st sequence aligned to background HMM
 * @param[in] pcRepresent1
 * (gapped) sequence representing background HMM aligned to 1st sequence
 * @param[in] pcPrealigned2
 * (gapped) 2nd sequence aligned to background HMM
 * @param[in] pcRepresent2
 * (gapped) sequence representing background HMM aligned to 2nd sequence
 * @param[in] rHhalignPara,
 * various parameters passed down from commandline
 * e.g., iMaxRamMB
 * @param[out] rHHscores
 * qualify goodness of alignment
 *
 * iFlag,zcAux,zcError are debugging arguments
 *
 * @return Non-zero on error
 */
extern "C" int 
hhalign(char **ppcFirstProf, int iFirstCnt, double *pdWeightsL,
        char **ppcSecndProf, int iSecndCnt, double *pdWeightsR,
        double *dScore_p, hmm_light *prHMM, hmm_light *prHMM2,
        char *pcPrealigned1, char *pcRepresent1, 
        char *pcPrealigned2, char *pcRepresent2, 
        hhalign_para rHhalignPara, hhalign_scores *rHHscores, 
        int iFlag, int iVerbosity, 
        char zcAux[], char zcError[]) {


#ifdef CLUSTALO
    int iRetVal = RETURN_OK;
    iAux_GLOBAL = iFlag;
#ifndef CLUSTALO_NOFILE
    int argc = 0;
    char **argv = NULL;
    argv = (char **)malloc(argc*sizeof(char *));
    for (int i = 0; i < argc; i++){
        argv[i] = (char *)malloc(100);
    }
    strcpy(argv[0], "./hhalign");
    strcpy(argv[1], "-t");
    strcpy(argv[2], "hhalign.C");
    strcpy(argv[3], "-i");
    strcpy(argv[4], "hhalign.C");
    strcpy(argv[5], "-ofas");
    strcpy(argv[6], "out");
#endif
#endif
    /*char* argv_conf[MAXOPT];*/ /* Input arguments from .hhdefaults file 
      (first=1: argv_conf[0] is not used) */
    /*int argc_conf;*/           // Number of arguments in argv_conf 
    /*char inext[IDLEN]="";*/    // Extension of query input file (hhm or a3m) 
    /*char text[IDLEN]="";*/     // Extension of template input file (hhm or a3m)
    /*int** ali=NULL;*/          // ali[i][j]=1 if (i,j) is part of an alignment
    /*int** alisto=NULL;*/       // ali[i][j]=1 if (i,j) is part of an alignment
    int Nali;                    /* number of normally backtraced alignments 
                                    in dot plot */
    
    SetDefaults(rHhalignPara); 
    strcpy(par.tfile,""); /* FS, Initialise Argument */
    strcpy(par.alnfile,""); 
    par.p=0.0 ;                  /* minimum threshold for inclusion in hit list 
                                    and alignment listing */
    par.E=1e6;                   /* maximum threshold for inclusion in hit list 
                                    and alignment listing */
    par.b=1;                     // min number of alignments
    par.B=100;                   // max number of alignments
    par.z=1;                     // min number of lines in hit list
    par.Z=100;                   // max number of lines in hit list
    par.append=0;                /* append alignment to output file 
                                    with -a option */
    par.altali=1;                /* find only ONE (possibly overlapping) 
                                    subalignment  */
    par.hitrank=0;               /* rank of hit to be printed as a3m alignment 
                                    (default=0) */
    par.outformat=3;             // default output format for alignment is a3m
    hit.self=0;                  // no self-alignment
    par.forward=2;               /* 0: Viterbi algorithm; 
                                    1: Viterbi+stochastic sampling; 
                                    2:Maximum Accuracy (MAC) algorithm */
    
    // Make command line input globally available
#ifndef CLUSTALO_NOFILE
    par.argv=argv; 
    par.argc=argc;
    RemovePathAndExtension(program_name,argv[0]); 
#endif
    
#ifndef CLUSTALO_NOFILE
    /* Enable changing verbose mode before defaults file 
       and command line are processed */
    for (int i=1; i<argc; i++)
        { 
            if (!strcmp(argv[i],"-def")) par.readdefaultsfile=1;
            else if (argc>1 && !strcmp(argv[i],"-v0")) v=0;
            else if (argc>1 && !strcmp(argv[i],"-v1")) v=1;
            else if (argc>2 && !strcmp(argv[i],"-v")) v=atoi(argv[i+1]);
        }
    
    // Read .hhdefaults file?
    if (par.readdefaultsfile) 
        {
            // Process default otpions from .hhconfig file
            ReadDefaultsFile(argc_conf,argv_conf);
            ProcessArguments(argc_conf,argv_conf);
        }
#endif
    /* Process command line options 
       (they override defaults from .hhdefaults file) */
#ifndef CLUSTALO_NOFILE
    ProcessArguments(argc,argv);
#endif
    
    
#ifdef CLUSTALO
    int iAuxLen1 = strlen(ppcFirstProf[0]);
    int iAuxLen2 = strlen(ppcSecndProf[0]); 
    if ( (0 == iAuxLen1) || (0 == iAuxLen2) ){ /* problem with empty profiles, FS, r249 -> r250 */
        sprintf(zcError, "%s:%s:%d: strlen(prof1)=%d, strlen(prof2)=%d -- Nothing to align!\n",
                __FUNCTION__, __FILE__, __LINE__, iAuxLen1, iAuxLen2);
        iRetVal = RETURN_UNKNOWN;
        /* Note: at this stage cannot do   'goto this_is_the_end;' 
           because would cross initialisation of several variables   */
        return iRetVal;
    }
    par.maxResLen  = iAuxLen1 > iAuxLen2 ? iAuxLen1 : iAuxLen2;
    const int ciGoodMeasureSeq = 10;
    par.maxResLen += ciGoodMeasureSeq;
    par.maxColCnt  = iAuxLen1 + iAuxLen2 + ciGoodMeasureSeq;
    //par.max_seqid=iFirstCnt+iSecndCnt+3; /* -id     */
    par.max_seqid=DEFAULT_FILTER;        /* -id     */
    par.loc=0; par.mact=0;              /* -glob   */
    par.nseqdis=iFirstCnt+iSecndCnt;   /* -seq    */
    par.showcons=0;                   /* -nocons */
    par.showdssp=0;                  /* -nodssp */
    par.Mgaps=100;                  /* -M      */
    par.M=2;                       /* -M      */
    par.pdWg1=pdWeightsL;         /* tree wg */
    par.pdWg2=pdWeightsR;        /* tree wg */
    v = 0;                      /* -v0     */
    /* NOTE: *qali declared globally but only created here,
       pass in number of sequences to get rid of statically 
       defined MAXSEQ (FS)
    */
    Alignment qali(iFirstCnt+iSecndCnt);
    HMM q(iFirstCnt+iSecndCnt);
    HMM t(iFirstCnt+iSecndCnt);
#endif
    
    
#ifndef CLUSTALO_NOFILE
    // Check command line input and default values
    if (!*par.infile) 
        {help(); cerr<<endl<<"Error in "<<program_name<<": no query alignment file given (-i file)\n"; exit(4);}
    if (par.forward==2 && par.loc==0) 
        {
            if (par.mact<0.301 || par.mact>0.300) 
                if (v>=1) fprintf(stderr,"REMARK: in -mac -global mode -mact is forced to 0\n");
            par.mact=0;
            //       par.loc=1; // global forward-backward algorithm seems to work fine! use it in place of local version.
        }
    
    // Get rootname (no directory path, no extension) and extension of infile
    RemoveExtension(q.file,par.infile);
    RemoveExtension(t.file,par.tfile); /* FS, NOFILE2 (commented out) */
    Extension(inext,par.infile); 
    Extension(text,par.tfile); /*FS, NOFILE2 (commented out) */
    
    // Check option compatibilities
    if (par.nseqdis>MAXSEQDIS-3-par.showcons) par.nseqdis=MAXSEQDIS-3-par.showcons; //3 reserved for secondary structure
    if (par.aliwidth<20) par.aliwidth=20; 
    if (par.pca<0.001) par.pca=0.001; // to avoid log(0)
    if (par.b>par.B) par.B=par.b;
    if (par.z>par.Z) par.Z=par.z;
    if (par.hitrank>0) par.altali=0;
    
    // Input parameters
    if (v>=3) 
        {
            cout<<"query file : "<<par.infile<<"\n";
            cout<<"template file: "<<par.tfile<<"\n";/*FS, NOFILE2 (commented out) */
            cout<<"Output file:  "<<par.outfile<<"\n";
            cout<<"Alignment file:  "<<par.alnfile<<"\n";
        }
#endif /* NOFILE2 */
    
    
    // Set (global variable) substitution matrix and derived matrices
	// DD: new experimental matrices and params for nucleotides
	if(rHhalignPara.bIsDna)
	{
		nucleomode = true;
		SetDnaDefaults(rHhalignPara);
		SetDnaSubstitutionMatrix();
	}
	else if(rHhalignPara.bIsRna)
	{
		nucleomode = true;
		SetRnaDefaults(rHhalignPara);
		SetRnaSubstitutionMatrix();
	}	
	else
		SetSubstitutionMatrix();
    
    // Set secondary structure substitution matrix
    SetSecStrucSubstitutionMatrix();


    /* moved Viterbi switch from after RnP() to here, 
       switch after RnP() ineffectual as RnP decides log/lin of transition, 
       however log/lin of transitions depends on MAC/Viterbi,
       FS, r228 -> r229 */
    int qL, tL; 
    if (iFirstCnt > 0) {
        qL = strlen(ppcFirstProf[0]);
    }
    else {
        qL = prHMM->L;
    }
    if (iSecndCnt > 0) {
        tL = strlen(ppcSecndProf[0]);
    }
    else {
        tL = prHMM2->L;
    }
    

    //  const float MEMSPACE_DYNPROG = 512*1024*1024;
    /* determine amount of memory available for MAC on command-line; FS, r240 -> r241 */
    const float MEMSPACE_DYNPROG = (double)1024*1024*rHhalignPara.iMacRamMB;
    // longest allowable length of database HMM
    int Lmaxmem=(int)((float)MEMSPACE_DYNPROG/qL/6/8); 
    if (par.forward==2 && tL+2>=Lmaxmem) {
        sprintf(zcError, "%s:%s:%d: switch to Viterbi (qL=%d, tL=%d, MAC-RAM=%d\n)", 
                __FUNCTION__, __FILE__, __LINE__, qL, tL, rHhalignPara.iMacRamMB);
        if (v>=1)
            cerr<<"WARNING: Not sufficient memory to realign with MAC algorithm. Using Viterbi algorithm."<<endl;
        par.forward=0;

        /* use different 'fudge' parameters for Viterbi, FS, r228 -> r229 */
        par.pca = par.pcaV;
        par.pcb = par.pcbV;
        par.pcc = par.pccV;
        par.pcw = par.pcwV;
        
        par.gapb = par.gapbV;
        par.gapd = par.gapdV;
        par.gape = par.gapeV;
        par.gapf = par.gapfV;
        par.gapg = par.gapgV;
        par.gaph = par.gaphV;
        par.gapi = par.gapiV;

    }
    
    // Read input file (HMM, HHM, or alignment format), and add pseudocounts etc.
    q.cQT = 'q';
    if (OK != ReadAndPrepare(INTERN_ALN_2_HMM,
                             ppcFirstProf, iFirstCnt, prHMM, 
                             pcPrealigned1, pcRepresent1, pdWeightsL, 
                             (char*)(""), q, &qali)) {
        sprintf(zcError, "%s:%s:%d: Problem Reading/Preparing q-profile (len=%d)\n",
                __FUNCTION__, __FILE__, __LINE__, qL);
        iRetVal = RETURN_FROM_RNP;
        goto this_is_the_end;
    }
    // Set query columns in His-tags etc to Null model distribution
    if (par.notags) q.NeutralizeTags();
    
    // Do self-comparison?
    if (0 /* !*par.tfile // FS, 2010-03-10 */) 
        {
            // Deep-copy q into t
            t = q;
            
            // Find overlapping alternative alignments
            hit.self=1;
        } 
    // Read template alignment/HMM t and add pseudocounts
    else 
        {
            char infile[] = "";
            /* Read input file (HMM, HHM, or alignment format), 
               and add pseudocounts etc. */
            t.cQT = 't';
            if (OK != ReadAndPrepare(INTERN_ALN_2_HMM,
                                     ppcSecndProf, iSecndCnt, prHMM2, 
                                     pcPrealigned2, pcRepresent2, pdWeightsR, 
                                     infile, t)) {
                sprintf(zcError, "%s:%s:%d: Problem Reading/Preparing t-profile (len=%d)\n",
                        __FUNCTION__, __FILE__, __LINE__, tL);
                iRetVal = RETURN_FROM_RNP;
                goto this_is_the_end;
            }
        }
    
  // Factor Null model into HMM t
  t.IncludeNullModelInHMM(q,t); 

  /* alignment will fail if one profile does not contain useful characters, FS, r259 -> r260 */
  if ( (q.L <= 0) || (t.L <= 0) ){
      sprintf(zcError, "%s:%s:%d: Problem Reading/Preparing profiles (len(q)=%d/len(t)=%d)\n",
              __FUNCTION__, __FILE__, __LINE__, q.L, t.L);
      iRetVal = RETURN_FROM_RNP;
      goto this_is_the_end;
  }

  /* switch at this stage is ineffectual as log/lin already decided in RnP(). 
     FS, r228 -> r229 */
  /*const float MEMSPACE_DYNPROG = 512*1024*1024;
  // longest allowable length of database HMM
  int Lmaxmem=(int)((float)MEMSPACE_DYNPROG/q.L/6/8); 
  if (par.forward==2 && t.L+2>=Lmaxmem) 
  {
  if (v>=1)
  cerr<<"WARNING: Not sufficient memory to realign with MAC algorithm. Using Viterbi algorithm."<<endl;
  par.forward=0;
  }*/

  // Allocate memory for dynamic programming matrix
  hit.AllocateBacktraceMatrix(q.L+2,t.L+2); // ...with a separate dynamic programming matrix (memory!!)
  if (par.forward>=1 || Nstochali) 
    hit.AllocateForwardMatrix(q.L+2,t.L+2);
  if (par.forward==2)     
    hit.AllocateBackwardMatrix(q.L+2,t.L+2);

#ifndef CLUSTALO_NOFILE
  // Read structure file for Forward() function?
  if (strucfile && par.wstruc>0) 
    {
      float PMIN=1E-20;
      Pstruc = new(float*[q.L+2]);
      for (int i=0; i<q.L+2; i++) Pstruc[i] = new(float[t.L+2]);
      Sstruc = new(float*[q.L+2]);
      for (int i=0; i<q.L+2; i++) Sstruc[i] = new(float[t.L+2]);
      FILE* strucf=NULL;
      if (strcmp(strucfile,"stdin"))
	{
	  strucf = fopen(strucfile, "r");
	  if (!strucf) OpenFileError(strucfile);
	}
      else 
	{
	  strucf = stdin;
	  if (v>=2) printf("Reading structure matrix from standard input ... (for UNIX use ^D for 'end-of-file')\n");
	}
      for (int i=1; i<=q.L; i++)
	{
	  for (int j=1; j<=t.L; j++)
	    {
	      float f;
	      if (fscanf(strucf,"%f",&f) <=0 )
	      {
		fprintf(stderr,"Error: too few numbers in file %s while reading line %i, column %i\n",strucfile,i,j);
		exit(1);
	      } 
	      if (par.wstruc==1)
		Pstruc[i][j]=fmax(f,PMIN);
	      else 
		Pstruc[i][j]=fmax(pow(f,par.wstruc),PMIN);
// 	      printf("%10.2E ",f);
	      Sstruc[i][j] = par.wstruc * log2(f);
	    }
// 	  printf("\n");
	}
      fclose(strucf);
    } /* (strucfile && par.wstruc>0) */
#endif
  
  /* Do (self-)comparison, store results if score>SMIN, 
     and try next best alignment */
  /* FIXME very ambigous and possibly faulty if-else */
  if (v>=2)
  {
    if (par.forward==2) {
        printf("Using maximum accuracy (MAC) alignment algorithm ...\n");
    }
    else if (par.forward==0) {
        printf("Using Viterbi algorithm ...\n");
    }
    else if (par.forward==1) {
        printf("Using stochastic sampling algorithm ...\n");
    }
    else {
        printf("\nWhat alignment algorithm are we using??\n");
    }
  }


  for (hit.irep=1; hit.irep<=imax(par.hitrank,par.altali); hit.irep++)
      {
          if (par.forward==0)        
              {
                  // generate Viterbi alignment
                  hit.Viterbi(q,t,Sstruc);
                  hit.Backtrace(q,t);
              } 
          else if (par.forward==1)   
              {
                  // generate a single stochastically sampled alignment
                  hit.Forward(q,t,Pstruc); 
                  srand( time(NULL) );   // initialize random generator 
                  hit.StochasticBacktrace(q,t);
                  hitlist.Push(hit);      /* insert hit at beginning of list 
                                             (last repeats first!) */
                  (hit.irep)++;
                  break;
              }
          else if (par.forward==2)   
              {
                  // generate forward alignment
                  if (OK != hit.Forward(q,t,Pstruc)){
                      fprintf(stderr, "%s:%s:%d: cannot complete hit.Forward\n", 
                              __FUNCTION__, __FILE__, __LINE__);
                      iRetVal = RETURN_FROM_MAC; /* spot double overflow in Forward(). FS, r241 -> r243 */
                      goto this_is_the_end;
                  }
                  if (OK != hit.Backward(q,t)){
                      fprintf(stderr, "%s:%s:%d: cannot complete hit.backward\n", 
                              __FUNCTION__, __FILE__, __LINE__);
                      iRetVal = RETURN_FROM_MAC; /* spot double overflow in hit.Backward(). FS, r241 -> r243 */
                      goto this_is_the_end;
                  }
                  hit.MACAlignment(q,t);
                  if ((isnan(hit.score)) || (isnan(hit.Pforward))){
                      printf("nan after MAC\n");
                  }
                  hit.BacktraceMAC(q,t);
                  if ((isnan(hit.score)) || (isnan(hit.Pforward))){
                      printf("nan after backtrace\n");
                  } 
              } /* use MAC algorithm */
          //       printf ("%-12.12s  %-12.12s   irep=%-2i  score=%6.2f hit.Pvalt=%.2g\n",hit.name,hit.fam,hit.irep,hit.score,hit.Pvalt);
          *dScore_p = hit.score;
          
          if (hit.irep<=par.hitrank || hit.score>SMIN || (hit.Pvalt<pself && hit.score>0 )) 
              hitlist.Push(hit);      // insert hit at beginning of list (last repeats first!) and do next alignment
          else 
              {
                  if (hit.irep==1) hitlist.Push(hit); // first hit will be inserted into hitlist anyway, even if not significant
                  hit = hitlist.ReadLast(); // last alignment was not significant => read last (significant) hit from list
                  break;
                  /* FIXME: memory leak during normal alignment: both push'es above are not freed:
                   * valgrind --leak-check=full --show-reachable=yes ./src/clustalo -i ~/void/protein.fa -o /dev/null -v
                   ==9506== 1,456 bytes in 1 blocks are still reachable in loss record 4 of 5
                   ==9506==    at 0x4C2726C: operator new(unsigned long) (vg_replace_malloc.c:230)
                   ==9506==    by 0x4618C5: List<Hit>::Push(Hit) (list-C.h:134)
                   ==9506==    by 0x45E87A: hhalign (hhalign.cpp:914)
                   ==9506==    by 0x413F8C: HhalignWrapper (hhalign_wrapper.c:405)
                   ==9506==    by 0x41A7B5: MyMain (mymain.c:676)
                   ==9506==    by 0x4031A6: main (main.cpp:37)
                   ==9506== 
                   ==9506== 
                   ==9506== 4,442 (4,368 direct, 74 indirect) bytes in 3 blocks are definitely lost in loss record 5 of 5
                   ==9506==    at 0x4C2726C: operator new(unsigned long) (vg_replace_malloc.c:230)
                   ==9506==    by 0x4618C5: List<Hit>::Push(Hit) (list-C.h:134)
                   ==9506==    by 0x45E7C0: hhalign (hhalign.cpp:911)
                   ==9506==    by 0x413F8C: HhalignWrapper (hhalign_wrapper.c:405)
                   ==9506==    by 0x41A7B5: MyMain (mymain.c:676)
                   ==9506==    by 0x4031A6: main (main.cpp:37)
                  */
              }
      } /* 1 <= hit.irep <= imax(par.hitrank,par.altali) */
  Nali = hit.irep;
  
  
#ifndef CLUSTALO_NOFILE
  // Write posterior probability matrix as TCoffee library file
  if (tcfile) 
      {
          if (v>=2) printf("Writing TCoffee library file to %s\n",tcfile);
          int i,j; 
          FILE* tcf=NULL;
          if (strcmp(tcfile,"stdout")) tcf = fopen(tcfile, "w"); else tcf = stdout;
          if (!tcf) OpenFileError(tcfile);
          fprintf(tcf,"! TC_LIB_FORMAT_01\n");
          fprintf(tcf,"%i\n",2); // two sequences in library file
          fprintf(tcf,"%s %i %s\n",q.name,q.L,q.seq[q.nfirst]+1);
          fprintf(tcf,"%s %i %s\n",hit.name,hit.L,hit.seq[hit.nfirst]+1);
          fprintf(tcf,"#1 2\n");
          for (i=1; i<=q.L; i++)  // print all pairs (i,j) with probability above PROBTCMIN
              for (j=1; j<=t.L; j++)
                  if (hit.B_MM[i][j]>probmin_tc) 
                      fprintf(tcf,"%5i %5i %5i\n",i,j,iround(100.0*hit.B_MM[i][j]));
          for (int step=hit.nsteps; step>=1; step--)  // print all pairs on MAC alignment which were not yet printed
              {
                  i=hit.i[step]; j=hit.j[step];
                  // 	  printf("%5i %5i %5i  %i\n",i,j,iround(100.0*hit.B_MM[i][j]),hit.states[step]);
                  if (hit.states[step]>=MM && hit.B_MM[i][j]<=probmin_tc) 
                      fprintf(tcf,"%5i %5i %5i\n",i,j,iround(100.0*hit.B_MM[i][j]));
              }
          
          
          fprintf(tcf,"! SEQ_1_TO_N\n");
          fclose(tcf);
          //       for (i=1; i<=q.L; i++)
          //        	{
          //        	  double sum=0.0;
          //        	  for (j=1; j<=t.L; j++) sum+=hit.B_MM[i][j];
          // 	  printf("i=%-3i sum=%7.4f\n",i,sum);
          //        	}
          //        printf("\n");
      } /* if (tcfile) */

  // Write last alignment into alitabfile
  if (alitabfile) 
      {
          FILE* alitabf=NULL;
          if (strcmp(alitabfile,"stdout")) alitabf = fopen(alitabfile, "w"); else alitabf = stdout;
          if (!alitabf) OpenFileError(alitabfile);
          if (par.forward==2) 
              {
                  fprintf(alitabf,"    i     j  score     SS  probab\n");
                  for (int step=hit.nsteps; step>=1; step--)
                      if (hit.states[step]>=MM) 
                          fprintf(alitabf,"%5i %5i %6.2f %6.2f %7.4f\n",hit.i[step],hit.j[step],hit.S[step],hit.S_ss[step],hit.P_posterior[step]);
              } 
          else 
              {
                  fprintf(alitabf,"    i     j  score     SS\n");
                  for (int step=hit.nsteps; step>=1; step--)
                      if (hit.states[step]>=MM) 
                          fprintf(alitabf,"%5i %5i %6.2f %6.2f\n",hit.i[step],hit.j[step],hit.S[step],hit.S_ss[step]);
              }
          fclose(alitabf);
      } /* if (alitabfile) */
#endif
  
  // Do Stochastic backtracing?
  if (par.forward==1)
      for (int i=1; i<Nstochali; i++) 
          {
              hit.StochasticBacktrace(q,t);
              hitlist.Push(hit);      //insert hit at beginning of list (last repeats first!)
              (hit.irep)++;
          }
  else // Set P-value, E-value and probability
      {
          if (q.mu) hitlist.GetPvalsFromCalibration(q);
          else if (t.mu) hitlist.GetPvalsFromCalibration(t);
      }
  
  // Print FASTA or A2M alignments?
#ifndef CLUSTALO_NOFILE
  if (*par.pairwisealisfile) 
#endif
      {
          if (v>=2) {
              cout<<"Printing alignments in " <<
                  (par.outformat==1? "FASTA" : par.outformat==2?"A2M" :"A3M") <<
                  " format to "<<par.pairwisealisfile<<"\n"; 
          }
          int iPrAliRtn = hitlist.PrintAlignments(
#ifdef CLUSTALO
                                                  ppcFirstProf, ppcSecndProf, zcAux, zcError, 
#endif
                                                  q, par.pairwisealisfile, par.outformat);
          if ( (OK != iPrAliRtn) ){
              sprintf(zcAux, "%s:%s:%d: Could not print alignments\n",
                      __FUNCTION__, __FILE__, __LINE__);
              strcat(zcError, zcAux);
              iRetVal = RETURN_FROM_PRINT_ALI; /* this is where mis-alignment was originally spotted, 
                                    hope to trap it now earlier. FS, r241 -> r243 */
              goto this_is_the_end;
          }
      } /* if (*par.pairwisealisfile) */
  
#ifndef CLUSTALO_NOFILE
  // Print hit list and  alignments
  if (*par.outfile) 
      {
          hitlist.PrintHitList(q, par.outfile); 
          hitlist.PrintAlignments(
#ifdef CLUSTALO
                                  ppcFirstProf, ppcSecndProf, zcAux, zcError, 
#endif
                                  q, par.outfile);
          if (v>=2) {
              WriteToScreen(par.outfile,1000); //write only hit list to screen
          }
      }
#endif
  

  ///////////////////////////////////////////////////////////////////////
#ifndef CLUSTALO_NOFILE
  // Show results for hit with rank par.hitrank
  if (par.hitrank==0) hit=hitlist.Read(1); else hit=hitlist.Read(par.hitrank);
  
  // Generate output alignment or HMM file?
  if (*par.alnfile || *par.psifile || *par.hhmfile) 
      {
          if (par.append==0) 
              {
                  if (v>=2 && *par.alnfile) printf("Merging template to query alignment and writing resulting alignment in A3M format to %s...\n",par.alnfile);
                  if (v>=2 && *par.psifile) printf("Merging template to query alignment and writing resulting alignment in PSI format to %s...\n",par.psifile);
              }
          else 
              {
                  if (v>=2 && *par.alnfile) printf("Merging template to query alignment and appending template alignment in A3M format to %s...\n",par.alnfile);
                  if (v>=2 && *par.psifile) printf("Merging template to query alignment and appending template alignment in PSI format to %s...\n",par.psifile);
              }
          
          // Read query alignment into Qali
          Alignment Qali;  // output A3M generated by merging A3M alignments for significant hits to the query alignment
          char qa3mfile[NAMELEN];
          RemoveExtension(qa3mfile,par.infile); // directory??
          strcat(qa3mfile,".a3m");
          FILE* qa3mf=fopen(qa3mfile,"r");
          if (!qa3mf) OpenFileError(qa3mfile);
          Qali.Read(qa3mf,qa3mfile);
          fclose(qa3mf);
          
          // Align query with template in master-slave mode 
          Qali.MergeMasterSlave(hit,par.tfile); /*FS, NOFILE2 (commented out) */
          
          // Write output A3M alignment?
          if (*par.alnfile) {
              Qali.WriteToFile(par.alnfile,"a3m");
          }
          
          if (*par.psifile) 
              {
                  /* Convert ASCII to int (0-20),throw out all insert states, 
                     record their number in I[k][i]  */
                  Qali.Compress("merged A3M file");
                  
                  // Write output PSI-BLAST-formatted alignment?
                  Qali.WriteToFile(par.psifile,"psi");
              }
      } /* if (*par.alnfile || *par.psifile || *par.hhmfile) */
#endif  /* NOFILE2 */
  //////////////////////////////////////////////////////////////////////////
  
  
  
  
  //   double log2Pvalue;
  //   if (par.ssm && (par.ssm1 || par.ssm2))
  //     {
  //       log2Pvalue=hit.logPval/0.693147181+0.45*(4.0*par.ssw/0.15-hit.score_ss);
  //       if (v>=2) 
  // 	printf("Aligned %s with %s:\nApproximate P-value INCLUDING SS SCORE = %7.2g\n",q.name,t.name,pow(2.0,log2Pvalue));
  //     } else {
  //       if (v>=2) 
  // 	printf("Aligned %s with %s:\nApproximate P-value (without SS score) = %7.2g\n",q.name,t.name,hit.Pval);
  //    }
  
  if (v>=2)
      {
          if (par.hitrank==0)
              printf("Aligned %s with %s: Score = %-7.2f  P-value = %-7.2g\n",
                     q.name,t.name,hit.score,hit.Pval);
          else
              printf("Aligned %s with %s (rank %i): Score = %-7.2f  P-value = %-7.2g\n",
                     q.name,t.name,par.hitrank,hit.score,hit.Pval);
      }
  
  rHHscores->hhScore     = hit.score; 
  rHHscores->forwardProb = hit.Pforward;
  rHHscores->sumPP       = hit.sum_of_probs;

  /* next few lines commented out as they caused segfaults in RNA mode */
  /*  rHHscores->PP = (double *)malloc((hit.L+GOOD_MEASURE) * sizeof(double));*/
  rHHscores->L  = hit.L;
/*  for (int i = 0; i < hit.L; i++){
      cout << "rHHscores->PP[i]" << rHHscores->PP[i] << endl;
      cout << "hit.P_posterior[i+1]" << hit.P_posterior[i+1] << endl;
      rHHscores->PP[i] = hit.P_posterior[i+1];
  } */


 this_is_the_end: 

  // Delete memory for dynamic programming matrix
  hit.DeleteBacktraceMatrix(q.L+2);
  if (par.forward>=1 || Nstochali) 
      hit.DeleteForwardMatrix(q.L+2);
  if (par.forward==2) 
      hit.DeleteBackwardMatrix(q.L+2);
  /*   if (Pstruc) { 
       for (int i=0; i<q.L+2; i++) delete[](Pstruc[i]); delete[](Pstruc);}
  */
  
  // Delete content of hits in hitlist
  hitlist.Reset();
  while (!hitlist.End()) 
      hitlist.ReadNext().Delete(); // Delete content of hit object

  if (strucfile && par.wstruc>0) 
      {
          for (int i=0; i<q.L+2; i++){
              delete[] Pstruc[i]; Pstruc[i] = NULL;
          }
          delete[] Pstruc; Pstruc = NULL;
          for (int i=0; i<q.L+2; i++){
              delete[] Sstruc[i]; Sstruc[i] = NULL;
          }
          delete[] Sstruc; Sstruc = NULL;
          delete[] strucfile; strucfile = NULL;
      }
  
  if (pngfile){
      delete[] pngfile; pngfile = NULL;
  }
  if (alitabfile){
      delete[] alitabfile; alitabfile = NULL;
  }
  if (tcfile){
      delete[] tcfile; tcfile = NULL;
  }
  if (par.exclstr){
      delete[] par.exclstr; par.exclstr = NULL;
  }
  
#ifndef CLUSTALO_NOFILE
  // Print 'Done!'
  FILE* outf=NULL;
  if (!strcmp(par.outfile,"stdout")) printf("Done!\n");
  else
      {
          if (*par.outfile)
              {
                  outf=fopen(par.outfile,"a"); //open for append
                  fprintf(outf,"Done!\n");
                  fclose(outf);
              }
          if (v>=2) printf("Done\n");
      }
#endif
  
  //qali.ClobberGlobal();
  hit.ClobberGlobal();
  if (iFirstCnt > 0){
      q.ClobberGlobal();
  }
  if (iSecndCnt > 0){
      t.ClobberGlobal();
  }
  hitlist.ClobberGlobal();

  return iRetVal;

} /* this is the end of hhalign() //end main */

//////////////////////////////////////////////////////////////////////////////
// END OF MAIN
//////////////////////////////////////////////////////////////////////////////

/*
 * EOF hhalign.C
 */

