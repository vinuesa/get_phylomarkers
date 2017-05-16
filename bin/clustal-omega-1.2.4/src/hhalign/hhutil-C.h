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
 * RCS $Id: hhutil-C.h 315 2016-12-15 17:18:30Z fabian $
 */

#define PARAMETERSFROMFILE 0

extern bool nucleomode;
 
//////////////////////////////////////////////////////////////////////////////
// Transform a character to lower case and '.' to '-' and vice versa
//////////////////////////////////////////////////////////////////////////////

inline char 
MatchChr(char c)  {return ((c>='a' && c<='z')? c-'a'+'A' : (c=='.'? '-':c) );}

inline char 
InsertChr(char c) {return ((c>='A' && c<='Z')? c+'a'-'A' : ((c>='0' && c<='9') || c=='-')? '.':c );}

inline int  
WordChr(char c) {return (int)((c>='A' && c<='Z') || (c>='a' && c<='z'));}


//////////////////////////////////////////////////////////////////////////////
/**
 * @brief Transforms the one-letter amino acid code into an integer between 0 and 22
 */
inline char 
aa2i(char c)
{
  if (c>='a' && c<='z') c+='A'-'a';
  if(!nucleomode)
  {	  
	  // proteins!
	  // A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V
	  switch (c)
		{
		case 'A': return 0;
		case 'R': return 1;
		case 'N': return 2;
		case 'D': return 3;
		case 'C': return 4;
		case 'Q': return 5;
		case 'E': return 6;
		case 'G': return 7;
		case 'H': return 8;
		case 'I': return 9;
		case 'L': return 10;
		case 'K': return 11;
		case 'M': return 12;
		case 'F': return 13;
		case 'P': return 14;
		case 'S': return 15;
		case 'T': return 16;
		case 'W': return 17;
		case 'Y': return 18;
		case 'V': return 19;
		case 'X': return ANY;
		case 'J': return ANY;
		case 'O': return ANY;
		case 'U': return 4;  //Selenocystein -> Cystein
		case 'B': return 3;  //D (or N)
		case 'Z': return 6;  //E (or Q)
		case '-': return GAP;
		case '.': return GAP;
		case '_': return GAP;
		}
  }
  else
  {
		// nucleotides!
	switch(c)
	{
    case 'A': return 0;
    case 'C': return 1;
    case 'G': return 2;
    case 'T': return 3;
    case 'U': return 4;
    case 'N': return ANY;
    case '-': return GAP;
    case '.': return GAP;
    case '_': return GAP;
	default: return ANY;
	}
  }
  if (c>=0 && c<=32) return -1; // white space and control characters
  return -2;
}

///////////////////////////////////////////////////////////////////////////////
/**
 * @brief Transforms integers between 0 and 22 into the one-letter amino acid code
 */
inline char 
i2aa(char c)
{
  //A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V
  if(!nucleomode)
  {
  switch (c)
    {
    case 0: return 'A';
    case 1: return 'R';
    case 2: return 'N';
    case 3: return 'D';
    case 4: return 'C';
    case 5: return 'Q';
    case 6: return 'E';
    case 7: return 'G';
    case 8: return 'H';
    case 9: return 'I';
    case 10: return 'L';
    case 11: return 'K';
    case 12: return 'M';
    case 13: return 'F';
    case 14: return 'P';
    case 15: return 'S';
    case 16: return 'T';
    case 17: return 'W';
    case 18: return 'Y';
    case 19: return 'V';
    case ANY: return 'X';
    case GAP: return '-';
    case ENDGAP: return '-';
    }
  }
  else
  {
	// DNA/RNA
	switch(c)
	{
		case 0: return 'A';
		case 1: return 'C';
		case 2: return 'G';
		case 3: return 'T';
		case 4: return 'U';
		case ANY: return 'N';
		case GAP: return '-';
		case ENDGAP: return '-';
	}
  }
  return '?';
}

//////////////////////////////////////////////////////////////////////////////
/**
 * @brief Transforms the dssp/psipred secondary structure code into an integer number
 */
inline char 
ss2i(char c)
{
  //- H E C S T G B
  if (c>='a' && c<='z') c+='A'-'a';
  switch (c)
    {
    case '.': return 0;
    case '-': return 0;
    case 'X': return 0;
    case 'H': return 1;
    case 'E': return 2;
    case 'C': return 3;
    case '~': return 3;
    case 'S': return 4;
    case 'T': return 5;
    case 'G': return 6;
    case 'B': return 7;
    case 'I': return 3;
    case ' ': return -1;
    case '\t': return -1;
    case '\n': return -1;
    }
  return -2;
}

//////////////////////////////////////////////////////////////////////////////
/**
 * @brief Transforms integers between 0 and 8 into the dssp/psipred secondary structure code
 */
inline char 
i2ss(int c)
{
  //- H E C S T G B
  switch (c)
    {
    case 0: return '-';
    case 1: return 'H';
    case 2: return 'E';
    case 3: return 'C';
    case 4: return 'S'; 
    case 5: return 'T';
    case 6: return 'G';
    case 7: return 'B';
    case 8: return 'I';
    }
  return '?';
}


//////////////////////////////////////////////////////////////////////////////
/**
 * @brief Transforms the solvend accessiblity code into an integer number
 */
inline char 
sa2i(char c)
{
  //- A B C D E 
  if (c>='a' && c<='z') c+='A'-'a';
  switch (c)
    {
    case '.': return 0;
    case '-': return 0;
    case 'A': return 1;
    case 'B': return 2;
    case 'C': return 3;
    case 'D': return 4;
    case 'E': return 5;
    case 'F': return 6;
    case ' ': return -1;
    case '\t': return -1;
    case '\n': return -1;
    }
  return -2;
}

//////////////////////////////////////////////////////////////////////////////
/**
 * @brief Transforms integers between 0 and 5 into the solvent accessibility code
 */
inline char i2sa(int c)
{
  //- H E C S T G B
  switch (c)
    {
    case 0: return '-';
    case 1: return 'A';
    case 2: return 'B';
    case 3: return 'C';
    case 4: return 'D'; 
    case 5: return 'E';
    case 6: return 'F';
    }
  return '?';
}


//////////////////////////////////////////////////////////////////////////////
/**
 * @brief Transforms alternative secondary structure symbols into symbols 
 */
inline char 
ss2ss(char c)
{
  //- H E C S T G B
  switch (c)
    {
    case '~': return 'C';
    case 'I': return 'C';
    case 'i': return 'c';      
    case 'H': 
    case 'E': 
    case 'C': 
    case 'S': 
    case 'T': 
    case 'G': 
    case 'B': 
    case 'h': 
    case 'e': 
    case 'c': 
    case 's': 
    case 't': 
    case 'g': 
    case 'b': 
    case '.': 
      return c;
    }
  return '-';
}

//////////////////////////////////////////////////////////////////////////////
/**
 * @brief Transforms confidence values of psipred into internal code
 */
inline char 
cf2i(char c)
{
  switch (c)
    {
    case '-': return 0;
    case '.': return 0;
    case '0': return 1;
    case '1': return 2;
    case '2': return 3;
    case '3': return 4;
    case '4': return 5;
    case '5': return 6;
    case '6': return 7;
    case '7': return 8;
    case '8': return 9;
    case '9': return 10;
    }
  return 0;
}

//////////////////////////////////////////////////////////////////////////////
/**
 * @brief Transforms internal representation of psipred confidence values into printable chars
 */
inline char 
i2cf(char c)
{
  switch (c)
    {
    case 0: return '-';
    case 1: return '0';
    case 2: return '1';
    case 3: return '2';
    case 4: return '3';
    case 5: return '4';
    case 6: return '5';
    case 7: return '6';
    case 8: return '7';
    case 9: return '8';
    case 10: return '9';
    }
  return '-';
}


//////////////////////////////////////////////////////////////////////////////
/**
 * @brief Fast lookup of log2(1+2^(-x)) for x>=0 (precision < 0.35%)
 */
inline float 
fast_addscore(float x) 
{
  static float val[2001];         // val[i]=log2(1+2^(-x))
  static char initialized;
  if (x>20) return 0.0;
  if (x<0) 
    {
      fprintf(stderr,"Error in function fast_addscore: argument %g is negative\n",x);
      exit(7);
    }
  if (!initialized)   //First fill in the log2-vector
    {
      for (int i=0; i<=2000; i++) val[i]=log2(1.0+pow(2,-0.01*(i+0.5)));
      initialized=1;
    }  
  return val[(int)(100.0*x)];
}



//////////////////////////////////////////////////////////////////////////////
/**
 * @brief Little utilities for output
 */
inline void 
fout(FILE* outf, int d)
{
  if (d>=99999) fprintf(outf,"*\t"); else fprintf(outf,"%i\t",d); 
  return;
}

//////////////////////////////////////////////////////////////////////////////
/**
 * @brief Errors
 */
int 
FormatError(const char infile[], const char details[]="")
{
  cerr<<"Error in "<</*par.argv[0],FS*/__FILE__<<": wrong format while reading file \'"<<infile<<". "<<details<<"\n"; 
  exit(1);
}

int 
OpenFileError(const char outfile[])
{
  cerr<<endl<<"Error in "<</*par.argv[0],FS*/__FILE__<<": could not open file \'"<<outfile<<"\'\n"; 
  exit(2);
}

int 
MemoryError(const char arrayname[])
{
  cerr<<"Error in "<</*par.argv[0],FS*/__FILE__<<": Memory overflow while creating \'"<<arrayname<<"\'. Please report this bug to developers\n"; 
  exit(3);
}

int 
InternalError(const char errstr[])
{
  cerr<<"Error in "<</*par.argv[0],FS*/__FILE__<<":  "<<errstr<<". Please report this bug to developers\n"; 
  exit(6);
}


//////////////////////////////////////////////////////////////////////////////
/**
 * @brief Takes family code (eg. a.1.2.3) and returns strings 'a', 'a.1', and 'a.1.2'
 */
inline void  
ScopID(char cl[], char fold[], char sfam[], const char fam[])
{
  char* ptr;

  //get scop class ID 
  strcpy(cl,fam); 
  ptr = strchr(cl,'.');               //return adress of next '.' in name
  if(ptr) ptr[0]='\0';  

  //get scop fold ID
  strcpy(fold,fam); 
  ptr = strchr(fold,'.');             //return adress of next '.' in name
  if(ptr) ptr = strchr(ptr+1,'.');    //return adress of next '.' in name
  if(ptr) ptr[0]='\0';

  //get scop superfamily ID
  strcpy(sfam,fam); 
  ptr = strchr(sfam,'.');            //return adress of next '.' in name
  if(ptr) ptr = strchr(ptr+1,'.');   //return adress of next '.' in name
  if(ptr) ptr = strchr(ptr+1,'.');   //return adress of next '.' in name
  if(ptr) ptr[0]='\0';
  return;
}

//////////////////////////////////////////////////////////////////////////////
/**
 * @brief Read up to n lines of outfile and write to screen (STDERR)
 */
void 
WriteToScreen(char* outfile, int n)
{
  char line[LINELEN]="";
  ifstream outf;
  outf.open(outfile, ios::in);
  if (!outf) {OpenFileError(outfile);}
  cout<<"\n";
  for(; n>0 && outf.getline(line,LINELEN); n--) cout<<line<<"\n";
  outf.close();
  cout<<"\n";
}
  
inline void 
WriteToScreen(char* outfile) {WriteToScreen(outfile,INT_MAX);}



/////////////////////////////////////////////////////////////////////////////////////
/**
 * @brief Read .hhdefaults file into array argv_conf (beginning at argv_conf[1])
 */
void 
ReadDefaultsFile(int& argc_conf, char** argv_conf)
{
  char line[LINELEN]="";
  char filename[NAMELEN];
  char* c_first;   //pointer to first character of argument string
  char* c;         //pointer to scan line read in for end of argument
  //  ifstream configf;
  FILE* configf=NULL;
  argc_conf=1;     //counts number of arguments read in 

  // Open config file
  strcpy(filename,"./.hhdefaults");
  configf = fopen(filename,"r");
  if (!configf && getenv("HOME"))  
    {
      strcpy(filename,getenv("HOME"));
      strcat(filename,"/.hhdefaults");
      configf = fopen(filename,"r");
      if (!configf)
	{
	  if (v>=3) cerr<<"Warning: could not find ./.hhdefaults or "<<filename<<"\n";
	  return;
	}
    }
  else if (!configf) return; // only webserver has no home directory => need no warning

  // Scan file until line 'program_nameANYTHING'
  while (fgets(line,LINELEN,configf))
    if (!strncmp(line,program_name,6)) break;
  // Found line 'program_nameANYTHING'?
  if (!strncmp(line,program_name,6))
    {
      // Read in options until end-of-file or empty line
	while (fgets(line,LINELEN,configf) && strcmp(line,"\n"))
	{
	  // Analyze line
	  c=line;
	  do
	    {
	      // Find next word
	      while (*c==' ' || *c=='\t') c++; //Advance until next non-white space
	      if (*c=='\0' || *c=='\n' || *c=='#') break;  //Is next word empty string?
	      c_first=c;
	      while (*c!=' ' && *c!='\t'  && *c!='#' && *c!='\0' && *c!='\n' ) c++; //Advance until next white space or '#'
	      if (*c=='\0' || *c=='\n' || *c=='#')         //Is end of line reached?
		{
		  *c='\0'; 
		  argv_conf[argc_conf]=new char[strlen(c_first)+1];
		  strcpy(argv_conf[argc_conf++],c_first); 
		  break;
		}
	      *c='\0'; 
	      argv_conf[argc_conf]=new char[strlen(c_first)+1];
	      strcpy(argv_conf[argc_conf++],c_first);
	      printf("Argument: %s\n",c_first);
	      c++;
	    } while (1);
	} //end read line
     if (v>=3) 
	{
	  cout<<"Arguments read in from .hhdefaults:";
	  for (int argc=1; argc<argc_conf; argc++) cout<<(argv_conf[argc][0]=='-'? " ":"")<<argv_conf[argc]<<" ";
	  cout<<"\n";
	}
     else if (v>=3) cout<<"Read in "<<argc_conf<<" default arguments for "<<program_name<<" from "<<filename<<"\n";
     }
  else //found no line 'program_name   anything"
    {
      if (v>=3) cerr<<endl<<"Warning: no default options for \'"<<program_name<<"\' found in "<<filename<<"\n";
      return; //no line 'program_name   anything' found
    }
//   configf.close();
  fclose(configf);
}


/////////////////////////////////////////////////////////////////////////////////////
/**
 * @brief Set default parameter values
 */
void 
SetDefaults(hhalign_para rHhalignPara)
{
    
  par.append=0;                // overwrite output file
  par.outformat=0;             // 0: hhr  1: FASTA  2:A2M   3:A3M 
  par.p=20.0f;                 // minimum threshold for inclusion in hit list and alignment listing
  par.E=1e6f;                  // maximum threshold for inclusion in hit list and alignment listing
  par.b=10;                    // min number of alignments
  par.B=500;                   // max number of alignments
  par.z=10;                    // min number of lines in hit list
  par.Z=500;                   // max number of lines in hit list
  par.e=1e-3f;                 // maximum E-value for inclusion in output alignment, output HMM, and PSI-BLAST checkpoint model
  par.showcons=1;              // show consensus sequence
  par.showdssp=1;              // show predicted secondary structure
  par.showpred=1;              // show dssp secondary structure
  par.cons=0;                  // show first non-SS sequence as main representative sequence (not consensus)
  par.nseqdis=1;               // maximum number of query sequences for output alignment
  par.mark=0;                  // 1: only marked sequences (or first) get displayed; 0: most divergent ones get displayed
  par.aliwidth=80;             // number of characters per line in output alignments for HMM search
  par.max_seqid=90;            // default for maximum sequence identity threshold
  par.qid=0;                   // default for minimum sequence identity with query
  par.qsc=-20.0f;              // default for minimum score per column with query
  par.coverage=0;              // default for minimum coverage threshold
  par.Ndiff=100;               // pick Ndiff most different sequences from alignment
  par.coverage_core=80;        // Minimum coverage for sequences in core alignment
  par.qsc_core=0.3f;           // Minimum score per column of core sequence with query
  par.coresc=-20.0f;           // Minimum score per column with core alignment (HMM)

  par.M=1;                     // match state assignment is by A2M/A3M
  par.Mgaps=50;                // Above this percentage of gaps, columns are assigned to insert states (for par.M=2)
  par.calibrate=0;             // default: no calibration
  par.calm=0;                  // derive P-values from: 0:query calibration  1:template calibration  2:both
  par.mode=0;                  // 
  
  par.wg=0;                    // 0: use local sequence weights   1: use local ones

  par.matrix=0;                // Subst.matrix 0: Gonnet, 1: HSDM, 2: BLOSUM50 3: BLOSUM62 
  par.pcm=2;                   // pseudocount mode: default=divergence-dependent (but not column-specific)
#if 1 /* Nelder-Meade on Baliscore*/
  par.pca=1.712190f;                // default values for substitution matrix pseudocounts 
  par.pcb=1.039640f;                // significant reduction of pcs by Neff_M starts around Neff_M-1=pcb
  par.pcc=0.878067f;                // pcs are reduced prop. to 1/Neff^pcc
  par.pcw=0.0f;                // wc>0 weighs columns according to their intra-clomun similarity

  par.gapb=1.405220;                // default values for transition pseudocounts
  par.gapd=1.316760;               // gap open penalty pseudocount; 0.25 corresponds to 7.1*gapf bits
  par.gape=1.793780;                // gap extension penalty pseudocount
  par.gapf=1.034710;                // factor for increasing gap open penalty for deletes
  par.gapg=0.894772;                // factor for increasing gap open penalty for inserts
  par.gaph=0.544072;                // factor for increasing gap extension penalty for deletes
  par.gapi=0.862559;                // factor for increasing gap extension penalty for inserts
#else  /* Soeding's default*/
  par.pca=1.0f;                // default values for substitution matrix pseudocounts 
  par.pcb=1.5f;                // significant reduction of pcs by Neff_M starts around Neff_M-1=pcb
  par.pcc=1.0f;                // pcs are reduced prop. to 1/Neff^pcc
  par.pcw=0.0f;                // wc>0 weighs columns according to their intra-clomun similarity

  par.gapb=1.0;                // default values for transition pseudocounts
  par.gapd=0.15;               // gap open penalty pseudocount; 0.25 corresponds to 7.1*gapf bits
  par.gape=1.0;                // gap extension penalty pseudocount
  par.gapf=0.6;                // factor for increasing gap open penalty for deletes
  par.gapg=0.6;                // factor for increasing gap open penalty for inserts
  par.gaph=0.6;                // factor for increasing gap extension penalty for deletes
  par.gapi=0.6;                // factor for increasing gap extension penalty for inserts
#endif

#if 0
  /* Viterbi parameters optimised on Sabre (R228), FS, r228 -> r229 */
  par.pcaV=1.245150f;                // default values for substitution matrix pseudocounts 
  par.pcbV=1.682110f;                // significant reduction of pcs by Neff_M starts around Neff_M-1=pcb
  par.pccV=1.483840f;                // pcs are reduced prop. to 1/Neff^pcc
  par.pcwV=0.0f;                // wc>0 weighs columns according to their intra-clomun similarity

  par.gapbV=0.818625;                // default values for transition pseudocounts
  par.gapdV=0.666110;               // gap open penalty pseudocount; 0.25 corresponds to 7.1*gapf bits
  par.gapeV=1.028050;                // gap extension penalty pseudocount
  par.gapfV=0.710760;                // factor for increasing gap open penalty for deletes
  par.gapgV=1.649800;                // factor for increasing gap open penalty for inserts
  par.gaphV=0.470604;                // factor for increasing gap extension penalty for deletes
  par.gapiV=0.829479;                // factor for increasing gap extension penalty for inserts
#elif 1
  /* Viterbi parameters optimised on Balibase, r244 -> r245 */
  par.pcaV=1.333860f;                // default values for substitution matrix pseudocounts 
  par.pcbV=1.934480f;                // significant reduction of pcs by Neff_M starts around Neff_M-1=pcb
  par.pccV=1.655610f;                // pcs are reduced prop. to 1/Neff^pcc
  par.pcwV=0.0f;                // wc>0 weighs columns according to their intra-clomun similarity

  par.gapbV=0.334525;                // default values for transition pseudocounts
  par.gapdV=0.074534;               // gap open penalty pseudocount; 0.25 corresponds to 7.1*gapf bits
  par.gapeV=0.320336;                // gap extension penalty pseudocount
  par.gapfV=0.151634;                // factor for increasing gap open penalty for deletes
  par.gapgV=0.641516;                // factor for increasing gap open penalty for inserts
  par.gaphV=0.266434;                // factor for increasing gap extension penalty for deletes
  par.gapiV=0.598414;                // factor for increasing gap extension penalty for inserts 
#else /* Soeding default*/ 
  par.pcaV=1.0f;                // default values for substitution matrix pseudocounts 
  par.pcbV=1.5f;                // significant reduction of pcs by Neff_M starts around Neff_M-1=pcb
  par.pccV=1.0f;                // pcs are reduced prop. to 1/Neff^pcc
  par.pcwV=0.0f;                // wc>0 weighs columns according to their intra-clomun similarity

  par.gapbV=1.0;                // default values for transition pseudocounts
  par.gapdV=0.15;               // gap open penalty pseudocount; 0.25 corresponds to 7.1*gapf bits
  par.gapeV=1.0;                // gap extension penalty pseudocount
  par.gapfV=0.6;                // factor for increasing gap open penalty for deletes
  par.gapgV=0.6;                // factor for increasing gap open penalty for inserts
  par.gaphV=0.6;                // factor for increasing gap extension penalty for deletes
  par.gapiV=0.6;                // factor for increasing gap extension penalty for inserts
#endif

  par.ssm=2;                   // ss scoring mode: 0:no ss score  1:score after alignment  2:score during alignment
  par.ssw=0.11f;               // weight of ss scoring
  par.ssa=1.0f;                // weight of ss evolution matrix
  par.shift=-0.01f;            // Shift match score up 
  par.mact=0.3001f;            // Score threshold for MAC alignment in local mode (set to 0.5001 to track user modification)
  par.corr=0.1f;               // Weight of correlations of scores for |i-j|<=4
  par.wstruc=1.0f;             // Weight of structure scores

  par.egq=0.0f;                // no charge for end gaps as default
  par.egt=0.0f;                // no charge for end gaps as default

  par.trans=0;                 // no transitive scoring as default
  par.Emax_trans=100.0f;       // use intermediate HMMs with E-values up to 100 between query and database HMM
  par.Emax_trans=100.0f;       // use intermediate HMMs with E-values up to 100 between query and database HMM
  par.wtrans=1.0f;             // Ztot[k] = Zq[k] + wtrans * (Zforward[k]+Zreverse[k])
  par.ssgap=0;                 // 1: add secondary structure-dependent gap penalties  0:off
  par.ssgapd=1.0f;             // secondary structure-dependent gap-opening penalty (per residue)
  par.ssgape=0.0f;             // secondary structure-dependent gap-extension penalty (per residue)
  par.ssgapi=4;                // max. number of inside-integer(ii); gap-open-penalty= -ii*ssgapd

  par.loc=1;                   // local vs. global alignment as default
  par.altali=2;                // find up to two (possibly overlapping) subalignments 
  par.forward=0;               // 0: Viterbi algorithm; 1: Viterbi+stochastic sampling; 3:Maximum Accuracy (MAC) algorithm
  par.realign=1;               // realign with MAC algorithm

  par.repmode=0;               // repeats score independently of one another
  par.columnscore=1;           // Default column score is 1: null model pnul = 1/2 * (q_av(a)+p_av(a))
  par.min_overlap=0;           // automatic minimum overlap used
  par.opt=0;                   // Default = optimization mode off
  par.readdefaultsfile=0;      // Default = do not read a defaults file ./.hhdefaults or HOME/.hhdefaults
  par.maxdbstrlen=200;         // maximum length of database string to be printed in 'Command' line of hhr file
  par.mode=0;
  par.idummy=par.jdummy=0;     //

  par.notags=1;                // neutralize His-tags, FLAG-tags, C-myc-tags

  // Initialize strings
  strcpy(par.infile,"stdin");
  strcpy(par.outfile,"");
  strcpy(par. pairwisealisfile,"");
  strcpy(par.buffer,"buffer.txt"); 
  strcpy(par.scorefile,""); 
  strcpy(par.wfile,""); 
  strcpy(par.alnfile,""); 
  strcpy(par.hhmfile,""); 
  strcpy(par.psifile,""); 
  par.exclstr=NULL; 


  //#if PARAMETERSFROMFILE  /* read parameter file from home-dir */
  //#include "hhutil-C-help.h"
  //#endif /* read parameter file from home-dir */
  if (rHhalignPara.pca  >= 0.00){ par.pca  = rHhalignPara.pca;  }
  if (rHhalignPara.pcb  >= 0.00){ par.pcb  = rHhalignPara.pcb;  }
  if (rHhalignPara.pcc  >= 0.00){ par.pcc  = rHhalignPara.pcc;  }
  if (rHhalignPara.pcw  >= 0.00){ par.pcw  = rHhalignPara.pcw;  }
  if (rHhalignPara.gapb >= 0.00){ par.gapb = rHhalignPara.gapb; }
  if (rHhalignPara.gapd >= 0.00){ par.gapd = rHhalignPara.gapd; }
  if (rHhalignPara.gape >= 0.00){ par.gape = rHhalignPara.gape; }
  if (rHhalignPara.gapf >= 0.00){ par.gapf = rHhalignPara.gapf; }
  if (rHhalignPara.gapg >= 0.00){ par.gapg = rHhalignPara.gapg; }
  if (rHhalignPara.gaph >= 0.00){ par.gaph = rHhalignPara.gaph; }
  if (rHhalignPara.gapi >= 0.00){ par.gapi = rHhalignPara.gapi; }
  if (rHhalignPara.pcaV  >= 0.00){ par.pcaV  = rHhalignPara.pcaV;  }
  if (rHhalignPara.pcbV  >= 0.00){ par.pcbV  = rHhalignPara.pcbV;  }
  if (rHhalignPara.pccV  >= 0.00){ par.pccV  = rHhalignPara.pccV;  }
  if (rHhalignPara.pcwV  >= 0.00){ par.pcwV  = rHhalignPara.pcwV;  }
  if (rHhalignPara.gapbV >= 0.00){ par.gapbV = rHhalignPara.gapbV; }
  if (rHhalignPara.gapdV >= 0.00){ par.gapdV = rHhalignPara.gapdV; }
  if (rHhalignPara.gapeV >= 0.00){ par.gapeV = rHhalignPara.gapeV; }
  if (rHhalignPara.gapfV >= 0.00){ par.gapfV = rHhalignPara.gapfV; }
  if (rHhalignPara.gapgV >= 0.00){ par.gapgV = rHhalignPara.gapgV; }
  if (rHhalignPara.gaphV >= 0.00){ par.gaphV = rHhalignPara.gaphV; }
  if (rHhalignPara.gapiV >= 0.00){ par.gapiV = rHhalignPara.gapiV; }

   return;
} /** this is the end of SetDefaults() **/

void SetRnaDefaults(hhalign_para rHhalignPara)
{
  par.pca=1.28f;                // default values for substitution matrix pseudocounts 
  par.pcb=1.75f;                // significant reduction of pcs by Neff_M starts around Neff_M-1=pcb
  par.pcc=0.87f;                // pcs are reduced prop. to 1/Neff^pcc
  par.pcw=0.0f;                // wc>0 weighs columns according to their intra-clomun similarity

  par.gapb=0.80;                // default values for transition pseudocounts
  par.gapd=0.34;               // gap open penalty pseudocount; 0.25 corresponds to 7.1*gapf bits
  par.gape=2.25;                // gap extension penalty pseudocount
  par.gapf=0.51;                // factor for increasing gap open penalty for deletes
  par.gapg=1.01;                // factor for increasing gap open penalty for inserts
  par.gaph=1.24;                // factor for increasing gap extension penalty for deletes
  par.gapi=0.45;                // factor for increasing gap extension penalty for inserts

#if 0  /* these are the parameters determined by Dave (pre r274) */
  par.pcaV=2.57f;                // default values for substitution matrix pseudocounts 
  par.pcbV=2.34f;                // significant reduction of pcs by Neff_M starts around Neff_M-1=pcb
  par.pccV=0.88f;                // pcs are reduced prop. to 1/Neff^pcc
  par.pcwV=0.0f;                // wc>0 weighs columns according to their intra-clomun similarity

  par.gapbV=1.41;                // default values for transition pseudocounts
  par.gapdV=1.8;               // gap open penalty pseudocount; 0.25 corresponds to 7.1*gapf bits
  par.gapeV=1.5;                // gap extension penalty pseudocount
  par.gapfV=1.03;                // factor for increasing gap open penalty for deletes
  par.gapgV=0.89;                // factor for increasing gap open penalty for inserts
  par.gaphV=0.54;                // factor for increasing gap extension penalty for deletes
  par.gapiV=0.86;                // factor for increasing gap extension penalty for inserts  
#else /* parameters determined for r274, using Bralibase, FS */
  par.pcaV=1.655620f;                // default values for substitution matrix pseudocounts 
  par.pcbV=0.438399f;                // significant reduction of pcs by Neff_M starts around Neff_M-1=pcb
  par.pccV=0.371491f;                // pcs are reduced prop. to 1/Neff^pcc
  par.pcwV=0.0f;                // wc>0 weighs columns according to their intra-clomun similarity

  par.gapbV=1.914490;                // default values for transition pseudocounts
  par.gapdV=0.104278;               // gap open penalty pseudocount; 0.25 corresponds to 7.1*gapf bits
  par.gapeV=1.100210;                // gap extension penalty pseudocount
  par.gapfV=0.335152;                // factor for increasing gap open penalty for deletes
  par.gapgV=0.786688;                // factor for increasing gap open penalty for inserts
  par.gaphV=0.667143;                // factor for increasing gap extension penalty for deletes
  par.gapiV=0.711993;                // factor for increasing gap extension penalty for inserts  
#endif

  //#if PARAMETERSFROMFILE  /* read parameter file from home-dir */
  //#include "hhutil-C-help.h"
  //#endif /* read parameter file from home-dir */
  if (rHhalignPara.pca  >= 0.00){ par.pca  = rHhalignPara.pca;  }
  if (rHhalignPara.pcb  >= 0.00){ par.pcb  = rHhalignPara.pcb;  }
  if (rHhalignPara.pcc  >= 0.00){ par.pcc  = rHhalignPara.pcc;  }
  if (rHhalignPara.pcw  >= 0.00){ par.pcw  = rHhalignPara.pcw;  }
  if (rHhalignPara.gapb >= 0.00){ par.gapb = rHhalignPara.gapb; }
  if (rHhalignPara.gapd >= 0.00){ par.gapd = rHhalignPara.gapd; }
  if (rHhalignPara.gape >= 0.00){ par.gape = rHhalignPara.gape; }
  if (rHhalignPara.gapf >= 0.00){ par.gapf = rHhalignPara.gapf; }
  if (rHhalignPara.gapg >= 0.00){ par.gapg = rHhalignPara.gapg; }
  if (rHhalignPara.gaph >= 0.00){ par.gaph = rHhalignPara.gaph; }
  if (rHhalignPara.gapi >= 0.00){ par.gapi = rHhalignPara.gapi; }
  if (rHhalignPara.pcaV  >= 0.00){ par.pcaV  = rHhalignPara.pcaV;  }
  if (rHhalignPara.pcbV  >= 0.00){ par.pcbV  = rHhalignPara.pcbV;  }
  if (rHhalignPara.pccV  >= 0.00){ par.pccV  = rHhalignPara.pccV;  }
  if (rHhalignPara.pcwV  >= 0.00){ par.pcwV  = rHhalignPara.pcwV;  }
  if (rHhalignPara.gapbV >= 0.00){ par.gapbV = rHhalignPara.gapbV; }
  if (rHhalignPara.gapdV >= 0.00){ par.gapdV = rHhalignPara.gapdV; }
  if (rHhalignPara.gapeV >= 0.00){ par.gapeV = rHhalignPara.gapeV; }
  if (rHhalignPara.gapfV >= 0.00){ par.gapfV = rHhalignPara.gapfV; }
  if (rHhalignPara.gapgV >= 0.00){ par.gapgV = rHhalignPara.gapgV; }
  if (rHhalignPara.gaphV >= 0.00){ par.gaphV = rHhalignPara.gaphV; }
  if (rHhalignPara.gapiV >= 0.00){ par.gapiV = rHhalignPara.gapiV; }

} /* this is the end of SetRnaDefaults() */

void SetDnaDefaults(hhalign_para rHhalignPara)
{
  par.pca=2.89f;                // default values for substitution matrix pseudocounts 
  par.pcb=1.17f;                // significant reduction of pcs by Neff_M starts around Neff_M-1=pcb
  par.pcc=0.88f;                // pcs are reduced prop. to 1/Neff^pcc
  par.pcw=0.0f;                // wc>0 weighs columns according to their intra-clomun similarity

  par.gapb=0.80;                // default values for transition pseudocounts
  par.gapd=0.34;               // gap open penalty pseudocount; 0.25 corresponds to 7.1*gapf bits
  par.gape=2.25;                // gap extension penalty pseudocount
  par.gapf=0.78;                // factor for increasing gap open penalty for deletes
  par.gapg=1.01;                // factor for increasing gap open penalty for inserts
  par.gaph=1.24;                // factor for increasing gap extension penalty for deletes
  par.gapi=0.45;                // factor for increasing gap extension penalty for inserts
  
#if 0 /* these are the parameters determined by Dave (pre r274) */
  par.pcaV=1.712f;                // default values for substitution matrix pseudocounts 
  par.pcbV=1.039f;                // significant reduction of pcs by Neff_M starts around Neff_M-1=pcb
  par.pccV=0.266f;                // pcs are reduced prop. to 1/Neff^pcc
  par.pcwV=0.0f;                // wc>0 weighs columns according to their intra-clomun similarity

  par.gapbV=1.405;                // default values for transition pseudocounts
  par.gapdV=1.8;               // gap open penalty pseudocount; 0.25 corresponds to 7.1*gapf bits
  par.gapeV=2.25;                // gap extension penalty pseudocount
  par.gapfV=1.034;                // factor for increasing gap open penalty for deletes
  par.gapgV=2.025;                // factor for increasing gap open penalty for inserts
  par.gaphV=0.544;                // factor for increasing gap extension penalty for deletes
  par.gapiV=1.35;                // factor for increasing gap extension penalty for inserts }
#else /* parameters determined for r274, using mdsa, FS */
  par.pcaV=2.196;                // default values for substitution matrix pseudocounts 
  par.pcbV=0.329;                // significant reduction of pcs by Neff_M starts around Neff_M-1=pcb
  par.pccV=0.393;                // pcs are reduced prop. to 1/Neff^pcc
  par.pcwV=0.0f;                // wc>0 weighs columns according to their intra-clomun similarity

  par.gapbV=0.570;                // default values for transition pseudocounts
  par.gapdV=0.048;               // gap open penalty pseudocount; 0.25 corresponds to 7.1*gapf bits
  par.gapeV=1.692;                // gap extension penalty pseudocount
  par.gapfV=0.398;                // factor for increasing gap open penalty for deletes
  par.gapgV=0.121;                // factor for increasing gap open penalty for inserts
  par.gaphV=0.012;                // factor for increasing gap extension penalty for deletes
  par.gapiV=0.645;                // factor for increasing gap extension penalty for inserts }
#endif

  //#if PARAMETERSFROMFILE  /* read parameter file from home-dir */
  //#include "hhutil-C-help.h"
  //#endif /* read parameter file from home-dir */
  if (rHhalignPara.pca  >= 0.00){ par.pca  = rHhalignPara.pca;  }
  if (rHhalignPara.pcb  >= 0.00){ par.pcb  = rHhalignPara.pcb;  }
  if (rHhalignPara.pcc  >= 0.00){ par.pcc  = rHhalignPara.pcc;  }
  if (rHhalignPara.pcw  >= 0.00){ par.pcw  = rHhalignPara.pcw;  }
  if (rHhalignPara.gapb >= 0.00){ par.gapb = rHhalignPara.gapb; }
  if (rHhalignPara.gapd >= 0.00){ par.gapd = rHhalignPara.gapd; }
  if (rHhalignPara.gape >= 0.00){ par.gape = rHhalignPara.gape; }
  if (rHhalignPara.gapf >= 0.00){ par.gapf = rHhalignPara.gapf; }
  if (rHhalignPara.gapg >= 0.00){ par.gapg = rHhalignPara.gapg; }
  if (rHhalignPara.gaph >= 0.00){ par.gaph = rHhalignPara.gaph; }
  if (rHhalignPara.gapi >= 0.00){ par.gapi = rHhalignPara.gapi; }
  if (rHhalignPara.pcaV  >= 0.00){ par.pcaV  = rHhalignPara.pcaV;  }
  if (rHhalignPara.pcbV  >= 0.00){ par.pcbV  = rHhalignPara.pcbV;  }
  if (rHhalignPara.pccV  >= 0.00){ par.pccV  = rHhalignPara.pccV;  }
  if (rHhalignPara.pcwV  >= 0.00){ par.pcwV  = rHhalignPara.pcwV;  }
  if (rHhalignPara.gapbV >= 0.00){ par.gapbV = rHhalignPara.gapbV; }
  if (rHhalignPara.gapdV >= 0.00){ par.gapdV = rHhalignPara.gapdV; }
  if (rHhalignPara.gapeV >= 0.00){ par.gapeV = rHhalignPara.gapeV; }
  if (rHhalignPara.gapfV >= 0.00){ par.gapfV = rHhalignPara.gapfV; }
  if (rHhalignPara.gapgV >= 0.00){ par.gapgV = rHhalignPara.gapgV; }
  if (rHhalignPara.gaphV >= 0.00){ par.gaphV = rHhalignPara.gaphV; }
  if (rHhalignPara.gapiV >= 0.00){ par.gapiV = rHhalignPara.gapiV; }


} /* this is the end of SetDnaDefaults() */


/*
 * EOF hhutil-C.h
 */
