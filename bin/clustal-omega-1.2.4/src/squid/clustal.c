/*****************************************************************
 * SQUID - a library of functions for biological sequence analysis
 * Copyright (C) 1992-2002 Washington University School of Medicine
 * 
 *     This source code is freely distributed under the terms of the
 *     GNU General Public License. See the files COPYRIGHT and LICENSE
 *     for details.
 *****************************************************************/

/* clustal.c
 * SRE, Sun Jun  6 17:50:45 1999 [bus from Madison, 1999 worm mtg]
 * 
 * Import/export of ClustalV/W multiple sequence alignment
 * formatted files. Derivative of msf.c; MSF is a pretty
 * generic interleaved format.   
 * 
 * RCS $Id: clustal.c 315 2016-12-15 17:18:30Z fabian $ (Original squid RCS Id: clustal.c,v 1.1 1999/07/15 22:26:53 eddy Exp)
 */

#include <stdbool.h>
#include <unistd.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include "squid.h"
#include "msa.h"

#ifdef CLUSTALO
/* needed for PACKAGE_VERSION */
#include "../config.h"

/* DD: isnumber is in BSD/OSX but not GCC/Linux */
#ifndef isnumber
	#define isnumber(c) ( (c>='0') && (c<='9'))
#endif

#ifndef min
	#define min( a, b ) ( ((a) < (b)) ? (a) : (b) )
#endif

/*These are all the positively scoring groups that occur in the Gonnet Pam250
matrix. There are strong and weak groups, defined as strong score >0.5 and
weak score =<0.5. Strong matching columns to be assigned ':' and weak matches
assigned '.' in the clustal output format.
amino_strong = res_cat1
amino_weak = res_cat2
*/

char *amino_strong[] = {"STA", "NEQK", "NHQK", "NDEQ", "QHRK", "MILV",
    "MILF", "HY", "FYW", NULL};

char *amino_weak[] = {"CSA", "ATV", "SAG", "STNK", "STPA", "SGND",
    "SNDEQK", "NDEQHK", "NEQHRK", "FVLIM", "HFY", NULL};

#endif

#ifdef TESTDRIVE_CLUSTAL
/*****************************************************************
 * msf.c test driver: 
 * cc -DTESTDRIVE_CLUSTAL -g -O2 -Wall -o test clustal.c msa.c gki.c sqerror.c sre_string.c file.c hsregex.c sre_math.c sre_ctype.c -lm
 * 
 */
int
main(int argc, char **argv)
{
  MSAFILE *afp;
  MSA     *msa;
  char    *file;
  
  file = argv[1];

  if ((afp = MSAFileOpen(file, MSAFILE_CLUSTAL, NULL)) == NULL)
    Die("Couldn't open %s\n", file);

  while ((msa = ReadClustal(afp)) != NULL)
    {
#ifdef CLUSTALO
#define LINE_WRAP 60
      WriteClustal(stdout, msa, LINE_WRAP);
#else
      WriteClustal(stdout, msa);
#endif
      MSAFree(msa); 
    }
  
  MSAFileClose(afp);
  exit(0);
}
/******************************************************************/
#endif /* testdrive_clustal */


/* Function: ReadClustal()
 * Date:     SRE, Sun Jun  6 17:53:49 1999 [bus from Madison, 1999 worm mtg]
 *
 * Purpose:  Parse an alignment read from an open Clustal format
 *           alignment file. Clustal is a single-alignment format.
 *           Return the alignment, or NULL if we have no data.
 *           
 * Args:     afp  - open alignment file
 *
 * Returns:  MSA * - an alignment object
 *                   caller responsible for an MSAFree()
 *           NULL if no more alignments
 *
 * Diagnostics: 
 *           Will Die() here with a (potentially) useful message
 *           if a parsing error occurs.
 */
MSA *
ReadClustal(MSAFILE *afp)
{
  MSA    *msa;
  char   *s;
  int     slen;
  int     sqidx;
  char   *name;
  char   *seq;
  char   *s2;

  if (feof(afp->f)) return NULL;

  /* Skip until we see the CLUSTAL header
   */
  while ((s = MSAFileGetLine(afp)) != NULL)
    {
      if (strncmp(s, "CLUSTAL", 7) == 0 &&
	  strstr(s, "multiple sequence alignment") != NULL)
	break;
    }
  if (s == NULL) return NULL;

  msa = MSAAlloc(10, 0);

  /* Now we're in the sequence section. 
   * As discussed above, if we haven't seen a sequence name, then we
   * don't include the sequence in the alignment.
   * Watch out for conservation markup lines that contain *.: chars
   */
  while ((s = MSAFileGetLine(afp)) != NULL) 
    {
      if ((name = sre_strtok(&s, WHITESPACE, NULL))  == NULL) continue;
      if ((seq  = sre_strtok(&s, WHITESPACE, &slen)) == NULL) continue;
      s2 = sre_strtok(&s, "\n", NULL);

      /* The test for a conservation markup line
       */
      if (strpbrk(name, ".*:") != NULL && strpbrk(seq, ".*:") != NULL)
	continue;
#ifdef CLUSTALO
	  /* extra bit at the end of a line might be the unaligned residue
		 count */
      if (s2 != NULL) {
	    int i;
		for (i=0; i<strlen(s2); i++) {
		  if (! isnumber(s2[i])) {
			Die("Parse failed at line %d, file %s: possibly using spaces as gaps",
				afp->linenumber, afp->fname);
          }
		}
	  }
#else
      if (s2 != NULL)
			  Die("Parse failed at line %d, file %s: possibly using spaces as gaps",
				  afp->linenumber, afp->fname);
#endif
      /* It's not blank, and it's not a coord line: must be sequence
       */
      sqidx = MSAGetSeqidx(msa, name, msa->lastidx+1);
      msa->lastidx = sqidx;
      msa->sqlen[sqidx] = sre_strcat(&(msa->aseq[sqidx]), msa->sqlen[sqidx], seq, slen); 
    }

  MSAVerifyParse(msa);		/* verifies, and also sets alen and wgt. */
  return msa;
}

size_t utf8len(char *s)
{
  size_t len = 0;
  for (; *s; ++s) if ((*s & 0xC0) != 0x80) ++len;
  return len;
}


/* Function: WriteClustal()
 * Date:     SRE, Sun Jun  6 18:12:47 1999 [bus from Madison, worm mtg 1999]
 *
 * Purpose:  Write an alignment in Clustal format to an open file.
 *
 * Args:     fp    - file that's open for writing.
 *           msa   - alignment to write. 
 *
 * Returns:  (void)
 */
void
#ifdef CLUSTALO
WriteClustal(FILE *fp, MSA *msa, int iWrap, int bResno, int iSeqType)
#else
WriteClustal(FILE *fp, MSA *msa)
#endif
{
  int    idx;			/* counter for sequences         */
  int    len;			/* tmp variable for name lengths */
  int    namelen;		/* maximum name length used      */
  int    pos;			/* position counter              */
#ifdef CLUSTALO
  char  *buf;    	        /* buffer for writing seq        */
  int    cpl = msa->alen<iWrap ? msa->alen+10 : (iWrap > 0 ? iWrap : 60);		/* char per line (< 64)          */
#else
  char   buf[80];	        /* buffer for writing seq        */
  int    cpl = 60;		/* char per line (< 64)          */
#endif

  /* consensus line stuff */
  int subpos;
  char first;
  int bail;
  int strong_bins[9];
  int weak_bins[11];
  /*int cons;*/
  int bin;

#ifdef CLUSTALO
  int *piResCnt = NULL;
#endif

#ifdef CLUSTALO
  if (1 == bResno){

    if (NULL == (piResCnt = (int *)malloc(msa->nseq * sizeof(int)))){
      Die("%s:%s:%d: could not malloc %d int for residue count", 
	  __FUNCTION__, __FILE__, __LINE__, msa->nseq);
    }
    else {
      memset(piResCnt, 0, msa->nseq * sizeof(int));
    }
  } /* do print residue numbers */

  if (NULL == (buf = (char *)malloc(cpl+20))){
    Die("%s:%s:%d: could not malloc %d char for buffer",
        __FUNCTION__, __FILE__, __LINE__, cpl+20);
  }
  else {
    memset(buf, 0, cpl+20);
  }

#endif


  /* calculate max namelen used */
  namelen = 0;
  for (idx = 0; idx < msa->nseq; idx++){
    /*if ((len = strlen(msa->sqname[idx])) > namelen) */ /* strlen() gives problems for unicode, FS, -> 290 */
    if ((len = utf8len(msa->sqname[idx])) > namelen)
      {
	namelen = len; 
      }
  }

#ifdef CLUSTALO
  fprintf(fp, "CLUSTAL O(%s) multiple sequence alignment\n", PACKAGE_VERSION);
#else
  fprintf(fp, "CLUSTAL W(1.5) multiple sequence alignment\n");
#endif
  
  /*****************************************************
   * Write the sequences
   *****************************************************/

#ifdef CLUSTALO
    fprintf(fp, "\n");	/* original had two blank lines */
#endif

  for (pos = 0; pos < msa->alen; pos += cpl)
    {
      fprintf(fp, "\n");	/* Blank line between sequence blocks */
      for (idx = 0; idx < msa->nseq; idx++)
      {
        strncpy(buf, msa->aseq[idx] + pos, cpl);
	    buf[cpl] = '\0';
#ifdef CLUSTALO
	    if (1 == bResno){
	      char *pc = NULL;
	      for (pc = buf; *pc != '\0'; pc++){
		if ( ( (*pc >= 'a') && (*pc <= 'z') ) || ( (*pc >= 'A') && (*pc <= 'Z') ) ){
		  piResCnt[idx]++;
		}
	      }
	      /* printf("%*s") gives problems for unicode, FS, -> 290 */
	      /*fprintf(fp, "%-*s\t%s\t%d\n", namelen+5, msa->sqname[idx], buf, piResCnt[idx]);*/
	      fprintf(fp, "%s%*s %s\t%d\n", msa->sqname[idx], (int)(namelen+5-utf8len(msa->sqname[idx])), "", buf, piResCnt[idx]);
	    }
	    else {
	      /* printf("%*s") gives problems for unicode, FS, -> 290 */
	      /*fprintf(fp, "%-*s\t%s\n", namelen+5, msa->sqname[idx], buf);*/
	      fprintf(fp, "%s%*s %s\n", msa->sqname[idx], (int)(namelen+5-utf8len(msa->sqname[idx])), "", buf);
	    }
#else
	    fprintf(fp, "%*s %s\n", namelen, msa->sqname[idx], buf);
#endif
	  }
#ifdef CLUSTALO
      /* do consensus dots */

      /* print namelen+5 spaces */
      for(subpos = 0; subpos <= namelen+5; subpos++)
        fprintf(fp, " ");

      for(subpos = pos; subpos < min(pos + cpl, msa->alen); subpos++)
      {
          /* see if 100% conservation */
          first = msa->aseq[0][subpos];
          bail = 0;
          for (idx = 1; idx < msa->nseq; idx++)
          {
            if(toupper(msa->aseq[idx][subpos]) != toupper(first)) /* toupper makes consensus case-insensitive, FS, r290 -> */
            {
              bail = 1;
              break;
            }
          }
          if(!bail)
            fprintf(fp, "*");
          else
          {
            /* if not then check strong */
            for(bin = 0; bin < 9; bin++)
              strong_bins[bin] = 0; /* clear the bins */

            for (idx = 0; (iSeqType == kAmino) && (idx < msa->nseq); idx++)
            { /* do this only for amino acids, no strong/weak consensus for nucleotide, FS, r290 -> */
              switch(toupper(msa->aseq[idx][subpos])) /* toupper makes consensus case-insensitive, FS, r290 -> */
              {
                case 'S': strong_bins[0]++; break;
                case 'T': strong_bins[0]++; break;
                case 'A': strong_bins[0]++; break;
                case 'N': strong_bins[1]++; strong_bins[2]++; strong_bins[3]++; break;
                case 'E': strong_bins[1]++; strong_bins[3]++; break;
                case 'Q': strong_bins[1]++; strong_bins[2]++; strong_bins[3]++; strong_bins[4]++; break;
                case 'K': strong_bins[1]++; strong_bins[2]++; strong_bins[4]++; break;
                case 'D': strong_bins[3]++; break;
                case 'R': strong_bins[4]++; break;
	        case 'H': strong_bins[2]++; strong_bins[4]++; strong_bins[7]++; break; /* added bin-2 (NHQK), FS 2016-07-14 */
                case 'M': strong_bins[5]++; strong_bins[6]++; break;
                case 'I': strong_bins[5]++; strong_bins[6]++; break;
                case 'L': strong_bins[5]++; strong_bins[6]++; break;
                case 'V': strong_bins[5]++; break;
                case 'F': strong_bins[6]++; strong_bins[8]++; break;
                case 'Y': strong_bins[7]++; strong_bins[8]++; break;
                case 'W': strong_bins[8]++; break;
              }
            }
            bail = 0;
            for(bin = 0; bin < 9; bin++)
              if(strong_bins[bin] == msa->nseq)
              {
                  bail = 1;
                  break;
              }
            if(bail)
              fprintf(fp, ":");
            else
            {
              /* check weak */
              for(bin = 0; bin < 11; bin++)
                weak_bins[bin] = 0; /* clear the bins */

              for(idx = 0; (iSeqType == kAmino) && (idx < msa->nseq); idx++)
	      { /* do this only for amino acids, no strong/weak consensus for nucleotide, FS, r290 -> */
                switch(toupper(msa->aseq[idx][subpos])) /* toupper makes consensus case-insensitive, FS, r290 -> */
                {
                  case 'C': weak_bins[0]++; break;
                  case 'S': weak_bins[0]++; weak_bins[2]++; weak_bins[3]++; weak_bins[4]++; weak_bins[5]++; weak_bins[6]++; break;
                  case 'A': weak_bins[0]++; weak_bins[1]++; weak_bins[2]++; weak_bins[4]++; break;
                  case 'T': weak_bins[1]++; weak_bins[3]++; weak_bins[4]++; break;
                  case 'V': weak_bins[1]++; weak_bins[9]++; break;
		  case 'G': weak_bins[2]++; weak_bins[5]++; break; /* Added bin-5, FS 2016-07-14 */
                  case 'N': weak_bins[3]++; weak_bins[5]++; weak_bins[6]++; weak_bins[7]++; weak_bins[8]++; break;
                  case 'K': weak_bins[3]++; weak_bins[6]++; weak_bins[7]++; weak_bins[8]++; break;
                  case 'D': weak_bins[5]++; weak_bins[6]++; weak_bins[7]++; break;
                  case 'E': weak_bins[6]++; weak_bins[7]++; weak_bins[8]++; break;
                  case 'Q': weak_bins[6]++; weak_bins[7]++; weak_bins[8]++; break;
                  case 'H': weak_bins[7]++; weak_bins[8]++; weak_bins[10]++; break;
                  case 'R': weak_bins[8]++; break;
                  case 'F': weak_bins[9]++; weak_bins[10]++; break;
                  case 'L': weak_bins[9]++; break;
                  case 'I': weak_bins[9]++; break;
                  case 'M': weak_bins[9]++; break;
                  case 'Y': weak_bins[10]++; break;
                }
              }
              bail = 0;
              for(bin = 0; bin < 11; bin++)
                if(weak_bins[bin] == msa->nseq)
                {
                    bail = 1;
                    break;
                }
              if(bail)
                fprintf(fp, ".");
              else
                fprintf(fp, " ");
            }
          }
      }
      fprintf(fp,"\n");
#endif
    }

#ifdef CLUSTAL
  free(piResCnt); piResCnt = NULL;
#endif

  return;
}



