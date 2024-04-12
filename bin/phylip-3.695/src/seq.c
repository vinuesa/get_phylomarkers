
#include "phylip.h"
#include "seq.h"

/* version 3.6. (c) Copyright 1993-2004 by the University of Washington.
   Written by Joseph Felsenstein, Akiko Fuseki, Sean Lamont, and Andrew Keeffe.
   Permission is granted to copy and use this program provided no fee is
   charged for it and provided that this copyright notice is not removed. */

long nonodes, endsite, outgrno, nextree, which;
boolean interleaved, printdata, outgropt, treeprint, dotdiff, transvp;
steptr weight, category, alias, location, ally;
sequence y;


void fix_x(node* p,long site, double maxx, long rcategs)
{ /* dnaml dnamlk */
  long i,j;
  p->underflows[site] += log(maxx);

  for ( i = 0 ; i < rcategs ; i++ ) {
    for ( j = 0 ; j < ((long)T - (long)A + 1) ; j++)
      p->x[site][i][j] /= maxx;
  }
} /* fix_x */


void fix_protx(node* p,long site, double maxx, long rcategs) 
{ /* proml promlk */
  long i,m;

  p->underflows[site] += log(maxx);

  for ( i = 0 ; i < rcategs  ; i++ ) 
    for (m = 0; m <= 19; m++)
      p->protx[site][i][m] /= maxx;
} /* fix_protx */


void alloctemp(node **temp, long *zeros, long endsite)
{
  /*used in dnacomp and dnapenny */
  *temp = (node *)Malloc(sizeof(node));
  (*temp)->numsteps = (steptr)Malloc(endsite*sizeof(long));
  (*temp)->base = (baseptr)Malloc(endsite*sizeof(long));
  (*temp)->numnuc = (nucarray *)Malloc(endsite*sizeof(nucarray));
  memcpy((*temp)->base, zeros, endsite*sizeof(long));
  memcpy((*temp)->numsteps, zeros, endsite*sizeof(long));
  zeronumnuc(*temp, endsite);
}  /* alloctemp */


void freetemp(node **temp)
{
  /* used in dnacomp, dnapars, & dnapenny */
  free((*temp)->numsteps);
  free((*temp)->base);
  free((*temp)->numnuc);
  free(*temp);
}  /* freetemp */


void freetree2 (pointarray treenode, long nonodes)
{
  /* The natural complement to alloctree2.  Free all elements of all
  the rings (normally triads) in treenode */
  long i;
  node *p, *q;

  /* The first spp elements are just nodes, not rings */
  for (i = 0; i < spp; i++)
    free (treenode[i]);

  /* The rest are rings */
  for (i = spp; i < nonodes; i++) {
    p = treenode[i]->next;
    while (p != treenode[i]) {
      q = p->next;
      free (p);
      p = q;
    }
    /* p should now point to treenode[i], which has yet to be freed */
    free (p);
  }
  free (treenode);
}  /* freetree2 */


void inputdata(long chars)
{
  /* input the names and sequences for each species */
  /* used by dnacomp, dnadist, dnainvar, dnaml, dnamlk, dnapars, & dnapenny */
  long i, j, k, l, basesread, basesnew=0;
  Char charstate;
  boolean allread, done;

  if (printdata)
    headings(chars, "Sequences", "---------");
  basesread = 0;
  allread = false;
  while (!(allread)) {
    /* eat white space -- if the separator line has spaces on it*/
    do {
      charstate = gettc(infile);
    } while (charstate == ' ' || charstate == '\t');
    ungetc(charstate, infile);
    if (eoln(infile))
      scan_eoln(infile);
    i = 1;
    while (i <= spp) {
      if ((interleaved && basesread == 0) || !interleaved)
        initname(i-1);
      j = (interleaved) ? basesread : 0;
      done = false;
      while (!done && !eoff(infile)) {
        if (interleaved)
          done = true;
        while (j < chars && !(eoln(infile) || eoff(infile))) {
          charstate = gettc(infile);
          if (charstate == '\n' || charstate == '\t')
            charstate = ' ';
          if (charstate == ' ' || (charstate >= '0' && charstate <= '9'))
            continue;
          uppercase(&charstate);
          if ((strchr("ABCDGHKMNRSTUVWXY?O-",charstate)) == NULL){
            printf("ERROR: bad base: %c at site %5ld of species %3ld\n",
                   charstate, j+1, i);
            if (charstate == '.') {
              printf("       Periods (.) may not be used as gap characters.\n");
              printf("       The correct gap character is (-)\n");
            }
            exxit(-1);
          }
          j++;
          y[i - 1][j - 1] = charstate;
        }
        if (interleaved)
          continue;
        if (j < chars) 
          scan_eoln(infile);
        else if (j == chars)
          done = true;
      }
      if (interleaved && i == 1)
        basesnew = j;

      scan_eoln(infile);
    
      if ((interleaved && j != basesnew) ||
          (!interleaved && j != chars)) {
        printf("\nERROR: sequences out of alignment at position %ld", j+1);
        printf(" of species %ld\n\n", i);
        exxit(-1);
      }
      i++;
    }
    
    if (interleaved) {
      basesread = basesnew;
      allread = (basesread == chars);
    } else
      allread = (i > spp);
  }
  if (!printdata)
    return;
  for (i = 1; i <= ((chars - 1) / 60 + 1); i++) {
    for (j = 1; j <= spp; j++) {
      for (k = 0; k < nmlngth; k++)
        putc(nayme[j - 1][k], outfile);
      fprintf(outfile, "   ");
      l = i * 60;
      if (l > chars)
        l = chars;
      for (k = (i - 1) * 60 + 1; k <= l; k++) {
        if (dotdiff && (j > 1 && y[j - 1][k - 1] == y[0][k - 1]))
          charstate = '.';
        else
          charstate = y[j - 1][k - 1];
        putc(charstate, outfile);
        if (k % 10 == 0 && k % 60 != 0)
          putc(' ', outfile);
      }
      putc('\n', outfile);
    }
    putc('\n', outfile);
  }
  putc('\n', outfile);
}  /* inputdata */


void alloctree(pointarray *treenode, long nonodes, boolean usertree)
{
  /* allocate treenode dynamically */
  /* used in dnapars, dnacomp, dnapenny & dnamove */
  long i, j;
  node *p, *q;

  *treenode = (pointarray)Malloc(nonodes*sizeof(node *));
  for (i = 0; i < spp; i++) {
    (*treenode)[i] = (node *)Malloc(sizeof(node));
    (*treenode)[i]->tip = true;
    (*treenode)[i]->index = i+1;
    (*treenode)[i]->iter = true;
    (*treenode)[i]->branchnum = 0;
    (*treenode)[i]->initialized = true;
  }
  if (!usertree)
    for (i = spp; i < nonodes; i++) {
      q = NULL;
      for (j = 1; j <= 3; j++) {
        p = (node *)Malloc(sizeof(node));
        p->tip = false;
        p->index = i+1;
        p->iter = true;
        p->branchnum = 0;
        p->initialized = false;
        p->next = q;
        q = p;
      }
      p->next->next->next = p;
      (*treenode)[i] = p;
    }
} /* alloctree */


void allocx(long nonodes, long rcategs, pointarray treenode, boolean usertree)
{
  /* allocate x dynamically */
  /* used in dnaml & dnamlk */
  long i, j, k;
  node *p;

  for (i = 0; i < spp; i++){
    treenode[i]->x = (phenotype)Malloc(endsite*sizeof(ratelike));
    treenode[i]->underflows = (double *)Malloc(endsite * sizeof (double));
    for (j = 0; j < endsite; j++)
      treenode[i]->x[j]  = (ratelike)Malloc(rcategs*sizeof(sitelike));
  }
  if (!usertree) {
    for (i = spp; i < nonodes; i++) {
      p = treenode[i];
      for (j = 1; j <= 3; j++) {
        p->underflows = (double *)Malloc (endsite * sizeof (double));
        p->x = (phenotype)Malloc(endsite*sizeof(ratelike));
        for (k = 0; k < endsite; k++)
          p->x[k] = (ratelike)Malloc(rcategs*sizeof(sitelike));
        p = p->next;
      }
    }
  }
}  /* allocx */


void prot_allocx(long nonodes, long rcategs, pointarray treenode, 
                        boolean usertree)
{
  /* allocate x dynamically */
  /* used in proml          */
  long i, j, k;
  node *p;

  for (i = 0; i < spp; i++){
    treenode[i]->protx = (pphenotype)Malloc(endsite*sizeof(pratelike));
    treenode[i]->underflows = (double *)Malloc(endsite*sizeof(double));
    for (j = 0; j < endsite; j++)
      treenode[i]->protx[j]  = (pratelike)Malloc(rcategs*sizeof(psitelike));
  }  
  if (!usertree) {
    for (i = spp; i < nonodes; i++) {
      p = treenode[i];
      for (j = 1; j <= 3; j++) {
        p->protx = (pphenotype)Malloc(endsite*sizeof(pratelike));
        p->underflows = (double *)Malloc(endsite*sizeof(double));
        for (k = 0; k < endsite; k++)
          p->protx[k] = (pratelike)Malloc(rcategs*sizeof(psitelike));
        p = p->next;
      } 
    }  
  } 
}  /* prot_allocx */




void setuptree(pointarray treenode, long nonodes, boolean usertree)
{
  /* initialize treenodes */
  long i;
  node *p;

  for (i = 1; i <= nonodes; i++) {
    if (i <= spp || !usertree) {
      treenode[i-1]->back = NULL;
      treenode[i-1]->tip = (i <= spp);
      treenode[i-1]->index = i;
      treenode[i-1]->numdesc = 0;
      treenode[i-1]->iter = true;
      treenode[i-1]->initialized = true;
      treenode[i-1]->tyme =  0.0;
    }
  }
  if (!usertree) {
    for (i = spp + 1; i <= nonodes; i++) {
      p = treenode[i-1]->next;
      while (p != treenode[i-1]) {
        p->back = NULL;
        p->tip = false;
        p->index = i;
        p->numdesc = 0;
        p->iter = true;
        p->initialized = false;
        p->tyme = 0.0;
        p = p->next;
      }
    }
  }
} /* setuptree */


void setuptree2(tree *a)
{
  /* initialize a tree */
  /* used in dnaml, dnamlk, & restml */

  a->likelihood = -999999.0;
  a->start = a->nodep[0]->back;
  a->root = NULL;
} /* setuptree2 */


void alloctip(node *p, long *zeros)
{ /* allocate a tip node */
  /* used by dnacomp, dnapars, & dnapenny */

  p->numsteps = (steptr)Malloc(endsite*sizeof(long));
  p->oldnumsteps = (steptr)Malloc(endsite*sizeof(long));
  p->base = (baseptr)Malloc(endsite*sizeof(long));
  p->oldbase = (baseptr)Malloc(endsite*sizeof(long));
  memcpy(p->base, zeros, endsite*sizeof(long));
  memcpy(p->numsteps, zeros, endsite*sizeof(long));
  memcpy(p->oldbase, zeros, endsite*sizeof(long));
  memcpy(p->oldnumsteps, zeros, endsite*sizeof(long));
}  /* alloctip */




void getbasefreqs(double freqa, double freqc, double freqg, double freqt,
            double *freqr, double *freqy, double *freqar, double *freqcy,
            double *freqgr, double *freqty, double *ttratio, double *xi,
            double *xv, double *fracchange, boolean freqsfrom,
            boolean printdata)
{
  /* used by dnadist, dnaml, & dnamlk */
  double aa, bb;

  if (printdata) {
    putc('\n', outfile);
    if (freqsfrom)
      fprintf(outfile, "Empirical ");
    fprintf(outfile, "Base Frequencies:\n\n");
    fprintf(outfile, "   A    %10.5f\n", freqa);
    fprintf(outfile, "   C    %10.5f\n", freqc);
    fprintf(outfile, "   G    %10.5f\n", freqg);
    fprintf(outfile, "  T(U)  %10.5f\n", freqt);
    fprintf(outfile, "\n");
  }
  *freqr = freqa + freqg;
  *freqy = freqc + freqt;
  *freqar = freqa / *freqr;
  *freqcy = freqc / *freqy;
  *freqgr = freqg / *freqr;
  *freqty = freqt / *freqy;
  aa = *ttratio * (*freqr) * (*freqy) - freqa * freqg - freqc * freqt;
  bb = freqa * (*freqgr) + freqc * (*freqty);
  *xi = aa / (aa + bb);
  *xv = 1.0 - *xi;
  if (*xi < 0.0) {
    printf("\n WARNING: This transition/transversion ratio\n");
    printf(" is impossible with these base frequencies!\n");
    *xi = 0.0;
    *xv = 1.0;
    (*ttratio) = (freqa*freqg+freqc*freqt)/((*freqr)*(*freqy));
    printf(" Transition/transversion parameter reset\n");
    printf("  so transition/transversion ratio is %10.6f\n\n", (*ttratio));
  }
  if (freqa <= 0.0)
    freqa = 0.000001;
  if (freqc <= 0.0)
    freqc = 0.000001;
  if (freqg <= 0.0)
    freqg = 0.000001;
  if (freqt <= 0.0)
    freqt = 0.000001;
  *fracchange = (*xi) * (2 * freqa * (*freqgr) + 2 * freqc * (*freqty)) +
      (*xv) * (1.0 - freqa * freqa - freqc * freqc - freqg * freqg
      - freqt * freqt);
}  /* getbasefreqs */


void empiricalfreqs(double *freqa, double *freqc, double *freqg,
                        double *freqt, steptr weight, pointarray treenode)
{
  /* Get empirical base frequencies from the data */
  /* used in dnaml & dnamlk */
  long i, j, k;
  double sum, suma, sumc, sumg, sumt, w;

  *freqa = 0.25;
  *freqc = 0.25;
  *freqg = 0.25;
  *freqt = 0.25;
  for (k = 1; k <= 8; k++) {
    suma = 0.0;
    sumc = 0.0;
    sumg = 0.0;
    sumt = 0.0;
    for (i = 0; i < spp; i++) {
      for (j = 0; j < endsite; j++) {
        w = weight[j];
        sum = (*freqa) * treenode[i]->x[j][0][0];
        sum += (*freqc) * treenode[i]->x[j][0][(long)C - (long)A];
        sum += (*freqg) * treenode[i]->x[j][0][(long)G - (long)A];
        sum += (*freqt) * treenode[i]->x[j][0][(long)T - (long)A];
        suma += w * (*freqa) * treenode[i]->x[j][0][0] / sum;
        sumc += w * (*freqc) * treenode[i]->x[j][0][(long)C - (long)A] / sum;
        sumg += w * (*freqg) * treenode[i]->x[j][0][(long)G - (long)A] / sum;
        sumt += w * (*freqt) * treenode[i]->x[j][0][(long)T - (long)A] / sum;
      }
    }
    sum = suma + sumc + sumg + sumt;
    *freqa = suma / sum;
    *freqc = sumc / sum;
    *freqg = sumg / sum;
    *freqt = sumt / sum;
  }
  if (*freqa <= 0.0)
    *freqa = 0.000001;
  if (*freqc <= 0.0)
    *freqc = 0.000001;
  if (*freqg <= 0.0)
    *freqg = 0.000001;
  if (*freqt <= 0.0)
    *freqt = 0.000001;
}  /* empiricalfreqs */


void sitesort(long chars, steptr weight)
{
  /* Shell sort keeping sites, weights in same order */
  /* used in dnainvar, dnapars, dnacomp & dnapenny */
  long gap, i, j, jj, jg, k, itemp;
  boolean flip, tied;

  gap = chars / 2;
  while (gap > 0) {
    for (i = gap + 1; i <= chars; i++) {
      j = i - gap;
      flip = true;
      while (j > 0 && flip) {
        jj = alias[j - 1];
        jg = alias[j + gap - 1];
        tied = true;
        k = 1;
        while (k <= spp && tied) {
          flip = (y[k - 1][jj - 1] > y[k - 1][jg - 1]);
          tied = (tied && y[k - 1][jj - 1] == y[k - 1][jg - 1]);
          k++;
        }
        if (!flip)
          break;
        itemp = alias[j - 1];
        alias[j - 1] = alias[j + gap - 1];
        alias[j + gap - 1] = itemp;
        itemp = weight[j - 1];
        weight[j - 1] = weight[j + gap - 1];
        weight[j + gap - 1] = itemp;
        j -= gap;
      }
    }
    gap /= 2;
  }
}  /* sitesort */


void sitecombine(long chars)
{
  /* combine sites that have identical patterns */
  /* used in dnapars, dnapenny, & dnacomp */
  long i, j, k;
  boolean tied;

  i = 1;
  while (i < chars) {
    j = i + 1;
    tied = true;
    while (j <= chars && tied) {
      k = 1;
      while (k <= spp && tied) {
        tied = (tied &&
            y[k - 1][alias[i - 1] - 1] == y[k - 1][alias[j - 1] - 1]);
        k++;
      }
      if (tied) {
        weight[i - 1] += weight[j - 1];
        weight[j - 1] = 0;
        ally[alias[j - 1] - 1] = alias[i - 1];
      }
      j++;
    }
    i = j - 1;
  }
}  /* sitecombine */


void sitescrunch(long chars)
{
  /* move so one representative of each pattern of
     sites comes first */
  /* used in dnapars & dnacomp */
  long i, j, itemp;
  boolean done, found;

  done = false;
  i = 1;
  j = 2;
  while (!done) {
    if (ally[alias[i - 1] - 1] != alias[i - 1]) {
      if (j <= i)
        j = i + 1;
      if (j <= chars) {
        do {
          found = (ally[alias[j - 1] - 1] == alias[j - 1]);
          j++;
        } while (!(found || j > chars));
        if (found) {
          j--;
          itemp = alias[i - 1];
          alias[i - 1] = alias[j - 1];
          alias[j - 1] = itemp;
          itemp = weight[i - 1];
          weight[i - 1] = weight[j - 1];
          weight[j - 1] = itemp;
        } else
          done = true;
      } else
        done = true;
    }
    i++;
    done = (done || i >= chars);
  }
}  /* sitescrunch */


void sitesort2(long sites, steptr aliasweight)
{
  /* Shell sort keeping sites, weights in same order */
  /* used in dnaml & dnamnlk */
  long gap, i, j, jj, jg, k, itemp;
  boolean flip, tied;

  gap = sites / 2;
  while (gap > 0) {
    for (i = gap + 1; i <= sites; i++) {
      j = i - gap;
      flip = true;
      while (j > 0 && flip) {
        jj = alias[j - 1];
        jg = alias[j + gap - 1];
        tied = (category[jj - 1] == category[jg - 1]);
        flip = (category[jj - 1] > category[jg - 1]);
        k = 1;
        while (k <= spp && tied) {
          flip = (y[k - 1][jj - 1] > y[k - 1][jg - 1]);
          tied = (tied && y[k - 1][jj - 1] == y[k - 1][jg - 1]);
          k++;
        }
        if (!flip)
          break;
        itemp = alias[j - 1];
        alias[j - 1] = alias[j + gap - 1];
        alias[j + gap - 1] = itemp;
        itemp = aliasweight[j - 1];
        aliasweight[j - 1] = aliasweight[j + gap - 1];
        aliasweight[j + gap - 1] = itemp;
        j -= gap;
      }
    }
    gap /= 2;
  }
}  /* sitesort2 */


void sitecombine2(long sites, steptr aliasweight)
{
  /* combine sites that have identical patterns */
  /* used in dnaml & dnamlk */
  long i, j, k;
  boolean tied;

  i = 1;
  while (i < sites) {
    j = i + 1;
    tied = true;
    while (j <= sites && tied) {
      tied = (category[alias[i - 1] - 1] == category[alias[j - 1] - 1]);
      k = 1;
      while (k <= spp && tied) {
        tied = (tied &&
            y[k - 1][alias[i - 1] - 1] == y[k - 1][alias[j - 1] - 1]);
        k++;
      }
      if (!tied)
        break;
      aliasweight[i - 1] += aliasweight[j - 1];
      aliasweight[j - 1] = 0;
      ally[alias[j - 1] - 1] = alias[i - 1];
      j++;
    }
    i = j;
  }
}  /* sitecombine2 */


void sitescrunch2(long sites, long i, long j, steptr aliasweight)
{
  /* move so positively weighted sites come first */
  /* used by dnainvar, dnaml, dnamlk, & restml */
  long itemp;
  boolean done, found;

  done = false;
  while (!done) {
    if (aliasweight[i - 1] > 0)
      i++;
    else {
      if (j <= i)
        j = i + 1;
      if (j <= sites) {
        do {
          found = (aliasweight[j - 1] > 0);
          j++;
        } while (!(found || j > sites));
        if (found) {
          j--;
          itemp = alias[i - 1];
          alias[i - 1] = alias[j - 1];
          alias[j - 1] = itemp;
          itemp = aliasweight[i - 1];
          aliasweight[i - 1] = aliasweight[j - 1];
          aliasweight[j - 1] = itemp;
        } else
          done = true;
      } else
        done = true;
    }
    done = (done || i >= sites);
  }
}  /* sitescrunch2 */


void makevalues(pointarray treenode, long *zeros, boolean usertree)
{
  /* set up fractional likelihoods at tips */
  /* used by dnacomp, dnapars, & dnapenny */
  long i, j;
  char ns = 0;
  node *p;

  setuptree(treenode, nonodes, usertree);
  for (i = 0; i < spp; i++)
    alloctip(treenode[i], zeros);
  if (!usertree) {
    for (i = spp; i < nonodes; i++) {
      p = treenode[i];
      do {
        allocnontip(p, zeros, endsite);
        p = p->next;
      } while (p != treenode[i]);
    }
  }
  for (j = 0; j < endsite; j++) {
    for (i = 0; i < spp; i++) {
      switch (y[i][alias[j] - 1]) {

      case 'A':
        ns = 1 << A;
        break;

      case 'C':
        ns = 1 << C;
        break;

      case 'G':
        ns = 1 << G;
        break;

      case 'U':
        ns = 1 << T;
        break;

      case 'T':
        ns = 1 << T;
        break;

      case 'M':
        ns = (1 << A) | (1 << C);
        break;

      case 'R':
        ns = (1 << A) | (1 << G);
        break;

      case 'W':
        ns = (1 << A) | (1 << T);
        break;

      case 'S':
        ns = (1 << C) | (1 << G);
        break;

      case 'Y':
        ns = (1 << C) | (1 << T);
        break;

      case 'K':
        ns = (1 << G) | (1 << T);
        break;

      case 'B':
        ns = (1 << C) | (1 << G) | (1 << T);
        break;

      case 'D':
        ns = (1 << A) | (1 << G) | (1 << T);
        break;

      case 'H':
        ns = (1 << A) | (1 << C) | (1 << T);
        break;

      case 'V':
        ns = (1 << A) | (1 << C) | (1 << G);
        break;

      case 'N':
        ns = (1 << A) | (1 << C) | (1 << G) | (1 << T);
        break;

      case 'X':
        ns = (1 << A) | (1 << C) | (1 << G) | (1 << T);
        break;

      case '?':
        ns = (1 << A) | (1 << C) | (1 << G) | (1 << T) | (1 << O);
        break;

      case 'O':
        ns = 1 << O;
        break;

      case '-':
        ns = 1 << O;
        break;
      }
      treenode[i]->base[j] = ns;
      treenode[i]->numsteps[j] = 0;
    }
  }
}  /* makevalues */


void makevalues2(long categs, pointarray treenode, long endsite,
                        long spp, sequence y, steptr alias)
{
  /* set up fractional likelihoods at tips */
  /* used by dnaml & dnamlk */
  long i, j, k, l;
  bases b;

  for (k = 0; k < endsite; k++) {
    j = alias[k];
    for (i = 0; i < spp; i++) {
      for (l = 0; l < categs; l++) {
        for (b = A; (long)b <= (long)T; b = (bases)((long)b + 1))
          treenode[i]->x[k][l][(long)b - (long)A] = 0.0;
        switch (y[i][j - 1]) {

        case 'A':
          treenode[i]->x[k][l][0] = 1.0;
          break;

        case 'C':
          treenode[i]->x[k][l][(long)C - (long)A] = 1.0;
          break;

        case 'G':
          treenode[i]->x[k][l][(long)G - (long)A] = 1.0;
          break;

        case 'T':
          treenode[i]->x[k][l][(long)T - (long)A] = 1.0;
          break;

        case 'U':
          treenode[i]->x[k][l][(long)T - (long)A] = 1.0;
          break;

        case 'M':
          treenode[i]->x[k][l][0] = 1.0;
          treenode[i]->x[k][l][(long)C - (long)A] = 1.0;
          break;

        case 'R':
          treenode[i]->x[k][l][0] = 1.0;
          treenode[i]->x[k][l][(long)G - (long)A] = 1.0;
          break;

        case 'W':
          treenode[i]->x[k][l][0] = 1.0;
          treenode[i]->x[k][l][(long)T - (long)A] = 1.0;
          break;

        case 'S':
          treenode[i]->x[k][l][(long)C - (long)A] = 1.0;
          treenode[i]->x[k][l][(long)G - (long)A] = 1.0;
          break;

        case 'Y':
          treenode[i]->x[k][l][(long)C - (long)A] = 1.0;
          treenode[i]->x[k][l][(long)T - (long)A] = 1.0;
          break;

        case 'K':
          treenode[i]->x[k][l][(long)G - (long)A] = 1.0;
          treenode[i]->x[k][l][(long)T - (long)A] = 1.0;
          break;

        case 'B':
          treenode[i]->x[k][l][(long)C - (long)A] = 1.0;
          treenode[i]->x[k][l][(long)G - (long)A] = 1.0;
          treenode[i]->x[k][l][(long)T - (long)A] = 1.0;
          break;

        case 'D':
          treenode[i]->x[k][l][0] = 1.0;
          treenode[i]->x[k][l][(long)G - (long)A] = 1.0;
          treenode[i]->x[k][l][(long)T - (long)A] = 1.0;
          break;

        case 'H':
          treenode[i]->x[k][l][0] = 1.0;
          treenode[i]->x[k][l][(long)C - (long)A] = 1.0;
          treenode[i]->x[k][l][(long)T - (long)A] = 1.0;
          break;

        case 'V':
          treenode[i]->x[k][l][0] = 1.0;
          treenode[i]->x[k][l][(long)C - (long)A] = 1.0;
          treenode[i]->x[k][l][(long)G - (long)A] = 1.0;
          break;

        case 'N':
          for (b = A; (long)b <= (long)T; b = (bases)((long)b + 1))
            treenode[i]->x[k][l][(long)b - (long)A] = 1.0;
          break;

        case 'X':
          for (b = A; (long)b <= (long)T; b = (bases)((long)b + 1))
            treenode[i]->x[k][l][(long)b - (long)A] = 1.0;
          break;

        case '?':
          for (b = A; (long)b <= (long)T; b = (bases)((long)b + 1))
            treenode[i]->x[k][l][(long)b - (long)A] = 1.0;
          break;

        case 'O':
          for (b = A; (long)b <= (long)T; b = (bases)((long)b + 1))
            treenode[i]->x[k][l][(long)b - (long)A] = 1.0;
          break;

        case '-':
          for (b = A; (long)b <= (long)T; b = (bases)((long)b + 1))
            treenode[i]->x[k][l][(long)b - (long)A] = 1.0;
          break;
        }
      }
    }
  }
}  /* makevalues2 */


void fillin(node *p, node *left, node *rt)
{
  /* sets up for each node in the tree the base sequence
     at that point and counts the changes.  */
  long i, j, k, n, purset, pyrset;
  node *q;

  purset = (1 << (long)A) + (1 << (long)G);
  pyrset = (1 << (long)C) + (1 << (long)T);
  if (!left) {
    memcpy(p->base, rt->base, endsite*sizeof(long));
    memcpy(p->numsteps, rt->numsteps, endsite*sizeof(long));
    q = rt;
  } else if (!rt) {
    memcpy(p->base, left->base, endsite*sizeof(long));
    memcpy(p->numsteps, left->numsteps, endsite*sizeof(long));
    q = left;
  } else {
    for (i = 0; i < endsite; i++) {
      p->base[i] = left->base[i] & rt->base[i];
      p->numsteps[i] = left->numsteps[i] + rt->numsteps[i];
      if (p->base[i] == 0) {
        p->base[i] = left->base[i] | rt->base[i];
        if (transvp) {
          if (!((p->base[i] == purset) || (p->base[i] == pyrset)))
            p->numsteps[i] += weight[i];
        }
        else p->numsteps[i] += weight[i];
      }
    }
    q = rt;
  }
  if (left && rt) n = 2;
  else n = 1;
  for (i = 0; i < endsite; i++)
    for (j = (long)A; j <= (long)O; j++)
      p->numnuc[i][j] = 0;
  for (k = 1; k <= n; k++) {
    if (k == 2) q = left;
    for (i = 0; i < endsite; i++) {
      for (j = (long)A; j <= (long)O; j++) {
        if (q->base[i] & (1 << j))
          p->numnuc[i][j]++;
      }
    }
  }
}  /* fillin */


long getlargest(long *numnuc)
{
  /* find the largest in array numnuc */
  long i, largest;

  largest = 0;
  for (i = (long)A; i <= (long)O; i++)
    if (numnuc[i] > largest)
      largest = numnuc[i];
  return largest;
} /* getlargest */


void multifillin(node *p, node *q, long dnumdesc)
{
  /* sets up for each node in the tree the base sequence
     at that point and counts the changes according to the
     changes in q's base */
  long i, j, b, largest, descsteps, purset, pyrset;

  memcpy(p->oldbase, p->base, endsite*sizeof(long));
  memcpy(p->oldnumsteps, p->numsteps, endsite*sizeof(long));
  purset = (1 << (long)A) + (1 << (long)G);
  pyrset = (1 << (long)C) + (1 << (long)T);
  for (i = 0; i < endsite; i++) {
    descsteps = 0;
    for (j = (long)A; j <= (long)O; j++) {
      b = 1 << j;
      if ((descsteps == 0) && (p->base[i] & b)) 
        descsteps = p->numsteps[i] 
                      - (p->numdesc - dnumdesc - p->numnuc[i][j]) * weight[i];
    }
    if (dnumdesc == -1)
      descsteps -= q->oldnumsteps[i];
    else if (dnumdesc == 0)
      descsteps += (q->numsteps[i] - q->oldnumsteps[i]);
    else
      descsteps += q->numsteps[i];
    if (q->oldbase[i] != q->base[i]) {
      for (j = (long)A; j <= (long)O; j++) {
        b = 1 << j;
        if (transvp) {
          if (b & purset) b = purset;
          if (b & pyrset) b = pyrset;
        }
        if ((q->oldbase[i] & b) && !(q->base[i] & b))
          p->numnuc[i][j]--;
        else if (!(q->oldbase[i] & b) && (q->base[i] & b))
          p->numnuc[i][j]++;
      }
    }
    largest = getlargest(p->numnuc[i]);
    if (q->oldbase[i] != q->base[i]) {
      p->base[i] = 0;
      for (j = (long)A; j <= (long)O; j++) {
        if (p->numnuc[i][j] == largest)
            p->base[i] |= (1 << j);
      }
    }
    p->numsteps[i] = (p->numdesc - largest) * weight[i] + descsteps;
  }
} /* multifillin */


void sumnsteps(node *p, node *left, node *rt, long a, long b)
{
  /* sets up for each node in the tree the base sequence
     at that point and counts the changes. */
  long i;
  long ns, rs, ls, purset, pyrset;

  if (!left) {
    memcpy(p->numsteps, rt->numsteps, endsite*sizeof(long));
    memcpy(p->base, rt->base, endsite*sizeof(long));
  } else if (!rt) {
    memcpy(p->numsteps, left->numsteps, endsite*sizeof(long));
    memcpy(p->base, left->base, endsite*sizeof(long));
  } else  {
    purset = (1 << (long)A) + (1 << (long)G);
    pyrset = (1 << (long)C) + (1 << (long)T);
    for (i = a; i < b; i++) {
      ls = left->base[i];
      rs = rt->base[i];
      ns = ls & rs;
      p->numsteps[i] = left->numsteps[i] + rt->numsteps[i];
      if (ns == 0) {
        ns = ls | rs;
        if (transvp) {
          if (!((ns == purset) || (ns == pyrset)))
            p->numsteps[i] += weight[i];
        }
        else p->numsteps[i] += weight[i];
      }
      p->base[i] = ns;
      }
    }
}  /* sumnsteps */


void sumnsteps2(node *p,node *left,node *rt,long a,long b,long *threshwt)
{
  /* counts the changes at each node.  */
  long i, steps;
  long ns, rs, ls, purset, pyrset;
  long term;

  if (a == 0) p->sumsteps = 0.0;
  if (!left)
    memcpy(p->numsteps, rt->numsteps, endsite*sizeof(long));
  else if (!rt)
    memcpy(p->numsteps, left->numsteps, endsite*sizeof(long));
  else {
    purset = (1 << (long)A) + (1 << (long)G);
    pyrset = (1 << (long)C) + (1 << (long)T);
    for (i = a; i < b; i++) {
      ls = left->base[i];
      rs = rt->base[i];
      ns = ls & rs;
      p->numsteps[i] = left->numsteps[i] + rt->numsteps[i];
      if (ns == 0) {
        ns = ls | rs;
        if (transvp) {
          if (!((ns == purset) || (ns == pyrset)))
            p->numsteps[i] += weight[i];
        }
        else p->numsteps[i] += weight[i];
      }
    }
  }
  for (i = a; i < b; i++) {
    steps = p->numsteps[i];
    if ((long)steps <= threshwt[i])
      term = steps;
    else
      term = threshwt[i];
    p->sumsteps += (double)term;
  }
}  /* sumnsteps2 */


void multisumnsteps(node *p, node *q, long a, long b, long *threshwt)
{
  /* computes the number of steps between p and q */
  long i, j, steps, largest, descsteps, purset, pyrset, b1;
  long term;

  if (a == 0) p->sumsteps = 0.0;
  purset = (1 << (long)A) + (1 << (long)G);
  pyrset = (1 << (long)C) + (1 << (long)T);
  for (i = a; i < b; i++) {
    descsteps = 0;
    for (j = (long)A; j <= (long)O; j++) {
      if ((descsteps == 0) && (p->base[i] & (1 << j))) 
        descsteps = p->numsteps[i] -
                        (p->numdesc - 1 - p->numnuc[i][j]) * weight[i];
    }
    descsteps += q->numsteps[i];
    largest = 0;
    for (j = (long)A; j <= (long)O; j++) {
      b1 = (1 << j);
      if (transvp) {
        if (b1 & purset) b1 = purset;
        if (b1 & pyrset) b1 = pyrset;
      }
      if (q->base[i] & b1)
        p->numnuc[i][j]++;
      if (p->numnuc[i][j] > largest)
        largest = p->numnuc[i][j];
    }
    steps = (p->numdesc - largest) * weight[i] + descsteps;
    if ((long)steps <= threshwt[i])
      term = steps;
    else
      term = threshwt[i];
    p->sumsteps += (double)term;
  }
} /* multisumnsteps */


void multisumnsteps2(node *p)
{
  /* counts the changes at each multi-way node. Sums up
     steps of all descendants */
  long i, j, largest, purset, pyrset, b1;
  node *q;
  baseptr b;

  purset = (1 << (long)A) + (1 << (long)G);
  pyrset = (1 << (long)C) + (1 << (long)T);
  for (i = 0; i < endsite; i++) {
    p->numsteps[i] = 0;
    q = p->next;
    while (q != p) {
      if (q->back) {
        p->numsteps[i] += q->back->numsteps[i];
        b = q->back->base;
        for (j = (long)A; j <= (long)O; j++) {
          b1 = (1 << j);   
          if (transvp) {
            if (b1 & purset) b1 = purset;
            if (b1 & pyrset) b1 = pyrset;
          }
          if (b[i] & b1) p->numnuc[i][j]++;
        }
      }
      q = q->next;
    }
    largest = getlargest(p->numnuc[i]);
    p->base[i] = 0;
    for (j = (long)A; j <= (long)O; j++) {
      if (p->numnuc[i][j] == largest)
        p->base[i] |= (1 << j);
    }
    p->numsteps[i] += ((p->numdesc - largest) * weight[i]);
  }
}  /* multisumnsteps2 */

boolean alltips(node *forknode, node *p)
{
  /* returns true if all descendants of forknode except p are tips; 
     false otherwise.  */
  node *q, *r;
  boolean tips;

  tips = true;
  r = forknode;
  q = forknode->next;
  do {
    if (q->back && q->back != p && !q->back->tip)
      tips = false;
    q = q->next;
  } while (tips && q != r);
  return tips;
} /* alltips */


void gdispose(node *p, node **grbg, pointarray treenode)
{
  /* go through tree throwing away nodes */
  node *q, *r;

  p->back = NULL;
  if (p->tip)
    return;
  treenode[p->index - 1] = NULL;
  q = p->next;
  while (q != p) {
    gdispose(q->back, grbg, treenode);
    q->back = NULL;
    r = q;
    q = q->next;
    chuck(grbg, r);
  }
  chuck(grbg, q);
}  /* gdispose */


void preorder(node *p, node *r, node *root, node *removing, node *adding,
                        node *changing, long dnumdesc)
{
  /* recompute number of steps in preorder taking both ancestoral and
     descendent steps into account. removing points to a node being 
     removed, if any */
  node *q, *p1, *p2;

  if (p && !p->tip && p != adding) {
    q = p;
    do {
      if (p->back != r) {
        if (p->numdesc > 2) {
          if (changing)
            multifillin (p, r, dnumdesc);
          else
            multifillin (p, r, 0);
        } else {
          p1 = p->next;
          if (!removing)
            while (!p1->back)
              p1 = p1->next;
          else
            while (!p1->back || p1->back == removing)
              p1 = p1->next;
          p2 = p1->next;
          if (!removing)
            while (!p2->back)
              p2 = p2->next;
          else
            while (!p2->back || p2->back == removing)
              p2 = p2->next;
          p1 = p1->back;
          p2 = p2->back;
          if (p->back == p1) p1 = NULL;
          else if (p->back == p2) p2 = NULL;
          memcpy(p->oldbase, p->base, endsite*sizeof(long));
          memcpy(p->oldnumsteps, p->numsteps, endsite*sizeof(long));
          fillin(p, p1, p2);
        }
      }
      p = p->next;
    } while (p != q);
    q = p;
    do {
      preorder(p->next->back, p->next, root, removing, adding, NULL, 0);
      p = p->next;
    } while (p->next != q);
  }
} /* preorder */


void updatenumdesc(node *p, node *root, long n)
{
  /* set p's numdesc to n.  If p is the root, numdesc of p's
  descendants are set to n-1. */
  node *q;

  q = p;
  if (p == root && n > 0) {
    p->numdesc = n;
    n--;
    q = q->next;
  }
  do {
    q->numdesc = n;
    q = q->next;
  } while (q != p);
}  /* updatenumdesc */


void add(node *below,node *newtip,node *newfork,node **root,
        boolean recompute,pointarray treenode,node **grbg,long *zeros)
{
  /* inserts the nodes newfork and its left descendant, newtip,
     to the tree.  below becomes newfork's right descendant.
     if newfork is NULL, newtip is added as below's sibling */
  /* used in dnacomp & dnapars */
  node *p;

  if (below != treenode[below->index - 1])
    below = treenode[below->index - 1];
  if (newfork) {
    if (below->back != NULL)
      below->back->back = newfork;
    newfork->back = below->back;
    below->back = newfork->next->next;
    newfork->next->next->back = below;
    newfork->next->back = newtip;
    newtip->back = newfork->next;
    if (*root == below)
      *root = newfork;
    updatenumdesc(newfork, *root, 2);
  } else {
    gnutreenode(grbg, &p, below->index, endsite, zeros);
    p->back = newtip;
    newtip->back = p;
    p->next = below->next;
    below->next = p;
    updatenumdesc(below, *root, below->numdesc + 1);
  }
  if (!newtip->tip)
    updatenumdesc(newtip, *root, newtip->numdesc);
  (*root)->back = NULL;
  if (!recompute)
    return;
  if (!newfork) {
    memcpy(newtip->back->base, below->base, endsite*sizeof(long));
    memcpy(newtip->back->numsteps, below->numsteps, endsite*sizeof(long));
    memcpy(newtip->back->numnuc, below->numnuc, endsite*sizeof(nucarray));
    if (below != *root) {
      memcpy(below->back->oldbase, zeros, endsite*sizeof(long));
      memcpy(below->back->oldnumsteps, zeros, endsite*sizeof(long));
      multifillin(newtip->back, below->back, 1);
    }
    if (!newtip->tip) {
      memcpy(newtip->back->oldbase, zeros, endsite*sizeof(long));
      memcpy(newtip->back->oldnumsteps, zeros, endsite*sizeof(long));
      preorder(newtip, newtip->back, *root, NULL, NULL, below, 1);
    }
    memcpy(newtip->oldbase, zeros, endsite*sizeof(long));
    memcpy(newtip->oldnumsteps, zeros, endsite*sizeof(long));
    preorder(below, newtip, *root, NULL, newtip, below, 1);
    if (below != *root)
      preorder(below->back, below, *root, NULL, NULL, NULL, 0);
  } else {
    fillin(newtip->back, newtip->back->next->back,
             newtip->back->next->next->back);
    if (!newtip->tip) {
      memcpy(newtip->back->oldbase, zeros, endsite*sizeof(long));
      memcpy(newtip->back->oldnumsteps, zeros, endsite*sizeof(long));
      preorder(newtip, newtip->back, *root, NULL, NULL, newfork, 1);
    }
    if (newfork != *root) {
      memcpy(below->back->base, newfork->back->base, endsite*sizeof(long));
      memcpy(below->back->numsteps, newfork->back->numsteps, endsite*sizeof(long));
      preorder(newfork, newtip, *root, NULL, newtip, NULL, 0);
    } else {
      fillin(below->back, newtip, NULL);
      fillin(newfork, newtip, below);
      memcpy(below->back->oldbase, zeros, endsite*sizeof(long));
      memcpy(below->back->oldnumsteps, zeros, endsite*sizeof(long));
      preorder(below, below->back, *root, NULL, NULL, newfork, 1);
    }
    if (newfork != *root) {
      memcpy(newfork->oldbase, below->base, endsite*sizeof(long));
      memcpy(newfork->oldnumsteps, below->numsteps, endsite*sizeof(long));
      preorder(newfork->back, newfork, *root, NULL, NULL, NULL, 0);
    }
  }
}  /* add */


void findbelow(node **below, node *item, node *fork)
{
  /* decide which of fork's binary children is below */

  if (fork->next->back == item)
    *below = fork->next->next->back;
  else
    *below = fork->next->back;
} /* findbelow */


void re_move(node *item, node **fork, node **root, boolean recompute,
                        pointarray treenode, node **grbg, long *zeros)
{
  /* removes nodes item and its ancestor, fork, from the tree.
     the new descendant of fork's ancestor is made to be
     fork's second descendant (other than item).  Also
     returns pointers to the deleted nodes, item and fork.
     If item belongs to a node with more than 2 descendants,
     fork will not be deleted */
  /* used in dnacomp & dnapars */
  node *p, *q, *other = NULL, *otherback = NULL;

  if (item->back == NULL) {
    *fork = NULL;
    return;
  }
  *fork = treenode[item->back->index - 1];
  if ((*fork)->numdesc == 2) {
    updatenumdesc(*fork, *root, 0);
    findbelow(&other, item, *fork);
    otherback = other->back;
    if (*root == *fork) {
      *root = other;
      if (!other->tip)
        updatenumdesc(other, *root, other->numdesc);
    }
    p = item->back->next->back;
    q = item->back->next->next->back;
    if (p != NULL)
      p->back = q;
    if (q != NULL)
      q->back = p;
    (*fork)->back = NULL;
    p = (*fork)->next;
    while (p != *fork) {
      p->back = NULL;
      p = p->next;
    }
  } else {
    updatenumdesc(*fork, *root, (*fork)->numdesc - 1);
    p = *fork;
    while (p->next != item->back)
      p = p->next;
    p->next = item->back->next;
  }
  if (!item->tip) {
    updatenumdesc(item, item, item->numdesc);
    if (recompute) {
      memcpy(item->back->oldbase, item->back->base, endsite*sizeof(long));
      memcpy(item->back->oldnumsteps, item->back->numsteps, endsite*sizeof(long));
      memcpy(item->back->base, zeros, endsite*sizeof(long));
      memcpy(item->back->numsteps, zeros, endsite*sizeof(long));
      preorder(item, item->back, *root, item->back, NULL, item, -1);
    }
  }
  if ((*fork)->numdesc >= 2)
    chuck(grbg, item->back);
  item->back = NULL;
  if (!recompute)
    return;
  if ((*fork)->numdesc == 0) {
    memcpy(otherback->oldbase, otherback->base, endsite*sizeof(long));  
    memcpy(otherback->oldnumsteps, otherback->numsteps, endsite*sizeof(long));
    if (other == *root) {
      memcpy(otherback->base, zeros, endsite*sizeof(long));
      memcpy(otherback->numsteps, zeros, endsite*sizeof(long));
    } else {
      memcpy(otherback->base, other->back->base, endsite*sizeof(long));
      memcpy(otherback->numsteps, other->back->numsteps, endsite*sizeof(long));
    }
    p = other->back;
    other->back = otherback;
    if (other == *root)
      preorder(other, otherback, *root, otherback, NULL, other, -1);
    else
      preorder(other, otherback, *root, NULL, NULL, NULL, 0);
    other->back = p;
    if (other != *root) {
      memcpy(other->oldbase,(*fork)->base, endsite*sizeof(long));
      memcpy(other->oldnumsteps,(*fork)->numsteps, endsite*sizeof(long));
      preorder(other->back, other, *root, NULL, NULL, NULL, 0);
    }
  } else {
    memcpy(item->oldbase, item->base, endsite*sizeof(long));
    memcpy(item->oldnumsteps, item->numsteps, endsite*sizeof(long));
    memcpy(item->base, zeros, endsite*sizeof(long));
    memcpy(item->numsteps, zeros, endsite*sizeof(long));
    preorder(*fork, item, *root, NULL, NULL, *fork, -1);
    if (*fork != *root)
      preorder((*fork)->back, *fork, *root, NULL, NULL, NULL, 0);
    memcpy(item->base, item->oldbase, endsite*sizeof(long));
    memcpy(item->numsteps, item->oldnumsteps, endsite*sizeof(long));
  }
}  /* remove */


void postorder(node *p)
{
  /* traverses an n-ary tree, suming up steps at a node's descendants */
  /* used in dnacomp, dnapars, & dnapenny */
  node *q;

  if (p->tip)
    return;
  q = p->next;
  while (q != p) {
    postorder(q->back);
    q = q->next;
  }
  zeronumnuc(p, endsite);
  if (p->numdesc > 2)
    multisumnsteps2(p);
  else
    fillin(p, p->next->back, p->next->next->back);
}  /* postorder */


void getnufork(node **nufork,node **grbg,pointarray treenode,long *zeros)
{
  /* find a fork not used currently */
  long i;

  i = spp;
  while (treenode[i] && treenode[i]->numdesc > 0) i++;
  if (!treenode[i])
    gnutreenode(grbg, &treenode[i], i, endsite, zeros);
  *nufork = treenode[i];
} /* getnufork */


void reroot(node *outgroup, node *root)
{
  /* reorients tree, putting outgroup in desired position. used if
     the root is binary. */
  /* used in dnacomp & dnapars */
  node *p, *q;

  if (outgroup->back->index == root->index)
    return;
  p = root->next;
  q = root->next->next;
  p->back->back = q->back;
  q->back->back = p->back;
  p->back = outgroup;
  q->back = outgroup->back;
  outgroup->back->back = q;
  outgroup->back = p;
}  /* reroot */


void reroot2(node *outgroup, node *root)
{
  /* reorients tree, putting outgroup in desired position. */
  /* used in dnacomp & dnapars */
  node *p;

  p = outgroup->back->next;
  while (p->next != outgroup->back)
    p = p->next;
  root->next = outgroup->back;
  p->next = root;
}  /* reroot2 */


void reroot3(node *outgroup, node *root, node *root2, node *lastdesc,
                        node **grbg)
{
  /* reorients tree, putting back outgroup in original position. */
  /* used in dnacomp & dnapars */
  node *p;

  p = root->next;
  while (p->next != root)
    p = p->next;
  chuck(grbg, root);
  p->next = outgroup->back;
  root2->next = lastdesc->next;
  lastdesc->next = root2;
}  /* reroot3 */


void savetraverse(node *p)
{
  /* sets BOOLEANs that indicate which way is down */
  node *q;

  p->bottom = true;
  if (p->tip)
    return;
  q = p->next;
  while (q != p) {
    q->bottom = false;
    savetraverse(q->back);
    q = q->next;
  }
}  /* savetraverse */


void newindex(long i, node *p)
{
  /* assigns index i to node p */

  while (p->index != i) {
    p->index = i;
    p = p->next;
  }
} /* newindex */


void flipindexes(long nextnode, pointarray treenode)
{
  /* flips index of nodes between nextnode and last node.  */
  long last;
  node *temp;

  last = nonodes;
  while (treenode[last - 1]->numdesc == 0)
    last--;
  if (last > nextnode) {
    temp = treenode[nextnode - 1];
    treenode[nextnode - 1] = treenode[last - 1];
    treenode[last - 1] = temp;
    newindex(nextnode, treenode[nextnode - 1]);
    newindex(last, treenode[last - 1]);
  }
} /* flipindexes */  


boolean parentinmulti(node *anode)
{
  /* sees if anode's parent has more than 2 children */
  node *p;

  while (!anode->bottom) anode = anode->next;
  p = anode->back;
  while (!p->bottom)
    p = p->next;
  return (p->numdesc > 2);
} /* parentinmulti */


long sibsvisited(node *anode, long *place)
{
  /* computes the number of nodes which are visited earlier than anode among 
  its siblings */
  node *p;
  long nvisited;

  while (!anode->bottom) anode = anode->next;
  p = anode->back->next;
  nvisited = 0;
  do {
    if (!p->bottom && place[p->back->index - 1] != 0)
      nvisited++;
    p = p->next;
  } while (p != anode->back);
  return nvisited;
}  /* sibsvisited */


long smallest(node *anode, long *place)
{
  /* finds the smallest index of sibling of anode */
  node *p;
  long min;

  while (!anode->bottom) anode = anode->next;
  p = anode->back->next;
  if (p->bottom) p = p->next;
  min = nonodes;
  do {
    if (p->back && place[p->back->index - 1] != 0) {
      if (p->back->index <= spp) {
        if (p->back->index < min)
          min = p->back->index;
      } else {
        if (place[p->back->index - 1] < min)
          min = place[p->back->index - 1];
      }
    }
    p = p->next;
    if (p->bottom) p = p->next;
  } while (p != anode->back);
  return min;
}  /* smallest */


void bintomulti(node **root, node **binroot, node **grbg, long *zeros)
{  /* attaches root's left child to its right child and makes
      the right child new root */
  node *left, *right, *newnode, *temp;

  right = (*root)->next->next->back;
  left = (*root)->next->back;
  if (right->tip) {
    (*root)->next = right->back;
    (*root)->next->next = left->back;
    temp = left;
    left = right;
    right = temp;
    right->back->next = *root;
  }
  gnutreenode(grbg, &newnode, right->index, endsite, zeros);
  newnode->next = right->next;
  newnode->back = left;
  left->back = newnode;
  right->next = newnode;
  (*root)->next->back = (*root)->next->next->back = NULL;
  *binroot = *root;
  (*binroot)->numdesc = 0;
  *root = right;
  (*root)->numdesc++;
  (*root)->back = NULL;
} /* bintomulti */


void backtobinary(node **root, node *binroot, node **grbg)
{ /* restores binary root */
  node *p;

  binroot->next->back = (*root)->next->back;
  (*root)->next->back->back = binroot->next;
  p = (*root)->next;
  (*root)->next = p->next;
  binroot->next->next->back = *root;
  (*root)->back = binroot->next->next;
  chuck(grbg, p);
  (*root)->numdesc--;
  *root = binroot;
  (*root)->numdesc = 2;
} /* backtobinary */


boolean outgrin(node *root, node *outgrnode)
{ /* checks if outgroup node is a child of root */
  node *p;

  p = root->next;
  while (p != root) {
    if (p->back == outgrnode)
      return true;
    p = p->next;
  }
  return false;
} /* outgrin */


void flipnodes(node *nodea, node *nodeb)
{ /* flip nodes */
  node *backa, *backb;

  backa = nodea->back;
  backb = nodeb->back;
  backa->back = nodeb;
  backb->back = nodea;
  nodea->back = backb;
  nodeb->back = backa;
} /* flipnodes */


void moveleft(node *root, node *outgrnode, node **flipback)
{ /* makes outgroup node to leftmost child of root */
  node *p;
  boolean done;

  p = root->next;
  done = false;
  while (p != root && !done) {
    if (p->back == outgrnode) {
      *flipback = p;
      flipnodes(root->next->back, p->back);
      done = true;
    }
    p = p->next;
  }
} /* moveleft */


void savetree(node *root,  long *place, pointarray treenode,
                        node **grbg, long *zeros)
{ /* record in place where each species has to be
     added to reconstruct this tree */
  /* used by dnacomp & dnapars */
  long i, j, nextnode, nvisited;
  node *p, *q, *r = NULL, *root2, *lastdesc, 
                *outgrnode, *binroot, *flipback;
  boolean done, newfork;

  binroot = NULL;
  lastdesc = NULL;
  root2 = NULL;
  flipback = NULL;
  outgrnode = treenode[outgrno - 1];
  if (root->numdesc == 2)
    bintomulti(&root, &binroot, grbg, zeros);
  if (outgrin(root, outgrnode)) {
    if (outgrnode != root->next->back)
      moveleft(root, outgrnode, &flipback);
  } else {
    root2 = root;
    lastdesc = root->next;
    while (lastdesc->next != root)
      lastdesc = lastdesc->next;
    lastdesc->next = root->next;
    gnutreenode(grbg, &root, outgrnode->back->index, endsite, zeros);
    root->numdesc = root2->numdesc;
    reroot2(outgrnode, root);
  }
  savetraverse(root);
  nextnode = spp + 1;
  for (i = nextnode; i <= nonodes; i++)
    if (treenode[i - 1]->numdesc == 0)
      flipindexes(i, treenode);
  for (i = 0; i < nonodes; i++)
    place[i] = 0;
  place[root->index - 1] = 1;
  for (i = 1; i <= spp; i++) {
    p = treenode[i - 1];
    while (place[p->index - 1] == 0) {
      place[p->index - 1] = i;
      while (!p->bottom)
        p = p->next;
      r = p;
      p = p->back;
    }
    if (i > 1) {
      q = treenode[i - 1]; 
      newfork = true;
      nvisited = sibsvisited(q, place);
      if (nvisited == 0) {
        if (parentinmulti(r)) {
          nvisited = sibsvisited(r, place);
          if (nvisited == 0)
            place[i - 1] = place[p->index - 1];
          else if (nvisited == 1)
            place[i - 1] = smallest(r, place);
          else {
            place[i - 1] = -smallest(r, place);
            newfork = false;
          }
        } else
          place[i - 1] = place[p->index - 1];
      } else if (nvisited == 1) {
        place[i - 1] = place[p->index - 1];
      } else {
        place[i - 1] = -smallest(q, place);
        newfork = false;
      }
      if (newfork) {
        j = place[p->index - 1];
        done = false;
        while (!done) {
          place[p->index - 1] = nextnode;
          while (!p->bottom)
            p = p->next;
          p = p->back;
          done = (p == NULL);
          if (!done)
            done = (place[p->index - 1] != j);
          if (done) {
            nextnode++;
          }
        }
      }
    }
  }
  if (flipback)
    flipnodes(outgrnode, flipback->back);
  else {
    if (root2) {
      reroot3(outgrnode, root, root2, lastdesc, grbg);
      root = root2;
    }
  }
  if (binroot)
    backtobinary(&root, binroot, grbg);
}  /* savetree */ 


void addnsave(node *p, node *item, node *nufork, node **root, node **grbg,
                boolean multf, pointarray treenode, long *place, long *zeros)
{  /* adds item to tree and save it.  Then removes item. */
  node *dummy;

  if (!multf)
    add(p, item, nufork, root, false, treenode, grbg, zeros);
  else
    add(p, item, NULL, root, false, treenode, grbg, zeros);
  savetree(*root, place, treenode, grbg, zeros);
  if (!multf)
    re_move(item, &nufork, root, false, treenode, grbg, zeros);
  else
    re_move(item, &dummy, root, false, treenode, grbg, zeros);
} /* addnsave */


void addbestever(long *pos, long *nextree, long maxtrees, boolean collapse,
                        long *place, bestelm *bestrees)
{ /* adds first best tree */

  *pos = 1;
  *nextree = 1;
  initbestrees(bestrees, maxtrees, true);
  initbestrees(bestrees, maxtrees, false);
  addtree(*pos, nextree, collapse, place, bestrees);
} /* addbestever */


void addtiedtree(long pos, long *nextree, long maxtrees, boolean collapse,
                        long *place, bestelm *bestrees)
{ /* add tied tree */

  if (*nextree <= maxtrees)
    addtree(pos, nextree, collapse, place, bestrees);
} /* addtiedtree */


void clearcollapse(pointarray treenode)
{
  /* clears collapse status at a node */
  long i;
  node *p;

  for (i = 0; i < nonodes; i++) {
    treenode[i]->collapse = undefined;
    if (!treenode[i]->tip) {
      p = treenode[i]->next;
      while (p != treenode[i]) {
        p->collapse = undefined;
        p = p->next;
      }
    }
  }
} /* clearcollapse */


void clearbottom(pointarray treenode)
{
  /* clears boolean bottom at a node */
  long i;
  node *p;

  for (i = 0; i < nonodes; i++) {
    treenode[i]->bottom = false;
    if (!treenode[i]->tip) {
      p = treenode[i]->next;
      while (p != treenode[i]) {
        p->bottom = false;
        p = p->next;
      }
    }
  }
} /* clearbottom */


void collabranch(node *collapfrom, node *tempfrom, node *tempto)
{ /* collapse branch from collapfrom */
  long i, j, b, largest, descsteps;
  boolean done;

  for (i = 0; i < endsite; i++) {
    descsteps = 0;
    for (j = (long)A; j <= (long)O; j++) {
      b = 1 << j;
      if ((descsteps == 0) && (collapfrom->base[i] & b)) 
        descsteps = tempfrom->oldnumsteps[i] 
                     - (collapfrom->numdesc - collapfrom->numnuc[i][j])
                       * weight[i];
    }
    done = false;
    for (j = (long)A; j <= (long)O; j++) {
      b = 1 << j;
      if (!done && (tempto->base[i] & b)) {
        descsteps += (tempto->numsteps[i] 
                      - (tempto->numdesc - collapfrom->numdesc
                        - tempto->numnuc[i][j]) * weight[i]);
        done = true;
      }
    }
    for (j = (long)A; j <= (long)O; j++)
      tempto->numnuc[i][j] += collapfrom->numnuc[i][j];
    largest = getlargest(tempto->numnuc[i]);
    tempto->base[i] = 0;
    for (j = (long)A; j <= (long)O; j++) {
      if (tempto->numnuc[i][j] == largest)
        tempto->base[i] |= (1 << j);
    }
    tempto->numsteps[i] = (tempto->numdesc - largest) * weight[i] + descsteps;
  }
} /* collabranch */


boolean allcommonbases(node *a, node *b, boolean *allsame)
{  /* see if bases are common at all sites for nodes a and b */    
  long i;
  boolean allcommon;

  allcommon = true;
  *allsame = true;
  for (i = 0; i < endsite; i++) {
    if ((a->base[i] & b->base[i]) == 0)
      allcommon = false;
    else if (a->base[i] != b->base[i])
      *allsame = false;
  }
  return allcommon;
} /* allcommonbases */


void findbottom(node *p, node **bottom)
{ /* find a node with field bottom set at node p */
  node *q;

  if (p->bottom)
    *bottom = p;
  else {
    q = p->next;
    while(!q->bottom && q != p)
      q = q->next;
    *bottom = q;
  }
} /* findbottom */


boolean moresteps(node *a, node *b)
{  /* see if numsteps of node a exceeds those of node b */    
  long i;

  for (i = 0; i < endsite; i++)
    if (a->numsteps[i] > b->numsteps[i])
      return true;
  return false;
} /* moresteps */


boolean passdown(node *desc, node *parent, node *start, node *below,
                        node *item, node *added, node *total, node *tempdsc,
            node *tempprt, boolean multf)
{ /* track down to node start to see if an ancestor branch can be collapsed */
  node *temp;
  boolean done, allsame;

  done = (parent == start);
  while (!done) {
    desc = parent;
    findbottom(parent->back, &parent);
    if (multf && start == below && parent == below)
      parent = added;
    memcpy(tempdsc->base, tempprt->base, endsite*sizeof(long));
    memcpy(tempdsc->numsteps, tempprt->numsteps, endsite*sizeof(long));
    memcpy(tempdsc->oldbase, desc->base, endsite*sizeof(long));
    memcpy(tempdsc->oldnumsteps, desc->numsteps, endsite*sizeof(long));
    memcpy(tempprt->base, parent->base, endsite*sizeof(long));
    memcpy(tempprt->numsteps, parent->numsteps, endsite*sizeof(long));
    memcpy(tempprt->numnuc, parent->numnuc, endsite*sizeof(nucarray));
    tempprt->numdesc = parent->numdesc;
    multifillin(tempprt, tempdsc, 0);
    if (!allcommonbases(tempprt, parent, &allsame))
      return false;
    else if (moresteps(tempprt, parent))
      return false;
    else if (allsame)
      return true;
    if (parent == added)
      parent = below;
    done = (parent == start);
    if (done && ((start == item) || (!multf && start == below))) {
      memcpy(tempdsc->base, tempprt->base, endsite*sizeof(long));
      memcpy(tempdsc->numsteps, tempprt->numsteps, endsite*sizeof(long));
      memcpy(tempdsc->oldbase, start->base, endsite*sizeof(long));
      memcpy(tempdsc->oldnumsteps, start->numsteps, endsite*sizeof(long));
      multifillin(added, tempdsc, 0);
      tempprt = added;
    }
  } 
  temp = tempdsc;
  if (start == below || start == item)
    fillin(temp, tempprt, below->back);
  else
    fillin(temp, tempprt, added);
  return !moresteps(temp, total);
} /* passdown */


boolean trycollapdesc(node *desc, node *parent, node *start,
                        node *below, node *item, node *added, node *total,
            node *tempdsc, node *tempprt, boolean multf, long *zeros)
  { /* see if branch between nodes desc and parent can be collapsed */
  boolean allsame;

  if (desc->numdesc == 1)
    return true;
  if (multf && start == below && parent == below)
    parent = added;
  memcpy(tempdsc->base, zeros, endsite*sizeof(long));
  memcpy(tempdsc->numsteps, zeros, endsite*sizeof(long));
  memcpy(tempdsc->oldbase, desc->base, endsite*sizeof(long));
  memcpy(tempdsc->oldnumsteps, desc->numsteps, endsite*sizeof(long));
  memcpy(tempprt->base, parent->base, endsite*sizeof(long));
  memcpy(tempprt->numsteps, parent->numsteps, endsite*sizeof(long));
  memcpy(tempprt->numnuc, parent->numnuc, endsite*sizeof(nucarray));
  tempprt->numdesc = parent->numdesc - 1;
  multifillin(tempprt, tempdsc, -1);
  tempprt->numdesc += desc->numdesc;
  collabranch(desc, tempdsc, tempprt);
  if (!allcommonbases(tempprt, parent, &allsame) || 
        moresteps(tempprt, parent)) {
    if (parent != added) {
      desc->collapse = nocollap;
      parent->collapse = nocollap;
    }
    return false;
  } else if (allsame) {
    if (parent != added) {
      desc->collapse = tocollap;
      parent->collapse = tocollap;
    }
    return true;
  }
  if (parent == added)
    parent = below;
  if ((start == item && parent == item) ||
        (!multf && start == below && parent == below)) {
    memcpy(tempdsc->base, tempprt->base, endsite*sizeof(long));
    memcpy(tempdsc->numsteps, tempprt->numsteps, endsite*sizeof(long));
    memcpy(tempdsc->oldbase, start->base, endsite*sizeof(long));
    memcpy(tempdsc->oldnumsteps, start->numsteps, endsite*sizeof(long));
    memcpy(tempprt->base, added->base, endsite*sizeof(long));
    memcpy(tempprt->numsteps, added->numsteps, endsite*sizeof(long));
    memcpy(tempprt->numnuc, added->numnuc, endsite*sizeof(nucarray));
    tempprt->numdesc = added->numdesc;
    multifillin(tempprt, tempdsc, 0);
    if (!allcommonbases(tempprt, added, &allsame))
      return false;
    else if (moresteps(tempprt, added))
      return false;
    else if (allsame)
      return true;
  }
  return passdown(desc, parent, start, below, item, added, total, tempdsc,
                    tempprt, multf);
} /* trycollapdesc */


void setbottom(node *p)
{ /* set field bottom at node p */
  node *q;

  p->bottom = true;
  q = p->next;
  do {
    q->bottom = false;
    q = q->next;
  } while (q != p);
} /* setbottom */

boolean zeroinsubtree(node *subtree, node *start, node *below, node *item,
                        node *added, node *total, node *tempdsc, node *tempprt,
                        boolean multf, node* root, long *zeros)
{ /* sees if subtree contains a zero length branch */
  node *p;

  if (!subtree->tip) {
    setbottom(subtree);
    p = subtree->next;
    do {
      if (p->back && !p->back->tip && 
         !((p->back->collapse == nocollap) && (subtree->collapse == nocollap))
           && (subtree->numdesc != 1)) {
        if ((p->back->collapse == tocollap) && (subtree->collapse == tocollap)
            && multf && (subtree != below))
          return true;
        /* when root->numdesc == 2
         * there is no mandatory step at the root, 
         * instead of checking at the root we check around it 
         * we only need to check p because the first if 
         * statement already gets rid of it for the subtree */
        else if ((p->back->index != root->index || root->numdesc > 2) && 
            trycollapdesc(p->back, subtree, start, below, item, added, total, 
                tempdsc, tempprt, multf, zeros))
          return true;
        else if ((p->back->index == root->index && root->numdesc == 2) && 
            !(root->next->back->tip) && !(root->next->next->back->tip) && 
            trycollapdesc(root->next->back, root->next->next->back, start, 
                below, item,added, total, tempdsc, tempprt, multf, zeros))
          return true;
      }
      p = p->next;
    } while (p != subtree);
    p = subtree->next;
    do {
      if (p->back && !p->back->tip) {
        if (zeroinsubtree(p->back, start, below, item, added, total, 
                            tempdsc, tempprt, multf, root, zeros))
          return true;
      }
      p = p->next;
    } while (p != subtree);
  }
  return false;
} /* zeroinsubtree */


boolean collapsible(node *item, node *below, node *temp, node *temp1,
                        node *tempdsc, node *tempprt, node *added, node *total,
            boolean multf, node *root, long *zeros, pointarray treenode)
{
  /* sees if any branch can be collapsed */
  node *belowbk;
  boolean allsame;

  if (multf) {
    memcpy(tempdsc->base, item->base, endsite*sizeof(long));
    memcpy(tempdsc->numsteps, item->numsteps, endsite*sizeof(long));
    memcpy(tempdsc->oldbase, zeros, endsite*sizeof(long));
    memcpy(tempdsc->oldnumsteps, zeros, endsite*sizeof(long));
    memcpy(added->base, below->base, endsite*sizeof(long));
    memcpy(added->numsteps, below->numsteps, endsite*sizeof(long));
    memcpy(added->numnuc, below->numnuc, endsite*sizeof(nucarray));
    added->numdesc = below->numdesc + 1;
    multifillin(added, tempdsc, 1);
  } else {
    fillin(added, item, below);
    added->numdesc = 2;
  }
  fillin(total, added, below->back);
  clearbottom(treenode);
  if (below->back) {
    if (zeroinsubtree(below->back, below->back, below, item, added, total,
                        tempdsc, tempprt, multf, root, zeros))
      return true;
  }
  if (multf) {
    if (zeroinsubtree(below, below, below, item, added, total,
                       tempdsc, tempprt, multf, root, zeros))
      return true;
  } else if (!below->tip) {
    if (zeroinsubtree(below, below, below, item, added, total,
                        tempdsc, tempprt, multf, root, zeros))
      return true;
  }
  if (!item->tip) {
    if (zeroinsubtree(item, item, below, item, added, total,
                        tempdsc, tempprt, multf, root, zeros))
      return true;
  }
  if (multf && below->back && !below->back->tip) {
    memcpy(tempdsc->base, zeros, endsite*sizeof(long));
    memcpy(tempdsc->numsteps, zeros, endsite*sizeof(long));
    memcpy(tempdsc->oldbase, added->base, endsite*sizeof(long));
    memcpy(tempdsc->oldnumsteps, added->numsteps, endsite*sizeof(long));
    if (below->back == treenode[below->back->index - 1])
      belowbk = below->back->next;
    else
      belowbk = treenode[below->back->index - 1];
    memcpy(tempprt->base, belowbk->base, endsite*sizeof(long));
    memcpy(tempprt->numsteps, belowbk->numsteps, endsite*sizeof(long));
    memcpy(tempprt->numnuc, belowbk->numnuc, endsite*sizeof(nucarray));
    tempprt->numdesc = belowbk->numdesc - 1;
    multifillin(tempprt, tempdsc, -1);
    tempprt->numdesc += added->numdesc;
    collabranch(added, tempdsc, tempprt);
    if (!allcommonbases(tempprt, belowbk, &allsame))
      return false;
    else if (allsame && !moresteps(tempprt, belowbk))
      return true;
    else if (belowbk->back) {
      fillin(temp, tempprt, belowbk->back);
      fillin(temp1, belowbk, belowbk->back);
      return !moresteps(temp, temp1);
    }
  }
  return false;
} /* collapsible */


void replaceback(node **oldback, node *item, node *forknode,
                        node **grbg, long *zeros)
{ /* replaces back node of item with another */
  node *p;

  p = forknode;
  while (p->next->back != item)
    p = p->next;
  *oldback = p->next;
  gnutreenode(grbg, &p->next, forknode->index, endsite, zeros);
  p->next->next = (*oldback)->next;
  p->next->back = (*oldback)->back;
  p->next->back->back = p->next;
  (*oldback)->next = (*oldback)->back = NULL;
} /* replaceback */


void putback(node *oldback, node *item, node *forknode, node **grbg)
{ /* restores node to back of item */
  node *p, *q;

  p = forknode;
  while (p->next != item->back)
    p = p->next;
  q = p->next;
  oldback->next = p->next->next;
  p->next = oldback;
  oldback->back = item;
  item->back = oldback;
  oldback->index = forknode->index;
  chuck(grbg, q);
} /* putback */


void savelocrearr(node *item, node *forknode, node *below, node *tmp,
        node *tmp1, node *tmp2, node *tmp3, node *tmprm, node *tmpadd,
        node **root, long maxtrees, long *nextree, boolean multf,
        boolean bestever, boolean *saved, long *place,
        bestelm *bestrees, pointarray treenode, node **grbg,
        long *zeros)
{ /* saves tied or better trees during local rearrangements by removing
     item from forknode and adding to below */
  node *other, *otherback = NULL, *oldfork, *nufork, *oldback;
  long pos;
  boolean found, collapse;

  if (forknode->numdesc == 2) {
    findbelow(&other, item, forknode);
    otherback = other->back;
    oldback = NULL;
  } else {
    other = NULL;
    replaceback(&oldback, item, forknode, grbg, zeros);
  }
  re_move(item, &oldfork, root, false, treenode, grbg, zeros);
  if (!multf)
    getnufork(&nufork, grbg, treenode, zeros);
  else
    nufork = NULL;
  addnsave(below, item, nufork, root, grbg, multf, treenode, place, zeros);
  pos = 0;
  findtree(&found, &pos, *nextree, place, bestrees);
  if (other) {
    add(other, item, oldfork, root, false, treenode, grbg, zeros);
    if (otherback->back != other)
      flipnodes(item, other);
  } else
    add(forknode, item, NULL, root, false, treenode, grbg, zeros);
  *saved = false;
  if (found) {
    if (oldback)
      putback(oldback, item, forknode, grbg);
  } else {
    if (oldback)
      chuck(grbg, oldback);
    re_move(item, &oldfork, root, true, treenode, grbg, zeros);
    collapse = collapsible(item, below, tmp, tmp1, tmp2, tmp3, tmprm,
                     tmpadd, multf, *root, zeros, treenode);
    if (!collapse) {
      if (bestever)
        addbestever(&pos, nextree, maxtrees, collapse, place, bestrees);
      else
        addtiedtree(pos, nextree, maxtrees, collapse, place, bestrees);
    }
    if (other)
      add(other, item, oldfork, root, true, treenode, grbg, zeros);
    else
      add(forknode, item, NULL, root, true, treenode, grbg, zeros);
    *saved = !collapse;
  }
} /* savelocrearr */


void clearvisited(pointarray treenode)
{
  /* clears boolean visited at a node */
  long i;
  node *p;

  for (i = 0; i < nonodes; i++) {
    treenode[i]->visited = false;
    if (!treenode[i]->tip) {
      p = treenode[i]->next;
      while (p != treenode[i]) {
        p->visited = false;
        p = p->next;
      }
    }
  }
} /* clearvisited */


void hyprint(long b1, long b2, struct LOC_hyptrav *htrav,
                        pointarray treenode, Char *basechar)
{
  /* print out states in sites b1 through b2 at node */
  long i, j, k, n;
  boolean dot;
  bases b;

  if (htrav->bottom) {
    if (!outgropt)
      fprintf(outfile, "       ");
    else
      fprintf(outfile, "root   ");
  } else
    fprintf(outfile, "%4ld   ", htrav->r->back->index - spp);
  if (htrav->r->tip) {
    for (i = 0; i < nmlngth; i++)
      putc(nayme[htrav->r->index - 1][i], outfile);
  } else
    fprintf(outfile, "%4ld      ", htrav->r->index - spp);
  if (htrav->bottom)
    fprintf(outfile, "          ");
  else if (htrav->nonzero)
    fprintf(outfile, "   yes    ");
  else if (htrav->maybe)
    fprintf(outfile, "  maybe   ");
  else
    fprintf(outfile, "   no     ");
  for (i = b1; i <= b2; i++) {
    j = location[ally[i - 1] - 1];
    htrav->tempset = htrav->r->base[j - 1];
    htrav->anc = htrav->hypset[j - 1];
    if (!htrav->bottom)
      htrav->anc = treenode[htrav->r->back->index - 1]->base[j - 1];
    dot = dotdiff && (htrav->tempset == htrav->anc && !htrav->bottom);
    if (dot)
      putc('.', outfile); 
    else if (htrav->tempset == (1 << A))
      putc('A', outfile);
    else if (htrav->tempset == (1 << C))
      putc('C', outfile);
    else if (htrav->tempset == (1 << G))
      putc('G', outfile);
    else if (htrav->tempset == (1 << T))
      putc('T', outfile);
    else if (htrav->tempset == (1 << O))
      putc('-', outfile);
    else {
      k = 1;
      n = 0;
      for (b = A; b <= O; b = b + 1) {
        if (((1 << b) & htrav->tempset) != 0)
          n += k;
        k += k;
      }
      putc(basechar[n - 1], outfile);
    }
    if (i % 10 == 0)
      putc(' ', outfile);
  }
  putc('\n', outfile);
}  /* hyprint */


void gnubase(gbases **p, gbases **garbage, long endsite)
{
  /* this and the following are do-it-yourself garbage collectors.
     Make a new node or pull one off the garbage list */
  if (*garbage != NULL) {
    *p = *garbage;
    *garbage = (*garbage)->next;
  } else {
    *p = (gbases *)Malloc(sizeof(gbases));
    (*p)->base = (baseptr)Malloc(endsite*sizeof(long));
  }
  (*p)->next = NULL;
}  /* gnubase */


void chuckbase(gbases *p, gbases **garbage)
{
  /* collect garbage on p -- put it on front of garbage list */
  p->next = *garbage;
  *garbage = p;
}  /* chuckbase */


void hyptrav(node *r_, long *hypset_, long b1, long b2, boolean bottom_,
                        pointarray treenode, gbases **garbage, Char *basechar)
{
  /*  compute, print out states at one interior node */
  struct LOC_hyptrav Vars;
  long i, j, k;
  long largest;
  gbases *ancset;
  nucarray *tempnuc;
  node *p, *q;

  Vars.bottom = bottom_;
  Vars.r = r_;
  Vars.hypset = hypset_;
  gnubase(&ancset, garbage, endsite);
  tempnuc = (nucarray *)Malloc(endsite*sizeof(nucarray));
  Vars.maybe = false;
  Vars.nonzero = false;
  if (!Vars.r->tip)
    zeronumnuc(Vars.r, endsite);
  for (i = b1 - 1; i < b2; i++) {
    j = location[ally[i] - 1];
    Vars.anc = Vars.hypset[j - 1];
    if (!Vars.r->tip) {
      p = Vars.r->next;
      for (k = (long)A; k <= (long)O; k++)
        if (Vars.anc & (1 << k))
          Vars.r->numnuc[j - 1][k]++;
      do {
        for (k = (long)A; k <= (long)O; k++)
          if (p->back->base[j - 1] & (1 << k))
            Vars.r->numnuc[j - 1][k]++;
        p = p->next;
      } while (p != Vars.r);
      largest = getlargest(Vars.r->numnuc[j - 1]);
      Vars.tempset = 0;
      for (k = (long)A; k <= (long)O; k++) {
        if (Vars.r->numnuc[j - 1][k] == largest)
          Vars.tempset |= (1 << k);
      }
      Vars.r->base[j - 1] = Vars.tempset;
    }
    if (!Vars.bottom)
      Vars.anc = treenode[Vars.r->back->index - 1]->base[j - 1];
    Vars.nonzero = (Vars.nonzero || (Vars.r->base[j - 1] & Vars.anc) == 0);
    Vars.maybe = (Vars.maybe || Vars.r->base[j - 1] != Vars.anc);
  }
  hyprint(b1, b2, &Vars, treenode, basechar);
  Vars.bottom = false;
  if (!Vars.r->tip) {
    memcpy(tempnuc, Vars.r->numnuc, endsite*sizeof(nucarray));
    q = Vars.r->next;
    do {
      memcpy(Vars.r->numnuc, tempnuc, endsite*sizeof(nucarray));
      for (i = b1 - 1; i < b2; i++) {
        j = location[ally[i] - 1];
        for (k = (long)A; k <= (long)O; k++)
          if (q->back->base[j - 1] & (1 << k))
            Vars.r->numnuc[j - 1][k]--;
        largest = getlargest(Vars.r->numnuc[j - 1]);
        ancset->base[j - 1] = 0;
        for (k = (long)A; k <= (long)O; k++)
          if (Vars.r->numnuc[j - 1][k] == largest)
            ancset->base[j - 1] |= (1 << k);
        if (!Vars.bottom)
          Vars.anc = ancset->base[j - 1];
      }
      hyptrav(q->back, ancset->base, b1, b2, Vars.bottom,
                treenode, garbage, basechar);
      q = q->next;
    } while (q != Vars.r);
  }
  chuckbase(ancset, garbage);
}  /* hyptrav */


void hypstates(long chars, node *root, pointarray treenode,
                        gbases **garbage, Char *basechar)
{
  /* fill in and describe states at interior nodes */
  /* used in dnacomp, dnapars, & dnapenny */
  long i, n;
  baseptr nothing;

  fprintf(outfile, "\nFrom    To     Any Steps?    State at upper node\n");
  fprintf(outfile, "                            ");
  if (dotdiff)
    fprintf(outfile, " ( . means same as in the node below it on tree)\n");
  nothing = (baseptr)Malloc(endsite*sizeof(long));
  for (i = 0; i < endsite; i++)
    nothing[i] = 0;
  for (i = 1; i <= ((chars - 1) / 40 + 1); i++) {
    putc('\n', outfile);
    n = i * 40;
    if (n > chars)
      n = chars;
    hyptrav(root, nothing, i * 40 - 39, n, true, treenode, garbage, basechar);
  }
  free(nothing);
}  /* hypstates */


void initbranchlen(node *p)
{
  node *q;

  p->v = 0.0;
  if (p->back)
    p->back->v = 0.0;
  if (p->tip)
    return;
  q = p->next;
  while (q != p) {
    initbranchlen(q->back);
    q = q->next;
  }
  q = p->next;
  while (q != p) {
    q->v = 0.0;
    q = q->next;
  }
} /* initbranchlen */


void initmin(node *p, long sitei, boolean internal)
{
  long i;

  if (internal) {
    for (i = (long)A; i <= (long)O; i++) {
      p->cumlengths[i] = 0;
      p->numreconst[i] = 1;
    }
  } else {
    for (i = (long)A; i <= (long)O; i++) {
      if (p->base[sitei - 1] & (1 << i)) {
        p->cumlengths[i] = 0;
        p->numreconst[i] = 1;
      } else {
        p->cumlengths[i] = -1;
        p->numreconst[i] = 0;
      }
    }
  }
} /* initmin */


void initbase(node *p, long sitei)
{
  /* traverse tree to initialize base at internal nodes */
  node *q;
  long i, largest;

  if (p->tip)
    return;
  q = p->next;
  while (q != p) {
    if (q->back) {
      memcpy(q->numnuc, p->numnuc, endsite*sizeof(nucarray));
      for (i = (long)A; i <= (long)O; i++) {
        if (q->back->base[sitei - 1] & (1 << i))
          q->numnuc[sitei - 1][i]--;
      }
      if (p->back) {
        for (i = (long)A; i <= (long)O; i++) {
          if (p->back->base[sitei - 1] & (1 << i))
            q->numnuc[sitei - 1][i]++;
        }
      }
      largest = getlargest(q->numnuc[sitei - 1]);
      q->base[sitei - 1] = 0;
      for (i = (long)A; i <= (long)O; i++) {
        if (q->numnuc[sitei - 1][i] == largest)
          q->base[sitei - 1] |= (1 << i);
      }
    }
    q = q->next;
  }
  q = p->next;
  while (q != p) {
    initbase(q->back, sitei);
    q = q->next;
  }
} /* initbase */


void inittreetrav(node *p, long sitei)
{
  /* traverse tree to clear boolean initialized and set up base */
  node *q;

  if (p->tip) {
    initmin(p, sitei, false);
    p->initialized = true;
    return;
  }
  q = p->next;
  while (q != p) {
    inittreetrav(q->back, sitei);
    q = q->next;
  }
  initmin(p, sitei, true);
  p->initialized = false;
  q = p->next;
  while (q != p) {
    initmin(q, sitei, true);
    q->initialized = false;
    q = q->next;
  }
} /* inittreetrav */


void compmin(node *p, node *desc)
{
  /* computes minimum lengths up to p */
  long i, j, minn, cost, desclen, descrecon=0, maxx;

  maxx = 10 * spp;
  for (i = (long)A; i <= (long)O; i++) {
    minn = maxx;
    for (j = (long)A; j <= (long)O; j++) {
      if (transvp) {
        if (
               (
                ((i == (long)A) || (i == (long)G))
             && ((j == (long)A) || (j == (long)G))
               )
            || (
                ((j == (long)C) || (j == (long)T))
             && ((i == (long)C) || (i == (long)T))
               )
           )
          cost = 0;
        else
          cost = 1;
      } else {
        if (i == j)
          cost = 0;
        else
          cost = 1;
      }
      if (desc->cumlengths[j] == -1) {
        desclen = maxx;
      } else {
        desclen = desc->cumlengths[j];
      }
      if (minn > cost + desclen) {
        minn = cost + desclen;
        descrecon = 0;
      }
      if (minn == cost + desclen) {
        descrecon += desc->numreconst[j];
      }
    }
    p->cumlengths[i] += minn;
    p->numreconst[i] *= descrecon;
  }
  p->initialized = true;
} /* compmin */


void minpostorder(node *p, pointarray treenode)
{
  /* traverses an n-ary tree, computing minimum steps at each node */
  node *q;

  if (p->tip) {
    return;
  }
  q = p->next;
  while (q != p) {
    if (q->back)
      minpostorder(q->back, treenode);
    q = q->next;
  }
  if (!p->initialized) {
    q = p->next;
    while (q != p) {
      if (q->back)
        compmin(p, q->back);
      q = q->next;
    }
  }
}  /* minpostorder */


void branchlength(node *subtr1, node *subtr2, double *brlen,
                        pointarray treenode)
{
  /* computes a branch length between two subtrees for a given site */
  long i, j, minn, cost, nom, denom;
  node *temp;

  if (subtr1->tip) {
    temp = subtr1;
    subtr1 = subtr2;
    subtr2 = temp;
  }
  if (subtr1->index == outgrno) {
    temp = subtr1;
    subtr1 = subtr2;
    subtr2 = temp;
  }
  minpostorder(subtr1, treenode);
  minpostorder(subtr2, treenode);
  minn = 10 * spp;
  nom = 0;
  denom = 0;
  for (i = (long)A; i <= (long)O; i++) {
    for (j = (long)A; j <= (long)O; j++) {
      if (transvp) {
        if (
               (
                ((i == (long)A) || (i == (long)G))
             && ((j == (long)A) || (j == (long)G))
               )
            || (
                ((j == (long)C) || (j == (long)T))
             && ((i == (long)C) || (i == (long)T))
               )
           )
          cost = 0;
        else
          cost = 1;
      } else {
        if (i == j)
          cost = 0;
        else
          cost = 1;
      }
      if (subtr1->cumlengths[i] != -1 && (subtr2->cumlengths[j] != -1)) {
        if (subtr1->cumlengths[i] + cost + subtr2->cumlengths[j] < minn) {
          minn = subtr1->cumlengths[i] + cost + subtr2->cumlengths[j];
          nom = 0;
          denom = 0;
        }
        if (subtr1->cumlengths[i] + cost + subtr2->cumlengths[j] == minn) {
          nom += subtr1->numreconst[i] * subtr2->numreconst[j] * cost;
          denom += subtr1->numreconst[i] * subtr2->numreconst[j];
        }
      }
    }
  }
  *brlen = (double)nom/(double)denom;
} /* branchlength */  


void printbranchlengths(node *p)
{
  node *q;
  long i;

  if (p->tip)
    return;
  q = p->next;
  do {
    fprintf(outfile, "%6ld      ",q->index - spp);
    if (q->back->tip) {
      for (i = 0; i < nmlngth; i++)
        putc(nayme[q->back->index - 1][i], outfile);
    } else
      fprintf(outfile, "%6ld    ", q->back->index - spp);
    fprintf(outfile, "   %f\n",q->v);
    if (q->back)
      printbranchlengths(q->back);
    q = q->next;
  } while (q != p);
} /* printbranchlengths */


void branchlentrav(node *p, node *root, long sitei, long chars,
                        double *brlen, pointarray treenode)
  {
  /*  traverses the tree computing tree length at each branch */
  node *q;

  if (p->tip)
    return;
  if (p->index == outgrno)
    p = p->back;
  q = p->next;
  do {
    if (q->back) {
      branchlength(q, q->back, brlen, treenode);
      q->v += ((weight[sitei - 1] / 10.0) * (*brlen)/chars);
      q->back->v += ((weight[sitei - 1] / 10.0) * (*brlen)/chars);
      if (!q->back->tip)
        branchlentrav(q->back, root, sitei, chars, brlen, treenode);
    }
    q = q->next;
  } while (q != p);
}  /* branchlentrav */


void treelength(node *root, long chars, pointarray treenode)
  {
  /*  calls branchlentrav at each site */
  long sitei;
  double trlen;

  initbranchlen(root);
  for (sitei = 1; sitei <= endsite; sitei++) {
    trlen = 0.0;
    initbase(root, sitei);
    inittreetrav(root, sitei);
    branchlentrav(root, root, sitei, chars, &trlen, treenode);
  }
} /* treelength */


void coordinates(node *p, long *tipy, double f, long *fartemp)
{
  /* establishes coordinates of nodes for display without lengths */
  node *q, *first, *last;
  node *mid1 = NULL, *mid2 = NULL;
  long numbranches, numb2;

  if (p->tip) {
    p->xcoord = 0;
    p->ycoord = *tipy;
    p->ymin = *tipy;
    p->ymax = *tipy;
    (*tipy) += down;
    return;
  }
  numbranches = 0;
  q = p->next;
  do {
    coordinates(q->back, tipy, f, fartemp);
    numbranches += 1;
    q = q->next;
  } while (p != q);
  first = p->next->back;
  q = p->next;
  while (q->next != p) 
    q = q->next;
  last = q->back;
  numb2 = 1;
  q = p->next;
  while (q != p) {
    if (numb2 == (long)(numbranches + 1)/2)
      mid1 = q->back;
    if (numb2 == (long)(numbranches/2 + 1))
      mid2 = q->back;
    numb2 += 1;
    q = q->next;
  }
  p->xcoord = (long)((double)(last->ymax - first->ymin) * f);
  p->ycoord = (long)((mid1->ycoord + mid2->ycoord) / 2);
  p->ymin = first->ymin;
  p->ymax = last->ymax;
  if (p->xcoord > *fartemp)
    *fartemp = p->xcoord;
}  /* coordinates */


void drawline(long i, double scale, node *root)
{
  /* draws one row of the tree diagram by moving up tree */
  node *p, *q, *r, *first =NULL, *last =NULL;
  long n, j;
  boolean extra, done, noplus;

  p = root;
  q = root;
  extra = false;
  noplus = false;
  if (i == (long)p->ycoord && p == root) {
    if (p->index - spp >= 10)
      fprintf(outfile, " %2ld", p->index - spp);
    else
      fprintf(outfile, "  %ld", p->index - spp);
    extra = true;
    noplus = true;
  } else
    fprintf(outfile, "  ");
  do {
    if (!p->tip) {
      r = p->next;
      done = false;
      do {
        if (i >= r->back->ymin && i <= r->back->ymax) {
          q = r->back;
          done = true;
        }
        r = r->next;
      } while (!(done || r == p));
      first = p->next->back;
      r = p->next;
      while (r->next != p)
        r = r->next;
      last = r->back;
    }
    done = (p == q);
    n = (long)(scale * (p->xcoord - q->xcoord) + 0.5);
    if (n < 3 && !q->tip)
      n = 3;
    if (extra) {
      n--;
      extra = false;
    }
    if ((long)q->ycoord == i && !done) {
      if (noplus) {
        putc('-', outfile);
        noplus = false;
      }
      else
        putc('+', outfile);
      if (!q->tip) {
        for (j = 1; j <= n - 2; j++)
          putc('-', outfile);
        if (q->index - spp >= 10)
          fprintf(outfile, "%2ld", q->index - spp);
        else
          fprintf(outfile, "-%ld", q->index - spp);
        extra = true;
        noplus = true;
      } else {
        for (j = 1; j < n; j++)
          putc('-', outfile);
      }
    } else if (!p->tip) {
      if ((long)last->ycoord > i && (long)first->ycoord < i
            && i != (long)p->ycoord) {
        putc('!', outfile);
        for (j = 1; j < n; j++)
          putc(' ', outfile);
      } else {
        for (j = 1; j <= n; j++)
          putc(' ', outfile);
      }
      noplus = false;
    } else {
      for (j = 1; j <= n; j++)
        putc(' ', outfile);
      noplus = false;
    }
    if (p != q)
      p = q;
  } while (!done);
  if ((long)p->ycoord == i && p->tip) {
    for (j = 0; j < nmlngth; j++)
      putc(nayme[p->index - 1][j], outfile);
  }
  putc('\n', outfile);
}  /* drawline */


void printree(node *root, double f)
{
  /* prints out diagram of the tree */
  /* used in dnacomp, dnapars, & dnapenny */
  long i, tipy, dummy;
  double scale;

  putc('\n', outfile);
  if (!treeprint)
    return;
  putc('\n', outfile);
  tipy = 1;
  dummy = 0;
  coordinates(root, &tipy, f, &dummy);
  scale = 1.5;
  putc('\n', outfile);
  for (i = 1; i <= (tipy - down); i++)
    drawline(i, scale, root);
  fprintf(outfile, "\n  remember:");
  if (outgropt)
    fprintf(outfile, " (although rooted by outgroup)");
  fprintf(outfile, " this is an unrooted tree!\n\n");
}  /* printree */


void writesteps(long chars, boolean weights, steptr oldweight, node *root)
{
  /* used in dnacomp, dnapars, & dnapenny */
  long i, j, k, l;

  putc('\n', outfile);
  if (weights)
    fprintf(outfile, "weighted ");
  fprintf(outfile, "steps in each site:\n");
  fprintf(outfile, "      ");
  for (i = 0; i <= 9; i++)
    fprintf(outfile, "%4ld", i);
  fprintf(outfile, "\n     *------------------------------------");
  fprintf(outfile, "-----\n");
  for (i = 0; i <= (chars / 10); i++) {
    fprintf(outfile, "%5ld", i * 10);
    putc('|', outfile);
    for (j = 0; j <= 9; j++) {
      k = i * 10 + j;
      if (k == 0 || k > chars)
        fprintf(outfile, "    ");
      else {
        l = location[ally[k - 1] - 1];
        if (oldweight[k - 1] > 0)
          fprintf(outfile, "%4ld",
                  oldweight[k - 1] *
                  (root->numsteps[l - 1] / weight[l - 1]));
        else
          fprintf(outfile, "   0");
      }
    }
    putc('\n', outfile);
  }
} /* writesteps */


void treeout(node *p, long nextree, long *col, node *root)
{
  /* write out file with representation of final tree */
  /* used in dnacomp, dnamove, dnapars, & dnapenny */
  node *q;
  long i, n;
  Char c;

  if (p->tip) {
    n = 0;
    for (i = 1; i <= nmlngth; i++) {
      if (nayme[p->index - 1][i - 1] != ' ')
        n = i;
    }
    for (i = 0; i < n; i++) {
      c = nayme[p->index - 1][i];
      if (c == ' ')
        c = '_';
      putc(c, outtree);
    }
    *col += n;
  } else {
    putc('(', outtree);
    (*col)++;
    q = p->next;
    while (q != p) {
      treeout(q->back, nextree, col, root);
      q = q->next;
      if (q == p)
        break;
      putc(',', outtree);
      (*col)++;
      if (*col > 60) {
        putc('\n', outtree);
        *col = 0;
      }
    }
    putc(')', outtree);
    (*col)++;
  }
  if (p != root)
    return;
  if (nextree > 2)
    fprintf(outtree, "[%6.4f];\n", 1.0 / (nextree - 1));
  else
    fprintf(outtree, ";\n");
}  /* treeout */


void treeout3(node *p, long nextree, long *col, node *root)
{
  /* write out file with representation of final tree */
  /* used in dnapars -- writes branch lengths */
  node *q;
  long i, n, w;
  double x;
  Char c;

  if (p->tip) {
    n = 0;
    for (i = 1; i <= nmlngth; i++) {
      if (nayme[p->index - 1][i - 1] != ' ')
        n = i;
    }
    for (i = 0; i < n; i++) {
      c = nayme[p->index - 1][i];
      if (c == ' ')
        c = '_';
      putc(c, outtree);
    }
    *col += n;
  } else {
    putc('(', outtree);
    (*col)++;
    q = p->next;
    while (q != p) {
      treeout3(q->back, nextree, col, root);
      q = q->next;
      if (q == p)
        break;
      putc(',', outtree);
      (*col)++;
      if (*col > 60) {
        putc('\n', outtree);
        *col = 0;
      }
    }
    putc(')', outtree);
    (*col)++;
  }
  x = p->v;
  if (x > 0.0)
    w = (long)(0.43429448222 * log(x));
  else if (x == 0.0)
    w = 0;
  else
    w = (long)(0.43429448222 * log(-x)) + 1;
  if (w < 0)
    w = 0;
  if (p != root) {
    fprintf(outtree, ":%*.5f", (int)(w + 7), x);
    *col += w + 8; 
  }
  if (p != root)
    return;
  if (nextree > 2)
    fprintf(outtree, "[%6.4f];\n", 1.0 / (nextree - 1));
  else
    fprintf(outtree, ";\n");
}  /* treeout3 */


/* FIXME curtree should probably be passed by reference */
void drawline2(long i, double scale, tree curtree)
{
  fdrawline2(outfile, i, scale, &curtree);
}

void fdrawline2(FILE *fp, long i, double scale, tree *curtree)
{
  /* draws one row of the tree diagram by moving up tree */
  /* used in dnaml, proml, & restml */
  node *p, *q;
  long n, j;
  boolean extra;
  node *r, *first =NULL, *last =NULL;
  boolean done;

  p = curtree->start;
  q = curtree->start;
  extra = false;
  if (i == (long)p->ycoord && p == curtree->start) {
    if (p->index - spp >= 10)
      fprintf(fp, " %2ld", p->index - spp);
    else
      fprintf(fp, "  %ld", p->index - spp);
    extra = true;
  } else
    fprintf(fp, "  ");
  do {
    if (!p->tip) {
      r = p->next;
      done = false;
      do {
        if (i >= r->back->ymin && i <= r->back->ymax) {
          q = r->back;
          done = true;
        }
        r = r->next;
      } while (!(done || (p != curtree->start && r == p) ||
                 (p == curtree->start && r == p->next)));
      first = p->next->back;
      r = p;
      while (r->next != p)
        r = r->next;
      last = r->back;
      if (p == curtree->start)
        last = p->back;
    }
    done = (p->tip || p == q);
    n = (long)(scale * (q->xcoord - p->xcoord) + 0.5);
    if (n < 3 && !q->tip)
      n = 3;
    if (extra) {
      n--;
      extra = false;
    }
    if ((long)q->ycoord == i && !done) {
      if ((long)p->ycoord != (long)q->ycoord)
        putc('+', fp);
      else
        putc('-', fp);
      if (!q->tip) {
        for (j = 1; j <= n - 2; j++)
          putc('-', fp);
        if (q->index - spp >= 10)
          fprintf(fp, "%2ld", q->index - spp);
        else
          fprintf(fp, "-%ld", q->index - spp);
        extra = true;
      } else {
        for (j = 1; j < n; j++)
          putc('-', fp);
      }
    } else if (!p->tip) {
      if ((long)last->ycoord > i && (long)first->ycoord < i &&
          (i != (long)p->ycoord || p == curtree->start)) {
        putc('|', fp);
        for (j = 1; j < n; j++)
          putc(' ', fp);
      } else {
        for (j = 1; j <= n; j++)
          putc(' ', fp);
      }
    } else {
      for (j = 1; j <= n; j++)
        putc(' ', fp);
    }
    if (q != p)
      p = q;
  } while (!done);
  if ((long)p->ycoord == i && p->tip) {
    for (j = 0; j < nmlngth; j++)
      putc(nayme[p->index-1][j], fp);
  }
  putc('\n', fp);
}  /* drawline2 */


void drawline3(long i, double scale, node *start)
{
  /* draws one row of the tree diagram by moving up tree */
  /* used in dnapars */
  node *p, *q;
  long n, j;
  boolean extra;
  node *r, *first =NULL, *last =NULL;
  boolean done;

  p = start;
  q = start;
  extra = false;
  if (i == (long)p->ycoord) {
    if (p->index - spp >= 10)
      fprintf(outfile, " %2ld", p->index - spp);
    else
      fprintf(outfile, "  %ld", p->index - spp);
    extra = true;
  } else
    fprintf(outfile, "  ");
  do {
    if (!p->tip) {
      r = p->next;
      done = false;
      do {
        if (i >= r->back->ymin && i <= r->back->ymax) {
          q = r->back;
          done = true;
        }
        r = r->next;
      } while (!(done || (r == p))); 
      first = p->next->back;
      r = p;
      while (r->next != p)
        r = r->next;
      last = r->back;
    }
    done = (p->tip || p == q);
    n = (long)(scale * (q->xcoord - p->xcoord) + 0.5);
    if (n < 3 && !q->tip)
      n = 3;
    if (extra) {
      n--;
      extra = false;
    }
    if ((long)q->ycoord == i && !done) {
      if ((long)p->ycoord != (long)q->ycoord)
        putc('+', outfile);
      else
        putc('-', outfile);
      if (!q->tip) {
        for (j = 1; j <= n - 2; j++)
          putc('-', outfile);
        if (q->index - spp >= 10)
          fprintf(outfile, "%2ld", q->index - spp);
        else
          fprintf(outfile, "-%ld", q->index - spp);
        extra = true;
      } else {
        for (j = 1; j < n; j++)
          putc('-', outfile);
      }
    } else if (!p->tip) {
      if ((long)last->ycoord > i && (long)first->ycoord < i &&
          (i != (long)p->ycoord || p == start)) {
        putc('|', outfile);
        for (j = 1; j < n; j++)
          putc(' ', outfile);
      } else {
        for (j = 1; j <= n; j++)
          putc(' ', outfile);
      }
    } else {
      for (j = 1; j <= n; j++)
        putc(' ', outfile);
    }
    if (q != p)
      p = q;
  } while (!done);
  if ((long)p->ycoord == i && p->tip) {
    for (j = 0; j < nmlngth; j++)
      putc(nayme[p->index-1][j], outfile);
  }
  putc('\n', outfile);
}  /* drawline3 */


void copynode(node *c, node *d, long categs)
{
  long i, j;

  for (i = 0; i < endsite; i++)
    for (j = 0; j < categs; j++)
      memcpy(d->x[i][j], c->x[i][j], sizeof(sitelike));
  memcpy(d->underflows,c->underflows,sizeof(double) * endsite);
  d->tyme = c->tyme;
  d->v = c->v;
  d->xcoord = c->xcoord;
  d->ycoord = c->ycoord;
  d->ymin = c->ymin;
  d->ymax = c->ymax;
  d->iter = c->iter;                   /* iter used in dnaml only */
  d->haslength = c->haslength;         /* haslength used in dnamlk only */
  d->initialized = c->initialized;     /* initialized used in dnamlk only */
}  /* copynode */


void prot_copynode(node *c, node *d, long categs)
{
  /* a version of copynode for proml */
  long i, j;

  for (i = 0; i < endsite; i++)
    for (j = 0; j < categs; j++)
      memcpy(d->protx[i][j], c->protx[i][j], sizeof(psitelike));
  memcpy(d->underflows,c->underflows,sizeof(double) * endsite);
  d->tyme = c->tyme;
  d->v = c->v;
  d->xcoord = c->xcoord;
  d->ycoord = c->ycoord;
  d->ymin = c->ymin;
  d->ymax = c->ymax;
  d->iter = c->iter;                   /* iter used in dnaml only */
  d->haslength = c->haslength;         /* haslength used in dnamlk only */
  d->initialized = c->initialized;     /* initialized used in dnamlk only */
}  /* prot_copynode */


void copy_(tree *a, tree *b, long nonodes, long categs)
{
  /* used in dnamlk */
  long i;
  node *p, *q, *r, *s, *t;

  for (i = 0; i < spp; i++) {
    copynode(a->nodep[i], b->nodep[i], categs);
    if (a->nodep[i]->back) {
      if (a->nodep[i]->back == a->nodep[a->nodep[i]->back->index - 1])
        b->nodep[i]->back = b->nodep[a->nodep[i]->back->index - 1];
      else if (a->nodep[i]->back == a->nodep[a->nodep[i]->back->index - 1]->next)
        b->nodep[i]->back = b->nodep[a->nodep[i]->back->index - 1]->next;
      else
        b->nodep[i]->back = b->nodep[a->nodep[i]->back->index - 1]->next->next;
    }
    else b->nodep[i]->back = NULL;
  }
  for (i = spp; i < nonodes; i++) {
    if (a->nodep[i]) {
      p = a->nodep[i];
      q = b->nodep[i];
          r = p;
      do {
        copynode(p, q, categs);
        if (p->back) {
          s = a->nodep[p->back->index - 1];
          t = b->nodep[p->back->index - 1];
          if (s->tip) {
            if(p->back == s)
              q->back = t;
          } else {
            do {
              if (p->back == s)
                q->back = t;
              s = s->next;
              t = t->next;
            } while (s != a->nodep[p->back->index - 1]);
          }
        }
        else
          q->back = NULL;
        p = p->next;
        q = q->next;
      } while (p != r);
    }
  }
  b->likelihood = a->likelihood;
  b->start = a->start;               /* start used in dnaml only */
  b->root = a->root;                 /* root used in dnamlk only */
}  /* copy_ */


void prot_copy_(tree *a, tree *b, long nonodes, long categs)
{
  /* used in promlk */
  /* identical to copy_() except for calls to prot_copynode rather */
  /* than copynode.                                                */
  long i;
  node *p, *q, *r, *s, *t;

  for (i = 0; i < spp; i++) {
    prot_copynode(a->nodep[i], b->nodep[i], categs);
    if (a->nodep[i]->back) {
      if (a->nodep[i]->back == a->nodep[a->nodep[i]->back->index - 1])
        b->nodep[i]->back = b->nodep[a->nodep[i]->back->index - 1];
      else if (a->nodep[i]->back == a->nodep[a->nodep[i]->back->index - 1]->next
) 
        b->nodep[i]->back = b->nodep[a->nodep[i]->back->index - 1]->next;
      else
        b->nodep[i]->back = b->nodep[a->nodep[i]->back->index - 1]->next->next;
    }
    else b->nodep[i]->back = NULL;
  }
  for (i = spp; i < nonodes; i++) {
    if (a->nodep[i]) {
      p = a->nodep[i];
      q = b->nodep[i];
          r = p;
      do {
        prot_copynode(p, q, categs);
        if (p->back) {
          s = a->nodep[p->back->index - 1];
          t = b->nodep[p->back->index - 1];
          if (s->tip)
            {
                if(p->back == s)
                  q->back = t;
          } else {
            do {
              if (p->back == s)
                q->back = t;
              s = s->next;
              t = t->next;
            } while (s != a->nodep[p->back->index - 1]);
          }
        }
        else
          q->back = NULL;
        p = p->next;
        q = q->next;
      } while (p != r);
    }
  }
  b->likelihood = a->likelihood;
  b->start = a->start;               /* start used in dnaml only */
  b->root = a->root;                 /* root used in dnamlk only */
}  /* prot_copy_ */


void standev(long chars, long numtrees, long minwhich, double minsteps,
                        double *nsteps, long **fsteps, longer seed)
{  /* do paired sites test (KHT or SH test) on user-defined trees */
   /* used in dnapars & protpars */
  long i, j, k;
  double wt, sumw, sum, sum2, sd;
  double temp;
  double **covar, *P, *f, *r;

#define SAMPLES 1000
  if (numtrees == 2) {
    fprintf(outfile, "Kishino-Hasegawa-Templeton test\n\n");
    fprintf(outfile, "Tree    Steps   Diff Steps   Its S.D.");
    fprintf(outfile, "   Significantly worse?\n\n");
    which = 1;
    while (which <= numtrees) {
      fprintf(outfile, "%3ld%10.1f", which, nsteps[which - 1] / 10);
      if (minwhich == which)
        fprintf(outfile, "  <------ best\n");
      else {
        sumw = 0.0;
        sum = 0.0;
        sum2 = 0.0;
        for (i = 0; i < endsite; i++) {
          if (weight[i] > 0) {
            wt = weight[i] / 10.0;
            sumw += wt;
            temp = (fsteps[which - 1][i] - fsteps[minwhich - 1][i]) / 10.0;
            sum += wt * temp;
            sum2 += wt * temp * temp;
          }
        }
        sd = sqrt(sumw / (sumw - 1.0) * (sum2 - sum * sum / sumw));
        fprintf(outfile, "%10.1f%12.4f",
                (nsteps[which - 1] - minsteps) / 10, sd);
        if ((sum > 0.0) && (sum > 1.95996 * sd))
          fprintf(outfile, "           Yes\n");
        else
          fprintf(outfile, "           No\n");
      }
      which++;
    }
    fprintf(outfile, "\n\n");
  } else {           /* Shimodaira-Hasegawa test using normal approximation */
    if(numtrees > MAXSHIMOTREES){
      fprintf(outfile, "Shimodaira-Hasegawa test on first %d of %ld trees\n\n"
              , MAXSHIMOTREES, numtrees);
      numtrees = MAXSHIMOTREES;
    } else {
      fprintf(outfile, "Shimodaira-Hasegawa test\n\n");
    }
    covar = (double **)Malloc(numtrees*sizeof(double *));  
    sumw = 0.0;
    for (i = 0; i < endsite; i++)
      sumw += weight[i] / 10.0;
    for (i = 0; i < numtrees; i++)
      covar[i] = (double *)Malloc(numtrees*sizeof(double));  
    for (i = 0; i < numtrees; i++) {        /* compute covariances of trees */
      sum = nsteps[i]/(10.0*sumw);
      for (j = 0; j <=i; j++) {
        sum2 = nsteps[j]/(10.0*sumw);
        temp = 0.0;
        for (k = 0; k < endsite; k++) {
          if (weight[k] > 0) {
            wt = weight[k]/10.0;
            temp = temp + wt*(fsteps[i][k]/10.0-sum)
                            *(fsteps[j][k]/10.0-sum2);
          }
        }
        covar[i][j] = temp;
        if (i != j)
          covar[j][i] = temp;
      }
    }
    for (i = 0; i < numtrees; i++) { /* in-place Cholesky decomposition
                                        of trees x trees covariance matrix */
      sum = 0.0;
      for (j = 0; j <= i-1; j++)
        sum = sum + covar[i][j] * covar[i][j];
      if (covar[i][i] <= sum)
        temp = 0.0;
      else
        temp = sqrt(covar[i][i] - sum);
      covar[i][i] = temp;
      for (j = i+1; j < numtrees; j++) {
        sum = 0.0;
        for (k = 0; k < i; k++)
          sum = sum + covar[i][k] * covar[j][k];
        if (fabs(temp) < 1.0E-12)
          covar[j][i] = 0.0;
        else
          covar[j][i] = (covar[j][i] - sum)/temp;
      }
    }
    f = (double *)Malloc(numtrees*sizeof(double)); /* resampled sums */
    P = (double *)Malloc(numtrees*sizeof(double)); /* vector of P's of trees */
    r = (double *)Malloc(numtrees*sizeof(double)); /* store Normal variates */
    for (i = 0; i < numtrees; i++)
      P[i] = 0.0;
    sum2 = nsteps[0]/10.0;               /* sum2 will be smallest # of steps */
    for (i = 1; i < numtrees; i++)
      if (sum2 > nsteps[i]/10.0)
        sum2 = nsteps[i]/10.0;
    for (i = 1; i <= SAMPLES; i++) {          /* loop over resampled trees */
      for (j = 0; j < numtrees; j++)          /* draw Normal variates */
        r[j] = normrand(seed);
      for (j = 0; j < numtrees; j++) {        /* compute vectors */
        sum = 0.0;
        for (k = 0; k <= j; k++)
          sum += covar[j][k]*r[k];
        f[j] = sum;
      }
      sum = f[1];
      for (j = 1; j < numtrees; j++)          /* get min of vector */
        if (f[j] < sum)
          sum = f[j];
      for (j = 0; j < numtrees; j++)          /* accumulate P's */
        if (nsteps[j]/10.0-sum2 <= f[j] - sum)
          P[j] += 1.0/SAMPLES;
    }
    fprintf(outfile, "Tree    Steps   Diff Steps   P value");
    fprintf(outfile, "   Significantly worse?\n\n");
    for (i = 0; i < numtrees; i++) {
      fprintf(outfile, "%3ld%10.1f", i+1, nsteps[i]/10);
      if ((minwhich-1) == i)
        fprintf(outfile, "  <------ best\n");
      else {
        fprintf(outfile, " %9.1f  %10.3f", nsteps[i]/10.0-sum2, P[i]);
        if (P[i] < 0.05)
          fprintf(outfile, "           Yes\n");
        else
          fprintf(outfile, "           No\n");
      }
    }
  fprintf(outfile, "\n");
  free(P);             /* free the variables we Malloc'ed */
  free(f);
  free(r);
  for (i = 0; i < numtrees; i++)
    free(covar[i]);
  free(covar);
  }
}  /* standev */


void standev2(long numtrees, long maxwhich, long a, long b, double maxlogl,
              double *l0gl, double **l0gf, steptr aliasweight, longer seed)
{  /* do paired sites test (KHT or SH) for user-defined trees */
  /* used in dnaml, dnamlk, proml, promlk, and restml */
  double **covar, *P, *f, *r;
  long i, j, k;
  double wt, sumw, sum, sum2, sd;
  double temp;

#define SAMPLES 1000
  if (numtrees == 2) {
    fprintf(outfile, "Kishino-Hasegawa-Templeton test\n\n");
    fprintf(outfile, "Tree    logL    Diff logL    Its S.D.");
    fprintf(outfile, "   Significantly worse?\n\n");
    which = 1;
    while (which <= numtrees) {
      fprintf(outfile, "%3ld %9.1f", which, l0gl[which - 1]);
      if (maxwhich == which)
        fprintf(outfile, "  <------ best\n");
      else {
        sumw = 0.0;
        sum = 0.0;
        sum2 = 0.0;
        for (i = a; i <= b; i++) {
          if (aliasweight[i] > 0) {
            wt = aliasweight[i];
            sumw += wt;
            temp = l0gf[which - 1][i] - l0gf[maxwhich - 1][i];
            sum += temp * wt;
            sum2 += wt * temp * temp;
          }
        }
        temp = sum / sumw;
        sd = sqrt(sumw / (sumw - 1.0) * (sum2 - sum * sum / sumw ));
        fprintf(outfile, "%10.1f %11.4f", (l0gl[which - 1])-maxlogl, sd);
        if ((sum < 0.0) && ((-sum) > 1.95996 * sd))
          fprintf(outfile, "           Yes\n");
        else
          fprintf(outfile, "           No\n");
      }
      which++;
    }
    fprintf(outfile, "\n\n");
  } else {           /* Shimodaira-Hasegawa test using normal approximation */
    if(numtrees > MAXSHIMOTREES){
      fprintf(outfile, "Shimodaira-Hasegawa test on first %d of %ld trees\n\n"
              , MAXSHIMOTREES, numtrees);
      numtrees = MAXSHIMOTREES;
    } else {
      fprintf(outfile, "Shimodaira-Hasegawa test\n\n");
    }
    covar = (double **)Malloc(numtrees*sizeof(double *));  
    sumw = 0.0;
    for (i = a; i <= b; i++)
      sumw += aliasweight[i];
    for (i = 0; i < numtrees; i++)
      covar[i] = (double *)Malloc(numtrees*sizeof(double));  
    for (i = 0; i < numtrees; i++) {        /* compute covariances of trees */
      sum = l0gl[i]/sumw;
      for (j = 0; j <=i; j++) {
        sum2 = l0gl[j]/sumw;
        temp = 0.0;
        for (k = a; k <= b ; k++) {
          if (aliasweight[k] > 0) {
            wt = aliasweight[k];
            temp = temp + wt*(l0gf[i][k]-sum)
                            *(l0gf[j][k]-sum2);
          }
        }
        covar[i][j] = temp;
        if (i != j)
          covar[j][i] = temp;
      }
    }
    for (i = 0; i < numtrees; i++) { /* in-place Cholesky decomposition
                                        of trees x trees covariance matrix */
      sum = 0.0;
      for (j = 0; j <= i-1; j++)
        sum = sum + covar[i][j] * covar[i][j];
      if (covar[i][i] <= sum)
        temp = 0.0;
      else
        temp = sqrt(covar[i][i] - sum);
      covar[i][i] = temp;
      for (j = i+1; j < numtrees; j++) {
        sum = 0.0;
        for (k = 0; k < i; k++)
          sum = sum + covar[i][k] * covar[j][k];
        if (fabs(temp) < 1.0E-12)
          covar[j][i] = 0.0;
        else
          covar[j][i] = (covar[j][i] - sum)/temp;
      }
    }
    f = (double *)Malloc(numtrees*sizeof(double)); /* resampled likelihoods */
    P = (double *)Malloc(numtrees*sizeof(double)); /* vector of P's of trees */
    r = (double *)Malloc(numtrees*sizeof(double)); /* store Normal variates */
    for (i = 0; i < numtrees; i++)
      P[i] = 0.0;
    for (i = 1; i <= SAMPLES; i++) {          /* loop over resampled trees */
      for (j = 0; j < numtrees; j++)          /* draw Normal variates */
        r[j] = normrand(seed);
      for (j = 0; j < numtrees; j++) {        /* compute vectors */
        sum = 0.0;
        for (k = 0; k <= j; k++)
          sum += covar[j][k]*r[k];
        f[j] = sum;
      }
      sum = f[1];
      for (j = 1; j < numtrees; j++)          /* get max of vector */
        if (f[j] > sum)
          sum = f[j];
      for (j = 0; j < numtrees; j++)          /* accumulate P's */
        if (maxlogl-l0gl[j] <= sum-f[j])
          P[j] += 1.0/SAMPLES;
    }
    fprintf(outfile, "Tree    logL    Diff logL    P value");
    fprintf(outfile, "   Significantly worse?\n\n");
    for (i = 0; i < numtrees; i++) {
      fprintf(outfile, "%3ld%10.1f", i+1, l0gl[i]);
      if ((maxwhich-1) == i)
        fprintf(outfile, "  <------ best\n");
      else {
        fprintf(outfile, " %9.1f  %10.3f", l0gl[i]-maxlogl, P[i]);
        if (P[i] < 0.05)
          fprintf(outfile, "           Yes\n");
        else
          fprintf(outfile, "           No\n");
      }
    }
  fprintf(outfile, "\n");
  free(P);             /* free the variables we Malloc'ed */
  free(f);
  free(r);
  for (i = 0; i < numtrees; i++)
    free(covar[i]);
  free(covar);
  }
}  /* standev */


void freetip(node *anode)
{
  /* used in dnacomp, dnapars, & dnapenny */

  free(anode->numsteps);
  free(anode->oldnumsteps);
  free(anode->base);
  free(anode->oldbase);
}  /* freetip */


void freenontip(node *anode)
{
  /* used in dnacomp, dnapars, & dnapenny */

  free(anode->numsteps);
  free(anode->oldnumsteps);
  free(anode->base);
  free(anode->oldbase);
  free(anode->numnuc);
}  /* freenontip */


void freenodes(long nonodes, pointarray treenode)
{
  /* used in dnacomp, dnapars, & dnapenny */
  long i;
  node *p;

  for (i = 0; i < spp; i++)
    freetip(treenode[i]);
  for (i = spp; i < nonodes; i++) {
    if (treenode[i] != NULL) {
      p = treenode[i]->next;
      do {
        freenontip(p);
        p = p->next;
      } while (p != treenode[i]);
      freenontip(p);
    }
  }
}  /* freenodes */


void freenode(node **anode)
{
  /* used in dnacomp, dnapars, & dnapenny */

  freenontip(*anode);
  free(*anode);
}  /* freenode */


void freetree(long nonodes, pointarray treenode)
{
  /* used in dnacomp, dnapars, & dnapenny */
  long i;
  node *p, *q;

  for (i = 0; i < spp; i++)
    free(treenode[i]);
  for (i = spp; i < nonodes; i++) {
    if (treenode[i] != NULL) {
      p = treenode[i]->next;
      do {
        q = p->next;
        free(p);
        p = q;
      } while (p != treenode[i]);
      free(p);
    }
  }
  free(treenode);
}  /* freetree */


void prot_freex_notip(long nonodes, pointarray treenode)
{
  /* used in proml */
  long i, j;
  node *p;

  for (i = spp; i < nonodes; i++) {
    p = treenode[i];
    if ( p == NULL ) continue;
    do {
      for (j = 0; j < endsite; j++){
        free(p->protx[j]);
        p->protx[j] = NULL;
      }
      free(p->underflows);
      p->underflows = NULL;
      free(p->protx);
      p->protx = NULL;
      p = p->next;
    } while (p != treenode[i]);
  }
}  /* prot_freex_notip */


void prot_freex(long nonodes, pointarray treenode)
{
  /* used in proml */
  long i, j;
  node *p;

  for (i = 0; i < spp; i++) {
    for (j = 0; j < endsite; j++)
      free(treenode[i]->protx[j]);
    free(treenode[i]->protx);
    free(treenode[i]->underflows);
  }
  for (i = spp; i < nonodes; i++) {
    p = treenode[i];
    do {
      for (j = 0; j < endsite; j++)
        free(p->protx[j]);
      free(p->protx);
      free(p->underflows);
      p = p->next;
    } while (p != treenode[i]);
  }
}  /* prot_freex */


void freex_notip(long nonodes, pointarray treenode)
{
  /* used in dnaml & dnamlk */
  long i, j;
  node *p;

  for (i = spp; i < nonodes; i++) {
    p = treenode[i];
    if ( p == NULL ) continue;
    do {
      for (j = 0; j < endsite; j++)
        free(p->x[j]);
      free(p->underflows);
      free(p->x);
      p = p->next;
    } while (p != treenode[i]);
  }
}  /* freex_notip */


void freex(long nonodes, pointarray treenode)
{
  /* used in dnaml & dnamlk */
  long i, j;
  node *p;

  for (i = 0; i < spp; i++) {
    for (j = 0; j < endsite; j++)
      free(treenode[i]->x[j]);
    free(treenode[i]->x);
    free(treenode[i]->underflows);
  }
  for (i = spp; i < nonodes; i++) {
    if(treenode[i]){
      p = treenode[i];
      do {
        for (j = 0; j < endsite; j++)
          free(p->x[j]);
        free(p->x);
        free(p->underflows);
        p = p->next;
      } while (p != treenode[i]);
    }
  }
}  /* freex */



void freegarbage(gbases **garbage)
{
  /* used in dnacomp, dnapars, & dnapenny */
  gbases *p;

  while (*garbage) {
    p = *garbage;
    *garbage = (*garbage)->next;
    free(p->base);
    free(p);
  }
}  /*freegarbage */


void freegrbg(node **grbg)
{
  /* used in dnacomp, dnapars, & dnapenny */
  node *p;

  while (*grbg) {
    p = *grbg;
    *grbg = (*grbg)->next;
    freenontip(p);
    free(p);
  }
} /*freegrbg */


void collapsetree(node *p, node *root, node **grbg, pointarray treenode, 
                  long *zeros)
{
  /*  Recurse through tree searching for zero length brances between */
  /*  nodes (not to tips).  If one exists, collapse the nodes together, */
  /*  removing the branch. */
  node *q, *x1, *y1, *x2, *y2;
  long i, j, index, index2, numd;
  if (p->tip)
    return;
  q = p->next;
  do {
    if (!q->back->tip && q->v == 0.000000) {
      /* merge the two nodes. */
      x1 = y2 = q->next;
      x2 = y1 = q->back->next;
      while(x1->next != q)
        x1 = x1-> next;
      while(y1->next != q->back)
        y1 = y1-> next;
      x1->next = x2;
      y1->next = y2;

      index = q->index;
      index2 = q->back->index;
      numd = treenode[index-1]->numdesc + q->back->numdesc -1;
      chuck(grbg, q->back);
      chuck(grbg, q);
      q = x2;

      /* update the indicies around the node circle */
      do{
        if(q->index != index){
          q->index = index;
        }
        q = q-> next;
      }while(x2 != q);
      updatenumdesc(treenode[index-1], root, numd);
       
      /* Alter treenode to point to real nodes, and update indicies */
      /* acordingly. */
       j = 0; i=0;
      for(i = (index2-1); i < nonodes-1 && treenode[i+1]; i++){ 
        treenode[i]=treenode[i+1];
        treenode[i+1] = NULL;
        x1=x2=treenode[i]; 
        do{ 
          x1->index = i+1; 
          x1 = x1 -> next; 
        } while(x1 != x2); 
      }

      /* Create a new empty fork in the blank spot of treenode */
      x1=NULL;
      for(i=1; i <=3 ; i++){
        gnutreenode(grbg, &x2, index2, endsite, zeros);
        x2->next = x1;
        x1 = x2;
      }
      x2->next->next->next = x2;
      treenode[nonodes-1]=x2;
      if (q->back)
        collapsetree(q->back, root, grbg, treenode, zeros);
    } else {
      if (q->back)
        collapsetree(q->back, root, grbg, treenode, zeros);
      q = q->next;
    }
  } while (q != p);
} /* collapsetree */


void collapsebestrees(node **root, node **grbg, pointarray treenode, 
                      bestelm *bestrees, long *place, long *zeros,
                      long chars, boolean recompute, boolean progress)
     
{
  /* Goes through all best trees, collapsing trees where possible, and  */
  /* deleting trees that are not unique.    */
  long i,j, k, pos, nextnode, oldnextree;
  boolean found;
  node *dummy;

  oldnextree = nextree;
  for(i = 0 ; i < (oldnextree - 1) ; i++){
    bestrees[i].collapse = true;
  }

  if(progress)
    printf("Collapsing best trees\n   ");
  k = 0;
  for(i = 0 ; i < (oldnextree - 1) ; i++){
    if(progress){
      if(i % (((oldnextree-1) / 72) + 1) == 0)
        putchar('.');
      fflush(stdout);
    }
    while(!bestrees[k].collapse)
      k++;
    /* Reconstruct tree. */
    *root = treenode[0];
    add(treenode[0], treenode[1], treenode[spp], root, recompute,
        treenode, grbg, zeros);
    nextnode = spp + 2;
    for (j = 3; j <= spp; j++) {
      if (bestrees[k].btree[j - 1] > 0)
        add(treenode[bestrees[k].btree[j - 1] - 1], treenode[j - 1],
            treenode[nextnode++ - 1], root, recompute, treenode, grbg,
            zeros);
      else
          add(treenode[treenode[-bestrees[k].btree[j - 1]-1]->back->index-1],
              treenode[j - 1], NULL, root, recompute, treenode, grbg, zeros);
    }
    reroot(treenode[outgrno - 1], *root);

    treelength(*root, chars, treenode);
    collapsetree(*root, *root, grbg, treenode, zeros);
    savetree(*root, place, treenode, grbg, zeros);
    /* move everything down in the bestree list */
    for(j = k ; j < (nextree - 2) ; j++){
      memcpy(bestrees[j].btree, bestrees[j + 1].btree, spp * sizeof(long));
      bestrees[j].gloreange = bestrees[j + 1].gloreange;
      bestrees[j + 1].gloreange = false;
      bestrees[j].locreange = bestrees[j + 1].locreange;
      bestrees[j + 1].locreange = false;
      bestrees[j].collapse = bestrees[j + 1].collapse;
    }
    pos=0;
    findtree(&found, &pos, nextree-1, place, bestrees);    

    /* put the new tree at the end of the list if it wasn't found */
    nextree--;
    if(!found)
      addtree(pos, &nextree, false, place, bestrees);

    /* Deconstruct the tree */
    for (j = 1; j < spp; j++){
      re_move(treenode[j], &dummy, root, recompute, treenode,
              grbg, zeros);
    }
  }
  if (progress) {
    putchar('\n');
#ifdef WIN32
    phyFillScreenColor();
#endif
  }
}
