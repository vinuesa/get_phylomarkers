
#include "phylip.h"
#include "discrete.h"

/* version 3.6 (c) Copyright 1993-2004 by the University of Washington.
   Written by Joseph Felsenstein, Akiko Fuseki, Sean Lamont, and Andrew Keeffe.
   Permission is granted to copy and use this program provided no fee is
   charged for it and provided that this copyright notice is not removed. */

#define MAXNUMTREES     1000000  /* bigger than number of user trees can be */

#ifndef OLDC
/* function prototypes */
void   getoptions(void);
void   allocrest(void);
void   doinit(void);
void   makeweights(void);
void   doinput(void);
void   initparsnode(node **, node **, node *, long, long, long *, long *,
                initops, pointarray, pointarray, Char *, Char *, FILE *);
void   evaluate(node *);
void   tryadd(node *, node *, node *);
void   addpreorder(node *, node *, node *);
void   trydescendants(node *, node *, node *, node *, boolean);

void   trylocal(node *, node *);
void   trylocal2(node *, node *, node *);
void   tryrearr(node *p, boolean *);
void   repreorder(node *p, boolean *);
void   rearrange(node **);
void   describe(void);
void   pars_coordinates(node *, double, long *, double *);
void   pars_printree(void);
void   globrearrange(void);
void   grandrearr(void);

void   maketree(void);
void   freerest(void);
void   load_tree(long treei);
void   reallocchars(void);
/* function prototypes */
#endif


Char infilename[FNMLNGTH], outfilename[FNMLNGTH], intreename[FNMLNGTH], outtreename[FNMLNGTH],
     weightfilename[FNMLNGTH];
node *root;
long chars, col, msets, ith, njumble, jumb, maxtrees;
/*   chars = number of sites in actual sequences */
long inseed, inseed0;
double threshold;
boolean jumble, usertree, thresh, weights, thorough, rearrfirst, trout,
          progress, stepbox, ancseq, mulsets, justwts, firstset, mulf, multf;
steptr oldweight;
longer seed;
pointarray treenode;            /* pointers to all nodes in tree */
long *enterorder;
char *progname;
long *zeros;
unsigned char *zeros2;

/* local variables for Pascal maketree, propagated globally for C version: */

long minwhich;
double like, minsteps, bestyet, bestlike, bstlike2;
boolean lastrearr, recompute;
double nsteps[maxuser];
long **fsteps;
node *there, *oldnufork;
long *place;
bestelm *bestrees;
long *threshwt;
discbaseptr nothing;
gbases *garbage;
node *temp, *temp1, *temp2, *tempsum, *temprm, *tempadd, *tempf, *tmp, *tmp1,
       *tmp2, *tmp3, *tmprm, *tmpadd;
boolean *names;
node *grbg;


void getoptions()
{
  /* interactively set options */
  long inseed0, loopcount, loopcount2;
  Char ch, ch2;

  fprintf(outfile, "\nDiscrete character parsimony algorithm, version %s\n\n",
            VERSION);
  jumble = false;
  njumble = 1;
  outgrno = 1;
  outgropt = false;
  thresh = false;
  thorough = true;
  rearrfirst = false;
  maxtrees = 100;
  trout = true;
  usertree = false;
  weights = false;
  mulsets = false;
  printdata = false;
  progress = true;
  treeprint = true;
  stepbox = false;
  ancseq = false;
  dotdiff = true;
  interleaved = true;
  loopcount = 0;
  for (;;) {
    cleerhome();
    printf("\nDiscrete character parsimony algorithm, version %s\n\n",VERSION);
    printf("Setting for this run:\n");
    printf("  U                 Search for best tree?  %s\n",
           (usertree ? "No, use user trees in input file" : "Yes"));
    if (!usertree) {
      printf("  S                        Search option?  ");
      if (thorough)
        printf("More thorough search\n");
      else if (rearrfirst)
        printf("Rearrange on one best tree\n");
      else
        printf("Less thorough\n");
      printf("  V              Number of trees to save?  %ld\n", maxtrees);
      printf("  J     Randomize input order of species?");
      if (jumble)
        printf("  Yes (seed =%8ld,%3ld times)\n", inseed0, njumble);
      else
        printf("  No. Use input order\n");
    }
    printf("  O                        Outgroup root?");
    if (outgropt)
      printf("  Yes, at species number %ld\n", outgrno);
    else
      printf("  No, use as outgroup species %ld\n", outgrno);
    printf("  T              Use Threshold parsimony?");
    if (thresh)
      printf("  Yes, count steps up to%4.1f per site\n", threshold);
    else
      printf("  No, use ordinary parsimony\n");
    printf("  W                       Sites weighted?  %s\n",
           (weights ? "Yes" : "No"));
    printf("  M           Analyze multiple data sets?");
    if (mulsets)
      printf("  Yes, %2ld %s\n", msets,
               (justwts ? "sets of weights" : "data sets"));
    else
      printf("  No\n");
    printf("  I            Input species interleaved?  %s\n",
           (interleaved ? "Yes" : "No, sequential"));
    printf("  0   Terminal type (IBM PC, ANSI, none)?  %s\n",
           ibmpc ? "IBM PC" : ansi  ? "ANSI"  : "(none)");
    printf("  1    Print out the data at start of run  %s\n",
           (printdata ? "Yes" : "No"));
    printf("  2  Print indications of progress of run  %s\n",
           progress ? "Yes" : "No");
    printf("  3                        Print out tree  %s\n",
           treeprint ? "Yes" : "No");
    printf("  4          Print out steps in each site  %s\n",
           stepbox ? "Yes" : "No");
    printf("  5  Print character at all nodes of tree  %s\n",
           ancseq ? "Yes" : "No");
    if (ancseq || printdata)
      printf("  .  Use dot-differencing to display them  %s\n",
           dotdiff ? "Yes" : "No");
    printf("  6       Write out trees onto tree file?  %s\n",
           trout ? "Yes" : "No");
    printf("\n  Y to accept these or type the letter for one to change\n");
#ifdef WIN32
    phyFillScreenColor();
#endif
    fflush(stdout);
    scanf("%c%*[^\n]", &ch);
    getchar();
    if (ch == '\n')
      ch = ' ';
    uppercase(&ch);
    if (ch == 'Y')
      break;
    if (((!usertree) && (strchr("WSVJOTUMI12345.60", ch) != NULL))
          || (usertree && (strchr("WSVOTUMI12345.60", ch) != NULL))){
      switch (ch) {
        
      case 'J':
        jumble = !jumble;
        if (jumble)
          initjumble(&inseed, &inseed0, seed, &njumble);
        else njumble = 1;
        break;
        
      case 'O':
        outgropt = !outgropt;
        if (outgropt)
          initoutgroup(&outgrno, spp);
        break;
        
      case 'T':
        thresh = !thresh;
        if (thresh)
          initthreshold(&threshold);
        break;
        
      case 'W':
        weights = !weights;
        break;

      case 'M':
        mulsets = !mulsets;
        if (mulsets) {
          printf("Multiple data sets or multiple weights?");
          loopcount2 = 0;
          do {
            printf(" (type D or W)\n");
#ifdef WIN32
            phyFillScreenColor();
#endif
            fflush(stdout);
            scanf("%c%*[^\n]", &ch2);
            getchar();
            if (ch2 == '\n')
              ch2 = ' ';
            uppercase(&ch2);
            countup(&loopcount2, 10);
          } while ((ch2 != 'W') && (ch2 != 'D'));
          justwts = (ch2 == 'W');
          if (justwts)
            justweights(&msets);
          else
            initdatasets(&msets);
          if (!jumble) {
            jumble = true;
            initjumble(&inseed, &inseed0, seed, &njumble);
          }
        }
        break;
        
      case 'U':
        usertree = !usertree;
        break;
        
      case 'S':
        thorough = !thorough;
        if (!thorough) {
          printf("Rearrange on just one best tree?");
          loopcount2 = 0;
          do {
            printf(" (type Y or N)\n");
#ifdef WIN32
            phyFillScreenColor();
#endif
            fflush(stdout);
            scanf("%c%*[^\n]", &ch2);
            getchar();
            if (ch2 == '\n')
              ch2 = ' ';
            uppercase(&ch2);
            countup(&loopcount2, 10);
          } while ((ch2 != 'Y') && (ch2 != 'N'));
          rearrfirst = (ch2 == 'Y');
        }
        break;

      case 'V':
        loopcount2 = 0;
        do {
          printf("type the number of trees to save\n");
#ifdef WIN32
          phyFillScreenColor();
#endif
          fflush(stdout);
          scanf("%ld%*[^\n]", &maxtrees);
          getchar();
          if (maxtrees  > MAXNUMTREES)
            maxtrees = MAXNUMTREES;
          countup(&loopcount2, 10);
        } while (maxtrees < 1);
        break;

      case 'I':
        interleaved = !interleaved;
        break;
        
      case '0':
        initterminal(&ibmpc, &ansi);
        break;
        
      case '1':
        printdata = !printdata;
        break;
        
      case '2':
        progress = !progress;
        break;
        
      case '3':
        treeprint = !treeprint;
        break;
        
      case '4':
        stepbox = !stepbox;
        break;
        
      case '5':
        ancseq = !ancseq;
        break;
        
      case '.':
        dotdiff = !dotdiff;
        break;
        
      case '6':
        trout = !trout;
        break;
      }
    } else
      printf("Not a possible option!\n");
    countup(&loopcount, 100);
  }
}  /* getoptions */


void reallocchars() 
{
  long i;

  for (i = 0; i < spp; i++) {
    free(y[i]);
    y[i] = (Char *)Malloc(chars*sizeof(Char));
  }
  for (i = 0; i < spp; i++){
    free(convtab[i]);
    convtab[i] = (Char *)Malloc(chars*sizeof(Char));
  }
  
  free(weight);
  free(oldweight);
  free(alias);
  free(ally);
  free(location);

  weight = (long *)Malloc(chars*sizeof(long));
  oldweight = (long *)Malloc(chars*sizeof(long));
  alias = (long *)Malloc(chars*sizeof(long));
  ally = (long *)Malloc(chars*sizeof(long));
  location = (long *)Malloc(chars*sizeof(long));
}

void allocrest()
{
  long i;

  y = (Char **)Malloc(spp*sizeof(Char *));
  for (i = 0; i < spp; i++)
    y[i] = (Char *)Malloc(chars*sizeof(Char));
  convtab = (Char **)Malloc(spp*sizeof(Char *));
  for (i = 0; i < spp; i++)
    convtab[i] = (Char *)Malloc(chars*sizeof(Char));
  bestrees = (bestelm *)Malloc(maxtrees*sizeof(bestelm));
  for (i = 1; i <= maxtrees; i++)
    bestrees[i - 1].btree = (long *)Malloc(nonodes*sizeof(long));
  nayme = (naym *)Malloc(spp*sizeof(naym));
  enterorder = (long *)Malloc(spp*sizeof(long));
  place = (long *)Malloc(nonodes*sizeof(long));
  weight = (long *)Malloc(chars*sizeof(long));
  oldweight = (long *)Malloc(chars*sizeof(long));
  alias = (long *)Malloc(chars*sizeof(long));
  ally = (long *)Malloc(chars*sizeof(long));
  location = (long *)Malloc(chars*sizeof(long));
}  /* alocrest */


void doinit()
{
  /* initializes variables */

  inputnumbers(&spp, &chars, &nonodes, 1);
  getoptions();
  if (printdata)
    fprintf(outfile, "%2ld species, %3ld  sites\n\n", spp, chars);
  alloctree(&treenode, nonodes, usertree);
  allocrest();
}  /* doinit */


void makeweights()
{
  /* make up weights vector to avoid duplicate computations */
  long i;

  for (i = 1; i <= chars; i++) {
    alias[i - 1] = i;
    oldweight[i - 1] = weight[i - 1];
    ally[i - 1] = i;
  }
  sitesort(chars, weight);
  sitecombine(chars);
  sitescrunch(chars);
  endsite = 0;
  for (i = 1; i <= chars; i++) {
    if (ally[i - 1] == i)
      endsite++;
  }
  for (i = 1; i <= endsite; i++)
    location[alias[i - 1] - 1] = i;
  if (!thresh)
    threshold = spp;
  threshwt = (long *)Malloc(endsite*sizeof(long));
  for (i = 0; i < endsite; i++) {
    weight[i] *= 10;
    threshwt[i] = (long)(threshold * weight[i] + 0.5);
  }
  zeros = (long *)Malloc(endsite*sizeof(long));
  for (i = 0; i < endsite; i++)
    zeros[i] = 0;
  zeros2 = (unsigned char *)Malloc(endsite*sizeof(unsigned char));
  for (i = 0; i < endsite; i++)
    zeros2[i] = 0;
}  /* makeweights */


void doinput()
{
  /* reads the input data */
  long i;

  if (justwts) {
    if (firstset)
      inputdata(chars);
    for (i = 0; i < chars; i++)
      weight[i] = 1;
    inputweights(chars, weight, &weights);
    if (justwts) {
      fprintf(outfile, "\n\nWeights set # %ld:\n\n", ith);
      if (progress)
        printf("\nWeights set # %ld:\n\n", ith);
    }
    if (printdata)
      printweights(outfile, 0, chars, weight, "Sites");
  } else {
    if (!firstset) {
      samenumsp(&chars, ith);
      reallocchars();
    }
    inputdata(chars);
    for (i = 0; i < chars; i++)
      weight[i] = 1;
    if (weights) {
      inputweights(chars, weight, &weights);
      if (printdata)
        printweights(outfile, 0, chars, weight, "Sites");
    }
  }
  makeweights();
  makevalues(treenode, zeros, zeros2, usertree);
  if (!usertree) {
    allocdiscnode(&temp, zeros, zeros2, endsite);
    allocdiscnode(&temp1, zeros, zeros2, endsite);
    allocdiscnode(&temp2, zeros, zeros2, endsite);
    allocdiscnode(&tempsum, zeros, zeros2, endsite);
    allocdiscnode(&temprm, zeros, zeros2, endsite);
    allocdiscnode(&tempadd, zeros, zeros2, endsite);
    allocdiscnode(&tempf, zeros, zeros2, endsite);
    allocdiscnode(&tmp, zeros, zeros2, endsite);
    allocdiscnode(&tmp1, zeros, zeros2, endsite);
    allocdiscnode(&tmp2, zeros, zeros2, endsite);
    allocdiscnode(&tmp3, zeros, zeros2, endsite);
    allocdiscnode(&tmprm, zeros, zeros2, endsite);
    allocdiscnode(&tmpadd, zeros, zeros2, endsite);
  }
}  /* doinput */


void initparsnode(node **p, node **grbg, node *q, long len, long nodei,
                        long *ntips, long *parens, initops whichinit,
                        pointarray treenode, pointarray nodep, Char *str, Char *ch,
                        FILE *intree)
{
  /* initializes a node */
  boolean minusread;
  double valyew, divisor;

  switch (whichinit) {
  case bottom:
    gnudisctreenode(grbg, p, nodei, endsite, zeros, zeros2);
    treenode[nodei - 1] = *p;
    break;
  case nonbottom:
    gnudisctreenode(grbg, p, nodei, endsite, zeros, zeros2);
    break;
  case tip:
    match_names_to_data (str, treenode, p, spp);
    break;
  case length:         /* if there is a length, read it and discard value */
    processlength(&valyew, &divisor, ch, &minusread, intree, parens);
    break;
  default:      /*cases hslength,hsnolength,treewt,unittrwt,iter,*/
    break;      /*length should never occur                         */
  }
} /* initparsnode */


void evaluate(node *r)
{
  /* determines the number of steps needed for a tree. this is
     the minimum number of steps needed to evolve sequences on
     this tree */
  long i, steps;
  long term;
  double sum;

  sum = 0.0;
  for (i = 0; i < endsite; i++) {
    steps = r->numsteps[i];
    if ((long)steps <= threshwt[i])
      term = steps;
    else
      term = threshwt[i];
    sum += (double)term;
    if (usertree && which <= maxuser)
      fsteps[which - 1][i] = term;
  }
  if (usertree && which <= maxuser) {
    nsteps[which - 1] = sum;
    if (which == 1) {
      minwhich = 1;
      minsteps = sum;
    } else if (sum < minsteps) {
      minwhich = which;
      minsteps = sum;
    }
  }
  like = -sum;
}  /* evaluate */


void tryadd(node *p, node *item, node *nufork)
{
  /* temporarily adds one fork and one tip to the tree.
     if the location where they are added yields greater
     "likelihood" than other locations tested up to that
     time, then keeps that location as there */
  long pos;
  double belowsum, parentsum;
  boolean found, collapse, changethere, trysave;

  if (!p->tip) {
    memcpy(temp->discbase, p->discbase, endsite*sizeof(unsigned char));
    memcpy(temp->numsteps, p->numsteps, endsite*sizeof(long));
    memcpy(temp->discnumnuc, p->discnumnuc, endsite*sizeof(discnucarray));
    temp->numdesc = p->numdesc + 1;
    if (p->back) {
      multifillin(temp, tempadd, 1);
      sumnsteps2(tempsum, temp, p->back, 0, endsite, threshwt);
    } else {
      multisumnsteps(temp, tempadd, 0, endsite, threshwt);
      tempsum->sumsteps = temp->sumsteps;
    }
    if (tempsum->sumsteps <= -bestyet) {
      if (p->back)
        sumnsteps2(tempsum, temp, p->back, endsite+1, endsite, threshwt);
      else {
        multisumnsteps(temp, temp1, endsite+1, endsite, threshwt);
        tempsum->sumsteps = temp->sumsteps;
      }
    }
    p->sumsteps = tempsum->sumsteps;
  }
  if (p == root)
    sumnsteps2(temp, item, p, 0, endsite, threshwt);
  else {
    sumnsteps(temp1, item, p, 0, endsite);
    sumnsteps2(temp, temp1, p->back, 0, endsite, threshwt);
  }
  if (temp->sumsteps <= -bestyet) {
    if (p == root)
      sumnsteps2(temp, item, p, endsite+1, endsite, threshwt);
    else {
      sumnsteps(temp1, item, p, endsite+1, endsite);
      sumnsteps2(temp, temp1, p->back, endsite+1, endsite, threshwt);
    }
  }
  belowsum = temp->sumsteps;
  multf = false;
  like = -belowsum;
  if (!p->tip && belowsum >= p->sumsteps) {
    multf = true;
    like = -p->sumsteps;
  }
  trysave = true;
  if (!multf && p != root) {
    parentsum = treenode[p->back->index - 1]->sumsteps;
    if (belowsum >= parentsum)
      trysave = false;
  }
  if (lastrearr) {
    changethere = true;
    if (like >= bstlike2 && trysave) {
      if (like > bstlike2)
        found = false;
      else {
        addnsave(p, item, nufork, &root, &grbg, multf, treenode, place,
                   zeros, zeros2);
        pos = 0;
        findtree(&found, &pos, nextree, place, bestrees);
      }
      if (!found) {
        collapse = collapsible(item, p, temp, temp1, temp2, tempsum, temprm,
                     tmpadd, multf, root, zeros, zeros2, treenode);
        if (!thorough)
          changethere = !collapse;
        if (thorough || !collapse || like > bstlike2 || (nextree == 1)) {
          if (like > bstlike2) {
            addnsave(p, item, nufork, &root, &grbg, multf, treenode, 
                       place, zeros, zeros2);
            bestlike = bstlike2 = like;
            addbestever(&pos, &nextree, maxtrees, collapse, place, bestrees);
          } else
            addtiedtree(pos, &nextree, maxtrees, collapse, place, bestrees);
        }
      }
    }
    if (like >= bestyet) {
      if (like > bstlike2)
        bstlike2 = like;
      if (changethere && trysave) {
        bestyet = like;
        there = p;
        mulf = multf;
      }
    }
  } else if ((like > bestyet) || (like >= bestyet && trysave)) {
    bestyet = like;
    there = p;
    mulf = multf;
  }
}  /* tryadd */


void addpreorder(node *p, node *item, node *nufork)
{
  /* traverses a n-ary tree, calling PROCEDURE tryadd
     at a node before calling tryadd at its descendants */
  node *q;

  if (p == NULL)
    return;
  tryadd(p, item, nufork);
  if (!p->tip) {
    q = p->next;
    while (q != p) {
      addpreorder(q->back, item, nufork);
      q = q->next;
    }
  }
}  /* addpreorder */


void trydescendants(node *item, node *forknode, node *parent,
                        node *parentback, boolean trybelow)
{
  /* tries rearrangements at parent and below parent's descendants */
  node *q, *tempblw;
  boolean bestever=0, belowbetter, multf=0, saved, trysave;
  double parentsum=0, belowsum;

  memcpy(temp->discbase, parent->discbase, endsite*sizeof(unsigned char));
  memcpy(temp->numsteps, parent->numsteps, endsite*sizeof(long));
  memcpy(temp->discnumnuc, parent->discnumnuc, endsite*sizeof(discnucarray));
  temp->numdesc = parent->numdesc + 1;
  multifillin(temp, tempadd, 1);
  sumnsteps2(tempsum, parentback, temp, 0, endsite, threshwt);
  belowbetter = true;
  if (lastrearr) {
    parentsum = tempsum->sumsteps;
    if (-tempsum->sumsteps >= bstlike2) {
      belowbetter = false;
      bestever = false;
      multf = true;
      if (-tempsum->sumsteps > bstlike2)
        bestever = true;
      savelocrearr(item, forknode, parent, tmp, tmp1, tmp2, tmp3, tmprm,
                    tmpadd, &root, maxtrees, &nextree, multf, bestever,
                    &saved, place, bestrees, treenode, &grbg, zeros, zeros2);
      if (saved) {
        like = bstlike2 = -tempsum->sumsteps;
        there = parent;
        mulf = true;
      }
    }
  } else if (-tempsum->sumsteps >= like) {
    there = parent;
    mulf = true;
    like = -tempsum->sumsteps;
  }
  if (trybelow) {
    sumnsteps(temp, parent, tempadd, 0, endsite);
    sumnsteps2(tempsum, temp, parentback, 0, endsite, threshwt);
    if (lastrearr) {
      belowsum = tempsum->sumsteps;
      if (-tempsum->sumsteps >= bstlike2 && belowbetter && 
            (forknode->numdesc > 2 ||
               (forknode->numdesc == 2 && 
                 parent->back->index != forknode->index))) {
        trysave = false;
        memcpy(temp->discbase, parentback->discbase, endsite*sizeof(unsigned char));
        memcpy(temp->numsteps, parentback->numsteps, endsite*sizeof(long));
        memcpy(temp->discnumnuc, parentback->discnumnuc, endsite*sizeof(discnucarray));
        temp->numdesc = parentback->numdesc + 1;
        multifillin(temp, tempadd, 1);
        sumnsteps2(tempsum, parent, temp, 0, endsite, threshwt);
        if (-tempsum->sumsteps < bstlike2) {
          multf = false;
          bestever = false;
          trysave = true;
        }
        if (-belowsum > bstlike2) {
          multf = false;
          bestever = true;
          trysave = true;
        }
        if (trysave) {
          if (treenode[parent->index - 1] != parent)
            tempblw = parent->back;
          else
            tempblw = parent;
          savelocrearr(item, forknode, tempblw, tmp, tmp1, tmp2, tmp3, tmprm,
                         tmpadd, &root, maxtrees, &nextree, multf, bestever,
                         &saved, place, bestrees, treenode, &grbg,
                         zeros, zeros2);
          if (saved) {
            like = bstlike2 = -belowsum;
            there = tempblw;
            mulf = false;
          }
        }
      }
    } else if (-tempsum->sumsteps > like) {
      like = -tempsum->sumsteps;
      if (-tempsum->sumsteps > bestyet) {
        if (treenode[parent->index - 1] != parent)
          tempblw = parent->back;
        else
          tempblw = parent;
        there = tempblw;
        mulf = false;
      }
    }
  }
  q = parent->next;
  while (q != parent) {
    if (q->back && q->back != item) {
      memcpy(temp1->discbase, q->discbase, endsite*sizeof(unsigned char));
      memcpy(temp1->numsteps, q->numsteps, endsite*sizeof(long));
      memcpy(temp1->discnumnuc, q->discnumnuc, endsite*sizeof(discnucarray));
      temp1->numdesc = q->numdesc;
      multifillin(temp1, parentback, 0);
      if (lastrearr)
        belowbetter = (-parentsum < bstlike2);
      if (!q->back->tip) {
        memcpy(temp->discbase, q->back->discbase, endsite*sizeof(unsigned char));
        memcpy(temp->numsteps, q->back->numsteps, endsite*sizeof(long));
        memcpy(temp->discnumnuc, q->back->discnumnuc, endsite*sizeof(discnucarray));
        temp->numdesc = q->back->numdesc + 1;
        multifillin(temp, tempadd, 1);
        sumnsteps2(tempsum, temp1, temp, 0, endsite, threshwt);
        if (lastrearr) {
          if (-tempsum->sumsteps >= bstlike2) {
            belowbetter = false;
            bestever = false;
            multf = true;
            if (-tempsum->sumsteps > bstlike2)
              bestever = true;
            savelocrearr(item, forknode, q->back, tmp, tmp1, tmp2, tmp3, tmprm,
                          tmpadd, &root, maxtrees, &nextree, multf, bestever,
                          &saved, place, bestrees, treenode, &grbg,
                          zeros, zeros2);
            if (saved) {
              like = bstlike2 = -tempsum->sumsteps;
              there = q->back;
              mulf = true;
            }
          }
        } else if (-tempsum->sumsteps >= like) {
          like = -tempsum->sumsteps;
          there = q->back;
          mulf = true;
        }
      }
      sumnsteps(temp, q->back, tempadd, 0, endsite);
      sumnsteps2(tempsum, temp, temp1, 0, endsite, threshwt);
      if (lastrearr) {
        if (-tempsum->sumsteps >= bstlike2) {
          trysave = false;
          multf = false;
          if (belowbetter) {
            bestever = false;
            trysave = true;
          }
          if (-tempsum->sumsteps > bstlike2) {
            bestever = true;
            trysave = true;
          }
          if (trysave) {
            if (treenode[q->back->index - 1] != q->back)
              tempblw = q;
            else
              tempblw = q->back;
            savelocrearr(item, forknode, tempblw, tmp, tmp1, tmp2, tmp3, tmprm,
                        tmpadd, &root, maxtrees, &nextree, multf, bestever,
                        &saved, place, bestrees, treenode, &grbg,
                        zeros, zeros2);
            if (saved) {
              like = bstlike2 = -tempsum->sumsteps;
              there = tempblw;
              mulf = false;
            }
          }
        }
      } else if (-tempsum->sumsteps > like) {
        like = -tempsum->sumsteps;
        if (-tempsum->sumsteps > bestyet) {
          if (treenode[q->back->index - 1] != q->back)
            tempblw = q;
          else
            tempblw = q->back;
          there = tempblw;
          mulf = false;
        }
      }
    }
    q = q->next;
  }
} /* trydescendants */


void trylocal(node *item, node *forknode)
{
  /* rearranges below forknode, below descendants of forknode when
     there are more than 2 descendants, then unroots the back of
     forknode and rearranges on its descendants */
  node *q;
  boolean bestever, multf, saved;

  memcpy(temprm->discbase, zeros2, endsite*sizeof(unsigned char));
  memcpy(temprm->numsteps, zeros, endsite*sizeof(long));
  memcpy(temprm->olddiscbase, item->discbase, endsite*sizeof(unsigned char));
  memcpy(temprm->oldnumsteps, item->numsteps, endsite*sizeof(long));
  memcpy(tempf->discbase, forknode->discbase, endsite*sizeof(unsigned char));
  memcpy(tempf->numsteps, forknode->numsteps, endsite*sizeof(long));
  memcpy(tempf->discnumnuc, forknode->discnumnuc, endsite*sizeof(discnucarray));
  tempf->numdesc = forknode->numdesc - 1;
  multifillin(tempf, temprm, -1);
  if (!forknode->back) {
    sumnsteps2(tempsum, tempf, tempadd, 0, endsite, threshwt);
    if (lastrearr) {
      if (-tempsum->sumsteps > bstlike2) {
        bestever = true;
        multf = false;
        savelocrearr(item, forknode, forknode, tmp, tmp1, tmp2, tmp3, tmprm,
                       tmpadd, &root, maxtrees, &nextree, multf, bestever,
                       &saved, place, bestrees, treenode, &grbg,
                       zeros, zeros2);
        if (saved) {
          like = bstlike2 = -tempsum->sumsteps;
          there = forknode;
          mulf = false;
        }
      }
    } else if (-tempsum->sumsteps > like) {
      like = -tempsum->sumsteps;
      if (-tempsum->sumsteps > bestyet) {
        there = forknode;
        mulf = false;
      }
    }
  } else {
    sumnsteps(temp, tempf, tempadd, 0, endsite);
    sumnsteps2(tempsum, temp, forknode->back, 0, endsite, threshwt);
    if (lastrearr) {
      if (-tempsum->sumsteps > bstlike2) {
        bestever = true;
        multf = false;
        savelocrearr(item, forknode, forknode, tmp, tmp1, tmp2, tmp3, tmprm,
                       tmpadd, &root, maxtrees, &nextree, multf, bestever,
                       &saved, place, bestrees, treenode, &grbg,
                       zeros, zeros2);
        if (saved) {
          like = bstlike2 = -tempsum->sumsteps;
          there = forknode;
          mulf = false;
        }
      }
    } else if (-tempsum->sumsteps > like) {
      like = -tempsum->sumsteps;
      if (-tempsum->sumsteps > bestyet) {
        there = forknode;
        mulf = false;
      }
    }
    trydescendants(item, forknode, forknode->back, tempf, false);
  }
  q = forknode->next;
  while (q != forknode) {
    if (q->back != item) {
      memcpy(temp2->discbase, q->discbase, endsite*sizeof(unsigned char));
      memcpy(temp2->numsteps, q->numsteps, endsite*sizeof(long));
      memcpy(temp2->discnumnuc, q->discnumnuc, endsite*sizeof(discnucarray));
      temp2->numdesc = q->numdesc - 1;
      multifillin(temp2, temprm, -1);
      if (!q->back->tip) {
        trydescendants(item, forknode, q->back, temp2, true);
      } else {
        sumnsteps(temp1, q->back, tempadd, 0, endsite);
        sumnsteps2(tempsum, temp1, temp2, 0, endsite, threshwt);
        if (lastrearr) {
          if (-tempsum->sumsteps > bstlike2) {
            multf = false;
            bestever = true;
            savelocrearr(item, forknode, q->back, tmp, tmp1, tmp2, tmp3,
                         tmprm, tmpadd, &root, maxtrees, &nextree, multf,
                         bestever, &saved, place, bestrees, treenode,
                         &grbg, zeros, zeros2);
            if (saved) {
              like = bstlike2 = -tempsum->sumsteps;
              there = q->back;
              mulf = false;
            }
          }
        } else if (-tempsum->sumsteps > like) {
          like = -tempsum->sumsteps;
          if (-tempsum->sumsteps > bestyet) {
            there = q->back;
            mulf = false;
          }
        }
      }
    }
    q = q->next;
  }
} /* trylocal */


void trylocal2(node *item, node *forknode, node *other)
{
  /* rearranges below forknode, below descendants of forknode when
     there are more than 2 descendants, then unroots the back of
     forknode and rearranges on its descendants.  Used if forknode
     has binary descendants */
  node *q;
  boolean bestever=0, multf, saved, belowbetter, trysave;

  memcpy(tempf->discbase, other->discbase, endsite*sizeof(unsigned char));
  memcpy(tempf->numsteps, other->numsteps, endsite*sizeof(long));
  memcpy(tempf->olddiscbase, forknode->discbase, endsite*sizeof(unsigned char));
  memcpy(tempf->oldnumsteps, forknode->numsteps, endsite*sizeof(long));
  tempf->numdesc = other->numdesc;
  if (forknode->back)
    trydescendants(item, forknode, forknode->back, tempf, false);
  if (!other->tip) {
    memcpy(temp->discbase, other->discbase, endsite*sizeof(unsigned char));
    memcpy(temp->numsteps, other->numsteps, endsite*sizeof(long));
    memcpy(temp->discnumnuc, other->discnumnuc, endsite*sizeof(discnucarray));
    temp->numdesc = other->numdesc + 1;
    multifillin(temp, tempadd, 1);
    if (forknode->back)
      sumnsteps2(tempsum, forknode->back, temp, 0, endsite, threshwt);
    else
      sumnsteps2(tempsum, NULL, temp, 0, endsite, threshwt);
    belowbetter = true;
    if (lastrearr) {
      if (-tempsum->sumsteps >= bstlike2) {
        belowbetter = false;
        bestever = false;
        multf = true;
        if (-tempsum->sumsteps > bstlike2)
          bestever = true;
        savelocrearr(item, forknode, other, tmp, tmp1, tmp2, tmp3, tmprm,
                       tmpadd, &root, maxtrees, &nextree, multf, bestever,
                       &saved, place, bestrees, treenode, &grbg,
                       zeros, zeros2);
        if (saved) {
          like = bstlike2 = -tempsum->sumsteps;
          there = other;
          mulf = true;
        }
      }
    } else if (-tempsum->sumsteps >= like) {
      there = other;
      mulf = true;
      like = -tempsum->sumsteps;
    }
    if (forknode->back) {
      memcpy(temprm->discbase, forknode->back->discbase, endsite*sizeof(unsigned char));
      memcpy(temprm->numsteps, forknode->back->numsteps, endsite*sizeof(long));
    } else {
      memcpy(temprm->discbase, zeros2, endsite*sizeof(unsigned char));
      memcpy(temprm->numsteps, zeros, endsite*sizeof(long));
    }
    memcpy(temprm->olddiscbase, other->back->discbase, endsite*sizeof(unsigned char));
    memcpy(temprm->oldnumsteps, other->back->numsteps, endsite*sizeof(long));
    q = other->next;
    while (q != other) {
      memcpy(temp2->discbase, q->discbase, endsite*sizeof(unsigned char));
      memcpy(temp2->numsteps, q->numsteps, endsite*sizeof(long));
      memcpy(temp2->discnumnuc, q->discnumnuc, endsite*sizeof(discnucarray));
      if (forknode->back) {
        temp2->numdesc = q->numdesc;
        multifillin(temp2, temprm, 0);
      } else {
        temp2->numdesc = q->numdesc - 1;
        multifillin(temp2, temprm, -1);
      }
      if (!q->back->tip)
        trydescendants(item, forknode, q->back, temp2, true);
      else {
        sumnsteps(temp1, q->back, tempadd, 0, endsite);
        sumnsteps2(tempsum, temp1, temp2, 0, endsite, threshwt);
        if (lastrearr) {
          if (-tempsum->sumsteps >= bstlike2) {
            trysave = false;
            multf = false;
            if (belowbetter) {
              bestever = false;
              trysave = true;
            }
            if (-tempsum->sumsteps > bstlike2) {
              bestever = true;
              trysave = true;
            }
            if (trysave) {
              savelocrearr(item, forknode, q->back, tmp, tmp1, tmp2, tmp3,
                           tmprm, tmpadd, &root, maxtrees, &nextree, multf,
                           bestever, &saved, place, bestrees, treenode,
                           &grbg, zeros, zeros2);
              if (saved) {
                like = bstlike2 = -tempsum->sumsteps;
                there = q->back;
                mulf = false;
              }
            }
          }
        } else if (-tempsum->sumsteps > like) {
          like = -tempsum->sumsteps;
          if (-tempsum->sumsteps > bestyet) {
            there = q->back;
            mulf = false;
          }
        }
      }
      q = q->next;
    }
  }
} /* trylocal2 */


void tryrearr(node *p, boolean *success)
{
  /* evaluates one rearrangement of the tree.
     if the new tree has greater "likelihood" than the old
     one sets success = TRUE and keeps the new tree.
     otherwise, restores the old tree */
  node *forknode, *newfork, *other, *oldthere;
  double oldlike;
  boolean oldmulf;

  if (p->back == NULL)
    return;
  forknode = treenode[p->back->index - 1]; 
  if (!forknode->back && forknode->numdesc <= 2 && alltips(forknode, p))
    return;
  oldlike = bestyet;
  like = -10.0 * spp * chars;
  memcpy(tempadd->discbase, p->discbase, endsite*sizeof(unsigned char));
  memcpy(tempadd->numsteps, p->numsteps, endsite*sizeof(long));
  memcpy(tempadd->olddiscbase, zeros2, endsite*sizeof(unsigned char));
  memcpy(tempadd->oldnumsteps, zeros, endsite*sizeof(long));
  if (forknode->numdesc > 2) {
    oldthere = there = forknode;
    oldmulf = mulf = true;
    trylocal(p, forknode);
  } else {
    findbelow(&other, p, forknode);
    oldthere = there = other;
    oldmulf = mulf = false;
    trylocal2(p, forknode, other);
  }
  if ((like <= oldlike) || (there == oldthere && mulf == oldmulf))
    return;
  recompute = true;
  re_move(p, &forknode, &root, recompute, treenode, &grbg, zeros, zeros2);
  if (mulf)
    add(there, p, NULL, &root, recompute, treenode, &grbg, zeros, zeros2);
  else {
    if (forknode->numdesc > 0)
      getnufork(&newfork, &grbg, treenode, zeros, zeros2);
    else
      newfork = forknode;
    add(there, p, newfork, &root, recompute, treenode, &grbg, zeros, zeros2);
  } 
  if (like - oldlike > LIKE_EPSILON) {
    *success = true;
    bestyet = like;
  }
}  /* tryrearr */


void repreorder(node *p, boolean *success)
{
  /* traverses a binary tree, calling PROCEDURE tryrearr
     at a node before calling tryrearr at its descendants */
  node *q, *this;

  if (p == NULL)
    return;
  if (!p->visited) {
    tryrearr(p, success);
    p->visited = true;
  }
  if (!p->tip) {
    q = p;
    while (q->next != p) {
      this = q->next->back;
      repreorder(q->next->back,success);
      if (q->next->back == this)
        q = q->next;
    }
  }
}  /* repreorder */


void rearrange(node **r)
{
  /* traverses the tree (preorder), finding any local
     rearrangement which decreases the number of steps.
     if traversal succeeds in increasing the tree's
     "likelihood", PROCEDURE rearrange runs traversal again */
  boolean success=true;

  while (success) {
    success = false;
    clearvisited(treenode);
    repreorder(*r, &success);
  }
}  /* rearrange */


void describe()
{
  /* prints ancestors, steps and table of numbers of steps in
     each site */

  if (treeprint) {
    fprintf(outfile, "\nrequires a total of %10.3f\n", like / -10.0);
    fprintf(outfile, "\n  between      and       length\n");
    fprintf(outfile, "  -------      ---       ------\n");
    printbranchlengths(root);
  }
  if (stepbox)
    writesteps(chars, weights, oldweight, root);
  if (ancseq) {
    hypstates(chars, root, treenode, &garbage);
    putc('\n', outfile);
  }
  putc('\n', outfile);
  if (trout) {
    col = 0;
    treeout3(root, nextree, &col, root);
  }
}  /* describe */


void pars_coordinates(node *p, double lengthsum, long *tipy,
                double *tipmax)
{
  /* establishes coordinates of nodes */
  node *q, *first, *last;
  double xx;

  if (p == NULL)
    return;
  if (p->tip) {
    p->xcoord = (long)(over * lengthsum + 0.5);
    p->ycoord = (*tipy);
    p->ymin = (*tipy);
    p->ymax = (*tipy);
    (*tipy) += down;
    if (lengthsum > (*tipmax))
      (*tipmax) = lengthsum;
    return;
  }
  q = p->next;
  do {
    xx = q->v;
    if (xx > 100.0)
      xx = 100.0;
    pars_coordinates(q->back, lengthsum + xx, tipy,tipmax);
    q = q->next;
  } while (p != q);
  first = p->next->back;
  q = p;
  while (q->next != p)
    q = q->next;
  last = q->back;
  p->xcoord = (long)(over * lengthsum + 0.5);
  if ((p == root) || count_sibs(p) > 2)
    p->ycoord = p->next->next->back->ycoord;
  else
    p->ycoord = (first->ycoord + last->ycoord) / 2;
  p->ymin = first->ymin;
  p->ymax = last->ymax;
}  /* pars_coordinates */


void pars_printree()
{
  /* prints out diagram of the tree2 */
  long tipy;
  double scale, tipmax;
  long i;
 
  if (!treeprint)
    return;
  putc('\n', outfile);
  tipy = 1;
  tipmax = 0.0;
  pars_coordinates(root, 0.0, &tipy, &tipmax);
  scale = 1.0 / (long)(tipmax + 1.000);
  for (i = 1; i <= (tipy - down); i++)
    drawline3(i, scale, root);
  putc('\n', outfile);
}  /* pars_printree */


void globrearrange()
{
  /* does global rearrangements */
  long j;
  double gotlike;
  boolean frommulti;
  node *item, *nufork;

  recompute = true;
  do {
    printf("   ");
    gotlike = bestlike;
    for (j = 0; j < nonodes; j++) {
      bestyet = -10.0 * spp * chars;
      if (j < spp)
        item = treenode[enterorder[j] -1];
      else 
        item = treenode[j];
      if ((item != root) && 
           ((j < spp) || ((j >= spp) && (item->numdesc > 0))) &&
           !((item->back->index == root->index) && (root->numdesc == 2)
              && alltips(root, item))) {
        re_move(item, &nufork, &root, recompute, treenode, &grbg, zeros, zeros2);
        frommulti = (nufork->numdesc > 0);
        clearcollapse(treenode);
        there = root;
        memcpy(tempadd->discbase, item->discbase, endsite*sizeof(unsigned char));
        memcpy(tempadd->numsteps, item->numsteps, endsite*sizeof(long));
        memcpy(tempadd->olddiscbase, zeros2, endsite*sizeof(unsigned char));
        memcpy(tempadd->oldnumsteps, zeros, endsite*sizeof(long));
        if (frommulti){
          oldnufork = nufork;
          getnufork(&nufork, &grbg, treenode, zeros, zeros2);
        }
        addpreorder(root, item, nufork);
        if (frommulti)
          oldnufork = NULL;
        if (!mulf)
          add(there, item, nufork, &root, recompute, treenode, &grbg,
                zeros, zeros2);
        else
          add(there, item, NULL, &root, recompute, treenode, &grbg,
                zeros, zeros2);
      }
      if (progress) {
        if (j % ((nonodes / 72) + 1) == 0)
          putchar('.');
        fflush(stdout);
      }
    }
    if (progress)
      putchar('\n');
  } while (bestlike > gotlike);
} /* globrearrange */


void load_tree(long treei)
{ /* restores a tree from bestrees */
  long j, nextnode;
  boolean recompute = false;
  node *dummy;

  for (j = spp - 1; j >= 1; j--)
    re_move(treenode[j], &dummy, &root, recompute, treenode, &grbg,
              zeros, zeros2);
  root = treenode[0];
  recompute = true;
  add(treenode[0], treenode[1], treenode[spp], &root, recompute,
    treenode, &grbg, zeros, zeros2);
  nextnode = spp + 2;
  for (j = 3; j <= spp; j++) {
    if (bestrees[treei].btree[j - 1] > 0)
      add(treenode[bestrees[treei].btree[j - 1] - 1], treenode[j - 1],
            treenode[nextnode++ - 1], &root, recompute, treenode, &grbg,
            zeros, zeros2);
    else
      add(treenode[treenode[-bestrees[treei].btree[j-1]-1]->back->index-1],
            treenode[j - 1], NULL, &root, recompute, treenode, &grbg,
            zeros, zeros2);
  }
} /* load_tree */


void grandrearr()
{
  /* calls either global rearrangement or local rearrangement on best trees */
  long treei;
  boolean done;

  done = false;
  do {
    treei = findunrearranged(bestrees, nextree, true);
    if (treei < 0)
      done = true;
    else 
        bestrees[treei].gloreange = true;
    
    if (!done) {
      load_tree(treei);
      globrearrange();
      done = rearrfirst;
    }
  } while (!done);
} /* grandrearr */


void maketree()
{
  /* constructs a binary tree from the pointers in treenode.
     adds each node at location which yields highest "likelihood"
  then rearranges the tree for greatest "likelihood" */
  long i, j, numtrees, nextnode;
  boolean done, firsttree, goteof, haslengths;
  node *item, *nufork, *dummy;
  pointarray nodep;

  if (!usertree) {
    for (i = 1; i <= spp; i++)
      enterorder[i - 1] = i;
    if (jumble)
      randumize(seed, enterorder);
    recompute = true;
    root = treenode[enterorder[0] - 1];
    add(treenode[enterorder[0] - 1], treenode[enterorder[1] - 1],
        treenode[spp], &root, recompute, treenode, &grbg, zeros, zeros2);
    if (progress) {
      printf("Adding species:\n");
      writename(0, 2, enterorder);
#ifdef WIN32
      phyFillScreenColor();
#endif
    }
    lastrearr = false;
    oldnufork = NULL;
    for (i = 3; i <= spp; i++) {
      bestyet = -10.0 * spp * chars;
      item = treenode[enterorder[i - 1] - 1];
      getnufork(&nufork, &grbg, treenode, zeros, zeros2);
      there = root;
      memcpy(tempadd->discbase, item->discbase, endsite*sizeof(unsigned char));
      memcpy(tempadd->numsteps, item->numsteps, endsite*sizeof(long));
      memcpy(tempadd->olddiscbase, zeros2, endsite*sizeof(unsigned char));
      memcpy(tempadd->oldnumsteps, zeros, endsite*sizeof(long));
      addpreorder(root, item, nufork);
      if (!mulf)
        add(there, item, nufork, &root, recompute, treenode, &grbg,
              zeros, zeros2);
      else
        add(there, item, NULL, &root, recompute, treenode, &grbg,
              zeros, zeros2);
      like = bestyet;
      rearrange(&root);
      if (progress) {
        writename(i - 1, 1, enterorder);
#ifdef WIN32
        phyFillScreenColor();
#endif
      }
      lastrearr = (i == spp);
      if (lastrearr) {
        bestlike = bestyet;
        if (jumb == 1) {
          bstlike2 = bestlike;
          nextree = 1;
          initbestrees(bestrees, maxtrees, true);
          initbestrees(bestrees, maxtrees, false);
        }
        if (progress) {
          printf("\nDoing global rearrangements");
          if (rearrfirst)
            printf(" on the first of the trees tied for best\n");
          else
            printf(" on all trees tied for best\n");
          printf("  !");
          for (j = 0; j < nonodes; j++)
            if (j % ((nonodes / 72) + 1) == 0)
              putchar('-');
          printf("!\n");
#ifdef WIN32
          phyFillScreenColor();
#endif
    }
        globrearrange();
        rearrange(&root);
      }
    }
    done = false;
    while (!done && findunrearranged(bestrees, nextree, true) >= 0) {
      grandrearr();
      done = rearrfirst;
    }
    if (progress) {
      putchar('\n');
#ifdef WIN32
      phyFillScreenColor();
#endif
    }
    recompute = false;
    for (i = spp - 1; i >= 1; i--)
      re_move(treenode[i], &dummy, &root, recompute, treenode, &grbg,
                zeros, zeros2);
    if (jumb == njumble) {
      collapsebestrees(&root, &grbg, treenode, bestrees, place, zeros,
                       zeros2, chars, recompute, progress);
      if (treeprint) {
        putc('\n', outfile);
        if (nextree == 2)
          fprintf(outfile, "One most parsimonious tree found:\n");
        else
          fprintf(outfile, "%6ld trees in all found\n", nextree - 1);
      }
      if (nextree > maxtrees + 1) {
        if (treeprint)
          fprintf(outfile, "here are the first %4ld of them\n", (long)maxtrees);
        nextree = maxtrees + 1;
      }
      if (treeprint)
        putc('\n', outfile);
      for (i = 0; i <= (nextree - 2); i++) {
        root = treenode[0];
        add(treenode[0], treenode[1], treenode[spp], &root, recompute,
                   treenode, &grbg, zeros, zeros2);
        nextnode = spp + 2;
        for (j = 3; j <= spp; j++) {
          if (bestrees[i].btree[j - 1] > 0)
            add(treenode[bestrees[i].btree[j - 1] - 1], treenode[j - 1],
              treenode[nextnode++ - 1], &root, recompute, treenode, &grbg,
              zeros, zeros2);
          else
            add(treenode[treenode[-bestrees[i].btree[j - 1]-1]->back->index-1],
              treenode[j - 1], NULL, &root, recompute, treenode, &grbg,
              zeros, zeros2);
        }
        reroot(treenode[outgrno - 1], root);
        postorder(root);
        evaluate(root);
        treelength(root, chars, treenode);
        pars_printree();
        describe();
        for (j = 1; j < spp; j++)
          re_move(treenode[j], &dummy, &root, recompute, treenode,
                    &grbg, zeros, zeros2);
      }
    }
  } else {
    /* Open in binary: ftell() is broken for UNIX line-endings under WIN32 */
    openfile(&intree,INTREE,"input tree", "rb",progname,intreename);
    numtrees = countsemic(&intree);
    if (numtrees > MAXNUMTREES) {
      printf(
            "\n\nERROR: number of input trees is read incorrectly from %s\n\n",
             intreename);
      exxit(-1);
    }
    if (numtrees > 2)
      initseed(&inseed, &inseed0, seed);
    if (treeprint) {
      fprintf(outfile, "User-defined tree");
      if (numtrees > 1)
        putc('s', outfile);
      fprintf(outfile, ":\n");
    }
    fsteps = (long **)Malloc(maxuser*sizeof(long *));
    for (j = 1; j <= maxuser; j++)
      fsteps[j - 1] = (long *)Malloc(endsite*sizeof(long));
    nodep = NULL;
    which = 1;
    while (which <= numtrees) {
      firsttree = true;
      nextnode = 0;
      haslengths = true;
      treeread(intree, &root, treenode, &goteof, &firsttree, nodep, &nextnode,
                &haslengths, &grbg, initparsnode,false,nonodes);
      if (treeprint)
        fprintf(outfile, "\n\n");
      if (outgropt)
        reroot(treenode[outgrno - 1], root);
      postorder(root);
      evaluate(root);
      treelength(root, chars, treenode);
      pars_printree();
      describe();
      if (which < numtrees)
        gdispose(root, &grbg, treenode);
      which++;
    }
    FClose(intree);
    putc('\n', outfile);
    if (numtrees > 1 && chars > 1 )
      standev(chars, numtrees, minwhich, minsteps, nsteps, fsteps, seed);
    for (j = 1; j <= maxuser; j++)
      free(fsteps[j - 1]);
    free(fsteps);
  }
  if (jumb == njumble) {
    if (progress) {
      printf("\nOutput written to file \"%s\"\n", outfilename);
      if (trout) {
        printf("\nTree");
        if ((usertree && numtrees > 1) || (!usertree && nextree != 2))
          printf("s");
        printf(" also written onto file \"%s\"\n", outtreename);
      }
    }
  }
}  /* maketree */


void freerest()
{
  if (!usertree) {
    freenode(&temp);
    freenode(&temp1);
    freenode(&temp2);
    freenode(&tempsum);
    freenode(&temprm);
    freenode(&tempadd);
    freenode(&tempf);
    freenode(&tmp);
    freenode(&tmp1);
    freenode(&tmp2);
    freenode(&tmp3);
    freenode(&tmprm);
    freenode(&tmpadd);
  }
  freegrbg(&grbg);
  if (ancseq)
    freegarbage(&garbage);
  free(threshwt);
  free(zeros);
  free(zeros2);
  freenodes(nonodes, treenode);
}  /* freerest*/


int main(int argc, Char *argv[])
{  /* Discrete character parsimony by uphill search */

  /* reads in spp, chars, and the data. Then calls maketree to
     construct the tree */
#ifdef MAC
   argc = 1;                /* macsetup("Pars","");                */
   argv[0]="Pars";
#endif
  init(argc, argv);
  progname = argv[0];
  openfile(&infile,INFILE,"input file", "r",argv[0],infilename);
  openfile(&outfile,OUTFILE,"output file", "w",argv[0],outfilename);

  ibmpc = IBMCRT;
  ansi = ANSICRT;
  msets = 1;
  firstset = true;
  garbage = NULL;
  grbg = NULL;
  doinit();
  if (weights || justwts)
    openfile(&weightfile,WEIGHTFILE,"weights file","r",argv[0],weightfilename);
  if (trout)
    openfile(&outtree,OUTTREE,"output tree file", "w",argv[0],outtreename);
  for (ith = 1; ith <= msets; ith++) {
    if (msets > 1 && !justwts) {
      fprintf(outfile, "\nData set # %ld:\n\n", ith);
      if (progress)
        printf("\nData set # %ld:\n\n", ith);
    }
    doinput();
    if (ith == 1)
      firstset = false;
    for (jumb = 1; jumb <= njumble; jumb++)
      maketree();
  freerest();    
  }
  FClose(infile);
  FClose(outfile);
  if (weights || justwts)
    FClose(weightfile);
  if (trout)
    FClose(outtree);
  if (usertree)
    FClose(intree);
#ifdef MAC
  fixmacfile(outfilename);
  fixmacfile(outtreename);
#endif
  if (progress)
    printf("\nDone.\n\n");
#ifdef WIN32
  phyRestoreConsoleAttributes();
#endif
  return 0;
}  /* Discrete character parsimony by uphill search */

