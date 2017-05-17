#include "phylip.h"
#include "cons.h"

/* version 3.6. (c) Copyright 1993-2008 by the University of Washington.
   Written by Joseph Felsenstein, Hisashi Horino,
   Akiko Fuseki, Dan Fineman, Sean Lamont, and Andrew Keeffe.
   Permission is granted
   to copy and use this program provided no fee is charged for it and
   provided that this copyright notice is not removed. */


/* The following extern's refer to things declared in cons.c */

extern int tree_pairing;

extern Char outfilename[FNMLNGTH], intreename[FNMLNGTH], intree2name[FNMLNGTH], outtreename[FNMLNGTH];
extern node *root;

extern long numopts, outgrno, col;
extern long maxgrp;               /* max. no. of groups in all trees found  */

extern boolean trout, firsttree, noroot, outgropt, didreroot, prntsets,
          progress, treeprint, goteof, strict, mr, mre, ml;
extern pointarray nodep;                 /* pointers to all nodes in tree */
extern group_type **grouping, **grping2, **group2;/* to store groups found  */
extern long **order, **order2, lasti;
extern group_type *fullset;
extern long tipy;

extern double trweight, ntrees, mlfrac;

#ifndef OLDC
/* function prototypes */
void   getoptions(void);
void   count_siblings(node **p);
void   treeout(node *);
/* function prototypes */
#endif


void getoptions()
{
  /* interactively set options */
  long loopcount, loopcount2;
  Char ch;
  boolean done, done1;

  /* Initial settings */
  ibmpc          = IBMCRT;
  ansi           = ANSICRT;
  didreroot      = false;
  firsttree      = true;
  spp            = 0 ;
  col            = 0 ;

  /* This is needed so functions in cons.c work */
  tree_pairing   = NO_PAIRING ;  

  fprintf(outfile, "\nConsensus tree");
  fprintf(outfile, " program, version %s\n\n", VERSION);
  putchar('\n');
  strict = false;
  mr = false;
  mre = true;
  ml = false;
  mlfrac = 0.5;
  noroot = true;
  numopts = 0;
  outgrno = 1;
  outgropt = false;
  trout = true;
  prntsets = true;
  progress = true;
  treeprint = true;
  loopcount = 0;
  do {
    cleerhome();
    printf("\nConsensus tree");
    printf(" program, version %s\n\n", VERSION);
    printf("Settings for this run:\n");
    printf(" C         Consensus type (MRe, strict, MR, Ml):");
    if (strict)
      printf("  strict\n");
    else if (mr)
        printf("  Majority rule\n");
      else if (mre)
          printf("  Majority rule (extended)\n");
        else if (ml)
            printf("  Ml\n");
          else printf("  Adams\n");
    if (noroot) {
      printf(" O                                Outgroup root:");
      if (outgropt)
        printf("  Yes, at species number%3ld\n", outgrno);
      else
        printf("  No, use as outgroup species%3ld\n", outgrno);
      }
    printf(" R                Trees to be treated as Rooted:");
    if (noroot)
      printf("  No\n");
    else
      printf("  Yes\n");
    printf(" T           Terminal type (IBM PC, ANSI, none):");
    if (ibmpc)
      printf("  IBM PC\n");
    if (ansi)
      printf("  ANSI\n");
    if (!(ibmpc || ansi))
      printf("  (none)\n");
    printf(" 1                Print out the sets of species:");
    if (prntsets)
      printf("  Yes\n");
    else
      printf("  No\n");
    printf(" 2         Print indications of progress of run:  %s\n",
           (progress ? "Yes" : "No"));
    printf(" 3                               Print out tree:");
    if (treeprint)
      printf("  Yes\n");
    else
      printf("  No\n");
    printf(" 4               Write out trees onto tree file:");
    if (trout)
      printf("  Yes\n");
    else
      printf("  No\n");

    printf("\nAre these settings correct? (type Y or the letter for one to change)\n");
#ifdef WIN32
    phyFillScreenColor();
#endif
    fflush(stdout);
    scanf("%c%*[^\n]", &ch);
    getchar();
    uppercase(&ch);
    done = (ch == 'Y');
    if (!done) {
      if ((noroot && (ch == 'O')) || strchr("CRT1234",ch) != NULL) {
        switch (ch) {

        case 'C':
          if (strict) {
            strict = false;
            mr = true;
          } else {
            if (ml) {
              ml = false;
              mre = true;
            } else {
              if (mre) {
                mre = false;
                strict = true;
              } else {
                if (mr) {
                  mr = false;
                  ml = true;
                }
              }
            }
          }
          break;

        case 'O':
          outgropt = !outgropt;
          if (outgropt) {
            numopts++;
            loopcount2 = 0;
            do {
              printf("Type number of the outgroup:\n");
#ifdef WIN32
              phyFillScreenColor();
#endif
              fflush(stdout);
              scanf("%ld%*[^\n]", &outgrno);
              getchar();
              done1 = (outgrno >= 1);
              if (!done1) {
                printf("ERROR: Bad outgroup number: %ld\n", outgrno);
                printf("  Must be greater than zero\n");
              }
            countup(&loopcount2, 10);
            } while (done1 != true);
          }
          break;

        case 'R':
          noroot = !noroot;
          break;

        case 'T':
          initterminal(&ibmpc, &ansi);
          break;

        case '1':
          prntsets = !prntsets;
          break;

        case '2':
          progress = !progress;
          break;
        
        case '3':
          treeprint = !treeprint;
          break;

        case '4':
          trout = !trout;
          break;

        }
      } else
        printf("Not a possible option!\n");
    }
    countup(&loopcount, 100);
  } while (!done);
  if (ml) {
    do {
      printf("\nFraction (l) of times a branch must appear\n");
      fflush(stdout);
      scanf("%lf%*[^\n]", &mlfrac);
      getchar();
    } while ((mlfrac < 0.5) || (mlfrac > 1.0));
  }
}  /* getoptions */


void count_siblings(node **p)
{
  node *tmp_node;
  int i;

  if (!(*p)) {
    /* This is a leaf, */
    return;
  } else {
    tmp_node = (*p)->next;
  }

  for (i = 0 ; i < 1000; i++) {
    if (tmp_node == (*p)) {
      /* When we've gone through all the siblings, */
      break;
    } else if (tmp_node) {
      tmp_node = tmp_node->next;
    } else  {
      /* Should this be executed? */
      return ;
    }
  }
} /* count_siblings */


void treeout(node *p)
{
  /* write out file with representation of final tree */
  long i, n = 0;
  Char c;
  node *q;
  double x;

  count_siblings (&p);  

  if (p->tip) {
    /* If we're at a node which is a leaf, figure out how long the
       name is and print it out. */
    for (i = 1; i <= MAXNCH; i++) {
      if (p->nayme[i - 1] != '\0')
        n = i;
    }
    for (i = 0; i < n; i++) {
      c = p->nayme[i];
      if (c == ' ')
        c = '_';
      putc(c, outtree);
    }
    col += n;
  } else {
    /* If we're at a furcation, print out the proper formatting, loop
       through all the children, calling the procedure recursively. */
    putc('(', outtree);
    col++;
    q = p->next;
    while (q != p) {
      /* This should terminate when we've gone through all the
         siblings, */
      treeout(q->back);
      q = q->next;
      if (q == p)
        break;
      putc(',', outtree);
      col++;
      if (col > 60) {
        putc('\n', outtree);
        col = 0;
      }
    }
    putc(')', outtree);
    col++;
  }

  if (p->tip)
    x = ntrees;
  else
    x = (double)p->deltav;

  if (p == root) {
    /* When we're all done with this tree, */
    fprintf(outtree, ";\n");
    return;
  }

  /* Figure out how many characters the branch length requires: */
  else {
    if (!strict) {
      if (x >= 100.0) { 
        fprintf(outtree, ":%5.1f", x);
        col += 4;
      } else if (x >= 10.0) {
          fprintf(outtree, ":%4.1f", x); 
          col += 3;
        } else if (x >= 1.00) {
            fprintf(outtree, ":%4.2f", x); 
            col += 3;
          }
    }
  }
}  /* treeout */


int main(int argc, Char *argv[])
{  
  /* Local variables added by Dan F. */
  pattern_elm  ***pattern_array;
  long trees_in = 0;
  long i, j;
  long tip_count = 0;
  node *p, *q;

#ifdef MAC
  argc = 1;                /* macsetup("Consense", "");        */
  argv[0] = "Consense";
#endif
  init(argc, argv);

  /* Open in binary: ftell() is broken for UNIX line-endings under WIN32 */
  openfile(&intree, INTREE, "input tree file", "rb", argv[0], intreename);
  openfile(&outfile, OUTFILE, "output file", "w", argv[0], outfilename);

  /* Initialize option-based variables, then ask for changes regarding
     their values. */
  getoptions();

  ntrees = 0.0;
  maxgrp = 32767;   /* initial size of set hash table */
  lasti  = -1;

  if (trout) 
    openfile(&outtree, OUTTREE, "output tree file", "w", argv[0], outtreename);
  if (prntsets)
    fprintf(outfile, "Species in order: \n\n");

  trees_in = countsemic(&intree);
  countcomma(&intree,&tip_count);
  tip_count++; /* countcomma does a raw comma count, tips is one greater */



  /* Read the tree file and put together grouping, order, and timesseen */
  read_groups (&pattern_array, trees_in, tip_count, intree);
  /* Compute the consensus tree. */
  putc('\n', outfile);
  nodep      = (pointarray)Malloc(2*(1+spp)*sizeof(node *));
  for (i = 0; i < spp; i++) {
    nodep[i] = (node *)Malloc(sizeof(node));
    for (j = 0; j < MAXNCH; j++)
      nodep[i]->nayme[j] = '\0';
    strncpy(nodep[i]->nayme, nayme[i], MAXNCH);
  }
  for (i = spp; i < 2*(1+spp); i++)
    nodep[i] = NULL;
  consensus(pattern_array, trees_in);
  printf("\n");
  if (trout) {
    treeout(root);
    if (progress)
      printf("Consensus tree written to file \"%s\"\n\n", outtreename);
  }
  if (progress)
    printf("Output written to file \"%s\"\n\n", outfilename);
  for (i = 0; i < spp; i++)
    free(nodep[i]);
  for (i = spp; i < 2*(1 + spp); i++) {
    if (nodep[i] != NULL) {
      p = nodep[i]->next;
      do {
        q = p->next;
        free(p);
        p = q;
      } while (p != nodep[i]);
      free(p);
    }
  }
  free(nodep);
  FClose(outtree);
  FClose(intree);
  FClose(outfile);

#ifdef MAC
  fixmacfile(outfilename);
  fixmacfile(outtreename);
#endif
printf("Done.\n\n");

#ifdef WIN32
  phyRestoreConsoleAttributes();
#endif

return 0;
}  /* main */

