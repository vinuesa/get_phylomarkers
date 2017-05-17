#include "phylip.h"
#include "cons.h"

int tree_pairing;

Char outfilename[FNMLNGTH], intreename[FNMLNGTH], intree2name[FNMLNGTH], outtreename[FNMLNGTH];
node *root;

long numopts, outgrno, col, setsz;
long maxgrp;               /* max. no. of groups in all trees found  */

boolean trout, firsttree, noroot, outgropt, didreroot, prntsets,
               progress, treeprint, goteof, strict, mr=false, mre=false,
               ml=false; /* initialized all false for Treedist */
pointarray nodep;
pointarray treenode;
group_type **grouping, **grping2, **group2;/* to store groups found  */
double *lengths, *lengths2;
long **order, **order2, lasti;
group_type *fullset;
node *grbg;
long tipy;

double **timesseen, **tmseen2, **times2 ;
double *tchange2;
double trweight, ntrees, mlfrac;

/* prototypes */
void censor(void);
boolean compatible(long, long);
void elimboth(long);
void enterpartition (group_type*, long*);
void reorient(node* n);

/* begin hash table code */

#define NUM_BUCKETS 100

typedef struct namenode {
  struct namenode *next;
  plotstring naym;
  int hitCount;
} namenode;

typedef namenode **hashtype;

hashtype hashp;

long namesGetBucket(plotstring);
void namesAdd(plotstring);
boolean namesSearch(plotstring);
void namesDelete(plotstring);
void namesClearTable(void);
void namesCheckTable(void);
void missingnameRecurs(node *p);

/**
 * namesGetBucket - return the bucket for a given name
 */
long namesGetBucket(plotstring searchname) {
  long i;
  long sum = 0;

  for (i = 0; (i < MAXNCH) && (searchname[i] != '\0'); i++) {
    sum += searchname[i];
  }
  return (sum % NUM_BUCKETS);
}


/**
 * namesAdd - add a name to the hash table
 *
 * The argument is added at the head of the appropriate linked list.  No 
 * checking is done for duplicates.  The caller can call 
 * namesSearch to check for an existing name prior to calling
 * namesAdd.
 */
void namesAdd(plotstring addname) {
  long bucket = namesGetBucket(addname);
  namenode *hp, *temp;

  temp = hashp[bucket];
  hashp[bucket] = (namenode *)Malloc(sizeof(namenode));
  hp = hashp[bucket];
  strcpy(hp->naym, addname);
  hp->next = temp;
  hp->hitCount = 0;
}

/**
 * namesSearch - search for a name in the hash table
 *
 * Return true if the name is found, else false.
 */
boolean namesSearch(plotstring searchname) {
  long i = namesGetBucket(searchname);
  namenode *p;

  p = hashp[i];
  if (p == NULL) {
    return false;
  }
  do {
    if (strcmp(searchname,p->naym) == 0) {
      p->hitCount++;
      return true;
    }
    p = p->next;
  } while (p != NULL);

  return false;  
}

/**
 * Go through hash table and check that the hit count on all entries is one.
 * If it is zero, then a species was missed, if it is two, then there is a
 * duplicate species.
 */

void namesCheckTable(void) {
  namenode *p;  
  long i;

  for (i=0; i< NUM_BUCKETS; i++) {
    p = hashp[i];
    while (p != NULL){
      if(p->hitCount >1){
        printf("\n\nERROR in user tree: duplicate name found: ");
        puts(p->naym);
        printf("\n\n");
        exxit(-1);
      } else if(p->hitCount == 0){
        printf("\n\nERROR in user tree: name %s not found\n\n\n",
               p->naym);
        exxit(-1);
      }
      p->hitCount = 0;
      p = p->next;
    } 
  }  
}

/**
 * namesClearTable - empty names out of the table and 
 *                   return allocated memory
 */
void namesClearTable(void) {
  long i;
  namenode *p, *temp;

  for (i=0; i< NUM_BUCKETS; i++) {
    p = hashp[i];
    if (p != NULL) {
      do {
        temp = p;
        p = p->next;
        free(temp);
      } while (p != NULL);
    hashp[i] = NULL;
    }
  }
}
/* end hash table code */

void initconsnode(node **p, node **grbg, node *q, long len, long nodei,
                        long *ntips, long *parens, initops whichinit,
                        pointarray treenode, pointarray nodep, Char *str,
                        Char *ch, FILE *intree)
{
  /* initializes a node */
  long i;
  char c;
  boolean minusread;
  double valyew, divisor, fracchange;

  switch (whichinit) {
  case bottom:
    gnu(grbg, p);
    (*p)->index = nodei;
    (*p)->tip = false;
    for (i=0; i<MAXNCH; i++)
      (*p)->nayme[i] = '\0';
    nodep[(*p)->index - 1] = (*p);
    (*p)->v = 0;
    break;
  case nonbottom:
    gnu(grbg, p);
    (*p)->index = nodei;
    (*p)->v = 0;
    break;
  case tip:
    (*ntips)++;
    gnu(grbg, p);
    nodep[(*ntips) - 1] = *p;
    setupnode(*p, *ntips);
    (*p)->tip = true;
    strncpy ((*p)->nayme, str, MAXNCH);
    if (firsttree && prntsets) {
      fprintf(outfile, "  %ld. ", *ntips);
      for (i = 0; i < len; i++)
        putc(str[i], outfile);
      putc('\n', outfile);
      if ((*ntips > 0) && (((*ntips) % 10) == 0))
        putc('\n', outfile);
    }
    (*p)->v = 0;
    break;
  case length:
    processlength(&valyew, &divisor, ch, &minusread, intree, parens);
    fracchange = 1.0;
    (*p)->v = valyew / divisor / fracchange;
    break;
  case treewt:
    if (!eoln(intree)) {
      if (fscanf(intree, "%lf", &trweight) == 1) {
        getch(ch, parens, intree);
        if (*ch != ']') {
          printf("\n\nERROR: Missing right square bracket\n\n");
          exxit(-1);
        } else {
          getch(ch, parens, intree);
          if (*ch != ';') {
            printf("\n\nERROR: Missing semicolon after square brackets\n\n");
            exxit(-1);
          }
        }
      }
      else {
        printf("\n\nERROR: Expecting tree weight in last comment field\n\n");
        exxit(-1);
      }
    }
    break;
  case unittrwt:
    /* This comes not only when setting trweight but also at the end of
     * any tree. The following code saves the current position in a 
     * file and reads to a new line. If there is a new line then we're 
     * at the end of tree, otherwise warn the user. This function should
     * really leave the file alone, so once we're done with 'intree' 
     * we seek the position back so that it doesn't look like we did 
     * anything */
    trweight = 1.0 ;
    i = ftell (intree);
    c = ' ';
    while (c == ' ')  {
      if (eoff(intree)) {
        fseek(intree,i,SEEK_SET);
        return;
      }
      c = gettc(intree);
    }
    fseek(intree,i,SEEK_SET);
    if ( c != '\n' && c!= '\r')
      printf("WARNING: Tree weight set to 1.0\n");
    if ( c == '\r' )
      if ( (c == gettc(intree)) != '\n')
        ungetc(c, intree);
    break;
  case hsnolength:
    (*p)->v = -1;         /* signal value that a length is missing */
    break;
  default:                /* cases hslength, iter, hsnolength      */ 
    break;                /* should there be an error message here?*/
  }
} /* initconsnode */


void censor(void)
{
  /* delete groups that are too rare to be in the consensus tree */
  long i;

  i = 1;
  do {
    if (timesseen[i-1])
      if (!(mre || (mr && (2*(*timesseen[i-1]) > ntrees))
                || (ml && ((*timesseen[i-1]) > mlfrac*ntrees))
                || (strict && ((*timesseen[i-1]) == ntrees)))) {
        free(grouping[i - 1]);
        free(timesseen[i - 1]);
        grouping[i - 1] = NULL;
        timesseen[i - 1] = NULL;
    }
    i++;
  } while (i < maxgrp);
} /* censor */


void compress(long *n)
{
  /* push all the nonempty subsets to the front end of their array */
  long i, j;

  i = 1;
  j = 1;
  do {
    while (grouping[i - 1] != NULL)
      i++;
    if (j <= i)
      j = i + 1;
    while ((grouping[j - 1] == NULL) && (j < maxgrp))
      j++;
    if (j < maxgrp) {
      grouping[i - 1] = (group_type *)Malloc(setsz * sizeof(group_type));
      timesseen[i - 1] = (double *)Malloc(sizeof(double));
      memcpy(grouping[i - 1], grouping[j - 1], setsz * sizeof(group_type));
      *timesseen[i - 1] = *timesseen[j - 1];
      free(grouping[j - 1]);
      free(timesseen[j - 1]);
      grouping[j - 1] = NULL;
      timesseen[j - 1] = NULL;
    }
  } while (j != maxgrp);
  (*n) = i - 1;
}  /* compress */


void sort(long n)
{
  /* Shell sort keeping grouping, timesseen in same order */
  long gap, i, j;
  group_type *stemp;
  double rtemp;

  gap = n / 2;
  stemp = (group_type *)Malloc(setsz * sizeof(group_type));
  while (gap > 0) {
    for (i = gap + 1; i <= n; i++) {
      j = i - gap;
      while (j > 0) {
        if (*timesseen[j - 1] < *timesseen[j + gap - 1]) {
          memcpy(stemp, grouping[j - 1], setsz * sizeof(group_type));
          memcpy(grouping[j - 1], grouping[j + gap - 1], setsz * sizeof(group_type));
          memcpy(grouping[j + gap - 1], stemp, setsz * sizeof(group_type));
          rtemp = *timesseen[j - 1];
          *timesseen[j - 1] = *timesseen[j + gap - 1];
          *timesseen[j + gap - 1] = rtemp;
        }
        j -= gap;
      }
    }
    gap /= 2;
  }
  free(stemp);
}  /* sort */


boolean compatible(long i, long j)
{
  /* are groups i and j compatible? */
  boolean comp;
  long k;

  comp = true;
  for (k = 0; k < setsz; k++)
    if ((grouping[i][k] & grouping[j][k]) != 0)
      comp = false;
  if (!comp) {
    comp = true;
    for (k = 0; k < setsz; k++)
      if ((grouping[i][k] & ~grouping[j][k]) != 0)
        comp = false;
    if (!comp) {
      comp = true;
      for (k = 0; k < setsz; k++)
        if ((grouping[j][k] & ~grouping[i][k]) != 0)
          comp = false;
      if (!comp) {
        comp = noroot;
        if (comp) {
          for (k = 0; k < setsz; k++)
            if ((fullset[k] & ~grouping[i][k] & ~grouping[j][k]) != 0)
              comp = false;
        }
      }
    }
  }
  return comp;
} /* compatible */


void eliminate(long *n, long *n2)
{
  /* eliminate groups incompatible with preceding ones */
  long i, j, k;
  boolean comp;

  for (i = 2; i <= (*n); i++) {
    comp = true;
    for (j = 0; comp && (j <= i - 2); j++) {
      if ((timesseen[j] != NULL) && *timesseen[j] > 0) {
        comp = compatible(i-1,j);
        if (!comp) {
          (*n2)++;
          times2[(*n2) - 1] = (double *)Malloc(sizeof(double));
          group2[(*n2) - 1] = (group_type *)Malloc(setsz * sizeof(group_type));
          *times2[(*n2) - 1] = *timesseen[i - 1];
          memcpy(group2[(*n2) - 1], grouping[i - 1], setsz * sizeof(group_type));
          *timesseen[i - 1] = 0.0;
          for (k = 0; k < setsz; k++)
            grouping[i - 1][k] = 0;
        }
      }
    }
    if (*timesseen[i - 1] == 0.0) {
      free(grouping[i - 1]);
      free(timesseen[i -  1]);
      timesseen[i - 1] = NULL;
      grouping[i - 1] = NULL;
    }
  }
}  /* eliminate */


void printset(long n)
{
  /* print out the n sets of species */
  long i, j, k, size;
  boolean noneprinted;

  fprintf(outfile, "\nSet (species in order)   ");
  for (i = 1; i <= spp - 25; i++)
    putc(' ', outfile);
  fprintf(outfile, "  How many times out of %7.2f\n\n", ntrees);
  noneprinted = true;
  for (i = 0; i < n; i++) {
    if ((timesseen[i] != NULL) && (*timesseen[i] > 0)) {
      size = 0;
      k = 0;
      for (j = 1; j <= spp; j++) {
        if (j == ((k+1)*SETBITS+1)) k++;
        if (((1L << (j - 1 - k*SETBITS)) & grouping[i][k]) != 0)
          size++;
      }
      if (size != 1 && !(noroot && size >= (spp-1))) {
        noneprinted = false;
        k = 0;
        for (j = 1; j <= spp; j++) {
          if (j == ((k+1)*SETBITS+1)) k++;
          if (((1L << (j - 1 - k*SETBITS)) & grouping[i][k]) != 0)
            putc('*', outfile);
          else
            putc('.', outfile);
          if (j % 10 == 0)
            putc(' ', outfile);
        }
        for (j = 1; j <= 23 - spp; j++)
          putc(' ', outfile);
        fprintf(outfile, "    %5.2f\n", *timesseen[i]);
      }
    }
  }
  if (noneprinted)
    fprintf(outfile, " NONE\n");
}  /* printset */


void bigsubset(group_type *st, long n)
{
  /* Find a maximal subset of st among the n groupings,
     to be the set at the base of the tree.  */
  long i, j;
  group_type *su;
  boolean max, same;

  su = (group_type *)Malloc(setsz * sizeof(group_type));
  for (i = 0; i < setsz; i++)
    su[i] = 0;
  for (i = 0; i < n; i++) {
    max = true;
    for (j = 0; j < setsz; j++)
      if ((grouping[i][j] & ~st[j]) != 0)
        max = false;
    if (max) {
      same = true;
      for (j = 0; j < setsz; j++)
        if (grouping[i][j] != st[j])
          same = false;
      max = !same;
    }
    if (max) {
      for (j = 0; j < setsz; j ++)
        if ((su[j] & ~grouping[i][j]) != 0)
          max = false;
      if (max) {
        same = true;
        for (j = 0; j < setsz; j ++)
          if (su[j] != grouping[i][j])
            same = false;
        max = !same;
      }
      if (max)
        memcpy(su, grouping[i], setsz * sizeof(group_type));
    }
  }
  memcpy(st, su, setsz * sizeof(group_type));
  free(su);
}  /* bigsubset */


void recontraverse(node **p, group_type *st, long n, long *nextnode)
{
  /* traverse to add next node to consensus tree */
  long i, j = 0, k = 0, l = 0;

  boolean found, same = 0, zero, zero2;
  group_type *tempset, *st2;
  node *q, *r;

  for (i = 1; i <= spp; i++) {  /* count species in set */
    if (i == ((l+1)*SETBITS+1)) l++;
    if (((1L << (i - 1 - l*SETBITS)) & st[l]) != 0) {
      k++;               /* k  is the number of species in the set */
      j = i;             /* j  is set to last species in the set */
    }
  }
  if (k == 1) {           /* if only 1, set up that tip */
    *p = nodep[j - 1];
    (*p)->tip = true;
    (*p)->index = j;
    return;
  }
  gnu(&grbg, p);          /* otherwise make interior node */
  (*p)->tip = false;
  (*p)->index = *nextnode;
  nodep[*nextnode - 1] = *p;
  (*nextnode)++;
  (*p)->deltav = 0.0;
  for (i = 0; i < n; i++) { /* go through all sets */
    same = true;            /* to find one which is this one */
    for (j = 0; j < setsz; j++)
      if (grouping[i][j] != st[j])
        same = false;
    if (same)
      (*p)->deltav = *timesseen[i];
  }
  tempset = (group_type *)Malloc(setsz * sizeof(group_type));
  memcpy(tempset, st, setsz * sizeof(group_type));
  q = *p;
  st2 = (group_type *)Malloc(setsz * sizeof(group_type));
  memcpy(st2, st, setsz * sizeof(group_type));
  zero = true;      /* having made two copies of the set ... */
  for (j = 0; j < setsz; j++)       /* see if they are empty */
    if (tempset[j] != 0)
      zero = false;
  if (!zero)
    bigsubset(tempset, n);        /* find biggest set within it */
  zero = zero2 = false;           /* ... tempset is that subset */
  while (!zero && !zero2) {
    zero = zero2 = true;
    for (j = 0; j < setsz; j++) {
      if (st2[j] != 0)
        zero = false;
      if (tempset[j] != 0)
        zero2 = false;
    }
    if (!zero && !zero2) {
      gnu(&grbg, &q->next);
      q->next->index = q->index;
      q = q->next;
      q->tip = false;
      r = *p;
      recontraverse(&q->back, tempset, n, nextnode); /* put it on tree */
      *p = r;
      q->back->back = q;
      for (j = 0; j < setsz; j++)
        st2[j] &= ~tempset[j];     /* remove that subset from the set */
      memcpy(tempset, st2, setsz * sizeof(group_type));  /* that becomes set */
      found = false;
      i = 1;
      while (!found && i <= n) {
        if (grouping[i - 1] != 0) {
          same = true;
          for (j = 0; j < setsz; j++)
            if (grouping[i - 1][j] != tempset[j])
              same = false;
        }
        if ((grouping[i - 1] != 0) && same)
          found = true;
        else
          i++;
      }
      zero = true;
      for (j = 0; j < setsz; j++)
        if (tempset[j] != 0)
          zero = false;
      if (!zero && !found)
        bigsubset(tempset, n);
    }
  }
  q->next = *p;
  free(tempset);
  free(st2);
}  /* recontraverse */


void reconstruct(long n)
{
  /* reconstruct tree from the subsets */
  long nextnode;
  group_type *s;

  nextnode = spp + 1;
  s = (group_type *)Malloc(setsz * sizeof(group_type));
  memcpy(s, fullset, setsz * sizeof(group_type));
  recontraverse(&root, s, n, &nextnode);
  free(s);
}  /* reconstruct */


void coordinates(node *p, long *tipy)
{
  /* establishes coordinates of nodes */
  node *q, *first, *last;
  long maxx;

  if (p->tip) {
    p->xcoord = 0;
    p->ycoord = *tipy;
    p->ymin = *tipy;
    p->ymax = *tipy;
    (*tipy) += down;
    return;
  }
  q = p->next;
  maxx = 0;
  while (q != p) {
    coordinates(q->back, tipy);
    if (!q->back->tip) {
      if (q->back->xcoord > maxx)
        maxx = q->back->xcoord;
    }
    q = q->next;
  }
  first = p->next->back;
  q = p;
  while (q->next != p)
    q = q->next;
  last = q->back;
  p->xcoord = maxx + OVER;
  p->ycoord = (long)((first->ycoord + last->ycoord) / 2);
  p->ymin = first->ymin;
  p->ymax = last->ymax;
}  /* coordinates */


void drawline(long i)
{
  /* draws one row of the tree diagram by moving up tree */
  node *p, *q;
  long n, j;
  boolean extra, done, trif;
  node *r,  *first = NULL, *last = NULL;
  boolean found;

  p = root;
  q = root;
  fprintf(outfile, "  ");
  extra = false;
  trif = false;
  do {
    if (!p->tip) {
      found = false;
      r = p->next;
      while (r != p && !found) {
        if (i >= r->back->ymin && i <= r->back->ymax) {
          q = r->back;
          found = true;
        } else
          r = r->next;
      }
      first = p->next->back;
      r = p;
      while (r->next != p)
        r = r->next;
      last = r->back;
    }
    done = (p->tip || p == q);
    n = p->xcoord - q->xcoord;
    if (extra) {
      n--;
      extra = false;
    }
    if (q->ycoord == i && !done) {
      if (trif)
        putc('-', outfile);
      else
        putc('+', outfile);
      trif = false;
      if (!q->tip) {
        for (j = 1; j <= n - 8; j++)
          putc('-', outfile);
        if (noroot && (root->next->next->next == root) &&
            (((root->next->back == q) && root->next->next->back->tip)
             || ((root->next->next->back == q) && root->next->back->tip)))
          fprintf(outfile, "-------|");
        else {
          if (!strict) {   /* write number of times seen */
            if (q->deltav >= 10000)
              fprintf(outfile, "-%5.0f-|", (double)q->deltav);
            else if (q->deltav >= 1000)
              fprintf(outfile, "--%4.0f-|", (double)q->deltav);
            else if (q->deltav >= 100)
              fprintf(outfile, "-%5.1f-|", (double)q->deltav);
            else if (q->deltav >= 10)
              fprintf(outfile, "--%4.1f-|", (double)q->deltav);
            else
              fprintf(outfile, "--%4.2f-|", (double)q->deltav);
          } else
            fprintf(outfile, "-------|");
        }
        extra = true;
        trif = true;
      } else {
        for (j = 1; j < n; j++)
          putc('-', outfile);
      }
    } else if (!p->tip && last->ycoord > i && first->ycoord < i &&
               (i != p->ycoord || p == root)) {
      putc('|', outfile);
      for (j = 1; j < n; j++)
        putc(' ', outfile);
    } else {
      for (j = 1; j <= n; j++)
        putc(' ', outfile);
      if (trif)
        trif = false;
    }
    if (q != p)
      p = q;
  } while (!done);
  if (p->ycoord == i && p->tip) {
    for (j = 0; (j < MAXNCH) && (p->nayme[j] != '\0'); j++)
      putc(p->nayme[j], outfile);
  }
  putc('\n', outfile);
}  /* drawline */


void printree()
{
  /* prints out diagram of the tree */
  long i;
  long tipy;

  if (treeprint) {
    fprintf(outfile, "\nCONSENSUS TREE:\n");
    if (mr || mre || ml) {
      if (noroot) {
        fprintf(outfile, "the numbers on the branches indicate the number\n");
 fprintf(outfile, "of times the partition of the species into the two sets\n");
        fprintf(outfile, "which are separated by that branch occurred\n");
      } else {
        fprintf(outfile, "the numbers forks indicate the number\n");
        fprintf(outfile, "of times the group consisting of the species\n");
        fprintf(outfile, "which are to the right of that fork occurred\n");
      }
      fprintf(outfile, "among the trees, out of %6.2f trees\n", ntrees);
      if (ntrees <= 1.001)
        fprintf(outfile, "(trees had fractional weights)\n");
    }
    tipy = 1;
    coordinates(root, &tipy);
    putc('\n', outfile);
    for (i = 1; i <= tipy - down; i++)
      drawline(i);
    putc('\n', outfile);
  }
  if (noroot) {
    fprintf(outfile, "\n  remember:");
    if (didreroot)
      fprintf(outfile, " (though rerooted by outgroup)");
    fprintf(outfile, " this is an unrooted tree!\n");
  }
  putc('\n', outfile);
}  /* printree */


void enterpartition (group_type *s1, long *n)
{
  /* try to put this partition in list of partitions.  If implied by others,
     don't bother.  If others implied by it, replace them.  If this one
     vacuous because only one element in s1, forget it */
  long i, j;
  boolean found;

/* this stuff all to be rewritten but left here so pieces can be used */
  found = false;
  for (i = 0; i < (*n); i++) {  /* go through looking whether it is there */
    found = true;
    for (j = 0; j < setsz; j++) {  /* check both parts of partition */
      found = found && (grouping[i][j] == s1[j]);
      found = found && (group2[i][j] == (fullset[j] & (~s1[j])));
    }
    if (found)
      break;
  }
  if (!found) {    /* if not, add it to the slot after the end,
                      which must be empty */
    grouping[i] = (group_type *)Malloc(setsz * sizeof(group_type));
    timesseen[i] = (double *)Malloc(sizeof(double));
    group2[i] = (group_type *)Malloc(setsz * sizeof(group_type));
    for (j = 0; j < setsz; j++)
      grouping[i][j] = s1[j];
    *timesseen[i] = 1;
    (*n)++;
  }
} /* enterpartition */


void elimboth(long n)
{
  /* for Adams case: eliminate pairs of groups incompatible with each other */
  long i, j;
  boolean comp;

  for (i = 0; i < n-1; i++) {
    for (j = i+1; j < n; j++) {
      comp = compatible(i,j);
      if (!comp) {
        *timesseen[i] = 0.0;
        *timesseen[j] = 0.0;
      }
    }
    if (*timesseen[i] == 0.0) {
      free(grouping[i]);
      free(timesseen[i]);
      timesseen[i] = NULL;
      grouping[i] = NULL;
    }
  }
  if (*timesseen[n-1] == 0.0) {
    free(grouping[n-1]);
    free(timesseen[n-1]);
    timesseen[n-1] = NULL;
    grouping[n-1] = NULL;
  }
}  /* elimboth */


void consensus(pattern_elm ***pattern_array, long trees_in)
{
  long i, n, n2, tipy;

  group2 = (group_type **)  Malloc(maxgrp*sizeof(group_type *));
  for (i = 0; i < maxgrp; i++)
    group2[i] = NULL;
  times2 = (double **)Malloc(maxgrp*sizeof(double *));
  for (i = 0; i < maxgrp; i++)
    times2[i] = NULL;
  n2 = 0;
  censor();                /* drop groups that are too rare */
  compress(&n);            /* push everybody to front of array */
  if (!strict) {           /* drop those incompatible, if any */
    sort(n);
    eliminate(&n, &n2);
    compress(&n);
    }
  reconstruct(n);
  tipy = 1;
  coordinates(root, &tipy);
  if (prntsets) {
    fprintf(outfile, "\nSets included in the consensus tree\n");
    printset(n);
    for (i = 0; i < n2; i++) {
      if (!grouping[i]) {
        grouping[i] = (group_type *)Malloc(setsz * sizeof(group_type));
        timesseen[i] = (double *)Malloc(sizeof(double));
        }
      memcpy(grouping[i], group2[i], setsz * sizeof(group_type));
      *timesseen[i] = *times2[i];
    }
    n = n2;
    fprintf(outfile, "\n\nSets NOT included in consensus tree:");
    if (n2 == 0)
      fprintf(outfile, " NONE\n");
    else {
      putc('\n', outfile);
      printset(n);
    }
  }
  putc('\n', outfile);
  if (strict)
    fprintf(outfile, "\nStrict consensus tree\n");
  if (mre)
    fprintf(outfile, "\nExtended majority rule consensus tree\n");
  if (ml) {
    fprintf(outfile, "\nM  consensus tree (l = %4.2f)\n", mlfrac);
    fprintf(outfile, " l\n");
  }
  if (mr)
    fprintf(outfile, "\nMajority rule consensus tree\n");
  printree();
  free(nayme);  
  for (i = 0; i < maxgrp; i++)
    free(grouping[i]);
  free(grouping);
  for (i = 0; i < maxgrp; i++)
    free(order[i]);
  free(order);
  for (i = 0; i < maxgrp; i++)
    if (timesseen[i] != NULL)
      free(timesseen[i]);
  free(timesseen);
}  /* consensus */


void rehash()
{
  group_type *s;
  long i, j;
  double temp, ss, smult;
  boolean done;

  long old_maxgrp = maxgrp;
  long new_maxgrp = maxgrp*2;

  tmseen2 =     (double **)Malloc(new_maxgrp*sizeof(double *));
  grping2 = (group_type **)Malloc(new_maxgrp*sizeof(group_type *));
  order2 =        (long **)Malloc(new_maxgrp*sizeof(long *));
  lengths2 =     (double *)Malloc(new_maxgrp*sizeof(double));
  tchange2 =     (double *)Malloc(new_maxgrp*sizeof(double));

  for (i = 0; i < new_maxgrp; i++)
  {
    tmseen2[i] = NULL;
    grping2[i] = NULL;
    order2[i] = NULL;
    lengths2[i] = 0.0;
    tchange2[i] = 0.0;
  }


  smult = (sqrt(5.0) - 1) / 2;
  s = (group_type *)Malloc(setsz * sizeof(group_type));

  for (i = 0; i < old_maxgrp; i++) {
    long old_index = *order[i];
    long new_index = -1;
    memcpy(s, grouping[old_index], setsz * sizeof(group_type));
    ss = 0.0;
    for (j = 0; j < setsz; j++)
      ss += s[j] /* pow(2, SETBITS*j)*/;
    temp = ss * smult;
    new_index = (long)(new_maxgrp * (temp - floor(temp)));
    done = false;
    while (!done) {
      if (!grping2[new_index])
      {

        grping2[new_index] = (group_type *)Malloc(setsz * sizeof(group_type));
        memcpy(grping2[new_index], grouping[old_index], setsz * sizeof(group_type));

        order2[i] = (long *)Malloc(sizeof(long));
        *order2[i] = new_index;

        tmseen2[new_index] = (double *)Malloc(sizeof(double));
        *tmseen2[new_index] = *timesseen[old_index];

        lengths2[new_index] = lengths[old_index];

        free(grouping[old_index]);
        free(timesseen[old_index]);
        free(order[i]);

        grouping[old_index] = NULL;
        timesseen[old_index] = NULL;
        order[i] = NULL;

        done = true; /* successfully found place for this item */

      } else {
        new_index++;
        if (new_index >= new_maxgrp) new_index -= new_maxgrp;
      }
    }
  }

  free(lengths);
  free(timesseen);
  free(grouping);
  free(order);

  free(s);

  timesseen = tmseen2;
  grouping = grping2;
  lengths = lengths2;
  order = order2;

  maxgrp = new_maxgrp;

}  /* rehash */


void enternodeset(node* r)
{ /* enter a set of species into the hash table */
  long i, j, start;
  double ss, n;
  boolean done, same;
  double times ;
  group_type *s;

  s = r->nodeset;


  /* do not enter full sets */
  same = true;
  for (i = 0; i < setsz; i++)
    if (s[i] != fullset[i])
      same = false;
  if (same) 
    return;

  times = trweight;
  ss = 0.0;                        /* compute the hashcode for the set */
  n = ((sqrt(5.0) - 1.0) / 2.0);   /* use an irrational multiplier */
  for (i = 0; i < setsz; i++)
    ss += s[i] * n;
  i = (long)(maxgrp * (ss - floor(ss))) + 1; /* use fractional part of code */
  start = i;
  done = false;                   /* go through seeing if it is there */

  while (!done) {
    if (grouping[i - 1]) {        /* ... i.e. if group is absent, or  */
      same = false;               /* (will be false if timesseen = 0) */
      if (!(timesseen[i-1] == 0)) {
        same = true;
        for (j = 0; j < setsz; j++) {
          if (s[j] != grouping[i - 1][j])
            same = false;
        }
      }
      else {                       /* if group is present but timessen = 0 */
        for (j = 0; j < setsz; j++)   /* replace by correct group */
          grouping[i - 1][j] = s[j];
        *timesseen[i-1] = 1;
      }
    }
    if (grouping[i - 1] && same) {  /* if it is there, increment timesseen */
      *timesseen[i - 1] += times;
      lengths[i - 1] = nodep[r->index - 1]->v;
      done = true;
    } else if (!grouping[i - 1]) {  /* if not there and slot empty ... */
      grouping[i - 1] = (group_type *)Malloc(setsz * sizeof(group_type));
      lasti++;
      order[lasti] = (long *)Malloc(sizeof(long));
      timesseen[i - 1] = (double *)Malloc(sizeof(double));
      memcpy(grouping[i - 1], s, setsz * sizeof(group_type));
      *timesseen[i - 1] = times;
      *order[lasti] = i - 1;
      done = true;
      lengths[i - 1] = nodep[r->index -1]->v;
    } else {  /* otherwise look to put it in next slot ... */
        i++;
        if (i > maxgrp) i -= maxgrp;
      }
    if (!done && i == start) {  /* if no place to put it, expand hash table */

      rehash();

      done = true;
      enternodeset(r); /* calls this procedure again, but now there
                          should be space */
    }
  }
}  /* enternodeset */


/* recursively crawls through tree, setting nodeset values to be the
 * bitwise OR of bits from downstream nodes
 */
void accumulate(node *r)
{
  node *q;
  long i;

  /* zero out nodeset values. since we are re-using tree nodes,
   * the malloc only happens the first time we encounter a node. */
  if (!r->nodeset)
  {
    r->nodeset = (group_type *)Malloc(setsz * sizeof(group_type));
  }
  for (i = 0; i < setsz; i++)
  {
    r->nodeset[i] = 0L;
  }

  if (r->tip) {
    /* tip nodes should have a single bit set corresponding to index-1 */
    i = (r->index-1) / (long)SETBITS;
    r->nodeset[i] = 1L << (r->index - 1 - i*SETBITS);
  }
  else {
    /* for loop should not visit r->back -- we've likely come from there */
    for (q = r->next; q != r; q = q->next) {

      /* recursive call to this function */
      accumulate(q->back);

      /* bitwise OR of bits from downstream nodes */
      for (i = 0; i < setsz; i++)
        r->nodeset[i] |= q->back->nodeset[i];
    }
  }

  if ((!r->tip && (r->next->next != r)) || r->tip) 
    enternodeset(r);
}  /* accumulate */


void dupname2(Char *name, node *p, node *this)
{
  /* search for a duplicate name recursively */
  node *q;

  if (p->tip) {
    if (p != this) {
      if (namesSearch(p->nayme)) {
        printf("\n\nERROR in user tree: duplicate name found: ");
        puts(p->nayme);
        printf("\n\n");
        exxit(-1);
      } else {
        namesAdd(p->nayme);
      }
    }
  } else {
    q = p;
    while (p->next != q) {
      dupname2(name, p->next->back, this);
      p = p->next;
    }
  }
}  /* dupname2 */


void dupname(node *p)
{
  /* Recursively searches tree, starting at p, to verify that
   * each tip name occurs only once. When called with root as
   * its argument, at final recusive exit, all tip names should 
   * be in the hash "hashp".
   */
  node *q;

  if (p->tip) {
    if (namesSearch(p->nayme)) {
      printf("\n\nERROR in user tree: duplicate name found: ");
      puts(p->nayme);
      printf("\n\n");
      exxit(-1);
    } else {
      namesAdd(p->nayme);
    }
  } else {
    q = p;
    while (p->next != q) {
      dupname(p->next->back);
      p = p->next;
    }
  }
}  /* dupname */


void missingnameRecurs(node *p)
{
  /* search for missing names in first tree */
  node *q;

  if (p->tip) {
    if (!namesSearch(p->nayme)) {
      printf("\n\nERROR in user tree: name %s not found in first tree\n\n\n",
             p->nayme);
      exxit(-1);
    }
  } else {
    q = p;
    while (p->next != q) {
      missingnameRecurs(p->next->back);
      p = p->next;
    }
  }
}  /* missingnameRecurs */

/**
 * wrapper for recursive missingname function 
 */
void missingname(node *p){
  missingnameRecurs(p);
  namesCheckTable();
} /* missingname */


void gdispose(node *p)
{
  /* go through tree throwing away nodes */
  node *q, *r;

  if (p->tip) {
    chuck(&grbg, p);
    return;
  }
  q = p->next;
  while (q != p) {
    gdispose(q->back);
    r = q;
    q = q->next;
    chuck(&grbg, r);
  }
  chuck(&grbg, p);
}  /* gdispose */


void initreenode(node *p)
{
  /* traverse tree and assign species names to tip nodes */
  node *q;

  if (p->tip) {
    memcpy(nayme[p->index - 1], p->nayme, MAXNCH);
  } else {
    q = p->next;
    while (q && q != p) {
      initreenode(q->back);
      q = q->next;
    }
  }
} /* initreenode */


void reroot(node *outgroup, long *nextnode)
{
  /* reroots and reorients tree, placing root at outgroup  */
  long i;
  node *p, *q;
  double newv;

  /* count root's children & find last */
  p = root;
  i = 0;
  while (p->next != root) {
    p = p->next;
    i++;
  }
  if (i == 2) {
    /* 2 children: */
    q = root->next;

    newv = q->back->v + p->back->v;
    
    /* if outgroup is already here, just move
     * its length to the other branch and finish */
    if (outgroup == p->back) {
      /* flip branch order at root so that outgroup 
       * is first, just to be consistent */
      root->next = p;
      p->next = q;
      q->next = root;
      
      q->back->v = newv;
      p->back->v = 0;
      return;
    }
    if (outgroup == q) {
      p->back->v = newv;
      q->back->v = 0;
      return;
    }
   
    /* detach root by linking child nodes */
    q->back->back = p->back;
    p->back->back = q->back;
    p->back->v = newv;
    q->back->v = newv;
  } else { /* 3+ children */
    p->next = root->next;              /* join old root nodes */
    nodep[root->index-1] = root->next; /* make root->next the primary node */
    
    /* create new root nodes */
    gnu(&grbg, &root->next);
    q = root->next;
    gnu(&grbg, &q->next);
    p = q->next;
    p->next = root;
    q->tip = false;
    p->tip = false;
    nodep[*nextnode] = root;
    (*nextnode)++;
    root->index = *nextnode;
    root->next->index = root->index;
    root->next->next->index = root->index;
  }
  newv = outgroup->v;
  /* root is 3 "floating" nodes */
  /* q == root->next */
  /* p == root->next->next */

  /* attach root at outgroup */
  q->back = outgroup;
  p->back = outgroup->back;
  outgroup->back->back = p;
  outgroup->back = q;
  outgroup->v = 0;
  outgroup->back->v = 0;
  root->v = 0;
  p->v = newv;
  p->back->v = newv;
  reorient(root);
}  /* reroot */


void reorient(node* n) {
  node* p;
  
  if ( n->tip ) return;
  if ( nodep[n->index - 1] != n )  {
    nodep[n->index - 1] = n;
    if ( n->back )
      n->v = n->back->v;
  }
  
  for ( p = n->next ; p != n ; p = p->next) 
    reorient(p->back);
}


void store_pattern (pattern_elm ***pattern_array, int trees_in_file)
{ 
  /* put a tree's groups into a pattern array.
     Don't forget that when not Adams, grouping[] is not compressed. . . */
  long i, total_groups=0, j=0, k;

  /* First, find out how many groups exist in the given tree. */
  for (i = 0 ; i < maxgrp ; i++)
    if ((grouping[i] != NULL) &&
		  (*timesseen[i] > 0))
      /* If this is group exists and is present in the current tree, */
      total_groups++ ;

  /* Then allocate a space to store the bit patterns. . . */
  for (i = 0 ; i < setsz ; i++) {
    pattern_array[i][trees_in_file]
      = (pattern_elm *) Malloc(sizeof(pattern_elm)) ;
    pattern_array[i][trees_in_file]->apattern = 
      (group_type *) Malloc (total_groups * sizeof (group_type)) ;
    pattern_array[i][trees_in_file]->length = 
      (double *) Malloc (maxgrp * sizeof (double)) ;
      for ( j = 0 ; j < maxgrp ; j++ ) {
        pattern_array[i][trees_in_file]->length[j] = -1;
      }
    pattern_array[i][trees_in_file]->patternsize = (long *)Malloc(sizeof(long));
  }
  j = 0;
  /* Then go through groupings again, and copy in each element
     appropriately. */
  for (i = 0 ; i < maxgrp ; i++)
    if (grouping[i] != NULL) {
      if (*timesseen[i] > 0) {
        for (k = 0 ; k < setsz ; k++)
          pattern_array[k][trees_in_file]->apattern[j] = grouping[i][k] ;  
        pattern_array[0][trees_in_file]->length[j] = lengths[i];
        j++ ;

          /* 
             EWFIX.BUG.756

             updates timesseen_changes to the current value
             pointed to by timesseen

             treedist uses this to determine if group i has been seen
             by comparing timesseen_changes[i] (the count now) with 
             timesseen[i] (the count after reading next tree)

             We could make treedist more efficient by not keeping
             timesseen (and groupings, etc) around, but doing it
             this way allows us to share code between treedist and
             consense.
             
          */
      *timesseen[i] = 0;
      }
    }
  *pattern_array[0][trees_in_file]->patternsize = total_groups;

}  /* store_pattern */


boolean samename(naym name1, plotstring name2)
{
  return !(strncmp(name1, name2, MAXNCH)); 
}  /* samename */


void reordertips()
{
  /* Reorders nodep[] and indexing to match species order from first tree */
  /* Assumes tree has spp tips and nayme[] has spp elements, and that there is a
   * one-to-one mapping between tip names and the names in nayme[].
   */

  long i, j;
  node *t;

  for (i = 0; i < spp-1; i++) {
    for (j = i + 1; j < spp; j++) {
      if (samename(nayme[i], nodep[j]->nayme)) {
        /* switch the pointers in
         * nodep[] and set index accordingly for each node. */
        t = nodep[i];
        
        nodep[i] = nodep[j];
        nodep[i]->index = i+1;
        
        nodep[j] = t;
        nodep[j]->index = j+1;
        
        break;  /* next i */
      }
    }
  }
}  /* reordertips */

void read_groups (pattern_elm ****pattern_array, 
        long total_trees, long tip_count, FILE *intree)
{
  /* read the trees.  Accumulate sets. */
  int i, j, k;
  boolean haslengths, initial;
  long nextnode, trees_read = 0;

  /* do allocation first *****************************************/
  grouping  = (group_type **)  Malloc(maxgrp*sizeof(group_type *));
  lengths  = (double *)  Malloc(maxgrp*sizeof(double));

  for (i = 0; i < maxgrp; i++)
    grouping[i] = NULL;

  order     = (long **) Malloc(maxgrp*sizeof(long *));
  for (i = 0; i < maxgrp; i++)
    order[i] = NULL;

  timesseen = (double **)Malloc(maxgrp*sizeof(double *));
  for (i = 0; i < maxgrp; i++)
    timesseen[i] = NULL;

  nayme = (naym *)Malloc(tip_count*sizeof(naym));
  hashp = (hashtype)Malloc(sizeof(namenode) * NUM_BUCKETS);
  for (i=0;i<NUM_BUCKETS;i++) {
      hashp[i] = NULL;
  }
  setsz = (long)ceil((double)tip_count/(double)SETBITS);
  if (tree_pairing != NO_PAIRING)
    {
      /* Now that we know setsz, we can malloc pattern_array
      and pattern_array[n] accordingly. */
      (*pattern_array) =
        (pattern_elm ***)Malloc(setsz * sizeof(pattern_elm **));

      /* For this assignment, let's assume that there will be no
         more than maxtrees. */
      for (j = 0 ; j < setsz ; j++) 
      {
          (*pattern_array)[j] = 
            (pattern_elm **)Malloc(total_trees * sizeof(pattern_elm *));
        for(k = 0 ; k < total_trees ; k++ )
        {
            (*pattern_array)[j][k] = NULL;
        }
      }
    }

  fullset = (group_type *)Malloc(setsz * sizeof(group_type));
  for (j = 0; j < setsz; j++)
    fullset[j] = 0L;
  k = 0;
  for (j = 1; j <= tip_count; j++) {
    if (j == ((k+1)*SETBITS+1)) k++;
    fullset[k] |= 1L << (j - k*SETBITS - 1);
  }
  /* end allocation **********************************************/

  firsttree = true;
  grbg = NULL;
  initial = true;
  while (!eoff(intree)) {          /* go till end of input tree file */
    for (i = 0; i < maxgrp; i++) {
      lengths[i] = -1;
    }
    goteof = false;
    nextnode = 0;
    haslengths = true;
    allocate_nodep(&nodep, &intree, &spp);
    assert(spp == tip_count);
    treeread(intree, &root, treenode, &goteof, &firsttree, nodep, 
              &nextnode, &haslengths, &grbg, initconsnode,true,-1);
    if (!initial) { 
      missingname(root);
      reordertips();
    } else {
      initial = false;
      dupname(root);
      initreenode(root);
    }
    if (goteof)
      continue;
    ntrees += trweight;
    if (noroot) {
      reroot(nodep[outgrno - 1], &nextnode);
      didreroot = outgropt;
    }
    accumulate(root);
    gdispose(root);
    for (j = 0; j < 2*(1+spp); j++)
      nodep[j] = NULL;
    free(nodep);
    /* Added by Dan F. */
    if (tree_pairing != NO_PAIRING) {
        /* If we're computing pairing or need separate tree sets, store the
           current pattern as an element of it's trees array. */
      store_pattern ((*pattern_array), trees_read) ;
      trees_read++ ;
      for (i = 0; i < maxgrp; i++)
        if (grouping[i] != NULL) {
          *timesseen[i] = 0;
        }
    }
  }
} /* read_groups */


void clean_up_final()
{
    long i;
    for(i=0;i<maxgrp;i++)
    {
        if(grouping[i] != NULL) {
            free(grouping[i]);
        }
        if(order[i] != NULL) {
            free(order[i]);
        }
        if(timesseen[i] != NULL) {
            free(timesseen[i]);
        }
    }
    free(grouping);
    free(nayme);
    free(order);
    free(timesseen);
    free(fullset);
    free(lengths);

    namesClearTable();
    free(hashp);
}
