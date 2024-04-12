#include "phylip.h"
#include "discrete.h"

/* version 3.6. (c) Copyright 1993-2004 by the University of Washington.
   Written by Joseph Felsenstein, Akiko Fuseki, Sean Lamont, and Andrew Keeffe.
   Permission is granted to copy and use this program provided no fee is
   charged for it and provided that this copyright notice is not removed. */

long nonodes, endsite, outgrno, nextree, which;
boolean interleaved, printdata, outgropt, treeprint, dotdiff;
steptr weight, category, alias, location, ally;
sequence y, convtab;

void inputdata(long chars)
{
  /* input the names and sequences for each species */
  /* used by pars */
  long i, j, k, l;
  long basesread=0, basesnew=0, nsymbol=0, convsymboli=0;
  Char charstate;
  boolean allread, done, found;

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
        initname(i - 1);
      j = (interleaved) ? basesread : 0;
      done = false;
      while (!done && !eoff(infile)) {
        if (interleaved)
          done = true;
        while (j < chars && !(eoln(infile) || eoff(infile))) {
          charstate = gettc(infile);
          if (charstate == '\n' || charstate == '\t')
            charstate = ' ';
          if (charstate == ' ')
            continue;
          if ((strchr("!\"#$%&'()*+,-./0123456789:;<=>?@\
                       ABCDEFGHIJKLMNOPQRSTUVWXYZ[\\]^_`\
                       abcdefghijklmnopqrstuvwxyz{|}~",charstate)) == NULL){
            printf(
                "\n\nERROR: Bad symbol: %c at position %ld of species %ld\n\n",
                   charstate, j+1, i);
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
        printf("\n\nERROR: Sequences out of alignment at position %ld\n\n", j);
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
  if (printdata) {
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
  }
  for (i = 1; i <= chars; i++) {
    nsymbol = 0;
    for (j = 1; j <= spp; j++) {
      if ((nsymbol == 0) && (y[j - 1][i - 1] != '?')) {
        nsymbol = 1;
        convsymboli = 1;
        convtab[0][i-1] = y[j-1][i-1];
      } else if (y[j - 1][i - 1] != '?'){
        found = false;
        for (k = 1; k <= nsymbol; k++) {
          if (convtab[k - 1][i - 1] == y[j - 1][i - 1]) {
            found = true;
            convsymboli = k;
          }
        }
        if (!found) {
          nsymbol++;
          convtab[nsymbol-1][i - 1] = y[j - 1][i - 1];
          convsymboli = nsymbol;
        } 
      }
      if (nsymbol <= 8) {
        if (y[j - 1][i - 1] != '?')
          y[j - 1][i - 1] = (Char)('0' + (convsymboli - 1));
      } else {
        printf(
          "\n\nERROR: More than maximum of 8 symbols in column %ld\n\n", i);
        exxit(-1);
      }
    }
  }
}  /* inputdata */


void alloctree(pointarray *treenode, long nonodes, boolean usertree)
{
  /* allocate treenode dynamically */
  /* used in pars */
  long i, j;
  node *p, *q;

  *treenode = (pointarray)Malloc(nonodes*sizeof(node *));
  for (i = 0; i < spp; i++) {
    (*treenode)[i] = (node *)Malloc(sizeof(node));
    (*treenode)[i]->tip = true;
    (*treenode)[i]->index = i+1;
    (*treenode)[i]->iter = true;
    (*treenode)[i]->branchnum = i+1;
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
        p->branchnum = i+1;
        p->initialized = false;
        p->next = q;
        q = p;
      }
      p->next->next->next = p;
      (*treenode)[i] = p;
    }
} /* alloctree */


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
        p = p->next;
      }
    }
  }
} /* setuptree */


void alloctip(node *p, long *zeros, unsigned char *zeros2)
{ /* allocate a tip node */
  /* used by pars */

  p->numsteps = (steptr)Malloc(endsite*sizeof(long));
  p->oldnumsteps = (steptr)Malloc(endsite*sizeof(long));
  p->discbase = (discbaseptr)Malloc(endsite*sizeof(unsigned char));
  p->olddiscbase = (discbaseptr)Malloc(endsite*sizeof(unsigned char));
  memcpy(p->discbase, zeros2, endsite*sizeof(unsigned char));
  memcpy(p->numsteps, zeros, endsite*sizeof(long));
  memcpy(p->olddiscbase, zeros2, endsite*sizeof(unsigned char));
  memcpy(p->oldnumsteps, zeros, endsite*sizeof(long));
}  /* alloctip */


void sitesort(long chars, steptr weight)
{
  /* Shell sort keeping sites, weights in same order */
  /* used in pars */
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
  /* used in pars */
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
  /* used in pars */
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


void makevalues(pointarray treenode, long *zeros, unsigned char *zeros2,
                        boolean usertree)
{
  /* set up fractional likelihoods at tips */
  /* used by pars */
  long i, j;
  unsigned char ns=0;
  node *p;

  setuptree(treenode, nonodes, usertree);
  for (i = 0; i < spp; i++)
    alloctip(treenode[i], zeros, zeros2);
  if (!usertree) {
    for (i = spp; i < nonodes; i++) {
      p = treenode[i];
      do {
        allocdiscnontip(p, zeros, zeros2, endsite);
        p = p->next;
      } while (p != treenode[i]);
    }
  }
  for (j = 0; j < endsite; j++) {
    for (i = 0; i < spp; i++) {
      switch (y[i][alias[j] - 1]) {

      case '0':
        ns = 1 << zero;
        break;

      case '1':
        ns = 1 << one;
        break;

      case '2':
        ns = 1 << two;
        break;

      case '3':
        ns = 1 << three;
        break;

      case '4':
        ns = 1 << four;
        break;

      case '5':
        ns = 1 << five;
        break;

      case '6':
        ns = 1 << six;
        break;

      case '7':
        ns = 1 << seven;
        break;

      case '?':
        ns = (1 << zero) | (1 << one) | (1 << two) | (1 << three) |
               (1 << four) | (1 << five) | (1 << six) | (1 << seven);
        break;
      }
      treenode[i]->discbase[j] = ns;
      treenode[i]->numsteps[j] = 0;
    }
  }
}  /* makevalues */


void fillin(node *p, node *left, node *rt)
{
  /* sets up for each node in the tree the base sequence
     at that point and counts the changes.  */
  long i, j, k, n;
  node *q;

  if (!left) {
    memcpy(p->discbase, rt->discbase, endsite*sizeof(unsigned char));
    memcpy(p->numsteps, rt->numsteps, endsite*sizeof(long));
    q = rt;
  } else if (!rt) {
    memcpy(p->discbase, left->discbase, endsite*sizeof(unsigned char));
    memcpy(p->numsteps, left->numsteps, endsite*sizeof(long));
    q = left;
  } else {
    for (i = 0; i < endsite; i++) {
      p->discbase[i] = left->discbase[i] & rt->discbase[i];
      p->numsteps[i] = left->numsteps[i] + rt->numsteps[i];
      if (p->discbase[i] == 0) {
        p->discbase[i] = left->discbase[i] | rt->discbase[i];
        p->numsteps[i] += weight[i];
      }
    }
    q = rt;
  }
  if (left && rt) n = 2;
  else n = 1;
  for (i = 0; i < endsite; i++)
    for (j = (long)zero; j <= (long)seven; j++)
      p->discnumnuc[i][j] = 0;
  for (k = 1; k <= n; k++) {
    if (k == 2) q = left;
    for (i = 0; i < endsite; i++) {
      for (j = (long)zero; j <= (long)seven; j++) {
        if (q->discbase[i] & (1 << j))
          p->discnumnuc[i][j]++;
      }
    }
  }
}  /* fillin */


long getlargest(long *discnumnuc)
{
  /* find the largest in array numnuc */
  long i, largest;

  largest = 0;
  for (i = (long)zero; i <= (long)seven; i++)
    if (discnumnuc[i] > largest)
      largest = discnumnuc[i];
  return largest;
} /* getlargest */


void multifillin(node *p, node *q, long dnumdesc)
{
  /* sets up for each node in the tree the base sequence
     at that point and counts the changes according to the
     changes in q's base */
  long i, j, largest, descsteps;
  unsigned char b;

  memcpy(p->olddiscbase, p->discbase, endsite*sizeof(unsigned char));
  memcpy(p->oldnumsteps, p->numsteps, endsite*sizeof(long));
  for (i = 0; i < endsite; i++) {
    descsteps = 0;
    for (j = (long)zero; j <= (long)seven; j++) {
      b = 1 << j;
      if ((descsteps == 0) && (p->discbase[i] & b)) 
        descsteps = p->numsteps[i] 
                      - (p->numdesc - dnumdesc - p->discnumnuc[i][j])
                        * weight[i];
    }
    if (dnumdesc == -1)
      descsteps -= q->oldnumsteps[i];
    else if (dnumdesc == 0)
      descsteps += (q->numsteps[i] - q->oldnumsteps[i]);
    else
      descsteps += q->numsteps[i];
    if (q->olddiscbase[i] != q->discbase[i]) {
      for (j = (long)zero; j <= (long)seven; j++) {
        b = 1 << j;
        if ((q->olddiscbase[i] & b) && !(q->discbase[i] & b))
          p->discnumnuc[i][j]--;
        else if (!(q->olddiscbase[i] & b) && (q->discbase[i] & b))
          p->discnumnuc[i][j]++;
      }
    }
    largest = getlargest(p->discnumnuc[i]);
    if (q->olddiscbase[i] != q->discbase[i]) {
      p->discbase[i] = 0;
      for (j = (long)zero; j <= (long)seven; j++) {
        if (p->discnumnuc[i][j] == largest)
            p->discbase[i] |= (1 << j);
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
  unsigned char ns, rs, ls;

  if (!left) {
    memcpy(p->numsteps, rt->numsteps, endsite*sizeof(long));
    memcpy(p->discbase, rt->discbase, endsite*sizeof(unsigned char));
  } else if (!rt) {
    memcpy(p->numsteps, left->numsteps, endsite*sizeof(long));
    memcpy(p->discbase, left->discbase, endsite*sizeof(unsigned char));
  } else 
    for (i = a; i < b; i++) {
      ls = left->discbase[i];
      rs = rt->discbase[i];
      ns = ls & rs;
      p->numsteps[i] = left->numsteps[i] + rt->numsteps[i];
      if (ns == 0) {
        ns = ls | rs;
        p->numsteps[i] += weight[i];
      }
      p->discbase[i] = ns;
    }
}  /* sumnsteps */


void sumnsteps2(node *p, node *left, node *rt, long a, long b, long *threshwt)
{
  /* counts the changes at each node.  */
  long i, steps;
  unsigned char ns, rs, ls;
  long term;

  if (a == 0) p->sumsteps = 0.0;
  if (!left)
    memcpy(p->numsteps, rt->numsteps, endsite*sizeof(long));
  else if (!rt)
    memcpy(p->numsteps, left->numsteps, endsite*sizeof(long));
  else 
    for (i = a; i < b; i++) {
      ls = left->discbase[i];
      rs = rt->discbase[i];
      ns = ls & rs;
      p->numsteps[i] = left->numsteps[i] + rt->numsteps[i];
      if (ns == 0)
        p->numsteps[i] += weight[i];
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
  /* sets up for each node in the tree the base sequence
     at that point and counts the changes according to the
     changes in q's base */
  long i, j, steps, largest, descsteps;
  long term;

  if (a == 0) p->sumsteps = 0.0;
  for (i = a; i < b; i++) {
    descsteps = 0;
    for (j = (long)zero; j <= (long)seven; j++) {
      if ((descsteps == 0) && (p->discbase[i] & (1 << j))) 
        descsteps = p->numsteps[i] -
                        (p->numdesc - 1 - p->discnumnuc[i][j]) * weight[i];
    }
    descsteps += q->numsteps[i];
    largest = 0;
    for (j = (long)zero; j <= (long)seven; j++) {
      if (q->discbase[i] & (1 << j))
        p->discnumnuc[i][j]++;
      if (p->discnumnuc[i][j] > largest)
        largest = p->discnumnuc[i][j];
    }
    steps = ((p->numdesc - largest) * weight[i] + descsteps);
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
  long i, j, largest;
  node *q;
  discbaseptr b;

  for (i = 0; i < endsite; i++) {
    p->numsteps[i] = 0;
    q = p->next;
    while (q != p) {
      if (q->back) {
        p->numsteps[i] += q->back->numsteps[i];
        b = q->back->discbase;
        for (j = (long)zero; j <= (long)seven; j++)
          if (b[i] & (1 << j))
            p->discnumnuc[i][j]++;
      }
      q = q->next;
    }
    largest = getlargest(p->discnumnuc[i]);
    p->numsteps[i] += ((p->numdesc - largest) * weight[i]);
    p->discbase[i] = 0;
    for (j = (long)zero; j <= (long)seven; j++) {
      if (p->discnumnuc[i][j] == largest)
       p->discbase[i] |= (1 << j);
    }
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
          memcpy(p->olddiscbase, p->discbase, endsite*sizeof(unsigned char));
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
}


void add(node *below, node *newtip, node *newfork, node **root, boolean recompute,
                        pointarray treenode, node **grbg, long *zeros, unsigned char *zeros2)
{
  /* inserts the nodes newfork and its left descendant, newtip,
     to the tree.  below becomes newfork's right descendant.
     if newfork is NULL, newtip is added as below's sibling */
  /* used in pars */
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
    gnudisctreenode(grbg, &p, below->index, endsite, zeros, zeros2);
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
    memcpy(newtip->back->discbase, below->discbase, endsite*sizeof(unsigned char));
    memcpy(newtip->back->numsteps, below->numsteps, endsite*sizeof(long));
    memcpy(newtip->back->discnumnuc, below->discnumnuc, endsite*sizeof(discnucarray));
    if (below != *root) {
      memcpy(below->back->olddiscbase, zeros2, endsite*sizeof(unsigned char));
      memcpy(below->back->oldnumsteps, zeros, endsite*sizeof(long));
      multifillin(newtip->back, below->back, 1);
    }
    if (!newtip->tip) {
      memcpy(newtip->back->olddiscbase, zeros2, endsite*sizeof(unsigned char));
      memcpy(newtip->back->oldnumsteps, zeros, endsite*sizeof(long));
      preorder(newtip, newtip->back, *root, NULL, NULL, below, 1);
    }
    memcpy(newtip->olddiscbase, zeros2, endsite*sizeof(unsigned char));
    memcpy(newtip->oldnumsteps, zeros, endsite*sizeof(long));
    preorder(below, newtip, *root, NULL, newtip, below, 1);
    if (below != *root)
      preorder(below->back, below, *root, NULL, NULL, NULL, 0);
  } else {
    fillin(newtip->back, newtip->back->next->back,
             newtip->back->next->next->back);
    if (!newtip->tip) {
      memcpy(newtip->back->olddiscbase, zeros2, endsite*sizeof(unsigned char));
      memcpy(newtip->back->oldnumsteps, zeros, endsite*sizeof(long));
      preorder(newtip, newtip->back, *root, NULL, NULL, newfork, 1);
    }
    if (newfork != *root) {
      memcpy(below->back->discbase, newfork->back->discbase, endsite*sizeof(unsigned char));
      memcpy(below->back->numsteps, newfork->back->numsteps, endsite*sizeof(long));
      preorder(newfork, newtip, *root, NULL, newtip, NULL, 0);
    } else {
      fillin(below->back, newtip, NULL);
      fillin(newfork, newtip, below);
      memcpy(below->back->olddiscbase, zeros2, endsite*sizeof(unsigned char));
      memcpy(below->back->oldnumsteps, zeros, endsite*sizeof(long));
      preorder(below, below->back, *root, NULL, NULL, newfork, 1);
    }
    if (newfork != *root) {
      memcpy(newfork->olddiscbase, below->discbase, endsite*sizeof(unsigned char));
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
                        pointarray treenode, node **grbg, long *zeros, unsigned char *zeros2)
{
  /* removes nodes item and its ancestor, fork, from the tree.
     the new descendant of fork's ancestor is made to be
     fork's second descendant (other than item).  Also
     returns pointers to the deleted nodes, item and fork.
     If item belongs to a node with more than 2 descendants,
     fork will not be deleted */
  /* used in pars */
  node *p, *q, *other, *otherback = NULL;

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
      if (other->tip)
        *root = NULL;
      else {
        *root = other;
        updatenumdesc(other, *root, other->numdesc);
      }
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
      memcpy(item->back->olddiscbase, item->back->discbase,
               endsite*sizeof(unsigned char));
      memcpy(item->back->oldnumsteps, item->back->numsteps, endsite*sizeof(long));
      memcpy(item->back->discbase, zeros2, endsite*sizeof(unsigned char));
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
    memcpy(otherback->olddiscbase, otherback->discbase,
             endsite*sizeof(unsigned char));  
    memcpy(otherback->oldnumsteps, otherback->numsteps, endsite*sizeof(long));
    if (other == *root) {
      memcpy(otherback->discbase, zeros2, endsite*sizeof(unsigned char));
      memcpy(otherback->numsteps, zeros, endsite*sizeof(long));
    } else {
      memcpy(otherback->discbase, other->back->discbase,
               endsite*sizeof(unsigned char));
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
      memcpy(other->olddiscbase,(*fork)->discbase, endsite*sizeof(unsigned char));
      memcpy(other->oldnumsteps,(*fork)->numsteps, endsite*sizeof(long));
      preorder(other->back, other, *root, NULL, NULL, NULL, 0);
    }
  } else {
    memcpy(item->olddiscbase, item->discbase, endsite*sizeof(unsigned char));
    memcpy(item->oldnumsteps, item->numsteps, endsite*sizeof(long));
    memcpy(item->discbase, zeros2, endsite*sizeof(unsigned char));
    memcpy(item->numsteps, zeros, endsite*sizeof(long));
    preorder(*fork, item, *root, NULL, NULL, *fork, -1);
    if (*fork != *root)
      preorder((*fork)->back, *fork, *root, NULL, NULL, NULL, 0);
    memcpy(item->discbase, item->olddiscbase, endsite*sizeof(unsigned char));
    memcpy(item->numsteps, item->oldnumsteps, endsite*sizeof(long));
  }
}  /* re_move */


void postorder(node *p)
{
  /* traverses an n-ary tree, suming up steps at a node's descendants */
  /* used in pars */
  node *q;

  if (p->tip)
    return;
  q = p->next;
  while (q != p) {
    postorder(q->back);
    q = q->next;
  }
  zerodiscnumnuc(p, endsite);
  if (p->numdesc > 2)
    multisumnsteps2(p);
  else
    fillin(p, p->next->back, p->next->next->back);
}  /* postorder */


void getnufork(node **nufork, node **grbg, pointarray treenode, long *zeros,
                        unsigned char *zeros2)
{
  /* find a fork not used currently */
  long i;

  i = spp;
  while (treenode[i] && treenode[i]->numdesc > 0) i++;
  if (!treenode[i])
    gnudisctreenode(grbg, &treenode[i], i, endsite, zeros, zeros2);
  *nufork = treenode[i];
} /* getnufork */


void reroot(node *outgroup, node *root)
{
  /* reorients tree, putting outgroup in desired position. used if
     the root is binary. */
  /* used in pars */
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
  /* used in pars */
  node *p;

  p = outgroup->back->next;
  while (p->next != outgroup->back)
    p = p->next;
  root->next = outgroup->back;
  p->next = root;
}  /* reroot2 */


void reroot3(node *outgroup,node *root,node *root2,node *lastdesc,node **grbg)
{
  /* reorients tree, putting back outgroup in original position. */
  /* used in pars */
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


void bintomulti(node **root, node **binroot, node **grbg, long *zeros,
                        unsigned char *zeros2)
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
  gnudisctreenode(grbg, &newnode, right->index, endsite, zeros, zeros2);
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


void savetree(node *root, long *place, pointarray treenode,
                        node **grbg, long *zeros, unsigned char *zeros2)
{  /* record in place where each species has to be
     added to reconstruct this tree */
  /* used by pars */
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
    bintomulti(&root, &binroot, grbg, zeros, zeros2);
  if (outgrin(root, outgrnode)) {
    if (outgrnode != root->next->back)
      moveleft(root, outgrnode, &flipback);
  } else {
    root2 = root;
    lastdesc = root->next;
    while (lastdesc->next != root)
      lastdesc = lastdesc->next;
    lastdesc->next = root->next;
    gnudisctreenode(grbg, &root, outgrnode->back->index, endsite, zeros, zeros2);
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
                        boolean multf, pointarray treenode, long *place, long *zeros,
                        unsigned char *zeros2)
{  /* adds item to tree and save it.  Then removes item. */
  node *dummy;

  if (!multf)
    add(p, item, nufork, root, false, treenode, grbg, zeros, zeros2);
  else
    add(p, item, NULL, root, false, treenode, grbg, zeros, zeros2);
  savetree(*root, place, treenode, grbg, zeros, zeros2);
  if (!multf)
    re_move(item, &nufork, root, false, treenode, grbg, zeros, zeros2);
  else
    re_move(item, &dummy, root, false, treenode, grbg, zeros, zeros2);
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
  long i, j, largest, descsteps;
  boolean done;
  unsigned char b;

  for (i = 0; i < endsite; i++) {
    descsteps = 0;
    for (j = (long)zero; j <= (long)seven; j++) {
      b = 1 << j;
      if ((descsteps == 0) && (collapfrom->discbase[i] & b)) 
        descsteps = tempfrom->oldnumsteps[i] 
                     - (collapfrom->numdesc - collapfrom->discnumnuc[i][j])
                       * weight[i];
    }
    done = false;
    for (j = (long)zero; j <= (long)seven; j++) {
      b = 1 << j;
      if (!done && (tempto->discbase[i] & b)) {
        descsteps += (tempto->numsteps[i] 
                      - (tempto->numdesc - collapfrom->numdesc
                        - tempto->discnumnuc[i][j]) * weight[i]);
        done = true;
      }
    }
    for (j = (long)zero; j <= (long)seven; j++)
      tempto->discnumnuc[i][j] += collapfrom->discnumnuc[i][j];
    largest = getlargest(tempto->discnumnuc[i]);
    tempto->discbase[i] = 0;
    for (j = (long)zero; j <= (long)seven; j++) {
      if (tempto->discnumnuc[i][j] == largest)
        tempto->discbase[i] |= (1 << j);
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
    if ((a->discbase[i] & b->discbase[i]) == 0)
      allcommon = false;
    else if (a->discbase[i] != b->discbase[i])
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


boolean passdown(node *desc, node *parent, node *start, node *below, node *item,
                        node *added, node *total, node *tempdsc, node *tempprt,
                        boolean multf)
{ /* track down to node start to see if an ancestor branch can be collapsed */
  node *temp;
  boolean done, allsame;

  done = (parent == start);
  while (!done) {
    desc = parent;
    findbottom(parent->back, &parent);
    if (multf && start == below && parent == below)
      parent = added;
    memcpy(tempdsc->discbase, tempprt->discbase, endsite*sizeof(unsigned char));
    memcpy(tempdsc->numsteps, tempprt->numsteps, endsite*sizeof(long));
    memcpy(tempdsc->olddiscbase, desc->discbase, endsite*sizeof(unsigned char));
    memcpy(tempdsc->oldnumsteps, desc->numsteps, endsite*sizeof(long));
    memcpy(tempprt->discbase, parent->discbase, endsite*sizeof(unsigned char));
    memcpy(tempprt->numsteps, parent->numsteps, endsite*sizeof(long));
    memcpy(tempprt->discnumnuc, parent->discnumnuc, endsite*sizeof(discnucarray));
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
      memcpy(tempdsc->discbase, tempprt->discbase, endsite*sizeof(unsigned char));
      memcpy(tempdsc->numsteps, tempprt->numsteps, endsite*sizeof(long));
      memcpy(tempdsc->olddiscbase, start->discbase, endsite*sizeof(unsigned char));
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


boolean trycollapdesc(node *desc, node *parent, node *start, node *below,
                        node *item, node *added, node *total, node *tempdsc, node *tempprt,
                        boolean multf,long *zeros, unsigned char *zeros2)
{ /* see if branch between nodes desc and parent can be collapsed */
  boolean allsame;

  if (desc->numdesc == 1)
    return true;
  if (multf && start == below && parent == below)
    parent = added;
  memcpy(tempdsc->discbase, zeros2, endsite*sizeof(unsigned char));
  memcpy(tempdsc->numsteps, zeros, endsite*sizeof(long));
  memcpy(tempdsc->olddiscbase, desc->discbase, endsite*sizeof(unsigned char));
  memcpy(tempdsc->oldnumsteps, desc->numsteps, endsite*sizeof(long));
  memcpy(tempprt->discbase, parent->discbase, endsite*sizeof(unsigned char));
  memcpy(tempprt->numsteps, parent->numsteps, endsite*sizeof(long));
  memcpy(tempprt->discnumnuc, parent->discnumnuc, endsite*sizeof(discnucarray));
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
    memcpy(tempdsc->discbase, tempprt->discbase, endsite*sizeof(unsigned char));
    memcpy(tempdsc->numsteps, tempprt->numsteps, endsite*sizeof(long));
    memcpy(tempdsc->olddiscbase, start->discbase, endsite*sizeof(unsigned char));
    memcpy(tempdsc->oldnumsteps, start->numsteps, endsite*sizeof(long));
    memcpy(tempprt->discbase, added->discbase, endsite*sizeof(unsigned char));
    memcpy(tempprt->numsteps, added->numsteps, endsite*sizeof(long));
    memcpy(tempprt->discnumnuc, added->discnumnuc, endsite*sizeof(discnucarray));
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
                        boolean multf, node* root, long *zeros, 
                        unsigned char *zeros2)
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
                tempdsc, tempprt, multf, zeros, zeros2))
          return true;
        else if ((p->back->index == root->index && root->numdesc == 2) &&
            !(root->next->back->tip) && !(root->next->next->back->tip) &&
            trycollapdesc(root->next->back, root->next->next->back, start,
                below, item, added, total, tempdsc, tempprt, multf, zeros,
                zeros2))
          return true;
      }
      p = p->next;
    } while (p != subtree);
    p = subtree->next;
    do {
      if (p->back && !p->back->tip) {
        if (zeroinsubtree(p->back, start, below, item, added, total, 
                            tempdsc, tempprt, multf, root, zeros, zeros2))
          return true;
      }
      p = p->next;
    } while (p != subtree);
  }
  return false;
} /* zeroinsubtree */


boolean collapsible(node *item, node *below, node *temp, node *temp1, 
                node *tempdsc, node *tempprt, node *added, node *total, 
                boolean multf, node *root, long *zeros, unsigned char *zeros2, 
                pointarray treenode)
{
  /* sees if any branch can be collapsed */
  node *belowbk;
  boolean allsame;

  if (multf) {
    memcpy(tempdsc->discbase, item->discbase, endsite*sizeof(unsigned char));
    memcpy(tempdsc->numsteps, item->numsteps, endsite*sizeof(long));
    memcpy(tempdsc->olddiscbase, zeros2, endsite*sizeof(unsigned char));
    memcpy(tempdsc->oldnumsteps, zeros, endsite*sizeof(long));
    memcpy(added->discbase, below->discbase, endsite*sizeof(unsigned char));
    memcpy(added->numsteps, below->numsteps, endsite*sizeof(long));
    memcpy(added->discnumnuc, below->discnumnuc, endsite*sizeof(discnucarray));
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
                        tempdsc, tempprt, multf, root, zeros, zeros2))
      return true;
  }
  if (multf) {
    if (zeroinsubtree(below, below, below, item, added, total,
                       tempdsc, tempprt, multf, root, zeros, zeros2))
      return true;
  } else if (!below->tip) {
    if (zeroinsubtree(below, below, below, item, added, total,
                        tempdsc, tempprt, multf, root, zeros, zeros2))
      return true;
  }
  if (!item->tip) {
    if (zeroinsubtree(item, item, below, item, added, total,
                        tempdsc, tempprt, multf, root, zeros, zeros2))
      return true;
  }
  if (multf && below->back && !below->back->tip) {
    memcpy(tempdsc->discbase, zeros2, endsite*sizeof(unsigned char));
    memcpy(tempdsc->numsteps, zeros, endsite*sizeof(long));
    memcpy(tempdsc->olddiscbase, added->discbase, endsite*sizeof(unsigned char));
    memcpy(tempdsc->oldnumsteps, added->numsteps, endsite*sizeof(long));
    if (below->back == treenode[below->back->index - 1])
      belowbk = below->back->next;
    else
      belowbk = treenode[below->back->index - 1];
    memcpy(tempprt->discbase, belowbk->discbase, endsite*sizeof(unsigned char));
    memcpy(tempprt->numsteps, belowbk->numsteps, endsite*sizeof(long));
    memcpy(tempprt->discnumnuc, belowbk->discnumnuc, endsite*sizeof(discnucarray));
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


void replaceback(node **oldback, node *item, node *forknode, node **grbg,
                        long *zeros, unsigned char *zeros2)
{ /* replaces back node of item with another */
  node *p;

  p = forknode;
  while (p->next->back != item)
    p = p->next;
  *oldback = p->next;
  gnudisctreenode(grbg, &p->next, forknode->index, endsite, zeros, zeros2);
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


void savelocrearr(node *item, node *forknode, node *below, node *tmp, node *tmp1,
                        node *tmp2, node *tmp3, node *tmprm, node *tmpadd, node **root,
                        long maxtrees, long *nextree, boolean multf, boolean bestever,
                        boolean *saved, long *place, bestelm *bestrees, pointarray treenode,
                        node **grbg, long *zeros, unsigned char *zeros2)
{ /* saves tied or better trees during local rearrangements by removing
     item from forknode and adding to below */
  node *other, *otherback=NULL, *oldfork, *nufork, *oldback;
  long pos;
  boolean found, collapse;

  if (forknode->numdesc == 2) {
    findbelow(&other, item, forknode);
    otherback = other->back;
    oldback = NULL;
  } else {
    other = NULL;
    replaceback(&oldback, item, forknode, grbg, zeros, zeros2);
  }
  re_move(item, &oldfork, root, false, treenode, grbg, zeros, zeros2);
  if (!multf)
    getnufork(&nufork, grbg, treenode, zeros, zeros2);
  else
    nufork = NULL;
  addnsave(below, item, nufork, root, grbg, multf, treenode, place,
             zeros, zeros2);
  pos = 0;
  findtree(&found, &pos, *nextree, place, bestrees);
  if (other) {
    add(other, item, oldfork, root, false, treenode, grbg, zeros, zeros2);
    if (otherback->back != other)
      flipnodes(item, other);
  } else
    add(forknode, item, NULL, root, false, treenode, grbg, zeros, zeros2);
  *saved = false;
  if (found) {
    if (oldback)
      putback(oldback, item, forknode, grbg);
  } else {
    if (oldback)
      chuck(grbg, oldback);
    re_move(item, &oldfork, root, true, treenode, grbg, zeros, zeros2);
    collapse = collapsible(item, below, tmp, tmp1, tmp2, tmp3, tmprm,
                     tmpadd, multf, *root, zeros, zeros2, treenode);
    if (!collapse) {
      if (bestever)
        addbestever(&pos, nextree, maxtrees, collapse, place, bestrees);
      else
        addtiedtree(pos, nextree, maxtrees, collapse, place, bestrees);
    }
    if (other)
      add(other, item, oldfork, root, true, treenode, grbg, zeros, zeros2);
    else
      add(forknode, item, NULL, root, true, treenode, grbg, zeros, zeros2);
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


void hyprint(long b1,long b2,struct LOC_hyptrav *htrav,pointarray treenode)
{
  /* print out states in sites b1 through b2 at node */
  long i, j, k;
  boolean dot, found;

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
    htrav->tempset = htrav->r->discbase[j - 1];
    htrav->anc = htrav->hypset[j - 1];
    if (!htrav->bottom)
      htrav->anc = treenode[htrav->r->back->index - 1]->discbase[j - 1];
    dot = dotdiff && (htrav->tempset == htrav->anc && !htrav->bottom);
    if (dot)
      putc('.', outfile); 
    else {
      found = false;
      k = (long)zero;
      do {
        if (htrav->tempset == (1 << k)) {
          putc(convtab[k][i - 1], outfile);
          found = true;
        }
        k++;
      } while (!found && k <= (long)seven);
      if (!found) 
        putc('?', outfile);
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
    (*p)->discbase = (discbaseptr)Malloc(endsite*sizeof(unsigned char));
  }
  (*p)->next = NULL;
}  /* gnubase */


void chuckbase(gbases *p, gbases **garbage)
{
  /* collect garbage on p -- put it on front of garbage list */
  p->next = *garbage;
  *garbage = p;
}  /* chuckbase */


void hyptrav(node *r_, discbaseptr hypset_, long b1, long b2, boolean bottom_,
                        pointarray treenode, gbases **garbage)
{
  /*  compute, print out states at one interior node */
  struct LOC_hyptrav Vars;
  long i, j, k;
  long largest;
  gbases *ancset;
  discnucarray *tempnuc;
  node *p, *q;

  Vars.bottom = bottom_;
  Vars.r = r_;
  Vars.hypset = hypset_;
  gnubase(&ancset, garbage, endsite);
  tempnuc = (discnucarray *)Malloc(endsite*sizeof(discnucarray));
  Vars.maybe = false;
  Vars.nonzero = false;
  if (!Vars.r->tip)
    zerodiscnumnuc(Vars.r, endsite);
  for (i = b1 - 1; i < b2; i++) {
    j = location[ally[i] - 1];
    Vars.anc = Vars.hypset[j - 1];
    if (!Vars.r->tip) {
      p = Vars.r->next;
      for (k = (long)zero; k <= (long)seven; k++)
        if (Vars.anc & (1 << k))
          Vars.r->discnumnuc[j - 1][k]++;
      do {
        for (k = (long)zero; k <= (long)seven; k++)
          if (p->back->discbase[j - 1] & (1 << k))
            Vars.r->discnumnuc[j - 1][k]++;
        p = p->next;
      } while (p != Vars.r);
      largest = getlargest(Vars.r->discnumnuc[j - 1]);
      Vars.tempset = 0;
      for (k = (long)zero; k <= (long)seven; k++) {
        if (Vars.r->discnumnuc[j - 1][k] == largest)
          Vars.tempset |= (1 << k);
      }
      Vars.r->discbase[j - 1] = Vars.tempset;
    }
    if (!Vars.bottom)
      Vars.anc = treenode[Vars.r->back->index - 1]->discbase[j - 1];
    Vars.nonzero = (Vars.nonzero || (Vars.r->discbase[j - 1] & Vars.anc) == 0);
    Vars.maybe = (Vars.maybe || Vars.r->discbase[j - 1] != Vars.anc);
  }
  hyprint(b1, b2, &Vars, treenode);
  Vars.bottom = false;
  if (!Vars.r->tip) {
    memcpy(tempnuc, Vars.r->discnumnuc, endsite*sizeof(discnucarray));
    q = Vars.r->next;
    do {
      memcpy(Vars.r->discnumnuc, tempnuc, endsite*sizeof(discnucarray));
      for (i = b1 - 1; i < b2; i++) {
        j = location[ally[i] - 1];
        for (k = (long)zero; k <= (long)seven; k++)
          if (q->back->discbase[j - 1] & (1 << k))
            Vars.r->discnumnuc[j - 1][k]--;
        largest = getlargest(Vars.r->discnumnuc[j - 1]);
        ancset->discbase[j - 1] = 0;
        for (k = (long)zero; k <= (long)seven; k++)
          if (Vars.r->discnumnuc[j - 1][k] == largest)
            ancset->discbase[j - 1] |= (1 << k);
        if (!Vars.bottom)
          Vars.anc = ancset->discbase[j - 1];
      }
      hyptrav(q->back, ancset->discbase, b1, b2, Vars.bottom,
                treenode, garbage);
      q = q->next;
    } while (q != Vars.r);
  }
  chuckbase(ancset, garbage);
}  /* hyptrav */


void hypstates(long chars, node *root, pointarray treenode, gbases **garbage)
{
  /* fill in and describe states at interior nodes */
  /* used in pars */
  long i, n;
  discbaseptr nothing;

  fprintf(outfile, "\nFrom    To     Any Steps?    State at upper node\n");
  fprintf(outfile, "                            ");
  if (dotdiff)
    fprintf(outfile, " ( . means same as in the node below it on tree)\n");
  nothing = (discbaseptr)Malloc(endsite*sizeof(unsigned char));
  for (i = 0; i < endsite; i++)
    nothing[i] = 0;
  for (i = 1; i <= ((chars - 1) / 40 + 1); i++) {
    putc('\n', outfile);
    n = i * 40;
    if (n > chars)
      n = chars;
    hyptrav(root, nothing, i * 40 - 39, n, true, treenode, garbage);
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
    for (i = (long)zero; i <= (long)seven; i++) {
      p->disccumlengths[i] = 0;
      p->discnumreconst[i] = 1;
    }
  } else {
    for (i = (long)zero; i <= (long)seven; i++) {
      if (p->discbase[sitei - 1] & (1 << i)) {
        p->disccumlengths[i] = 0;
        p->discnumreconst[i] = 1;
      } else {
        p->disccumlengths[i] = -1;
        p->discnumreconst[i] = 0;
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
      memcpy(q->discnumnuc, p->discnumnuc, endsite*sizeof(discnucarray));
      for (i = (long)zero; i <= (long)seven; i++) {
        if (q->back->discbase[sitei - 1] & (1 << i))
          q->discnumnuc[sitei - 1][i]--;
      }
      if (p->back) {
        for (i = (long)zero; i <= (long)seven; i++) {
          if (p->back->discbase[sitei - 1] & (1 << i))
            q->discnumnuc[sitei - 1][i]++;
        }
      }
      largest = getlargest(q->discnumnuc[sitei - 1]);
      q->discbase[sitei - 1] = 0;
      for (i = (long)zero; i <= (long)seven; i++) {
        if (q->discnumnuc[sitei - 1][i] == largest)
          q->discbase[sitei - 1] |= (1 << i);
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
  for (i = (long)zero; i <= (long)seven; i++) {
    minn = maxx;
    for (j = (long)zero; j <= (long)seven; j++) {
      if (i == j) 
        cost = 0;
      else
        cost = 1;
      if (desc->disccumlengths[j] == -1) {
        desclen = maxx;
      } else {
        desclen = desc->disccumlengths[j];
      }
      if (minn > cost + desclen) {
        minn = cost + desclen;
        descrecon = 0;
      }
      if (minn == cost + desclen) {
        descrecon += desc->discnumreconst[j];
      }
    }
    p->disccumlengths[i] += minn;
    p->discnumreconst[i] *= descrecon;
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
  for (i = (long)zero; i <= (long)seven; i++) {
    for (j = (long)zero; j <= (long)seven; j++) {
      if (i == j)
        cost = 0;
      else
        cost = 1;
      if (subtr1->disccumlengths[i] != -1 && (subtr2->disccumlengths[j] != -1)) {
        if (subtr1->disccumlengths[i] + cost + subtr2->disccumlengths[j] < minn) {
          minn = subtr1->disccumlengths[i] + cost + subtr2->disccumlengths[j];
          nom = 0;
          denom = 0;
        }
        if (subtr1->disccumlengths[i] + cost + subtr2->disccumlengths[j] == minn) {
          nom += subtr1->discnumreconst[i] * subtr2->discnumreconst[j] * cost;
          denom += subtr1->discnumreconst[i] * subtr2->discnumreconst[j];
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
    fprintf(outfile, "     %.2f\n",q->v);
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
      q->v += ((weight[sitei - 1] / 10.0) * (*brlen));
      q->back->v += ((weight[sitei - 1] / 10.0) * (*brlen));
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
  node *q, *first, *last, *mid1 = NULL, *mid2 = NULL;
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
    if (numb2 == (numbranches + 1)/2)
      mid1 = q->back;
    if (numb2 == (numbranches/2 + 1))
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
  node *p, *q, *r, *first = NULL, *last = NULL;
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
  /* used in pars */
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
  /* used in pars */
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
  /* used in pars */
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
    fprintf(outtree, ":%*.2f", (int)(w + 4), x);
  }
  if (p != root)
    return;
  if (nextree > 2)
    fprintf(outtree, "[%6.4f];\n", 1.0 / (nextree - 1));
  else
    fprintf(outtree, ";\n");
}  /* treeout3 */


void drawline3(long i, double scale, node *start)
{
  /* draws one row of the tree diagram by moving up tree */
  /* used in pars */
  node *p, *q;
  long n, j;
  boolean extra;
  node *r, *first = NULL, *last = NULL;
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


void standev(long chars, long numtrees, long minwhich, double minsteps,
                        double *nsteps, long **fsteps, longer seed)
{  /* do paired sites test (KHT or SH) on user trees */
   /* used in pars */
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
            sum += temp;
            sum2 += temp * temp / wt;
          }
        }
        temp = sum / sumw;
        sd = sqrt(sumw / (sumw - 1.0) * (sum2 - sum * sum / sumw));
        fprintf(outfile, "%9.1f %12.4f",
                (nsteps[which - 1] - minsteps) / 10, sd);
        if (sum > 1.95996 * sd)
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
    for (i = 0; i < numtrees; i++)
      covar[i] = (double *)Malloc(numtrees*sizeof(double));  
    sumw = 0.0;
    for (i = 0; i < endsite; i++)
      sumw += weight[i];
    for (i = 0; i < numtrees; i++) {        /* compute covariances of trees */
      sum = nsteps[i]/sumw;
      for (j = 0; j <=i; j++) {
        sum2 = nsteps[j]/sumw;
        temp = 0.0;
        for (k = 0; k < endsite; k++) {
          wt = weight[k]/10.0;
          if (weight[k] > 0) {
            temp = temp + wt*(fsteps[i][k]/(wt*10.0)-sum)
                            *(fsteps[j][k]/(wt*10.0)-sum2);
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
        if (fabs(temp) < 1.0E-23)
          covar[j][i] = 0.0;
        else
          covar[j][i] = (covar[j][i] - sum)/temp;
      }
    }
    f = (double *)Malloc(numtrees*sizeof(double)); /* resampled sum of steps */
    P = (double *)Malloc(numtrees*sizeof(double)); /* vector of P's of trees */
    r = (double *)Malloc(numtrees*sizeof(double)); /* store normal variates */
    for (i = 0; i < numtrees; i++)
      P[i] = 0.0;
    sum2 = nsteps[0]/10.0;               /* sum2 will be smallest # of steps */
    for (i = 1; i < numtrees; i++)
      if (sum2 > nsteps[i]/10.0)
        sum2 = nsteps[i]/10.0;
    for (i = 1; i <= SAMPLES; i++) {          /* loop over resampled trees */
      for (j = 0; j < numtrees; j++)          /* draw normal variates */
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
        fprintf(outfile, "  %9.1f %10.3f", nsteps[i]/10.0-sum2, P[i]);
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
  /* used in pars */

  free(anode->numsteps);
  free(anode->oldnumsteps);
  free(anode->discbase);
  free(anode->olddiscbase);
}  /* freetip */


void freenontip(node *anode)
{
  /* used in pars */
  free(anode->numsteps);
  free(anode->oldnumsteps);
  free(anode->discbase);
  free(anode->olddiscbase);
  free(anode->discnumnuc);
}  /* freenontip */


void freenodes(long nonodes, pointarray treenode)
{
  /* used in pars */
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
  /* used in pars */

  freenontip(*anode);
  free(*anode);
}  /* freenode */


void freetree(long nonodes, pointarray treenode)
{
  /* used in pars */
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


void freegarbage(gbases **garbage)
{
  /* used in pars */
  gbases *p;

  while (*garbage) {
    p = *garbage;
    *garbage = (*garbage)->next;
    free(p->discbase);
    free(p);
  }
}  /* freegarbage */


void freegrbg(node **grbg)
{
  /* used in pars */
  node *p;

  while (*grbg) {
    p = *grbg;
    *grbg = (*grbg)->next;
    freenontip(p);
    free(p);
  }
} /*freegrbg */


void collapsetree(node *p, node *root, node **grbg, pointarray treenode, 
                  long *zeros, unsigned char *zeros2)
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
        gnudisctreenode(grbg, &x2, index2, endsite, zeros, zeros2);
        x2->next = x1;
        x1 = x2;
      }
      x2->next->next->next = x2;
      treenode[nonodes-1]=x2;
      if (q->back)
        collapsetree(q->back, root, grbg, treenode, zeros, zeros2);
    } else {
      if (q->back)
        collapsetree(q->back, root, grbg, treenode, zeros, zeros2);
      q = q->next;
    }
  } while (q != p);
} /* collapsetree */


void collapsebestrees(node **root, node **grbg, pointarray treenode, 
                      bestelm *bestrees, long *place, long *zeros, 
                      unsigned char *zeros2, long chars, boolean recompute, 
                      boolean progress)
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
        treenode, grbg, zeros, zeros2);
    nextnode = spp + 2;
    for (j = 3; j <= spp; j++) {
      if (bestrees[k].btree[j - 1] > 0)
        add(treenode[bestrees[k].btree[j - 1] - 1], treenode[j - 1],
            treenode[nextnode++ - 1], root, recompute, treenode, grbg,
            zeros, zeros2);
      else
          add(treenode[treenode[-bestrees[k].btree[j - 1]-1]->back->index-1],
              treenode[j - 1], NULL, root, recompute, treenode, grbg, zeros, zeros2);
    }
    reroot(treenode[outgrno - 1], *root);

    treelength(*root, chars, treenode);
    collapsetree(*root, *root, grbg, treenode, zeros, zeros2);
    savetree(*root, place, treenode, grbg, zeros, zeros2);
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
              grbg, zeros, zeros2);
    }
  }
  if (progress) {
    putchar('\n');
#ifdef WIN32
    phyFillScreenColor();
#endif
  }
}
