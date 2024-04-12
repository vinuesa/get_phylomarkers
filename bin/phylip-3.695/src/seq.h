/* version 3.6. (c) Copyright 1993-2000 by the University of Washington.
   Written by Joseph Felsenstein, Akiko Fuseki, Sean Lamont, and Andrew Keeffe.
   Permission is granted to copy and use this program provided no fee is
   charged for it and provided that this copyright notice is not removed. */

/*
    seq.h:  included in dnacomp, dnadist, dnainvar, dnaml, dnamlk, dnamove,
            dnapars, dnapenny, protdist, protpars, restdist & restml
*/

#ifndef SEQ_H
#define SEQ_H

#define ebcdic          EBCDIC
#define MAXNCH          20

/* All of this came over from cons.h    -plc*/ 
#define OVER              7
#define ADJACENT_PAIRS    1
#define CORR_IN_1_AND_2   2
#define ALL_IN_1_AND_2    3
#define NO_PAIRING        4
#define ALL_IN_FIRST      5
#define TREE1             8
#define TREE2             9

#define FULL_MATRIX       11
#define VERBOSE           22
#define SPARSE            33

/* Number of columns per block in a matrix output */
#define COLUMNS_PER_BLOCK 10


/*end move*/


typedef struct gbases {
  baseptr base;
  struct gbases *next;
} gbases;

typedef struct nuview_data {
  /* A big 'ol collection of pointers used in nuview */
  double *yy, *wwzz, *vvzz, *vzsumr, *vzsumy, *sum, *sumr, *sumy;
  sitelike *xx;
} nuview_data;

struct LOC_hyptrav {
  boolean bottom;
  node *r;
  long *hypset;
  boolean maybe, nonzero;
  long tempset, anc;
} ;


extern long nonodes, endsite, outgrno, nextree, which;

extern boolean interleaved, printdata, outgropt, treeprint, dotdiff, transvp;
extern steptr weight, category, alias, location, ally;
extern sequence y;

#ifndef OLDC
/* function prototypes */
void   alloctemp(node **, long *, long);
void   freetemp(node **);
void   freetree2 (pointarray, long);
void   inputdata(long);
void   alloctree(pointarray *, long, boolean);
void   allocx(long, long, pointarray, boolean);

void   prot_allocx(long, long, pointarray, boolean);
void   setuptree(pointarray, long, boolean);
void   setuptree2(tree *);
void   alloctip(node *, long *);
void   getbasefreqs(double, double, double, double, double *, double *,
                        double *, double *, double *, double *, double *,
            double *xi, double *, double *, boolean, boolean);
void   empiricalfreqs(double *,double *,double *,double *,steptr,pointarray);
void   sitesort(long, steptr);
void   sitecombine(long);

void   sitescrunch(long);
void   sitesort2(long, steptr);
void   sitecombine2(long, steptr);
void   sitescrunch2(long, long, long, steptr);
void   makevalues(pointarray, long *, boolean);
void   makevalues2(long, pointarray, long, long, sequence, steptr);
void   fillin(node *, node *, node *);
long   getlargest(long *);
void   multifillin(node *, node *, long);
void   sumnsteps(node *, node *, node *, long, long);

void   sumnsteps2(node *, node *, node *, long, long, long *);
void   multisumnsteps(node *, node *, long, long, long *);
void   multisumnsteps2(node *);
boolean alltips(node *, node *);
void   gdispose(node *, node **, pointarray);
void   preorder(node *, node *, node *, node *, node *, node *, long);
void   updatenumdesc(node *, node *, long);
void   add(node *,node *,node *,node **,boolean,pointarray,node **,long *);
void   findbelow(node **below, node *item, node *fork);

void   re_move(node *item, node **fork, node **root, boolean recompute,
                pointarray, node **, long *);
void   postorder(node *p);
void   getnufork(node **, node **, pointarray, long *);
void   reroot(node *, node *);
void   reroot2(node *, node *);
void   reroot3(node *, node *, node *, node *, node **);
void   savetraverse(node *);
void   newindex(long, node *);
void   flipindexes(long, pointarray);
boolean parentinmulti(node *);

long   sibsvisited(node *, long *);
long   smallest(node *, long *);
void   bintomulti(node **, node **, node **, long *);
void   backtobinary(node **, node *, node **);
boolean outgrin(node *, node *);
void   flipnodes(node *, node *);
void   moveleft(node *, node *, node **);
void   savetree(node *, long *, pointarray, node **, long *);
void   addnsave(node *, node *, node *, node **, node **,boolean,
                pointarray, long *, long *);
void   addbestever(long *, long *, long, boolean, long *, bestelm *);

void   addtiedtree(long, long *, long, boolean,long *, bestelm *);
void   clearcollapse(pointarray);
void   clearbottom(pointarray);
void   collabranch(node *, node *, node *);
boolean allcommonbases(node *, node *, boolean *);
void   findbottom(node *, node **);
boolean moresteps(node *, node *);
boolean passdown(node *, node *, node *, node *, node *, node *,
                node *, node *, node *, boolean);
boolean trycollapdesc(node *, node *, node *, node *, node *,
                node *, node *, node *, node *, boolean , long *);
void   setbottom(node *);

boolean zeroinsubtree(node *, node *, node *, node *, node *,
                node *, node *, node *, boolean, node *, long *);
boolean collapsible(node *, node *, node *, node *, node *,
                node *, node *, node *, boolean, node *, long *, pointarray);
void   replaceback(node **, node *, node *, node **, long *);
void   putback(node *, node *, node *, node **);
void   savelocrearr(node *, node *, node *, node *, node *, node *,
                node *, node *, node *, node **, long, long *, boolean,
                boolean , boolean *, long *, bestelm *, pointarray ,
                node **, long *);
void   clearvisited(pointarray);
void   hyprint(long, long, struct LOC_hyptrav *,pointarray, Char *);
void   gnubase(gbases **, gbases **, long);
void   chuckbase(gbases *, gbases **);
void   hyptrav(node *, long *, long, long, boolean,pointarray,
                gbases **, Char *);

void   hypstates(long , node *, pointarray, gbases **, Char *);
void   initbranchlen(node *p);
void   initmin(node *, long, boolean);
void   initbase(node *, long);
void   inittreetrav(node *, long);
void   compmin(node *, node *);
void   minpostorder(node *, pointarray);
void   branchlength(node *,node *,double *,pointarray);
void   printbranchlengths(node *);
void   branchlentrav(node *,node *,long,long,double *,pointarray);

void   treelength(node *, long, pointarray);
void   coordinates(node *, long *, double, long *);
void   drawline(long, double, node *);
void   printree(node *, double);
void   writesteps(long, boolean, steptr, node *);
void   treeout(node *, long, long *, node *);
void   treeout3(node *, long, long *, node *);
void   fdrawline2(FILE *fp, long i, double scale, tree *curtree);
void   drawline2(long, double, tree);
void   drawline3(long, double, node *);
void   copynode(node *, node *, long);

void   prot_copynode(node *, node *, long);
void   copy_(tree *, tree *, long, long);
void   prot_copy_(tree *, tree *, long, long);
void   standev(long, long, long, double, double *, long **, longer);
void   standev2(long, long, long, long, double, double *, double **,
              steptr, longer);
void   freetip(node *);
void   freenontip(node *);
void   freenodes(long, pointarray);
void   freenode(node **);
void   freetree(long, pointarray);

void   freex(long, pointarray);
void   freex_notip(long, pointarray);
void   prot_freex_notip(long nonodes, pointarray treenode);
void   prot_freex(long nonodes, pointarray treenode);
void   freegarbage(gbases **);
void   freegrbg(node **);

void   collapsetree(node *, node *, node **, pointarray, long *);
void   collapsebestrees(node **, node **, pointarray, bestelm *, long *,
                      long *, long, boolean, boolean);
void   fix_x(node* p,long site, double maxx, long rcategs);
void   fix_protx(node* p,long site,double maxx, long rcategs);
/*function prototypes*/
#endif

#endif /* SEQ_H */
