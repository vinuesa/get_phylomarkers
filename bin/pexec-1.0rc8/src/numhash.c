/*****************************************************************************/
/* numhash.c								     */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/* Hash table for numerical (integer) keys				     */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/* Copyright (C) 2007; Pal, A. (apal@szofi.elte.hu)			     */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*  This library is free software: you can redistribute it and/or modify     */
/*  it under the terms of the GNU General Public License as published by     */
/*  the Free Software Foundation, either version 3 of the License, or	     */
/*  (at your option) any later version.					     */
/*									     */
/*  This program is distributed in the hope that it will be useful,	     */
/*  but WITHOUT ANY WARRANTY; without even the implied warranty of	     */
/*  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the	     */
/*  GNU General Public License for more details.			     */
/*									     */
/*  You should have received a copy of the GNU General Public License	     */
/*  along with the program.  If not, see <http://www.gnu.org/licenses/>.     */
/*****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "numhash.h"

#ifdef HAVE_LIBINTL_H
#include <libintl.h>
#define _(string)	gettext(string)
#else
#define _(string)	string
#endif

#define	NUMHASH_ABORT_ON_MALLOC
#ifdef	NUMHASH_ABORT_ON_MALLOC
#define	realloc_check(ptr,size) \
	do \
	 {	if ( ptr==NULL && size>0 ) \
		 {	fprintf(stderr,"numhash.c: %s.\n",_("memory exhausted"));\
			abort(); \
		 } \
	 } while(0)
#else
#define	realloc_check(ptr,size) 
#endif
#define	malloc_check(ptr)	realloc_check(ptr,1)

/*****************************************************************************/

#define	NUMHASH_BIT(bitsize,depth,off,value) \
	(((value)>>((depth)*((((bitsize)-1)/(depth))-(off))))&((1<<(depth))-1))

#define	NUMHASH_IOF(bitsize,depth,off) \
	((off)<((bitsize)+(depth)-1)/(depth))

#define	NUMHASH_LST(bitsize,depth,off) \
	((off)==((bitsize)-1)/(depth))

#define	NUMHASH_RLF(bitsize,depth,off) \
	((off)<((bitsize)-1)/(depth))

#define	NUMHASH_NTN(bitsize,depth,off) \
	(1<<((off)?(depth):((bitsize)-(depth)*(((bitsize)-1)/(depth)))))

#define	NUMHASH_LCB(bitsize,depth,off) \
	(1<<((depth)*(((bitsize)-1)/(depth)-(off))))

/*****************************************************************************/

int numhash_search(numhashtable *nt,int value,void **ret)
{
 int		off,bit;
 numhashnode	*nh;

 if ( nt==NULL )	return(0);

 off=0;
 nh=&nt->table;

 if ( nt->bitsize<sizeof(int)*8 ) 
	value=value & ((1<<nt->bitsize)-1);

 do
  {	/*shift=depth*(((bitsize-1)/depth)-off);
	bit=(value >> shift) & ((1<<depth)-1);*/
	bit=NUMHASH_BIT(nt->bitsize,nt->depth,off,value);

	if ( nh->node.leaves==NULL )
		return(0);

	off++;
	nh=&nh->node.leaves[bit];

  } while ( off*nt->depth < nt->bitsize );

 if ( ret != NULL )	*ret=nh->node.data;
 return(nh->nchild);
}

/*****************************************************************************/

int numhash_init(numhashtable *nt,int bitsize,int depth)
{
 nt->table.node.leaves=NULL;
 nt->table.node.data=NULL;
 nt->table.nchild=0;

 nt->depth=depth;
 nt->bitsize=bitsize;

 return(0);
}

/*****************************************************************************/

int numhash_add(numhashtable *nt,int value,void *data)
{
 int		off,bit,bitsize,depth,shift,tnode,nodesize;
 numhashnode	*nh,*tree[64];

 if ( nt==NULL )
	return(-1);

 off=0;
 nh=&nt->table;
 depth=nt->depth;
 bitsize=nt->bitsize;

 tnode=0;

 if ( bitsize<sizeof(int)*8 )
	value=value & ((1<<bitsize)-1);

 do
  {	shift=depth*(((bitsize-1)/depth)-off);
	bit=(value >> shift) & ((1<<depth)-1);

	if ( nh->node.leaves==NULL )
	 {	nodesize=(1<<depth);
		nh->node.leaves=malloc(sizeof(numhashnode)*nodesize);
		malloc_check(nh->node.leaves);
		memset(nh->node.leaves,0,sizeof(numhashnode)*nodesize);
	 }
	off++;
	tree[tnode]=nh;
	tnode++;
	nh=&nh->node.leaves[bit];

  } while ( off*depth < bitsize );

 nh->node.data=data;

 if ( nh->nchild )
 	return(0);
 else
  {	nh->nchild=1;
	while ( tnode>0 )
	 {	tnode--;
		(tree[tnode]->nchild)++;
	 };
	return(1);
  }

}

/*****************************************************************************/

int numhash_total(numhashtable *nt)
{
 return(nt->table.nchild);
}

/*****************************************************************************/

int numhash_remove(numhashtable *nt,int value)
{
 int		off,bit,bitsize,depth,shift,tnode;
 numhashnode	*nh,*tree[64];

 if ( nt==NULL )	return(0);

 off=0;
 nh=&nt->table;
 depth=nt->depth;
 bitsize=nt->bitsize;

 tnode=0;

 if ( bitsize<sizeof(int)*8 )
	value=value & ((1<<bitsize)-1);

 do
  {	shift=depth*(((bitsize-1)/depth)-off);
	bit=(value >> shift) & ((1<<depth)-1);

	if ( nh->node.leaves==NULL )
		return(0);

	off++;
	tree[tnode]=nh;
	tnode++;
	nh=&nh->node.leaves[bit];

  } while ( off*depth < bitsize );

 if ( nh->nchild <= 0 )
	return(0);

 nh->nchild=0;

 while ( tnode>0 )
  {	tnode--;
	nh=tree[tnode];
	nh->nchild--;
	if ( nh->nchild <= 0 )
	 {	free(nh->node.leaves);
		nh->node.leaves=NULL;
	 }
  };

 return(1);
}

/*****************************************************************************/

static int numhash_get_terminal_free(numhashtable *nt,int dir)
{
 int		capacity,cdepth,i,k,off,shift;
 numhashnode	*nh;

 if ( (1<<nt->bitsize) <= nt->table.nchild ) 
	return(-1);	/* no free node */

 off=0;
 nh=&nt->table;

 if ( nh->node.leaves==NULL )
	return(0);	/* no root node: 0 is the smallest free */

 for ( off=0,k=0 ; nh != NULL ; )
  {	shift=nt->depth*(((nt->bitsize-1)/nt->depth)-off);
	cdepth=(nt->bitsize-shift<nt->depth?nt->bitsize-shift:nt->depth);
	capacity=(1<<shift);
	if ( dir>=0 )	/* forward */
	 {	for ( i=0 ; i<(1<<cdepth) ; i++ )
		 {	if ( nh->node.leaves[i].nchild <= 0 )
				return(k);
			else if ( nh->node.leaves[i].nchild < capacity )
				break;
			else
				k+=capacity;
		 }
		if ( i==(1<<cdepth) )
			return(-1);	/* theoretically unreachable point */
	 }
	else		/* backward */
	 {	for ( i=(1<<cdepth)-1 ; i>=0 ; i-- )
		 {	if ( nh->node.leaves[i].nchild <= 0 )
				return((1<<nt->bitsize)-1-k);
			else if ( nh->node.leaves[i].nchild < capacity )
				break;
			else
				k+=capacity;
		 }
		if ( i<0 )
			return(-1);	/* theoretically unreachable point */
	 }
	nh=&nh->node.leaves[i];
	off++;

  };
 
 return(-1);	/* theoretically unreachable point */
}

int numhash_get_smallest_free(numhashtable *nt)
{
 return(numhash_get_terminal_free(nt,+1));
}
int numhash_get_largest_free(numhashtable *nt)
{
 return(numhash_get_terminal_free(nt,-1));
}

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

static int numhash_get_terminal_used(numhashtable *nt,int dir)
{
 int		capacity,cdepth,i,k,off,shift,tnd;
 numhashnode	*nh;

 if ( nt->table.nchild <= 0 )
	return(-1);	/* no used node */

 off=0;
 nh=&nt->table;

 if ( nh->node.leaves==NULL )
	return(-1);	/* no root node, so no used node */

 tnd=(nt->bitsize-1)/nt->depth;

 for ( off=0,k=0 ; nh != NULL ; )
  {	shift=nt->depth*(((nt->bitsize-1)/nt->depth)-off);
	cdepth=(nt->bitsize-shift<nt->depth?nt->bitsize-shift:nt->depth);
	capacity=(1<<shift);
	if ( dir>=0 )	/* forward */
	 {	for ( i=0 ; i<(1<<cdepth) ; i++ )
		 {	if ( nh->node.leaves[i].nchild && off<tnd )
				break;
			else if ( nh->node.leaves[i].nchild )
				return(k);
			else
				k+=capacity;
		 }
		if ( i==(1<<cdepth) )
			return(-1);	/* theoretically unreachable point */
	 }
	else		/* backward */
	 {	for ( i=(1<<cdepth)-1 ; i>=0 ; i-- )
		 {	if ( nh->node.leaves[i].nchild && off<tnd )
				break;
			else if ( nh->node.leaves[i].nchild )
				return((1<<nt->bitsize)-1-k);
			else
				k+=capacity;
		 }
		if ( i<0 )
			return(-1);	/* theoretically unreachable point */
	 }
	nh=&nh->node.leaves[i];
	off++;

  };
 
 return(-1);	/* theoretically unreachable point */
}

int numhash_get_smallest_used(numhashtable *nt)
{
 return(numhash_get_terminal_used(nt,+1));
}
int numhash_get_largest_used(numhashtable *nt)
{
 return(numhash_get_terminal_used(nt,-1));
}

/*****************************************************************************/

static	int numhash_local_walk(numhashnode *nh,int off,int num,
	int depth,int bitsize,int (*callback)(int,void *,void *),void *param)
{
 int	shift,cdepth,i,ret;

 if ( nh->node.leaves==NULL )
	return(0);

 shift=depth*(((bitsize-1)/depth)-off);
 cdepth=(bitsize-shift<depth?bitsize-shift:depth);

 ret=0;

 if ( shift>0 )
  {	for ( i=0 ; i<(1<<cdepth) ; i++ )
	 {	if ( nh->node.leaves[i].node.leaves != NULL )
			ret+=numhash_local_walk(&nh->node.leaves[i],off+1,
				num|(i<<shift),depth,bitsize,callback,param);
	 }
  }
 else if ( callback != NULL )
  {	for ( i=0 ; i<(1<<cdepth) ; i++ )
	 {	if ( nh->node.leaves[i].nchild>0 )
		 {	callback(num+i,nh->node.leaves[i].node.data,param);
			ret++;
		 }
	 }
  }
 else
  {	for ( i=0 ; i<(1<<cdepth) ; i++ )
	 {	if ( nh->node.leaves[i].nchild>0 )	ret++;		}
  }

 return(ret);
}

static	int numhash_local_walk_desc(numhashnode *nh,int off,int num,
	int depth,int bitsize,int (*callback)(int,void *,void *),void *param)
{
 int	shift,cdepth,i,ret;

 if ( nh->node.leaves==NULL )
	return(0);

 shift=depth*(((bitsize-1)/depth)-off);
 cdepth=(bitsize-shift<depth?bitsize-shift:depth);

 ret=0;

 if ( shift>0 )
  {	for ( i=(1<<cdepth)-1 ; i>=0 ; i-- )
	 {	if ( nh->node.leaves[i].node.leaves != NULL )
			ret+=numhash_local_walk_desc(&nh->node.leaves[i],off+1,
				num|(i<<shift),depth,bitsize,callback,param);
	 }
  }
 else if ( callback != NULL )
  {	for ( i=(1<<cdepth)-1 ; i>=0 ; i-- )
	 {	if ( nh->node.leaves[i].nchild>0 )
		 {	callback(num+i,nh->node.leaves[i].node.data,param);
			ret++;
		 }
	 }
  }
 else
  {	for ( i=0 ; i<(1<<cdepth) ; i++ )
	 {	if ( nh->node.leaves[i].nchild>0 )	ret++;		}
  }

 return(ret);
}

int numhash_walk(numhashtable *nt,
	int (*callback)(int,void *,void *),void *param)
{
 int	ret;
 ret=numhash_local_walk(&nt->table,0,0,nt->depth,nt->bitsize,callback,param);
 return(ret);
}

int numhash_walk_asc(numhashtable *nt,
	int (*callback)(int,void *,void *),void *param)
{
 return(numhash_walk(nt,callback,param));
}

int numhash_walk_desc(numhashtable *nt,
	int (*callback)(int,void *,void *),void *param)
{
 int	ret;
 ret=numhash_local_walk_desc(&nt->table,0,0,nt->depth,nt->bitsize,callback,param);
 return(ret);
}

int numhash_walk_dir(numhashtable *nt,
	int (*callback)(int,void *,void *),void *param,int dir)
{
 int	ret;

 if ( dir>=0 )
	ret=numhash_local_walk(&nt->table,0,0,nt->depth,nt->bitsize,callback,param);
 else 
	ret=numhash_local_walk_desc(&nt->table,0,0,nt->depth,nt->bitsize,callback,param);

 return(ret);
}

/*****************************************************************************/

int numhash_first(numhashtable *nt)
{
 return(numhash_get_terminal_used(nt,+1));
}
int numhash_first_wdata(numhashtable *nt,void **dret)
{
 int	c;
 c=numhash_get_terminal_used(nt,+1);
 if ( c>=0 && dret != NULL )	numhash_search(nt,c,dret);
 return(c);
}
int numhash_last(numhashtable *nt)
{
 return(numhash_get_terminal_used(nt,-1));
}
int numhash_last_wdata(numhashtable *nt,void **dret)
{
 int	c;
 c=numhash_get_terminal_used(nt,-1);
 if ( c>=0 && dret != NULL )	numhash_search(nt,c,dret);
 return(c);
}

static int numhash_local_first(numhashtable *nt,numhashnode *nh,
	int fcc,int off,void **dret)
{
 int	i,ntn,lcb,ret;

 if ( nh->node.leaves==NULL )
	return(-1);

 ntn=NUMHASH_NTN(nt->bitsize,nt->depth,off);
 lcb=NUMHASH_LCB(nt->bitsize,nt->depth,off);
 ret=fcc;
 for ( i=0 ; i<ntn ; )
  {	if ( nh->node.leaves[i].nchild )
	 {	if ( NUMHASH_RLF(nt->bitsize,nt->depth,off) )
		 {	nh=&nh->node.leaves[i];
			off++;
			ntn=NUMHASH_NTN(nt->bitsize,nt->depth,off);
			i=0;
			lcb=lcb>>nt->depth;
		 }
		else
		 {	if ( dret != NULL )
				*dret=nh->node.leaves[i].node.data;
			return(ret);
		 }
	 }
	else
	 {	i++;
		ret+=lcb;
	 }
  }
 return(-1);
}

static int numhash_local_last(numhashtable *nt,numhashnode *nh,
	int fcc,int off,void **dret)
{
 int	i,ntn,lcb,ret;

 if ( nh->node.leaves==NULL )
	return(-1);

 ntn=NUMHASH_NTN(nt->bitsize,nt->depth,off);
 lcb=NUMHASH_LCB(nt->bitsize,nt->depth,off);
 ret=fcc+lcb*(ntn-1);
 for ( i=ntn-1 ; i>=0 ; )
  {	if ( nh->node.leaves[i].nchild )
	 {	if ( NUMHASH_RLF(nt->bitsize,nt->depth,off) )
		 {	nh=&nh->node.leaves[i];
			off++;
			ntn=NUMHASH_NTN(nt->bitsize,nt->depth,off);
			i=ntn-1;
			lcb=lcb>>nt->depth;
			ret+=lcb*(ntn-1);
		 }
		else
		 {	if ( dret != NULL )
				*dret=nh->node.leaves[i].node.data;
			return(ret);
		 }
	 }
	else
	 {	i--;
		ret-=lcb;
	 }
  }
 return(-1);
}

static int numhash_local_nextprev(numhashtable *nt,numhashnode *nh0,
	int value,int off,int fcc,int *ret,int dir,void **dret)
{
 numhashnode	*nh;
 int		c,i,k,is_real_leaf,ntn,lcb;

 if ( nh0->node.leaves==NULL )
	return(0);

 k  =NUMHASH_BIT(nt->bitsize,nt->depth,off,value);
 nh=&nh0->node.leaves[k];
 lcb=NUMHASH_LCB(nt->bitsize,nt->depth,off);

 is_real_leaf=NUMHASH_RLF(nt->bitsize,nt->depth,off);
 if ( nh->node.leaves != NULL && is_real_leaf )
  {	c=numhash_local_nextprev(nt,nh,value,off+1,fcc+k*lcb,ret,dir,dret);
	if ( c )	return(1);
  }

 ntn=NUMHASH_NTN(nt->bitsize,nt->depth,off);

 if ( dir>=0 )	/* next */
  {	for ( i=k+1 ; i<ntn ; i++ )
	 {	nh=&nh0->node.leaves[i];
		if ( ! is_real_leaf && nh->nchild )
		 {	*ret=fcc+i;
			if ( dret != NULL )	*dret=nh->node.data;
			return(1);
		 }
		else if ( is_real_leaf && nh->nchild )
		 {	*ret=numhash_local_first(nt,nh,fcc+i*lcb,off+1,dret);
			return(1);
		 }
	 }
  }
 else		/* prev */
  {	for ( i=k-1 ; i>=0 ; i-- )
	 {	nh=&nh0->node.leaves[i];
		if ( ! is_real_leaf && nh->nchild )
		 {	*ret=fcc+i;
			if ( dret != NULL )	*dret=nh->node.data;
			return(1);
		 }
		else if ( is_real_leaf && nh->nchild )
		 {	*ret=numhash_local_last(nt,nh,fcc+i*lcb,off+1,dret);
			return(1);
		 }
	 }
  }

 return(0);
}

int numhash_next_wdata(numhashtable *nt,int value,void **dret)
{
 int	ret,c;
 c=numhash_local_nextprev(nt,&nt->table,value,0,0,&ret,+1,dret);
 if ( c	)	return(ret);
 else		return(-1);
}

int numhash_prev_wdata(numhashtable *nt,int value,void **dret)
{
 int	ret,c;
 c=numhash_local_nextprev(nt,&nt->table,value,0,0,&ret,-1,dret);
 if ( c	)	return(ret);
 else		return(-1);
}

int numhash_next(numhashtable *nt,int value)
{
 int	ret,c;
 c=numhash_local_nextprev(nt,&nt->table,value,0,0,&ret,+1,NULL);
 if ( c	)	return(ret);
 else		return(-1);
}

int numhash_prev(numhashtable *nt,int value)
{
 int	ret,c;
 c=numhash_local_nextprev(nt,&nt->table,value,0,0,&ret,-1,NULL);
 if ( c	)	return(ret);
 else		return(-1);
}

int numhash_loop_start(numhashtable *nt,int dir,void **ret)
{
 if ( dir>=0 )	return(numhash_first_wdata(nt,ret));
 else		return(numhash_last_wdata(nt,ret));
}
int numhash_loop_next(numhashtable *nt,int dir,int key,void **ret)
{
 if ( dir>=0 )	return(numhash_next_wdata(nt,key,ret));
 else		return(numhash_prev_wdata(nt,key,ret));
}
 
/*****************************************************************************/

static	int numhash_local_free(numhashnode *nh,int off,int depth,int bitsize)
{
 int	shift,cdepth,i;

 if ( nh->node.leaves==NULL )
	return(0);

 shift=depth*(((bitsize-1)/depth)-off);
 cdepth=(bitsize-shift<depth?bitsize-shift:depth);

 if ( shift>0 )
  {	for ( i=0 ; i<(1<<cdepth) ; i++ )
	 {	if ( nh->node.leaves[i].node.leaves != NULL )
			numhash_local_free(&nh->node.leaves[i],
				off+1,depth,bitsize);
	 }
  }

 free(nh->node.leaves);

 return(0);
}

int numhash_free(numhashtable *nt)
{
 numhash_local_free(&nt->table,0,nt->depth,nt->bitsize);
 numhash_init(nt,nt->bitsize,nt->depth);
 return(0);
}

/*****************************************************************************/
                         
