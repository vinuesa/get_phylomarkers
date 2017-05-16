/*****************************************************************************/
/* numhash.h								     */
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

#ifndef	__NUMHASH_H_INCLUDED
#define	__NUMHASH_H_INCLUDED	1

/*****************************************************************************/

typedef struct	numhashnode	numhashnode;

struct numhashnode
 {	union
	 {	numhashnode	*leaves;
		void		*data;
	 } node;
	int		nchild;
 };

typedef struct
 {	numhashnode	table;
	int		depth;
	int		bitsize;
 } numhashtable;

/*****************************************************************************/

/* numhash_init():
   Initializes the hash table 'nt' to be capable to store 1<<bitsize numbers
   (between 0 and 1<<bitsize-1). The size of the subnodes of the hash tree
   is going to be 1<<depth. 						     */
int	numhash_init(numhashtable *nt,int bitsize,int depth);

/* numhash_add():
   Adds a new key 'key' to the hash table with the optional data 'data'. 
   'data' either can be NULL, it is a valid entry. The function returns
   0 if the key 'key' has been found in the table. In this case, the data
   is updated to 'data'. Otherwise (if 'key' is a new key) the function
   returns a positive value.						     */
int	numhash_add(numhashtable *nt,int key,void *data);

/* numhash_search():
   Searches for the key 'key' in the hash table 'nt'. If the key is found
   in the table, the function returns a non-zero key and if 'ret' is not 
   NULL, stores the optional user-supplied data in *ret. If the key is not
   found in the hash table, the function returns zero and *ret is unchanged. */
int	numhash_search(numhashtable *nt,int key,void **ret);

/* numhash_remove():
   Removes the key 'key' from the hash table 'nt'. If the key was found
   in the table, the function returns a positive value. Otherwise, if the
   key is not found in the hash table, the function returns zero. In this
   case the table is unchanged.						     */
int	numhash_remove(numhashtable *nt,int key);

/* numhash_total():
   Returns the total number of keys stored in the hash table 'nt'.	     */
int	numhash_total(numhashtable *nt);

/* numhash_get_{smallest,largest}_{free,used}():
   These functions return the smallest or largest free or used key in the
   hash table. If the table is empty (i.e.numhash_total() is zero), the 
   "used" functions return a negative value. If the table is full
   (i.e. numhash_total() equals to 1<<bitsize), the "free" functions
   return a negative value.						     */
int	numhash_get_smallest_free(numhashtable *nt);
int	numhash_get_largest_free(numhashtable *nt);
int	numhash_get_smallest_used(numhashtable *nt);
int	numhash_get_largest_used(numhashtable *nt);

/* numhash_walk(), numhash_walk_asc(), numhash_walk_desc(), numhash_walk_dir():
   These functions walk through the used keys in the hash table 'nt' in
   ascending or descending order. If 'callback' is not NULL, it is called 
   each time when a key is used: the first parameter is the numeric key itself,
   the second is the associated (and user-supplied) data and the last
   is the optional parameter 'param'. The function numhash_walk_dir() walks
   through in ascending order (forward) if 'dir' is positive or zero, otherwise
   it walks through in descending order (backward).			     */
int	numhash_walk(numhashtable *nt,
	int (*callback)(int,void *,void *),void *param);
int	numhash_walk_asc(numhashtable *nt,
	int (*callback)(int,void *,void *),void *param);
int	numhash_walk_desc(numhashtable *nt,
	int (*callback)(int,void *,void *),void *param);
int	numhash_walk_dir(numhashtable *nt,
	int (*callback)(int,void *,void *),void *param,int dir);

/* numhash_first(), numhash_last():
   These functions do the same what numhash_get_smallest_used() and 
   numhash_get_largest_used() do.					     */
int	numhash_first(numhashtable *nt);
int	numhash_last(numhashtable *nt);
/* numhash_next(), numhash_prev():
   These functions returns the leaf with the smallest (or largest) key after 
   (before) the key 'key'. If there is no keys after (before) the specified 
   'key', or the hash table is empty, these functions return a negative 
   value. These functions can also be used for classic loop control:
	for ( v=numhash_first(table) ; 0 <= v ; v=numhash_next(table,v) );
   or like so for descending order:
	for ( v=numhash_last(table)  ; 0 <= v ; v=numhash_prev(table,v) );   */
int	numhash_next(numhashtable *nt,int key);
int	numhash_prev(numhashtable *nt,int key);

/* numhash_{first|last}_wdata(), numhash_{next|prev}_wdata():
   These functions do the same what numhash_{frist|last}() and 
   numhash_{next|prev}() but also return the associated data of the key if
   the pointer 'ret' is not NULL.					     */
int	numhash_first_wdata(numhashtable *nt,void **ret);
int	numhash_last_wdata(numhashtable *nt,void **ret);
int	numhash_next_wdata(numhashtable *nt,int key,void **ret);
int	numhash_prev_wdata(numhashtable *nt,int key,void **ret);

/* numhash_loop_start(), numhash_loop_next():
   These functions unify the forward and backward loop control functions. 
   The direction can be specified with 'dir'. If it is positive or zero, the
   loop is going to be a forward loop, otherwise the loop is a backward loop.*/
int	numhash_loop_start(numhashtable *nt,int dir,void **ret);
int	numhash_loop_next(numhashtable *nt,int dir,int key,void **ret);

/* numhash_free():
   Releases the hash table 'nt'.					     */
int	numhash_free(numhashtable *nt);

/*****************************************************************************/

#endif

/*****************************************************************************/
            
