/*****************************************************************************/
/* list.h								     */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/* Simple set of macros for manipulating linked and doubly linked lists.     */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/* Copyright (C) 2006, 2008; Pal, A. (apal@szofi.elte.hu)		     */
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

#ifndef	__LIST_H_INCLUDED
#define	__LIST_H_INCLUDED	1

/* 
   Each object used in the linked list shoud contain two fields with the 
   type of a pointer to the object itself. These pointers are then used to 
   define the links between the list elements. Both in the case of the linked
   list and doubly-linked list, the object must contain a 'prev' and 'next'
   field. To traverse on the list, one can use the for ( ; ; ) loop directly:

   for ( o=firstpointer ; o != NULL ; o=o->next ) { ... };
   for ( o=lastpointer  ; o != NULL ; o=o->prev ) { ... };

   (the latter one can only be used in doubly-linked lists).
*/

#define	__LIST_PREV		prev
#define	__LIST_NEXT		next

#ifndef	__LIST_MALLOC
#define	__LIST_MALLOC		malloc
#endif

#ifndef	__LIST_FREE
#define	__LIST_FREE		free
#endif

/*****************************************************************************/

/* list_new():
   Simple wrapper to malloc/calloc. */
#define		list_new(type)		((type *)__LIST_MALLOC(sizeof(type)))

/* list_insert_first():
   Inserts the object 'object' to the beginning of the list of which first
   element is 'firstpointer'. 						     */
#define		list_insert_first(firstpointer,object) \
 do 							\
  {	(object)->__LIST_NEXT=(firstpointer); 		\
	(object)->__LIST_PREV=NULL; 			\
	if ( (firstpointer) != NULL ) 			\
		(firstpointer)->__LIST_PREV=(object); 	\
	(firstpointer)=(object);			\
  } while ( 0 )

/* list_remove():
   Deletes the object 'object' from the list of which first element is
   'firstpointer'.							     */
#define		list_remove(firstpointer,object) \
 do 									  \
  {	if ( (object)->__LIST_NEXT != NULL )				  \
		(object)->__LIST_NEXT->__LIST_PREV=(object)->__LIST_PREV; \
	if ( (object)->__LIST_PREV != NULL )				  \
		(object)->__LIST_PREV->__LIST_NEXT=(object)->__LIST_NEXT; \
	else								  \
		(firstpointer)=(object)->__LIST_NEXT; 			  \
  } while ( 0 )

/* list_delete():
   Like list_remove(), followed by a free() to the ``object'' pointer.       */
#define		list_delete(firstpointer,object) \
 do 									  \
  {	list_remove((firstpointer),(object));				  \
	__LIST_FREE(object);						  \
  } while ( 0 )								  \

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/* dlist_new():
   Simple wrapper to malloc/calloc. */
#define		dlist_new(type)		((type *)__LIST_MALLOC(sizeof(type)))

/* dlist_insert_first(): 
   Inserts the object 'object' to the beginning of the list of which first
   element is 'firstpointer' and last element is 'lastpointer':		     */
#define		dlist_insert_first(firstpointer,lastpointer,object) \
 do 							\
  {	(object)->__LIST_NEXT=(firstpointer); 		\
	(object)->__LIST_PREV=NULL; 			\
	if ( (firstpointer) != NULL )			\
		(firstpointer)->__LIST_PREV=(object);	\
	(firstpointer)=(object); 			\
	if ( (lastpointer)  == NULL )			\
		(lastpointer)=(object); 		\
  } while ( 0 )

/* dlist_insert_last(): 
   Inserts the object 'object' to the end of the list of which first
   element is 'firstpointer' and last element is 'lastpointer':		     */
#define		dlist_insert_last(firstpointer,lastpointer,object) \
 do 							\
  {	(object)->__LIST_PREV=(lastpointer); 		\
	(object)->__LIST_NEXT=NULL; 			\
	if ( (lastpointer)  != NULL )			\
		(lastpointer)->__LIST_NEXT=(object); 	\
	(lastpointer)=(object); 			\
	if ( (firstpointer) == NULL )			\
		(firstpointer)=(object); 		\
  } while ( 0 ) 

/* dlist_remove():
   Deletes the object 'object' from the list of which first element is
   'firstpointer' and last element is 'lastpointer':		    	     */
#define		dlist_remove(firstpointer,lastpointer,object) \
 do									  \
  {	if ( (object)->__LIST_NEXT != NULL )				  \
		(object)->__LIST_NEXT->__LIST_PREV=(object)->__LIST_PREV; \
	else								  \
		(lastpointer) =(object)->__LIST_PREV; 			  \
	if ( (object)->__LIST_PREV != NULL )				  \
		(object)->__LIST_PREV->__LIST_NEXT=(object)->__LIST_NEXT; \
	else								  \
		(firstpointer)=(object)->__LIST_NEXT; 			  \
  } while ( 0 )

/* dlist_delete():
   Like dlist_remove(), followed by a free() to the ``object'' pointer.      */
#define		dlist_delete(firstpointer,lastpointer,object) \
 do 									  \
  {	list_remove((firstpointer),(lastpointer),(object));		  \
	__LIST_FREE(object);						  \
  } while ( 0 )								  \

/*****************************************************************************/

#endif

/*****************************************************************************/
                 
