/*********************************************************************
 * Clustal Omega - Multiple sequence alignment
 *
 * Copyright (C) 2010 University College Dublin
 *
 * Clustal-Omega is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License as
 * published by the Free Software Foundation; either version 2 of the
 * License, or (at your option) any later version.
 *
 * This file is part of Clustal-Omega.
 *
 ********************************************************************/

/*
 * RCS $Id: list.h 193 2011-02-07 15:45:21Z andreas $
 *
 * Generic single linked list storing pointers to data
 *
 */

#ifndef CLUSTALO_LIST_H
#define CLUSTALO_LIST_H

#include <stdlib.h>

typedef struct list_elem_s {
    void  *data;
    struct list_elem_s *next;
} list_elem_t;

typedef struct {
    /* size of list */
    int size;
    /* user defined function for freeing data */
    void (*destroy)(void *data);
    list_elem_t *head;
    list_elem_t *tail;
} list_t;

void ListInit(list_t *prList, void (*destroy)(void *data));

void ListDestroy(list_t *prList);

int ListInsertNext(list_t *prList, list_elem_t *prElement, const void *data);

#define LIST_APPEND(prList, data)  ListInsertNext((prList), LIST_TAIL(prList), (data))

#define LIST_PREPEND(prList, data)  ListInsertNext((prList), CLUSTALO_LIST_HEAD(prList), (data))

int ListRemoveNext(list_t *prList, list_elem_t *prElement, void **data);

#define LIST_SIZE(prList) ((prList)->size)

#define CLUSTALO_LIST_HEAD(prList) ((prList)->head)

#define LIST_TAIL(prList) ((prList)->tail)

#define LIST_IS_HEAD(prList, prElement) ((prElement) == (prList)->head ? 1 : 0)

#define LIST_IS_TAIL(prElement) ((prElement)->next == NULL ? 1 : 0)

#define LIST_DATA(prElement) ((prElement)->data)

#define LIST_NEXT(prElement) ((prElement)->next)





/* special int list: stores ints by copying them (instead of storing
 * pointers as generic list)
 *
 */

typedef list_t int_list_t;

#define INT_LIST_INIT(prList)  ListInit((prList), free)

#define INT_LIST_DESTROY(prList) ListDestroy((prList));

int IntListInsertNext(list_t *prList, list_elem_t *prElement, const int data);

#define INT_LIST_APPEND(prList, data) IntListInsertNext((prList), LIST_TAIL(prList), (data))

#define INT_LIST_PREPEND(prList, data) IntListInsertNext((prList), CLUSTALO_LIST_HEAD(prList), (data))

int IntListRemoveNext(list_t *prList, list_elem_t *prElement, int *data);

#define INT_LIST_SIZE(prList)   LIST_SIZE(prList)

#define INT_CLUSTALO_LIST_HEAD(prList) CLUSTALO_LIST_HEAD_INT((prList))

#define INT_LIST_TAIL(prList)  LIST_TAIL_INT((prList) )

#define INT_LIST_IS_HEAD(prList, prElement) LIST_IS_HEAD(prList, prElement) 

#define INT_LIST_IS_TAIL(prElement) LIST_IS_TAIL((prElement))

#define INT_LIST_DATA(prElement) LIST_DATA((prElement))

#define INT_LIST_NEXT(prElement) LIST_NEXT((prElement))


#endif
