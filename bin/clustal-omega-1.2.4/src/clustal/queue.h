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
 * RCS $Id: queue.h 193 2011-02-07 15:45:21Z andreas $
 *
 * Functions/Macros for FIFOs/Queues
 *
 */

#ifndef CLUSTALO_QUEUE_H
#define CLUSTALO_QUEUE_H

#include "list.h"



/*  FIFO/Queue as list_t, storing data pointers
 * 
 */

typedef list_t queue_t;

/* setup queue */
#define QUEUE_INIT(prQueue, destroy_func) ListInit((prQueue), (destroy_func))

/* free all elements from queue */
#define QUEUE_DESTROY(prQueue) ListDestroy((prQueue))

/* enqueue */
#define QUEUE_PUSH(prQueue, data) LIST_APPEND((prQueue), (data))

/* dequeue */
#define QUEUE_POP(prQueue, data) ListRemoveNext((prQueue), NULL, (data))

/* is queue empty ? */
#define QUEUE_EMPTY(prQueue) (0==LIST_SIZE((prQueue)))



/* Special int FIF/Queue, storing ints by copying them instead of
 * keeping pointers only
 */

typedef queue_t int_queue_t;

/* setup queue */
#define INT_QUEUE_INIT(prQueue) INT_LIST_INIT((prQueue))

/* free all elements from queue */
#define INT_QUEUE_DESTROY(prQueue) INT_LIST_DESTROY((prQueue))

/* enqueue */
#define INT_QUEUE_PUSH(prQueue, data) INT_LIST_APPEND((prQueue), (data))

/* dequeue */
#define INT_QUEUE_POP(prQueue, data) IntListRemoveNext((prQueue), NULL, (data))

/* is queue empty ? */
#define INT_QUEUE_EMPTY(prQueue) (0==INT_LIST_SIZE((prQueue)))



#endif
