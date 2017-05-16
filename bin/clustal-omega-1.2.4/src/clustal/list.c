/* -*- mode: c; tab-width: 4; c-basic-offset: 4;  indent-tabs-mode: nil -*- */

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
 *  RCS $Id: list.c 156 2010-11-18 10:52:40Z andreas $
 *
 * Single linked list and FIFO/Queue
 *
 * Based on Kyle Loudon's Mastering Algorithms with C
 * http://oreilly.com/catalog/9781565924536
 *
 * Allows generic data types by using void* pointers, which works fine
 * as long as your data is just pointers.
 *
 */


#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <assert.h>
#include <string.h> /* for memset */

#include "list.h"

#ifdef LIST_TEST
#include <stdio.h>
#include "queue.h"
#include "time.h"
#endif

/**
 * @brief Initialise data members of a list
 *
 * @param[in] prList
 * List to initialise
 * @param[in] destroy
 * A function to be called with pointer to data when destroying the
 * list. NULL if in doubt, free in most other cases.
 * Note: doxygen will always fail to parse this...
 */
void
ListInit(list_t *prList, void (*destroy)(void *data)) {
    prList->size = 0;
    prList->destroy = destroy;
    prList->head = NULL;
    prList->tail = NULL;
    
    return;
}
/* end of ListInit() */



/**
 * @brief Calls user defined function to free data in list and resets
 * the list to NULL. Call even if your destroy function is NULL.
 *
 * @param[in] prList
 * The list to destroy
 *
 */
void
ListDestroy(list_t *prList) {
    void *pvData;

    while (LIST_SIZE(prList) > 0) {
        if (0 == ListRemoveNext(prList, NULL, (void **)&pvData)
            &&
            NULL != prList->destroy) {
            prList->destroy(pvData);
        }
    }
    /* precaution */
    memset(prList, 0, sizeof(list_t));
    
    return;
}
/* end of ListDestroy() */



/**
 *
 * @brief Insert data next to given element
 *
 * @param[in] prList
 * List into which to insert
 * @param[in] prElement
 * Current position/element. Element after which to insert. If NULL
 * head is used.
 * @param[in] pvData
 * Pointer to data to store
 *
 * @return Non-zero on failure
 */
int
ListInsertNext(list_t *prList, list_elem_t *prElement, const void *pvData)
{
    list_elem_t *prNewElement;

    if (NULL == (prNewElement = (list_elem_t *) malloc(sizeof(list_elem_t))))
        return -1;

    prNewElement->data = (void *)pvData;
 
    if (NULL == prElement) {
        /* insertion at head */
        if (LIST_SIZE(prList) == 0)
            prList->tail = prNewElement;
        prNewElement->next = prList->head;
        prList->head = prNewElement;
        
    } else {
        /* insert somewhere other than at head */
        if (NULL == prElement->next)
            prList->tail = prNewElement;
        prNewElement->next = prElement->next;
        prElement->next = prNewElement;
    }

    prList->size++;

    return 0;
}
/* end of ListInsertNext() */



/**
 * @brief Remove next element from current element/position.
 *
 * @param[in] prList
 * List from which an element is to be removed.
 * @param[in] prElement
 * Current element/position. Next item will be removed. If NULL head
 * is used.
 * @param[out] pvData_p
 * Will be pointed to removed elements data.
 *
 * @return Non-zero on failure
 */
int
ListRemoveNext(list_t *prList, list_elem_t *prElement, void **pvData_p)
{
    list_elem_t *prOldElement;

    if (0 == LIST_SIZE(prList))
        return -1;

    if (NULL == prElement) {
        /* handle removal from the head of the list */
        *pvData_p = prList->head->data;
        prOldElement = prList->head;
        prList->head = prList->head->next;
        if (1 == LIST_SIZE(prList))
            prList->tail = NULL;
        
    } else {
        /* handle removal from somewhere other than the head */        
        if (NULL == prElement->next)
            return -1;

        *pvData_p = prElement->next->data;
        prOldElement = prElement->next;
        prElement->next = prElement->next->next;
        
        if (NULL == prElement->next)
            prList->tail = prElement;
    }

    free(prOldElement);

    prList->size--;

    return 0;
}
/* end of ListRemoveNext() */



/**
 * @brief Insert int next to given element
 *
 * @param[in] prList
 * List into which to insert
 * @param[in] prElement
 * Current position/element. Element after which to insert. If NULL
 * head is used.
 * @param[in] data
 * int to store
 *
 * @return Non-zero on failure
 *
 */
int
IntListInsertNext(list_t *prList, list_elem_t *prElement, const int data)
{
    int *piInt;

    if (NULL == (piInt = malloc(sizeof(int)))) {
        return -1;
    }
    *piInt = data;
    
    return ListInsertNext(prList, prElement, piInt);
 }
/* end if IntListInsertNext() */


/**
 * @brief Remove next element from current element/position.
 *
 * @param[in] prList
 * List from which an element is to be removed.
 * @param[in] prElement
 * Current element/position. Next item will be removed. If NULL head
 * is used.
 * @param[out] iData_p
 * Will be pointed to removed elements data.
 *
 * @return Non-zero on failure
 *
 */
int
IntListRemoveNext(list_t *prList, list_elem_t *prElement, int *iData_p)
{
    int *piData;
    int res;
    res = ListRemoveNext(prList, prElement, (void **)&piData);
    *iData_p = *piData;
    prList->destroy(piData);
    return res;
}
/* end of IntListRemoveNext */




#ifdef LIST_TEST
/* gcc list.c -o list_test -ansi -Wall -DLIST_TEST */
   
int main(int argc, char **argv)
{

    int i;
    list_t *mylist;
    int_list_t *myintlist;
    queue_t *myqueue;
    int_queue_t *myintqueue;
    int res;

    int iSeed = (int)time(NULL);
    srand((unsigned int)iSeed);
    
    printf("%s", "list test! also delete #include!\n");
        

    mylist = malloc(sizeof(list_t));
    ListInit(mylist, NULL);

    printf("LIST test\n");

    for (i=0; i<argc; i++) {
        res = LIST_APPEND(mylist, argv[i]);
        printf("LIST Result for appending '%s' was %d\n", argv[i], res);
    }

    while (LIST_SIZE(mylist)) {
        char *argv_ptr;
        if (ListRemoveNext(mylist, NULL, (void **)&argv_ptr))
            perror("ListRemoveNext() failed");
        printf("LIST Popped %s\n", argv_ptr);
    }
    printf("LIST %s", "No more elements to pop");
        
    /* could become list_free */
    ListDestroy(mylist);
    free(mylist);

        


    myintlist = malloc(sizeof(list_t));
    INT_LIST_INIT(myintlist);

    printf("\n");
    printf("%s", "INT_LIST test");

    for (i=0; i<argc; i++) {
        int data = 666-i;
        res = INT_LIST_APPEND(myintlist, data);
        printf("INT_LIST Result for appending '%d' was %d\n", data, res);
    }

    while (INT_LIST_SIZE(myintlist)) {
        int data;
        if (IntListRemoveNext(myintlist, NULL, &data))
            perror("ListRemoveNext() failed\n");
        printf("INT_LIST Popped %d\n", data);
    }
    printf("INT_LIST %s\n", "No more elements to pop");
        
    /* could become list_free */
    INT_LIST_DESTROY(myintlist);
    free(myintlist);




    
    myqueue = malloc(sizeof(queue_t));
    QUEUE_INIT(myqueue, NULL);

    printf("\n");
    printf("%s", "QUEUE test\n");

    for (i=0; i<argc; i++) {
        res = QUEUE_PUSH(myqueue, argv[i]);
        printf("QUEUE Result for pushing '%s' was %d\n", argv[i], res);
    }

    while (! QUEUE_EMPTY(myqueue)) {
        char *argv_ptr;
        if (QUEUE_POP(myqueue, (void **)&argv_ptr))
            perror("QUEUE_POP() failed\n");
        printf("QUEUE Popped %s\n", argv_ptr);
    }
    printf("QUEUE %s\n", "QUEUE No more elements to pop");
        
    /* could become list_free */
    QUEUE_DESTROY(myqueue);
    free(myqueue);




    myintqueue = malloc(sizeof(queue_t));
    INT_QUEUE_INIT(myintqueue);

    printf("\n");
    printf("%s\n", "INT_QUEUE test");

    for (i=0; i<argc; i++) {
        res = INT_QUEUE_PUSH(myintqueue, i);
        printf("INT_QUEUE Result for appending '%d' was %d\n", i, res);
    }

    while (! INT_QUEUE_EMPTY(myintqueue)) {
        int data;
        int rand_data;
        if (INT_QUEUE_POP(myintqueue, &data))
            perror("INT_QUEUE_POP() failed\n");
        printf("INT_QUEUE Popped %d\n", data);

        rand_data = (int)( 10.0 * rand() / ( RAND_MAX+1.0));
        if (! (rand_data%3)) {
            res = INT_QUEUE_PUSH(myintqueue, rand_data);
            printf("INT_QUEUE Result for pushing random number '%d' was %d\n", rand_data, res);            
        }
    }
    printf("INT_QUEUE %s\n", "INT_QUEUE No more elements to pop");
        
    /* could become list_free */
    INT_QUEUE_DESTROY(myintqueue);
    free(myintqueue);
    

    exit(0);
}

#endif



