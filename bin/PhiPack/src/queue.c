/* This file is part of PhiPack.
 Copyright (c)2005, Trevor Bruen
 Foobar is free software: you can redistribute it and/or modify
 it under the terms of the GNU Lesser General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.
 
 Foobar is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU Lesser General Public License for more details.
 
 You should have received a copy of the GNU General Public License
 along with PhiPack.  If not, see <http://www.gnu.org/licenses/>.
 */
/*
   Copyright (c)2005, Trevor Bruen
   All rights reserved.                          

   Any feedback is very welcome.
   email: david.bryant@otago.ac.nz
*/

#include <stdlib.h>
#include <stdio.h>
#include "queue.h"
#include "mem.h"

void init_queue(queue *q)
{
  (*q).elements=(queue_elem *)mcalloc( (q->max_size),sizeof(queue_elem));
  q->front=0;
  q->back=0;
  q->cur_size=0;
}

void destroy(queue *q)
{
  free((*q).elements);
}

queue *enqueue(queue *q, queue_elem elem)
{
  if(q->cur_size == q->max_size)
    {
      return NULL;
    }
  else
    {
      (*q).elements[q->back]=elem;
      (q->cur_size)=(q->cur_size)+1;
      /* Loop around storage space */
      q->back=(q->back+1) % (q->max_size);
    }
  return q;
}

queue *element_at(queue *q, queue_elem *elem, int index)
{
  int mapped_index=0;
  if(index > q->cur_size)
    return NULL;
  else
    {
      mapped_index=((q->front)+(index)) % (q->max_size);
      *elem=(*q).elements[mapped_index];
      
    }
  return q;
}

void clear_queue(queue *q)
{
  q->front=0;
  q->back=0;
  q->cur_size=0;
}
void print_queue(queue *q)
{
  int cur=q->cur_size,i;
  queue_elem elem=0;
  fprintf(stdout,"Queue: ");
  for(i=0;i<cur;i++)
    {
      element_at(q,&elem,i);
      fprintf(stdout,"%d ",elem);

    }
  fprintf(stdout,"\n");
}
queue *dequeue_front(queue *q, queue_elem *elem)
{
  if(q->cur_size == 0)
    {
      //fprintf(stderr,"Cannot dequeue\n");
      return NULL;
    }
  else
    {
      *elem=(*q).elements[q->front];
      (q->cur_size)=(q->cur_size)-1;
      /* Loop around storage space */
      q->front=(q->front+1) % (q->max_size);
    
    }
  return q;
}
