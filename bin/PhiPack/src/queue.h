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
#ifndef QUEUE
#define QUEUE 1

typedef int queue_elem;

typedef struct {
  int max_size;
  int front;
  int back;
  int cur_size;
  queue_elem *elements;
} queue;



void init_queue(queue *q);

void destroy(queue *q);

queue *enqueue(queue *q, queue_elem elem);

queue *element_at(queue *q, queue_elem *elem, int index);

void clear_queue(queue *q);

void print_queue(queue *q);

queue *dequeue_front(queue *q, queue_elem *elem);

#endif 
