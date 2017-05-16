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

#include "global.h"
#include "mem.h"

static size_t allocated=0;
static int files_open=0;

size_t alloc_size()
{
  return allocated;
}


void *mcalloc(size_t num, size_t size)
{
  void *value = calloc (num,size);
  if (value == NULL)
    {
    error("Memory Exhausted");
    }
  else
    {
      /*allocated=allocated+num*size; */
    }
  return value;
}



void *mmalloc(size_t size)
{
  void *value = malloc (size);
  if (value == NULL)
    {
      error("Memory Exhausted");
    }
  else
    {
      /* allocated=allocated+size; */
    }
  return value;
}

void mfree(void *ptr, size_t size)
{
  free(ptr);

  /* allocated=allocated-size; */

}

void *mrealloc(void *ptr, size_t old_size,size_t new_size)
{
  void *value = realloc (ptr,new_size);
  if (value == NULL)
    {
    error("Memory Exhausted");
    }
  else
    {
      /*allocated=allocated + new_size-old_size */
    }
  return value;
}


FILE *ffopen(char *name,char *mod)
{
  char s[MAX_SIZE+1];
  FILE *f;
  f=fopen(name,mod);

  if(f == NULL)
    {
      sprintf(s,"Could not open file %s",name);
      error(s);
    }
  files_open++;
  return f;
}

char ffclose(FILE *handle)
{
  char f;
  char s[MAX_SIZE+1];

  f=fclose(handle);
  
  if(f == EOF)
    {
      sprintf(s,"Could not close file hanlde ");
      error(s);
    }
  files_open--;
  return f;
}

