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
#ifndef MEM
#define MEM


size_t alloc_size();

void *mcalloc(size_t num, size_t size);

void *mmalloc(size_t size);


void *mrealloc(void *ptr, size_t old_size,size_t new_size);

void mfree(void *ptr, size_t size);

void error(const char * message); 


char ffclose(FILE *handle);


FILE *ffopen(char *name,char *mod);

#endif
