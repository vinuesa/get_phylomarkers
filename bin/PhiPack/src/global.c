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
const char DNA_ALPHA[] = {'A','C','G','T'};
const int DNA_ALPHA_SIZE=4;

const char DNA_MISSING[] = {'X','N','?'};
const int DNA_MISSING_SIZE=3;

const char DNA_AMBIG[] = {'R','Y','M','K','S','W','H','B','V','D'};
const int DNA_AMBIG_SIZE = 10;

const char AA_ALPHA[] = {'A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y'};
const int AA_ALPHA_SIZE=20;
 
const char AA_MISSING[] ={'X','?'};
const int AA_MISSING_SIZE =2;

const char AA_AMBIG[] = {'B','Z'};
const int AA_AMBIG_SIZE =2;

const char GAP[] = {'-'};
const int GAP_SIZE=2;

cbool memberOf(const char *set, const int num, char ch)
{
  int i;
 
  for(i=0;i<num;i++)
    if(ch == (set[i]-CHAR_START))
      return TRUE;
  return FALSE;
}

cbool validState(alignmentClass alignKind, char ch)
{

  switch(alignKind)
    {
    case DNA:
      //  fprintf(stdout,"checking based %c -have %d\n",ch,memberOf(DNA_ALPHA,DNA_ALPHA_SIZE,ch));
      return memberOf(DNA_ALPHA,DNA_ALPHA_SIZE,ch);
      break;
    case AA:
      return memberOf(AA_ALPHA,AA_ALPHA_SIZE,ch);
      break;
    case OTHER:
      if((ch != (GLOBAL_GAP-CHAR_START)) && (ch != (GLOBAL_MISSING-CHAR_START)))
	return TRUE;
      else
	return FALSE;
      
    }

  return FALSE;
}

cbool valid_gap(alignmentClass alignKind, char ch)
{
  return (memberOf(GAP,GAP_SIZE,ch));
}
cbool missing_ambig_State(alignmentClass alignKind, char ch)
{

  switch(alignKind)
    {
    case DNA:
      return (memberOf(DNA_MISSING,DNA_MISSING_SIZE,ch) || memberOf(DNA_AMBIG,DNA_AMBIG_SIZE,ch));
      break;
    case AA:
      return (memberOf(AA_MISSING,AA_MISSING_SIZE,ch) || memberOf(AA_AMBIG,AA_AMBIG_SIZE,ch));
      break;
    case OTHER:
      if(ch == GLOBAL_MISSING)
	return TRUE;
      else
	return FALSE;
    }
  
  return FALSE;
}
