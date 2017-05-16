/*  
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
 

   Any feedback is very welcome.
   email: david.bryant@otago.ac.nz
*/

#include <stdlib.h>
#include "global.h"
#include "mem.h"
#include "math.h"


#define RED 0
#define BLUE 1
#define GREEN 2


#define MIN(a,b) ((a<b) ? (a):(b))

typedef unsigned char image_type;
/* Convert a color from hue/saturation/value space to rgb space */
/* hsv between 0 (red) and (360) */
void hsv2rgb(double hsv_h, double hsv_s,double hsv_v,  double **rgb_col, int index)
{
  
   if (hsv_h < 120) 
     {
       rgb_col[index][RED] = (120 - hsv_h) / 60.0;
       rgb_col[index][GREEN] = hsv_h / 60.0;
       rgb_col[index][BLUE] = 0;
     } 
   else if (hsv_h < 240) 
     {
       rgb_col[index][RED] = 0;
       rgb_col[index][GREEN] = (240 - hsv_h) / 60.0;
       rgb_col[index][BLUE] = (hsv_h - 120) / 60.0;
     } 
   else 
     {
       
       rgb_col[index][RED] = (hsv_h - 240) / 60.0;
       rgb_col[index][GREEN] = 0;
       rgb_col[index][BLUE] = (360 - hsv_h) / 60.0;
     }
   rgb_col[index][RED] = MIN(rgb_col[index][RED],1);
   rgb_col[index][GREEN] = MIN(rgb_col[index][GREEN],1);
   rgb_col[index][BLUE] = MIN(rgb_col[index][BLUE],1);
   
   rgb_col[index][RED] = 255*(1 - hsv_s + hsv_s * rgb_col[index][RED]) * hsv_v;
   rgb_col[index][GREEN] = 255*(1 - hsv_s + hsv_s * rgb_col[index][GREEN]) * hsv_v;
   rgb_col[index][BLUE] = 255*(1 - hsv_s + hsv_s * rgb_col[index][BLUE]) * hsv_v;
}


void write_PPM(FILE *out_file,image_type **r_array, image_type **g_array, image_type **b_array,int y_size,int x_size)
{
  int i, j;
  fprintf(out_file,"P3\n%d %d\n255\n",x_size,y_size); 
  
  for(i=0;i<y_size;i++)
    {
      for(j=0;j<x_size;j++)
	{
	  fprintf(out_file,"%d %d %d",r_array[i][j],g_array[i][j],b_array[i][j]);
	  if( ((j+1)% 4 == 0))
	    fprintf(out_file,"\n");
	  else
	    fprintf(out_file," ");
	}
      fprintf(out_file,"\n");
    }
}


void draw_line(image_type ***r_array,image_type ***g_array,image_type ***b_array,int y0,int y1,int x0,int x1,int r_0,int g_0, int b_0)
{
  cbool steep=FALSE;
  int deltax,deltay,deltaerr,x,y,xstep=-1,ystep=-1,temp,error=0;

  if(abs(y1-y0) > abs (x1 -x0))
    {
      steep=TRUE;
      /* Swap x0 and y0 */
      temp=x0;
      x0=y0;
      y0=temp;
      /* Swap x1 and y1 */
      temp=x1;
      x1=y1;
      y1=temp;
    }

  deltax=abs(x1-x0);
  deltay=abs(y1-y0);
  deltaerr=deltay;
  x=x0;
  y=y0;
  if(x0 < x1)
    xstep=1;
  if(y0 < y1)
    ystep=1;
  
 
  if(steep==TRUE)
    {
      (*r_array)[x][y]=r_0;
      (*g_array)[x][y]=g_0;
      (*b_array)[x][y]=b_0;

    }
  else
    {
      (*r_array)[y][x]=r_0;
      (*g_array)[y][x]=g_0;
      (*b_array)[y][x]=b_0;
    }

  while(x != x1)
    {
      x=x+xstep;
      error=error+deltaerr;
      if(2*error >= deltax)
	{
	  y=y+ystep;
	  error=error - deltax;
	}
    
      if(steep==TRUE)
	{
	  (*r_array)[x][y]=r_0;
	  (*g_array)[x][y]=g_0;
	  (*b_array)[x][y]=b_0;
	  
	}
      else
	{
	  (*r_array)[y][x]=r_0;
	  (*g_array)[y][x]=g_0;
	  (*b_array)[y][x]=b_0;
	}
      
    }
  
}

 

void draw_box(image_type ***r_array,image_type ***g_array,image_type ***b_array,int y1,int y2,int x1,int x2,int r_0,int g_0, int b_0)
{
  int i,j;
  
  for(i=y1;i<y2;i++)
    for(j=x1;j<x2;j++)
      {
	(*r_array)[i][j]=r_0;
	(*g_array)[i][j]=g_0;
	(*b_array)[i][j]=b_0;
      }
}


void create_profile_pic(int num_sites,int k,int num_inf, double *break_values,site *site_desc)
{
  FILE *out_file;
  int cur_value,next_value;
  int i,j,temp;
  double max_value=0,ratio;
  double peak_in=6.63;
  /* Size of graph - total number of states available */
  
  int header=100,margin=100;
  int height=400+2*header,width=num_sites+2*margin;
  
  /* Array of actual image */
  image_type **r_array,**g_array,**b_array;

  r_array = (image_type **)mmalloc(height * sizeof(int *));
  g_array = (image_type **)mmalloc(height * sizeof(int *));
  b_array = (image_type **)mmalloc(height * sizeof(int *));
  
  for(i=0;i<height;i++)
    {      
      r_array[i] = (image_type *)mmalloc(width * sizeof(int));
      g_array[i] = (image_type *)mmalloc(width * sizeof(int));
      b_array[i] = (image_type *)mmalloc(width * sizeof(int));
    }
  for(i=0;i<height;i++)
    {
      for(j=0;j<width;j++)
	{
	  r_array[i][j] = 255;
	  g_array[i][j] = 255;
	  b_array[i][j] = 255;
	}
    }
  for(i=0;i<num_inf;i++)
    {
      if(break_values[i]> max_value)
	max_value=break_values[i];
    }

  temp=(int)(max_value/10) + 1;
  ratio=40.0/((double)temp);
  /* Draw borders */
  draw_line(&r_array,&g_array,&b_array,height-header+1,height-header+1,0+margin,num_sites+margin,0,0,0);
  
  draw_line(&r_array,&g_array,&b_array,height-header+1,header,margin-1,margin-1,0,0,0);
  
  /* Draw ticks */
  /* First height ticks */
  
  for(i=1;i<2*temp;i++)
    {
      cur_value=(ratio*5*i);
      draw_line(&r_array,&g_array,&b_array,height-((int)cur_value+header),height-((int)cur_value+header),margin-11,margin-1,0,0,0);
    }
  /* Draw line at sig <0.01 */
  cur_value=(ratio*peak_in);
  draw_line(&r_array,&g_array,&b_array,height-((int)cur_value+header),height-((int)cur_value+header),margin,margin+num_sites,255,0,0);
  
  /* Next draw site ticks */
  for(i=1;i<num_sites;i++)
    {
      if(i % 500 == 0)
	{
	  draw_line(&r_array,&g_array,&b_array,height-header+1,height-header+21,margin+i,margin+i,0,0,0);
	}
      else if(i % 100 == 0)
	{
	  draw_line(&r_array,&g_array,&b_array,height-header+1,height-header+11,margin+i,margin+i,0,0,0);
	}
    }
  /* Add profile */
  draw_line(&r_array,&g_array,&b_array,height-header,height-header,0+margin,site_desc[k].orig_index+margin,0,0,255);
  for(i=0;i<(num_inf-2*k);i++)
    {
      cur_value=(ratio*break_values[i]);
      next_value=(ratio*break_values[i+1]);
             
      draw_line(&r_array,&g_array,&b_array,height-((int)cur_value+header),height-((int)next_value+header),site_desc[i+k].orig_index+margin,site_desc[i+1+k].orig_index+margin,0,0,255);
      
     
    }
  draw_line(&r_array,&g_array,&b_array,height-header,height-header,margin+site_desc[i+k].orig_index,num_sites+margin,0,0,255);
 
  /* Print image */
  out_file=fopen("profile.ppm","w");
  write_PPM(out_file,r_array,g_array,b_array,height,width);

}

void create_incompat_pic(int num_states,int max_inc,inc_type **inc_matrix, cbool embellish, int num_chars)
{

  FILE *out_file;
  double hue, saturation, value;

  int i,j,i2,j2,cur;

  /* Size of graph - total number of states available */
  int size;
  
  /* Embellish sizes */
 
  int tick_height_small=12;
  int freq_small_tick=10;
  int freq_big_tick=50;
 
  int tick_height_big=24;
  int tick_width_small=1,tick_width_large=3;
 
  int header=tick_height_big+5;
  int block=5;
  int spacer=1;
  int b_size=block+spacer;
  int border=5;
  int offset=header+border;
  int cur_a,cur_b;

  /* 60 degrees - i.e. yellow */
  int start_color = 60;

  /* 0 degrees - i.e. red */
  int end_color = 0;
  
  double ratio=0.5;
  /* Array of actual image */
  image_type **r_array,**g_array,**b_array;


  /* Map to colors */
  double **palette;
  
  /* Size of graph */
  if(embellish)
    size=(num_chars)*(b_size) - spacer+2*offset;
  else
    size=num_chars;
    
  /* Allocate memory */
  palette = (double **)mmalloc(max_inc * sizeof(double *)); 
  
  for(i=0;i<max_inc;i++) 
    palette[i]=(double *)mmalloc( 3 * sizeof(double));
  
  r_array = (image_type **)mmalloc(size * sizeof(int *));
  g_array = (image_type **)mmalloc(size * sizeof(int *));
  b_array = (image_type **)mmalloc(size * sizeof(int *));
  
  for(i=0;i<size;i++)
    {
      r_array[i] = (image_type *)mmalloc(size * sizeof(image_type));
      g_array[i] = (image_type *)mmalloc(size * sizeof(image_type));
      b_array[i] = (image_type *)mmalloc(size * sizeof(image_type));
    }
  
  
  for(i=0;i<size;i++)
    for(j=0;j<size;j++)
      {
	r_array[i][j] = 32;
	g_array[i][j] = 32;
	b_array[i][j] = 32;
      }
  
  /* Create colors */
  for(i=0;i<max_inc;i++)
    {
      saturation =1.0;
      value=1.0;
      hue=end_color-((double)(end_color - start_color)) * pow(ratio,i);
           
      hsv2rgb(hue,saturation,value,palette,i);
    }
  /* Add border if embellish */
  if(embellish)
    {
      /* Fill background with white */
      for(i=0;i<size;i++)
	for(j=0;j<size;j++)
	  {
	    if((i<offset) || (i>= (size-offset)) || (j<offset) || (j>= (size-offset)))
	      {
		r_array[i][j] = 255;
		g_array[i][j] = 255;
		b_array[i][j] = 255;
	      }
	  }
      
      /* Add borders */
      draw_box(&r_array,&g_array,&b_array,header,offset,offset,size-offset,0,0,0);
      draw_box(&r_array,&g_array,&b_array,size-offset,size-header,offset,size-offset,0,0,0);
      draw_box(&r_array,&g_array,&b_array,header,size-header,header,offset,0,0,0);
      draw_box(&r_array,&g_array,&b_array,header,size-header,size-offset,size-header,0,0,0);
      
      /* Add ticks */
      for(i=0;i<(num_chars-1);i++)
	{
	  /* Draw big tick */
	  if((i+1) % freq_big_tick == 0)
	    {
	      draw_box(&r_array,&g_array,&b_array,offset+b_size*i+(block-tick_width_large)/2,offset+b_size*i+(block+tick_width_large)/2,header-tick_height_big,header,0,0,0);
	      draw_box(&r_array,&g_array,&b_array,header-tick_height_big,header,offset+b_size*i+(block-tick_width_large)/2,offset+b_size*i+(block+tick_width_large)/2,0,0,0);
	      draw_box(&r_array,&g_array,&b_array,offset+b_size*i+(block-tick_width_large)/2,offset+b_size*i+(block+tick_width_large)/2,size-header,size-header+tick_height_big,0,0,0);
	      draw_box(&r_array,&g_array,&b_array,size-header,size-header+tick_height_big,offset+b_size*i+(block-tick_width_large)/2,offset+b_size*i+(block+tick_width_large)/2,0,0,0);
	      
	      
	    }
	  else if ((i+1) % freq_small_tick == 0)
	    {
	      draw_box(&r_array,&g_array,&b_array,offset+b_size*i+(block-tick_width_small)/2,offset+b_size*i+(block+tick_width_small)/2,header-tick_height_small,header,0,0,0);
	      
	      draw_box(&r_array,&g_array,&b_array,header-tick_height_small,header,offset+b_size*i+(block-tick_width_small)/2,offset+b_size*i+(block+tick_width_small)/2,0,0,0);
	      draw_box(&r_array,&g_array,&b_array,offset+b_size*i+(block-tick_width_small)/2,offset+b_size*i+(block+tick_width_small)/2,size-header,size-header+tick_height_small,0,0,0);
	      draw_box(&r_array,&g_array,&b_array,size-header,size-header+tick_height_small,offset+b_size*i+(block-tick_width_small)/2,offset+b_size*i+(block+tick_width_small)/2,0,0,0);
	      
	    }
	  
	}
    }
  /* Fill image with values from matrix */
  if(embellish)
    {
      for(i=0;i<num_chars;i++)
	{
	  for(i2=0;i2<block;i2++)
	    {
	      cur_a=b_size*i+offset+i2;
	      
	      for(j=0;j<num_chars;j++)
		{
		  for(j2=0;j2<block;j2++)
		    {
		      cur_b=b_size*j+offset+j2;
		      cur=inc_matrix[i][j];
		      
		      r_array[cur_a][cur_b] = (image_type)palette[cur][RED];
		      g_array[cur_a][cur_b] =  (image_type)palette[cur][GREEN];
		      b_array[cur_a][cur_b] =  (image_type)palette[cur][BLUE];
		      
		      
		      /* Black down diagonal except every x times */
		      if(i==j)
			{
			  
			  if((i+1) % freq_small_tick == 0)
			    {
			      r_array[cur_a][cur_b] = 255;
			      g_array[cur_a][cur_b] = 255;
			      b_array[cur_a][cur_b] = 255;
			    }
			  
			  else
			    {
			      r_array[cur_a][cur_b] = 32;
			      g_array[cur_a][cur_b] = 32;
			      b_array[cur_a][cur_b] = 32;
			    }
			  
			}
		    }
		}
	    }
	}
    }
  else
    {
      for(i=0;i<size;i++)
	for(j=0;j<size;j++)
	  {
	    cur=inc_matrix[i][j];
	    r_array[i][j]= (image_type)palette[cur][RED];
	    g_array[i][j]= (image_type)palette[cur][GREEN];
	    b_array[i][j]= (image_type)palette[cur][BLUE];
	  }
    }
  
  /* Print image */
  out_file=fopen("matrix.ppm","w");
  write_PPM(out_file,r_array,g_array,b_array,size,size);
}
