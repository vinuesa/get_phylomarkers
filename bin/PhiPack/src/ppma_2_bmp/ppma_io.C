# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <fstream>
# include <cmath>
# include <ctime>
# include <cstring>

using namespace std;

# include "ppma_io.H"

# define MAX(a,b) ( (a)>(b) ? (a) : (b) ) 

//****************************************************************************

bool ppma_check_data ( int xsize, int ysize, int maxrgb, int *r,
  int *g, int *b )

//****************************************************************************
//
//  Purpose:
//
//    PPMA_CHECK_DATA checks the data for an ASCII portable pixel map file.
//
//  Discussion:
//
//    XSIZE and YSIZE must be positive, the pointers must not be null,
//    and the data must be nonnegative and no greater than MAXRGB.
//
//  Example:
//
//    P3
//    # feep.ppm
//    4 4
//    15
//     0  0  0    0  0  0    0  0  0   15  0 15
//     0  0  0    0 15  7    0  0  0    0  0  0
//     0  0  0    0  0  0    0 15  7    0  0  0
//    15  0 15    0  0  0    0  0  0    0  0  0
//
//  Modified:
//
//    28 February 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int XSIZE, YSIZE, the number of rows and columns of data.
//
//    Input, int MAXRGB, the maximum RGB value.
//
//    Input, int *R, *G, *B, the arrays of XSIZE by YSIZE data values.
//
//    Output, bool PPMA_CHECK_DATA, is
//    true, if an error was detected, or
//    false, if the data was legal.
//
{
  char c;
  int i;
  int *index;
  int j;
  int k;

  if ( xsize <= 0 )
  {
    cout << "\n";
    cout << "PPMA_CHECK_DATA: Error!\n";
    cout << "  xsize <= 0.\n";
    cout << "  xsize = " << xsize << "\n";
    return true;
  }

  if ( ysize <= 0 )
  {
    cout << "\n";
    cout << "PPMA_CHECK_DATA: Error!\n";
    cout << "  ysize <= 0.\n";
    cout << "  ysize = " << ysize << "\n";
    return true;
  }

  if ( r == NULL )
  {
    cout << "\n";
    cout << "PPMA_CHECK_DATA: Error!\n";
    cout << "  Null pointer to R.\n";
    return true;
  }

  if ( g == NULL )
  {
    cout << "\n";
    cout << "PPMA_CHECK_DATA: Error!\n";
    cout << "  Null pointer to G.\n";
    return true;
  }

  if ( b == NULL )
  {
    cout << "\n";
    cout << "PPMA_CHECK_DATA: Error!\n";
    cout << "  Null pointer to B.\n";
    return true;
  }

  for ( k = 0; k < 3; k++ )
  {

    if ( k == 0 )
    {
      index = r;
      c = 'R';
    }
    else if ( k == 1 )
    {
      index = g;
      c = 'G';
    }
    else if ( k == 2 )
    {
      index = b;
      c = 'B';
    }

    for ( j = 0; j < ysize; j++ )
    {
      for ( i = 0; i < xsize; i++ )
      {
        if ( *index < 0 )
        {
          cout << "\n";
          cout << "PPMA_CHECK_DATA - Fatal error!\n";
          cout << "  Negative data.\n";
          cout << "  " << c << "(" << i << "," << j << ")=" << *index << "\n";
          return true;
        }
        else if ( maxrgb < *index )
        {
          cout << "\n";
          cout << "PPMA_CHECK_DATA - Fatal error!\n";
          cout << "  Data exceeds MAXRGB = " << maxrgb << "\n";
          cout << "  " << c << "(" << i << "," << j << ")=" << *index << "\n";
          return true;
        }

        index = index + 1;
      }
    } 
  }

  return false;
}
//****************************************************************************

bool ppma_example ( int xsize, int ysize, int *r, int *g, int *b )

//****************************************************************************
//
//  Purpose:
//
//    PPMA_EXAMPLE sets up some RGB data.
//
//  Modified:
//
//    28 February 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int XSIZE, YSIZE, the number of rows and columns of data.
//
//    Output, int *R, *G, *B, the arrays of XSIZE by YSIZE RGB values.
//
//    Output, bool PPMA_EXAMPLE, is
//    false, if no error occurred,
//    true, if an error occurred.
//
{
  float f1;
  float f2;
  float f3;
  int i;
  int *indexr;
  int *indexg;
  int *indexb;
  int j;
  float x;
  float y;

  indexr = r;
  indexg = g;
  indexb = b;

  for ( i = 0; i < ysize; i++ )
  {
    y = ( float ) ( ysize + 1 - i ) / ( float ) ( ysize - 1 );
    for ( j = 0; j < xsize; j++ )
    {
      x = ( float ) ( j ) / ( float ) ( xsize - 1 );

      f1 = 4.0 * ( x - 0.5 ) * ( x - 0.5 );
      f2 = sin ( 3.14159265 * x );
      f3 = x;

      if ( y <= f1 )
      {
        *indexr = ( int ) ( 255.0 * f1 );
      }
      else
      {
        *indexr = 50;
      }

      if ( y <= f2 )
      {
        *indexg = ( int ) ( 255.0 * f2 );
      }
      else
      {
        *indexg = 150;
      }

      if ( y <= f3 )
      {
        *indexb = ( int ) ( 255.0 * f3 );
      }
      else
      {
        *indexb = 250;
      }

      indexr = indexr + 1;
      indexg = indexg + 1;
      indexb = indexb + 1;
    }
  }

  return false;
}
//****************************************************************************

bool ppma_read ( char *file_in_name, int *xsize, int *ysize, int *maxrgb,
  int **r, int **g, int **b )

//****************************************************************************
//
//  Purpose:
//
//    PPMA_READ reads the header and data from an ASCII portable pixel map file.
// 
//  Modified:
// 
//    28 February 2003
// 
//  Author:
// 
//    John Burkardt
//
//  Parameters:
//
//    Input, char *FILE_IN_NAME, the name of the file containing the ASCII
//    portable pixel map data.
//
//    Output, int *XSIZE, *YSIZE, the number of rows and columns of data.
//
//    Output, int *MAXRGB, the maximum RGB value.
//
//    Output, int **R, **G, **B, the arrays of XSIZE by YSIZE data values.
//
//    Output, bool PPMA_READ, is
//    true, if an error was detected, or
//    false, if the file was read.
//
{
  bool error;
  ifstream file_in;
  int numbytes;

  file_in.open ( file_in_name );

  if ( !file_in )
  {
    cout << "\n";
    cout << "PPMA_READ - Fatal error!\n";
    cout << "  Cannot open the input file \"" << file_in_name << "\".\n";
    return true;
  }
//
//  Read the header.
//
  error = ppma_read_header ( file_in, xsize, ysize, maxrgb );

  if ( error )
  {
    cout << "\n";
    cout << "PPMA_READ - Fatal error!\n";
    cout << "  PPMA_READ_HEADER failed.\n";
    return true;
  }
//
//  Allocate storage for the data.
//
  numbytes = (*xsize) * (*ysize) * sizeof ( int );

  *r = new int[numbytes];
  *g = new int[numbytes];
  *b = new int[numbytes];
//
//  Read the data.
//
  error = ppma_read_data ( file_in, *xsize, *ysize, *r, *g, *b );

  if ( error )
  {
    cout << "\n";
    cout << "PPMA_READ - Fatal error!\n";
    cout << "  PPMA_READ_DATA failed.\n";
    return true;
  }
//
//  Close the file.
//
  file_in.close ( );

  return false;
}
//****************************************************************************

bool ppma_read_data ( ifstream &file_in, int xsize, int ysize, int *r,
  int *g, int *b )

//****************************************************************************
//
//  Purpose:
//
//    PPMA_READ_DATA reads the data in an ASCII portable pixel map file.
//
//  Modified:
//
//    07 August 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, ifstream &FILE_IN, a pointer to the file containing the ASCII
//    portable pixel map data.
//
//    Input, int XSIZE, YSIZE, the number of rows and columns of data.
//
//    Output, int *R, *G, *B, the arrays of XSIZE by YSIZE data values.
//
//    Output, bool PPMA_READ_DATA, is
//    true, if an error was detected, or
//    false, if the data was read.
//
{
  int i;
  int j;

  for ( j = 0; j < ysize; j++ )
  {
    for ( i = 0; i < xsize; i++ )
    {
      file_in >> *r;
      if ( file_in.eof() )
      {
        return true;
      }
      r = r + 1;

      file_in >> *g;
      if ( file_in.eof() )
      {
        return true;
      }
      g = g + 1;

      file_in >> *b;
      if ( file_in.eof() )
      {
        return true;
      }
      b = b + 1;
    }
  }

  return false;
}
//****************************************************************************

bool ppma_read_header ( ifstream &file_in, int *xsize, int *ysize, int *maxrgb )

//****************************************************************************
//
//  Purpose:
//
//    PPMA_READ_HEADER reads the header of an ASCII portable pixel map file.
//
//  Modified:
//
//    28 February 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, ifstream &FILE_IN, a pointer to the file containing the ASCII
//    portable pixel map data.
//
//    Output, int *XSIZE, *YSIZE, the number of rows and columns of data.
//
//    Output, int *MAXRGB, the maximum RGB value.
//
//    Output, bool PPMA_READ_HEADER, is
//    true, if an error was detected, or
//    false, if the header was read.
//
{
  int count;
  char line[255];
  char *next;
  int step;
  int width;
  char word[255];

  step = 0;

  while ( 1 )
  {
    file_in.getline ( line, sizeof ( line ) );

    if ( file_in.eof() )
    {
      cout << "\n";
      cout << "PPMA_READ_HEADER - Fatal error!\n";
      cout << "  End of file.\n";
      return true;
    }

    next = line;

    if ( line[0] == '#' )
    {
      continue;
    }

    if ( step == 0 )
    {
      count = sscanf ( next, "%s%n", word, &width );
      if ( count == EOF || count ==0 )
      {
        continue;
      }
      next = next + width;
      if ( strcmp ( word, "P3" ) != 0 && strcmp ( word, "p3" ) != 0 )
      {
        cout << "\n";
        cout << "PPMA_READ_HEADER - Fatal error.\n";
        cout << "  Bad magic number = \"" << word << "\".\n";
        return true;
      }
      step = 1;
    }

    if ( step == 1 )
    {

      count = sscanf ( next, "%d%n", xsize, &width );
      next = next + width;
      if ( count == EOF  || count ==0 )
      {
        continue;
      }
      step = 2;
    }

    if ( step == 2 )
    {
      count = sscanf ( next, "%d%n", ysize, &width );
      next = next + width;
      if ( count == EOF || count ==0  )
      {
        continue;
      }
      step = 3;
    }

    if ( step == 3 )
    {
      count = sscanf ( next, "%d%n", maxrgb, &width );
      next = next + width;
      if ( count == EOF  || count ==0 )
      {
        continue;
      }
      break;
    }

  }

  return false;
}
//****************************************************************************

bool ppma_read_test ( char *file_in_name )

//****************************************************************************
//
//  Purpose:
//
//    PPMA_READ_TEST tests the ASCII portable pixel map read routines.
//
//  Modified:
//
//    07 August 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, char *FILE_IN_NAME, the name of the file containing the ASCII
//    portable pixel map data.
//
//    Output, bool PPMA_READ_TEST, is
//    true, if an error was detected, or
//    false, if the test was carried out.
//
{
  int *b;
  bool error;
  int *g;
  int maxrgb;
  int *r;
  int xsize;
  int ysize;

  r = NULL;
  g = NULL;
  b = NULL;
//
//  Read the data.
//
  error = ppma_read ( file_in_name, &xsize, &ysize, &maxrgb, &r, &g, &b );

  if ( error )
  {
    cout << "\n";
    cout << "PPMA_READ_TEST - Fatal error!\n";
    cout << "  PPMA_READ failed.\n";

    delete [] r;
    delete [] g;
    delete [] b;

    return true;
  }

  cout << "\n";
  cout << "PPMA_READ_TEST:\n";
  cout << "  PPMA_READ was able to read \"" << file_in_name << "\".\n";
//
//  Check the data.
//
  error = ppma_check_data ( xsize, ysize, maxrgb, r, g, b );

  delete [] r;
  delete [] g;
  delete [] b;

  if ( error )
  {
    cout << "\n";
    cout << "PPMA_READ_TEST - Fatal error!\n";
    cout << "  PPMA_CHECK_DATA reports bad data in the file.\n";
    return true;
  }

  cout << "\n";
  cout << "PPMA_READ_TEST:\n";
  cout << "  PPMA_CHECK_DATA has approved the data from the file.\n";

  return false;
}
//****************************************************************************

bool ppma_write ( char *file_out_name, int xsize, int ysize, int *r, 
  int *g, int *b )

//****************************************************************************
//
//  Purpose:
//
//    PPMA_WRITE writes the header and data for an ASCII portable pixel map file.
// 
//  Example:
//
//    P3
//    # feep.ppm
//    4 4
//    15
//     0  0  0    0  0  0    0  0  0   15  0 15
//     0  0  0    0 15  7    0  0  0    0  0  0
//     0  0  0    0  0  0    0 15  7    0  0  0
//    15  0 15    0  0  0    0  0  0    0  0  0
//
//  Modified:
// 
//    28 February 2003
// 
//  Author:
// 
//    John Burkardt
//
//  Parameters:
//
//    Input, char *FILE_OUT_NAME, the name of the file to contain the ASCII
//    portable pixel map data.
//
//    Input, int XSIZE, YSIZE, the number of rows and columns of data.
//
//    Input, int *R, *G, *B, the arrays of XSIZE by YSIZE data values.
//
//    Output, bool PPMA_WRITE, is
//    true, if an error was detected, or
//    false, if the file was written.
//
{
  bool error;
  ofstream file_out;
  int i;
  int *indexb;
  int *indexg;
  int *indexr;
  int j;
  int maxrgb;
//
//  Open the output file.
//
  file_out.open ( file_out_name );

  if ( !file_out )
  {
    cout << "\n";
    cout << "PPMA_WRITE - Fatal error!\n";
    cout << "  Cannot open the output file \"" << file_out_name << "\".\n";
    return true;
  }
//
//  Compute the maximum.
//
  maxrgb = 0;
  indexr = r;
  indexg = g;
  indexb = b;

  for ( j = 0; j < ysize; j++ )
  {
    for ( i = 0; i < xsize; i++ )
    {
      maxrgb = MAX ( maxrgb, *indexr );
      maxrgb = MAX ( maxrgb, *indexg );
      maxrgb = MAX ( maxrgb, *indexb );
      indexr = indexr + 1;
      indexg = indexg + 1;
      indexb = indexb + 1;
    }
  }
//
//  Write the header.
//
  error = ppma_write_header ( file_out, file_out_name, xsize, ysize, maxrgb );

  if ( error )
  {
    cout << "\n";
    cout << "PPMA_WRITE - Fatal error!\n";
    cout << "  PPMA_WRITE_HEADER failed.\n";
    return true;
  }
//
//  Write the data.
//
  error = ppma_write_data ( file_out, xsize, ysize, r, g, b );

  if ( error )
  {
    cout << "\n";
    cout << "PPMA_WRITE - Fatal error!\n";
    cout << "  PPMA_WRITE_DATA failed.\n";
    return true;
  }
//
//  Close the file.
//
  file_out.close ( );

  return false;
}
//****************************************************************************

bool ppma_write_data ( ofstream &file_out, int xsize, int ysize, int *r,
  int *g, int *b )

//****************************************************************************
//
//  Purpose:
//
//    PPMA_WRITE_DATA writes the data for an ASCII portable pixel map file.
//
//  Modified:
//
//    28 February 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, ofstream &FILE_OUT, a pointer to the file to contain the ASCII
//    portable pixel map data.
//
//    Input, int XSIZE, YSIZE, the number of rows and columns of data.
//
//    Input, int *R, *G, *B, the arrays of XSIZE by YSIZE data values.
//
//    Output, bool PPMA_WRITE_DATA, is
//    true, if an error was detected, or
//    false, if the data was written.
//
{
  int i;
  int *indexb;
  int *indexg;
  int *indexr;
  int j;
  int numval;

  indexr = r;
  indexg = g;
  indexb = b;
  numval = 0;

  for ( j = 0; j < ysize; j++ )
  {
    for ( i = 0; i < xsize; i++ )
    {
      file_out << *indexr << " " << *indexg << " " << *indexb;
      numval = numval + 3;
      indexr = indexr + 1;
      indexg = indexg + 1;
      indexb = indexb + 1;

      if ( numval % 12 == 0 || i == xsize - 1 || numval == 3 * xsize * ysize )
      {
        file_out << "\n";
      }
      else
      {
        file_out << " ";
      }

    }
  }
  return false;
}
//****************************************************************************

bool ppma_write_header ( ofstream &file_out, char *file_out_name, int xsize, 
  int ysize, int maxrgb )

//****************************************************************************
//
//  Purpose:
//
//    PPMA_WRITE_HEADER writes the header of an ASCII portable pixel map file.
//
//  Modified:
//
//    28 February 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, ofstream &FILE_OUT, a pointer to the file to contain the ASCII
//    portable pixel map data.
//
//    Input, char *FILE_OUT_NAME, the name of the file.
//
//    Input, int XSIZE, YSIZE, the number of rows and columns of data.
//
//    Input, int MAXRGB, the maximum RGB value.
//
//    Output, bool PPMA_WRITE_HEADER, is
//    true, if an error was detected, or
//    false, if the header was written.
//
{
  file_out << "P3\n";
  file_out << "# " << file_out_name << " created by PPMA_IO::PPMA_WRITE.\n";
  file_out << xsize << "  " << ysize << "\n";
  file_out << maxrgb << "\n";

  return false;
}
//****************************************************************************

bool ppma_write_test ( char *file_out_name )

//****************************************************************************
//
//  Purpose:
//
//    PPMA_WRITE_TEST tests the ASCII portable pixel map write routines.
//
//  Modified:
//
//    28 February 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, char *FILE_OUT_NAME, the name of the file to contain the ASCII
//    portable pixel map data.
//
//    Output, bool PPMA_WRITE_TEST, equals
//    true, if the test could not be carried out,
//    false, if the test was carried out.
//
{
  int *b;
  bool error;
  int *g;
  int maxrgb;
  int *r;
  int xsize;
  int ysize;

  xsize = 300;
  ysize = 300;
//
//  Allocate memory.
//
  r = new int[xsize * ysize];
  g = new int[xsize * ysize];
  b = new int[xsize * ysize];
//
//  Set the data.
//
  error = ppma_example ( xsize, ysize, r, g, b );

  if ( error )
  {
    cout << "\n";
    cout << "PPMA_WRITE_TEST - Fatal error!\n";
    cout << "  PPMA_EXAMPLE failed.\n";
    return true;
  }
//
//  Write the data to the file.
//
  error = ppma_write ( file_out_name, xsize, ysize, r, g, b );

  delete [] r;
  delete [] g;
  delete [] b;
 
  if ( error )
  {
    cout << "\n";
    cout << "PPMA_WRITE_TEST - Fatal error!\n";
    cout << "  PPMA_WRITE failed.\n";
    return true;
  }

  return false;
}
