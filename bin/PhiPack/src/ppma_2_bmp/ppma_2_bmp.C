# include <cstdlib> 
# include <iostream>
# include <iomanip>
# include <fstream>
# include <ctime>
# include <cstring>

using namespace std;

# include "bmp_io.H"
# include "ppma_io.H"

int main ( int argc, char *argv[] );
void imat_vert_flip ( int xsize, int ysize, int *array );
void int_2_u_char ( int xsize, int ysize, int *a_int, unsigned char *a_u );
bool ppma_2_bmp ( char *file_in_name, char *file_out_name );
void timestamp ( void );

//****************************************************************************

int main ( int argc, char *argv[] )

//****************************************************************************
//
//  Purpose:
//
//    PPMA_2_BMP converts an ASCII PPM file to a BMP file.
//
//  Discussion:
//
//    The program requires the BMP_IO and PPMA_IO libraries.
//
//  Modified:
//
//    09 April 2003
//
//  Author:
//
//    John Burkardt
//
//  Usage:
//
//    ppma_2_bmp file.ppma file.bmp
//
//  Parameters:
//
//    FILE.PPMA is the name of the input ASCII PPM file to be read.
//
//    FILE.BMP is the name of the output BMP file to be created.
//
{
  bool error;
  char file_in_name[80];
  char file_out_name[80];

  timestamp ( );

  cout << "\n";
  cout << "PPMA_2_BMP:\n";
  cout << "\n";
  cout << "  Compiled on " << __DATE__ << " at " << __TIME__ << ".\n";
  cout << "\n";
  cout << "  Convert a NetPBM ASCII PPM file to Microsoft BMP format.\n";
//
//  Get the specification for the input file.
//
  if ( argc < 2 ) 
  {
    cout << "\n";
    cout << "PPMA_2_BMP:\n";
    cout << "  Please enter the input PPMA file name:\n";
    
    cin.getline ( file_in_name, sizeof ( file_in_name ) );

    error = cin.eof();

    if ( error )
    {
      exit ( 1 );
    }
  }
  else 
  {
    strcpy ( file_in_name, argv[1] );
  }
//
//  Get the specification for the output file.
//
  if ( argc < 3 ) 
  {
    cout << "\n";
    cout << "PPMA_2_BMP:\n";
    cout << "  Please enter the output BMP file name:\n";

    cin.getline ( file_out_name, sizeof ( file_out_name ) );

    error = cin.eof();

    if ( error )
    {
      exit ( 1 );
    }
  }
  else 
  {
    strcpy ( file_out_name, argv[2] );
  }

  error = ppma_2_bmp ( file_in_name, file_out_name );

  if ( error )
  {
    cout << "\n";
    cout << "PPMA_2_BMP - Fatal error!\n";
    cout << "  The conversion was not successful.\n";
    exit ( 1 );
  }

  cout << "\n";
  cout << "PPMA_2_BMP:\n";
  cout << "  Normal end of execution.\n";

  cout << "\n";
  timestamp ( );

  return 0;
}
//****************************************************************************

bool ppma_2_bmp ( char *file_in_name, char *file_out_name )

//****************************************************************************
//
//  Purpose:
//
//    PPMA_2_BMP reads a PPMA (ASCII PPM) file and writes a BMP file.
//
//  Modified:
//
//    06 April 2001
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, char *FILE_IN_NAME, the name of the input PPMA file.
//
//    Input, char *FILE_OUT_NAME, the name of the output BMP file.
//
//    Output, bool ppma_2_bmp, is true if an error occurred.
//
{
  int *b_int;
  unsigned char *b_u_char;
  bool error;
  int *g_int;
  unsigned char *g_u_char;
  long int height;
  int maxrgb;
  int *r_int;
  unsigned char *r_u_char;
  unsigned long int width;
  int xsize;
  int ysize;

  r_int = NULL;
  g_int = NULL;
  b_int = NULL;
//
//  Read the integer data from the PPMA file.
//
  error = ppma_read ( file_in_name, &xsize, &ysize, &maxrgb,
    &r_int, &g_int, &b_int );

  if ( error )
  {
    cout << "\n";
    cout << "PPMA_2_BMP - Fatal error!\n";
    cout << "  PPMA_READ failed.\n";
    return error;
  }

  cout << "\n";
  cout << "PPMA_2_BMP:\n";
  cout << "  XSIZE = " << xsize << ".\n";
  cout << "  YSIZE = " << ysize << ".\n";
//
//  The BMP up-down orientation is the opposite of what it is in PPM.
//
  cout << "\n";
  cout << "PPMA_2_BMP:\n";
  cout << "  Flipping data orientation.\n";

  imat_vert_flip ( xsize, ysize, r_int );
  imat_vert_flip ( xsize, ysize, g_int );
  imat_vert_flip ( xsize, ysize, b_int );
//
//  BMP expects unsigned characters, not integers.
//
  r_u_char = new unsigned char [xsize*ysize];
  g_u_char = new unsigned char [xsize*ysize];
  b_u_char = new unsigned char [xsize*ysize];

  int_2_u_char ( xsize, ysize, r_int, r_u_char );
  int_2_u_char ( xsize, ysize, g_int, g_u_char );
  int_2_u_char ( xsize, ysize, b_int, b_u_char );

  delete [] r_int;
  delete [] g_int;
  delete [] b_int;

  height = ( long int ) ysize;
  width = ( unsigned long int ) xsize;

  error = bmp_24_write ( file_out_name, width, height, 
    r_u_char, g_u_char, b_u_char );

  if ( error )
  {
    cout << "\n";
    cout << "PPMA_2_BMP - Fatal error!\n";
    cout << "  BMP_24_WRITE failed.\n";
    return error;
  }

  delete [] r_u_char;
  delete [] g_u_char;
  delete [] b_u_char;

  cout << "\n";
  cout << "PPMA_2_BMP:\n";
  cout << "  Normal end of translation.\n";

  return false;
}
//****************************************************************************

void imat_vert_flip ( int xsize, int ysize, int *array )

//****************************************************************************

//
//  Purpose:
//
//    IMAT_VERT_FLIP swaps rows of a 2D array, to flip it vertically.
//
//  Modified:
//
//    07 April 2001
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int XSIZE, int YSIZE, the number of columns and rows.
//
//    Input, int *ARRAY, the address of the first element of the array.
//
{
  int i;
  int j;
  int k1;
  int k2;
  int temp;

  for ( j = 0; j <= ( ysize / 2 ); j++ ) 
  {
    k1 = xsize * j;
    k2 = xsize * ( ysize - 1 - j );

    for ( i = 0; i < xsize; i++ ) 
    {
      temp        = array[k1+i];
      array[k1+i] = array[k2+i];
      array[k2+i] = temp;
    }
  }

  return;
}
//****************************************************************************

void int_2_u_char ( int xsize, int ysize, int *a_int, unsigned char *a_u )

//****************************************************************************
//
//  Purpose:
//
//    INT_2_U_CHAR converts a 2D array of ints to unsigned chars.
//
//  Modified:
//
//    09 April 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int XSIZE, int YSIZE, the number of columns and rows.
//
//    Input, int char *A_INT, the array of ints.
//
//    Output, unsigned char *A_U, the array of unsigned chars.
//
{
  int *b_int;
  unsigned char *b_u;
  int i;
  int j;

  b_u = a_u;
  b_int = a_int;
 
  for ( i = 0; i < xsize; i++ )
  {
    for ( j = 0; j < ysize; j++ )
    {
      *b_u = ( unsigned char ) *b_int;
      b_int = b_int + 1;
      b_u = b_u + 1;
    }
  }

  return;
}
//**********************************************************************

void timestamp ( void )

//**********************************************************************
//
//  Purpose:
//
//    TIMESTAMP prints the current YMDHMS date as a time stamp.
//
//  Example:
//
//    May 31 2001 09:45:54 AM
//
//  Modified:
//
//    24 September 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    None
//
{
#define TIME_SIZE 40

  static char time_buffer[TIME_SIZE];
  const struct tm *tm;
  size_t len;
  time_t now;

  now = time ( NULL );
  tm = localtime ( &now );

  len = strftime ( time_buffer, TIME_SIZE, "%d %B %Y %I:%M:%S %p", tm );

  cout << time_buffer << "\n";

  return;
#undef TIME_SIZE
}
