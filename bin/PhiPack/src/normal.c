#include "global.h"
#include "math.h"
double normal_01_cdf ( double x )

//******************************************************************************
//
//  Purpose:
//
//    NORMAL_01_CDF evaluates the Normal 01 CDF.
//
//  Modified:
//
//    10 February 1999
//
//  Author:
//
//    John Burkardt
//
//  Reference: 
//
//    A G Adams,
//    Areas Under the Normal Curve,
//    Algorithm 39, 
//    Computer j., 
//    Volume 12, pages 197-198, 1969.
//
//  Parameters:
//
//    Input, double X, the argument of the CDF.
//
//    Output, double CDF, the value of the CDF.
//
{
  double a1 = 0.398942280444E+00;
  double a2 = 0.399903438504E+00;
  double a3 = 5.75885480458E+00;
  double a4 = 29.8213557808E+00;
  double a5 = 2.62433121679E+00;
  double a6 = 48.6959930692E+00;
  double a7 = 5.92885724438E+00;
  double b0 = 0.398942280385E+00;
  double b1 = 3.8052E-08;
  double b2 = 1.00000615302E+00;
  double b3 = 3.98064794E-04;
  double b4 = 1.98615381364E+00;
  double b5 = 0.151679116635E+00;
  double b6 = 5.29330324926E+00;
  double b7 = 4.8385912808E+00;
  double b8 = 15.1508972451E+00;
  double b9 = 0.742380924027E+00;
  double b10 = 30.789933034E+00;
  double b11 = 3.99019417011E+00;
  double cdf;
  double q;
  double y;
//
//  |X| <= 1.28.
//
  if ( fabs ( x ) <= 1.28 )
  {
    y = 0.5 * x * x;

    q = 0.5 - fabs ( x ) * ( a1 - a2 * y / ( y + a3 - a4 / ( y + a5 
      + a6 / ( y + a7 ) ) ) );
//
//  1.28 < |X| <= 12.7
//
  }
  else if ( fabs ( x ) <= 12.7 )
  {
    y = 0.5 * x * x;

    q = exp ( - y ) * b0 / ( fabs ( x ) - b1 
      + b2 / ( fabs ( x ) + b3 
      + b4 / ( fabs ( x ) - b5 
      + b6 / ( fabs ( x ) + b7 
      - b8 / ( fabs ( x ) + b9 
      + b10 / ( fabs ( x ) + b11 ) ) ) ) ) );
//
//  12.7 < |X|
//
  }
  else
  {
    q = 0.0;
  }
//
//  Take account of negative X.
//
  if ( x < 0.0 )
  {
    cdf = q;
  }
  else
  {
    cdf = 1.0 - q;
  }

  return cdf;
}
//**********************************************************************
