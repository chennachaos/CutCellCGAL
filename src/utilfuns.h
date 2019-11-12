
#include "headersBasic.h"

#include <algorithm>

using namespace std;


inline double AreaTriangle(double* pt1, double* pt2, double *pt3)
{
  return 0.5*(pt1[0]*(pt2[1]-pt3[1]) + pt2[0]*(pt3[1]-pt1[1]) + pt3[0]*(pt1[1]-pt2[1]));
}


inline double AreaQuad(double* pt1, double* pt2, double *pt3, double* pt4)
{
   return  0.5*( (pt1[0]-pt4[0])*(pt2[1]-pt3[1]) + (pt2[0]-pt3[0])*(pt4[1]-pt1[1]) );
}



