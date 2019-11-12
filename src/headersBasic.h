//#define EIGEN_YES_I_KNOW_SPARSE_MODULE_IS_NOT_STABLE_YET
//#define EIGEN_SUPERLU_SUPPORT


#ifndef incl_headersBasic_h
#define incl_headersBasic_h


#include <iostream>
#include <limits.h>
#include <float.h>
#include <vector>
#include <assert.h>
#include <math.h>
#include <cmath>
#include <stdio.h>
#include <set>
#include <fstream>
#include <string.h>
#include <iomanip>
#include <stdlib.h>
#include <algorithm>
#include <numeric>
#include <iterator>
#include <functional>
#include <list>
#include <limits>


struct myGridData
{
    int elType;
    int nNode, nElem;
    int nNx, nNy, nNz;
    int nEx, nEy, nEz;
    double x0, x1, y0, y1, z0, z1;

    printInfo()
    {
      cout << " Grid data... " << endl;
      cout << " ............ " << endl;
      cout << " Dimensions " << endl;
      cout << " ............ " << endl;
      cout << x0 << '\t' << x1 << endl;
      cout << y0 << '\t' << y1 << endl;
      cout << z0 << '\t' << z1 << endl;
      cout << " ............ " << endl;
      cout << " nEx = " << nEx << endl;
      cout << " nEy = " << nEy << endl;
      cout << " nEz = " << nEz << endl;
      cout << " ............ " << endl;
      cout << " nNx = " << nNx << endl;
      cout << " nNy = " << nNy << endl;
      cout << " nNz = " << nNz << endl;
      cout << " ............ " << endl;
      cout << " nElem = " << nElem << endl;
      cout << " nNode = " << nNode << endl;
      cout << " ............ " << endl;
      cout << endl; cout << endl; cout << endl;
    }
};



#endif
