
#include "headersVTK.h"
#include "headersBasic.h"
#include "utilfuns.h"
#include <algorithm>
#include <iostream>
#include <list>

using namespace std;

typedef double REAL;



int main(int argc, char* argv[])
{
    char infile_coords[]="platethickshort-nodes-coords.dat";

    ifstream  infile(infile_coords);

    if(infile.fail())
    {
       cout << " Could not open the input file" << endl;
      exit(1);
    }

    double  val[10];

    vector<vector<double> > vecCoords;
    vector<double> point3d(3);

    while(infile >> val[0] >> val[1] >> val[2] >> val[3])
    {
       //printf("%12.6f \t %12.6f \t %12.6f \n", val[0], val[1], val[2]);

       point3d[0] = val[1];
       point3d[1] = val[2];
       point3d[2] = val[3];

       vecCoords.push_back(point3d);
    }

    int nNode = vecCoords.size();

    cout << " size = " << nNode << endl;

    char infile_fixed[]="platethickshort-nodes-fixed.dat";

    ifstream  infile2(infile_fixed);

    if(infile2.fail())
    {
      cout << " Could not open the input file" << endl;
      exit(1);
    }

    int  valint[10], nn, ii;

    vector<int> vectemp(nNode,0);

    while(infile2 >> val[0] >> val[1] >> val[2] >> val[3] >> val[4])
    {
       //printf("%12.6f \t %12.6f \t %12.6f \n", val[0], val[1], val[2]);

       //vectemp.push_back(valint[0]);

       nn = (int) val[0];
       nn -= 1;

       vectemp[nn] = 1;
    }

    ofstream  outfile("nodes-new.dat");

    if(outfile.fail())
    {
      cout << " Could not open the input file" << endl;
      exit(1);
    }

    if(outfile.is_open())
    {
      for(ii=0; ii<nNode; ii++)
      {
        if(vectemp[ii] == 1)
        {
          outfile << ii << '\t' << 1 << '\t' << 1 << '\t' << 1 << '\t' << 1 << '\t' << 1 << '\t' << 1 << '\t' << 0 ;
          outfile << '\t' << vecCoords[ii][0] << '\t' << vecCoords[ii][1] << '\t' << vecCoords[ii][2] << endl;
        }
        else
        {
          outfile << ii << '\t' << 0 << '\t' << 0 << '\t' << 0 << '\t' << 0 << '\t' << 0 << '\t' << 0 << '\t' << 0;
          outfile << '\t' << vecCoords[ii][0] << '\t' << vecCoords[ii][1] << '\t' << vecCoords[ii][2] << endl;
        }
      }
    }

    outfile.close();

    //sort(vectemp.begin(), vectemp.end());
    //vectemp.erase(unique(vectemp.begin(), vectemp.end()), vectemp.end());

  return 0;
}










