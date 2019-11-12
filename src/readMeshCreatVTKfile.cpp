
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

    char infile_trias[]="plate-tria-faces.dat";

    ifstream  infile2(infile_trias);

    if(infile2.fail())
    {
      cout << " Could not open the input file" << endl;
      exit(1);
    }

    int  valint[10], nn, ii, jj;

    vector<int> vectemp, nodenums_new(nNode,-1);

    while(infile2 >> val[0] >> val[1] >> val[2] >> val[3] >> val[4] >> val[5] >> val[6] >> val[7])
    {
       //printf("%12.6f \t %12.6f \t %12.6f \n", val[0], val[1], val[2]);

       vectemp.push_back(val[5]-1);
       vectemp.push_back(val[6]-1);
       vectemp.push_back(val[7]-1);
    }

    sort(vectemp.begin(), vectemp.end());
    vectemp.erase(unique(vectemp.begin(), vectemp.end()), vectemp.end());

    nn=vectemp.size();

    for(ii=0; ii<nn; ii++)
    {
      nodenums_new[vectemp[ii]] = 1;
    }

    jj=0;
    for(ii=0; ii<nNode; ii++)
    {
      if(nodenums_new[ii] == 1)
      {
        nodenums_new[ii] = jj++;
      }
    }

    cout << " jj = " << jj << endl;

    vtkSmartPointer<vtkPoints>              pointsPoly  =  vtkSmartPointer<vtkPoints>::New();
    vtkSmartPointer<vtkPoints>               pointsGrid  =  vtkSmartPointer<vtkPoints>::New();
    vtkSmartPointer<vtkTriangle>                triaVTK  =  vtkSmartPointer<vtkTriangle>::New();
    vtkSmartPointer<vtkPolygon>              polygonVTK  =  vtkSmartPointer<vtkPolygon>::New();


    vtkSmartPointer<vtkUnstructuredGrid>         uGrid   =  vtkSmartPointer<vtkUnstructuredGrid>::New();
    vtkSmartPointer<vtkXMLUnstructuredGridWriter>  writeruGrid   =  vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
    vtkSmartPointer<vtkCellArray> polyList  =   vtkSmartPointer<vtkCellArray>::New();
    vtkSmartPointer<vtkPolyData> polyData =   vtkSmartPointer<vtkPolyData>::New();
    vtkSmartPointer<vtkXMLPolyDataWriter>  writerPolyData =     vtkSmartPointer<vtkXMLPolyDataWriter>::New();

    vtkIdType pts[10];

    ofstream  outfile("nodes-new.dat");

    if(outfile.fail())
    {
      cout << " Could not open the input file" << endl;
      exit(1);
    }

    if(outfile.is_open())
    {
      for(ii=0; ii<nn; ii++)
      {
        jj = vectemp[ii];
        outfile << (ii+1) << '\t' << vecCoords[jj][0] << '\t' << vecCoords[jj][1] << '\t' << vecCoords[jj][2] << endl;
        pts[0] = pointsPoly->InsertNextPoint(vecCoords[jj][0], vecCoords[jj][1], vecCoords[jj][2]);
      }
    }

    outfile.close();


    char infile_trias2[]="plate-tria-faces.dat";

    ifstream  infile3(infile_trias2);

    if(infile3.fail())
    {
      cout << " Could not open the input file" << endl;
      exit(1);
    }

    while(infile3 >> val[0] >> val[1] >> val[2] >> val[3] >> val[4] >> val[5] >> val[6] >> val[7])
    {
       printf("%6d \t %6d \t %6d \n", nodenums_new[val[5]-1]+1, nodenums_new[val[6]-1]+1, nodenums_new[val[7]-1]+1);

        polygonVTK->GetPointIds()->SetNumberOfIds(3); //make a triangle

        for(ii=0; ii<3; ii++)
          polygonVTK->GetPointIds()->SetId(ii, nodenums_new[val[5+ii]-1]);

        polyList->InsertNextCell(polygonVTK);
    }

    polyData->SetPoints(pointsPoly);
    polyData->SetPolys(polyList);

    writerPolyData->SetFileName("immersedpoly.vtp");
    writerPolyData->SetInputData(polyData);
    writerPolyData->Write();


  return 0;
}










