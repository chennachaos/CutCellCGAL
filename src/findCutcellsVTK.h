/*
*
* Find cutcells using subroutines in VTK library
*
* Copyright (c) 2017, Chennakesava Kadapa (c.kadapa@swansea.ac.uk).
* All rights reserved.
* Date: 12-Nov-2017
*
* This library is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
* The author is not responsible for any damage.
*
*
* =========================================================================

* Create triangulation for the base grid and perform 
* subtriangulation for the cut-cells

* distVTK
*       type    - vtkFloatArray
*       details - stores distance function value at a node of the background mesh

* nodeInOutVTK
*       type    - vtkIntArray
*       details - stores true/false depending upon whether a node is inside/outside
                these values are computed based on 'distVTK' data
* cutCellTypeVTK
*       type    - vtkIntArray
*       details - stores 0/1/2 depending upon whether an element is cut or not. 
                If it is not then which domain (#1 or #2) the element belongs to.
*/

#ifndef CUTCELLS_VTK_H
#define CUTCELLS_VTK_H

#include "headersVTK.h"
#include "headersBasic.h"
#include <iostream>

using namespace std;



inline int findCutcellsVTK(myGridData& gridData, char* fileNamePoints, char* fileNameFaces)
{
    //////////////////////////////////////////////
    //
    // declare vtk variables
    //
    //////////////////////////////////////////////

    time_t tstart, tend;

    tstart = time(0);


    vtkSmartPointer<vtkDataSetMapper>        mapperVTK    =  vtkSmartPointer<vtkDataSetMapper>::New();
    vtkSmartPointer<vtkUnstructuredGrid>     uGridVTK     =  vtkSmartPointer<vtkUnstructuredGrid>::New();
    vtkSmartPointer<vtkPoints>               pointsGrid   =  vtkSmartPointer<vtkPoints>::New();
    vtkSmartPointer<vtkPoints>               pointsVTK2   =  vtkSmartPointer<vtkPoints>::New();
    vtkSmartPointer<vtkVertex>               vertexVTK    =  vtkSmartPointer<vtkVertex>::New();
    vtkSmartPointer<vtkLine>                 lineVTK      =  vtkSmartPointer<vtkLine>::New();
    vtkSmartPointer<vtkQuad>                 quadVTK      =  vtkSmartPointer<vtkQuad>::New();
    vtkSmartPointer<vtkHexahedron>           hexVTK       =  vtkSmartPointer<vtkHexahedron>::New();
    vtkSmartPointer<vtkTriangle>             triaVTK      =  vtkSmartPointer<vtkTriangle>::New();
    vtkSmartPointer<vtkPolygon>              polygonVTK   =  vtkSmartPointer<vtkPolygon>::New();
    vtkSmartPointer<vtkTetra>                tetraVTK     =  vtkSmartPointer<vtkTetra>::New();
    vtkSmartPointer<vtkPyramid>              pyramidVTK   =  vtkSmartPointer<vtkPyramid>::New();
    vtkSmartPointer<vtkWedge>                wedgeVTK     =  vtkSmartPointer<vtkWedge>::New();

    vtkSmartPointer<vtkIntArray>          nodeInOutVTK       =  vtkSmartPointer<vtkIntArray>::New();
    vtkSmartPointer<vtkIntArray>          cutCellTypeVTK       =  vtkSmartPointer<vtkIntArray>::New();
    vtkSmartPointer<vtkIntArray>          cellOrientationVTK    =  vtkSmartPointer<vtkIntArray>::New();
    vtkSmartPointer<vtkIntArray>          cellOrientationVTK2    =  vtkSmartPointer<vtkIntArray>::New();
    vtkSmartPointer<vtkFloatArray>          vecVTK       =  vtkSmartPointer<vtkFloatArray>::New();
    vtkSmartPointer<vtkFloatArray>          vecVTK2       =  vtkSmartPointer<vtkFloatArray>::New();
    vtkSmartPointer<vtkFloatArray>          distVTK       =  vtkSmartPointer<vtkFloatArray>::New();
    vtkSmartPointer<vtkFloatArray>          scaVTK2       =  vtkSmartPointer<vtkFloatArray>::New();
    vtkSmartPointer<vtkFloatArray>          cellDataVTK       =  vtkSmartPointer<vtkFloatArray>::New();
    vtkSmartPointer<vtkFloatArray>          cellDataVTK2       =  vtkSmartPointer<vtkFloatArray>::New();
    vtkSmartPointer<vtkExtractEdges>     extractEdgesVTK  =    vtkSmartPointer<vtkExtractEdges>::New();

    vtkSmartPointer<vtkUnstructuredGrid> uGrid   =  vtkSmartPointer<vtkUnstructuredGrid>::New();
    vtkSmartPointer<vtkUnstructuredGrid> uGrid2  =  vtkSmartPointer<vtkUnstructuredGrid>::New();

    vtkSmartPointer<vtkPoints>              pointsPoly  =  vtkSmartPointer<vtkPoints>::New();

    // Add the polygon to a list of polygons
    vtkSmartPointer<vtkCellArray> polyList  =   vtkSmartPointer<vtkCellArray>::New();

    // Create a PolyData
    vtkSmartPointer<vtkPolyData> polyData =   vtkSmartPointer<vtkPolyData>::New();
    vtkSmartPointer<vtkPolyData> polyData2 =   vtkSmartPointer<vtkPolyData>::New();

    vtkSmartPointer<vtkXMLUnstructuredGridWriter>  writerUGridVTK =  vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
    vtkSmartPointer<vtkXMLPolyDataWriter>  writerPolyData =     vtkSmartPointer<vtkXMLPolyDataWriter>::New();
    vtkSmartPointer<vtkXMLUnstructuredGridWriter>  writeruGrid   =  vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
    vtkSmartPointer<vtkXMLUnstructuredGridWriter>  writeruGrid2  =  vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();

    vtkSmartPointer<vtkMergePoints> mergePoints =     vtkSmartPointer<vtkMergePoints>::New();


    //////////////////////////////////////////////
    //
    // generate base triangulation
    //
    //////////////////////////////////////////////

    uGridVTK->Reset();
    pointsGrid->Reset();
    cellDataVTK->Reset();
    cellDataVTK->Reset();
    cutCellTypeVTK->Reset();
    nodeInOutVTK->Reset();

    int  ii, jj, kk, ll, n1, n2, n3, n4, nn;
    int  ind, ind1, ind2, ind3, ind4, ind5, ind6;

    int nEx = gridData.nEx;
    int nEy = gridData.nEy;
    int nEz = gridData.nEz;

    int nNx = gridData.nNx;
    int nNy = gridData.nNy;
    int nNz = gridData.nNz;

    int nNode = gridData.nNode;
    int nElem = gridData.nElem;

    double x0 = gridData.x0;
    double x1 = gridData.x1;
    double y0 = gridData.y0;
    double y1 = gridData.y1;
    double z0 = gridData.z0;
    double z1 = gridData.z1;

    double dx = (x1-x0)/nEx;
    double dy = (y1-y0)/nEy;
    double dz = (z1-z0)/nEz;

    double xx, yy, zz, fact;

    //gridData.printInfo();

    vtkIdType pts[10];

    distVTK->SetName("dist");
    distVTK->SetNumberOfTuples(nNode);
    nodeInOutVTK->SetName("InOut");
    nodeInOutVTK->SetNumberOfTuples(nNode);
    cutCellTypeVTK->SetName("cellType");
    cutCellTypeVTK->SetNumberOfTuples(nElem);
    cellOrientationVTK->SetName("elType");
    cellOrientationVTK->SetNumberOfTuples(nElem);
    cellOrientationVTK2->SetName("elType");
    cellOrientationVTK2->SetNumberOfTuples(nElem);

    //////////////////////////////////////////////
    //
    // create nodes/points
    //
    //////////////////////////////////////////////

    ind = 0;
    zz = z0;
    for(kk=0; kk<nNz; kk++)
    {
      yy = y0;
      for(jj=0;jj<nNy;jj++)
      {
        xx = x0;
        for(ii=0;ii<nNx;ii++)
        {
          pts[0] = pointsGrid->InsertNextPoint(xx, yy, zz);

          //distVTK->SetTuple1(ind, fact);
          //nodeInOutVTK->SetTuple1(ind, (fact >= 0.0));

          xx += dx;
          ind++;
        }
        yy += dy;
      }
      zz += dz;
    }


    ////////////////////////////////////////////////
    ////////////////////////////////////////////////
    //
    // read the surface and generate polydata
    //
    ////////////////////////////////////////////////
    ////////////////////////////////////////////////

    ifstream  infileNodes(fileNamePoints);
    ifstream  infileFaces(fileNameFaces);

    if(infileNodes.fail())
    {
      cout << " Could not open the input file" << endl;
      exit(1);
    }

    if(infileFaces.fail())
    {
      cout << " Could not open the input file" << endl;
      exit(1);
    }


    double  val[10];
    vtkIdType pt[12];

    while(infileNodes >> val[0] >> val[1] >> val[2] >> val[3])
    {
       //printf("%12.6f \t %12.6f \t %12.6f \n", val[0], val[1], val[2]);

       pt[0] = pointsPoly->InsertNextPoint(val[1], val[2], val[3]);
    }


    int  valInt[10];

    if(gridData.elType == 3)
    {
      while(infileFaces >> valInt[0] >> valInt[1] >> valInt[2] >> valInt[3] >> valInt[4])
      {
        //printf("%12.6f \t %12.6f \t %12.6f \n", valInt[0], valInt[1], valInt[2]);
       
        polygonVTK->GetPointIds()->SetNumberOfIds(3); //make a triangle

        polygonVTK->GetPointIds()->SetId(0, valInt[2]-1);
        polygonVTK->GetPointIds()->SetId(1, valInt[3]-1);
        polygonVTK->GetPointIds()->SetId(2, valInt[4]-1);

        polyList->InsertNextCell(polygonVTK);
      }
    }
    else
    {
      while(infileFaces >> valInt[0] >> valInt[1] >> valInt[2] >> valInt[3] >> valInt[4] >> valInt[5])
      {
        //printf("%12.6f \t %12.6f \t %12.6f \n", valInt[0], valInt[1], valInt[2]);

        polygonVTK->GetPointIds()->SetNumberOfIds(4); //make a quad

        polygonVTK->GetPointIds()->SetId(0, valInt[2]-1);
        polygonVTK->GetPointIds()->SetId(1, valInt[3]-1);
        polygonVTK->GetPointIds()->SetId(2, valInt[5]-1);
        polygonVTK->GetPointIds()->SetId(3, valInt[4]-1);

        polyList->InsertNextCell(polygonVTK);

        /*
        polygonVTK->GetPointIds()->SetNumberOfIds(3); //make a triangle

        polygonVTK->GetPointIds()->SetId(0, valInt[2]-1);
        polygonVTK->GetPointIds()->SetId(1, valInt[3]-1);
        polygonVTK->GetPointIds()->SetId(2, valInt[5]-1);

        polyList->InsertNextCell(polygonVTK);

        polygonVTK->GetPointIds()->SetNumberOfIds(3); //make a triangle

        polygonVTK->GetPointIds()->SetId(0, valInt[2]-1);
        polygonVTK->GetPointIds()->SetId(1, valInt[5]-1);
        polygonVTK->GetPointIds()->SetId(2, valInt[4]-1);

        polyList->InsertNextCell(polygonVTK);
        */
      }
    }

    polyData->SetPoints(pointsPoly);
    polyData->SetPolys(polyList);

    writerPolyData->SetFileName("immersedpoly-vtk.vtp");
    writerPolyData->SetInputData(polyData);
    writerPolyData->Write();

    /////////////////////////////////////////////////////
    /////////////////////////////////////////////////////
    //
    // compute using VTK
    //
    /////////////////////////////////////////////////////
    /////////////////////////////////////////////////////

    vtkSmartPointer<vtkPolyData> pointsPolydataVTK   =   vtkSmartPointer<vtkPolyData>::New();
    pointsPolydataVTK->SetPoints(pointsGrid);

    vtkSmartPointer<vtkSelectEnclosedPoints> selectEnclosedPoints = vtkSmartPointer<vtkSelectEnclosedPoints>::New();

    selectEnclosedPoints->SetInputData(pointsPolydataVTK);
    selectEnclosedPoints->SetTolerance(0.0000001);
    selectEnclosedPoints->Initialize(polyData);


    //////////////////////////////////////////////
    //
    // create cells/elements
    //
    //////////////////////////////////////////////

    double  bnds[6];
    bool flag;
    double ptTemp[3];

    vector<int>  ff(8);

    vtkIdType npts=3, cellId, cellId2, type, id;
    vtkCell  *cellVTK, *cellVTK2;

    cout << " AAAAAAAAAAAAA " << endl;

    nn = nNx*nNy;

    cellId = 0;
    ind=0;
    for(kk=0; kk<nEz; kk++)
    {
      ind5 = nn*kk;
      ind6 = nn*(kk+1);

      for(jj=0; jj<nEy; jj++)
      {
        ind1 = ind5 + nNx*jj;
        ind2 = ind5 + nNx*(jj+1);

        ind3 = ind6 + nNx*jj;
        ind4 = ind6 + nNx*(jj+1);

        for(ii=0;ii<nEx;ii++)
        {
          pts[0] = ind1+ii;          pts[4] = ind3+ii;
          pts[1] = pts[0]+1;         pts[5] = pts[4]+1;
          pts[3] = ind2+ii;          pts[7] = ind4+ii;
          pts[2] = pts[3]+1;         pts[6] = pts[7]+1;

          ff.assign(8,0);

          for(ll=0;ll<8;ll++)
          {
            hexVTK->GetPointIds()->SetId(ll, pts[ll]);
            //cout << ll << '\t' << pts[ll] << endl;
            pointsGrid->GetPoint(pts[ll], ptTemp);
            //ff[ll] = selectEnclosedPoints->IsInside(pts[ll]);
            ff[ll] = selectEnclosedPoints->IsInsideSurface(ptTemp[0], ptTemp[1], ptTemp[2]);
            //cout << ptTemp[0] << '\t' << ptTemp[1] << '\t' << ptTemp[2] << '\t' << ff[ll] << endl;
          }

          flag = std::equal( ff.begin()+1, ff.end(), ff.begin() );

          /*
          uGrid->InsertNextCell(hexVTK->GetCellType(), hexVTK->GetPointIds());
          if( flag && (ff[0] == 0) )
          {
            cellOrientationVTK->SetTuple1(cellId, ff[0]);
          }
          else
          {
            cellOrientationVTK->SetTuple1(cellId, -1);
          }
          cellId++;
          */

          if( flag )
          {
            if(ff[0] == 0)
            {
              uGrid->InsertNextCell(hexVTK->GetCellType(), hexVTK->GetPointIds());
              cellOrientationVTK->SetTuple1(cellId, ff[0]);
              cellId++;
            }
          }
          else
          {
            uGrid->InsertNextCell(hexVTK->GetCellType(), hexVTK->GetPointIds());
            cellOrientationVTK->SetTuple1(cellId, -1);
            cellId++;
          }
        }
      }
    }

    //////////////////////////////////////////////
    //
    // setup and write grid data
    //
    //////////////////////////////////////////////

    char fname1[50];
    sprintf(fname1,"%s%d%s", "cutFEMhex-vtk-",nEx,".vtu");

    uGrid->SetPoints(pointsGrid);
//    uGrid->GetPointData()->SetScalars(distVTK);
//    uGrid->GetPointData()->AddArray(nodeInOutVTK);

    cout << " uGrid->GetNumberOfPoints() = " <<  uGrid->GetNumberOfPoints() << endl;
    cout << " uGrid->GetNumberOfCells()  = " <<  uGrid->GetNumberOfCells() << endl;

    uGrid->GetCellData()->SetScalars(cellOrientationVTK);

    writeruGrid->SetFileName(fname1);
    writeruGrid->SetInputData(uGrid);
    writeruGrid->Write();

    tend = time(0); 
    cout << "Time taken for VTK = "<< difftime(tend, tstart) <<" second(s)."<< endl;

    return 0;
}


#endif

