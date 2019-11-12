/*
*
* Find cutcells using subroutines in CGAL library
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
*/

#ifndef CUTCELLS_CGAL_H
#define CUTCELLS_CGAL_H

#include "headersVTK.h"
#include "headersBasic.h"
#include <iostream>
#include <boost/foreach.hpp>

#include <CGAL/Simple_cartesian.h>
#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/AABB_triangle_primitive.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>

#include <CGAL/Polyhedron_3.h>
#include <CGAL/boost/graph/graph_traits_Polyhedron_3.h>
#include <CGAL/AABB_face_graph_triangle_primitive.h>
#include <CGAL/Side_of_triangle_mesh.h>
 
#include<CGAL/Polyhedron_incremental_builder_3.h>

typedef CGAL::Simple_cartesian<double> CGAL_K;

typedef CGAL_K::FT CGAL_FT;
typedef CGAL_K::Ray_3 CGAL_Ray;
typedef CGAL_K::Line_3 CGAL_Line;
typedef CGAL_K::Point_3 CGAL_Point;
typedef CGAL_K::Segment_3 CGAL_Segment;
typedef CGAL_K::Triangle_3 CGAL_Triangle;
//typedef std::list<CGAL_Triangle>::iterator CGAL_Iterator;
//typedef CGAL::AABB_triangle_primitive<K, CGAL_Iterator> CGAL_Primitive;
//typedef CGAL::AABB_traits<K, CGAL_Primitive> AABB_triangle_traits;
//typedef CGAL::AABB_tree<AABB_triangle_traits> CGAL_Tree;

typedef CGAL::Polyhedron_3<CGAL_K> CGAL_Polyhedron;
typedef CGAL::AABB_face_graph_triangle_primitive<CGAL_Polyhedron> CGAL_Primitive;
typedef CGAL::AABB_traits<CGAL_K, CGAL_Primitive> CGAL_Traits;
typedef CGAL::AABB_tree<CGAL_Traits> CGAL_Tree;

typedef CGAL::Side_of_triangle_mesh<CGAL_Polyhedron, CGAL_K> CGAL_Point_inside;

typedef CGAL_Polyhedron::HalfedgeDS  CGAL_HalfedgeDS;
using namespace std;



// A modifier creating a triangle with the incremental builder.
template<class HDS>
class polyhedron_builder : public CGAL::Modifier_base<HDS>
{
  public:
    std::vector<CGAL_Point> &vertices;
    std::vector<int>    &faces;
    int  elType;

    polyhedron_builder( std::vector<CGAL_Point> &_vertices, std::vector<int> &_faces, int eT=3 ) 
                : vertices(_vertices), faces(_faces), elType(eT) {}

    void operator()( HDS& hds)
    {
        typedef typename HDS::Vertex    Vertex;
        typedef typename Vertex::Point  Point;
 
        // create a cgal incremental builder
        CGAL::Polyhedron_incremental_builder_3<HDS> B( hds, true);
        
        if(elType == 3)
          B.begin_surface( vertices.size(), faces.size()/3 );
        else
          B.begin_surface( vertices.size(), faces.size()/2 );

       // add the polyhedron vertices
       for(int i=0; i<vertices.size(); i++ )
       {
         B.add_vertex( vertices[i] );
       }

       // add the polyhedron triangles
       if(elType == 3)
       {
         for( int i=0; i<(int)faces.size(); i+=3 )
         {
           B.begin_facet();
           B.add_vertex_to_facet( faces[i+0] );
           B.add_vertex_to_facet( faces[i+1] );
           B.add_vertex_to_facet( faces[i+2] );
           B.end_facet();
         }
       }
       else
       {
         for( int i=0; i<(int)faces.size(); i+=4 )
         {
           B.begin_facet();
           B.add_vertex_to_facet( faces[i+0] );
           B.add_vertex_to_facet( faces[i+1] );
           B.add_vertex_to_facet( faces[i+3] );
           B.end_facet();

           B.begin_facet();
           B.add_vertex_to_facet( faces[i+0] );
           B.add_vertex_to_facet( faces[i+3] );
           B.add_vertex_to_facet( faces[i+2] );
           B.end_facet();
         }
       }

       // finish up the surface
       B.end_surface();
    }
};


/*
Point_inside*  createPointer1(std::vector<Point> &vertices, std::vector<int> &faces, int eT=3)
{
    // build a polyhedron from the loaded arrays
    Polyhedron  polyhedronTemp;
    polyhedron_builder<HalfedgeDS>  poly_builder( vertices, faces, eT );
    polyhedronTemp.delegate( poly_builder );

    Tree tree;
    tree.insert(polyhedronTemp.facets_begin(), polyhedronTemp.facets_end(), polyhedronTemp);

    Point_inside*  inside_tester = new Point_inside(tree);
    
    return inside_tester;
}
*/




inline int findCutcellsCGAL(myGridData& gridData, char* fileNamePoints, char* fileNameFaces)
{
    time_t tstart, tend;

    tstart = time(0); 

    //////////////////////////////////////////////
    //
    // declare vtk variables
    //
    //////////////////////////////////////////////

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

    vtkIdType pts[10];
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

          xx += dx;
          ind++;
        }
        yy += dy;
      }
      zz += dz;
    }


    /*
    ifstream  infile4(infile_trias);

    if(infile4.fail())
    {
      cout << " Could not open the input file" << endl;
      exit(1);
    }

    std::list<Triangle> triangles;

    if(gridData.elType==3)
    {
      while(infile4 >> valInt[0] >> valInt[1] >> valInt[2] >> valInt[3] >> valInt[4])
      {
        triangles.push_back(Triangle(pointsCGAL[valInt[2]-1], pointsCGAL[valInt[3]-1], pointsCGAL[valInt[4]-1]));
      }
    }
    else
    {
      while(infile4 >> valInt[0] >> valInt[1] >> valInt[2] >> valInt[3] >> valInt[4] >> valInt[5])
      {
        triangles.push_back(Triangle(pointsCGAL[valInt[2]-1], pointsCGAL[valInt[3]-1], pointsCGAL[valInt[5]-1]));
        triangles.push_back(Triangle(pointsCGAL[valInt[2]-1], pointsCGAL[valInt[5]-1], pointsCGAL[valInt[4]-1]));
      }
    }
    */

    // create a cgal incremental builder

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
    //vector<double>  points;
    vector<CGAL_Point> points;

    // add the polyhedron points
    while(infileNodes >> val[0] >> val[1] >> val[2] >> val[3])
    {
      //cout << val[0] << '\t' << val[1] << '\t' << val[2] << '\t' << val[3] << endl;
      //points.push_back(val[1]);
      //points.push_back(val[2]);
      //points.push_back(val[3]);

      points.push_back(CGAL_Point(val[1], val[2], val[3]));
    }
   
    int  valInt[10];
    vector<int>  faces;

    // add the polyhedron triangles
    if(gridData.elType == 3)
    {
      while(infileFaces >> valInt[0] >> valInt[1] >> valInt[2] >> valInt[3] >> valInt[4])
      {
        //cout << valInt[0] << '\t' << valInt[1] << '\t' << valInt[2] << '\t' << valInt[3] << endl;
        faces.push_back(valInt[2]-1);
        faces.push_back(valInt[3]-1);
        faces.push_back(valInt[4]-1);
      }
    }
    else
    {
      while(infileFaces >> valInt[0] >> valInt[1] >> valInt[2] >> valInt[3] >> valInt[4] >> valInt[5])
      {
        faces.push_back(valInt[2]-1);
        faces.push_back(valInt[3]-1);
        faces.push_back(valInt[4]-1);
        faces.push_back(valInt[5]-1);
      }
    }

    cout << " points.size()  = " << points.size() << endl;
    cout << " faces.size()   = " << faces.size() << endl;
    
    // build a polyhedron from the loaded arrays
    CGAL_Polyhedron  polyhedronTemp;
    polyhedron_builder<CGAL_HalfedgeDS>  poly_builder( points, faces, gridData.elType );
    polyhedronTemp.delegate( poly_builder );

    // construct AABB tree
    CGAL_Tree treeCGAL;
    treeCGAL.insert(polyhedronTemp.facets_begin(), polyhedronTemp.facets_end(), polyhedronTemp);
    treeCGAL.accelerate_distance_queries();


    //Point_inside inside_tester(treeCGAL);
    CGAL_Point_inside  *point_inside_tester;
    point_inside_tester = new CGAL_Point_inside(treeCGAL);


    //////////////////////////////////////////////
    //
    // create cells/elements
    //
    //////////////////////////////////////////////

    bool flag;
    double ptTemp[3];
    vector<int>  ff(8);
    vtkIdType cellId;

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

            pointsGrid->GetPoint(pts[ll], ptTemp);

            CGAL_Point  ray_begin(ptTemp[0], ptTemp[1], ptTemp[2]);

            //cout << ptTemp[0] << '\t' << ptTemp[1] << '\t' << ptTemp[2] << endl;

            //CGAL::Bounded_side res = inside_tester(ray_begin);
            CGAL::Bounded_side res = (*point_inside_tester)(ray_begin);

            if( (res == CGAL::ON_BOUNDED_SIDE) || (res == CGAL::ON_BOUNDARY) )
            {
              ff[ll] = 1;
            }
            else
            {
              ff[ll] = 0;
            }
          }

          flag = std::equal( ff.begin()+1, ff.end(), ff.begin() );

          if( flag )
          {
            if(ff[0] == 0)
            {
              uGrid2->InsertNextCell(hexVTK->GetCellType(), hexVTK->GetPointIds());
              cellOrientationVTK2->SetTuple1(cellId, ff[0]);
              cellId++;
            }
          }
          else
          {
            uGrid2->InsertNextCell(hexVTK->GetCellType(), hexVTK->GetPointIds());
            cellOrientationVTK2->SetTuple1(cellId, -1);
            cellId++;
          }
        }
      }
    }

    //////////////////////////////////////////////
    //
    // setup and write polyData
    //
    //////////////////////////////////////////////

    char fname2[50];
    sprintf(fname2,"%s%d%s", "cutFEMhex-cgal-",nEx,".vtu");

    uGrid2->SetPoints(pointsGrid);
//    uGrid->GetPointData()->SetScalars(distVTK);
//    uGrid->GetPointData()->AddArray(nodeInOutVTK);

    cout << " uGrid2->GetNumberOfPoints() = " <<  uGrid2->GetNumberOfPoints() << endl;
    cout << " uGrid2->GetNumberOfCells()  = " <<  uGrid2->GetNumberOfCells() << endl;

    uGrid2->GetCellData()->SetScalars(cellOrientationVTK2);

    writeruGrid->SetFileName(fname2);
    writeruGrid->SetInputData(uGrid2);
    writeruGrid->Write();

    tend = time(0); 
    cout << "Time taken for CGAL = "<< difftime(tend, tstart) <<" second(s)."<< endl;

    return 0;
}



#endif


