/*
* Copyright (c) 2017, Chennakesava Kadapa (c.kadapa@swansea.ac.uk).
* All rights reserved.
* Date: 17-July-2017
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

* =========================================================================
*/

#include "headersVTK.h"
#include "headersBasic.h"
#include "utilfuns.h"

#include <algorithm>
#include <iostream>
#include <list>
#include <vector>
#include <fstream>
#include <limits>
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
 


typedef CGAL::Simple_cartesian<double> K;
//typedef CGAL::Exact_predicates_exact_constructions_kernel K;

typedef K::FT FT;
typedef K::Ray_3 Ray;
typedef K::Line_3 Line;
typedef K::Point_3 Point;
typedef K::Segment_3 Segment;
typedef K::Triangle_3 Triangle;
//typedef std::list<Triangle>::iterator Iterator;
//typedef CGAL::AABB_triangle_primitive<K, Iterator> Primitive;
//typedef CGAL::AABB_traits<K, Primitive> AABB_triangle_traits;
//typedef CGAL::AABB_tree<AABB_triangle_traits> Tree;

typedef CGAL::Polyhedron_3<K> Polyhedron;
typedef CGAL::AABB_face_graph_triangle_primitive<Polyhedron> Primitive;
typedef CGAL::AABB_traits<K, Primitive> Traits;
typedef CGAL::AABB_tree<Traits> Tree;

typedef CGAL::Side_of_triangle_mesh<Polyhedron, K> Point_inside;

typedef Polyhedron::HalfedgeDS  HalfedgeDS;


// A modifier creating a triangle with the incremental builder.
template<class HDS>
class polyhedron_builder : public CGAL::Modifier_base<HDS>
{
  public:
    std::vector<Point> &vertices;
    std::vector<int>    &faces;
    int  elType;

    polyhedron_builder( std::vector<Point> &_coords, std::vector<int> &_tris, int eT=3 ) : vertices(_coords), faces(_tris), elType(eT) {}
    
    void operator()( HDS& hds)
    {
        typedef typename HDS::Vertex   Vertex;
        typedef typename Vertex::Point Point;
 
        // create a cgal incremental builder
        CGAL::Polyhedron_incremental_builder_3<HDS> B( hds, true);
        
        if(elType == 3)
          B.begin_surface( vertices.size(), faces.size()/3 );
        else
          B.begin_surface( vertices.size(), faces.size()/4 );

       // add the polyhedron vertices
       for( int i=0; i<(int)vertices.size(); i++ )
       {
         B.add_vertex(vertices[i]);
         //B.add_vertex( Point( coords[i+0], coords[i+1], coords[i+2] ) );
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


Point_inside*  createPointer2(Polyhedron& polyhedronTemp)
{
    Tree tree;
    tree.insert(polyhedronTemp.facets_begin(), polyhedronTemp.facets_end(), polyhedronTemp);

    Point_inside*  inside_tester = new Point_inside(tree);
    
    return inside_tester;
}


inline Point_inside*  createPointer3(Tree& tree)
{
    Point_inside*  inside_tester = new Point_inside(tree);
    
    return inside_tester;
}



using namespace std;


int main(int argc, char* argv[])
{
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

    double  x0, x1, y0, y1, z0, z1;

    int nEx, nEy, nEz, nNx, nNy, nNz, nNode, nElem;

    x0 = -1.6;    x1 =  1.6;    nEx = 50;
    y0 = -1.6;    y1 =  1.6;    nEy = 50;
    z0 = -1.6;    z1 =  1.6;    nEz = 50;

/*
    if(argc == 0)
      cerr << " Input data " << endl;
    else
    {
       x0  = atof(argv[1]);
       x1  = atof(argv[2]);
       nEx = atoi(argv[3]);

       y0  = atof(argv[4]);
       y1  = atof(argv[5]);
       nEy = atoi(argv[6]);

       z0  = atof(argv[7]);
       z1  = atof(argv[8]);
       nEz = atoi(argv[9]);
    }
*/
    nElem = nEx*nEy*nEz;

    nNx = nEx+1;
    nNy = nEy+1;
    nNz = nEz+1;

    nNode = nNx*nNy*nNz;


    uGridVTK->Reset();
    pointsGrid->Reset();
    cellDataVTK->Reset();
    cellDataVTK->Reset();
    cutCellTypeVTK->Reset();
    nodeInOutVTK->Reset();

    int  ii, jj, kk, ll, n1, n2, n3, n4, nn;
    int  ind, ind1, ind2, ind3, ind4, ind5, ind6;

    double dx = (x1-x0)/nEx;
    double dy = (y1-y0)/nEy;
    double dz = (z1-z0)/nEz;
    double xx, yy, zz, fact;

    cout << x0 << '\t' << x1 << '\t' << dx << endl;
    cout << y0 << '\t' << y1 << '\t' << dy << endl;
    cout << z0 << '\t' << z1 << '\t' << dz << endl;


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

    //char infile_nodes[]="sphere-nodes.dat";
    //char infile_trias[]="sphere-trias.dat";

    //char infile_nodes[]="reliefvalve3Dblock-nodes.dat";
    //char infile_trias[]="reliefvalve3Dblock-trias.dat";

    //char infile_nodes[]="turekbeam3d-nodes.dat";
    //char infile_trias[]="turekbeam3d-quads.dat";

    char infile_nodes[]="thickplate-nodes.dat";
    char infile_trias[]="thickplate-elems.dat";


    ifstream  infile(infile_nodes);

    if(infile.fail())
    {
      cout << " Could not open the input file" << endl;
      exit(1);
    }


    double  val[10];

    vtkIdType pt[12];

    while(infile >> val[0] >> val[1] >> val[2] >> val[3])
    {
       //printf("%12.6f \t %12.6f \t %12.6f \n", val[0], val[1], val[2]);

       pt[0] = pointsPoly->InsertNextPoint(val[1], val[2], val[3]);
    }

    ifstream  infile2(infile_trias);

    if(infile2.fail())
    {
       cout << " Could not open the input file" << endl;
      exit(1);
    }

    int  valInt[10], ETYPE=4;

    if(ETYPE==3)
    {
      while(infile2 >> valInt[0] >> valInt[1] >> valInt[2] >> valInt[3] >> valInt[4])
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
      while(infile2 >> valInt[0] >> valInt[1] >> valInt[2] >> valInt[3] >> valInt[4] >> valInt[5])
      {
        //printf("%12.6f \t %12.6f \t %12.6f \n", valInt[0], valInt[1], valInt[2]);
       
        //polygonVTK->GetPointIds()->SetNumberOfIds(4); //make a quad

        //polygonVTK->GetPointIds()->SetId(0, valInt[2]-1);
        //polygonVTK->GetPointIds()->SetId(1, valInt[3]-1);
        //polygonVTK->GetPointIds()->SetId(2, valInt[5]-1);
        //polygonVTK->GetPointIds()->SetId(3, valInt[4]-1);

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
      }
    }

    polyData->SetPoints(pointsPoly);
    polyData->SetPolys(polyList);

    writerPolyData->SetFileName("immersedpoly.vtp");
    writerPolyData->SetInputData(polyData);
    writerPolyData->Write();

    /////////////////////////////////////////////////////
    /////////////////////////////////////////////////////
    //
    // compute using VTK
    //
    /////////////////////////////////////////////////////
    /////////////////////////////////////////////////////

    time_t tstart, tend;

    tstart = time(0);


    vtkSmartPointer<vtkPolyData> pointsPolydataVTK   =   vtkSmartPointer<vtkPolyData>::New();
    pointsPolydataVTK->SetPoints(pointsGrid);

    vtkSmartPointer<vtkSelectEnclosedPoints> selectEnclosedPoints = vtkSmartPointer<vtkSelectEnclosedPoints>::New();

    selectEnclosedPoints->SetInputData(pointsPolydataVTK);
    selectEnclosedPoints->SetTolerance(0.00001);
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

          uGrid->InsertNextCell(hexVTK->GetCellType(), hexVTK->GetPointIds());

          flag = std::equal( ff.begin()+1, ff.end(), ff.begin() );

          if( flag )
          {
            cellOrientationVTK->SetTuple1(cellId, ff[0]);
          }
          else
          {
            cellOrientationVTK->SetTuple1(cellId, -1);
          }

          cellId++;
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
    cout << "Time taken = "<< difftime(tend, tstart) <<" second(s)."<< endl;

    //////////////////////////////////////////////////////
    //////////////////////////////////////////////////////
    //
    // compute using CGAL
    //
    //////////////////////////////////////////////////////
    //////////////////////////////////////////////////////

    tstart = time(0); 

    /*
    ifstream  infile3(infile_nodes);
    
    if(infile3.fail())
    {
      cout << " Could not open the input file" << endl;
      exit(1);
    }

    std::vector<Point> pointsCGAL;

    while(infile3 >> val[0] >> val[1] >> val[2] >> val[3])
    {
      Point a(val[1], val[2], val[3]);
      pointsCGAL.push_back(a);
    }

    ifstream  infile4(infile_trias);

    if(infile4.fail())
    {
      cout << " Could not open the input file" << endl;
      exit(1);
    }

    std::list<Triangle> triangles;

    if(ETYPE==3)
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

    ifstream  infile3(infile_nodes);
    
    if(infile3.fail())
    {
      cout << " Could not open the input file" << endl;
      exit(1);
    }

    // create a cgal incremental builder

    //vector<double>  vertices;
    vector<Point> vertices;
   
    // add the polyhedron vertices
    while(infile3 >> val[0] >> val[1] >> val[2] >> val[3])
    {
      //vertices.push_back(val[1]);
      //vertices.push_back(val[2]);
      //vertices.push_back(val[3]);

      vertices.push_back(Point(val[1], val[2], val[3]));
    }
   
    ifstream  infile4(infile_trias);

    if(infile4.fail())
    {
      cout << " Could not open the input file" << endl;
      exit(1);
    }
    
    vector<int>  faces;

    // add the polyhedron triangles
    if(ETYPE==3)
    {
      while(infile4 >> valInt[0] >> valInt[1] >> valInt[2] >> valInt[3] >> valInt[4])
      {
        faces.push_back(valInt[2]-1);
        faces.push_back(valInt[3]-1);
        faces.push_back(valInt[4]-1);
      }
    }
    else
    {
      while(infile4 >> valInt[0] >> valInt[1] >> valInt[2] >> valInt[3] >> valInt[4] >> valInt[5])
      {
        faces.push_back(valInt[2]-1);
        faces.push_back(valInt[3]-1);
        faces.push_back(valInt[4]-1);
        faces.push_back(valInt[5]-1);
      }
    }
    
    /*
    // build a polyhedron from the loaded arrays
    Polyhedron  polyhedronTemp;
    polyhedron_builder<HalfedgeDS>  poly_builder( vertices, faces, ETYPE );
    polyhedronTemp.delegate( poly_builder );


    // constructs AABB tree
    //Tree tree(polyhedronTemp.facets_begin(), polyhedronTemp.facets_end(), polyhedronTemp);

    Tree tree;
    tree.insert(polyhedronTemp.facets_begin(), polyhedronTemp.facets_end(), polyhedronTemp);

    //Point_inside inside_tester(tree);

    Point_inside  *inside_tester;
    inside_tester = new Point_inside(tree);
    */

    Point_inside  *inside_tester;

    
    // function #1 ... does not work
    //inside_tester = createPointer1(vertices, faces, ETYPE);


    // function #2 ... does not work
    Polyhedron  polyhedronTemp;
    polyhedron_builder<HalfedgeDS>  poly_builder( vertices, faces, ETYPE );
    polyhedronTemp.delegate( poly_builder );

    //inside_tester = createPointer2(polyhedronTemp);

    // function #3

    Tree tree;
    tree.insert(polyhedronTemp.facets_begin(), polyhedronTemp.facets_end(), polyhedronTemp);

    inside_tester = createPointer3(tree);



    int nints;
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

            //Point  ray_begin(ptTemp[0], ptTemp[1], ptTemp[2]);
            //Point  ray_end(10.0, ptTemp[1], ptTemp[2]);

            //FT sqd = tree.squared_distance(ray_begin);
            //cout << " sqd = " << sqd << endl;
            //ff[ll] = (sqd >= 0.0);

            //count number of intersections
            //Ray ray_query(ray_begin, ray_end);
            //Segment segment_query(ray_begin, ray_end);

            //ff[ll] = tree.number_of_intersected_primitives(Ray(ray_begin, ray_end)) % 2;
            //ff[ll] = tree.number_of_intersected_primitives(segment_query) % 2;

            //ff[ll] = tree.number_of_intersected_primitives(Ray(
            //              Point(ptTemp[0], ptTemp[1], ptTemp[2]), 
            //              Point(100.0, ptTemp[1], ptTemp[2]))) % 2;


            CGAL::Bounded_side res = (*inside_tester)(Point(ptTemp[0], ptTemp[1], ptTemp[2]));

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
    cout << "Time taken = "<< difftime(tend, tstart) <<" second(s)."<< endl;

  return 0;
}










