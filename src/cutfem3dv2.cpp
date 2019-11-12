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
*/

#include "headersVTK.h"
#include "headersBasic.h"
#include "utilfuns.h"

#include <boost/foreach.hpp>

#include "findCutcellsVTK.h"
//#include "findCutcellsCGAL.h"
#include "findCutcellsCGAL2.h"


int main(int argc, char* argv[])
{
    //char infile_nodes[]="sphere-nodes.dat";
    //char infile_trias[]="sphere-trias.dat";

    //char infile_nodes[]="reliefvalve3Dblock-nodes.dat";
    //char infile_trias[]="reliefvalve3Dblock-trias.dat";

    char fileNamePoints[]="turekbeam3d-nodes.dat";
    char fileNameFaces[]="turekbeam3d-quads.dat";

    // char infile_nodes[]="thickplate-nodes.dat";
    // char infile_trias[]="thickplate-elems.dat";

    // char fileNamePoints[]="cylinder3d-nodes.dat";
    // char fileNameFaces[]="cylinder3d-quads.dat";


    //////////////////////////////////////////
    //////////////////////////////////////////
    //////////////////////////////////////////

    myGridData  gridData;

    gridData.elType=4;


    gridData.x0 = -1.6;    gridData.x1 =  1.6;    gridData.nEx = 5;
    gridData.y0 = -1.6;    gridData.y1 =  1.6;    gridData.nEy = 5;
    gridData.z0 = -1.6;    gridData.z1 =  1.6;    gridData.nEz = 5;


    if(argc == 0)
      cerr << " Input data " << endl;
    else
    {
       gridData.x0  = atof(argv[1]);
       gridData.x1  = atof(argv[2]);
       gridData.nEx = atoi(argv[3]);

       gridData.y0  = atof(argv[4]);
       gridData.y1  = atof(argv[5]);
       gridData.nEy = atoi(argv[6]);

       gridData.z0  = atof(argv[7]);
       gridData.z1  = atof(argv[8]);
       gridData.nEz = atoi(argv[9]);
    }

    gridData.nElem = gridData.nEx*gridData.nEy*gridData.nEz;

    gridData.nNx = gridData.nEx+1;
    gridData.nNy = gridData.nEy+1;
    gridData.nNz = gridData.nEz+1;

    gridData.nNode = gridData.nNx*gridData.nNy*gridData.nNz;

    gridData.printInfo();

    /////////////////////////////////////////
    //
    // VTK algorithm

    findCutcellsVTK(gridData, fileNamePoints, fileNameFaces);

    /////////////////////////////////////////
    //
    // CGAL algorithm

    //findCutcellsCGAL(gridData, fileNamePoints, fileNameFaces);

    /////////////////////////////////////////
    //
    // CGAL algorithm

    findCutcellsCGAL2(gridData, fileNamePoints, fileNameFaces);

    return 0;
}


