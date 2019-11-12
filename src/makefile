#!smake

CC    = gcc
COPT  = -std=c++11 -O3 -w -Wno-write-strings -Wno-uninitialized -Wno-deprecated  -Wno-format -fpermissive -frounding-math

OBJFILES := $(patsubst %.cpp,%.o,$(wildcard *.cpp))


VTKLIB58 = /usr/local/lib -lvtkWidgets -lvtkRendering -lvtkImaging -lvtkGraphics -lvtkFiltering \
       -lvtkGenericFiltering -lvtkHybrid -lvtkIO -lvtkftgl -lvtkVolumeRendering -lvtksys -lvtkInfovis\
       -lvtkDICOMParser -lvtkCommon -lvtkexoIIc -lGL -lGLU

VTKLIB62 = /usr/lib/x86_64-linux-gnu  -lvtkCommonCore-6.2  -lvtkIOCore-6.2  -lvtkIOGeometry-6.2  -lvtkIOXML-6.2 -lvtkIOImage-6.2 -lvtkIOExport-6.2 \
         -lvtkRenderingCore-6.2   -lvtkFiltersCore-6.2  -lvtkFiltersGeneral-6.2    -lvtkFiltersModeling-6.2  -lvtkFiltersSources-6.2  \
         -lvtkRenderingAnnotation-6.2  -lvtkRenderingVolume-6.2  -lvtkRenderingVolumeOpenGL-6.2  -lvtkRenderingFreeTypeOpenGL-6.2 \
         -lvtkInteractionStyle-6.2  -lvtkCommonExecutionModel-6.2 -lvtkCommonDataModel-6.2 \
		 -lvtkFiltersExtraction-6.2


CINCL = -I /usr/include/c++/5.4.0 \
        -I /usr/include/vtk-6.2 \
        -I /usr/include/boost \

#        -I /usr/include/CGAL \
#        -I /home/chenna/Documents/IGACode/mpap3/eigen-devel \
#        -I /home/chenna/Documents/IGACode/mpap3/eigen-devel/unsupported \
#        -I /home/chenna/Documents/IGACode/mpap3/src/HBsplines \


LIB   = -L/usr/lib/gcc/x86_64-linux-gnu/5.4.0 -lstdc++  -lgmp -lmpfr \
        -L$(VTKLIB62) \
		-L/usr/lib/x86_64-linux-gnu -lboost_thread -lboost_system \
    	-L/usr/lib/x86_64-linux-gnu -lCGAL -lCGAL_Core -lCGAL_ImageIO \



##################################################################################


vtk1: $(OBJFILES)
	$(CC) $(OBJFILES)  $(LIB) -o $@

%.o: %.cpp
	$(CC) -c $(COPT)  $(CINCL) $< -o $@

ex1: ex1.o ${OBJECTS}
	${CXX} -o ex1 ex1.o ${OBJECTS} ${LIB}
	rm  ex1.o

ex2: ex2.o ${OBJECTS}
	${CXX} -o ex2 ex2.o ${OBJECTS} ${LIB}
	rm  ex2.o

cutfem3d: cutfem3d.o ${OBJECTS}
	${CXX} -o cutfem3d cutfem3d.o ${OBJECTS} ${LIB}
	rm  cutfem3d.o

cutfem3dv1: cutfem3dv1.o ${OBJECTS}
	${CXX} -o cutfem3dv1 cutfem3dv1.o ${OBJECTS} ${LIB}
	rm  cutfem3dv1.o

cutfem3dv2: cutfem3dv2.o ${OBJECTS}
	${CXX} -o cutfem3dv2 cutfem3dv2.o ${OBJECTS} ${LIB}
	rm  cutfem3dv2.o


polygonmesh: polygonmesh.o ${OBJECTS}
	${CXX} -o polygonmesh polygonmesh.o ${OBJECTS} ${LIB}
	rm  polygonmesh.o


.PHONEY: clean
clean:
	rm *.o
	rm vtk1.*


