cmake_minimum_required (VERSION 2.8)

project (cgal)

enable_language(CXX)

##################################
# project directories
##################################


#PROJECT_BINARY_DIR
#PROJECT_SOURCE_DIR

##################################
# compiler options
##################################

# C++ compiler

set(CMAKE_CXX_COMPILER "/usr/bin/gcc-5")

if(CXX_COMPILER_NAME MATCHES "icc*")
    message ("Using INTEL C++ compiler")
    set(COPT "-O3 -Wall -DSELF_CHECK -DUNIX")
endif()

if(CXX_COMPILER_NAME MATCHES "gcc*")
    message ("Using GNU C++ compiler")
    set(COPT "-std=c++11 -O3 -w -Wno-write-strings -Wno-uninitialized -Wno-deprecated  -Wno-format -fpermissive -frounding-math")
endif()

set(COPT "-std=c++11 -O3 -w -Wno-write-strings -Wno-uninitialized -Wno-deprecated  -Wno-format -fpermissive -frounding-math")


set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${COPT}")


##################################
# include directories
##################################
include_directories(/usr/include/c++/5.4.0  /usr/include/vtk-6.2  /usr/include/boost)


##################################
# lib directories
##################################
link_directories(/usr/lib/gcc/x86_64-linux-gnu/5.4.0  /usr/lib/x86_64-linux-gnu)



##################################
# source files
##################################


file(GLOB sources src/*h)

##################################
##################################

add_executable(cutfem3dv2  src/cutfem3dv2.cpp ${sources})

set(VTKLIB62 "vtkCommonCore-6.2 vtkIOCore-6.2  vtkIOGeometry-6.2  vtkIOXML-6.2  vtkIOImage-6.2  vtkIOExport-6.2   vtkRenderingCore-6.2  vtkFiltersCore-6.2  vtkFiltersGeneral-6.2   vtkFiltersModeling-6.2  vtkFiltersSources-6.2  vtkRenderingAnnotation-6.2  vtkRenderingVolume-6.2  vtkRenderingVolumeOpenGL-6.2  vtkRenderingFreeTypeOpenGL-6.2  vtkInteractionStyle-6.2  vtkCommonExecutionModel-6.2  vtkCommonDataModel-6.2  vtkFiltersExtraction-6.2")


target_link_libraries(cutfem3dv2 stdc++  gmp  mpfr  vtkCommonCore-6.2  vtkIOCore-6.2  vtkIOGeometry-6.2  vtkIOXML-6.2  vtkIOImage-6.2  vtkIOExport-6.2   vtkRenderingCore-6.2  vtkFiltersCore-6.2  vtkFiltersGeneral-6.2   vtkFiltersModeling-6.2  vtkFiltersSources-6.2  vtkRenderingAnnotation-6.2  vtkRenderingVolume-6.2  vtkRenderingVolumeOpenGL-6.2  vtkRenderingFreeTypeOpenGL-6.2  vtkInteractionStyle-6.2  vtkCommonExecutionModel-6.2  vtkCommonDataModel-6.2  vtkFiltersExtraction-6.2  boost_thread boost_system  CGAL  CGAL_Core  CGAL_ImageIO m)

#target_link_libraries(cutfem3dv2 stdc++  gmp  mpfr  "${VTKLIB62}"  boost_thread boost_system  CGAL  CGAL_Core  CGAL_ImageIO m)

