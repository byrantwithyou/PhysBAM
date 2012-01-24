//#####################################################################
// Copyright 2012, Alexey Stomakhin.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Data_Structures/HASHTABLE.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/ARRAYS_ND.h>
#include <PhysBAM_Tools/Read_Write/Grids_Uniform/READ_WRITE_GRID.h>
#include <PhysBAM_Tools/Read_Write/Utilities/FILE_UTILITIES.h>
#include <PhysBAM_Geometry/Geometry_Particles/REGISTER_GEOMETRY_READ_WRITE.h>
#include <PhysBAM_Geometry/Implicit_Objects_Uniform/LEVELSET_IMPLICIT_OBJECT.h>
#include <PhysBAM_Geometry/Read_Write/Geometry/READ_WRITE_BOX.h>
#include <PhysBAM_Geometry/Read_Write/Geometry/READ_WRITE_HEXAHEDRALIZED_VOLUME.h>
#include <PhysBAM_Geometry/Read_Write/Geometry/READ_WRITE_SEGMENTED_CURVE.h>
#include <PhysBAM_Geometry/Read_Write/Geometry/READ_WRITE_SEGMENTED_CURVE_2D.h>
#include <PhysBAM_Geometry/Read_Write/Geometry/READ_WRITE_TETRAHEDRALIZED_VOLUME.h>
#include <PhysBAM_Geometry/Read_Write/Geometry/READ_WRITE_TRIANGULATED_AREA.h>
#include <PhysBAM_Geometry/Read_Write/Geometry/READ_WRITE_TRIANGULATED_SURFACE.h>
#include <PhysBAM_Geometry/Read_Write/Implicit_Objects_Uniform/READ_WRITE_LEVELSET_IMPLICIT_OBJECT.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/HEXAHEDRALIZED_VOLUME.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TETRAHEDRALIZED_VOLUME.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TRIANGULATED_AREA.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TRIANGULATED_SURFACE.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Particles/PARTICLES.h>
#include <PhysBAM_Dynamics/Geometry/GENERAL_GEOMETRY_FORWARD.h>
#include <PhysBAM_Dynamics/Particles/PARTICLES_FORWARD.h>
#include <cstring>
#include <fstream>
using namespace PhysBAM;
using namespace FILE_UTILITIES;
template<class T> void Add_File(const std::string& filename,int number);
template<class T> void Add_Tri2D_File(const std::string& filename,int number);
template<class T> void Add_Tri_File(const std::string& filename,int number);
template<class T> void Add_Phi_File(const std::string& filename,int number);
template<class T> void Add_Curve_File(const std::string& filename,int number);
template<class T> void Add_Curve2D_File(const std::string& filename,int number);
template<class T> void Add_Tet_File(const std::string& filename,int number);
//#################################################################
static bool triangulated_surface_highlight_boundary=false;
static bool print_statistics=true;
//#################################################################
// Function main
//#################################################################
int main(int argc,char *argv[])
{
    bool type_double=false;    
    Initialize_Geometry_Particle();
    Initialize_Read_Write_General_Structures();

    ARRAY<std::string> files;
    for(int i=1;i<argc;i++){
        if(!strcmp(argv[i],"-float")) type_double=false;
        else if (!strcmp(argv[i],"-double")) type_double=true;
        else if (!strcmp(argv[i],"-show_boundary")) triangulated_surface_highlight_boundary = true;
        else if (!strcmp(argv[i],"-nostat")) print_statistics = false;
        else files.Append(argv[i]);}

    for(int i=0;i<files.m;i++){
        if(!type_double) Add_File<float>(files(i),i);
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
        else Add_File<double>(files(i),i);
#else
        else{LOG::cerr<<"Double support not enabled."<<std::endl;exit(1);}
#endif
        }
    return 0;
}
//#################################################################
// Function Add_File
//#################################################################
template<class T> void Add_File(const std::string& filename,int number)
{
    FILE_TYPE type=Get_File_Type(filename);
    switch(type){
        case TRI_FILE: Add_Tri_File<T>(filename,number);break;
        case TRI2D_FILE: Add_Tri2D_File<T>(filename,number);break;
        case PHI_FILE: Add_Phi_File<T>(filename,number);break;
        case CURVE_FILE: Add_Curve_File<T>(filename,number);break;
        case CURVE2D_FILE: Add_Curve2D_File<T>(filename,number);break;
        case TET_FILE: Add_Tet_File<T>(filename,number);break;
        default: LOG::cerr<<"Unrecognized file "<<filename<<std::endl;}
    LOG::cout<<std::flush;
}
//#################################################################
// Function Add_Tri2D_File
//#################################################################
template<class T> void Add_Tri2D_File(const std::string& filename,int number)
{
    try{
        TRIANGULATED_AREA<T>* area;
        FILE_UTILITIES::Create_From_File<T>(filename,area);}
    catch(FILESYSTEM_ERROR&){}
}
//#################################################################
// Function Add_Tri_File
//#################################################################
template<class T> void Add_Tri_File(const std::string& filename,int number)
{
    try{
        TRIANGULATED_SURFACE<T>* surface;
        FILE_UTILITIES::Create_From_File<T>(filename,surface);}
    catch(FILESYSTEM_ERROR&){}
}
//#################################################################
// Function Add_Phi_File
//#################################################################
template<class T> void Add_Phi_File(const std::string& filename,int number)
{
    try{
        LEVELSET_IMPLICIT_OBJECT<VECTOR<T,3> >* surface;
        FILE_UTILITIES::Create_From_File<T>(filename,surface);}
    catch(FILESYSTEM_ERROR&){}
}
//#################################################################
// Function Add_Curve_File
//#################################################################
template<class T> void Add_Curve_File(const std::string& filename,int number)
{
    try{
        SEGMENTED_CURVE<VECTOR<T,3> >* curve;
        FILE_UTILITIES::Create_From_File<T>(filename,curve);}
    catch(FILESYSTEM_ERROR&){}
}
//#################################################################
// Function Add_Curve2D_File
//#################################################################
template<class T> void Add_Curve2D_File(const std::string& filename,int number)
{
    try{
        SEGMENTED_CURVE_2D<T>* curve;
        FILE_UTILITIES::Create_From_File<T>(filename,curve);}
    catch(FILESYSTEM_ERROR&){}
}
//#################################################################
// Function Add_Tet_File
//#################################################################
template<class T> void Add_Tet_File(const std::string& filename,int number)
{
    try{
        TETRAHEDRALIZED_VOLUME<T>* tetrahedralized_volume;
        FILE_UTILITIES::Create_From_File<T>(filename,tetrahedralized_volume);}
    catch(FILESYSTEM_ERROR&){}
}
//#################################################################
