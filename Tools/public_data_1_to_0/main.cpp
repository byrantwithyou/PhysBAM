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
#include <PhysBAM_Solids/PhysBAM_Deformables/Particles/DEFORMABLE_PARTICLES.h>
#include <PhysBAM_Dynamics/Geometry/GENERAL_GEOMETRY_FORWARD.h>
#include <PhysBAM_Dynamics/Particles/PARTICLES_FORWARD.h>
#include <cstring>
#include <fstream>
using namespace PhysBAM;
using namespace FILE_UTILITIES;
template<class T> void Convert_File(const std::string& ifilename,const std::string& ofilename);
template<class T> void Convert_Tri_File(const std::string& ifilename,const std::string& ofilename);
template<class T> void Convert_Tri2D_File(const std::string& ifilename,const std::string& ofilename);
template<class T> void Convert_Phi_File(const std::string& ifilename,const std::string& ofilename);
template<class T> void Convert_Phi2D_File(const std::string& ifilename,const std::string& ofilename);
template<class T> void Convert_Tet_File(const std::string& ifilename,const std::string& ofilename);
template<class T> void Convert_Curve_File(const std::string& ifilename,const std::string& ofilename);
template<class T> void Convert_Curve2D_File(const std::string& ifilename,const std::string& ofilename);
//#################################################################
// Function main
//#################################################################
int main(int argc,char *argv[])
{
#ifndef COMPILE_WITH_READ_ONE_BASED_DATA
    LOG::cerr<<"1-indexed data reading functionality is not enabled."<<std::endl;exit(1);
#endif
    bool type_double=false;    
    Initialize_Geometry_Particle();
    Initialize_Read_Write_General_Structures();

    ARRAY<std::string> files;
    for(int i=1;i<argc;i++){
        if(!strcmp(argv[i],"-float")) type_double=false;
        else if (!strcmp(argv[i],"-double")) type_double=true;
        else files.Append(argv[i]);}

    if(files.m!=2){LOG::cerr<<"Input and output file names expected."<<std::endl;exit(1);}

    if(!type_double) Convert_File<float>(files(0),files(1));
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
    else Convert_File<double>(files(0),files(1));
#else
    else{LOG::cerr<<"Double support not enabled."<<std::endl;exit(1);}
#endif
    return 0;
}
//#################################################################
// Function Convert_File
//#################################################################
template<class T> void Convert_File(const std::string& ifilename,const std::string& ofilename)
{
    FILE_TYPE itype=Get_File_Type(ifilename);
    FILE_TYPE otype=Get_File_Type(ofilename);

    if(itype!=otype){LOG::cerr<<"Input and output file types do not match."<<std::endl;exit(1);}

    switch(itype){
        case TRI_FILE: Convert_Tri_File<T>(ifilename,ofilename);break;
        case TRI2D_FILE: Convert_Tri2D_File<T>(ifilename,ofilename);break;
        case PHI_FILE: Convert_Phi_File<T>(ifilename,ofilename);break;
        case PHI2D_FILE: Convert_Phi2D_File<T>(ifilename,ofilename);break;
        case TET_FILE: Convert_Tet_File<T>(ifilename,ofilename);break;
        case CURVE_FILE: Convert_Curve_File<T>(ifilename,ofilename);break;
        case CURVE2D_FILE: Convert_Curve2D_File<T>(ifilename,ofilename);break;
        default: LOG::cerr<<"Unrecognized file type "<<ifilename<<std::endl;}
    LOG::cout<<std::flush;
}
//#################################################################
// Function Convert_Tri_File
//#################################################################
template<class T> void Convert_Tri_File(const std::string& ifilename,const std::string& ofilename)
{
    LOG::cout<<"TRI: "<<ifilename<<" -> "<<ofilename<<std::endl;
    try{
        TRIANGULATED_SURFACE<T>* surface;
        FILE_UTILITIES::Create_From_File<T>(ifilename,surface);
        for(int i=0; i<surface->mesh.elements.m;i++){
            VECTOR<int,3>& element=surface->mesh.elements(i);
            element-=1;
            if (element.Min()<0){
                LOG::cerr<<"Negative vertex index"<<std::endl; PHYSBAM_FATAL_ERROR();}}
        FILE_UTILITIES::Write_To_File<T>(ofilename,*surface);
    }
    catch(FILESYSTEM_ERROR&){}
}
//#################################################################
// Function Convert_Tri2D_File
//#################################################################
template<class T> void Convert_Tri2D_File(const std::string& ifilename,const std::string& ofilename)
{
    LOG::cout<<"TRI2D: "<<ifilename<<" -> "<<ofilename<<std::endl;
    try{
        TRIANGULATED_AREA<T>* area;
        FILE_UTILITIES::Create_From_File<T>(ifilename,area);
        for(int i=0; i<area->mesh.elements.m;i++){
            VECTOR<int,3>& element=area->mesh.elements(i);
            element-=1;
            if (element.Min()<0){
                LOG::cerr<<"Negative vertex index"<<std::endl; PHYSBAM_FATAL_ERROR();}}
        FILE_UTILITIES::Write_To_File<T>(ofilename,*area);
    }
    catch(FILESYSTEM_ERROR&){}
}
//#################################################################
// Function Convert_Phi_File
//#################################################################
template<class T> void Convert_Phi_File(const std::string& ifilename,const std::string& ofilename)
{
    LOG::cout<<"PHI: "<<ifilename<<" -> "<<ofilename<<std::endl;
    try{
        LEVELSET_IMPLICIT_OBJECT<VECTOR<T,3> >* surface;
        FILE_UTILITIES::Create_From_File<T>(ifilename,surface);
        VECTOR<int,3>& min_corner = surface->levelset.phi.domain.min_corner;
        if (min_corner.Min()<0){
            LOG::cerr<<"Negative vertex index"<<std::endl; PHYSBAM_FATAL_ERROR();}
        FILE_UTILITIES::Write_To_File<T>(ofilename,*surface);
    }
    catch(FILESYSTEM_ERROR&){}
}
//#################################################################
// Function Convert_Phi2D_File
//#################################################################
template<class T> void Convert_Phi2D_File(const std::string& ifilename,const std::string& ofilename)
{
    LOG::cout<<"PHI2: "<<ifilename<<" -> "<<ofilename<<std::endl;
    try{
        LEVELSET_IMPLICIT_OBJECT<VECTOR<T,2> >* area;
        FILE_UTILITIES::Create_From_File<T>(ifilename,area);
        VECTOR<int,2>& min_corner = area->levelset.phi.domain.min_corner;
        if (min_corner.Min()<0){
            LOG::cerr<<"Negative vertex index"<<std::endl; PHYSBAM_FATAL_ERROR();}
        FILE_UTILITIES::Write_To_File<T>(ofilename,*area);
    }
    catch(FILESYSTEM_ERROR&){}
}
//#################################################################
// Function Convert_Tet_File
//#################################################################
template<class T> void Convert_Tet_File(const std::string& ifilename,const std::string& ofilename)
{
    LOG::cout<<"TET: "<<ifilename<<" -> "<<ofilename<<std::endl;
    try{
        TETRAHEDRALIZED_VOLUME<T>* tetrahedralized_volume;
        FILE_UTILITIES::Create_From_File<T>(ifilename,tetrahedralized_volume);
        for(int i=0; i<tetrahedralized_volume->mesh.elements.m;i++){
            VECTOR<int,4>& element=tetrahedralized_volume->mesh.elements(i);
            element-=1;
            if (element.Min()<0){
                LOG::cerr<<"Negative vertex index"<<std::endl; PHYSBAM_FATAL_ERROR();}}
        FILE_UTILITIES::Write_To_File<T>(ofilename,*tetrahedralized_volume);
    }
    catch(FILESYSTEM_ERROR&){}
}
//#################################################################
// Function Convert_Curve_File
//#################################################################
template<class T> void Convert_Curve_File(const std::string& ifilename,const std::string& ofilename)
{
    LOG::cout<<"CUR: "<<ifilename<<" -> "<<ofilename<<std::endl;
    try{
        SEGMENTED_CURVE<VECTOR<T,3> >* curve;
        FILE_UTILITIES::Create_From_File<T>(ifilename,curve);
        for(int i=0; i<curve->mesh.elements.m;i++){
            VECTOR<int,2>& element=curve->mesh.elements(i);
            element-=1;
            if (element.Min()<0){
                LOG::cerr<<"Negative vertex index"<<std::endl; PHYSBAM_FATAL_ERROR();}}
        FILE_UTILITIES::Write_To_File<T>(ofilename,*curve);
    }
    catch(FILESYSTEM_ERROR&){}
}
//#####################################################################
// Function Convert_Curve2_File
//#####################################################################
template<class T> void Convert_Curve2D_File(const std::string& ifilename,const std::string& ofilename)
{
    LOG::cout<<"CUR2: "<<ifilename<<" -> "<<ofilename<<std::endl;
    try{
        SEGMENTED_CURVE_2D<T>* curve;
        FILE_UTILITIES::Create_From_File<T>(ifilename,curve);
        for(int i=0; i<curve->mesh.elements.m;i++){
            VECTOR<int,2>& element=curve->mesh.elements(i);
            element-=1;
            if (element.Min()<0){
                LOG::cerr<<"Negative vertex index"<<std::endl; PHYSBAM_FATAL_ERROR();}}
        FILE_UTILITIES::Write_To_File<T>(ofilename,*curve);
    }
    catch(FILESYSTEM_ERROR&){}
}
//#################################################################
