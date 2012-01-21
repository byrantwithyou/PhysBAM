//#####################################################################
// Copyright 2009, Avi Robinson-Mosher.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Parsing/PARSE_ARGS.h>
#include <PhysBAM_Tools/Read_Write/Utilities/FILE_UTILITIES.h>
#include <PhysBAM_Geometry/Read_Write/Geometry/READ_WRITE_TRIANGULATED_SURFACE.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TRIANGLE_SUBDIVISION.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TRIANGULATED_SURFACE.h>
#include <PhysBAM_Geometry/Basic_Geometry/BOX.h>

using namespace PhysBAM;
template<class T>
void Fractal_Midpoint_Subdivide(TRIANGULATED_SURFACE<T>& surface,const int steps_of_subdivision=1,const T fractal_dimension=(T).5,const T scale_factor=(T).5)
{
    T ratio=pow((T)2,-fractal_dimension);
    surface.Update_Bounding_Box();
    T scale=scale_factor*surface.bounding_box->Edge_Lengths().Max();
    T std=scale*ratio;
    
    for(int i=0;i<steps_of_subdivision;i++){
        std::cout<<"Subdivision step "<<i<<std::endl;
        TRIANGLE_SUBDIVISION subdivision(surface.mesh);
        TRIANGLE_MESH refined_mesh;
        subdivision.Refine_Mesh(refined_mesh);
        int number_new_particles=refined_mesh.number_nodes-surface.particles.array_collection->Size();
        ARRAY<VECTOR<T,3> > X_save(surface.particles.X);
        surface.particles.array_collection->Add_Elements(number_new_particles);
        subdivision.Apply_Fractal_Subdivision(X_save,surface.particles.X,std);
        surface.mesh.Initialize_Mesh(refined_mesh);
        std*=ratio;
    }
    surface.Loop_Subdivide(); // one final subdivision to smooth it a bit
}

int main(int argc,char *argv[])
{
    typedef float T;
    PARSE_ARGS parse_args;
    parse_args.Add_String_Argument("-i",""); // input surface
    parse_args.Add_String_Argument("-o",""); // output surface
    parse_args.Add_Integer_Argument("-depth",1);
    parse_args.Add_Double_Argument("-power",(T).25);
    parse_args.Add_Double_Argument("-scale",(T).25);
    parse_args.Parse(argc,argv);
    TRIANGULATED_SURFACE<T>& triangulated_surface=*TRIANGULATED_SURFACE<T>::Create();
    FILE_UTILITIES::Read_From_File<T>(parse_args.Get_String_Value("-i"),triangulated_surface);
    Fractal_Midpoint_Subdivide(triangulated_surface,parse_args.Get_Integer_Value("-depth"),(T)parse_args.Get_Double_Value("-power"),(T)parse_args.Get_Double_Value("-scale"));
    FILE_UTILITIES::Write_To_File<T>(parse_args.Get_String_Value("-o"),triangulated_surface);
    return 0;
}
