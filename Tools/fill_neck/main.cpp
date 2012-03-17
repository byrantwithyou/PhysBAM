//#####################################################################
// Copyright 2005, Melody Wu, Jerry Talton.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Parsing/PARSE_ARGS.h>
#include <PhysBAM_Tools/Random_Numbers/RANDOM_NUMBERS.h>
#include <PhysBAM_Tools/Read_Write/FILE_UTILITIES.h>
#include <PhysBAM_Geometry/Topology/SEGMENT_MESH.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TRIANGULATED_SURFACE.h>
#include <cstring>
#include <fstream>
using namespace PhysBAM;
using namespace std;

//#################################################################
// Function FillNeck
//#################################################################
template<class T> void FillNeck(std::string tri_filename,std::string new_filename)
{
    TRIANGULATED_SURFACE<T>* triangulated_surface=0;
    FILE_UTILITIES::Create_From_File<T>(tri_filename,triangulated_surface);
    std::cout<<"done reading from file"<<std::endl;
    triangulated_surface->Print_Mesh_Diagnostics();
    triangulated_surface->triangle_mesh.Initialize_Boundary_Mesh();
    std::cout<<"1"<<std::endl;
    triangulated_surface->triangle_mesh.boundary_mesh->Initialize_Connected_Segments();
    std::cout<<"2"<<std::endl;
    SEGMENT_MESH* boundary_mesh = triangulated_surface->triangle_mesh.boundary_mesh;
    std::cout<<"3"<<std::endl;
    ARRAY<ARRAYS<int> >* connected_segments = boundary_mesh->connected_segments;
    std::cout<<"number of holes to fill = "<<connected_segments->m<<std::endl;
    ARRAYS<int> largest_segment_array;
    connected_segments->Get(1,largest_segment_array);
    std::cout<<"4"<<std::endl;
    int triangle_count = triangulated_surface->triangle_mesh.triangles.m;
    triangulated_surface->triangle_mesh.triangles.Exact_Resize(3,triangle_count+largest_segment_array.m);
    VECTOR_3D<T> center_point(0,-.2,0);
    triangulated_surface->particles.Increase_Array_Size(1);
    int center_index = triangulated_surface->particles.Add_Particle();
    triangulated_surface->particles.X(center_index) = center_point;
    for(int i=0;i<largest_segment_array.m;i++){
        int j,k;
        largest_segment_array.Get(i,j,k);
        triangle_count++;
        triangulated_surface->triangle_mesh.triangles.Set(triangle_count,j,center_index,k);
    }
    triangulated_surface->triangle_mesh.number_nodes = triangulated_surface->particles.number;

    std::cout<<"closing surface"<<std::endl;
    triangulated_surface->Close_Surface(false,0.0,true,true);
    std::cout<<"surface closed"<<std::endl;
    std::cout<<"writing to file"<<std::endl;
    FILE_UTILITIES::Write_To_File<double>(new_filename+".tri",*triangulated_surface);
}

//#################################################################
// Function main
//#################################################################
int main(int argc,char *argv[])
{
    PARSE_ARGS args;
    args.Set_Extra_Arguments(2,"<tri filename> <new filename>");
    args.Parse(argc,argv);
    std::cout<<"num extra "<<args.Num_Extra_Args()<<std::endl;
    std::string tri_filename = args.Extra_Arg(1);
    std::string new_filename = args.Extra_Arg(2);
    FillNeck<float>(tri_filename,new_filename);
    return 0;
}
//#################################################################
