//#####################################################################
// Copyright 2005, Melody Wu.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Parsing/PARSE_ARGS.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TRIANGULATED_SURFACE.h>
#include <cstring>
#include <fstream>
using namespace PhysBAM;
using namespace std;

//#################################################################
// Function PatchTogether
//#################################################################
template<class T> void PatchTogether(std::string old_filename,std::string old_color_filename,std::string new_filename,std::string new_color_filename,std::string condensation_mapping_filename)
{
    ARRAY<int> condensation_mapping;
    ARRAY<VECTOR_3D<T> > old_vertex_color_array, new_vertex_color_array;
    TRIANGULATED_SURFACE<T>* old_triangulated_surface=0;
    TRIANGULATED_SURFACE<T>* new_triangulated_surface=0;
    
    FILE_UTILITIES::Create_From_File<T>(old_filename,old_triangulated_surface);
    FILE_UTILITIES::Create_From_File<T>(new_filename,new_triangulated_surface);
    FILE_UTILITIES::Read_From_File<T>(old_color_filename,old_vertex_color_array);
    FILE_UTILITIES::Read_From_File<T>(new_color_filename,new_vertex_color_array);
    FILE_UTILITIES::Read_From_File<T>(condensation_mapping_filename,condensation_mapping);
    //std::cout<<"# of vertex colors = "<<new_vertex_color_array.m<<std::endl;

    int original_size = old_triangulated_surface->particles.number;
    //int original_size = condensation_mapping.m;
    std::cout<<"number of particles in old = "<<original_size<<std::endl;
    int mapping_count = 0, number_new_vertices_needed = 0, number_vertices_in_new;
    number_vertices_in_new = new_triangulated_surface->particles.number;
    std::cout<<"number of new particles = "<<number_vertices_in_new<<std::endl;
    /*std::cout<<"old --- new "<<std::endl;
    for(int i=1; i<=old_size; i++){
        std::cout<<"i = "<<i<<" ";
        std::cout<<old_triangulated_surface->particles.X(i)<<" --- "<<new_triangulated_surface->particles.X(i)<<std::endl;
    }*/
    for(int i=1; i<=condensation_mapping.m; i++){
        if(condensation_mapping(i)>mapping_count){
            mapping_count = condensation_mapping(i);}}
    std::cout<<std::endl;
    std::cout<<"mapping count = "<<mapping_count<<std::endl;
    number_new_vertices_needed = number_vertices_in_new - mapping_count;
    std::cout<<"number of new vertices needed = "<<number_new_vertices_needed<<std::endl;
    
    if(number_new_vertices_needed!=0){
        condensation_mapping.Resize(original_size+number_new_vertices_needed);
        old_vertex_color_array.Resize(original_size+number_new_vertices_needed);
        old_triangulated_surface->particles.Add_Particles(number_new_vertices_needed);
        old_triangulated_surface->triangle_mesh.number_nodes = old_triangulated_surface->particles.number;
        for(int j=1; j<=number_new_vertices_needed; j++){
            //std::cout<<" j = "<<j<<std::endl;
            condensation_mapping(original_size+j) = mapping_count+j;
            old_vertex_color_array(original_size+j) = new_vertex_color_array(mapping_count+j);}}
    for(int i=1; i<=condensation_mapping.m; i++){
        int new_vertex_index = condensation_mapping(i);
        if(new_vertex_index!=0){
            old_triangulated_surface->particles.X(i) = new_triangulated_surface->particles.X(new_vertex_index);}}
    
    int triangle_count = old_triangulated_surface->triangle_mesh.triangles.m;
    std::cout<<"number of triangles = " <<triangle_count<<std::endl;
    int number_of_new_triangles = 0;
    for(int i=1; i<=new_triangulated_surface->triangle_mesh.triangles.m; i++){
        int v1,v2,v3;
        new_triangulated_surface->triangle_mesh.triangles.Get(i,v1,v2,v3);
        if((v1>mapping_count)||(v2>mapping_count)||(v3>mapping_count)){number_of_new_triangles++;}}
    std::cout<<"number of new triangles = "<<number_of_new_triangles<<std::endl;
    old_triangulated_surface->triangle_mesh.triangles.Exact_Resize(3,triangle_count+number_of_new_triangles);
    triangle_count++;
    
    for(int i=1; i<=new_triangulated_surface->triangle_mesh.triangles.m; i++){
        int v1,v2,v3;
        new_triangulated_surface->triangle_mesh.triangles.Get(i,v1,v2,v3);
        if((v1>mapping_count)||(v2>mapping_count)||(v3>mapping_count)){
            int p1,p2,p3;
            bool temp = condensation_mapping.Find(v1,p1);
            temp = condensation_mapping.Find(v2,p2);
            temp = condensation_mapping.Find(v3,p3);
            std::cout<<p1<<" "<<p2<<" "<<p3<<std::endl;
            old_triangulated_surface->triangle_mesh.triangles.Set(triangle_count,p1,p2,p3);
            triangle_count++;}}

    FILE_UTILITIES::Write_To_File<T>("patched.tri",*old_triangulated_surface);
    FILE_UTILITIES::Write_To_File<T>("patched.col",old_vertex_color_array);
    
}
//#################################################################
// Function main
//#################################################################
int main(int argc,char *argv[])
{
    PARSE_ARGS args;
    args.Set_Extra_Arguments(5,"<old filename> <old color filename> <new filename> <new color filename> <condensation mapping file>");
    args.Parse(argc,argv);
    std::cout<<"num extra "<<args.Num_Extra_Args()<<std::endl;
    std::string old_filename = args.Extra_Arg(1);
    std::string old_color_filename = args.Extra_Arg(2);
    std::string new_filename = args.Extra_Arg(3);
    std::string new_color_filename = args.Extra_Arg(4);
    std::string condensation_mapping_filename = args.Extra_Arg(5);
    PatchTogether<float>(old_filename,old_color_filename,new_filename,new_color_filename,condensation_mapping_filename);
    return 0;
}
//#################################################################