//#####################################################################
// Copyright 2010, Jon Gretarsson.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Namespace CONVERSION_TOOLS
//#####################################################################
#ifndef __READ_AEROF_FLUID__
#define __READ_AEROF_FLUID__

#include <PhysBAM_Tools/Vectors/VECTOR_3D.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TETRAHEDRALIZED_VOLUME.h>
#include <cstdio>

namespace PhysBAM{
namespace CONVERSION_TOOLS{

template<class T> void
Read_AeroF_Fluid_Mesh(const std::string& fluid_mesh_filename,TETRAHEDRALIZED_VOLUME<T>& ale_fluid_mesh)
{
    char* read_result;
    char input_buffer[1024];

    LOG::Time("Counting Elements");

    FILE* file=fopen(fluid_mesh_filename.c_str(),"r");
    int num_nodes=0,num_mesh_elements=0;

    read_result=fgets(input_buffer,1024,file);
    while(read_result != NULL){
        if(strncmp(input_buffer,"Nodes FluidNodes",16)==0){
            int particle_index;float x,y,z;
            read_result=fgets(input_buffer,1024,file);
            while(read_result != NULL && sscanf(input_buffer,"%d %f %f %f",&particle_index,&x,&y,&z)==4){
                num_nodes=max(num_nodes,particle_index);
                read_result=fgets(input_buffer,1024,file);}}
        else if(strncmp(input_buffer,"Elements FluidMesh_0 using FluidNodes",37)==0){
            int mesh_element_index,count,node1,node2,node3,node4;
            read_result=fgets(input_buffer,1024,file);
            while(read_result != NULL && sscanf(input_buffer,"%d %d %d %d %d %d",&mesh_element_index,&count,&node1,&node2,&node3,&node4)==6){
                num_mesh_elements=max(num_mesh_elements,mesh_element_index);
                read_result=fgets(input_buffer,1024,file);}}
        else read_result=fgets(input_buffer,1024,file);}
    fclose(file);
    LOG::cout<<"Reading in "<<num_nodes<<" nodes, and "<<num_mesh_elements<<" tets"<<std::endl;
    ale_fluid_mesh.particles.array_collection->Resize(num_nodes);
    ale_fluid_mesh.mesh.elements.Exact_Resize(num_mesh_elements,false);

    LOG::Time("Reading in Mesh Elements");
    file=fopen(fluid_mesh_filename.c_str(),"r");
    read_result=fgets(input_buffer,1024,file);
    while(read_result != NULL){
        if(strncmp(input_buffer,"Nodes FluidNodes",16)==0){
            int particle_index;float x,y,z;
            read_result=fgets(input_buffer,1024,file);
            while(read_result != NULL && sscanf(input_buffer,"%d %f %f %f",&particle_index,&x,&y,&z)==4){
                ale_fluid_mesh.particles.X(particle_index)=VECTOR<T,3>(x,y,z);
                read_result=fgets(input_buffer,1024,file);}}
        else if(strncmp(input_buffer,"Elements FluidMesh_0 using FluidNodes",37)==0){
            int mesh_element_index,count,node1,node2,node3,node4;
            read_result=fgets(input_buffer,1024,file);
            while(read_result != NULL && sscanf(input_buffer,"%d %d %d %d %d %d",&mesh_element_index,&count,&node1,&node2,&node3,&node4)==6){
                ale_fluid_mesh.mesh.elements(mesh_element_index)=VECTOR<int,4>(node1,node2,node3,node4);
                read_result=fgets(input_buffer,1024,file);}}
        else read_result=fgets(input_buffer,1024,file);}
    fclose(file);

    LOG::Stop_Time();
    LOG::cout<<std::endl;
}

template<class T> void
Read_AeroF_Fluid_Data(const std::string& prefix,COMPRESSIBLE_FLUID_PARTICLES<VECTOR<T,3> >& ale_fluid_data)
{
    return;
}
}
}
#endif
