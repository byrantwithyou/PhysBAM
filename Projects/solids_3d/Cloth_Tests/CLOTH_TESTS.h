//#####################################################################
// Copyright 2006, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __CLOTH_TESTS__
#define __CLOTH_TESTS__

#include <PhysBAM_Solids/PhysBAM_Deformables/Forces/MULTILINEAR_SPRINGS.h>

namespace PhysBAM{

template<class T,class RW>
class CLOTH_TESTS:public SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<TV>,RW>
{
public:
    typedef SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<TV>,RW> BASE;
    using BASE::last_frame;using BASE::restart;using BASE::restart_frame;using BASE::frame_rate;
    using BASE::output_directory;using BASE::data_directory;using BASE::fluids_parameters;
    using BASE::solids_parameters;

    int number_side_panels;
    T aspect_ratio,side_length;
    int example_number;
    int resolution;
    GRID<TV> cloth_grid;

    CLOTH_TESTS(int example_number_input,int resolution_input)
        :BASE(0,fluids_parameters.NONE),aspect_ratio((T)1.5),side_length(1),example_number(example_number_input),resolution(resolution_input)
    {
        std::cout<<"We have example = "<<example_number<<" and resolution="<<resolution<<std::endl;
        
        solids_parameters.cfl=(T).9;solids_parameters.implicit_solve_parameters.cg_tolerance=(T)1e-3;
        number_side_panels=resolution*15;
        output_directory=STRING_UTILITIES::string_sprintf("Cloth_Tests/%d_%d",example_number,resolution);

        last_frame=(int)((T)7*frame_rate);
    }

//#####################################################################
// Function Get_Initial_Data
//#####################################################################
void Get_Initial_Data()
{
    int index=solids_parameters.deformable_body_parameters.list.Add_Deformable_Triangulated_Surface();
    TRIANGULATED_SURFACE<T>& triangulated_surface=*solids_parameters.deformable_body_parameters.list(index).triangulated_surface;
    TRIANGLE_MESH& triangle_mesh=triangulated_surface.triangle_mesh;
    PARTICLES<T,VECTOR_3D<T> >& particles=triangulated_surface.particles;
    cloth_grid.Initialize(number_side_panels+1,(int)(aspect_ratio*number_side_panels)+1,0,1,0,aspect_ratio);
    triangle_mesh.Initialize_Herring_Bone_Mesh(cloth_grid.m,cloth_grid.n);
    particles.array_collection->Add_Elements(triangle_mesh.number_nodes);
    for(int i=1;i<=cloth_grid.m;i++) for(int j=1;j<=cloth_grid.n;j++){
        int node=i+cloth_grid.m*(j-1);
        particles.X(node)=VECTOR_3D<T>(cloth_grid.x(i),(T)1,cloth_grid.y(j));
        particles.V(node)=VECTOR_3D<T>();}
    T mass_node=aspect_ratio*sqr(side_length)/(cloth_grid.m*cloth_grid.n);ARRAY<T>::copy(mass_node,particles.mass.array);
}
//#####################################################################
// Function Initialize_Bodies
//#####################################################################
void Initialize_Bodies() PHYSBAM_OVERRIDE
{
    Get_Initial_Data();
    solids_parameters.deformable_body_parameters.list(1).Add_Body_Forces(*solids_parameters.deformable_body_parameters.list(1).triangulated_surface);
    // choose springs to use based on example
    bool use_multilinear_edge_springs=false,use_edge_springs=true,use_altitude_springs=true,use_bending_elements=true;
    T edge_spring_constant=2/(1+sqrt((T)2));;
    if(example_number==1){use_bending_elements=use_altitude_springs=false;}
    else if(example_number==2){use_bending_elements=use_altitude_springs=use_edge_springs=false;use_multilinear_edge_springs=true;}
    else if(example_number==3){use_bending_elements=use_altitude_springs=false;edge_spring_constant*=10;}
    // add springs
    if(use_multilinear_edge_springs){
        ARRAY<T> compression(2,1);ARRAY<T> extension(2,1);
        compression(1,1)=-.1;compression(2,1)=(T).1;extension(1,1)=.1;extension(2,1)=(T)10;
        solids_parameters.deformable_body_parameters.list(1).Add_Multilinear_Edge_Springs(solids_parameters.deformable_body_parameters.list(1).triangulated_surface->triangle_mesh,
                                                                                          compression,extension,edge_spring_constant,2);}
    if(use_edge_springs) solids_parameters.deformable_body_parameters.list(1).Add_Edge_Springs(solids_parameters.deformable_body_parameters.list(1).triangulated_surface->triangle_mesh,edge_spring_constant,2);
    if(use_altitude_springs) solids_parameters.deformable_body_parameters.list(1).Add_Altitude_Springs(solids_parameters.deformable_body_parameters.list(1).triangulated_surface->triangle_mesh,2*4/(1+sqrt((T)2)),4);
    if(use_bending_elements) solids_parameters.deformable_body_parameters.list(1).Add_Bending_Elements(solids_parameters.deformable_body_parameters.list(1).triangulated_surface->triangle_mesh);
}
//#####################################################################
// Function Set_External_Velocities
//#####################################################################
// for external forces and velocities
void Set_External_Velocities(ARRAY<VECTOR_3D<T> >& V,const T time)
{
    switch(id_number){
    case 1:
        {int i,j;i=1;j=1;V(i+cloth_grid.m*(j-1))=VECTOR_3D<T>(0,0,0);i=cloth_grid.m;j=1;V(i+cloth_grid.m*(j-1))=VECTOR_3D<T>(0,0,0);}
      break;
    default:std::cout<<"Unrecognized deformable object id number"<<std::endl;exit(1);}
}
//#####################################################################
// Function Zero_Out_Enslaved_Velocity_Nodes
//#####################################################################
// for external forces and velocities
void Zero_Out_Enslaved_Velocity_Nodes(ARRAY<VECTOR_3D<T> >& V,const T time)
{
    switch(id_number){
    case 1:
        {int i,j;i=1;j=1;V(i+cloth_grid.m*(j-1))=VECTOR_3D<T>(0,0,0);i=cloth_grid.m;j=1;V(i+cloth_grid.m*(j-1))=VECTOR_3D<T>(0,0,0);}
      break;
    default:std::cout<<"Unrecognized deformable object id number"<<std::endl;exit(1);}
}
//#####################################################################
// Function Write_Output_Files
//#####################################################################
void Write_Output_Files(const int frame) const PHYSBAM_OVERRIDE
{
    BASE::Write_Output_Files(frame);
    if(example_number==2){
        MULTILINEAR_SPRINGS<T,VECTOR_3D<T> >& springs=*solids_parameters.deformable_body_parameters.list(1).multilinear_springs(1);
        ARRAY<T> deformation(springs.youngs_modulus.m);
        if(springs.optimization_current_length.m==0) ARRAY<T>::copy((T)1,deformation);
        else for(int i=1;i<=deformation.m;i++) deformation(i)=(springs.optimization_current_length(i)-springs.visual_restlength(i))/springs.restlength(i);
        FILE_UTILITIES::Write_To_File<RW>(STRING_UTILITIES::string_sprintf("%s/deformation.%d",output_directory.c_str(),frame),deformation);
        FILE_UTILITIES::Write_To_File<RW>(STRING_UTILITIES::string_sprintf("%s/youngs_modulus.%d",output_directory.c_str(),frame),springs.youngs_modulus);}
}
//#####################################################################
};
}
#endif
