//#####################################################################
// Copyright 2002, 2003, Robert Bridson, Ronald Fedkiw, Eran Guendelman, Neil Molino.
// This file is part of PhysBAM whose distribution is governed by the license 
// contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class HANGING_CURTAIN_TEST
//##################################################################### 
#ifndef __HANGING_CURTAIN_TEST__
#define __HANGING_CURTAIN_TEST__

#include "../CLOTH_EXAMPLE.h"
namespace PhysBAM{

template<class T>
class HANGING_CURTAIN_TEST:public CLOTH_EXAMPLE<T>
{
public:
    int number_side_panels;
    T aspect_ratio,side_length;
    T cloth_density;
    bool cloth_mass_area_weighted;

    HANGING_CURTAIN_TEST()
        :number_side_panels(100),aspect_ratio((T)1.5),side_length(10),cloth_density(0.04),cloth_mass_area_weighted(true)
    {
        T edge_stiffness_scaling=80;
        T altitude_stiffness_scaling=80;

        final_time=200;
        restart_step_number=0;
        cfl_number=(T)5.9;
        cg_tolerance=(T)1e-3;
        use_masses_and_springs=true;use_altitude_springs=true;
        edge_spring_stiffness=(T)(edge_stiffness_scaling/(1+root_two));edge_spring_overdamping_fraction=2;
        use_shortest_altitude_spring_only=true;
        use_altitude_springs_compressed_by_threshold_only=true;altitude_spring_fraction_compression=(T).1;
        altitude_spring_stiffness=altitude_stiffness_scaling*4/(1+(T)root_two);altitude_spring_overdamping_fraction=4;
        bending_stiffness=(T).001;bending_damping=(T).001;
        output_directory="Hanging_Curtain_Test/output";
        check_initial_mesh_for_self_intersection=false;
        perform_self_collision=true;
    }

    ~HANGING_CURTAIN_TEST()
    {}

//#####################################################################
// Function Get_Initial_Data
//#####################################################################
void Get_Initial_Data(TRIANGULATED_SURFACE<T>& triangulated_surface)
{
    TRIANGLE_MESH& triangle_mesh=triangulated_surface.triangle_mesh;
    PARTICLES<T,VECTOR_3D<T> >& particles=triangulated_surface.particles;
    
    particles.Update_Position_And_Velocity();particles.Store_Mass();

    int m=(int)(aspect_ratio*number_side_panels)+1,n=number_side_panels+1;
    triangle_mesh.Initialize_Herring_Bone_Mesh(m,n);
    for(int k=0;k<triangle_mesh.number_nodes;k++) particles.array_collection->Add_Element();
    T mass_node=aspect_ratio*sqr(side_length)/(m*n);ARRAY<T,VECTOR<int,1> >::copy(mass_node,particles.mass); 
    T dx=aspect_ratio*side_length/(m-1),dy=side_length/(n-1);

    VECTOR_3D<T> corner(10,19,5);
    VECTOR_3D<T> cloth_x_direction(0,-1,0);
    VECTOR_3D<T> cloth_y_direction(0,0,1);

    for(int i=0;i<m;i++) for(int j=0;j<n;j++){int node=i+m*(j-1);
        particles.X(node)=corner + (i-1)*dx*cloth_x_direction + (j-1)*dy*cloth_y_direction;
        particles.V(node)=VECTOR_3D<T>(0,0,0);}

    // set the mass
    triangulated_surface.Set_Density(cloth_density);
    triangulated_surface.Set_Mass_Of_Particles(!cloth_mass_area_weighted);
}
//#####################################################################
// Function Set_External_Velocities
//#####################################################################
// for external forces and velocities
void Set_External_Velocities(ARRAY<VECTOR_3D<T> ,VECTOR<int,1> >& V,const T time)
{
    int i,j;
    int m=(int)(aspect_ratio*number_side_panels)+1,n=number_side_panels+1;
    i=1;j=1;V(i+m*(j-1))=VECTOR_3D<T>(0,0,0);
    i=1;j=n;V(i+m*(j-1))=VECTOR_3D<T>(0,0,0);
}
//#####################################################################
// Function Zero_Out_Enslaved_Velocity_Nodes
//#####################################################################
// for external forces and velocities
void Zero_Out_Enslaved_Velocity_Nodes(ARRAY<VECTOR_3D<T> ,VECTOR<int,1> >& V,const T time)
{
    int i,j;
    int m=(int)(aspect_ratio*number_side_panels)+1,n=number_side_panels+1;
    i=1;j=1;V(i+m*(j-1))=VECTOR_3D<T>(0,0,0);
    i=1;j=n;V(i+m*(j-1))=VECTOR_3D<T>(0,0,0);
}
//#####################################################################
};
}
#endif

