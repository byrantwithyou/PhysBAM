//#####################################################################
// Copyright 2003, Neil Molino, Joseph Teran.
// This file is part of PhysBAM whose distribution is governed by the license
// contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class BALL_AND_WALL_EXAMPLE
//#####################################################################
#ifndef __BALL_AND_WALL_EXAMPLE__
#define __BALL_AND_WALL_EXAMPLE__

#include <PhysBAM_Tools/Math_Tools/constants.h>
#include <PhysBAM_Tools/Matrices/QUATERNION.h>
#include "../TET_SIM_EMBEDDED_EXAMPLE.h"
using namespace PhysBAM;

template <class T>
class BALL_AND_WALL_EXAMPLE:public TET_SIM_EMBEDDED_EXAMPLE<T>
{
public:
    T initial_height;
    QUATERNION<T> initial_orientation;
    VECTOR_3D<T> initial_velocity,initial_angular_velocity;
    int m_input,n_input,mn_input;
    T xmin_input,xmax_input,ymin_input,ymax_input,zmin_input,zmax_input; 
    GRID<TV> mattress_grid;
    ARRAY<int> constrained_nodes;
    bool cube_mesh;
    
    BALL_AND_WALL_EXAMPLE()
        :initial_height(0),initial_orientation(0,VECTOR_3D<T>(1,0,0)),initial_velocity(0,0,0),initial_angular_velocity(0,0,0),
        m_input(2),n_input(9),mn_input(17),xmin_input((T)-.1),xmax_input((T)-xmin_input),ymin_input((T)-.72),ymax_input((T)-ymin_input),zmin_input((T)-1.5),zmax_input((T)-zmin_input),
        mattress_grid(m_input,n_input,mn_input,xmin_input,xmax_input,ymin_input,ymax_input,zmin_input,zmax_input),
        cube_mesh(true)
    {
        restart_step_number=0;
        final_time=5;                                         
        frame_rate=24;
        cfl_number=(T).5;
        cg_iterations=250;
        youngs_modulus=(T)7.9e5;poissons_ratio=(T).4;Rayleigh_coefficient=(T).1;
        strcpy(output_directory,"Ball_And_Wall/output");
        
        perform_fracture=true;
        fracture_threshold_1=(T)1e4;fracture_threshold_2=(T)1e5;fracture_threshold_3=(T)1e5;
        compressive_threshold_1=(T)-fracture_threshold_1*1000;compressive_threshold_2=(T)-1000*fracture_threshold_1;compressive_threshold_3=(T)-1000*fracture_threshold_1;

        interpolation_fraction_threshhold=(T).1;
        max_number_of_cuts=4;
        number_of_fracture_initiation_points=2000;
        fracture_bias_propagation=0;
        perturb_amount_for_collision_freeness=(T)1e-5;
        perform_self_collision=false;
        push_out=false;
        floor_friction=0;
    }

    virtual ~BALL_AND_WALL_EXAMPLE()
    {}

//#####################################################################
// Function Initialize_Tetrahedralized_Volume
//#####################################################################
void Initialize_Tetrahedralized_Volume(TETRAHEDRALIZED_VOLUME<T>& tetrahedralized_volume)
{
    tetrahedralized_volume.particles.Update_Position();
    if(!cube_mesh) tetrahedralized_volume.Initialize_Octahedron_Mesh_And_Particles(mattress_grid);
    else tetrahedralized_volume.Initialize_Cube_Mesh_And_Particles(mattress_grid);

    for(int i=0;i<mattress_grid.m;i++)for(int k=0;k<mattress_grid.mn;k++){
        constrained_nodes.Append(1+i+k*mattress_grid.m*mattress_grid.n);
        constrained_nodes.Append(tetrahedralized_volume.particles.Size()+1-(1+i+k*mattress_grid.m*mattress_grid.n));}
   
    tetrahedralized_volume.particles.Delete_Velocity_And_Acceleration();
    tetrahedralized_volume.particles.Delete_Mass(); // in case they're accidently stored in the .tet file
    tetrahedralized_volume.particles.Update_Position_And_Velocity();
    tetrahedralized_volume.particles.Store_Mass();

    tetrahedralized_volume.tetrahedron_mesh.Initialize_Incident_Tetrahedrons();

    tetrahedralized_volume.Update_Bounding_Box();
    VECTOR_3D<T> center(tetrahedralized_volume.bounding_box->Center());T bottom=tetrahedralized_volume.bounding_box->ymin;
    for(int i=0;i<tetrahedralized_volume.particles.array_size;i++){
        tetrahedralized_volume.particles.V(i)=initial_velocity+VECTOR_3D<T>::Cross_Product(initial_angular_velocity,tetrahedralized_volume.particles.X(i)-center);
        tetrahedralized_volume.particles.X(i)=center+initial_orientation.Rotate(tetrahedralized_volume.particles.X(i)-center);
        tetrahedralized_volume.particles.X(i).y+=initial_height-bottom;}
    std::cout << "total vertices = " << tetrahedralized_volume.particles.Size() << std::endl;
    std::cout << "total tets = " << tetrahedralized_volume.tetrahedron_mesh.tetrahedrons.m << std::endl;
    tetrahedralized_volume.Set_Density(1000);
    tetrahedralized_volume.Set_Mass_Of_Particles(solids_parameters.use_constant_mass);
    tetrahedralized_volume.Update_Bounding_Box();
}
//#####################################################################
// Function Set_External_Velocities
//#####################################################################
virtual void Set_External_Velocities(ARRAY<VECTOR_3D<T> ,VECTOR<int,1> >& V,const T time)
{
    for(int i=0;i<constrained_nodes.m;i++) 
        if(i%2 == 1) V(constrained_nodes(i))=VECTOR_3D<T>(0,0,0);
        else V(constrained_nodes(i))=VECTOR_3D<T>(0,.1,0);
}
//#####################################################################
// Function Zero_Out_Enslaved_Velocity_Nodes
//#####################################################################
virtual void Zero_Out_Enslaved_Velocity_Nodes(ARRAY<VECTOR_3D<T> ,VECTOR<int,1> >& V,const T time)
{
    for(int i=0;i<constrained_nodes.m;i++) V(constrained_nodes(i))=VECTOR_3D<T>(0,0,0);
} 
//#####################################################################
// Function Initialize_Rigid_Bodies
//#####################################################################
virtual void Initialize_Rigid_Bodies(ARRAY<RIGID_BODY<TV>*,VECTOR<int,1> >& rigid_bodies)
{    
    std::fstream input;char filename[256];
    
    // plane
    rigid_bodies.Resize(1,rigid_bodies.m+1);rigid_bodies(rigid_bodies.m)=new RIGID_BODY<TV>;
    // triangulated surface
    int index=triangulated_surface_list.Add_Triangulated_Surface();
    sprintf(filename,"../../Public_Data/Rigid_Bodies/ground.tri");input.open(filename,std::ios::in|std::ios::binary);
    assert(input.is_open());triangulated_surface_list.triangulated_surface(index)->Read_Float(input);input.close();
    rigid_bodies(rigid_bodies.m)->Initialize_Triangulated_Surface(*triangulated_surface_list.triangulated_surface(index));
    // implicit surface
    index=implicit_surface_list.Add_Levelset_Implicit_Surface();
    sprintf(filename,"../../Public_Data/Rigid_Bodies/ground.phi");input.open(filename,std::ios::in|std::ios::binary);
    implicit_surface_list.implicit_surface(index)->Read_Float(input);input.close();
    rigid_bodies(rigid_bodies.m)->Initialize_Implicit_Surface(*implicit_surface_list.implicit_surface(index));
    // rigid body
    sprintf(filename,"../../Public_Data/Rigid_Bodies/ground.rgd");input.open(filename,std::ios::in|std::ios::binary);
    rigid_bodies(rigid_bodies.m)->Read_Float(input);input.close();
    rigid_bodies(rigid_bodies.m)->coefficient_of_friction=floor_friction;

    // moving sphere
    rigid_bodies.Resize(1,rigid_bodies.m+1);rigid_bodies(rigid_bodies.m)=new RIGID_BODY<TV>;
    // triangulated surface
    index=triangulated_surface_list.Add_Triangulated_Surface();
    sprintf(filename,"../../Public_Data/Rigid_Bodies/sphere.tri");input.open(filename,std::ios::in|std::ios::binary);
    triangulated_surface_list.triangulated_surface(index)->Read_Float(input);input.close();
    triangulated_surface_list.triangulated_surface(index)->Rescale(.25); // resize radius from the default of 1 to .25
    rigid_bodies(rigid_bodies.m)->Initialize_Triangulated_Surface(*triangulated_surface_list.triangulated_surface(index));
    // implicit surface
    index=implicit_surface_list.Add_Levelset_Implicit_Surface();
    sprintf(filename,"../../Public_Data/Rigid_Bodies/sphere.phi");input.open(filename,std::ios::in|std::ios::binary);
    implicit_surface_list.implicit_surface(index)->Read_Float(input);input.close();
    implicit_surface_list.implicit_surface(index)->Rescale(.25); // resize radius from the default of 1 to .25
    rigid_bodies(rigid_bodies.m)->Initialize_Implicit_Surface(*implicit_surface_list.implicit_surface(index));
    // rigid body
    sprintf(filename,"../../Public_Data/Rigid_Bodies/sphere.rgd");input.open(filename,std::ios::in|std::ios::binary);
    rigid_bodies(rigid_bodies.m)->Read_Float(input);input.close();
    rigid_bodies(rigid_bodies.m)->coefficient_of_friction=0;

    // put it in starting position
    Update_Rigid_Body_Positions(rigid_bodies,0);
    Update_Rigid_Body_Velocities(rigid_bodies,0);
}
//#####################################################################
// Function Update_Rigid_Body_Positions
//#####################################################################
void Update_Rigid_Body_Positions(ARRAY<RIGID_BODY<TV>*,VECTOR<int,1> >& rigid_bodies,const T time)
{
    rigid_bodies(2)->position=VECTOR_3D<T>((T)(5-6*(time)),(T).6,0);
}
//#####################################################################
// Function Update_Rigid_Body_Velocities
//#####################################################################
void Update_Rigid_Body_Velocities(ARRAY<RIGID_BODY<TV>*,VECTOR<int,1> >& rigid_bodies,const T time)
{
    rigid_bodies(2)->velocity=VECTOR_3D<T>(-6,0,0);
}
//##################################################################
};
#endif
