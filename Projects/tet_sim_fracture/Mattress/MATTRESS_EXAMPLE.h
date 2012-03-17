//#####################################################################
// Copyright 2003, Neil Molino, Joseph Teran.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class MATTRESS_EXAMPLE
//#####################################################################
#ifndef __MATTRESS_EXAMPLE__
#define __MATTRESS_EXAMPLE__

#include <PhysBAM_Tools/Math_Tools/constants.h>
#include <PhysBAM_Tools/Matrices/QUATERNION.h>
#include "../TET_SIM_EMBEDDED_EXAMPLE.h"
using namespace PhysBAM;

template <class T>
class MATTRESS_EXAMPLE:public TET_SIM_EMBEDDED_EXAMPLE<T>
{
public:
    T initial_height;
    QUATERNION<T> initial_orientation;
    VECTOR_3D<T> initial_velocity,initial_angular_velocity;
    int m_input,n_input,mn_input;
    T xmin_input,xmax_input,ymin_input,ymax_input,zmin_input,zmax_input; 
    GRID<TV> mattress_grid;
    bool cube_mesh;
    
    MATTRESS_EXAMPLE()
        :initial_height((T)1.2),initial_orientation(0,VECTOR_3D<T>(0,1,0)),
        initial_velocity(0,(T)-sqrt((T)280),0),initial_angular_velocity(0,0,0),
        m_input(8),n_input(8),mn_input(2),           
        xmin_input((T)-.5),xmax_input((T)-xmin_input),ymin_input((T)-.5),ymax_input((T)-ymin_input),zmin_input((T)-.125),zmax_input((T)-zmin_input),
        mattress_grid(m_input,n_input,mn_input,xmin_input,xmax_input,ymin_input,ymax_input,zmin_input,zmax_input),
        cube_mesh(true)
    {
        restart_step_number=0;
        final_time=10;
        frame_rate=24;
        cfl_number=.5;
        solids_parameters.implicit_solve_parameters.cg_iterations=250;
        youngs_modulus=(T)1.9e6;poissons_ratio=(T).4;Rayleigh_coefficient=(T).01;
        strcpy(output_directory,"MATTRESS/output");
        gravity=10;
        solids_parameters.implicit_solve_parameters.cg_tolerance=(T)1e-9;
        
        // Fracture Parameters
        perform_fracture=true;
        fracture_threshold_1=(T)1e3;fracture_threshold_2=(T)fracture_threshold_1;fracture_threshold_3=(T)fracture_threshold_1;
        compressive_threshold_1=(T)-1e12;compressive_threshold_2=(T)compressive_threshold_1;compressive_threshold_3=(T)compressive_threshold_1;
        interpolation_fraction_threshhold=(T).001;
        number_of_fracture_initiation_points=1;
        fracture_bias_propagation=5e6;

        // Collision Parameters
        solids_parameters.perform_self_collision=false;
        push_out=false;
        perturb_amount_for_collision_freeness=(T)1e-5;
        solids_parameters.collisions_nonrigid_collision_attempts=0;//22000;
        floor_friction=1;
    }

    virtual ~MATTRESS_EXAMPLE()
    {}

//#####################################################################
// Function Initialize_Tetrahedralized_Volume
//#####################################################################
void Initialize_Tetrahedralized_Volume(TETRAHEDRALIZED_VOLUME<T>& tetrahedralized_volume)
{
    tetrahedralized_volume.particles.Update_Position();
    if(!cube_mesh) tetrahedralized_volume.Initialize_Octahedron_Mesh_And_Particles(mattress_grid);
    else tetrahedralized_volume.Initialize_Cube_Mesh_And_Particles(mattress_grid);

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
}
//##################################################################
};
#endif
