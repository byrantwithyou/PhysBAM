//#####################################################################
// Copyright 2003-2006, Geoffrey Irving, Neil Molino, Eftychios Sifakis, Joseph Teran.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class GRAIN_MATTRESS_EXAMPLE
//#####################################################################
#ifndef __GRAIN_MATTRESS_EXAMPLE__
#define __GRAIN_MATTRESS_EXAMPLE__

#include <PhysBAM_Tools/Math_Tools/constants.h>
#include <PhysBAM_Tools/Matrices/QUATERNION.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Constitutive_Models/ROTATED_LINEAR.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Forces/FINITE_VOLUME.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Forces_And_Torques/GRAVITY.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Fracture/GRAIN_BOUNDARIES.h>
#include "../TET_SIM_FRACTURE_EXAMPLE.h"
using namespace PhysBAM;

template<class T_input>
class GRAIN_MATTRESS_EXAMPLE:public TET_SIM_FRACTURE_EXAMPLE<T_input>
{
    typedef T_input T;
    typedef VECTOR<T,3> TV;
public:
    typedef TET_SIM_FRACTURE_EXAMPLE<T> BASE;
    using BASE::last_frame;using BASE::frame_rate;using BASE::output_directory;using BASE::data_directory;using BASE::fluids_parameters;
    using BASE::solids_parameters;using BASE::perform_fracture;using BASE::perturb_amount_for_collision_freeness;using BASE::fracture_object;
    using BASE::stream_type;using BASE::tests;

    T initial_height;
    ROTATION<TV> initial_orientation;
    TV initial_velocity,initial_angular_velocity;
    GRID<TV> mattress_grid;
    bool cube_mesh;   
    
    GRAIN_MATTRESS_EXAMPLE(const STREAM_TYPE stream_type)
        :BASE(stream_type),initial_height((T)1.2),initial_orientation((T).5*(T)pi,TV(1,0,0)),initial_velocity(0,(T)-10,0),initial_angular_velocity(0,0,0),
        mattress_grid(8,8,2,RANGE<TV>(-(T).5,(T).5,-(T).5,(T).5,-(T)1/14,(T)1/14)),cube_mesh(true)
    {
        frame_rate=50;
        last_frame=(int)(5*frame_rate);
        solids_parameters.cfl=(T).5;
        solids_parameters.implicit_solve_parameters.cg_iterations=250;
        solids_parameters.implicit_solve_parameters.cg_tolerance=(T)1e-3;
        output_directory="Grain_Mattress/output";        
        
        // Fracture Parameters
        perform_fracture=true;

        // Collision Parameters
        solids_parameters.perform_self_collision=false;
        perturb_amount_for_collision_freeness=(T)1e-5;
        solids_parameters.allow_intersections=true;solids_parameters.allow_intersections_tolerance=(T)1e-7;
        solids_parameters.collisions_nonrigid_collision_attempts=3;
    }

    virtual ~GRAIN_MATTRESS_EXAMPLE()
    {}

//#####################################################################
// Function Initialize_Tetrahedralized_Volume
//#####################################################################
void Get_Initial_Data()
{
    DEFORMABLE_OBJECT<TV>& deformable_object=solid_body_collection.deformable_object;
    PARTICLES<TV>& particles=deformable_object.particles;
    EMBEDDED_TETRAHEDRALIZED_VOLUME_BOUNDARY_SURFACE<T>& embedding=*EMBEDDED_TETRAHEDRALIZED_VOLUME_BOUNDARY_SURFACE<T>::Create(particles);
    deformable_object.Add_Structure(&embedding);
    EMBEDDED_TETRAHEDRALIZED_VOLUME<T>& embedded_object=embedding.embedded_object;
    TETRAHEDRALIZED_VOLUME<T>& tetrahedralized_volume=embedded_object.simplicial_object;

    if(!cube_mesh) tetrahedralized_volume.Initialize_Octahedron_Mesh_And_Particles(mattress_grid);
    else tetrahedralized_volume.Initialize_Cube_Mesh_And_Particles(mattress_grid);
    embedded_object.embedded_particles.Update_Subset_Index_From_Element_Index(); // TODO: make this unnecessary
    embedded_object.embedded_mesh.number_nodes=particles.array_collection->Size();
    particles.Store_Velocity();

    tetrahedralized_volume.Update_Bounding_Box();
    TV center(tetrahedralized_volume.bounding_box->Center());T bottom=tetrahedralized_volume.bounding_box->min_corner.y;
    for(int i=0;i<particles.array_size;i++){
        particles.V(i)=initial_velocity+TV::Cross_Product(initial_angular_velocity,particles.X(i)-center);
        particles.X(i)=center+initial_orientation.Rotate(particles.X(i)-center);
        particles.X(i).y+=initial_height-bottom;}
    LOG::cout << "total vertices = " << particles.array_collection->Size() << std::endl;
    LOG::cout << "total tets = " << tetrahedralized_volume.mesh.elements.m << std::endl;
    tetrahedralized_volume.Set_Density(1000);
    tetrahedralized_volume.Set_Mass_Of_Particles();

    // rigid bodies
    tests.Add_Ground();
    tests.Add_Rigid_Body("sphere",(T).5,0).is_static=true;;
    
    solids_parameters.collision_body_list.Add_Bodies(solid_body_collection.deformable_object.rigid_body_particles);
}
//#####################################################################
// Function Initialize_Bodies
//#####################################################################
void Initialize_Bodies() PHYSBAM_OVERRIDE
{
    Get_Initial_Data();

    // initialize forces

    DEFORMABLE_OBJECT<TV>& deformable_object=solid_body_collection.deformable_object;
    EMBEDDED_TETRAHEDRALIZED_VOLUME_BOUNDARY_SURFACE<T>& embedding=deformable_object.template Find_Structure<EMBEDDED_TETRAHEDRALIZED_VOLUME_BOUNDARY_SURFACE<T>&>();
    EMBEDDED_TETRAHEDRALIZED_VOLUME<T>& embedded_object=embedding.embedded_object;
    TETRAHEDRALIZED_VOLUME<T>& tetrahedralized_volume=embedded_object.simplicial_object;

    solid_body_collection.Add_Force(new GRAVITY<TV>(deformable_object.particles,deformable_object.rigid_body_particles,true,true));
    solid_body_collection.Add_Force(Create_Finite_Volume(deformable_object.particles,deformable_object.rigid_body_particles,tetrahedralized_volume,
        new ROTATED_LINEAR<T,3>((T)1e6,(T).4,(T).01),true,(T).1,true,true));

    // initialize_fracture_tetrahedralized_volume
    fracture_object=new FRACTURE_TETRAHEDRALIZED_VOLUME<T>(embedded_object);
    fracture_object->fracture_threshold=TV((T)1e3,(T)5e3,(T)1e5);
    fracture_object->compressive_threshold=TV(-(T)1e4,-(T)5e4,-(T)1e6);
    fracture_object->embedded_object.interpolation_fraction_threshold=(T).2;
    fracture_object->number_of_fracture_initiation_points=1;

    // Initialize grain boundaries
    fracture_object->bias_stress=true;
    fracture_object->fracture_bias_direction_coefficient=1;
    fracture_object->eigenvector_coefficient=0;
    fracture_object->fracture_bias_propagation_coefficient=0;
    fracture_object->Set_Fracture_Bias_Stress_Scaling(1);
    fracture_object->Set_Fracture_Bias_Magnitude((T)1e5);
    fracture_object->Set_Fracture_Bias_Propagation(0);

    GRAIN_BOUNDARIES<TV,3> grain_boundaries(*fracture_object,tetrahedralized_volume.mesh);
    int number_of_regions=2;
    grain_boundaries.Fill_Node_Regions_With_Uniform_Vectors(number_of_regions);
    //grain_boundaries.Print_Number_In_Regions(number_of_regions);
    for(int t=0;t<fracture_object->fracture_bias_direction.m;t++) LOG::cout<<"fracture_bias_direction("<<t<<")="<<fracture_object->fracture_bias_direction(t)<<std::endl;
    //grain_boundaries.Smooth_Fracture_Bias_Directions(9);

    deformable_object.collisions.collision_structures.Append_Elements(deformable_object.structures);
    // initialize rest of fracture
    TET_SIM_FRACTURE_EXAMPLE<T>::Initialize_Bodies();
}
//##################################################################
};
#endif
