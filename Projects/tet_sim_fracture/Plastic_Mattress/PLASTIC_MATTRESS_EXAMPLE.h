//#####################################################################
// Copyright 2003-2006, Neil Molino, Eftychios Sifakis, Joseph Teran.
// This file is part of PhysBAM whose distribution is governed by the license
// contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class PLASTIC_MATTRESS_EXAMPLE
//#####################################################################
#ifndef __PLASTIC_MATTRESS_EXAMPLE__
#define __PLASTIC_MATTRESS_EXAMPLE__

#include <PhysBAM_Tools/Grids_Uniform/GRID.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/FACE_ARRAYS.h>
#include <PhysBAM_Tools/Math_Tools/constants.h>
#include <PhysBAM_Tools/Matrices/ROTATION.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Constitutive_Models/ROTATED_LINEAR.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Constitutive_Models/SIMPLE_PLASTICITY.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Deformable_Objects/DEFORMABLE_BODY_COLLECTION.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Forces/FINITE_VOLUME.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Forces_And_Torques/GRAVITY.h>
#include "../TET_SIM_FRACTURE_EXAMPLE.h"
using namespace PhysBAM;

template<class T_input>
class PLASTIC_MATTRESS_EXAMPLE:public TET_SIM_FRACTURE_EXAMPLE<T_input>
{
    typedef T_input T;
    typedef VECTOR<T,3> TV;
public:
    typedef TET_SIM_FRACTURE_EXAMPLE<T> BASE;
    using BASE::last_frame;using BASE::frame_rate;using BASE::restart;using BASE::restart_frame;using BASE::output_directory;using BASE::data_directory;using BASE::fluids_parameters;
    using BASE::solids_parameters;using BASE::perform_fracture;using BASE::perturb_amount_for_collision_freeness;using BASE::fracture_object;
    using BASE::plasticity_model;using BASE::tests;

    T initial_height;
    ROTATION<TV> initial_orientation;
    TV initial_velocity,initial_angular_velocity;
    GRID<TV> mattress_grid;
    bool cube_mesh;

    PLASTIC_MATTRESS_EXAMPLE(const STREAM_TYPE stream_type)
        :BASE(stream_type),initial_height((T)1.2),initial_orientation(1,TV(1,1,1)),initial_velocity(0,0,0),initial_angular_velocity(0,0,0),
        mattress_grid(2,2,2,RANGE<TV>(-.5,.5,-.5,.5,-.5,.5)),cube_mesh(true)
    {
        last_frame=240;
        output_directory="Plastic_Mattress/output";

        // Fracture Parameters
        solids_parameters.fracture_evolution->perform_fracture=true;

        // Collision Parameters
        solids_parameters.perform_self_collision=false;
        perturb_amount_for_collision_freeness=(T)1e-5;
        solids_parameters.allow_intersections=true;solids_parameters.allow_intersections_tolerance=(T)1e-7;
        solids_parameters.collisions_nonrigid_collision_attempts=3;
        solids_parameters.perform_collision_body_collisions=true;
        solids_parameters.use_post_cg_constraints=false;
    }

    virtual ~PLASTIC_MATTRESS_EXAMPLE()
    {}

//#####################################################################
// Function Initialize_Tetrahedralized_Volume
//#####################################################################
void Get_Initial_Data()
{
    DEFORMABLE_OBJECT<TV>& deformable_object=solid_body_collection.deformable_object;
    DEFORMABLE_PARTICLES<TV>& particles=deformable_object.particles;
    EMBEDDED_TETRAHEDRALIZED_VOLUME_BOUNDARY_SURFACE<T>& embedding=*EMBEDDED_TETRAHEDRALIZED_VOLUME_BOUNDARY_SURFACE<T>::Create(particles);
    deformable_object.Add_Structure(&embedding);
    deformable_object.collisions.collision_structures.Append(&embedding);
    EMBEDDED_TETRAHEDRALIZED_VOLUME<T>& embedded_object=embedding.embedded_object;
    TETRAHEDRALIZED_VOLUME<T>& tetrahedralized_volume=embedded_object.simplicial_object;

    if(!cube_mesh) tetrahedralized_volume.Initialize_Octahedron_Mesh_And_Particles(mattress_grid);
    else tetrahedralized_volume.Initialize_Cube_Mesh_And_Particles(mattress_grid);
    tetrahedralized_volume.Update_Number_Nodes();
    embedded_object.embedded_particles.Update_Subset_Index_From_Element_Index(); // TODO: make this unnecessary
    particles.Store_Velocity();

    tetrahedralized_volume.Update_Bounding_Box();
    TV center(tetrahedralized_volume.bounding_box->Center());T bottom=tetrahedralized_volume.bounding_box->min_corner.y;
    for(int i=0;i<particles.Size();i++){
        particles.V(i)=initial_velocity+TV::Cross_Product(initial_angular_velocity,particles.X(i)-center);
        particles.X(i)=center+initial_orientation.Rotate(particles.X(i)-center);
        particles.X(i).y+=initial_height-bottom;}
    std::cout << "total vertices = " << particles.Size() << std::endl;
    std::cout << "total tets = " << tetrahedralized_volume.mesh.elements.m << std::endl;
    tetrahedralized_volume.Set_Density(1000);
    tetrahedralized_volume.Set_Mass_Of_Particles(false);

    // rigid bodies
    tests.Add_Ground();
    solids_parameters.collision_body_list.Add_Bodies(deformable_object.rigid_body_particles);
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
    tests.Add_Gravity();
    if(false) plasticity_model=new SIMPLE_PLASTICITY<T,3>(tetrahedralized_volume.mesh.elements.m,2,2);
    solid_body_collection.Add_Force(Create_Finite_Volume(tetrahedralized_volume,new ROTATED_LINEAR<T,3>((T)1e6,(T).3,(T).02),true,(T).1,true,true,false,plasticity_model));

    // initialize_fracture_tetrahedralized_volume
    fracture_object=new FRACTURE_TETRAHEDRALIZED_VOLUME<T>(embedded_object);
    fracture_object->fracture_threshold=TV((T)1e1,(T)1e2,(T)1e3);
    fracture_object->compressive_threshold=TV((T)-1e12,(T)-5e12,(T)-1e12);
    fracture_object->number_of_fracture_initiation_points=4;
    embedded_object.Set_Interpolation_Fraction_Threshold((T)1e-1);

    // initialize rest of fracture
    TET_SIM_FRACTURE_EXAMPLE<T>::Initialize_Bodies();
}
//##################################################################
};
#endif
