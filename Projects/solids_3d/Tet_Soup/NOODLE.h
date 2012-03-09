//#####################################################################
// Copyright 2007-2008, Tamar Shinar.
// This file is part of PhysBAM whose distribution is governed by the license 
// contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class NOODLE
//##################################################################### 
#ifndef __NOODLE__
#define __NOODLE__

#include <PhysBAM_Tools/Parsing/PARSE_ARGS.h>
#include <PhysBAM_Geometry/Basic_Geometry/TETRAHEDRON.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Forces/TRIANGLE_BENDING_ELEMENTS.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Forces/TRIANGLE_BENDING_SPRINGS.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Forces_And_Torques/ETHER_DRAG.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Standard_Tests/SOLIDS_STANDARD_TESTS.h>
#include <PhysBAM_Dynamics/Solids_And_Fluids/SOLIDS_FLUIDS_EXAMPLE_UNIFORM.h>
#include <Constitutive_Models/DIAGONALIZED_NEO_HOOKEAN_3D.h>
#include <Forces_And_Torques/DIAGONALIZED_FINITE_VOLUME_3D.h>
#include <Solids_And_Fluids/SOLIDS_PARAMETERS_3D.h>
namespace PhysBAM{

template<class T,class RW>
class NOODLE:public SOLIDS_FLUIDS_EXAMPLE_UNIFORM<T,GRID<TV>,RW>
{
    typedef VECTOR<T,3> TV;
public:
    typedef SOLIDS_FLUIDS_EXAMPLE_UNIFORM<T,GRID<TV>,RW> BASE;
    using BASE::first_frame;using BASE::last_frame;using BASE::frame_rate;using BASE::restart;using BASE::restart_frame;using BASE::output_directory;using BASE::solids_parameters;
    using BASE::write_last_frame;using BASE::data_directory;using BASE::solids_parameters;

    SOLIDS_STANDARD_TESTS<TV,RW> tests;
    T ribbon_length,ribbon_width;
    T bowl_scale;

    NOODLE(const PARSE_ARGS& parse_args):
        SOLIDS_FLUIDS_EXAMPLE_UNIFORM<T,GRID<TV>,RW>(0,FLUIDS_PARAMETERS_UNIFORM<T,GRID<TV> >::NONE),tests(*this,solids_parameters),ribbon_length(250),ribbon_width(2.5),bowl_scale(10)
    {
        solids_parameters.rigid_body_evolution_parameters.simulate_rigid_bodies=true;
        solids_parameters.cfl=(T).2;
        last_frame=10000;
        frame_rate=24;
        output_directory="Tet_Soup/noodle/output";
        std::cout << "Frame rate: "<<frame_rate<<std::endl;
        solids_parameters.implicit_solve_parameters.cg_tolerance=1e-6;
        solids_parameters.implicit_solve_parameters.cg_iterations=500;
        solids_parameters.perform_self_collision=true;
        solids_parameters.collide_with_interior=false;
    }

    // unused callbacks
    void Update_Solids_Parameters(const T time) PHYSBAM_OVERRIDE {}
    void Limit_Solids_Dt(T& dt,const T time) PHYSBAM_OVERRIDE {}
    void Postprocess_Solids_Substep(const T time,const int substep) PHYSBAM_OVERRIDE {}
    void Apply_Constraints(const T dt,const T time) PHYSBAM_OVERRIDE {}
    void Set_External_Velocities(ARRAY_VIEW<TWIST<TV> > twist,const T velocity_time,const T current_position_time) PHYSBAM_OVERRIDE {}
    void Zero_Out_Enslaved_Velocity_Nodes(ARRAY_VIEW<TWIST<TV> > twist,const T velocity_time,const T current_position_time) PHYSBAM_OVERRIDE {}

//#####################################################################
// Initialize_Bodies
//#####################################################################
void Initialize_Bodies() PHYSBAM_OVERRIDE
{
    DEFORMABLE_OBJECT<T,TV>& deformable_object=solid_body_collection.deformable_object;
    DEFORMABLE_PARTICLES<T,TV>& particles=deformable_object.particles;

    RIGID_BODY<TV>& rigid_body=tests.Add_Rigid_Body("bowl",bowl_scale,0,true);
    rigid_body.frame.r=QUATERNION<T>::Rotation_Quaternion(-TV::Axis_Vector(3),-TV::Axis_Vector(2));
    rigid_body.frame.t=TV(0,-.5,0);

    // create the cloth ribbon
    RIGID_BODY_STATE<TV> state(FRAME_3D<T>(TV(0,.7*ribbon_length,0),QUATERNION<T>::From_Euler_Angles(-.5*pi,0,0)));
    TRIANGULATED_SURFACE<T>& cloth_panel=tests.Create_Cloth_Panel((int)ribbon_length*2,ribbon_length,ribbon_width/ribbon_length,&state);
    TV X=particles.X(1);for(int p=0;p<particles.array_collection->Size();p++){particles.X(p)-=X;particles.V(p)=TV(0,-10,0);}

    // correct number nodes
    for(int i=0;i<deformable_object.structures.m;i++) deformable_object.structures(i)->Update_Number_Nodes();

    // collisions
    deformable_object.collisions.collision_structures.Append_Elements(deformable_object.structures);
    solid_body_collection.deformable_body_collection.triangle_repulsions_and_collisions_geometry.structures.Append_Elements(deformable_object.structures);
    solids_parameters.collision_body_list.Add_Bodies(solids_parameters.rigid_body_parameters.list);

    // correct mass
    particles.Compute_Auxiliary_Attributes(solid_body_collection.deformable_body_collection.soft_bindings);

    // add forces
    solid_body_collection.Add_Force(Create_Edge_Springs<T>(cloth_panel,(T)1e3,1));
    TRIANGLE_BENDING_ELEMENTS<T>* bend=Create_Bending_Elements(cloth_panel,(T)1e3);
    bend->Set_Area_Cutoff_From_Triangulated_Surface(cloth_panel,.1);
    bend->use_force_differential=false;
    solid_body_collection.Add_Force(bend);

    solid_body_collection.Add_Force(new GRAVITY<T,TV>(cloth_panel,2));
    ETHER_DRAG<T,GRID<TV> >& ether_drag=*new ETHER_DRAG<T,GRID<TV> >(particles,.1);
    solid_body_collection.Add_Force(&ether_drag);

    solid_body_collection.Update_Fragments();

    LOG::cout<<"NUMBER of FRAGMENTS="<<solid_body_collection.deformable_object.particles_of_fragment.m<<std::endl;
    SOLIDS_FLUIDS_EXAMPLE_UNIFORM<T,GRID<TV>,RW>::Initialize_Bodies();
}
//#####################################################################
};
}
#endif
