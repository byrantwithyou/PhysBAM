//#####################################################################
// Copyright 2009, Avi Robinson-Mosher.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class WATERFALL 
//##################################################################### 
#ifndef __WATERFALL__
#define __WATERFALL__

#include <Core/Math_Tools/cube.h>
#include <Core/Random_Numbers/NOISE.h>
#include <Grid_Tools/Grids/NODE_ITERATOR.h>
#include <Geometry/Basic_Geometry/CYLINDER.h>
#include <Incompressible/Boundaries/BOUNDARY_OPEN_CALLBACKS.h>
#include <Dynamics/Incompressible_Flows/SPH_EVOLUTION_UNIFORM.h>
#include <Dynamics/Solids_And_Fluids/SOLIDS_FLUIDS_EXAMPLE_UNIFORM.h>
namespace PhysBAM{

template<class T_input,class TV=VECTOR<T_input,3> >
class WATERFALL:public SOLIDS_FLUIDS_EXAMPLE_UNIFORM<TV>,public BOUNDARY_OPEN_CALLBACKS<TV>
{
    typedef T_input T;
    typedef VECTOR<int,TV::m> TV_INT;
    typedef ARRAY<PARTICLE_LEVELSET_REMOVED_PARTICLES<TV>*,TV_INT> T_ARRAYS_PARTICLE_LEVELSET_REMOVED_PARTICLES;
public:
    typedef SOLIDS_FLUIDS_EXAMPLE_UNIFORM<TV> BASE;
    using BASE::first_frame;using BASE::last_frame;using BASE::frame_rate;using BASE::restart;using BASE::restart_frame;using BASE::output_directory;using BASE::data_directory;
    using BASE::fluids_parameters;using BASE::solids_parameters;using BASE::stream_type;using BASE::solid_body_collection;using BASE::test_number;using BASE::resolution;
    using BASE::user_last_frame;

    RIGID_BODY_COLLECTION<TV>& rigid_body_collection;
    int waterfall;

    FLUID_COLLISION_BODY_INACCURATE_UNION<TV> inaccurate_union;

    RANGE<TV> source_box;
    TV source_velocity;
    MATRIX<T,4> world_to_source;

    WATERFALL(const STREAM_TYPE stream_type_input,PARSE_ARGS& parse_args)
        :SOLIDS_FLUIDS_EXAMPLE_UNIFORM<TV>(stream_type_input,parse_args,1,fluids_parameters.WATER),
        rigid_body_collection(solid_body_collection.rigid_body_collection),inaccurate_union(*fluids_parameters.grid)
    {
        parse_args.Parse();

        first_frame=0;
        if(!user_last_frame) last_frame=1000;
        if(!this->user_frame_rate) frame_rate=36;
        restart=false;restart_frame=18;
        int cells=1*resolution;
        fluids_parameters.grid->Initialize(TV_INT(2*cells+1,3*cells+1,2*cells+1),RANGE<TV>(TV(-4,-8,0),TV(4,4,8)));
        fluids_parameters.domain_walls[0][0]=false;fluids_parameters.domain_walls[0][1]=true;fluids_parameters.domain_walls[1][0]=true;
        fluids_parameters.domain_walls[1][1]=false;fluids_parameters.domain_walls[2][0]=true;fluids_parameters.domain_walls[2][1]=true;
        fluids_parameters.number_particles_per_cell=32;
        fluids_parameters.viscosity=(T)0;fluids_parameters.implicit_viscosity=false;
        fluids_parameters.write_levelset=true;fluids_parameters.write_velocity=true;fluids_parameters.write_particles=true;
        fluids_parameters.use_removed_positive_particles=true;fluids_parameters.use_removed_negative_particles=true;
        fluids_parameters.write_removed_positive_particles=true;fluids_parameters.write_removed_negative_particles=true;
        fluids_parameters.write_debug_data=true;
        if(!this->user_output_directory)
            output_directory=LOG::sprintf("Waterfall/Test_%d_Waterfall_Resolution_%d_%d_%d",test_number,(fluids_parameters.grid->counts.x-1),(fluids_parameters.grid->counts.y-1),
            (fluids_parameters.grid->counts.z-1));
        fluids_parameters.delete_fluid_inside_objects=true;
        fluids_parameters.enforce_divergence_free_extrapolation=false;
        fluids_parameters.store_particle_ids=true;
        fluids_parameters.second_order_cut_cell_method=true;
        fluids_parameters.cfl=(T)1.9;
        fluids_parameters.incompressible_iterations=40;
        fluids_parameters.solid_affects_fluid=true;

        // MacCormack parameters
        fluids_parameters.use_maccormack_semi_lagrangian_advection=true;
        fluids_parameters.use_maccormack_compute_mask=true;
        fluids_parameters.use_maccormack_for_incompressible=true;
        fluids_parameters.bandwidth_without_maccormack_near_interface=1;

        // Source parameters
        source_velocity=TV(0,-1,1);
        source_box=RANGE<TV>(TV(-2,2,-1),TV(2,3.9,1.5));
        world_to_source=MATRIX<T,4>::Identity_Matrix();
    }

    ~WATERFALL() 
    {}

//#####################################################################
// Function Get_Source_Velocities
//#####################################################################
void Get_Source_Velocities(ARRAY<T,FACE_INDEX<3> >& face_velocities,ARRAY<bool,FACE_INDEX<3> >& psi_N,const T time) override
{
    BASE::Get_Source_Velocities(source_box,world_to_source,source_velocity);
}
//#####################################################################
// Function Initialize_Phi
//#####################################################################
void Initialize_Phi() override
{
    ARRAY<T,VECTOR<int,3> >& phi=fluids_parameters.particle_levelset_evolution->phi;
    phi.Fill(1);
}
//#####################################################################
// Function Adjust_Phi_With_Sources
//#####################################################################
bool Adjust_Phi_With_Sources(const T time) override
{
    BASE::Adjust_Phi_With_Source(source_box,world_to_source);
    return true;
}
//#####################################################################
// Function Initialize_Advection
//#####################################################################
void Initialize_Advection() override
{
    fluids_parameters.Use_Fluid_Coupling_Defaults();
}
//#####################################################################
// Function Initial_Phi_Object
//#####################################################################
T Initial_Phi_Object(const VECTOR<T,3>& X) const
{
    return rigid_body_collection.Rigid_Body(waterfall).Implicit_Geometry_Extended_Value(X);
}
//#####################################################################
// Function Initialize_Bodies
//#####################################################################
void Initialize_Bodies() override
{    
    waterfall=rigid_body_collection.Add_Rigid_Body(data_directory+"/Rigid_Bodies/fractal_box_mesh",(T)2,true,true,false);
    rigid_body_collection.rigid_body_particles.frame(waterfall).t=VECTOR<T,3>((T)25,(T)-2,(T)30);rigid_body_collection.Rigid_Body(waterfall).Is_Kinematic()=true;
    inaccurate_union.collision_bodies.Add_Bodies(rigid_body_collection);
    fluids_parameters.collision_bodies_affecting_fluid->collision_geometry_collection.Add_Body(&inaccurate_union,0,false);
}
//#####################################################################
};
}
#endif
