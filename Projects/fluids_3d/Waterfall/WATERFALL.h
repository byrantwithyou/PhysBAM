//#####################################################################
// Copyright 2009, Avi Robinson-Mosher.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class WATERFALL 
//##################################################################### 
#ifndef __WATERFALL__
#define __WATERFALL__

#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_NODE.h>
#include <PhysBAM_Tools/Math_Tools/cube.h>
#include <PhysBAM_Tools/Random_Numbers/NOISE.h>
#include <PhysBAM_Tools/Read_Write/Data_Structures/READ_WRITE_PAIR.h>
#include <PhysBAM_Geometry/Basic_Geometry/CYLINDER.h>
#include <PhysBAM_Dynamics/Boundaries/BOUNDARY_OPEN_CALLBACKS.h>
#include <PhysBAM_Dynamics/Boundaries/BOUNDARY_OPEN_WATER.h>
#include <PhysBAM_Dynamics/Incompressible_Flows/SPH_EVOLUTION_UNIFORM.h>
#include <PhysBAM_Dynamics/Solids_And_Fluids/SOLIDS_FLUIDS_EXAMPLE_UNIFORM.h>
namespace PhysBAM{

template<class T_input,class T_GRID=GRID<VECTOR<T_input,3> > >
class WATERFALL:public SOLIDS_FLUIDS_EXAMPLE_UNIFORM<T_GRID>,public BOUNDARY_OPEN_CALLBACKS<T_GRID>
{
    typedef T_input T;
    typedef typename GRID_ARRAYS_POLICY<T_GRID>::FACE_ARRAYS T_FACE_ARRAYS_SCALAR;typedef typename T_FACE_ARRAYS_SCALAR::template REBIND<bool>::TYPE T_FACE_ARRAYS_BOOL;
    typedef typename T_GRID::FACE_ITERATOR FACE_ITERATOR;typedef typename T_GRID::VECTOR_INT TV_INT;typedef VECTOR<T,3> TV;
    typedef typename GRID_ARRAYS_POLICY<GRID<TV> >::ARRAYS_SCALAR::template REBIND<PARTICLE_LEVELSET_REMOVED_PARTICLES<TV>*>::TYPE T_ARRAYS_PARTICLE_LEVELSET_REMOVED_PARTICLES;
    typedef typename GRID<TV>::CELL_ITERATOR CELL_ITERATOR;typedef typename GRID<TV>::NODE_ITERATOR NODE_ITERATOR;
public:
    typedef SOLIDS_FLUIDS_EXAMPLE_UNIFORM<T_GRID> BASE;
    using BASE::first_frame;using BASE::last_frame;using BASE::frame_rate;using BASE::restart;using BASE::restart_frame;using BASE::output_directory;using BASE::data_directory;
    using BASE::fluids_parameters;using BASE::solids_parameters;using BASE::stream_type;using BASE::solid_body_collection;using BASE::test_number;using BASE::resolution;

    RIGID_BODY_COLLECTION<TV>& rigid_body_collection;
    int waterfall;

    FLUID_COLLISION_BODY_INACCURATE_UNION<T_GRID> inaccurate_union;

    RANGE<TV> source_box;
    TV source_velocity;
    MATRIX<T,4> world_to_source;

    WATERFALL(const STREAM_TYPE stream_type)
        :SOLIDS_FLUIDS_EXAMPLE_UNIFORM<T_GRID>(stream_type,1,fluids_parameters.WATER),
        rigid_body_collection(solid_body_collection.rigid_body_collection),inaccurate_union(*fluids_parameters.grid)
    {
    }

    ~WATERFALL() 
    {}

    // Unused callbacks
    void Preprocess_Frame(const int frame) PHYSBAM_OVERRIDE {}
    void Limit_Dt(T& dt,const T time) PHYSBAM_OVERRIDE {}
    void Postprocess_Phi(const T time) PHYSBAM_OVERRIDE {}
    void Apply_Constraints(const T dt,const T time) PHYSBAM_OVERRIDE {}
    void Postprocess_Solids_Substep(const T time,const int substep) PHYSBAM_OVERRIDE {}
    void Limit_Solids_Dt(T& dt,const T time) PHYSBAM_OVERRIDE {}
    void Get_Source_Reseed_Mask(ARRAY<bool,VECTOR<int,3> >*& cell_centered_mask,const T time) PHYSBAM_OVERRIDE {}

//#####################################################################
// Function Register_Options
//#####################################################################
void Register_Options() PHYSBAM_OVERRIDE
{
    BASE::Register_Options();
}
//#####################################################################
// Function Parse_Options
//#####################################################################
void Parse_Options() PHYSBAM_OVERRIDE
{
    BASE::Parse_Options();
    first_frame=0;last_frame=1000;
    frame_rate=36;
    restart=false;restart_frame=18;
    int cells=1*resolution;
    fluids_parameters.grid->Initialize(2*cells+1,3*cells+1,2*cells+1,-4,4,-8,4,0,8);
    fluids_parameters.domain_walls[0][0]=false;fluids_parameters.domain_walls[0][1]=true;fluids_parameters.domain_walls[1][0]=true;
    fluids_parameters.domain_walls[1][1]=false;fluids_parameters.domain_walls[2][0]=true;fluids_parameters.domain_walls[2][1]=true;
    fluids_parameters.number_particles_per_cell=32;
    fluids_parameters.viscosity=(T)0;fluids_parameters.implicit_viscosity=false;
    fluids_parameters.write_levelset=true;fluids_parameters.write_velocity=true;fluids_parameters.write_particles=true;
    fluids_parameters.use_removed_positive_particles=true;fluids_parameters.use_removed_negative_particles=true;
    fluids_parameters.write_removed_positive_particles=true;fluids_parameters.write_removed_negative_particles=true;
    fluids_parameters.write_debug_data=true;
    output_directory=STRING_UTILITIES::string_sprintf("Waterfall/Test_%d_Waterfall_Resolution_%d_%d_%d",test_number,(fluids_parameters.grid->counts.x-1),(fluids_parameters.grid->counts.y-1),
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
void Parse_Late_Options() PHYSBAM_OVERRIDE {BASE::Parse_Late_Options();}
//#####################################################################
// Function Get_Source_Velocities
//#####################################################################
void Get_Source_Velocities(ARRAY<T,FACE_INDEX<3> >& face_velocities,ARRAY<bool,FACE_INDEX<3> >& psi_N,const T time) PHYSBAM_OVERRIDE
{
    BASE::Get_Source_Velocities(source_box,world_to_source,source_velocity);
}
//#####################################################################
// Function Initialize_Phi
//#####################################################################
void Initialize_Phi() PHYSBAM_OVERRIDE
{
    ARRAY<T,VECTOR<int,3> >& phi=fluids_parameters.particle_levelset_evolution->phi;
    phi.Fill(1);
}
//#####################################################################
// Function Adjust_Phi_With_Sources
//#####################################################################
bool Adjust_Phi_With_Sources(const T time) PHYSBAM_OVERRIDE
{
    BASE::Adjust_Phi_With_Source(source_box,world_to_source);
    return true;
}
//#####################################################################
// Function Initialize_Advection
//#####################################################################
void Initialize_Advection() PHYSBAM_OVERRIDE
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
void Initialize_Bodies() PHYSBAM_OVERRIDE
{    
    waterfall=rigid_body_collection.Add_Rigid_Body(stream_type,data_directory+"/Rigid_Bodies/fractal_box_mesh",(T)2,true,true,false);
    rigid_body_collection.rigid_body_particle.frame(waterfall).t=VECTOR<T,3>((T)25,(T)-2,(T)30);rigid_body_collection.Rigid_Body(waterfall).Is_Kinematic()=true;
    inaccurate_union.collision_bodies.Add_Bodies(rigid_body_collection.rigid_geometry_collection);
    fluids_parameters.collision_bodies_affecting_fluid->collision_geometry_collection.Add_Body(&inaccurate_union,0,false);
}
//#####################################################################
};
}
#endif
