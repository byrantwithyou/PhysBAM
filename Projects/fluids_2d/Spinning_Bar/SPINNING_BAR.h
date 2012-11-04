//#####################################################################
// Copyright 2004-2006, Geoffrey Irving, Frank Losasso.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class SPINNING_BAR
//#####################################################################
#ifndef __SPINNING_BAR__
#define __SPINNING_BAR__

#include <PhysBAM_Dynamics/Solids_And_Fluids/SOLIDS_FLUIDS_EXAMPLE_UNIFORM.h>
namespace PhysBAM{

template<class T_input>
class SPINNING_BAR:public SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<VECTOR<T_input,2> > >
{
    typedef T_input T;
    typedef VECTOR<T,2> TV;typedef VECTOR<int,2> TV_INT;
    typedef typename GRID<TV>::CELL_ITERATOR CELL_ITERATOR;typedef typename GRID<TV>::FACE_ITERATOR FACE_ITERATOR;
public:
    typedef SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<TV> > BASE;
    using BASE::first_frame;using BASE::last_frame;using BASE::frame_rate;using BASE::restart;using BASE::restart_frame;using BASE::output_directory;using BASE::Adjust_Phi_With_Sources;
    using BASE::Get_Source_Reseed_Mask;using BASE::Get_Source_Velocities;using BASE::fluids_parameters;using BASE::solids_parameters;using BASE::data_directory;
    using BASE::fluid_collection;using BASE::parse_args;using BASE::resolution;

    T angular_velocity;

    SPINNING_BAR(const STREAM_TYPE stream_type)
        :SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<TV> >(stream_type,1,fluids_parameters.WATER)
    {
    }

    // Unused callbacks
    void Preprocess_Frame(const int frame) PHYSBAM_OVERRIDE {}
    void Adjust_Velocity_With_Objects(const T time){}
    void Postprocess_Solids_Substep(const T time,const int substep) PHYSBAM_OVERRIDE {}
    void Postprocess_Frame(const int frame) PHYSBAM_OVERRIDE {}
    void Postprocess_Phi(const T time) PHYSBAM_OVERRIDE {}
    bool Adjust_Phi_With_Sources(const T time) PHYSBAM_OVERRIDE {return false;}
    void Get_Source_Velocities(const T time) PHYSBAM_OVERRIDE {}
    void Get_Source_Reseed_Mask(ARRAY<bool,VECTOR<int,2> >*& cell_centered_mask,const T time) PHYSBAM_OVERRIDE {}
    void Extrapolate_Phi_Into_Objects(const T time) PHYSBAM_OVERRIDE {}
    void Limit_Dt(T& dt,const T time) PHYSBAM_OVERRIDE {}
    void Align_Deformable_Bodies_With_Rigid_Bodies() PHYSBAM_OVERRIDE {}
    void Preprocess_Solids_Substep(const T time,const int substep) PHYSBAM_OVERRIDE {}
    void Apply_Constraints(const T dt,const T time) PHYSBAM_OVERRIDE {}

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
    fluids_parameters.grid->Initialize(TV_INT(resolution*20,resolution*20),RANGE<TV>(TV(0,0),TV(1,1)));
    first_frame=0;last_frame=2000;
    frame_rate=24;
    restart=false;restart_frame=18;
    output_directory=STRING_UTILITIES::string_sprintf("Spinning_Bar/output_%d",resolution);
    fluids_parameters.domain_walls[0][0]=true;fluids_parameters.domain_walls[0][1]=true;fluids_parameters.domain_walls[1][0]=true;fluids_parameters.domain_walls[1][1]=false; 
    fluids_parameters.number_particles_per_cell=32;
    fluids_parameters.particle_half_bandwidth=1; // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    fluids_parameters.reseeding_frame_rate=1;
    fluids_parameters.bias_towards_negative_particles=true;
    fluids_parameters.incompressible_iterations=20;
    fluids_parameters.write_levelset=true;fluids_parameters.write_velocity=true;fluids_parameters.write_particles=true;
    fluids_parameters.write_removed_positive_particles=false;fluids_parameters.write_removed_negative_particles=false;
    fluids_parameters.write_debug_data=true;fluids_parameters.write_ghost_values=true;
    fluids_parameters.delete_fluid_inside_objects=true;
    //fluids_parameters.enforce_divergence_free_extrapolation=true;

    //fluids_parameters.levelset_substeps=.9;

    fluids_parameters.use_strain=true;fluids_parameters.write_strain=true;
    fluids_parameters.elastic_modulus=1000;
    fluids_parameters.plasticity_alpha=0;
    fluids_parameters.cfl/=10;
    fluids_parameters.gravity=0;
    frame_rate*=5;

    angular_velocity=5;
    //use_external_velocity=true;

    fluids_parameters.viscosity=0;//(T)25;
    if(fluids_parameters.viscosity) output_directory+=STRING_UTILITIES::string_sprintf("_v%g",fluids_parameters.viscosity);
}
void Parse_Late_Options() PHYSBAM_OVERRIDE {BASE::Parse_Late_Options();}
//#####################################################################
// Function Initialize_Advection
//#####################################################################
void Initialize_Advection()    PHYSBAM_OVERRIDE
{
    fluids_parameters.Use_No_Fluid_Coupling_Defaults();
}
//#####################################################################
// Function Initialize_Phi
//#####################################################################
void Initialize_Phi() PHYSBAM_OVERRIDE
{
    GRID<TV>& grid=*fluids_parameters.grid;
    ARRAY<T,VECTOR<int,2> >& phi=fluids_parameters.particle_levelset_evolution->phi;
    RANGE<TV> box((T).25,(T).75,(T).4,(T).6);

    for(int i=0;i<grid.counts.x;i++) for(int j=0;j<grid.counts.y;j++)
        phi(i,j)=box.Signed_Distance(grid.X(TV_INT(i,j)));
}
//#####################################################################
// Function Initialize_Velocities
//#####################################################################
void Initialize_Velocities() PHYSBAM_OVERRIDE
{
    GRID<TV>& grid=*fluids_parameters.grid;
    TV center=grid.Domain().Center();
    for(FACE_ITERATOR iterator(grid);iterator.Valid();iterator.Next()){
        int axis=iterator.Axis();
        fluid_collection.incompressible_fluid_collection.face_velocities.Component(axis)(iterator.Face_Index())=
            TV::Cross_Product(VECTOR<T,1>(angular_velocity),iterator.Location()-center)[axis];}
    //ARRAY<SYMMETRIC_MATRIX<T,2> ,VECTOR<int,2> >::copy(SYMMETRIC_MATRIX<T,2>(1,0,-1),e);
}
//#####################################################################
// Function Get_External_Velocity
//#####################################################################
void Get_External_Velocity(ARRAY<TV,VECTOR<int,2> >& V_blend,ARRAY<T,VECTOR<int,2> >& blend,const T time) PHYSBAM_OVERRIDE
{
    //TV center=grid.Domain().Center();
    //for(int i=0;i<grid.m;i++)for(int j=0;j<grid.n;j++)
    //    V(i,j)=angular_velocity*(grid.X(TV_INT(i,j))-center).Rotate_Counterclockwise_90();
} 
//#####################################################################
};      
}
#endif
