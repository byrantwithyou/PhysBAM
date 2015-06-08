//#####################################################################
// Copyright 2006-2007, Tamar Shinar.
// This file is part of PhysBAM whose distribution is governed by the license 
// contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class DSD_FIRE_BALL_EXAMPLE  
//##################################################################### 
#ifndef __DSD_FIRE_BALL_EXAMPLE__
#define __DSD_FIRE_BALL_EXAMPLE__

#include <Tools/Grids_Uniform_Advection/ADVECTION_MACCORMACK_UNIFORM.h>
#include <Tools/Interpolation/INTERPOLATION_CURVE.h>
#include <Geometry/Basic_Geometry/SPHERE.h>
#include <Dynamics/Incompressible_Flows/DETONATION_SHOCK_DYNAMICS.h>
#include <Dynamics/Solids_And_Fluids/SOLIDS_FLUIDS_EXAMPLE_UNIFORM.h>
using namespace PhysBAM;

template<class T_input>
class DSD_FIRE_BALL_EXAMPLE:public SOLIDS_FLUIDS_EXAMPLE_UNIFORM<VECTOR<T_input,3> >
{
    typedef T_input T;
public:
    typedef VECTOR<T,3> TV;typedef VECTOR<int,3> TV_INT;

    typedef SOLIDS_FLUIDS_EXAMPLE_UNIFORM<TV> BASE;
    using BASE::fluids_parameters;using BASE::fluid_collection;using BASE::solids_parameters;using BASE::first_frame;using BASE::last_frame;using BASE::frame_rate;using BASE::write_output_files;
    using BASE::output_directory;using BASE::restart;using BASE::restart_frame;using BASE::data_directory;

    SPHERE<TV> source_sphere;
    T source_end_time;
    T normal_velocity;
    T Dn_initial;

    DSD_FIRE_BALL_EXAMPLE(const STREAM_TYPE stream_type_input,PARSE_ARGS& parse_args)
        :SOLIDS_FLUIDS_EXAMPLE_UNIFORM<TV>(stream_type_input,parse_args,2,fluids_parameters.FIRE),source_sphere(TV(1,(T)0.5,1),(T)0.2),source_end_time((T)3),normal_velocity(4),Dn_initial((T)0.1)
    {
        frame_rate=96;
        fluids_parameters.grid->Initialize(TV_INT(151,226,151),RANGE<TV>(TV(),TV(2,3,2)));
        fluids_parameters.domain_walls[2][1]=fluids_parameters.domain_walls[2][0]=fluids_parameters.domain_walls[0][0]=fluids_parameters.domain_walls[0][1]=fluids_parameters.domain_walls[1][1]=false;fluids_parameters.domain_walls[1][0]=true;
        last_frame=int(T(20)*frame_rate);
        fluids_parameters.incompressible_iterations=200;
        fluids_parameters.gravity=TV();
        fluids_parameters.density_container.Set_Ambient_Density(0);
        fluids_parameters.temperature_container.Set_Cooling_Constant((T)5000);
        fluids_parameters.temperature_container.Set_Ambient_Temperature((T)283.15);
        fluids_parameters.temperature_products=2800;fluids_parameters.temperature_fuel=2500;
        fluids_parameters.temperature_container.Use_Maccormack_Advection(fluids_parameters.maccormack_cell_mask);
        fluids_parameters.density_container.Use_Maccormack_Advection(fluids_parameters.maccormack_cell_mask);
        fluids_parameters.density_buoyancy_constant=fluids_parameters.temperature_buoyancy_constant=(T)0.055;
        fluids_parameters.use_body_force=true;
//        fluids_parameters.confinement_parameters(1)=.2;
//        fluids_parameters.confinement_parameters(2)=.1;
        fluids_parameters.use_reacting_flow=true;
        fluids_parameters.use_dsd=true;
        fluids_parameters.densities(1)=(T)3;
        fluids_parameters.densities(2)=(T)1;
        fluids_parameters.fuel_region(1)=true;
        fluids_parameters.fuel_region(2)=false;
        fluids_parameters.normal_flame_speeds(1,2)=fluids_parameters.normal_flame_speeds(2,1)=1;
        fluids_parameters.delete_fluid_inside_objects=true;
        write_output_files=true;
        fluids_parameters.write_debug_data=true;
        fluids_parameters.write_velocity=true;
        fluids_parameters.write_particles=true;
        fluids_parameters.number_particles_per_cell=0;
        output_directory="DSD_Fire_Ball/output";
        fluids_parameters.use_maccormack_semi_lagrangian_advection=true;
        fluids_parameters.use_maccormack_for_level_set=true;
        fluids_parameters.use_maccormack_compute_mask=true;
        parse_args.Parse();
    }

    virtual ~DSD_FIRE_BALL_EXAMPLE()
    {}

    // unused callbacks
    void Limit_Solids_Dt(T& dt,const T time) PHYSBAM_OVERRIDE {}
    void Postprocess_Solids_Substep(const T time,const int substep) PHYSBAM_OVERRIDE {}
    void Apply_Constraints(const T dt,const T time) PHYSBAM_OVERRIDE {}
    void Preprocess_Frame(const int frame) PHYSBAM_OVERRIDE {}
    void Extrapolate_Phi_Into_Objects(const T time) PHYSBAM_OVERRIDE {}
    void Get_Source_Reseed_Mask(ARRAY<bool,TV_INT>*& cell_centered_mask,const T time) PHYSBAM_OVERRIDE {}
    void Postprocess_Phi(const T time) PHYSBAM_OVERRIDE {}
    void Postprocess_Frame(const int frame) PHYSBAM_OVERRIDE {}
    void Adjust_Density_And_Temperature_With_Sources(const T time) PHYSBAM_OVERRIDE {}
    void Limit_Dt(T& dt,const T time) PHYSBAM_OVERRIDE {}

void After_Initialization() PHYSBAM_OVERRIDE {BASE::After_Initialization();}
//#####################################################################
// Function Get_Flame_Speed_Multiplier
//#####################################################################
void Get_Flame_Speed_Multiplier(const T dt,const T time) PHYSBAM_OVERRIDE
{
}
//#####################################################################
// Function Update_Fluid_Parameters
//#####################################################################
void Update_Fluid_Parameters(const T dt,const T time) PHYSBAM_OVERRIDE
{
    SOLIDS_FLUIDS_EXAMPLE_UNIFORM<TV>::Update_Fluid_Parameters(dt,time);
    DETONATION_SHOCK_DYNAMICS<TV>& dsd=*fluids_parameters.incompressible_multiphase->projection.dsd;
    dsd.order=3;
    dsd.Dcj=(T)0.2;
    dsd.Dcj_min_clamp=(T)0.0;
    dsd.Dcj_max_clamp=1000;
    dsd.A_coeff=10;
    dsd.B_coeff=(T).1;
    dsd.C_coeff=100;
    dsd.D_coeff=0;
    dsd.mutheta=(T)2;
    dsd.dtheta=(T)2.5;
    dsd.use_log_Lcj=true;
    dsd.nb_width=5;
}
//#####################################################################
// Function Initialize_Advection
//#####################################################################
void Initialize_Advection()    PHYSBAM_OVERRIDE
{
    fluids_parameters.Use_No_Fluid_Coupling_Defaults();
}
//#####################################################################
// Function Initialize_Velocities
//#####################################################################
void Initialize_Velocities() PHYSBAM_OVERRIDE
{
    GRID<TV>& grid=*fluids_parameters.grid;
    for(FACE_ITERATOR<TV> iterator(grid);iterator.Valid();iterator.Next()){
        const TV location=iterator.Location();
        if(source_sphere.Lazy_Inside(location)){
            const int axis=iterator.Axis();const TV_INT face_index=iterator.Face_Index();
            fluid_collection.incompressible_fluid_collection.face_velocities.Component(axis)(face_index)=source_sphere.Normal(location)[axis]*normal_velocity;}}

    DETONATION_SHOCK_DYNAMICS<TV>& dsd=*fluids_parameters.incompressible_multiphase->projection.dsd;
    fluids_parameters.incompressible_multiphase->projection.dsd->Dn.array.Fill(Dn_initial);
    dsd.Dn.boundary->Fill_Ghost_Cells(dsd.Dn.grid,dsd.Dn.array,dsd.Dn.array,0,0,3); // TODO: use real time/dt
    LEVELSET<TV>& levelset=*fluids_parameters.particle_levelset_evolution_multiple->particle_levelset_multiple.levelset_multiple.levelsets(1);
    levelset.Compute_Curvature();dsd.curvature_old.array=*levelset.curvature;
}
//#####################################################################
// Function Get_Source_Velocities
//#####################################################################
void Get_Source_Velocities(ARRAY<T,FACE_INDEX<3> >& face_velocities,ARRAY<bool,FACE_INDEX<3> >& psi_N,const T time) PHYSBAM_OVERRIDE
{
    GRID<TV>& grid=*fluids_parameters.grid;
    for(FACE_ITERATOR<TV> iterator(grid);iterator.Valid();iterator.Next()){
        const TV location=iterator.Location();
        if(source_sphere.Lazy_Inside(location)){
            const int axis=iterator.Axis();const TV_INT face_index=iterator.Face_Index();
            psi_N.Component(axis)(face_index)=true;
            face_velocities.Component(axis)(face_index)=source_sphere.Normal(location)[axis]*normal_velocity;}}
}
//#####################################################################
// Function Initialize_Phi
//#####################################################################
void Initialize_Phi() PHYSBAM_OVERRIDE
{
    GRID<TV>& grid=*fluids_parameters.grid;
    ARRAY<ARRAY<T,TV_INT>>& phis=fluids_parameters.particle_levelset_evolution_multiple->phis;
    for(CELL_ITERATOR<TV> iterator(grid);iterator.Valid();iterator.Next()) phis(1)(iterator.Cell_Index())=source_sphere.Signed_Distance(iterator.Location());
    phis(2).Copy(-1,phis(1));
}
//#####################################################################
// Function Adjust_Phi_With_Sources
//#####################################################################
bool Adjust_Phi_With_Sources(const T time) PHYSBAM_OVERRIDE
{
    GRID<TV>& grid=*fluids_parameters.grid;
    ARRAY<ARRAY<T,TV_INT>>& phis=fluids_parameters.particle_levelset_evolution_multiple->phis;
    for(CELL_ITERATOR<TV> iterator(grid);iterator.Valid();iterator.Next())
        if(source_sphere.Lazy_Inside(iterator.Location())) phis(1)(iterator.Cell_Index())=-grid.dX.x;
    return false;
}
//##################################################################### 
};      
#endif
