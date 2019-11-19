//#####################################################################
// Copyright 2007, Nipun Kwatra, Jonathan Su
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class DSD_NO_NAVIER_STOKES
//#####################################################################
#ifndef __DSD_NO_NAVIER_STOKES__
#define __DSD_NO_NAVIER_STOKES__

#include <Grid_Tools/Grids/CELL_ITERATOR.h>
#include <Geometry/Basic_Geometry/CYLINDER.h>
#include <Geometry/Basic_Geometry/SPHERE.h>
#include <Incompressible/Collisions_And_Interactions/GRID_BASED_COLLISION_GEOMETRY_UNIFORM.h>
#include <Dynamics/Incompressible_Flows/SPH_EVOLUTION_UNIFORM.h>
#include <Dynamics/Solids_And_Fluids/FLUIDS_PARAMETERS_UNIFORM.h>
#include <Dynamics/Solids_And_Fluids/SOLIDS_FLUIDS_EXAMPLE_UNIFORM.h>
#include <Dynamics/Standard_Tests/WATER_STANDARD_TESTS_3D.h>
namespace PhysBAM{

template<class T_input>
class DSD_NO_NAVIER_STOKES:public SOLIDS_FLUIDS_EXAMPLE_UNIFORM<VECTOR<T_input,3> >
{
    typedef T_input T;
public:
    typedef VECTOR<T,3> TV;typedef VECTOR<int,3> TV_INT;
    typedef ARRAY<PARTICLE_LEVELSET_REMOVED_PARTICLES<TV>*,TV_INT> T_ARRAYS_PARTICLE_LEVELSET_REMOVED_PARTICLES;

    typedef SOLIDS_FLUIDS_EXAMPLE_UNIFORM<TV> BASE;
    using BASE::last_frame;using BASE::frame_rate;using BASE::restart;using BASE::restart_frame;using BASE::output_directory;using BASE::Adjust_Phi_With_Sources;
    using BASE::Get_Source_Reseed_Mask;using BASE::Get_Source_Velocities;using BASE::fluids_parameters;using BASE::fluid_collection;using BASE::solids_parameters;using BASE::data_directory;
    using BASE::test_number;using BASE::resolution;using BASE::user_last_frame;

    ARRAY<SPHERE<TV> > sources;
    T source_end_time;
    T normal_velocity;
    T Dn_initial;

    DSD_NO_NAVIER_STOKES(const STREAM_TYPE stream_type_input,PARSE_ARGS& parse_args)
        :SOLIDS_FLUIDS_EXAMPLE_UNIFORM<TV>(stream_type_input,parse_args,1,fluids_parameters.FIRE),source_end_time(3.0),normal_velocity(4),Dn_initial((T)0.2)
    {
        parse_args.Parse();

        // set up the standard fluid environment
        if(!this->user_frame_rate) frame_rate=30;
        restart=false;restart_frame=0;
        fluids_parameters.domain_walls[0][0]=true;fluids_parameters.domain_walls[0][1]=true;fluids_parameters.domain_walls[1][0]=true;
        fluids_parameters.domain_walls[1][1]=false;fluids_parameters.domain_walls[2][0]=true;fluids_parameters.domain_walls[2][1]=true;
        fluids_parameters.number_particles_per_cell=16;
        fluids_parameters.viscosity=(T)0;fluids_parameters.implicit_viscosity=false;
        fluids_parameters.write_levelset=true;fluids_parameters.write_velocity=true;fluids_parameters.write_particles=true;fluids_parameters.write_debug_data=true;
        fluids_parameters.write_ghost_values=true;fluids_parameters.write_removed_positive_particles=true;fluids_parameters.write_removed_negative_particles=true;
        fluids_parameters.delete_fluid_inside_objects=true;
        fluids_parameters.incompressible_iterations=40;
        fluids_parameters.use_removed_positive_particles=true;fluids_parameters.use_removed_negative_particles=true;
        fluids_parameters.second_order_cut_cell_method=true;
        fluids_parameters.store_particle_ids=true;
        fluids_parameters.use_vorticity_confinement=false;
        fluids_parameters.use_vorticity_confinement_fuel=false;
        fluids_parameters.cfl=(T)1.9;
        fluids_parameters.gravity=TV();
        fluids_parameters.fuel_region(1)=true;

        //DSD parameters
        fluids_parameters.analytic_test=true;
        fluids_parameters.use_dsd=true;

        int cells=1*resolution;
        if(!user_last_frame) last_frame=10;
        GRID<TV>& grid=*fluids_parameters.grid;
        grid.Initialize(TV_INT(10*cells+1,10*cells+1,10*cells+1),RANGE<TV>(TV(0,0,0),TV(8,8,8)));

        if(!user_last_frame) last_frame=1000;
        
        if(!this->user_output_directory)
            output_directory=LOG::sprintf("DSD_No_Navier_Stokes/DSD_No_Navier_Stokes_%d__Resolution_%d_%d_%d",test_number,(grid.counts.x-1),(grid.counts.y-1),(grid.counts.z-1));
        LOG::cout<<"Running DSD simulation to "<<output_directory<<std::endl;

        //sources
        sources.Append(SPHERE<TV>(TV(3.75,4,4),(T)0.2));
        sources.Append(SPHERE<TV>(TV(4.25,4,4),(T)0.2));
        sources.Append(SPHERE<TV>(TV(4,3.75,4),(T)0.2));
        sources.Append(SPHERE<TV>(TV(4,4,3.75),(T)0.2));
        sources.Append(SPHERE<TV>(TV(4,4,4.25),(T)0.2));
        sources.Append(SPHERE<TV>(TV(4,3.75,4),(T)0.2));
        sources.Append(SPHERE<TV>(TV(4,4,4),(T)0.2));
        //sources.Append(SPHERE<TV>(TV(4,4.25,4),(T)0.2));
    }

    ~DSD_NO_NAVIER_STOKES()
    {}

//#####################################################################
// Function Initialize_Advection
//#####################################################################
void Initialize_Advection() override
{
    fluids_parameters.Use_No_Fluid_Coupling_Defaults();
}
//#####################################################################
// Function Initialize_Velocities
//#####################################################################
void Initialize_Velocities() override
{
    GRID<TV>& grid=*fluids_parameters.grid;
    for(FACE_ITERATOR<TV> iterator(grid);iterator.Valid();iterator.Next()){
        const TV location=iterator.Location();
        for(int s=0;s<sources.m;s++) if(sources(s).Lazy_Inside(location)){
            const int axis=iterator.Axis();const TV_INT face_index=iterator.Face_Index();
            fluid_collection.incompressible_fluid_collection.face_velocities.Component(axis)(face_index)=sources(s).Normal(location)[axis]*normal_velocity;}}

    DETONATION_SHOCK_DYNAMICS<TV>& dsd=*fluids_parameters.incompressible->projection.dsd;
    fluids_parameters.incompressible->projection.dsd->Dn.array.Fill(Dn_initial);
    dsd.Dn.boundary->Fill_Ghost_Cells(dsd.Dn.grid,dsd.Dn.array,dsd.Dn.array,0,0,3); // TODO: use real dt/time
    LEVELSET<TV>& levelset=fluids_parameters.particle_levelset_evolution->Particle_Levelset(0).levelset;
    levelset.Compute_Curvature();dsd.curvature_old.array=*levelset.curvature;
}
//#####################################################################
// Function Initialize_Phi
//#####################################################################
void Initialize_Phi() override
{
    PHYSBAM_DEBUG_WRITE_SUBSTEP("After Initialize_Phi",0);
    GRID<TV>& grid=*fluids_parameters.grid;
    ARRAY<T,VECTOR<int,3> >& phi=fluids_parameters.particle_levelset_evolution->phi;
    T initial_phi;
    for(CELL_ITERATOR<TV> iterator(grid);iterator.Valid();iterator.Next()){
        TV X=iterator.Location();
        initial_phi=1;
        for(int s=0;s<sources.m;s++) initial_phi=min(initial_phi,sources(s).Signed_Distance(X));
        phi(iterator.Cell_Index())=initial_phi;}
    PHYSBAM_DEBUG_WRITE_SUBSTEP("After Initialize_Phi",0);
}
//#####################################################################
// Function Update_Fluid_Parameters
//#####################################################################
void Update_Fluid_Parameters(const T dt,const T time) override
{
    SOLIDS_FLUIDS_EXAMPLE_UNIFORM<TV>::Update_Fluid_Parameters(dt,time);
    DETONATION_SHOCK_DYNAMICS<TV>& dsd=*fluids_parameters.incompressible->projection.dsd;
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
// Function Get_Analytic_Velocities
//#####################################################################
void Get_Analytic_Velocities(const T time) const override
{
    PHYSBAM_FATAL_ERROR("broken");
#if 0
    GRID<TV>& grid=*fluids_parameters.grid;
    DETONATION_SHOCK_DYNAMICS<TV>& dsd=*fluids_parameters.incompressible->projection.dsd;
    dsd.Dn.boundary->Fill_Ghost_Cells(dsd.Dn.grid,dsd.Dn.array,dsd.Dn.array,0,time);
    LEVELSET<TV>& levelset=fluids_parameters.particle_levelset_evolution->Particle_Levelset(0).levelset;
    ARRAY<T,FACE_INDEX<TV::m> >& face_velocities=fluid_collection.incompressible_fluid_collection.face_velocities;
    for(FACE_ITERATOR<TV> iterator(grid);iterator.Valid();iterator.Next()){
        int axis=iterator.Axis();TV_INT face=iterator.Face_Index();
        TV_INT cell1,cell2;grid.Cells_Touching_Face(axis,face,cell1,cell2);
        face_velocities.Component(axis)(face)=dsd.Normal_Flame_Speed(axis,face)*(levelset.phi(cell2)-levelset.phi(cell1))*grid.One_Over_DX()[axis];}
#endif
}
//#####################################################################
};
}
#endif
