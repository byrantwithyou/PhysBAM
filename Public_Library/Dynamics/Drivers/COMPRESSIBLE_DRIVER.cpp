//#####################################################################
// Copyright 2009, Nipun Kwatra.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Core/Log/LOG.h>
#include <Core/Log/SCOPE.h>
#include <Tools/Ordinary_Differential_Equations/RUNGEKUTTA.h>
#include <Grid_Tools/Grids/CELL_ITERATOR.h>
#include <Grid_Tools/Grids/GRID.h>
#include <Rigids/Collisions/COLLISION_GEOMETRY.h>
#include <Incompressible/Collisions_And_Interactions/GRID_BASED_COLLISION_GEOMETRY_UNIFORM.h>
#include <Incompressible/Grids_Uniform_PDE_Linear/POISSON_COLLIDABLE_UNIFORM.h>
#include <Compressible/Euler_Equations/EULER_LAPLACE.h>
#include <Compressible/Euler_Equations/EULER_UNIFORM.h>
#include <Dynamics/Coupled_Evolution/SOLID_COMPRESSIBLE_FLUID_COUPLING_UTILITIES.h>
#include <Dynamics/Drivers/COMPRESSIBLE_DRIVER.h>
#include <Dynamics/Drivers/COMPRESSIBLE_EXAMPLE.h>
using namespace PhysBAM;
//#####################################################################
// Initialize
//#####################################################################
template<class TV> COMPRESSIBLE_DRIVER<TV>::
COMPRESSIBLE_DRIVER(COMPRESSIBLE_EXAMPLE<TV>& example)
    :DRIVER<TV>(example),example(example)
{}
//#####################################################################
// Initialize
//#####################################################################
template<class TV> COMPRESSIBLE_DRIVER<TV>::
~COMPRESSIBLE_DRIVER()
{}
//#####################################################################
// Initialize
//#####################################################################
template<class TV> void COMPRESSIBLE_DRIVER<TV>::
Initialize()
{
    BASE::Initialize();

    EULER_UNIFORM<TV>& euler=example.euler;
    SOLID_COMPRESSIBLE_FLUID_COUPLING_UTILITIES<TV>& euler_solid_fluid_coupling_utilities=example.euler_solid_fluid_coupling_utilities;

    // setup grids
    euler.Initialize_Domain(example.mac_grid);
    example.compressible_fluid_collection.Initialize_Grids();
    example.compressible_fluid_collection.Set_Equation_Of_State(new EOS_GAMMA<T>());
    if(example.conservation_method) euler.Set_Custom_Conservation(*example.conservation_method);
    if(example.boundary) euler.Set_Custom_Boundary(*example.boundary);
    if(example.pressure_boundary) euler.euler_projection.Set_Custom_Pressure_Boundary(*example.pressure_boundary);
    if(example.set_max_time_step) euler.Set_Max_Time_Step(example.max_time_step);
    euler.conservation->Set_Order(example.spatial_order);
    euler.timesplit=example.timesplit;
    euler_solid_fluid_coupling_utilities.fluid_affects_solid=example.fluid_affects_solid;
    euler.perform_rungekutta_for_implicit_part=example.perform_rungekutta_for_implicit_part;
    euler.use_sound_speed_for_cfl=example.use_sound_speed_for_cfl;
    euler.use_sound_speed_based_dt_multiple_for_cfl=example.use_sound_speed_based_dt_multiple_for_cfl;
    euler.multiplication_factor_for_sound_speed_based_dt=example.multiplication_factor_for_sound_speed_based_dt;

    // setup elliptic solver
    euler.euler_projection.elliptic_solver->solve_single_cell_neumann_regions=example.solve_single_cell_neumann_regions;
    euler.euler_projection.elliptic_solver->Set_Relative_Tolerance(example.tolerance);
    euler.euler_projection.elliptic_solver->pcg.Set_Maximum_Iterations(example.iterations);
    euler.euler_projection.elliptic_solver->pcg.Show_Results();
    euler.euler_projection.elliptic_solver->pcg.Show_Residuals();

    // initialize bodies
    example.Initialize_Bodies();
    example.Initialize_Solid_Fluid_Coupling();

    // initialize state variables
    example.Initialize_Euler_State();

    example.After_Initialization();

    // Boundary conditions
    example.Set_Boundary_Conditions(time); // get so CFL is correct
    euler_solid_fluid_coupling_utilities.Fill_Solid_Cells();

    // Logging
    if(example.monitor_conservation_error) euler.Compute_Total_Conserved_Quantity(false,(T)0,euler.initial_total_conserved_quantity);

    Write_Output_Files(example.first_frame);
}
//#####################################################################
// Advance_To_Target_Time
//#####################################################################
template<class TV> void COMPRESSIBLE_DRIVER<TV>::
Advance_To_Target_Time(const T target_time)
{
    bool done=false;for(int substep=1;!done;substep++){
        LOG::SCOPE scope("SUBSTEP","substep %d",substep);

        example.collision_bodies_affecting_fluid->Save_State(COLLISION_GEOMETRY<TV>::FLUID_COLLISION_GEOMETRY_OLD_STATE,time);
        Setup_Fluids(time);
        T dt=Compute_Dt(time,target_time,done);

        example.Advance_Kinematic_Collision_Bodies(dt,time);
        Restore_Solids_To_Time_N(time+dt);

        PHYSBAM_DEBUG_WRITE_SUBSTEP("object compatibility",1);
        if(example.solid_affects_fluid) example.euler_solid_fluid_coupling_utilities.Fill_Solid_Cells();
        PHYSBAM_DEBUG_WRITE_SUBSTEP("advect fluid",1);
        Advect_Fluid(dt,substep);

        Restore_Solids_To_Time_N_Plus_One();

        PHYSBAM_DEBUG_WRITE_SUBSTEP("project fluid at end of substep",1);
        Advance_Fluid_One_Time_Step_Implicit_Part(dt,time,substep);

        time+=dt;

        PHYSBAM_DEBUG_WRITE_SUBSTEP("END Substep %d",0,substep);}
}
//#####################################################################
// Setup_Fluids
//#####################################################################
template<class TV> void COMPRESSIBLE_DRIVER<TV>::
Setup_Fluids(const T time)
{
    example.euler.Clamp_Internal_Energy((T)1./example.frame_rate,time); // fictitious dt
}
//#####################################################################
// Restore_Solids_To_Time_N
//#####################################################################
template<class TV> void COMPRESSIBLE_DRIVER<TV>::
Restore_Solids_To_Time_N(const T time_n_plus_one)
{
    GRID_BASED_COLLISION_GEOMETRY_UNIFORM<TV>& collision_bodies_affecting_fluid=*example.collision_bodies_affecting_fluid;

    collision_bodies_affecting_fluid.Save_State(COLLISION_GEOMETRY<TV>::FLUID_COLLISION_GEOMETRY_NEW_STATE,time_n_plus_one);
    collision_bodies_affecting_fluid.Restore_State(COLLISION_GEOMETRY<TV>::FLUID_COLLISION_GEOMETRY_OLD_STATE);
    collision_bodies_affecting_fluid.Update_Intersection_Acceleration_Structures(true,COLLISION_GEOMETRY<TV>::FLUID_COLLISION_GEOMETRY_NEW_STATE,
        COLLISION_GEOMETRY<TV>::FLUID_COLLISION_GEOMETRY_OLD_STATE);
    collision_bodies_affecting_fluid.Rasterize_Objects(); // non-swept
}
//#####################################################################
// Restore_Solids_To_Time_N_Plus_One
//#####################################################################
template<class TV> void COMPRESSIBLE_DRIVER<TV>::
Restore_Solids_To_Time_N_Plus_One()
{
    GRID_BASED_COLLISION_GEOMETRY_UNIFORM<TV>& collision_bodies_affecting_fluid=*example.collision_bodies_affecting_fluid;

    collision_bodies_affecting_fluid.Restore_State(COLLISION_GEOMETRY<TV>::FLUID_COLLISION_GEOMETRY_NEW_STATE);
    collision_bodies_affecting_fluid.Update_Intersection_Acceleration_Structures(false); // NON-swept acceleration structures
    collision_bodies_affecting_fluid.Rasterize_Objects(); // non-swept
}
//#####################################################################
// Advect_Fluid
//#####################################################################
template<class TV> void COMPRESSIBLE_DRIVER<TV>::
Advect_Fluid(const T dt,const int substep)
{
    EULER_UNIFORM<TV>& euler=example.euler;
    SOLID_COMPRESSIBLE_FLUID_COUPLING_UTILITIES<TV>& euler_solid_fluid_coupling_utilities=example.euler_solid_fluid_coupling_utilities;
    COMPRESSIBLE_FLUID_COLLECTION<TV>& compressible_fluid_collection=example.compressible_fluid_collection;

    PHYSBAM_DEBUG_WRITE_SUBSTEP("before compressible explicit solve",1);
    if(euler_solid_fluid_coupling_utilities.thinshell){
        euler_solid_fluid_coupling_utilities.uncovered_cells.Fill(false);
        for(CELL_ITERATOR<TV> iterator(euler.grid);iterator.Valid();iterator.Next())
            if(example.collision_bodies_affecting_fluid->Any_Simplex_Crossover(iterator.Location(),iterator.Location(),dt)) euler_solid_fluid_coupling_utilities.uncovered_cells(iterator.Cell_Index())=true;}

    euler_solid_fluid_coupling_utilities.Update_Cut_Out_Grid();
    example.Set_Boundary_Conditions(time); // Should be computed at time 'n' for explicit part
    LOG::Time("compressible explicit update");

    // initialize p_advected to current eos pressure
    if(euler.timesplit && !euler.perform_rungekutta_for_implicit_part){
        for(CELL_ITERATOR<TV> iterator(euler.grid);iterator.Valid();iterator.Next()){TV_INT cell_index=iterator.Cell_Index();
            euler.euler_projection.p_advected(cell_index)=euler.eos->p(compressible_fluid_collection.U(cell_index)(1),EULER<TV>::e(compressible_fluid_collection.U(cell_index)));}}

    RUNGEKUTTA<T_ARRAYS_DIMENSION_SCALAR> rungekutta_u(compressible_fluid_collection.U,example.rungekutta_order,dt,time);
    RUNGEKUTTA<ARRAY<T,TV_INT>> rungekutta_p_advected(euler.euler_projection.p_advected,example.rungekutta_order,dt,time);
    for(int rk_substep=0;rk_substep<rungekutta_u.order;rk_substep++){
        euler.Advance_One_Time_Step_Explicit_Part(dt,rungekutta_u.time,rk_substep,rungekutta_u.order);
        if(euler.timesplit && euler.perform_rungekutta_for_implicit_part){
            euler.Get_Dirichlet_Boundary_Conditions(dt,time);
            example.Set_Boundary_Conditions(rungekutta_u.time+dt);
            euler.Advance_One_Time_Step_Implicit_Part(dt,rungekutta_u.time);}
        if(!euler.timesplit || euler.perform_rungekutta_for_implicit_part) example.Apply_Isobaric_Fix(dt,rungekutta_u.time);
        euler.Remove_Added_Internal_Energy(dt,rungekutta_u.time);
        rungekutta_u.Next();
        euler.boundary->Apply_Boundary_Condition(euler.grid,compressible_fluid_collection.U,rungekutta_u.time);
        if(euler.timesplit && !euler.perform_rungekutta_for_implicit_part) rungekutta_p_advected.Next();
        if(rk_substep!=rungekutta_u.order) euler.Clamp_Internal_Energy(dt,rungekutta_u.time);}
    if(euler.timesplit && !euler.perform_rungekutta_for_implicit_part) euler.Get_Dirichlet_Boundary_Conditions(dt,time);
    LOG::Stop_Time();
    PHYSBAM_DEBUG_WRITE_SUBSTEP("after compressible explicit solve",1);
}
//#####################################################################
// Advance_Fluid_One_Time_Step_Implicit_Part
//#####################################################################
template<class TV> void COMPRESSIBLE_DRIVER<TV>::
Advance_Fluid_One_Time_Step_Implicit_Part(const T dt_projection,const T time_projection,const int substep)
{
    EULER_UNIFORM<TV>& euler=example.euler;

    if(euler.timesplit && !euler.perform_rungekutta_for_implicit_part){
        LOG::Time("compressible implicit update");
        euler.Clamp_Internal_Energy(dt_projection,time_projection+dt_projection);
        example.Set_Boundary_Conditions(time_projection+dt_projection);
        euler.Advance_One_Time_Step_Implicit_Part(dt_projection,time_projection);
        example.Apply_Isobaric_Fix(dt_projection,time_projection);
        euler.Remove_Added_Internal_Energy(dt_projection,time_projection);
        PHYSBAM_DEBUG_WRITE_SUBSTEP("after compressible implicit solve",1);
        LOG::Stop_Time();}
}
//#####################################################################
// Function Compute_Dt
//#####################################################################
template<class TV> typename TV::SCALAR COMPRESSIBLE_DRIVER<TV>::
Compute_Dt(const T time,const T target_time,bool& done)
{
    T dt=example.cfl*example.euler.CFL(time);
    example.Limit_Dt(dt,time);
    done=false;
    example.Clamp_Time_Step_With_Target_Time(time,target_time,dt,done,0);
    return dt;
}
//#####################################################################
namespace PhysBAM{
template class COMPRESSIBLE_DRIVER<VECTOR<float,1> >;
template class COMPRESSIBLE_DRIVER<VECTOR<float,2> >;
template class COMPRESSIBLE_DRIVER<VECTOR<float,3> >;
template class COMPRESSIBLE_DRIVER<VECTOR<double,1> >;
template class COMPRESSIBLE_DRIVER<VECTOR<double,2> >;
template class COMPRESSIBLE_DRIVER<VECTOR<double,3> >;
}
