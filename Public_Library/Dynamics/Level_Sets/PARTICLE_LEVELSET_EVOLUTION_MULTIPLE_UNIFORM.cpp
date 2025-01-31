//#####################################################################
// Copyright 2005, Ron Fedkiw, Eran Guendelman, Frank Losasso, Tamar Shinar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Core/Log/LOG.h>
#include <Tools/Ordinary_Differential_Equations/RUNGEKUTTA.h>
#include <Grid_Tools/Grids/CELL_ITERATOR.h>
#include <Grid_PDE/Boundaries/BOUNDARY.h>
#include <Geometry/Level_Sets/LEVELSET.h>
#include <Incompressible/Level_Sets/LEVELSET_MULTIPLE.h>
#include <Dynamics/Level_Sets/LEVELSET_ADVECTION_MULTIPLE.h>
#include <Dynamics/Level_Sets/LEVELSET_CALLBACKS.h>
#include <Dynamics/Level_Sets/PARTICLE_LEVELSET_EVOLUTION_MULTIPLE_UNIFORM.h>
#include <Dynamics/Particles/PARTICLE_LEVELSET_REMOVED_PARTICLES.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class TV> PARTICLE_LEVELSET_EVOLUTION_MULTIPLE_UNIFORM<TV>::
PARTICLE_LEVELSET_EVOLUTION_MULTIPLE_UNIFORM(const GRID<TV>& grid_input,GRID_BASED_COLLISION_GEOMETRY_UNIFORM<TV>& collision_body_list_input,const int number_of_ghost_cells_input)
    :PARTICLE_LEVELSET_EVOLUTION_UNIFORM<TV>(grid_input,collision_body_list_input,number_of_ghost_cells_input,true),
    particle_levelset_multiple(*new PARTICLE_LEVELSET_MULTIPLE_UNIFORM<TV>(grid,phis,number_of_ghost_cells_input)),
    levelset_advection_multiple(*new LEVELSET_ADVECTION_MULTIPLE<TV>(particle_levelset_multiple.levelset_multiple))
{
    Use_Semi_Lagrangian_Advection();
    Track_Mass(false);
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> PARTICLE_LEVELSET_EVOLUTION_MULTIPLE_UNIFORM<TV>::
~PARTICLE_LEVELSET_EVOLUTION_MULTIPLE_UNIFORM()
{
    delete &particle_levelset_multiple;
    delete &levelset_advection_multiple;
}
//#####################################################################
// Function Time_Step
//#####################################################################
template<class TV> typename TV::SCALAR PARTICLE_LEVELSET_EVOLUTION_MULTIPLE_UNIFORM<TV>::
Time_Step(const T stopping_time,bool& limited_by_stopping_time)
{
    T dt=CFL();limited_by_stopping_time=false;
    if(time+dt>=stopping_time){dt=stopping_time-time;limited_by_stopping_time=true;}
    else if(time+2*dt>=stopping_time) dt=(T).51*(stopping_time-time);
    return dt;
}
//#####################################################################
// Function CFL
//#####################################################################
template<class TV> typename TV::SCALAR PARTICLE_LEVELSET_EVOLUTION_MULTIPLE_UNIFORM<TV>::
CFL(const bool need_to_get_velocity,const bool analytic_test)
{
    if(need_to_get_velocity) particle_levelset_multiple.levelset_multiple.levelset_callbacks->Get_Levelset_Velocity(grid,particle_levelset_multiple.levelset_multiple,V,time);
    if(analytic_test){T max_time_step=particle_levelset_multiple.levelset_multiple.levelsets(0)->max_time_step;
        for(int i=1;i<particle_levelset_multiple.levelset_multiple.levelsets.m;i++) max_time_step=min(max_time_step,particle_levelset_multiple.levelset_multiple.levelsets(i)->max_time_step);
        return cfl_number/max(V.Max_Abs().Max(),1/max_time_step);}
    return cfl_number*particle_levelset_multiple.levelset_multiple.CFL(V);
}
//#####################################################################
// Function Advance_To_Time
//#####################################################################
template<class TV> void PARTICLE_LEVELSET_EVOLUTION_MULTIPLE_UNIFORM<TV>::
Advance_To_Time(ARRAY<T,FACE_INDEX<TV::m> >* face_velocities,const T stopping_time,const bool verbose)
{
    int substep=0;bool done=false;
    while(!done){substep++;
        T dt=Time_Step(stopping_time,done);
        if(verbose) LOG::cout<<"substep = "<<substep<<", dt = "<<dt<<std::endl;
        Advance_One_Time_Step(face_velocities,dt);}
    if(track_mass)for(int i=0;i<particle_levelset_multiple.particle_levelsets.m;i++){
        T mass=Levelset_Advection(i).Approximate_Negative_Material();
        LOG::cout<<"negative material("<<i<<") = "<<mass<<" - change = "<<(mass-initial_mass(i))/initial_mass(i)*100<<"%"<<std::endl;}
}
//#####################################################################
// Function Advance_One_Time_Step
//#####################################################################
template<class TV> void PARTICLE_LEVELSET_EVOLUTION_MULTIPLE_UNIFORM<TV>::
Advance_One_Time_Step(ARRAY<T,FACE_INDEX<TV::m> >* face_velocities,const T dt)
{
    LOG::Time("advancing levelset");
    Advance_Levelset(dt);
    LOG::Time("advancing particles");
    Advance_Particles(*face_velocities,dt);
    Modify_Levelset_And_Particles(face_velocities);
}
//#####################################################################
// Function Advance_Levelset
//#####################################################################
template<class TV> void PARTICLE_LEVELSET_EVOLUTION_MULTIPLE_UNIFORM<TV>::
Advance_Levelset(const T dt)
{
    ARRAY<RUNGEKUTTA<ARRAY<T,TV_INT> >*> rungekutta_phis(phis.m);
    for(int i=0;i<phis.m;i++) rungekutta_phis(i)=new RUNGEKUTTA<ARRAY<T,TV_INT> >(phis(i),runge_kutta_order_levelset,dt,time);
    for(int k=0;k<runge_kutta_order_levelset;k++){
        if(k==0 || !use_frozen_velocity) particle_levelset_multiple.levelset_multiple.levelset_callbacks->Get_Levelset_Velocity(grid,particle_levelset_multiple.levelset_multiple,V,time);
        levelset_advection_multiple.Euler_Step(V,dt,time,particle_levelset_multiple.number_of_ghost_cells);
        for(int i=0;i<rungekutta_phis.m;i++){
            rungekutta_phis(i)->Next();
            time=rungekutta_phis(i)->time;
            particle_levelset_multiple.particle_levelsets(i)->levelset.boundary->Apply_Boundary_Condition(grid,phis(i),time);}}
    rungekutta_phis.Delete_Pointers_And_Clean_Memory();
}
//#####################################################################
// Function Advance_Particles
//#####################################################################
template<class TV> void PARTICLE_LEVELSET_EVOLUTION_MULTIPLE_UNIFORM<TV>::
Advance_Particles(const ARRAY<T,FACE_INDEX<TV::m> >& face_velocities,const T dt,const bool analytic_test)
{
    if(analytic_test) PHYSBAM_NOT_IMPLEMENTED("analytic_test");
    if(use_particle_levelset){
        time-=dt; // to fix up time advancement due to Advance_Levelset()
        if(runge_kutta_order_particles==1 || (runge_kutta_order_particles==2 && use_frozen_velocity)){
            // TODO: still needed?
            //particle_levelset_multiple.levelset_multiple.levelset_callbacks->Get_Levelset_Velocity(grid,particle_levelset_multiple.levelset_multiple,V,time);
            particle_levelset_multiple.Euler_Step_Particles(face_velocities,dt,time,runge_kutta_order_particles==2);time+=dt;}
        else if(runge_kutta_order_particles==2 || runge_kutta_order_particles==3){
            T start_time=time;
            for(int i=0;i<particle_levelset_multiple.particle_levelsets.m;i++)
                time=Advance_Particles(particle_levelset_multiple.particle_levelsets(i)->negative_particles,PARTICLE_LEVELSET_NEGATIVE,dt,start_time);
            particle_levelset_multiple.Euler_Step_Removed_Particles(dt,start_time,true);}} // can only Euler Step removed particles
}
//#####################################################################
// Function Advance_Particles
//#####################################################################
template<class TV> typename TV::SCALAR PARTICLE_LEVELSET_EVOLUTION_MULTIPLE_UNIFORM<TV>::
Advance_Particles(ARRAY<PARTICLE_LEVELSET_PARTICLES<TV>*,TV_INT>& particles,const PARTICLE_LEVELSET_PARTICLE_TYPE particle_type,const T dt,const T input_time)
{
    PHYSBAM_NOT_IMPLEMENTED();
}
//#####################################################################
// Function Advance_Particles
//#####################################################################
template<class TV> typename TV::SCALAR PARTICLE_LEVELSET_EVOLUTION_MULTIPLE_UNIFORM<TV>::
Advance_Particles(ARRAY<PARTICLE_LEVELSET_REMOVED_PARTICLES<TV>*,TV_INT>& particles,const PARTICLE_LEVELSET_PARTICLE_TYPE particle_type,const T dt,const T input_time)
{
    PHYSBAM_NOT_IMPLEMENTED();
}
//#####################################################################
// Function Modify_Levelset_And_Particles
//#####################################################################
template<class TV> void PARTICLE_LEVELSET_EVOLUTION_MULTIPLE_UNIFORM<TV>::
Modify_Levelset_And_Particles(ARRAY<T,FACE_INDEX<TV::m> >* face_velocities)
{
    if(use_particle_levelset){
        LOG::Time("modifying levelset");
        particle_levelset_multiple.Modify_Levelset_Using_Escaped_Particles(face_velocities);}
    LOG::Time("reinitializing and projecting levelset");
    particle_levelset_multiple.levelset_multiple.Project_Levelset();
    Make_Signed_Distance();
    if(use_particle_levelset){
        LOG::Time("modifying levelset");
        particle_levelset_multiple.Modify_Levelset_Using_Escaped_Particles(face_velocities);
        LOG::Time("adjusting particle radii");
        particle_levelset_multiple.Adjust_Particle_Radii();
        LOG::Stop_Time();}
    particle_levelset_multiple.levelset_multiple.Project_Levelset();
}
//#####################################################################
// Function Reseed_Particles
//#####################################################################
template<class TV> void PARTICLE_LEVELSET_EVOLUTION_MULTIPLE_UNIFORM<TV>::
Reseed_Particles(const T time,const int time_step,ARRAY<bool,TV_INT>* cell_centered_mask,const bool verbose)
{
    if((use_particle_levelset && cell_centered_mask) || (!time_step || (reseeding_frequency && time_step%reseeding_frequency==0))){
        int new_particles=particle_levelset_multiple.Reseed_Particles(time,cell_centered_mask);
        if(verbose) LOG::cout<<"Reseeding... "<<new_particles<<" new particles"<<std::endl;} // need to reset based on new number of particles
}
//#####################################################################
// Function Fill_Levelset_Ghost_Cells
//#####################################################################
template<class TV> void PARTICLE_LEVELSET_EVOLUTION_MULTIPLE_UNIFORM<TV>::
Fill_Levelset_Ghost_Cells(const T time)
{
    for(int i=0;i<phis.m;i++) particle_levelset_multiple.levelset_multiple.levelsets(i)->boundary->Fill_Ghost_Cells(grid,phis(i),phis(i),0,time,particle_levelset_multiple.number_of_ghost_cells);
}
template<class TV> void PARTICLE_LEVELSET_EVOLUTION_MULTIPLE_UNIFORM<TV>::
Use_Semi_Lagrangian_Advection()
{
    for(int i=0;i<phis.m;i++) levelset_advection_multiple.levelset_advections(i).Use_Local_Semi_Lagrangian_Advection();
}
template<class TV> void PARTICLE_LEVELSET_EVOLUTION_MULTIPLE_UNIFORM<TV>::
Use_Hamilton_Jacobi_Weno_Advection()
{
    for(int i=0;i<phis.m;i++) levelset_advection_multiple.levelset_advections(i).Use_Local_WENO_For_Advection();
}
template<class TV> void PARTICLE_LEVELSET_EVOLUTION_MULTIPLE_UNIFORM<TV>::
Use_Hamilton_Jacobi_Eno_Advection(const int order)
{
    assert(order>=1 && order<=3);
    for(int i=0;i<phis.m;i++) levelset_advection_multiple.levelset_advections(i).Use_Local_ENO_For_Advection(order);
}
template<class TV> void PARTICLE_LEVELSET_EVOLUTION_MULTIPLE_UNIFORM<TV>::
Track_Mass(const bool track_mass_input)
{
    track_mass=track_mass_input;
    if(track_mass) for(int i=0;i<phis.m;i++){
        initial_mass(i)=levelset_advection_multiple.levelset_advections(i).Approximate_Negative_Material();
        LOG::cout<<"negative material("<<i<<") = "<<initial_mass(i)<<std::endl;}
}
template<class TV> void PARTICLE_LEVELSET_EVOLUTION_MULTIPLE_UNIFORM<TV>::
Initialize_Domain(const GRID<TV>& grid_input,GRID_BASED_COLLISION_GEOMETRY_UNIFORM<TV>& collision_body_list_input,const int number_of_regions,const bool use_only_negative_particles)  // don't call up to the base class here because we don't need those variables initialized OVERRIDE PROBLEM
{
    assert(grid_input.Is_MAC_Grid());
    grid=grid_input;
    phis.Resize(number_of_regions);
    for(int i=0;i<phis.m;i++) phis(i).Resize(grid.Domain_Indices(particle_levelset_multiple.number_of_ghost_cells));
    V.Resize(grid);
    particle_levelset_multiple.Initialize_Particle_Levelsets_And_Grid_Values(grid,phis,collision_body_list_input,number_of_regions,use_only_negative_particles);
    levelset_advection_multiple.Initialize();
    for(int i=0;i<phis.m;i++)
        if(levelset_advection_multiple.levelset_advections(i).semi_lagrangian_collidable)
            particle_levelset_multiple.particle_levelsets(i)->levelset.Initialize_Valid_Masks(grid);
    initial_mass.Resize(number_of_regions);
}
template<class TV> void PARTICLE_LEVELSET_EVOLUTION_MULTIPLE_UNIFORM<TV>::
Initialize_Domain(const GRID<TV>& grid_input)
{
    PHYSBAM_FATAL_ERROR();
}
template<class TV> void PARTICLE_LEVELSET_EVOLUTION_MULTIPLE_UNIFORM<TV>::
Make_Signed_Distance()
{
    if(use_fmm){
        ARRAY<int> local_advection_spatial_orders(phis.m);
        for(int i=0;i<phis.m;i++)
            local_advection_spatial_orders(i)=levelset_advection_multiple.levelset_advections(i).local_advection_spatial_order;
        particle_levelset_multiple.levelset_multiple.Fast_Marching_Method(local_advection_spatial_orders);}
    else if(use_reinitialization) levelset_advection_multiple.Reinitialize();
}
template<class TV> void PARTICLE_LEVELSET_EVOLUTION_MULTIPLE_UNIFORM<TV>::
Set_Number_Particles_Per_Cell(const int number_particles_per_cell)
{
    for(int i=0;i<phis.m;i++) particle_levelset_multiple.particle_levelsets(i)->Set_Number_Particles_Per_Cell(number_particles_per_cell);
}
template<class TV> void PARTICLE_LEVELSET_EVOLUTION_MULTIPLE_UNIFORM<TV>::
Set_Levelset_Callbacks(LEVELSET_CALLBACKS<TV>& levelset_callbacks)
{
    particle_levelset_multiple.levelset_multiple.Set_Levelset_Callbacks(levelset_callbacks);
}
template<class TV> void PARTICLE_LEVELSET_EVOLUTION_MULTIPLE_UNIFORM<TV>::
Bias_Towards_Negative_Particles(const bool bias_towards_negative_particles)
{
    for(int i=0;i<phis.m;i++) particle_levelset_multiple.particle_levelsets(i)->Bias_Towards_Negative_Particles(bias_towards_negative_particles);
}
template<class TV> void PARTICLE_LEVELSET_EVOLUTION_MULTIPLE_UNIFORM<TV>::
Set_Seed(const int seed)
{
    for(int i=0;i<phis.m;i++) particle_levelset_multiple.particle_levelsets(i)->random.Set_Seed(seed);
}
template<class TV> void PARTICLE_LEVELSET_EVOLUTION_MULTIPLE_UNIFORM<TV>::
Seed_Particles(const T time)
{
    for(int i=0;i<phis.m;i++) particle_levelset_multiple.particle_levelsets(i)->Seed_Particles(time);
}
template<class TV> void PARTICLE_LEVELSET_EVOLUTION_MULTIPLE_UNIFORM<TV>::
Delete_Particles_Outside_Grid()
{
    for(int i=0;i<phis.m;i++) particle_levelset_multiple.particle_levelsets(i)->Delete_Particles_Outside_Grid();
}
template<class TV> void PARTICLE_LEVELSET_EVOLUTION_MULTIPLE_UNIFORM<TV>::
Set_CFL_Number(const T cfl_number_input)
{
    PARTICLE_LEVELSET_EVOLUTION<T>::Set_CFL_Number(cfl_number_input);
    for(int i=0;i<phis.m;i++) particle_levelset_multiple.particle_levelsets(i)->cfl_number=cfl_number_input;
}
template<class TV> PARTICLE_LEVELSET_MULTIPLE_UNIFORM<TV>& PARTICLE_LEVELSET_EVOLUTION_MULTIPLE_UNIFORM<TV>::
Particle_Levelset_Multiple()
{
    return particle_levelset_multiple;
}
template<class TV> LEVELSET_MULTIPLE<TV>& PARTICLE_LEVELSET_EVOLUTION_MULTIPLE_UNIFORM<TV>::
Levelset_Multiple()
{
    return particle_levelset_multiple.levelset_multiple;
}
template<class TV> PARTICLE_LEVELSET_UNIFORM<TV>& PARTICLE_LEVELSET_EVOLUTION_MULTIPLE_UNIFORM<TV>::
Particle_Levelset(const int i)
{
    return *particle_levelset_multiple.particle_levelsets(i);
}
template<class TV> LEVELSET<TV>& PARTICLE_LEVELSET_EVOLUTION_MULTIPLE_UNIFORM<TV>::
Levelset(const int i)
{
    return *particle_levelset_multiple.levelset_multiple.levelsets(i);
}
template<class TV> LEVELSET_ADVECTION<TV>& PARTICLE_LEVELSET_EVOLUTION_MULTIPLE_UNIFORM<TV>::
Levelset_Advection(const int i)
{
    return levelset_advection_multiple.levelset_advections(i);
}
//#####################################################################
namespace PhysBAM{
template class PARTICLE_LEVELSET_EVOLUTION_MULTIPLE_UNIFORM<VECTOR<float,1> >;
template class PARTICLE_LEVELSET_EVOLUTION_MULTIPLE_UNIFORM<VECTOR<float,2> >;
template class PARTICLE_LEVELSET_EVOLUTION_MULTIPLE_UNIFORM<VECTOR<float,3> >;
template class PARTICLE_LEVELSET_EVOLUTION_MULTIPLE_UNIFORM<VECTOR<double,1> >;
template class PARTICLE_LEVELSET_EVOLUTION_MULTIPLE_UNIFORM<VECTOR<double,2> >;
template class PARTICLE_LEVELSET_EVOLUTION_MULTIPLE_UNIFORM<VECTOR<double,3> >;
}
