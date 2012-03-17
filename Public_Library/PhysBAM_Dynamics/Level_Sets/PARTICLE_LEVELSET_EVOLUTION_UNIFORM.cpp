//#####################################################################
// Copyright 2002-2009, Ronald Fedkiw, Frederic Gibou, Geoffrey Irving, Frank Losasso, Neil Molino, Avi Robinson-Mosher, Tamar Shinar, Elliot English.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Grids_Uniform/GRID.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_CELL.h>
#include <PhysBAM_Tools/Grids_Uniform_Boundaries/BOUNDARY_UNIFORM.h>
#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Tools/Ordinary_Differential_Equations/RUNGEKUTTA.h>
#include <PhysBAM_Geometry/Grids_Uniform_Level_Sets/FAST_LEVELSET.h>
#include <PhysBAM_Dynamics/Level_Sets/LEVELSET_CALLBACKS.h>
#include <PhysBAM_Dynamics/Level_Sets/PARTICLE_LEVELSET_EVOLUTION_UNIFORM.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class T_GRID> PARTICLE_LEVELSET_EVOLUTION_UNIFORM<T_GRID>::
PARTICLE_LEVELSET_EVOLUTION_UNIFORM(const T_GRID& grid_input,const int number_of_ghost_cells_input)
    :grid(grid_input),particle_levelset(grid,phi,number_of_ghost_cells_input),rungekutta_phi(0),levelset_advection(&particle_levelset.levelset)
{
    Use_Semi_Lagrangian_Advection();
    Track_Mass(false);
}
//#####################################################################
// Destructor
//#####################################################################
template<class T_GRID> PARTICLE_LEVELSET_EVOLUTION_UNIFORM<T_GRID>::
~PARTICLE_LEVELSET_EVOLUTION_UNIFORM()
{
    delete rungekutta_phi;
}
//#####################################################################
// Function Initialize_Runge_Kutta
//#####################################################################
template<class T_GRID> void PARTICLE_LEVELSET_EVOLUTION_UNIFORM<T_GRID>::
Initialize_Runge_Kutta()
{
    if(runge_kutta_order_levelset > 1){
        delete rungekutta_phi;rungekutta_phi=new RUNGEKUTTA<T_ARRAYS_SCALAR>(phi);
        rungekutta_phi->Set_Order(runge_kutta_order_levelset);rungekutta_phi->Set_Time(time);
        rungekutta_phi->Set_Grid_And_Boundary_Condition(grid,*particle_levelset.levelset.boundary);}
}
//#####################################################################
// Function Time_Step
//#####################################################################
template<class T_GRID> typename T_GRID::SCALAR PARTICLE_LEVELSET_EVOLUTION_UNIFORM<T_GRID>::
Time_Step(const T stopping_time,bool& limited_by_stopping_time)
{
    T dt=CFL();limited_by_stopping_time=false;
    if(time+dt >= stopping_time){dt=stopping_time-time;limited_by_stopping_time=true;}else if(time+2*dt >= stopping_time) dt=(T).51*(stopping_time-time);
    return dt;
}
//#####################################################################
// Function CFL
//#####################################################################
template<class T_GRID> typename T_GRID::SCALAR PARTICLE_LEVELSET_EVOLUTION_UNIFORM<T_GRID>::
CFL(const bool need_to_get_velocity,const bool analytic_test)
{
    if(need_to_get_velocity) particle_levelset.levelset.levelset_callbacks->Get_Levelset_Velocity(grid,particle_levelset.levelset,V,time);
    if(analytic_test) return cfl_number/max((V.Maxabs()/grid.dX).Max(),1/particle_levelset.levelset.max_time_step);
    return cfl_number*particle_levelset.levelset.CFL(V);
}
//#####################################################################
// Function Advance_To_Time
//#####################################################################
template<class T_GRID> void PARTICLE_LEVELSET_EVOLUTION_UNIFORM<T_GRID>::
Advance_To_Time(T_FACE_ARRAYS_SCALAR* face_velocities,const T stopping_time,const bool verbose)
{
    int substep=0;bool done=false;
    while(!done){substep++;
        T dt=Time_Step(stopping_time,done);
        if(verbose) LOG::cout << "substep = " << substep << ", dt = " << dt << std::endl;
        Advance_One_Time_Step(face_velocities,dt);}
    if(track_mass){
        T mass=levelset_advection.Approximate_Negative_Material();
        LOG::cout << "negative material = " << mass << " - change = " << (mass-initial_mass)/initial_mass*100 << "%" << std::endl;}
}
//#####################################################################
// Function Advance_One_Time_Step
//#####################################################################
template<class T_GRID> void PARTICLE_LEVELSET_EVOLUTION_UNIFORM<T_GRID>::
Advance_One_Time_Step(T_FACE_ARRAYS_SCALAR* face_velocities,const T dt)
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
// only does an Euler step for the level set
template<class T_GRID> void PARTICLE_LEVELSET_EVOLUTION_UNIFORM<T_GRID>::
Advance_Levelset(const T dt)
{
    if(runge_kutta_order_levelset == 1){
        particle_levelset.levelset.levelset_callbacks->Get_Levelset_Velocity(grid,particle_levelset.levelset,V,time);
        levelset_advection.Euler_Step(V,dt,time,particle_levelset.number_of_ghost_cells);time+=dt;}
    else{
        if(!rungekutta_phi) Initialize_Runge_Kutta();
        rungekutta_phi->Start(dt);
        for(int k=0;k<runge_kutta_order_levelset;k++){
            if(k == 1 || !use_frozen_velocity) particle_levelset.levelset.levelset_callbacks->Get_Levelset_Velocity(grid,particle_levelset.levelset,V,time);
            levelset_advection.Euler_Step(V,dt,time,particle_levelset.number_of_ghost_cells);
            time=rungekutta_phi->Main();}}
}
//#####################################################################
// Function Advance_Particles
//#####################################################################
template<class T_GRID> void PARTICLE_LEVELSET_EVOLUTION_UNIFORM<T_GRID>::
Advance_Particles(const T_FACE_ARRAYS_SCALAR& face_velocities,const T dt,const bool analytic_test)
{
    if(use_particle_levelset){
        time-=dt; // to fix up time advancement due to Advance_Levelset()
        if(runge_kutta_order_particles == 1 || (runge_kutta_order_particles == 2 && use_frozen_velocity)){
            // TODO: still needed?
            //particle_levelset.levelset.levelset_callbacks->Get_Levelset_Velocity(grid,particle_levelset.levelset,V,time);
            particle_levelset.Euler_Step_Particles(face_velocities,dt,time,runge_kutta_order_particles==2,true,analytic_test);time+=dt;}
        else if(runge_kutta_order_particles == 2 || runge_kutta_order_particles == 3){
            T start_time=time;
            time=Advance_Particles(particle_levelset.positive_particles,PARTICLE_LEVELSET_POSITIVE,dt,start_time);
            time=Advance_Particles(particle_levelset.negative_particles,PARTICLE_LEVELSET_NEGATIVE,dt,start_time);
            if(analytic_test){
                time=Advance_Particles(particle_levelset.removed_positive_particles,PARTICLE_LEVELSET_POSITIVE,dt,start_time);
                time=Advance_Particles(particle_levelset.removed_negative_particles,PARTICLE_LEVELSET_NEGATIVE,dt,start_time);}
            else particle_levelset.Euler_Step_Removed_Particles(dt,start_time,true);}} // can only Euler Step removed particles
}
//#####################################################################
// Function Advance_Particles
//#####################################################################
template<class T_GRID> typename T_GRID::SCALAR PARTICLE_LEVELSET_EVOLUTION_UNIFORM<T_GRID>::
Advance_Particles(T_ARRAYS_PARTICLE_LEVELSET_PARTICLES& particles,const PARTICLE_LEVELSET_PARTICLE_TYPE particle_type,const T dt,const T input_time)
{
    T_GRID mac_grid=grid.Get_MAC_Grid();
    T current_time=input_time;
    T_ARRAYS_RUNGEKUTTA rungekutta_particles;rungekutta_particles.Resize(mac_grid.Domain_Indices(1));
    for(CELL_ITERATOR iterator(grid,1);iterator.Valid();iterator.Next()) if(particles(iterator.Cell_Index()))
        rungekutta_particles(iterator.Cell_Index())=RUNGEKUTTA<ARRAY_VIEW<TV> >::Create(particles(iterator.Cell_Index())->X,runge_kutta_order_particles,dt,current_time);
    for(int k=0;k<runge_kutta_order_particles;k++){
        if(k == 1 || !use_frozen_velocity) particle_levelset.levelset.levelset_callbacks->Get_Levelset_Velocity(grid,particle_levelset.levelset,V,current_time);
        particle_levelset.Euler_Step_Particles_Wrapper(V,particles,particle_type,dt,current_time,false,k==0,false);
        for(CELL_ITERATOR iterator(grid,1);iterator.Valid();iterator.Next()) if(particles(iterator.Cell_Index())) current_time=rungekutta_particles(iterator.Cell_Index())->Main();}
    for(CELL_ITERATOR iterator(grid,1);iterator.Valid();iterator.Next()) if(particles(iterator.Cell_Index())){
        PARTICLE_LEVELSET_PARTICLES<TV>& cell_particles=*particles(iterator.Cell_Index());
        for(int k=0;k<cell_particles.Size();k++){
            TV velocity;
            particle_levelset.levelset.levelset_callbacks->Adjust_Particle_For_Domain_Boundaries(cell_particles,k,velocity,particle_type,(T)1,current_time);
            cell_particles.X(k)+=velocity;}}
    particle_levelset.Update_Particle_Cells(particles);
    rungekutta_particles.Delete_Pointers_And_Clean_Memory();
    return current_time;
}
//#####################################################################
// Function Advance_Particles
//#####################################################################
template<class T_GRID> typename T_GRID::SCALAR PARTICLE_LEVELSET_EVOLUTION_UNIFORM<T_GRID>::
Advance_Particles(T_ARRAYS_PARTICLE_LEVELSET_REMOVED_PARTICLES& particles,const PARTICLE_LEVELSET_PARTICLE_TYPE particle_type,const T dt,const T input_time)
{
    T_GRID mac_grid=grid.Get_MAC_Grid();
    T current_time=input_time;
    T_ARRAYS_RUNGEKUTTA rungekutta_particles;rungekutta_particles.Resize(mac_grid.Domain_Indices(1));
    for(CELL_ITERATOR iterator(grid,1);iterator.Valid();iterator.Next()) if(particles(iterator.Cell_Index()))
        rungekutta_particles(iterator.Cell_Index())=RUNGEKUTTA<ARRAY_VIEW<TV> >::Create(particles(iterator.Cell_Index())->X,runge_kutta_order_particles,dt,current_time);
    for(int k=0;k<runge_kutta_order_particles;k++){
        if(k == 1 || !use_frozen_velocity) particle_levelset.levelset.levelset_callbacks->Get_Levelset_Velocity(grid,particle_levelset.levelset,V,current_time);
        particle_levelset.Euler_Step_Particles_Wrapper(V,particles,particle_type,dt,current_time,false,k==0,false);
        for(CELL_ITERATOR iterator(grid,1);iterator.Valid();iterator.Next()) if(particles(iterator.Cell_Index())) current_time=rungekutta_particles(iterator.Cell_Index())->Main();}
    for(CELL_ITERATOR iterator(grid,1);iterator.Valid();iterator.Next()) if(particles(iterator.Cell_Index())){
        PARTICLE_LEVELSET_PARTICLES<TV>& cell_particles=*particles(iterator.Cell_Index());
        for(int k=0;k<cell_particles.Size();k++){
            TV velocity;
            particle_levelset.levelset.levelset_callbacks->Adjust_Particle_For_Domain_Boundaries(cell_particles,k,velocity,particle_type,(T)1,current_time);
            cell_particles.X(k)+=velocity;}}
    particle_levelset.Update_Particle_Cells(particles);
    rungekutta_particles.Delete_Pointers_And_Clean_Memory();
    return input_time+dt; // there may be no removed particles
}
//#####################################################################
// Function Modify_Levelset_And_Particles
//#####################################################################
template<class T_GRID> void PARTICLE_LEVELSET_EVOLUTION_UNIFORM<T_GRID>::
Modify_Levelset_And_Particles(T_FACE_ARRAYS_SCALAR* face_velocities)
{
    // TODO: a call for creating particles from the geometry if necessary
    if(use_particle_levelset){
        LOG::Time("modifying levelset");
        particle_levelset.Modify_Levelset_Using_Escaped_Particles(face_velocities);}
    LOG::Time("reinitializing levelset");
    Make_Signed_Distance();
    if(use_particle_levelset){
        LOG::Time("modifying levelset");        
        particle_levelset.Modify_Levelset_Using_Escaped_Particles(face_velocities);
        LOG::Time("adjusting particle radii");
        particle_levelset.Adjust_Particle_Radii();
        LOG::Stop_Time();}
}
//#####################################################################
// Function Reseed_Particles
//#####################################################################
template<class T_GRID> void PARTICLE_LEVELSET_EVOLUTION_UNIFORM<T_GRID>::
Reseed_Particles(const T time,const int time_step,T_ARRAYS_BOOL* cell_centered_mask,const bool verbose)
{
    if((use_particle_levelset && cell_centered_mask) || (!time_step || (reseeding_frequency && time_step%reseeding_frequency == 0))){
        int new_particles=particle_levelset.Reseed_Particles(time,cell_centered_mask);
        if(verbose) LOG::cout << "Reseeding... " << new_particles << " new particles" << std::endl;
        Initialize_Runge_Kutta();} // need to reset based on new number of particles
}
//#####################################################################
template class PARTICLE_LEVELSET_EVOLUTION_UNIFORM<GRID<VECTOR<float,1> > >;
template class PARTICLE_LEVELSET_EVOLUTION_UNIFORM<GRID<VECTOR<float,2> > >;
template class PARTICLE_LEVELSET_EVOLUTION_UNIFORM<GRID<VECTOR<float,3> > >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class PARTICLE_LEVELSET_EVOLUTION_UNIFORM<GRID<VECTOR<double,1> > >;
template class PARTICLE_LEVELSET_EVOLUTION_UNIFORM<GRID<VECTOR<double,2> > >;
template class PARTICLE_LEVELSET_EVOLUTION_UNIFORM<GRID<VECTOR<double,3> > >;
#endif
