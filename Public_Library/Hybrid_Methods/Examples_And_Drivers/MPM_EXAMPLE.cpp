//#####################################################################
// Copyright 2015, Craig Schroeder.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Tools/Log/LOG.h>
#include <Geometry/Geometry_Particles/DEBUG_PARTICLES.h>
#include <Deformables/Forces/DEFORMABLES_FORCES.h>
#include <Hybrid_Methods/Collisions/MPM_COLLISION_OBJECT.h>
#include <Hybrid_Methods/Examples_And_Drivers/MPM_EXAMPLE.h>
#include <Hybrid_Methods/Examples_And_Drivers/MPM_PARTICLES.h>
#include <Hybrid_Methods/Forces/PARTICLE_GRID_FORCES.h>
#include <Hybrid_Methods/Iterators/GATHER_SCATTER.h>
#include <Hybrid_Methods/Iterators/PARTICLE_GRID_WEIGHTS.h>
#include <Hybrid_Methods/System/MPM_KRYLOV_VECTOR.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class TV> MPM_EXAMPLE<TV>::
MPM_EXAMPLE(const STREAM_TYPE stream_type)
    :stream_type(stream_type),particles(*new MPM_PARTICLES<TV>),debug_particles(*new DEBUG_PARTICLES<TV>),
    weights(0),gather_scatter(*new GATHER_SCATTER<TV>(simulated_particles)),initial_time(0),last_frame(100),
    write_substeps_level(-1),substeps_delay_frame(-1),output_directory("output"),
    restart(0),dt(0),time(0),frame_dt((T)1/24),min_dt(0),max_dt(frame_dt),ghost(3),
    use_reduced_rasterization(false),use_affine(false),use_midpoint(false),
    use_particle_collision(false),print_stats(false),flip(0),cfl(1),newton_tolerance(1),
    newton_iterations(10),solver_tolerance(1e-4),solver_iterations(1000),test_diff(false),threads(1)
{
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> MPM_EXAMPLE<TV>::
~MPM_EXAMPLE()
{
    delete &particles;
    delete &debug_particles;
    delete weights;
    delete &gather_scatter;
    collision_objects.Delete_Pointers_And_Clean_Memory();
    forces.Delete_Pointers_And_Clean_Memory();
    lagrangian_forces.Delete_Pointers_And_Clean_Memory();
    av.Delete_Pointers_And_Clean_Memory();
}
//#####################################################################
// Function Write_Output_Files
//#####################################################################
template<class TV> void MPM_EXAMPLE<TV>::
Write_Output_Files(const int frame)
{
    std::string f=STRING_UTILITIES::string_sprintf("%d",frame);

    FILE_UTILITIES::Write_To_File(stream_type,output_directory+"/common/grid",grid);
    FILE_UTILITIES::Write_To_File(stream_type,STRING_UTILITIES::string_sprintf("%s/%d/mpm_particles",output_directory.c_str(),frame),particles);
    FILE_UTILITIES::Write_To_File(stream_type,STRING_UTILITIES::string_sprintf("%s/%d/centered_velocities",output_directory.c_str(),frame),velocity);
    FILE_UTILITIES::Write_To_File(stream_type,STRING_UTILITIES::string_sprintf("%s/%d/density",output_directory.c_str(),frame),mass);
    FILE_UTILITIES::Write_To_File(stream_type,STRING_UTILITIES::string_sprintf("%s/%d/restart_data",output_directory.c_str(),frame),time);

    for(int i=0;i<particles.X.m;i++){
        Add_Debug_Particle(particles.X(i),VECTOR<T,3>(1,1,1));
        Debug_Particle_Set_Attribute<TV>(ATTRIBUTE_ID_V,particles.V(i));}
    debug_particles.Write_Debug_Particles(stream_type,output_directory,frame);
}
//#####################################################################
// Function Read_Output_Files
//#####################################################################
template<class TV> void MPM_EXAMPLE<TV>::
Read_Output_Files(const int frame)
{
    std::string f=STRING_UTILITIES::string_sprintf("%d",frame);
    FILE_UTILITIES::Read_From_File(stream_type,STRING_UTILITIES::string_sprintf("%s/%d/mpm_particles",output_directory.c_str(),frame),particles);
    FILE_UTILITIES::Read_From_File(stream_type,STRING_UTILITIES::string_sprintf("%s/%d/restart_data",output_directory.c_str(),frame),time);
}
//#####################################################################
// Function Precompute_Forces
//#####################################################################
template<class TV> void MPM_EXAMPLE<TV>::
Capture_Stress()
{
    for(int i=0;i<forces.m;i++)
        forces(i)->Capture_Stress();
}
//#####################################################################
// Function Precompute_Forces
//#####################################################################
template<class TV> void MPM_EXAMPLE<TV>::
Precompute_Forces(const T time)
{
    for(int i=0;i<forces.m;i++)
        forces(i)->Precompute(time);
    for(int i=0;i<lagrangian_forces.m;i++)
        lagrangian_forces(i)->Update_Position_Based_State(time,false);
}
//#####################################################################
// Function Potential_Energy
//#####################################################################
template<class TV> typename TV::SCALAR MPM_EXAMPLE<TV>::
Potential_Energy(const T time) const
{
    typename TV::SCALAR pe=0;
    for(int i=0;i<forces.m;i++)
        pe+=forces(i)->Potential_Energy(time);
    for(int i=0;i<lagrangian_forces.m;i++)
        pe+=lagrangian_forces(i)->Potential_Energy(time);
    return pe;
}
//#####################################################################
// Function Add_Forces
//#####################################################################
template<class TV> void MPM_EXAMPLE<TV>::
Add_Forces(ARRAY<TV,TV_INT>& F,const T time) const
{
    for(int i=0;i<forces.m;i++)
        forces(i)->Add_Forces(F,time);

    if(!lagrangian_forces.m) return;
    lagrangian_forces_F.Remove_All();
    lagrangian_forces_F.Resize(particles.X.m);
    for(int i=0;i<lagrangian_forces.m;i++)
        lagrangian_forces(i)->Add_Velocity_Independent_Forces(lagrangian_forces_F,time);
    gather_scatter.Scatter(
        [this,&F](int p,const PARTICLE_GRID_ITERATOR<TV>& it,int tid){
            F(it.Index())+=it.Weight()*lagrangian_forces_F(p);},false);
}
//#####################################################################
// Function Add_Hessian_Times
//#####################################################################
template<class TV> void MPM_EXAMPLE<TV>::
Add_Hessian_Times(ARRAY<TV,TV_INT>& F,const ARRAY<TV,TV_INT>& V,const T time) const
{
    for(int i=0;i<forces.m;i++)
        forces(i)->Add_Hessian_Times(F,V,time);

    if(!lagrangian_forces.m) return;
    lagrangian_forces_F.Remove_All();
    lagrangian_forces_V.Remove_All();
    lagrangian_forces_F.Resize(particles.X.m);
    lagrangian_forces_V.Resize(particles.X.m);
    gather_scatter.Gather(
        [this,&V](int p,const PARTICLE_GRID_ITERATOR<TV>& it,int tid){
            lagrangian_forces_V(p)+=it.Weight()*V(it.Index());},false);
    for(int i=0;i<lagrangian_forces.m;i++)
        lagrangian_forces(i)->Add_Implicit_Velocity_Independent_Forces(lagrangian_forces_V,lagrangian_forces_F,1,time);
    gather_scatter.Scatter(
        [this,&F](int p,const PARTICLE_GRID_ITERATOR<TV>& it,int tid){
            F(it.Index())+=it.Weight()*lagrangian_forces_F(p);},false);
}
//#####################################################################
// Function Add_Force
//#####################################################################
template<class TV> int MPM_EXAMPLE<TV>::
Add_Force(PARTICLE_GRID_FORCES<TV>& force)
{
    return forces.Append(&force);
}
//#####################################################################
// Function Add_Force
//#####################################################################
template<class TV> int MPM_EXAMPLE<TV>::
Add_Force(DEFORMABLES_FORCES<TV>& force)
{
    return lagrangian_forces.Append(&force);
}
//#####################################################################
// Function Set_Weights
//#####################################################################
template<class TV> void MPM_EXAMPLE<TV>::
Set_Weights(PARTICLE_GRID_WEIGHTS<TV>* weights_input)
{
    weights=weights_input;
    gather_scatter.Set_Weights(weights);
}
//#####################################################################
// Function Total_Particle_Linear_Momentum
//#####################################################################
template<class TV> TV MPM_EXAMPLE<TV>::
Total_Particle_Linear_Momentum() const
{
    TV result;
    for(int k=0;k<simulated_particles.m;k++){
        int p=simulated_particles(k);
        result+=particles.mass(p)*particles.V(p);}
    return result;
}
//#####################################################################
// Function Total_Grid_Linear_Momentum
//#####################################################################
template<class TV> TV MPM_EXAMPLE<TV>::
Total_Grid_Linear_Momentum(const ARRAY<TV,TV_INT>& u) const
{
    TV result;
    for(int i=0;i<valid_grid_indices.m;i++){
        int j=valid_grid_indices(i);
        result+=mass.array(j)*u.array(j);}
    return result;
}
//#####################################################################
// Function Total_Grid_Angular_Momentum
//#####################################################################
template<class TV> typename TV::SPIN MPM_EXAMPLE<TV>::
Total_Grid_Angular_Momentum(T dt,const ARRAY<TV,TV_INT>& u,const ARRAY<TV,TV_INT>* u0) const
{
    typename TV::SPIN result;
    for(int k=0;k<valid_grid_indices.m;k++){
        int i=valid_grid_indices(k);
        TV X=location.array(i);
        if(use_midpoint && u0) X+=dt/2*u0->array(i);
        result+=mass.array(i)*TV::Cross_Product(X,u.array(i));}
    return result;
}
//#####################################################################
// Function Total_Grid_Kinetic_Energy
//#####################################################################
template<class TV> typename TV::SCALAR MPM_EXAMPLE<TV>::
Total_Grid_Kinetic_Energy(const ARRAY<TV,TV_INT>& u) const
{
    T result=0;
    for(int i=0;i<valid_grid_indices.m;i++){
        int j=valid_grid_indices(i);
        result+=(T).5*mass.array(j)*u.array(j).Magnitude_Squared();}
    return result;
}
//#####################################################################
namespace PhysBAM{
template class MPM_EXAMPLE<VECTOR<float,2> >;
template class MPM_EXAMPLE<VECTOR<float,3> >;
template class MPM_EXAMPLE<VECTOR<double,2> >;
template class MPM_EXAMPLE<VECTOR<double,3> >;
}
