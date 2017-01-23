//#####################################################################
// Copyright 2015, Craig Schroeder.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Tools/Read_Write/OCTAVE_OUTPUT.h>
#include <Grid_PDE/Poisson/LAPLACE_UNIFORM.h>
#include <Grid_PDE/Poisson/PROJECTION_UNIFORM.h>
#include <Hybrid_Methods/Collisions/MPM_COLLISION_IMPLICIT_OBJECT.h>
#include <Hybrid_Methods/Examples_And_Drivers/MPM_MAC_EXAMPLE.h>
#include <Hybrid_Methods/Examples_And_Drivers/MPM_PARTICLES.h>
#include <Hybrid_Methods/Forces/MPM_FORCE_HELPER.h>
#include <Hybrid_Methods/Forces/MPM_PLASTICITY_MODEL.h>
#include <Hybrid_Methods/Iterators/GATHER_SCATTER.h>
#include <Hybrid_Methods/Iterators/PARTICLE_GRID_WEIGHTS_SPLINE.h>
#include <Hybrid_Methods/Projection/MPM_PROJECTION_SYSTEM.h>
#include <Hybrid_Methods/Projection/MPM_PROJECTION_VECTOR.h>
#include <Hybrid_Methods/System/MPM_KRYLOV_VECTOR.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class TV> MPM_MAC_EXAMPLE<TV>::
MPM_MAC_EXAMPLE(const STREAM_TYPE stream_type)
    :stream_type(stream_type),particles(*new MPM_PARTICLES<TV>),
    projection_system(*new MPM_PROJECTION_SYSTEM<TV>),
    sol(*new MPM_PROJECTION_VECTOR<TV>),rhs(*new MPM_PROJECTION_VECTOR<TV>),
    ghost(3),gather_scatter(*new GATHER_SCATTER<TV>(grid,simulated_particles)),
    use_affine(true),use_early_gradient_transfer(false),initial_time(0),
    last_frame(100),write_substeps_level(-1),substeps_delay_frame(-1),
    output_directory("output"),data_directory("../../Public_Data"),use_test_output(false),
    restart(0),dt(0),time(0),frame_dt((T)1/24),min_dt(0),max_dt(frame_dt),
    only_write_particles(false),cfl(1),
    solver_tolerance(std::numeric_limits<T>::epsilon()*10),solver_iterations(1000),
    threads(1),use_particle_volumes(false),debug_particles(*new DEBUG_PARTICLES<TV>),
    print_stats(false),last_te(0),last_grid_ke(0),test_system(false)
{
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> MPM_MAC_EXAMPLE<TV>::
~MPM_MAC_EXAMPLE()
{
    delete &particles;
    delete &debug_particles;
    for(int i=0;i<TV::m;i++) delete weights(i);
    delete &gather_scatter;
    collision_objects.Delete_Pointers_And_Clean_Memory();
    fluid_walls.Delete_Pointers_And_Clean_Memory();
}
//#####################################################################
// Function Write_Output_Files
//#####################################################################
template<class TV> void MPM_MAC_EXAMPLE<TV>::
Write_Output_Files(const int frame)
{
    std::string f=LOG::sprintf("%d",frame);
    if(this->use_test_output){
        std::string file=LOG::sprintf("%s/%s-%03d.txt",output_directory.c_str(),test_output_prefix.c_str(),frame);
        OCTAVE_OUTPUT<T> oo(file.c_str());
        oo.Write("X",particles.X.Flattened());
        oo.Write("V",particles.V.Flattened());
        oo.Write("u",velocity.array);}

#pragma omp parallel
#pragma omp single
    {
#pragma omp task
        FILE_UTILITIES::Write_To_File(stream_type,output_directory+"/common/grid",grid);
#pragma omp task
        FILE_UTILITIES::Write_To_File(stream_type,LOG::sprintf("%s/%d/restart_data",output_directory.c_str(),frame),time);
#pragma omp task
        FILE_UTILITIES::Write_To_File(stream_type,LOG::sprintf("%s/%d/mpm_particles",output_directory.c_str(),frame),particles);

        if(!only_write_particles){
#pragma omp task
            FILE_UTILITIES::Write_To_File(stream_type,LOG::sprintf("%s/%d/mac_velocities",output_directory.c_str(),frame),velocity);
#pragma omp task
            {
                GRID<TV> ghost_grid(grid.numbers_of_cells+2*ghost,grid.Ghost_Domain(ghost),true);
                for(int i=0;i<collision_objects.m;i++)
                    if(IMPLICIT_OBJECT<TV>* io=collision_objects(i)->Get_Implicit_Object(time))
                        Dump_Levelset(ghost_grid,*io,VECTOR<T,3>(0.7,0.3,0.3));
                debug_particles.Write_Debug_Particles(stream_type,output_directory,frame);
            }
        }
    }
}
//#####################################################################
// Function Read_Output_Files
//#####################################################################
template<class TV> void MPM_MAC_EXAMPLE<TV>::
Read_Output_Files(const int frame)
{
    std::string f=LOG::sprintf("%d",frame);
    FILE_UTILITIES::Read_From_File(stream_type,LOG::sprintf("%s/%d/mpm_particles",output_directory.c_str(),frame),particles);
    FILE_UTILITIES::Read_From_File(stream_type,LOG::sprintf("%s/%d/restart_data",output_directory.c_str(),frame),time);
}
//#####################################################################
// Function Potential_Energy
//#####################################################################
template<class TV> typename TV::SCALAR MPM_MAC_EXAMPLE<TV>::
Potential_Energy(const T time) const
{
    // TODO:
    return 0;
}
//#####################################################################
// Function Add_Forces
//#####################################################################
template<class TV> void MPM_MAC_EXAMPLE<TV>::
Add_Forces(ARRAY<TV,TV_INT>& F,const T time) const
{
    // TODO:
}
//#####################################################################
// Function Set_Weights
//#####################################################################
template<class TV> void MPM_MAC_EXAMPLE<TV>::
Set_Weights(int order)
{
    for(int i=0;i<TV::m;i++){
        GRID<TV> face_grid=grid.Get_Face_MAC_Grid(i);
        if(order==1)
            weights(i)=new PARTICLE_GRID_WEIGHTS_SPLINE<TV,1>(face_grid);
        else if(order==2)
            weights(i)=new PARTICLE_GRID_WEIGHTS_SPLINE<TV,2>(face_grid);
        else if(order==3)
            weights(i)=new PARTICLE_GRID_WEIGHTS_SPLINE<TV,3>(face_grid);
        else PHYSBAM_FATAL_ERROR("Unrecognized interpolation order");}
    gather_scatter.face_weights=weights;
}
//#####################################################################
// Function Total_Particle_Linear_Momentum
//#####################################################################
template<class TV> TV MPM_MAC_EXAMPLE<TV>::
Total_Particle_Linear_Momentum() const
{
    TV result;
#pragma omp parallel for
    for(int t=0;t<threads;t++){
        int a=t*simulated_particles.m/threads;
        int b=(t+1)*simulated_particles.m/threads;
        TV result_local;
        for(int k=a;k<b;k++){
            int p=simulated_particles(k);
            result_local+=particles.mass(p)*particles.V(p);}
#pragma omp critical
        result+=result_local;}
    return result;
}
//#####################################################################
// Function Total_Grid_Linear_Momentum
//#####################################################################
template<class TV> TV MPM_MAC_EXAMPLE<TV>::
Total_Grid_Linear_Momentum() const
{
    TV result;
#pragma omp parallel for
    for(int t=0;t<threads;t++){
        int a=t*valid_flat_indices.m/threads;
        int b=(t+1)*valid_flat_indices.m/threads;
        TV result_local;
        for(int k=a;k<b;k++){
            int j=valid_flat_indices(k);
            result_local(valid_indices(k).axis)+=mass.array(j)*velocity.array(j);}
#pragma omp critical
        result+=result_local;}
    return result;
}
//#####################################################################
// Function Total_Grid_Angular_Momentum
//#####################################################################
template<class TV> typename TV::SPIN MPM_MAC_EXAMPLE<TV>::
Total_Particle_Angular_Momentum() const
{
    typename TV::SPIN result;
#pragma omp parallel for
    for(int t=0;t<threads;t++){
        int a=t*simulated_particles.m/threads;
        int b=(t+1)*simulated_particles.m/threads;
        typename TV::SPIN result_local;
        for(int k=a;k<b;k++){
            int p=simulated_particles(k);
            result_local+=particles.mass(p)*particles.X(p).Cross(particles.V(p));
            if(particles.store_B) result_local-=particles.mass(p)*particles.B(p).Contract_Permutation_Tensor();}
#pragma omp critical
        result+=result_local;}
    return result;
}
//#####################################################################
// Function Total_Grid_Angular_Momentum
//#####################################################################
template<class TV> typename TV::SPIN MPM_MAC_EXAMPLE<TV>::
Total_Grid_Angular_Momentum(T dt) const
{
    typename TV::SPIN result;
#pragma omp parallel for
    for(int t=0;t<threads;t++){
        int a=t*valid_flat_indices.m/threads;
        int b=(t+1)*valid_flat_indices.m/threads;
        typename TV::SPIN result_local;
        for(int k=a;k<b;k++){
            int i=valid_flat_indices(k);
            TV X=location.array(i);
            result_local+=mass.array(i)*TV::Cross_Product(X,velocity.array(i)*TV::Axis_Vector(valid_indices(k).axis));}
#pragma omp critical
        result+=result_local;}
    return result;
}
//#####################################################################
// Function Total_Grid_Kinetic_Energy
//#####################################################################
template<class TV> typename TV::SCALAR MPM_MAC_EXAMPLE<TV>::
Total_Grid_Kinetic_Energy() const
{
    T result=0;
#pragma omp parallel for reduction(+:result)
    for(int i=0;i<valid_flat_indices.m;i++){
        int j=valid_flat_indices(i);
        result+=(T).5*mass.array(j)*sqr(velocity.array(j));}
    return result;
}
//#####################################################################
// Function Total_Particle_Kinetic_Energy
//#####################################################################
template<class TV> typename TV::SCALAR MPM_MAC_EXAMPLE<TV>::
Total_Particle_Kinetic_Energy() const
{
    T result=0,Dp_inverse=0;
    if(use_affine && weights(0)->constant_scalar_inertia_tensor)
        Dp_inverse=weights(0)->Constant_Scalar_Inverse_Dp();
#pragma omp parallel for reduction(+:result)
    for(int k=0;k<simulated_particles.m;k++){
        int p=simulated_particles(k);
        T result_local=particles.mass(p)/2*particles.V(p).Magnitude_Squared();
        if(particles.store_B) result_local+=particles.mass(p)/2*Dp_inverse*particles.B(p).Frobenius_Norm_Squared();
        result+=result_local;}
    return result;
}
//#####################################################################
// Function Average_Particle_Mass
//#####################################################################
template<class TV> typename TV::SCALAR MPM_MAC_EXAMPLE<TV>::
Average_Particle_Mass() const
{
    T result=0;
#pragma omp parallel for reduction(+:result)
    for(int k=0;k<simulated_particles.m;k++)
        result+=particles.mass(simulated_particles(k));
    return result/(T)particles.number;
}
//#####################################################################
// Function Add_Collision_Object
//#####################################################################
template<class TV> void MPM_MAC_EXAMPLE<TV>::
Add_Collision_Object(IMPLICIT_OBJECT<TV>* io,COLLISION_TYPE type,T friction,std::function<FRAME<TV>(T)> func_frame,std::function<TWIST<TV>(T)> func_twist)
{
    collision_objects.Append(new MPM_COLLISION_IMPLICIT_OBJECT<TV>(io,type,friction,func_frame,func_twist));
}
//#####################################################################
// Function Add_Fluid_Wall
//#####################################################################
template<class TV> void MPM_MAC_EXAMPLE<TV>::
Add_Fluid_Wall(IMPLICIT_OBJECT<TV>* io)
{
    fluid_walls.Append(io);
}
//#####################################################################
namespace PhysBAM{
template class MPM_MAC_EXAMPLE<VECTOR<float,2> >;
template class MPM_MAC_EXAMPLE<VECTOR<float,3> >;
template class MPM_MAC_EXAMPLE<VECTOR<double,2> >;
template class MPM_MAC_EXAMPLE<VECTOR<double,3> >;
}
