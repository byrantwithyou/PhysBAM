//#####################################################################
// Copyright 2015, Craig Schroeder.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Tools/Read_Write/OCTAVE_OUTPUT.h>
#include <Deformables/Deformable_Objects/DEFORMABLE_BODY_COLLECTION.h>
#include <Deformables/Forces/LAGGED_FORCE.h>
#include <Solids/Solids/SOLID_BODY_COLLECTION.h>
#include <Hybrid_Methods/Collisions/MPM_COLLISION_OBJECT.h>
#include <Hybrid_Methods/Examples_And_Drivers/MPM_EXAMPLE.h>
#include <Hybrid_Methods/Examples_And_Drivers/MPM_PARTICLES.h>
#include <Hybrid_Methods/Forces/MPM_FORCE_HELPER.h>
#include <Hybrid_Methods/Forces/MPM_PLASTICITY_MODEL.h>
#include <Hybrid_Methods/Iterators/GATHER_SCATTER.h>
#include <Hybrid_Methods/Iterators/PARTICLE_GRID_WEIGHTS_SPLINE.h>
#include <Hybrid_Methods/System/MPM_KRYLOV_VECTOR.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class TV> MPM_EXAMPLE<TV>::
MPM_EXAMPLE(const STREAM_TYPE stream_type)
    :stream_type(stream_type),particles(*new MPM_PARTICLES<TV>),
    solid_body_collection(*new SOLID_BODY_COLLECTION<TV>(&particles)),
    debug_particles(*new DEBUG_PARTICLES<TV>),
    lagrangian_forces(solid_body_collection.deformable_body_collection.deformables_forces),
    gather_scatter(*new GATHER_SCATTER<TV>(grid,simulated_particles)),
    force_helper(*new MPM_FORCE_HELPER<TV>(particles,quad_F_coeff))
{
    debug_particles.debug_particles.template Add_Array<T>("display_size");
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> MPM_EXAMPLE<TV>::
~MPM_EXAMPLE()
{
    delete &solid_body_collection;
    delete &particles;
    delete &debug_particles;
    delete weights;
    delete &gather_scatter;
    delete &force_helper;
    plasticity_models.Delete_Pointers_And_Clean_Memory();
    collision_objects.Delete_Pointers_And_Clean_Memory();
    forces.Delete_Pointers_And_Clean_Memory();
    av.Delete_Pointers_And_Clean_Memory();
}
//#####################################################################
// Function Write_Output_Files
//#####################################################################
template<class TV> void MPM_EXAMPLE<TV>::
Write_Output_Files()
{
    if(use_test_output){
        std::string file=LOG::sprintf("%s/%s-%03d.txt",viewer_dir.output_directory.c_str(),test_output_prefix.c_str(),viewer_dir.frame_stack(0));
        OCTAVE_OUTPUT<T> oo(file.c_str());
        oo.Write("X",particles.X.Flattened());
        oo.Write("V",particles.V.Flattened());
        oo.Write("u",velocity.array.Flattened());}

    for(int i=0;i<write_output_files.m;i++) write_output_files(i)();
#pragma omp parallel
#pragma omp single
    {
#pragma omp task
        Write_To_File(stream_type,viewer_dir.output_directory+"/common/grid",grid);
#pragma omp task
        if(!system(LOG::sprintf("rm -f %s/mpm_particles.gz ;  ln -s ./deformable_object_particles.gz %s/mpm_particles.gz",viewer_dir.current_directory,viewer_dir.current_directory).c_str())){}
#pragma omp task
        Write_To_File(stream_type,viewer_dir.current_directory+"/restart_data",time);
#pragma omp task
        {
            solid_body_collection.Write(stream_type,viewer_dir,write_structures_every_frame);
        }

        if(!only_write_particles){
#pragma omp task
            Write_To_File(stream_type,viewer_dir.current_directory+"/centered_velocities",velocity);
#pragma omp task
            {
                GRID<TV> ghost_grid(grid.numbers_of_cells+2*ghost,grid.Ghost_Domain(ghost),true);
                for(int i=0;i<collision_objects.m;i++)
                    if(IMPLICIT_OBJECT<TV>* io=collision_objects(i)->Get_Implicit_Object(time))
                        Dump_Levelset(ghost_grid,*io,VECTOR<T,3>(0.7,0.3,0.3));
                if(mass_contour>=0)
                    Dump_Levelset(grid,mass,VECTOR<T,3>(0.2,0.6,0.2),mass_contour*Average_Particle_Mass());
                debug_particles.Write_Debug_Particles(stream_type,viewer_dir);
            }
        }
    }
}
//#####################################################################
// Function Read_Output_Files
//#####################################################################
template<class TV> void MPM_EXAMPLE<TV>::
Read_Output_Files()
{
    Read_From_File(viewer_dir.current_directory+"/mpm_particles",particles);
    Read_From_File(viewer_dir.current_directory+"/restart_data",time);
    for(int i=0;i<read_output_files.m;i++) read_output_files(i)();
}
//#####################################################################
// Function Capture_Stress
//#####################################################################
template<class TV> void MPM_EXAMPLE<TV>::
Capture_Stress()
{
    force_helper.Fn=particles.F;
    if(particles.store_S) force_helper.Sn=particles.S;
}
//#####################################################################
// Function Precompute_Forces
//#####################################################################
template<class TV> void MPM_EXAMPLE<TV>::
Precompute_Forces(const T time,const T dt,const bool update_hessian)
{
    for(int i=0;i<forces.m;i++)
        forces(i)->Precompute(time,dt,true,update_hessian);
    solid_body_collection.Update_Position_Based_State(time,false,update_hessian);
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
    lagrangian_forces_F.Resize(particles.X.m,no_init);
#pragma omp parallel for
    for(int i=0;i<lagrangian_forces_F.m;i++)
        lagrangian_forces_F(i)=TV();
    solid_body_collection.deformable_body_collection.Add_Velocity_Independent_Forces(lagrangian_forces_F,time);
    gather_scatter.template Scatter<int>(false,
        [this,&F](int p,const PARTICLE_GRID_ITERATOR<TV>& it,int tid){
            F(it.Index())+=it.Weight()*lagrangian_forces_F(p);});
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
    lagrangian_forces_F.Resize(particles.X.m,no_init);
    lagrangian_forces_V.Resize(particles.X.m,no_init);
#pragma omp parallel for
    for(int i=0;i<lagrangian_forces_F.m;i++){
        lagrangian_forces_F(i)=TV();
        lagrangian_forces_V(i)=TV();}
    gather_scatter.template Gather<int>(false,
        [this,&V](int p,const PARTICLE_GRID_ITERATOR<TV>& it,int tid){
            lagrangian_forces_V(p)+=it.Weight()*V(it.Index());});
    solid_body_collection.deformable_body_collection.Add_Implicit_Velocity_Independent_Forces(lagrangian_forces_V,lagrangian_forces_F,time);
    gather_scatter.template Scatter<int>(false,
        [this,&F](int p,const PARTICLE_GRID_ITERATOR<TV>& it,int tid){
            F(it.Index())-=it.Weight()*lagrangian_forces_F(p);});
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
    return solid_body_collection.deformable_body_collection.Add_Force(&force);
}
//#####################################################################
// Function Set_Weights
//#####################################################################
template<class TV> void MPM_EXAMPLE<TV>::
Set_Weights(int order)
{
    delete weights;
    if(order==1) gather_scatter.weights=weights=new PARTICLE_GRID_WEIGHTS_SPLINE<TV,1>(grid,threads);
    else if(order==2) gather_scatter.weights=weights=new PARTICLE_GRID_WEIGHTS_SPLINE<TV,2>(grid,threads);
    else if(order==3) gather_scatter.weights=weights=new PARTICLE_GRID_WEIGHTS_SPLINE<TV,3>(grid,threads);
    else PHYSBAM_FATAL_ERROR("Unrecognized interpolation order");
}
//#####################################################################
// Function Total_Particle_Linear_Momentum
//#####################################################################
template<class TV> TV MPM_EXAMPLE<TV>::
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
template<class TV> TV MPM_EXAMPLE<TV>::
Total_Grid_Linear_Momentum(const ARRAY<TV,TV_INT>& u) const
{
    TV result;
#pragma omp parallel for
    for(int t=0;t<threads;t++){
        int a=t*valid_grid_indices.m/threads;
        int b=(t+1)*valid_grid_indices.m/threads;
        TV result_local;
        for(int k=a;k<b;k++){
            int j=valid_grid_indices(k);
            result_local+=mass.array(j)*u.array(j);}
#pragma omp critical
        result+=result_local;}
    return result;
}
//#####################################################################
// Function Total_Grid_Angular_Momentum
//#####################################################################
template<class TV> typename TV::SPIN MPM_EXAMPLE<TV>::
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
template<class TV> typename TV::SPIN MPM_EXAMPLE<TV>::
Total_Grid_Angular_Momentum(T dt,const ARRAY<TV,TV_INT>& u,const ARRAY<TV,TV_INT>* u0) const
{
    typename TV::SPIN result;
#pragma omp parallel for
    for(int t=0;t<threads;t++){
        int a=t*valid_grid_indices.m/threads;
        int b=(t+1)*valid_grid_indices.m/threads;
        typename TV::SPIN result_local;
        for(int k=a;k<b;k++){
            int i=valid_grid_indices(k);
            TV X=location.array(i);
            if(use_midpoint && u0) X+=dt/2*u0->array(i);
            result_local+=mass.array(i)*TV::Cross_Product(X,u.array(i));}
#pragma omp critical
        result+=result_local;}
    return result;
}
//#####################################################################
// Function Total_Grid_Kinetic_Energy
//#####################################################################
template<class TV> typename TV::SCALAR MPM_EXAMPLE<TV>::
Total_Grid_Kinetic_Energy(const ARRAY<TV,TV_INT>& u) const
{
    T result=0;
#pragma omp parallel for reduction(+:result)
    for(int i=0;i<valid_grid_indices.m;i++){
        int j=valid_grid_indices(i);
        result+=(T).5*mass.array(j)*u.array(j).Magnitude_Squared();}
    return result;
}
//#####################################################################
// Function Total_Particle_Kinetic_Energy
//#####################################################################
template<class TV> typename TV::SCALAR MPM_EXAMPLE<TV>::
Total_Particle_Kinetic_Energy() const
{
    T result=0;
#pragma omp parallel for reduction(+:result)
    for(int k=0;k<simulated_particles.m;k++){
        int p=simulated_particles(k);
        T result_local=particles.mass(p)/2*particles.V(p).Magnitude_Squared();
        if(particles.store_B){
            SYMMETRIC_MATRIX<T,TV::m> D=lag_Dp?Dp_inv(p).Inverse():Dp_inv(p);
            T a=(particles.B(p)*D).Double_Contract(particles.B(p));
            result_local+=particles.mass(p)/2*a;}
        result+=result_local;}
    return result;
}
//#####################################################################
// Function Average_Particle_Mass
//#####################################################################
template<class TV> typename TV::SCALAR MPM_EXAMPLE<TV>::
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
template<class TV> void MPM_EXAMPLE<TV>::
Add_Collision_Object(IMPLICIT_OBJECT<TV>* io,COLLISION_TYPE type,T friction,std::function<FRAME<TV>(T)> func_frame,std::function<TWIST<TV>(T)> func_twist)
{
    auto* co=new MPM_COLLISION_OBJECT<TV>(io,type,friction);
    co->func_frame=func_frame;
    co->func_twist=func_twist;
    collision_objects.Append(co);
}
//#####################################################################
// Function Update_Lagged_Forces
//#####################################################################
template<class TV> void MPM_EXAMPLE<TV>::
Update_Lagged_Forces(const T time) const
{
    for(int i=0;i<solid_body_collection.deformable_body_collection.deformables_forces.m;i++)
        if(LAGGED_FORCE<TV>* lf=dynamic_cast<LAGGED_FORCE<TV>*>(solid_body_collection.deformable_body_collection.deformables_forces(i)))
            lf->Lagged_Update_Position_Based_State(time);
}
//#####################################################################
// Function Add_Callbacks
//#####################################################################
template<class TV> void MPM_EXAMPLE<TV>::
Add_Callbacks(bool is_begin,const char* func_name,std::function<void()> func)
{
    time_step_callbacks.Get_Or_Insert(func_name).y(is_begin).Append(func);
}
namespace PhysBAM{
template class MPM_EXAMPLE<VECTOR<float,2> >;
template class MPM_EXAMPLE<VECTOR<float,3> >;
template class MPM_EXAMPLE<VECTOR<double,2> >;
template class MPM_EXAMPLE<VECTOR<double,3> >;
}
