//#####################################################################
// Copyright 2015, Craig Schroeder, apt.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Tools/Log/LOG.h>
#include <Geometry/Geometry_Particles/DEBUG_PARTICLES.h>
#include <Geometry/Implicit_Objects/IMPLICIT_OBJECT.h>
#include <Geometry/Implicit_Objects/LEVELSET_IMPLICIT_OBJECT.h>
#include <Deformables/Deformable_Objects/DEFORMABLE_BODY_COLLECTION.h>
#include <Deformables/Forces/DEFORMABLES_FORCES.h>
#include <Deformables/Forces/LAGGED_FORCE.h>
#include <Hybrid_Methods/Collisions/MPM_COLLISION_IMPLICIT_OBJECT.h>
#include <Hybrid_Methods/Examples_And_Drivers/MPM_KKT_EXAMPLE.h>
#include <Hybrid_Methods/Examples_And_Drivers/MPM_PARTICLES.h>
#include <Hybrid_Methods/Forces/MPM_FORCE_HELPER.h>
#include <Hybrid_Methods/Forces/PARTICLE_GRID_FORCES.h>
#include <Hybrid_Methods/Iterators/GATHER_SCATTER.h>
#include <Hybrid_Methods/Iterators/PARTICLE_GRID_WEIGHTS.h>
#include <Hybrid_Methods/Iterators/PARTICLE_GRID_WEIGHTS_SPLINE.h>
#include <Hybrid_Methods/System/MPM_KRYLOV_VECTOR.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class TV> MPM_KKT_EXAMPLE<TV>::
MPM_KKT_EXAMPLE(const STREAM_TYPE stream_type)
    :stream_type(stream_type),particles(*new MPM_PARTICLES<TV>),
    deformable_body_collection(*new DEFORMABLE_BODY_COLLECTION<TV>(&particles,0)),
    debug_particles(*new DEBUG_PARTICLES<TV>),current_velocity(0),
    lagrangian_forces(deformable_body_collection.deformables_forces),
    weights(0),gather_scatter(*new GATHER_SCATTER<TV>(grid,simulated_particles)),
    gather_scatter_coarse(*new GATHER_SCATTER<TV>(coarse_grid,simulated_particles)),
    force_helper(*new MPM_FORCE_HELPER<TV>(particles,quad_F_coeff)),
    initial_time(0),last_frame(100),
    write_substeps_level(-1),substeps_delay_frame(-1),output_directory("output"),data_directory("../../Public_Data"),
    mass_contour(-1),restart(0),dt(0),time(0),frame_dt((T)1/24),min_dt(0),max_dt(frame_dt),ghost(5),
    use_reduced_rasterization(false),use_affine(false),use_f2p(false),use_midpoint(false),use_symplectic_euler(false),
    use_particle_collision(false),use_early_gradient_transfer(false),use_oldroyd(false),print_stats(false),flip(0),
    cfl(1),inv_Wi(0),newton_tolerance(1),newton_iterations(100),solver_tolerance(.5),solver_iterations(1000),
    test_diff(false),threads(1),last_te(0),last_grid_ke(0),output_structures_each_frame(false),
    quad_F_coeff(0),use_FEM_mass(0)
{
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> MPM_KKT_EXAMPLE<TV>::
~MPM_KKT_EXAMPLE()
{
    delete &deformable_body_collection;
    delete &particles;
    delete &debug_particles;
    delete weights;
    delete &gather_scatter;
    delete &gather_scatter_coarse;
    delete &force_helper;
    collision_objects.Delete_Pointers_And_Clean_Memory();
    fluid_walls.Delete_Pointers_And_Clean_Memory();
    forces.Delete_Pointers_And_Clean_Memory();
    av.Delete_Pointers_And_Clean_Memory();
}
//#####################################################################
// Function Write_Output_Files
//#####################################################################
template<class TV> void MPM_KKT_EXAMPLE<TV>::
Write_Output_Files(const int frame)
{
    std::string f=LOG::sprintf("%d",frame);

    FILE_UTILITIES::Write_To_File(stream_type,output_directory+"/common/grid",grid);
    FILE_UTILITIES::Write_To_File(stream_type,LOG::sprintf("%s/%d/mpm_particles",output_directory.c_str(),frame),particles);
    FILE_UTILITIES::Write_To_File(stream_type,LOG::sprintf("%s/%d/centered_velocities",output_directory.c_str(),frame),*current_velocity);
    FILE_UTILITIES::Write_To_File(stream_type,LOG::sprintf("%s/%d/density",output_directory.c_str(),frame),mass);
    FILE_UTILITIES::Write_To_File(stream_type,LOG::sprintf("%s/%d/restart_data",output_directory.c_str(),frame),time);
    int static_frame=output_structures_each_frame?frame:-1;
    bool write_structures=(frame==0 || output_structures_each_frame);
    deformable_body_collection.Write(stream_type,output_directory,output_directory,frame,static_frame,write_structures,false);

    ARRAY_VIEW<VECTOR<T,3> >* color_attribute=particles.template Get_Array<VECTOR<T,3> >(ATTRIBUTE_ID_COLOR);
    for(int i=0;i<particles.X.m;i++){
        Add_Debug_Particle(particles.X(i),particles.valid(i)?(*color_attribute)(i):VECTOR<T,3>(1,0,1));
        Debug_Particle_Set_Attribute<TV>(ATTRIBUTE_ID_V,particles.V(i));}
    for(int i=0;i<valid_pressure_indices.m;i++){
        int id=valid_pressure_indices(i);
        Add_Debug_Particle(coarse_location.array(id),VECTOR<T,3>(0,0,1));}
    for(int i=0;i<valid_velocity_indices.m;i++){
        int id=valid_velocity_indices(i);
        Add_Debug_Particle(location.array(id),VECTOR<T,3>(1,0,0));
        //Debug_Particle_Set_Attribute<VECTOR<T,1>>(ATTRIBUTE_ID_EFFECTIVE_MASS,VECTOR<T,1>(mass.array(id)));
        Debug_Particle_Set_Attribute<TV>(ATTRIBUTE_ID_V,velocity.array(id));}
    GRID<TV> ghost_grid(grid.numbers_of_cells+2*ghost,grid.Ghost_Domain(ghost),true);
    for(int i=0;i<collision_objects.m;i++)
        if(IMPLICIT_OBJECT<TV>* io=collision_objects(i)->Get_Implicit_Object(time))
            Dump_Levelset(ghost_grid,*io,VECTOR<T,3>(0.7,0.3,0.3));
    if(mass_contour>=0)
        Dump_Levelset(grid,mass,VECTOR<T,3>(0.2,0.6,0.2),mass_contour*Average_Particle_Mass());
    debug_particles.Write_Debug_Particles(stream_type,output_directory,frame);
}
//#####################################################################
// Function Read_Output_Files
//#####################################################################
template<class TV> void MPM_KKT_EXAMPLE<TV>::
Read_Output_Files(const int frame)
{
    std::string f=LOG::sprintf("%d",frame);
    FILE_UTILITIES::Read_From_File(stream_type,LOG::sprintf("%s/%d/mpm_particles",output_directory.c_str(),frame),particles);
    FILE_UTILITIES::Read_From_File(stream_type,LOG::sprintf("%s/%d/restart_data",output_directory.c_str(),frame),time);
}
//#####################################################################
// Function Capture_Stress
//#####################################################################
template<class TV> void MPM_KKT_EXAMPLE<TV>::
Capture_Stress()
{
    force_helper.Fn=particles.F;
}
//#####################################################################
// Function Precompute_Forces
//#####################################################################
template<class TV> void MPM_KKT_EXAMPLE<TV>::
Precompute_Forces(const T time,const T dt,const bool update_hessian)
{
    for(int i=0;i<forces.m;i++)
        forces(i)->Precompute(time,dt,true,update_hessian);
    deformable_body_collection.Update_Position_Based_State(time,false,update_hessian);
}
//#####################################################################
// Function Potential_Energy
//#####################################################################
template<class TV> typename TV::SCALAR MPM_KKT_EXAMPLE<TV>::
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
template<class TV> void MPM_KKT_EXAMPLE<TV>::
Add_Forces(ARRAY<TV,TV_INT>& F,const T time) const
{
    for(int i=0;i<forces.m;i++)
        forces(i)->Add_Forces(F,time);

    if(!lagrangian_forces.m) return;
    lagrangian_forces_F.Remove_All();
    lagrangian_forces_F.Resize(particles.X.m);
    deformable_body_collection.Add_Velocity_Independent_Forces(lagrangian_forces_F,time);
    gather_scatter.template Scatter<int>(
        [this,&F](int p,const PARTICLE_GRID_ITERATOR<TV>& it,int tid){
            F(it.Index())+=it.Weight()*lagrangian_forces_F(p);},false);
}
//#####################################################################
// Function Add_Hessian_Times
//#####################################################################
template<class TV> void MPM_KKT_EXAMPLE<TV>::
Add_Hessian_Times(ARRAY<TV,TV_INT>& F,const ARRAY<TV,TV_INT>& V,const T time) const
{
    for(int i=0;i<forces.m;i++)
        forces(i)->Add_Hessian_Times(F,V,time);

    if(!lagrangian_forces.m) return;
    lagrangian_forces_F.Remove_All();
    lagrangian_forces_V.Remove_All();
    lagrangian_forces_F.Resize(particles.X.m);
    lagrangian_forces_V.Resize(particles.X.m);
    gather_scatter.template Gather<int>(
        [this,&V](int p,const PARTICLE_GRID_ITERATOR<TV>& it,int tid){
            lagrangian_forces_V(p)+=it.Weight()*V(it.Index());},false);
    deformable_body_collection.Add_Implicit_Velocity_Independent_Forces(lagrangian_forces_V,lagrangian_forces_F,time);
    gather_scatter.template Scatter<int>(
        [this,&F](int p,const PARTICLE_GRID_ITERATOR<TV>& it,int tid){
            F(it.Index())-=it.Weight()*lagrangian_forces_F(p);},false);
}
//#####################################################################
// Function Add_Force
//#####################################################################
template<class TV> int MPM_KKT_EXAMPLE<TV>::
Add_Force(PARTICLE_GRID_FORCES<TV>& force)
{
    return forces.Append(&force);
}
//#####################################################################
// Function Add_Force
//#####################################################################
template<class TV> int MPM_KKT_EXAMPLE<TV>::
Add_Force(DEFORMABLES_FORCES<TV>& force)
{
    return deformable_body_collection.Add_Force(&force);
}
//#####################################################################
// Function Set_Weights
//#####################################################################
template<class TV> void MPM_KKT_EXAMPLE<TV>::
Set_Weights(PARTICLE_GRID_WEIGHTS<TV>* weights_input)
{
    weights=weights_input;
    gather_scatter.weights=weights;
    if(weights->Order()==2)
        coarse_weights=new PARTICLE_GRID_WEIGHTS_SPLINE<TV,1>(coarse_grid,threads);
    else if(weights->Order()==3)
        coarse_weights=new PARTICLE_GRID_WEIGHTS_SPLINE<TV,2>(coarse_grid,threads);
    else PHYSBAM_FATAL_ERROR("Unrecognized interpolation order");
    gather_scatter_coarse.weights=coarse_weights;
}
//#####################################################################
// Function Total_Particle_Linear_Momentum
//#####################################################################
template<class TV> TV MPM_KKT_EXAMPLE<TV>::
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
template<class TV> TV MPM_KKT_EXAMPLE<TV>::
Total_Grid_Linear_Momentum(const ARRAY<TV,TV_INT>& u) const
{
    TV result;
#pragma omp parallel for
    for(int t=0;t<threads;t++){
        int a=t*valid_velocity_indices.m/threads;
        int b=(t+1)*valid_velocity_indices.m/threads;
        TV result_local;
        for(int k=a;k<b;k++){
            int j=valid_velocity_indices(k);
            result_local+=mass.array(j)*u.array(j);}
#pragma omp critical
        result+=result_local;}
    return result;
}
//#####################################################################
// Function Total_Grid_Angular_Momentum
//#####################################################################
template<class TV> typename TV::SPIN MPM_KKT_EXAMPLE<TV>::
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
template<class TV> typename TV::SPIN MPM_KKT_EXAMPLE<TV>::
Total_Grid_Angular_Momentum(T dt,const ARRAY<TV,TV_INT>& u,const ARRAY<TV,TV_INT>* u0) const
{
    typename TV::SPIN result;
#pragma omp parallel for
    for(int t=0;t<threads;t++){
        int a=t*valid_velocity_indices.m/threads;
        int b=(t+1)*valid_velocity_indices.m/threads;
        typename TV::SPIN result_local;
        for(int k=a;k<b;k++){
            int i=valid_velocity_indices(k);
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
template<class TV> typename TV::SCALAR MPM_KKT_EXAMPLE<TV>::
Total_Grid_Kinetic_Energy(const ARRAY<TV,TV_INT>& u) const
{
    T result=0;
#pragma omp parallel for reduction(+:result)
    for(int i=0;i<valid_velocity_indices.m;i++){
        int j=valid_velocity_indices(i);
        result+=(T).5*mass.array(j)*u.array(j).Magnitude_Squared();}
    return result;
}
//#####################################################################
// Function Total_Particle_Kinetic_Energy
//#####################################################################
template<class TV> typename TV::SCALAR MPM_KKT_EXAMPLE<TV>::
Total_Particle_Kinetic_Energy() const
{
    T result=0,Dp_inverse=0;
    if(use_affine && weights->constant_scalar_inertia_tensor)
        Dp_inverse=weights->Constant_Scalar_Inverse_Dp();
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
template<class TV> typename TV::SCALAR MPM_KKT_EXAMPLE<TV>::
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
template<class TV> void MPM_KKT_EXAMPLE<TV>::
Add_Collision_Object(IMPLICIT_OBJECT<TV>* io,COLLISION_TYPE type,T friction)
{
    collision_objects.Append(new MPM_COLLISION_IMPLICIT_OBJECT<TV>(io,type,friction));
}
//#####################################################################
// Function Add_Fluid_Wall
//#####################################################################
template<class TV> void MPM_KKT_EXAMPLE<TV>::
Add_Fluid_Wall(IMPLICIT_OBJECT<TV>* io)
{
    fluid_walls.Append(io);
}
//#####################################################################
// Function Update_Lagged_Forces
//#####################################################################
template<class TV> void MPM_KKT_EXAMPLE<TV>::
Update_Lagged_Forces(const T time) const
{
    for(int i=0;i<deformable_body_collection.deformables_forces.m;i++)
        if(LAGGED_FORCE<TV>* lf=dynamic_cast<LAGGED_FORCE<TV>*>(deformable_body_collection.deformables_forces(i)))
            lf->Lagged_Update_Position_Based_State(time);
}
//#####################################################################
namespace PhysBAM{
template class MPM_KKT_EXAMPLE<VECTOR<float,1> >;
template class MPM_KKT_EXAMPLE<VECTOR<float,2> >;
template class MPM_KKT_EXAMPLE<VECTOR<float,3> >;
template class MPM_KKT_EXAMPLE<VECTOR<double,1> >;
template class MPM_KKT_EXAMPLE<VECTOR<double,2> >;
template class MPM_KKT_EXAMPLE<VECTOR<double,3> >;
}
