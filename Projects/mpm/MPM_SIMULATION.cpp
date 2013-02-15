//#####################################################################
// Copyright 2013, Chenfanfu Jiang
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Vectors/VECTOR.h>
#include <PhysBAM_Tools/Matrices/MATRIX.h>
#include <PhysBAM_Tools/Matrices/DIAGONAL_MATRIX.h>
#include <PhysBAM_Tools/Matrices/SYMMETRIC_MATRIX.h>
#include <PhysBAM_Tools/Krylov_Solvers/CONJUGATE_GRADIENT.h>
#include <PhysBAM_Tools/Utilities/DEBUG_CAST.h>
#include <omp.h>
#include "TIMING.h"
#include "MPM_SYSTEM.h"
#include "MPM_VECTOR.h"
#include "MPM_SIMULATION.h"
namespace PhysBAM{
using ::std::exp;
//#####################################################################
// Constructor
//#####################################################################
template<class TV> MPM_SIMULATION<TV>::
MPM_SIMULATION()
    :min_mass(1e-5),min_pho(1e-5)
{}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> MPM_SIMULATION<TV>::
~MPM_SIMULATION()
{}
//#####################################################################
// Function Initialize
//#####################################################################
template<class TV> void MPM_SIMULATION<TV>::
Initialize()
{
// #pragma omp parallel for    
//     for(int t=0;t<1000;t++){
//         LOG::cout<<omp_get_thread_num()<<" ";}
//     LOG::cout<<std::endl;
    LOG::cout<<"Allocating Memories for Simulation..."<<std::endl;
    TIMING_START;
    gravity_constant=TV();gravity_constant(1)=-(T)9.8;
    N_particles=particles.X.m;
    mu.Resize(N_particles);
    lambda.Resize(N_particles);
    Je.Resize(N_particles);
    Jp.Resize(N_particles);
    Ue.Resize(N_particles);
    Ve.Resize(N_particles);
    Re.Resize(N_particles);
    Se.Resize(N_particles);
    SIGMAe.Resize(N_particles);
    influence_corner.Resize(N_particles);
    weight.Resize(N_particles);
    grad_weight.Resize(N_particles);
    for(int p=0;p<N_particles;p++){
        weight(p).Resize(RANGE<TV_INT>(TV_INT(),TV_INT()+IN));
        grad_weight(p).Resize(RANGE<TV_INT>(TV_INT(),TV_INT()+IN));}
    node_mass.Resize(RANGE<TV_INT>(TV_INT(),grid.counts));
    node_V.Resize(RANGE<TV_INT>(TV_INT(),grid.counts));
    node_V_star.Resize(RANGE<TV_INT>(TV_INT(),grid.counts));
    node_V_old.Resize(RANGE<TV_INT>(TV_INT(),grid.counts));
    node_force.Resize(RANGE<TV_INT>(TV_INT(),grid.counts));
    node_external_force.Resize(RANGE<TV_INT>(TV_INT(),grid.counts));
    frame=0;
    if(PROFILING) TIMING_END("");
}
//#####################################################################
// Function Advance_One_Time_Step_Forward_Euler
//#####################################################################
template<class TV> void MPM_SIMULATION<TV>::
Advance_One_Time_Step_Forward_Euler()
{
    Build_Weights_And_Grad_Weights();
    Build_Helper_Structures_For_Constitutive_Model();
    Rasterize_Particle_Data_To_The_Grid();
    if(frame==0) Compute_Particle_Volumes_And_Densities();
    Compute_Grid_Forces();
    if(use_gravity) Apply_Gravity_To_Grid_Forces();
    Update_Velocities_On_Grid();
    Grid_Based_Body_Collisions();
    node_V_old=node_V;node_V=node_V_star;
    Update_Deformation_Gradient();
    Update_Particle_Velocities();
    Particle_Based_Body_Collisions();
    Update_Particle_Positions();
    frame++;
}
//#####################################################################
// Function Advance_One_Time_Step_Backward_Euler
//#####################################################################
template<class TV> void MPM_SIMULATION<TV>::
Advance_One_Time_Step_Backward_Euler()
{
    Build_Weights_And_Grad_Weights();
    Build_Helper_Structures_For_Constitutive_Model();
    Rasterize_Particle_Data_To_The_Grid();
    if(frame==0) Compute_Particle_Volumes_And_Densities();
    Compute_Grid_Forces();
    if(use_gravity) Apply_Gravity_To_Grid_Forces();
    Update_Velocities_On_Grid();
    Grid_Based_Body_Collisions();
    Solve_The_Linear_System();
    Update_Deformation_Gradient();
    Update_Particle_Velocities();
    Particle_Based_Body_Collisions();
    Update_Particle_Positions();
    frame++;
}
//#####################################################################
// Function Build_Weights_And_Grad_Weights
//#####################################################################
template<class TV> void MPM_SIMULATION<TV>::
Build_Weights_And_Grad_Weights()
{
    TIMING_START;
    for(int p=0;p<N_particles;p++){
        grid_basis_function.Build_Weights_And_Grad_Weights_Exact(particles.X(p),grid,influence_corner(p),weight(p),grad_weight(p));}
    if(PROFILING) TIMING_END("Build_Weights_And_Grad_Weights");
}
//#####################################################################
// Function Build_Helper_Structures_For_Constitutive_Model
//#####################################################################
template<class TV> void MPM_SIMULATION<TV>::
Build_Helper_Structures_For_Constitutive_Model()
{
    TIMING_START;
    for(int p=0;p<N_particles;p++){
        constitutive_model.Compute_Helper_Quantities_Using_F(particles.Fe(p),particles.Fp(p),Je(p),Jp(p),Ue(p),SIGMAe(p),Ve(p),Re(p),Se(p));
        T lame_scale=exp(xi*((T)1-Jp(p)));
        mu(p)=mu0*lame_scale;
        lambda(p)=lambda0*lame_scale;}
    if(PROFILING) TIMING_END("Build_Helper_Structures_For_Constitutive_Model");
}
//#####################################################################
// Function Rasterize_Particle_Data_To_The_Grid
//#####################################################################
template<class TV> void MPM_SIMULATION<TV>::
Rasterize_Particle_Data_To_The_Grid()
{
    TIMING_START;
    TV_INT TV_INT_IN=TV_INT()+IN;
    RANGE<TV_INT> range(TV_INT(),TV_INT_IN);
    node_mass.Fill(T(0));
    for(int p=0;p<N_particles;p++){
        for(RANGE_ITERATOR<TV::m> it(range);it.Valid();it.Next())
            node_mass(influence_corner(p)+it.index)+=particles.mass(p)*weight(p)(it.index);}
    node_V.Fill(TV());
    for(int p=0;p<N_particles;p++){
        for(RANGE_ITERATOR<TV::m> it(range);it.Valid();it.Next()){
            TV_INT ind=influence_corner(p)+it.index;
            if(node_mass(ind)>min_mass) node_V(ind)+=particles.V(p)*particles.mass(p)*weight(p)(it.index)/node_mass(ind);}}
    if(PROFILING) TIMING_END("Rasterize_Particle_Data_To_The_Grid");
}
//#####################################################################
// Function Compute_Particle_Volumes_And_Densities
//#####################################################################
template<class TV> void MPM_SIMULATION<TV>::
Compute_Particle_Volumes_And_Densities()
{
    TIMING_START;
    TV_INT TV_INT_IN=TV_INT()+IN;
    T one_over_cell_volume=grid.one_over_dX.Product();
    RANGE<TV_INT> range(TV_INT(),TV_INT_IN);
    particles.density.Fill(T(0));
    particles.volume.Fill(T(0));
    for(int p=0;p<N_particles;p++){
        for(RANGE_ITERATOR<TV::m> it(range);it.Valid();it.Next()){
            TV_INT ind=influence_corner(p)+it.index;
            particles.density(p)+=node_mass(ind)*weight(p)(it.index)*one_over_cell_volume;}
        if(particles.density(p)>min_pho) particles.volume(p)=particles.mass(p)/particles.density(p);}
    if(PROFILING) TIMING_END("Compute_Particle_Volumes_And_Densities");
}
//#####################################################################
// Function Compute_Grid_Forces
//#####################################################################
template<class TV> void MPM_SIMULATION<TV>::
Compute_Grid_Forces()
{
    TIMING_START;
    TV_INT TV_INT_IN=TV_INT()+IN;
    RANGE<TV_INT> range(TV_INT(),TV_INT_IN);
    node_force.Fill(TV());
    for(int p=0;p<N_particles;p++){
        MATRIX<T,TV::m> B=particles.volume(p)*constitutive_model.Compute_dPsi_dFe(mu(p),lambda(p),particles.Fe(p),Re(p),Je(p))*(particles.Fe(p).Transposed());
        for(RANGE_ITERATOR<TV::m> it(range);it.Valid();it.Next()){
            TV_INT ind=influence_corner(p)+it.index;
            node_force(ind)-=B*grad_weight(p)(it.index);}}
    if(PROFILING) TIMING_END("Compute_Grid_Forces");
}
//#####################################################################
// Function Apply_Gravity_To_Grid_Forces
//#####################################################################
template<class TV> void MPM_SIMULATION<TV>::
Apply_Gravity_To_Grid_Forces()
{
    TIMING_START;
    TV_INT TV_INT_IN=TV_INT()+IN;
    RANGE<TV_INT> range(TV_INT(),TV_INT_IN);
    node_external_force.Fill(TV());
    for(int p=0;p<N_particles;p++){
        for(RANGE_ITERATOR<TV::m> it(range);it.Valid();it.Next()){
            TV_INT ind=influence_corner(p)+it.index;
            node_external_force(ind)=node_mass(ind)*gravity_constant;}}
    node_force+=node_external_force;
    if(PROFILING) TIMING_END("Apply_Gravity_To_Grid_Forces");
}
//#####################################################################
// Function Update_Velocities_On_Grid
//#####################################################################
template<class TV> void MPM_SIMULATION<TV>::
Update_Velocities_On_Grid()
{
    TIMING_START;
    static T eps=1e-5;
    for(int i=0;i<node_V.array.m;i++){
        node_V_star.array(i)=node_V.array(i);
        if(node_mass.array(i)>eps) node_V_star.array(i)+=dt/node_mass.array(i)*node_force.array(i);}
    if(PROFILING) TIMING_END("Update_Velocities_On_Grid");
}
//#####################################################################
// Function Grid_Based_Body_Collisions
//#####################################################################
template<class TV> void MPM_SIMULATION<TV>::
Grid_Based_Body_Collisions()
{
    TIMING_START;
    static T eps=1e-5;
    RANGE<TV_INT> range(TV_INT(),grid.counts);
    for(RANGE_ITERATOR<TV::m> it(range);it.Valid();it.Next()){
        if(grid.Node(it.index)(1)<=ground_level && node_V_star(it.index)(1)<=(T)0){
            TV vt(node_V_star(it.index));vt(1)=(T)0;
            T vn=node_V_star(it.index)(1);
            T vt_mag=vt.Magnitude();
            if(vt_mag>eps) node_V_star(it.index)=vt+friction_coefficient*vn*vt/vt_mag;
            else node_V_star(it.index)=vt;}}
    if(PROFILING) TIMING_END("Grid_Based_Body_Collisions");
}
//#####################################################################
// Function Solve_The_Linear_System
//#####################################################################
template<class TV> void MPM_SIMULATION<TV>::
Solve_The_Linear_System()
{
    TIMING_START;
    MPM_SYSTEM<TV> system(debug_cast<MPM_SIMULATION<TV>&>(*this));
    MPM_VECTOR<TV> rhs,x;
    ARRAY<KRYLOV_VECTOR_BASE<T>*> vectors;
    rhs.v.Resize(RANGE<TV_INT>(TV_INT(),grid.counts));
    x.v.Resize(RANGE<TV_INT>(TV_INT(),grid.counts));
    KRYLOV_SOLVER<T>::Ensure_Size(vectors,x,3);
    rhs.v=node_V_star;
    // system.Test_System(*vectors(0),*vectors(1),*vectors(2));
    CONJUGATE_GRADIENT<T> cg;
    KRYLOV_SOLVER<T>* solver=&cg;
    solver->print_residuals=true;
    solver->Solve(system,x,rhs,vectors,(T)1e-6,0,1000);
    node_V_old=node_V;
    node_V=x.v;
    if(PROFILING) TIMING_END("Solve_The_Linear_System");
}
//#####################################################################
// Function Update_Deformation_Gradient
//#####################################################################
template<class TV> void MPM_SIMULATION<TV>::
Update_Deformation_Gradient()
{
    TIMING_START;
    TV_INT TV_INT_IN=TV_INT()+IN;
    RANGE<TV_INT> range(TV_INT(),TV_INT_IN);
    if(!use_plasticity){
        for(int p=0;p<N_particles;p++){
            MATRIX<T,TV::m> grad_vp;
            for(RANGE_ITERATOR<TV::m> it(range);it.Valid();it.Next()){
                TV_INT ind=influence_corner(p)+it.index;
                grad_vp+=MATRIX<T,TV::m>::Outer_Product(node_V(ind),grad_weight(p)(it.index));}
            particles.Fe(p)=particles.Fe(p)+dt*grad_vp*particles.Fe(p);}}
    else{
        for(int p=0;p<N_particles;p++){
            MATRIX<T,TV::m> grad_vp;
            for(RANGE_ITERATOR<TV::m> it(range);it.Valid();it.Next()){
                TV_INT ind=influence_corner(p)+it.index;
                grad_vp+=MATRIX<T,TV::m>::Outer_Product(node_V(ind),grad_weight(p)(it.index));}
            MATRIX<T,TV::m> Fe_hat=particles.Fe(p)+dt*grad_vp*particles.Fe(p);
            MATRIX<T,TV::m> Fp_hat=particles.Fp(p);
            MATRIX<T,TV::m> F=Fe_hat*Fp_hat;
            MATRIX<T,TV::m> U_hat,V_hat;
            DIAGONAL_MATRIX<T,TV::m> SIGMA_hat;
            Fe_hat.Fast_Singular_Value_Decomposition(U_hat,SIGMA_hat,V_hat);
            SIGMA_hat=SIGMA_hat.Clamp_Min((T)1-theta_c);
            SIGMA_hat=SIGMA_hat.Clamp_Max((T)1+theta_s);
            particles.Fe(p)=U_hat*SIGMA_hat*(V_hat.Transposed());
            particles.Fp(p)=V_hat*(SIGMA_hat.Inverse())*(U_hat.Transposed())*F;}}
    if(PROFILING) TIMING_END("Update_Deformation_Gradient")
}
//#####################################################################
// Function Update_Particle_Velocities
//#####################################################################
template<class TV> void MPM_SIMULATION<TV>::
Update_Particle_Velocities()
{
    TIMING_START;
    TV_INT TV_INT_IN=TV_INT()+IN;
    RANGE<TV_INT> range(TV_INT(),TV_INT_IN);
    for(int p=0;p<N_particles;p++){
        TV V_PIC;
        TV V_FLIP=particles.V(p);
        for(RANGE_ITERATOR<TV::m> it(range);it.Valid();it.Next()){
            TV_INT ind=influence_corner(p)+it.index;
            T w=weight(p)(it.index);
            V_PIC+=node_V(ind)*w;
            V_FLIP+=(node_V(ind)-node_V_old(ind))*w;}
        particles.V(p)=((T)1-FLIP_alpha)*V_PIC+FLIP_alpha*V_FLIP;}
    if(PROFILING) TIMING_END("Update_Particle_Velocities");
}
//#####################################################################
// Function Particle_Based_Body_Collisions
//#####################################################################
template<class TV> void MPM_SIMULATION<TV>::
Particle_Based_Body_Collisions()
{
    TIMING_START;
    static T eps=1e-5;
    for(int p=0;p<N_particles;p++){
        if(particles.X(p)(1)<=ground_level && particles.V(p)(1)<=(T)0){
            TV vt(particles.V(p));vt(1)=(T)0;
            T vn=particles.V(p)(1);
            T vt_mag=vt.Magnitude();
            if(vt_mag>eps) particles.V(p)=vt+friction_coefficient*vn*vt/vt_mag;
            else particles.V(p)=vt;}}
    if(PROFILING) TIMING_END("Particle_Based_Body_Collisions");
}
//#####################################################################
// Function Update_Particle_Positions
//#####################################################################
template<class TV> void MPM_SIMULATION<TV>::
Update_Particle_Positions()
{
    TIMING_START;
    for(int p=0;p<N_particles;p++) particles.X(p)+=dt*particles.V(p);
    if(PROFILING) TIMING_END("Update_Particle_Positions");
}
//#####################################################################
template class MPM_SIMULATION<VECTOR<float,2> >;
template class MPM_SIMULATION<VECTOR<float,3> >;
template class MPM_SIMULATION<VECTOR<double,2> >;
template class MPM_SIMULATION<VECTOR<double,3> >;
}
