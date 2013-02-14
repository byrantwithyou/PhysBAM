//#####################################################################
// Copyright 2013, Chenfanfu Jiang
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Vectors/VECTOR.h>
#include <PhysBAM_Tools/Matrices/MATRIX.h>
#include <PhysBAM_Tools/Matrices/DIAGONAL_MATRIX.h>
#include <PhysBAM_Tools/Matrices/SYMMETRIC_MATRIX.h>
#include <omp.h>
#include "Fanfu_Utilities/FLATTEN_INDEX.h"
#include "MPM_SIMULATION.h"
namespace PhysBAM{
using ::std::exp;
//#####################################################################
// Constructor
//#####################################################################
template<class TV> MPM_SIMULATION<TV>::
MPM_SIMULATION()
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
//#pragma omp parallel for    
//    for(int t=0;t<1000000;t++){
//        LOG::cout<<t<<std::endl;}

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
    IN_POWER=1;
    for(int d=0;d<TV::m;d++) IN_POWER*=IN;
    for(int p=0;p<N_particles;p++){
        weight(p).Resize(IN_POWER);
        grad_weight(p).Resize(IN_POWER);}
    N_nodes=grid.counts.Product();
    node_mass.Resize(N_nodes);
    node_V.Resize(N_nodes);
    node_V_star.Resize(N_nodes);
    node_V_old.Resize(N_nodes);
    node_force.Resize(N_nodes);
    node_external_force.Resize(N_nodes);
    frame=0;
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
    for(int p=0;p<N_particles;p++){
        grid_basis_function.Build_Weights_And_Grad_Weights_Exact(particles.X(p),grid,influence_corner(p),weight(p),grad_weight(p));}
}
//#####################################################################
// Function Build_Helper_Structures_For_Constitutive_Model
//#####################################################################
template<class TV> void MPM_SIMULATION<TV>::
Build_Helper_Structures_For_Constitutive_Model()
{
    for(int p=0;p<N_particles;p++){
        constitutive_model.Compute_Helper_Quantities_Using_F(particles.Fe(p),particles.Fp(p),Je(p),Jp(p),Ue(p),SIGMAe(p),Ve(p),Re(p),Se(p));
        T lame_scale=exp(xi*((T)1-Jp(p)));
        mu(p)=mu0*lame_scale;
        lambda(p)=lambda0*lame_scale;}
}
//#####################################################################
// Function Rasterize_Particle_Data_To_The_Grid
//#####################################################################
template<class TV> void MPM_SIMULATION<TV>::
Rasterize_Particle_Data_To_The_Grid()
{
    static T eps=1e-5;
    TV_INT TV_INT_IN=TV_INT()+IN;
    RANGE<TV_INT> range(TV_INT(),TV_INT_IN);
    node_mass.Fill(T(0));
    for(int p=0;p<N_particles;p++){
        for(RANGE_ITERATOR<TV::m> it(range);it.Valid();it.Next()){
            int ind=Flatten_Index(influence_corner(p)+it.index,grid.counts);
            node_mass(ind)+=particles.mass(p)*weight(p)(Flatten_Index(it.index,TV_INT_IN));}}
    node_V.Fill(TV());
    for(int p=0;p<N_particles;p++){
        for(RANGE_ITERATOR<TV::m> it(range);it.Valid();it.Next()){
            int ind=Flatten_Index(influence_corner(p)+it.index,grid.counts);
            if(node_mass(ind)>eps) node_V(ind)+=particles.V(p)*particles.mass(p)*weight(p)(Flatten_Index(it.index,TV_INT_IN))/node_mass(ind);}}
}
//#####################################################################
// Function Compute_Particle_Volumes_And_Densities
//#####################################################################
template<class TV> void MPM_SIMULATION<TV>::
Compute_Particle_Volumes_And_Densities()
{
    static T eps=1e-5; 
    TV_INT TV_INT_IN=TV_INT()+IN;
    T one_over_cell_volume=grid.one_over_dX.Product();
    RANGE<TV_INT> range(TV_INT(),TV_INT_IN);
    particles.density.Fill(T(0));
    particles.volume.Fill(T(0));
#pragma omp parallel for
    for(int p=0;p<N_particles;p++){
        for(RANGE_ITERATOR<TV::m> it(range);it.Valid();it.Next()){
            int ind=Flatten_Index(influence_corner(p)+it.index,grid.counts);
            particles.density(p)+=node_mass(ind)*weight(p)(Flatten_Index(it.index,TV_INT_IN))*one_over_cell_volume;}
        if(particles.density(p)>eps) particles.volume(p)=particles.mass(p)/particles.density(p);}
}
//#####################################################################
// Function Compute_Grid_Forces
//#####################################################################
template<class TV> void MPM_SIMULATION<TV>::
Compute_Grid_Forces()
{
    TV_INT TV_INT_IN=TV_INT()+IN;
    RANGE<TV_INT> range(TV_INT(),TV_INT_IN);
    node_force.Fill(TV());
    for(int p=0;p<N_particles;p++){
        MATRIX<T,TV::m> B=particles.volume(p)*constitutive_model.Compute_dPsi_dFe(mu(p),lambda(p),particles.Fe(p),Re(p),Je(p))*(particles.Fe(p).Transposed());
        for(RANGE_ITERATOR<TV::m> it(range);it.Valid();it.Next()){
            int ind=Flatten_Index(influence_corner(p)+it.index,grid.counts);
            node_force(ind)-=B*grad_weight(p)(Flatten_Index(it.index,TV_INT_IN));}}
}
//#####################################################################
// Function Apply_Gravity_To_Grid_Forces
//#####################################################################
template<class TV> void MPM_SIMULATION<TV>::
Apply_Gravity_To_Grid_Forces()
{
    TV_INT TV_INT_IN=TV_INT()+IN;
    RANGE<TV_INT> range(TV_INT(),TV_INT_IN);
    node_external_force.Fill(TV());
    for(int p=0;p<N_particles;p++){
        for(RANGE_ITERATOR<TV::m> it(range);it.Valid();it.Next()){
            int ind=Flatten_Index(influence_corner(p)+it.index,grid.counts);
            node_external_force(ind)=node_mass(ind)*gravity_constant;}}
    node_force+=node_external_force;
}
//#####################################################################
// Function Update_Velocities_On_Grid
//#####################################################################
template<class TV> void MPM_SIMULATION<TV>::
Update_Velocities_On_Grid()
{
    static T eps=1e-5;
    for(int i=0;i<N_nodes;i++){
        node_V_star(i)=node_V(i);
        if(node_mass(i)>eps) node_V_star(i)+=dt/node_mass(i)*node_force(i);}
}
//#####################################################################
// Function Grid_Based_Body_Collisions
//#####################################################################
template<class TV> void MPM_SIMULATION<TV>::
Grid_Based_Body_Collisions()
{
    static T eps=1e-5;
    RANGE<TV_INT> range(TV_INT(),grid.counts);
    for(RANGE_ITERATOR<TV::m> it(range);it.Valid();it.Next()){
        int ind=Flatten_Index(it.index,grid.counts);
        if(grid.Node(it.index)(1)<=ground_level && node_V_star(ind)(1)<=(T)0){
            TV vt(node_V_star(ind));vt(1)=(T)0;
            T vn=node_V_star(ind)(1);
            T vt_mag=vt.Magnitude();
            if(vt_mag>eps) node_V_star(ind)=vt+friction_coefficient*vn*vt/vt_mag;
            else node_V_star(ind)=vt;}}
}
//#####################################################################
// Function Solve_The_Linear_System
//#####################################################################
template<class TV> void MPM_SIMULATION<TV>::
Solve_The_Linear_System()
{
    //TODO
}
//#####################################################################
// Function Update_Deformation_Gradient
//#####################################################################
template<class TV> void MPM_SIMULATION<TV>::
Update_Deformation_Gradient()
{
    TV_INT TV_INT_IN=TV_INT()+IN;
    RANGE<TV_INT> range(TV_INT(),TV_INT_IN);
    if(!use_plasticity){
        for(int p=0;p<N_particles;p++){
            MATRIX<T,TV::m> grad_vp;
            for(RANGE_ITERATOR<TV::m> it(range);it.Valid();it.Next()){
                int ind=Flatten_Index(influence_corner(p)+it.index,grid.counts);
                grad_vp+=MATRIX<T,TV::m>::Outer_Product(node_V(ind),grad_weight(p)(Flatten_Index(it.index,TV_INT_IN)));}
            particles.Fe(p)=particles.Fe(p)+dt*grad_vp*particles.Fe(p);}}
    else{
        for(int p=0;p<N_particles;p++){
            MATRIX<T,TV::m> grad_vp;
            for(RANGE_ITERATOR<TV::m> it(range);it.Valid();it.Next()){
                int ind=Flatten_Index(influence_corner(p)+it.index,grid.counts);
                grad_vp+=MATRIX<T,TV::m>::Outer_Product(node_V(ind),grad_weight(p)(Flatten_Index(it.index,TV_INT_IN)));}
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
}
//#####################################################################
// Function Update_Particle_Velocities
//#####################################################################
template<class TV> void MPM_SIMULATION<TV>::
Update_Particle_Velocities()
{
    TV_INT TV_INT_IN=TV_INT()+IN;
    RANGE<TV_INT> range(TV_INT(),TV_INT_IN);
    for(int p=0;p<N_particles;p++){
        TV V_PIC;
        TV V_FLIP=particles.V(p);
        for(RANGE_ITERATOR<TV::m> it(range);it.Valid();it.Next()){
            int ind=Flatten_Index(influence_corner(p)+it.index,grid.counts);
            T w=weight(p)(Flatten_Index(it.index,TV_INT_IN));
            V_PIC+=node_V(ind)*w;
            V_FLIP+=(node_V(ind)-node_V_old(ind))*w;}
        particles.V(p)=((T)1-FLIP_alpha)*V_PIC+FLIP_alpha*V_FLIP;}
}
//#####################################################################
// Function Particle_Based_Body_Collisions
//#####################################################################
template<class TV> void MPM_SIMULATION<TV>::
Particle_Based_Body_Collisions()
{
    static T eps=1e-5;
    for(int p=0;p<N_particles;p++){
        if(particles.X(p)(1)<=ground_level && particles.V(p)(1)<=(T)0){
            TV vt(particles.V(p));vt(1)=(T)0;
            T vn=particles.V(p)(1);
            T vt_mag=vt.Magnitude();
            if(vt_mag>eps) particles.V(p)=vt+friction_coefficient*vn*vt/vt_mag;
            else particles.V(p)=vt;}}
}
//#####################################################################
// Function Update_Particle_Positions
//#####################################################################
template<class TV> void MPM_SIMULATION<TV>::
Update_Particle_Positions()
{
    for(int p=0;p<N_particles;p++) particles.X(p)+=dt*particles.V(p);
}
//#####################################################################
// Function Compute_df
//#####################################################################
template<class TV> void MPM_SIMULATION<TV>::
Compute_df(const ARRAY<TV>& du,ARRAY<TV>& df)
{
    TV_INT TV_INT_IN=TV_INT()+IN;
    RANGE<TV_INT> range(TV_INT(),TV_INT_IN);
    df.Resize(N_nodes);
    df.Fill(TV());
    for(int p=0;p<N_particles;p++){
        MATRIX<T,TV::m> Cp;
        for(RANGE_ITERATOR<TV::m> it(range);it.Valid();it.Next()){
            int j=Flatten_Index(influence_corner(p)+it.index,grid.counts);
            Cp+=MATRIX<T,TV::m>::Outer_Product(du(j),grad_weight(p)(Flatten_Index(it.index,TV_INT_IN)));}
        MATRIX<T,TV::m> Ep=Cp*particles.Fe(p);
        MATRIX<T,TV::m> Ap=constitutive_model.Compute_d2Psi_dFe_dFe_Action_dF(mu(p),lambda(p),particles.Fe(p),Je(p),Re(p),Se(p),Ep);
        MATRIX<T,TV::m> Gp=particles.volume(p)*Ap*(particles.Fe(p).Transposed());
        for(RANGE_ITERATOR<TV::m> it(range);it.Valid();it.Next()){
            int i=Flatten_Index(influence_corner(p)+it.index,grid.counts);
            df(i)-=Gp*grad_weight(p)(Flatten_Index(it.index,TV_INT_IN));}}
}
//#####################################################################
template class MPM_SIMULATION<VECTOR<float,2> >;
template class MPM_SIMULATION<VECTOR<float,3> >;
template class MPM_SIMULATION<VECTOR<double,2> >;
template class MPM_SIMULATION<VECTOR<double,3> >;
}
