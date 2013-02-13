//#####################################################################
// Copyright 2013, Chenfanfu Jiang
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Vectors/VECTOR.h>
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
Initialize(/*to be determined*/)
{
    //TODO: fill in all members in particles, grid, dt, mu0, lambda0, xi
    //      Fe and Fp should be Identity
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
    node_V_old.Resize(N_nodes);
    node_force.Resize(N_nodes);
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
    Update_Grid_Velocity_Forward_Euler_Time_Integration();
}
//#####################################################################
// Function Build_Weights_And_Grad_Weights
//#####################################################################
template<class TV> void MPM_SIMULATION<TV>::
Build_Weights_And_Grad_Weights()
{
    for(int p=0;p<N_particles;p++) grid_basis_function.Build_Weights_And_Grad_Weights_Exact(particles.X(p),grid,influence_corner(p),weight(p),grad_weight(p));
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
    TV_INT TV_INT_IN=TV_INT()+IN;
    T one_over_cell_volume=grid.one_over_dX.Product();
    RANGE<TV_INT> range(TV_INT(),TV_INT_IN);
    particles.density.Fill(T(0));
    particles.volume.Fill(T(0));
    for(int p=0;p<N_particles;p++){
        for(RANGE_ITERATOR<TV::m> it(range);it.Valid();it.Next()){
            int ind=Flatten_Index(influence_corner(p)+it.index,grid.counts);
            particles.density(p)+=node_mass(ind)*weight(p)(Flatten_Index(it.index,TV_INT_IN))*one_over_cell_volume;}
        particles.volume(p)=particles.mass(p)/particles.density(p);}
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
// Function Update_Grid_Velocity_Forward_Euler_Time_Integration
//#####################################################################
template<class TV> void MPM_SIMULATION<TV>::
Update_Grid_Velocity_Forward_Euler_Time_Integration()
{
    static T eps=1e-5;
    for(int i=0;i<N_nodes;i++){
        node_V_old(i)=node_V(i);
        if(node_mass(i)>eps) node_V(i)+=dt/node_mass(i)*node_force(i);}
}
//#####################################################################
template class MPM_SIMULATION<VECTOR<float,2> >;
template class MPM_SIMULATION<VECTOR<float,3> >;
template class MPM_SIMULATION<VECTOR<double,2> >;
template class MPM_SIMULATION<VECTOR<double,3> >;
}
