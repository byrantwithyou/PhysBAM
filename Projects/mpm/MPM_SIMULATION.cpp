//#####################################################################
// Copyright 2013, Chenfanfu Jiang
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include "Fanfu_Utilities/FLATTEN_INDEX.h"
#include "MPM_SIMULATION.cpp"
namespace PhysBAM{
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
    //TODO: fill in all members in particles, constitutive_model, grid, dt
    //      Fe and Fp should be Identity

    N_particles=particles.X.m;
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
        weight(p).Resize(influenced_node_N);
        grad_weight(p).Resize(influenced_node_N);}
    N_nodes=grid.counts.Product();
    nodes_mass.Resize(N_nodes);
    nodes_V.Resize(N_nodes);
    nodes_force.Resize(N_nodes);

    frame=0;
}
//#####################################################################
// Function Advance_One_Time_Step_Forward_Euler
//#####################################################################
template<class TV> void MPM_SIMULATION<TV>::
Advance_One_Time_Step_Forward_Euler()
{
    Build_Weights_And_Grad_Weights();
    Rasterize_Particle_Data_To_Grid();
    if(frame==0) Compute_Particle_Volumes();
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
// Function Rasterize_Particle_Data_To_Grid
//#####################################################################
template<class TV> void MPM_SIMULATION<TV>::
Rasterize_Particle_Data_To_Grid()
{
    static T eps=1e-5;
    TV_INT TV_INT_IN=TV_INT()+IN;
    RANGE<TV_INT> range(TV_INT(),TV_INT_IN);
    node_mass.Fill(T(0));
    for(int p=0;p<N_particles;p++){
        for(RANGE_ITERATOR<TV::m> it(range);it.Valid();it.Next()){
            int ind=Flatten_Index(influence_corner(p)+it.index,grid.counts);
            node_mass(ind)+=particles.mass(p)*weight(p,Flatten_Index(it.index,TV_INT_IN));}}
    node_V.Fill(TV());
    for(int p=0;p<N_particles;p++){
        for(RANGE_ITERATOR<TV::m> it(range);it.Valid();it.Next()){
            int ind=Flatten_Index(influence_corner(p)+it.index,grid.counts);
            if(node_mass(ind)>eps) node_V(ind)+=particles.V(p)*particles.mass(p)*weight(p,Flatten_Index(it.index,TV_INT_IN))/node_mass(ind);}}
}
//#####################################################################
// Function Compute_Particle_Volumes
//#####################################################################
template<class TV> void MPM_SIMULATION<TV>::
Compute_Particle_Volumes()
{
    TV_INT TV_INT_IN=TV_INT()+IN;
    T one_over_cell_volume=grid.one_over_dX.Product();
    RANGE<TV_INT> range(TV_INT(),TV_INT_IN);
    particles.density.Fill(T(0));
    particles.volume.Fill(T(0));
    for(int p=0;p<N_particles;p++){
        for(RANGE_ITERATOR<TV::m> it(range);it.Valid();it.Next()){
            int ind=Flatten_Index(influence_corner(p)+it.index,grid.counts);
            particles.density(p)+=node_mass(ind)*weight(p,Flatten_Index(it.index,TV_INT_IN))*one_over_cell_volume;}
        particles.volume(p)=particles.mass(p)/particles.density(p);}
}
//#####################################################################
template class MPM_SIMULATION<VECTOR<float,1> >;
template class MPM_SIMULATION<VECTOR<float,2> >;
template class MPM_SIMULATION<VECTOR<float,3> >;
template class MPM_SIMULATION<VECTOR<double,1> >;
template class MPM_SIMULATION<VECTOR<double,2> >;
template class MPM_SIMULATION<VECTOR<double,3> >;
}
