//#####################################################################
// Copyright 2015, Craig Schroeder.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Core/Arrays_Nd/ARRAYS_ND.h>
#include <Core/Log/LOG.h>
#include <Core/Matrices/DIAGONAL_MATRIX.h>
#include <Core/Utilities/NONCOPYABLE.h>
#include <Deformables/Constitutive_Models/DIAGONALIZED_ISOTROPIC_STRESS_DERIVATIVE.h>
#include <Deformables/Constitutive_Models/ISOTROPIC_CONSTITUTIVE_MODEL.h>
#include <Hybrid_Methods/Examples_And_Drivers/MPM_PARTICLES.h>
#include <Hybrid_Methods/Forces/MPM_FORCE_HELPER.h>
#include <Hybrid_Methods/Forces/MPM_GRAVITY.h>
#include <Hybrid_Methods/Iterators/GATHER_SCATTER.h>
namespace PhysBAM{
//#####################################################################
// Constructor
//#####################################################################
template<class TV> MPM_GRAVITY<TV>::
MPM_GRAVITY(const MPM_FORCE_HELPER<TV>& force_helper,const TV& gravity,
    GATHER_SCATTER<TV>& gather_scatter_input,ARRAY<int>* affected_particles)
    :BASE(force_helper),gravity(gravity),affect_all(!affected_particles),
    gather_scatter(affect_all?gather_scatter_input:*new GATHER_SCATTER<TV>(gather_scatter_input.grid,*new ARRAY<int>(*affected_particles)))
{
    if(!affect_all){
        gather_scatter.weights=gather_scatter_input.weights;
        gather_scatter.threads=gather_scatter_input.threads;
        gather_scatter.Prepare_Scatter(particles);}
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> MPM_GRAVITY<TV>:: 
~MPM_GRAVITY()
{
    if(!affect_all){
        delete &gather_scatter.simulated_particles;
        delete &gather_scatter;}
}
//#####################################################################
// Function Precompute
//#####################################################################
template<class TV> void MPM_GRAVITY<TV>:: 
Precompute(const T time,const T dt,bool want_dE,bool want_ddE)
{
}
//#####################################################################
// Function Potential_Energy
//#####################################################################
template<class TV> typename TV::SCALAR MPM_GRAVITY<TV>:: 
Potential_Energy(const T time) const
{
    T pe=0;
#pragma omp parallel for reduction(+:pe)
    for(int k=0;k<gather_scatter.simulated_particles.m;k++){
        int p=gather_scatter.simulated_particles(k);
        pe-=particles.mass(p)*particles.X(p).Dot(gravity);}
    return pe;
}
//#####################################################################
// Function Add_Forces
//#####################################################################
template<class TV> void MPM_GRAVITY<TV>:: 
Add_Forces(ARRAY<TV,TV_INT>& F,const T time) const
{
    gather_scatter.template Scatter<int>(false,0,
        [this,&F](int p,const PARTICLE_GRID_ITERATOR<TV>& it,int tid){
            F(it.Index())+=it.Weight()*particles.mass(p)*gravity;});
}
//#####################################################################
// Function Add_Hessian_Times
//#####################################################################
template<class TV> void MPM_GRAVITY<TV>:: 
Add_Hessian_Times(ARRAY<TV,TV_INT>& F,const ARRAY<TV,TV_INT>& Z,const T time) const
{
}
template class MPM_GRAVITY<VECTOR<float,1> >;
template class MPM_GRAVITY<VECTOR<float,2> >;
template class MPM_GRAVITY<VECTOR<float,3> >;
template class MPM_GRAVITY<VECTOR<double,1> >;
template class MPM_GRAVITY<VECTOR<double,2> >;
template class MPM_GRAVITY<VECTOR<double,3> >;
}
