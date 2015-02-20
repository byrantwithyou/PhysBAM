//#####################################################################
// Copyright 2015, Craig Schroeder.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Tools/Grids_Uniform_Arrays/ARRAYS_ND.h>
#include <Tools/Log/LOG.h>
#include <Tools/Matrices/DIAGONAL_MATRIX.h>
#include <Tools/Utilities/NONCOPYABLE.h>
#include <Tools/Utilities/PHYSBAM_OVERRIDE.h>
#include <Deformables/Constitutive_Models/DIAGONALIZED_ISOTROPIC_STRESS_DERIVATIVE.h>
#include <Deformables/Constitutive_Models/ISOTROPIC_CONSTITUTIVE_MODEL.h>
#include <Hybrid_Methods/Examples_And_Drivers/MPM_PARTICLES.h>
#include <Hybrid_Methods/Forces/MPM_FINITE_ELEMENTS.h>
#include <Hybrid_Methods/Iterators/GATHER_SCATTER.h>
namespace PhysBAM{
//#####################################################################
// Constructor
//#####################################################################
template<class TV> MPM_FINITE_ELEMENTS<TV>::
MPM_FINITE_ELEMENTS(MPM_PARTICLES<TV>& particles,ISOTROPIC_CONSTITUTIVE_MODEL<T,TV::m>& constitutive_model,
    GATHER_SCATTER<TV>& gather_scatter_input,ARRAY<int>* affected_particles)
    :BASE(particles),constitutive_model(constitutive_model),affect_all(!affected_particles),
    gather_scatter(affect_all?gather_scatter_input:*new GATHER_SCATTER<TV>(*new ARRAY<int>(*affected_particles)))
{
    if(!affect_all) gather_scatter.Set_Weights(gather_scatter_input.weights);
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> MPM_FINITE_ELEMENTS<TV>:: 
~MPM_FINITE_ELEMENTS()
{
    if(!affect_all){
        delete &gather_scatter.simulated_particles;
        delete &gather_scatter;}
    delete &constitutive_model;
}
//#####################################################################
// Function Precompute
//#####################################################################
template<class TV> void MPM_FINITE_ELEMENTS<TV>:: 
Precompute(const T time)
{
    U.Resize(particles.X.m);
    sigma.Resize(particles.X.m);
    dPi_dF.Resize(particles.X.m);
    for(int k=0;k<gather_scatter.simulated_particles.m;k++){
        int p=gather_scatter.simulated_particles(k);
        MATRIX<T,TV::m> V_local;
        particles.F(p).Fast_Singular_Value_Decomposition(U(p),sigma(p),V_local);
        constitutive_model.Isotropic_Stress_Derivative(sigma(p),dPi_dF(p),p);}
}
//#####################################################################
// Function Potential_Energy
//#####################################################################
template<class TV> typename TV::SCALAR MPM_FINITE_ELEMENTS<TV>:: 
Potential_Energy(const T time) const
{
    T pe=0;
    for(int k=0;k<gather_scatter.simulated_particles.m;k++){
        int p=gather_scatter.simulated_particles(k);
        pe+=particles.volume(p)*constitutive_model.Energy_Density(sigma(p),p);}
     LOG::printf("PE: %g\n",pe);
    return pe;
}
//#####################################################################
// Function Add_Forces
//#####################################################################
template<class TV> void MPM_FINITE_ELEMENTS<TV>:: 
Add_Forces(ARRAY<TV,TV_INT>& F,const T time) const
{
    ARRAY<MATRIX<T,TV::m> > V0_P_FT(gather_scatter.threads);
    gather_scatter.Scatter(
        [this,&V0_P_FT](int p,int tid)
        {V0_P_FT(tid)=U(p)*(constitutive_model.P_From_Strain(sigma(p),particles.volume(p),p)*sigma(p)).Times_Transpose(U(p));},
        [&V0_P_FT,&F](int p,const PARTICLE_GRID_ITERATOR<TV>& it,int tid)
        {F(it.Index())-=V0_P_FT(tid)*it.Gradient();},
        [](int p,int tid){},true);
}
//#####################################################################
// Function Add_Hessian_Times
//#####################################################################
template<class TV> void MPM_FINITE_ELEMENTS<TV>:: 
Add_Hessian_Times(ARRAY<TV,TV_INT>& F,const ARRAY<TV,TV_INT>& V,const T time) const
{
    tmp.Resize(particles.X.m);

    gather_scatter.Gather(
        [this](int p,int tid){tmp(p)=MATRIX<T,TV::m>();},
        [this,&V](int p,const PARTICLE_GRID_ITERATOR<TV>& it,int tid)
        {tmp(p)+=MATRIX<T,TV::m>::Outer_Product(it.Gradient(),V(it.Index()));},
        [this](int p,int tid){
            MATRIX<T,TV::m> dFh=(tmp(p)*U(p)).Transpose_Times(U(p))*sigma(p);
            MATRIX<T,TV::m> dPh=dPi_dF(p).Differential(dFh);
            tmp(p)=(U(p)*(sigma(p)*particles.volume(p))).Times_Transpose(U(p)*dPh);},true);

    gather_scatter.Scatter(
        [this,&F](int p,const PARTICLE_GRID_ITERATOR<TV>& it,int tid)
        {F(it.Index())+=tmp(p).Transpose_Times(it.Gradient());},true);
}
template class MPM_FINITE_ELEMENTS<VECTOR<float,2> >;
template class MPM_FINITE_ELEMENTS<VECTOR<float,3> >;
template class MPM_FINITE_ELEMENTS<VECTOR<double,2> >;
template class MPM_FINITE_ELEMENTS<VECTOR<double,3> >;
}
