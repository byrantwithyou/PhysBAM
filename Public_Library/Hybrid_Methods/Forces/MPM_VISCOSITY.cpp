//#####################################################################
// Copyright 2015, Craig Schroeder.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Tools/Grids_Uniform_Arrays/ARRAYS_ND.h>
#include <Tools/Log/LOG.h>
#include <Tools/Matrices/DIAGONAL_MATRIX.h>
#include <Tools/Utilities/NONCOPYABLE.h>
#include <Deformables/Constitutive_Models/DIAGONALIZED_ISOTROPIC_STRESS_DERIVATIVE.h>
#include <Deformables/Constitutive_Models/ISOTROPIC_CONSTITUTIVE_MODEL.h>
#include <Hybrid_Methods/Examples_And_Drivers/MPM_PARTICLES.h>
#include <Hybrid_Methods/Forces/MPM_FORCE_HELPER.h>
#include <Hybrid_Methods/Forces/MPM_VISCOSITY.h>
#include <Hybrid_Methods/Iterators/GATHER_SCATTER.h>
namespace PhysBAM{
//#####################################################################
// Constructor
//#####################################################################
template<class TV> MPM_VISCOSITY<TV>::
MPM_VISCOSITY(MPM_FORCE_HELPER<TV>& force_helper,GATHER_SCATTER<TV>& gather_scatter_input,
    ARRAY<int>* affected_particles,const T constant_viscosity)
    :BASE(force_helper),affect_all(!affected_particles),
    gather_scatter(affect_all?gather_scatter_input:*new GATHER_SCATTER<TV>(gather_scatter_input.grid,*new ARRAY<int>(*affected_particles))),
    stored_dt(0),constant_viscosity(constant_viscosity)
{
    if(!affect_all){
        gather_scatter.weights=gather_scatter_input.weights;
        gather_scatter.threads=gather_scatter_input.threads;
        gather_scatter.Prepare_Scatter(particles);}
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> MPM_VISCOSITY<TV>:: 
~MPM_VISCOSITY()
{
    if(!affect_all){
        delete &gather_scatter.simulated_particles;
        delete &gather_scatter;}
}
//#####################################################################
// Function Precompute
//#####################################################################
template<class TV> void MPM_VISCOSITY<TV>:: 
Precompute(const T time,const T dt,bool want_dE,bool want_ddE)
{
    stored_dt=dt;
}
//#####################################################################
// Function Potential_Energy
//#####################################################################
template<class TV> typename TV::SCALAR MPM_VISCOSITY<TV>:: 
Potential_Energy(const T time) const
{
    T pe=0;
#pragma omp parallel for reduction(+:pe)
    for(int k=0;k<gather_scatter.simulated_particles.m;k++){
        int p=gather_scatter.simulated_particles(k);
        T c=(viscosity.m?viscosity(p):constant_viscosity)/sqr(stored_dt);
        MATRIX<T,TV::m> B=force_helper.B(p);
        SYMMETRIC_MATRIX<T,TV::m> S=B.Twice_Symmetric_Part();
        pe+=particles.volume(p)*force_helper.Fn(p).Determinant()/2*c*S.Frobenius_Norm_Squared();}
    return pe;
}
//#####################################################################
// Function Add_Forces
//#####################################################################
template<class TV> void MPM_VISCOSITY<TV>:: 
Add_Forces(ARRAY<TV,TV_INT>& F,const T time) const
{
    gather_scatter.template Scatter<MATRIX<T,TV::m> >(true,
        [this](int p,MATRIX<T,TV::m>& A)
        {
            T c=(viscosity.m?viscosity(p):constant_viscosity)/sqr(stored_dt);
            MATRIX<T,TV::m> B=force_helper.B(p);
            SYMMETRIC_MATRIX<T,TV::m> S=B.Twice_Symmetric_Part();
            A=2*c*particles.volume(p)*force_helper.Fn(p).Determinant()*S;
        },
        [this,&F](int p,const PARTICLE_GRID_ITERATOR<TV>& it,const MATRIX<T,TV::m>& A)
        {F(it.Index())-=A*it.Gradient();});
}
//#####################################################################
// Function Add_Hessian_Times
//#####################################################################
template<class TV> void MPM_VISCOSITY<TV>:: 
Add_Hessian_Times(ARRAY<TV,TV_INT>& F,const ARRAY<TV,TV_INT>& V,const T time) const
{
    tmp.Resize(particles.X.m);

    gather_scatter.template Gather<int>(true,
        [this](int p,int data){tmp(p)=MATRIX<T,TV::m>();},
        [this,&V](int p,const PARTICLE_GRID_ITERATOR<TV>& it,int data)
        {tmp(p)+=SYMMETRIC_MATRIX<T,TV::m>::Symmetric_Outer_Product(V(it.Index()),it.Gradient());},
        [this](int p,int data)
        {
            T c=(viscosity.m?viscosity(p):constant_viscosity)/sqr(stored_dt);
            tmp(p)*=2*c*force_helper.Fn(p).Determinant()*particles.volume(p);
        });
    gather_scatter.template Scatter<int>(true,0,
        [this,&F](int p,const PARTICLE_GRID_ITERATOR<TV>& it,int data)
        {F(it.Index())+=tmp(p)*it.Gradient();});
}
//#####################################################################
// Function Use_Variable_Viscosity
//#####################################################################
template<class TV> void MPM_VISCOSITY<TV>:: 
Use_Variable_Viscosity()
{
    const_cast<MPM_PARTICLES<TV>&>(particles).Add_Array(ATTRIBUTE_ID_VISCOSITY,&viscosity);
    viscosity.Fill(constant_viscosity);
}
template class MPM_VISCOSITY<VECTOR<float,2> >;
template class MPM_VISCOSITY<VECTOR<float,3> >;
template class MPM_VISCOSITY<VECTOR<double,2> >;
template class MPM_VISCOSITY<VECTOR<double,3> >;
}
