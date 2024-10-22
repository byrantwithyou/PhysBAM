//#####################################################################
// Copyright 2002-2003, Robert Bridson, Ron Fedkiw, Eran Guendelman.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Core/Arrays_Nd/ARRAYS_ND.h>
#include <Core/Math_Tools/constants.h>
#include <Core/Math_Tools/sqr.h>
#include <Grid_Tools/Fourier_Transforms/FFT.h>
#include <Grid_Tools/Fourier_Transforms/FFT_UTILITIES.h>
#include <Grid_Tools/Grids/GRID.h>
#include <Dynamics/Incompressible_Flows/TURBULENCE.h>
using namespace PhysBAM;
//#####################################################################
// Function Generate_Random_Turbulence
//#####################################################################
template<class TV> void TURBULENCE<TV>::
Generate_Random_Turbulence(const GRID<TV>& grid,VECTOR<ARRAY<T,TV_INT>,TV::m>& u) const
{
    FFT<TV> fft;
    TV_INT size(grid.counts-1);
    size(TV::m-1)=grid.counts(TV::m-1)/2;
    RANGE<TV_INT> range(TV_INT(),size);
    VECTOR<ARRAY<std::complex<T>,TV_INT>,TV::m> u_hat;
    for(int i=0;i<TV::m;i++) u_hat(i).Resize(range);

    T coeff=(T)pi*(1<<(TV::m-1));
    TV coefficients=(T)(2*pi)/grid.domain.Edge_Lengths();
    for(RANGE_ITERATOR<TV::m> it(range);it.Valid();it.Next()){
        TV k;
        for(int i=0;i<TV::m;i++) k(i)=coefficients(i)*(it.index(i)<grid.counts(i)/2?i:i-grid.counts(i));
        T k2=k.Magnitude_Squared(),km=sqrt(k2);
        T area=coeff*k2;
        T energy=0;
        if(km > k_inertial) energy=constant*pow(epsilon,((T)2/3))*pow((T)km,-(T)5/3)/area;
        T sqrt_energy_over_two=sqrt((T).5*energy);
        for(int i=0;i<TV::m;i++){
            T r=random->Get_Gaussian(),theta=(T)pi*random->Get_Uniform_Number((T)0,(T)1);
            u_hat(i)(it.index)=std::polar(sqrt_energy_over_two*r,theta);}}

    for(int i=0;i<TV::m;i++) Enforce_Real_Valued_Symmetry(u_hat(i));
    if(incompressible) Make_Divergence_Free(u_hat,grid.domain.Edge_Lengths());
    for(int i=0;i<TV::m;i++) fft.Inverse_Transform(u_hat(i),u(i));

    // rescale the final velocity
    if(rescaled_average_velocity){
        T average_velocity=0;
        for(RANGE_ITERATOR<TV::m> it(grid.Domain_Indices());it.Valid();it.Next()){
            T s=0;
            for(int i=0;i<TV::m;i++) s+=sqr(u(i)(it.index));
            average_velocity+=sqrt(s);}
        average_velocity/=(grid.counts.Product());
        T scaling=rescaled_average_velocity/average_velocity;
        for(int i=0;i<TV::m;i++) u(i)*=scaling;}
}
//#####################################################################
namespace PhysBAM{
template class TURBULENCE<VECTOR<float,1> >;
template class TURBULENCE<VECTOR<float,2> >;
template class TURBULENCE<VECTOR<float,3> >;
template class TURBULENCE<VECTOR<double,1> >;
template class TURBULENCE<VECTOR<double,2> >;
template class TURBULENCE<VECTOR<double,3> >;
}
