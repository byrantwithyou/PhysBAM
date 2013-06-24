//#####################################################################
// Copyright 2002-2003, Robert Bridson, Ron Fedkiw, Eran Guendelman.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Tools/Fourier_Transforms/FFT_1D.h>
#include <Tools/Fourier_Transforms/FFT_2D.h>
#include <Tools/Fourier_Transforms/FFT_3D.h>
#include <Tools/Fourier_Transforms/FFT_POLICY.h>
#include <Tools/Fourier_Transforms_Calculations/TURBULENCE.h>
#include <Tools/Grids_Uniform/GRID.h>
#include <Tools/Grids_Uniform_Arrays/ARRAYS_ND.h>
#include <Tools/Math_Tools/constants.h>
#include <Tools/Math_Tools/sqr.h>
#include <Tools/Vectors/COMPLEX.h>
using namespace PhysBAM;
//#####################################################################
// Function Generate_Random_Turbulence
//#####################################################################
template<class TV> void TURBULENCE<TV>::
Generate_Random_Turbulence(const GRID<TV>& grid,VECTOR<ARRAY<T,TV_INT>,TV::m>& u) const
{
    typename FFT_POLICY<TV>::FFT fft(grid);
    TV_INT size(grid.counts-1);
    size(TV::m-1)=grid.counts(TV::m-1)/2;
    RANGE<TV_INT> range(TV_INT(),size);
    VECTOR<ARRAY<COMPLEX<T>,TV_INT>,TV::m> u_hat;
    for(int i=0;i<TV::m;i++) u_hat(i).Resize(range);

    T coeff=(T)pi*(1<<(TV::m-1));
    TV coefficients=(T)(2*pi)/grid.domain.Edge_Lengths();
    for(RANGE_ITERATOR<TV::m> it(range);it.Valid();it.Next()){
        TV k;
        for(int i=0;i<TV::m;i++) k(i)=coefficients(i)*(it.index(i)<grid.counts(i)/2?i:i-grid.counts(i));
        T k2=k.Magnitude_Squared(),km=sqrt(k2);
        T area=coeff*k2;
        T energy=0;
        if(km > k_inertial) energy=constant*pow(epsilon,(T)two_thirds)*pow((T)km,(T)-five_thirds)/area;
        T sqrt_energy_over_two=sqrt((T).5*energy);
        for(int i=0;i<TV::m;i++){
            T r=random->Get_Gaussian(),theta=(T)pi*random->Get_Uniform_Number((T)0,(T)1);
            u_hat(i)(it.index)=COMPLEX<T>::Polar(sqrt_energy_over_two*r,theta);}}

    for(int i=0;i<TV::m;i++) fft.Enforce_Real_Valued_Symmetry(u_hat(i));
    if(incompressible) fft.Make_Divergence_Free(u_hat);
    for(int i=0;i<TV::m;i++) fft.Inverse_Transform(u_hat(i),u(i),false,false);

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
