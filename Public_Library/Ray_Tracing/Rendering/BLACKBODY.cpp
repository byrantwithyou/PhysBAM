//#####################################################################
// Copyright 2002, Ronald Fedkiw, Geoffrey Irving.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//##################################################################### 
#include <Core/Math_Tools/constants.h>
#include <Core/Math_Tools/cube.h>
#include <Core/Math_Tools/sqr.h>
#include <Ray_Tracing/Rendering/BLACKBODY.h>
#include <cmath>
using namespace PhysBAM;
//#####################################################################
// Function Calculate_Blackbody_Spectrum
//#####################################################################
// radiance_spectrum is in units of Watts/(steradian*m^2)/m 
// input grid of wavelengths in nanometers
template<class T> void BLACKBODY<T>::
Calculate_Radiance_Spectrum(const T temperature,const GRID<VECTOR<T,1> >& grid,ARRAY<T,VECTOR<int,1> >& radiance_spectrum) const
{
    T constant_1=T(plancks_constant*sqr(speed_of_light)/2),constant_2=T(plancks_constant*speed_of_light/boltzmanns_constant);
    for(int i=0;i<grid.counts.x;i++){
        T lambda=grid.X(VECTOR<int,1>(i)).x;
        radiance_spectrum(i)=constant_1/(cube(lambda)*sqr(lambda)*(exp(constant_2/(lambda*temperature))-1));}
}
//#####################################################################
// Function Calculate_XYZ
//#####################################################################
// assumes f is defined on the grid
template<class T> VECTOR<T,3> BLACKBODY<T>::
Calculate_XYZ(const T temperature) const
{
    ARRAY<T,VECTOR<int,1> > radiance_spectrum(cie.grid.counts);
    Calculate_Radiance_Spectrum(temperature,cie.grid,radiance_spectrum);
    return cie.Calculate_XYZ(radiance_spectrum);
}
//#####################################################################
namespace PhysBAM{
template class BLACKBODY<float>;
template class BLACKBODY<double>;
}
