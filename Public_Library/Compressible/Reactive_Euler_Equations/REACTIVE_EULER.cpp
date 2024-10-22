//#####################################################################
// Copyright 2002, Ronald Fedkiw.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Core/Math_Tools/sqr.h>
#include <Grid_Tools/Grids/GRID.h>
#include <Compressible/Reactive_Euler_Equations/REACTIVE_EULER.h>
#include <cmath>
namespace PhysBAM{
//#####################################################################
// Function e
//#####################################################################
// internal energy
template<class TV> typename TV::SCALAR REACTIVE_EULER<TV>::
e(const T rho,const T rho_u,const T E)
{
    return E/rho-sqr(rho_u/rho)/2;
}
//#####################################################################
// Function e
//#####################################################################
// internal energy
template<class TV> typename TV::SCALAR REACTIVE_EULER<TV>::
e(const T rho,const T rho_u,const T rho_v,const T E)
{
    return E/rho-(sqr(rho_u/rho)+sqr(rho_v/rho))/2;
}
//#####################################################################
// Function e
//#####################################################################
// internal energy - 3D
template<class TV> typename TV::SCALAR REACTIVE_EULER<TV>::
e(const T rho,const T rho_u,const T rho_v,const T rho_w,const T E)
{
    return E/rho-(sqr(rho_u/rho)+sqr(rho_v/rho)+sqr(rho_w/rho))/2;
}
//#####################################################################
template class REACTIVE_EULER<VECTOR<float,1> >;
template class REACTIVE_EULER<VECTOR<float,2> >;
template class REACTIVE_EULER<VECTOR<float,3> >;
template class REACTIVE_EULER<VECTOR<double,1> >;
template class REACTIVE_EULER<VECTOR<double,2> >;
template class REACTIVE_EULER<VECTOR<double,3> >;
}
