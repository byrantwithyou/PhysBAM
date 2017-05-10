//#####################################################################
// Copyright 2017, Ounan Ding, Craig Schroeder.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Core/Math_Tools/sqr.h>
namespace PhysBAM{

// From C. Macdonald and S. Ruuth, "Level set equations on surfaces via the Closest Point Method"
// Interpolate the data at location x.
// If x=0, then interpolates at z2; x=1 interpolates at z3
// Assumes x in [0,1].
// eps is used to prevent division by zero; the paper uses eps=1e-6.
template<class T> T
WENO_Interpolation(T x,T z0,T z1,T z2,T z3,T z4,T z5,T eps)
{
    T A0=z2,A1=((T)1/6)*(3*z2+2*z3-6*z1+z0),A2=(T).5*(-2*z2+z3+z1),A3=((T)1/6)*(z3-3*z2+3*z1-z0);
    T B0=z2,B1=((T)1/6)*(-3*z2+6*z3-2*z1-z4),B2=(T).5*(-2*z2+z3+z1),B3=((T)1/6)*(z4-3*z3+3*z2-z1);
    T C0=z2,C1=((T)1/6)*(-11*z2+18*z3-9*z4+2*z5),C2=(T).5*(2*z2-5*z3+4*z4-z5),C3=((T)1/6)*(-z2+3*z3-3*z4+z5);
    T PA=A0+(A1+(A2+A3*x)*x)*x,PB=B0+(B1+(B2+B3*x)*x)*x,PC=C0+(C1+(C2+C3*x)*x)*x;
    T cA=(T)1/20*(x-2)*(x-3),cB=-(T)1/10*(x+2)*(x-3),cC=(T)1/20*(x+2)*(x+1);
    T ISA=((T)1/180)*(244*z0*z0-1659*z0*z1+1854*z0*z2-683*z0*z3+2976*z1*z1-6927*z1*z2+2634*z1*z3+4326*z2*z2-3579*z2*z3+814*z3*z3);
    T ISB=((T)1/180)*(244*z1*z1-1269*z1*z2+1074*z1*z3-293*z1*z4+1986*z2*z2-3777*z2*z3+1074*z2*z4+1986*z3*z3-1269*z3*z4+244*z4*z4);
    T ISC=((T)1/180)*(814*z2*z2-3579*z2*z3+2634*z2*z4-683*z2*z5+4326*z3*z3-6927*z3*z4+1854*z3*z5+2976*z4*z4-1659*z4*z5+244*z5*z5);
    T aA=cA/sqr(eps+ISA),aB=cB/sqr(eps+ISB),aC=cC/sqr(eps+ISC);
    T ww=(T)1/(aA+aB+aC),wA=aA*ww,wB=aB*ww,wC=aC*ww;
    return PA*wA+PB*wB+PC*wC;
}

template float WENO_Interpolation(float x,float z0,float z1,float z2,float z3,float z4,float z5,float eps);
template double WENO_Interpolation(double x,double z0,double z1,double z2,double z3,double z4,double z5,double eps);
}
