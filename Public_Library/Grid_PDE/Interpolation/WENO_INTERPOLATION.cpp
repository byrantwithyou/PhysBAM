//#####################################################################
// Copyright 2017, Ounan Ding, Craig Schroeder.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Core/Arrays_Nd/ARRAYS_ND.h>
#include <Core/Arrays_Nd/ARRAYS_ND_VIEW.h>
#include <Core/Log/LOG.h>
#include <Core/Math_Tools/sqr.h>
#include <Core/Matrices/SYMMETRIC_MATRIX.h>
#include <Core/Vectors/VECTOR.h>
#include <Grid_PDE/Interpolation/WENO_INTERPOLATION.h>
namespace PhysBAM{

// From C. Macdonald and S. Ruuth, "Level set equations on surfaces via the Closest Point Method"
// Interpolate the data at location x.
// If x=0, then interpolates at z2; x=1 interpolates at z3
// Assumes x in [0,1].
// eps is used to prevent division by zero; the paper uses eps=1e-6.
template<class T> T
WENO_Interpolation_Helper(T x,T z0,T z1,T z2,T z3,T z4,T z5,T eps)
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

template<class T,int d> static VECTOR<T,d>
WENO_Interpolation_Helper(T x,const VECTOR<T,d>& z0,const VECTOR<T,d>& z1,const VECTOR<T,d>& z2,
    const VECTOR<T,d>& z3,const VECTOR<T,d>& z4,const VECTOR<T,d>& z5,T eps)
{
    VECTOR<T,d> r;
    for(int i=0;i<d;i++) r(i)=WENO_Interpolation(x,z0(i),z1(i),z2(i),z3(i),z4(i),z5(i),eps);
    return r;
}

template<class T,int d> static SYMMETRIC_MATRIX<T,d>
WENO_Interpolation_Helper(T x,const SYMMETRIC_MATRIX<T,d>& z0,const SYMMETRIC_MATRIX<T,d>& z1,const SYMMETRIC_MATRIX<T,d>& z2,
    const SYMMETRIC_MATRIX<T,d>& z3,const SYMMETRIC_MATRIX<T,d>& z4,const SYMMETRIC_MATRIX<T,d>& z5,T eps)
{
    SYMMETRIC_MATRIX<T,d> r;
    for(int i=0;i<d;i++)
        for(int j=0;j<=i;j++)
            r(i,j)=WENO_Interpolation(x,z0(i,j),z1(i,j),z2(i,j),z3(i,j),z4(i,j),z5(i,j),eps);
    return r;
}

template<class T,class U> U
WENO_Interpolation(T x,U z0,U z1,U z2,U z3,U z4,U z5,T eps)
{
    return WENO_Interpolation_Helper(x,z0,z1,z2,z3,z4,z5,eps);
}

template<class T,class U> U
WENO_Interpolation_Helper(const VECTOR<T,1>& X,const ARRAYS_ND_BASE<U,VECTOR<int,1> >& z,const VECTOR<int,1>& index,T eps)
{
    VECTOR<int,1> r;
    U u[6];
    int i;
    for(i=0,r.x=index.x-2;i<6;i++,r.x++)
         u[i]=z(r);
     return WENO_Interpolation(X.x,u,eps);
}

template<class T,class U> U
WENO_Interpolation_Helper(const VECTOR<T,2>& X,const ARRAYS_ND_BASE<U,VECTOR<int,2> >& z,const VECTOR<int,2>& index,T eps)
{
    VECTOR<int,2> r;
    U u[6],v[6];
    int i,j;
    for(i=0,r.x=index.x-2;i<6;i++,r.x++){
        for(j=0,r.y=index.y-2;j<6;j++,r.y++)
            v[j]=z(r);
        u[i]=WENO_Interpolation(X.y,v,eps);}
    return WENO_Interpolation(X.x,u,eps);
}

template<class T,class U> U
WENO_Interpolation_Helper(const VECTOR<T,3>& X,const ARRAYS_ND_BASE<U,VECTOR<int,3> >& z,const VECTOR<int,3>& index,T eps)
{
    VECTOR<int,3> r;
    U u[6],v[6],w[6];
    int i,j,k;
    for(i=0,r.x=index.x-2;i<6;i++,r.x++){
        for(j=0,r.y=index.y-2;j<6;j++,r.y++){
            for(k=0,r.z=index.z-2;k<6;k++,r.z++)
                w[k]=z(r);
            v[j]=WENO_Interpolation(X.z,w,eps);}
        u[i]=WENO_Interpolation(X.y,v,eps);}
    return WENO_Interpolation(X.x,u,eps);
}

template<class T,class U,int d> U
WENO_Interpolation(const VECTOR<T,d>& x,const ARRAYS_ND_BASE<U,VECTOR<int,d> >& z,const VECTOR<int,d>& index,T eps)
{
    return WENO_Interpolation_Helper(x,z,index,eps);
}

template float WENO_Interpolation(float x,float z0,float z1,float z2,float z3,float z4,float z5,float eps);
template double WENO_Interpolation(double x,double z0,double z1,double z2,double z3,double z4,double z5,double eps);

// Instantiation for field
// T^d -> T(e.g. scalar field),
// T^d -> T^d(e.g. gradient),
// T^d -> T^{dxd}(e.g. Hessian)
#define INSTV(T,d)\
    template T\
        WENO_Interpolation<T,T,d>(VECTOR<T,d> const&,\
            ARRAYS_ND_BASE<T,VECTOR<int,d> > const&,VECTOR<int,d> const&,T);\
    template VECTOR<T,d>\
        WENO_Interpolation<T,VECTOR<T,d>,d>(VECTOR<T,d> const&,\
            ARRAYS_ND_BASE<VECTOR<T,d>,VECTOR<int,d> > const&,VECTOR<int,d> const&,T);\
    template SYMMETRIC_MATRIX<T,d>\
        WENO_Interpolation<T,SYMMETRIC_MATRIX<T,d>,d>(VECTOR<T,d> const&,\
            ARRAYS_ND_BASE<SYMMETRIC_MATRIX<T,d>,VECTOR<int,d> > const&,VECTOR<int,d> const&,T)

INSTV(float,1);
INSTV(float,2);
INSTV(float,3);
INSTV(double,1);
INSTV(double,2);
INSTV(double,3);
}
