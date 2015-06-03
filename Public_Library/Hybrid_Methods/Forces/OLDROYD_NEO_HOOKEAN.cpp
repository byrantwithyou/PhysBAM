//#####################################################################
// Copyright 2015, Craig Schroeder.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Hybrid_Methods/Forces/OLDROYD_NEO_HOOKEAN.h>
namespace PhysBAM{
//#####################################################################
// Constructor
//#####################################################################
template<class TV> OLDROYD_NEO_HOOKEAN<TV>::
OLDROYD_NEO_HOOKEAN()
    :mu(0),lambda(0)
{
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> OLDROYD_NEO_HOOKEAN<TV>::
~OLDROYD_NEO_HOOKEAN()
{}
//#####################################################################
// Function Resize
//#####################################################################
template<class TV> void OLDROYD_NEO_HOOKEAN<TV>::
Resize(int n)
{
    psi.Resize(n);
    H.Resize(n);
    P.Resize(n);
    Q.Resize(n);
    b.Resize(n);
    c.Resize(n);
}
//#####################################################################
// Function Precompute
//#####################################################################
template<class TV> void OLDROYD_NEO_HOOKEAN<TV>::
Precompute(const MATRIX<T,TV::m>& F,const SYMMETRIC_MATRIX<T,TV::m>& S,int p)
{
    T J=F.Determinant();
    psi(p)=mu/2*(S.Trace()-TV::m)-mu*log(J)+lambda/2*sqr(J-1);
    MATRIX<T,TV::m> HH=F.Cofactor_Matrix();
    H(p)=HH;
    T a=lambda*(J-1)-mu/J;
    P(p)=a*HH;
    Q(p)=SYMMETRIC_MATRIX<T,TV::m>()+mu/2;
    b(p)=a/J;
    c(p)=lambda*(2*J-1)/J;
}
//#####################################################################
// Function Energy_Density
//#####################################################################
template<class TV> typename TV::SCALAR OLDROYD_NEO_HOOKEAN<TV>::
Energy_Density(int p) const
{
    return psi(p);
}
//#####################################################################
// Function Gradient
//#####################################################################
template<class TV> void OLDROYD_NEO_HOOKEAN<TV>::
Gradient(MATRIX<T,TV::m>& dF,SYMMETRIC_MATRIX<T,TV::m>& dS,int p) const
{
    dF=P(p);
    dS=Q(p);
}
//#####################################################################
// Function Hessian
//#####################################################################
template<class TV> void OLDROYD_NEO_HOOKEAN<TV>::
Hessian(const MATRIX<T,TV::m>& F,const SYMMETRIC_MATRIX<T,TV::m>& S,
    MATRIX<T,TV::m>& dF,SYMMETRIC_MATRIX<T,TV::m>& dS,int p) const
{
    MATRIX<T,TV::m> A=H(p).Transpose_Times(F);
    dF=c(p)*A.Trace()*H(p)-b(p)*H(p).Times_Transpose(A);
    dS=SYMMETRIC_MATRIX<T,TV::m>();
}
template class OLDROYD_NEO_HOOKEAN<VECTOR<float,2> >;
template class OLDROYD_NEO_HOOKEAN<VECTOR<float,3> >;
template class OLDROYD_NEO_HOOKEAN<VECTOR<double,2> >;
template class OLDROYD_NEO_HOOKEAN<VECTOR<double,3> >;
}
