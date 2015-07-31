//#####################################################################
// Copyright 2015
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Hybrid_Methods/Forces/VOLUME_PRESERVING_OB_NEO_HOOKEAN.h>
namespace PhysBAM{
//#####################################################################
// Constructor
//#####################################################################
template<class TV> VOLUME_PRESERVING_OB_NEO_HOOKEAN<TV>::
VOLUME_PRESERVING_OB_NEO_HOOKEAN()
    :mu(0),lambda(0)
{
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> VOLUME_PRESERVING_OB_NEO_HOOKEAN<TV>::
~VOLUME_PRESERVING_OB_NEO_HOOKEAN()
{}
//#####################################################################
// Function Resize
//#####################################################################
template<class TV> void VOLUME_PRESERVING_OB_NEO_HOOKEAN<TV>::
Resize(int n)
{
    psi.Resize(n);
    H.Resize(n);
    N.Resize(n);
    G.Resize(n);
    M.Resize(n);
    c.Resize(n);
    S_mat.Resize(n);
    F_mat.Resize(n);
}
//#####################################################################
// Function Precompute
//#####################################################################
template<class TV> void VOLUME_PRESERVING_OB_NEO_HOOKEAN<TV>::
Precompute(const MATRIX<T,TV::m>& F,const SYMMETRIC_MATRIX<T,TV::m>& S,int p)
{
    T J=F.Determinant();
    T K=S.Determinant();
    T e=1./(TV::m);
    if(K<=0) LOG::cout <<"WARNING: K Det <= 0" << std::endl;
    if(J<=0) LOG::cout <<"WARNING: J Det <= 0" << std::endl;
    
    psi(p)=(mu/2.0)*pow(sqr(J),e)*pow(K,-e)*S.Trace()-mu*log(J)+(lambda/2)*sqr(J-1);
    G(p)=F.Inverse();
    M(p)=S.Inverse();
    S_mat(p)=S;
    F_mat(p)=F;
    H(p)=G(p).Transposed()*(mu*e*pow(sqr(J),e)*pow(K,-e)*S.Trace()-mu+lambda*(J-1)*J); //F part
    N(p)=(M(p)*-e*S.Trace()+IDENTITY_MATRIX<T,TV::m>())*(mu/2)*pow(sqr(J),e)*pow(K,-e); //S part
    c(p)=mu*e*pow(sqr(J),e)*pow(K,-e);   
}
//#####################################################################
// Function Energy_Density
//#####################################################################
template<class TV> typename TV::SCALAR VOLUME_PRESERVING_OB_NEO_HOOKEAN<TV>::
Energy_Density(int p) const
{
    return psi(p);
}
//#####################################################################
// Function Gradient
//#####################################################################
template<class TV> void VOLUME_PRESERVING_OB_NEO_HOOKEAN<TV>::
Gradient(MATRIX<T,TV::m>& dF,SYMMETRIC_MATRIX<T,TV::m>& dS,int p) const
{
    dF=H(p);
    dS=N(p);
}
//#####################################################################
// Function Hessian
//#####################################################################
template<class TV> void VOLUME_PRESERVING_OB_NEO_HOOKEAN<TV>::
Hessian(const MATRIX<T,TV::m>& F,const SYMMETRIC_MATRIX<T,TV::m>& S,
    MATRIX<T,TV::m>& dF,MATRIX<T,TV::m>& dS,int p) const
{
    MATRIX<T,TV::m> A=G(p).Transposed()*(G(p).Transposed().Double_Contract(F));
    MATRIX<T,TV::m> B=(G(p)*F*G(p)).Transposed();
    MATRIX<T,TV::m> D=M(p)*S*M(p);
    MATRIX<T,TV::m> E=M(p)*S.Trace();
    T e = 1./(TV::m);
    T d=G(p).Transposed().Double_Contract(F);
    T g=M(p).Double_Contract(S);
    T h=S.Trace();
    T m=(G(p).Transposed()*S.Transposed()).Trace();
    T n=(M(p)*F.Transposed()).Trace();
    T k=F.Trace();
    T J=F_mat(p).Determinant();
    
    MATRIX<T,TV::m> first= (M(p)*(-e)*(S_mat(p).Trace()) + IDENTITY_MATRIX<T,TV::m>())*g;
    MATRIX<T,TV::m> secon= M(p)*h;
    MATRIX<T,TV::m> third= -D*S_mat(p).Trace();
    MATRIX<T,TV::m> fourt= M(p)*(-e)*(S_mat(p).Trace()) + IDENTITY_MATRIX<T,TV::m>();

    dS=(first+secon+third)*(-c(p)/2.0) + (fourt)*c(p)*d;
    
    dF=(2*e*A-B)*c(p)*S_mat(p).Trace() + (A*J+(A-B)*(J-1))*lambda*J + mu*B + G(p).Transposed()*(-e*g*S_mat(p).Trace()+h)*c(p) ;    
}
template class VOLUME_PRESERVING_OB_NEO_HOOKEAN<VECTOR<float,2> >;
template class VOLUME_PRESERVING_OB_NEO_HOOKEAN<VECTOR<float,3> >;
template class VOLUME_PRESERVING_OB_NEO_HOOKEAN<VECTOR<double,2> >;
template class VOLUME_PRESERVING_OB_NEO_HOOKEAN<VECTOR<double,3> >;
}
