//#####################################################################
// Copyright 2015.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Tools/Auto_Diff/AUTO_DIFF_EXT.h>
#include <Tools/Auto_Diff/AUTO_HESS_EXT.h>
#include <Deformables/Constitutive_Models/MOONEY_RIVLIN_CURVATURE.h>
namespace PhysBAM{

//#####################################################################
// Constructor
//#####################################################################
template<class T> MOONEY_RIVLIN_CURVATURE<T>::
MOONEY_RIVLIN_CURVATURE(T c1,T c2,T thickness)
:c1(c1),c2(c2),thickness(thickness)
{}
//#####################################################################
// Function Extract
//#####################################################################
template<class T> template<class T1> void MOONEY_RIVLIN_CURVATURE<T>::
Extract_Helper(const T1& in,TT* he) const
{
    TT h;
    Extract(h,in.ddx);
    *he+=h;
}
//#####################################################################
// Function Extract
//#####################################################################
template<class T> template<class T1> void MOONEY_RIVLIN_CURVATURE<T>::
Extract_Helper(const T1& in,int he) const
{
}
//#####################################################################
// Function Helper
//#####################################################################
template<class T> template<class T1,class T2,class T3,class T4,class T5,class T6> T MOONEY_RIVLIN_CURVATURE<T>::
Helper(const T1& a1,const T2& a2,const T3& lambda_n,const T4& dlambda_n1,const T5& dlambda_n2,const TM& S,TM2& ge,T6 he,T w,T scale) const
{
    auto g1=a1+w*dlambda_n1;
    auto g2=a2+w*dlambda_n2;
    auto GtG00=g1.Magnitude_Squared();
    auto GtG01=g1.Dot(g2);
    auto GtG02=g1.Dot(lambda_n);
    auto GtG11=g2.Magnitude_Squared();
    auto GtG12=g2.Dot(lambda_n);
    auto GtG22=lambda_n.Magnitude_Squared();
    auto c00=sqr(S(0,0))*GtG00+2*S(0,0)*S(1,0)*GtG01+
        2*S(0,0)*S(2,0)*GtG02+sqr(S(1,0))*GtG11+
        2*S(1,0)*S(2,0)*GtG12+sqr(S(2,0))*GtG22; 
    auto c01=S(0,0)*S(0,1)*GtG00+
        (S(0,1)*S(1,0)+S(0,0)*S(1,1))*GtG01+
        (S(0,1)*S(2,0)+S(0,0)*S(2,1))*GtG02+
        S(1,0)*S(1,1)*GtG11+(S(1,1)*S(2,0)+
            S(1,0)*S(2,1))*GtG12+S(2,0)*S(2,1)*GtG22; 
    auto c02=S(0,0)*S(0,2)*GtG00+
        (S(0,2)*S(1,0)+S(0,0)*S(1,2))*GtG01+
        (S(0,2)*S(2,0)+S(0,0)*S(2,2))*GtG02+
        S(1,0)*S(1,2)*GtG11+(S(1,2)*S(2,0)+
            S(1,0)*S(2,2))*GtG12+S(2,0)*S(2,2)*GtG22;
    auto c11=sqr(S(0,1))*GtG00+2*S(0,1)*S(1,1)*GtG01+
        2*S(0,1)*S(2,1)*GtG02+sqr(S(1,1))*GtG11+
        2*S(1,1)*S(2,1)*GtG12+sqr(S(2,1))*GtG22; 
    auto c12=S(0,1)*S(0,2)*GtG00+
        (S(0,2)*S(1,1)+S(0,1)*S(1,2))*GtG01+
        (S(0,2)*S(2,1)+S(0,1)*S(2,2))*GtG02+
        S(1,1)*S(1,2)*GtG11+(S(1,2)*S(2,1)+
            S(1,1)*S(2,2))*GtG12+S(2,1)*S(2,2)*GtG22;  
    auto c22=sqr(S(0,2))*GtG00+2*S(0,2)*S(1,2)*GtG01+
        2*S(0,2)*S(2,2)*GtG02+sqr(S(1,2))*GtG11+
        2*S(1,2)*S(2,2)*GtG12+sqr(S(2,2))*GtG22;
    auto s=c11+c22;
    auto i1=c00+s;
    auto i2=c00*s+c11*c22-sqr(c01)-sqr(c02)-sqr(c12);
    T C1=scale*c1;
    T C2=scale*c2;
    auto psi=C1*(i1-3)+C2*(i1-3);
    TM2 g;
    Extract(g,psi.dx);
    ge+=g;
    Extract_Helper(psi,he);
    return psi.x;
}
//#####################################################################
// Function Potential_Energy_Helper
//#####################################################################
template<class T> template<class T1,class T2,class T3,class T4,class T5,class T6> T MOONEY_RIVLIN_CURVATURE<T>::
Potential_Energy_Helper(const T1& a1,const T2& a2,const T3& da11,const T4& da12,const T5& da22,const VECTOR<TM,3>& G0_inv,const VECTOR<T,3>& G0_det,TM2& ge,T6 he,T weight) const
{
   auto a3=a1.Cross(a2);
   auto da31=da11.Cross(a2)+a1.Cross(da12);
   auto da32=da12.Cross(a2)+a1.Cross(da22);

   auto m_sqr=a3.Magnitude_Squared();
   auto dm_sqr1=(T)2*da31.Dot(a3);
   auto dm_sqr2=(T)2*da32.Dot(a3);

   T J0=G0_det(1);
   auto beta=J0/m_sqr;

   auto lambda_n=beta*a3;
   auto dlambda_n1=beta*(da31-(dm_sqr1/m_sqr)*a3);
   auto dlambda_n2=beta*(da31-(dm_sqr2/m_sqr)*a3);
   T scale=weight*thickness/6;
   T e=0;
   ge=TM2();
   e+=Helper(a1,a2,lambda_n,dlambda_n1,dlambda_n2,G0_inv(1),ge,he,0,4*scale*G0_det(1));
   e+=Helper(a1,a2,lambda_n,dlambda_n1,dlambda_n2,G0_inv(0),ge,he,-thickness/2,scale*G0_det(0));
   e+=Helper(a1,a2,lambda_n,dlambda_n1,dlambda_n2,G0_inv(2),ge,he,thickness/2,scale*G0_det(2));

   return e;
}
//#####################################################################
// Function Potential_Energy
//#####################################################################
template<class T> T MOONEY_RIVLIN_CURVATURE<T>::
Potential_Energy(TV A1,TV A2,TV dA11,TV dA12,TV dA22,const VECTOR<TM,3>& G0_inv,const VECTOR<T,3>& G0_det,TM2& ge,TT& he,T weight) const
{
    typedef DIFF_LAYOUT<T,TV::m,TV::m,TV::m,TV::m,TV::m> LAYOUT;
    auto a1=Hess_From_Var<LAYOUT,0>(A1);
    auto a2=Hess_From_Var<LAYOUT,1>(A2);
    auto da11=Hess_From_Var<LAYOUT,2>(dA11);
    auto da12=Hess_From_Var<LAYOUT,3>(dA12);
    auto da22=Hess_From_Var<LAYOUT,4>(dA22);
    he=TT();
    return Potential_Energy_Helper(a1,a2,da11,da12,da22,G0_inv,G0_det,ge,&he,weight);
}
//#####################################################################
// Function Potential_Energy
//#####################################################################
template<class T> T MOONEY_RIVLIN_CURVATURE<T>::
Potential_Energy(TV A1,TV A2,TV dA11,TV dA12,TV dA22,const VECTOR<TM,3>& G0_inv,const VECTOR<T,3>& G0_det,TM2& ge,T weight) const
{
    typedef DIFF_LAYOUT<T,TV::m,TV::m,TV::m,TV::m,TV::m> LAYOUT;
    auto a1=Diff_From_Var<LAYOUT,0>(A1);
    auto a2=Diff_From_Var<LAYOUT,1>(A2);
    auto da11=Diff_From_Var<LAYOUT,2>(dA11);
    auto da12=Diff_From_Var<LAYOUT,3>(dA12);
    auto da22=Diff_From_Var<LAYOUT,4>(dA22);
    return Potential_Energy_Helper(a1,a2,da11,da12,da22,G0_inv,G0_det,ge,0,weight);
}
template class MOONEY_RIVLIN_CURVATURE<float>;
template class MOONEY_RIVLIN_CURVATURE<double>;
}
