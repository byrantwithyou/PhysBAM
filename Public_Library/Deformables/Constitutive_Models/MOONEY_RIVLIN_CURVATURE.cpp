//#####################################################################
// Copyright 2015.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Deformables/Constitutive_Models/MOONEY_RIVLIN_CURVATURE.h>
#include <Tools/Auto_Diff/AUTO_DIFF_EXT.h>
#include <Tools/Auto_Diff/AUTO_HESS_EXT.h>
namespace PhysBAM{

//#####################################################################
// Constructor
//#####################################################################
template<class T>
MOONEY_RIVLIN_CURVATURE<T>::
MOONEY_RIVLIN_CURVATURE(T c1,T c2,T thickness)
:c1(c1),c2(c2),thickness(thickness)
{}

template<class T> T MOONEY_RIVLIN_CURVATURE<T>::
Potential_Energy(TV A1,TV A2,TV dA11,TV dA12,TV dA22,const VECTOR<TM,3>& G0_inv,const VECTOR<T,3>& G0_det,TM2& ge,TT& he) const
{
   auto a1=Hess_From_Var<5,0>(A1);
   auto a2=Hess_From_Var<5,1>(A2);
   auto da11=Hess_From_Var<5,2>(dA11);
   auto da12=Hess_From_Var<5,3>(dA12);
   auto da22=Hess_From_Var<5,4>(dA22);

   auto a3=a1.Cross(a2);
   auto da31=da11.Cross(a2)+a1.Cross(da12);
   auto da32=da12.Cross(a2)+a1.Cross(da22);

   auto m_sqr=a3.Magnitude_Squared();
   auto dm_sqr1=2*da31.Dot(a3);
   auto dm_sqr2=2*da32.Dot(a3);

   T J0=G0_det(1);
   auto beta=J0/m_sqr;

   auto lambda_n=beta*a3;
   auto dlambda_n1=beta*(da31-(dm_sqr1/m_sqr)*a3);
   auto dlambda_n2=beta*(da31-(dm_sqr2/m_sqr)*a3);

   // G=[a1|a2|lambda_n]
   T e;
   {
      auto GtG00=a1.Magnitude_Squared();
      auto GtG01=a1.Dot(a2);
      auto GtG02=a1.Dot(lambda_n);
      auto GtG11=a2.Magnitude_Squared();
      auto GtG12=a2.Dot(lambda_n);
      auto GtG22=lambda_n.Magnitude_Squared();
      const TM& S=G0_inv(1);
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
      T C1=(2*thickness*G0_det(1)/3)*c1;
      T C2=(2*thickness*G0_det(1)/3)*c2;
      auto psi=C1*(i1-3)+C2*(i1-3);
      e=psi.x;
      Extract(ge,psi.dx);
      Extract(he,psi.ddx);
   }

   TM2 g;
   TT h;
   // G=[g1|g2|lambda_n]
   {
      T w=-thickness/2;
      auto g1=a1+w*dlambda_n1;
      auto g2=a2+w*dlambda_n2;
      auto GtG00=g1.Magnitude_Squared();
      auto GtG01=g1.Dot(g2);
      auto GtG02=g1.Dot(lambda_n);
      auto GtG11=g2.Magnitude_Squared();
      auto GtG12=g2.Dot(lambda_n);
      auto GtG22=lambda_n.Magnitude_Squared();
      const TM& S=G0_inv(0);
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
      T C1=(thickness*G0_det(0)/6)*c1;
      T C2=(thickness*G0_det(0)/6)*c2;
      auto psi=C1*(i1-3)+C2*(i1-3);
      e+=psi.x;
      Extract(g,psi.dx);
      Extract(h,psi.ddx);
      ge+=g;
      he+=h;
   }

   // G=[g1|g2|lambda_n]
   {
      T w=thickness/2;
      auto g1=a1+w*dlambda_n1;
      auto g2=a2+w*dlambda_n2;
      auto GtG00=g1.Magnitude_Squared();
      auto GtG01=g1.Dot(g2);
      auto GtG02=g1.Dot(lambda_n);
      auto GtG11=g2.Magnitude_Squared();
      auto GtG12=g2.Dot(lambda_n);
      auto GtG22=lambda_n.Magnitude_Squared();
      const TM& S=G0_inv(2);
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
      T C1=(thickness*G0_det(2)/6)*c1;
      T C2=(thickness*G0_det(2)/6)*c2;
      auto psi=C1*(i1-3)+C2*(i1-3);
      e+=psi.x;
      Extract(g,psi.dx);
      Extract(h,psi.ddx);
      ge+=g;
      he+=h;
   }

   return e;
}

template<class T> T MOONEY_RIVLIN_CURVATURE<T>::
Potential_Energy(TV A1,TV A2,TV dA11,TV dA12,TV dA22,const VECTOR<TM,3>& G0_inv,const VECTOR<T,3>& G0_det,TM2& ge) const
{
   auto a1=Diff_From_Var<5,0>(A1);
   auto a2=Diff_From_Var<5,1>(A2);
   auto da11=Diff_From_Var<5,2>(dA11);
   auto da12=Diff_From_Var<5,3>(dA12);
   auto da22=Diff_From_Var<5,4>(dA22);

   auto a3=a1.Cross(a2);
   auto da31=da11.Cross(a2)+a1.Cross(da12);
   auto da32=da12.Cross(a2)+a1.Cross(da22);

   auto m_sqr=a3.Magnitude_Squared();
   auto dm_sqr1=2*da31.Dot(a3);
   auto dm_sqr2=2*da32.Dot(a3);

   T J0=G0_det(1);
   auto beta=J0/m_sqr;

   auto lambda_n=beta*a3;
   auto dlambda_n1=beta*(da31-(dm_sqr1/m_sqr)*a3);
   auto dlambda_n2=beta*(da31-(dm_sqr2/m_sqr)*a3);

   // G=[a1|a2|lambda_n]
   T e;
   {
      auto GtG00=a1.Magnitude_Squared();
      auto GtG01=a1.Dot(a2);
      auto GtG02=a1.Dot(lambda_n);
      auto GtG11=a2.Magnitude_Squared();
      auto GtG12=a2.Dot(lambda_n);
      auto GtG22=lambda_n.Magnitude_Squared();
      const TM& S=G0_inv(1);
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
      T C1=(2*thickness*G0_det(1)/3)*c1;
      T C2=(2*thickness*G0_det(1)/3)*c2;
      auto psi=C1*(i1-3)+C2*(i1-3);
      e=psi.x;
      Extract(ge,psi.dx);
   }

   TM2 g;
   // G=[g1|g2|lambda_n]
   {
      T w=-thickness/2;
      auto g1=a1+w*dlambda_n1;
      auto g2=a2+w*dlambda_n2;
      auto GtG00=g1.Magnitude_Squared();
      auto GtG01=g1.Dot(g2);
      auto GtG02=g1.Dot(lambda_n);
      auto GtG11=g2.Magnitude_Squared();
      auto GtG12=g2.Dot(lambda_n);
      auto GtG22=lambda_n.Magnitude_Squared();
      const TM& S=G0_inv(0);
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
      T C1=(thickness*G0_det(0)/6)*c1;
      T C2=(thickness*G0_det(0)/6)*c2;
      auto psi=C1*(i1-3)+C2*(i1-3);
      e+=psi.x;
      Extract(g,psi.dx);
      ge+=g;
   }

   // G=[g1|g2|lambda_n]
   {
      T w=thickness/2;
      auto g1=a1+w*dlambda_n1;
      auto g2=a2+w*dlambda_n2;
      auto GtG00=g1.Magnitude_Squared();
      auto GtG01=g1.Dot(g2);
      auto GtG02=g1.Dot(lambda_n);
      auto GtG11=g2.Magnitude_Squared();
      auto GtG12=g2.Dot(lambda_n);
      auto GtG22=lambda_n.Magnitude_Squared();
      const TM& S=G0_inv(2);
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
      T C1=(thickness*G0_det(2)/6)*c1;
      T C2=(thickness*G0_det(2)/6)*c2;
      auto psi=C1*(i1-3)+C2*(i1-3);
      e+=psi.x;
      Extract(g,psi.dx);
      ge+=g;
   }

   return e;
}

template class MOONEY_RIVLIN_CURVATURE<float>;
template class MOONEY_RIVLIN_CURVATURE<double>;

}
