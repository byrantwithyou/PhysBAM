//#####################################################################
// Copyright 2013, Chenfanfu Jiang
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Tools/Random_Numbers/RANDOM_NUMBERS.h>
#include "MPM_CONSTITUTIVE_MODEL.h"
namespace PhysBAM{
//#####################################################################
// Constructor
//#####################################################################
template<class TV> MPM_CONSTITUTIVE_MODEL<TV>::
MPM_CONSTITUTIVE_MODEL()
{}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> MPM_CONSTITUTIVE_MODEL<TV>::
~MPM_CONSTITUTIVE_MODEL()
{}
//#####################################################################
// Function Compute_Helper_Quantities_Using_F
//#####################################################################
template<class TV> void MPM_CONSTITUTIVE_MODEL<TV>::
Compute_Helper_Quantities_Using_F(const MATRIX<T,TV::m>& Fe,const MATRIX<T,TV::m>& Fp,T& Je,T& Jp,MATRIX<T,TV::m>& Ue,DIAGONAL_MATRIX<T,TV::m>& SIGMAe,MATRIX<T,TV::m>& Ve,MATRIX<T,TV::m>& Re,MATRIX<T,TV::m>& Se) const
{
    Je=Fe.Determinant();
    Jp=Fp.Determinant();
    Fe.Fast_Singular_Value_Decomposition(Ue,SIGMAe,Ve);
    Re=Ue*Ve.Transposed();
    Se=Re.Transposed()*Fe;
}
//#####################################################################
// Function Compute_Energy_Density_Psi
//#####################################################################
template<class TV> typename TV::SCALAR MPM_CONSTITUTIVE_MODEL<TV>::
Compute_Elastic_Energy_Density_Psi(const T& mu,const T& lambda,const MATRIX<T,TV::m>& Fe,const MATRIX<T,TV::m>& Re,const T& Je) const
{
    return (Fe-Re).Frobenius_Norm_Squared()*mu+0.5*lambda*sqr(Je-1);
}
//#####################################################################
// Function Compute_dPsi_dFe
//#####################################################################
template<class TV> MATRIX<typename TV::SCALAR,TV::m> MPM_CONSTITUTIVE_MODEL<TV>::
Compute_dPsi_dFe(const T& mu,const T& lambda,const MATRIX<T,TV::m>& Fe,const MATRIX<T,TV::m>& Re,const T& Je) const
{
    return (Fe-Re)*2.0*mu+Fe.Inverse_Transposed()*Je*lambda*(Je-1);
}
//#####################################################################
// Helper Function Compute_dJFinvT // 2d
//#####################################################################
template<class T> static MATRIX<T,2> Compute_dJFinvT(const MATRIX<T,2>& F,const MATRIX<T,2>& dF)
{
    return MATRIX<T,2>(dF(1,1),-dF(0,1),-dF(1,0),dF(0,0));
}
//#####################################################################
// Helper Function Compute_dJFinvT // 3d
//#####################################################################
template<class T> static MATRIX<T,3> Compute_dJFinvT(const MATRIX<T,3>&F,const MATRIX<T,3>& dF)
{
    MATRIX<T,3> result;
    result(0,0)=dF(1,1)*F(2,2)+F(1,1)*dF(2,2)-dF(2,1)*F(1,2)-F(2,1)*dF(1,2);
    result(0,1)=dF(2,0)*F(1,2)+F(2,0)*dF(1,2)-dF(1,0)*F(2,2)-F(1,0)*dF(2,2);
    result(0,2)=dF(1,0)*F(2,1)+F(1,0)*dF(2,1)-dF(2,0)*F(1,1)-F(2,0)*dF(1,1);
    result(1,0)=dF(2,1)*F(0,2)+F(2,1)*dF(0,2)-dF(0,1)*F(2,2)-F(0,1)*dF(2,2);
    result(1,1)=dF(0,0)*F(2,2)+F(0,0)*dF(2,2)-dF(2,0)*F(0,2)-F(2,0)*dF(0,2);
    result(1,2)=dF(2,0)*F(0,1)+F(2,0)*dF(0,1)-dF(0,0)*F(2,1)-F(0,0)*dF(2,1);
    result(2,0)=dF(0,1)*F(1,2)+F(0,1)*dF(1,2)-dF(1,1)*F(0,2)-F(1,1)*dF(0,2);
    result(2,1)=dF(1,0)*F(0,2)+F(1,0)*dF(0,2)-dF(0,0)*F(1,2)-F(0,0)*dF(1,2);
    result(2,2)=dF(0,0)*F(1,1)+F(0,0)*dF(1,1)-dF(1,0)*F(0,1)-F(1,0)*dF(0,1);
    return result;
}
//#####################################################################
// Helper Function Compute_dR // 2d
//#####################################################################
template<class T> static MATRIX<T,2> Compute_dR(const MATRIX<T,2>& R,const MATRIX<T,2>& S,const MATRIX<T,2>& dF)
{
    T y=((R.Transposed()*dF)-(dF.Transposed()*R))(0,1)/(S(0,0)+S(1,1));
    return R*MATRIX<T,2>(0,-y,y,0);
}
//#####################################################################
// Helper Function Compute_dR // 3d
//#####################################################################
template<class T> static MATRIX<T,3> Compute_dR(const MATRIX<T,3>& R,const MATRIX<T,3>& S,const MATRIX<T,3>& dF)
{
    MATRIX<T,3> H=(R.Transposed()*dF)-(dF.Transposed()*R);
    VECTOR<T,3> b(H(0,1),H(0,2),H(1,2));
    MATRIX<T,3> A(S(0,0)+S(1,1),S(1,2),-S(0,2),S(1,2),S(0,0)+S(2,2),S(1,0),-S(0,2),S(0,1),S(1,1)+S(2,2));
    VECTOR<T,3> x=A.Solve_Linear_System(b);
    return R*MATRIX<T,3>(0,-x(0),-x(1),x(0),0,-x(2),x(1),x(2),0);
}
//#####################################################################
// Function Compute_d2Psi_dFe_dFe_Action_dF
//#####################################################################
template<class TV> MATRIX<typename TV::SCALAR,TV::m> MPM_CONSTITUTIVE_MODEL<TV>::
Compute_d2Psi_dFe_dFe_Action_dF(const T& mu,const T& lambda,const MATRIX<T,TV::m>& Fe,const T& Je,const MATRIX<T,TV::m>& Re,const MATRIX<T,TV::m>& Se,const MATRIX<T,TV::m>& dF) const
{
    MATRIX<T,TV::m> JFinvT=Fe.Inverse_Transposed()*Je;
    MATRIX<T,TV::m> dJFinvT=Compute_dJFinvT(Fe,dF);
    MATRIX<T,TV::m> dR=Compute_dR(Re,Se,dF);
    T contraction=MATRIX<T,TV::m>::Inner_Product(JFinvT,dF);
    return 2*mu*(dF-dR)+lambda*JFinvT*contraction+lambda*(Je-1)*dJFinvT;
}
//#####################################################################
// Function Derivative_Test
//#####################################################################
template<class TV> void MPM_CONSTITUTIVE_MODEL<TV>::
Derivative_Test() const
{
    MATRIX<T,TV::m> Fe=MATRIX<T,TV::m>::Identity_Matrix();
    MATRIX<T,TV::m> Fp=MATRIX<T,TV::m>::Identity_Matrix();
    MATRIX<T,TV::m> Ue,Ve,Re,Se;
    DIAGONAL_MATRIX<T,TV::m> SIGMAe;
    T Je,Jp;
    T youngs_modulus=3000,poisson_ratio=0.3;
    T mu=youngs_modulus/(2.0*(1.0+poisson_ratio));
    T lambda=youngs_modulus*poisson_ratio/((1.0+poisson_ratio)*(1.0-2.0*poisson_ratio));

    RANDOM_NUMBERS<T> rand_generator;
    MATRIX<T,TV::m> F1,F2,dF;
    T Psi1,Psi2;
    MATRIX<T,TV::m> P1,P2,dP1,dP2;
    T eps=1e-5;
    for(int i=0;i<10000;i++){
        rand_generator.Fill_Uniform(F1,-1,1);
        rand_generator.Fill_Uniform(dF,-eps,eps);
        Fe=F1;
        Compute_Helper_Quantities_Using_F(Fe,Fp,Je,Jp,Ue,SIGMAe,Ve,Re,Se);
        Psi1=Compute_Elastic_Energy_Density_Psi(mu,lambda,Fe,Re,Je);
        P1=Compute_dPsi_dFe(mu,lambda,Fe,Re,Je);
        dP1=Compute_d2Psi_dFe_dFe_Action_dF(mu,lambda,Fe,Je,Re,Se,dF);
        F2=F1+dF;
        Fe=F2;
        Compute_Helper_Quantities_Using_F(Fe,Fp,Je,Jp,Ue,SIGMAe,Ve,Re,Se);
        Psi2=Compute_Elastic_Energy_Density_Psi(mu,lambda,Fe,Re,Je);
        P2=Compute_dPsi_dFe(mu,lambda,Fe,Re,Je);
        dP2=Compute_d2Psi_dFe_dFe_Action_dF(mu,lambda,Fe,Je,Re,Se,dF);

        T energy_test_result=((Psi2-Psi1)-(P2+P1).Times_Transpose(dF).Trace()/2)/eps;
        T force_test_result=((P2-P1)-(dP2+dP1)/2).Frobenius_Norm()/eps;
   
        LOG::cout<<"Energy Test: "<<energy_test_result<<std::endl;
        LOG::cout<<"P test     : "<<force_test_result<<std::endl;}
}
//#####################################################################
template class MPM_CONSTITUTIVE_MODEL<VECTOR<float,2> >;
template class MPM_CONSTITUTIVE_MODEL<VECTOR<float,3> >;
template class MPM_CONSTITUTIVE_MODEL<VECTOR<double,2> >;
template class MPM_CONSTITUTIVE_MODEL<VECTOR<double,3> >;
}
