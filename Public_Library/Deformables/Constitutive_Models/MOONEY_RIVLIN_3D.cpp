//#####################################################################
// Copyright 2004-2005, Geoffrey Irving, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Deformables/Constitutive_Models/DIAGONALIZED_ISOTROPIC_STRESS_DERIVATIVE.h>
#include <Deformables/Constitutive_Models/MOONEY_RIVLIN_3D.h>
namespace PhysBAM{
//#####################################################################
// Constructor
//#####################################################################
template<class T> MOONEY_RIVLIN_3D<T>::
MOONEY_RIVLIN_3D(const T mu_10_input,const T mu_01_input,const T kappa_input,const T Rayleigh_coefficient,const T failure_threshold_input)
    :mu_10(mu_10_input),mu_01(mu_01_input),kappa(kappa_input),failure_threshold(failure_threshold_input)
{
    constant_mu=2*(mu_10+mu_01);
    constant_lambda=((T)1/3)*(8*mu_01-4*mu_10)+kappa;
    constant_alpha=Rayleigh_coefficient*constant_lambda;
    constant_beta=Rayleigh_coefficient*constant_mu;
}
//#####################################################################
// Function P_From_Strain
//#####################################################################
template<class T> DIAGONAL_MATRIX<T,3> MOONEY_RIVLIN_3D<T>::
P_From_Strain(const DIAGONAL_MATRIX<T,3>& F,const int id) const
{
    DIAGONAL_MATRIX<T,3> F_threshold=F.Clamp_Min(failure_threshold),C=F_threshold*F_threshold,F_cube=C*F_threshold,F_inverse=F_threshold.Inverse();
    T I_C=C.Trace(),II_C=(C*C).Trace(),J=F_threshold.Determinant(),Jcc=pow(J,-((T)2/3));
    return (2*Jcc*(mu_10+Jcc*mu_01*I_C))*F_threshold-(2*Jcc*Jcc*mu_01)*F_cube+((kappa*log(J)-((T)2/3)*Jcc*(mu_10*I_C+Jcc*mu_01*(I_C*I_C-II_C))))*F_inverse;
}
//#####################################################################
// Function P_From_Strain_Rate
//#####################################################################
template<class T> MATRIX<T,3> MOONEY_RIVLIN_3D<T>::
P_From_Strain_Rate(const DIAGONAL_MATRIX<T,3>& F,const MATRIX<T,3>& F_dot,const int id) const
{
    T id_alpha=Alpha(id),id_beta=Beta(id);
    SYMMETRIC_MATRIX<T,3> strain_rate=F_dot.Symmetric_Part(); 
    return 2*id_beta*strain_rate+id_alpha*strain_rate.Trace();
}
//#####################################################################
// Function Isotropic_Stress_Derivative
//#####################################################################
template<class T> void MOONEY_RIVLIN_3D<T>::
Isotropic_Stress_Derivative(const DIAGONAL_MATRIX<T,3>& F,DIAGONALIZED_ISOTROPIC_STRESS_DERIVATIVE<VECTOR<T,3> >& dPi_dF,const int id) const
{
    DIAGONAL_MATRIX<T,3> F_threshold=F.Clamp_Min(failure_threshold),C=F_threshold*F_threshold,F_cube=C*F_threshold,F_inverse=F_threshold.Inverse();
    T I_C=C.Trace(),II_C=(C*C).Trace(),J=F_threshold.Determinant(),Jcc=pow(J,-((T)2/3));
    SYMMETRIC_MATRIX<T,3> alpha;
    alpha.x00=2*Jcc*(mu_10+Jcc*mu_01*(I_C-C.x.x-C.x.x));alpha.x10=2*Jcc*(mu_10+Jcc*mu_01*(I_C-C.x.x-C.x.x));alpha.x20=2*Jcc*(mu_10+Jcc*mu_01*(I_C-C.x.z-C.x.x));
    alpha.x11=2*Jcc*(mu_10+Jcc*mu_01*(I_C-C.x.y-C.x.y));alpha.x21=2*Jcc*(mu_10+Jcc*mu_01*(I_C-C.x.z-C.x.y));alpha.x22=2*Jcc*(mu_10+Jcc*mu_01*(I_C-C.x.z-C.x.z));
    SYMMETRIC_MATRIX<T,3> beta;
    beta.x00=(Jcc*((T)2/3)*(mu_10*I_C+Jcc*mu_01*(I_C*I_C-II_C))-kappa*log(J))*F_inverse.x.x*F_inverse.x.x-2*Jcc*Jcc*mu_01*F_threshold.x.x*F_threshold.x.x;
    beta.x10=(Jcc*((T)2/3)*(mu_10*I_C+Jcc*mu_01*(I_C*I_C-II_C))-kappa*log(J))*F_inverse.x.y*F_inverse.x.x-2*Jcc*Jcc*mu_01*F_threshold.x.y*F_threshold.x.x;
    beta.x20=(Jcc*((T)2/3)*(mu_10*I_C+Jcc*mu_01*(I_C*I_C-II_C))-kappa*log(J))*F_inverse.x.z*F_inverse.x.x-2*Jcc*Jcc*mu_01*F_threshold.x.z*F_threshold.x.x;
    beta.x11=(Jcc*((T)2/3)*(mu_10*I_C+Jcc*mu_01*(I_C*I_C-II_C))-kappa*log(J))*F_inverse.x.y*F_inverse.x.y-2*Jcc*Jcc*mu_01*F_threshold.x.y*F_threshold.x.y;
    beta.x21=(Jcc*((T)2/3)*(mu_10*I_C+Jcc*mu_01*(I_C*I_C-II_C))-kappa*log(J))*F_inverse.x.z*F_inverse.x.y-2*Jcc*Jcc*mu_01*F_threshold.x.z*F_threshold.x.y;
    beta.x22=(Jcc*((T)2/3)*(mu_10*I_C+Jcc*mu_01*(I_C*I_C-II_C))-kappa*log(J))*F_inverse.x.z*F_inverse.x.z-2*Jcc*Jcc*mu_01*F_threshold.x.z*F_threshold.x.z;
    SYMMETRIC_MATRIX<T,3> eta;
    eta.x00=4*Jcc*Jcc*mu_01;
    eta.x20=-Jcc*((T)1/3)*(4*mu_10+8*Jcc*mu_01*I_C);
    eta.x21=8*((T)1/3)*Jcc*Jcc*mu_01;
    eta.x22=Jcc*((T)1/9)*(4*mu_10*I_C+Jcc*8*mu_01*(I_C*I_C-II_C))+kappa;
    MATRIX<T,3> F_base(F_threshold.x,F_cube.x,F_inverse.x);
    SYMMETRIC_MATRIX<T,3> gamma=SYMMETRIC_MATRIX<T,3>::Conjugate(F_base,eta);
    dPi_dF.H(0,0)=alpha.x00+beta.x00+gamma.x00;dPi_dF.H(1,1)=alpha.x11+beta.x11+gamma.x11;dPi_dF.H(2,2)=alpha.x22+beta.x22+gamma.x22;
    dPi_dF.H(1,0)=gamma.x10;dPi_dF.H(2,0)=gamma.x20;dPi_dF.H(2,1)=gamma.x21;
    dPi_dF.B(2)=alpha.x10;dPi_dF.B(1)=alpha.x20;dPi_dF.B(0)=alpha.x21;
    dPi_dF.C(2)=beta.x10;dPi_dF.C(1)=beta.x20;dPi_dF.C(0)=beta.x20;
    if(enforce_definiteness) dPi_dF.Enforce_Definiteness();
}
//#####################################################################
// Function Energy_Density
//#####################################################################
template<class T> T MOONEY_RIVLIN_3D<T>::
Energy_Density(const DIAGONAL_MATRIX<T,3>& F,const int id) const
{
    DIAGONAL_MATRIX<T,3> F_threshold=F.Clamp_Min(failure_threshold),C=F_threshold*F_threshold;
    T I_C=C.Trace(),II_C=(C*C).Trace(),J=F_threshold.Determinant(),Jcc=pow(J,-((T)2/3));
    return (T)0.5*kappa*sqr(log(J))+(mu_10*I_C+mu_01*II_C*Jcc)*Jcc;
}
template class MOONEY_RIVLIN_3D<float>;
template class MOONEY_RIVLIN_3D<double>;
}
