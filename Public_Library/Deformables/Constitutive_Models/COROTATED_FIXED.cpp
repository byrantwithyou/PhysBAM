//#####################################################################
// Copyright 2012, Alexey Stomakhin.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Core/Log/LOG.h>
#include <Core/Math_Tools/cube.h>
#include <Core/Matrices/DIAGONAL_MATRIX.h>
#include <Core/Matrices/MATRIX_1X1.h>
#include <Core/Matrices/MATRIX_2X2.h>
#include <Core/Matrices/MATRIX_3X3.h>
#include <Core/Matrices/SYMMETRIC_MATRIX_1X1.h>
#include <Core/Matrices/SYMMETRIC_MATRIX_2X2.h>
#include <Core/Matrices/SYMMETRIC_MATRIX_3X3.h>
#include <Deformables/Constitutive_Models/COROTATED_FIXED.h>
#include <Deformables/Constitutive_Models/DIAGONALIZED_ISOTROPIC_STRESS_DERIVATIVE.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class T,int d> COROTATED_FIXED<T,d>::
COROTATED_FIXED(const T youngs_modulus_input,const T poissons_ratio_input,const T Rayleigh_coefficient,const T failure_threshold_input)
    :panic_threshold((T)1e-6)
{
    Set_Parameters(youngs_modulus_input,poissons_ratio_input,Rayleigh_coefficient);
}
//#####################################################################
// Destructor
//#####################################################################
template<class T,int d> COROTATED_FIXED<T,d>::
~COROTATED_FIXED()
{
}
//#####################################################################
// Function Set_Parameters
//#####################################################################
template<class T,int d> void COROTATED_FIXED<T,d>::
Set_Parameters(const T youngs_modulus_input,const T poissons_ratio_input,const T Rayleigh_coefficient)
{
    assert(poissons_ratio_input>-1&&poissons_ratio_input<.5);
    youngs_modulus=youngs_modulus_input;
    poissons_ratio=poissons_ratio_input;
    constant_lambda=youngs_modulus*poissons_ratio/((1+poissons_ratio)*(1-2*poissons_ratio));
    constant_mu=youngs_modulus/(2*(1+poissons_ratio));
    constant_alpha=Rayleigh_coefficient*constant_lambda;
    constant_beta=Rayleigh_coefficient*constant_mu;
}
//#####################################################################
// Function Zero_Out_Mu
//#####################################################################
template<class T,int d> void COROTATED_FIXED<T,d>::
Zero_Out_Mu()
{
    constant_mu=0;
    constant_beta=0;
}
//#####################################################################
// Function P_From_Strain
//#####################################################################
template<class T,int d> DIAGONAL_MATRIX<T,d> COROTATED_FIXED<T,d>::
P_From_Strain(const DIAGONAL_MATRIX<T,d>& F,const int id) const
{
    T id_mu=Mu(id),id_lambda=Lambda(id);
    DIAGONAL_MATRIX<T,d> Fm1=F-1;
    return 2*id_mu*Fm1+id_lambda*(F.Determinant()-1)*F.Cofactor_Matrix();
}
//#####################################################################
// Function P_From_Strain_Rate
//#####################################################################
template<class T,int d> MATRIX<T,d> COROTATED_FIXED<T,d>::
P_From_Strain_Rate(const DIAGONAL_MATRIX<T,d>& F,const MATRIX<T,d>& F_dot,const int id) const
{
    T id_alpha=Alpha(id),id_beta=Beta(id);
    SYMMETRIC_MATRIX<T,d> strain_rate=F_dot.Symmetric_Part(); // use linear damping because of problems with inverting elements...
    return 2*id_beta*strain_rate+id_alpha*strain_rate.Trace();
}
//#####################################################################
// Function P_From_Strain_Rate_Forces_Size
//#####################################################################
template<class T,int d> int COROTATED_FIXED<T,d>::
P_From_Strain_Rate_Forces_Size() const
{
    return sizeof(MATRIX<T,d>)/sizeof(T);
}
//#####################################################################
// Function P_From_Strain_Rate_First_Half
//#####################################################################
template<class T,int d> void COROTATED_FIXED<T,d>::
P_From_Strain_Rate_First_Half(const DIAGONAL_MATRIX<T,d>& F,ARRAY_VIEW<T> aggregate,const MATRIX<T,d>& F_dot,const int id) const
{
    T id_alpha=Alpha(id),id_beta=Beta(id);
    SYMMETRIC_MATRIX<T,d> strain_rate=F_dot.Symmetric_Part(); // use linear damping because of problems with inverting elements...
    T sb=sqrt(2*id_beta);
    T dd=sb/TV::m;
    T sa=sqrt(id_alpha/TV::m+dd*dd)-dd;
    SYMMETRIC_MATRIX<T,d> s=sb*strain_rate+sa*strain_rate.Trace();
    *(MATRIX<T,d>*)aggregate.Get_Array_Pointer()+=s;
}
//#####################################################################
// Function P_From_Strain_Rate_Second_Half
//#####################################################################
template<class T,int d> MATRIX<T,d> COROTATED_FIXED<T,d>::
P_From_Strain_Rate_Second_Half(const DIAGONAL_MATRIX<T,d>& F,ARRAY_VIEW<const T> aggregate,const int id) const
{
    T id_alpha=Alpha(id),id_beta=Beta(id);
    SYMMETRIC_MATRIX<T,d> strain_rate=(*(const MATRIX<T,d>*)aggregate.Get_Array_Pointer()).Symmetric_Part(); // use linear damping because of problems with inverting elements...
    T sb=sqrt(2*id_beta);
    T dd=sb/TV::m;
    T sa=sqrt(id_alpha/TV::m+dd*dd)-dd;
    return sb*strain_rate+sa*strain_rate.Trace();
}
//#####################################################################
// Function Isotropic_Stress_Derivative_Helper
//#####################################################################
template<class T,int d> void COROTATED_FIXED<T,d>::
Isotropic_Stress_Derivative_Helper(const DIAGONAL_MATRIX<T,1>& F,DIAGONALIZED_ISOTROPIC_STRESS_DERIVATIVE<VECTOR<T,1> >& dP_dF,const int id) const
{
    T id_mu=Mu(id),id_lambda=Lambda(id);
    T mu=id_mu,la=id_lambda, mu2=2*mu;
    
    dP_dF.H(0,0)=mu2+la;
    if(enforce_definiteness) dP_dF.Enforce_Definiteness();
}
//#####################################################################
// Function Isotropic_Stress_Derivative_Helper
//#####################################################################
template<class T,int d> void COROTATED_FIXED<T,d>::
Isotropic_Stress_Derivative_Helper(const DIAGONAL_MATRIX<T,2>& F,DIAGONALIZED_ISOTROPIC_STRESS_DERIVATIVE<VECTOR<T,2> >& dP_dF,const int id) const
{
    T id_mu=Mu(id),id_lambda=Lambda(id);
    T mu=id_mu,la=id_lambda, mu2=2*mu;
    T J=F.Determinant();
    T d01=F.x.x+F.x.y;if(fabs(d01)<panic_threshold) d01=d01<0?-panic_threshold:panic_threshold;
    
    dP_dF.H(0,0)=mu2+la*sqr(F.x.y);
    dP_dF.H(1,1)=mu2+la*sqr(F.x.x);
    dP_dF.H(1,0)=la*(2*J-1);
    dP_dF.C(0)=mu2/d01-la*(J-1);
    dP_dF.B(0)=mu2*(1-1/d01);
    if(enforce_definiteness) dP_dF.Enforce_Definiteness();
}
//#####################################################################
// Function Isotropic_Stress_Derivative
//#####################################################################
template<class T,int d> void COROTATED_FIXED<T,d>::
Isotropic_Stress_Derivative_Helper(const DIAGONAL_MATRIX<T,3>& F,DIAGONALIZED_ISOTROPIC_STRESS_DERIVATIVE<VECTOR<T,3> >& dPi_dF,const int id) const
{
    T id_mu=Mu(id),id_lambda=Lambda(id);
    T mu=id_mu,la=id_lambda,mu2=2*mu,J=F.Determinant();
    DIAGONAL_MATRIX<T,3> F_cofactor=F.Cofactor_Matrix();
    T d01=F.x.x+F.x.y;if(fabs(d01)<panic_threshold) d01=d01<0?-panic_threshold:panic_threshold;
    T d02=F.x.x+F.x.z;if(fabs(d02)<panic_threshold) d02=d02<0?-panic_threshold:panic_threshold;
    T d12=F.x.y+F.x.z;if(fabs(d12)<panic_threshold) d12=d12<0?-panic_threshold:panic_threshold;
 
    dPi_dF.H(0,0)=mu2+la*sqr(F_cofactor.x.x);
    dPi_dF.H(1,1)=mu2+la*sqr(F_cofactor.x.y);
    dPi_dF.H(2,2)=mu2+la*sqr(F_cofactor.x.z);
    dPi_dF.H(1,0)=la*F.x.z*(2*J-1);
    dPi_dF.H(2,0)=la*F.x.y*(2*J-1);
    dPi_dF.H(2,1)=la*F.x.x*(2*J-1);
    dPi_dF.C(2)=mu2/d01+la*F.x.z*(1-J);
    dPi_dF.C(1)=mu2/d02+la*F.x.y*(1-J);
    dPi_dF.C(0)=mu2/d12+la*F.x.x*(1-J);
    dPi_dF.B(2)=mu2*(1-1/d01);
    dPi_dF.B(1)=mu2*(1-1/d02);
    dPi_dF.B(0)=mu2*(1-1/d12);
    if(enforce_definiteness) dPi_dF.Enforce_Definiteness();
}
//#####################################################################
// Function Isotropic_Stress_Derivative
//#####################################################################
template<class T,int d> void COROTATED_FIXED<T,d>::
Isotropic_Stress_Derivative(const DIAGONAL_MATRIX<T,d>& F,DIAGONALIZED_ISOTROPIC_STRESS_DERIVATIVE<TV>& dPi_dF,const int id) const
{
    Isotropic_Stress_Derivative_Helper(F,dPi_dF,id);
}
//#####################################################################
// Function Energy_Density
//#####################################################################
template<class T,int d> T COROTATED_FIXED<T,d>::
Energy_Density(const DIAGONAL_MATRIX<T,d>& F,const int id) const
{
    T id_mu=Mu(id),id_lambda=Lambda(id);
    DIAGONAL_MATRIX<T,d> Fm1=F-1;
    return id_mu*(Fm1*Fm1).Trace()+(T).5*id_lambda*sqr(F.Determinant()-1);
}
//#####################################################################
// Function Robust_Divided_Pressure
//#####################################################################
template<class T,int d> T COROTATED_FIXED<T,d>::
Robust_Divided_Pressure(T J,const int id) const
{
    // Computes dpsi_dJ/(J-I) robustly. (should fail if not a pressure-based model)
    T id_mu=Mu(id),id_lambda=Lambda(id);
    PHYSBAM_ASSERT(!id_mu);
    return id_lambda;
}
//#####################################################################
// Function Robust_Divided_Pressure
//#####################################################################
template<class T,int d> T COROTATED_FIXED<T,d>::
Pressure_Bound(T J,const int id) const
{
    T id_mu=Mu(id),id_lambda=Lambda(id);
    PHYSBAM_ASSERT(!id_mu);
    return id_lambda;
}
namespace PhysBAM{
template class COROTATED_FIXED<float,1>;
template class COROTATED_FIXED<float,2>;
template class COROTATED_FIXED<float,3>;
template class COROTATED_FIXED<double,1>;
template class COROTATED_FIXED<double,2>;
template class COROTATED_FIXED<double,3>;
}
