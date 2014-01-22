//#####################################################################
// Copyright 2012, Alexey Stomakhin.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Tools/Log/LOG.h>
#include <Tools/Math_Tools/cube.h>
#include <Tools/Math_Tools/pow.h>
#include <Tools/Matrices/DIAGONAL_MATRIX.h>
#include <Tools/Matrices/MATRIX_2X2.h>
#include <Tools/Matrices/MATRIX_3X3.h>
#include <Tools/Matrices/SYMMETRIC_MATRIX_2X2.h>
#include <Tools/Matrices/SYMMETRIC_MATRIX_3X3.h>
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
// Function P_From_Strain
//#####################################################################
template<class T,int d> DIAGONAL_MATRIX<T,d> COROTATED_FIXED<T,d>::
P_From_Strain(const DIAGONAL_MATRIX<T,d>& F,const T scale,const int simplex) const
{
    T scale_mu=scale*constant_mu,scale_lambda=scale*constant_lambda;
    DIAGONAL_MATRIX<T,d> Fm1=F-1;
    return 2*scale_mu*Fm1+scale_lambda*(F.Determinant()-1)*F.Cofactor_Matrix();
}
//#####################################################################
// Function P_From_Strain_Rate
//#####################################################################
template<class T,int d> MATRIX<T,d> COROTATED_FIXED<T,d>::
P_From_Strain_Rate(const DIAGONAL_MATRIX<T,d>& F,const MATRIX<T,d>& F_dot,const T scale,const int simplex) const
{
    SYMMETRIC_MATRIX<T,d> strain_rate=F_dot.Symmetric_Part(); // use linear damping because of problems with inverting elements...
    return 2*scale*constant_beta*strain_rate+scale*constant_alpha*strain_rate.Trace();
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
P_From_Strain_Rate_First_Half(const DIAGONAL_MATRIX<T,d>& F,ARRAY_VIEW<T> aggregate,const MATRIX<T,d>& F_dot,const T scale,const int simplex) const
{
    SYMMETRIC_MATRIX<T,d> strain_rate=scale*F_dot.Symmetric_Part(); // use linear damping because of problems with inverting elements...
    T sb=sqrt(2*constant_beta);
    T dd=sb/TV::dimension;
    T sa=sqrt(constant_alpha/TV::dimension+dd*dd)-dd;
    SYMMETRIC_MATRIX<T,d> s=sb*strain_rate+sa*strain_rate.Trace();
    *(MATRIX<T,d>*)aggregate.Get_Array_Pointer()+=s;
}
//#####################################################################
// Function P_From_Strain_Rate_Second_Half
//#####################################################################
template<class T,int d> MATRIX<T,d> COROTATED_FIXED<T,d>::
P_From_Strain_Rate_Second_Half(const DIAGONAL_MATRIX<T,d>& F,ARRAY_VIEW<const T> aggregate,const T scale,const int simplex) const
{
    SYMMETRIC_MATRIX<T,d> strain_rate=scale*(*(const MATRIX<T,d>*)aggregate.Get_Array_Pointer()).Symmetric_Part(); // use linear damping because of problems with inverting elements...
    T sb=sqrt(2*constant_beta);
    T dd=sb/TV::dimension;
    T sa=sqrt(constant_alpha/TV::dimension+dd*dd)-dd;
    return sb*strain_rate+sa*strain_rate.Trace();
}
//#####################################################################
// Function Isotropic_Stress_Derivative
//#####################################################################
template<class T,int d> void COROTATED_FIXED<T,d>::
Isotropic_Stress_Derivative(const DIAGONAL_MATRIX<T,2>& F,DIAGONALIZED_ISOTROPIC_STRESS_DERIVATIVE<T,2>& dP_dF,const int triangle) const
{
    T mu=constant_mu,la=constant_lambda, mu2=2*mu;
    T J=F.Determinant();
    T d01=F.x.x+F.x.y;if(fabs(d01)<panic_threshold) d01=d01<0?-panic_threshold:panic_threshold;
    
    dP_dF.x0000=mu2+la*sqr(F.x.y);
    dP_dF.x1111=mu2+la*sqr(F.x.x);
    dP_dF.x1100=la*(2*J-1);
    dP_dF.x1001=mu2/d01-la*(J-1);
    dP_dF.x1010=mu2*(1-1/d01);
    if(enforce_definiteness) dP_dF.Enforce_Definiteness();
}
//#####################################################################
// Function Isotropic_Stress_Derivative
//#####################################################################
template<class T,int d> void COROTATED_FIXED<T,d>::
Isotropic_Stress_Derivative(const DIAGONAL_MATRIX<T,3>& F,DIAGONALIZED_ISOTROPIC_STRESS_DERIVATIVE<T,3>& dPi_dF,const int tetrahedron) const
{
    T mu=constant_mu,la=constant_lambda,mu2=2*mu,J=F.Determinant();
    DIAGONAL_MATRIX<T,3> F_cofactor=F.Cofactor_Matrix();
    T d01=F.x.x+F.x.y;if(fabs(d01)<panic_threshold) d01=d01<0?-panic_threshold:panic_threshold;
    T d02=F.x.x+F.x.z;if(fabs(d02)<panic_threshold) d02=d02<0?-panic_threshold:panic_threshold;
    T d12=F.x.y+F.x.z;if(fabs(d12)<panic_threshold) d12=d12<0?-panic_threshold:panic_threshold;
 
    dPi_dF.x0000=mu2+la*sqr(F_cofactor.x.x);
    dPi_dF.x1111=mu2+la*sqr(F_cofactor.x.y);
    dPi_dF.x2222=mu2+la*sqr(F_cofactor.x.z);
    dPi_dF.x1100=la*F.x.z*(2*J-1);
    dPi_dF.x2200=la*F.x.y*(2*J-1);
    dPi_dF.x2211=la*F.x.x*(2*J-1);
    dPi_dF.x1001=mu2/d01+la*F.x.z*(1-J);
    dPi_dF.x2002=mu2/d02+la*F.x.y*(1-J);
    dPi_dF.x2112=mu2/d12+la*F.x.x*(1-J);
    dPi_dF.x1010=mu2*(1-1/d01);
    dPi_dF.x2020=mu2*(1-1/d02);
    dPi_dF.x2121=mu2*(1-1/d12);
    if(enforce_definiteness) dPi_dF.Enforce_Definiteness();
}
//#####################################################################
// Function Energy_Density
//#####################################################################
template<class T,int d> T COROTATED_FIXED<T,d>::
Energy_Density(const DIAGONAL_MATRIX<T,d>& F,const int simplex) const
{
    DIAGONAL_MATRIX<T,d> Fm1=F-1;
    return constant_mu*(Fm1*Fm1).Trace()+(T).5*constant_lambda*sqr(F.Determinant()-1);
}
namespace PhysBAM{
template class COROTATED_FIXED<float,2>;
template class COROTATED_FIXED<float,3>;
template class COROTATED_FIXED<double,2>;
template class COROTATED_FIXED<double,3>;
}
