//#####################################################################
// Copyright 2003-2007, Ron Fedkiw, Geoffrey Irving, Igor Neverov, Eftychios Sifakis, Joseph Teran.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Core/Log/LOG.h>
#include <Core/Math_Tools/cube.h>
#include <Core/Math_Tools/pow.h>
#include <Core/Matrices/DIAGONAL_MATRIX.h>
#include <Core/Matrices/MATRIX_2X2.h>
#include <Core/Matrices/MATRIX_3X3.h>
#include <Core/Matrices/SYMMETRIC_MATRIX_2X2.h>
#include <Core/Matrices/SYMMETRIC_MATRIX_3X3.h>
#include <Deformables/Constitutive_Models/COROTATED.h>
#include <Deformables/Constitutive_Models/DIAGONALIZED_ISOTROPIC_STRESS_DERIVATIVE.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class T,int d> COROTATED<T,d>::
COROTATED(const T youngs_modulus_input,const T poissons_ratio_input,const T Rayleigh_coefficient,const T failure_threshold_input)
    :panic_threshold((T)1e-6)
{
    Set_Parameters(youngs_modulus_input,poissons_ratio_input,Rayleigh_coefficient);
}
//#####################################################################
// Destructor
//#####################################################################
template<class T,int d> COROTATED<T,d>::
~COROTATED()
{
}
//#####################################################################
// Function Set_Parameters
//#####################################################################
template<class T,int d> void COROTATED<T,d>::
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
// clamp to hyperbola to avoid indefiniteness "automatically"
template<class T,int d> DIAGONAL_MATRIX<T,d> COROTATED<T,d>::
P_From_Strain(const DIAGONAL_MATRIX<T,d>& F,const int id) const
{
    T id_mu=Mu(id),id_lambda=Lambda(id);
    DIAGONAL_MATRIX<T,d> Fm1=F-1;
    return 2*id_mu*Fm1+id_lambda*Fm1.Trace();
}
//#####################################################################
// Function P_From_Strain_Rate
//#####################################################################
template<class T,int d> MATRIX<T,d> COROTATED<T,d>::
P_From_Strain_Rate(const DIAGONAL_MATRIX<T,d>& F,const MATRIX<T,d>& F_dot,const int id) const
{
    T id_alpha=Alpha(id),id_beta=Beta(id);
    SYMMETRIC_MATRIX<T,d> strain_rate=F_dot.Symmetric_Part(); // use linear damping because of problems with inverting elements...
    return 2*id_beta*strain_rate+id_alpha*strain_rate.Trace();
}
//#####################################################################
// Function P_From_Strain_Rate_Forces_Size
//#####################################################################
template<class T,int d> int COROTATED<T,d>::
P_From_Strain_Rate_Forces_Size() const
{
    return sizeof(MATRIX<T,d>)/sizeof(T);
}
//#####################################################################
// Function P_From_Strain_Rate_First_Half
//#####################################################################
template<class T,int d> void COROTATED<T,d>::
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
template<class T,int d> MATRIX<T,d> COROTATED<T,d>::
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
template<class T,int d> void COROTATED<T,d>::
Isotropic_Stress_Derivative_Helper(const DIAGONAL_MATRIX<T,2>& F,DIAGONALIZED_ISOTROPIC_STRESS_DERIVATIVE<T,2>& dP_dF,const int id) const
{
    T id_mu=Mu(id),id_lambda=Lambda(id);
    T mu=id_mu,la=id_lambda,mu2la=2*mu+la,la2mu2=2*la+2*mu;
    T d01=F.x.x+F.x.y;if(fabs(d01)<panic_threshold) d01=d01<0?-panic_threshold:panic_threshold;
    T i01=la2mu2/d01;
    dP_dF.x0000=mu2la;
    dP_dF.x1001=i01-la;
    dP_dF.x1010=mu2la-i01;
    dP_dF.x1100=la;
    dP_dF.x1111=mu2la;
    if(enforce_definiteness) dP_dF.Enforce_Definiteness();
}
//#####################################################################
// Function Isotropic_Stress_Derivative_Helper
//#####################################################################
template<class T,int d> void COROTATED<T,d>::
Isotropic_Stress_Derivative_Helper(const DIAGONAL_MATRIX<T,3>& F,DIAGONALIZED_ISOTROPIC_STRESS_DERIVATIVE<T,3>& dPi_dF,const int id) const
{
    T id_mu=Mu(id),id_lambda=Lambda(id);
    T mu=id_mu,la=id_lambda,mu2la=2*mu+la,la3mu2=3*la+2*mu;
    T d01=F.x.x+F.x.y;if(fabs(d01)<panic_threshold) d01=d01<0?-panic_threshold:panic_threshold;
    T d02=F.x.x+F.x.z;if(fabs(d02)<panic_threshold) d02=d02<0?-panic_threshold:panic_threshold;
    T d12=F.x.y+F.x.z;if(fabs(d12)<panic_threshold) d12=d12<0?-panic_threshold:panic_threshold;
    T i01=(la3mu2-la*F.x.z)/d01,i02=(la3mu2-la*F.x.y)/d02,i12=(la3mu2-la*F.x.x)/d12;
    dPi_dF.x0000=mu2la;
    dPi_dF.x1111=mu2la;
    dPi_dF.x2222=mu2la;
    dPi_dF.x1100=la;
    dPi_dF.x2200=la;
    dPi_dF.x2211=la;
    dPi_dF.x1001=i01-la;
    dPi_dF.x2002=i02-la;
    dPi_dF.x2112=i12-la;
    dPi_dF.x1010=mu2la-i01;
    dPi_dF.x2020=mu2la-i02;
    dPi_dF.x2121=mu2la-i12;
    if(enforce_definiteness) dPi_dF.Enforce_Definiteness();
}
//#####################################################################
// Function Isotropic_Stress_Derivative
//#####################################################################
template<class T,int d> void COROTATED<T,d>::
Isotropic_Stress_Derivative(const DIAGONAL_MATRIX<T,d>& F,DIAGONALIZED_ISOTROPIC_STRESS_DERIVATIVE<T,d>& dP_dF,const int id) const
{
    Isotropic_Stress_Derivative_Helper(F,dP_dF,id);
}
//#####################################################################
// Function Energy_Density
//#####################################################################
template<class T,int d> T COROTATED<T,d>::
Energy_Density(const DIAGONAL_MATRIX<T,d>& F,const int id) const
{
    T id_mu=Mu(id),id_lambda=Lambda(id);
    DIAGONAL_MATRIX<T,d> Fm1=F-1;
    return id_mu*(Fm1*Fm1).Trace()+(T).5*id_lambda*sqr(Fm1.Trace());
}
namespace PhysBAM{
template class COROTATED<float,2>;
template class COROTATED<float,3>;
template class COROTATED<double,2>;
template class COROTATED<double,3>;
}
