//#####################################################################
// Copyright 2003-2007, Ron Fedkiw, Geoffrey Irving, Igor Neverov, Eftychios Sifakis, Joseph Teran.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Tools/Math_Tools/cube.h>
#include <PhysBAM_Tools/Math_Tools/pow.h>
#include <PhysBAM_Tools/Matrices/DIAGONAL_MATRIX_2X2.h>
#include <PhysBAM_Tools/Matrices/DIAGONAL_MATRIX_3X3.h>
#include <PhysBAM_Tools/Matrices/MATRIX_2X2.h>
#include <PhysBAM_Tools/Matrices/MATRIX_3X3.h>
#include <PhysBAM_Tools/Matrices/SYMMETRIC_MATRIX_2X2.h>
#include <PhysBAM_Tools/Matrices/SYMMETRIC_MATRIX_3X3.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Constitutive_Models/COROTATED.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Constitutive_Models/DIAGONALIZED_ISOTROPIC_STRESS_DERIVATIVE.h>
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
    assert(poissons_ratio>-1&&poissons_ratio<.5);
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
P_From_Strain(const DIAGONAL_MATRIX<T,d>& F,const T scale,const int simplex) const
{
    T scale_mu=scale*constant_mu,scale_lambda=scale*constant_lambda;
    DIAGONAL_MATRIX<T,d> Fm1=F-1;
    return 2*scale_mu*Fm1+scale_lambda*Fm1.Trace();
}
//#####################################################################
// Function P_From_Strain_Rate
//#####################################################################
template<class T,int d> MATRIX<T,d> COROTATED<T,d>::
P_From_Strain_Rate(const DIAGONAL_MATRIX<T,d>& F,const MATRIX<T,d>& F_dot,const T scale,const int simplex) const
{
    SYMMETRIC_MATRIX<T,d> strain_rate=F_dot.Symmetric_Part(); // use linear damping because of problems with inverting elements...
    return 2*scale*constant_beta*strain_rate+scale*constant_alpha*strain_rate.Trace();
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
template<class T,int d> MATRIX<T,d> COROTATED<T,d>::
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
template<class T,int d> void COROTATED<T,d>::
Isotropic_Stress_Derivative(const DIAGONAL_MATRIX<T,2>& F,DIAGONALIZED_ISOTROPIC_STRESS_DERIVATIVE<T,2>& dP_dF,const int triangle) const
{
    T mu=constant_mu,la=constant_lambda,mu2la=2*mu+la,la2mu2=2*la+2*mu;
    T d12=F.x11+F.x22;if(fabs(d12)<panic_threshold) d12=d12<0?-panic_threshold:panic_threshold;
    T i12=la2mu2/d12;
    dP_dF.x1111=mu2la;
    dP_dF.x2112=i12-la;
    dP_dF.x2121=mu2la-i12;
    dP_dF.x2211=la;
    dP_dF.x2222=mu2la;
    if(enforce_definiteness) dP_dF.Enforce_Definiteness();
}
//#####################################################################
// Function Isotropic_Stress_Derivative
//#####################################################################
template<class T,int d> void COROTATED<T,d>::
Isotropic_Stress_Derivative(const DIAGONAL_MATRIX<T,3>& F,DIAGONALIZED_ISOTROPIC_STRESS_DERIVATIVE<T,3>& dPi_dF,const int tetrahedron) const
{
    T mu=constant_mu,la=constant_lambda,mu2la=2*mu+la,la3mu2=3*la+2*mu;
    T d12=F.x11+F.x22;if(fabs(d12)<panic_threshold) d12=d12<0?-panic_threshold:panic_threshold;
    T d13=F.x11+F.x33;if(fabs(d13)<panic_threshold) d13=d13<0?-panic_threshold:panic_threshold;
    T d23=F.x22+F.x33;if(fabs(d23)<panic_threshold) d23=d23<0?-panic_threshold:panic_threshold;
    T i12=(la3mu2-la*F.x33)/d12,i13=(la3mu2-la*F.x22)/d13,i23=(la3mu2-la*F.x11)/d23;
    dPi_dF.x1111=mu2la;
    dPi_dF.x2222=mu2la;
    dPi_dF.x3333=mu2la;
    dPi_dF.x2211=la;
    dPi_dF.x3311=la;
    dPi_dF.x3322=la;
    dPi_dF.x2112=i12-la;
    dPi_dF.x3113=i13-la;
    dPi_dF.x3223=i23-la;
    dPi_dF.x2121=mu2la-i12;
    dPi_dF.x3131=mu2la-i13;
    dPi_dF.x3232=mu2la-i23;
    if(enforce_definiteness) dPi_dF.Enforce_Definiteness();
}
//#####################################################################
// Function Energy_Density
//#####################################################################
template<class T,int d> T COROTATED<T,d>::
Energy_Density(const DIAGONAL_MATRIX<T,d>& F,const int simplex) const
{
    DIAGONAL_MATRIX<T,d> Fm1=F-1;
    return constant_mu*(Fm1*Fm1).Trace()+(T).5*constant_lambda*sqr(Fm1.Trace());
}
template class COROTATED<float,2>;
template class COROTATED<float,3>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class COROTATED<double,2>;
template class COROTATED<double,3>;
#endif
