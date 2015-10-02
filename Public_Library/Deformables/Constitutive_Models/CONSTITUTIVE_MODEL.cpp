//#####################################################################
// Copyright 2003-2007, Ron Fedkiw, Geoffrey Irving, Igor Neverov, Andrew Selle, Eftychios Sifakis, Russell Howes, Joseph Teran.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Tools/Log/DEBUG_UTILITIES.h>
#include <Tools/Matrices/MATRIX_1X1.h>
#include <Tools/Matrices/MATRIX_2X2.h>
#include <Tools/Matrices/MATRIX_3X3.h>
#include <Deformables/Constitutive_Models/CONSTITUTIVE_MODEL.h>
using namespace PhysBAM;
template<class T,int d> CONSTITUTIVE_MODEL<T,d>::
CONSTITUTIVE_MODEL()
    :enforce_definiteness(false),constant_lambda(0),constant_mu(0),constant_alpha(0),constant_beta(0)
{
}
template<class T,int d> CONSTITUTIVE_MODEL<T,d>::
~CONSTITUTIVE_MODEL()
{
}
template<class T,int d> T CONSTITUTIVE_MODEL<T,d>::
Maximum_Elastic_Stiffness(const int simplex) const // for elastic CFL computation
{
    return lambda.m?lambda(simplex)+2*mu(simplex):constant_lambda+2*constant_mu;
}
template<class T,int d> T CONSTITUTIVE_MODEL<T,d>::
Maximum_Damping_Stiffness(const int simplex) const // for damping CFL computation
{
    return alpha.m?alpha(simplex)+2*beta(simplex):constant_alpha+2*constant_beta;
}
template<class T,int d> void CONSTITUTIVE_MODEL<T,d>::
Isotropic_Stress_Derivative(const DIAGONAL_MATRIX<T,d>& F,DIAGONALIZED_ISOTROPIC_STRESS_DERIVATIVE<T,d>& dPi_dF,const int simplex) const
{
    PHYSBAM_FUNCTION_IS_NOT_DEFINED();
}
template<class T,int d> int CONSTITUTIVE_MODEL<T,d>::
P_From_Strain_Rate_Forces_Size() const
{
    PHYSBAM_FUNCTION_IS_NOT_DEFINED();
}
template<class T,int d> void CONSTITUTIVE_MODEL<T,d>::
P_From_Strain_Rate_First_Half(const DIAGONAL_MATRIX<T,d>& F,ARRAY_VIEW<T> aggregate,const MATRIX<T,d>& F_dot,const T scale,const int simplex) const
{
    PHYSBAM_FUNCTION_IS_NOT_DEFINED();
}
template<class T,int d> MATRIX<T,d> CONSTITUTIVE_MODEL<T,d>::
P_From_Strain_Rate_Second_Half(const DIAGONAL_MATRIX<T,d>& F,const ARRAY_VIEW<const T> aggregate,const T scale,const int simplex) const
{
    PHYSBAM_FUNCTION_IS_NOT_DEFINED();
}
template<class T,int d> void CONSTITUTIVE_MODEL<T,d>::
Update_Lame_Constants(const T youngs_modulus_input, const T poissons_ratio_input,const T Rayleigh_coefficient_input)
{
    constant_lambda=youngs_modulus_input*poissons_ratio_input/((1+poissons_ratio_input)*(1-2*poissons_ratio_input));
    constant_mu=youngs_modulus_input/(2*(1+poissons_ratio_input));
    constant_alpha=Rayleigh_coefficient_input*constant_lambda;
    constant_beta=Rayleigh_coefficient_input*constant_mu;
}
namespace PhysBAM{
template class CONSTITUTIVE_MODEL<float,1>;
template class CONSTITUTIVE_MODEL<float,2>;
template class CONSTITUTIVE_MODEL<float,3>;
template class CONSTITUTIVE_MODEL<double,1>;
template class CONSTITUTIVE_MODEL<double,2>;
template class CONSTITUTIVE_MODEL<double,3>;
}
