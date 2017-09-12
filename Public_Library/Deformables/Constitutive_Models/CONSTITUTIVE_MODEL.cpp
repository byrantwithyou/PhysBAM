//#####################################################################
// Copyright 2003-2007, Ron Fedkiw, Geoffrey Irving, Igor Neverov, Andrew Selle, Eftychios Sifakis, Russell Howes, Joseph Teran.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Core/Log/DEBUG_UTILITIES.h>
#include <Core/Matrices/MATRIX_1X1.h>
#include <Core/Matrices/MATRIX_2X2.h>
#include <Core/Matrices/MATRIX_3X3.h>
#include <Deformables/Constitutive_Models/CONSTITUTIVE_MODEL.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class T,int d> CONSTITUTIVE_MODEL<T,d>::
CONSTITUTIVE_MODEL()
    :enforce_definiteness(false),constant_lambda(0),constant_mu(0),constant_alpha(0),constant_beta(0),
    lambda(0),mu(0),alpha(0),beta(0)
{
}
//#####################################################################
// Destructor
//#####################################################################
template<class T,int d> CONSTITUTIVE_MODEL<T,d>::
~CONSTITUTIVE_MODEL()
{
}
//#####################################################################
// Function Maximum_Elastic_Stiffness
//#####################################################################
template<class T,int d> T CONSTITUTIVE_MODEL<T,d>::
Maximum_Elastic_Stiffness(const int id) const // for elastic CFL computation
{
    T id_mu=Mu(id),id_lambda=Lambda(id);
    return id_lambda+2*id_mu;
}
//#####################################################################
// Function Maximum_Damping_Stiffness
//#####################################################################
template<class T,int d> T CONSTITUTIVE_MODEL<T,d>::
Maximum_Damping_Stiffness(const int id) const // for damping CFL computation
{
    T id_alpha=Alpha(id),id_beta=Beta(id);
    return id_alpha+2*id_beta;
}
//#####################################################################
// Function Isotropic_Stress_Derivative
//#####################################################################
template<class T,int d> void CONSTITUTIVE_MODEL<T,d>::
Isotropic_Stress_Derivative(const DIAGONAL_MATRIX<T,d>& F,DIAGONALIZED_ISOTROPIC_STRESS_DERIVATIVE<TV>& dPi_dF,const int id) const
{
    PHYSBAM_FUNCTION_IS_NOT_DEFINED();
}
//#####################################################################
// Function P_From_Strain_Rate_Forces_Size
//#####################################################################
template<class T,int d> int CONSTITUTIVE_MODEL<T,d>::
P_From_Strain_Rate_Forces_Size() const
{
    PHYSBAM_FUNCTION_IS_NOT_DEFINED();
}
//#####################################################################
// Function P_From_Strain_Rate_First_Half
//#####################################################################
template<class T,int d> void CONSTITUTIVE_MODEL<T,d>::
P_From_Strain_Rate_First_Half(const DIAGONAL_MATRIX<T,d>& F,ARRAY_VIEW<T> aggregate,const MATRIX<T,d>& F_dot,const int id) const
{
    PHYSBAM_FUNCTION_IS_NOT_DEFINED();
}
//#####################################################################
// Function P_From_Strain_Rate_Second_Half
//#####################################################################
template<class T,int d> MATRIX<T,d> CONSTITUTIVE_MODEL<T,d>::
P_From_Strain_Rate_Second_Half(const DIAGONAL_MATRIX<T,d>& F,const ARRAY_VIEW<const T> aggregate,const int id) const
{
    PHYSBAM_FUNCTION_IS_NOT_DEFINED();
}
//#####################################################################
// Function Update_Lame_Constants
//#####################################################################
template<class T,int d> void CONSTITUTIVE_MODEL<T,d>::
Update_Lame_Constants(const T youngs_modulus_input, const T poissons_ratio_input,const T Rayleigh_coefficient_input)
{
    constant_lambda=youngs_modulus_input*poissons_ratio_input/((1+poissons_ratio_input)*(1-2*poissons_ratio_input));
    constant_mu=youngs_modulus_input/(2*(1+poissons_ratio_input));
    constant_alpha=Rayleigh_coefficient_input*constant_lambda;
    constant_beta=Rayleigh_coefficient_input*constant_mu;
}
//#####################################################################
// Function P_From_Strain_Rate
//#####################################################################
template<class T,int d> MATRIX<T,d> CONSTITUTIVE_MODEL<T,d>::
P_From_Strain_Rate(const DIAGONAL_MATRIX<T,d>& F,const MATRIX<T,d>& F_dot,const int id) const
{
    PHYSBAM_FUNCTION_IS_NOT_DEFINED();
}
namespace PhysBAM{
template class CONSTITUTIVE_MODEL<float,1>;
template class CONSTITUTIVE_MODEL<float,2>;
template class CONSTITUTIVE_MODEL<float,3>;
template class CONSTITUTIVE_MODEL<double,1>;
template class CONSTITUTIVE_MODEL<double,2>;
template class CONSTITUTIVE_MODEL<double,3>;
}
