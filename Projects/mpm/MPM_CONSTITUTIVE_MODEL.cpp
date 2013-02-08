//#####################################################################
// Copyright 2013, Chenfanfu Jiang
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################

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
// Function Initialize
//#####################################################################
template<class TV> void MPM_CONSTITUTIVE_MODEL<TV>::
Initialize(const T youngs_modulus,const T poisson_ratio,const T hardening_coefficient)
{
    mu0=youngs_modulus/(2.0*(1.0+poisson_ratio));
    lambda0=youngs_modulus*poisson_ratio/((1.0+poisson_ratio)*(1.0-2.0*poisson_ratio));
    xi=hardening_coefficient;
    Fe=MTX::Identity_Matrix();
    Fp=MTX::Identity_Matrix();
}

//#####################################################################
// Function Update_Quantities_Using_Current_Deformation_Gradient
//#####################################################################
template<class TV> void MPM_CONSTITUTIVE_MODEL<TV>::
Update_Quantities_Using_Current_Deformation_Gradient()
{
    Je=Fe.Determinant();
    Jp=Fp.Determinant();
    Fast_Singular_Value_Decomposition(Ue,SIGMAe,Ve);
    Re=Ue*Ve.Transposed();
    Se=Re.Transposed()*Fe;
    T lame_factor=exp(xi*(1.0-Jp));
    mu=mu0*lame_factor;
    lambda=lambda0*lame_factor;
}

//#####################################################################
// Function Energy_Density_Psi
//#####################################################################
template<class TV> T MPM_CONSTITUTIVE_MODEL<TV>::
Energy_Density_Psi()
{
    return (Fe-Re).Frobenius_Norm_Squared()*mu+0.5*lambda*sqr(Je-1);
}

//#####################################################################
// Function dPsi_dFe
//#####################################################################
template<class TV> MTX MPM_CONSTITUTIVE_MODEL<TV>::
dPsi_dFe()
{
    return (Fe-Re)*2.0*mu+Fe.Inverse_Transposed()*Je*lambda*(Je-1);
}

//#####################################################################
// Function d2Psi_dFe_dFe_Action_dF
//#####################################################################
template<class TV> MTX MPM_CONSTITUTIVE_MODEL<TV>::
d2Psi_dFe_dFe_Action_dF(const MTX& dF)
{
    MTX JFinvT=Fe.Inverse_Transposed()*Je;
    //TODO
}


//#####################################################################
template class MPM_CONSTITUTIVE_MODEL<VECTOR<float,1> >;
template class MPM_CONSTITUTIVE_MODEL<VECTOR<float,2> >;
template class MPM_CONSTITUTIVE_MODEL<VECTOR<float,3> >;
template class MPM_CONSTITUTIVE_MODEL<VECTOR<double,1> >;
template class MPM_CONSTITUTIVE_MODEL<VECTOR<double,2> >;
template class MPM_CONSTITUTIVE_MODEL<VECTOR<double,3> >;
template class MPM_CONSTITUTIVE_MODEL<VECTOR<float,1> >;
template class MPM_CONSTITUTIVE_MODEL<VECTOR<float,2> >;
template class MPM_CONSTITUTIVE_MODEL<VECTOR<float,3> >;
template class MPM_CONSTITUTIVE_MODEL<VECTOR<double,1> >;
template class MPM_CONSTITUTIVE_MODEL<VECTOR<double,2> >;
template class MPM_CONSTITUTIVE_MODEL<VECTOR<double,3> >;
}
