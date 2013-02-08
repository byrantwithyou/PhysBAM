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
    Fe=MATRIX<T,TV::m>::Identity_Matrix();
    Fp=MATRIX<T,TV::m>::Identity_Matrix();
}

//#####################################################################
// Function Update_Quantities_Using_Current_Deformation_Gradient
//#####################################################################
template<class TV> void MPM_CONSTITUTIVE_MODEL<TV>::
Update_Quantities_Using_Current_Deformation_Gradient()
{
    Je=Fe.Determinant();
    Jp=Fp.Determinant();
    Fe.Fast_Singular_Value_Decomposition(Ue,SIGMAe,Ve);
    Re=Ue*Ve.Transposed();
    Se=Re.Transposed()*Fe;
    T lame_factor=exp(xi*(1.0-Jp));
    mu=mu0*lame_factor;
    lambda=lambda0*lame_factor;
}

//#####################################################################
// Function Energy_Density_Psi
//#####################################################################
template<class TV> typename TV::SCALAR MPM_CONSTITUTIVE_MODEL<TV>::
Energy_Density_Psi()
{
    return (Fe-Re).Frobenius_Norm_Squared()*mu+0.5*lambda*sqr(Je-1);
}

//#####################################################################
// Function dPsi_dFe
//#####################################################################
template<class TV> MATRIX<typename TV::SCALAR,TV::m> MPM_CONSTITUTIVE_MODEL<TV>::
dPsi_dFe()
{
    return (Fe-Re)*2.0*mu+Fe.Inverse_Transposed()*Je*lambda*(Je-1);
}

//#####################################################################
// Helper Function Compute_dJFinvT // 2d
//#####################################################################
template<class T> MATRIX<T,2> Compute_dJFinvT(const MATRIX<T,2>& F,const MATRIX<T,2>& dF)
{
    return MATRIX<T,2>(dF(1,1),-dF(0,1),-dF(1,0),dF(0,0));
}

//#####################################################################
// Helper Function Compute_dJFinvT // 3d
//#####################################################################
template<class T> MATRIX<T,3> Compute_dJFinvT(const MATRIX<T,3>&F,const MATRIX<T,3>& dF)
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
template<class T> MATRIX<T,2> Compute_dR(const MATRIX<T,2>& R,const MATRIX<T,2>& S,const MATRIX<T,2>& dF)
{
    T y=((R.Transposed()*dF)-(dF.Transposed()*R))(0,1)/(S(0,0)+S(1,1));
    return R*MATRIX<T,2>(0,-y,y,0);
}

//#####################################################################
// Helper Function Compute_dR // 3d
//#####################################################################
template<class T> MATRIX<T,3> Compute_dR(const MATRIX<T,3>& R,const MATRIX<T,3>& S,const MATRIX<T,3>& dF)
{
    MATRIX<T,3> H=(R.Transposed()*dF)-(dF.Transposed()*R);
    VECTOR<T,3> b(H(0,1),H(0,2),H(1,2));
    MATRIX<T,3> A(S(0,0)+S(1,1),S(1,2),-S(0,2),S(1,2),S(0,0)+S(2,2),S(1,0),-S(0,2),S(0,1),S(1,1)+S(2,2));
    VECTOR<T,3> x=A.Solve_Linear_System(b);
    return R*MATRIX<T,3>(0,-x(0),-x(1),x(0),0,-x(2),x(1),x(2),0);
}

//#####################################################################
// Function d2Psi_dFe_dFe_Action_dF
//#####################################################################
template<class TV> MATRIX<typename TV::SCALAR,TV::m> MPM_CONSTITUTIVE_MODEL<TV>::
d2Psi_dFe_dFe_Action_dF(const MATRIX<T,TV::m>& dF)
{
    MATRIX<T,TV::m> JFinvT=Fe.Inverse_Transposed()*Je;
    MATRIX<T,TV::m> dJFinvT=Compute_dJFinvT(Fe,dF);
    MATRIX<T,TV::m> dR=Compute_dR(Re,Se,dF);
    T contraction=MATRIX<T,TV::m>::Inner_Product(JFinvT,dF);
    return 2*mu*(dF-dR)+lambda*JFinvT*contraction+lambda*(Je-1)*dJFinvT;
}

//#####################################################################
template class MPM_CONSTITUTIVE_MODEL<VECTOR<float,2> >;
template class MPM_CONSTITUTIVE_MODEL<VECTOR<float,3> >;
template class MPM_CONSTITUTIVE_MODEL<VECTOR<double,2> >;
template class MPM_CONSTITUTIVE_MODEL<VECTOR<double,3> >;
}
