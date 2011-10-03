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
    :youngs_modulus(youngs_modulus_input),poissons_ratio(poissons_ratio_input),panic_threshold((T)1e-6)
{
    assert(poissons_ratio>-1&&poissons_ratio<.5);
    constant_lambda=youngs_modulus*poissons_ratio/((1+poissons_ratio)*(1-2*poissons_ratio));
    constant_mu=youngs_modulus/(2*(1+poissons_ratio));
    constant_alpha=Rayleigh_coefficient*constant_lambda;
    constant_beta=Rayleigh_coefficient*constant_mu;
}
//#####################################################################
// Destructor
//#####################################################################
template<class T,int d> COROTATED<T,d>::
~COROTATED()
{
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
#if 0
    DIAGONAL_MATRIX<T,2> F_inverse=F.Clamp_Min(failure_threshold).Inverse();
    T mu_minus_lambda_logJ=constant_mu+constant_lambda*log(F_inverse.Determinant());
    SYMMETRIC_MATRIX<T,2> F_inverse_outer=SYMMETRIC_MATRIX<T,2>::Outer_Product(F_inverse.To_Vector());
    dP_dF.x1111=constant_mu+(constant_lambda+mu_minus_lambda_logJ)*F_inverse_outer.x11;//alpha+beta+gamma
    dP_dF.x2222=constant_mu+(constant_lambda+mu_minus_lambda_logJ)*F_inverse_outer.x22;
    dP_dF.x2211=constant_lambda*F_inverse_outer.x21;//gamma
    dP_dF.x2121=constant_mu;//alpha
    dP_dF.x2112=mu_minus_lambda_logJ*F_inverse_outer.x21;//beta
    if(enforce_definiteness) dP_dF.Enforce_Definiteness();
#endif
}
//#####################################################################
// Function dI_dsigma_Helper
//#####################################################################
namespace{
template<class T> MATRIX<T,3> dI_dsigma_Helper(const DIAGONAL_MATRIX<T,3>& F)
{
    typedef VECTOR<T,3> TV;
    TV FV=F.To_Vector();
    TV row1=(T)2*FV;
    TV row2=(T)4*FV*FV*FV;
    TV row3=(T)2*FV.Product()*F.Cofactor_Matrix().To_Vector();
    return MATRIX<T,3>(row1,row2,row3).Transposed();
}
template<class T> MATRIX<T,2> dI_dsigma_Helper(const DIAGONAL_MATRIX<T,2>& F)
{
    typedef VECTOR<T,2> TV;
    TV FV=F.To_Vector();
    TV row1=(T)2*FV;
    TV row2=(T)2*FV.Product()*F.Cofactor_Matrix().To_Vector();
    return MATRIX<T,2>(row1,row2).Transposed();
}
}
//#####################################################################
// Function Isotropic_Stress_Derivative
//#####################################################################
template<class T,int d> void COROTATED<T,d>::
Isotropic_Stress_Derivative(const DIAGONAL_MATRIX<T,3>& F,DIAGONALIZED_ISOTROPIC_STRESS_DERIVATIVE<T,3>& dPi_dF,const int tetrahedron) const
{
    typedef VECTOR<T,3> TV3;
    DIAGONAL_MATRIX<T,3> Fm1=F-1,F2=F*F,CF=F.Cofactor_Matrix();
    TV3 FV=F.To_Vector(),FV2=F2.To_Vector(),CFV=CF.To_Vector();
    TV3 dE_ds=(T)2*constant_mu*Fm1.To_Vector()+constant_lambda*Fm1.Trace();
    MATRIX<T,3> dE_dsds;for(int i=1;i<=d;i++) for(int j=1;j<=d;j++) dE_dsds(i,j)=constant_lambda;dE_dsds+=2*constant_mu;
    MATRIX<T,3> dI_ds=dI_dsigma_Helper(F);
    TV3 dE_dI;
    if(fabs(dI_ds.Determinant())>1e-15) dE_dI=dI_ds.Solve_Linear_System(dE_ds);
    else{
        MATRIX<T,3> U,V;
        DIAGONAL_MATRIX<T,3> singular_values;
        dI_ds.Fast_Singular_Value_Decomposition(U,singular_values,V);
        singular_values=DIAGONAL_MATRIX<T,3>(clamp_min(singular_values.To_Vector(),panic_threshold));
        dE_dI=V*singular_values.Solve_Linear_System(U.Transpose_Times(dE_ds));}

    MATRIX<T,3> cross_term_III=(T)2*dE_dI.z*MATRIX<T,3>::Outer_Product(CFV,CFV);
    MATRIX<T,3> diff_III=(T)2*cross_term_III-(T)2*CF*CF;
    MATRIX<T,3> a;for(int i=1;i<=d;i++) for(int j=1;j<=d;j++) a(i,j)=(T)2*dE_dI.x+(T)4*(FV2(i)+FV2(j))*dE_dI.y;
    MATRIX<T,3> b=(T)4*dE_dI.y*MATRIX<T,3>::Outer_Product(FV,FV)-cross_term_III;
    MATRIX<T,3> g=dE_dsds-(T)2*dE_dI.x-(T)12*dE_dI.y*F2+(T)2*CF*CF*dE_dI.z;
    MATRIX<T,3> A=g+a.Diagonal_Part()+b.Diagonal_Part();

    dPi_dF.x1111=A(1,1);
    dPi_dF.x2222=A(2,2);
    dPi_dF.x3333=A(3,3);
    dPi_dF.x2211=A(2,1);
    dPi_dF.x3311=A(3,1);
    dPi_dF.x3322=A(3,2);
    dPi_dF.x2121=a(1,2);
    dPi_dF.x3131=a(1,3);
    dPi_dF.x3232=a(2,3);
    dPi_dF.x2112=b(1,2);
    dPi_dF.x3113=b(1,3);
    dPi_dF.x3223=b(2,3);
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
