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
#include <Deformables/Constitutive_Models/TAIT_PRESSURE_FORCE.h>
#include <Deformables/Constitutive_Models/DIAGONALIZED_ISOTROPIC_STRESS_DERIVATIVE.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class TV> TAIT_PRESSURE_FORCE<TV>::
TAIT_PRESSURE_FORCE(const T stiffness_input,const T tait_const_input)
    :stiffness(stiffness_input),tait_const(tait_const_input)
{
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> TAIT_PRESSURE_FORCE<TV>::
~TAIT_PRESSURE_FORCE()
{
}
//#####################################################################
// Function P_From_Strain
//#####################################################################
// clamp to hyperbola to avoid indefiniteness "automatically"
template<class TV> DIAGONAL_MATRIX<typename TV::SCALAR,TV::m> TAIT_PRESSURE_FORCE<TV>::
P_From_Strain(const DIAGONAL_MATRIX<typename TV::SCALAR,TV::m>& F,const int id) const
{
    T J=F.Determinant();
    return -tait_const*stiffness*(exp((1-J)/tait_const)-1)*F.Cofactor_Matrix();
}
template<class T> static VECTOR<T,3> Helper(const VECTOR<T,3>& f){return f;}
template<class T> static VECTOR<T,1> Helper(const VECTOR<T,2>& f){return VECTOR<T,1>(1);}
template<class T> static VECTOR<T,0> Helper(const VECTOR<T,1>& f){return VECTOR<T,0>();}
//#####################################################################
// Function Isotropic_Stress_Derivative
//#####################################################################
template<class TV> void TAIT_PRESSURE_FORCE<TV>::
Isotropic_Stress_Derivative(const DIAGONAL_MATRIX<T,TV::m>& F,DIAGONALIZED_ISOTROPIC_STRESS_DERIVATIVE<TV>& dP_dF,const int id) const
{
    T J=F.Determinant(),e=exp((1-J)/tait_const);
    T g=stiffness*(e-1)*tait_const;
    T f=stiffness*e*J-g;
    DIAGONAL_MATRIX<T,TV::m> CF=F.Cofactor_Matrix();
    SYMMETRIC_MATRIX<T,TV::m> hess(stiffness*e*CF*CF);
    typename TV::SPIN Fh=Helper(F.x),div_diff=g*Fh,div_sum=-div_diff;
    for(int i=0;i<TV::SPIN::m;i++){
        int j=(i+1)%TV::m,k=(j+1)%TV::m;
        hess(j,k)=f*Fh(i);}
    dP_dF.Compute_From_Derivatives_And_Differences(hess,div_diff,div_sum);
}
//#####################################################################
// Function Energy_Density
//#####################################################################
template<class TV> typename TV::SCALAR TAIT_PRESSURE_FORCE<TV>::
Energy_Density(const DIAGONAL_MATRIX<T,TV::m>& F,const int id) const
{
    T J=F.Determinant(),e=exp((1-J)/tait_const);
    return tait_const*stiffness*(tait_const*(e-1)+J-1);
}
namespace PhysBAM{
template class TAIT_PRESSURE_FORCE<VECTOR<float,1> >;
template class TAIT_PRESSURE_FORCE<VECTOR<float,2> >;
template class TAIT_PRESSURE_FORCE<VECTOR<float,3> >;
template class TAIT_PRESSURE_FORCE<VECTOR<double,1> >;
template class TAIT_PRESSURE_FORCE<VECTOR<double,2> >;
template class TAIT_PRESSURE_FORCE<VECTOR<double,3> >;
}
