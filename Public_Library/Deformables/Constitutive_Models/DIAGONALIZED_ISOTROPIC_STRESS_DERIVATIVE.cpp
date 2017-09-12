//#####################################################################
// Copyright 2004-2007, Geoffrey Irving, Igor Neverov, Eftychios Sifakis, Joseph Teran.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class DIAGONALIZED_ISOTROPIC_STRESS_DERIVATIVE
//#####################################################################
#include <Core/Log/LOG.h>
#include <Core/Matrices/DIAGONAL_MATRIX.h>
#include <Core/Matrices/MATRIX_1X1.h>
#include <Core/Matrices/MATRIX_2X2.h>
#include <Core/Matrices/MATRIX_3X3.h>
#include <Core/Matrices/SYMMETRIC_MATRIX_1X1.h>
#include <Core/Matrices/SYMMETRIC_MATRIX_2X2.h>
#include <Core/Matrices/SYMMETRIC_MATRIX_3X3.h>
#include <Core/Vectors/VECTOR.h>
#include <Deformables/Constitutive_Models/DIAGONALIZED_ISOTROPIC_STRESS_DERIVATIVE.h>
using namespace PhysBAM;
//#####################################################################
// Function Enforce_Definiteness
//#####################################################################
template<class TV> void DIAGONALIZED_ISOTROPIC_STRESS_DERIVATIVE<TV>::
Enforce_Definiteness()
{
    DIAGONAL_MATRIX<T,TV::m> D;
    MATRIX<T,TV::m> V;
    H.Fast_Solve_Eigenproblem(D,V);
    H=SYMMETRIC_MATRIX<T,TV::m>::Conjugate(V,D.Clamp_Min(0));
    for(int i=0;i<B.m;i++){
        SYMMETRIC_MATRIX<T,2> A(B(i),C(i),B(i));
        DIAGONAL_MATRIX<T,2> D;
        MATRIX<T,2> V;
        A.Solve_Eigenproblem(D,V);
        A=SYMMETRIC_MATRIX<T,2>::Conjugate(V,D.Clamp_Min(0));
        B(i)=A(0,0);
        C(i)=A(0,1);}
}
//#####################################################################
// Function Compute_From_Singular_Value_Derivatives
//#####################################################################
template<class TV> void DIAGONALIZED_ISOTROPIC_STRESS_DERIVATIVE<TV>::
Compute_From_Singular_Value_Derivatives(const DIAGONAL_MATRIX<T,TV::m>& F,const TV& dE_ds,const SYMMETRIC_MATRIX<T,TV::m>& dE_dsds)
{
    H=dE_dsds;
    for(int i=0;i<B.m;i++){
        int j=(i+1)%B.m,k=(j+1)%B.m;
        T d=1/(sqr(F.x(j))-sqr(F.x(k)));
        B(i)=(-dE_ds(k)*F.x(k)+dE_ds(j)*F.x(j))*d;
        C(i)=(-dE_ds(k)*F.x(j)+dE_ds(j)*F.x(k))*d;}
}
//#####################################################################
namespace PhysBAM{
template class DIAGONALIZED_ISOTROPIC_STRESS_DERIVATIVE<VECTOR<float,1> >;
template class DIAGONALIZED_ISOTROPIC_STRESS_DERIVATIVE<VECTOR<float,2> >;
template class DIAGONALIZED_ISOTROPIC_STRESS_DERIVATIVE<VECTOR<float,3> >;
template class DIAGONALIZED_ISOTROPIC_STRESS_DERIVATIVE<VECTOR<double,1> >;
template class DIAGONALIZED_ISOTROPIC_STRESS_DERIVATIVE<VECTOR<double,2> >;
template class DIAGONALIZED_ISOTROPIC_STRESS_DERIVATIVE<VECTOR<double,3> >;
}
