//#####################################################################
// Copyright 2004-2007, Geoffrey Irving, Igor Neverov, Eftychios Sifakis, Joseph Teran.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class DIAGONALIZED_ISOTROPIC_STRESS_DERIVATIVE
//#####################################################################
#include <PhysBAM_Tools/Arrays/SORT.h>
#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Tools/Matrices/DIAGONAL_MATRIX_2X2.h>
#include <PhysBAM_Tools/Matrices/DIAGONAL_MATRIX_3X3.h>
#include <PhysBAM_Tools/Matrices/MATRIX_2X2.h>
#include <PhysBAM_Tools/Matrices/MATRIX_3X3.h>
#include <PhysBAM_Tools/Matrices/SYMMETRIC_MATRIX_2X2.h>
#include <PhysBAM_Tools/Matrices/SYMMETRIC_MATRIX_3X3.h>
#include <PhysBAM_Tools/Vectors/VECTOR.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Constitutive_Models/DIAGONALIZED_ISOTROPIC_STRESS_DERIVATIVE.h>
using namespace PhysBAM;
//#####################################################################
// Function Enforce_Definiteness
//#####################################################################
template<class T> void DIAGONALIZED_ISOTROPIC_STRESS_DERIVATIVE<T,2>::
Enforce_Definiteness()
{   LOG::cout << "WARNING: Definiteness Enforcement called!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << std::endl;
    SYMMETRIC_MATRIX<T,2> A1(x1111,x2211,x2222);DIAGONAL_MATRIX<T,2> D1;MATRIX<T,2> V1;
    A1.Solve_Eigenproblem(D1,V1);D1=D1.Clamp_Min(0);A1=SYMMETRIC_MATRIX<T,2>::Conjugate(V1,D1);
    x1111=A1.x11;x2211=A1.x21;x2222=A1.x22;
    SYMMETRIC_MATRIX<T,2> A2(x2121,x2112,x2121);DIAGONAL_MATRIX<T,2> D2;MATRIX<T,2> V2;
    A2.Solve_Eigenproblem(D2,V2);D2=D2.Clamp_Min(0);A2=SYMMETRIC_MATRIX<T,2>::Conjugate(V2,D2);
    x2121=A2.x11;x2112=A2.x21;
}
//#####################################################################
// Function Enforce_Definiteness
//#####################################################################
template<class T> void DIAGONALIZED_ISOTROPIC_STRESS_DERIVATIVE<T,3>::
Enforce_Definiteness(const T eigenvalue_clamp_percentage,const T epsilon)
{
    SYMMETRIC_MATRIX<T,3> A1(x1111,x2211,x3311,x2222,x3322,x3333);DIAGONAL_MATRIX<T,3> D1;MATRIX<T,3> V1;A1.Fast_Solve_Eigenproblem(D1,V1);
    SYMMETRIC_MATRIX<T,2> A2(x2121,x2112,x2121);DIAGONAL_MATRIX<T,2> D2;MATRIX<T,2> V2;A2.Solve_Eigenproblem(D2,V2);
    SYMMETRIC_MATRIX<T,2> A3(x3131,x3113,x3131);DIAGONAL_MATRIX<T,2> D3;MATRIX<T,2> V3;A3.Solve_Eigenproblem(D3,V3);
    SYMMETRIC_MATRIX<T,2> A4(x3232,x3223,x3232);DIAGONAL_MATRIX<T,2> D4;MATRIX<T,2> V4;A4.Solve_Eigenproblem(D4,V4);
    VECTOR<T,9> eigenvalues;
    eigenvalues(0)=abs(D1.x11);eigenvalues(1)=abs(D1.x22);eigenvalues(2)=abs(D1.x33);
    eigenvalues(3)=abs(D2.x11);eigenvalues(4)=abs(D2.x22);eigenvalues(5)=abs(D3.x11);eigenvalues(6)=abs(D3.x22);eigenvalues(7)=abs(D4.x11);eigenvalues(8)=abs(D4.x22);
    Sort(eigenvalues);T min_nonzero_absolute_eigenvalue=epsilon;for(int i=0;i<9;i++) if(min_nonzero_absolute_eigenvalue<eigenvalues(i)){min_nonzero_absolute_eigenvalue=eigenvalues(i);break;}
    if(D1.x11<-epsilon) D1.x11=eigenvalue_clamp_percentage*min_nonzero_absolute_eigenvalue;else if (D1.x11<(T)0) D1.x11=(T)0;
    if(D1.x22<-epsilon) D1.x22=eigenvalue_clamp_percentage*min_nonzero_absolute_eigenvalue;else if (D1.x22<(T)0) D1.x22=(T)0;
    if(D1.x33<-epsilon) D1.x33=eigenvalue_clamp_percentage*min_nonzero_absolute_eigenvalue;else if (D1.x33<(T)0) D1.x33=(T)0;
    if(D2.x11<-epsilon) D2.x11=eigenvalue_clamp_percentage*min_nonzero_absolute_eigenvalue;else if (D2.x11<(T)0) D2.x11=(T)0;
    if(D2.x22<-epsilon) D2.x22=eigenvalue_clamp_percentage*min_nonzero_absolute_eigenvalue;else if (D2.x22<(T)0) D2.x22=(T)0;
    if(D3.x11<-epsilon) D3.x11=eigenvalue_clamp_percentage*min_nonzero_absolute_eigenvalue;else if (D3.x11<(T)0) D3.x11=(T)0;
    if(D3.x22<-epsilon) D3.x22=eigenvalue_clamp_percentage*min_nonzero_absolute_eigenvalue;else if (D3.x22<(T)0) D3.x22=(T)0;
    if(D4.x11<-epsilon) D4.x11=eigenvalue_clamp_percentage*min_nonzero_absolute_eigenvalue;else if (D4.x11<(T)0) D4.x11=(T)0;
    if(D4.x22<-epsilon) D4.x22=eigenvalue_clamp_percentage*min_nonzero_absolute_eigenvalue;else if (D4.x22<(T)0) D4.x22=(T)0;
    A1=SYMMETRIC_MATRIX<T,3>::Conjugate(V1,D1);x1111=A1.x11;x2211=A1.x21;x3311=A1.x31;x2222=A1.x22;x3322=A1.x32;x3333=A1.x33;
    A2=SYMMETRIC_MATRIX<T,2>::Conjugate(V2,D2);x2121=A2.x11;x2112=A2.x21;
    A3=SYMMETRIC_MATRIX<T,2>::Conjugate(V3,D3);x3131=A3.x11;x3113=A3.x21;
    A4=SYMMETRIC_MATRIX<T,2>::Conjugate(V4,D4);x3232=A4.x11;x3223=A4.x21;
}
//#####################################################################
// Function Compute_From_Singular_Value_Derivatives
//#####################################################################
template<class T> void DIAGONALIZED_ISOTROPIC_STRESS_DERIVATIVE<T,2>::
Compute_From_Singular_Value_Derivatives(const DIAGONAL_MATRIX<T,2>& F,const VECTOR<T,2>& dE_ds,const SYMMETRIC_MATRIX<T,2>& dE_dsds)
{
    T ss1=sqr(F.x11),ss2=sqr(F.x22);
    T s12=1/(ss1-ss2);
    
    x1111=dE_dsds.x11;
    x2211=dE_dsds.x21;
    x2222=dE_dsds.x22;
    x2112=(-dE_ds.y*F.x11+dE_ds.x*F.x22)*s12;
    x2121=(-dE_ds.y*F.x22+dE_ds.x*F.x11)*s12;
}
//#####################################################################
// Function Compute_From_Singular_Value_Derivatives
//#####################################################################
template<class T> void DIAGONALIZED_ISOTROPIC_STRESS_DERIVATIVE<T,3>::
Compute_From_Singular_Value_Derivatives(const DIAGONAL_MATRIX<T,3>& F,const VECTOR<T,3>& dE_ds,const SYMMETRIC_MATRIX<T,3>& dE_dsds)
{
    T ss1=sqr(F.x11),ss2=sqr(F.x22),ss3=sqr(F.x33);
    T s12=1/(ss1-ss2),s13=1/(ss1-ss3),s23=1/(ss2-ss3);
    
    x1111=dE_dsds.x11;
    x2211=dE_dsds.x21;
    x2222=dE_dsds.x22;
    x3311=dE_dsds.x31;
    x3322=dE_dsds.x32;
    x3333=dE_dsds.x33;
    x2112=(-dE_ds.y*F.x11+dE_ds.x*F.x22)*s12;
    x2121=(-dE_ds.y*F.x22+dE_ds.x*F.x11)*s12;
    x3113=(-dE_ds.z*F.x11+dE_ds.x*F.x33)*s13;
    x3131=(-dE_ds.z*F.x33+dE_ds.x*F.x11)*s13;
    x3223=(-dE_ds.z*F.x22+dE_ds.y*F.x33)*s23;
    x3232=(-dE_ds.z*F.x33+dE_ds.y*F.x22)*s23;
}
//#####################################################################
template class DIAGONALIZED_ISOTROPIC_STRESS_DERIVATIVE<float,2>;
template class DIAGONALIZED_ISOTROPIC_STRESS_DERIVATIVE<float,3>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class DIAGONALIZED_ISOTROPIC_STRESS_DERIVATIVE<double,2>;
template class DIAGONALIZED_ISOTROPIC_STRESS_DERIVATIVE<double,3>;
#endif
