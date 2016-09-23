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
template<class T> void DIAGONALIZED_ISOTROPIC_STRESS_DERIVATIVE<T,1>::
Enforce_Definiteness()
{
    x0000=clamp_min(x0000,(T)0);
}
//#####################################################################
// Function Enforce_Definiteness
//#####################################################################
template<class T> void DIAGONALIZED_ISOTROPIC_STRESS_DERIVATIVE<T,2>::
Enforce_Definiteness()
{   LOG::cout << "WARNING: Definiteness Enforcement called!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << std::endl;
    SYMMETRIC_MATRIX<T,2> A0(x0000,x1100,x1111);DIAGONAL_MATRIX<T,2> D0;MATRIX<T,2> V0;
    A0.Solve_Eigenproblem(D0,V0);D0=D0.Clamp_Min(0);A0=SYMMETRIC_MATRIX<T,2>::Conjugate(V0,D0);
    x0000=A0.x00;x1100=A0.x10;x1111=A0.x11;
    SYMMETRIC_MATRIX<T,2> A1(x1010,x1001,x1010);DIAGONAL_MATRIX<T,2> D1;MATRIX<T,2> V1;
    A1.Solve_Eigenproblem(D1,V1);D1=D1.Clamp_Min(0);A1=SYMMETRIC_MATRIX<T,2>::Conjugate(V1,D1);
    x1010=A1.x00;x1001=A1.x10;
}
//#####################################################################
// Function Enforce_Definiteness
//#####################################################################
template<class T> void DIAGONALIZED_ISOTROPIC_STRESS_DERIVATIVE<T,3>::
Enforce_Definiteness(const T eigenvalue_clamp_percentage,const T epsilon)
{
    SYMMETRIC_MATRIX<T,3> A0(x0000,x1100,x2200,x1111,x2211,x2222);DIAGONAL_MATRIX<T,3> D0;MATRIX<T,3> V0;A0.Fast_Solve_Eigenproblem(D0,V0);
    SYMMETRIC_MATRIX<T,2> A1(x1010,x1001,x1010);DIAGONAL_MATRIX<T,2> D1;MATRIX<T,2> V1;A1.Solve_Eigenproblem(D1,V1);
    SYMMETRIC_MATRIX<T,2> A2(x2020,x2002,x2020);DIAGONAL_MATRIX<T,2> D2;MATRIX<T,2> V2;A2.Solve_Eigenproblem(D2,V2);
    SYMMETRIC_MATRIX<T,2> A3(x2121,x2112,x2121);DIAGONAL_MATRIX<T,2> D3;MATRIX<T,2> V3;A3.Solve_Eigenproblem(D3,V3);
    VECTOR<T,9> eigenvalues;
    eigenvalues(0)=abs(D0.x.x);eigenvalues(1)=abs(D0.x.y);eigenvalues(2)=abs(D0.x.z);
    eigenvalues(3)=abs(D1.x.x);eigenvalues(4)=abs(D1.x.y);eigenvalues(5)=abs(D2.x.x);eigenvalues(6)=abs(D2.x.y);eigenvalues(7)=abs(D3.x.x);eigenvalues(8)=abs(D3.x.y);
    eigenvalues.Sort();
    T min_nonzero_absolute_eigenvalue=epsilon;for(int i=0;i<9;i++) if(min_nonzero_absolute_eigenvalue<eigenvalues(i)){min_nonzero_absolute_eigenvalue=eigenvalues(i);break;}
    if(D0.x.x<-epsilon) D0.x.x=eigenvalue_clamp_percentage*min_nonzero_absolute_eigenvalue;else if(D0.x.x<(T)0) D0.x.x=(T)0;
    if(D0.x.y<-epsilon) D0.x.y=eigenvalue_clamp_percentage*min_nonzero_absolute_eigenvalue;else if(D0.x.y<(T)0) D0.x.y=(T)0;
    if(D0.x.z<-epsilon) D0.x.z=eigenvalue_clamp_percentage*min_nonzero_absolute_eigenvalue;else if(D0.x.z<(T)0) D0.x.z=(T)0;
    if(D1.x.x<-epsilon) D1.x.x=eigenvalue_clamp_percentage*min_nonzero_absolute_eigenvalue;else if(D1.x.x<(T)0) D1.x.x=(T)0;
    if(D1.x.y<-epsilon) D1.x.y=eigenvalue_clamp_percentage*min_nonzero_absolute_eigenvalue;else if(D1.x.y<(T)0) D1.x.y=(T)0;
    if(D2.x.x<-epsilon) D2.x.x=eigenvalue_clamp_percentage*min_nonzero_absolute_eigenvalue;else if(D2.x.x<(T)0) D2.x.x=(T)0;
    if(D2.x.y<-epsilon) D2.x.y=eigenvalue_clamp_percentage*min_nonzero_absolute_eigenvalue;else if(D2.x.y<(T)0) D2.x.y=(T)0;
    if(D3.x.x<-epsilon) D3.x.x=eigenvalue_clamp_percentage*min_nonzero_absolute_eigenvalue;else if(D3.x.x<(T)0) D3.x.x=(T)0;
    if(D3.x.y<-epsilon) D3.x.y=eigenvalue_clamp_percentage*min_nonzero_absolute_eigenvalue;else if(D3.x.y<(T)0) D3.x.y=(T)0;
    A0=SYMMETRIC_MATRIX<T,3>::Conjugate(V0,D0);x0000=A0.x00;x1100=A0.x10;x2200=A0.x20;x1111=A0.x11;x2211=A0.x21;x2222=A0.x22;
    A1=SYMMETRIC_MATRIX<T,2>::Conjugate(V1,D1);x1010=A1.x00;x1001=A1.x10;
    A2=SYMMETRIC_MATRIX<T,2>::Conjugate(V2,D2);x2020=A2.x00;x2002=A2.x10;
    A3=SYMMETRIC_MATRIX<T,2>::Conjugate(V3,D3);x2121=A3.x00;x2112=A3.x10;
}
//#####################################################################
// Function Compute_From_Singular_Value_Derivatives
//#####################################################################
template<class T> void DIAGONALIZED_ISOTROPIC_STRESS_DERIVATIVE<T,1>::
Compute_From_Singular_Value_Derivatives(const DIAGONAL_MATRIX<T,1>& F,const VECTOR<T,1>& dE_ds,const SYMMETRIC_MATRIX<T,1>& dE_dsds)
{
    x0000=dE_dsds.x00;
}
//#####################################################################
// Function Compute_From_Singular_Value_Derivatives
//#####################################################################
template<class T> void DIAGONALIZED_ISOTROPIC_STRESS_DERIVATIVE<T,2>::
Compute_From_Singular_Value_Derivatives(const DIAGONAL_MATRIX<T,2>& F,const VECTOR<T,2>& dE_ds,const SYMMETRIC_MATRIX<T,2>& dE_dsds)
{
    T ss1=sqr(F.x.x),ss2=sqr(F.x.y);
    T s01=1/(ss1-ss2);
    
    x0000=dE_dsds.x00;
    x1100=dE_dsds.x10;
    x1111=dE_dsds.x11;
    x1001=(-dE_ds.y*F.x.x+dE_ds.x*F.x.y)*s01;
    x1010=(-dE_ds.y*F.x.y+dE_ds.x*F.x.x)*s01;
}
//#####################################################################
// Function Compute_From_Singular_Value_Derivatives
//#####################################################################
template<class T> void DIAGONALIZED_ISOTROPIC_STRESS_DERIVATIVE<T,3>::
Compute_From_Singular_Value_Derivatives(const DIAGONAL_MATRIX<T,3>& F,const VECTOR<T,3>& dE_ds,const SYMMETRIC_MATRIX<T,3>& dE_dsds)
{
    T ss1=sqr(F.x.x),ss2=sqr(F.x.y),ss3=sqr(F.x.z);
    T s01=1/(ss1-ss2),s02=1/(ss1-ss3),s12=1/(ss2-ss3);
    
    x0000=dE_dsds.x00;
    x1100=dE_dsds.x10;
    x1111=dE_dsds.x11;
    x2200=dE_dsds.x20;
    x2211=dE_dsds.x21;
    x2222=dE_dsds.x22;
    x1001=(-dE_ds.y*F.x.x+dE_ds.x*F.x.y)*s01;
    x1010=(-dE_ds.y*F.x.y+dE_ds.x*F.x.x)*s01;
    x2002=(-dE_ds.z*F.x.x+dE_ds.x*F.x.z)*s02;
    x2020=(-dE_ds.z*F.x.z+dE_ds.x*F.x.x)*s02;
    x2112=(-dE_ds.z*F.x.y+dE_ds.y*F.x.z)*s12;
    x2121=(-dE_ds.z*F.x.z+dE_ds.y*F.x.y)*s12;
}
//#####################################################################
namespace PhysBAM{
template class DIAGONALIZED_ISOTROPIC_STRESS_DERIVATIVE<float,1>;
template class DIAGONALIZED_ISOTROPIC_STRESS_DERIVATIVE<float,2>;
template class DIAGONALIZED_ISOTROPIC_STRESS_DERIVATIVE<float,3>;
template class DIAGONALIZED_ISOTROPIC_STRESS_DERIVATIVE<double,1>;
template class DIAGONALIZED_ISOTROPIC_STRESS_DERIVATIVE<double,2>;
template class DIAGONALIZED_ISOTROPIC_STRESS_DERIVATIVE<double,3>;
}
