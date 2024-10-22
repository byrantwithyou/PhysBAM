//#####################################################################
// Copyright 2002-2004, Ronald Fedkiw, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class BURGERS_1D_EIGENSYSTEM_F  
//##################################################################### 
#include <Core/Arrays_Nd/ARRAYS_ND.h>
#include <Core/Matrices/MATRIX.h>
#include <Compressible/Burgers_Equation/BURGERS_1D_EIGENSYSTEM_F.h>
using namespace PhysBAM;
//#####################################################################
// Function Flux
//#####################################################################
// F(U) for i in (-2,m+3) - 3 ghost cells
template<class T> void BURGERS_1D_EIGENSYSTEM_F<T>::
Flux(const int m,const ARRAY<TV_DIMENSION,VECTOR<int,1> >& U,ARRAY<TV_DIMENSION,VECTOR<int,1> >& F,ARRAY<TV_DIMENSION,VECTOR<int,1> >* U_clamped)       
{
    for(int i=-3;i<m+3;i++) F(i)(0)=sqr(U(i)(0)); // u^2
}
//#####################################################################
// Function Eigenvalues
//#####################################################################
// eigenvalues for F(U) at flux i and and at points i and i+1
template<class T> bool BURGERS_1D_EIGENSYSTEM_F<T>::
Eigenvalues(const ARRAY<TV_DIMENSION,VECTOR<int,1> >& U,const int i,VECTOR<T,d>& lambda,VECTOR<T,d>& lambda_left,VECTOR<T,d>& lambda_right)
{
    // eigenvalues on the left - at point i
    lambda_left(0)=U(i)(0);      
    // eigenvalues on the right - at point i+1
    lambda_right(0)=U(i+1)(0);
    // eigenvalues in the center - at flux i
    lambda(0)=(U(i)(0)+U(i+1)(0))/2;

    return true; // eigensystem is well defined
}  
//#####################################################################
// Function Eigenvectors
//#####################################################################
// no eigensystem for the scalar
template<class T> void BURGERS_1D_EIGENSYSTEM_F<T>::
Eigenvectors(const ARRAY<TV_DIMENSION,VECTOR<int,1> >& U,const int i,MATRIX<T,d,d>& L,MATRIX<T,d,d>& R)
{
    L(0,0)=1;
    R(0,0)=1;
}  
//#####################################################################
namespace PhysBAM{
template class BURGERS_1D_EIGENSYSTEM_F<float>;
template class BURGERS_1D_EIGENSYSTEM_F<double>;
}
