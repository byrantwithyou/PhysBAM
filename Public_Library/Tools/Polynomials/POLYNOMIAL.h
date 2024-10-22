//#####################################################################
// Copyright 2002-2005, Ronald Fedkiw, Geoffrey Irving.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class POLYNOMIAL
//#####################################################################
#ifndef __POLYNOMIAL__
#define __POLYNOMIAL__

#include <Core/Arrays/ARRAY.h>
#include <Core/Matrices/MATRIX_MXN.h>
#include <Tools/Nonlinear_Equations/NONLINEAR_FUNCTION.h>
#include <cassert>
namespace PhysBAM{

template<class T>
class POLYNOMIAL:public NONLINEAR_FUNCTION<T(T)>
{
public:
    int degree;
    ARRAY<T> c; // coefficients

    POLYNOMIAL(const int degree_input)
        :degree(degree_input),c(degree+1)
    {}

    T operator()(const T x_input) const override
    {T y=c(degree);for(int i=degree-1;i>=0;i--) y=c(i)+x_input*y;return y;}

    void Compute_Coefficients(const ARRAY<T>& x,const ARRAY<T>& y)
    {assert(x.n == degree+1 && y.n == degree+1);MATRIX_MXN<T> A(degree+1,degree+1);
    for(int i=0;i<degree+1;i++){A(i,1)=1;for(int j=1;j<degree+1;j++)A(i,j)=x(i)*A(i,j-1);}
    ARRAY<T> c_temp=A.PLU_Solve(y);for(int i=0;i<degree+1;i++)c(i)=c_temp(i);}

//#####################################################################
};
}
#endif
