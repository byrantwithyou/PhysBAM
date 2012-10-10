//#####################################################################
// Copyright 2002, Ronald Fedkiw.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class NEWTON_POLYNOMIAL
//#####################################################################
#ifndef __NEWTON_POLYNOMIAL__
#define __NEWTON_POLYNOMIAL__

#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Nonlinear_Equations/NONLINEAR_FUNCTION.h>
namespace PhysBAM{

template<class T>
class NEWTON_POLYNOMIAL:public NONLINEAR_FUNCTION<T(T)>
{
public:
    int degree;
    ARRAY<T> c; // coefficients
    ARRAY<T> x; // data points

    NEWTON_POLYNOMIAL(const int degree_input)
        :degree(degree_input),c(degree+1),x(degree+1)
    {}

    T operator()(const T x_input) const PHYSBAM_OVERRIDE
    {T y=c(degree);for(int i=degree-1;i>=0;i--) y=c(i)+(x_input-x(i))*y;return y;}

    void Compute_Coefficients(const ARRAY<T>& x_input,const ARRAY<T>& y_input)
    {int i,j;assert(x_input.n == degree+1 && y_input.n == degree+1);for(i=0;i<degree+1;i++) x(i)=x_input(i);
    ARRAY<T> f(degree+1);for(i=0;i<degree+1;i++) f(i)=y_input(i); c(0)=f(0);
    for(j=0;j<degree;j++){for(i=0;i<degree+1-j;i++) f(i)=(f(i+1)-f(i))/(x(i+1+j-1)-x(i)); c(j)=f(0);}}

//#####################################################################
};
}
#endif
