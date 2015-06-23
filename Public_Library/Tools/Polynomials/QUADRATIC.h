//#####################################################################
// Copyright 2002-2006, Robert Bridson, Ronald Fedkiw, Frederic Gibou, Geoffrey Irving.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class QUADRATIC
//#####################################################################
#ifndef __QUADRATIC__
#define __QUADRATIC__

#include <Tools/Math_Tools/exchange.h>
#include <Tools/Math_Tools/sqr.h>
#include <Tools/Nonlinear_Equations/NONLINEAR_FUNCTION.h>
namespace PhysBAM{

template<class T>
class QUADRATIC:public NONLINEAR_FUNCTION<T(T)>
{
public:
    T a,b,c; // coefficients
    int roots; // number of roots, -1 indicates a=b=c=0 - always a root!
    T root1,root2; // root1 < root2

    QUADRATIC();
    QUADRATIC(const T a_input,const T b_input,const T c_input);
    ~QUADRATIC();

    T Value(const T x) const
    {return (a*x+b)*x+c;}

    T Discriminant() const
    {return sqr(b)-4*a*c;}

    void Compute(const T x,T* ddf,T* df,T* f) const override;
    void Coefficients_From_Interpolation(T x0,T y0,T x1,T y1,T x2,T y2);
    void Compute_Roots();
    void Compute_Roots_In_Interval(const T xmin,const T xmax);
//#####################################################################
};
}
#endif
