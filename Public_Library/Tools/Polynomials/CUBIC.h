//#####################################################################
// Copyright 2002, 2003, Zhaosheng Bao, Ronald Fedkiw.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class CUBIC
//##################################################################### 
#ifndef __CUBIC__
#define __CUBIC__
#include <Tools/Nonlinear_Equations/NONLINEAR_FUNCTION.h>

namespace PhysBAM{

template<class T> class INTERVAL;
template<class T>
class CUBIC
{
public:
    T c3,c2,c1,c0; // coefficients 
    int roots; // number of roots, -1 indicates a=b=c=0 - always a root!
    T root[3]; // root[0] < root[1] < root[2]
    T error_tolerance; // there will be errors in the iterative solver

    T Value(const T x) const
    {return ((c3*x+c2)*x+c1)*x+c0;}

    T operator()(const T x) const
    {return Value(x);}

    T Prime(const T x) const
    {return (3*c3*x+2*c2)*x+c1;}
    
    CUBIC(const T c3_input,const T c2_input,const T c1_input,const T c0_input);
    void Coefficients_From_Interpolation(T x0,T y1,T x1,T y2,T x2,T y3,T x3,T y4);
    void Compute_Roots_Noniterative_In_Interval(const T xmin,const T xmax);
    void Compute_Roots_Noniterative();
    void Compute_Roots();
    void Compute_Roots_In_Interval(const T xmin,const T xmax);
    void Compute_Intervals(const T& a,const T& b,int& intervals,INTERVAL<T> interval[3]);
    void Compute_Intervals(int& intervals,INTERVAL<T> interval[3]);
    void Insert_Root_In_Extrema_Interval(const T xmin,const T xmax);
//#####################################################################
};   
}
#endif
