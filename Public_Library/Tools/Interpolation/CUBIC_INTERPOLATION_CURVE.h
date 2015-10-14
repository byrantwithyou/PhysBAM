//#####################################################################
// Copyright 2007, Craig Schroeder.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class CUBIC_INTERPOLATION_CURVE
//#####################################################################
#ifndef __CUBIC_INTERPOLATION_CURVE__
#define __CUBIC_INTERPOLATION_CURVE__

#include <Tools/Arrays/ARRAY.h>
#include <Tools/Vectors/VECTOR_FORWARD.h>
namespace PhysBAM{

template<class T,class T2>
class CUBIC_INTERPOLATION_CURVE
{
    struct CONTROL_POINT
    {
        CONTROL_POINT()
        {}

        CONTROL_POINT(T time,T2 value)
        {t=time;x=value;}

        T t;T2 x;
        mutable T2 b,c,d;
    };
    ARRAY<CONTROL_POINT> control_points;
    mutable bool need_compute_coefficients;
public:

    CUBIC_INTERPOLATION_CURVE();
    int Locate_Interval(const T t) const;
    T2 Value(T t) const;
    T2 Derivative(T t) const;
    T2 Second_Derivative(T t) const;
    void Add_Control_Point(T t,const T2& value);
    void Compute_Coefficients() const;
//#####################################################################
};
}
#endif
