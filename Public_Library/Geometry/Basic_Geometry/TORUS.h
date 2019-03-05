//#####################################################################
// Copyright 2006-2007, Geoffrey Irving, Craig Schroeder, Tamar Shinar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class TORUS
//##################################################################### 
#ifndef __TORUS__
#define __TORUS__

#include <Core/Arrays/ARRAYS_FORWARD.h>
#include <Core/Math_Tools/RANGE.h>
#include <Core/Matrices/SYMMETRIC_MATRIX.h>
#include <Core/Vectors/Dot_Product.h>
#include <string>
namespace PhysBAM{

template<class T>
class TORUS
{
    typedef VECTOR<T,3> TV;
public:
    typedef int HAS_UNTYPED_READ_WRITE;
    typedef TV VECTOR_T;
    TV center,axis;
    T inner_radius,outer_radius;

    TORUS()
        :axis(TV(0,0,1)),inner_radius(1),outer_radius(2)
    {}

    TORUS(const TV& center_input,const TV& axis_input,const T inner_radius_input,
        const T outer_radius_input)
        :center(center_input),axis(axis_input.Normalized()),
        inner_radius(inner_radius_input),outer_radius(outer_radius_input)
    {}

    static std::string Name()
    {return "TORUS<T>";}

    template<class RW>
    void Read(std::istream& input)
    {Read_Binary<RW>(input,center,axis,inner_radius,outer_radius);}

    template<class RW>
    void Write(std::ostream& output) const
    {Write_Binary<RW>(output,center,axis,inner_radius,outer_radius);}

    T Signed_Distance(const TV& X) const;
    TV Normal(const TV& X) const;
    SYMMETRIC_MATRIX<T,3> Hessian(const TV& X) const PHYSBAM_FLATTEN;
    VECTOR<T,2> Principal_Curvatures(const TV& X) const;
    RANGE<TV> Bounding_Box() const;
//#####################################################################
};   
}
#endif
