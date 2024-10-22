//#####################################################################
// Copyright 2003-2007, Zhaosheng Bao, Geoffrey Irving, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class LINE_2D
//##################################################################### 
#ifndef __LINE_2D__
#define __LINE_2D__

#include <Core/Math_Tools/RANGE.h>
#include <Core/Matrices/SYMMETRIC_MATRIX.h>
#include <Core/Vectors/VECTOR_2D.h>
#include <cfloat>
namespace PhysBAM{

template<class T>
class LINE_2D
{
    typedef VECTOR<T,2> TV;
public:
    typedef int HAS_UNTYPED_READ_WRITE;
    typedef TV VECTOR_T;

    TV normal;
    TV x0; // point on the line

    LINE_2D()
        :normal(0,1)
    {}

    LINE_2D(const TV& normal_input,const TV& x1_input)
        :normal(normal_input),x0(x1_input)
    {}

    TV Normal(const TV&) const
    {return normal;}

    RANGE<TV> Bounding_Box() const
    {return RANGE<TV>::Full_Box();}

    T Signed_Distance(const TV& location) const
    {return TV::Dot_Product(normal,location-x0);}

    SYMMETRIC_MATRIX<T,2> Hessian(const TV& X) const
    {return SYMMETRIC_MATRIX<T,2>();}

    // inside is the half space behind the normal
    bool Inside(const TV& location,const T thickness_over_two) const
    {return Signed_Distance(location)<=-thickness_over_two;}

    bool Lazy_Inside(const TV& location) const
    {return Signed_Distance(location)<=0;}

    bool Outside(const TV& location,const T thickness_over_two) const
    {return !Inside(location,-thickness_over_two);}

    bool Lazy_Outside(const TV& location) const
    {return !Lazy_Inside(location);}

    bool Segment_Intersection(const TV& endpoint1,const TV& endpoint2,T& interpolation_fraction) const
    {return Segment_Line_Intersection(endpoint1,endpoint2,interpolation_fraction);}

    bool Segment_Line_Intersection(const TV& endpoint1,const TV& endpoint2,T& interpolation_fraction) const 
    {T denominator=TV::Dot_Product(endpoint2-endpoint1,normal);
    if(!denominator){interpolation_fraction=FLT_MAX;return false;} // parallel
    interpolation_fraction=TV::Dot_Product(x0-endpoint1,normal)/denominator;
    return interpolation_fraction>=0 && interpolation_fraction<=1;}

    VECTOR<T,1> Principal_Curvatures(const TV& X) const
    {return VECTOR<T,1>();}

    static std::string Name()
    {return "LINE_2D<T>";}

    template<class RW> void Read(std::istream& input)
    {Read_Binary<RW>(input,normal,x0);}

    template<class RW> void Write(std::ostream& output) const
    {Write_Binary<RW>(output,normal,x0);}
//##################################################################### 
};
}
#endif
