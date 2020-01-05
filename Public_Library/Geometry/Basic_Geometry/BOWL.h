//#####################################################################
// Copyright 2011, Alexey Stomakhin.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class BOWL
//##################################################################### 
#ifndef __BOWL__
#define __BOWL__
#include <Core/Matrices/FRAME.h>
namespace PhysBAM{

template<class T>
class BOWL
{
    typedef VECTOR<T,3> TV;

public:
    typedef int HAS_UNTYPED_READ_WRITE;
    typedef TV VECTOR_T;

    FRAME<TV> frame;
    T hole_radius,depth,thickness;
    T height,inner_radius,outer_radius;
    
    BOWL()
        :BOWL(TV(),TV(0,1,0),2,1,(T).5)
    {}

    BOWL(const TV& location,const TV& axis,const T hole_radius,const T depth,const T thickness)
        :hole_radius(hole_radius),depth(depth),thickness(thickness)
    {
        frame.r=ROTATION<TV>::From_Rotated_Vector(axis,TV(0,-1,0));
        frame.t=TV(0,depth,0)-frame.r.Rotate(location);
        LOG::printf("%P\n",frame);
        height=depth+thickness;
        inner_radius=hole_radius+depth;
        outer_radius=inner_radius+thickness;
    }

    template<class RW> void Read(std::istream& input)
    {Read_Binary<RW>(input,frame,hole_radius,depth,thickness,height,inner_radius,outer_radius);}

    template<class RW> void Write(std::ostream& output) const
    {Write_Binary<RW>(output,frame,hole_radius,depth,thickness,height,inner_radius,outer_radius);}

//#####################################################################
    RANGE<TV> Bounding_Box() const;
    T Signed_Distance(const TV& X) const;
    TV Surface(const TV& X) const;
    TV Normal(const TV& X) const;
    TV Normal(const TV& X,const int aggregate) const;
    SYMMETRIC_MATRIX<T,3> Hessian(const TV& X) const PHYSBAM_FLATTEN;
    VECTOR<T,2> Principal_Curvatures(const TV& X) const;
    bool Lazy_Inside(const TV& X) const;
    bool Lazy_Outside(const TV& X) const;
    bool Inside(const TV& X,const T thickness_over_two=0) const;
    bool Outside(const TV& X,const T thickness_over_two=0) const;
    bool Boundary(const TV& X,const T thickness_over_two) const;
    static std::string Name();
//#####################################################################
};   
}
#endif
