//#####################################################################
// Copyright 2011, Alexey Stomakhin.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class BOWL
//##################################################################### 
#ifndef __BOWL__
#define __BOWL__
#include <PhysBAM_Geometry/Basic_Geometry/PLANE.h>
namespace PhysBAM{

template<class T>
class BOWL
{
    typedef VECTOR<T,3> TV;

public:
    typedef int HAS_UNTYPED_READ_WRITE;
    typedef TV VECTOR_T;

    T hole_radius,depth,thickness;
    T height,inner_radius,outer_radius;

    BOWL()
        :hole_radius(2),depth(1),thickness(0.5)
    {}

    BOWL(const T hole_radius,const T depth,const T thickness)
        :hole_radius(hole_radius),depth(depth),thickness(thickness)
    {
        height=depth+thickness;
        inner_radius=hole_radius+depth;
        outer_radius=inner_radius+thickness;
    }

    struct HELPER
    {
        T radius;            // from cylindric coordinates 
        TV radial;           // from cylindric coordinates 
        VECTOR<T,2> dX;      // position on vertical slice *relative* to the concentric quarter-circles center
        T dr;
        T signed_distance,c1,c2;
    };

    template<class RW> void Read(std::istream& input)
    {Read_Binary<RW>(input,hole_radius,depth,thickness,height,inner_radius,outer_radius);}

    template<class RW> void Write(std::ostream& output) const
    {Write_Binary<RW>(output,hole_radius,depth,thickness,height,inner_radius,outer_radius);}

//#####################################################################
    RANGE<TV> Bounding_Box() const;
    void Compute_Helper(const TV& X,HELPER& h) const;
    T Signed_Distance(const TV& X) const;
    TV Surface(const TV& X) const;
    TV Surface(const TV& X,const HELPER& h) const;
    TV Normal(const TV& X) const;
    TV Normal(const TV& X,const HELPER& h) const;
    TV Normal(const TV& X,const int aggregate) const;
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
