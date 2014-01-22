//#####################################################################
// Copyright 2003-2008, Eran Guendelman, Geoffrey Irving, Sergey Koltakov.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class RING
//##################################################################### 
#ifndef __RING__
#define __RING__

#include <Geometry/Basic_Geometry/PLANE.h>
namespace PhysBAM{

template<class T>
class RING
{
    typedef VECTOR<T,3> TV;
public:
    typedef int HAS_UNTYPED_READ_WRITE;
    typedef TV VECTOR_T;

    T height,outer_radius,inner_radius;
    PLANE<T> plane1,plane2; // plane2 is height units behind circle

    RING()
        :height(1),outer_radius(2),inner_radius(1),plane2(plane1.x0-height*plane1.normal,-plane1.normal)
    {}

    RING(const TV& X0,const TV& X1,const T outer_radius,const T inner_radius)
        :outer_radius(outer_radius),inner_radius(inner_radius)
    {
        Set_Endpoints(X0,X1);
    }

    void Set_Height(const T height_input)
    {height=height_input;plane2.x0=plane1.x0-height*plane1.normal;}

    void Set_Endpoints(const TV& X0,const TV& X1)
    {TV plane_normal=X0-X1;height=plane_normal.Normalize();
    plane1.x0=X0;plane1.normal=plane_normal;
    plane2.x0=X1;plane2.normal=-plane_normal;}

    template<class RW>
    void Read(std::istream& input)
    {Read_Binary<RW>(input,plane1.x0,plane2.x0,outer_radius,inner_radius);
    Set_Endpoints(plane1.x0,plane2.x0);}

    template<class RW>
    void Write(std::ostream& output) const
    {Write_Binary<RW>(output,plane1.x0,plane2.x0,outer_radius,inner_radius);}

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
