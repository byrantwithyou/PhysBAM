//#####################################################################
// Copyright 2003-2008, Eran Guendelman, Geoffrey Irving, Sergey Koltakov.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class CONE
//##################################################################### 
#ifndef __CONE__
#define __CONE__

#include <Geometry/Basic_Geometry/PLANE.h>
namespace PhysBAM{

template<class T>
class CONE
{
    typedef VECTOR<T,3> TV;
public:
    typedef int HAS_UNTYPED_READ_WRITE;
    typedef TV VECTOR_T;

    TV base,dir;
    T radius,height;

    CONE()
        :dir(1,0,0),radius(1),height(1)
    {}

    CONE(const TV& bottom,const TV& top,const T bottom_radius)
        :base(bottom),dir(top-bottom),radius(bottom_radius),height(dir.Normalize())
    {}

    template<class RW>
    void Read(std::istream& input)
    {Read_Binary<RW>(input,base,dir,radius,height);}

    template<class RW>
    void Write(std::ostream& output) const
    {Write_Binary<RW>(output,base,dir,radius,height);}

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
