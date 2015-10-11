//#####################################################################
// Copyright 2002-2007, Geoffrey Irving.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class SPHERE_IMPLICIT_SURFACE
//#####################################################################
#ifndef __SPHERE_IMPLICIT_SURFACE__
#define __SPHERE_IMPLICIT_SURFACE__

#include <Geometry/Basic_Geometry/SPHERE.h>
#include <Geometry/Implicit_Objects/IMPLICIT_OBJECT.h>
namespace PhysBAM{

template<class T_input>
class SPHERE_IMPLICIT_SURFACE:public IMPLICIT_OBJECT<VECTOR<T_input,3> >
{
    typedef T_input T;typedef VECTOR<T,3> TV;
public:
    typedef IMPLICIT_OBJECT<TV> BASE;
    using BASE::box;

    SPHERE<TV> sphere;
protected:
    T minimum_cell_size;
public:

    SPHERE_IMPLICIT_SURFACE(const T minimum_cell_size=0)
        :minimum_cell_size(minimum_cell_size)
    {
        Update_Box();
    }

    SPHERE_IMPLICIT_SURFACE(const TV& center,const T radius,const T minimum_cell_size=0)
        :sphere(center,radius),minimum_cell_size(minimum_cell_size)
    {
        Update_Box();
    }

    void Update_Box() override
    {box.Reset_Bounds(sphere.center);box.Change_Size(sphere.radius);}

    T Minimum_Cell_Size() const override
    {return minimum_cell_size;}

    T Minimum_Cell_Size_Within_Box(const RANGE<TV>& box) const override
    {return minimum_cell_size;}

    T operator()(const TV& location) const override
    {return Extended_Phi(location);}

    T Extended_Phi(const TV& location) const override
    {return sphere.Signed_Distance(location);}

    TV Normal(const TV& location,const int aggregate=-1) const override
    {return sphere.Normal(location);}

    TV Extended_Normal(const TV& location,const int aggregate=-1) const override
    {return sphere.Normal(location);}

    VECTOR<T,2> Principal_Curvatures(const TV& X) const override
    {return sphere.Principal_Curvatures(X);}

    bool Lazy_Inside(const TV& location,const T contour_value=0) const override
    {return Extended_Phi(location)<=contour_value;}

    bool Lazy_Inside_And_Value(const TV& location,T& phi_value,const T contour_value=0) const override
    {phi_value=Extended_Phi(location);return phi_value<=contour_value;}

    bool Lazy_Outside(const TV& location,const T contour_value=0) const override
    {return Extended_Phi(location)>=contour_value;}

    T Min_Phi() const override
    {return -sphere.radius;}

    virtual void Read(TYPED_ISTREAM& input) override {PHYSBAM_FATAL_ERROR();}
    virtual void Write(TYPED_OSTREAM& output) const override {PHYSBAM_FATAL_ERROR();}

//#####################################################################
};
}
#endif
