//#####################################################################
// Copyright 2002-2005, Ronald Fedkiw, Geoffrey Irving.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class RENDERING_CYLINDER
//#####################################################################
#ifndef __RENDERING_CYLINDER__
#define __RENDERING_CYLINDER__

#include <Geometry/Basic_Geometry/CYLINDER.h>
#include <Geometry/Intersections/RAY_CYLINDER_INTERSECTION.h>
#include <Ray_Tracing/Rendering_Objects/RENDERING_OBJECT.h>
namespace PhysBAM{

template<class T>
class RENDERING_CYLINDER:public RENDERING_OBJECT<T>
{
    typedef VECTOR<T,3> TV;
public:
    using RENDERING_OBJECT<T>::small_number;using RENDERING_OBJECT<T>::Object_Space_Point;using RENDERING_OBJECT<T>::Object_Space_Ray;
    using RENDERING_OBJECT<T>::Inside;using RENDERING_OBJECT<T>::Intersection;using RENDERING_OBJECT<T>::World_Space_Vector;
    using RENDERING_OBJECT<T>::World_Space_Point;

    CYLINDER<T> cylinder;

    RENDERING_CYLINDER()
    {}

    RENDERING_CYLINDER(const TV& point1_input,const TV& point2_input,const T radius_input)
        :cylinder(point1_input,point2_input,radius_input)
    {}

    virtual ~RENDERING_CYLINDER()
    {}

    bool Intersection(RAY<TV>& ray) const override
    {RAY<TV> object_space_ray=Object_Space_Ray(ray);
    if(INTERSECTION::Intersects(object_space_ray,cylinder,small_number)){
        ray.semi_infinite=false;ray.t_max=object_space_ray.t_max;ray.aggregate_id=object_space_ray.aggregate_id;return true;}
    else return false;}

    TV Normal(const TV& location,const int aggregate=0) const override
    {return World_Space_Vector(cylinder.Normal(Object_Space_Point(location),aggregate));}

    bool Inside(const TV& location) const override
    {return cylinder.Inside(Object_Space_Point(location),small_number);}

    bool Lazy_Inside(const TV& location) const override
    {return cylinder.Lazy_Inside(Object_Space_Point(location));}
    
    bool Outside(const TV& location) const override
    {return cylinder.Outside(Object_Space_Point(location),small_number);}
    
    bool Lazy_Outside(const TV& location) const override
    {return cylinder.Lazy_Outside(Object_Space_Point(location));}

    bool Boundary(const TV& location) const override
    {return cylinder.Boundary(Object_Space_Point(location),small_number);}

    TV Surface(const TV& location) const override
    {return World_Space_Point(cylinder.Surface(Object_Space_Point(location)));}

    T Signed_Distance(const TV& location) const override
    {return cylinder.Signed_Distance(Object_Space_Point(location));}

//#####################################################################
};   
}
#endif

