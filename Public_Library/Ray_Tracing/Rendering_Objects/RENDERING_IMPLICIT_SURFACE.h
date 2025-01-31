//#####################################################################
// Copyright 2002-2006, Ronald Fedkiw, Geoffrey Irving.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class RENDERING_IMPLICIT_SURFACE
//#####################################################################
#ifndef __RENDERING_IMPLICIT_SURFACE__
#define __RENDERING_IMPLICIT_SURFACE__

#include <Geometry/Implicit_Objects/IMPLICIT_OBJECT.h>
#include <Geometry/Level_Sets/LEVELSET_POLICY.h>
#include <Geometry/Tessellation/IMPLICIT_OBJECT_TESSELLATION.h>
#include <Ray_Tracing/Rendering_Objects/RENDERING_OBJECT.h>
namespace PhysBAM{

template<class TV> struct GRID_ARRAYS_POLICY;
template<class TV> class GRID;

template<class T>
class RENDERING_IMPLICIT_SURFACE:public RENDERING_OBJECT<T>
{
    typedef VECTOR<T,3> TV;typedef VECTOR<int,TV::m> TV_INT;
public:
    using RENDERING_OBJECT<T>::small_number;using RENDERING_OBJECT<T>::inverse_transform;using RENDERING_OBJECT<T>::Object_Space_Ray;
    using RENDERING_OBJECT<T>::Inside;using RENDERING_OBJECT<T>::Intersection;using RENDERING_OBJECT<T>::Object_Space_Point;
    using RENDERING_OBJECT<T>::World_Space_Vector;

    IMPLICIT_OBJECT<TV>* implicit_surface;

    RENDERING_IMPLICIT_SURFACE(IMPLICIT_OBJECT<TV>* implicit_surface_input)
        :implicit_surface(implicit_surface_input)
    {}

    template<class TV>
    RENDERING_IMPLICIT_SURFACE(GRID<TV>& grid_input,ARRAY<T,TV_INT>& phi_input)
    {
        implicit_surface=new LEVELSET_IMPLICIT_OBJECT<TV>(grid_input,phi_input);
    }

    virtual ~RENDERING_IMPLICIT_SURFACE()
    {delete implicit_surface;}

    bool Intersection(RAY<TV>& ray) const override
    {RAY<TV> object_space_ray=Object_Space_Ray(ray);
    if(implicit_surface->Intersection(object_space_ray,small_number)){
        ray.semi_infinite=false;ray.t_max=object_space_ray.t_max;ray.aggregate_id=object_space_ray.aggregate_id;return true;}
    else return false;}

    TV Normal(const TV& location,const int aggregate=0) const override
    {return World_Space_Vector(implicit_surface->Normal(Object_Space_Point(location),aggregate));}

    bool Inside(const TV& location) const override
    {return implicit_surface->Inside(Object_Space_Point(location),small_number);}

    bool Outside(const TV& location) const override
    {return implicit_surface->Outside(Object_Space_Point(location),small_number);}

    bool Boundary(const TV& location) const override
    {return implicit_surface->Boundary(Object_Space_Point(location),small_number);}

    T Signed_Distance(const TV& location) const override
    {return (*implicit_surface)(Object_Space_Point(location));}

    TRIANGULATED_SURFACE<T>* Generate_Triangles() const override
    {TRIANGULATED_SURFACE<T>* surface=TESSELLATION::Generate_Triangles(*implicit_surface);surface->Update_Triangle_List();return surface;}

    RANGE<TV> Object_Space_Bounding_Box() const override
    {return implicit_surface->box;}

//#####################################################################
};
}
#endif
