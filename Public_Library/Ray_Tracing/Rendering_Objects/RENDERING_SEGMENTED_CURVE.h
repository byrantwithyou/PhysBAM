//#####################################################################
// Copyright 2007, Andrew Selle, Justin Solomon.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class RENDERING_SEGMENTED_CURVE
//#####################################################################
#ifndef __RENDERING_SEGMENTED_CURVE__
#define __RENDERING_SEGMENTED_CURVE__

#include <Tools/Arrays/PROJECTED_ARRAY.h>
#include <Geometry/Basic_Geometry/CYLINDER.h>
#include <Geometry/Intersections/RAY_CYLINDER_INTERSECTION.h>
#include <Geometry/Topology_Based_Geometry/SEGMENTED_CURVE.h>
#include <Ray_Tracing/Rendering_Objects/RENDERING_OBJECT.h>
namespace PhysBAM{

template<class T>
class RENDERING_SEGMENTED_CURVE:public RENDERING_OBJECT<T>
{
    typedef VECTOR<T,3> TV;
public:
    using RENDERING_OBJECT<T>::small_number;using RENDERING_OBJECT<T>::World_Space_Vector;using RENDERING_OBJECT<T>::World_Space_Point;
    using RENDERING_OBJECT<T>::World_Space_Bounding_Box;using RENDERING_OBJECT<T>::Object_Space_Ray;using RENDERING_OBJECT<T>::Object_Space_Point;
    using RENDERING_OBJECT<T>::Inside;using RENDERING_OBJECT<T>::Intersection; // silence -Woverloaded-virtual

    SEGMENTED_CURVE<TV>& segmented_curve;
    bool add_segments_to_acceleration_structure;
    T thickness;

    RENDERING_SEGMENTED_CURVE(SEGMENTED_CURVE<TV>& segmented_curve_input, T thickness_input=1e-4);
    virtual ~RENDERING_SEGMENTED_CURVE();

    bool Intersection(RAY<TV> &ray) const override;
    TV Normal(const TV& location,const int aggregate=0) const override;
    bool Inside(const TV& location) const override; // segmented curves have no inside
    bool Outside(const TV& location) const override;
    bool Boundary(const TV& location) const override;
    bool Intersection(RAY<TV>& ray,const int aggregate) const override;
    void Get_Aggregate_World_Space_Bounding_Boxes(ARRAY<RENDERING_OBJECT_ACCELERATION_PRIMITIVE<T> >& primitives) const override;
    bool Has_Bounding_Box() const override;
    RANGE<TV> Object_Space_Bounding_Box() const override;
    void Get_Object_Space_Tangent_And_Bitangent(const TV& object_space_point,const TV& object_space_normal,const int aggregate,TV& object_tangent,TV& object_bitangent) const override;
    void Get_World_Space_Tangent_And_Bitangent(const TV& world_space_point,const TV& world_space_normal,const int aggregate,TV& world_tangent,TV& world_bitangent) const override;
//#####################################################################
};
}
#endif
