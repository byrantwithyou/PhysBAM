//#####################################################################
// Copyright 2007, Andrew Selle, Justin Solomon.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Core/Arrays/PROJECTED_ARRAY.h>
#include <Geometry/Basic_Geometry/CYLINDER.h>
#include <Geometry/Intersections/RAY_CYLINDER_INTERSECTION.h>
#include <Geometry/Topology_Based_Geometry/SEGMENTED_CURVE.h>
#include <Ray_Tracing/Rendering_Objects/RENDERING_OBJECT.h>
#include <Ray_Tracing/Rendering_Objects/RENDERING_SEGMENTED_CURVE.h>
namespace PhysBAM{
//#####################################################################
// Constructor
//#####################################################################
template<class T> RENDERING_SEGMENTED_CURVE<T>::
RENDERING_SEGMENTED_CURVE(SEGMENTED_CURVE<TV>& segmented_curve_input, T thickness_input)
    :segmented_curve(segmented_curve_input),add_segments_to_acceleration_structure(true),thickness(thickness_input)
{
    segmented_curve.Update_Bounding_Box();
    segmented_curve.bounding_box->Change_Size(thickness); // thickness accounts for cylinder size
    if(!segmented_curve.hierarchy) segmented_curve.Initialize_Hierarchy(true); // TODO: when spatial partition, do we need hierachy?
}
//#####################################################################
// Destructor
//#####################################################################
template<class T> RENDERING_SEGMENTED_CURVE<T>::
~RENDERING_SEGMENTED_CURVE()
{
}
//#####################################################################
// Function Intersection
//#####################################################################
template<class T> bool RENDERING_SEGMENTED_CURVE<T>::
Intersection(RAY<TV> &ray) const
{
    RAY<TV> object_space_ray=Object_Space_Ray(ray);
    bool intersected=false;
    for(int t=0;t<segmented_curve.mesh.elements.m;t++){ // TODO: make this more efficient using hierarchy.
        const VECTOR<int,2>& nodes=segmented_curve.mesh.elements(t);
        if(INTERSECTION::Intersects(ray,CYLINDER<T>(segmented_curve.particles.X(nodes[0]),segmented_curve.particles.X(nodes[1]),thickness))) {
            ray.aggregate_id=t;intersected=true;}}
    return intersected;
}
//#####################################################################
// Function Normal
//#####################################################################
template<class T> auto RENDERING_SEGMENTED_CURVE<T>::
Normal(const TV& location,const int aggregate) const -> TV
{
    const VECTOR<int,2>& nodes=segmented_curve.mesh.elements(aggregate);
    return World_Space_Vector(CYLINDER<T>(segmented_curve.particles.X(nodes[0]),segmented_curve.particles.X(nodes[1]),thickness).Normal(Object_Space_Point(location),1));
}
//#####################################################################
// Function Inside
//#####################################################################
template<class T> bool RENDERING_SEGMENTED_CURVE<T>::
Inside(const TV& location) const
{
    return false;
} // segmented curves have no inside
//#####################################################################
// Function Outside
//#####################################################################
template<class T> bool RENDERING_SEGMENTED_CURVE<T>::
Outside(const TV& location) const
{
    return true;
}
//#####################################################################
// Function Boundary
//#####################################################################
template<class T> bool RENDERING_SEGMENTED_CURVE<T>::
Boundary(const TV& location) const
{
    T distance;segmented_curve.Closest_Point_On_Curve(Object_Space_Point(location),0,NULL,&distance);
    return (abs(distance-thickness) < small_number);
}
//#####################################################################
// Function Intersection
//#####################################################################
template<class T> bool RENDERING_SEGMENTED_CURVE<T>::
Intersection(RAY<TV>& ray,const int aggregate) const
{
    if(!add_segments_to_acceleration_structure) return Intersection(ray);
    const VECTOR<int,2>& nodes=segmented_curve.mesh.elements(aggregate);
    if(INTERSECTION::Intersects(ray,CYLINDER<T>(segmented_curve.particles.X(nodes[0]),segmented_curve.particles.X(nodes[1]),thickness))){
        ray.aggregate_id=aggregate;return true;}
    return false;
}
//#####################################################################
// Function Get_Aggregate_World_Space_Bounding_Boxes
//#####################################################################
template<class T> void RENDERING_SEGMENTED_CURVE<T>::
Get_Aggregate_World_Space_Bounding_Boxes(ARRAY<RENDERING_OBJECT_ACCELERATION_PRIMITIVE<T> >& primitives) const
{
    if(add_segments_to_acceleration_structure){
        for(int t=0;t<segmented_curve.mesh.elements.m;t++){
            const VECTOR<int,2>& nodes=segmented_curve.mesh.elements(t);
            RANGE<TV> bounding_box=RANGE<TV>::Bounding_Box(World_Space_Point(segmented_curve.particles.X(nodes[0])),World_Space_Point(segmented_curve.particles.X(nodes[1])));
            bounding_box.Change_Size(thickness);
            primitives.Append(RENDERING_OBJECT_ACCELERATION_PRIMITIVE<T>(bounding_box,this,t));}}
    else{
        primitives.Append(RENDERING_OBJECT_ACCELERATION_PRIMITIVE<T>(Object_Space_Bounding_Box(),this,1));}
}
//#####################################################################
// Function Has_Bounding_Box
//#####################################################################
template<class T> bool RENDERING_SEGMENTED_CURVE<T>::
Has_Bounding_Box() const
{
    return true;
}
//#####################################################################
// Function Object_Space_Bounding_Box
//#####################################################################
template<class T> auto RENDERING_SEGMENTED_CURVE<T>::
Object_Space_Bounding_Box() const -> RANGE<TV>
{
    if(!segmented_curve.bounding_box){
        segmented_curve.Update_Bounding_Box();
        segmented_curve.bounding_box->Change_Size(thickness);} // thickness accounts for cylinder size
    return *segmented_curve.bounding_box;
}
//#####################################################################
// Function Get_Object_Space_Tangent_And_Bitangent
//#####################################################################
template<class T> void RENDERING_SEGMENTED_CURVE<T>::
Get_Object_Space_Tangent_And_Bitangent(const TV& object_space_point,const TV& object_space_normal,const int aggregate,
    TV& object_tangent,TV& object_bitangent) const
{
    const VECTOR<int,2>& nodes=segmented_curve.mesh.elements(aggregate);
    object_tangent=segmented_curve.particles.X(nodes[0])-segmented_curve.particles.X(nodes[1]).Normalized();
    object_bitangent=TV::Cross_Product(object_tangent,object_space_normal).Normalized();
}
//#####################################################################
// Function Get_World_Space_Tangent_And_Bitangent
//#####################################################################
template<class T> void RENDERING_SEGMENTED_CURVE<T>::
Get_World_Space_Tangent_And_Bitangent(const TV& world_space_point,const TV& world_space_normal,const int aggregate,
    TV& world_tangent,TV& world_bitangent) const
{
    const VECTOR<int,2>& nodes=segmented_curve.mesh.elements(aggregate);
    world_tangent=World_Space_Vector(segmented_curve.particles.X(nodes[0])-segmented_curve.particles.X(nodes[1])).Normalized();
    world_bitangent=TV::Cross_Product(world_tangent,world_space_normal).Normalized();
}
template class RENDERING_SEGMENTED_CURVE<double>;
template class RENDERING_SEGMENTED_CURVE<float>;
}
