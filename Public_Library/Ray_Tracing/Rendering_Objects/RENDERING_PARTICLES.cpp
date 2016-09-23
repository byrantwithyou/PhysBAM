//#####################################################################
// Copyright 2002-2005, Tamar Shinar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################/
#include <Grid_Tools/Grids/GRID.h>
#include <Geometry/Basic_Geometry/SPHERE.h>
#include <Geometry/Geometry_Particles/GEOMETRY_PARTICLES.h>
#include <Geometry/Intersections/RAY_SPHERE_INTERSECTION.h>
#include <Ray_Tracing/Rendering_Objects/RENDERING_OBJECT.h>
#include <Ray_Tracing/Rendering_Objects/RENDERING_PARTICLES.h>
namespace PhysBAM {
//#####################################################################
// Constructor
//#####################################################################
template<class T> RENDERING_PARTICLES<T>::
RENDERING_PARTICLES(const ARRAY<GEOMETRY_PARTICLES<TV>*,VECTOR<int,3> >& p,const GRID<TV>& g,T s)
    :particles_array(p),grid(g),scale(s),
    particle_to_aggregate_id(particles_array.domain.min_corner.x,particles_array.domain.max_corner.x,particles_array.domain.min_corner.y,particles_array.domain.max_corner.y,particles_array.domain.min_corner.z,particles_array.domain.max_corner.z,true),
    aggregate_id_to_particle()
{
    Create_Aggregate_Ids();
}
//#####################################################################
// Destructor
//#####################################################################
template<class T> RENDERING_PARTICLES<T>::
~RENDERING_PARTICLES()
{
}
//#####################################################################
// Function Intersection
//#####################################################################
template<class T> bool RENDERING_PARTICLES<T>::
Intersection(RAY<TV>& ray) const
{
    RAY<TV> object_space_ray=Object_Space_Ray(ray);
    bool intersection=false;
    for(int i=particles_array.domain.min_corner.x;i<particles_array.domain.max_corner.x;i++)
        for(int j=particles_array.domain.min_corner.y;j<particles_array.domain.max_corner.y;j++)
            for(int ij=particles_array.domain.min_corner.z;ij<particles_array.domain.max_corner.z;ij++){
                GEOMETRY_PARTICLES<TV>* particles=particles_array(i,j,ij);
                if(particles){
                    ARRAY_VIEW<T> radius=*particles->template Get_Array<T>(ATTRIBUTE_ID(15)); // radius is attribute 15
                    for(int p=0;p<particles->Size();p++)
                        if(INTERSECTION::Intersects(object_space_ray,SPHERE<TV>(particles->X(p),scale*radius(p)),small_number)){
                            object_space_ray.aggregate_id=particle_to_aggregate_id(i,j,ij)(p);intersection=true;}}}
    if(intersection){ray.t_max=object_space_ray.t_max;ray.semi_infinite=object_space_ray.semi_infinite;ray.aggregate_id=object_space_ray.aggregate_id;}
    return intersection;
}
//#####################################################################
// Function Intersection
//#####################################################################
template<class T> bool RENDERING_PARTICLES<T>::
Intersection(RAY<TV>& ray,const int aggregate) const
{
    RAY<TV> object_space_ray=Object_Space_Ray(ray);bool intersection=false;
    GEOMETRY_PARTICLES<TV>* particles=particles_array(aggregate_id_to_particle(aggregate).x);int p=aggregate_id_to_particle(aggregate).y;
    ARRAY_VIEW<T> radius=*particles->template Get_Array<T>(ATTRIBUTE_ID(15)); // radius is attribute 15
    if(INTERSECTION::Intersects(object_space_ray,(SPHERE<TV>(particles->X(p),scale*radius(p))),small_number)){
        ray.t_max=object_space_ray.t_max;ray.semi_infinite=object_space_ray.semi_infinite;ray.aggregate_id=aggregate;intersection=true;}
    return intersection;
}
//#####################################################################
// Function Inside
//#####################################################################
template<class T> bool RENDERING_PARTICLES<T>::
Inside(const TV& location) const
{
    for(int i=particles_array.domain.min_corner.x;i<particles_array.domain.max_corner.x;i++)
        for(int j=particles_array.domain.min_corner.y;j<particles_array.domain.max_corner.y;j++)
            for(int ij=particles_array.domain.min_corner.z;ij<particles_array.domain.max_corner.z;ij++){
                GEOMETRY_PARTICLES<TV>* particles=particles_array(i,j,ij);
                if(particles){
                    ARRAY_VIEW<T> radius=*particles->template Get_Array<T>(ATTRIBUTE_ID(15)); // radius is attribute 15
                    for(int p=0;p<particles->Size();p++)
                        if(SPHERE<TV>(particles->X(p),scale*radius(p)).Inside(Object_Space_Point(location),small_number))
                            return true;}}
    return false;
}
//#####################################################################
// Function Normal
//#####################################################################
template<class T> auto RENDERING_PARTICLES<T>::
Normal(const TV& location,const int aggregate) const -> TV
{
    GEOMETRY_PARTICLES<TV>* particles=particles_array(aggregate_id_to_particle(aggregate).x);int p=aggregate_id_to_particle(aggregate).y;
    ARRAY_VIEW<T> radius=*particles->template Get_Array<T>(ATTRIBUTE_ID(15)); // radius is attribute 15
    return (SPHERE<TV>(particles->X(p),scale*radius(p))).Normal(location);
}
//#####################################################################
// Function Get_Aggregate_World_Space_Bounding_Boxes
//#####################################################################
template<class T> void RENDERING_PARTICLES<T>::
Get_Aggregate_World_Space_Bounding_Boxes(ARRAY<RENDERING_OBJECT_ACCELERATION_PRIMITIVE<T> >& primitives) const
{
    for(int i=particles_array.domain.min_corner.x;i<particles_array.domain.max_corner.x;i++)
        for(int j=particles_array.domain.min_corner.y;j<particles_array.domain.max_corner.y;j++)
            for(int ij=particles_array.domain.min_corner.z;ij<particles_array.domain.max_corner.z;ij++){
                GEOMETRY_PARTICLES<TV>* particles=particles_array(i,j,ij);
                if(particles){
                    ARRAY_VIEW<T> radius=*particles->template Get_Array<T>(ATTRIBUTE_ID(15)); // radius is attribute 15
                    for(int p=0;p<particles->Size();p++){
                        T radius_scalar=scale*radius(p);
                        TV radius_vector(radius_scalar,radius_scalar,radius_scalar),world_center=transform.Homogeneous_Times(particles->X(p));
                        RANGE<TV> world_space_bounding_box(world_center-radius_vector,world_center+radius_vector);
                        int aggregate_id=particle_to_aggregate_id(i,j,ij)(p);
                        primitives.Append(RENDERING_OBJECT_ACCELERATION_PRIMITIVE<T>(world_space_bounding_box,this,aggregate_id));}}}
}
//#####################################################################
// Function Create_Aggregate_Ids
//#####################################################################
template<class T> void RENDERING_PARTICLES<T>::
Create_Aggregate_Ids()
{
    int index=0;
    for(int i=particles_array.domain.min_corner.x;i<particles_array.domain.max_corner.x;i++)
        for(int j=particles_array.domain.min_corner.y;j<particles_array.domain.max_corner.y;j++)
            for(int ij=particles_array.domain.min_corner.z;ij<particles_array.domain.max_corner.z;ij++){
                GEOMETRY_PARTICLES<TV>* particles=particles_array(i,j,ij);
                if(particles)
                    for(int p=0;p<particles->Size();p++){
                        particle_to_aggregate_id(i,j,ij).Append(index);aggregate_id_to_particle.Append(PAIR<VECTOR<int,3>,int>(VECTOR<int,3>(i,j,ij),p));index++;}}
}
template class RENDERING_PARTICLES<double>;
template class RENDERING_PARTICLES<float>;
}
