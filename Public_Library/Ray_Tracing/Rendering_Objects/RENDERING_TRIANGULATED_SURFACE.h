//#####################################################################
// Copyright 2002-2006, Robert Bridson, Ronald Fedkiw, Geoffrey Irving, Sergey Koltakov, Andrew Selle, Michael Turitzin, Jiayi Chong.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class RENDERING_TRIANGLULATED_SURFACE
//#####################################################################
#ifndef __RENDERING_TRIANGULATED_SURFACE__
#define __RENDERING_TRIANGULATED_SURFACE__

#include <Tools/Arrays/PROJECTED_ARRAY.h>
#include <Tools/Grids_Uniform_Interpolation/LINEAR_INTERPOLATION_UNIFORM.h>
#include <Tools/Images/BMP_FILE.h>
#include <Tools/Interpolation/LINEAR_INTERPOLATION.h>
#include <Tools/Log/LOG.h>
#include <Tools/Log/PROGRESS_INDICATOR.h>
#include <Geometry/Intersections/RAY_TRIANGLE_3D_INTERSECTION.h>
#include <Geometry/Intersections/RAY_TRIANGULATED_SURFACE_INTERSECTION.h>
#include <Geometry/Topology_Based_Geometry/TRIANGULATED_SURFACE.h>
#include <Ray_Tracing/Rendering_Objects/RENDERING_OBJECT.h>
namespace PhysBAM{

template<class T> class RENDERING_BSSRDF_SHADER;

template<class T>
class RENDERING_TRIANGULATED_SURFACE:public RENDERING_OBJECT<T>
{
    typedef VECTOR<T,3> TV;
public:
    using RENDERING_OBJECT<T>::small_number;using RENDERING_OBJECT<T>::transform;using RENDERING_OBJECT<T>::material_shader;
    using RENDERING_OBJECT<T>::flip_normal;using RENDERING_OBJECT<T>::World_Space_Bounding_Box;using RENDERING_OBJECT<T>::name;using RENDERING_OBJECT<T>::two_sided;
    using RENDERING_OBJECT<T>::Inside;using RENDERING_OBJECT<T>::Intersection;using RENDERING_OBJECT<T>::Object_Space_Ray;using RENDERING_OBJECT<T>::Object_Space_Point;
    using RENDERING_OBJECT<T>::Object_Space_Vector;using RENDERING_OBJECT<T>::World_Space_Vector;using RENDERING_OBJECT<T>::World_Space_Point;

    TRIANGULATED_SURFACE<T>& triangulated_surface;
    mutable TRIANGULATED_SURFACE<T>* base_triangulated_surface;
    mutable ARRAY<TRIANGLE_3D<T> > world_space_triangles;
    bool closed_volume;
    ARRAY<TV,VECTOR<int,2> > bump_map_pixels;
    LINEAR_INTERPOLATION_UNIFORM<VECTOR<T,2>,TV> interpolation;
    GRID<VECTOR<T,2> > grid;
    ARRAY<VECTOR<T,2> >* texture_coordinates;
    ARRAY<VECTOR<int,3> >* triangle_texture_coordinates;
    ARRAY<TV>* tangent_vectors;
    std::string sample_locations_file;
    bool add_triangles_to_acceleration_structure;

    RENDERING_TRIANGULATED_SURFACE(TRIANGULATED_SURFACE<T>& triangulated_surface_input,
        const int triangles_per_hierarchy_group=0);
    virtual ~RENDERING_TRIANGULATED_SURFACE();
    bool Intersection(RAY<TV>& ray) const override;
    TV Normal(const TV& location,const int aggregate=0) const override;
    bool Inside(const TV& location) const override;
    bool Outside(const TV& location) const override;
    bool Boundary(const TV& location) const override;
    TV Surface(const TV& location) const override;
    bool Has_Bounding_Box() const  override;
    RANGE<TV> Object_Space_Bounding_Box() const override;
    bool Closed_Volume() const override;
    bool Close_To_Open_Surface(const TV& location,const T threshold_distance) const override;
    bool Intersection(RAY<TV>& ray,const int aggregate) const override;
    void Get_Aggregate_World_Space_Bounding_Boxes(
        ARRAY<RENDERING_OBJECT_ACCELERATION_PRIMITIVE<T> >& primitives) const override;
    TRIANGULATED_SURFACE<T>* Generate_Triangles() const override;
    void Get_Texture_Coordinates(const TV& object_space_point,const int aggregate,T& s,T& t) const override;
    void Get_Object_Space_Tangent_And_Bitangent(const TV& object_space_point,const TV& object_space_normal,
        const int aggregate,TV& object_tangent,TV& object_bitangent) const override;
    void Get_World_Space_Tangent_And_Bitangent(const TV& world_space_point,const TV& world_space_normal,
        const int aggregate,TV& world_tangent,TV& world_bitangent) const override;
    void Compute_Per_Vertex_Tangent_Vectors();
    template<class RW>
    void Read_Texture_Coordinates(const std::string& filename);
    void Rescale_Texture_Coordinates(T scale);
    void Initialize_Bump_Map(const std::string& filename);
    void Do_Displacement_Map_Per_Vertex(const T perturb_factor,const T perturb_power);
//#####################################################################
};
}
#endif
