//#####################################################################
// Copyright 2002-2006, Robert Bridson, Ronald Fedkiw, Geoffrey Irving, Sergey Koltakov, Andrew Selle, Michael Turitzin, Jiayi Chong.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Tools/Arrays/PROJECTED_ARRAY.h>
#include <Tools/Grids_Uniform_Interpolation/LINEAR_INTERPOLATION_UNIFORM.h>
#include <Tools/Images/BMP_FILE.h>
#include <Tools/Images/IMAGE.h>
#include <Tools/Interpolation/LINEAR_INTERPOLATION.h>
#include <Tools/Log/LOG.h>
#include <Tools/Log/PROGRESS_INDICATOR.h>
#include <Geometry/Intersections/RAY_TRIANGLE_3D_INTERSECTION.h>
#include <Geometry/Intersections/RAY_TRIANGULATED_SURFACE_INTERSECTION.h>
#include <Geometry/Topology_Based_Geometry/TRIANGULATED_SURFACE.h>
#include <Ray_Tracing/Rendering_Objects/RENDERING_OBJECT.h>
#include <Ray_Tracing/Rendering_Objects/RENDERING_TRIANGULATED_SURFACE.h>
namespace PhysBAM{
//#####################################################################
// Constructor
//#####################################################################
template<class T> RENDERING_TRIANGULATED_SURFACE<T>::
RENDERING_TRIANGULATED_SURFACE(TRIANGULATED_SURFACE<T>& triangulated_surface_input,
    const int triangles_per_hierarchy_group)
    :triangulated_surface(triangulated_surface_input),base_triangulated_surface(&triangulated_surface),
    closed_volume(true),texture_coordinates(0),triangle_texture_coordinates(0),tangent_vectors(0),
    add_triangles_to_acceleration_structure(true)
{
    triangulated_surface.Update_Bounding_Box();
    if(!triangulated_surface.hierarchy)
        triangulated_surface.Initialize_Hierarchy(true,triangles_per_hierarchy_group);
    if(triangulated_surface.use_vertex_normals && !triangulated_surface.face_vertex_normals)
        triangulated_surface.Update_Vertex_Normals();
}
//#####################################################################
// Destructor
//#####################################################################
template<class T> RENDERING_TRIANGULATED_SURFACE<T>::
~RENDERING_TRIANGULATED_SURFACE()
{
    delete texture_coordinates;
    delete triangle_texture_coordinates;
    delete tangent_vectors;
}
//#####################################################################
// Function Intersection
//#####################################################################
template<class T> bool RENDERING_TRIANGULATED_SURFACE<T>::
Intersection(RAY<TV>& ray) const
{
    RAY<TV> object_space_ray=Object_Space_Ray(ray);
    if(INTERSECTION::Intersects(object_space_ray,*base_triangulated_surface,small_number)){
        ray.semi_infinite=false;ray.t_max=object_space_ray.t_max;ray.aggregate_id=object_space_ray.aggregate_id;return true;}
    else return false;
}
//#####################################################################
// Function Normal
//#####################################################################
template<class T> VECTOR<T,3> RENDERING_TRIANGULATED_SURFACE<T>::
Normal(const TV& location,const int aggregate) const
{
    return World_Space_Vector(triangulated_surface.Normal(Object_Space_Point(location),aggregate));
}
//#####################################################################
// Function Inside
//#####################################################################
template<class T> bool RENDERING_TRIANGULATED_SURFACE<T>::
Inside(const TV& location) const
{
    return closed_volume && triangulated_surface.Inside(Object_Space_Point(location),small_number);
}
//#####################################################################
// Function Outside
//#####################################################################
template<class T> bool RENDERING_TRIANGULATED_SURFACE<T>::
Outside(const TV& location) const
{
    return !closed_volume || triangulated_surface.Outside(Object_Space_Point(location),small_number);
}
//#####################################################################
// Function Boundary
//#####################################################################
template<class T> bool RENDERING_TRIANGULATED_SURFACE<T>::
Boundary(const TV& location) const
{
    return triangulated_surface.Boundary(Object_Space_Point(location),small_number);
}
//#####################################################################
// Constructor
//#####################################################################
template<class T> VECTOR<T,3> RENDERING_TRIANGULATED_SURFACE<T>::
Surface(const TV& location) const
{
    return World_Space_Point(triangulated_surface.Surface(Object_Space_Point(location)));
}
//#####################################################################
// Function Has_Bounding_Box
//#####################################################################
template<class T> bool RENDERING_TRIANGULATED_SURFACE<T>::
Has_Bounding_Box() const 
{
    return true;
}
//#####################################################################
// Function Object_Space_Bounding_Box
//#####################################################################
template<class T> auto RENDERING_TRIANGULATED_SURFACE<T>::
Object_Space_Bounding_Box() const -> RANGE<TV>
{
    if(!triangulated_surface.bounding_box) triangulated_surface.Update_Bounding_Box();
    return *triangulated_surface.bounding_box;
}
//#####################################################################
// Function Closed_Volume
//#####################################################################
template<class T> bool RENDERING_TRIANGULATED_SURFACE<T>::
Closed_Volume() const
{
    return closed_volume;
}
//#####################################################################
// Function Close_To_Open_Surface
//#####################################################################
template<class T> bool RENDERING_TRIANGULATED_SURFACE<T>::
Close_To_Open_Surface(const TV& location,const T threshold_distance) const
{
    assert(!Closed_Volume());int closest_triangle;T distance;
    triangulated_surface.Surface(Object_Space_Point(location),threshold_distance,small_number,&closest_triangle,&distance);
    return distance<=threshold_distance;
}
//#####################################################################
// Function Intersection
//#####################################################################
template<class T> bool RENDERING_TRIANGULATED_SURFACE<T>::
Intersection(RAY<TV>& ray,const int aggregate) const
{
    if(!add_triangles_to_acceleration_structure) return Intersection(ray);
    if(INTERSECTION::Intersects(ray,(world_space_triangles)(aggregate),small_number)){
        ray.aggregate_id=aggregate;return true;}
    else return false;
}
//#####################################################################
// Function Get_Aggregate_World_Space_Bounding_Boxes
//#####################################################################
template<class T> void RENDERING_TRIANGULATED_SURFACE<T>::
Get_Aggregate_World_Space_Bounding_Boxes(ARRAY<RENDERING_OBJECT_ACCELERATION_PRIMITIVE<T> >& primitives) const
{
    if(add_triangles_to_acceleration_structure){
        world_space_triangles.Remove_All();
        for(int i=0;i<triangulated_surface.mesh.elements.m;i++){
            int node1,node2,node3;triangulated_surface.mesh.elements(i).Get(node1,node2,node3);
            TRIANGLE_3D<T> world_space_triangle(transform.Homogeneous_Times(triangulated_surface.particles.X(node1)),transform.Homogeneous_Times(triangulated_surface.particles.X(node2)),transform.Homogeneous_Times(triangulated_surface.particles.X(node3)));
            world_space_triangles.Append(world_space_triangle);
            primitives.Append(RENDERING_OBJECT_ACCELERATION_PRIMITIVE<T>(world_space_triangle.Bounding_Box(),this,i));}}
    else{
        triangulated_surface.Update_Bounding_Box();
        primitives.Append(RENDERING_OBJECT_ACCELERATION_PRIMITIVE<T>(World_Space_Bounding_Box(),this,1));}
}
//#####################################################################
// Function Generate_Triangles
//#####################################################################
template<class T> TRIANGULATED_SURFACE<T>* RENDERING_TRIANGULATED_SURFACE<T>::
Generate_Triangles() const
{
    return triangulated_surface.Create_Compact_Copy();
}
//#####################################################################
// Function Get_Texture_Coordinates
//#####################################################################
template<class T> void RENDERING_TRIANGULATED_SURFACE<T>::
Get_Texture_Coordinates(const TV& object_space_point,const int aggregate,T& s,T& t) const
{
    assert(texture_coordinates);
    int node1,node2,node3;triangulated_surface.mesh.elements(aggregate).Get(node1,node2,node3);
    int uv1,uv2,uv3;(*triangle_texture_coordinates)(aggregate).Get(uv1,uv2,uv3);
    assert(triangulated_surface.mesh.elements.m==triangle_texture_coordinates->m);
    VECTOR<T,3> weights=TRIANGLE_3D<T>::Barycentric_Coordinates(object_space_point,triangulated_surface.particles.X(node1),triangulated_surface.particles.X(node2),triangulated_surface.particles.X(node3));
    VECTOR<T,2> coordinates=(*texture_coordinates)(uv1)*weights.x+(*texture_coordinates)(uv2)*weights.y+(*texture_coordinates)(uv3)*weights.z;
    s=coordinates.x;
    t=coordinates.y;
}
//#####################################################################
// Function Get_Object_Space_Tangent_And_Bitangent
//#####################################################################
template<class T> void RENDERING_TRIANGULATED_SURFACE<T>::
Get_Object_Space_Tangent_And_Bitangent(const TV& object_space_point,const TV& object_space_normal,const int aggregate,TV& object_tangent,TV& object_bitangent) const
{
    assert(tangent_vectors);
    int node1,node2,node3;triangulated_surface.mesh.elements(aggregate).Get(node1,node2,node3);
    VECTOR<T,3> weights=TRIANGLE_3D<T>::Barycentric_Coordinates(object_space_point,triangulated_surface.particles.X(node1),triangulated_surface.particles.X(node2),triangulated_surface.particles.X(node3));
    VECTOR<T,3> interpolated_tangent=(*tangent_vectors)(node1)*weights.x+(*tangent_vectors)(node2)*weights.y+(*tangent_vectors)(node3)*weights.z;
    object_tangent=(interpolated_tangent-object_space_normal*TV::Dot_Product(object_space_normal,interpolated_tangent)).Normalized();
    object_bitangent=-TV::Cross_Product(object_space_normal,object_tangent);
}
//#####################################################################
// Function Get_World_Space_Tangent_And_Bitangent
//#####################################################################
template<class T> void RENDERING_TRIANGULATED_SURFACE<T>::
Get_World_Space_Tangent_And_Bitangent(const TV& world_space_point,const TV& world_space_normal,const int aggregate,TV& world_tangent,TV& world_bitangent) const
{
    TV object_tangent,object_bitangent;
    Get_Object_Space_Tangent_And_Bitangent(Object_Space_Point(world_space_point),Object_Space_Vector(world_space_normal),aggregate,object_tangent,object_bitangent);
    world_tangent=World_Space_Vector(object_tangent);world_bitangent=World_Space_Vector(object_bitangent);
}
//#####################################################################
// Function Compute_Per_Vertex_Tangent_Vectors
//#####################################################################
template<class T> void RENDERING_TRIANGULATED_SURFACE<T>::
Compute_Per_Vertex_Tangent_Vectors()
{
    assert(texture_coordinates&&triangulated_surface.vertex_normals);assert(!tangent_vectors);
    tangent_vectors=new ARRAY<TV>;
    tangent_vectors->Resize(triangulated_surface.mesh.number_nodes);
    for(int i=0;i<triangulated_surface.mesh.elements.m;i++){
        int index1,index2,index3;triangulated_surface.mesh.elements(i).Get(index1,index2,index3);
        int uv1,uv2,uv3;(*triangle_texture_coordinates)(i).Get(uv1,uv2,uv3);
        VECTOR<T,3> p1=triangulated_surface.particles.X(index1),p2=triangulated_surface.particles.X(index2),p3=triangulated_surface.particles.X(index3);
        VECTOR<T,2> st1=(*texture_coordinates)(uv1),st2=(*texture_coordinates)(uv2),st3=(*texture_coordinates)(uv3);
        T denominator=((st2.y-st1.y)*(st3.x-st1.x)-(st3.y-st1.y)*(st2.x-st1.x));
        VECTOR<T,3> dpds=((st2.y-st1.y)*(p3-p1)-(st3.y-st1.y)*(p2-p1))/denominator;
        (*tangent_vectors)(index1)+=dpds;(*tangent_vectors)(index2)+=dpds;(*tangent_vectors)(index3)+=dpds;}
    for(int i=0;i<tangent_vectors->m;i++){
        VECTOR<T,3> cur_normal=(*triangulated_surface.vertex_normals)(i);
        // Gram-Schmidt orthogonalize and normalize each tangent vector
        (*tangent_vectors)(i)=((*tangent_vectors)(i)-cur_normal*TV::Dot_Product(cur_normal,(*tangent_vectors)(i))).Normalized();}
}
//#####################################################################
// Function Read_Texture_Coordinates
//#####################################################################
template<class T> template<class RW> void RENDERING_TRIANGULATED_SURFACE<T>::
Read_Texture_Coordinates(const std::string& filename)
{
    assert(!texture_coordinates);
    texture_coordinates=new ARRAY<VECTOR<T,2> >;
    triangle_texture_coordinates=new ARRAY<VECTOR<int,3> >;
    int backward_compatible;
    FILE_UTILITIES::Read_From_File<RW>(filename,*texture_coordinates,backward_compatible,*triangle_texture_coordinates);
    Compute_Per_Vertex_Tangent_Vectors();
}
//#####################################################################
// Function Rescale_Texture_Coordinates
//#####################################################################
template<class T> void RENDERING_TRIANGULATED_SURFACE<T>::
Rescale_Texture_Coordinates(T scale)
{
    assert(texture_coordinates);
    for(int i=0;i<texture_coordinates->m;i++)
    {T x_rescaled=(*texture_coordinates)(i).x*scale,y_rescaled=(*texture_coordinates)(i).y*scale;
        (*texture_coordinates)(i).x=x_rescaled-floor(x_rescaled);
        (*texture_coordinates)(i).y=y_rescaled-floor(y_rescaled);}
}
//#####################################################################
// Function Initialize_Bump_Map
//#####################################################################
template<class T> void RENDERING_TRIANGULATED_SURFACE<T>::
Initialize_Bump_Map(const std::string& filename)
{
    IMAGE<T>::Read(filename,bump_map_pixels);
    grid.Initialize(bump_map_pixels.Size(),RANGE<VECTOR<T,2> >::Unit_Box());
}
//#####################################################################
// Function Do_Displacement_Map_Per_Vertex
//#####################################################################
template<class T> void RENDERING_TRIANGULATED_SURFACE<T>::
Do_Displacement_Map_Per_Vertex(const T perturb_factor,const T perturb_power)
{
    LOG::cerr<<"Doing per vertex displacement..." << std::endl;
    // TODO: fix to average all vertex uv's
    assert(texture_coordinates->m==triangulated_surface.particles.Size());
    for(int i=0;i<triangulated_surface.mesh.number_nodes;i++){
        VECTOR<T,3> current_perturbation=interpolation.Clamped_To_Array_Cell(grid,bump_map_pixels,VECTOR<T,2>((*texture_coordinates)(i).x,(*texture_coordinates)(i).y));
        current_perturbation=(current_perturbation-TV(0.5,0.5,0.5))*(T)2;
        T current_perturbation_sum=current_perturbation.x+current_perturbation.y+current_perturbation.z;
        triangulated_surface.particles.X(i)+=perturb_factor*std::pow(current_perturbation_sum,perturb_power)*(*triangulated_surface.vertex_normals)(i);}
    triangulated_surface.Initialize_Hierarchy();triangulated_surface.Update_Vertex_Normals();
}
template class RENDERING_TRIANGULATED_SURFACE<double>;
template class RENDERING_TRIANGULATED_SURFACE<float>;
template void RENDERING_TRIANGULATED_SURFACE<double>::Read_Texture_Coordinates<double>(std::string const&);
template void RENDERING_TRIANGULATED_SURFACE<double>::Read_Texture_Coordinates<float>(std::string const&);
template void RENDERING_TRIANGULATED_SURFACE<float>::Read_Texture_Coordinates<float>(std::string const&);
template void RENDERING_TRIANGULATED_SURFACE<float>::Read_Texture_Coordinates<double>(std::string const&);
}
