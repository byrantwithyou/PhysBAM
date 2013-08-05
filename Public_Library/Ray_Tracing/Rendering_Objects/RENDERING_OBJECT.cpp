//#####################################################################
// Copyright 2002-2007, Jiayi Chong, Ronald Fedkiw, Eran Guendelman, Igor Neverov, Andrew Selle, Mike Turitzin.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Ray_Tracing/Rendering_Objects/RENDERING_OBJECT.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class T> RENDERING_OBJECT<T>::
RENDERING_OBJECT()
    :add_to_spatial_partition(false),index_of_refraction(1),priority(0),support_transparent_overlapping_objects(false),material_shader(0),volumetric_shader(0),two_sided(true),
    flip_normal(false),
    bssrdf_shader(0)
{
    small_number=Default_Small_Number();
    transform=MATRIX<T,4>::Identity_Matrix();inverse_transform=MATRIX<T,4>::Identity_Matrix();solid_texture_transform=MATRIX<T,4>::Identity_Matrix();
}
//#####################################################################
// Destructor
//#####################################################################
template<class T> RENDERING_OBJECT<T>::
~RENDERING_OBJECT()
{
}
//#####################################################################
// Function Index_Of_Refraction
//#####################################################################
template<class T> T RENDERING_OBJECT<T>::
Index_Of_Refraction(const TV& world_space_location) const
{
    return index_of_refraction;
}
//#####################################################################
// Function Get_Texture_Coordinates
//#####################################################################
template<class T> void RENDERING_OBJECT<T>::
Get_Texture_Coordinates(const TV& object_space_point,const int aggregate,T& s,T& t) const
{
    PHYSBAM_FUNCTION_IS_NOT_DEFINED();
}
//#####################################################################
// Function Get_Solid_Texture_Coordinates
//#####################################################################
template<class T> VECTOR<T,3> RENDERING_OBJECT<T>::
Get_Solid_Texture_Coordinates(const TV& object_space_point,const int aggregate) const
{
    return solid_texture_transform.Homogeneous_Times(object_space_point);
}
//#####################################################################
// Function Get_World_Space_Tangent_And_Bitangent
//#####################################################################
template<class T> void RENDERING_OBJECT<T>::
Get_World_Space_Tangent_And_Bitangent(const TV& world_space_point,const TV& world_space_normal,const int aggregate,TV& world_tangent,TV& world_bitangent) const
{
    PHYSBAM_FUNCTION_IS_NOT_DEFINED();
}
//#####################################################################
// Function Get_Object_Space_Tangent_And_Bitangent
//#####################################################################
template<class T> void RENDERING_OBJECT<T>::
Get_Object_Space_Tangent_And_Bitangent(const TV& object_space_point,const TV& object_space_normal,const int aggregate,TV& object_tangent,TV& object_bitangent) const
{
    PHYSBAM_FUNCTION_IS_NOT_DEFINED();
}
//#####################################################################
// Function Intersection
//#####################################################################
template<class T> bool RENDERING_OBJECT<T>::
Intersection(RAY<TV>& ray,const int lowest_priority,const RENDERING_OBJECT<T>** intersected_object) const
{
    if(priority>=lowest_priority&&Intersection(ray)){*intersected_object=this;return true;}
    else return false;
}
//#####################################################################
// Function Inside
//#####################################################################
template<class T> bool RENDERING_OBJECT<T>::
Inside(const TV& location,RENDERING_OBJECT<T>** intersected_object) const
{
    if(support_transparent_overlapping_objects && Inside(location)){
        *intersected_object=const_cast<RENDERING_OBJECT<T>*>(this);
        return true;}
    return false;
}
//#####################################################################
// Function Closed_Volume
//#####################################################################
template<class T> bool RENDERING_OBJECT<T>::
Closed_Volume() const // indicates whether Inside/Outside are meaningful for this object
{
    return true;
}
//#####################################################################
// Function Intersection
//#####################################################################
template<class T> bool RENDERING_OBJECT<T>::
Intersection(RAY<TV>& ray) const
{
    PHYSBAM_FUNCTION_IS_NOT_DEFINED();
}
//#####################################################################
// Function Intersection
//#####################################################################
template<class T> bool RENDERING_OBJECT<T>::
Intersection(RAY<TV>& ray,const int aggregate) const
{
    PHYSBAM_FUNCTION_IS_NOT_DEFINED();
}
//#####################################################################
// Function Preprocess_Efficiency_Structures
//#####################################################################
template<class T> void RENDERING_OBJECT<T>::
Preprocess_Efficiency_Structures(RENDER_WORLD<T>& world)
{
}
//#####################################################################
// Function Get_Aggregate_World_Space_Bounding_Boxes
//#####################################################################
template<class T> void RENDERING_OBJECT<T>::
Get_Aggregate_World_Space_Bounding_Boxes(ARRAY<RENDERING_OBJECT_ACCELERATION_PRIMITIVE<T> >& primitives) const
{
    PHYSBAM_FUNCTION_IS_NOT_DEFINED();
}
//#####################################################################
// Function Normal
//#####################################################################
template<class T> VECTOR<T,3> RENDERING_OBJECT<T>::
Normal(const TV& location,const int aggregate) const
{
    PHYSBAM_FUNCTION_IS_NOT_DEFINED();
}
//#####################################################################
// Function Inside
//#####################################################################
template<class T> bool RENDERING_OBJECT<T>::
Inside(const TV& location) const
{
    PHYSBAM_FUNCTION_IS_NOT_DEFINED();
}
//#####################################################################
// Function Lazy_Inside
//#####################################################################
template<class T> bool RENDERING_OBJECT<T>::
Lazy_Inside(const TV& location) const
{
    PHYSBAM_FUNCTION_IS_NOT_DEFINED();
}
//#####################################################################
// Function Outside
//#####################################################################
template<class T> bool RENDERING_OBJECT<T>::
Outside(const TV& location) const
{
    PHYSBAM_FUNCTION_IS_NOT_DEFINED();
}
//#####################################################################
// Function Lazy_Outside
//#####################################################################
template<class T> bool RENDERING_OBJECT<T>::
Lazy_Outside(const TV& location) const
{
    PHYSBAM_FUNCTION_IS_NOT_DEFINED();
}
//#####################################################################
// Function Boundary
//#####################################################################
template<class T> bool RENDERING_OBJECT<T>::
Boundary(const TV& location) const
{
    PHYSBAM_FUNCTION_IS_NOT_DEFINED();
}
//#####################################################################
// Function Surface
//#####################################################################
template<class T> VECTOR<T,3> RENDERING_OBJECT<T>::
Surface(const TV& location) const
{
    PHYSBAM_FUNCTION_IS_NOT_DEFINED();
}
//#####################################################################
// Function Signed_Distance
//#####################################################################
template<class T> T RENDERING_OBJECT<T>::
Signed_Distance(const TV& location) const
{
    PHYSBAM_FUNCTION_IS_NOT_DEFINED();
}
//#####################################################################
// Function Has_Bounding_Box
//#####################################################################
template<class T> bool RENDERING_OBJECT<T>::
Has_Bounding_Box() const
{
    return false;
}
//#####################################################################
// Function Object_Space_Bounding_Box
//#####################################################################
template<class T> RANGE<VECTOR<T,3> > RENDERING_OBJECT<T>::
Object_Space_Bounding_Box() const
{
    PHYSBAM_FUNCTION_IS_NOT_DEFINED();
}
//#####################################################################
// Function Close_To_Open_Surface
//#####################################################################
template<class T> bool RENDERING_OBJECT<T>::
Close_To_Open_Surface(const TV& location,const T threshold_distance) const
{
    PHYSBAM_FUNCTION_IS_NOT_DEFINED();
}
//#####################################################################
// Function Get_Intersection_Range
//#####################################################################
template<class T> bool RENDERING_OBJECT<T>::
Get_Intersection_Range(const RAY<TV>& ray,T& start_t,T& end_t) const
{
    return false;
}
//#####################################################################
// Function Generate_BSSRDF_Tree
//#####################################################################
template<class T> void RENDERING_OBJECT<T>::
Generate_BSSRDF_Tree(RENDER_WORLD<T>& world)
{
}
//#####################################################################
// Function Generate_Triangles
//#####################################################################
template<class T> TRIANGULATED_SURFACE<T>* RENDERING_OBJECT<T>::
Generate_Triangles() const
{
    PHYSBAM_FUNCTION_IS_NOT_DEFINED();
    return 0;
};
namespace PhysBAM{
template class RENDERING_OBJECT<double>;
template class RENDERING_OBJECT<float>;
}
