//#####################################################################
// Copyright 2009, Jon Gretarsson, Nipun Kwatra.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Geometry/Registry/STRUCTURE_REGISTRY.h>
#include <Geometry/Topology_Based_Geometry/POINT_SIMPLICES_1D.h>
namespace PhysBAM{
//#####################################################################
// Constructor
//#####################################################################
template<class T> POINT_SIMPLICES_1D<T>::
POINT_SIMPLICES_1D()
    :POINT_SIMPLICES_1D(*new POINT_SIMPLEX_MESH,*new GEOMETRY_PARTICLES<TV>)
{
    this->need_destroy_mesh=this->need_destroy_particles=true;
}
//#####################################################################
// Constructor
//#####################################################################
template<class T> POINT_SIMPLICES_1D<T>::
POINT_SIMPLICES_1D(POINT_SIMPLEX_MESH& point_simplex_mesh_input,GEOMETRY_PARTICLES<TV>& particles_input)
    :MESH_OBJECT<TV,POINT_SIMPLEX_MESH>(point_simplex_mesh_input,particles_input),point_simplex_list(0),particle_partition(0),hierarchy(0),number_point_simplices(0)
{
}
//#####################################################################
// Function Closest_Point_On_Boundary
//#####################################################################
template<class T> VECTOR<T,1> POINT_SIMPLICES_1D<T>::
Closest_Point_On_Boundary(const TV& location,const T max_depth,const T thickness_over_two,int* closest_point_simplex,T* distance) const 
{
    T min_distance_squared=FLT_MAX;TV closest_point;
    for(int i=0;i<mesh.elements.m;i++){
        POINT_SIMPLEX_1D<T> point_simplex=point_simplex_list?(*point_simplex_list)(i):Get_Element(i);
        TV new_point=point_simplex.X.x;
        T distance_squared=(new_point-location).Magnitude_Squared();
        if(distance_squared < min_distance_squared){
            min_distance_squared=distance_squared;closest_point=new_point;
            if(closest_point_simplex) *closest_point_simplex=i;}}
    if(distance) *distance=sqrt(min_distance_squared);
    return closest_point;
}
//#####################################################################
// Function Calculate_Signed_Distance
//#####################################################################
template<class T> T POINT_SIMPLICES_1D<T>:: 
Calculate_Signed_Distance(const TV& location,T thickness_over_two) const
{
    T distance;
    Closest_Point_On_Boundary(location,0,thickness_over_two,0,&distance);
    if(Inside(location,thickness_over_two)) distance*=-1;
    return distance;
}
//#####################################################################
// Function Inside
//#####################################################################
template<class T> bool POINT_SIMPLICES_1D<T>::
Inside(const TV& location,T thickness_over_two) const
{
    for(int i=0;i<mesh.elements.m;i++){
        POINT_SIMPLEX_1D<T> point_simplex=point_simplex_list?(*point_simplex_list)(i):Get_Element(i);
        T direction=(T)(mesh.directions(i)?1:-1);
        T robust_point_simplex_location=point_simplex.X.x.x+direction*thickness_over_two;
        if(direction*(location.x-robust_point_simplex_location)>0) return false;}
    return true;
}
//#####################################################################
// Function Outside
//#####################################################################
template<class T> bool POINT_SIMPLICES_1D<T>::
Outside(const TV& location,T thickness_over_two) const
{
    return !Inside(location,thickness_over_two);
}
//#####################################################################
// Function Inside_Any_Simplex
//#####################################################################
template<class T> bool POINT_SIMPLICES_1D<T>::
Inside_Any_Simplex(const TV& point,int& point_simplex_id,const T thickness_over_two) const
{
    for(int i=0;i<mesh.elements.m;i++){
        POINT_SIMPLEX_1D<T> point_simplex=point_simplex_list?(*point_simplex_list)(i):Get_Element(i);
        if(point_simplex.X.x.x-thickness_over_two <= point.x && point_simplex.X.x.x+thickness_over_two >= point.x) return true;}
    return false;
}
//#####################################################################
template class POINT_SIMPLICES_1D<float>;
template class POINT_SIMPLICES_1D<double>;
}
