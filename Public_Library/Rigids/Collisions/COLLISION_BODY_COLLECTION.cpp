//#####################################################################
// Copyright 2004-2006, Ron Fedkiw, Geoffrey Irving, Frank Losasso, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Core/Arrays_Nd/ARRAYS_ND.h>
#include <Core/Read_Write/FILE_UTILITIES.h>
#include <Geometry/Basic_Geometry/RAY.h>
#include <Geometry/Basic_Geometry/SEGMENT_2D.h>
#include <Geometry/Basic_Geometry/TRIANGLE_3D.h>
#include <Rigids/Collisions/COLLISION_BODY_COLLECTION.h>
#include <Rigids/Collisions/COLLISION_GEOMETRY.h>
#include <Rigids/Collisions/COLLISION_GEOMETRY_SPATIAL_PARTITION.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class TV> COLLISION_BODY_COLLECTION<TV>::
COLLISION_BODY_COLLECTION()
    :spatial_partition(0),collision_body_thickness(0)
{}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> COLLISION_BODY_COLLECTION<TV>::
~COLLISION_BODY_COLLECTION()
{
    delete spatial_partition;
    for(COLLISION_GEOMETRY_ID id(0);id<bodies.m;id++) if(owns_collision_geometry(id)) delete bodies(id);
}
//##################################################################### 
// Function Add_Body
//##################################################################### 
template<class TV> COLLISION_GEOMETRY_ID COLLISION_BODY_COLLECTION<TV>::
Add_Body(COLLISION_GEOMETRY<TV>* body,const int geometry_id,bool owns_collision_geometry_input)
{
    COLLISION_GEOMETRY_ID id;
    if(deletion_list.m){
        id=deletion_list.Pop_Value();bodies(id)=body;
        owns_collision_geometry(id)=owns_collision_geometry_input;}
    else{
        id=bodies.Append(body);
        owns_collision_geometry.Append(owns_collision_geometry_input);}
    body->collision_geometry_id=id;
    owns_collision_geometry(id)=owns_collision_geometry_input;
    geometry_id_to_collision_geometry_id.Set(geometry_id,id);
    collision_geometry_id_to_geometry_id.Set(id,geometry_id);
    return id;
}
//##################################################################### 
// Function Add_Bodies
//##################################################################### 
template<class TV> void COLLISION_BODY_COLLECTION<TV>::
Add_Bodies(COLLISION_BODY_COLLECTION<TV>& collision_body_list)
{
    for(COLLISION_GEOMETRY_ID i(0);i<collision_body_list.bodies.m;i++)
        if(collision_body_list.Is_Active(i))
            Add_Body(collision_body_list.bodies(i),collision_body_list.collision_geometry_id_to_geometry_id.Get(i),false);
}
//#####################################################################
// Function Remove_Body
//#####################################################################
template<class TV> void COLLISION_BODY_COLLECTION<TV>::
Remove_Body(COLLISION_GEOMETRY_ID id)
{
    PHYSBAM_ASSERT(Is_Active(id));
    bodies(id)->collision_geometry_id=COLLISION_GEOMETRY_ID(0);
    if(owns_collision_geometry(id)){delete bodies(id);owns_collision_geometry(id)=false;}
    geometry_id_to_collision_geometry_id.Delete_If_Present(collision_geometry_id_to_geometry_id.Get(id));
    collision_geometry_id_to_geometry_id.Delete(id);
    bodies(id)=0;
    deletion_list.Append(id);
}
//##################################################################### 
// Function Intersection_Between_Points
//##################################################################### 
template<class TV> bool COLLISION_BODY_COLLECTION<TV>::
Intersection_Between_Points(const TV& x0,const TV& x1,COLLISION_GEOMETRY_ID& body_id,int& triangle_id,TV& intersection_point,const ARRAY<COLLISION_GEOMETRY_ID>* objects) const
{
    TV ray_vector=x1-x0;T ray_vector_length_squared=ray_vector.Magnitude_Squared();
    if(ray_vector_length_squared==0){
        if(Inside_Any_Simplex_Of_Any_Body(x0,body_id,triangle_id,objects)){intersection_point=x0;return true;}
        else return false;}
    else{
        T ray_t_max=sqrt(ray_vector_length_squared);RAY<TV> ray(x0,ray_vector/ray_t_max,true);ray.semi_infinite=false;ray.t_max=ray_t_max;
        if(Intersection_With_Any_Simplicial_Object(ray,body_id,objects)){triangle_id=ray.aggregate_id;intersection_point=ray.Point(ray.t_max);return true;}
        else return false;}
}
//##################################################################### 
// Function Intersection_Between_Points
//##################################################################### 
template<class TV> bool COLLISION_BODY_COLLECTION<TV>::
Intersection_Between_Points(const TV& x0,const TV& x1,const ARRAY<COLLISION_GEOMETRY_ID>* objects) const
{
    TV ray_vector=x1-x0;T ray_vector_length_squared=ray_vector.Magnitude_Squared();COLLISION_GEOMETRY_ID body_id;int triangle_id;
    if(ray_vector_length_squared==0)return Inside_Any_Simplex_Of_Any_Body(x0,body_id,triangle_id,objects);
    else{
        T ray_t_max=sqrt(ray_vector_length_squared);RAY<TV> ray(x0,ray_vector/ray_t_max,true);ray.semi_infinite=false;ray.t_max=ray_t_max;
        return Intersection_With_Any_Simplicial_Object(ray,body_id,objects);}
}
//##################################################################### 
// Function Intersection_With_Any_Simplicial_Object
//##################################################################### 
template<class TV> bool COLLISION_BODY_COLLECTION<TV>::
Intersection_With_Any_Simplicial_Object(RAY<TV>& ray,COLLISION_GEOMETRY_ID& body_id,const ARRAY<COLLISION_GEOMETRY_ID>* objects) const 
{
    bool hit=false;
    if(!objects){for(COLLISION_GEOMETRY_ID i(0);i<bodies.m;i++) if(Is_Active(i) && bodies(i)->active && bodies(i)->Simplex_Intersection(ray)){body_id=i;hit=true;}}
    else for(int k=0;k<objects->m;k++){COLLISION_GEOMETRY_ID i=(*objects)(k);if(Is_Active(i) && bodies(i)->active && bodies(i)->Simplex_Intersection(ray)){body_id=i;hit=true;}}
    return hit;
}
//##################################################################### 
// Function Closest_Non_Intersecting_Point_Of_Any_Simplicial_Object
//##################################################################### 
template<class TV> bool COLLISION_BODY_COLLECTION<TV>::
Closest_Non_Intersecting_Point_Of_Any_Simplicial_Object(RAY<TV>& ray,COLLISION_GEOMETRY_ID& body_id,const ARRAY<COLLISION_GEOMETRY_ID>* objects) const 
{
    bool hit=false;
    if(!objects){for(COLLISION_GEOMETRY_ID i(0);i<bodies.m;i++)if(Is_Active(i) && bodies(i)->active && bodies(i)->Simplex_Closest_Non_Intersecting_Point(ray)){body_id=i;hit=true;}}
    else for(int k=0;k<objects->m;k++){COLLISION_GEOMETRY_ID i=(*objects)(k);if(Is_Active(i) && bodies(i)->active && bodies(i)->Simplex_Closest_Non_Intersecting_Point(ray)){body_id=i;hit=true;}}
    return hit;
}
//##################################################################### 
// Function Inside_Any_Simplex_Of_Any_Body
//##################################################################### 
template<class TV> bool COLLISION_BODY_COLLECTION<TV>::
Inside_Any_Simplex_Of_Any_Body(const TV& location,COLLISION_GEOMETRY_ID& body_id,int& simplex_id,const ARRAY<COLLISION_GEOMETRY_ID>* objects) const
{
    if(!objects){
        for(COLLISION_GEOMETRY_ID i(0);i<bodies.m;i++)
            if(Is_Active(i) && bodies(i)->active && bodies(i)->Inside_Any_Simplex(location,simplex_id)){
                body_id=i;
                return true;}}
    else
        for(int k=0;k<objects->m;k++){
            COLLISION_GEOMETRY_ID i=(*objects)(k);
            if(Is_Active(i) && bodies(i)->active && bodies(i)->Inside_Any_Simplex(location,simplex_id)){
                body_id=i;
                return true;}}
    return false;
}
//##################################################################### 
// Function Inside_Any_Body
//##################################################################### 
template<class TV> bool COLLISION_BODY_COLLECTION<TV>::
Inside_Any_Body(const TV& location,const T thickness_over_two,COLLISION_GEOMETRY_ID& body_id,const ARRAY<COLLISION_GEOMETRY_ID>* objects) const
{
    if(!objects){
        for(COLLISION_GEOMETRY_ID i(0);i<bodies.m;i++)
            if(Is_Active(i) && bodies(i)->active && bodies(i)->Inside(location,thickness_over_two)){
                body_id=i;
                return true;}}
    else
        for(int k=0;k<objects->m;k++){
            COLLISION_GEOMETRY_ID i=(*objects)(k);
            if(Is_Active(i) && bodies(i)->active && bodies(i)->Inside(location,thickness_over_two)){
                body_id=i;
                return true;}}
    return false;
}
//##################################################################### 
// Function Implicit_Geometry_Lazy_Inside_Any_Body
//##################################################################### 
template<class TV> bool COLLISION_BODY_COLLECTION<TV>::
Implicit_Geometry_Lazy_Inside_Any_Body(const TV& location,COLLISION_GEOMETRY_ID& body_id,const ARRAY<COLLISION_GEOMETRY_ID>* objects) const
{
    if(!objects){
        for(COLLISION_GEOMETRY_ID i(0);i<bodies.m;i++){
            if(Is_Active(i) && bodies(i)->active && bodies(i)->Implicit_Geometry_Lazy_Inside(location)){body_id=i;return true;}}}
    else{
        for(int k=0;k<objects->m;k++){COLLISION_GEOMETRY_ID i=(*objects)(k);
            if(Is_Active(i) && bodies(i)->active && bodies(i)->Implicit_Geometry_Lazy_Inside(location)){body_id=i;return true;}}}
    return false;
}
//##################################################################### 
// Function Closest_Boundary_Point
//##################################################################### 
template<class TV> TV COLLISION_BODY_COLLECTION<TV>::
Closest_Boundary_Point(const TV& location,const T max_distance,T& distance,COLLISION_GEOMETRY_ID& body_id,int& simplex_id,const ARRAY<COLLISION_GEOMETRY_ID>* objects) const
{
    distance=FLT_MAX;TV point,current_point;T current_distance;int current_simplex_id;
    if(!objects){
        for(COLLISION_GEOMETRY_ID i(0);i<bodies.m;i++)
            if(Is_Active(i) && bodies(i)->active){
// TODO The results from this give really bad garbage...?
//         current_point=bodies(i)->Closest_Point_On_Boundary(location,max_distance,(T)0,1,&current_simplex_id,&current_distance);
                current_point=bodies(i)->Simplex_Closest_Point_On_Boundary(location,max_distance,(T)0,&current_simplex_id,&current_distance);
                if(current_distance<distance){distance=current_distance;
                    point=current_point;
                    simplex_id=current_simplex_id;
                    body_id=i;}}}
    else
        for(int k=0;k<objects->m;k++){
            COLLISION_GEOMETRY_ID i=(*objects)(k);
            if(Is_Active(i) && bodies(i)->active){
//         current_point=bodies(i)->Closest_Point_On_Boundary(location,max_distance,(T)0,1,&current_simplex_id,&current_distance);
                current_point=bodies(i)->Simplex_Closest_Point_On_Boundary(location,max_distance,(T)0,&current_simplex_id,&current_distance);
                if(current_distance<distance){
                    distance=current_distance;
                    point=current_point;
                    simplex_id=current_simplex_id;
                    body_id=i;}}}
    return point;
}
//#####################################################################
// Function Update_Spatial_Partition
//#####################################################################
template<class TV> void COLLISION_BODY_COLLECTION<TV>::
Update_Spatial_Partition(const SPATIAL_PARTITION_VOXEL_SIZE_HEURISTIC heuristic,const int number_of_boxes,const T voxel_size_scale_factor)
{
    if(!spatial_partition){
        spatial_partition=new COLLISION_GEOMETRY_SPATIAL_PARTITION<TV,const ARRAY<COLLISION_GEOMETRY<TV>*,COLLISION_GEOMETRY_ID>,COLLISION_GEOMETRY_ID>(bodies,collision_body_thickness);
        spatial_partition->Compute_Voxel_Size(heuristic,number_of_boxes,voxel_size_scale_factor);}
    spatial_partition->Reinitialize();
}
//##################################################################### 
// Function Update_Bounding_Boxes
//##################################################################### 
template<class TV> void COLLISION_BODY_COLLECTION<TV>::
Update_Bounding_Boxes()
{
    for(COLLISION_GEOMETRY_ID i(0);i<bodies.m;i++) if(Is_Active(i)) bodies(i)->Update_Bounding_Box();
}
//#####################################################################
namespace PhysBAM{
template class COLLISION_BODY_COLLECTION<VECTOR<float,3> >;
template class COLLISION_BODY_COLLECTION<VECTOR<float,2> >;
template class COLLISION_BODY_COLLECTION<VECTOR<float,1> >;
template class COLLISION_BODY_COLLECTION<VECTOR<double,3> >;
template class COLLISION_BODY_COLLECTION<VECTOR<double,2> >;
template class COLLISION_BODY_COLLECTION<VECTOR<double,1> >;
}
