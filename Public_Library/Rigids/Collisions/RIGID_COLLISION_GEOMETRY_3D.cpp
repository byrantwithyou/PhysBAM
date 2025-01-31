//#####################################################################
// Copyright 2003-2008, Zhaosheng Bao, Ronald Fedkiw, Eran Guendelman, Sergey Koltakov, Frank Losasso, Avi Robinson-Mosher, Craig Schroeder, Andrew Selle, Tamar Shinar, Eftychios Sifakis, Rachel Weinstein.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Geometry/Implicit_Objects/IMPLICIT_OBJECT_TRANSFORMED.h>
#include <Geometry/Spatial_Acceleration/TRIANGLE_HIERARCHY.h>
#include <Geometry/Topology_Based_Geometry/TRIANGULATED_SURFACE.h>
#include <Rigids/Collisions/RIGID_COLLISION_GEOMETRY_3D.h>
#include <Rigids/Rigid_Bodies/RIGID_BODY.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class T> RIGID_COLLISION_GEOMETRY<VECTOR<T,3> >::
RIGID_COLLISION_GEOMETRY(RIGID_BODY<TV>& rigid_body_input)
    :RIGID_COLLISION_GEOMETRY_BASE<TV>(rigid_body_input)
{}
//#####################################################################
// Destructor
//#####################################################################
template<class T> RIGID_COLLISION_GEOMETRY<VECTOR<T,3> >::
~RIGID_COLLISION_GEOMETRY()
{}
//##################################################################### 
// Function Earliest_Simplex_Crossover
//##################################################################### 
template<class T> bool RIGID_COLLISION_GEOMETRY<VECTOR<T,3> >::
Earliest_Simplex_Crossover(const TV& start_X,const TV& end_X,const T dt,T& hit_time,TV& weights,int& triangle_id) const
{
    const T collision_thickness_over_two=(T).5*collision_thickness;
    T min_time=FLT_MAX;bool collision=false;T current_hit_time;TV current_normal;VECTOR<T,TV::m+1> current_weights;
    if(rigid_body.moving_simplex_hierarchy){
        ARRAY<int> triangles_to_check;
        if(start_X==end_X){rigid_body.moving_simplex_hierarchy->Intersection_List(start_X,triangles_to_check,collision_thickness);}
        else rigid_body.moving_simplex_hierarchy->Intersection_List(RANGE<TV>::Bounding_Box(start_X,end_X),triangles_to_check,collision_thickness);
        for(int i=0;i<triangles_to_check.m;i++){
            int t=triangles_to_check(i);
            TRIANGLE_3D<T> initial_triangle=World_Space_Simplex(t),final_triangle=World_Space_Simplex(t,saved_states(0).x);
            POINT_SIMPLEX_COLLISION_TYPE collision_type=TRIANGLE_3D<T>::Robust_Point_Face_Collision(initial_triangle,final_triangle,start_X,end_X,dt,collision_thickness_over_two,current_hit_time,
                current_normal,current_weights);
            if(collision_type!=POINT_SIMPLEX_NO_COLLISION && current_hit_time < min_time){
                min_time=hit_time=current_hit_time;weights=current_weights.Remove_Index(0);triangle_id=t;collision=true;}}}
    else{
        LOG::cerr<<"Earliest_Simplex_Crossover: no moving triangle hierarchy"<<std::endl;
        for(int t=0;t<rigid_body.simplicial_object->mesh.elements.m;t++){   
            TRIANGLE_3D<T> initial_triangle=World_Space_Simplex(t),final_triangle=World_Space_Simplex(t,saved_states(0).x);
            POINT_SIMPLEX_COLLISION_TYPE collision_type=TRIANGLE_3D<T>::Robust_Point_Face_Collision(initial_triangle,final_triangle,start_X,end_X,dt,collision_thickness_over_two,current_hit_time,
                current_normal,current_weights);
            if(collision_type!=POINT_SIMPLEX_NO_COLLISION && current_hit_time < min_time){
                min_time=hit_time=current_hit_time;weights=current_weights.Remove_Index(0);triangle_id=t;collision=true;}}}
    return collision;
}
//##################################################################### 
// Function Latest_Simplex_Crossover
//##################################################################### 
template<class T> bool RIGID_COLLISION_GEOMETRY<VECTOR<T,3> >::
Latest_Simplex_Crossover(const TV& start_X,const TV& end_X,const T dt,T& hit_time,TV& weights,int& triangle_id,
    POINT_SIMPLEX_COLLISION_TYPE& returned_collision_type) const
{
    const T collision_thickness_over_two=(T).5*collision_thickness;
    returned_collision_type=POINT_SIMPLEX_NO_COLLISION;
    T max_time=-FLT_MAX;bool collision=false;T current_hit_time;TV current_normal;VECTOR<T,TV::m+1> current_weights;
    if(rigid_body.moving_simplex_hierarchy){
        ARRAY<int> triangles_to_check;
        if(start_X==end_X){rigid_body.moving_simplex_hierarchy->Intersection_List(start_X,triangles_to_check,collision_thickness);}
        else rigid_body.moving_simplex_hierarchy->Intersection_List(RANGE<TV>::Bounding_Box(start_X,end_X),triangles_to_check,collision_thickness);
        for(int i=0;i<triangles_to_check.m;i++){
            int t=triangles_to_check(i);
            TRIANGLE_3D<T> initial_triangle=World_Space_Simplex(t),final_triangle=World_Space_Simplex(t,saved_states(0).x);
            POINT_SIMPLEX_COLLISION_TYPE collision_type=TRIANGLE_3D<T>::Robust_Point_Face_Collision(initial_triangle,final_triangle,start_X,end_X,dt,collision_thickness_over_two,current_hit_time,
                current_normal,current_weights);
            if(collision_type!=POINT_SIMPLEX_NO_COLLISION && current_hit_time > max_time){
                max_time=hit_time=current_hit_time;weights=current_weights.Remove_Index(0);triangle_id=t;collision=true;returned_collision_type=collision_type;}}}
    else{
        LOG::cerr<<"Latest_Simplex_Crossover: no moving triangle hierarchy"<<std::endl;
        for(int t=0;t<rigid_body.simplicial_object->mesh.elements.m;t++){   
            TRIANGLE_3D<T> initial_triangle=World_Space_Simplex(t),final_triangle=World_Space_Simplex(t,saved_states(0).x);
            POINT_SIMPLEX_COLLISION_TYPE collision_type=TRIANGLE_3D<T>::Robust_Point_Face_Collision(initial_triangle,final_triangle,start_X,end_X,dt,collision_thickness_over_two,current_hit_time,
                current_normal,current_weights);
            if(collision_type!=POINT_SIMPLEX_NO_COLLISION && current_hit_time > max_time){
                max_time=hit_time=current_hit_time;weights=current_weights.Remove_Index(0);triangle_id=t;collision=true;returned_collision_type=collision_type;}}}
    return collision;
}
//##################################################################### 
// Function Any_Simplex_Crossover
//##################################################################### 
template<class T> bool RIGID_COLLISION_GEOMETRY<VECTOR<T,3> >::
Any_Simplex_Crossover(const TV& start_X,const TV& end_X,const T dt) const
{
    const T collision_thickness_over_two=(T).5*collision_thickness;
    T hit_time;TV normal;VECTOR<T,TV::m+1> weights;
    if(rigid_body.moving_simplex_hierarchy){
        ARRAY<int> triangles_to_check;
        if(start_X==end_X){rigid_body.moving_simplex_hierarchy->Intersection_List(start_X,triangles_to_check,collision_thickness);}
        else rigid_body.moving_simplex_hierarchy->Intersection_List(RANGE<TV>::Bounding_Box(start_X,end_X),triangles_to_check,collision_thickness);
        for(int i=0;i<triangles_to_check.m;i++){
            int t=triangles_to_check(i);
            TRIANGLE_3D<T> initial_triangle=World_Space_Simplex(t),final_triangle=World_Space_Simplex(t,saved_states(0).x);
            POINT_SIMPLEX_COLLISION_TYPE collision_type=TRIANGLE_3D<T>::Robust_Point_Face_Collision(initial_triangle,final_triangle,start_X,end_X,dt,collision_thickness_over_two,hit_time,normal,
                weights);
            if(collision_type!=POINT_SIMPLEX_NO_COLLISION) return true;}}
    else{
        LOG::cerr<<"Any_Simplex_Crossover: no moving triangle hierarchy"<<std::endl;
        for(int t=0;t<rigid_body.simplicial_object->mesh.elements.m;t++){   
            TRIANGLE_3D<T> initial_triangle=World_Space_Simplex(t),final_triangle=World_Space_Simplex(t,saved_states(0).x);
            POINT_SIMPLEX_COLLISION_TYPE collision_type=TRIANGLE_3D<T>::Robust_Point_Face_Collision(initial_triangle,final_triangle,start_X,end_X,dt,collision_thickness_over_two,hit_time,normal,
                weights);
            if(collision_type!=POINT_SIMPLEX_NO_COLLISION) return true;}}
    return false;
}
//#####################################################################
// Function Get_Simplex_Bounding_Boxes
//#####################################################################
template<class T> void RIGID_COLLISION_GEOMETRY<VECTOR<T,3> >::
Get_Simplex_Bounding_Boxes(ARRAY<RANGE<TV> >& bounding_boxes,const bool with_body_motion,const T extra_thickness,const T body_thickness_factor) const
{
    if(!rigid_body.simplicial_object->triangle_list) rigid_body.simplicial_object->Update_Triangle_List();
    for(int t=0;t<rigid_body.simplicial_object->mesh.elements.m;t++){
        RANGE<TV> box=rigid_body.World_Space_Simplex_Bounding_Box(t);
        if(with_body_motion) box.Enlarge_To_Include_Box(rigid_body.World_Space_Simplex_Bounding_Box(t,saved_states(0).x));
        box.Change_Size(extra_thickness+body_thickness_factor*collision_thickness);
        bounding_boxes.Append(box);}
}
//##################################################################### 
// Function Update_Intersection_Acceleration_Structures
//##################################################################### 
template<class T> void RIGID_COLLISION_GEOMETRY<VECTOR<T,3> >::
Update_Intersection_Acceleration_Structures(const bool use_swept_triangle_hierarchy,const int state1,const int state2)
{
    if(!use_swept_triangle_hierarchy) return;
    if(!rigid_body.moving_simplex_hierarchy) rigid_body.moving_simplex_hierarchy=new TRIANGLE_HIERARCHY<T>(rigid_body.simplicial_object->mesh,rigid_body.simplicial_object->particles,false);
    rigid_body.moving_simplex_hierarchy->Update_Boxes(saved_states(state1).x,saved_states(state2).x);
}
//#####################################################################
// Function World_Space_Simplex
//#####################################################################
template<class T> TRIANGLE_3D<T> RIGID_COLLISION_GEOMETRY<VECTOR<T,3> >::
World_Space_Simplex(const int triangle_id,const bool use_saved_state) const
{
    if(use_saved_state) return World_Space_Simplex(triangle_id,saved_states(0).x);
    return rigid_body.World_Space_Simplex(triangle_id);
}
//#####################################################################
// Function World_Space_Simplex
//#####################################################################
template<class T> TRIANGLE_3D<T> RIGID_COLLISION_GEOMETRY<VECTOR<T,3> >::
World_Space_Simplex(const int triangle_id,const FRAME<TV>& frame) const
{
    return TRIANGLE_3D<T>(frame*rigid_body.simplicial_object->particles.X(rigid_body.simplicial_object->mesh.elements(triangle_id)(0)),
        frame*rigid_body.simplicial_object->particles.X(rigid_body.simplicial_object->mesh.elements(triangle_id)(1)),
        frame*rigid_body.simplicial_object->particles.X(rigid_body.simplicial_object->mesh.elements(triangle_id)(2)));
}
//#####################################################################
namespace PhysBAM{
template class RIGID_COLLISION_GEOMETRY<VECTOR<float,3> >;
template class RIGID_COLLISION_GEOMETRY<VECTOR<double,3> >;
}
