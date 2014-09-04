//#####################################################################
// Copyright 2003-2007, Zhaosheng Bao, Eran Guendelman, Geoffrey Irving, Sergey Koltakov, Frank Losasso, Avi Robinson-Mosher, Craig Schroeder, Tamar Shinar, Rachel Weinstein.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Geometry/Basic_Geometry/SEGMENT_2D.h>
#include <Geometry/Implicit_Objects/IMPLICIT_OBJECT_TRANSFORMED.h>
#include <Geometry/Implicit_Objects_Uniform/LEVELSET_IMPLICIT_OBJECT.h>
#include <Geometry/Topology_Based_Geometry/SEGMENTED_CURVE_2D.h>
#include <Geometry/Topology_Based_Geometry/TRIANGULATED_AREA.h>
#include <Rigids/Collisions/RIGID_COLLISION_GEOMETRY_2D.h>
#include <Rigids/Rigid_Bodies/RIGID_BODY.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class T> RIGID_COLLISION_GEOMETRY<VECTOR<T,2> >::
RIGID_COLLISION_GEOMETRY(RIGID_BODY<TV>& rigid_body_input)
    :RIGID_COLLISION_GEOMETRY_BASE<TV>(rigid_body_input)
{}
//#####################################################################
// Destructor
//#####################################################################
template<class T> RIGID_COLLISION_GEOMETRY<VECTOR<T,2> >::
~RIGID_COLLISION_GEOMETRY()
{}
//##################################################################### 
// Function Earliest_Simplex_Crossover
//##################################################################### 
template<class T> bool RIGID_COLLISION_GEOMETRY<VECTOR<T,2> >::
Earliest_Simplex_Crossover(const VECTOR<T,2>& start_X,const VECTOR<T,2>& end_X,const T dt,T& hit_time,TV& weights,int& segment_id) const
{
    const T collision_thickness_over_two=(T).5*collision_thickness;
    T min_time=FLT_MAX;bool collision=false;T current_hit_time;VECTOR<T,TV::m+1> current_weights;
    TV normal;
    for(int segment_number=0;segment_number<rigid_body.simplicial_object->mesh.elements.m;segment_number++){
        SEGMENT_2D<T> initial_segment=World_Space_Simplex(segment_number),final_segment=World_Space_Simplex(segment_number,saved_states(0).x);
        POINT_SIMPLEX_COLLISION_TYPE collision_type=
            SEGMENT_2D<T>::Robust_Point_Face_Collision(initial_segment,final_segment,start_X,end_X,dt,collision_thickness_over_two,current_hit_time,normal,current_weights);
        if(collision_type!=POINT_SIMPLEX_NO_COLLISION && current_hit_time<min_time){
            min_time=hit_time=current_hit_time;weights=current_weights.Remove_Index(0);segment_id=segment_number;collision=true;}}
    return collision;
}
//##################################################################### 
// Function Latest_Simplex_Crossover
//##################################################################### 
template<class T> bool RIGID_COLLISION_GEOMETRY<VECTOR<T,2> >::
Latest_Simplex_Crossover(const VECTOR<T,2>& start_X,const VECTOR<T,2>& end_X,const T dt,T& hit_time,TV& weights,int& segment_id,POINT_SIMPLEX_COLLISION_TYPE& returned_collision_type) const
{
    const T collision_thickness_over_two=(T).5*collision_thickness;
    returned_collision_type=POINT_SIMPLEX_NO_COLLISION;
    T max_time=-FLT_MAX;bool collision=false;T current_hit_time;VECTOR<T,TV::m+1> current_weights;
    TV normal;
    for(int segment_number=0;segment_number<rigid_body.simplicial_object->mesh.elements.m;segment_number++){
        SEGMENT_2D<T> initial_segment=World_Space_Simplex(segment_number),final_segment=World_Space_Simplex(segment_number,saved_states(0).x);
        POINT_SIMPLEX_COLLISION_TYPE collision_type=
            SEGMENT_2D<T>::Robust_Point_Face_Collision(initial_segment,final_segment,start_X,end_X,dt,collision_thickness_over_two,current_hit_time,normal,current_weights);
        if(collision_type!=POINT_SIMPLEX_NO_COLLISION && current_hit_time>max_time){
            max_time=hit_time=current_hit_time;weights=current_weights.Remove_Index(0);segment_id=segment_number;collision=true;returned_collision_type=collision_type;}}
    return collision;
}
//#####################################################################
// Function Any_Simplex_Crossover
//#####################################################################
template<class T> bool RIGID_COLLISION_GEOMETRY<VECTOR<T,2> >::
Any_Simplex_Crossover(const VECTOR<T,2>& start_X,const VECTOR<T,2>& end_X,const T dt) const
{
    const T collision_thickness_over_two=(T).5*collision_thickness;
    T current_hit_time;VECTOR<T,TV::m+1> current_weights;
    TV normal;
    for(int segment_number=0;segment_number<rigid_body.simplicial_object->mesh.elements.m;segment_number++){
        SEGMENT_2D<T> initial_segment=World_Space_Simplex(segment_number),final_segment=World_Space_Simplex(segment_number,saved_states(0).x);
        POINT_SIMPLEX_COLLISION_TYPE collision_type=
            SEGMENT_2D<T>::Robust_Point_Face_Collision(initial_segment,final_segment,start_X,end_X,dt,collision_thickness_over_two,current_hit_time,normal,current_weights);
        if(collision_type!=POINT_SIMPLEX_NO_COLLISION) return true;}
    return false;
}
//#####################################################################
// Function Get_Simplex_Bounding_Boxes
//#####################################################################
template<class T> void RIGID_COLLISION_GEOMETRY<VECTOR<T,2> >::
Get_Simplex_Bounding_Boxes(ARRAY<RANGE<TV> >& bounding_boxes,const bool with_body_motion,const T extra_thickness,const T body_thickness_factor) const
{
    if(!rigid_body.simplicial_object->segment_list) rigid_body.simplicial_object->Update_Segment_List();
    for(int t=0;t<rigid_body.simplicial_object->mesh.elements.m;t++){
        RANGE<TV> box=rigid_body.World_Space_Simplex_Bounding_Box(t);
        if(with_body_motion) box.Enlarge_To_Include_Box(rigid_body.World_Space_Simplex_Bounding_Box(t,saved_states(0).x));
        box.Change_Size(extra_thickness+body_thickness_factor*collision_thickness);
        bounding_boxes.Append(box);}
}
//#####################################################################
// Function Update_Intersection_Acceleration_Structures
//#####################################################################
template<class T> void RIGID_COLLISION_GEOMETRY<VECTOR<T,2> >::
Update_Intersection_Acceleration_Structures(const bool use_swept_triangle_hierarchy,const int state1,const int state2)
{}
//##################################################################### 
// Function World_Space_Simplex
//##################################################################### 
template<class T> SEGMENT_2D<T> RIGID_COLLISION_GEOMETRY<VECTOR<T,2> >::
World_Space_Simplex(const int segment_id,const bool use_saved_state) const
{
    if(use_saved_state) return World_Space_Simplex(segment_id,saved_states(0).x);
    return rigid_body.World_Space_Simplex(segment_id);
}
//##################################################################### 
// Function World_Space_Simplex
//##################################################################### 
template<class T> SEGMENT_2D<T> RIGID_COLLISION_GEOMETRY<VECTOR<T,2> >::
World_Space_Simplex(const int segment_id,const FRAME<TV>& frame) const
{
    SEGMENT_2D<T>& object_space_segment=(*rigid_body.simplicial_object->segment_list)(segment_id);
    return SEGMENT_2D<T>(frame*object_space_segment.X.x,frame*object_space_segment.X.y);
}
//##################################################################### 
namespace PhysBAM{
template class RIGID_COLLISION_GEOMETRY<VECTOR<float,2> >;
template class RIGID_COLLISION_GEOMETRY<VECTOR<double,2> >;
}
