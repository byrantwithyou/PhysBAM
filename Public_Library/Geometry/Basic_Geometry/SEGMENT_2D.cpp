//#####################################################################
// Copyright 2003-2007, Zhaosheng Bao, Ronald Fedkiw, Eran Guendelman, Nipun Kwatra, Avi Robinson-Mosher, Andrew Selle, Jonathan Su.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class SEGMENT_2D  
//##################################################################### 
#include <Core/Arrays/ARRAY.h>
#include <Core/Arrays/INDIRECT_ARRAY.h>
#include <Core/Log/LOG.h>
#include <Tools/Interpolation/LINEAR_INTERPOLATION.h>
#include <Tools/Polynomials/CUBIC.h>
#include <Tools/Polynomials/QUADRATIC.h>
#include <Geometry/Basic_Geometry/LINE_2D.h>
#include <Geometry/Basic_Geometry/ORIENTED_BOX.h>
#include <Geometry/Basic_Geometry/SEGMENT_2D.h>
#include <Geometry/Level_Sets/LEVELSET_UTILITIES.h>
using namespace PhysBAM;
//#####################################################################
// Segment_Line_Intersection
//#####################################################################
template<class T> bool SEGMENT_2D<T>::
Segment_Line_Intersection(const TV& point_on_line,const TV& normal_of_line,T &interpolation_fraction) const
{
    T denominator=TV::Dot_Product(X.y-X.x,normal_of_line);
    if(!denominator){interpolation_fraction=FLT_MAX;return false;} // parallel
    interpolation_fraction=TV::Dot_Product(point_on_line-X.x,normal_of_line)/denominator;
    return interpolation_fraction<=1 && interpolation_fraction>=0;
}
//#####################################################################
// Function Closest_Point_On_Segment
//#####################################################################
template<class T> VECTOR<T,2> SEGMENT_2D<T>::
Closest_Point_On_Segment(const TV& point) const
{                  
    TV v=X.y-X.x;
    T denominator=TV::Dot_Product(v,v);
    if(denominator == 0) return X.x; // X.x and X.y are a single point
    else{
        T t=TV::Dot_Product(point-X.x,v)/denominator;
        if(t <= 0) return X.x;
        else if(t >= 1) return X.y;
        else{v=X.x+(X.y-X.x)*t;return v;}}
}
//#####################################################################
// Function Distance_From_Point_To_Segment
//#####################################################################
template<class T> T SEGMENT_2D<T>::
Distance_From_Point_To_Segment(const TV& point) const
{                  
    TV v=Closest_Point_On_Segment(point),d=v-point;
    return d.Magnitude();
}
//#####################################################################
// Function Closest_Point_On_Line
//#####################################################################
template<class T> VECTOR<T,2> SEGMENT_2D<T>::
Closest_Point_On_Line(const TV& point) const
{                  
    TV v=X.y-X.x;
    T denominator=TV::Dot_Product(v,v);
    if(denominator == 0) return X.x; // X.x and X.y are a single point
    else{
        T t=TV::Dot_Product(point-X.x,v)/denominator;
        v=X.x+(X.y-X.x)*t;return v;}
}
//#####################################################################
// Function Distance_From_Point_To_Line
//#####################################################################
template<class T> T SEGMENT_2D<T>::
Distance_From_Point_To_Line(const TV& point) const
{                  
    TV v=Closest_Point_On_Line(point),d=v-point;
    return d.Magnitude();
}
//#####################################################################
// Function Shortest_Vector_Between_Segments
//#####################################################################
// vector points from input segment to this segment
// not accurate as the segments become parallel
template<class T> VECTOR<T,2> SEGMENT_2D<T>::
Shortest_Vector_Between_Segments(const SEGMENT_2D<T>& segment,T& a,T& b) const
{
    TV u=X.y-X.x,v=segment.X.y-segment.X.x,w=segment.X.x-X.x;
    T u_magnitude=u.Magnitude(),v_magnitude=v.Magnitude();
    T u_dot_u=sqr(u_magnitude),v_dot_v=sqr(v_magnitude),u_dot_v=TV::Dot_Product(u,v),
               u_dot_w=TV::Dot_Product(u,w),v_dot_w=TV::Dot_Product(v,w);
    T denominator=u_dot_u*v_dot_v-sqr(u_dot_v);
    T rhs1=v_dot_v*u_dot_w-u_dot_v*v_dot_w,rhs2=u_dot_v*u_dot_w-u_dot_u*v_dot_w;
    if(rhs1 <= 0) a=0;else if(denominator < rhs1) a=1;else a=rhs1/denominator;
    if(rhs2 <= 0) b=0;else if(denominator < rhs2) b=1;else b=rhs2/denominator;
    if(a > 0 && a < 1){
        if(b == 0){a=u_dot_w/u_dot_u;if(a < 0) a=0;else if(a > 1) a=1;}
        else if(b == 1){a=TV::Dot_Product(segment.X.y-X.x,u)/u_dot_u;if(a < 0) a=0;else if(a > 1) a=1;}}
    else if(b > 0 && b < 1){
        if(a == 0){b=-v_dot_w/v_dot_v;if(b < 0) b=0;else if(b > 1) b=1;}
        else if(a == 1){b=TV::Dot_Product(X.y-segment.X.x,v)/v_dot_v;if(b < 0) b=0;else if(b > 1) b=1;}}
    else{
        T a_distance=(a == 0 ? -rhs1:rhs1-denominator)*u_magnitude,b_distance=(b == 0 ? -rhs2:rhs2-denominator)*v_magnitude;
        if(a_distance > b_distance){
            if(a == 0) b=-v_dot_w/v_dot_v;else b=TV::Dot_Product(X.y-segment.X.x,v)/v_dot_v;if(b < 0) b=0;else if(b > 1) b=1;}
        else{if(b == 0) a=u_dot_w/u_dot_u;else a=TV::Dot_Product(segment.X.y-X.x,u)/u_dot_u;if(a < 0) a=0;else if(a > 1) a=1;}}
    return a*u-w-b*v;
}
//#####################################################################
// Function Segment_Segment_Interaction
//#####################################################################
template<class T> int SEGMENT_2D<T>::
Segment_Segment_Interaction(const SEGMENT_2D<T>& segment,const TV& v1,const TV& v2,const TV& v3,const TV& v4,const T interaction_distance,
                            T& distance,TV& normal,T& a,T& b,const T small_number) const
{
    normal=Shortest_Vector_Between_Segments(segment,a,b);
    distance=normal.Magnitude();if(distance > interaction_distance) return 0; // no interaction
    TV velocity0=(1-a)*v1+a*v2,velocity1=(1-b)*v3+b*v4;
    if(distance > small_number) normal/=distance;
    else{ // set normal based on relative velocity perpendicular to the two points
        TV relative_velocity=velocity0-velocity1;
        TV u=X.y-X.x;
        normal=relative_velocity-TV::Dot_Product(relative_velocity,u)/TV::Dot_Product(u,u)*u;
        T normal_magnitude=normal.Magnitude();
        if(normal_magnitude > small_number) normal/=normal_magnitude;
        else{ // relative velocity perpendicular to the segment is 0, pick any direction perpendicular to the segment
            if(abs(u.x) > abs(u.y) && abs(u.x)) normal=TV(0,1);
            else if(abs(u.y) > abs(u.x) && abs(u.y)) normal=TV(1,0);
            else normal=TV(0,1);
            normal=normal-TV::Dot_Product(normal,u)/TV::Dot_Product(u,u)*u;normal.Normalize();
            LOG::cout << "                                            PICKING RANDOM NORMAL !!!!!!!!!!!!!!!!!!!!!!!" <<  std::endl;
    }}
    return 1;
}
//#####################################################################
// Function Segment_Segment_Collision
//#####################################################################
// Needs to be fixed
#if 0
template<class T> int SEGMENT_2D<T>::
Segment_Segment_Collision(const SEGMENT_2D<T>& segment,const TV& v1,const TV& v2,const TV& v3,const TV& v4,const T dt,
                          const T collision_thickness,T& collision_time,TV& normal,T& a,T& b,const T small_number) const
{
    // find cubic and compute the roots as possible collision times
    TV ABo=X.y-X.x,ABv=dt*(v2-v1),ACo=segment.X.y-segment.X.x,ACv=dt*(v4-v3);
    VECTOR<T,3> No_3D=TV::Cross_Product(ABo,ACo),
                             Nv_3D=TV::Cross_Product(ABo,ACv)+TV::Cross_Product(ABv,ACo),
                             Na_3D=TV::Cross_Product(ABv,ACv);
    TV No(No_3D.x,No_3D.y),Nv(Nv_3D.x,Nv_3D.y), Na(Na_3D.x,Na_3D.y);
    TV APo=segment.X.x-X.x,APv=dt*(v3-v1);
    CUBIC<T> cubic(TV::Dot_Product(Na,APv),TV::Dot_Product(Nv,APv)+TV::Dot_Product(Na,APo),
                              TV::Dot_Product(No,APv)+TV::Dot_Product(Nv,APo),TV::Dot_Product(No,APo));
    cubic.Compute_Roots_In_Interval(0,1+1e-6);

    // check the collision times
    T distance;
    for(int roots=0;roots<cubic.roots;roots++){
        collision_time=dt*cubic.root[roots-1];
        SEGMENT_2D segment2(X.x+collision_time*v1,X.y+collision_time*v2);
        if(segment2.Segment_Segment_Interaction(SEGMENT_2D(segment.X.x+collision_time*v3,segment.X.y+collision_time*v4),v1,v2,v3,v4,collision_thickness,distance,normal,a,b,
            small_number)) return 1;}
    return 0;
}
#endif
//#####################################################################
// Function Thickened_Box 
//#####################################################################
template<class T> ORIENTED_BOX<VECTOR<T,2> > SEGMENT_2D<T>::
Thickened_Oriented_Box(const T thickness_over_two) const 
{
    // make norm and tangent direction of thickness_over_two length
    TV segment_vector=X.y-X.x,segment_tangent=segment_vector.Normalized()*thickness_over_two,segment_normal(-segment_tangent.y,segment_tangent.x);
    // Form box point and intersect
    return ORIENTED_BOX<TV>(X.x-segment_tangent-segment_normal,MATRIX<T,2>(segment_vector+segment_tangent*2,segment_normal*2));
}
//#####################################################################
// Function Inside
//#####################################################################
template<class T> bool SEGMENT_2D<T>::
Inside(const TV& point,const T thickness_over_two) const 
{
    ORIENTED_BOX<TV> thickened_oriented_box=Thickened_Oriented_Box(thickness_over_two);
    return thickened_oriented_box.Lazy_Inside(point);
}
//#####################################################################
// Function Linear_Point_Inside_Segment
//#####################################################################
template<class T> bool SEGMENT_2D<T>::
Linear_Point_Inside_Segment(const TV& X,const T thickness_over_2) const
{
    TV weights=Barycentric_Coordinates(X);
    return weights.x>=-thickness_over_2 && weights.y>=-thickness_over_2;
}
//#####################################################################
// Function Point_Face_Interaction
//#####################################################################
// outputs unsigned distance
template<class T> bool SEGMENT_2D<T>::
Point_Face_Interaction(const TV& x,const T interaction_distance,const bool allow_negative_weights,T& distance) const
{      
    distance=TV::Dot_Product(x-X.x,Normal());
    return abs(distance)<=interaction_distance && Linear_Point_Inside_Segment(x,allow_negative_weights?interaction_distance:0);
}
//#####################################################################
// Function Point_Face_Interaction_Data 
//#####################################################################
template<class T> void SEGMENT_2D<T>::
Point_Face_Interaction_Data(const TV& x,T& distance,TV& interaction_normal,VECTOR<T,TV::m+1>& weights,const bool perform_attractions) const
{
    interaction_normal=Normal();weights=Barycentric_Coordinates(x).Insert(-1,0);
    if(!perform_attractions && distance<0){distance*=-1;interaction_normal*=-1;} // distance > 0, interaction_normal points from the triangle to the point
}
//#####################################################################
// Function Point_Face_Interaction
//#####################################################################
template<class T> bool SEGMENT_2D<T>::
Point_Face_Interaction(const TV& x,const TV& v,const TV& v1,const TV& v2,const T interaction_distance,T& distance,
    TV& interaction_normal,VECTOR<T,TV::m+1>& weights,const bool allow_negative_weights,const bool exit_early) const
{
    if(!Point_Face_Interaction(x,interaction_distance,allow_negative_weights,distance)) return false;
    if(!exit_early) Point_Face_Interaction_Data(x,distance,interaction_normal,weights,false);
    return true;
}
//#####################################################################
// Function Robust_Point_Face_Collision
//#####################################################################
template<class T> POINT_SIMPLEX_COLLISION_TYPE SEGMENT_2D<T>::
Robust_Point_Face_Collision(const SEGMENT_2D<T>& initial_segment,const SEGMENT_2D<T>& final_segment,const TV &x,const TV &final_x,const T dt, const T collision_thickness,
    T& collision_time,TV& normal,VECTOR<T,TV::m+1>& weights)
{
    if(final_segment.Thickened_Oriented_Box(collision_thickness).Lazy_Inside(final_x)){
        collision_time=dt;
        weights=final_segment.Barycentric_Coordinates(final_x).Insert(-1,0);
        return POINT_SIMPLEX_COLLISION_ENDS_INSIDE;}
    if(initial_segment.Thickened_Oriented_Box(collision_thickness).Lazy_Inside(x)){
        collision_time=0;
        weights=initial_segment.Barycentric_Coordinates(x).Insert(-1,0);
        return POINT_SIMPLEX_COLLISION_ENDS_OUTSIDE;}
    TV v1=(final_segment.X.x-initial_segment.X.x)/dt,v2=(final_segment.X.y-initial_segment.X.y)/dt,v=(final_x-x)/dt;
    VECTOR<T,3> tweights;
    if(initial_segment.Point_Face_Collision(x,v,v1,v2,dt,collision_thickness,collision_time,normal,tweights,false)){
        weights=tweights;
        return POINT_SIMPLEX_COLLISION_ENDS_OUTSIDE;}
    return POINT_SIMPLEX_NO_COLLISION;
}
//#####################################################################
// Function Point_Face_Collision
//#####################################################################
template<class T> bool SEGMENT_2D<T>::
Point_Face_Collision(const TV& x,const TV& v,const TV& v1,const TV& v2,const T dt,const T collision_thickness,T& collision_time,TV& normal,VECTOR<T,TV::m+1>& weights,
    const bool exit_early) const 
{
    VECTOR<double,2> v_minus_v1(v-v1),v2_minus_v1(v2-v1),x_minus_x1(x-X.x),x2_minus_x1(X.y-X.x);
    QUADRATIC<double> quadratic(VECTOR<double,2>::Cross_Product(v_minus_v1,v2_minus_v1).x,
        VECTOR<double,2>::Cross_Product(x_minus_x1,v2_minus_v1).x+VECTOR<double,2>::Cross_Product(v_minus_v1,x2_minus_x1).x,
        VECTOR<double,2>::Cross_Product(x_minus_x1,x2_minus_x1).x);

    double collision_time_temp(0),distance(0);VECTOR<double,2> normal_temp;VECTOR<double,3> weights_temp;
    if(abs(1e3*quadratic.a)<abs(quadratic.b)){
        collision_time_temp=-quadratic.c/quadratic.b;
        if(collision_time_temp<0 || collision_time_temp>dt) return false;}
    else{
        quadratic.Compute_Roots_In_Interval(0,dt);
        if(quadratic.roots==0)return false;
        else if(quadratic.roots==-1){LOG::cout<<"VERY SINGULAR ON QUADRATIC SOLVE"<<std::endl;collision_time_temp=0;}
        else if(quadratic.roots==1)collision_time_temp=quadratic.root[0];
        else collision_time_temp=quadratic.root[0];}
    SEGMENT_2D<double> segment((VECTOR<double,2>)X.x+collision_time_temp*(VECTOR<double,2>)v1,(VECTOR<double,2>)X.y+collision_time_temp*(VECTOR<double,2>)v2);
    bool interaction=segment.Point_Face_Interaction((VECTOR<double,2>)x+collision_time_temp*(VECTOR<double,2>)v,(VECTOR<double,2>)v,(VECTOR<double,2>)v1,(VECTOR<double,2>)v2,
        (double)collision_thickness,distance,normal_temp,weights_temp,true,exit_early);
    collision_time=(T)collision_time_temp; 
    normal=(TV)normal_temp;
    weights=(VECTOR<T,3>)weights_temp;
    return interaction;
}
//#####################################################################
// Function Clip_To_Box
//#####################################################################
template<class T> void SEGMENT_2D<T>::
Clip_To_Box(const RANGE<TV>& box,ARRAY<SEGMENT_2D<T> >& clipped_simplices) const
{
    // cut with all sides of box
    clipped_simplices.Remove_All();
    clipped_simplices.Append(*this);
    for(int axis=0;axis<TV::m;axis++){
        for(int i=clipped_simplices.m-1;i>=0;i--){
            SEGMENT_2D<T> triangle=clipped_simplices(i);clipped_simplices.Remove_Index_Lazy(i);
            Cut_With_Hyperplane_And_Discard_Outside_Simplices(triangle,LINE_2D<T>(-TV::Axis_Vector(axis),box.min_corner),clipped_simplices);}
        for(int i=clipped_simplices.m-1;i>=0;i--){
            SEGMENT_2D<T> triangle=clipped_simplices(i);clipped_simplices.Remove_Index_Lazy(i);
            Cut_With_Hyperplane_And_Discard_Outside_Simplices(triangle,LINE_2D<T>(TV::Axis_Vector(axis),box.max_corner),clipped_simplices);}}
}
//#####################################################################
// Function Cut_With_Hyperplane
//#####################################################################
template<class T> void SEGMENT_2D<T>::
Cut_With_Hyperplane_And_Discard_Outside_Simplices(const SEGMENT_2D<T>& segment,const LINE_2D<T>& cutting_plane,ARRAY<SEGMENT_2D<T> >& negative_segments)
{
    TV phi_nodes;
    VECTOR<TV,2> X_nodes;
    X_nodes(0)=segment.X.x;
    X_nodes(1)=segment.X.y;
    for(int i=0;i<2;i++){phi_nodes[i]=cutting_plane.Signed_Distance(X_nodes[i]);}
    int positive_count=0;
    for(int i=0;i<2;i++) if(phi_nodes[i]>0) positive_count++;
    switch(positive_count){
        case 0: // in negative halfspace
            negative_segments.Append(segment);break;
        case 1:{
            // draw positive triangle. has correct positive/negative area based on whether triangle is backwards or not
            TV interface_location=LINEAR_INTERPOLATION<T,TV>::Linear(X_nodes[0],X_nodes[1],LEVELSET_UTILITIES<T>::Theta(phi_nodes[0],phi_nodes[1]));
            if(phi_nodes[0]>0) negative_segments.Append(SEGMENT_2D<T>(interface_location,segment.X.y));
            else negative_segments.Append(SEGMENT_2D<T>(segment.X.x,interface_location));
            break;}
        case 2: // in positive halfspace
            break;}
}
//#####################################################################
// Function Cut_With_Hyperplane
//#####################################################################
template<class T> void SEGMENT_2D<T>::
Cut_With_Hyperplane(const SEGMENT_2D<T>& segment,const LINE_2D<T>& cutting_plane,ARRAY<SEGMENT_2D<T> >& negative_segments,ARRAY<SEGMENT_2D<T> >& positive_segments,T tol)
{
    TV phi_nodes;
    VECTOR<TV,2> X_nodes;
    X_nodes(0)=segment.X.x;
    X_nodes(1)=segment.X.y;
    for(int i=0;i<2;i++){phi_nodes[i]=cutting_plane.Signed_Distance(X_nodes[i]);}
    int positive_count=0;
    for(int i=0;i<2;i++) if(phi_nodes[i]>0) positive_count++;
    switch(positive_count){
        case 0: // in negative halfspace
            negative_segments.Append(segment);break;
        case 1:{
            // draw positive triangle. has correct positive/negative area based on whether triangle is backwards or not
            T theta=LEVELSET_UTILITIES<T>::Theta(phi_nodes[0],phi_nodes[1]);
            TV interface_location=LINEAR_INTERPOLATION<T,TV>::Linear(X_nodes[0],X_nodes[1],theta);
            if(phi_nodes[0]>0){
                if(1-theta>tol) negative_segments.Append(SEGMENT_2D<T>(interface_location,segment.X.y));
                if(theta>tol) positive_segments.Append(SEGMENT_2D<T>(segment.X.x,interface_location));}
            else{
                if(theta>tol) negative_segments.Append(SEGMENT_2D<T>(segment.X.x,interface_location));
                if(1-theta>tol) positive_segments.Append(SEGMENT_2D<T>(interface_location,segment.X.y));}
            break;}
        case 2: // in positive halfspace
            positive_segments.Append(segment);break;
            break;}
}
//#####################################################################
// Function Clip_To_Box_Helper
//#####################################################################
template<class T> bool
Clip_To_Box_Helper(T z0,T z1,T& a,T& b)
{
    if(z0>z1){z0=1-z0;z1=1-z1;}
    if(z0<0){
        if(z1<0) return false;
        a=max(a,z0/(z0-z1));}
    if(z1>1){
        if(z0>1) return false;
        b=min(b,(z0-1)/(z0-z1));}
    return true;
}
//#####################################################################
// Function Clip_To_Box
//#####################################################################
template<class T> bool SEGMENT_2D<T>::
Clip_To_Box(const RANGE<TV>& box,T& a,T& b) const
{
    TV z0=(X.x-box.min_corner)/box.Edge_Lengths();
    TV z1=(X.y-box.min_corner)/box.Edge_Lengths();
    a=0;
    b=1;
    return Clip_To_Box_Helper(z0.x,z1.x,a,b) && Clip_To_Box_Helper(z0.y,z1.y,a,b) && a<b;
}
//#####################################################################
namespace PhysBAM{
template class SEGMENT_2D<float>;
template class SEGMENT_2D<double>;
}
