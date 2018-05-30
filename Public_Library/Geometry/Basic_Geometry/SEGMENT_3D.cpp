//#####################################################################
// Copyright 2002-2007, Robert Bridson, Ronald Fedkiw, Sergey Koltakov, Andrew Selle, Joseph Teran.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class SEGMENT_3D  
//##################################################################### 
#include <Core/Arrays/ARRAY.h>
#include <Core/Log/LOG.h>
#include <Core/Math_Tools/clamp.h>
#include <Core/Math_Tools/INTERVAL.h>
#include <Core/Vectors/VECTOR.h>
#include <Tools/Nonlinear_Equations/ITERATIVE_SOLVER.h>
#include <Tools/Polynomials/CUBIC.h>
#include <Geometry/Basic_Geometry/SEGMENT_3D.h>
using namespace PhysBAM;
//#####################################################################
// Function Closest_Point_On_Segment
//#####################################################################
template<class T> VECTOR<T,3> SEGMENT_3D<T>::
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
template<class T> T SEGMENT_3D<T>::
Distance_From_Point_To_Segment(const TV& point) const
{                  
    TV v=Closest_Point_On_Segment(point),d=v-point;
    return d.Magnitude();
}
//#####################################################################
// Function Closest_Point_On_Line
//#####################################################################
template<class T> VECTOR<T,3> SEGMENT_3D<T>::
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
template<class T> T SEGMENT_3D<T>::
Distance_From_Point_To_Line(const TV& point) const
{                  
    TV v=Closest_Point_On_Line(point),d=v-point;
    return d.Magnitude();
}
//#####################################################################
// Function Shortest_Vector_Between_Lines
//#####################################################################
// vector points from argument segment to class segment; not accurate as the segments become parallel
template<class T> VECTOR<T,3> SEGMENT_3D<T>::
Shortest_Vector_Between_Lines(const SEGMENT_3D<T>& segment,VECTOR<T,2>& weights) const
{
    TV u=X.y-X.x,v=segment.X.y-segment.X.x,w=segment.X.x-X.x;
    T u_magnitude_squared=u.Magnitude_Squared(),v_magnitude_squared=v.Magnitude_Squared(),u_dot_u=u_magnitude_squared,v_dot_v=v_magnitude_squared,u_dot_v=TV::Dot_Product(u,v),
        u_dot_w=TV::Dot_Product(u,w),v_dot_w=TV::Dot_Product(v,w);
    T denominator=u_dot_u*v_dot_v-sqr(u_dot_v),rhs1=v_dot_v*u_dot_w-u_dot_v*v_dot_w,rhs2=u_dot_v*u_dot_w-u_dot_u*v_dot_w;
    weights.x=rhs1/denominator;weights.y=rhs2/denominator;
    return weights.x*u-w-weights.y*v;
}
//#####################################################################
// Function Shortest_Vector_Between_Segments
//#####################################################################
// vector points from argument segment to class segment; not accurate as the segments become parallel
template<class T> VECTOR<T,3> SEGMENT_3D<T>::
Shortest_Vector_Between_Segments(const SEGMENT_3D<T>& segment,VECTOR<T,2>& weights) const
{
    TV u=X.y-X.x,v=segment.X.y-segment.X.x,w=segment.X.x-X.x;
    T u_magnitude_squared=u.Magnitude_Squared(),v_magnitude_squared=v.Magnitude_Squared(),u_dot_u=u_magnitude_squared,v_dot_v=v_magnitude_squared,u_dot_v=TV::Dot_Product(u,v),
        u_dot_w=TV::Dot_Product(u,w),v_dot_w=TV::Dot_Product(v,w);
    T denominator=u_dot_u*v_dot_v-sqr(u_dot_v),rhs1=v_dot_v*u_dot_w-u_dot_v*v_dot_w,rhs2=u_dot_v*u_dot_w-u_dot_u*v_dot_w;
    bool check_boundary=false;
    if(rhs1 <= 0 || denominator <= rhs1) check_boundary=true;else weights.x=rhs1/denominator; 
    if(rhs2 <= 0 || denominator <= rhs2) check_boundary=true;else weights.y=rhs2/denominator; 
    if(check_boundary){ // check boundaries of [0,1]x[0,1] weights domain
        T v_plus_w_dot_u=u_dot_v+u_dot_w,u_minus_w_dot_v=u_dot_v-v_dot_w,distance_squared_minus_w_dot_w;
        weights.x=0; // check weights.x=0 side
        if(v_dot_w>=0){distance_squared_minus_w_dot_w=0;weights.y=0;}
        else if(v_dot_v<=-v_dot_w){distance_squared_minus_w_dot_w=v_dot_v+2*v_dot_w;weights.y=1;}
        else{weights.y=-v_dot_w/v_dot_v;distance_squared_minus_w_dot_w=weights.y*v_dot_w;}
        // check weights.x=1 side
        if(u_minus_w_dot_v<=0){T new_distance_squared=u_dot_u-2*u_dot_w;
            if(new_distance_squared<distance_squared_minus_w_dot_w){distance_squared_minus_w_dot_w=new_distance_squared;weights=VECTOR<T,2>(1,0);}}
        else if(v_dot_v<=u_minus_w_dot_v){T new_distance_squared=v_dot_v+2*(v_dot_w-v_plus_w_dot_u)+u_dot_u;
            if(new_distance_squared<distance_squared_minus_w_dot_w){distance_squared_minus_w_dot_w=new_distance_squared;weights=VECTOR<T,2>(1,1);}}
        else{T weights_y_temp=u_minus_w_dot_v/v_dot_v,new_distance_squared=u_dot_u-2*u_dot_w-weights_y_temp*u_minus_w_dot_v;
            if(new_distance_squared<distance_squared_minus_w_dot_w){distance_squared_minus_w_dot_w=new_distance_squared;weights=VECTOR<T,2>(1,weights_y_temp);}}
        // check weights.y=0 side ignoring corners (already handled above)
        if(u_dot_w>0 && u_dot_u>u_dot_w){T weights_x_temp=u_dot_w/u_dot_u,new_distance_squared=-weights_x_temp*u_dot_w;
            if(new_distance_squared<distance_squared_minus_w_dot_w){distance_squared_minus_w_dot_w=new_distance_squared;weights=VECTOR<T,2>(weights_x_temp,0);}}
        // check weights.y=1 side ignoring corners (already handled above)
        if(v_plus_w_dot_u>0 && u_dot_u>v_plus_w_dot_u){T weights_x_temp=v_plus_w_dot_u/u_dot_u,new_distance_squared=v_dot_v+2*v_dot_w-weights_x_temp*v_plus_w_dot_u;
            if(new_distance_squared<distance_squared_minus_w_dot_w){distance_squared_minus_w_dot_w=new_distance_squared;weights=VECTOR<T,2>(weights_x_temp,1);}}}
    return weights.x*u-w-weights.y*v;
}
//#####################################################################
// Function Edge_Edge_Interaction
//#####################################################################
// returns the distance, normal and weights
template<class T> bool SEGMENT_3D<T>::
Edge_Edge_Interaction(const SEGMENT_3D<T>& segment,const T interaction_distance,T& distance,TV& normal,VECTOR<T,TV::m+1>& weights,bool allow_negative_weights) const
{
    VECTOR<T,2> w;
    normal=Shortest_Vector_Between_Segments(segment,w);
    weights=VECTOR<T,TV::m+1>(-w.x,w.x-1,w.y,1-w.y);
    distance=normal.Magnitude();
    if(!allow_negative_weights && (weights.Contains(0) || weights.Contains(1))) return false;
    return distance<=interaction_distance;
}
//#####################################################################
// Function template<cl
//#####################################################################
// input the distance, normal and weights
template<class T> void SEGMENT_3D<T>::
Edge_Edge_Interaction_Data(const SEGMENT_3D<T>& segment,const TV& v1,const TV& v2,const TV& v3,const TV& v4,const T& distance,
    TV& normal,const VECTOR<T,TV::m+1>& weights,const T small_number,const bool verbose) const
{
    if(distance > small_number) normal/=distance;
    else{ // set normal based on relative velocity perpendicular to the two points
        TV relative_velocity=-weights(0)*v1-weights(1)*v2-weights(2)*v3-weights(3)*v4,u=X.y-X.x;
        normal=relative_velocity-TV::Dot_Product(relative_velocity,u)/TV::Dot_Product(u,u)*u;
        T normal_magnitude=normal.Magnitude();
        if(normal_magnitude > small_number) normal/=normal_magnitude;
        else{ // relative velocity perpendicular to the segment is 0, pick any direction perpendicular to the segment
            if(abs(u.x) > abs(u.y) && abs(u.x) > abs(u.z)) normal=TV(0,1,1);
            else if(abs(u.y) > abs(u.x) && abs(u.y) > abs(u.z)) normal=TV(1,0,1);
            else normal=TV(1,1,0);
            normal=normal-TV::Dot_Product(normal,u)/TV::Dot_Product(u,u)*u;normal.Normalize();
            if(verbose) LOG::cout << "                                            PICKING RANDOM NORMAL !!!!!!!!!!!!!!!!!!!!!!!" << std::endl;
    }}
}
//#####################################################################
// Function Edge_Edge_Interaction
//#####################################################################
template<class T> bool SEGMENT_3D<T>::
Edge_Edge_Interaction(const SEGMENT_3D<T>& segment,const TV& v1,const TV& v2,const TV& v3,const TV& v4,const T interaction_distance,
    T& distance,TV& normal,VECTOR<T,TV::m+1>& weights,bool allow_negative_weights,const T small_number,const bool exit_early,const bool verbose) const
{
    if(!Edge_Edge_Interaction(segment,interaction_distance,distance,normal,weights,allow_negative_weights)) return false;
    if(!exit_early) Edge_Edge_Interaction_Data(segment,v1,v2,v3,v4,distance,normal,weights,small_number,verbose);
    return true;
}
//#####################################################################
// Function Edge_Edge_Interaction_Velocity
//#####################################################################
template<class T> bool SEGMENT_3D<T>::
Edge_Edge_Interaction_Velocity(const SEGMENT_3D<T>& segment,const TV& v1,const TV& v2,const TV& v3,const TV& v4,const T interaction_distance,
    T& distance,TV& normal,VECTOR<T,TV::m+1>& weights,T& relative_speed,const T small_number,const bool exit_early,const bool verbose) const
{
    Edge_Edge_Interaction(segment,interaction_distance,distance,normal,weights,false);
    Edge_Edge_Interaction_Data(segment,v1,v2,v3,v4,distance,normal,weights,small_number,verbose);
    relative_speed=TV::Dot_Product(weights(0)*v1+weights(1)*v2-weights(2)*v3+weights(3)*v4,normal); // relative speed is in the normal direction
    return true;
}
//#####################################################################
// Function Edge_Edge_Collision
//#####################################################################
template<class T> bool SEGMENT_3D<T>::
Edge_Edge_Collision(const SEGMENT_3D<T>& segment,const TV& v1,const TV& v2,const TV& v3,const TV& v4,const T dt,const T collision_thickness,T& collision_time,TV& normal,
    VECTOR<T,TV::m+1>& weights,const T small_number,const bool exit_early) const
{
    // find cubic and compute the roots as possible collision times
    TV ABo=X.y-X.x,ABv=dt*(v2-v1),ACo=segment.X.y-segment.X.x,ACv=dt*(v4-v3);
    TV No=TV::Cross_Product(ABo,ACo),Nv=TV::Cross_Product(ABo,ACv)+TV::Cross_Product(ABv,ACo),Na=TV::Cross_Product(ABv,ACv);
    TV APo=segment.X.x-X.x,APv=dt*(v3-v1);
    
    CUBIC<double> cubic((double)TV::Dot_Product(Na,APv),(double)TV::Dot_Product(Nv,APv)+TV::Dot_Product(Na,APo),
                                       (double)TV::Dot_Product(No,APv)+TV::Dot_Product(Nv,APo),(double)TV::Dot_Product(No,APo));
    double xmin=0,xmax=1.000001;
    int num_intervals=0;
    INTERVAL<double> intervals[3];
    cubic.Compute_Intervals(xmin,xmax,num_intervals,intervals);
    if(!num_intervals) return false;

    // find and check roots
    T distance;
    for(int k=0;k<num_intervals;k++){
        collision_time=dt*(T)Bisection_Secant_Root<double>(cubic,intervals[k].min_corner,intervals[k].max_corner);
        SEGMENT_3D<T> segment2(X.x+collision_time*v1,X.y+collision_time*v2);
        if(segment2.Edge_Edge_Interaction(SEGMENT_3D<T>(segment.X.x+collision_time*v3,segment.X.y+collision_time*v4),v1,v2,v3,v4,collision_thickness,distance,normal,weights,
                false,small_number,exit_early)) return true;}

    return false;
}
//#####################################################################
// Function Interpolation_Fraction
//#####################################################################
template<class T> T SEGMENT_3D<T>::
Interpolation_Fraction(const TV& location) const
{  
    return SEGMENT_3D::Interpolation_Fraction(location,X.x,X.y);
}
//#####################################################################
// Function Barycentric_Coordinates
//#####################################################################
template<class T> VECTOR<T,2> SEGMENT_3D<T>::
Barycentric_Coordinates(const TV& location) const
{  
    return SEGMENT_3D::Barycentric_Coordinates(location,X.x,X.y);
}
//#####################################################################
// Function Clamped_Barycentric_Coordinates
//#####################################################################
template<class T> VECTOR<T,2> SEGMENT_3D<T>::
Clamped_Barycentric_Coordinates(const TV& location,const T tolerance) const
{  
    return SEGMENT_3D::Clamped_Barycentric_Coordinates(location,X.x,X.y);
}
//#####################################################################
// Function Interpolation_Fraction
//#####################################################################
template<class T> T SEGMENT_3D<T>::
Interpolation_Fraction(const TV& location,const TV& x0,const TV& x1) 
{  
    TV v=x1-x0;
    T denominator=TV::Dot_Product(v,v);
    if(denominator==0) return 0; // x0 and x1 are a single point
    else return TV::Dot_Product(location-x0,v)/denominator;
}
//#####################################################################
// Function Barycentric_Coordinates
//#####################################################################
template<class T> VECTOR<T,2> SEGMENT_3D<T>::
Barycentric_Coordinates(const TV& location,const TV& x0,const TV& x1)
{  
    T t=Interpolation_Fraction(location,x0,x1);
    return VECTOR<T,2>(1-t,t);
}
//#####################################################################
// Function Clamped_Barycentric_Coordinates
//#####################################################################
template<class T> VECTOR<T,2> SEGMENT_3D<T>::
Clamped_Barycentric_Coordinates(const TV& location,const TV& x0,const TV& x1,const T tolerance)
{  
    TV v=x1-x0;
    T denominator=TV::Dot_Product(v,v);
    if(abs(denominator)<tolerance) return VECTOR<T,2>((T).5,(T).5);
    T a=clamp(TV::Dot_Product(location-x0,v)/denominator,(T)0,(T)1);
    return VECTOR<T,2>((T)1-a,a);
}
//#####################################################################
// Function Inside
//#####################################################################
template<class T> bool SEGMENT_3D<T>::
Inside(const TV& point,const T thickness_over_two) const
{
    return Distance_From_Point_To_Segment(point)<thickness_over_two;
}
//#####################################################################
namespace PhysBAM{
template class SEGMENT_3D<float>;
template class SEGMENT_3D<double>;
}
