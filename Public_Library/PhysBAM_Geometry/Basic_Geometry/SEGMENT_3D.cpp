//#####################################################################
// Copyright 2002-2007, Robert Bridson, Ronald Fedkiw, Sergey Koltakov, Andrew Selle, Joseph Teran.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class SEGMENT_3D  
//##################################################################### 
#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Tools/Math_Tools/clamp.h>
#include <PhysBAM_Tools/Math_Tools/INTERVAL.h>
#include <PhysBAM_Tools/Nonlinear_Equations/ITERATIVE_SOLVER.h>
#include <PhysBAM_Tools/Polynomials/CUBIC.h>
#include <PhysBAM_Tools/Vectors/VECTOR.h>
#include <PhysBAM_Geometry/Basic_Geometry/SEGMENT_3D.h>
using namespace PhysBAM;
//#####################################################################
// Function Closest_Point_On_Segment
//#####################################################################
template<class T> VECTOR<T,3> SEGMENT_3D<T>::
Closest_Point_On_Segment(const VECTOR<T,3>& point) const
{                  
    VECTOR<T,3> v=X.y-X.x;
    T denominator=VECTOR<T,3>::Dot_Product(v,v);
    if(denominator == 0) return X.x; // X.x and X.y are a single point
    else{
        T t=VECTOR<T,3>::Dot_Product(point-X.x,v)/denominator;
        if(t <= 0) return X.x;
        else if(t >= 1) return X.y;
        else{v=X.x+(X.y-X.x)*t;return v;}}
}
//#####################################################################
// Function Distance_From_Point_To_Segment
//#####################################################################
template<class T> T SEGMENT_3D<T>::
Distance_From_Point_To_Segment(const VECTOR<T,3>& point) const
{                  
    VECTOR<T,3> v=Closest_Point_On_Segment(point),d=v-point;
    return d.Magnitude();
}
//#####################################################################
// Function Closest_Point_On_Line
//#####################################################################
template<class T> VECTOR<T,3> SEGMENT_3D<T>::
Closest_Point_On_Line(const VECTOR<T,3>& point) const
{                  
    VECTOR<T,3> v=X.y-X.x;
    T denominator=VECTOR<T,3>::Dot_Product(v,v);
    if(denominator == 0) return X.x; // X.x and X.y are a single point
    else{
        T t=VECTOR<T,3>::Dot_Product(point-X.x,v)/denominator;
        v=X.x+(X.y-X.x)*t;return v;}
}
//#####################################################################
// Function Distance_From_Point_To_Line
//#####################################################################
template<class T> T SEGMENT_3D<T>::
Distance_From_Point_To_Line(const VECTOR<T,3>& point) const
{                  
    VECTOR<T,3> v=Closest_Point_On_Line(point),d=v-point;
    return d.Magnitude();
}
//#####################################################################
// Function Shortest_Vector_Between_Lines
//#####################################################################
// vector points from argument segment to class segment; not accurate as the segments become parallel
template<class T> VECTOR<T,3> SEGMENT_3D<T>::
Shortest_Vector_Between_Lines(const SEGMENT_3D<T>& segment,VECTOR<T,2>& weights) const
{
    VECTOR<T,3> u=X.y-X.x,v=segment.X.y-segment.X.x,w=segment.X.x-X.x;
    T u_magnitude_squared=u.Magnitude_Squared(),v_magnitude_squared=v.Magnitude_Squared(),u_dot_u=u_magnitude_squared,v_dot_v=v_magnitude_squared,u_dot_v=VECTOR<T,3>::Dot_Product(u,v),
        u_dot_w=VECTOR<T,3>::Dot_Product(u,w),v_dot_w=VECTOR<T,3>::Dot_Product(v,w);
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
    VECTOR<T,3> u=X.y-X.x,v=segment.X.y-segment.X.x,w=segment.X.x-X.x;
    T u_magnitude_squared=u.Magnitude_Squared(),v_magnitude_squared=v.Magnitude_Squared(),u_dot_u=u_magnitude_squared,v_dot_v=v_magnitude_squared,u_dot_v=VECTOR<T,3>::Dot_Product(u,v),
        u_dot_w=VECTOR<T,3>::Dot_Product(u,w),v_dot_w=VECTOR<T,3>::Dot_Product(v,w);
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
Edge_Edge_Interaction(const SEGMENT_3D<T>& segment,const T interaction_distance,T& distance,VECTOR<T,3>& normal,VECTOR<T,TV::m+1>& weights,bool allow_negative_weights) const
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
        VECTOR<T,3> relative_velocity=-weights(0)*v1-weights(1)*v2-weights(2)*v3-weights(3)*v4,u=X.y-X.x;
        normal=relative_velocity-VECTOR<T,3>::Dot_Product(relative_velocity,u)/VECTOR<T,3>::Dot_Product(u,u)*u;
        T normal_magnitude=normal.Magnitude();
        if(normal_magnitude > small_number) normal/=normal_magnitude;
        else{ // relative velocity perpendicular to the segment is 0, pick any direction perpendicular to the segment
            if(abs(u.x) > abs(u.y) && abs(u.x) > abs(u.z)) normal=VECTOR<T,3>(0,1,1);
            else if(abs(u.y) > abs(u.x) && abs(u.y) > abs(u.z)) normal=VECTOR<T,3>(1,0,1);
            else normal=VECTOR<T,3>(1,1,0);
            normal=normal-VECTOR<T,3>::Dot_Product(normal,u)/VECTOR<T,3>::Dot_Product(u,u)*u;normal.Normalize();
            if(verbose) LOG::cout << "                                            PICKING RANDOM NORMAL !!!!!!!!!!!!!!!!!!!!!!!" << std::endl;
    }}
}
//#####################################################################
// Function Edge_Edge_Interaction
//#####################################################################
template<class T> bool SEGMENT_3D<T>::
Edge_Edge_Interaction(const SEGMENT_3D<T>& segment,const TV& v1,const TV& v2,const TV& v3,const TV& v4,const T interaction_distance,
    T& distance,VECTOR<T,3>& normal,VECTOR<T,TV::m+1>& weights,bool allow_negative_weights,const T small_number,const bool exit_early,const bool verbose) const
{
    if(!Edge_Edge_Interaction(segment,interaction_distance,distance,normal,weights,allow_negative_weights)) return false;
    if(!exit_early) Edge_Edge_Interaction_Data(segment,v1,v2,v3,v4,distance,normal,weights,small_number,verbose);
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
    VECTOR<T,3> ABo=X.y-X.x,ABv=dt*(v2-v1),ACo=segment.X.y-segment.X.x,ACv=dt*(v4-v3);
    VECTOR<T,3> No=VECTOR<T,3>::Cross_Product(ABo,ACo),Nv=VECTOR<T,3>::Cross_Product(ABo,ACv)+VECTOR<T,3>::Cross_Product(ABv,ACo),Na=VECTOR<T,3>::Cross_Product(ABv,ACv);
    VECTOR<T,3> APo=segment.X.x-X.x,APv=dt*(v3-v1);
    
    CUBIC<double> cubic((double)VECTOR<T,3>::Dot_Product(Na,APv),(double)VECTOR<T,3>::Dot_Product(Nv,APv)+VECTOR<T,3>::Dot_Product(Na,APo),
                                       (double)VECTOR<T,3>::Dot_Product(No,APv)+VECTOR<T,3>::Dot_Product(Nv,APo),(double)VECTOR<T,3>::Dot_Product(No,APo));
    double xmin=0,xmax=1.000001;
    int num_intervals=0;VECTOR<INTERVAL<double>,3> intervals;
    cubic.Compute_Intervals(xmin,xmax,num_intervals,intervals(0),intervals(1),intervals(2));
    if(!num_intervals) return false;

    // find and check roots
    T distance;
    ITERATIVE_SOLVER<double> iterative_solver;iterative_solver.tolerance=1e-14;
    for(int k=0;k<num_intervals;k++){
        T collision_time=dt*(T)iterative_solver.Bisection_Secant_Root(cubic,intervals(k).min_corner,intervals(k).max_corner);
        SEGMENT_3D<T> segment2(X.x+collision_time*v1,X.y+collision_time*v2);
        if(segment2.Edge_Edge_Interaction(SEGMENT_3D<T>(segment.X.x+collision_time*v3,segment.X.y+collision_time*v4),v1,v2,v3,v4,collision_thickness,distance,normal,weights,
                false,small_number,exit_early)) return true;}

    return false;
}
//#####################################################################
// Function Interpolation_Fraction
//#####################################################################
template<class T> T SEGMENT_3D<T>::
Interpolation_Fraction(const VECTOR<T,3>& location) const
{  
    return SEGMENT_3D::Interpolation_Fraction(location,X.x,X.y);
}
//#####################################################################
// Function Barycentric_Coordinates
//#####################################################################
template<class T> VECTOR<T,2> SEGMENT_3D<T>::
Barycentric_Coordinates(const VECTOR<T,3>& location) const
{  
    return SEGMENT_3D::Barycentric_Coordinates(location,X.x,X.y);
}
//#####################################################################
// Function Clamped_Barycentric_Coordinates
//#####################################################################
template<class T> VECTOR<T,2> SEGMENT_3D<T>::
Clamped_Barycentric_Coordinates(const VECTOR<T,3>& location,const T tolerance) const
{  
    return SEGMENT_3D::Clamped_Barycentric_Coordinates(location,X.x,X.y);
}
//#####################################################################
// Function Interpolation_Fraction
//#####################################################################
template<class T> T SEGMENT_3D<T>::
Interpolation_Fraction(const VECTOR<T,3>& location,const VECTOR<T,3>& x1,const VECTOR<T,3>& x2) 
{  
    VECTOR<T,3> v=x2-x1;
    T denominator=VECTOR<T,3>::Dot_Product(v,v);
    if(denominator==0) return 0; // x1 and x2 are a single point
    else return VECTOR<T,3>::Dot_Product(location-x1,v)/denominator;
}
//#####################################################################
// Function Barycentric_Coordinates
//#####################################################################
template<class T> VECTOR<T,2> SEGMENT_3D<T>::
Barycentric_Coordinates(const VECTOR<T,3>& location,const VECTOR<T,3>& x1,const VECTOR<T,3>& x2)
{  
    T t=Interpolation_Fraction(location,x1,x2);
    return VECTOR<T,2>(1-t,t);
}
//#####################################################################
// Function Clamped_Barycentric_Coordinates
//#####################################################################
template<class T> VECTOR<T,2> SEGMENT_3D<T>::
Clamped_Barycentric_Coordinates(const VECTOR<T,3>& location,const VECTOR<T,3>& x1,const VECTOR<T,3>& x2,const T tolerance)
{  
    VECTOR<T,3> v=x2-x1;
    T denominator=VECTOR<T,3>::Dot_Product(v,v);
    if(abs(denominator)<tolerance) return VECTOR<T,2>((T).5,(T).5);
    T a=clamp(VECTOR<T,3>::Dot_Product(location-x1,v)/denominator,(T)0,(T)1);
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
