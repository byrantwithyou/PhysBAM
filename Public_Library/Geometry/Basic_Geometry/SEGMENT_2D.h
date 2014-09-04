//#####################################################################
// Copyright 2003-2007, Zhaosheng Bao, Ronald Fedkiw, Eran Guendelman, Nipun Kwatra, Andrew Selle, Jonathan Su.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class SEGMENT_2D  
//##################################################################### 
#ifndef __SEGMENT_2D__
#define __SEGMENT_2D__

#include <Tools/Arrays/ARRAYS_FORWARD.h>
#include <Tools/Utilities/TYPE_UTILITIES.h>
#include <Tools/Vectors/VECTOR_2D.h>
#include <Geometry/Basic_Geometry/BASIC_GEOMETRY_FORWARD.h>
namespace PhysBAM{

template<class TV> class ORIENTED_BOX;

template<class T>
class SEGMENT_2D
{
    typedef VECTOR<T,2> TV;
public:
    VECTOR<TV,2> X;

    SEGMENT_2D()
        :X(TV(),TV(1,0))
    {}

    SEGMENT_2D(const TV& x1_input,const TV& x2_input)
        :X(x1_input,x2_input)
    {}

    template<class T_ARRAY>
    explicit SEGMENT_2D(const T_ARRAY& X_input)
        :X(X_input)
    {
        STATIC_ASSERT(T_ARRAY::m==2);
    }

    T Length() const
    {return (X.y-X.x).Magnitude();}

    T Size() const
    {return Length();}

    template<class T_ARRAY>
    static T Size(const T_ARRAY& X)
    {STATIC_ASSERT(T_ARRAY::m==2);return (X(1)-X(0)).Magnitude();}

    template<class T_ARRAY>
    static T Signed_Size(const T_ARRAY& X)
    {return Size(X);}

    TV Center() const
    {return (T).5*(X.x+X.y);}

    static TV Normal(const TV& x0,const TV& x1) 
    {return (x1-x0).Normalized().Rotate_Clockwise_90();}

    TV Normal() const
    {return SEGMENT_2D<T>::Normal(X.x,X.y);}

    T Signed_Distance(const TV& location) const
    {return Normal().Dot(location-X.x);}

    template<class T_ARRAY>
    static TV Normal(const T_ARRAY& X)
    {STATIC_ASSERT(T_ARRAY::m==2);return Normal(X(0),X(1));}

    static TV Barycentric_Coordinates(const TV& location,const TV& x0,const TV& x1) 
    {TV v=x1-x0;
    T denominator=TV::Dot_Product(v,v);
    if(denominator == 0) return TV(1,0); // x0 and x1 are a single point
    else{
        T t=TV::Dot_Product(location-x0,v)/denominator;
        return TV(1-t,t);}}

    static TV Clamped_Barycentric_Coordinates(const TV& location,const TV& x0,const TV& x1) 
    {TV v=x1-x0;
    T denominator=TV::Dot_Product(v,v);
    if(denominator == 0) return TV(1,0); // x0 and x1 are a single point
    else{
        T t=clamp(TV::Dot_Product(location-x0,v)/denominator,(T)0,(T)1);
        return TV(1-t,t);}}

    template<class T_ARRAY>
    static TV Clamped_Barycentric_Coordinates(const TV& location,const T_ARRAY& X)
    {STATIC_ASSERT(T_ARRAY::m==2);return Clamped_Barycentric_Coordinates(location,X(0),X(1));}

    TV Sum_Barycentric_Coordinates(const SEGMENT_2D<T>& embedded_segment) const
    {return Barycentric_Coordinates(embedded_segment.X.x)+Barycentric_Coordinates(embedded_segment.X.y);}

    TV Barycentric_Coordinates(const TV& location) const
    {return Barycentric_Coordinates(location,X.x,X.y);}

    static TV Point_From_Barycentric_Coordinates(const T alpha,const TV& x0,const TV& x1)
    {return (x1-x0)*alpha+x0;}

    template<class T_ARRAY>
    static TV Point_From_Barycentric_Coordinates(const TV& weights,const T_ARRAY& X)
    {STATIC_ASSERT(T_ARRAY::m==2);return weights.x*X(0)+weights.y*X(1);}

    TV Point_From_Barycentric_Coordinates(const T alpha) const 
    {return (X.y-X.x)*alpha+X.x;}

    TV Point_From_Barycentric_Coordinates(const VECTOR<T,2>& weights)
    {return weights.x*X.x+weights.y*X.y;}

    template<class T_ARRAY>
    static TV Point_From_Barycentric_Coordinates(const T alpha,const T_ARRAY& X)
    {STATIC_ASSERT(T_ARRAY::m==2);return Point_From_Barycentric_Coordinates(alpha,X(0),X(1));}

    bool Point_Face_Collision(const TV& x,const TV& v,const INDIRECT_ARRAY<ARRAY_VIEW<TV>,VECTOR<int,2>&> V_face,const T dt,const T collision_thickness,T& collision_time,TV& normal,
        VECTOR<T,TV::m+1>& weights,const bool exit_early=false) const
    {return Point_Face_Collision(x,v,V_face(0),V_face(1),dt,collision_thickness,collision_time,normal,weights,exit_early);}
    
    bool Point_Face_Interaction(const TV& x,const TV& v,const INDIRECT_ARRAY<ARRAY_VIEW<TV>,VECTOR<int,2>&> V_face,const T interaction_distance,T& distance,
            TV& interaction_normal,VECTOR<T,TV::m+1>& weights,const bool allow_negative_weights,const bool exit_early) const
    {return Point_Face_Interaction(x,v,V_face(0),V_face(1),interaction_distance,distance,interaction_normal,weights,allow_negative_weights,exit_early);}

    RANGE<TV> Bounding_Box() const
    {return RANGE<TV>::Bounding_Box(X.x,X.y);}

    T Distance_To_Element(const TV& location) const
    {return Distance_From_Point_To_Segment(location);}

    TV Closest_Point(const TV& point) const
    {return Closest_Point_On_Segment(point);}

//#####################################################################
    bool Segment_Line_Intersection(const TV& point_on_line,const TV& normal_of_line,T &interpolation_fraction) const;
    TV Closest_Point_On_Segment(const TV& point) const;
    T Distance_From_Point_To_Segment(const TV& point) const;
    TV Closest_Point_On_Line(const TV& point) const;
    T Distance_From_Point_To_Line(const TV& point) const;
    TV Shortest_Vector_Between_Segments(const SEGMENT_2D<T>& segment,T& a,T& b) const;
    int Segment_Segment_Interaction(const SEGMENT_2D<T>& segment,const TV& v1,const TV& v2,const TV& v3,const TV& v4,
        const T interaction_distance,T& distance,TV& normal,T& a,T& b,const T small_number=0) const;
//    int Segment_Segment_Collision(const SEGMENT_2D<T>& segment,const TV& v1,const TV& v2,const TV& v3,const TV& v4,const T dt,
//        const T collision_thickness,T& collision_time,TV& normal,T& a,T& b,const T small_number=0) const;
    ORIENTED_BOX<TV> Thickened_Oriented_Box(const T thickness_over_two=0) const;
    bool Inside(const TV& point,const T thickness_over_two=0) const;
    bool Linear_Point_Inside_Segment(const TV& X,const T thickness_over_2) const;
    static POINT_SIMPLEX_COLLISION_TYPE Robust_Point_Face_Collision(const SEGMENT_2D<T>& initial_segment,const SEGMENT_2D<T>& final_segment,const TV& x,
        const TV& final_x,const T dt,const T collision_thickness,T& collision_time,TV& normal,VECTOR<T,TV::m+1>& weights);
    bool Point_Face_Collision(const TV& x,const TV& v,const TV& v1,const TV& v2,const T dt,const T collision_thickness,T& collision_time,TV& normal,VECTOR<T,TV::m+1>& weights,
        const bool exit_early=false) const;
    bool Point_Face_Interaction(const TV& x,const T interaction_distance,const bool allow_negative_weights,T& distance) const;
    void Point_Face_Interaction_Data(const TV& x,T& distance,TV& interaction_normal,VECTOR<T,TV::m+1>& weights,const bool perform_attractions) const;
    bool Point_Face_Interaction(const TV& x,const TV& v,const TV& v1,const TV& v2,const T interaction_distance,T& distance,
        TV& interaction_normal,VECTOR<T,TV::m+1>& weights,const bool allow_negative_weights,const bool exit_early) const;
    void Clip_To_Box(const RANGE<TV>& box,ARRAY<SEGMENT_2D<T> >& clipped_simplices) const;
    static void Cut_With_Hyperplane_And_Discard_Outside_Simplices(const SEGMENT_2D<T>& segment,const LINE_2D<T>& cutting_plane,ARRAY<SEGMENT_2D<T> >& negative_segments);
    static void Cut_With_Hyperplane(const SEGMENT_2D<T>& segment,const LINE_2D<T>& cutting_plane,ARRAY<SEGMENT_2D<T> >& negative_segments,ARRAY<SEGMENT_2D<T> >& positive_segments,T tol=0);
    bool Clip_To_Box(const RANGE<TV>& box,T& a,T& b) const;
//#####################################################################
};

template<class T> std::ostream &operator<<(std::ostream &output,const SEGMENT_2D<T> &segment)
{output << segment.X.x << ", " << segment.X.y;return output;}

}
#endif

