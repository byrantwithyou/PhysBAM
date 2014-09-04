//#####################################################################
// Copyright 2002-2007, Robert Bridson, Ronald Fedkiw, Geoffrey Irving, Sergey Koltakov, Nipun Kwatra, Andrew Selle, Jonathan Su, Joseph Teran.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class SEGMENT_3D  
//##################################################################### 
#ifndef __SEGMENT_3D__
#define __SEGMENT_3D__

#include <Tools/Arrays/ARRAYS_FORWARD.h>
#include <Tools/Vectors/VECTOR_3D.h>
namespace PhysBAM{

template<class T>
class SEGMENT_3D
{
    typedef VECTOR<T,3> TV;
public:
    VECTOR<TV,2> X;

    SEGMENT_3D()
        :X(TV(),TV(1,0,0))
    {}

    SEGMENT_3D(const TV& x1_input,const TV& x2_input)
        :X(x1_input,x2_input)
    {}

    template<class T_ARRAY>
    SEGMENT_3D(const T_ARRAY& X_input)
        :X(X_input(0),X_input(1))
    {
        STATIC_ASSERT(T_ARRAY::m==2);
    }

    T Length() const
    {return (X.y-X.x).Magnitude();}

    T Size() const
    {return Length();}

    static T Size(const TV& x0,const TV& x1)
    {return (x1-x0).Magnitude();}

    template<class T_ARRAY>
    static T Signed_Size(const T_ARRAY& X)
    {STATIC_ASSERT(T_ARRAY::m==2);return Size(X(0),X(1));}

    template<class T_ARRAY>
    static T Size(const T_ARRAY& X)
    {return Signed_Size(X);}

    static TV Normal(const TV& x0,const TV& x1,const TV& x2) 
    {TV v=x1-x0,face_normal=TV::Cross_Product(v,x2-x0);
    return TV::Cross_Product(face_normal,v).Normalized();} // rotate by 90 degrees clockwise

    static TV Normal_Direction(const TV& x0,const TV& x1,const TV& x2) 
    {TV v=x1-x0,face_normal=TV::Cross_Product(v,x2-x0);
    return TV::Cross_Product(face_normal,v);} // can have any magnitude

    template<class T_ARRAY>
    void Edge_Edge_Interaction_Data(const SEGMENT_3D<T>& segment,const INDIRECT_ARRAY<T_ARRAY,VECTOR<int,4>&> V_edges,const T& distance,TV& normal,const VECTOR<T,TV::m+1>& weights,
        const T small_number=0,const bool verbose=true) const
    {Edge_Edge_Interaction_Data(segment,V_edges(0),V_edges(1),V_edges(2),V_edges(3),distance,normal,weights,small_number,verbose);}

    template<class T_ARRAY>
    bool Edge_Edge_Interaction(const SEGMENT_3D<T>& segment,const INDIRECT_ARRAY<T_ARRAY,VECTOR<int,4>&> V_edges,const T interaction_distance,T& distance,TV& normal,
        VECTOR<T,TV::m+1>& weights,bool allow_negative_weights,const T small_number=0,const bool exit_early=false,const bool verbose=true) const
    {return Edge_Edge_Interaction(segment,V_edges(0),V_edges(1),V_edges(2),V_edges(3),interaction_distance,distance,normal,weights,allow_negative_weights,small_number,exit_early,verbose);}

    template<class T_ARRAY>
    bool Edge_Edge_Collision(const SEGMENT_3D<T>& segment,const INDIRECT_ARRAY<T_ARRAY,VECTOR<int,4>&> V_edges,const T dt,const T collision_thickness,T& collision_time,
        TV& normal,VECTOR<T,TV::m+1>& weights,const T small_number=0,const bool exit_early=false) const
    {return Edge_Edge_Collision(segment,V_edges(0),V_edges(1),V_edges(2),V_edges(3),dt,collision_thickness,collision_time,normal,weights,small_number,exit_early);}

    T Distance_To_Element(const TV& location) const
    {return Distance_From_Point_To_Segment(location);}

//#####################################################################
    TV Closest_Point_On_Segment(const TV& point) const;
    T Distance_From_Point_To_Segment(const TV& point) const;
    TV Closest_Point_On_Line(const TV& point) const;
    T Distance_From_Point_To_Line(const TV& point) const;
    TV Shortest_Vector_Between_Lines(const SEGMENT_3D<T>& segment,VECTOR<T,2>& weights) const;
    TV Shortest_Vector_Between_Segments(const SEGMENT_3D<T>& segment,VECTOR<T,2>& weights) const;
    bool Edge_Edge_Interaction(const SEGMENT_3D<T>& segment,const T interaction_distance,T& distance,TV& normal,VECTOR<T,TV::m+1>& weights,bool allow_negative_weights) const;
    void Edge_Edge_Interaction_Data(const SEGMENT_3D<T>& segment,const TV& v1,const TV& v2,const TV& v3,const TV& v4,
        const T& distance,TV& normal,const VECTOR<T,TV::m+1>& weights,const T small_number=0,const bool verbose=true) const;
    bool Edge_Edge_Interaction(const SEGMENT_3D<T>& segment,const TV& v1,const TV& v2,const TV& v3,const TV& v4,
        const T interaction_distance,T& distance,TV& normal,VECTOR<T,TV::m+1>& weights,bool allow_negative_weights,const T small_number=0,
        const bool exit_early=false,const bool verbose=true) const;
    bool Edge_Edge_Collision(const SEGMENT_3D<T>& segment,const TV& v1,const TV& v2,const TV& v3,const TV& v4,const T dt,
        const T collision_thickness,T& collision_time,TV& normal,VECTOR<T,TV::m+1>& weights,const T small_number=0,const bool exit_early=false) const;
    T Interpolation_Fraction(const TV& location) const;
    VECTOR<T,2> Barycentric_Coordinates(const TV& location) const;
    VECTOR<T,2> Clamped_Barycentric_Coordinates(const TV& location,const T tolerance=1e-7) const;
    static T Interpolation_Fraction(const TV& location,const TV& x0,const TV& x1);
    static VECTOR<T,2> Barycentric_Coordinates(const TV& location,const TV& x0,const TV& x1);
    static VECTOR<T,2> Clamped_Barycentric_Coordinates(const TV& location,const TV& x0,const TV& x1,const T tolerance=1e-7);
    bool Inside(const TV& point,const T thickness_over_two=0) const;
//#####################################################################
};   
}
#endif
