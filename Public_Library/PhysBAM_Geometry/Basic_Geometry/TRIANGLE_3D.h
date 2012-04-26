//#####################################################################
// Copyright 2002-2007, Robert Bridson, Ronald Fedkiw, Eran Guendelman, Geoffrey Irving, Nipun Kwatra, Neil Molino, Andrew Selle, Jonathan Su, Joseph Teran.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class TRIANGLE_3D
//##################################################################### 
#ifndef __TRIANGLE_3D__
#define __TRIANGLE_3D__

#include <PhysBAM_Tools/Math_Tools/constants.h>
#include <PhysBAM_Tools/Math_Tools/RANGE.h>
#include <PhysBAM_Tools/Utilities/TYPE_UTILITIES.h>
#include <PhysBAM_Geometry/Basic_Geometry/PLANE.h>
#include <PhysBAM_Geometry/Basic_Geometry/SEGMENT_3D.h>
namespace PhysBAM{

template<class T>
class TRIANGLE_3D
{
    typedef VECTOR<T,3> TV;
public:
    VECTOR<TV,3> X; // clockwise order when looking at the plane

    TRIANGLE_3D()
        :X(TV(0,0,0),TV(0,1,0),TV(1,0,0))
    {
    }

    TRIANGLE_3D(const TV& x1_input,const TV& x2_input,const TV& x3_input)
        :X(x1_input,x2_input,x3_input)
    {
    }

    template<class T_ARRAY>
    TRIANGLE_3D(const T_ARRAY& X_input)
        :X(X_input)
    {
        STATIC_ASSERT(T_ARRAY::m==3);
    }

    T Area() const 
    {return Area(X.x,X.y,X.z);}
    
    static T Area(const TV& x1,const TV& x2,const TV& x3) // always positive for clockwise vertices: x1, x2, x3 
    {return (T).5*TV::Cross_Product(x2-x1,x3-x1).Magnitude();}
    
    static T Area_Squared(const TV& x1,const TV& x2,const TV& x3) // always positive for clockwise vertices: x1, x2, x3 
    {return (T).25*TV::Cross_Product(x2-x1,x3-x1).Magnitude_Squared();}

    T Size() const
    {return Area();}

    template<class T_ARRAY>
    static T Size(const T_ARRAY& X)
    {STATIC_ASSERT(T_ARRAY::m==3);return Area(X(0),X(1),X(2));}

    template<class T_ARRAY>
    static T Signed_Size(const T_ARRAY& X)
    {return Size(X);}

    T Aspect_Ratio() const
    {return Aspect_Ratio(X.x,X.y,X.z);}

    PLANE<T> Plane() const
    {return PLANE<T>(Normal(),X.x);}

    static T Aspect_Ratio(const TV& x1_input,const TV& x2_input,const TV& x3_input)
    {TV u=x1_input-x2_input,v=x2_input-x3_input,w=x3_input-x1_input;
    T u2=TV::Dot_Product(u,u),v2=TV::Dot_Product(v,v),w2=TV::Dot_Product(w,w);
    return max(u2,v2,w2)/sqrt(TV::Cross_Product(u,v).Magnitude_Squared());}

    static T Minimum_Edge_Length(const TV& x1,const TV& x2,const TV& x3)
    {return sqrt(min((x2-x1).Magnitude_Squared(),(x3-x1).Magnitude_Squared(),(x3-x2).Magnitude_Squared()));}

    static T Maximum_Edge_Length(const TV& x1,const TV& x2,const TV& x3)
    {return sqrt(max((x2-x1).Magnitude_Squared(),(x3-x1).Magnitude_Squared(),(x3-x2).Magnitude_Squared()));}

    static T Minimum_Altitude(const TV& x1,const TV& x2,const TV& x3)
    {return 2*Area(x1,x2,x3)/Maximum_Edge_Length(x1,x2,x3);}
    
    template<class T_ARRAY>
    static TV Barycentric_Coordinates(const TV& location,const T_ARRAY& X)
    {STATIC_ASSERT(T_ARRAY::m==3);return Barycentric_Coordinates(location,X(0),X(1),X(2));}

    template<class T_ARRAY>
    static TV Clamped_Barycentric_Coordinates(const TV& location,const T_ARRAY& X)
    {STATIC_ASSERT(T_ARRAY::m==3);return Clamped_Barycentric_Coordinates(location,X(0),X(1),X(2));}

    TV Sum_Barycentric_Coordinates(const TRIANGLE_3D<T>& embedded_triangle) const
    {return Barycentric_Coordinates(embedded_triangle.X.x)+Barycentric_Coordinates(embedded_triangle.X.y)+Barycentric_Coordinates(embedded_triangle.X.z);}

    TV Barycentric_Coordinates(const TV& location) const 
    {return Barycentric_Coordinates(location,X.x,X.y,X.z);}

    static TV Point_From_Barycentric_Coordinates(const TV& weights,const TV& x1,const TV& x2,const TV& x3) // clockwise vertices
    {return weights.x*x1+weights.y*x2+weights.z*x3;}

    template<class T_ARRAY>
    static TV Point_From_Barycentric_Coordinates(const TV& weights,const T_ARRAY& X) // clockwise vertices
    {STATIC_ASSERT(T_ARRAY::m==3);return weights.x*X(0)+weights.y*X(1)+weights.z*X(2);}

    TV Point_From_Barycentric_Coordinates(const TV& weights) const
    {return Point_From_Barycentric_Coordinates(weights,X.x,X.y,X.z);}

    template<class T_ARRAY>
    static TV Normal(const T_ARRAY& X)
    {return PLANE<T>::Normal(X);}

    TV Normal() const
    {return Normal(X);}

    TV Raw_Normal() const
    {return TV::Cross_Product(X.y-X.x,X.z-X.x);}

    static TV Center(const TV& x1,const TV& x2,const TV& x3) // centroid
    {return (T)one_third*(x1+x2+x3);}

    TV Center() const // centroid
    {return Center(X.x,X.y,X.z);}

    TV Incenter() const // intersection of angle bisectors
    {TV edge_lengths((X.z-X.y).Magnitude(),(X.x-X.z).Magnitude(),(X.y-X.x).Magnitude());T perimeter=edge_lengths.x+edge_lengths.y+edge_lengths.z;assert(perimeter>0);
    return Point_From_Barycentric_Coordinates(edge_lengths/perimeter);}

    bool Point_Face_Interaction(const TV& x,const TV& v,const INDIRECT_ARRAY<ARRAY_VIEW<TV>,VECTOR<int,3>&> V_face,const T interaction_distance,T& distance,
        TV& interaction_normal,VECTOR<T,TV::m+1>& weights,const bool allow_negative_weights,const bool exit_early) const
    {return Point_Face_Interaction(x,v,V_face(0),V_face(1),V_face(2),interaction_distance,distance,interaction_normal,weights,allow_negative_weights,exit_early);}

    bool Point_Face_Collision(const TV& x,const TV& v,const INDIRECT_ARRAY<ARRAY_VIEW<TV>,VECTOR<int,3>&> V_face,const T dt,const T collision_thickness,
        T& collision_time,TV& normal,VECTOR<T,TV::m+1>& weights,const bool exit_early) const
    {return Point_Face_Collision(x,v,V_face(0),V_face(1),V_face(2),dt,collision_thickness,collision_time,normal,weights,exit_early);}

    RANGE<TV> Bounding_Box() const
    {return RANGE<TV>::Bounding_Box(X.x,X.y,X.z);}

    T Signed_Distance(const TV& location) const
    {return Normal().Dot(location-X.x);}

    bool Lazy_Inside_Plane(const TV& location) const
    {return Raw_Normal().Dot(location-X.x)<=0;}  

    bool Lazy_Outside_Plane(const TV& location) const
    {return !Lazy_Inside_Plane(location);}

    bool Inside_Plane(const TV& location,const T thickness_over_two) const
    {return Signed_Distance(location)<=-thickness_over_two;}

    bool Outside_Plane(const TV& location,const T thickness_over_two) const
    {return !Inside_Plane(location,-thickness_over_two);}

//#####################################################################
    void Change_Size(const T delta);
    bool Inside(const TV& point,const T thickness_over_two=0) const;
    bool Point_Inside_Triangle(const TV& point,const T thickness_over_2=0) const;
    bool Planar_Point_Inside_Triangle(const TV& point,const T thickness_over_2=0) const;
    bool Lazy_Planar_Point_Inside_Triangle(const TV& point) const;
    T Minimum_Edge_Length() const;
    T Maximum_Edge_Length() const;
    int Region(const TV& location,int& region_id,const T tolerance) const;
    TV Closest_Point(const TV& location,TV& weights) const;
    T Distance_To_Triangle(const TV& location) const;
    T Minimum_Angle() const;
    T Maximum_Angle() const;
    T Signed_Solid_Angle(const TV& center) const;
    bool Point_Face_Interaction(const TV& x,const T interaction_distance,const bool allow_negative_weights,T& distance) const;
    void Point_Face_Interaction_Data(const TV& x,T& distance,TV& interaction_normal,VECTOR<T,TV::m+1>& weights,const bool perform_attractions) const;    
    bool Point_Face_Interaction(const TV& x,const TV& v,const TV& v1,const TV& v2,const TV& v3,const T interaction_distance,
        T& distance,TV& interaction_normal,VECTOR<T,TV::m+1>& weights,const bool allow_negative_weights,const bool exit_early) const;
    static POINT_SIMPLEX_COLLISION_TYPE Robust_Point_Face_Collision(const TRIANGLE_3D<T>& initial_triangle,const TRIANGLE_3D<T>& final_triangle,const TV& x,
        const TV& final_x,const T dt,const T collision_thickness,T& collision_time,TV& normal,VECTOR<T,TV::m+1>& weights);
    bool Point_Face_Collision(const TV& x,const TV& v,const TV& v1,const TV& v2,const TV& v3,const T dt,const T collision_thickness,
        T& collision_time,TV& normal,VECTOR<T,TV::m+1>& weights,const bool exit_early) const;
    void Clip_To_Box(const RANGE<TV>& box,ARRAY<TRIANGLE_3D<T> >& clipped_simplices) const;
    static void Cut_With_Hyperplane_And_Discard_Outside_Simplices(const TRIANGLE_3D<T>& triangle,const PLANE<T>& cutting_plane,ARRAY<TRIANGLE_3D<T> >& negative_triangles);
    static void Cut_With_Hyperplane(const TRIANGLE_3D<T>& triangle,const PLANE<T>& cutting_plane,ARRAY<TRIANGLE_3D<T> >& negative_triangles,
        ARRAY<TRIANGLE_3D<T> >& positive_triangles,T tol=0);
    struct INTERSECTS_HELPER
    {
        TV n,x,w;
        int pos,neg,is,i[2];
        T th[2],t[2];
    };
    bool Intersects(const TRIANGLE_3D<T>& triangle,T theta_tol=0,INTERSECTS_HELPER* ih=0) const;
    static TV Barycentric_Coordinates(const TV& location,const TV& x1,const TV& x2,const TV& x3); // clockwise vertices
    static TV Clamped_Barycentric_Coordinates(const TV& location,const TV& x1,const TV& x2,const TV& x3,const T tolerance=1e-7); // clockwise vertices
//#####################################################################
};

template<class T> std::ostream& operator<<(std::ostream& output,const TRIANGLE_3D<T>& triangle)
{output<<triangle.X.x<<", "<<triangle.X.y<<", "<<triangle.X.z;return output;}

}
#endif
