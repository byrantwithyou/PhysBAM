//#####################################################################
// Copyright 2003-2006, Ronald Fedkiw, Geoffrey Irving, Neil Molino, Duc Nguyen, Avi Robinson-Mosher, Andrew Selle, Joseph Teran.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class TRIANGLE_2D
//#####################################################################
#ifndef __TRIANGLE_2D__
#define __TRIANGLE_2D__

#include <Tools/Math_Tools/constants.h>
#include <Tools/Math_Tools/RANGE.h>
#include <Tools/Matrices/MATRIX_2X2.h>
#include <Tools/Vectors/VECTOR.h>
#include <Geometry/Basic_Geometry/BASIC_GEOMETRY_FORWARD.h>
namespace PhysBAM{

template<class T>
class TRIANGLE_2D
{
    typedef VECTOR<T,2> TV;
public:
    VECTOR<TV,3> X;

    TRIANGLE_2D()
    {
        Specify_Three_Points(TV(0,0),TV(0,1),TV(1,0));
    }

    TRIANGLE_2D(const TV& x1_input,const TV& x2_input,const TV& x3_input)
    {
        Specify_Three_Points(x1_input,x2_input,x3_input);
    }

    template<class T_ARRAY>
    TRIANGLE_2D(const T_ARRAY& X_input)
        :X(X_input)
    {
        STATIC_ASSERT(T_ARRAY::m==3);
    }

    void Specify_Three_Points(const TV& x1_input,const TV& x2_input,const TV& x3_input)
    {X[0]=x1_input;X[1]=x2_input;X[2]=x3_input;}

    static T Signed_Area(const TV& x0,const TV& x1,const TV& x2)
    {return (T).5*TV::Cross_Product(x1-x0,x2-x0).x;}

    T Signed_Area() const
    {return Signed_Area(X[0],X[1],X[2]);}

    static T Area(const TV& x0,const TV& x1,const TV& x2)
    {return abs(Signed_Area(x0,x1,x2));}

    T Area() const
    {return abs(Signed_Area());}

    T Size() const
    {return Area();}

    T Signed_Size() const
    {return Signed_Area();}

    template<class T_ARRAY>
    static T Signed_Size(const T_ARRAY& X)
    {STATIC_ASSERT(T_ARRAY::m==3);return Signed_Area(X(0),X(1),X(2));}

    template<class T_ARRAY>
    static T Size(const T_ARRAY& X)
    {return abs(Signed_Size(X));}

    static bool Check_Orientation(const TV& x0,const TV& x1,const TV& x2)
    {return Signed_Area(x0,x1,x2)>=0;}

    bool Check_Orientation() const
    {return Signed_Area()>=0;}

    bool Fix_Orientation()
    {if(Check_Orientation()) return false;exchange(X[1],X[2]);return true;}

    template<class T_ARRAY>
    static T Half_Boundary_Measure(const T_ARRAY& X)
    {STATIC_ASSERT(T_ARRAY::m==3);return (T).5*((X(0)-X(1)).Magnitude()+(X(1)-X(2)).Magnitude()+(X(2)-X(0)).Magnitude());}

    T Aspect_Ratio() const
    {return Aspect_Ratio(X[0],X[1],X[2]);}

    static T Aspect_Ratio(const TV& x1_input,const TV& x2_input,const TV& x3_input)
    {TV u=x1_input-x2_input,v=x2_input-x3_input,w=x3_input-x1_input;
    T u2=TV::Dot_Product(u,u),v2=TV::Dot_Product(v,v),w2=TV::Dot_Product(w,w);
    return max(u2,v2,w2)/abs(TV::Cross_Product(u,v).x);}

    static T Minimum_Edge_Length(const TV& x0,const TV& x1,const TV& x2)
    {return sqrt(Minimum_Edge_Length_Squared(x0,x1,x2));}

    static T Minimum_Edge_Length_Squared(const TV& x0,const TV& x1,const TV& x2)
    {return min((x1-x0).Magnitude_Squared(),(x2-x0).Magnitude_Squared(),(x2-x1).Magnitude_Squared());}

    static T Maximum_Edge_Length(const TV& x0,const TV& x1,const TV& x2)
    {return sqrt(Maximum_Edge_Length_Squared(x0,x1,x2));}

    static T Maximum_Edge_Length_Squared(const TV& x0,const TV& x1,const TV& x2)
    {return max((x1-x0).Magnitude_Squared(),(x2-x0).Magnitude_Squared(),(x2-x1).Magnitude_Squared());}

    T Minimum_Altitude() const
    {return Minimum_Altitude(X[0],X[1],X[2]);}

    static T Minimum_Altitude(const TV& x0,const TV& x1,const TV& x2)
    {return 2*Area(x0,x1,x2)/Maximum_Edge_Length(x0,x1,x2);}

    static TV First_Two_Barycentric_Coordinates(const TV& location,const TV& x0,const TV& x1,const TV& x2)
    {return MATRIX<T,2>(x0-x2,x1-x2).Robust_Solve_Linear_System(location-x2);}

    VECTOR<T,3> Barycentric_Coordinates(const TV& location) const
    {return Barycentric_Coordinates(location,X[0],X[1],X[2]);}

    static VECTOR<T,3> Barycentric_Coordinates(const TV& location,const TV& x0,const TV& x1,const TV& x2)
    {TV w=First_Two_Barycentric_Coordinates(location,x0,x1,x2);return VECTOR<T,3>(w.x,w.y,1-w.x-w.y);}

    template<class T_ARRAY>
    static VECTOR<T,3> Barycentric_Coordinates(const TV& location,const T_ARRAY& X)
    {STATIC_ASSERT(T_ARRAY::m==3);return Barycentric_Coordinates(location,X(0),X(1),X(2));}

    static TV Point_From_Barycentric_Coordinates(const VECTOR<T,3>& weights,const TV& x0,const TV& x1,const TV& x2)
    {return weights.x*x0+weights.y*x1+weights.z*x2;}

    static TV Point_From_Barycentric_Coordinates(const TV& weights,const TV& x0,const TV& x1,const TV& x2)
    {return weights.x*x0+weights.y*x1+(1-weights.x-weights.y)*x2;}

    TV Point_From_Barycentric_Coordinates(const VECTOR<T,3>& weights) const
    {return Point_From_Barycentric_Coordinates(weights,X[0],X[1],X[2]);}

    template<class T_ARRAY>
    static TV Point_From_Barycentric_Coordinates(const VECTOR<T,3>& weights,const T_ARRAY& X)
    {STATIC_ASSERT(T_ARRAY::m==3);return Point_From_Barycentric_Coordinates(weights,X(0),X(1),X(2));}

    VECTOR<T,3> Sum_Barycentric_Coordinates(const TRIANGLE_2D<T>& embedded_triangle) const
    {return Barycentric_Coordinates(embedded_triangle.X[0],X)+Barycentric_Coordinates(embedded_triangle.X[1],X)+Barycentric_Coordinates(embedded_triangle.X[2],X);}

    static TV Center(const TV& x0,const TV& x1,const TV& x2) // centroid
    {return (T)one_third*(x0+x1+x2);}

    TV Center() const // centroid
    {return Center(X[0],X[1],X[2]);}

    TV Incenter() const // intersection of angle bisectors
    {VECTOR<T,3> edge_lengths((X[2]-X[1]).Magnitude(),(X[0]-X[2]).Magnitude(),(X[1]-X[0]).Magnitude());T perimeter=edge_lengths.x+edge_lengths.y+edge_lengths.z;assert(perimeter>0);
    return Point_From_Barycentric_Coordinates(edge_lengths/perimeter);}

    static TV Circumcenter(const TV& x0,const TV& x1,const TV& x2)
    {TV x1x2=x1-x0,x1x3=x2-x0,m1=(T).5*(x0+x1),m2=(T).5*(x0+x2),m1m2=m2-m1,x1x2_perp(-x1x2.y,x1x2.x);
    return m1+x1x2_perp*(TV::Dot_Product(m1m2,x1x3)/TV::Dot_Product(x1x2_perp,x1x3));}

    static VECTOR<T,3> Circumcenter_Barycentric_Coordinates(const TV& x0,const TV& x1,const TV& x2)
    {TV a=x2-x1,b=x2-x0,c=x1-x0;T aa=a.Magnitude_Squared(),bb=b.Magnitude_Squared(),cc=c.Magnitude_Squared();
    VECTOR<T,3> weights(aa*(bb+cc-aa),bb*(cc+aa-bb),cc*(aa+bb-cc));return weights/(weights.x+weights.y+weights.z);}

    T Minimum_Angle() const
    {TV s1=(X[0]-X[1]).Normalized(),s2=(X[1]-X[2]).Normalized(),s3=(X[2]-X[0]).Normalized();
    return acos(max(TV::Dot_Product(s1,-s2),TV::Dot_Product(-s1,s3),TV::Dot_Product(s2,-s3)));}

    T Maximum_Angle() const
    {TV s1=(X[0]-X[1]).Normalized(),s2=(X[1]-X[2]).Normalized(),s3=(X[2]-X[0]).Normalized();
    return acos(min(TV::Dot_Product(s1,-s2),TV::Dot_Product(-s1,s3),TV::Dot_Product(s2,-s3)));}

    bool Outside(const TV& location,const T thickness_over_2=0) const
    {return Outside(location,X[0],X[1],X[2],thickness_over_2);}

    static bool Outside(const TV& location,const TV& x0,const TV& x1,const TV& x2,const T thickness_over_2=0)
    {assert(Check_Orientation(x0,x1,x2));TV location_minus_x1=location-x0;
    TV edge1=x1-x0;if(TV::Cross_Product(location_minus_x1,edge1).x>thickness_over_2*edge1.Magnitude()) return true;
    TV edge3=x0-x2;if(TV::Cross_Product(location-x2,edge3).x>thickness_over_2*edge3.Magnitude()) return true;
    TV edge2=x2-x1;if(TV::Cross_Product(location-x1,edge2).x>thickness_over_2*edge2.Magnitude()) return true;
    return false;}

    bool Inside(const TV& location,const T thickness_over_2=0) const
    {return !Outside(location,thickness_over_2);}

    RANGE<TV> Bounding_Box() const
    {return RANGE<TV>::Bounding_Box(X[0],X[1],X[2]);}

    static bool Check_Delaunay_Criterion(TV a,TV b,TV c,TV d)
    {assert(Check_Orientation(a,b,c) && Check_Orientation(d,c,b));b-=a;c-=a;d-=a;
    return VECTOR<T,3>::Triple_Product(b.Append(b.Magnitude_Squared()),c.Append(c.Magnitude_Squared()),d.Append(d.Magnitude_Squared()))>=0;}

//#####################################################################
    bool Intersects(const TRIANGLE_2D& tri) const;
    static T Negative_Material(const ARRAY<TV>& X,const ARRAY<T>& phis,const VECTOR<int,3>& indices);
    static void Cut_With_Hyperplane(ARRAY<TV>& X,const LINE_2D<T>& cutting_plane,const VECTOR<int,3>& indices,ARRAY<VECTOR<int,3> >& left_tris,
        ARRAY<VECTOR<int,3> >& right_tris);
    void Clip_To_Box(const RANGE<TV>& box,ARRAY<TRIANGLE_2D<T> >& clipped_simplices) const;
    void Cut_With_Hyperplane_And_Discard_Outside_Simplices(const TRIANGLE_2D<T>& triangle,const LINE_2D<T>& cutting_plane,ARRAY<TRIANGLE_2D<T> >& negative_triangles) const;
private:
    static void Cut_Simplex(ARRAY<TV>& X,const VECTOR<int,3>& indices,const VECTOR<TV,3>& X_nodes,const VECTOR<T,3>& phi_nodes,
        ARRAY<VECTOR<int,3> >& left_tris,ARRAY<VECTOR<int,3> >& right_tris);
//#####################################################################
};

template<class T> std::ostream& operator<<(std::ostream& output,const TRIANGLE_2D<T>& triangle)
{return output<<triangle.X;}

}
#endif

