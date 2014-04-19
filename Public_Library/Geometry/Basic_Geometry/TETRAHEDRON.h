//#####################################################################
// Copyright 2002-2006, Robert Bridson, Ronald Fedkiw, Eran Guendelman, Geoffrey Irving, Neil Molino, Avi Robinson-Mosher, Andrew Selle, Eftychios Sifakis, Joseph Teran.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class TETRAHEDRON
//#####################################################################
#ifndef __TETRAHEDRON__
#define __TETRAHEDRON__

#include <Tools/Matrices/MATRIX_3X3.h>
#include <Tools/Vectors/VECTOR.h>
#include <Geometry/Basic_Geometry/TRIANGLE_3D.h>
namespace PhysBAM{

template<class T>
class TETRAHEDRON
{
    typedef VECTOR<T,3> TV;
public:
    VECTOR<TV,4> X;
    VECTOR<TRIANGLE_3D<T>,4> triangle; // includes the planes and the normals

    TETRAHEDRON();

    TETRAHEDRON(const TV& x1_input,const TV& x2_input,const TV& x3_input,const TV& x4_input);

    template<class T_ARRAY>
    TETRAHEDRON(const T_ARRAY& X_input)
        :X(X_input)
    {
        STATIC_ASSERT(T_ARRAY::m==4);
        Create_Triangles();
    }

    static T Volume(const TV& x0,const TV& x1,const TV& x2,const TV& x3)
    {return abs(Signed_Volume(x0,x1,x2,x3));}

    static T Signed_Volume(const TV& x0,const TV& x1,const TV& x2,const TV& x3)
    {return (T)1/6*TV::Triple_Product(x1-x0,x2-x0,x3-x0);} 

    TV First_Three_Barycentric_Coordinates(const TV& location) const
    {return First_Three_Barycentric_Coordinates(location,X[0],X[1],X[2],X[3]);}

    VECTOR<T,4> Barycentric_Coordinates(const TV& location) const
    {return Barycentric_Coordinates(location,X[0],X[1],X[2],X[3]);}

    static TV First_Three_Barycentric_Coordinates(const TV& location,const TV& x0,const TV& x1,const TV& x2,const TV& x3)
    {return MATRIX<T,3>(x0-x3,x1-x3,x2-x3).Robust_Inverse_Times(location-x3);}

    static VECTOR<T,4> Barycentric_Coordinates(const TV& location,const TV& x0,const TV& x1,const TV& x2,const TV& x3)
    {TV w=First_Three_Barycentric_Coordinates(location,x0,x1,x2,x3);return VECTOR<T,4>(w.x,w.y,w.z,1-w.x-w.y-w.z);}

    template<class T_ARRAY>
    static TV First_Three_Barycentric_Coordinates(const TV& location,const T_ARRAY& X)
    {STATIC_ASSERT(T_ARRAY::m==4);return First_Three_Barycentric_Coordinates(location,X(0),X(1),X(2),X(3));}

    template<class T_ARRAY>
    static VECTOR<T,4> Barycentric_Coordinates(const TV& location,const T_ARRAY& X)
    {STATIC_ASSERT(T_ARRAY::m==4);return Barycentric_Coordinates(location,X(0),X(1),X(2),X(3));}

    TV Point_From_Barycentric_Coordinates(const TV& weights) const
    {return Point_From_Barycentric_Coordinates(weights,X[0],X[1],X[2],X[3]);}

    TV Point_From_Barycentric_Coordinates(const VECTOR<T,4>& weights) const
    {return Point_From_Barycentric_Coordinates(weights,X[0],X[1],X[2],X[3]);}

    static TV Point_From_Barycentric_Coordinates(const TV& weights,const TV& x0,const TV& x1,const TV& x2,const TV& x3)
    {return weights.x*x0+weights.y*x1+weights.z*x2+((T)1-weights.x-weights.y-weights.z)*x3;}

    static TV Point_From_Barycentric_Coordinates(const VECTOR<T,4>& weights,const TV& x0,const TV& x1,const TV& x2,const TV& x3)
    {return weights[0]*x0+weights[1]*x1+weights[2]*x2+weights[3]*x3;}

    template<class T_ARRAY>
    static TV Point_From_Barycentric_Coordinates(const VECTOR<T,4>& weights,const T_ARRAY& X)
    {STATIC_ASSERT(T_ARRAY::m==4);return Point_From_Barycentric_Coordinates(weights,X(0),X(1),X(2),X(3));}

    static bool Barycentric_Inside(const TV& location,const TV& x0,const TV& x1,const TV& x2,const TV& x3,const T tolerance=0)
    {TV w=First_Three_Barycentric_Coordinates(location,x0,x1,x2,x3);
    return w.x>=-tolerance && w.y>=-tolerance && w.z>=-tolerance && w.x+w.y+w.z<=1+tolerance;}

    VECTOR<T,4> Sum_Barycentric_Coordinates(const TETRAHEDRON<T>& embedded_tetrahedron) const
    {VECTOR<T,4> coordinates(Barycentric_Coordinates(embedded_tetrahedron.X[0])+Barycentric_Coordinates(embedded_tetrahedron.X[1])+Barycentric_Coordinates(embedded_tetrahedron.X[2])+Barycentric_Coordinates(embedded_tetrahedron.X[3]));
    return coordinates;}

    static TV Center(const TV& x0,const TV& x1,const TV& x2,const TV& x3) // centroid
    {return (T).25*(x0+x1+x2+x3);}

    TV Center() const // centroid
    {return Center(X[0],X[1],X[2],X[3]);}

    static TV Circumcenter(const TV& x0,const TV& x1,const TV& x2,const TV& x3)
    {TV u=x1-x0,v=x2-x0,w=x3-x0,cross_uv=TV::Cross_Product(u,v),cross_wu=TV::Cross_Product(w,u),
                               cross_vw=TV::Cross_Product(v,w);
    T determinant=TV::Dot_Product(cross_uv,w),uu=u.Magnitude_Squared(),vv=v.Magnitude_Squared(),ww=w.Magnitude_Squared();
    return x0+(ww*cross_uv+vv*cross_wu+uu*cross_vw)*((T).5/determinant);}

    static T Minimum_Edge_Length_Squared(const TV& x0,const TV& x1,const TV& x2,const TV& x3)
    {return min((x0-x1).Magnitude_Squared(),(x1-x2).Magnitude_Squared(),(x2-x0).Magnitude_Squared(),(x0-x3).Magnitude_Squared(),(x1-x3).Magnitude_Squared(),(x2-x3).Magnitude_Squared());}

    static T Minimum_Edge_Length(const TV& x0,const TV& x1,const TV& x2,const TV& x3)
    {return sqrt(Minimum_Edge_Length_Squared(x0,x1,x2,x3));}

    T Minimum_Edge_Length() const
    {return Minimum_Edge_Length(X[0],X[1],X[2],X[3]);}

    static T Maximum_Edge_Length_Squared(const TV& x0,const TV& x1,const TV& x2,const TV& x3)
    {return max((x0-x1).Magnitude_Squared(),(x1-x2).Magnitude_Squared(),(x2-x0).Magnitude_Squared(),(x0-x3).Magnitude_Squared(),(x1-x3).Magnitude_Squared(),(x2-x3).Magnitude_Squared());}

    static T Maximum_Edge_Length(const TV& x0,const TV& x1,const TV& x2,const TV& x3)
    {return sqrt(Maximum_Edge_Length_Squared(x0,x1,x2,x3));}

    T Maximum_Edge_Length() const
    {return Maximum_Edge_Length(X[0],X[1],X[2],X[3]);}

    T Size() const
    {return Volume();}

    T Signed_Size() const
    {return Signed_Volume();}

    template<class T_ARRAY>
    static T Signed_Size(const T_ARRAY& X)
    {STATIC_ASSERT(T_ARRAY::m==4);return Signed_Volume(X(0),X(1),X(2),X(3));}

    template<class T_ARRAY>
    static T Size(const T_ARRAY& X)
    {return abs(Signed_Size(X));}

    template<class T_ARRAY>
    static T Half_Boundary_Measure(const T_ARRAY& X)
    {STATIC_ASSERT(T_ARRAY::m==4);return (T).5*(TRIANGLE_3D<T>::Area(X(0),X(1),X(2))+TRIANGLE_3D<T>::Area(X(0),X(3),X(1))+TRIANGLE_3D<T>::Area(X(0),X(2),X(3))+TRIANGLE_3D<T>::Area(X(1),X(3),X(2)));}

//#####################################################################
    void Create_Triangles();
    TV Normal(const TV& location,const int aggregate=0) const;
    bool Inside(const TV& location,const T thickness_over_two=0) const;
    bool Outside(const TV& location,const T thickness_over_two=0) const;
    bool Boundary(const TV& location,const T thickness_over_two) const;
    TETRAHEDRON<T> Thickened(const T thickness_over_two) const;
    RANGE<TV> Bounding_Box() const;
    TV Surface(const TV& location) const;
    TV Closest_Point(const TV& location,TV& weights) const;
    T Volume() const;
    T Signed_Volume() const;
    T Minimum_Angle() const;
    T Maximum_Angle() const;
    T Minimum_Altitude() const;
    T Aspect_Ratio() const;  // largest_edge/smallest_altitude
    T Minimum_Dihedral_Angle() const;
    T Maximum_Dihedral_Angle() const;
    static T Signed_Reciprocal_Aspect_Ratio(const TV& x0,const TV& x1,const TV& x2,const TV& x3);
    static T Negative_Material(const ARRAY<TV>& X,const ARRAY<T>& phis,const VECTOR<int,4>& indices);
    static void Cut_With_Hyperplane(ARRAY<TV>& X,const PLANE<T>& cutting_plane,const VECTOR<int,4>& indices,ARRAY<VECTOR<int,4> >& left_simplices,
        ARRAY<VECTOR<int,4> >& right_simplices);
    void Clip_To_Box(const RANGE<TV>& box,ARRAY<TETRAHEDRON<T> >& clipped_simplices) const;
    void Cut_With_Hyperplane_And_Discard_Outside_Simplices(const TETRAHEDRON<T>& tetrahedron,const PLANE<T>& cutting_plane,ARRAY<TETRAHEDRON<T> >& negative_triangles) const;
private:
    static void Cut_Simplex(ARRAY<TV>& X,const VECTOR<int,4>& indices,const VECTOR<TV,4>& X_nodes,const VECTOR<T,4>& phi_nodes,
        ARRAY<VECTOR<int,4> >& left_simplices,ARRAY<VECTOR<int,4> >& right_simplices);
//#####################################################################
};
}
#endif
