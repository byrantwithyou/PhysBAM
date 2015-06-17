//#####################################################################
// Copyright 2007-2009, Jon Gretarsson, Geoffrey Irving, Nipun Kwatra, Andrew Selle, Jonathan Su.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class POINT_SIMPLEX_1D  
//##################################################################### 
#ifndef __POINT_SIMPLEX_1D__
#define __POINT_SIMPLEX_1D__

#include <Tools/Arrays/ARRAYS_FORWARD.h>
#include <Tools/Log/DEBUG_UTILITIES.h>
#include <Tools/Math_Tools/ONE.h>
#include <Tools/Matrices/SYMMETRIC_MATRIX.h>
#include <Tools/Utilities/TYPE_UTILITIES.h>
#include <Tools/Vectors/VECTOR_1D.h>
#include <Geometry/Basic_Geometry/BASIC_GEOMETRY_FORWARD.h>
namespace PhysBAM{

template<class TV> class RANGE;

template<class T>
class POINT_SIMPLEX_1D
{
    typedef VECTOR<T,1> TV;
public:
    typedef int HAS_UNTYPED_READ_WRITE;
    typedef TV VECTOR_T;

    VECTOR<TV,1> X;
    bool direction; // false for -1 and true for 1 direction
 
    POINT_SIMPLEX_1D()
        :direction(true)
    {}

    POINT_SIMPLEX_1D(const TV& x0,const bool direction)
        :X(x0),direction(direction)
    {}

    template<class T_ARRAY>
    explicit POINT_SIMPLEX_1D(const T_ARRAY& X_input)
        :X(X_input(0)),direction(true)
    {
        STATIC_ASSERT(T_ARRAY::m==1); // TODO: initialize directions properly
    }

    const TV& Center() const
    {return X.x;}

    T Size() const
    {return (T)1;}

    static TV Normal(const TV& x0,const int direction)
    {return TV((T)(direction?1:-1));}

    TV Normal() const
    {return POINT_SIMPLEX_1D<T>::Normal(X.x,direction);}

    TV Normal(const TV& location) const
    {return Normal();}

    SYMMETRIC_MATRIX<T,1> Hessian(const TV& X) const
    {return SYMMETRIC_MATRIX<T,1>();}

    static ONE Clamped_Barycentric_Coordinates(const TV& location,const TV& x0)
    {return ONE();}

    template<class T_ARRAY>
    static typename enable_if<is_same<typename T_ARRAY::ELEMENT,TV>::value && T_ARRAY::m==1,ONE>::type
    Clamped_Barycentric_Coordinates(const TV& location,const T_ARRAY& X)
    {return ONE();}

    TV Sum_Barycentric_Coordinates(const POINT_SIMPLEX_1D<T>& embedded_point_simplex) const
    {return TV::All_Ones_Vector();} // TODO

    ONE Barycentric_Coordinates(const TV& location) const
    {return ONE();}

    static TV Point_From_Barycentric_Coordinates(const ONE,const TV& x0)
    {return x0;}

    TV Point_From_Barycentric_Coordinates(const TV& weights)
    {return X.x;}

    template<class T_ARRAY>
    static typename enable_if<is_same<typename T_ARRAY::ELEMENT,TV>::value && T_ARRAY::m==1,TV>::type
    Point_From_Barycentric_Coordinates(const ONE,const T_ARRAY& X)
    {return X(1);}

    TV Point_From_Barycentric_Coordinates(const ONE) const
    {return X.x;}

    RANGE<TV> Bounding_Box() const
    {return RANGE<TV>(X.x);}

    T Signed_Distance(const TV& location) const
    {return (direction?(T)1:(T)-1)*(location.x-X.x.x);}

    void Clip_To_Box(const RANGE<TV>& box,ARRAY<POINT_SIMPLEX_1D<T> >& clipped_simplices) const
    {clipped_simplices.Remove_All();
    if(box.Lazy_Inside(X.x)) clipped_simplices.Append(*this);}

    bool Inside(const TV& point,const T thickness_over_two=0) const
    {if(X.x.x-thickness_over_two <= point.x && X.x.x+thickness_over_two >= point.x) return true;
    else return false;}

    VECTOR<T,0> Principal_Curvatures(const TV& X) const
    {return VECTOR<T,0>();}

    static std::string Name()
    {return "POINT_SIMPLEX_1D<T>";}

    template<class RW> void Read(std::istream& input)
    {Read_Binary<RW>(input,direction,X.x);}

    template<class RW> void Write(std::ostream& output) const
    {Write_Binary<RW>(output,direction,X.x);}

    static bool Point_Face_Collision(const POINT_SIMPLEX_1D<T>& initial_simplex,const VECTOR<T,1>& x,const VECTOR<T,1>& v,const VECTOR<T,1>& v1,const T dt,const T collision_thickness,
        T& collision_time,VECTOR<T,1>& normal,VECTOR<T,2>& weights,const bool exit_early)
    {
        weights.x=-1;
        weights.y=1;
        T distance=x.x-initial_simplex.X.x.x,relative_speed=-weights.x*v.x-weights.y*v1.x;
        if(distance*relative_speed>=0 || abs(distance)>abs(dt*relative_speed)) return false;
        collision_time=distance/relative_speed;
        if(distance>0) normal=VECTOR<T,1>(1);
        else normal=VECTOR<T,1>(-1);
        return true;
    }

    static POINT_SIMPLEX_COLLISION_TYPE Robust_Point_Face_Collision(const POINT_SIMPLEX_1D<T>& initial_simplex,const POINT_SIMPLEX_1D<T>& final_simplex,
        const VECTOR<T,1>& x,const VECTOR<T,1>& final_x,const T dt,const T collision_thickness,T& collision_time,VECTOR<T,1>& normal,VECTOR<T,2>& weights)
    {
        if(final_simplex.Bounding_Box().Thickened(collision_thickness).Lazy_Inside(final_x)){
            collision_time=dt;
            return POINT_SIMPLEX_COLLISION_ENDS_INSIDE;}
        if(initial_simplex.Bounding_Box().Thickened(collision_thickness).Lazy_Inside(x)){
            collision_time=0;
            return POINT_SIMPLEX_COLLISION_ENDS_OUTSIDE;}
        VECTOR<T,1> v1=(final_simplex.X.x-initial_simplex.X.x)/dt,v=(final_x-x)/dt;
        if(Point_Face_Collision(initial_simplex,x,v,v1,dt,collision_thickness,collision_time,normal,weights,false))
            return POINT_SIMPLEX_COLLISION_ENDS_OUTSIDE;
        return POINT_SIMPLEX_NO_COLLISION;
    }
};
}
#endif
