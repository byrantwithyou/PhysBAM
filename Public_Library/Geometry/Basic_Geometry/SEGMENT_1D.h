//#####################################################################
// Copyright 2007, Nipun Kwatra, Jonathan Su.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class SEGMENT_1D  
//##################################################################### 
#ifndef __SEGMENT_1D__
#define __SEGMENT_1D__

#include <Tools/Arrays/ARRAYS_FORWARD.h>
#include <Tools/Log/DEBUG_UTILITIES.h>
#include <Tools/Utilities/TYPE_UTILITIES.h>
#include <Tools/Vectors/VECTOR_2D.h>
namespace PhysBAM{

template<class TV> class RANGE;

template<class T>
class SEGMENT_1D
{
    typedef VECTOR<T,1> TV;
public:
    VECTOR<TV,2> X;

    SEGMENT_1D()
        :X(TV(),TV(1))
    {}

    SEGMENT_1D(const VECTOR<T,1>& x1_input,const VECTOR<T,1>& x2_input)
        :X(x1_input,x2_input)
    {}

    template<class T_ARRAY>
    explicit SEGMENT_1D(const T_ARRAY& X_input)
        :X(X_input)
    {
        STATIC_ASSERT(T_ARRAY::m==2);
    }

    T Length() const
    {return (X.y-X.x).Magnitude();}

    T Size() const
    {return Length();}

    VECTOR<T,2> Barycentric_Coordinates(const TV& location) const
    {return Barycentric_Coordinates(location,X.x,X.y);}

    VECTOR<T,2> Sum_Barycentric_Coordinates(const SEGMENT_1D<T>& embedded_segment) const
    {return Barycentric_Coordinates(embedded_segment.X.x)+Barycentric_Coordinates(embedded_segment.X.y);}

    static VECTOR<T,2> Barycentric_Coordinates(const TV& location,const TV& x0,const TV& x1)
    {T fraction=(location.x-x0.x)/(x1-x0).Magnitude();return VECTOR<T,2>((T)1-fraction,fraction);}

    template<class T_ARRAY>
    static VECTOR<T,2> Barycentric_Coordinates(const TV& location,const T_ARRAY& X)
    {STATIC_ASSERT(T_ARRAY::m==2);return Barycentric_Coordinates(location,X(0),X(1));}

    static TV Point_From_Barycentric_Coordinates(const VECTOR<T,2>& weights,const TV& x0,const TV& x1)
    {return weights.x*x0+weights.y*x1;}

    static TV Point_From_Barycentric_Coordinates(const TV& weights,const TV& x0,const TV& x1)
    {return weights.x*x0+(1-weights.x)*x1;}

    TV Point_From_Barycentric_Coordinates(const VECTOR<T,2>& weights) const
    {return Point_From_Barycentric_Coordinates(weights,X.x,X.y);}

    template<class T_ARRAY>
    static TV Point_From_Barycentric_Coordinates(const VECTOR<T,2>& weights,const T_ARRAY& X)
    {STATIC_ASSERT(T_ARRAY::m==2);return Point_From_Barycentric_Coordinates(weights,X(0),X(1));}

    VECTOR<T,1> Center() const
    {return (T).5*(X.x+X.y);}

    RANGE<TV> Bounding_Box() const
    {return RANGE<TV>::Bounding_Box(X.x,X.y);}

    bool Inside(const TV& location,const T thickness_over_2=0)
    {return X.x.x<=location.x && location.x<=X.y.x;}

//#####################################################################
    static T Negative_Material(const ARRAY<TV>& X,const ARRAY<T>& phis,const VECTOR<int,2>& indices){PHYSBAM_NOT_IMPLEMENTED();}
    void Clip_To_Box(const RANGE<TV>& box,ARRAY<SEGMENT_1D<T> >& clipped_simplices) const {PHYSBAM_NOT_IMPLEMENTED();}
    VECTOR<T,1> Closest_Point_On_Segment(const VECTOR<T,1>& point) const {PHYSBAM_NOT_IMPLEMENTED();}
//#####################################################################
};

template<class T> std::ostream &operator<<(std::ostream &output,const SEGMENT_1D<T> &segment)
{output << segment.X.x << ", " << segment.X.y;return output;}

}
#endif

