//#####################################################################
// Copyright 2007, Geoffrey Irving.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __PYTHON_UTILITIES__
#define __PYTHON_UTILITIES__

#include <Tools/Matrices/MATRIX_FORWARD.h>
#include <Tools/Vectors/VECTOR_FORWARD.h>
#include <boost/format.hpp>
#include <boost/python.hpp>
#include <PhysBAM_Geometry/Geometry/GEOMETRY_FORWARD.h>
namespace PhysBAM{

template<class T> class CYLINDER;
template<class T> class RING;
template<class TV> class GRID;

//#####################################################################
// Function Dereference_If_Pointer
//#####################################################################
template<class T> inline T& Dereference_If_Pointer(T& x)
{return x;}

template<class T> inline T& Dereference_If_Pointer(T* x)
{return *x;}

//#####################################################################
// Function Repr
//#####################################################################
template<class T>
std::string Class(const T& x)
{
    return boost::python::extract<std::string>(boost::python::object(x).attr("__class__").attr("__name__"));
}

// by default, defer to python __repr__
template<class T> std::string Repr_Helper(const T& x)
{return boost::python::extract<std::string>(boost::python::object(x).attr("__repr__")());}

inline std::string Repr_Helper(const float x)
{return (boost::format("%.9g")%x).str();}

inline std::string Repr_Helper(const double x)
{return (boost::format("%.17g")%x).str();}

inline std::string Repr_Helper(const long double x)
{return (boost::format("%.21lg")%x).str();}

template<class T,class TV> std::string Repr_Vector(const TV& x)
{std::ostringstream stream;stream<<Class(x)<<'(';
int d=sizeof(x)/sizeof(T);const T* data=(const T*)&x;
for(int i=0;i<d-1;i++) stream<<Repr(data[i])<<',';stream<<Repr(data[d-1])<<')';
return stream.str();}

template<class TV> std::string Repr_Grid(const GRID<TV>& grid)
{std::ostringstream stream;stream<<Class(grid)<<'(';
for(int i=0;i<TV::m;i++) stream<<Repr(grid.Counts()[i])<<',';
stream<<Repr(grid.Domain())<<(grid.Is_MAC_Grid()?",True)":",False)");
return stream.str();}

template<class T,int d> std::string Repr_Helper(const VECTOR<T,d>& x){return Repr_Vector<T>(x);}
template<class T,int m,int n> std::string Repr_Helper(const MATRIX<T,m,n>& x){return Repr_Vector<T>(x);}
template<class T,int d> std::string Repr_Helper(const SYMMETRIC_MATRIX<T,d>& x){return Repr_Vector<T>(x);}
template<class T,int d> std::string Repr_Helper(const DIAGONAL_MATRIX<T,d>& x){return Repr_Vector<T>(x);}
template<class T,int d> std::string Repr_Helper(const UPPER_TRIANGULAR_MATRIX<T,d>& x){return Repr_Vector<T>(x);}

template<class TV> std::string Repr_Helper(const RANGE<TV>& box)
{return (boost::format("%s(%s,%s)")%Class(box)%Repr(box.min_corner)%Repr(box.max_corner)).str();}

template<class TV> std::string Repr_Helper(const BOX<TV>& box)
{return (boost::format("%s(%s,%s)")%Class(box)%Repr(box.min_corner)%Repr(box.max_corner)).str();}

template<class TV> std::string Repr_Helper(const SPHERE<TV>& sphere)
{return (boost::format("%s(%s,%s)")%Class(sphere)%Repr(sphere.center)%Repr(sphere.radius)).str();}

template<class T> std::string Repr_Helper(const PLANE<T>& plane)
{return (boost::format("%s(%s,%s)")%Class(plane)%Repr(plane.normal)%Repr(plane.x1)).str();}

template<class T> std::string Repr_Helper(const CYLINDER<T>& cylinder)
{return (boost::format("%s(%s,%s,%s)")%Class(cylinder)%Repr(cylinder.plane1.x1)%Repr(cylinder.plane2.x1)%Repr(cylinder.radius)).str();}

template<class T> std::string Repr_Helper(const RING<T>& ring)
{return (boost::format("%s(%s,%s,%s,%s)")%Class(ring)%Repr(ring.plane1.x1)%Repr(ring.plane2.x1)%Repr(ring.outer_radius)%Repr(ring.inner_radius)).str();}

template<class T> std::string Repr_Helper(const ROTATION<VECTOR<T,3> >& rotation)
{T angle;VECTOR<T,3> axis;rotation.Get_Angle_Axis(angle,axis); // not exact, unfortunately
return (boost::format("%s(%s,%s)")%Class(rotation)%Repr(angle)%Repr(axis)).str();}

template<class TV> std::string Repr_Helper(const FRAME<TV>& frame)
{return (boost::format("%s(%s,%s)")%Class(frame)%Repr(frame.t)%Repr(frame.r)).str();}

template<class TV> std::string Repr_Helper(const TWIST<TV>& twist)
{return (boost::format("%s(%s,%s)")%Class(twist)%Repr(twist.linear)%Repr(twist.angular)).str();}

template<class TV> std::string Repr_Helper(const GRID<TV>& grid){return Repr_Grid(grid);}

template<class T> inline std::string Repr(const T& x)
{return Repr_Helper(x);}

}
#endif
