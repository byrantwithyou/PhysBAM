//#####################################################################
// Copyright 2012.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class HIDDEN_SURFACE_PRIMITIVES
//#####################################################################
#ifndef __HIDDEN_SURFACE_PRIMITIVES__
#define __HIDDEN_SURFACE_PRIMITIVES__

#ifdef USE_BOOST_GEOMETRY
#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Data_Structures/DATA_STRUCTURES_FORWARD.h>
#include <PhysBAM_Geometry/Basic_Geometry/PLANE.h>
#include <PhysBAM_Geometry/Basic_Geometry/TRIANGLE_3D.h>
#include <boost/geometry.hpp>
#include <boost/geometry/geometries/polygon.hpp>
#include <boost/geometry/geometries/register/point.hpp>
#include <boost/geometry/geometries/register/ring.hpp>
#include <boost/geometry/multi/geometries/multi_polygon.hpp>
#include <boost/geometry/multi/geometries/register/multi_point.hpp>
#include <boost/geometry/multi/geometries/register/multi_polygon.hpp>

typedef PhysBAM::VECTOR<float,2> VECTOR_float_2;
typedef PhysBAM::VECTOR<double,2> VECTOR_double_2;
typedef PhysBAM::VECTOR<float,3> VECTOR_float_3;
typedef PhysBAM::VECTOR<double,3> VECTOR_double_3;
BOOST_GEOMETRY_REGISTER_POINT_2D(VECTOR_double_2, double, cs::cartesian, x, y);
BOOST_GEOMETRY_REGISTER_POINT_2D(VECTOR_float_2, float, cs::cartesian, x, y);
BOOST_GEOMETRY_REGISTER_POINT_3D(VECTOR_double_3, double, cs::cartesian, x, y, z);
BOOST_GEOMETRY_REGISTER_POINT_3D(VECTOR_float_3, float, cs::cartesian, x, y, z);
template<class T,class A=void> struct ALLOCATOR_ARRAY:public PhysBAM::ARRAY<T>
{
    T& operator[](const int i)
    {return (*this)(i);}

    const T& operator[](const int i) const
    {return (*this)(i);}

    size_t size() const
    {return this->Size();}

    void push_back(const T& element)
    {this->Append(element);}

    void resize(size_t size)
    {this->Resize(size);}

    void clear()
    {this->Clean_Memory();}

    const T& back() const
    {return this->Last();}

    T& back()
    {return this->Last();}

    typedef const T* const_iterator;
    typedef T* iterator;
    typedef const T& const_reference;
};
typedef boost::geometry::model::polygon<VECTOR_float_2,true,true,ALLOCATOR_ARRAY,ALLOCATOR_ARRAY> POLYGON_float;
typedef boost::geometry::model::polygon<VECTOR_double_2,true,true,ALLOCATOR_ARRAY,ALLOCATOR_ARRAY> POLYGON_double;
BOOST_GEOMETRY_REGISTER_RING_TEMPLATED(ALLOCATOR_ARRAY);
BOOST_GEOMETRY_REGISTER_MULTI_POLYGON(ALLOCATOR_ARRAY<POLYGON_float>);
BOOST_GEOMETRY_REGISTER_MULTI_POLYGON(ALLOCATOR_ARRAY<POLYGON_double>);

namespace PhysBAM{

template<class T>
class SURFACE_PRIMITIVE
{
    typedef VECTOR<T,3> TV;
    typedef VECTOR<T,2> TV2;

public:
    typedef boost::geometry::model::polygon<VECTOR<T,2>,true,true,ALLOCATOR_ARRAY,ALLOCATOR_ARRAY> POLYGON;
    typedef boost::geometry::model::multi_polygon<POLYGON, ALLOCATOR_ARRAY> MULTI_POLYGON;

    MULTI_POLYGON projection;
    TRIANGLE_3D<T> triangle;
    int parent;

    SURFACE_PRIMITIVE(): parent(-1) {}
    SURFACE_PRIMITIVE(const TRIANGLE_3D<T> &t,const int pa);
    void Init_Projection();
};

template<class T>
class HIDDEN_SURFACE_PRIMITIVES
{
    typedef VECTOR<T,3> TV;
    typedef VECTOR<T,2> TV2;

public:

    ARRAY<SURFACE_PRIMITIVE<T> > primitives;
    ARRAY<int> order;

    int Add(const TRIANGLE_3D<T> &t,int pa=-1);
    void Initialize(DIRECTED_GRAPH<>& dg);
    int Divide_Primitive(int divide,int cutter); // return: divide=no inside, -1=no outside, else both
    bool Test_Edge(int a,int b) {return Projections_Intersect(a,b);} // 0=no edge, 1=a->b
    void Emit_Node(int a) {order.Append(a);}
    bool Projections_Intersect(int a,int b);
    void Handle_Intersection(int a,int b,ARRAY<ARRAY<int> >& adjacency_list,
        ARRAY<VECTOR<int,2> >& pairs,HASHTABLE<VECTOR<int,2> >& edges);
    void Add_Edge(HASHTABLE<VECTOR<int,2> >& edges,int a,int b);
    void Add_Edge(HASHTABLE<VECTOR<int,2> >& edges,int a,int b,TV2 pt);
};
}
#endif
#endif
