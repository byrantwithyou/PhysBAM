//#####################################################################
// Copyright 2011, Jeffrey Hellrung
// This file is part of PhysBAM whose distribution is governed by the
// license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################

#ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_GEOMETRY_SUB_CUBE_POLYTOPE_HPP
#define PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_GEOMETRY_SUB_CUBE_POLYTOPE_HPP

#include <cassert>

#include <PhysBAM_Tools/Vectors/VECTOR.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TOPOLOGY_BASED_GEOMETRY_FORWARD.h>
#include <Jeffrey_Utilities/BOUNDED_LIST.h>
#include <Jeffrey_Utilities/Geometry/Area_Weighted_Normal.h>
#include <Jeffrey_Utilities/RESULT_OF.h>

namespace PhysBAM
{

enum SUB_CUBE_POLYTOPE_ALIGNMENT
{
    SUB_CUBE_POLYTOPE_ALIGNMENT_NULL = 0,
    SUB_CUBE_POLYTOPE_ALIGNMENT_BOUNDARY_ALIGNED,
    SUB_CUBE_POLYTOPE_ALIGNMENT_CUBE_ALIGNED,
    SUB_CUBE_POLYTOPE_ALIGNMENT_SIMPLEX_ALIGNED,
    SUB_CUBE_POLYTOPE_ALIGNMENT_UNALIGNED
};

template< class T, int D >
struct SUB_CUBE_POLYTOPE;

template< class T >
struct SUB_CUBE_POLYTOPE<T,1>
{
    typedef T SCALAR_TYPE;
    static const int DIMENSION = 1;

    typedef SUB_CUBE_POLYTOPE_ALIGNMENT ALIGNMENT_TYPE;

    ALIGNMENT_TYPE alignment;
    T normal; // +1 or -1
    int vertex;

    SUB_CUBE_POLYTOPE();
    SUB_CUBE_POLYTOPE(
        const ALIGNMENT_TYPE alignment_,
        const T normal_);
    SUB_CUBE_POLYTOPE(
        const ALIGNMENT_TYPE alignment_,
        const T normal_,
        const int vertex_);

    template< class T_VERTICES_X, class T_F >
    typename RESULT_OF< T_F ( VECTOR<T,1> ) >::type
    Integrate(const T_VERTICES_X& vertices_x, T_F f) const;
    template< class T_QUADRATURE_RULE, class T_VERTICES_X, class T_F >
    typename RESULT_OF< T_F ( VECTOR<T,1> ) >::type
    Integrate(const T_VERTICES_X& vertices_x, const T_F& f) const;
    template< class T_VERTICES_X, class T_F >
    T Integrate_Flux(const T_VERTICES_X& vertices_x, T_F f) const;
    template< class T_QUADRATURE_RULE, class T_VERTICES_X, class T_F >
    T Integrate_Flux(const T_VERTICES_X& vertices_x, const T_F& f) const;
};

template< class T >
struct SUB_CUBE_POLYTOPE<T,2>
{
    typedef T SCALAR_TYPE;
    static const int DIMENSION = 2;

    typedef SUB_CUBE_POLYTOPE_ALIGNMENT ALIGNMENT_TYPE;

    ALIGNMENT_TYPE alignment;
    VECTOR<T,2> normal;
    VECTOR<int,2> vertices;

    SUB_CUBE_POLYTOPE();
    SUB_CUBE_POLYTOPE(
        const ALIGNMENT_TYPE alignment_,
        const VECTOR<T,2>& normal_);
    SUB_CUBE_POLYTOPE(
        const ALIGNMENT_TYPE alignment_,
        const VECTOR<T,2>& normal_,
        const int vertex1, const int vertex2);

    template< class T_QUADRATURE_RULE, class T_VERTICES_X, class T_F >
    typename RESULT_OF< T_F ( VECTOR<T,2> ) >::type
    Integrate(const T_VERTICES_X& vertices_x, T_F f) const;
    template< class T_QUADRATURE_RULE, class T_VERTICES_X, class T_F >
    T Integrate_Flux(const T_VERTICES_X& vertices_x, T_F f) const;

    template< class T_VERTICES_X >
    void Add_To(
        SEGMENTED_CURVE_2D<T>& curve,
        const VECTOR<T,2> x0,
        const T_VERTICES_X& vertices_x) const;
};

template< class T >
struct SUB_CUBE_POLYTOPE<T,3>
{
    typedef T SCALAR_TYPE;
    static const int DIMENSION = 3;

    typedef SUB_CUBE_POLYTOPE_ALIGNMENT ALIGNMENT_TYPE;

    ALIGNMENT_TYPE alignment;
    VECTOR<T,3> normal;
    BOUNDED_LIST<int,4> vertices;

    SUB_CUBE_POLYTOPE();
    SUB_CUBE_POLYTOPE(
        const ALIGNMENT_TYPE alignment_,
        const VECTOR<T,3>& normal_);
    SUB_CUBE_POLYTOPE(
        const ALIGNMENT_TYPE alignment_,
        const VECTOR<T,3>& normal_,
        const int vertex1, const int vertex2, const int vertex3);
    SUB_CUBE_POLYTOPE(
        const ALIGNMENT_TYPE alignment_,
        const VECTOR<T,3>& normal_,
        const int vertex1, const int vertex2, const int vertex3, const int vertex4);

    template< class T_QUADRATURE_RULE, class T_VERTICES_X, class T_F >
    typename RESULT_OF< T_F ( VECTOR<T,3> ) >::type
    Integrate(const T_VERTICES_X& vertices_x, T_F f) const;
    template< class T_QUADRATURE_RULE, class T_VERTICES_X, class T_F >
    T Integrate_Flux(const T_VERTICES_X& vertices_x, T_F f) const;

    template< class T_VERTICES_X >
    void Add_To(
        TRIANGULATED_SURFACE<T>& surface,
        const VECTOR<T,3> x0,
        const T_VERTICES_X& vertices_x) const;
};

//#####################################################################
//#####################################################################

template< class T >
inline
SUB_CUBE_POLYTOPE<T,1>::
SUB_CUBE_POLYTOPE()
    : alignment(SUB_CUBE_POLYTOPE_ALIGNMENT_NULL)
{ }

template< class T >
inline
SUB_CUBE_POLYTOPE<T,1>::
SUB_CUBE_POLYTOPE(
    const ALIGNMENT_TYPE alignment_,
    const T normal_)
    : alignment(alignment_),
      normal(normal_)
{ }

template< class T >
inline
SUB_CUBE_POLYTOPE<T,1>::
SUB_CUBE_POLYTOPE(
    const ALIGNMENT_TYPE alignment_,
    const T normal_,
    const int vertex_)
    : alignment(alignment_),
      normal(normal_),
      vertex(vertex_)
{ }

template< class T >
template< class T_VERTICES_X, class T_F >
inline typename RESULT_OF< T_F ( VECTOR<T,1> ) >::type
SUB_CUBE_POLYTOPE<T,1>::
Integrate(const T_VERTICES_X& vertices_x, T_F f) const
{ return f(vertices_x(vertex)); }

template< class T >
template< class T_QUADRATURE_RULE, class T_VERTICES_X, class T_F >
inline typename RESULT_OF< T_F ( VECTOR<T,1> ) >::type
SUB_CUBE_POLYTOPE<T,1>::
Integrate(const T_VERTICES_X& vertices_x, const T_F& f) const
{ return Integrate(f, vertices_x); }

template< class T >
template< class T_VERTICES_X, class T_F >
inline T
SUB_CUBE_POLYTOPE<T,1>::
Integrate_Flux(const T_VERTICES_X& vertices_x, T_F f) const
{ return f(vertices_x(vertex), normal); }

template< class T >
template< class T_QUADRATURE_RULE, class T_VERTICES_X, class T_F >
inline T
SUB_CUBE_POLYTOPE<T,1>::
Integrate_Flux(const T_VERTICES_X& vertices_x, const T_F& f) const
{ return Integrate_Flux(f, vertices_x); }

//#####################################################################
//#####################################################################

template< class T >
inline
SUB_CUBE_POLYTOPE<T,2>::
SUB_CUBE_POLYTOPE()
    : alignment(SUB_CUBE_POLYTOPE_ALIGNMENT_NULL)
{ }

template< class T >
inline
SUB_CUBE_POLYTOPE<T,2>::
SUB_CUBE_POLYTOPE(
    const ALIGNMENT_TYPE alignment_,
    const VECTOR<T,2>& normal_)
    : alignment(alignment_),
      normal(normal_)
{ }

template< class T >
inline
SUB_CUBE_POLYTOPE<T,2>::
SUB_CUBE_POLYTOPE(
    const ALIGNMENT_TYPE alignment_,
    const VECTOR<T,2>& normal_,
    const int vertex1, const int vertex2)
    : alignment(alignment_),
      normal(normal_),
      vertices(vertex1, vertex2)
{ }

template< class T >
template< class T_QUADRATURE_RULE, class T_VERTICES_X, class T_F >
inline typename RESULT_OF< T_F ( VECTOR<T,2> ) >::type
SUB_CUBE_POLYTOPE<T,2>::
Integrate(const T_VERTICES_X& vertices_x, T_F f) const
{
    typedef typename RESULT_OF< T_F ( VECTOR<T,2> ) >::type RESULT_TYPE;
    const VECTOR<T,2> x1 = vertices_x(vertices[1]);
    const VECTOR<T,2> x2 = vertices_x(vertices[2]);
    const T length = VECTOR<T,2>::Dot_Product(normal, Area_Weighted_Normal(x1, x2));
    assert(length >= 0);
    RESULT_TYPE unscaled_result = RESULT_TYPE();
    for(int i = 0; i != T_QUADRATURE_RULE::n_weight_x; ++i)
        unscaled_result += T_QUADRATURE_RULE::weight_x[i].weight *
            f(T_QUADRATURE_RULE::weight_x[i].x[0] * x1 +
              T_QUADRATURE_RULE::weight_x[i].x[1] * x2);
    return length * unscaled_result;
}

template< class T >
template< class T_QUADRATURE_RULE, class T_VERTICES_X, class T_F >
inline T
SUB_CUBE_POLYTOPE<T,2>::
Integrate_Flux(const T_VERTICES_X& vertices_x, T_F f) const
{
    const VECTOR<T,2> x1 = vertices_x(vertices[1]);
    const VECTOR<T,2> x2 = vertices_x(vertices[2]);
    const T length = VECTOR<T,2>::Dot_Product(normal, Area_Weighted_Normal(x1, x2));
    assert(length >= 0);
    T unscaled_result = 0;
    for(int i = 0; i != T_QUADRATURE_RULE::n_weight_x; ++i)
        unscaled_result += T_QUADRATURE_RULE::weight_x[i].weight *
            f(T_QUADRATURE_RULE::weight_x[i].x[0] * x1 +
              T_QUADRATURE_RULE::weight_x[i].x[1] * x2, normal);
    return length * unscaled_result;
}

template< class T >
template< class T_VERTICES_X >
inline void
SUB_CUBE_POLYTOPE<T,2>::
Add_To(
    SEGMENTED_CURVE_2D<T>& curve,
    const VECTOR<T,2> x0,
    const T_VERTICES_X& vertices_x) const
{
    const int particle_offset = curve.particles.array_collection->Size();
    curve.particles.array_collection->Add_Elements(2);
    curve.particles.X(particle_offset + 1) = x0 + vertices_x(vertices[1]);
    curve.particles.X(particle_offset + 2) = x0 + vertices_x(vertices[2]);
    curve.mesh.elements.Append(particle_offset + VECTOR<int,2>(1,2));
}

//#####################################################################
//#####################################################################

template< class T >
inline
SUB_CUBE_POLYTOPE<T,3>::
SUB_CUBE_POLYTOPE()
    : alignment(SUB_CUBE_POLYTOPE_ALIGNMENT_NULL)
{ }

template< class T >
inline
SUB_CUBE_POLYTOPE<T,3>::
SUB_CUBE_POLYTOPE(
    const ALIGNMENT_TYPE alignment_,
    const VECTOR<T,3>& normal_)
    : alignment(alignment_),
      normal(normal_)
{ }

template< class T >
inline
SUB_CUBE_POLYTOPE<T,3>::
SUB_CUBE_POLYTOPE(
    const ALIGNMENT_TYPE alignment_,
    const VECTOR<T,3>& normal_,
    const int vertex1, const int vertex2, const int vertex3)
    : alignment(alignment_),
      normal(normal_)
{
    vertices.Append(vertex1);
    vertices.Append(vertex2);
    vertices.Append(vertex3);
}

template< class T >
inline
SUB_CUBE_POLYTOPE<T,3>::
SUB_CUBE_POLYTOPE(
    const ALIGNMENT_TYPE alignment_,
    const VECTOR<T,3>& normal_,
    const int vertex1, const int vertex2, const int vertex3, const int vertex4)
    : alignment(alignment_),
      normal(normal_)
{
    vertices.Append(vertex1);
    vertices.Append(vertex2);
    vertices.Append(vertex3);
    vertices.Append(vertex4);
}

template< class T >
template< class T_QUADRATURE_RULE, class T_VERTICES_X, class T_F >
inline typename RESULT_OF< T_F ( VECTOR<T,3> ) >::type
SUB_CUBE_POLYTOPE<T,3>::
Integrate(const T_VERTICES_X& vertices_x, T_F f) const
{
    typedef typename RESULT_OF< T_F ( VECTOR<T,3> ) >::type RESULT_TYPE;
    assert(vertices.Size() >= 3);
    RESULT_TYPE result = RESULT_TYPE();
    const VECTOR<T,3> x1 = vertices_x(vertices(1));
    VECTOR<T,3> xv = vertices_x(vertices(2));
    for(int v = 3; v <= vertices.Size(); ++v){
        const VECTOR<T,3> xvm1 = xv;
        xv = vertices_x(vertices(v));
        const T partial_area = VECTOR<T,3>::Dot_Product(normal, Area_Weighted_Normal(x1, xvm1, xv)) / 2;
        assert(partial_area >= 0);
        RESULT_TYPE unscaled_partial_result = RESULT_TYPE();
        for(int i = 0; i != T_QUADRATURE_RULE::n_weight_x; ++i)
            unscaled_partial_result += T_QUADRATURE_RULE::weight_x[i].weight *
                f(T_QUADRATURE_RULE::weight_x[i].x[0] * x1 +
                  T_QUADRATURE_RULE::weight_x[i].x[1] * xvm1 +
                  T_QUADRATURE_RULE::weight_x[i].x[2] * xv);
        result += partial_area * unscaled_partial_result;
    }
    return result;
}

template< class T >
template< class T_QUADRATURE_RULE, class T_VERTICES_X, class T_F >
inline T
SUB_CUBE_POLYTOPE<T,3>::
Integrate_Flux(const T_VERTICES_X& vertices_x, T_F f) const
{
    T result = 0;
    const VECTOR<T,3> x1 = vertices_x(vertices(1));
    VECTOR<T,3> xv = vertices_x(vertices(2));
    for(int v = 3; v <= vertices.Size(); ++v){
        const VECTOR<T,3> xvm1 = xv;
        xv = vertices_x(vertices(v));
        const T partial_area = VECTOR<T,3>::Dot_Product(normal, Area_Weighted_Normal(x1, xvm1, xv)) / 2;
        assert(partial_area >= 0);
        T unscaled_partial_result = 0;
        for(int i = 0; i != T_QUADRATURE_RULE::n_weight_x; ++i)
            unscaled_partial_result += T_QUADRATURE_RULE::weight_x[i].weight *
                f(T_QUADRATURE_RULE::weight_x[i].x[0] * x1 +
                  T_QUADRATURE_RULE::weight_x[i].x[1] * xvm1 +
                  T_QUADRATURE_RULE::weight_x[i].x[2] * xv,
                  normal);
        result += partial_area * unscaled_partial_result;
    }
    return result;
}

template< class T >
template< class T_VERTICES_X >
inline void
SUB_CUBE_POLYTOPE<T,3>::
Add_To(
    TRIANGULATED_SURFACE<T>& surface,
    const VECTOR<T,3> x0,
    const T_VERTICES_X& vertices_x) const
{
    const int particle_offset = surface.particles.array_collection->Size();
    surface.particles.array_collection->Add_Elements(vertices.Size());
    for(int v = 1; v <= vertices.Size(); ++v)
        surface.particles.X(particle_offset + v) = x0 + vertices_x(vertices(v));
    for(int v = 3; v <= vertices.Size(); ++v)
        surface.mesh.elements.Append(particle_offset + VECTOR<int,3>(1,v-1,v));
}

} // namespace PhysBAM

#endif // #ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_GEOMETRY_SUB_CUBE_POLYTOPE_HPP
