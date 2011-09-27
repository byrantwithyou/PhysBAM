//#####################################################################
// Copyright 2011, Jeffrey Hellrung
// This file is part of PhysBAM whose distribution is governed by the
// license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################

#ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_STENCILS_CROSS_CONSTBETA_STENCIL_HPP
#define PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_STENCILS_CROSS_CONSTBETA_STENCIL_HPP

#include <cassert>

#include <algorithm>
#include <numeric>

#include <boost/concept/assert.hpp>
#include <boost/foreach.hpp>
#include <boost/mpl/assert.hpp>
#include <boost/type_traits/is_same.hpp>

#include <PhysBAM_Tools/Vectors/VECTOR.h>
#include <Jeffrey_Utilities/LOOP.h>
#include <Jeffrey_Utilities/Stencils/STENCIL_PROXY_CONCEPT.h>

namespace PhysBAM
{

template< int D >
struct CROSS_CONSTBETA_STENCIL
{
    typedef VECTOR<int,D> MULTI_INDEX_TYPE;

    typedef unsigned char VALUE_TYPE;
    static const int DIMENSION = D;
    static const int N_VALUES = 1 + 2*D;

    unsigned char values[N_VALUES];

    static CROSS_CONSTBETA_STENCIL Construct_Zero();

    void Zero();
    int N_Zero() const;
    int N_Nonzero() const;

    void Add(const MULTI_INDEX_TYPE cell_multi_offset);

    template< class T >
    void Assign(
        const VECTOR<T,D> beta_dv_over_dx_dx,
        const MULTI_INDEX_TYPE multi_offset,
        const T value);
    template< class T >
    void Add(
        const VECTOR<T,D> beta_dv_over_dx_dx,
        const MULTI_INDEX_TYPE multi_offset,
        const T value);

    template< class T, class T_STENCIL_PROXY >
    void Assign(
        const VECTOR<T,D>& beta_dv_over_dx_dx,
        const MULTI_INDEX_TYPE base_multi_index,
        const T_STENCIL_PROXY stencil_proxy);
    template< class T, class T_STENCIL_PROXY >
    void Add(
        const VECTOR<T,D>& beta_dv_over_dx_dx,
        const MULTI_INDEX_TYPE base_multi_index,
        const T_STENCIL_PROXY stencil_proxy);

    unsigned char& Center();
    const unsigned char& Center() const;
    template< class T >
    T Center(const VECTOR<T,D>& beta_dv_over_dx_dx) const;
    unsigned char& operator()(const MULTI_INDEX_TYPE multi_offset);
    const unsigned char& operator()(const MULTI_INDEX_TYPE multi_offset) const;
    unsigned char& operator()(const int d, const int s);
    const unsigned char& operator()(const int d, const int s) const;
    template< class T >
    T operator()(
        const MULTI_INDEX_TYPE multi_offset,
        const VECTOR<T,D>& beta_dv_over_dx_dx) const;
    template< class T >
    T operator()(
        const int d, const int s,
        const VECTOR<T,D>& beta_dv_over_dx_dx) const;

    unsigned char Sum() const;
    template< class T >
    T Sum(const VECTOR<T,D> beta_dv_over_dx_dx) const;

    template< class T, class T_VECTOR >
    T Apply(
        const VECTOR<T,D> beta_dv_over_dx_dx,
        const int base_index, const VECTOR<int,D> strides, const T_VECTOR& x) const;
    template< class T, class T_VECTOR >
    void Apply_Transpose(
        const VECTOR<T,D> beta_dv_over_dx_dx,
        const int base_index, const VECTOR<int,D> strides, T_VECTOR& y, const T x) const;
};

//#####################################################################
//#####################################################################

template< int D >
inline 
CROSS_CONSTBETA_STENCIL<D>
CROSS_CONSTBETA_STENCIL<D>::
Construct_Zero()
{
    CROSS_CONSTBETA_STENCIL result;
    result.Zero();
    return result;
}

template< int D >
inline void
CROSS_CONSTBETA_STENCIL<D>::
Zero()
{ std::fill(&values[0], &values[N_VALUES], static_cast< unsigned char >(0)); }

template< int D >
inline int
CROSS_CONSTBETA_STENCIL<D>::
N_Zero() const
{ return static_cast< int >(std::count(&values[0], &values[N_VALUES], static_cast< unsigned char >(0))); }

template< int D >
inline int
CROSS_CONSTBETA_STENCIL<D>::
N_Nonzero() const
{ return N_VALUES - N_Zero(); }

template< int D >
inline void
CROSS_CONSTBETA_STENCIL<D>::
Add(const MULTI_INDEX_TYPE cell_multi_offset)
{
    ++Center();
    for(int d = 1; d <= D; ++d) {
        assert(cell_multi_offset[d] == -1 || cell_multi_offset[d] == 0);
        const int s = 1 + 2 * cell_multi_offset[d];
        ++operator()(d,s);
    }
}

template< int D >
template< class T >
inline void
CROSS_CONSTBETA_STENCIL<D>::
Assign(
    const VECTOR<T,D> beta_dv_over_dx_dx,
    const MULTI_INDEX_TYPE multi_offset,
    const T value)
{
    assert(multi_offset.L1_Norm() <= 1);
    for(int d = 1; d != D; ++d) {
        if(const int s = multi_offset[d]) {
            assert(value <= 0);
            values[D + s * (D - d + 1)] = static_cast< unsigned char >(
                -value / beta_dv_over_dx_dx[d] + static_cast<T>(0.5)
            );
            return;
        }
    }
    assert(value >= 0);
    Center() = static_cast< unsigned char >(
        value / beta_dv_over_dx_dx.Sum() + static_cast<T>(0.5)
    );
}

template< int D >
template< class T >
inline void
CROSS_CONSTBETA_STENCIL<D>::
Add(
    const VECTOR<T,D> beta_dv_over_dx_dx,
    const MULTI_INDEX_TYPE multi_offset,
    const T value)
{
    assert(multi_offset.L1_Norm() <= 1);
    for(int d = 1; d != D; ++d) {
        if(const int s = multi_offset[d]) {
            assert(value <= 0);
            values[D + s * (D - d + 1)] += static_cast< unsigned char >(
                -value / beta_dv_over_dx_dx[d] + static_cast<T>(0.5)
            );
            return;
        }
    }
    assert(value >= 0);
    Center() += static_cast< unsigned char >(
        value / beta_dv_over_dx_dx.Sum() + static_cast<T>(0.5)
    );
}

template< int D >
template< class T, class T_STENCIL_PROXY >
inline void
CROSS_CONSTBETA_STENCIL<D>::
Assign(
    const VECTOR<T,D>& beta_dv_over_dx_dx,
    const MULTI_INDEX_TYPE base_multi_index,
    const T_STENCIL_PROXY stencil_proxy)
{
    BOOST_CONCEPT_ASSERT((CONST_STENCIL_PROXY_CONCEPT< T_STENCIL_PROXY >));
    BOOST_MPL_ASSERT((boost::is_same< typename T_STENCIL_PROXY::INDEX_TYPE, MULTI_INDEX_TYPE >));
    BOOST_MPL_ASSERT((boost::is_same< typename T_STENCIL_PROXY::SCALAR_TYPE, T >));
    Zero();
    BOOST_FOREACH( typename T_STENCIL_PROXY::reference index_value, stencil_proxy )
        if(index_value.value != 0)
            Assign(beta_dv_over_dx_dx, index_value.index - base_multi_index, index_value.value);
}

template< int D >
template< class T, class T_STENCIL_PROXY >
inline void
CROSS_CONSTBETA_STENCIL<D>::
Add(
    const VECTOR<T,D>& beta_dv_over_dx_dx,
    const MULTI_INDEX_TYPE base_multi_index,
    const T_STENCIL_PROXY stencil_proxy)
{
    BOOST_CONCEPT_ASSERT((CONST_STENCIL_PROXY_CONCEPT< T_STENCIL_PROXY >));
    BOOST_MPL_ASSERT((boost::is_same< typename T_STENCIL_PROXY::INDEX_TYPE, MULTI_INDEX_TYPE >));
    BOOST_MPL_ASSERT((boost::is_same< typename T_STENCIL_PROXY::SCALAR_TYPE, T >));
    BOOST_FOREACH( typename T_STENCIL_PROXY::reference index_value, stencil_proxy )
        if(index_value.value != 0)
            Add(beta_dv_over_dx_dx, index_value.index - base_multi_index, index_value.value);
}

template< int D >
inline unsigned char&
CROSS_CONSTBETA_STENCIL<D>::
Center()
{ return values[D]; }

template< int D >
inline const unsigned char&
CROSS_CONSTBETA_STENCIL<D>::
Center() const
{ return values[D]; }

template< int D >
template< class T >
inline T
CROSS_CONSTBETA_STENCIL<D>::
Center(const VECTOR<T,D>& beta_dv_over_dx_dx) const
{ return Center() * beta_dv_over_dx_dx.Sum(); }

template< int D >
inline unsigned char&
CROSS_CONSTBETA_STENCIL<D>::
operator()(const MULTI_INDEX_TYPE multi_offset)
{
    assert(multi_offset.L1_Norm() <= 1);
    for(int d = 1; d != D; ++d)
        if(const int s = multi_offset[d])
            return operator()(d,s);
    return Center();
}

template< int D >
inline const unsigned char&
CROSS_CONSTBETA_STENCIL<D>::
operator()(const MULTI_INDEX_TYPE multi_offset) const
{
    assert(multi_offset.L1_Norm() <= 1);
    for(int d = 1; d != D; ++d)
        if(const int s = multi_offset[d])
            return operator()(d,s);
    return Center();
}

template< int D >
inline unsigned char&
CROSS_CONSTBETA_STENCIL<D>::
operator()(const int d, const int s)
{
    assert(1 <= d && d <= D);
    assert(s == -1 || s == +1);
    return values[D + s * (D - d + 1)];
}

template< int D >
inline const unsigned char&
CROSS_CONSTBETA_STENCIL<D>::
operator()(const int d, const int s) const
{
    assert(1 <= d && d <= D);
    assert(s == -1 || s == +1);
    return values[D + s * (D - d + 1)];
}

template< int D >
template< class T >
inline T
CROSS_CONSTBETA_STENCIL<D>::
operator()(
    const MULTI_INDEX_TYPE multi_offset,
    const VECTOR<T,D>& beta_dv_over_dx_dx) const
{
    assert(multi_offset.L1_Norm() <= 1);
    for(int d = 1; d != D; ++d)
        if(const int s = multi_offset[d])
            return operator()(d, s, beta_dv_over_dx_dx);
    return Center(beta_dv_over_dx_dx);
}

template< int D >
template< class T >
inline T
CROSS_CONSTBETA_STENCIL<D>::
operator()(
    const int d, const int s,
    const VECTOR<T,D>& beta_dv_over_dx_dx) const
{ return operator()(d,s) * (-beta_dv_over_dx_dx[d]); }

template< int D >
inline unsigned char
CROSS_CONSTBETA_STENCIL<D>::
Sum() const
{ return 4 * Center() - std::accumulate(&values[0], &values[N_VALUES], static_cast< unsigned char >(0)); }

template< int D >
template< class T >
inline T
CROSS_CONSTBETA_STENCIL<D>::
Sum(const VECTOR<T,D> beta_dv_over_dx_dx) const
{
    T sum = Center() * beta_dv_over_dx_dx.Sum();
    for(int d = 1; d <= D; ++d)
        sum -= (values[d - 1] + values[2*D - d + 1]) * beta_dv_over_dx_dx[d];
    return sum;
}

template< int D >
template< class T, class T_VECTOR >
inline T
CROSS_CONSTBETA_STENCIL<D>::
Apply(
    const VECTOR<T,D> beta_dv_over_dx_dx,
    const int base_index, const VECTOR<int,D> strides, const T_VECTOR& x) const
{
    T y = Center() * beta_dv_over_dx_dx.Sum() * x(base_index);
    for(int d = 1; d <= D; ++d) {
        if(values[d - 1] != 0) {
            assert(1 <= base_index - strides[d]);
            y -= values[d - 1] * beta_dv_over_dx_dx[d] * x(base_index - strides[d]);
        }
        if(values[2*D - d + 1] != 0) {
            assert(base_index + strides[d] <= x.Size());
            y -= values[2*D - d + 1] * beta_dv_over_dx_dx[d] * x(base_index + strides[d]);
        }
    }
    return y;
}

template< int D >
template< class T, class T_VECTOR >
inline void
CROSS_CONSTBETA_STENCIL<D>::
Apply_Transpose(
    const VECTOR<T,D> beta_dv_over_dx_dx,
    const int base_index, const VECTOR<int,D> strides, T_VECTOR& y, const T x) const
{
    y(base_index) += Center() * beta_dv_over_dx_dx.Sum() * x;
    for(int d = 1; d <= D; ++d) {
        if(values[d - 1] != 0) {
            assert(1 <= base_index - strides[d]);
            y(base_index - strides[d]) -= values[d - 1] * beta_dv_over_dx_dx[d] * x;
        }
        if(values[2*D - d + 1] != 0) {
            assert(base_index + strides[d] <= y.Size());
            y(base_index + strides[d]) -= values[2*D - d + 1] * beta_dv_over_dx_dx[d] * x;
        }
    }
}

} // namespace PhysBAM

#endif // #ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_STENCILS_CROSS_CONSTBETA_STENCIL_HPP
